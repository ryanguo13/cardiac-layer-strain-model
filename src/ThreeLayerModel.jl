# ------------------------------------------------------------------
#  LV three-layer myofiber mechanics
#  ------------------------------------------------------------------
#  Models the left ventricle (LV) as three concentric transmural layers
#  of helical myofibers whose orientations rotate smoothly from a
#  right-handed helix at the endocardium to a left-handed helix at the
#  epicardium (Streeter & Bassett classic transmural rotation).
#
#  Layer 1 (inner, sub-endocardium)  → longitudinal, helix angle ≈ +60°
#  Layer 2 (mid, mesocardium)         → circumferential / oblique ≈  0°
#  Layer 3 (outer, sub-epicardium)   → radial-oblique, helix angle ≈ –60°
#
#  Each layer contracts simultaneously and produces three components of
#  deformation (narrowing, shortening, twisting).
#
#  We track the *health* (0..1) of each layer across ageing.  Healthy
#  tissue generates a contractile strain amplitude Amax.  During ageing
#  damage accumulates and surviving layers can compensate for failing
#  ones.  The compensation model follows the model brief:
#
#      longitudinal (L) layer fails first  →
#      radial (R) layer compensates        →
#      radial layer eventually de-compensates.
#
#  From the instantaneous strain tensors of the three layers we derive
#  the clinically important **Global Longitudinal Strain (GLS)** and
#  Global Radial Strain (GRS).
# ------------------------------------------------------------------

module ThreeLayerModel

using LinearAlgebra
using Statistics

export LayerState, HeartState, MyofiberParams, AgingParams,
       initialize_heart, step_age!, beat_strain!, gls_trace,
       layer_strain_contributions, helix_angle, transmural_rotation

# ------------------------------------------------------------------
#  Geometry / constitutive parameters
# ------------------------------------------------------------------

"""
    MyofiberParams

Bundle-level myofiber parameters.
- `A_max`         maximum strain amplitude (healthy, fully recruited)
- `τ_beat`        beat duration in seconds
- `τ_twist`       relative phase of twist vs. shortening
- `peak_phase`    phase within the beat at which shortening peaks
"""
Base.@kwdef struct MyofiberParams
    A_max::Float64          = 0.20   # 20 % peak shortening (peak systolic strain)
    τ_beat::Float64         = 0.85   # s, one cardiac cycle
    τ_twist::Float64        = 0.12   # phase shift, s, between shortening and twist
    peak_phase::Float64     = 0.45   # 0..1, fraction of the beat at which peak occurs
end

"""
    AgingParams

Ageing and compensation rules.
- `L_fail_age`     age at which the longitudinal layer starts to decline
- `R_fail_age`     age at which the radial layer also starts to decline
- `L_fail_rate`    per-year rate of damage to L once it has started
- `R_fail_rate`    per-year rate of damage to R once it has started
- `compensation`   0..1, how much a healthy layer can take over for a failing one
"""
Base.@kwdef struct AgingParams
    L_fail_age::Float64     = 35.0
    R_fail_age::Float64     = 60.0
    L_fail_rate::Float64    = 0.018
    R_fail_rate::Float64    = 0.030
    compensation::Float64   = 0.55
end

# ------------------------------------------------------------------
#  Helix-angle field
# ------------------------------------------------------------------

"""
    helix_angle(layer)

Returns the canonical Streeter–Bassett transmural helix angle (deg)
for a given layer index (1 = endocardial, 2 = mid, 3 = epicardial).
"""
helix_angle(layer::Int) = [60.0, 0.0, -60.0][layer]

"""
    transmural_rotation(t, layer)

Fraction of the twist component contributed by `layer` at phase `t`
(0..1).  Inner and outer layers twist in opposite senses, so this is
signed.  Used to combine the three layers into a single net twist.
"""
function transmural_rotation(t::Float64, layer::Int)
    θ = deg2rad(helix_angle(layer))
    return sign(θ) * abs(sin(π * t)) * cos(θ)
end

# ------------------------------------------------------------------
#  Layer / heart state
# ------------------------------------------------------------------

"""
    LayerState

One transmural myofiber layer.
- `health`  ∈ [0,1] remaining contractile capacity
- `helix`   helix angle (deg)
- `phase_offset` small per-layer phase offset (s) for physiological dyssynchrony
"""
Base.@kwdef mutable struct LayerState
    health::Float64 = 1.0
    helix::Float64  = 0.0
    phase_offset::Float64 = 0.0
end

"""
    HeartState

Three layers plus their parameters.
"""
mutable struct HeartState
    L::LayerState   # longitudinal  (inner,   helix ≈ +60°)
    O::LayerState   # oblique       (middle,  helix ≈   0°)
    R::LayerState   # radial-oblique(outer,   helix ≈ -60°)
    params::MyofiberParams
    age::Float64
end

function initialize_heart(; age::Float64 = 25.0,
                          params::MyofiberParams = MyofiberParams(),
                          aging::AgingParams = AgingParams())
    layers = (
        L = LayerState(health = 1.0, helix =  60.0, phase_offset = -0.04),
        O = LayerState(health = 1.0, helix =   0.0, phase_offset =  0.00),
        R = LayerState(health = 1.0, helix = -60.0, phase_offset =  0.04),
    )
    # Run a fast-forward ageing pass to set initial health for the
    # requested age.
    h = HeartState(layers.L, layers.O, layers.R, params, age)
    step_age!(h, aging, target_age = age)
    return h
end

# ------------------------------------------------------------------
#  Ageing dynamics — damage + compensation
# ------------------------------------------------------------------

"""
    step_age!(h, aging; target_age)

Advance the heart to `target_age`.  The longitudinal layer starts to
fail at `aging.L_fail_age`; once it does the outer (radial-oblique)
layer compensates.  The radial layer itself begins to fail at
`aging.R_fail_age`, after which there is no further reserve.

`compensation` is the *fraction* of lost longitudinal health that the
radial layer can take over (and vice-versa when R fails).
"""
function step_age!(h::HeartState, aging::AgingParams;
                   target_age::Float64 = h.age + 1.0)
    Δt = max(target_age - h.age, 0.0)
    h.age = target_age

    # --- intrinsic damage -----------------------------------------------
    L_damage = h.age ≥ aging.L_fail_age ?
               aging.L_fail_rate * (h.age - aging.L_fail_age) : 0.0
    R_damage = h.age ≥ aging.R_fail_age ?
               aging.R_fail_rate * (h.age - aging.R_fail_age) : 0.0

    new_L = clamp(1.0 - L_damage, 0.0, 1.0)
    new_R = clamp(1.0 - R_damage, 0.0, 1.0)
    new_O = 1.0                       # middle layer degrades only weakly

    # --- compensation: L → R then R → L -------------------------------
    # When L drops, R hypertrophies up to a ceiling set by `compensation`.
    if new_L < h.L.health
        deficit = (h.L.health - new_L) * aging.compensation
        new_R = clamp(new_R + deficit, 0.0, 1.0 + 0.5 * aging.compensation)
    end
    # When R also drops, no further compensation is available.
    if new_R < h.R.health
        new_R = clamp(new_R, 0.0, 1.0)
    end

    # gentle mid-wall decline at very old age
    if h.age > 70
        new_O = clamp(new_O - 0.004 * (h.age - 70), 0.4, 1.0)
    end

    h.L.health = new_L
    h.O.health = new_O
    h.R.health = new_R
    return h
end

# ------------------------------------------------------------------
#  Beat kinematics
# ------------------------------------------------------------------

"""
    beat_envelope(t, params)

Symmetric systolic strain envelope normalised to [0,1] within one
cardiac cycle (t in seconds).  Used by both shortening and narrowing.
"""
function beat_envelope(t::Float64, p::MyofiberParams)
    τ   = p.τ_beat
    φ   = mod(t, τ) / τ                      # 0..1
    x   = (φ - p.peak_phase) * π
    return exp(-x^2)                         # Gaussian peak
end

"""
    beat_twist_envelope(t, params)

Slightly phase-shifted twist envelope — physiological torsion
peaks ~30 ms after shortening.
"""
function beat_twist_envelope(t::Float64, p::MyofiberParams)
    τ   = p.τ_beat
    φ   = mod(t + p.τ_twist, τ) / τ
    x   = (φ - p.peak_phase) * π
    return exp(-x^2)
end

"""
    beat_strain!(t, h)

Populate the three scalar strain components for `h` at time `t`:
- `short[t]`  longitudinal shortening (GLS component, negative = shortening)
- `radial[t]` radial thickening (positive = thickening inward)
- `twist[t]`  net twist angle (rad)

Returns a NamedTuple of (short, radial, twist).
"""
function beat_strain!(t::Float64, h::HeartState)
    p = h.params
    env = beat_envelope(t, p)
    tenv = beat_twist_envelope(t, p)

    # each layer's intrinsic strain amplitude is scaled by its health
    A_L = -p.A_max * h.L.health   # negative: shortening
    A_O = -p.A_max * h.O.health
    A_R = -p.A_max * h.R.health

    # micro dyssynchrony
    function layer_short(layer::LayerState)
        env_local = beat_envelope(t - layer.phase_offset, p)
        return layer == h.L ? A_L * env_local :
               layer == h.O ? A_O * env_local :
                              A_R * env_local
    end
    s_L = layer_short(h.L)
    s_O = layer_short(h.O)
    s_R = layer_short(h.R)

    # longitudinal shortening = weighted sum, longitudinal layer dominates
    short = 0.50 * s_L + 0.30 * s_O + 0.20 * s_R

    # radial thickening is primarily the radial layer (helix ≈ ±60°)
    rad_L = 0.20 * (-s_L)
    rad_O = 0.30 * (-s_O)
    rad_R = 0.50 * (-s_R)
    radial = (rad_L + rad_O + rad_R) * env

    # net twist from opposite helicities: subtract endo and epi
    twist_layer = (transmural_rotation(mod(t, p.τ_beat) / p.τ_beat, 1) * h.L.health
                 - transmural_rotation(mod(t, p.τ_beat) / p.τ_beat, 3) * h.R.health) *
                  p.A_max * 1.8 * tenv

    return (short = short, radial = radial, twist = twist_layer)
end

# ------------------------------------------------------------------
#  Global longitudinal strain (GLS) over one beat
# ------------------------------------------------------------------

"""
    gls_trace(h; n = 240)

Sample the global longitudinal strain (most-negative peak) over one
beat at `n` time points.  Returns (t, ε_L) arrays.
"""
function gls_trace(h::HeartState; n::Int = 240)
    τ = h.params.τ_beat
    t = range(0.0, τ; length = n)
    ε = similar(t, Float64)
    for i in eachindex(t)
        ε[i] = beat_strain!(t[i], h).short
    end
    return collect(t), ε
end

"""
    layer_strain_contributions(h)

The per-layer contribution to GLS at peak systole — useful to plot
the compensatory shift from L-dominant to R-dominant strain.
"""
function layer_strain_contributions(h::HeartState)
    return (
        L = -h.params.A_max * h.L.health,
        O = -h.params.A_max * h.O.health,
        R = -h.params.A_max * h.R.health,
    )
end

end # module
