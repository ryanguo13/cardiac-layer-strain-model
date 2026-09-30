# ------------------------------------------------------------------
#  Beat-cycle animation — short-axis contraction + twist
#  Renders ./out/anim_beat.mp4 (or a series of PNGs if FFmpeg is missing)
# ------------------------------------------------------------------

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using CairoMakie
using GeometryBasics
using Colors
using Printf

include(joinpath(@__DIR__, "..", "src", "ThreeLayerModel.jl"))
using .ThreeLayerModel

const OUTDIR = joinpath(@__DIR__, "..", "out")
isdir(OUTDIR) || mkpath(OUTDIR)

const ANIMDIR = joinpath(OUTDIR, "anim_frames")
isdir(ANIMDIR) || mkpath(ANIMDIR)

# ------------------------------------------------------------------
#  helpers
# ------------------------------------------------------------------

function short_axis_strain_frame(t::Float64, h::HeartState;
                                 age::Float64 = h.age)
    s = beat_strain!(t, h)
    return s
end

function render_frame!(scene_obs, t::Float64, h::HeartState)
    s = short_axis_strain_frame(t, h)
    # s.short ∈ [-A_max, 0], s.radial ∈ [0, +A_max], s.twist ∈ rad
    short = clamp(-s.short / h.params.A_max, 0.0, 1.0)   # 0..1 progress
    thick = clamp( s.radial / h.params.A_max, 0.0, 1.0)
    twist_deg = rad2deg(s.twist)

    # geometry contracts radially and thickens the wall
    base_radii = (3.2, 2.7, 2.2, 1.5)
    r = base_radii .* (1.0 - 0.10 * short)
    wall_growth = (1.0 + 0.45 * thick)
    r = (r[1],
         r[2] + 0.05 * (wall_growth - 1.0),
         r[3] - 0.05 * (wall_growth - 1.0),
         r[4] - 0.10 * (wall_growth - 1.0))

    ax = scene_obs.ax
    empty!(ax)

    θ = range(0, 2π; length = 600)
    # layer fills (annular regions between consecutive radii)
    band_colors = (RGBA(0.20,0.70,0.30,0.30),  # R green (outer, epi)
                   RGBA(0.20,0.45,0.90,0.30),  # O blue  (mid, mid-wall)
                   RGBA(0.85,0.20,0.20,0.30))  # L red   (inner, sub-endo)
    # r = (epi-outer, epi-inner=mid-outer, mid-inner=endo-outer, endo-inner)
    # annulus k: r[k+1] (inner) .. r[k] (outer), for k=1..3
    for k in 1:3
        r_out = r[k]
        r_in  = r[k+1]
        xs = vcat(r_out .* cos.(θ), reverse(r_in .* cos.(θ)))
        ys = vcat(r_out .* sin.(θ), reverse(r_in .* sin.(θ)))
        poly!(ax, xs, ys; color = band_colors[k])
    end
    # outer/inner outlines
    for rr in r
        lines!(ax, rr .* cos.(θ), rr .* sin.(θ);
               color = (:black), linewidth = 0.8, linestyle = (:dash))
    end

    # helix glyphs rotated by twist
    cmap = (RGBA(0.85,0.20,0.20,0.85),   # k=1 L red (innermost, sub-endo)
            RGBA(0.20,0.45,0.90,0.85),   # k=2 O blue (mid, mid-wall)
            RGBA(0.20,0.70,0.30,0.85))   # k=3 R green (outermost, epi)
    helix_deg = (+60.0, 0.0, -60.0)
    glyph_radii = ((r[3]+r[4])/2, (r[2]+r[3])/2, (r[1]+r[2])/2)
    for k in 1:3
        for φ in range(0, 2π; length = 24)[1:end-1]
            φt = φ + deg2rad(twist_deg) * (k == 1 ? +1 : k == 3 ? -1 : 0)
            x, y = glyph_radii[k] * cos(φt), glyph_radii[k] * sin(φt)
            β = deg2rad(helix_deg[k])
            vx, vy = -sin(φt + β), cos(φt + β)
            L = 0.55
            lines!(ax, [x - vx*L/2, x + vx*L/2], [y - vy*L/2, y + vy*L/2];
                   color = cmap[k], linewidth = 2.0)
            scatter!(ax, [x + vx*L/2], [y + vy*L/2]; color = cmap[k], markersize = 4)
        end
    end

    # GLS read-out
    ε_pct = s.short * 100
    text!(ax, 0, 3.95;
          text = @sprintf("t = %4.2f s    GLS = %+5.2f %%    twist = %+5.2f°",
                          t, ε_pct, twist_deg),
          color = (:black), align = (:center, :center), fontsize = 14)

    return s
end

# ------------------------------------------------------------------
#  Build the animation
# ------------------------------------------------------------------

function make_animation(; ages = (25, 80),
                         frames_per_beat::Int = 120)
    out_mp4 = joinpath(OUTDIR, "anim_beat.mp4")

    fig = Figure(size = (900, 900), fontsize = 16)
    ax  = Axis(fig[1, 1]; aspect = DataAspect(),
               title = "LV short-axis beat — three-layer myofiber model")
    xlims!(ax, -4.4, 4.4); ylims!(ax, -4.4, 4.4)
    hidedecorations!(ax, grid = false, ticks = false)
    scene_obs = (fig = fig, ax = ax)

    h = initialize_heart(age = Float64(ages[1]))

    τ = h.params.τ_beat
    times = range(0, τ; length = frames_per_beat)

    println("Rendering $(frames_per_beat) frames …")
    @time begin
        record(fig, out_mp4, times; framerate = 30) do t
            render_frame!(scene_obs, t, h)
        end
    end
    println("Saved animation: $out_mp4")
    return out_mp4
end

# ------------------------------------------------------------------
#  Optional: also produce a side-by-side young-vs-old comparison
# ------------------------------------------------------------------

function make_comparison_animation(; ages = (25, 80),
                                    frames_per_beat::Int = 90)
    out_mp4 = joinpath(OUTDIR, "anim_compare_young_old.mp4")

    fig = Figure(size = (1700, 850), fontsize = 16)
    ax_young = Axis(fig[1, 1]; aspect = DataAspect(),
                    title = "Young (age $(ages[1]))")
    ax_old   = Axis(fig[2, 1]; aspect = DataAspect(),
                    title = "Aged  (age $(ages[2]))")
    for ax in (ax_young, ax_old)
        xlims!(ax, -4.4, 4.4); ylims!(ax, -4.4, 4.4)
        hidedecorations!(ax, grid = false, ticks = false)
    end
    obs = (young = ax_young, old = ax_old)

    h_young = initialize_heart(age = Float64(ages[1]))
    h_old   = initialize_heart(age = Float64(ages[2]))

    τ = h_young.params.τ_beat
    times = range(0, τ; length = frames_per_beat)

    println("Rendering $(frames_per_beat) comparison frames …")
    @time begin
        record(fig, out_mp4, times; framerate = 30) do t
            empty!(ax_young); empty!(ax_old)
            render_frame!((ax = ax_young,), t, h_young)
            render_frame!((ax = ax_old,),   t, h_old)
        end
    end
    println("Saved comparison animation: $out_mp4")
    return out_mp4
end

# ------------------------------------------------------------------
#  Driver
# ------------------------------------------------------------------

function main()
    p1 = make_animation()
    p2 = make_comparison_animation()
    return (p1, p2)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
