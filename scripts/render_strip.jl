# ------------------------------------------------------------------
#  Render a strip of still PNGs across one beat
#  (lets us inspect the animation visually without playing the MP4)
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

# ---- local copy of the per-frame renderer (same logic as in animate_beat.jl)
function render_local!(ax::Axis, t::Float64, h::HeartState)
    s = beat_strain!(t, h)
    short = clamp(-s.short / h.params.A_max, 0.0, 1.0)
    thick = clamp( s.radial / h.params.A_max, 0.0, 1.0)
    twist_deg = rad2deg(s.twist)

    base_radii = (3.2, 2.7, 2.2, 1.5)
    r = base_radii .* (1.0 - 0.10 * short)
    wall_growth = (1.0 + 0.45 * thick)
    r = (r[1],
         r[2] + 0.05 * (wall_growth - 1.0),
         r[3] - 0.05 * (wall_growth - 1.0),
         r[4] - 0.10 * (wall_growth - 1.0))

    θ = range(0, 2π; length = 600)
    band_colors = (RGBA(0.20,0.70,0.30,0.30),  # R green (outer, epi)
                   RGBA(0.20,0.45,0.90,0.30),  # O blue  (mid, mid-wall)
                   RGBA(0.85,0.20,0.20,0.30))  # L red   (inner, sub-endo)
    for k in 1:3
        r_out = r[k]; r_in = r[k+1]
        xs = vcat(r_out .* cos.(θ), reverse(r_in .* cos.(θ)))
        ys = vcat(r_out .* sin.(θ), reverse(r_in .* sin.(θ)))
        poly!(ax, xs, ys; color = band_colors[k])
    end
    for rr in r
        lines!(ax, rr .* cos.(θ), rr .* sin.(θ);
               color = (:black), linewidth = 0.8, linestyle = (:dash))
    end

    cmap = (RGBA(0.85,0.20,0.20,0.85),   # k=1 L red (innermost, sub-endo)
            RGBA(0.20,0.45,0.90,0.85),   # k=2 O blue (mid, mid-wall)
            RGBA(0.20,0.70,0.30,0.85))   # k=3 R green (outermost, epi)
    helix_deg = (+60.0, 0.0, -60.0)
    # glyph_radii: k=1 → innermost ring (between r[3] and r[4])
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

    ε_pct = s.short * 100
    text!(ax, 0, 3.95;
          text = @sprintf("t=%4.2fs  GLS=%+5.2f%%  twist=%+5.2f°",
                          t, ε_pct, twist_deg),
          color = (:black), align = (:center, :center), fontsize = 12)
end

function render_strip(; ages = (25, 55, 80), n_per_age::Int = 8)
    τ = MyofiberParams().τ_beat
    phases = range(0, τ; length = n_per_age)
    for age in ages
        h = initialize_heart(age = Float64(age))
        for (i, t) in enumerate(phases)
            fig = Figure(size = (640, 640), fontsize = 14)
            ax  = Axis(fig[1, 1]; aspect = DataAspect(),
                       title = "age $(age)   t = $(round(t, digits=2)) s")
            xlims!(ax, -4.4, 4.4); ylims!(ax, -4.4, 4.4)
            hidedecorations!(ax, grid = false, ticks = false)
            render_local!(ax, t, h)
            save(joinpath(OUTDIR, "strip_age$(age)_$(i).png"), fig, px_per_unit = 2)
        end
        println("strip for age $age rendered")
    end
end

function render_comparison_panel(; ages = (25, 80))
    # side-by-side at peak systole
    fig = Figure(size = (1300, 700), fontsize = 16)
    τ = MyofiberParams().τ_beat
    peak_t = τ * 0.45
    for (k, age) in enumerate(ages)
        ax = Axis(fig[1, k]; aspect = DataAspect(),
                  title = "age $(age) — peak systole")
        xlims!(ax, -4.4, 4.4); ylims!(ax, -4.4, 4.4)
        hidedecorations!(ax, grid = false, ticks = false)
        h = initialize_heart(age = Float64(age))
        render_local!(ax, peak_t, h)
    end
    save(joinpath(OUTDIR, "fig4_peak_compare.png"), fig, px_per_unit = 2)
    return joinpath(OUTDIR, "fig4_peak_compare.png")
end

if abspath(PROGRAM_FILE) == @__FILE__
    render_strip()
    p = render_comparison_panel()
    println("saved $p")
end
