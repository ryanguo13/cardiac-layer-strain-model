# ------------------------------------------------------------------
#  Plots for the LV three-layer myofiber model
#  Generates three static PNGs and one MP4 animation in ./out/
# ------------------------------------------------------------------

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using CairoMakie
using GeometryBasics
using Colors
using Statistics
using Printf

include(joinpath(@__DIR__, "..", "src", "ThreeLayerModel.jl"))
using .ThreeLayerModel

const OUTDIR = joinpath(@__DIR__, "..", "out")
isdir(OUTDIR) || mkpath(OUTDIR)

# ------------------------------------------------------------------
#  Fig 1 — short-axis cross-section of the three transmural layers
# ------------------------------------------------------------------

function fig1_helix(; out = joinpath(OUTDIR, "fig1_helix.png"))
    fig = Figure(size = (1100, 900), fontsize = 18)
    ax  = Axis(fig[1, 1];
               aspect = DataAspect(),
               title  = "LV short-axis — three transmural myofiber layers",
               xlabel = "x (cm)", ylabel = "y (cm)")

    # short-axis concentric circles (epi / mid / endo)
    θ = range(0, 2π; length = 600)
    radii = (3.2, 2.7, 2.2, 1.5)        # epi-outer, epi-inner, endo-outer, endo-inner
    # r = (epi-outer, epi-inner=mid-outer, mid-inner=endo-outer, endo-inner)
    # anatomy: endo = innermost (L, +60°), mid = O (0°), epi = outermost (R, –60°)
    # annulus k fills region r[k+1] .. r[k]:
    #   k=1 (outer): R (epi, green) — outer ring
    #   k=2 (mid) : O (mid-wall, blue) — middle ring
    #   k=3 (inner): L (sub-endo, red) — innermost ring
    layer_fills = (
        RGBA(0.20,0.70,0.30,0.30),  # R green (outermost, epi)
        RGBA(0.20,0.45,0.90,0.30),  # O blue  (middle, mid-wall)
        RGBA(0.85,0.20,0.20,0.30),  # L red   (innermost, sub-endo)
    )

    # fill the annular rings between consecutive radii
    for i in 1:length(radii)-1
        r_out = radii[i]
        r_in  = radii[i+1]
        xs = vcat(r_out .* cos.(θ), reverse(r_in .* cos.(θ)))
        ys = vcat(r_out .* sin.(θ), reverse(r_in .* sin.(θ)))
        poly!(ax, xs, ys; color = layer_fills[i])
    end
    # crisp ring outlines
    for r in radii
        lines!(ax, r .* cos.(θ), r .* sin.(θ);
               color = (:black), linewidth = 0.8, linestyle = (:dash))
    end

    # helix-angle glyphs around the ring
    cmap = (
        RGBA(0.85, 0.20, 0.20, 0.85),   # L red
        RGBA(0.20, 0.45, 0.90, 0.85),   # O blue
        RGBA(0.20, 0.70, 0.30, 0.85),   # R green
    )
    labels = ["L (longitudinal, sub-endo, +60°)",
              "O (oblique, mid-wall, 0°)",
              "R (radial-oblique, epi, –60°)"]
    radii_layers = (1.85, 2.45, 2.95)   # L (innermost), O (mid), R (outer)
    helix_deg    = (+60.0, 0.0, -60.0)

    # in-ring layer labels (positions match anatomy: L innermost, O middle, R outermost)
    text!(ax, 0.0,  1.80;
          text = "L\n(longitudinal,\nsub-endo, +60°)",
          color = cmap[1],
          align = (:center, :center), fontsize = 12)
    text!(ax, 0.0,  2.45;
          text = "O\n(oblique,\nmid-wall, 0°)",
          color = cmap[2],
          align = (:center, :center), fontsize = 12)
    text!(ax, 3.4,  0.0;
          text = "R\n(radial-oblique,\nepi, –60°)",
          color = cmap[3],
          align = (:left, :center), fontsize = 12)

    for (k, r0) in enumerate(radii_layers)
        for φ in range(0, 2π; length = 24)[1:end-1]
            x, y = r0 * cos(φ), r0 * sin(φ)
            β = deg2rad(helix_deg[k])            # helix pitch
            # short bar tangent to ring, rotated by helix pitch
            vx, vy = -sin(φ + β), cos(φ + β)
            L = 0.55
            lines!(ax, [x - vx*L/2, x + vx*L/2], [y - vy*L/2, y + vy*L/2];
                   color = cmap[k], linewidth = 2.4)
            # arrow tip
            arr = scatter!(ax, [x + vx*L/2], [y + vy*L/2];
                           color = cmap[k], markersize = 5)
        end
    end

    # legend
    legend_elements = [
        LineElement(color = cmap[1], linewidth = 4, linestyle = :solid),
        LineElement(color = cmap[2], linewidth = 4, linestyle = :solid),
        LineElement(color = cmap[3], linewidth = 4, linestyle = :solid),
    ]
    Legend(fig[1, 2], legend_elements, labels;
           framevisible = true)

    # annotations
    text!(ax, 0, 0; text = "LV cavity", color = (:black),
          align = (:center, :center), fontsize = 16)
    text!(ax, 0,  3.65; text = "epicardium", color = (:black),
          align = (:center, :center), fontsize = 13)

    xlims!(ax, -4.4, 4.4); ylims!(ax, -4.4, 4.4)
    hidedecorations!(ax, grid = false, ticks = false)
    save(out, fig, px_per_unit = 2)
    return out
end

# ------------------------------------------------------------------
#  Fig 2 — GLS curves at three representative ages
# ------------------------------------------------------------------

function fig2_gls_curves(; out = joinpath(OUTDIR, "fig2_gls_curves.png"))
    ages = (25, 55, 80)
    fig  = Figure(size = (1100, 700), fontsize = 18)
    ax   = Axis(fig[1, 1];
                xlabel = "time (s)",
                ylabel = "global longitudinal strain ε_L (%)",
                title  = "Global longitudinal strain (GLS) over one beat")
    ax.yreversed = true                      # clinical convention (shortening = negative)

    cmap = (RGBA(0.10,0.55,0.10,0.95),
            RGBA(0.95,0.55,0.10,0.95),
            RGBA(0.85,0.10,0.10,0.95))
    peak_vals = Float64[]

    for (k, a) in enumerate(ages)
        h = initialize_heart(age = Float64(a))
        t, ε = gls_trace(h; n = 400)
        ε_pct = ε .* 100
        lines!(ax, t, ε_pct; color = cmap[k], linewidth = 3,
               label = "age $(a)")
        pk, ipk = findmin(ε_pct)
        push!(peak_vals, pk)
        scatter!(ax, [t[ipk]], [pk]; color = cmap[k], markersize = 12)
        text!(ax, t[ipk]+0.04, pk+0.4;
              text = @sprintf("peak %.1f%%", pk),
              color = cmap[k], fontsize = 13, align = (:left, :center))
    end

    axislegend(ax; position = :rb)
    save(out, fig, px_per_unit = 2)
    return out
end

# ------------------------------------------------------------------
#  Fig 3 — layer contribution to GLS over age (stacked area)
# ------------------------------------------------------------------

function fig3_layer_dominance(; out = joinpath(OUTDIR, "fig3_layer_dominance.png"))
    ages = range(20, 90; length = 80)
    Lc = Float64[]; Oc = Float64[]; Rc = Float64[]
    healthL = Float64[]; healthR = Float64[]
    for a in ages
        h = initialize_heart(age = Float64(a))
        c = layer_strain_contributions(h)
        push!(Lc, -c.L * 100); push!(Oc, -c.O * 100); push!(Rc, -c.R * 100)
        push!(healthL, h.L.health); push!(healthR, h.R.health)
    end

    fig = Figure(size = (1200, 700), fontsize = 18)
    ax1 = Axis(fig[1, 1];
               xlabel = "age (yr)",
               ylabel = "GLS contribution (%)",
               title  = "Per-layer contribution to GLS — longitudinal\nfails first, radial compensates")
    ax1.yreversed = true

    # stacked area: positive numbers plotted vs age
    lower = zeros(length(ages))
    for (series, col, name) in ((Lc, RGBA(0.85,0.20,0.20,0.55), "L (longitudinal)"),
                                (Oc, RGBA(0.20,0.45,0.90,0.55), "O (oblique)"),
                                (Rc, RGBA(0.20,0.70,0.30,0.55), "R (radial-oblique)"))
        upper = lower .+ series
        band!(ax1, collect(ages), lower, upper; color = col, label = name)
        lower = upper
    end
    axislegend(ax1; position = :lb)

    # overlay GLS line (sum)
    ax2 = Axis(fig[1, 2];
               xlabel = "age (yr)",
               ylabel = "layer health",
               title  = "Layer health over age")
    lines!(ax2, collect(ages), healthL; color = RGBA(0.85,0.20,0.20,0.9),
           linewidth = 3, label = "L health")
    lines!(ax2, collect(ages), healthR; color = RGBA(0.20,0.70,0.30,0.9),
           linewidth = 3, label = "R health")
    lines!(ax2, collect(ages), fill(1.0, length(ages));
           color = (:grey, 0.4), linewidth = 1, linestyle = :dash)
    axislegend(ax2; position = :rb)

    save(out, fig, px_per_unit = 2)
    return out
end

# ------------------------------------------------------------------
#  Driver
# ------------------------------------------------------------------

function main()
    println("→ fig1_helix")
    p1 = fig1_helix()
    println("  saved $p1")
    println("→ fig2_gls_curves")
    p2 = fig2_gls_curves()
    println("  saved $p2")
    println("→ fig3_layer_dominance")
    p3 = fig3_layer_dominance()
    println("  saved $p3")
    return (p1, p2, p3)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
