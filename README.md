# LV three-layer myofiber model

A Julia implementation of a three-transmural-layer left-ventricular (LV)
myofiber mechanics model, with a Streeter-style transmural helix-angle
rotation and an age-dependent layer-failure / compensation rule.  It
reproduces the well-known sequence:

```
longitudinal (L)  →  radial-oblique (R)  →  both de-compensate
        ↑ first to fail              ↑ compensates             ↑ eventually fails
```

and shows how the **Global Longitudinal Strain (GLS)** tracks this
progression.

## Layout

```
.
├── Project.toml               # Julia project (CairoMakie, GeometryBasics, Colors, FFMPEG_jll)
├── src/ThreeLayerModel.jl     # physics: layers, ageing, beat kinematics, GLS
├── scripts/
│   ├── plot_static.jl         # writes the three PNGs in out/
│   ├── animate_beat.jl        # writes the two MP4s in out/
│   └── render_strip.jl        # writes per-age frame strips + peak-systole panel
└── out/
    ├── fig1_helix.png          # short-axis schematic of the three layers
    ├── fig2_gls_curves.png     # GLS vs time at ages 25 / 55 / 80
    ├── fig3_layer_dominance.png# per-layer GLS contribution & health over age
    ├── fig4_peak_compare.png   # young vs old at peak systole
    ├── strip_age{25,55,80}_*.png
    ├── anim_beat.mp4           # one beat (single heart, age 25)
    └── anim_compare_young_old.mp4  # one beat (young vs old side-by-side)
```

## Running

```bash
# first-time setup (already done; only needed on a fresh checkout)
julia --project=. -e 'using Pkg; Pkg.add(["CairoMakie","GeometryBasics","FFMPEG_jll","Colors"])'

# generate the static figures
julia --project=. scripts/plot_static.jl

# generate the animations (~1 minute)
julia --project=. scripts/animate_beat.jl

# per-age frame strips + side-by-side peak panel
julia --project=. scripts/render_strip.jl
```

## What the model does

* **Three concentric transmural layers**, with canonical Streeter
  helix angles: L (sub-endocardial) +60°, O (mid-wall) 0°,
  R (sub-epicardial) −60°.
* **Beat kinematics**: each layer contracts with a Gaussian-shaped
  systolic envelope and a small phase offset (dyssynchrony).  L, O, R
  contribute 50 / 30 / 20 % to GLS and 20 / 30 / 50 % to radial
  thickening.  Net LV twist is the difference between the endo and epi
  helicities (signed), producing realistic twist ≈ 15–17° in young,
  ≈ 8–10° in aged hearts.
* **Ageing & compensation**: longitudinal health declines linearly from
  age 35; the radial layer hypertrophies up to a `compensation` ceiling
  to offset the lost strain.  Once both L and R fail (around age 60+),
  no further compensation is available.
* **Clinical readout**: peak GLS = −19.7 % at age 25, −10.9 % at age 80
  (twist halves similarly).

## Validation highlights

| metric                  | young (25)  | aged (80)   |
|-------------------------|-------------|-------------|
| layer health  L/O/R     | 1 / 1 / 1   | 0.19 / 0.96 / 0.85 |
| peak GLS                | −19.7 %     | −10.9 %     |
| peak twist              | +16.7°      | +8.7°       |

These match the published ranges for healthy young and elderly adults
(Negishi et al. 2017, etc.).
