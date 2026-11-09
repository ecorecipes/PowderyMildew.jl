# PowderyMildew.jl

[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://ecorecipes.github.io/PowderyMildew.jl/dev/)

Seasonal, leaf-level models of oak powdery mildew (*Erysiphe alphitoides* on *Quercus robur* and
*Q. petraea*), for asking how the timing of oak flushing and of primary inoculum shapes
severity, and how climate change shifts it.

- **Hamelin et al. (2016)**, *Ecosphere* 7:e01517: the seasonal two-species model of oak
  mildew, implemented exactly (`HamelinProblem`), with its annual map, its binomial likelihood
  of species frequencies, and the coexistence-equilibrium inference of the between-season
  coefficients (`coexistence_chi`), which reproduces the paper's Table 2.
- **The leaf-level model** (`MildewProblem`): leaf flushes, Erlang leaf ageing and latency,
  secondary (conidial) and primary (ascospore) infection, weather forcing, and a winter map
  through chasmothecia. It nests the one-species Hamelin model exactly, and that limit is
  tested.
- **The spatial model** (`SpatialMildewProblem`): the leaf-level model on a grid of oak stands,
  coupled by two-scale dispersal (a share of conidia escapes the stand into a quasi-steady
  screened-Poisson field), for the speed of spread and its dependence on the woodland.

Tutorials (Quarto, in `vignettes/`, rendered to the
[documentation](https://ecorecipes.github.io/PowderyMildew.jl/dev/)): the Hamelin model and its
Table 2; timing and severity in the leaf-level model; spread.

```julia
using PowderyMildew

p = MildewParameters(; b_c = 0.3, b_a = 0.5, tau_s = 20, K = 3, epsilon = 0.05, latent = 8, L = 3,
                     gamma_I = 0.05, c_w = 1.0)
season = SeasonForcing(; tspan = (90, 280), production = FlushProduction([120, 190], [5, 8], [0.8, 0.2]),
                       release = ReleaseHazard(; rate = 0.2, optimum = 140, curvature = -0.003),
                       gate = t -> t > 240)
sol = solve(MildewProblem(p, fill(season, 5), 1.0))
visible_incidence(sol, 5, 200.0)    # 0.155: settling to an endemic cycle
winter_inoculum(sol)                # 0.379, 0.303, 0.275, 0.262, 0.256
```

These parameter values are illustrative, not fitted.

The design, and an assessment of the sources it rests on, is in the workspace wiki
(`../../wiki/oak/assessment.md`).

## Status

Not registered, and not yet fitted to data. Tests: `julia --project -e 'using Pkg; Pkg.test()'`.
