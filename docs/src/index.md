# PowderyMildew.jl

```@meta
CurrentModule = PowderyMildew
```

```@docs
PowderyMildew
```

PowderyMildew.jl models oak powdery mildew (*Erysiphe alphitoides* and relatives on
*Quercus robur* and *Q. petraea*): when epidemics start, how severe they become, and how they
spread. It has three layers, each tested against the one below.

1. **Hamelin et al. (2016)**, *Ecosphere* 7:e01517. The seasonal two-species model, exactly
   ([`HamelinProblem`](@ref)), with its annual map, the binomial likelihood of species
   frequencies, and the coexistence closed form for the between-season coefficients
   ([`coexistence_chi`](@ref)), which reproduces the paper's Table 2.
2. **The leaf-level model** ([`MildewProblem`](@ref)). Leaf flushes, Erlang leaf ageing and
   latency, secondary (conidial) and primary (ascospore) infection with a timed release,
   weather forcing, and a winter map through chasmothecia. It nests the one-species Hamelin
   model exactly.
3. **The spatial model** ([`SpatialMildewProblem`](@ref)). The leaf-level model on every cell of
   a grid of oak stands, coupled by two-scale dispersal: a share 1 − ρ of conidia infects within
   the stand and a share ρ escapes into a quasi-steady screened-Poisson field with length scale
   ℓ. A uniform landscape reduces to the leaf-level model.

## Tutorials

- [The model of Hamelin et al. (2016)](tutorials/01_hamelin.md): the model, Table 2, and why its
  between-season coefficients only restate the assumed initial incidence.
- [Timing and severity](tutorials/02_leaf_model.md): the leaf-level model, its nesting of
  Hamelin's, and how the timing of flushing against ascospore release sets severity.
- [Spread](tutorials/03_spread.md): the dispersal field, front speeds, and spread across a patchy
  landscape.

The tutorials use illustrative parameter values; the package has not yet been fitted to data.

## Quick start

```julia
using PowderyMildew

p = MildewParameters(; b_c = 0.25, b_a = 0.5, tau_s = 20, K = 3, epsilon = 0.05, latent = 8, L = 3,
                     gamma_I = 0.05, c_w = 2.0)
season = SeasonForcing(; tspan = (90, 280), production = FlushProduction([120, 190], [5, 8], [0.8, 0.2]),
                       release = ReleaseHazard(; rate = 0.2, optimum = 140, curvature = -0.003),
                       gate = t -> t > 240 ? 1.0 : 0.0, tstops = [240.0])
sol = solve(MildewProblem(p, fill(season, 5), 1.0))
visible_incidence(sol, 5, 280.0)

land = Landscape(ones(100, 1); cellsize = 1.0)
spread = solve(SpatialMildewProblem(p, land, Dispersal(; rho = 0.2, ell = 4.0), fill(season, 5),
                                    [1.0; zeros(99)]))
visible_incidence(spread, 5, 280.0)    # one value per stand
```
