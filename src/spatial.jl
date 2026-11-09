# Two-scale spatial coupling of the leaf-level model (wiki: oak/assessment.md).
#
# Every cell of a regular grid that holds oak runs the leaf-level model. A share 1 − ρ of the
# conidia a stand produces infects within the stand; the remaining share ρ escapes the canopy and
# is carried across the landscape. Conidia live for days, much less than an epidemic, so the
# airborne field is quasi-steady: it solves the screened Poisson equation
#   (1 − ℓ² Δ) φ = s,   s = host · I,
# with no-flux (Neumann) boundaries, whose kernel in the plane is K₀(r/ℓ)/(2πℓ²). The conidial
# pressure on cell i is then (1 − ρ) I_i + ρ φ_i. Primary inoculum (chasmothecia) stays local.
# The field conserves spores (Σφ = Σs) and leaves a uniform source uniform, so a landscape of
# uniform host density H reproduces the nonspatial model with b_c (1 − ρ + ρ H), whatever ℓ: a
# stand in sparse woodland receives fewer spores from the regional cloud.

"""
    Landscape(host::AbstractMatrix; cellsize)

A regular grid of square cells, `cellsize` km on a side, with `host[i, j]` the oak leaf area of
cell `(i, j)` relative to a fully stocked cell (usually between 0 and 1). Cells with `host = 0`
(open land, sea) carry airborne spores but no disease. The cells with oak, in column-major
order, are the model's stands; spatial results are vectors over them, and [`gridded`](@ref) puts
such a vector back on the grid.
"""
struct Landscape
    host::Matrix{Float64}
    cellsize::Float64
    active::Vector{Int}
    function Landscape(host::AbstractMatrix; cellsize::Real)
        isfinite(cellsize) && cellsize > 0 || throw(ArgumentError("cellsize must be finite and positive, got $cellsize"))
        all(h -> isfinite(h) && h ≥ 0, host) || throw(ArgumentError("host must be finite and nonnegative"))
        active = findall(>(0), vec(host))
        isempty(active) && throw(ArgumentError("host must be positive in at least one cell"))
        return new(Matrix{Float64}(host), Float64(cellsize), active)
    end
end

"""
    nstands(land::Landscape)

The number of cells with oak, which are the stands of the spatial model.
"""
nstands(land::Landscape) = length(land.active)

"""
    gridded(land::Landscape, values)

Put a vector over the stands of `land` (as returned by the spatial accessors) back on the grid,
with `NaN` in cells without oak.
"""
function gridded(land::Landscape, values::AbstractVector)
    length(values) == nstands(land) ||
        throw(ArgumentError("need one value per stand ($(nstands(land))), got $(length(values))"))
    out = fill(NaN, size(land.host))
    out[land.active] .= values
    return out
end

"""
    Dispersal(; rho, ell)

Two-scale conidial dispersal: a share `rho` of conidia escapes the stand, and travels a typical
distance `ell` (km), the length scale of the screened Poisson field (see
[`dispersal_field`](@ref)). `rho = 0` decouples the stands; `ell = 0` deposits escaping spores
in the cell they came from, weighted by its host density.
"""
struct Dispersal
    rho::Float64
    ell::Float64
    function Dispersal(rho::Real, ell::Real)
        0 ≤ rho ≤ 1 || throw(ArgumentError("rho must be in [0, 1], got $rho"))
        isfinite(ell) && ell ≥ 0 || throw(ArgumentError("ell must be finite and nonnegative, got $ell"))
        return new(Float64(rho), Float64(ell))
    end
end

Dispersal(; rho, ell) = Dispersal(rho, ell)

# The finite-volume operator 1 − ℓ² Δ on the grid with no-flux boundaries: the identity plus
# (ℓ/h)² times the grid graph's Laplacian. Its rows and columns sum to one.
function dispersal_operator(land::Landscape, ell::Real)
    nx, ny = size(land.host)
    c = (ell / land.cellsize)^2
    index = LinearIndices((nx, ny))
    rows, cols, vals = Int[], Int[], Float64[]
    for j in 1:ny, i in 1:nx
        k = index[i, j]
        degree = 0
        for (di, dj) in ((1, 0), (-1, 0), (0, 1), (0, -1))
            ii, jj = i + di, j + dj
            if 1 ≤ ii ≤ nx && 1 ≤ jj ≤ ny
                push!(rows, k)
                push!(cols, index[ii, jj])
                push!(vals, -c)
                degree += 1
            end
        end
        push!(rows, k)
        push!(cols, k)
        push!(vals, 1 + c * degree)
    end
    return sparse(rows, cols, vals, nx * ny, nx * ny)
end

"""
    dispersal_field(land::Landscape, ell, source::AbstractMatrix)

The quasi-steady airborne field `φ` for a source on the grid of `land`: the solution of
`(1 - ℓ² Δ) φ = source` with no-flux boundaries, discretised by finite volumes. It conserves the
source (`sum(φ) ≈ sum(source)`) and, in the plane, approaches the kernel `K₀(r/ℓ) / (2π ℓ²)`.
"""
function dispersal_field(land::Landscape, ell::Real, source::AbstractMatrix)
    size(source) == size(land.host) ||
        throw(ArgumentError("source must be $(size(land.host)), got $(size(source))"))
    isfinite(ell) && ell ≥ 0 || throw(ArgumentError("ell must be finite and nonnegative, got $ell"))
    φ = cholesky(dispersal_operator(land, ell)) \ vec(Matrix{Float64}(source))
    return reshape(φ, size(land.host))
end

"""
    SpatialMildewProblem(p::MildewParameters, land::Landscape, dispersal::Dispersal, forcings, C0)

The leaf-level model on every stand of `land`, coupled by [`Dispersal`](@ref), over consecutive
seasons. `forcings` is a [`SeasonForcing`](@ref) (one season, the same everywhere), a vector of
them (one per season, the same everywhere) or a matrix with one row per stand and one column per
season, so that phenology and weather can vary in space; all stands share each season's
`tspan`. `C0` is the primary inoculum at the start, a number for every stand or a vector with
one value per stand. Solve it with [`solve`](@ref).
"""
struct SpatialMildewProblem{F}
    p::MildewParameters
    land::Landscape
    dispersal::Dispersal
    forcings::Matrix{F}
    C0::Vector{Float64}
    function SpatialMildewProblem(p::MildewParameters, land::Landscape, dispersal::Dispersal,
                                  forcings::AbstractMatrix, C0::AbstractVector)
        n = nstands(land)
        size(forcings, 1) == n ||
            throw(ArgumentError("forcings needs one row per stand ($n), got $(size(forcings, 1))"))
        size(forcings, 2) ≥ 1 || throw(ArgumentError("need at least one season"))
        length(C0) == n || throw(ArgumentError("C0 needs one value per stand ($n), got $(length(C0))"))
        all(c -> isfinite(c) && c ≥ 0, C0) || throw(ArgumentError("C0 must be finite and nonnegative"))
        for s in axes(forcings, 2)
            spans = unique(f.tspan for f in forcings[:, s])
            length(spans) == 1 || throw(ArgumentError("season $s: all stands must share one tspan, got $spans"))
        end
        f = collect(forcings)
        return new{eltype(f)}(p, land, dispersal, f, collect(Float64, C0))
    end
end

function SpatialMildewProblem(p::MildewParameters, land::Landscape, dispersal::Dispersal, forcings, C0)
    n = nstands(land)
    F = forcings isa SeasonForcing ? [forcings] : forcings
    M = F isa AbstractMatrix ? F : permutedims(repeat(collect(F), 1, n))
    C = C0 isa Real ? fill(Float64(C0), n) : C0
    return SpatialMildewProblem(p, land, dispersal, M, C)
end

nseasons(prob::SpatialMildewProblem) = size(prob.forcings, 2)

struct SpatialRHS{F,FA}
    p::MildewParameters
    land::Landscape
    rho::Float64
    forcings::Vector{F}
    factor::FA
    source::Vector{Float64}
end

function (r::SpatialRHS)(du, u, _, t)
    p, land, rho = r.p, r.land, r.rho
    ns, iI = nstate(p), sporulating(p)
    U = reshape(u, ns, :)
    dU = reshape(du, ns, :)
    if rho > 0
        fill!(r.source, 0.0)
        for (j, k) in enumerate(land.active)
            r.source[k] = land.host[k] * U[iI, j]
        end
        φ = r.factor \ r.source
    else
        φ = r.source
    end
    for (j, k) in enumerate(land.active)
        conidia = (1 - rho) * U[iI, j] + (rho > 0 ? rho * φ[k] : 0.0)
        leaf_rhs!(view(dU, :, j), view(U, :, j), p, r.forcings[j], t, conidia)
    end
    return nothing
end

"""
    SpatialMildewSolution

The solution of a [`SpatialMildewProblem`](@ref): one dense ODE solution per season over all
stands. `sol(season, t)` is the state on day `t`, a matrix with one column per stand laid out as
in [`MildewSolution`](@ref). [`leaf_area`](@ref), [`incidence`](@ref) and
[`visible_incidence`](@ref) return one value per stand ([`gridded`](@ref) maps them back), and
[`regional_incidence`](@ref) the host-weighted mean.
"""
struct SpatialMildewSolution{P,S}
    prob::P
    seasons::Vector{S}
    winter::Matrix{Float64}
end

function (sol::SpatialMildewSolution)(season::Integer, t::Real)
    S = length(sol.seasons)
    1 ≤ season ≤ S || throw(ArgumentError("season must be in 1:$S, got $season"))
    t0, t1 = sol.prob.forcings[1, season].tspan
    t0 ≤ t ≤ t1 || throw(ArgumentError("t must be in [$t0, $t1] for season $season, got $t"))
    return reshape(sol.seasons[season](t), nstate(sol.prob.p), :)
end

"""
    solve(prob::SpatialMildewProblem; alg = Tsit5(), abstol = 1e-12, reltol = 1e-10)

Integrate each season of `prob` over all stands at once, with the airborne field solved at every
right-hand-side evaluation from one sparse Cholesky factorisation, and apply each stand's winter
map between seasons. Returns a [`SpatialMildewSolution`](@ref).
"""
function CommonSolve.solve(prob::SpatialMildewProblem; alg = Tsit5(), abstol = 1e-12, reltol = 1e-10)
    p, land = prob.p, prob.land
    ns, n = nstate(p), nstands(land)
    factor = cholesky(dispersal_operator(land, prob.dispersal.ell))
    C = copy(prob.C0)
    seasons = []
    winter = zeros(n, nseasons(prob))
    for s in 1:nseasons(prob)
        fs = prob.forcings[:, s]
        u0 = zeros(ns * n)
        U0 = reshape(u0, ns, n)
        for j in 1:n
            season_start!(view(U0, :, j), p, fs[j], C[j])
        end
        rhs = SpatialRHS(p, land, prob.dispersal.rho, fs, factor, zeros(length(land.host)))
        tstops = sort!(unique!(reduce(vcat, (f.tstops for f in fs))))
        ode = ODEProblem(rhs, u0, fs[1].tspan)
        sol = CommonSolve.solve(ode, alg; abstol, reltol, tstops)
        SciMLBase.successful_retcode(sol) ||
            throw(SolverError("season $s integration failed with retcode $(sol.retcode)"))
        push!(seasons, sol)
        Uend = reshape(sol.u[end], ns, n)
        for j in 1:n
            C[j] = winter_source(p, fs[j], view(Uend, :, j))
            winter[j, s] = C[j]
        end
    end
    return SpatialMildewSolution(prob, [s for s in seasons], winter)
end

function leaf_area(sol::SpatialMildewSolution, season::Integer, t::Real)
    p, U = sol.prob.p, sol(season, t)
    return [stand_leaves(p, view(U, :, j)) for j in axes(U, 2)]
end

function incidence(sol::SpatialMildewSolution, season::Integer, t::Real)
    p, U = sol.prob.p, sol(season, t)
    return [stand_incidence(p, view(U, :, j)) for j in axes(U, 2)]
end

function visible_incidence(sol::SpatialMildewSolution, season::Integer, t::Real; latent_visible::Real = 0.0)
    check_latent_visible(latent_visible)
    p, U = sol.prob.p, sol(season, t)
    return [stand_visible_incidence(p, view(U, :, j), latent_visible) for j in axes(U, 2)]
end

"""
    regional_incidence(sol::SpatialMildewSolution, season, t; visible = false, latent_visible = 0.0)

The share of all the region's leaves on day `t` of `season` that are infected, weighting each
stand by its host density; with `visible = true`, the share with visible mildew as in
[`visible_incidence`](@ref).
"""
function regional_incidence(sol::SpatialMildewSolution, season::Integer, t::Real; visible::Bool = false,
                            latent_visible::Real = 0.0)
    check_latent_visible(latent_visible)
    p, U, land = sol.prob.p, sol(season, t), sol.prob.land
    leaves = sick = 0.0
    for (j, k) in enumerate(land.active)
        u = view(U, :, j)
        leaves += land.host[k] * stand_leaves(p, u)
        sick += land.host[k] * (visible ? stand_visible(p, u, latent_visible) : stand_infected(p, u))
    end
    return leaves > 0 ? sick / leaves : 0.0
end

"""
    winter_inoculum(sol::SpatialMildewSolution)

The primary inoculum carried into the season after each season: a matrix with one row per stand
and one column per season.
"""
winter_inoculum(sol::SpatialMildewSolution) = copy(sol.winter)
