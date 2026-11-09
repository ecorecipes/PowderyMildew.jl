# The seasonal model of Hamelin, Bisson, Desprez-Loustau, Fabre & Mailleret (2016),
# "Temporal niche differentiation of parasites sharing the same plant host: oak powdery
# mildew as a case study", Ecosphere 7:e01517, doi:10.1002/ecs2.1517.
#
# Within a season of length τ, for species i = 1..n,
#     dS/dt = -S Σⱼ βⱼ Iⱼ,    dIᵢ/dt = Iᵢ (βᵢ S - α),    dR/dt = α Σⱼ Iⱼ,
# where R (removed leaves) is the bookkeeping completion of the displayed S-I system. Between
# seasons, with H = Σⱼ χⱼ Iⱼ(τ),
#     S(0) = S₀ e^{-H},    Iᵢ(0) = S₀ (1 - e^{-H}) χᵢ Iᵢ(τ) / H,    R(0) = 0.
# The fitted version takes α = 0 and S₀ = 1 (leaf densities normalised by the total), with
# time in days.

"""
    HamelinParameters(; beta, chi, alpha = 0.0, S0 = 1.0, season_length = 200.0)

Parameters of the seasonal model of Hamelin et al. (2016) for `n = length(beta)` parasite
species sharing one host: within-season transmission rates `beta` (per day, per unit leaf
density), composite between-season coefficients `chi` (production of survival forms, their
survival and primary-infection success together; not a survival probability), the removal
rate `alpha` (per day), the leaf density `S0` renewed each year before primary infection, and
the season length `season_length` (days). The paper's fit uses `alpha = 0` and `S0 = 1`, with
species 1 *Erysiphe alphitoides* and species 2 *E. quercicola*.
"""
struct HamelinParameters
    beta::Vector{Float64}
    chi::Vector{Float64}
    alpha::Float64
    S0::Float64
    season_length::Float64
    function HamelinParameters(beta::Vector{Float64}, chi::Vector{Float64}, alpha::Float64,
                               S0::Float64, season_length::Float64)
        n = length(beta)
        n ≥ 1 || throw(ArgumentError("beta must have at least one species"))
        length(chi) == n ||
            throw(ArgumentError("beta and chi must have the same length, got $n and $(length(chi))"))
        all(b -> isfinite(b) && b ≥ 0, beta) || throw(ArgumentError("beta must be finite and nonnegative, got $beta"))
        all(c -> isfinite(c) && c ≥ 0, chi) || throw(ArgumentError("chi must be finite and nonnegative, got $chi"))
        isfinite(alpha) && alpha ≥ 0 || throw(ArgumentError("alpha must be finite and nonnegative, got $alpha"))
        isfinite(S0) && S0 > 0 || throw(ArgumentError("S0 must be finite and positive, got $S0"))
        isfinite(season_length) && season_length > 0 ||
            throw(ArgumentError("season_length must be finite and positive, got $season_length"))
        return new(beta, chi, alpha, S0, season_length)
    end
end

function HamelinParameters(; beta, chi, alpha = 0.0, S0 = 1.0, season_length = 200.0)
    return HamelinParameters(collect(Float64, beta), collect(Float64, chi), Float64(alpha), Float64(S0),
                             Float64(season_length))
end

"""
    nspecies(p::HamelinParameters)

The number of parasite species in `p`.
"""
nspecies(p::HamelinParameters) = length(p.beta)

"""
    HamelinProblem(p::HamelinParameters, I0; S = p.S0 - sum(I0), seasons = 1)

The model of Hamelin et al. (2016) run for `seasons` growing seasons, starting from infected
leaf densities `I0` (one per species) and susceptible density `S` at the start of the first
season. The default `S` starts the first season as the annual map would, with `S + sum(I0) = S0`.
Solve it with [`solve`](@ref).
"""
struct HamelinProblem
    p::HamelinParameters
    u0::Vector{Float64}
    seasons::Int
end

function HamelinProblem(p::HamelinParameters, I0; S = p.S0 - sum(I0), seasons::Integer = 1)
    length(I0) == nspecies(p) ||
        throw(ArgumentError("I0 must have one entry per species ($(nspecies(p))), got $(length(I0))"))
    all(≥(0), I0) || throw(ArgumentError("I0 must be nonnegative, got $I0"))
    S ≥ 0 || throw(ArgumentError("S must be nonnegative, got $S"))
    seasons ≥ 1 || throw(ArgumentError("seasons must be at least 1, got $seasons"))
    return HamelinProblem(p, [Float64(S); collect(Float64, I0); 0.0], Int(seasons))
end

# state u = [S; I₁ … Iₙ; R]
function hamelin_rhs!(du, u, p::HamelinParameters, t)
    n = nspecies(p)
    S = u[1]
    force = zero(S)
    total = zero(S)
    @inbounds for i in 1:n
        force += p.beta[i] * u[1 + i]
        total += u[1 + i]
    end
    du[1] = -S * force
    @inbounds for i in 1:n
        du[1 + i] = u[1 + i] * (p.beta[i] * S - p.alpha)
    end
    du[n + 2] = p.alpha * total
    return nothing
end

"""
    between_season(p::HamelinParameters, u)

The annual map of Hamelin et al. (2016): from the state `u = [S; I₁ … Iₙ; R]` at the end of a
season to the state at the start of the next. With `H = Σ χᵢ Iᵢ`, susceptible leaves start at
`S0 exp(-H)` and species `i` infects `S0 (1 - exp(-H)) χᵢ Iᵢ / H`, the limit as `H → 0`
being no infection; `R` restarts at zero.
"""
function between_season(p::HamelinParameters, u::AbstractVector)
    n = nspecies(p)
    length(u) == n + 2 || throw(ArgumentError("state must have length $(n + 2), got $(length(u))"))
    weighted = [p.chi[i] * u[1 + i] for i in 1:n]
    H = sum(weighted)
    if H == 0
        return [p.S0; zeros(n); 0.0]
    end
    infected = -expm1(-H) * p.S0
    return [p.S0 * exp(-H); [infected * w / H for w in weighted]; 0.0]
end

"""
    HamelinSolution

The solution of a [`HamelinProblem`](@ref): one dense ODE solution per season, each over days
`0 … season_length` of that season. Call `sol(season, t)` for the state `[S; I₁ … Iₙ; R]` on day
`t` of `season`, and see [`frequency`](@ref) and [`incidence`](@ref).
"""
struct HamelinSolution{S}
    prob::HamelinProblem
    seasons::Vector{S}
end

function (sol::HamelinSolution)(season::Integer, t::Real)
    1 ≤ season ≤ length(sol.seasons) ||
        throw(ArgumentError("season must be in 1:$(length(sol.seasons)), got $season"))
    0 ≤ t ≤ sol.prob.p.season_length ||
        throw(ArgumentError("t must be in [0, $(sol.prob.p.season_length)], got $t"))
    return sol.seasons[season](t)
end

"""
    solve(prob::HamelinProblem; alg = Tsit5(), abstol = 1e-12, reltol = 1e-10)

Integrate each season of `prob` with the ODE solver `alg` and apply [`between_season`](@ref)
between seasons. Returns a [`HamelinSolution`](@ref).
"""
function CommonSolve.solve(prob::HamelinProblem; alg = Tsit5(), abstol = 1e-12, reltol = 1e-10)
    p = prob.p
    seasons = [solve_season(p, copy(prob.u0), alg, abstol, reltol)]
    for _ in 2:prob.seasons
        u = between_season(p, seasons[end].u[end])
        push!(seasons, solve_season(p, u, alg, abstol, reltol))
    end
    return HamelinSolution(prob, seasons)
end

function solve_season(p, u, alg, abstol, reltol)
    ode = ODEProblem(hamelin_rhs!, u, (0.0, p.season_length), p)
    sol = CommonSolve.solve(ode, alg; abstol, reltol)
    SciMLBase.successful_retcode(sol) || throw(SolverError("season integration failed with retcode $(sol.retcode)"))
    return sol
end

"""
    frequency(sol::HamelinSolution, season, t; species = 1)

The share of infected leaves infected by `species` on day `t` of `season`, `Iᵢ / Σⱼ Iⱼ`
(the paper's F for species 1). It is `NaN` when no leaf is infected.
"""
function frequency(sol::HamelinSolution, season::Integer, t::Real; species::Integer = 1)
    n = nspecies(sol.prob.p)
    1 ≤ species ≤ n || throw(ArgumentError("species must be in 1:$n, got $species"))
    u = sol(season, t)
    total = sum(@view u[2:(n + 1)])
    return total > 0 ? u[1 + species] / total : NaN
end

"""
    incidence(sol::HamelinSolution, season, t)

The share of leaves infected by any species on day `t` of `season`, `Σⱼ Iⱼ / (S + Σⱼ Iⱼ + R)`.
"""
function incidence(sol::HamelinSolution, season::Integer, t::Real)
    n = nspecies(sol.prob.p)
    u = sol(season, t)
    return sum(@view u[2:(n + 1)]) / sum(u)
end

"""
    coexistence_chi(p::HamelinParameters, I0)

The between-season coefficients χ for which the start-of-season state with infected densities
`I0` and `S = S0 - sum(I0)` is a fixed point of one season followed by
[`between_season`](@ref), given the within-season parameters of `p` (whose own `chi` is
ignored). Hamelin et al. (2016) inferred χ this way, under a coexistence-equilibrium
assumption, rather than from winter observations. With `H = -log(1 - Σ I0 / S0)` and `Iᵢ(τ)`
the end-of-season densities, `χᵢ = H I0ᵢ / (Σⱼ I0ⱼ Iᵢ(τ))`.
"""
function coexistence_chi(p::HamelinParameters, I0; alg = Vern9(), abstol = 1e-14, reltol = 1e-13)
    length(I0) == nspecies(p) ||
        throw(ArgumentError("I0 must have one entry per species ($(nspecies(p))), got $(length(I0))"))
    all(>(0), I0) || throw(ArgumentError("every species must be present for coexistence, got I0 = $I0"))
    total0 = sum(I0)
    total0 < p.S0 || throw(ArgumentError("sum(I0) must be below S0 = $(p.S0), got $total0"))
    season = solve_season(p, [p.S0 - total0; collect(Float64, I0); 0.0], alg, abstol, reltol)
    u = season.u[end]
    H = -log1p(-total0 / p.S0)
    return [H * I0[i] / (total0 * u[1 + i]) for i in eachindex(I0)]
end

"""
    frequency_loglikelihood(sol::HamelinSolution, season, t, sampled, alphitoides; species = 1)

The binomial log-likelihood of Hamelin et al. (2016): on observation `k`, `alphitoides[k]` of
the `sampled[k]` singly infected leaves sampled on day `t[k]` of `season[k]` carry `species`,
with probability [`frequency`](@ref). `sampled` counts identified single infections, not all
leaves examined. It is `-Inf` when an observed count is impossible, including when no leaf
is infected.
"""
function frequency_loglikelihood(sol::HamelinSolution, season::AbstractVector{<:Integer}, t::AbstractVector,
                                 sampled::AbstractVector{<:Integer}, alphitoides::AbstractVector{<:Integer};
                                 species::Integer = 1)
    K = length(season)
    length(t) == length(sampled) == length(alphitoides) == K ||
        throw(ArgumentError("season, t, sampled and alphitoides must have the same length"))
    ll = 0.0
    for k in 1:K
        0 ≤ alphitoides[k] ≤ sampled[k] ||
            throw(ArgumentError("observation $k: need 0 ≤ alphitoides ≤ sampled, got $(alphitoides[k]) and $(sampled[k])"))
        sampled[k] == 0 && continue
        F = frequency(sol, season[k], t[k]; species)
        isnan(F) && return -Inf
        ll += logpdf(Binomial(sampled[k], clamp(F, 0.0, 1.0)), alphitoides[k])
    end
    return ll
end
