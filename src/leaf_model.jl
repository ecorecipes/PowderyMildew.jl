# The leaf-level seasonal model proposed for oak powdery mildew (wiki: oak/assessment.md).
#
# Per stand, as fractions of the stand's maximal leaf number, with t in days:
#   U₁ … U_K   young susceptible leaves, aged by an Erlang chain (mean susceptible time τ_s)
#   M          mature leaves, susceptibility ε relative to young ones
#   E₁ … E_L   latent infected leaves (Erlang, mean latent period `latent`; none if L = 0)
#   I          sporulating leaves
#   Z          leaves that have stopped sporulating, still diseased
#   C          unreleased primary inoculum, in units of potential infections
#   W          ∫ η(t) I(t) dt, the late-season source of next year's C
#
#   λ = f(t) (b_c I + b_a g(t) C + b_ext)
#   U₁' = l(t) − (λ + K/τ_s) U₁,   U_k' = (K/τ_s)(U_{k−1} − U_k) − λ U_k,
#   M'  = (K/τ_s) U_K − ε λ M,
#   E₁' = λ (ΣU + ε M) − r E₁, …,   I' = r E_L − (γ_I + δ) I,   Z' = γ_I I − δ Z,   r = (L/latent) d(t),
#   C'  = −g(t) C,   W' = η(t) I,   with d(t) the relative rate of latent development (1 by default)
# and δ the rate at which diseased leaves are shed (0 by default; shed leaves leave the stand),
# and the winter map C_{a+1} = c_w W_a (or c_w I_a(end) with `gate = nothing`). With all
# leaves present at the start, τ_s = ∞, no latency, a pulse release and the end-of-season
# gate, it is the one-species model of Hamelin et al. (2016) with β = b_c, α = γ_I and
# χ = b_a c_w.

"""
    MildewParameters(; b_c, b_a, b_ext = 0.0, tau_s = Inf, K = 1, epsilon = 0.0, latent = 0.0,
                     L = 0, gamma_I = 0.0, c_w = 1.0, shed = 0.0)

Parameters of the leaf-level model:
- `b_c`: secondary (conidial) transmission per day per unit sporulating leaf fraction;
- `b_a`: infections per unit of released primary inoculum (ascospores);
- `b_ext`: background force of infection from outside the stand, per day;
- `tau_s`: mean time (days) a young leaf stays fully susceptible, Erlang with `K` stages (`Inf`
  for leaves that never mature);
- `epsilon`: susceptibility of mature leaves relative to young ones;
- `latent`: mean latent period (days) at unit development rate (see `development` in
  [`SeasonForcing`](@ref)), Erlang with `L` stages; `L = 0` makes infected leaves sporulate at once;
- `gamma_I`: rate (per day) at which sporulation stops;
- `c_w`: next year's primary inoculum per unit of gated sporulating leaf-days;
- `shed`: rate (per day) at which diseased leaves, sporulating or spent, fall prematurely. Shed
  leaves leave the stand, so the leaf area and the visibly diseased share can fall.
"""
struct MildewParameters
    b_c::Float64
    b_a::Float64
    b_ext::Float64
    tau_s::Float64
    K::Int
    epsilon::Float64
    latent::Float64
    L::Int
    gamma_I::Float64
    c_w::Float64
    shed::Float64
    function MildewParameters(b_c, b_a, b_ext, tau_s, K, epsilon, latent, L, gamma_I, c_w, shed = 0.0)
        for (name, v) in (("b_c", b_c), ("b_a", b_a), ("b_ext", b_ext), ("gamma_I", gamma_I), ("c_w", c_w),
                          ("shed", shed))
            isfinite(v) && v ≥ 0 || throw(ArgumentError("$name must be finite and nonnegative, got $v"))
        end
        tau_s > 0 || throw(ArgumentError("tau_s must be positive (Inf for no maturation), got $tau_s"))
        K ≥ 1 || throw(ArgumentError("K must be at least 1, got $K"))
        0 ≤ epsilon ≤ 1 || throw(ArgumentError("epsilon must be in [0, 1], got $epsilon"))
        L ≥ 0 || throw(ArgumentError("L must be nonnegative, got $L"))
        if L == 0
            latent == 0 || throw(ArgumentError("latent must be 0 when L = 0, got $latent"))
        else
            isfinite(latent) && latent > 0 || throw(ArgumentError("latent must be finite and positive when L > 0, got $latent"))
        end
        return new(b_c, b_a, b_ext, tau_s, K, epsilon, latent, L, gamma_I, c_w, shed)
    end
end

function MildewParameters(; b_c, b_a, b_ext = 0.0, tau_s = Inf, K = 1, epsilon = 0.0, latent = 0.0, L = 0,
                          gamma_I = 0.0, c_w = 1.0, shed = 0.0)
    return MildewParameters(Float64(b_c), Float64(b_a), Float64(b_ext), Float64(tau_s), Int(K), Float64(epsilon),
                            Float64(latent), Int(L), Float64(gamma_I), Float64(c_w), Float64(shed))
end

"""
    FlushProduction(dates, spreads, sizes)

Leaf production `l(t)` from flushes: flush `f` puts out `sizes[f]` leaves (as a fraction of the
stand's maximum) with normally distributed emergence dates, mean `dates[f]` and standard
deviation `spreads[f]` (days), so `l(t) = Σ_f sizes[f] φ((t - dates[f]) / spreads[f]) / spreads[f]`.
The spread stands for differences in flushing between trees.
"""
struct FlushProduction
    dates::Vector{Float64}
    spreads::Vector{Float64}
    sizes::Vector{Float64}
    function FlushProduction(dates, spreads, sizes)
        length(dates) == length(spreads) == length(sizes) ||
            throw(ArgumentError("dates, spreads and sizes must have the same length"))
        all(>(0), spreads) || throw(ArgumentError("spreads must be positive, got $spreads"))
        all(≥(0), sizes) || throw(ArgumentError("sizes must be nonnegative, got $sizes"))
        return new(collect(Float64, dates), collect(Float64, spreads), collect(Float64, sizes))
    end
end

function (f::FlushProduction)(t)
    s = 0.0
    for i in eachindex(f.dates)
        z = (t - f.dates[i]) / f.spreads[i]
        s += f.sizes[i] * exp(-z^2 / 2) / (sqrt(2π) * f.spreads[i])
    end
    return s
end

"""
    ReleaseHazard(; rate, optimum, curvature, intercept = 0.0, wet = t -> 1.0)

The release hazard of primary inoculum, `g(t) = rate ψ(t) wet(t)`, with the readiness of
chasmothecia `ψ(t) = logistic(intercept + curvature (t - optimum)²)`, the logit-quadratic in day
of year fitted by Marçais, Kavkova & Desprez-Loustau (2009); `curvature < 0`. `wet(t)` (between
0 and 1) gates release by rain or free moisture, for instance a [`DailySeries`](@ref).
"""
struct ReleaseHazard{W}
    rate::Float64
    optimum::Float64
    curvature::Float64
    intercept::Float64
    wet::W
end

function ReleaseHazard(; rate, optimum, curvature, intercept = 0.0, wet = t -> 1.0)
    rate ≥ 0 || throw(ArgumentError("rate must be nonnegative, got $rate"))
    curvature < 0 || throw(ArgumentError("curvature must be negative, got $curvature"))
    return ReleaseHazard(Float64(rate), Float64(optimum), Float64(curvature), Float64(intercept), wet)
end

(g::ReleaseHazard)(t) = g.rate * logistic(g.intercept + g.curvature * (t - g.optimum)^2) * g.wet(t)

"""
    DailySeries(values; first_day = 1)

A piecewise-constant function of time from daily values: day `first_day + k - 1` (the interval
`[first_day + k - 1, first_day + k)`) takes `values[k]`. Use it for daily weather.
"""
struct DailySeries
    values::Vector{Float64}
    first_day::Int
end

DailySeries(values; first_day::Integer = 1) = DailySeries(collect(Float64, values), Int(first_day))

function (d::DailySeries)(t)
    k = clamp(floor(Int, t) - d.first_day + 1, 1, length(d.values))
    return d.values[k]
end

"""
    beta_response(T; Tmin = 5.0, Topt = 21.5, Tmax = 30.0)

A temperature response between 0 and 1: zero outside `(Tmin, Tmax)` and 1 at `Topt`, the beta
function of Yin et al. (1995), `((T - Tmin)/(Topt - Tmin)) ((Tmax - T)/(Tmax - Topt))^((Tmax -
Topt)/(Topt - Tmin))`. The defaults follow the germination of *E. alphitoides* conidia reported by
Marçais & Desprez-Loustau (2014): high at 15–25 °C, optimum 20–23 °C, none at 30 °C.
"""
function beta_response(T; Tmin = 5.0, Topt = 21.5, Tmax = 30.0)
    Tmin < Topt < Tmax || throw(ArgumentError("need Tmin < Topt < Tmax, got $Tmin, $Topt, $Tmax"))
    (T ≤ Tmin || T ≥ Tmax) && return 0.0
    return ((T - Tmin) / (Topt - Tmin)) * ((Tmax - T) / (Tmax - Topt))^((Tmax - Topt) / (Topt - Tmin))
end

"""
    degree_day_rate(T; base, reference)

A linear thermal-time development rate, `max(T - base, 0) / (reference - base)`: 0 at or below
`base` and 1 at `reference` (°C). Used as `development` in [`SeasonForcing`](@ref), it makes the
latent period last `latent` days at `reference` and a fixed number of degree-days above `base`.
"""
function degree_day_rate(T; base, reference)
    base < reference || throw(ArgumentError("need base < reference, got $base and $reference"))
    return max(T - base, 0.0) / (reference - base)
end

"""
    FungicideSchedule(dates; efficacy, decay)

Protection by a fungicide applied on `dates` (days): a multiplier of the force of infection,
`1 - efficacy * exp(-(t - t_last) / decay)` with `t_last` the latest application on or before
`t`, and 1 before the first. `efficacy` (between 0 and 1) is the reduction just after an
application and `decay` (days) the e-folding time of its residual activity. Compose it with the
weather, e.g. `weather = t -> beta_response(T(t)) * schedule(t)`, and put the dates in `tstops`.
"""
struct FungicideSchedule
    dates::Vector{Float64}
    efficacy::Float64
    decay::Float64
    function FungicideSchedule(dates::AbstractVector{<:Real}; efficacy::Real, decay::Real)
        0 ≤ efficacy ≤ 1 || throw(ArgumentError("efficacy must be in [0, 1], got $efficacy"))
        isfinite(decay) && decay > 0 || throw(ArgumentError("decay must be finite and positive, got $decay"))
        return new(sort!(collect(Float64, dates)), Float64(efficacy), Float64(decay))
    end
end

function (s::FungicideSchedule)(t)
    k = searchsortedlast(s.dates, t)
    k == 0 && return 1.0
    return 1 - s.efficacy * exp(-(t - s.dates[k]) / s.decay)
end

"""
    SeasonForcing(; tspan, leaves0 = 0.0, production = t -> 0.0, release = t -> 0.0, pulse = false,
                  weather = t -> 1.0, development = t -> 1.0, gate = nothing, inoculations = [],
                  tstops = Float64[])

The forcing of one season, over `tspan = (t0, t1)` in days (day of year, say):
- `leaves0`: young susceptible leaves present at `t0`;
- `production`: leaf production `l(t)`, e.g. a [`FlushProduction`](@ref);
- `release`: the hazard `g(t)` of primary-inoculum release, e.g. a [`ReleaseHazard`](@ref);
  `pulse = true` instead releases all of it at `t0`, infecting a share `1 - exp(-b_a C)` of the
  leaves present (the fast-release limit);
- `weather`: the multiplier `f(t)` of the force of infection, e.g. [`beta_response`](@ref) of a
  daily temperature series;
- `development`: `d(t) ≥ 0`, the rate of latent development relative to the one at which the mean
  latent period is `latent`, e.g. [`degree_day_rate`](@ref) of a daily temperature series; the
  latent stages then advance at `(L/latent) d(t)`, so the latent period is shorter in warm weather;
- `gate`: `η(t)`, the weight of sporulating leaves in next year's inoculum; `nothing` uses the
  sporulating leaves at the end of the season, as Hamelin et al. (2016) do;
- `inoculations`: added inoculum, as `(day, dose)` pairs with `t0 < day < t1`. At each, a share
  `1 - exp(-dose f(day))` of the young susceptible leaves (`1 - exp(-ε dose f(day))` of the mature
  ones) is infected at once, as when spores are blown onto seedlings;
- `tstops`: times the solver must step to, such as day boundaries of daily forcings.
"""
struct SeasonForcing{P,R,F,D,G}
    tspan::Tuple{Float64,Float64}
    leaves0::Float64
    production::P
    release::R
    pulse::Bool
    weather::F
    development::D
    gate::G
    inoculations::Vector{Tuple{Float64,Float64}}
    tstops::Vector{Float64}
end

function SeasonForcing(; tspan, leaves0 = 0.0, production = t -> 0.0, release = t -> 0.0, pulse::Bool = false,
                       weather = t -> 1.0, development = t -> 1.0, gate = nothing,
                       inoculations = Tuple{Float64,Float64}[], tstops = Float64[])
    t0, t1 = Float64.(tspan)
    t0 < t1 || throw(ArgumentError("tspan must be increasing, got $tspan"))
    leaves0 ≥ 0 || throw(ArgumentError("leaves0 must be nonnegative, got $leaves0"))
    inoc = sort!([(Float64(d), Float64(x)) for (d, x) in inoculations])
    for (d, x) in inoc
        t0 < d < t1 || throw(ArgumentError("inoculation day $d must lie strictly inside tspan $tspan"))
        isfinite(x) && x ≥ 0 || throw(ArgumentError("inoculation dose must be finite and nonnegative, got $x"))
    end
    return SeasonForcing((t0, t1), Float64(leaves0), production, release, pulse, weather, development, gate, inoc,
                         sort!(unique!([collect(Float64, tstops); first.(inoc)])))
end

# Apply the inoculations of forcing f that fall on day t to one stand's state u.
function inoculate!(u, p::MildewParameters, f, t)
    for (d, dose) in f.inoculations
        d == t || continue
        x = dose * f.weather(t)
        target = p.L == 0 ? sporulating(p) : first(latent_stages(p))
        for k in young(p)
            infected = -expm1(-x) * u[k]
            u[k] -= infected
            u[target] += infected
        end
        infected = -expm1(-p.epsilon * x) * u[mature(p)]
        u[mature(p)] -= infected
        u[target] += infected
    end
    return u
end

function inoculation_callback(p::MildewParameters, f)
    days = first.(f.inoculations)
    condition(u, t, integrator) = t in days
    affect!(integrator) = inoculate!(integrator.u, p, f, integrator.t)
    return SciMLBase.DiscreteCallback(condition, affect!)
end

# state layout
nstate(p::MildewParameters) = p.K + 1 + p.L + 4
young(p) = 1:(p.K)
mature(p) = p.K + 1
latent_stages(p) = (p.K + 2):(p.K + 1 + p.L)
sporulating(p) = p.K + p.L + 2
spent(p) = p.K + p.L + 3
inoculum(p) = p.K + p.L + 4
gated(p) = p.K + p.L + 5

struct LeafRHS{F}
    p::MildewParameters
    forcing::F
end

(r::LeafRHS)(du, u, _, t) = leaf_rhs!(du, u, r.p, r.forcing, t, u[sporulating(r.p)])

# One stand's right-hand side, with `conidia` the conidial pressure that multiplies b_c: the
# stand's own sporulating leaves in the nonspatial model, a mix of local and dispersed spores in
# the spatial one.
function leaf_rhs!(du, u, p::MildewParameters, f, t, conidia)
    K, L = p.K, p.L
    iI, iC = sporulating(p), inoculum(p)
    g = f.pulse ? 0.0 : f.release(t)
    λ = f.weather(t) * (p.b_c * conidia + p.b_a * g * u[iC] + p.b_ext)
    age = isinf(p.tau_s) ? 0.0 : K / p.tau_s
    susceptible = zero(eltype(u))
    @inbounds for k in 1:K
        inflow = k == 1 ? f.production(t) : age * u[k - 1]
        du[k] = inflow - (λ + age) * u[k]
        susceptible += u[k]
    end
    iM = mature(p)
    @inbounds du[iM] = age * u[K] - p.epsilon * λ * u[iM]
    susceptible += p.epsilon * u[iM]
    new_infections = λ * susceptible
    if L == 0
        du[iI] = new_infections - (p.gamma_I + p.shed) * u[iI]
    else
        rate = L / p.latent * f.development(t)
        first = K + 2
        @inbounds du[first] = new_infections - rate * u[first]
        @inbounds for j in (first + 1):(K + 1 + L)
            du[j] = rate * (u[j - 1] - u[j])
        end
        @inbounds du[iI] = rate * u[K + 1 + L] - (p.gamma_I + p.shed) * u[iI]
    end
    du[spent(p)] = p.gamma_I * u[iI] - p.shed * u[spent(p)]
    du[iC] = -g * u[iC]
    du[gated(p)] = f.gate === nothing ? 0.0 : f.gate(t) * u[iI]
    return nothing
end

"""
    MildewProblem(p::MildewParameters, forcings, C0)

The leaf-level model run over consecutive seasons, one [`SeasonForcing`](@ref) each, starting with
primary inoculum `C0`. Leaves are renewed every season; the inoculum carries over through the
winter map `C = c_w W` (see [`MildewParameters`](@ref)). Solve it with [`solve`](@ref).
"""
struct MildewProblem{F}
    p::MildewParameters
    forcings::Vector{F}
    C0::Float64
    function MildewProblem(p::MildewParameters, forcings::AbstractVector, C0::Real)
        isempty(forcings) && throw(ArgumentError("need at least one season"))
        C0 ≥ 0 || throw(ArgumentError("C0 must be nonnegative, got $C0"))
        v = collect(forcings)
        return new{eltype(v)}(p, v, Float64(C0))
    end
end

MildewProblem(p::MildewParameters, forcing::SeasonForcing, C0::Real) = MildewProblem(p, [forcing], C0)

"""
    MildewSolution

The solution of a [`MildewProblem`](@ref): one dense ODE solution per season. `sol(season, t)` is
the state on day `t` of `season`, laid out as `[U₁ … U_K; M; E₁ … E_L; I; Z; C; W]`. See
[`leaf_area`](@ref), [`incidence`](@ref), [`visible_incidence`](@ref) and
[`winter_inoculum`](@ref).
"""
struct MildewSolution{P,S}
    prob::P
    seasons::Vector{S}
    winter::Vector{Float64}
end

function (sol::MildewSolution)(season::Integer, t::Real)
    1 ≤ season ≤ length(sol.seasons) ||
        throw(ArgumentError("season must be in 1:$(length(sol.seasons)), got $season"))
    t0, t1 = sol.prob.forcings[season].tspan
    t0 ≤ t ≤ t1 || throw(ArgumentError("t must be in [$t0, $t1] for season $season, got $t"))
    return sol.seasons[season](t)
end

"""
    solve(prob::MildewProblem; alg = Tsit5(), abstol = 1e-12, reltol = 1e-10)

Integrate each season of `prob` and apply the winter map between seasons. Returns a
[`MildewSolution`](@ref).
"""
function CommonSolve.solve(prob::MildewProblem; alg = Tsit5(), abstol = 1e-12, reltol = 1e-10)
    p = prob.p
    C = prob.C0
    seasons = []
    winter = Float64[]
    for f in prob.forcings
        u0 = zeros(nstate(p))
        season_start!(u0, p, f, C)
        ode = ODEProblem(LeafRHS(p, f), u0, f.tspan)
        s = if isempty(f.inoculations)
            CommonSolve.solve(ode, alg; abstol, reltol, tstops = f.tstops)
        else
            CommonSolve.solve(ode, alg; abstol, reltol, tstops = f.tstops, callback = inoculation_callback(p, f))
        end
        SciMLBase.successful_retcode(s) || throw(SolverError("season integration failed with retcode $(s.retcode)"))
        push!(seasons, s)
        C = winter_source(p, f, s.u[end])
        push!(winter, C)
    end
    return MildewSolution(prob, [s for s in seasons], winter)
end

# A stand's state at the start of a season with primary inoculum C (u0 must be zeroed).
function season_start!(u0, p::MildewParameters, f, C)
    u0[1] = f.leaves0
    if f.pulse
        # release everything at t0: each leaf present escapes with probability exp(-b_a C)
        infected = -expm1(-p.b_a * C * f.weather(f.tspan[1]))
        u0[1] = f.leaves0 * (1 - infected)
        u0[p.L == 0 ? sporulating(p) : first(latent_stages(p))] = f.leaves0 * infected
    else
        u0[inoculum(p)] = C
    end
    return u0
end

# The winter map: next season's primary inoculum from a stand's state at the end of this one.
winter_source(p::MildewParameters, f, uend) = p.c_w * (f.gate === nothing ? uend[sporulating(p)] : uend[gated(p)])

"""
    leaf_area(sol::MildewSolution, season, t)

The leaves present on day `t` of `season` (young, mature, latent, sporulating and spent), as a
fraction of the stand's maximum.
"""
leaf_area(sol::MildewSolution, season::Integer, t::Real) = stand_leaves(sol.prob.p, sol(season, t))

# per-stand summaries of a state vector u
stand_leaves(p::MildewParameters, u) = sum(@view u[1:spent(p)])
stand_infected(p::MildewParameters, u) = sum(@view u[(p.K + 2):spent(p)])
function stand_visible(p::MildewParameters, u, latent_visible)
    return u[sporulating(p)] + u[spent(p)] + latent_visible * sum(@view u[latent_stages(p)])
end
function stand_incidence(p::MildewParameters, u)
    A = stand_leaves(p, u)
    return A > 0 ? stand_infected(p, u) / A : 0.0
end
function stand_visible_incidence(p::MildewParameters, u, latent_visible)
    A = stand_leaves(p, u)
    return A > 0 ? stand_visible(p, u, latent_visible) / A : 0.0
end
function check_latent_visible(latent_visible)
    0 ≤ latent_visible ≤ 1 || throw(ArgumentError("latent_visible must be in [0, 1], got $latent_visible"))
    return nothing
end

"""
    incidence(sol::MildewSolution, season, t)

The share of leaves present on day `t` of `season` that are infected (latent, sporulating or
spent).
"""
incidence(sol::MildewSolution, season::Integer, t::Real) = stand_incidence(sol.prob.p, sol(season, t))

"""
    visible_incidence(sol::MildewSolution, season, t; latent_visible = 0.0)

The share of leaves present on day `t` of `season` with visible mildew: sporulating and spent
leaves, plus a share `latent_visible` of latent ones.
"""
function visible_incidence(sol::MildewSolution, season::Integer, t::Real; latent_visible::Real = 0.0)
    check_latent_visible(latent_visible)
    return stand_visible_incidence(sol.prob.p, sol(season, t), latent_visible)
end

"""
    winter_inoculum(sol::MildewSolution)

The primary inoculum carried into the season after each season, `c_w W`.
"""
winter_inoculum(sol::MildewSolution) = copy(sol.winter)
