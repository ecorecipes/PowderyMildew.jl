# Fit the leaf-level model to the six dated mildew scores of the Pierroton seedling experiment
# (Barrès et al. 2024; dates, treatments and flush timing from Desprez-Loustau et al. 2014).
#
#   julia --project=vignettes scripts/fit_pierroton.jl [outdir]
#
# Data: data/pierroton/ (plot means, budburst, daily weather; see its README). Each of the nine
# plots is a replicate of its treatment: L (fungicide), M (natural), H (natural plus spores
# blown on). The model is run once per treatment for 2009-2011; the plot means of the six dated
# scores (first flush and final, each year) are fitted on the logit scale with normal error.
# Writes outdir/fit.txt and outdir/fit.png (default runs/pierroton_fit, git-ignored).

ENV["GKSwstype"] = "100"   # headless GR: no plot-window helper holding the output open
using PowderyMildew
using Optim
using Plots
using Printf
using Statistics
using Dates: Date, dayofyear
using LinearAlgebra: Symmetric, inv, diag, eigvals

const DATA = joinpath(@__DIR__, "..", "data", "pierroton")
const OUT = length(ARGS) ≥ 1 ? ARGS[1] : joinpath(@__DIR__, "..", "runs", "pierroton_fit")
mkpath(OUT)

function readcsv(path)
    lines = readlines(path)
    header = split(lines[1], ',')
    return [Dict(zip(header, split(l, ','))) for l in lines[2:end]]
end

doy(y, m, d) = dayofyear(Date(y, m, d))

# ---- data ---------------------------------------------------------------------------------
const YEARS = (2009, 2010, 2011)
scores = filter(r -> r["date_known"] == "1" && r["kind"] == "percent", readcsv(joinpath(DATA, "plot_scores.csv")))
obs = [(year = parse(Int, r["year"]), day = dayofyear(Date(r["date"])), treatment = r["treatment"],
        plot = parse(Int, r["plot"]), mean = parse(Float64, r["mean"]) / 100, n = parse(Int, r["n"])) for r in scores]
@assert length(obs) == 54
weather_rows = readcsv(joinpath(DATA, "weather.csv"))
temperature = Dict(y => DailySeries([parse(Float64, r["tmean"]) for r in weather_rows if r["year"] == string(y)])
                   for y in YEARS)

# ---- forcing (published dates) ------------------------------------------------------------
# Flushes: 2009 seedling emergence (sown 21-27 April, 50 % emerged after a month); 2010 and 2011
# first flush when the mean budburst score crosses 2 (weekly scores from week 15, taken as
# Wednesdays), second flushes noted 15 June 2010 and 24 May 2011, third 20 June 2011.
phen = readcsv(joinpath(DATA, "plot_phenology.csv"))
function score2_day(year)
    orders = sort(unique(parse(Int, r["order"]) for r in phen if r["year"] == string(year)))
    m = [mean(parse(Float64, r["mean"]) for r in phen if r["year"] == string(year) && r["order"] == string(k))
         for k in orders]
    k = findfirst(≥(2), m)
    frac = (2 - m[k - 1]) / (m[k] - m[k - 1])
    wednesday15 = year == 2010 ? doy(2010, 4, 14) : doy(2011, 4, 13)
    return wednesday15 + 7 * (orders[k - 1] - 1 + frac)
end
const FLUSH = Dict(2009 => FlushProduction([145.0, 190.0], [10.0, 10.0], [0.6, 0.4]),
                   2010 => FlushProduction([score2_day(2010), 160.0, 200.0], [5.0, 8.0, 8.0], [0.5, 0.3, 0.2]),
                   2011 => FlushProduction([score2_day(2011), 140.0, 171.0], [5.0, 6.0, 8.0], [0.5, 0.3, 0.2]))
# Fungicide (L): first application and the weeks of all applications.
fungicide_days(first, w0, weeks) = [first + 7.0 * (w - w0) for w in weeks]
const FUNGICIDE = Dict(2009 => fungicide_days(doy(2009, 6, 22), 25, (25, 28, 30, 32, 35, 37)),
                       2010 => fungicide_days(doy(2010, 6, 1), 22, (22, 27, 30, 31, 34, 39)),
                       2011 => fungicide_days(doy(2011, 5, 27), 21, (21, 23, 25, 27, 30, 31, 33, 35, 39)))
# Added inoculum (H).
const INOCULATION = Dict(2009 => [doy(2009, 6, 19)], 2010 => [doy(2010, 5, 11), doy(2010, 6, 24)],
                         2011 => [doy(2011, 4, 20), doy(2011, 4, 28)])
# Latent development in degree-days above LATENT_BASE (°C, an assumption), with the 7-day latent
# period reached at 20 °C: about 10 days at 15 °C.
const LATENT_BASE = 5.0
const TSPAN = (90.0, 300.0)
const DAYS = collect(91.0:1.0:299.0)

# Ascospore release peaks a fitted lag after each year's budburst (score 2; 2009, emergence), so
# that it follows the same spring warmth as the host. The 2009 release has no source (C0 = 0).
const BUDBURST = Dict(y => FLUSH[y].dates[1] for y in YEARS)
# Fungicide efficacy is fixed at 1 (full protection just after a spray): fitted, it went to the
# bound, with or without a free decay (wiki: oak/barres-fit.md). The decay is fitted.
const NAMES = ["b_c", "b_ext", "c_w", "dose", "tau_s", "release_lag", "fungicide_decay", "obs_sd"]
function unpack(θ)
    return (b_c = exp(θ[1]), b_ext = exp(θ[2]), c_w = exp(θ[3]), efficacy = 1.0, dose = exp(θ[4]),
            tau_s = exp(θ[5]), lag = θ[6], decay = exp(θ[7]), sd = exp(θ[8]))
end

function simulate(q, treatment)
    p = MildewParameters(; b_c = q.b_c, b_a = 1.0, b_ext = q.b_ext, tau_s = q.tau_s, K = 3, epsilon = 0.05,
                         latent = 7.0, L = 3, gamma_I = 0.05, c_w = q.c_w)
    forcings = map(YEARS) do y
        T = temperature[y]
        release = ReleaseHazard(; rate = 0.2, optimum = BUDBURST[y] + q.lag, curvature = -0.003)
        spray = treatment == "L" ? FungicideSchedule(FUNGICIDE[y]; efficacy = q.efficacy, decay = q.decay) : t -> 1.0
        inoc = treatment == "H" ? [(d, q.dose) for d in INOCULATION[y]] : Tuple{Float64,Float64}[]
        SeasonForcing(; tspan = TSPAN, production = FLUSH[y], release, weather = t -> beta_response(T(t)) * spray(t),
                      development = t -> degree_day_rate(T(t); base = LATENT_BASE, reference = 20.0),
                      # myclobutanil also suppresses chasmothecium formation on treated leaves
                      gate = t -> t ≥ 240 ? spray(t) : 0.0, inoculations = inoc,
                      tstops = [DAYS; treatment == "L" ? FUNGICIDE[y] : Float64[]])
    end
    return solve(MildewProblem(p, collect(forcings), 0.0); abstol = 1e-9, reltol = 1e-7)
end

logit(x) = log(x / (1 - x))
# the score is the percentage of leaf area infected, i.e. the visible share of leaf area
severity(sol, season, day) = visible_incidence(sol, season, day)

function predictions(q)
    sols = Dict(t => simulate(q, t) for t in ("L", "M", "H"))
    return [severity(sols[o.treatment], o.year - 2008, Float64(o.day)) for o in obs], sols
end

# weak priors keep the optimiser in a sensible region; release near budburst (Marçais et al. 2009)
# the fungicide's residual activity: median 10 days, log-sd 0.7 (95 % within about 2.5 to 40 days)
logprior(θ, q) = -0.5 * (q.lag / 15)^2 - 0.5 * ((θ[7] - log(10.0)) / 0.7)^2 -
                 0.5 * sum(abs2, (θ[1:5] .- log.([0.2, 0.005, 1.0, 1.0, 20.0])) ./ 3)

function negloglik(θ)
    q = unpack(θ)
    pred = try
        first(predictions(q))
    catch err
        err isa PowderyMildewError || rethrow()
        return Inf
    end
    r = [logit(clamp(o.mean, 1e-4, 1 - 1e-4)) - logit(clamp(m, 1e-4, 1 - 1e-4)) for (o, m) in zip(obs, pred)]
    return length(r) * log(q.sd) + sum(abs2, r) / (2q.sd^2) - logprior(θ, q)
end

println("budburst (score 2) days: ", join((@sprintf("%d %.1f", y, BUDBURST[y]) for y in YEARS), ", "))
# start from the fit with efficacy and decay both free (runs/pierroton_fit5; efficacy 1 - 3e-7)
θ0 = [-0.05054, -6.19660, -1.75241, -1.29904, 2.23563, -12.19812, 2.42790, 0.07621]
@printf("start: %.3f\n", negloglik(θ0))
# Nelder-Mead from three starts, each restarted from its own optimum until it stops improving
function minimise(start)
    res = optimize(negloglik, start, NelderMead(), Optim.Options(; iterations = 3000, f_reltol = 1e-10))
    for _ in 1:6
        again = optimize(negloglik, Optim.minimizer(res), NelderMead(), Optim.Options(; iterations = 3000, f_reltol = 1e-10))
        improved = Optim.minimum(res) - Optim.minimum(again)
        res = again
        improved < 1e-6 && break
    end
    return res
end
best = nothing
for (k, start) in enumerate((θ0,))
    r = minimise(start)
    @printf("start %d: objective %.4f (converged: %s)\n", k, Optim.minimum(r), Optim.converged(r))
    global best = best === nothing || Optim.minimum(r) < Optim.minimum(best) ? r : best
end
res = best
θ̂ = Optim.minimizer(res)
q̂ = unpack(θ̂)

# Laplace approximation: central-difference Hessian on the unconstrained scale
function hessian(f, x; h = 1e-2)
    n = length(x)
    H = zeros(n, n)
    for i in 1:n, j in i:n
        e_i, e_j = (1:n .== i) .* h, (1:n .== j) .* h
        H[i, j] = H[j, i] = (f(x + e_i + e_j) - f(x + e_i - e_j) - f(x - e_i + e_j) + f(x - e_i - e_j)) / (4h^2)
    end
    return H
end
H = hessian(negloglik, θ̂)
Σ = inv(Symmetric(H))
se = sqrt.(max.(diag(Σ), 0))

pred, sols = predictions(q̂)
open(joinpath(OUT, "fit.txt"), "w") do io
    for out in (stdout, io)
        @printf(out, "objective %.4f; Hessian eigenvalues %.3g to %.3g\n", Optim.minimum(res), extrema(eigvals(Symmetric(H)))...)
        println(out, "optimum on the unconstrained scale: ", θ̂)
        println(out, "parameter          estimate   (95% interval from the Laplace approximation)")
        tr = [exp, exp, exp, exp, exp, identity, exp, exp]
        for (k, name) in enumerate(NAMES)
            @printf(out, "%-18s %9.4g   (%.4g, %.4g)\n", name, tr[k](θ̂[k]), tr[k](θ̂[k] - 1.96se[k]), tr[k](θ̂[k] + 1.96se[k]))
        end
        println(out, "\nyear day trt  observed (plot means, %)        model (%)")
        for y in YEARS, d in sort(unique(o.day for o in obs if o.year == y)), t in ("L", "M", "H")
            idx = findall(o -> o.year == y && o.day == d && o.treatment == t, obs)
            @printf(out, "%d %3d  %s   %-30s %6.1f\n", y, d, t, join((@sprintf("%5.1f", 100obs[i].mean) for i in idx), " "),
                    100pred[idx[1]])
        end
    end
end

colors = Dict("L" => 1, "M" => 2, "H" => 3)
panels = map(YEARS) do y
    plt = plot(; title = string(y), xlabel = "day of year", ylabel = y == 2009 ? "leaf area infected (%)" : "",
               legend = y == 2009 ? :topleft : false, ylims = (0, 75))
    for t in ("L", "M", "H")
        plot!(plt, DAYS, [100severity(sols[t], y - 2008, d) for d in DAYS]; lw = 2, color = colors[t],
              label = Dict("L" => "fungicide", "M" => "natural", "H" => "added spores")[t])
        o = filter(o -> o.year == y && o.treatment == t, obs)
        scatter!(plt, [x.day for x in o], [100x.mean for x in o]; color = colors[t], ms = 4, label = nothing)
    end
    plt
end
savefig(plot(panels...; layout = (1, 3), size = (1100, 360), left_margin = 4Plots.mm, bottom_margin = 5Plots.mm),
        joinpath(OUT, "fit.png"))
println("wrote ", OUT)
