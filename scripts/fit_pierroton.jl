# Fit the leaf-level model to the mildew scores of the Pierroton seedling experiment
# (Barrès et al. 2024; dates, treatments and flush timing from Desprez-Loustau et al. 2014).
#
#   julia --project=vignettes scripts/fit_pierroton.jl [outdir] [variant]
#
# variant (default "base") adds options joined by "+", e.g. "winter+dates":
#   winter  a separate survival factor k for the 2010/11 winter (the 2011 primary inoculum is
#           k c_w W instead of c_w W), for the early 2011 epidemic;
#   dates   also fits the six undated 2009-2011 scores (two a year), each pair at the ordered
#           whole days between that year's dated scores that fit best (profiled out; every
#           treatment scored on the same day; oid3_10 and oid4_11, which have second-shoot
#           scores, no earlier than that flush). They are plotted as diamonds;
#   shed    diseased leaves (sporulating or spent) fall at a fitted rate (`shed` in
#           MildewParameters), so the visibly diseased share can fall;
#   severity  diseased leaves fall at a rate proportional to the visibly diseased share
#           (`shed_severity`), so shedding is heaviest in severe epidemics;
#   drought diseased leaves fall at a rate proportional to a soil-water deficit between 0 and 1
#           (`shed` times the `shedding` forcing): a daily bucket of CAPACITY mm, filled by
#           rain and emptied by Hargreaves reference evapotranspiration scaled by the water
#           left, run from 1 January 2008. Not with "shed"; with "severity" it scales both.
#   power   with "severity", the rate is shed_severity D^p with p fitted (`shed_power`), so mild
#           epidemics shed little;
#   late    two leaf cohorts, the first flush and the later ones (`cohorts`), and only later-flush
#           leaves are shed (`cohort_shedding = [0, 1]`): mildew does not make first-flush leaves
#           fall (Desprez-Loustau et al. 2014). Each cohort's shedding uses its own D.
#   shoots  with "dates", the same two cohorts, and the per-shoot scores fitted to them on the day of
#           their whole-seedling companion: oid3P1_10 (first shoot) to cohort 1, oid3P2_10 and
#           oid4P2_11 (second shoot) to cohort 2. In 2010 the whole-seedling score oid3_10 is the
#           seedlings' mean of their shoot scores (P1 + f (P2 - P1)/2, f the share with a second
#           shoot, matches all nine plots to 0.5 points), not a leaf-area share, so it is replaced
#           by the shoot scores; oid4_11 in 2011 is not such a mean and stays. Not with a shedding
#           option unless "late".
#   sizes   the later flushes' share of the leaves fitted for 2010 and 2011 (logit scale, prior
#           median 0.5, log-odds sd 1.5), split 3:2 between the second and third flush; the
#           default shares are 0.5 (sizes 0.5, 0.3, 0.2), which were assumed, not measured;
#   mean2010  with two cohorts, every 2010 score is the seedlings' mean of their shoot scores,
#           V1 + f (V2 - V1)/2, with V the cohorts' visible shares and f the share of the later
#           flushes' leaves out (for the share of seedlings with a second shoot), as oid3_10 is.
#   epsilon the susceptibility of mature leaves relative to young ones fitted (logit scale, prior
#           median 0.05, the value fixed otherwise, log-odds sd 1.5);
#   A start for the added parameters (log scale) may follow the variant on the command line; the
#   values replace the last entries of the start vector, so a full vector replaces all of it.
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
using Dates: Date, Day, dayofyear
using LinearAlgebra: Symmetric, inv, diag, eigvals

const DATA = joinpath(@__DIR__, "..", "data", "pierroton")
const OUT = length(ARGS) ≥ 1 ? ARGS[1] : joinpath(@__DIR__, "..", "runs", "pierroton_fit")
mkpath(OUT)
const VARIANT = length(ARGS) ≥ 2 ? ARGS[2] : "base"
const OPTIONS = Set(split(VARIANT, '+'))
issubset(OPTIONS, Set(["base", "winter", "dates", "shed", "severity", "drought", "power", "late", "shoots", "sizes", "mean2010", "epsilon"])) ||
    error("unknown variant: $VARIANT")
const WINTER, UNDATED, SHED = "winter" in OPTIONS, "dates" in OPTIONS, "shed" in OPTIONS
const SEVERITY, DROUGHT = "severity" in OPTIONS, "drought" in OPTIONS
const POWER, LATE = "power" in OPTIONS, "late" in OPTIONS
SHED && DROUGHT && error("\"shed\" and \"drought\" both set the shedding rate; use one")
POWER && !SEVERITY && error("\"power\" needs \"severity\"")
LATE && !(SHED || SEVERITY || DROUGHT) && error("\"late\" needs a shedding option")
const SHOOTS = "shoots" in OPTIONS
SHOOTS && !UNDATED && error("\"shoots\" needs \"dates\"")
# two cohorts with both shed would compute D per cohort, a different model from one cohort
SHOOTS && (SHED || SEVERITY || DROUGHT) && !LATE && error("\"shoots\" with shedding needs \"late\"")
const TWO = LATE || SHOOTS     # two leaf cohorts: the first flush and the later ones
const SIZES, MEAN2010 = "sizes" in OPTIONS, "mean2010" in OPTIONS
MEAN2010 && !TWO && error("\"mean2010\" needs two cohorts (\"shoots\" or \"late\")")
const ISHED = 9 + WINTER       # positions in θ of log(shed), log(shed_severity), log(drought rate)
const ISEV = ISHED + SHED
const IDRY = ISEV + SEVERITY
const IPOW = IDRY + DROUGHT
const ISIZE = IPOW + POWER     # logit of the later flushes' share, 2010 then 2011
const EPSILON = "epsilon" in OPTIONS
const IEPS = ISIZE + 2SIZES     # logit of the mature leaves' relative susceptibility
logistic(x) = 1 / (1 + exp(-x))

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
# the undated whole-seedling scores, two a year between the dated first-flush and final scores
const UNDATED_COLUMNS = Dict(2009 => ("oid2_09", "oid3_09"), 2010 => ("oid3_10", "oid4_10"),
                             2011 => ("oid3_11", "oid4_11"))
const WINDOW = Dict(y => extrema(o.day for o in obs if o.year == y) for y in YEARS)
# scores with a second-shoot companion (oid3P2_10, oid4P2_11) come after that flush was noted
const EARLIEST = Dict("oid3_10" => doy(2010, 6, 15), "oid4_11" => doy(2011, 5, 24))
const undated = [(year = parse(Int, r["year"]), column = String(r["column"]), treatment = r["treatment"],
                  plot = parse(Int, r["plot"]), mean = parse(Float64, r["mean"]) / 100)
                 for r in readcsv(joinpath(DATA, "plot_scores.csv"))
                 if r["column"] in Iterators.flatten(values(UNDATED_COLUMNS))]
@assert length(undated) == 54
# per-shoot scores: column => (whole-seedling companion, cohort)
const SHOOT_COLUMNS = Dict("oid3P1_10" => ("oid3_10", 1), "oid3P2_10" => ("oid3_10", 2), "oid4P2_11" => ("oid4_11", 2))
const REPLACED = SHOOTS ? Set(["oid3_10"]) : Set{String}()
const shoot = [(year = parse(Int, r["year"]), column = String(r["column"]), paired = SHOOT_COLUMNS[r["column"]][1],
                cohort = SHOOT_COLUMNS[r["column"]][2], treatment = r["treatment"], plot = parse(Int, r["plot"]),
                mean = parse(Float64, r["mean"]) / 100)
               for r in readcsv(joinpath(DATA, "plot_scores.csv")) if SHOOTS && haskey(SHOOT_COLUMNS, r["column"])]
@assert length(shoot) == (SHOOTS ? 27 : 0)
weather_rows = readcsv(joinpath(DATA, "weather.csv"))
temperature = Dict(y => DailySeries([parse(Float64, r["tmean"]) for r in weather_rows if r["year"] == string(y)])
                   for y in YEARS)

# Soil-water deficit for "drought": 1 - S/CAPACITY, with S' = rain - ET0 S/CAPACITY and ET0 from
# Hargreaves (FAO-56 extraterrestrial radiation at the site's latitude). The capacity is an
# assumption (sandy Landes soil, young seedlings).
const CAPACITY = 100.0
function extraterrestrial(J; latitude = deg2rad(44.747))   # mm/day equivalent
    dr, dec = 1 + 0.033cos(2π * J / 365), 0.409sin(2π * J / 365 - 1.39)
    ωs = acos(-tan(latitude) * tan(dec))
    return 0.408 * 24 * 60 / π * 0.0820 * dr * (ωs * sin(latitude) * sin(dec) + cos(latitude) * cos(dec) * sin(ωs))
end
const DEFICIT = let S = CAPACITY, out = Dict(y => Float64[] for y in YEARS)
    for r in weather_rows
        y, J = parse(Int, r["year"]), parse(Int, r["doy"])
        T, tn, tx, rain = (parse(Float64, r[k]) for k in ("tmean", "tmin", "tmax", "precip"))
        et0 = 0.0023 * extraterrestrial(J) * (T + 17.8) * sqrt(max(tx - tn, 0.0))
        S = clamp(S + rain - et0 * S / CAPACITY, 0.0, CAPACITY)
        y in YEARS && push!(out[y], 1 - S / CAPACITY)
    end
    Dict(y => DailySeries(out[y]) for y in YEARS)
end

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
# for two cohorts: the first flush (2009, emergence) and the later flushes
cohorts_of(f) = (FlushProduction(f.dates[1:1], f.spreads[1:1], f.sizes[1:1]),
                 FlushProduction(f.dates[2:end], f.spreads[2:end], f.sizes[2:end]))
const COHORTS = Dict(y => cohorts_of(f) for (y, f) in FLUSH)
# the year's flushes, with "sizes" the fitted share s of the later flushes (sizes 1 - s, 0.6s, 0.4s)
function flushes(q, y)
    f = FLUSH[y]
    haskey(q.later, y) || return f
    s = q.later[y]
    return FlushProduction(f.dates, f.spreads, [1 - s, 0.6s, 0.4s])
end
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
const NAMES = ["b_c", "b_ext", "c_w", "dose", "tau_s", "release_lag", "fungicide_decay", "obs_sd",
               (WINTER ? ["winter_2011"] : String[])..., (SHED ? ["shed"] : String[])...,
               (SEVERITY ? ["shed_severity"] : String[])..., (DROUGHT ? ["shed_drought"] : String[])...,
               (POWER ? ["shed_power"] : String[])..., (SIZES ? ["later_2010", "later_2011"] : String[])...,
               (EPSILON ? ["epsilon"] : String[])...]
function unpack(θ)
    return (b_c = exp(θ[1]), b_ext = exp(θ[2]), c_w = exp(θ[3]), efficacy = 1.0, dose = exp(θ[4]),
            tau_s = exp(θ[5]), lag = θ[6], decay = exp(θ[7]), sd = exp(θ[8]),
            winter = WINTER ? exp(θ[9]) : 1.0,
            shed = SHED ? exp(θ[ISHED]) : DROUGHT ? exp(θ[IDRY]) : 0.0,
            shed_severity = SEVERITY ? exp(θ[ISEV]) : 0.0, shed_power = POWER ? exp(θ[IPOW]) : 1.0,
            later = SIZES ? Dict(2010 => logistic(θ[ISIZE]), 2011 => logistic(θ[ISIZE + 1])) : Dict{Int,Float64}(),
            epsilon = EPSILON ? logistic(θ[IEPS]) : 0.05)
end

function simulate(q, treatment)
    p = MildewParameters(; b_c = q.b_c, b_a = 1.0, b_ext = q.b_ext, tau_s = q.tau_s, K = 3, epsilon = q.epsilon,
                         latent = 7.0, L = 3, gamma_I = 0.05, c_w = q.c_w, shed = q.shed,
                         shed_severity = q.shed_severity, shed_power = q.shed_power, cohorts = TWO ? 2 : 1,
                         cohort_shedding = LATE ? [0.0, 1.0] : TWO ? [1.0, 1.0] : [1.0])
    forcings = map(YEARS) do y
        T = temperature[y]
        release = ReleaseHazard(; rate = 0.2, optimum = BUDBURST[y] + q.lag, curvature = -0.003)
        spray = treatment == "L" ? FungicideSchedule(FUNGICIDE[y]; efficacy = q.efficacy, decay = q.decay) : t -> 1.0
        inoc = treatment == "H" ? [(d, q.dose) for d in INOCULATION[y]] : Tuple{Float64,Float64}[]
        SeasonForcing(; tspan = TSPAN, production = TWO ? cohorts_of(flushes(q, y)) : flushes(q, y), release, weather = t -> beta_response(T(t)) * spray(t),
                      development = t -> degree_day_rate(T(t); base = LATENT_BASE, reference = 20.0),
                      shedding = DROUGHT ? DEFICIT[y] : t -> 1.0,
                      # myclobutanil also suppresses chasmothecium formation on treated leaves;
                      # the 2010 gate carries the 2010/11 winter factor (1 unless "winter")
                      gate = t -> t ≥ 240 ? (y == 2010 ? q.winter : 1.0) * spray(t) : 0.0, inoculations = inoc,
                      tstops = [DAYS; treatment == "L" ? FUNGICIDE[y] : Float64[]])
    end
    return solve(MildewProblem(p, collect(forcings), 0.0); abstol = 1e-9, reltol = 1e-7)
end

logit(x) = log(x / (1 - x))
# the score is the percentage of leaf area infected, i.e. the visible share of leaf area
# (with "mean2010", in 2010 the seedlings' mean of their shoot scores; f, the share of the later
# flushes' leaves out, uses the leaves present, so it assumes no shedding)
function severity(sol, season, day, q)
    (MEAN2010 && season == 2) || return visible_incidence(sol, season, day)
    f = min(leaf_area(sol, season, day; cohort = 2) / sum(flushes(q, 2010).sizes[2:end]), 1.0)
    v1, v2 = (visible_incidence(sol, season, day; cohort = k) for k in 1:2)
    return v1 + f * (v2 - v1) / 2
end

lg(x) = logit(clamp(x, 1e-4, 1 - 1e-4))

# Model predictions and logit residuals. With "dates", each year's pair of undated scores is put
# at the ordered pair of whole days strictly between the dated scores that minimises their
# squared residuals: the profile over the dates, since nothing else in the objective depends on
# them.
function fitted(q)
    sols = Dict(t => simulate(q, t) for t in ("L", "M", "H"))
    pred = [severity(sols[o.treatment], o.year - 2008, Float64(o.day), q) for o in obs]
    r = [lg(o.mean) - lg(m) for (o, m) in zip(obs, pred)]
    dates = Dict{String,Int}()
    if UNDATED
        for y in YEARS
            days = (WINDOW[y][1] + 1):(WINDOW[y][2] - 1)
            curve = Dict(t => [lg(severity(sols[t], y - 2008, Float64(d), q)) for d in days] for t in ("L", "M", "H"))
            c1, c2 = UNDATED_COLUMNS[y]
            # with "shoots", each cohort's visible share; the shoot scores join their companion's day
            ys = [o for o in shoot if o.year == y]
            cc = Dict((t, k) => [lg(visible_incidence(sols[t], y - 2008, Float64(d); cohort = k)) for d in days]
                      for t in ("L", "M", "H"), k in unique(o.cohort for o in ys))
            ssr(c) = [sum((abs2(lg(o.mean) - curve[o.treatment][i]) for o in undated
                           if o.column == c && !(c in REPLACED)); init = 0.0) +
                      sum((abs2(lg(o.mean) - cc[(o.treatment, o.cohort)][i]) for o in ys if o.paired == c); init = 0.0)
                      for i in eachindex(days)]
            s1, s2 = ssr(c1), ssr(c2)
            for (s, c) in ((s1, c1), (s2, c2)), i in eachindex(days)
                days[i] < get(EARLIEST, c, 0) && (s[i] = Inf)
            end
            best, i1, i2, low, at = Inf, 0, 0, Inf, 0
            for j in 2:length(days)        # best i1 < i2: running minimum of s1 before j
                s1[j - 1] < low && ((low, at) = (s1[j - 1], j - 1))
                low + s2[j] < best && ((best, i1, i2) = (low + s2[j], at, j))
            end
            dates[c1], dates[c2] = days[i1], days[i2]
            for o in undated
                o.year == y && !(o.column in REPLACED) && push!(r, lg(o.mean) - curve[o.treatment][o.column == c1 ? i1 : i2])
            end
            for o in ys
                push!(r, lg(o.mean) - cc[(o.treatment, o.cohort)][o.paired == c1 ? i1 : i2])
            end
        end
    end
    return (; pred, r, dates, sols)
end

# weak priors keep the optimiser in a sensible region; release near budburst (Marçais et al. 2009)
# the fungicide's residual activity: median 10 days, log-sd 0.7 (95 % within about 2.5 to 40 days);
# the 2010/11 winter factor: median 1, log-sd 2; shedding: median 0.002 per day (a diseased
# leaf's mean life before premature fall 500 days, so most survive the season), log-sd 2;
# shed_severity: median 0.005 per day per unit diseased share, log-sd 2; the drought rate:
# median 0.003 per day at full deficit, log-sd 2 (both about 0.002 at a typical summer value)
logprior(θ, q) = -0.5 * (q.lag / 15)^2 - 0.5 * ((θ[7] - log(10.0)) / 0.7)^2 -
                 0.5 * sum(abs2, (θ[1:5] .- log.([0.2, 0.005, 1.0, 1.0, 20.0])) ./ 3) -
                 (WINTER ? 0.5 * (θ[9] / 2)^2 : 0.0) - (SHED ? 0.5 * ((θ[ISHED] - log(0.002)) / 2)^2 : 0.0) -
                 (SEVERITY ? 0.5 * ((θ[ISEV] - log(0.005)) / 2)^2 : 0.0) -
                 (DROUGHT ? 0.5 * ((θ[IDRY] - log(0.003)) / 2)^2 : 0.0) - (POWER ? 0.5 * θ[IPOW]^2 : 0.0) -
                 (SIZES ? 0.5 * sum(abs2, θ[ISIZE:(ISIZE + 1)] ./ 1.5) : 0.0) -
                 (EPSILON ? 0.5 * ((θ[IEPS] - logit(0.05)) / 1.5)^2 : 0.0)

function negloglik(θ)
    q = unpack(θ)
    r = try
        fitted(q).r
    catch err
        err isa PowderyMildewError || rethrow()
        return Inf
    end
    return length(r) * log(q.sd) + sum(abs2, r) / (2q.sd^2) - logprior(θ, q)
end

println("budburst (score 2) days: ", join((@sprintf("%d %.1f", y, BUDBURST[y]) for y in YEARS), ", "))
# start from the base fit (runs/pierroton_fit6, objective 31.8754), or with "dates" from the
# fit-8 optimum of the same variant (runs/pierroton_fit8_*); "winter" starts at k = 3 and
# "shed" at 0.001 per day
θ0 = if UNDATED && WINTER
    [-0.3198769484862051, -6.262513650703709, -2.7638127948248825, -1.6610414898299364, 2.6194371379575245,
     -7.344399879501469, 2.6475266204658663, -0.11983558463385383, 1.9212111431795236]
elseif UNDATED
    [-0.05027322394419971, -6.241311132277026, -1.6766906990844392, -1.467156260990763, 2.245497117133406,
     -14.028218466298291, 2.235622150388756, -0.08804940452979772]
else
    [-0.050825674964122325, -6.195781043528346, -1.7499149460273897, -1.2983793059662951,
     2.2350727981292726, -12.223854149725843, 2.4289737677520136, 0.0750593311357487,
     (WINTER ? [log(3.0)] : Float64[])...]
end
if SHOOTS && UNDATED && WINTER && !(SHED || SEVERITY || DROUGHT)   # the fit-12 optimum (pierroton_fit12_shoots)
    θ0 = [0.1350585576810398, -6.556696462460295, -3.4911550441787256, -2.136159349058975, 2.1734619326934452,
          -9.326002792648845, 2.390308990529697, 0.07156404061388498, 1.899343854954856]
end
SHED && push!(θ0, log(0.001))
if SEVERITY && LATE && SHOOTS && UNDATED && WINTER && !SHED && !DROUGHT   # fit-11 late optimum
    θ0 = [0.08873906110075355, -6.359842842405908, -2.997394889880566, -1.2616773924224822, 2.445459271591088,
          -7.383578925333428, 2.738274340389261, -0.14535722188626152, 1.7094265593500304, -2.7862709062711155]
elseif SEVERITY && UNDATED && WINTER && !SHED && !DROUGHT   # the fit-10 optimum (runs/pierroton_fit10_*)
    θ0 = [0.31759195163279214, -6.591071744648859, -3.04487710692643, -1.072858612183583, 2.588070829008947,
          -7.725750928580431, 3.1560008605176697, -0.16682358634763997, 1.5342647908926692, -2.9567556372800534]
else
    SEVERITY && push!(θ0, log(0.002))
end
DROUGHT && push!(θ0, log(0.001))
POWER && push!(θ0, 0.0)
SIZES && push!(θ0, 0.0, 0.0)
if SIZES && SHOOTS && UNDATED && WINTER && !(SHED || SEVERITY || DROUGHT || MEAN2010)   # fit-13 sizes optimum
    θ0[1:11] = [-0.047388940120668376, -6.62623463164331, -3.3375413521401125, -2.504881154748258, 2.4567268288039066,
                -9.537665012058062, 2.461783358221106, 0.054063520809400925, 2.192505837616085, 1.0057833577270743,
                1.1583217813912703]
end
EPSILON && push!(θ0, logit(0.05))
# a start for the added rates on the command line: julia ... outdir variant log_rate...
length(ARGS) > 2 && (θ0[(end - length(ARGS) + 3):end] = parse.(Float64, ARGS[3:end]))
@printf("variant %s: start %.3f\n", VARIANT, negloglik(θ0))
flush(stdout)
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
    flush(stdout)
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

fit = fitted(q̂)
pred, sols = fit.pred, fit.sols
open(joinpath(OUT, "fit.txt"), "w") do io
    for out in (stdout, io)
        @printf(out, "variant %s: %d plot means\n", VARIANT, length(fit.r))
        @printf(out, "objective %.4f; Hessian eigenvalues %.3g to %.3g\n", Optim.minimum(res), extrema(eigvals(Symmetric(H)))...)
        println(out, "optimum on the unconstrained scale: ", θ̂)
        println(out, "parameter          estimate   (95% interval from the Laplace approximation)")
        tr = [k == 6 ? identity : startswith(NAMES[k], "later") || NAMES[k] == "epsilon" ? logistic : exp for k in eachindex(NAMES)]
        for (k, name) in enumerate(NAMES)
            @printf(out, "%-18s %9.4g   (%.4g, %.4g)\n", name, tr[k](θ̂[k]), tr[k](θ̂[k] - 1.96se[k]), tr[k](θ̂[k] + 1.96se[k]))
        end
        println(out, "\nyear day trt  observed (plot means, %)        model (%)")
        for y in YEARS, d in sort(unique(o.day for o in obs if o.year == y)), t in ("L", "M", "H")
            idx = findall(o -> o.year == y && o.day == d && o.treatment == t, obs)
            @printf(out, "%d %3d  %s   %-30s %6.1f\n", y, d, t, join((@sprintf("%5.1f", 100obs[i].mean) for i in idx), " "),
                    100pred[idx[1]])
        end
        if UNDATED
            println(out, "\nundated scores at their profiled days")
            for y in YEARS, c in UNDATED_COLUMNS[y], t in ("L", "M", "H")
                d = fit.dates[c]
                idx = findall(o -> o.column == c && o.treatment == t, undated)
                @printf(out, "%d %3d  %s   %-30s %6.1f   %s, %s%s\n", y, d, t,
                        join((@sprintf("%5.1f", 100undated[i].mean) for i in idx), " "),
                        100severity(sols[t], y - 2008, Float64(d), q̂), c, Date(y) + Day(d - 1),
                        c in REPLACED ? " (not fitted: replaced by its shoot scores)" : "")
            end
        end
        if SHOOTS
            println(out, "\nshoot scores at their companion's day (model: the cohort's visible share; for 2010 also",
                    " the model's mean of the two shoots against oid3_10)")
            for c in sort(collect(keys(SHOOT_COLUMNS))), t in ("L", "M", "H")
                idx = findall(o -> o.column == c && o.treatment == t, shoot)
                o1 = shoot[idx[1]]
                d = fit.dates[o1.paired]
                @printf(out, "%d %3d  %s   %-30s %6.1f   %s, cohort %d\n", o1.year, d, t,
                        join((@sprintf("%5.1f", 100shoot[i].mean) for i in idx), " "),
                        100visible_incidence(sols[t], o1.year - 2008, Float64(d); cohort = o1.cohort), c, o1.cohort)
            end
            for t in ("L", "M", "H")
                d = Float64(fit.dates["oid3_10"])
                @printf(out, "2010 %3d  %s   shoot mean %6.1f   leaf-area share %6.1f   oid3_10: %s\n", d, t,
                        50sum(visible_incidence(sols[t], 2, d; cohort = k) for k in 1:2), 100visible_incidence(sols[t], 2, d),
                        join((@sprintf("%5.1f", 100o.mean) for o in undated if o.column == "oid3_10" && o.treatment == t), " "))
            end
        end
        if SHED || SEVERITY || DROUGHT
            println(out, "\nleaves present at the final score, as a share of the leaves produced (%)",
                    LATE ? "; first flush, later flushes" : "")
            for y in YEARS, t in ("L", "M", "H")
                d = Float64(WINDOW[y][2])
                @printf(out, "%d %3d  %s   %5.1f", y, d, t, 100leaf_area(sols[t], y - 2008, d) / sum(flushes(q̂, y).sizes))
                if LATE
                    @printf(out, "   %5.1f  %5.1f", (100leaf_area(sols[t], y - 2008, d; cohort = c) / sum(cohorts_of(flushes(q̂, y))[c].sizes)
                                                  for c in 1:2)...)
                end
                println(out)
            end
        end
    end
end

colors = Dict("L" => 1, "M" => 2, "H" => 3)
panels = map(YEARS) do y
    plt = plot(; title = string(y), xlabel = "day of year", ylabel = y == 2009 ? "leaf area infected (%)" : "",
               legend = y == 2009 ? :topleft : false, ylims = (0, SHOOTS ? 100 : 75))
    for t in ("L", "M", "H")
        plot!(plt, DAYS, [100severity(sols[t], y - 2008, d, q̂) for d in DAYS]; lw = 2, color = colors[t],
              label = Dict("L" => "fungicide", "M" => "natural", "H" => "added spores")[t])
        o = filter(o -> o.year == y && o.treatment == t, obs)
        scatter!(plt, [x.day for x in o], [100x.mean for x in o]; color = colors[t], ms = 4, label = nothing)
        if UNDATED
            u = filter(o -> o.year == y && o.treatment == t, undated)
            scatter!(plt, [fit.dates[x.column] for x in u], [100x.mean for x in u]; color = colors[t], ms = 4,
                     markershape = :diamond, label = nothing)
        end
        if SHOOTS && y != 2009    # shoot scores: first shoot squares (dotted), second triangles (dashed)
            for (k, shape, style) in ((1, :rect, :dot), (2, :utriangle, :dash))
                s = filter(o -> o.year == y && o.treatment == t && o.cohort == k, shoot)
                isempty(s) && continue
                plot!(plt, DAYS, [100visible_incidence(sols[t], y - 2008, d; cohort = k) for d in DAYS]; lw = 1,
                      ls = style, color = colors[t], label = nothing)
                scatter!(plt, [fit.dates[x.paired] for x in s], [100x.mean for x in s]; color = colors[t], ms = 4,
                         markershape = shape, label = nothing)
            end
        end
    end
    plt
end
savefig(plot(panels...; layout = (1, 3), size = (1100, 360), left_margin = 4Plots.mm, bottom_margin = 5Plots.mm),
        joinpath(OUT, "fit.png"))
println("wrote ", OUT)
