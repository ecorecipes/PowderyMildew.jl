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
#   tauprior  an informative prior on tau_s, the mean time a young leaf stays fully susceptible:
#           log-normal, median 12 days, log-sd 0.35 (95 % within 6-24 days), for leaves that
#           resist mildew as they mature over two to three weeks; on top of the weak default;
#   spread  colonies spread on diseased leaves (`cover0`, `cover_rate`): a leaf is covered cover0
#           when it shows mildew and the rest at rate cover_rate, and every score is the covered
#           share of leaf area (`infected_area`) instead of the share of leaves showing mildew.
#           cover0 on the logit scale (prior median 0.3, log-odds sd 1.5), cover_rate on the log
#           scale (prior median 0.03 per day, log-sd 1.5).
#   efficacy  the fungicide's efficacy fitted (logit scale, prior median 0.9, log-odds sd 1.5):
#           protection just after a spray is partial, `1 - efficacy` of the force of infection
#           getting through; otherwise fixed at 1.
#   humidity  the force of infection times exp(-h (RH - 70)/10), RH the day's mean relative
#           humidity (%): h > 0 if dry air favours infection, h < 0 if humid air does; h normal,
#           prior mean 0, sd 1 (a factor of e per 10 points of RH at one sd).
#   humsummer  with "humidity", the factor applies from 1 June (day 152) only, 1 before: a
#           humidity effect on the summer epidemic but not on spring infection (fit 25: fitted
#           without 2011, h of -2 had no 2011 epidemic after the dry spring).
#   humcap  with "humidity", RH clamped to 60-80 % in the factor, so it lies within exp(+-h): no
#           single humid or dry day weighs more than one at 80 or 60 %.
#   rain    the force of infection times exp(-r) on wet days (precipitation >= 1 mm): r > 0 if rain
#           suppresses infection (conidia washed off or burst); r normal, prior mean 0, sd 1.5.
#   rain5   with "rain", a wet day needs at least 5 mm (also for "wetrelease").
#   rainlog with "rain", the factor is exp(-r log(1 + P)) with P the day's precipitation (mm), so
#           the effect grows with the amount instead of switching at a threshold.
#   wetrelease  the overwintered ascospores are released at a fraction w of the rate on dry days
#           (`ReleaseHazard(; wet)`: w + (1 - w) on wet days, so 1 on wet days), w fitted on the logit
#           scale, prior median 0.5, log-odds sd 1.5; w = 1 is the model without it. Chasmothecia
#           of powdery mildews generally need wetting to discharge.
#   gap     with "dates", the profiled undated scores are at least 7 days from each other and from
#           the dated scores either side (otherwise any whole days strictly between them): field
#           scorings a day apart cannot differ as much as consecutive scores sometimes do.
#   june2011  with "dates", both undated 2011 scores (oid3_11, oid4_11) no later than 20 June, when
#           the third flush was noted: an assumption (there is no third-shoot score), to test the
#           fits that move oid4_11 to late July.
#   presence  also fits the early presence scores oid1_10 and oid1_11 (the share p of seedlings
#           with any mildew; undated, before the first dated score): each year's day is profiled
#           between leaf emergence (mean budburst score 1) and that score. The model's whole-seedling score v maps to presence
#           through cloglog p = a + b log v, with a and b fixed by regressing cloglog(1 - zeros/n)
#           on log(mean) over every 2009-2011 percent scoring (one link for all years and
#           treatments). Residuals on the cloglog scale have their own sd (log-normal prior,
#           median 0.5, log-sd 1), the regression's error included.
#   ext2    the outside force of infection in 2010-11 is b_ext times a fitted ratio (the season's
#           `external` multiplier; 2009, when the seedlings were newly sown with no local
#           inoculum, keeps b_ext): log-normal prior, median 1, log-sd 2.
#   effyear the fungicide's efficacy fitted separately for each year (three logits, each with the
#           "efficacy" prior), for spraying that protected better in some years. Not with "efficacy".
#   decayyear  the fungicide's residual activity (`fungicide_decay`) in 2010 and 2011 is 2009's times
#           a fitted ratio (log-normal priors, median 1, log-sd 1).
#   hold2009, hold2010, hold2011  leave that year's scores (plot means, undated, shoot and
#           presence) out of the fit and report how well the fitted model predicts them: the
#           Gaussian -log density of the held-out residuals at the fitted sds. Every year is still
#           simulated (each season's inoculum comes from the one before), and the held-out
#           undated scores are still put at their best days.
#   fix_<name>=<value>  holds a parameter (a name from the fit.txt table) at a value on the
#           unconstrained scale (log, logit or identity, as in the printed optimum), e.g.
#           fix_humidity=-0.4, for profile likelihoods; it is not counted as a parameter.
#   cluster the plot-mean residuals of one year and treatment (nine cells) share a normal offset
#           of sd tau (`cluster_sd`; log-normal prior, median 0.3, log-sd 1.5) on top of their own
#           errors (obs_sd, now the within-cell sd), integrated out exactly. Misfits come in blocks
#           (L autumn 2010, M and H autumn 2009), so independent errors overstate the evidence
#           and held-out predictions' confidence. The undated days are profiled on squared
#           residuals, or with "lm" on the clustered likelihood.
#   plotyear  like "cluster", but the residuals of one plot in one year share the offset (27 groups,
#           4-6 scores each): a plot that runs high or low all season. `plot_sd`, same prior.
#   plot    the residuals of one plot share an offset in every year (9 groups): a persistent plot
#           (or block) effect. The held-out density is then conditional on the plot's fitted
#           years, the predictive density for a new year. One grouping at a time.
#   yearsd  each year's plot means have their own noise sd, obs_sd e^δ with δ ~ N(0, s²) (`year_sd`
#           = s; log-normal prior, median 0.3, log-sd 1.5), integrated out exactly per year
#           (adaptive Gauss-Hermite), so obs_sd is the median year's sd. A held-out year's density
#           integrates over its unknown sd. Presence scores keep their own sd. Not with a
#           grouping; the profiled days are unchanged (a year's term rises with its squared residuals).
#   lm      a numerical option, not a model: the seasons are integrated with Vern9 at tolerances
#           1e-14/1e-11 (the default Tsit5 at 1e-9/1e-7 leaves noise of about 1e-3 in the
#           objective, from the tiny early-season states on the logit scale, which no simplex can
#           converge on), and the objective is minimised by Levenberg-Marquardt on the residuals
#           with a finite-difference Jacobian, from the start and from LM_STARTS - 1 jittered
#           starts (environment variable, default 3). The noise sds are optimised in closed
#           form given the residuals, the profiled days are refreshed at every trial point, and
#           the Hessian is taken with the days held at their optimum. Objectives differ from
#           the default path's by about 1e-2 (its integration error).
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
using LinearAlgebra: Symmetric, inv, diag, eigvals, Diagonal, norm, eigen, SymTridiagonal
using Random: MersenneTwister

const DATA = joinpath(@__DIR__, "..", "data", "pierroton")
const OUT = length(ARGS) ≥ 1 ? ARGS[1] : joinpath(@__DIR__, "..", "runs", "pierroton_fit")
mkpath(OUT)
const VARIANT = length(ARGS) ≥ 2 ? ARGS[2] : "base"
const OPTIONS = Set(split(VARIANT, '+'))
issubset(filter(!startswith("fix_"), OPTIONS), Set(["base", "winter", "dates", "shed", "severity", "drought", "power", "late", "shoots", "sizes", "mean2010", "epsilon",
                       "tauprior", "spread", "efficacy", "humidity", "rain", "gap", "rain5", "rainlog",
                       "wetrelease", "june2011", "presence", "ext2", "effyear", "decayyear",
                       "hold2009", "hold2010", "hold2011", "humsummer", "humcap", "cluster", "plotyear", "plot", "yearsd", "lm"])) ||
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
const TAUPRIOR, SPREAD = "tauprior" in OPTIONS, "spread" in OPTIONS
const ISPREAD = IEPS + EPSILON  # logit cover0, then log cover_rate
const EFFICACY = "efficacy" in OPTIONS
const IEFF = ISPREAD + 2SPREAD  # logit of the fungicide's efficacy
const HUMID, RAIN = "humidity" in OPTIONS, "rain" in OPTIONS
const IHUM = IEFF + EFFICACY    # humidity coefficient h, then the wet-day coefficient r
const IRAIN = IHUM + HUMID
const HUMSUMMER, HUMCAP = "humsummer" in OPTIONS, "humcap" in OPTIONS
(HUMSUMMER || HUMCAP) && !HUMID && error("\"humsummer\" and \"humcap\" need \"humidity\"")
const RAIN5, RAINLOG, WETREL = "rain5" in OPTIONS, "rainlog" in OPTIONS, "wetrelease" in OPTIONS
(RAIN5 || RAINLOG) && !RAIN && error("\"rain5\" and \"rainlog\" need \"rain\"")
RAIN5 && RAINLOG && error("\"rain5\" and \"rainlog\" are alternatives; use one")
const WET_MM = RAIN5 ? 5.0 : 1.0   # the wet-day threshold
const IWREL = IRAIN + RAIN       # logit of the dry-day share of ascospore release
const PRESENCE = "presence" in OPTIONS
const IPRES = IWREL + WETREL     # log of the presence residuals' sd
const EXT2 = "ext2" in OPTIONS
const IEXT = IPRES + PRESENCE    # log of the 2010-11 outside input relative to 2009's
const EFFYEAR, DECAYYEAR = "effyear" in OPTIONS, "decayyear" in OPTIONS
EFFYEAR && EFFICACY && error("\"effyear\" and \"efficacy\" both set the efficacy; use one")
const IEFFY = IEXT + EXT2        # logits of the 2009, 2010 and 2011 efficacies
const IDECY = IEFFY + 3EFFYEAR   # logs of the 2010 and 2011 decay relative to 2009's
const CLUSTER = "cluster" in OPTIONS
const PLOTYEAR, PLOTRE = "plotyear" in OPTIONS, "plot" in OPTIONS
CLUSTER + PLOTYEAR + PLOTRE ≤ 1 || error("\"cluster\", \"plotyear\" and \"plot\" are alternatives; use one")
const GROUPED = CLUSTER || PLOTYEAR || PLOTRE   # residuals share normal offsets within groups
# a residual's group: year and treatment, plot and year, or plot
group_key(y, t, plot) = CLUSTER ? (y, t) : PLOTYEAR ? (y, string(plot)) : (0, string(plot))
const LM = "lm" in OPTIONS       # tight integration and Levenberg-Marquardt (numerical only)
const ICLU = IDECY + 2DECAYYEAR  # log of the group offsets' sd (year by treatment, or plot)
const YEARSD = "yearsd" in OPTIONS
YEARSD && GROUPED && error("\"yearsd\" is not combined with a grouping")
const IYSD = ICLU + GROUPED      # log of the sd of log(each year's noise sd)
# the held-out year (0: none)
const HOLD = let h = [y for y in 2009:2011 if "hold$y" in OPTIONS]
    length(h) ≤ 1 || error("hold out one year at a time")
    isempty(h) ? 0 : only(h)
end
# with "gap", the profiled undated scores are at least GAP days from each other and from the dated
# scores either side (whole days strictly between them otherwise, GAP = 1)
const GAP = "gap" in OPTIONS ? 7 : 1
# with "june2011", both undated 2011 scores no later than 20 June, when the third flush was noted
# (an assumption: there is no third-shoot score, as if scoring came before that flush)
const LATEST = "june2011" in OPTIONS ? Dict("oid3_11" => 171, "oid4_11" => 171) : Dict{String,Int}()   # 171 = 20 June 2011
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
# early presence scores ("presence"): the share of seedlings with any mildew, continuity-corrected
const presence = [(year = parse(Int, r["year"]), column = String(r["column"]), treatment = r["treatment"],
                   plot = parse(Int, r["plot"]), p = clamp(parse(Float64, r["mean"]), 0.5 / parse(Int, r["n"]),
                                                          1 - 0.5 / parse(Int, r["n"])))
                  for r in readcsv(joinpath(DATA, "plot_scores.csv"))
                  if PRESENCE && r["kind"] == "presence" && parse(Int, r["year"]) in YEARS]
@assert length(presence) == (PRESENCE ? 18 : 0)
cloglog(p) = log(-log1p(-p))
# the presence link cloglog p = a + b log v: least squares over every 2009-2011 percent scoring
# (whole seedling or one shoot) with 0 < p < 1, p = 1 - zeros/n and v the plot's mean score
const PRESENCE_LINK = let pts = [(log(parse(Float64, r["mean"]) / 100), cloglog(1 - parse(Int, r["zeros"]) / parse(Int, r["n"])))
                                 for r in readcsv(joinpath(DATA, "plot_scores.csv"))
                                 if parse(Int, r["year"]) in YEARS && startswith(r["kind"], "percent") &&
                                    0 < parse(Int, r["zeros"]) < parse(Int, r["n"]) && parse(Float64, r["mean"]) > 0]
    x, y = first.(pts), last.(pts)
    b = sum((x .- mean(x)) .* (y .- mean(y))) / sum(abs2, x .- mean(x))
    (a = mean(y) - b * mean(x), b = b, n = length(pts))
end
weather_rows = readcsv(joinpath(DATA, "weather.csv"))
temperature = Dict(y => DailySeries([parse(Float64, r["tmean"]) for r in weather_rows if r["year"] == string(y)])
                   for y in YEARS)
# daily mean relative humidity (%) and wet days (precipitation >= 1 mm), for "humidity" and "rain"
rel_humidity = Dict(y => DailySeries([parse(Float64, r["rh_mean"]) for r in weather_rows if r["year"] == string(y)])
                    for y in YEARS)
wet_day = Dict(y => DailySeries([parse(Float64, r["precip"]) ≥ WET_MM ? 1.0 : 0.0 for r in weather_rows if r["year"] == string(y)])
               for y in YEARS)
# the covariate of "rain": the wet-day indicator, or log(1 + P) with "rainlog"
rain_cov = RAINLOG ?
           Dict(y => DailySeries([log1p(parse(Float64, r["precip"])) for r in weather_rows if r["year"] == string(y)]) for y in YEARS) :
           wet_day

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
# the day the mean budburst score crosses `level`, interpolated between weekly scorings (the first
# scoring day if it is already reached then)
function score_day(year, level = 2)
    orders = sort(unique(parse(Int, r["order"]) for r in phen if r["year"] == string(year)))
    m = [mean(parse(Float64, r["mean"]) for r in phen if r["year"] == string(year) && r["order"] == string(k))
         for k in orders]
    k = findfirst(≥(level), m)
    wednesday15 = year == 2010 ? doy(2010, 4, 14) : doy(2011, 4, 13)
    k == 1 && return Float64(wednesday15)
    frac = (level - m[k - 1]) / (m[k] - m[k - 1])
    return wednesday15 + 7 * (orders[k - 1] - 1 + frac)
end
score2_day(year) = score_day(year, 2)
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
const SAVEDAYS = collect(TSPAN[1]:1.0:TSPAN[2])
const DAYS = collect(91.0:1.0:299.0)

# Ascospore release peaks a fitted lag after each year's budburst (score 2; 2009, emergence), so
# that it follows the same spring warmth as the host. The 2009 release has no source (C0 = 0).
const BUDBURST = Dict(y => FLUSH[y].dates[1] for y in YEARS)
# Fungicide efficacy is fixed at 1 (full protection just after a spray) unless "efficacy": fitted
# in fits 4-5, it went to the bound, with or without a free decay (wiki: oak/barres-fit.md). The
# decay is fitted.
const NAMES = ["b_c", "b_ext", "c_w", "dose", "tau_s", "release_lag", "fungicide_decay", "obs_sd",
               (WINTER ? ["winter_2011"] : String[])..., (SHED ? ["shed"] : String[])...,
               (SEVERITY ? ["shed_severity"] : String[])..., (DROUGHT ? ["shed_drought"] : String[])...,
               (POWER ? ["shed_power"] : String[])..., (SIZES ? ["later_2010", "later_2011"] : String[])...,
               (EPSILON ? ["epsilon"] : String[])..., (SPREAD ? ["cover0", "cover_rate"] : String[])...,
               (EFFICACY ? ["efficacy"] : String[])..., (HUMID ? ["humidity"] : String[])...,
               (RAIN ? ["rain"] : String[])..., (WETREL ? ["dry_release"] : String[])...,
               (PRESENCE ? ["presence_sd"] : String[])..., (EXT2 ? ["ext_2010_11"] : String[])...,
               (EFFYEAR ? ["efficacy_2009", "efficacy_2010", "efficacy_2011"] : String[])...,
               (DECAYYEAR ? ["decay_ratio_2010", "decay_ratio_2011"] : String[])...,
               (CLUSTER ? ["cluster_sd"] : String[])..., (PLOTYEAR || PLOTRE ? ["plot_sd"] : String[])...,
               (YEARSD ? ["year_sd"] : String[])...]
# "fix_<name>=<value>": the held parameters' positions in θ and values; the rest are optimised
const FIXED = Dict(begin
                       name, value = split(o[5:end], '=')
                       k = findfirst(==(name), NAMES)
                       k === nothing && error("fix: no parameter $name in this variant")
                       k => parse(Float64, value)
                   end for o in OPTIONS if startswith(o, "fix_"))
const FREE = [k for k in eachindex(NAMES) if !haskey(FIXED, k)]
function expand(φ)
    θ = zeros(eltype(φ), length(NAMES))
    θ[FREE] = φ
    for (k, v) in FIXED
        θ[k] = v
    end
    return θ
end
function unpack(θ)
    return (b_c = exp(θ[1]), b_ext = exp(θ[2]), c_w = exp(θ[3]), efficacy = EFFICACY ? logistic(θ[IEFF]) : 1.0,
            dose = exp(θ[4]),
            tau_s = exp(θ[5]), lag = θ[6], decay = exp(θ[7]), sd = exp(θ[8]),
            winter = WINTER ? exp(θ[9]) : 1.0,
            shed = SHED ? exp(θ[ISHED]) : DROUGHT ? exp(θ[IDRY]) : 0.0,
            shed_severity = SEVERITY ? exp(θ[ISEV]) : 0.0, shed_power = POWER ? exp(θ[IPOW]) : 1.0,
            later = SIZES ? Dict(2010 => logistic(θ[ISIZE]), 2011 => logistic(θ[ISIZE + 1])) : Dict{Int,Float64}(),
            epsilon = EPSILON ? logistic(θ[IEPS]) : 0.05,
            cover0 = SPREAD ? logistic(θ[ISPREAD]) : 1.0, cover_rate = SPREAD ? exp(θ[ISPREAD + 1]) : 0.0,
            humid = HUMID ? θ[IHUM] : 0.0, rain = RAIN ? θ[IRAIN] : 0.0,
            dry_release = WETREL ? logistic(θ[IWREL]) : 1.0, presence_sd = PRESENCE ? exp(θ[IPRES]) : 1.0,
            ext_late = EXT2 ? exp(θ[IEXT]) : 1.0,
            eff_year = EFFYEAR ? logistic.(θ[IEFFY:(IEFFY + 2)]) : nothing,
            decay_ratio = DECAYYEAR ? [1.0; exp.(θ[IDECY:(IDECY + 1)])] : [1.0, 1.0, 1.0],
            cluster = GROUPED ? exp(θ[ICLU]) : 0.0, year_sd = YEARSD ? exp(θ[IYSD]) : 0.0)
end
# the moisture factor on the force of infection ("humidity", "rain"; 1 otherwise)
# (the series are looked up once, outside the closure: the dictionaries are untyped globals, and a
# lookup on every right-hand-side call made each evaluation about 3 times slower)
function moisture(q, y)
    (HUMID || RAIN) || return t -> 1.0
    RH, W, h, r = rel_humidity[y], rain_cov[y], q.humid, q.rain
    # "humsummer": from 1 June (a whole day, so a tstop); "humcap": RH within 60-80 %. Without
    # them the arithmetic is exactly fit 18's, so earlier fits reproduce to the last digit
    from = HUMSUMMER ? 152.0 : -Inf
    return t -> begin
        x = HUMCAP ? clamp(RH(t), 60.0, 80.0) : RH(t)
        t ≥ from ? exp(-h * (x - 70) / 10 - r * W(t)) : exp(-r * W(t))
    end
end
# the wet gate on ascospore release ("wetrelease"): w on dry days, 1 on wet days
function release_gate(q, y)
    WETREL || return t -> 1.0
    W, w = wet_day[y], q.dry_release
    return t -> w + (1 - w) * W(t)
end

function simulate(q, treatment)
    p = MildewParameters(; b_c = q.b_c, b_a = 1.0, b_ext = q.b_ext, tau_s = q.tau_s, K = 3, epsilon = q.epsilon,
                         latent = 7.0, L = 3, gamma_I = 0.05, c_w = q.c_w, shed = q.shed,
                         shed_severity = q.shed_severity, shed_power = q.shed_power, cohorts = TWO ? 2 : 1,
                         cohort_shedding = LATE ? [0.0, 1.0] : TWO ? [1.0, 1.0] : [1.0], cover0 = q.cover0,
                         cover_rate = q.cover_rate)
    forcings = map(YEARS) do y
        T = temperature[y]
        release = ReleaseHazard(; rate = 0.2, optimum = BUDBURST[y] + q.lag, curvature = -0.003, wet = release_gate(q, y))
        spray = treatment == "L" ? FungicideSchedule(FUNGICIDE[y]; efficacy = EFFYEAR ? q.eff_year[y - 2008] : q.efficacy,
                                                     decay = q.decay * q.decay_ratio[y - 2008]) : t -> 1.0
        inoc = treatment == "H" ? [(d, q.dose) for d in INOCULATION[y]] : Tuple{Float64,Float64}[]
        wet = moisture(q, y)
        SeasonForcing(; tspan = TSPAN, production = TWO ? cohorts_of(flushes(q, y)) : flushes(q, y), release,
                      weather = t -> beta_response(T(t)) * spray(t) * wet(t),
                      development = t -> degree_day_rate(T(t); base = LATENT_BASE, reference = 20.0),
                      shedding = DROUGHT ? DEFICIT[y] : t -> 1.0,
                      # myclobutanil also suppresses chasmothecium formation on treated leaves;
                      # the 2010 gate carries the 2010/11 winter factor (1 unless "winter")
                      gate = t -> t ≥ 240 ? (y == 2010 ? q.winter : 1.0) * spray(t) : 0.0, inoculations = inoc,
                      tstops = [DAYS; treatment == "L" ? FUNGICIDE[y] : Float64[]],
                      external = y == 2009 ? 1.0 : q.ext_late)
    end
    # every score, profile and plot reads whole days only, which tstops hit exactly, so saving only
    # those leaves every value unchanged and is about twice as fast
    prob = MildewProblem(p, collect(forcings), 0.0)
    LM && return solve(prob; alg = PowderyMildew.Vern9(), abstol = 1e-14, reltol = 1e-11, saveat = SAVEDAYS, dense = false)
    return solve(prob; abstol = 1e-9, reltol = 1e-7, saveat = SAVEDAYS, dense = false)
end

logit(x) = log(x / (1 - x))
# the score is the percentage of leaf area infected, i.e. the visible share of leaf area
# (with "mean2010", in 2010 the seedlings' mean of their shoot scores; f, the share of the later
# flushes' leaves out, uses the leaves present, so it assumes no shedding)
# The leaf area infected: the covered share with "spread", otherwise the share of leaves showing
# mildew (infected_area is then exactly visible_incidence)
score(sol, season, day; cohort = nothing) = infected_area(sol, season, day; cohort)
function severity(sol, season, day, q)
    (MEAN2010 && season == 2) || return score(sol, season, day)
    f = min(leaf_area(sol, season, day; cohort = 2) / sum(flushes(q, 2010).sizes[2:end]), 1.0)
    v1, v2 = (score(sol, season, day; cohort = k) for k in 1:2)
    return v1 + f * (v2 - v1) / 2
end

lg(x) = logit(clamp(x, 1e-4, 1 - 1e-4))

# Model predictions and logit residuals. With "dates", each year's pair of undated scores is put
# at the ordered pair of whole days strictly between the dated scores that minimises their
# squared residuals: the profile over the dates, since nothing else in the objective depends on
# them. `pinned` (a column => day table, as returned in `dates`) holds the days instead.
function fitted(q; pinned = nothing)
    sols = Dict(t => simulate(q, t) for t in ("L", "M", "H"))
    pred = [severity(sols[o.treatment], o.year - 2008, Float64(o.day), q) for o in obs]
    r = [lg(o.mean) - lg(m) for (o, m) in zip(obs, pred)]
    ry = [o.year for o in obs]   # each residual's year, for "hold<year>"
    rt = [o.treatment for o in obs]   # and treatment, for "cluster"
    rk = [o.plot for o in obs]        # and plot, for "plotyear" and "plot"
    dates = Dict{String,Int}()
    # each year's undated candidates, for the grouped day search: the days, the feasible first and
    # second days, and per residual its position in r, column (1 or 2) and values on every day
    cand = Dict{Int,Any}()
    if UNDATED
        for y in YEARS
            days = (WINDOW[y][1] + GAP):(WINDOW[y][2] - GAP)
            curve = Dict(t => [lg(severity(sols[t], y - 2008, Float64(d), q)) for d in days] for t in ("L", "M", "H"))
            c1, c2 = UNDATED_COLUMNS[y]
            # with "shoots", each cohort's visible share; the shoot scores join their companion's day
            ys = [o for o in shoot if o.year == y]
            cc = Dict((t, k) => [lg(score(sols[t], y - 2008, Float64(d); cohort = k)) for d in days]
                      for t in ("L", "M", "H"), k in unique(o.cohort for o in ys))
            ssr(c) = [sum((abs2(lg(o.mean) - curve[o.treatment][i]) for o in undated
                           if o.column == c && !(c in REPLACED)); init = 0.0) +
                      sum((abs2(lg(o.mean) - cc[(o.treatment, o.cohort)][i]) for o in ys if o.paired == c); init = 0.0)
                      for i in eachindex(days)]
            s1, s2 = ssr(c1), ssr(c2)
            for (s, c) in ((s1, c1), (s2, c2)), i in eachindex(days)
                days[i] < get(EARLIEST, c, 0) && (s[i] = Inf)
                days[i] > get(LATEST, c, typemax(Int)) && (s[i] = Inf)
            end
            best, i1, i2, low, at = Inf, 0, 0, Inf, 0
            for j in (1 + GAP):length(days)   # best i1 ≤ i2 - GAP: running minimum of s1 up to j - GAP
                s1[j - GAP] < low && ((low, at) = (s1[j - GAP], j - GAP))
                low + s2[j] < best && ((best, i1, i2) = (low + s2[j], at, j))
            end
            pinned === nothing ||((i1, i2) = (findfirst(==(pinned[c1]), days), findfirst(==(pinned[c2]), days)))
            dates[c1], dates[c2] = days[i1], days[i2]
            entries = []
            for o in undated
                if o.year == y && !(o.column in REPLACED)
                    push!(r, lg(o.mean) - curve[o.treatment][o.column == c1 ? i1 : i2])
                    push!(rt, o.treatment)
                    push!(rk, o.plot)
                    push!(entries, (pos = length(r), col = o.column == c1 ? 1 : 2, v = lg(o.mean) .- curve[o.treatment]))
                end
            end
            for o in ys
                push!(r, lg(o.mean) - cc[(o.treatment, o.cohort)][o.paired == c1 ? i1 : i2])
                push!(rt, o.treatment)
                push!(rk, o.plot)
                push!(entries, (pos = length(r), col = o.paired == c1 ? 1 : 2, v = lg(o.mean) .- cc[(o.treatment, o.cohort)]))
            end
            append!(ry, fill(y, length(r) - length(ry)))
            cand[y] = (; days, c1, c2, ok1 = isfinite.(s1), ok2 = isfinite.(s2), entries)
        end
    end
    rg = [group_key(ry[i], rt[i], rk[i]) for i in eachindex(r)]
    # Grouped errors with "lm": each year's pair of days minimises the grouped -log likelihood
    # instead of the squared residuals, so that the profiled objective is the one being minimised
    # (the squared-residual pair makes it jump where the two disagree, which stalls the
    # optimiser). Each group's quadratic form is Q/(2σ²) - τ²S²/(2σ²(σ² + nτ²)), S and Q the sum
    # and sum of squares of its residuals; residuals of the held-out year count only in its own
    # search (its predictive density, given the fitted years). A persistent plot effect couples
    # the years, so the years are swept three times.
    if UNDATED && GROUPED && LM && pinned === nothing
        σ2, τ2 = q.sd^2, q.cluster^2
        for _ in 1:(PLOTRE ? 3 : 1), y in YEARS
            C = cand[y]
            mine = falses(length(r))
            foreach(e -> mine[e.pos] = true, C.entries)
            base = Dict{Any,Vector{Float64}}()   # group => [S, Q, n] of the other residuals
            for i in eachindex(r)
                (mine[i] || (y != HOLD && ry[i] == HOLD)) && continue
                b = get!(base, rg[i], zeros(3))
                b .+= (r[i], r[i]^2, 1.0)
            end
            gs = unique(rg[e.pos] for e in C.entries)
            nd = length(C.days)
            st = Dict((g, c) => (zeros(nd), zeros(nd), 0) for g in gs, c in 1:2)
            for e in C.entries
                S, Q, n = st[(rg[e.pos], e.col)]
                S .+= e.v
                Q .+= e.v .^ 2
                st[(rg[e.pos], e.col)] = (S, Q, n + 1)
            end
            best, i1, i2 = Inf, 0, 0
            for j in (1 + GAP):nd, i in 1:(j - GAP)
                (C.ok1[i] && C.ok2[j]) || continue
                f = 0.0
                for g in gs
                    b = get(base, g, zeros(3))
                    a1, a2 = st[(g, 1)], st[(g, 2)]
                    S = b[1] + a1[1][i] + a2[1][j]
                    Q = b[2] + a1[2][i] + a2[2][j]
                    n = b[3] + a1[3] + a2[3]
                    f += Q / (2σ2) - τ2 * S^2 / (2σ2 * (σ2 + n * τ2))
                end
                f < best && ((best, i1, i2) = (f, i, j))
            end
            dates[C.c1], dates[C.c2] = C.days[i1], C.days[i2]
            for e in C.entries
                r[e.pos] = e.v[e.col == 1 ? i1 : i2]
            end
        end
    end
    # presence scores: each year's day profiled between leaf emergence (mean budburst score 1) and
    # the first dated score
    rp, rpy = Float64[], Int[]
    for y in YEARS
        ps = [o for o in presence if o.year == y]
        isempty(ps) && continue
        days = ceil(Int, score_day(y, 1)):(WINDOW[y][1] - GAP)
        link = Dict(t => [PRESENCE_LINK.a + PRESENCE_LINK.b * log(max(severity(sols[t], y - 2008, Float64(d), q), 1e-12))
                          for d in days] for t in ("L", "M", "H"))
        s = [sum(abs2(cloglog(o.p) - link[o.treatment][i]) for o in ps) for i in eachindex(days)]
        i = pinned === nothing ? argmin(s) : findfirst(==(pinned[first(ps).column]), days)
        dates[first(ps).column] = days[i]
        append!(rp, [cloglog(o.p) - link[o.treatment][i] for o in ps])
        append!(rpy, fill(y, length(ps)))
    end
    return (; pred, r, rp, ry, rt, rk, rg, rpy, dates, sols)
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
                 (EPSILON ? 0.5 * ((θ[IEPS] - logit(0.05)) / 1.5)^2 : 0.0) -
                 (TAUPRIOR ? 0.5 * ((θ[5] - log(12.0)) / 0.35)^2 : 0.0) -
                 (SPREAD ? 0.5 * ((θ[ISPREAD] - logit(0.3)) / 1.5)^2 + 0.5 * ((θ[ISPREAD + 1] - log(0.03)) / 1.5)^2 : 0.0) -
                 (EFFICACY ? 0.5 * ((θ[IEFF] - logit(0.9)) / 1.5)^2 : 0.0) -
                 (HUMID ? 0.5 * θ[IHUM]^2 : 0.0) - (RAIN ? 0.5 * (θ[IRAIN] / 1.5)^2 : 0.0) -
                 (WETREL ? 0.5 * (θ[IWREL] / 1.5)^2 : 0.0) -
                 (PRESENCE ? 0.5 * (θ[IPRES] - log(0.5))^2 : 0.0) -
                 (EXT2 ? 0.5 * (θ[IEXT] / 2)^2 : 0.0) -
                 (EFFYEAR ? 0.5 * sum(abs2, (θ[IEFFY:(IEFFY + 2)] .- logit(0.9)) ./ 1.5) : 0.0) -
                 (DECAYYEAR ? 0.5 * sum(abs2, θ[IDECY:(IDECY + 1)]) : 0.0) -
                 (GROUPED ? 0.5 * ((θ[ICLU] - log(0.3)) / 1.5)^2 : 0.0) -
                 (YEARSD ? 0.5 * ((θ[IYSD] - log(0.3)) / 1.5)^2 : 0.0)

# the Gaussian -log likelihood of the plot means (logit scale, sd q.sd) and, with "presence", of
# the presence scores (cloglog scale, sd q.presence_sd); `constant` adds the 2π terms. Only the
# fitted years count ("hold<year>" leaves one out); `held = true` gives the held-out year's instead.
gauss(r, sd) = length(r) * log(sd) + sum(abs2, r; init = 0.0) / (2sd^2)
# "cluster": the residuals of one year and treatment share a normal offset of sd τ on top of
# their own error of sd σ (covariance σ²I + τ²11ᵀ), integrated out exactly: gauss when τ = 0
function cell_nll(r, σ, τ)
    n, v = length(r), σ^2 + length(r) * τ^2
    return (n - 1) * log(σ) + log(v) / 2 + (sum(abs2, r; init = 0.0) - τ^2 / v * sum(r; init = 0.0)^2) / (2σ^2)
end
# "yearsd": one year's plot-mean residuals r with sd σ e^δ, δ ~ N(0, s²) integrated out:
# -log ∫ N(δ; 0, s²) Π N(r_i; 0, σ² e^{2δ}) dδ without the 2π of the residuals (as gauss), by
# Gauss-Hermite around the integrand's mode (log-concave in δ; Newton converges). Also returns
# w = E[e^{-2δ} | r], the weight that makes r √w / σ whiten the gradient exactly, and the mode δ.
# As s → 0 it is gauss(r, σ).
const GH = let m = 30, E = eigen(SymTridiagonal(zeros(m), sqrt.((1:(m - 1)) ./ 2)))   # Golub-Welsch
    (x = E.values, w = sqrt(π) .* E.vectors[1, :] .^ 2)
end
function year_nll(r, σ, s)
    n, S = length(r), sum(abs2, r; init = 0.0)
    n == 0 && return (nll = 0.0, w = 1.0, δ = 0.0)
    g(δ) = n * (log(σ) + δ) + S * exp(-2δ) / (2σ^2) + δ^2 / (2s^2)
    δ = 0.0
    for _ in 1:500
        a = S * exp(-2δ) / σ^2
        step = (n - a + δ / s^2) / (2a + 1 / s^2)
        δ -= step
        abs(step) < 1e-13 && break
    end
    c = sqrt(2 / (2S * exp(-2δ) / σ^2 + 1 / s^2))
    l = [log(wk) + xk^2 - g(δ + c * xk) for (xk, wk) in zip(GH.x, GH.w)]
    m = maximum(l)
    p = exp.(l .- m)
    nll = log(s * sqrt(2π)) - log(c) - m - log(sum(p))
    w = sum(p .* exp.(-2 .* (δ .+ c .* GH.x))) / sum(p)
    return (; nll, w, δ)
end
# With grouped errors, the held-out year's density is conditional on the fitted years' residuals
# in the same groups (joint minus fitted): the predictive density. For year-local groups
# ("cluster", "plotyear") the fitted part is empty and it is the held-out residuals' own density.
function data_nll(fr, q; constant = false, held = false)
    keep(ys) = held ? ys .== HOLD : ys .!= HOLD
    k = keep(fr.ry)
    r, rp = fr.r[k], fr.rp[keep(fr.rpy)]
    n = length(r) + length(rp)
    means = if GROUPED && held
        fitd = fr.ry .!= HOLD
        sum(cell_nll(fr.r[fr.rg .== Ref(g)], q.sd, q.cluster) - cell_nll(fr.r[(fr.rg .== Ref(g)) .& fitd], q.sd, q.cluster)
            for g in unique(fr.rg[k]); init = 0.0)
    elseif GROUPED
        rg = fr.rg[k]
        sum(cell_nll(r[rg .== Ref(g)], q.sd, q.cluster) for g in unique(rg); init = 0.0)
    elseif YEARSD   # years are independent given σ and s, so a held-out year's term is its own
        ry = fr.ry[k]
        sum(year_nll(r[ry .== y], q.sd, q.year_sd).nll for y in unique(ry); init = 0.0)
    else
        gauss(r, q.sd)
    end
    return means + gauss(rp, q.presence_sd) + (constant ? n / 2 * log(2π) : 0.0)
end

function negloglik(θ)
    q = unpack(θ)
    fr = try
        fitted(q)
    catch err
        err isa PowderyMildewError || rethrow()
        return Inf
    end
    return data_nll(fr, q) - logprior(θ, q)
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
SPREAD && push!(θ0, logit(0.3), log(0.03))
EFFICACY && push!(θ0, logit(0.9))
HUMID && push!(θ0, 0.0)
RAIN && push!(θ0, 0.0)
WETREL && push!(θ0, logit(0.9))
PRESENCE && push!(θ0, log(0.5))
EXT2 && push!(θ0, 0.0)
EFFYEAR && push!(θ0, logit(0.95), logit(0.95), logit(0.95))
DECAYYEAR && push!(θ0, 0.0, 0.0)
GROUPED && push!(θ0, log(0.3))
YEARSD && push!(θ0, log(0.3))
# a start for the added rates on the command line: julia ... outdir variant log_rate...
length(ARGS) > 2 && (θ0[(end - length(ARGS) + 3):end] = parse.(Float64, ARGS[3:end]))
for (k, v) in FIXED
    θ0[k] = v
end
@printf("variant %s: start %.3f\n", VARIANT, negloglik(θ0))
flush(stdout)
# the objective over the free parameters ("fix_" holds the others)
objective(φ) = negloglik(expand(φ))
# Nelder-Mead from three starts, each restarted from its own optimum until it stops improving
function minimise(start)
    res = optimize(objective, start[FREE], NelderMead(), Optim.Options(; iterations = 3000, f_reltol = 1e-10))
    for _ in 1:6
        again = optimize(objective, Optim.minimizer(res), NelderMead(), Optim.Options(; iterations = 3000, f_reltol = 1e-10))
        improved = Optim.minimum(res) - Optim.minimum(again)
        res = again
        improved < 1e-6 && break
    end
    return res
end
# ---- "lm": Levenberg-Marquardt on the residuals ------------------------------------------------
# The residuals depend on the dynamic parameters only; the noise sds (obs_sd, presence_sd,
# cluster_sd) are optimised at fixed residuals, with no integration.
const NOISE = [8; PRESENCE ? [IPRES] : Int[]; GROUPED ? [ICLU] : Int[]; YEARSD ? [IYSD] : Int[]]
const DYN = [k for k in FREE if !(k in NOISE)]
const NOISEFREE = [k for k in FREE if k in NOISE]
total(θ, fr) = (q = unpack(θ); data_nll(fr, q) - logprior(θ, q))
function with_noise(θ, fr)
    isempty(NOISEFREE) && return θ
    f(x) = (ϑ = copy(θ); ϑ[NOISEFREE] = x; total(ϑ, fr))
    x = θ[NOISEFREE]
    for _ in 1:3
        x = Optim.minimizer(optimize(f, x, NelderMead(), Optim.Options(; iterations = 5000, g_abstol = 1e-13)))
    end
    ϑ = copy(θ)
    ϑ[NOISEFREE] = x
    return ϑ
end
# v with sum(abs2, v)/2 the misfit's quadratic form: x/σ for independent errors; with "cluster",
# (x - c Σx)/σ in each cell, c = (1 - sqrt(σ²/(σ² + nτ²)))/n. Linear in x, so it whitens the
# Jacobian's columns too.
function whiten(x, xp, fr, q)
    k = fr.ry .!= HOLD
    v = x[k] ./ q.sd
    if GROUPED
        rg = fr.rg[k]
        for g in unique(rg)
            m = rg .== Ref(g)
            n = count(m)
            c = (1 - sqrt(q.sd^2 / (q.sd^2 + n * q.cluster^2))) / n
            v[m] .-= c * sum(v[m])
        end
    elseif YEARSD   # each year's weight from the current residuals (not from x, a Jacobian column)
        ry = fr.ry[k]
        for y in unique(ry)
            v[ry .== y] .*= sqrt(year_nll(fr.r[k][ry .== y], q.sd, q.year_sd).w)
        end
    end
    return [v; xp[fr.rpy .!= HOLD] ./ q.presence_sd]
end
# central differences, with the profiled days held at the current point's
function jacobian(θ, fr, q; h = 1e-4)
    cols = map(DYN) do k
        e = zeros(length(θ))
        e[k] = h
        a, b = fitted(unpack(θ + e); pinned = fr.dates), fitted(unpack(θ - e); pinned = fr.dates)
        whiten((a.r - b.r) ./ 2h, (a.rp - b.rp) ./ 2h, fr, q)
    end
    return reduce(hcat, cols)
end
function hessian(f, x; h = 1e-2)
    n = length(x)
    H = zeros(n, n)
    for i in 1:n, j in i:n
        e_i, e_j = (1:n .== i) .* h, (1:n .== j) .* h
        H[i, j] = H[j, i] = (f(x + e_i + e_j) - f(x + e_i - e_j) - f(x - e_i + e_j) + f(x - e_i - e_j)) / (4h^2)
    end
    return H
end
# the priors are quadratic in θ, so central differences give their derivatives to rounding
function prior_derivatives(θ; h = 1e-3)
    f(x) = (ϑ = copy(θ); ϑ[DYN] = x; -logprior(ϑ, unpack(ϑ)))
    x = θ[DYN]
    g = [(f(x + (eachindex(x) .== i) .* h) - f(x - (eachindex(x) .== i) .* h)) / 2h for i in eachindex(x)]
    return g, hessian(f, x; h)
end
function lm(start; maxit = 300)
    θ = copy(start)
    fr = fitted(unpack(θ))
    θ = with_noise(θ, fr)
    f = total(θ, fr)
    μ, its, gmax = 1e-3, 0, Inf
    for _ in 1:maxit
        q = unpack(θ)
        J, v = jacobian(θ, fr, q), whiten(fr.r, fr.rp, fr, q)
        gp, Hp = prior_derivatives(θ)
        g, H = J' * v + gp, J' * J + Hp
        gmax = maximum(abs, g)
        its += 1
        gain, step = -Inf, Inf
        while μ < 1e10
            Δ = -((H + μ * Diagonal(diag(H))) \ g)
            θt = copy(θ)
            θt[DYN] += Δ
            frt = try
                fitted(unpack(θt))   # the days are profiled afresh
            catch err
                err isa PowderyMildewError || rethrow()
                nothing
            end
            if frt !== nothing
                θt = with_noise(θt, frt)
                ft = total(θt, frt)
                if ft < f
                    gain, step = f - ft, norm(Δ)
                    θ, fr, f = θt, frt, ft
                    μ = max(μ / 10, 1e-12)
                    break
                end
            end
            μ *= 10
        end
        (gain < 1e-9 || step < 1e-7) && break   # no step lowers it (gain = -Inf), or the last one was negligible
    end
    return (; θ, f, its, gmax)
end

if LM
    rng = MersenneTwister(1)
    found = nothing
    for k in 1:parse(Int, get(ENV, "LM_STARTS", "3"))
        start = copy(θ0)
        # jittered starts: log-sd 0.3 on the dynamic parameters (3 days on the release lag)
        k > 1 && (start[DYN] .+= [j == 6 ? 3.0 : 0.3 for j in DYN] .* randn(rng, length(DYN)))
        t = @elapsed out = try
            lm(start)
        catch err
            err isa PowderyMildewError || rethrow()
            nothing
        end
        out === nothing && (println("start $k: the model failed at the start"); continue)
        @printf("start %d: objective %.6f after %d iterations (%.0f s), largest gradient entry %.2g\n", k, out.f, out.its, t, out.gmax)
        flush(stdout)
        global found = found === nothing || out.f < found.f ? out : found
    end
    θ̂, fmin = found.θ, found.f
else
    best = nothing
    for (k, start) in enumerate((θ0,))
        r = minimise(start)
        @printf("start %d: objective %.4f (converged: %s)\n", k, Optim.minimum(r), Optim.converged(r))
        flush(stdout)
        global best = best === nothing || Optim.minimum(r) < Optim.minimum(best) ? r : best
    end
    θ̂, fmin = expand(Optim.minimizer(best)), Optim.minimum(best)
end
q̂ = unpack(θ̂)

# Laplace approximation: central-difference Hessian on the unconstrained scale (over the free
# parameters; a held one has no interval). With "lm" the profiled days are held at the optimum's,
# which makes the objective smooth, and the step is smaller.
H = if LM
    days = fitted(q̂).dates
    hessian(φ -> (θ = expand(φ); total(θ, fitted(unpack(θ); pinned = days))), θ̂[FREE]; h = 2e-3)
else
    hessian(objective, θ̂[FREE])
end
Σ = inv(Symmetric(H))
se = zeros(length(θ̂))
se[FREE] = sqrt.(max.(diag(Σ), 0))

fit = fitted(q̂)
pred, sols = fit.pred, fit.sols
# (year, treatment) of each residual, in the order fitted() builds them
labels = [(o.year, o.treatment) for o in obs]
if UNDATED
    for y in YEARS
        append!(labels, [(o.year, o.treatment) for o in undated if o.year == y && !(o.column in REPLACED)])
        append!(labels, [(o.year, o.treatment) for o in shoot if o.year == y])
    end
end
@assert length(labels) == length(fit.r)
# the likelihood without the priors, for comparing variants whose priors differ ("tauprior"),
# and information criteria with the Gaussian constant (the profiled days of the undated scores,
# the same in every variant with "dates", are not counted)
nobs = count(!=(HOLD), fit.ry) + count(!=(HOLD), fit.rpy)
npar = length(FREE)
nll = data_nll(fit, q̂; constant = true)
open(joinpath(OUT, "fit.txt"), "w") do io
    for out in (stdout, io)
        @printf(out, "variant %s: %d plot means%s\n", VARIANT, length(fit.r),
                PRESENCE ? @sprintf(" and %d presence scores (link cloglog p = %.3f + %.3f log v, from %d scorings)",
                                    length(fit.rp), PRESENCE_LINK.a, PRESENCE_LINK.b, PRESENCE_LINK.n) : "")
        @printf(out, "objective %.4f; Hessian eigenvalues %.3g to %.3g\n", fmin, extrema(eigvals(Symmetric(H)))...)
        @printf(out, "-log likelihood %.4f (with the Gaussian constant), -log prior %.4f; %d parameters, AIC %.2f, BIC %.2f\n",
                nll, -logprior(θ̂, q̂), npar, 2nll + 2npar, 2nll + npar * log(nobs))
        isempty(FIXED) || println(out, "held fixed: ", join((@sprintf("%s = %.4g", NAMES[k], v) for (k, v) in sort(collect(FIXED))), ", "))
        if HOLD != 0
            rh, rph = fit.r[fit.ry .== HOLD], fit.rp[fit.rpy .== HOLD]
            @printf(out, "held out %d (%d plot means, %d presence scores): -log predictive density %.4f (with the Gaussian constant), sum of squared logit residuals %.2f%s; fitted on %d scores\n",
                    HOLD, length(rh), length(rph), data_nll(fit, q̂; constant = true, held = true), sum(abs2, rh),
                    isempty(rph) ? "" : @sprintf(", of cloglog presence residuals %.2f", sum(abs2, rph)), nobs)
        end
        if GROUPED   # each group's offset, its posterior mean τ² Σr / (σ² + n τ²) from the fitted residuals, logit scale
            fitd = fit.ry .!= HOLD
            label(g) = CLUSTER ? @sprintf("%d %s", g...) : PLOTYEAR ? @sprintf("%d %s%s", g[1], only(unique(fit.rt[fit.rg .== Ref(g)])), g[2]) :
                       @sprintf("%s%s", only(unique(fit.rt[fit.rg .== Ref(g)])), g[2])
            println(out, CLUSTER ? "year-by-treatment" : PLOTYEAR ? "plot-by-year" : "plot",
                    " offsets (n, mean residual, estimated offset): ",
                    join((let c = fit.r[(fit.rg .== Ref(g)) .& fitd], n = length(c)
                              @sprintf("%s %d %.2f %.2f", label(g), n, sum(c) / n, q̂.cluster^2 * sum(c) / (q̂.sd^2 + n * q̂.cluster^2))
                          end for g in sort(unique(fit.rg[fitd]))), "; "))
        end
        if YEARSD   # each year's noise sd at the posterior mode of δ, given that year's own residuals
            println(out, "year noise sds (σ e^δ at the mode given the year's residuals; n, rms residual): ",
                    join((let c = fit.r[fit.ry .== y]
                              @sprintf("%d%s %.3f (%d, %.3f)", y, y == HOLD ? " held out" : "",
                                       q̂.sd * exp(year_nll(c, q̂.sd, q̂.year_sd).δ), length(c), sqrt(sum(abs2, c) / length(c)))
                          end for y in YEARS), "; "))
        end
        @printf(out, "sum of squared logit residuals %.2f; by treatment L %.2f (%d), M %.2f (%d), H %.2f (%d); by year %s\n",
                sum(abs2, fit.r),
                Iterators.flatten((sum(abs2, fit.r[i] for i in eachindex(labels) if labels[i][2] == t),
                                   count(l -> l[2] == t, labels)) for t in ("L", "M", "H"))...,
                join((@sprintf("%d %.2f (%d)", y, sum(abs2, fit.r[i] for i in eachindex(labels) if labels[i][1] == y),
                               count(l -> l[1] == y, labels)) for y in YEARS), ", "))
        println(out, "optimum on the unconstrained scale: ", θ̂)
        println(out, "parameter          estimate   (95% interval from the Laplace approximation)")
        tr = [k == 6 || NAMES[k] in ("humidity", "rain") ? identity :
              startswith(NAMES[k], "later") || startswith(NAMES[k], "efficacy") || NAMES[k] in ("epsilon", "cover0", "dry_release") ? logistic : exp
              for k in eachindex(NAMES)]
        for (k, name) in enumerate(NAMES)
            if haskey(FIXED, k)
                @printf(out, "%-18s %9.4g   (held fixed)\n", name, tr[k](θ̂[k]))
            else
                @printf(out, "%-18s %9.4g   (%.4g, %.4g)\n", name, tr[k](θ̂[k]), tr[k](θ̂[k] - 1.96se[k]), tr[k](θ̂[k] + 1.96se[k]))
            end
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
        if PRESENCE
            println(out, "\npresence scores at their profiled days (share of seedlings with mildew, %; model through the link)")
            for y in YEARS, t in ("L", "M", "H")
                ps = [o for o in presence if o.year == y && o.treatment == t]
                isempty(ps) && continue
                d = fit.dates[first(ps).column]
                v = severity(sols[t], y - 2008, Float64(d), q̂)
                @printf(out, "%d %3d  %s   %-30s %6.1f   %s, %s (model score %.3f %%)\n", y, d, t,
                        join((@sprintf("%5.1f", 100o.p) for o in ps), " "),
                        100(1 - exp(-exp(PRESENCE_LINK.a + PRESENCE_LINK.b * log(max(v, 1e-12))))),
                        first(ps).column, Date(y) + Day(d - 1), 100v)
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
                        100score(sols[t], o1.year - 2008, Float64(d); cohort = o1.cohort), c, o1.cohort)
            end
            for t in ("L", "M", "H")
                d = Float64(fit.dates["oid3_10"])
                @printf(out, "2010 %3d  %s   shoot mean %6.1f   leaf-area share %6.1f   oid3_10: %s\n", d, t,
                        50sum(score(sols[t], 2, d; cohort = k) for k in 1:2), 100score(sols[t], 2, d),
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
                plot!(plt, DAYS, [100score(sols[t], y - 2008, d; cohort = k) for d in DAYS]; lw = 1,
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
