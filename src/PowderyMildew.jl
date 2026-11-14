"""
    PowderyMildew

Seasonal, leaf-level models of oak powdery mildew (*Erysiphe alphitoides* and relatives on
*Quercus robur* and *Q. petraea*). It holds a faithful implementation of the seasonal model of
Hamelin et al. (2016, doi:10.1002/ecs2.1517), the baseline against which the package's
extensions (leaf flushes and ageing, a latent stage, weather forcing, a chasmothecial winter
reservoir and spatial inoculum) are checked.
"""
module PowderyMildew

using CommonSolve: CommonSolve, solve
using Distributions: Binomial, logpdf
using LinearAlgebra: cholesky
using LogExpFunctions: logistic
using OrdinaryDiffEqTsit5: Tsit5
using OrdinaryDiffEqVerner: Vern9
using SciMLBase: SciMLBase, ODEProblem
using SparseArrays: sparse

export solve
export PowderyMildewError, SolverError
export HamelinParameters, HamelinProblem, HamelinSolution, nspecies, between_season, frequency, incidence,
       coexistence_chi, frequency_loglikelihood
export MildewParameters, SeasonForcing, MildewProblem, MildewSolution, FlushProduction, ReleaseHazard, DailySeries,
       FungicideSchedule, degree_day_rate,
       beta_response, leaf_area, visible_incidence, infected_area, winter_inoculum
export Landscape, nstands, gridded, Dispersal, dispersal_field, SpatialMildewProblem, SpatialMildewSolution,
       regional_incidence

include("errors.jl")
include("hamelin.jl")
include("leaf_model.jl")
include("spatial.jl")

end
