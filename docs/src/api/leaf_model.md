# Leaf-level model

```@meta
CurrentModule = PowderyMildew
```

## Parameters and forcing

```@docs
MildewParameters
SeasonForcing
FlushProduction
ReleaseHazard
DailySeries
beta_response
FungicideSchedule
degree_day_rate
```

## Problem and solution

```@docs
MildewProblem
CommonSolve.solve
MildewSolution
leaf_area
visible_incidence
winter_inoculum
```

## Errors

```@docs
PowderyMildewError
SolverError
```
