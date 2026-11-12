@testset "fungicide and added inoculum" begin
    @testset "fungicide schedule" begin
        s = FungicideSchedule([20.0, 10.0]; efficacy = 0.8, decay = 5.0)
        @test s(9.99) == 1.0                                   # before the first application
        @test s(10.0) ≈ 0.2
        @test s(15.0) ≈ 1 - 0.8exp(-1)
        @test s(20.0) ≈ 0.2                                    # the latest application counts
        @test s(30.0) ≈ 1 - 0.8exp(-2)
        @test_throws ArgumentError FungicideSchedule([1.0]; efficacy = 1.5, decay = 5.0)
        @test_throws ArgumentError FungicideSchedule([1.0]; efficacy = 0.5, decay = 0.0)
        # it scales the force of infection: no infection at all under full efficacy forever
        p = MildewParameters(; b_c = 0.0, b_a = 0.0, b_ext = 0.1)
        full = FungicideSchedule([0.0]; efficacy = 1.0, decay = 1e9)
        sol = solve(MildewProblem(p, SeasonForcing(; tspan = (0.0, 30.0), leaves0 = 1.0, weather = full), 0.0))
        @test incidence(sol, 1, 30.0) < 1e-6
        late = FungicideSchedule([0.5]; efficacy = 1.0, decay = 1e9)   # unprotected for half a day
        sol = solve(MildewProblem(p, SeasonForcing(; tspan = (0.0, 30.0), leaves0 = 1.0, weather = late,
                                                   tstops = [0.5]), 0.0))
        @test incidence(sol, 1, 30.0) ≈ 1 - exp(-0.05) rtol = 1e-6
        none = solve(MildewProblem(p, SeasonForcing(; tspan = (0.0, 30.0), leaves0 = 1.0), 0.0))
        @test incidence(none, 1, 30.0) ≈ 1 - exp(-3.0) rtol = 1e-8
    end

    @testset "an inoculation infects a share 1 - exp(-dose f)" begin
        p = MildewParameters(; b_c = 0.0, b_a = 0.0)          # K = 1, L = 0: [U M I Z C W]
        f = SeasonForcing(; tspan = (0.0, 30.0), leaves0 = 1.0, weather = t -> 0.5, inoculations = [(10.0, 0.7)])
        sol = solve(MildewProblem(p, f, 0.0))
        @test sol(1, 5.0)[1] ≈ 1.0
        @test sol(1, 20.0)[1] ≈ exp(-0.35) rtol = 1e-10
        @test sol(1, 20.0)[3] ≈ 1 - exp(-0.35) rtol = 1e-10
        @test 10.0 in f.tstops
        two = solve(MildewProblem(p, SeasonForcing(; tspan = (0.0, 30.0), leaves0 = 1.0,
                                                   inoculations = [(15.0, 0.4), (10.0, 0.3)]), 0.0))
        @test two(1, 30.0)[1] ≈ exp(-0.7) rtol = 1e-10
    end

    @testset "inoculated leaves enter the latent stage; mature leaves at ε" begin
        p = MildewParameters(; b_c = 0.0, b_a = 0.0, tau_s = 5.0, K = 2, epsilon = 0.2, latent = 1e6, L = 2)
        f = SeasonForcing(; tspan = (0.0, 40.0), leaves0 = 1.0, inoculations = [(30.0, 1.0)])
        sol = solve(MildewProblem(p, f, 0.0))
        before, after = sol(1, 30.0 - 1e-9), sol(1, 30.0 + 1e-9)   # leaves age at 0.4 per day
        young_before, mature_before = before[1] + before[2], before[3]
        @test after[1] + after[2] ≈ young_before * exp(-1) rtol = 1e-6
        @test after[3] ≈ mature_before * exp(-0.2) rtol = 1e-6
        infected = young_before * (1 - exp(-1)) + mature_before * (1 - exp(-0.2))
        @test after[4] ≈ infected rtol = 1e-5                  # first latent stage
        @test leaf_area(sol, 1, 40.0) ≈ 1.0 rtol = 1e-9         # leaves are conserved
    end

    @testset "argument checks" begin
        @test_throws ArgumentError SeasonForcing(; tspan = (0.0, 10.0), inoculations = [(0.0, 1.0)])
        @test_throws ArgumentError SeasonForcing(; tspan = (0.0, 10.0), inoculations = [(10.0, 1.0)])
        @test_throws ArgumentError SeasonForcing(; tspan = (0.0, 10.0), inoculations = [(5.0, -1.0)])
        p = MildewParameters(; b_c = 0.1, b_a = 0.1)
        f = SeasonForcing(; tspan = (0.0, 10.0), leaves0 = 1.0, inoculations = [(5.0, 1.0)])
        @test_throws ArgumentError SpatialMildewProblem(p, Landscape(ones(2, 2); cellsize = 1.0),
                                                        Dispersal(; rho = 0.1, ell = 1.0), f, 0.1)
    end
end

@testset "temperature-dependent latent period" begin
    @test degree_day_rate(20.0; base = 5.0, reference = 20.0) == 1.0
    @test degree_day_rate(12.5; base = 5.0, reference = 20.0) ≈ 0.5
    @test degree_day_rate(3.0; base = 5.0, reference = 20.0) == 0.0
    @test_throws ArgumentError degree_day_rate(10.0; base = 20.0, reference = 20.0)

    # a constant rate d is the same model with latent period latent / d
    common = (; b_c = 0.3, b_a = 0.0, b_ext = 0.01, tau_s = 20.0, K = 2, epsilon = 0.1, L = 3, gamma_I = 0.05)
    slow = solve(MildewProblem(MildewParameters(; common..., latent = 8.0),
                               SeasonForcing(; tspan = (0.0, 60.0), leaves0 = 1.0, development = t -> 0.5), 0.0))
    ref = solve(MildewProblem(MildewParameters(; common..., latent = 16.0),
                              SeasonForcing(; tspan = (0.0, 60.0), leaves0 = 1.0), 0.0))
    @test slow(1, 60.0) ≈ ref(1, 60.0) rtol = 1e-6
    @test visible_incidence(slow, 1, 40.0) ≈ visible_incidence(ref, 1, 40.0) rtol = 1e-6
    # no development: inoculated leaves stay in the first latent stage until it resumes
    p = MildewParameters(; b_c = 0.0, b_a = 0.0, latent = 2.0, L = 2)   # [U M E1 E2 I Z C W]
    f = SeasonForcing(; tspan = (0.0, 30.0), leaves0 = 1.0, inoculations = [(5.0, 1.0)],
                      development = t -> t < 20 ? 0.0 : 1.0, tstops = [20.0])
    sol = solve(MildewProblem(p, f, 0.0))
    @test sol(1, 19.0)[3] ≈ 1 - exp(-1) rtol = 1e-9
    @test sol(1, 19.0)[5] == 0.0
    @test sol(1, 30.0)[5] > 0.9 * (1 - exp(-1))                      # sporulating soon after day 20
    # the spatial model uses the same development rate
    q = MildewParameters(; common..., latent = 8.0)
    g = SeasonForcing(; tspan = (0.0, 60.0), leaves0 = 1.0, development = t -> 0.5 + 0.01t)
    ssol = solve(SpatialMildewProblem(q, Landscape(ones(1, 2); cellsize = 1.0), Dispersal(; rho = 0.0, ell = 1.0), g, 0.0))
    nsol = solve(MildewProblem(q, g, 0.0))
    @test incidence(ssol, 1, 60.0) ≈ fill(incidence(nsol, 1, 60.0), 2) rtol = 1e-6
end

@testset "shedding of diseased leaves" begin
    @test MildewParameters(; b_c = 0.1, b_a = 1.0).shed == 0.0
    @test_throws ArgumentError MildewParameters(; b_c = 0.1, b_a = 1.0, shed = -0.01)
    @test_throws ArgumentError MildewParameters(; b_c = 0.1, b_a = 1.0, shed = Inf)

    # one inoculation on day 1 and no further infection: with s = t - 1 and I0 = 1 - e^-1,
    # I = I0 e^{-(γ+δ)s}, Z = I0 e^{-δs}(1 - e^{-γs}), leaves = 1 - I0 + I0 e^{-δs}
    γ, δ = 0.1, 0.05
    p = MildewParameters(; b_c = 0.0, b_a = 0.0, gamma_I = γ, shed = δ)   # [U M I Z C W]
    f = SeasonForcing(; tspan = (0.0, 40.0), leaves0 = 1.0, inoculations = [(1.0, 1.0)])
    sol = solve(MildewProblem(p, f, 0.0); abstol = 1e-12, reltol = 1e-10)
    I0 = 1 - exp(-1)
    for t in (5.0, 11.0, 30.0)
        s = t - 1
        @test sol(1, t)[3] ≈ I0 * exp(-(γ + δ) * s) rtol = 1e-7
        @test sol(1, t)[4] ≈ I0 * exp(-δ * s) * (1 - exp(-γ * s)) rtol = 1e-7
        @test leaf_area(sol, 1, t) ≈ 1 - I0 + I0 * exp(-δ * s) rtol = 1e-7
        @test visible_incidence(sol, 1, t) ≈ I0 * exp(-δ * s) / (1 - I0 + I0 * exp(-δ * s)) rtol = 1e-7
    end
    # so the visibly diseased share falls, where without shedding it stays at I0
    keep = solve(MildewProblem(MildewParameters(; b_c = 0.0, b_a = 0.0, gamma_I = γ), f, 0.0))
    @test visible_incidence(keep, 1, 30.0) ≈ I0 rtol = 1e-6
    @test visible_incidence(sol, 1, 30.0) < visible_incidence(sol, 1, 5.0)

    # the spatial model sheds the same way
    q = MildewParameters(; b_c = 0.3, b_a = 0.0, b_ext = 0.01, tau_s = 20.0, K = 2, epsilon = 0.1, latent = 8.0,
                         L = 3, gamma_I = 0.05, shed = 0.03)
    g = SeasonForcing(; tspan = (0.0, 60.0), leaves0 = 1.0)
    ssol = solve(SpatialMildewProblem(q, Landscape(ones(1, 2); cellsize = 1.0), Dispersal(; rho = 0.0, ell = 1.0), g, 0.0))
    nsol = solve(MildewProblem(q, g, 0.0))
    @test incidence(ssol, 1, 60.0) ≈ fill(incidence(nsol, 1, 60.0), 2) rtol = 1e-6
    @test visible_incidence(ssol, 1, 60.0) ≈ fill(visible_incidence(nsol, 1, 60.0), 2) rtol = 1e-6
end

@testset "shedding that rises with severity, and a shedding forcing" begin
    @test MildewParameters(; b_c = 0.1, b_a = 1.0).shed_severity == 0.0
    @test_throws ArgumentError MildewParameters(; b_c = 0.1, b_a = 1.0, shed_severity = -0.1)
    @test_throws ArgumentError MildewParameters(; b_c = 0.1, b_a = 1.0, shed_severity = NaN)

    # one inoculation on day 1, no further infection and no end to sporulation (γ = 0): with
    # healthy leaves h = 1 - I0 and δ = a I/(h + I), I' = -a I²/(h + I), so
    # h/I - log I = h/I0 - log I0 + a s with s = t - 1
    a = 0.2
    p = MildewParameters(; b_c = 0.0, b_a = 0.0, shed_severity = a)   # [U M I Z C W]
    f = SeasonForcing(; tspan = (0.0, 40.0), leaves0 = 1.0, inoculations = [(1.0, 1.0)])
    sol = solve(MildewProblem(p, f, 0.0); abstol = 1e-12, reltol = 1e-10)
    I0 = 1 - exp(-1)
    h = 1 - I0
    for t in (5.0, 20.0, 40.0)
        It = sol(1, t)[3]
        @test h / It - log(It) ≈ h / I0 - log(I0) + a * (t - 1) rtol = 1e-7
        @test leaf_area(sol, 1, t) ≈ h + It rtol = 1e-9
    end

    # the forcing multiplies both terms: s = 1/2 is the same as halving shed and shed_severity
    common = (; b_c = 0.3, b_a = 0.0, b_ext = 0.01, tau_s = 20.0, K = 2, epsilon = 0.1, latent = 8.0, L = 3,
              gamma_I = 0.05)
    g(s) = SeasonForcing(; tspan = (0.0, 60.0), leaves0 = 1.0, shedding = s)
    half = solve(MildewProblem(MildewParameters(; common..., shed = 0.02, shed_severity = 0.1), g(t -> 0.5), 0.0))
    ref = solve(MildewProblem(MildewParameters(; common..., shed = 0.01, shed_severity = 0.05), g(t -> 1.0), 0.0))
    @test half(1, 60.0) ≈ ref(1, 60.0) rtol = 1e-6
    # s = 0 switches shedding off; no shedding until day 20, then shedding
    none = solve(MildewProblem(MildewParameters(; common...), g(t -> 1.0), 0.0))
    off = solve(MildewProblem(MildewParameters(; common..., shed = 0.02, shed_severity = 0.1), g(t -> 0.0), 0.0))
    @test off(1, 60.0) ≈ none(1, 60.0) rtol = 1e-6
    late = solve(MildewProblem(MildewParameters(; common..., shed_severity = 0.5),
                               SeasonForcing(; tspan = (0.0, 60.0), leaves0 = 1.0, shedding = t -> t < 20 ? 0.0 : 1.0,
                                             tstops = [20.0]), 0.0))
    @test late(1, 19.0) ≈ none(1, 19.0) rtol = 1e-6
    @test leaf_area(late, 1, 60.0) < leaf_area(none, 1, 60.0)

    # the spatial model sheds the same way
    q = MildewParameters(; common..., shed = 0.01, shed_severity = 0.2)
    s = SeasonForcing(; tspan = (0.0, 60.0), leaves0 = 1.0, shedding = t -> 0.5 + 0.01t)
    ssol = solve(SpatialMildewProblem(q, Landscape(ones(1, 2); cellsize = 1.0), Dispersal(; rho = 0.0, ell = 1.0), s, 0.0))
    nsol = solve(MildewProblem(q, s, 0.0))
    @test visible_incidence(ssol, 1, 60.0) ≈ fill(visible_incidence(nsol, 1, 60.0), 2) rtol = 1e-6
    @test leaf_area(ssol, 1, 60.0) ≈ fill(leaf_area(nsol, 1, 60.0), 2) rtol = 1e-6
end
