using Distributions: Binomial, logpdf
using OrdinaryDiffEqVerner: Vern9
using SciMLBase: ODEProblem

# Reference for α = 0 and S + ΣI = 1: with u(t) = ∫₀ᵗ S, Iᵢ(t) = Iᵢ(0) exp(βᵢ u) and
# u' = 1 - Σᵢ Iᵢ(0) exp(βᵢ u), a scalar ODE independent of the package's formulation.
function reference_infected(beta, I0, t)
    f(u, p, s) = 1 - sum(I0[i] * exp(beta[i] * u) for i in eachindex(I0))
    u = solve(ODEProblem(f, 0.0, (0.0, t)), Vern9(); abstol = 1e-15, reltol = 1e-14).u[end]
    return [I0[i] * exp(beta[i] * u) for i in eachindex(I0)]
end

@testset "Hamelin et al. 2016" begin
    @testset "one species is logistic" begin
        b, I0 = 0.15, 1e-3
        p = HamelinParameters(; beta = [b], chi = [0.01], season_length = 120)
        sol = solve(HamelinProblem(p, [I0]))
        for t in (0.0, 10.0, 47.5, 120.0)
            logistic = I0 * exp(b * t) / (1 - I0 + I0 * exp(b * t))
            @test sol(1, t)[2] ≈ logistic rtol = 1e-7
            @test incidence(sol, 1, t) ≈ logistic rtol = 1e-7
        end
        @test frequency(sol, 1, 60.0) == 1.0
    end

    @testset "two species against the cumulative-exposure reference" begin
        beta, I0 = [0.196, 0.098], [0.024, 0.976] .* 4e-4
        p = HamelinParameters(; beta, chi = [1e-5, 3e-3])
        sol = solve(HamelinProblem(p, I0))
        for t in (5.0, 40.0, 80.0, 200.0)
            ref = reference_infected(beta, I0, t)
            @test sol(1, t)[2:3] ≈ ref rtol = 1e-7
            @test frequency(sol, 1, t) ≈ ref[1] / sum(ref) rtol = 1e-7
        end
        # the share of the faster species rises through the season
        @test frequency(sol, 1, 0.0) ≈ 0.024
        @test frequency(sol, 1, 200.0) > 0.85
    end

    @testset "equal transmission leaves the species share unchanged" begin
        p = HamelinParameters(; beta = [0.2, 0.2], chi = [0.01, 0.02], alpha = 0.03)
        sol = solve(HamelinProblem(p, [3e-4, 1e-4]; S = 0.9))
        @test all(frequency(sol, 1, t) ≈ 0.75 for t in 0:10:200)
    end

    @testset "conservation and positivity with removal" begin
        p = HamelinParameters(; beta = [0.3, 0.1], chi = [0.05, 0.5], alpha = 0.04, S0 = 2.5)
        sol = solve(HamelinProblem(p, [1e-3, 2e-3]; seasons = 4))
        for s in 1:4, t in 0:25:200
            u = sol(s, t)
            @test sum(u) ≈ 2.5 rtol = 1e-8
            @test all(≥(-1e-12), u)
        end
        # the annual map renews the host pool and bounds infection by S0
        for s in 1:3
            u = between_season(p, sol(s, 200.0))
            @test u ≈ sol(s + 1, 0.0)
            @test sum(u) ≈ 2.5
            @test u[end] == 0
        end
    end

    @testset "annual map" begin
        p = HamelinParameters(; beta = [0.2, 0.1], chi = [2.0, 0.5], S0 = 1.0)
        @test between_season(p, [0.4, 0.0, 0.0, 0.6]) == [1.0, 0.0, 0.0, 0.0]
        u = between_season(p, [0.1, 0.3, 0.6, 0.0])
        H = 2.0 * 0.3 + 0.5 * 0.6
        @test u[1] ≈ exp(-H)
        @test u[2:3] ≈ (1 - exp(-H)) .* [0.6, 0.3] ./ H
        # small pressure: no cancellation in 1 - exp(-H)
        tiny = between_season(p, [1.0, 1e-18, 0.0, 0.0])
        @test tiny[2] ≈ 2e-18 rtol = 1e-12
        @test tiny[1] + tiny[2] ≈ 1
        @test_throws ArgumentError between_season(p, [0.5, 0.5])
    end

    @testset "single-species limit of the two-species model" begin
        one = solve(HamelinProblem(HamelinParameters(; beta = [0.2], chi = [0.3], alpha = 0.01), [1e-3]; seasons = 3))
        two = solve(HamelinProblem(HamelinParameters(; beta = [0.2, 0.5], chi = [0.3, 0.9], alpha = 0.01), [1e-3, 0.0];
                                   seasons = 3))
        for s in 1:3, t in (0.0, 30.0, 200.0)
            @test two(s, t)[[1, 2, 4]] ≈ one(s, t) rtol = 1e-8
            @test two(s, t)[3] == 0
        end
    end

    @testset "coexistence equilibrium reproduces Table 2" begin
        # Table 2: (I(0), β1, β2, F0, Fτ, χ1, χ2), β per day. χ is a function of the rounded
        # β, F0 and I(0) alone; the 2009 fit with I(0) = 5.6e-5 ends so close to F = 1 that its
        # χ2 = H (1 - F0) / (1 - Fτ) is too sensitive to the three-figure β to check at 10 %.
        table2 = [(1e-3, 0.181, 0.083, 0.024, 0.878, 2.75e-5, 0.008, 0.02),
                  (4e-4, 0.196, 0.098, 0.024, 0.877, 1.12e-5, 0.003, 0.12),
                  (1.4e-4, 0.312, 0.089, 0.020, 0.994, 2.82e-6, 0.025, 0.05),
                  (5.6e-5, 0.316, 0.112, 0.020, 0.994, 1.59e-6, 0.010, 0.35)]
        for (I, b1, b2, F0, Fτ, c1, c2, tol) in table2
            pt = HamelinParameters(; beta = [b1, b2], chi = [0.0, 0.0])
            chi = coexistence_chi(pt, [F0, 1 - F0] .* I)
            st = solve(HamelinProblem(HamelinParameters(; beta = pt.beta, chi), [F0, 1 - F0] .* I))
            @test frequency(st, 1, 200.0) ≈ Fτ atol = 0.006
            @test chi ≈ [c1, c2] rtol = tol
            @test chi[2] ≈ -log1p(-I) * (1 - F0) / (1 - frequency(st, 1, 200.0)) rtol = 1e-6
        end
        # 2014, lower initial incidence, in detail
        p = HamelinParameters(; beta = [0.196, 0.098], chi = [0.0, 0.0])
        I0 = [0.024, 0.976] .* 4e-4
        chi = coexistence_chi(p, I0)
        @test chi[1] ≈ 1.12e-5 rtol = 0.05
        @test round(chi[2]; sigdigits = 1) == 0.003
        # and the start of the season is then a fixed point of the annual cycle
        q = HamelinParameters(; beta = p.beta, chi)
        sol = solve(HamelinProblem(q, I0; seasons = 3))
        for s in 2:3
            @test sol(s, 0.0)[2:3] ≈ I0 rtol = 1e-6
        end
        # chi only restates the assumed initial incidence: the season ends with every leaf infected
        @test incidence(sol, 1, 200.0) > 1 - 1e-6
        @test_throws ArgumentError coexistence_chi(p, [4e-4, 0.0])
        @test_throws ArgumentError coexistence_chi(p, [0.6, 0.6])
    end

    @testset "binomial likelihood of species frequencies" begin
        p = HamelinParameters(; beta = [0.196, 0.098], chi = [1.12e-5, 3e-3])
        sol = solve(HamelinProblem(p, [0.024, 0.976] .* 4e-4; seasons = 2))
        season, t, m, y = [1, 1, 2], [20.0, 90.0, 150.0], [30, 28, 0], [2, 20, 0]
        expected = logpdf(Binomial(30, frequency(sol, 1, 20.0)), 2) + logpdf(Binomial(28, frequency(sol, 1, 90.0)), 20)
        @test frequency_loglikelihood(sol, season, t, m, y) ≈ expected
        @test frequency_loglikelihood(sol, season, t, m, m .- y; species = 2) ≈ expected
        none = solve(HamelinProblem(p, [0.0, 0.0]))
        @test isnan(frequency(none, 1, 10.0))
        @test frequency_loglikelihood(none, [1], [10.0], [5], [1]) == -Inf
        @test_throws ArgumentError frequency_loglikelihood(sol, [1], [10.0], [5], [6])
        @test_throws ArgumentError frequency_loglikelihood(sol, [1, 1], [10.0], [5], [1])
    end

    @testset "argument checks" begin
        @test_throws ArgumentError HamelinParameters(; beta = Float64[], chi = Float64[])
        @test_throws ArgumentError HamelinParameters(; beta = [0.1, 0.2], chi = [0.1])
        @test_throws ArgumentError HamelinParameters(; beta = [-0.1], chi = [0.1])
        @test_throws ArgumentError HamelinParameters(; beta = [0.1], chi = [0.1], alpha = -1)
        @test_throws ArgumentError HamelinParameters(; beta = [0.1], chi = [0.1], S0 = 0)
        p = HamelinParameters(; beta = [0.1], chi = [0.1])
        @test_throws ArgumentError HamelinProblem(p, [0.1, 0.1])
        @test_throws ArgumentError HamelinProblem(p, [-0.1])
        @test_throws ArgumentError HamelinProblem(p, [0.1]; seasons = 0)
        sol = solve(HamelinProblem(p, [0.01]))
        @test_throws ArgumentError sol(2, 1.0)
        @test_throws ArgumentError sol(1, 250.0)
        @test_throws ArgumentError frequency(sol, 1, 1.0; species = 2)
        @test nspecies(p) == 1
    end
end
