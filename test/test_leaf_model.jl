using Distributions: Normal, cdf

@testset "leaf-level model" begin
    @testset "nests Hamelin et al. 2016 (one species)" begin
        β, α, χ, I0, τ = 0.2, 0.03, 1.5, 1e-3, 150.0
        hamelin = solve(HamelinProblem(HamelinParameters(; beta = [β], chi = [χ], alpha = α, season_length = τ), [I0];
                                       seasons = 3))
        b_a = 0.7
        p = MildewParameters(; b_c = β, b_a, gamma_I = α, c_w = χ / b_a)
        f = SeasonForcing(; tspan = (0.0, τ), leaves0 = 1.0, pulse = true)
        sol = solve(MildewProblem(p, fill(f, 3), -log1p(-I0) / b_a))
        iI, iZ = 3, 4    # [U; M; I; Z; C; W] with K = 1, L = 0
        for s in 1:3, t in (0.0, 15.0, 60.0, τ)
            h = hamelin(s, t)
            u = sol(s, t)
            @test [u[1] + u[2], u[iI], u[iZ]] ≈ h rtol = 1e-7
        end
        @test winter_inoculum(sol) ≈ [p.c_w * sol(s, τ)[iI] for s in 1:3]
    end

    @testset "a fast release converges to the pulse" begin
        p = MildewParameters(; b_c = 0.15, b_a = 2.0, gamma_I = 0.02)
        pulse = solve(MildewProblem(p, SeasonForcing(; tspan = (0.0, 60.0), leaves0 = 1.0, pulse = true), 0.01))
        err(r) = abs(incidence(solve(MildewProblem(p, SeasonForcing(; tspan = (0.0, 60.0), leaves0 = 1.0,
                                                                     release = t -> r), 0.01)), 1, 30.0) -
                     incidence(pulse, 1, 30.0))
        @test err(500.0) < err(50.0) / 5
        @test err(500.0) < 1e-3
    end

    @testset "Erlang ageing of uninfected leaves" begin
        p = MildewParameters(; b_c = 0.0, b_a = 0.0, tau_s = 20.0, K = 3)
        sol = solve(MildewProblem(p, SeasonForcing(; tspan = (0.0, 100.0), leaves0 = 1.0), 0.0))
        r = 3 / 20
        for t in (5.0, 20.0, 45.0)
            young = sum(exp(-r * t) * (r * t)^k / factorial(k) for k in 0:2)
            @test sum(sol(1, t)[1:3]) ≈ young rtol = 1e-7
            @test sol(1, t)[4] ≈ 1 - young rtol = 1e-7
        end
    end

    @testset "Erlang latent period" begin
        p = MildewParameters(; b_c = 0.0, b_a = 1.0, latent = 7.0, L = 2)
        sol = solve(MildewProblem(p, SeasonForcing(; tspan = (0.0, 40.0), leaves0 = 1.0, pulse = true), 0.1))
        infected = -expm1(-0.1)
        r = 2 / 7
        iI = 2 + 2 + 1    # [U; M; E₁; E₂; I; Z; C; W]
        for t in (3.0, 7.0, 20.0)
            @test sol(1, t)[iI] ≈ infected * (1 - exp(-r * t) * (1 + r * t)) rtol = 1e-7
            @test incidence(sol, 1, t) ≈ infected rtol = 1e-8
            @test visible_incidence(sol, 1, t) ≈ sol(1, t)[iI] / leaf_area(sol, 1, t) rtol = 1e-10
            @test visible_incidence(sol, 1, t; latent_visible = 1) ≈ incidence(sol, 1, t) rtol = 1e-10
        end
    end

    @testset "leaf production, conservation and positivity" begin
        flush = FlushProduction([120.0, 190.0], [5.0, 8.0], [0.8, 0.2])
        hazard = ReleaseHazard(; rate = 0.2, optimum = 140.0, curvature = -0.003, intercept = 1.0)
        temperature = DailySeries(fill(18.0, 365))
        weather = t -> beta_response(temperature(t))
        p = MildewParameters(; b_c = 0.3, b_a = 0.5, b_ext = 1e-5, tau_s = 20.0, K = 3, epsilon = 0.05, latent = 8.0,
                             L = 3, gamma_I = 0.05, c_w = 0.02)
        f = SeasonForcing(; tspan = (90.0, 280.0), production = flush, release = hazard, weather, gate = t -> t > 240,
                          tstops = 90:280)
        sol = solve(MildewProblem(p, fill(f, 3), 1.0))
        for s in 1:3, t in (100.0, 125.0, 150.0, 200.0, 280.0)
            expected = 0.8 * cdf(Normal(120, 5), t) + 0.2 * cdf(Normal(190, 8), t) -
                       (0.8 * cdf(Normal(120, 5), 90) + 0.2 * cdf(Normal(190, 8), 90))
            @test leaf_area(sol, s, t) ≈ expected rtol = 1e-6 atol = 1e-12
            @test all(≥(-1e-12), sol(s, t))
            @test 0 ≤ visible_incidence(sol, s, t) ≤ incidence(sol, s, t) ≤ 1
        end
        # the gate collects sporulating leaf-days after day 240; I is entry 8, W the last
        u = sol(1, 280.0)
        ts = range(240.0, 280.0; length = 4001)
        integral = (sum(sol(1, t)[8] for t in ts) - (sol(1, 240.0)[8] + u[8]) / 2) * step(ts)
        @test u[end] ≈ integral rtol = 1e-5
        @test winter_inoculum(sol)[1] ≈ 0.02 * u[end]
        @test sol(2, 90.0)[end - 1] ≈ winter_inoculum(sol)[1]
    end

    @testset "forcing helpers" begin
        d = DailySeries([1.0, 2.0, 3.0]; first_day = 10)
        @test d(10.0) == 1.0 && d(10.99) == 1.0 && d(11.0) == 2.0 && d(12.5) == 3.0
        @test d(5.0) == 1.0 && d(20.0) == 3.0
        @test beta_response(21.5) ≈ 1
        @test beta_response(5.0) == 0 && beta_response(30.0) == 0 && beta_response(35.0) == 0
        @test 0 < beta_response(15.0) < 1
        @test beta_response(25.0; Topt = 25.0, Tmax = 34.0) ≈ 1
        g = ReleaseHazard(; rate = 2.0, optimum = 140.0, curvature = -0.01, intercept = 0.5)
        @test g(140.0) ≈ 2.0 / (1 + exp(-0.5))
        @test g(120.0) < g(140.0) && g(160.0) ≈ g(120.0)
        @test FlushProduction([100.0], [4.0], [1.0])(100.0) ≈ 1 / (sqrt(2π) * 4)
    end

    @testset "argument checks" begin
        @test_throws ArgumentError MildewParameters(; b_c = -1, b_a = 0)
        @test_throws ArgumentError MildewParameters(; b_c = 0, b_a = 0, epsilon = 2)
        @test_throws ArgumentError MildewParameters(; b_c = 0, b_a = 0, K = 0)
        @test_throws ArgumentError MildewParameters(; b_c = 0, b_a = 0, latent = 5.0)
        @test_throws ArgumentError MildewParameters(; b_c = 0, b_a = 0, L = 2)
        @test_throws ArgumentError MildewParameters(; b_c = 0, b_a = 0, tau_s = 0)
        @test_throws ArgumentError SeasonForcing(; tspan = (10, 5))
        @test_throws ArgumentError FlushProduction([1.0], [0.0], [1.0])
        @test_throws ArgumentError ReleaseHazard(; rate = 1, optimum = 140, curvature = 0.01)
        @test_throws ArgumentError beta_response(20.0; Tmin = 25.0)
        p = MildewParameters(; b_c = 0.1, b_a = 1.0)
        @test_throws ArgumentError MildewProblem(p, SeasonForcing[], 1.0)
        @test_throws ArgumentError MildewProblem(p, SeasonForcing(; tspan = (0, 10)), -1.0)
        sol = solve(MildewProblem(p, SeasonForcing(; tspan = (0, 10), leaves0 = 1.0), 0.0))
        @test_throws ArgumentError sol(2, 1.0)
        @test_throws ArgumentError sol(1, 11.0)
        @test_throws ArgumentError visible_incidence(sol, 1, 1.0; latent_visible = 2)
    end
end
