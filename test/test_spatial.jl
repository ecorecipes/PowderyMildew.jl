using LinearAlgebra: I, Diagonal, kron

# K₀(x) = ∫₀^∞ exp(-x cosh t) dt; the trapezoid rule is spectrally accurate here because the
# integrand is even in t and decays doubly exponentially.
function besselk0(x; T = 8.0, n = 4000)
    h = T / n
    return h * (exp(-x) / 2 + sum(exp(-x * cosh(k * h)) for k in 1:n))
end

# the Neumann graph Laplacian of a path of n cells
function path_laplacian(n)
    L = zeros(n, n)
    for i in 1:(n - 1)
        L[i, i + 1] = L[i + 1, i] = -1.0
        L[i, i] += 1.0
        L[i + 1, i + 1] += 1.0
    end
    return L
end

with_bc(p, b) = MildewParameters(b, p.b_a, p.b_ext, p.tau_s, p.K, p.epsilon, p.latent, p.L, p.gamma_I, p.c_w)

@testset "spatial model" begin
    p = MildewParameters(; b_c = 0.15, b_a = 3.0, tau_s = 25.0, K = 2, epsilon = 0.1, latent = 8.0, L = 2,
                         gamma_I = 0.02, c_w = 0.6)
    forcing(date; w = 1.0) = SeasonForcing(; tspan = (90.0, 280.0), production = FlushProduction([date], [7.0], [1.0]),
                                           release = ReleaseHazard(; rate = 0.3, optimum = 125.0, curvature = -0.003),
                                           weather = t -> w)
    ts = (100.0, 150.0, 200.0, 280.0)
    function agrees(ssol, stand_sols, seasons; rtol = 1e-6, atol = 1e-10)
        return all(isapprox(ssol(s, t)[:, j], stand_sols[j](s, t); rtol, atol)
                   for j in eachindex(stand_sols), s in 1:seasons, t in ts)
    end

    @testset "dispersal field matches a dense inverse" begin
        nx, ny, h, ell = 5, 4, 0.5, 1.7
        land = Landscape(ones(nx, ny); cellsize = h)
        A = I + (ell / h)^2 * (kron(Matrix(1.0I, ny, ny), path_laplacian(nx)) + kron(path_laplacian(ny), Matrix(1.0I, nx, nx)))
        s = reshape(collect(range(0.1, 2.0; length = nx * ny)) .^ 2, nx, ny)
        φ = dispersal_field(land, ell, s)
        @test vec(φ) ≈ A \ vec(s) rtol = 1e-12
        @test sum(φ) ≈ sum(s) rtol = 1e-12                       # spores are conserved
        @test dispersal_field(land, ell, fill(0.4, nx, ny)) ≈ fill(0.4, nx, ny) rtol = 1e-12
        @test dispersal_field(land, 0.0, s) ≈ s rtol = 1e-14       # ℓ = 0 deposits locally
    end

    @testset "dispersal field approaches the K₀ kernel" begin
        n, ell = 101, 5.0
        land = Landscape(ones(n, n); cellsize = 1.0)
        source = zeros(n, n)
        source[51, 51] = 1.0
        φ = dispersal_field(land, ell, source)
        for r in (10, 15)                                         # 2ℓ and 3ℓ from the source
            @test φ[51 + r, 51] ≈ besselk0(r / ell) / (2π * ell^2) rtol = 0.02
            @test φ[51, 51 + r] ≈ φ[51 + r, 51] rtol = 1e-10
        end
        @test besselk0(1.0) ≈ 0.42102443824070834 rtol = 1e-12
    end

    @testset "ρ = 0 gives independent stands" begin
        host = [0.3 0.0; 1.0 0.6; 0.8 0.2]
        land = Landscape(host; cellsize = 2.0)
        @test nstands(land) == 5
        dates = [100.0, 104.0, 108.0, 112.0, 116.0]
        weights = [1.0, 0.8, 0.9, 0.7, 1.1]
        F = [forcing(dates[j] + 2s; w = weights[j]) for j in 1:5, s in 1:2]
        C0 = [0.05, 0.0, 0.2, 0.1, 0.01]
        ssol = solve(SpatialMildewProblem(p, land, Dispersal(; rho = 0.0, ell = 3.0), F, C0))
        stands = [solve(MildewProblem(p, F[j, :], C0[j])) for j in 1:5]
        @test agrees(ssol, stands, 2)
        @test winter_inoculum(ssol) ≈ reduce(vcat, permutedims.(winter_inoculum.(stands))) rtol = 1e-6
        @test incidence(ssol, 2, 200.0) ≈ [incidence(s, 2, 200.0) for s in stands] rtol = 1e-6
        @test all(iszero, incidence(ssol, 2, 200.0)[[2]])          # no inoculum, no coupling
        g = gridded(land, incidence(ssol, 1, 280.0))
        @test isnan(g[1, 2]) && g[3, 1] == incidence(ssol, 1, 280.0)[3]
    end

    @testset "a uniform landscape stays uniform" begin
        for H in (1.0, 0.7)
            rho = 0.8
            land = Landscape(fill(H, 4, 3); cellsize = 1.0)
            ssol = solve(SpatialMildewProblem(p, land, Dispersal(; rho, ell = 2.5), [forcing(105.0), forcing(108.0)], 0.1))
            ref = solve(MildewProblem(with_bc(p, p.b_c * (1 - rho + rho * H)), [forcing(105.0), forcing(108.0)], 0.1))
            @test agrees(ssol, fill(ref, 12), 2)
            @test regional_incidence(ssol, 2, 200.0) ≈ incidence(ref, 2, 200.0) rtol = 1e-6
            @test regional_incidence(ssol, 2, 200.0; visible = true) ≈ visible_incidence(ref, 2, 200.0) rtol = 1e-6
        end
    end

    @testset "ℓ = 0 gives the local source" begin
        host = [0.2 0.9; 0.5 1.0]
        rho = 0.6
        land = Landscape(host; cellsize = 1.0)
        ssol = solve(SpatialMildewProblem(p, land, Dispersal(; rho, ell = 0.0), forcing(106.0), [0.1, 0.02, 0.3, 0.0]))
        stands = [solve(MildewProblem(with_bc(p, p.b_c * (1 - rho + rho * host[k])), forcing(106.0), c))
                  for (k, c) in zip(1:4, [0.1, 0.02, 0.3, 0.0])]
        @test agrees(ssol, stands, 1)
    end

    @testset "spread along a strip" begin
        land = Landscape(ones(30, 1); cellsize = 1.0)
        C0 = [1.0; zeros(29)]
        seasons = [forcing(105.0) for _ in 1:4]
        run(rho, ell) = solve(SpatialMildewProblem(p, land, Dispersal(; rho, ell), seasons, C0))
        isolated, near, far = run(0.0, 2.0), run(0.3, 1.0), run(0.3, 3.0)
        @test all(iszero, incidence(isolated, 4, 280.0)[2:end])    # no escape, no spread
        inc = incidence(far, 1, 280.0)
        @test inc[2] > 0 && issorted(inc; rev = true)              # falls away from the source
        @test incidence(far, 4, 280.0)[15] > incidence(near, 4, 280.0)[15]
        @test incidence(far, 4, 280.0)[10] > incidence(far, 1, 280.0)[10]  # it advances
        @test all(s -> minimum(s(4, 280.0)) > -1e-10, (isolated, near, far))
        @test size(winter_inoculum(far)) == (30, 4)
        @test all(≥(0), leaf_area(far, 2, 150.0))
        @test visible_incidence(far, 2, 280.0; latent_visible = 1.0) ≈ incidence(far, 2, 280.0) rtol = 1e-10
    end

    @testset "argument checks" begin
        @test_throws ArgumentError Landscape([1.0 -0.1]; cellsize = 1.0)
        @test_throws ArgumentError Landscape([1.0 NaN]; cellsize = 1.0)
        @test_throws ArgumentError Landscape(zeros(2, 2); cellsize = 1.0)
        @test_throws ArgumentError Landscape(ones(2, 2); cellsize = 0.0)
        @test_throws ArgumentError Dispersal(; rho = 1.2, ell = 1.0)
        @test_throws ArgumentError Dispersal(; rho = 0.5, ell = -1.0)
        @test_throws ArgumentError Dispersal(; rho = 0.5, ell = Inf)
        land = Landscape(ones(2, 2); cellsize = 1.0)
        d = Dispersal(; rho = 0.5, ell = 1.0)
        @test_throws ArgumentError SpatialMildewProblem(p, land, d, forcing(100.0), [0.1, 0.1])
        @test_throws ArgumentError SpatialMildewProblem(p, land, d, forcing(100.0), -0.1)
        @test_throws ArgumentError SpatialMildewProblem(p, land, d, fill(forcing(100.0), 3, 1), 0.1)
        bad = [forcing(100.0); forcing(100.0); forcing(100.0);
               SeasonForcing(; tspan = (90.0, 270.0), production = FlushProduction([100.0], [7.0], [1.0]),
                             release = ReleaseHazard(; rate = 0.3, optimum = 125.0, curvature = -0.003))]
        @test_throws ArgumentError SpatialMildewProblem(p, land, d, reshape(bad, 4, 1), 0.1)
        @test_throws ArgumentError gridded(land, [1.0, 2.0])
        @test_throws ArgumentError dispersal_field(land, 1.0, ones(3, 3))
        sol = solve(SpatialMildewProblem(p, land, d, forcing(100.0), 0.1))
        @test_throws ArgumentError sol(2, 100.0)
        @test_throws ArgumentError sol(1, 300.0)
        @test_throws ArgumentError regional_incidence(sol, 1, 200.0; latent_visible = 2.0)
    end
end
