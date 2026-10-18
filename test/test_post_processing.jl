# Post-processing observables on finite element eigenfunctions.
#
# Reference: ground and first excited states of the 1D harmonic oscillator
# (V = x²/2, atomic units), whose observables are known in closed form:
#
#   ⟨x⟩ = 0,   ⟨x²⟩ₙ = n + 1/2,
#   Shannon entropy of ρ₀ = e^{-x²}/√π:   S = ½ log₂(π e)   [bits],
#   Rényi entropy of order 2 of ρ₀:       S₂ = -log₂ ∫ρ₀² = ½ log₂(2π).

@testset "Post-processing" begin
    model = CartesianDiscreteModel((-10.0, 10.0), (200,))
    ε, ϕ, Ω, dΩ, U = solve_on_mesh(model; p = electron_kinetic, q = harmonic_1d, nev = 2)

    @testset "Densities are normalized" begin
        ρ = density(ϕ)
        @test all(sum(∫(ρn) * dΩ) ≈ 1.0 for ρn in ρ)
    end

    @testset "Moments of the coordinate" begin
        first_moment = coord_first_moment(ϕ, U, Ω, dΩ, 1)
        second_moment = coord_second_moment(ϕ, U, Ω, dΩ, 1)
        @test all(abs.(first_moment) .< 1e-8)
        # discretization error of the moments on this mesh: ~1e-6 (relative)
        @test second_moment ≈ [0.5, 1.5] rtol = 1e-5
        @test coord_standar_deviation(ϕ, U, Ω, dΩ, 1) ≈ sqrt.([0.5, 1.5]) rtol = 1e-5
    end

    @testset "Shannon entropy" begin
        S = time_indep_entropy(ϕ[1:1], U, dΩ)
        @test S[1] ≈ 0.5 * log2(π * ℯ) atol = 1e-6
    end

    @testset "Rényi entropy (known issue in v0.1.0)" begin
        # time_indep_entropy evaluates ρ^(α/ln 2) instead of ρ^α (see CHANGELOG).
        # When the bug is fixed this test will report "unexpected pass": replace
        # @test_broken by @test.
        S2 = time_indep_entropy(ϕ[1:1], U, dΩ, 2.0)
        @test_broken S2[1] ≈ 0.5 * log2(2π) atol = 1e-6
    end
end
