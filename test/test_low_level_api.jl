# Low-level interface: correctness of the Galerkin eigensolver.
#
# Reference problem: 1D harmonic oscillator V = x²/2 in Ω = (-10, 10) with
# Dirichlet walls. Its low-lying levels coincide with εₙ = n + 1/2 to better
# than 1e-20, so they serve as exact values.
#
# Properties checked (see the "Theory" page of the documentation):
#   1. Variational upper bound: εₙʰ ≥ εₙ for every level (exact quadrature).
#   2. A priori convergence rate: the eigenvalue error decreases as h^(2ℓ) for
#      Lagrange elements of degree ℓ.
#   3. The eigenvectors are orthonormal in L².
#   4. `eigen_problem` + `solve` and `eigen_values_and_eigen_vectors` agree.

@testset "Low-level interface (1D harmonic oscillator)" begin
    domain = (-10.0, 10.0)
    exact = [n + 0.5 for n in 0:5]

    @testset "Upper bounds and h^(2ℓ) convergence, ℓ = $ℓ" for (ℓ, coarse_cells) in
                                                               ((1, 80), (2, 40), (3, 20))
        ε_coarse, = solve_on_mesh(CartesianDiscreteModel(domain, (coarse_cells,));
                                  p = electron_kinetic, q = harmonic_1d, order = ℓ)
        ε_fine, = solve_on_mesh(CartesianDiscreteModel(domain, (2coarse_cells,));
                                p = electron_kinetic, q = harmonic_1d, order = ℓ)
        error_coarse = ε_coarse[1] - exact[1]
        error_fine = ε_fine[1] - exact[1]

        @test all(ε_coarse .> exact)          # upper bounds on the coarse mesh
        @test all(ε_fine .> exact)            # ... and on the fine mesh
        @test ε_fine[1] < ε_coarse[1]         # monotone under refinement
        observed_rate = log2(error_coarse / error_fine)
        @test observed_rate ≈ 2ℓ atol = 0.3
    end

    @testset "Orthonormal eigenvectors" begin
        ε, ϕ, _, dΩ, _ = solve_on_mesh(CartesianDiscreteModel(domain, (200,));
                                       p = electron_kinetic, q = harmonic_1d)
        @test all(norm_l2(ψ, dΩ) ≈ 1.0 for ψ in ϕ)
        @test maximum(orthogonality_check(ϕ, dΩ)) < 1e-8
        @test eigenstates_normalization(ϕ, dΩ) ≈ ones(length(ϕ))
    end

    @testset "eigen_problem + solve agree with eigen_values_and_eigen_vectors" begin
        model = CartesianDiscreteModel(domain, (100,))
        _, dΩ, _, _ = measures(model, 6, "boundary")
        V, U = fe_spaces(model, ReferenceFE(lagrangian, Float64, 2); BC_type = "Dirichlet",
                         BC_data = (zero(ComplexF64), "boundary"))
        a, b = bilineal_forms(electron_kinetic, harmonic_1d, x -> 1.0, dΩ)
        problem = eigen_problem(a, b, U, V; nev = 4, sigma = 0.0, tol = 1e-12, maxiter = 1000)
        ε_two_steps, _ = FEMTISE.solve(problem)

        ε_one_call, = solve_on_mesh(model; p = electron_kinetic, q = harmonic_1d, nev = 4)
        @test sort(real.(ε_two_steps)) ≈ ε_one_call rtol = 1e-10
        A, B = problem.op.hamiltonian, problem.op.energy
        relative_asymmetry(M) = norm(M - transpose(M)) / norm(M)
        @test relative_asymmetry(A) < 1e-12                # A = Aᵀ (Hermitian form)
        @test relative_asymmetry(B) < 1e-12                # B = Bᵀ (mass matrix)
    end
end
