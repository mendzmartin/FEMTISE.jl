# Two particles of different masses in a common harmonic trap.
#
# On a 2D mesh with coordinates (x₁, x₂), the anisotropic kinetic form
# p₁ ∂₁u ∂₁v + p₂ ∂₂u ∂₂v with p₁ = 1/2 and p₂ = 1/(2M₂) describes a particle of
# mass 1 and one of mass M₂ in V = (x₁² + x₂²)/2. The problem is separable:
#
#     ε(n₁, n₂) = (n₁ + 1/2) + (n₂ + 1/2)/√M₂.
#
# M₂ = 3 makes the spectrum non-degenerate, so each level can be compared
# individually.

@testset "Two particles with different masses (2D)" begin
    heavy_mass = 3.0
    model = CartesianDiscreteModel((-7.0, 7.0, -7.0, 7.0), (40, 40))

    ε, = solve_on_mesh(model; p = electron_kinetic, p2 = x -> 0.5 / heavy_mass,
                       q = harmonic_2d, nev = 5)

    separable_levels = sort([(n1 + 0.5) + (n2 + 0.5) / sqrt(heavy_mass)
                             for n1 in 0:5 for n2 in 0:8])
    exact = separable_levels[1:5]

    @test all(ε .> exact)                              # upper bounds
    @test maximum((ε .- exact) ./ exact) < 1e-3

    # Equal masses: the isotropic method must recover the (n + 1)-fold degeneracy.
    ε_equal, = solve_on_mesh(model; p = electron_kinetic, q = harmonic_2d, nev = 3)
    @test ε_equal ≈ [1.0, 2.0, 2.0] rtol = 1e-3
    @test ε_equal[3] - ε_equal[2] < 1e-8              # x₁ ↔ x₂ symmetry is exact
end
