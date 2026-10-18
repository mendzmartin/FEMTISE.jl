# Position-dependent effective mass (BenDaniel–Duke Hamiltonian)
#
#     H = -d/dx [ 1/(2 m*(x)) d/dx ] + V(x)
#
# for a square quantum well of width w: m* = m_w, V = 0 inside; m* = m_b,
# V = V_b outside (GaAs/AlGaAs-like parameters, atomic units). The weak form
# imposes the interface conditions [ψ] = 0 and [ψ'/m*] = 0 automatically, so the
# finite element levels must converge to the exact levels given by
#
#     even states:  (k/m_w) tan(k w/2) =  κ/m_b
#     odd states :  (k/m_w) cot(k w/2) = -κ/m_b
#
# with k = √(2 m_w ε) and κ = √(2 m_b (V_b − ε)). The mesh is aligned with the
# interfaces (the element size divides w/2 and the half-width of the box).

@testset "Position-dependent mass (BenDaniel–Duke)" begin
    m_well, m_barrier = 0.067, 0.092
    barrier_height = 230.0 / 27211.386          # 230 meV in hartree
    well_width = 200.0                          # bohr
    box_length = 3000.0                         # bohr

    effective_mass(z) = abs(z) < well_width / 2 ? m_well : m_barrier
    kinetic_coefficient(x) = 0.5 / effective_mass(x[1])
    band_profile(x) = abs(x[1]) < well_width / 2 ? 0.0 : barrier_height

    k(ε) = sqrt(2m_well * ε)
    κ(ε) = sqrt(2m_barrier * (barrier_height - ε))
    even_condition(ε) = k(ε) / m_well * sin(k(ε) * well_width / 2) -
                        κ(ε) / m_barrier * cos(k(ε) * well_width / 2)
    # the ground state is the only root of the even condition below 60 meV
    exact_ground = bisection(even_condition, 1e-6, 60.0 / 27211.386)

    model = CartesianDiscreteModel((-box_length / 2, box_length / 2), (480,))   # h = 6.25 bohr
    ε, = solve_on_mesh(model; p = kinetic_coefficient, q = band_profile, nev = 3)

    @test ε[1] > exact_ground                                 # upper bound
    @test (ε[1] - exact_ground) / exact_ground < 2e-6         # ≈ 9e-7 for this mesh
    @test count(<(barrier_height), ε) == 3                    # three bound states

    # Neglecting the mass jump changes the ground level by ~7 %: the test fails if
    # the position-dependent coefficient is silently ignored.
    ε_constant_mass, = solve_on_mesh(model; p = x -> 0.5 / m_well, q = band_profile, nev = 1)
    @test (ε_constant_mass[1] - exact_ground) / exact_ground > 0.05
end
