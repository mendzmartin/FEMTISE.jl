# Features introduced in v0.2.0, each checked against an exact result.
#
#   * solve_tise with arbitrary masses: oscillator of mass M, εₙ = (n + ½)/√M.
#   * tabulated potentials: tabulated V = x²/2 (1D) and (x² + y²)/2 (2D).
#   * sweep_eigenvalues with one and two parameters: ε₀ = ω/2.
#   * autocorrelation: coherent state of the oscillator (displaced ground state),
#     whose weights |cₙ|² follow a Poisson distribution and whose autocorrelation
#     is C(t) = e^{-it/2} exp(α²(e^{-it} − 1)); the Fourier transform of C(t) must
#     reproduce the broadened stick spectrum.
#   * pauli_fierz_potential + 3D solver: quadratic nucleus–electron–photon model,
#     solved exactly by its normal modes.

@testset "Masses (solve_tise)" begin
    model = CartesianDiscreteModel((-10.0, 10.0), (200,))
    M = 4.0
    sol = solve_tise(model, harmonic_1d; masses = M, nev = 4)
    exact = [(n + 0.5) / sqrt(M) for n in 0:3]
    @test sol.energies ≈ exact rtol = 1e-5
    @test all(sol.energies .> exact)                       # upper bounds
    @test kinetic_coefficients((1.0, 2.0), 2) == (0.5, 0.25)
    @test_throws ArgumentError kinetic_coefficients((1.0, 2.0), 3)
    @test_throws ArgumentError kinetic_coefficients(-1.0, 1)
end

@testset "Tabulated potentials" begin
    model = CartesianDiscreteModel((-10.0, 10.0), (200,))
    x = collect(-10.0:0.1:10.0)
    exact = [n + 0.5 for n in 0:3]
    linear = solve_tise(model, tabulated_potential(x, 0.5 .* x .^ 2); nev = 4)
    cubic = solve_tise(model, tabulated_potential(x, 0.5 .* x .^ 2; interpolation = :cubic); nev = 4)
    @test maximum(abs.(linear.energies .- exact) ./ exact) < 3e-3     # error ∝ Δ²
    @test maximum(abs.(cubic.energies .- exact) ./ exact) < 1e-5      # error ∝ Δ⁴

    # the same table read from a data file
    mktempdir() do folder
        path = joinpath(folder, "potential.dat")
        open(path, "w") do io
            println(io, "# x   V(x)")
            for xi in x
                println(io, xi, " ", 0.5xi^2)
            end
        end
        from_file = solve_tise(model, tabulated_potential(path; interpolation = :cubic); nev = 4)
        @test from_file.energies ≈ cubic.energies rtol = 1e-12
    end

    # 2D bilinear table
    grid = collect(-6.0:0.1:6.0)
    table = [0.5 * (a^2 + b^2) for a in grid, b in grid]
    model2 = CartesianDiscreteModel((-6.0, 6.0, -6.0, 6.0), (40, 40))
    sol2 = solve_tise(model2, tabulated_potential(grid, grid, table); nev = 3)
    @test sol2.energies ≈ [1.0, 2.0, 2.0] rtol = 2e-3

    # evaluation outside the table is an error unless flat extrapolation is requested
    V = tabulated_potential(x, 0.5 .* x .^ 2)
    @test_throws DomainError V(11.0)
    @test tabulated_potential(x, 0.5 .* x .^ 2; extrapolation = :flat)(11.0) ≈ 50.0
    @test_throws ArgumentError tabulated_potential([0.0, 2.0, 1.0], [0.0, 1.0, 2.0])
end

@testset "Parameter sweeps (sweep_eigenvalues)" begin
    model = CartesianDiscreteModel((-10.0, 10.0), (200,))
    ω = [0.5, 1.0, 2.0]
    one_parameter = sweep_eigenvalues(model, (x, ω) -> 0.5ω^2 * x[1]^2, ω; nev = 2)
    @test size(one_parameter) == (2, 3)
    @test one_parameter[1, :] ≈ ω ./ 2 rtol = 1e-5

    centres = [0.0, 1.0]
    two_parameters = sweep_eigenvalues(model, (x, ω, x₀) -> 0.5ω^2 * (x[1] - x₀)^2, ω, centres; nev = 2)
    @test size(two_parameters) == (2, 3, 2)
    @test two_parameters[1, :, 1] ≈ two_parameters[1, :, 2] rtol = 1e-6   # shift does not matter
    @test two_parameters[2, :, 2] ≈ 1.5 .* ω rtol = 1e-5
end

@testset "Autocorrelation and spectrum" begin
    model = CartesianDiscreteModel((-10.0, 10.0), (200,))
    sol = solve_tise(model, harmonic_1d; nev = 30)
    x₀ = 2.0                                     # displacement; α² = x₀²/2
    α² = x₀^2 / 2
    ψ₀(x) = exp(-0.5 * (x[1] - x₀)^2) / π^0.25
    c = expansion_coefficients(ψ₀, sol.states, sol.dΩ)
    poisson = [Float64(exp(-α²) * α²^n / factorial(big(n))) for n in 0:29]
    @test sum(abs2, c) ≈ 1.0 atol = 1e-8         # 30 states span the wave packet
    @test abs2.(c) ≈ poisson atol = 1e-5

    t = collect(0.0:0.05:60.0)
    C = autocorrelation(c, sol.energies, t)
    exact = [cis(-τ / 2) * exp(α² * (cis(-τ) - 1)) for τ in t]
    @test maximum(abs.(C .- exact)) < 1e-3
    @test autocorrelation(c, sol.energies, 0.0) ≈ sum(abs2, c)

    # Fourier transform of the sampled C(t) = broadened stick spectrum
    ω = collect(0.0:0.02:6.0)
    for (window, lineshape) in ((:exponential, :lorentzian), (:gaussian, :gaussian))
        analytic = autocorrelation_spectrum(c, sol.energies, ω; width = 0.2, lineshape = lineshape)
        numeric = fourier_transform(t, C, ω; window = window, width = 0.2)
        @test maximum(abs.(numeric .- analytic)) < 1e-3 * maximum(analytic)
    end
    # the spectrum is normalized: ∫σ dω = Σ|cₙ|² (Gaussian lines, wide grid)
    wide = collect(-2.0:0.01:35.0)
    σ = autocorrelation_spectrum(c, sol.energies, wide; width = 0.2, lineshape = :gaussian)
    @test sum(σ) * 0.01 ≈ sum(abs2, c) rtol = 1e-6
    @test_throws ArgumentError fourier_transform([0.0, 1.0, 3.0], [1.0, 1.0, 1.0], ω)
end

@testset "Pauli–Fierz Hamiltonian, nucleus–electron–photon (3D)" begin
    # Quadratic model: harmonic nucleus (R, mass 4) and electron (z), dipole
    # d = R − z, one cavity mode q (ω_c = 1, g = 0.1). Exact levels: Σₖ Ωₖ(nₖ + ½),
    # with Ωₖ² the eigenvalues of the mass-weighted Hessian.
    kR, kz, M_R, ω_c, g = 1.0, 1.0, 4.0, 1.0, 0.1
    V = pauli_fierz_potential(r -> 0.5kR * r[1]^2 + 0.5kz * r[2]^2, r -> r[1] - r[2];
                              frequency = ω_c, coupling = g)
    shift = g * sqrt(2ω_c) / ω_c^2
    u, e_q = [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]
    hessian = Diagonal([kR, kz, ω_c^2]) + ω_c^2 * shift^2 * u * u' + ω_c^2 * shift * (u * e_q' + e_q * u')
    weights = Diagonal(1 ./ sqrt.([M_R, 1.0, 1.0]))
    Ω = sort(sqrt.(eigvals(Symmetric(weights * hessian * weights))))
    exact = [sum(Ω) / 2, sum(Ω) / 2 + Ω[1], sum(Ω) / 2 + Ω[2]]

    model = CartesianDiscreteModel((-6.0, 6.0, -6.0, 6.0, -6.0, 6.0), (16, 16, 16))
    sol = solve_tise(model, V; masses = (M_R, 1.0, 1.0), nev = 3)
    @test sol.energies ≈ exact rtol = 5e-3
    @test all(sol.energies .> exact)                       # upper bounds

    # without the dipole self-energy the potential differs by g²d²/ω
    V_no_dse = pauli_fierz_potential(r -> 0.0, r -> r[1] - r[2]; frequency = ω_c, coupling = g,
                                     dipole_self_energy = false)
    point = VectorValue(1.0, -0.5, 0.3)
    V_dse = pauli_fierz_potential(r -> 0.0, r -> r[1] - r[2]; frequency = ω_c, coupling = g)
    @test V_dse(point) - V_no_dse(point) ≈ g^2 / ω_c * 1.5^2
end

@testset "Plotting extension (Makie)" begin
    λ, μ = [0.5, 1.0, 2.0], [0.0, 1.0]
    energies = rand(3, 3, 2)
    # Makie is loaded with `import` in runtests.jl: `using Makie` would make names
    # exported by both packages (e.g. `density`) ambiguous.
    @test plot_energy_map(λ, μ, energies) isa Makie.Figure
    @test plot_energy_map(λ, μ, energies; gap = (1, 2)) isa Makie.Figure
    @test plot_energy_curves(λ, energies[:, :, 1]) isa Makie.Figure
end
