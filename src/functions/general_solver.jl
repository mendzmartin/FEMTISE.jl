# =============================================================================
# General solver: one interface for 1, 2 and 3 coordinates, arbitrary masses.
#
#   H = - Σₖ (1 / 2Mₖ) ∂²/∂xₖ² + V(x),        x = (x₁, ..., x_d),  d = 1, 2, 3
#
# in atomic units (ħ = 1, electron mass = 1). Each coordinate xₖ may belong to a
# different particle (e.g. nucleus, electron, photon mode), so each one has its
# own mass Mₖ. The weak form uses the anisotropic kinetic coefficients
# pₖ = 1/(2Mₖ):
#
#   a(u, v) = ∫ [ Σₖ pₖ ∂ₖu ∂ₖv + V u v ] dΩ ,     b(u, v) = ∫ u v dΩ .
# =============================================================================

"""
    kinetic_coefficients(masses, dimension)

Kinetic coefficients ``p_k = 1/(2M_k)`` (atomic units) for each of the `dimension`
coordinates. `masses` is a number (same mass for every coordinate) or a collection
with one mass per coordinate, in units of the electron mass.

```julia
kinetic_coefficients(1.0, 1)                 # (0.5,)          electron
kinetic_coefficients((1836.15, 1.0, 1.0), 3) # proton, electron, photon mode
```
"""
function kinetic_coefficients(masses, dimension::Integer)
    mass_vector = masses isa Number ? fill(Float64(masses), dimension) : Float64.(collect(masses))
    length(mass_vector) == dimension || throw(ArgumentError(
        "got $(length(mass_vector)) masses for $dimension coordinates; give one mass per coordinate"))
    all(>(0), mass_vector) || throw(ArgumentError("masses must be positive, got $mass_vector"))
    return Tuple(0.5 ./ mass_vector)
end

"""
    hamiltonian_forms(masses, potential, dΩ, dimension; weight = x -> 1.0)

Forms ``a(u,v)`` (Hamiltonian) and ``b(u,v)`` (overlap) for the kinetic coefficients
given by `masses` (see [`kinetic_coefficients`](@ref)) and the potential function
`potential(x)`. When all masses are equal the isotropic form of
[`bilineal_forms`](@ref) is used.
"""
function hamiltonian_forms(masses, potential, dΩ, dimension::Integer; weight = x -> 1.0)
    p = kinetic_coefficients(masses, dimension)
    if all(==(p[1]), p)
        return bilineal_forms(x -> p[1], potential, weight, dΩ)
    end
    unit_vectors = ntuple(k -> VectorValue(ntuple(i -> i == k ? 1.0 : 0.0, dimension)), dimension)
    kinetic(u, v) = sum(p[k] * ((∇(v) ⋅ unit_vectors[k]) * (∇(u) ⋅ unit_vectors[k])) for k in 1:dimension)
    a(u, v) = ∫(kinetic(u, v) + potential * v * u) * dΩ
    b(u, v) = ∫(weight * u * v) * dΩ
    return a, b
end

"""
    TISESolution

Result of [`solve_tise`](@ref).

# Fields
- `energies::Vector{Float64}`: eigenvalues in ascending order (a.u.).
- `states::Vector{CellField}`: eigenfunctions in the same order (finite element
  functions, normalized: ``\\int|\\psi|^2 d\\Omega = 1``).
- `Ω`, `dΩ`: triangulation and integration measure.
- `trial`, `test`: finite element spaces.
- `model`: the mesh.
"""
struct TISESolution{Ω_T,dΩ_T,U_T,V_T,M_T}
    energies::Vector{Float64}
    states::Vector{CellField}
    Ω::Ω_T
    dΩ::dΩ_T
    trial::U_T
    test::V_T
    model::M_T
end

"""
    discretization(model; order = 2, quadrature_degree = 2order + 2,
                   boundary_tags = "boundary")

Triangulation, measure and finite element spaces (Lagrange elements of degree
`order`, homogeneous Dirichlet conditions on `boundary_tags`) of a mesh. They do
not depend on the potential, so they can be reused for many potentials.

Returns `(Ω, dΩ, trial, test)`.
"""
function discretization(model; order::Integer = 2, quadrature_degree::Integer = 2order + 2,
                        boundary_tags = "boundary")
    Ω, dΩ, _, _ = measures(model, quadrature_degree, boundary_tags)
    test, trial = fe_spaces(model, ReferenceFE(lagrangian, Float64, order);
                            BC_type = "Dirichlet", BC_data = (zero(ComplexF64), boundary_tags))
    return Ω, dΩ, trial, test
end

"""
    nearest_sorted(ε, ϕ, nev, sigma)

The `nev` eigenpairs closest to `sigma`, sorted by increasing energy (real parts of
the eigenvalues).
"""
function nearest_sorted(ε, ϕ, nev::Integer, sigma::Real)
    closest = sortperm(abs.(real.(ε) .- sigma))[1:min(nev, length(ε))]
    order = closest[sortperm(real.(ε[closest]))]
    return real.(ε[order]), ϕ[order]
end

"""
    solve_tise(model, potential; masses = 1.0, order = 2,
               quadrature_degree = 2order + 2, nev = 10, sigma = 0.0,
               tol = 1e-10, maxiter = 1000, guard_states = 4,
               boundary_tags = "boundary")

Solve the time-independent Schrödinger equation

```math
\\Big[-\\sum_k \\frac{1}{2M_k}\\frac{\\partial^2}{\\partial x_k^2} + V(\\mathbf{x})\\Big]\\psi
= \\varepsilon\\,\\psi
```

on `model` (a Gridap mesh in 1, 2 or 3 dimensions) with ``\\psi = 0`` on the
boundary, and return the `nev` eigenpairs closest to `sigma` as a
[`TISESolution`](@ref), sorted by energy.

# Arguments
- `potential`: function of the mesh point `x` (`x[1]`, `x[2]`, `x[3]`), e.g. a
  Julia function or a [`tabulated_potential`](@ref).
- `masses`: mass of every coordinate, or one mass per coordinate, in units of the
  electron mass (default: electrons).
- `order`: degree of the Lagrange elements.
- `quadrature_degree`: degree of exactness of the Gauss quadrature. The default
  integrates exactly potentials up to degree 2 (eigenvalues are upper bounds).
- `nev`, `sigma`, `tol`, `maxiter`: settings of the shift-and-invert eigensolver.
- `guard_states`: number of additional eigenpairs computed and discarded. Krylov
  eigensolvers can miss members of a degenerate level at the edge of the requested
  window; the extra states make the `nev` returned eigenpairs reliable.

# Example
```julia
using FEMTISE, Gridap
model = CartesianDiscreteModel((-10.0, 10.0), (200,))
sol = solve_tise(model, x -> 0.5 * x[1]^2; nev = 5)
sol.energies          # ≈ [0.5, 1.5, 2.5, 3.5, 4.5]
```
"""
function solve_tise(model, potential; masses = 1.0, order::Integer = 2,
                    quadrature_degree::Integer = 2order + 2, nev::Integer = 10,
                    sigma::Real = 0.0, tol::Real = 1e-10, maxiter::Integer = 1000,
                    guard_states::Integer = 4, boundary_tags = "boundary")
    Ω, dΩ, trial, test = discretization(model; order, quadrature_degree, boundary_tags)
    a, b = hamiltonian_forms(masses, potential, dΩ, num_cell_dims(model))
    problem = eigen_problem(a, b, trial, test; nev = nev + guard_states, sigma = Float64(sigma),
                            tol = Float64(tol), maxiter = maxiter)
    ε, ϕ = nearest_sorted(solve(problem)..., nev, sigma)
    return TISESolution(ε, ϕ, Ω, dΩ, trial, test, model)
end

"""
    sweep_eigenvalues(model, potential, values1[, values2]; masses = 1.0, order = 2,
                      quadrature_degree = 2order + 2, nev = 10, sigma = 0.0,
                      tol = 1e-10, maxiter = 1000, guard_states = 4,
                      boundary_tags = "boundary")

Eigenvalues as functions of one or two potential parameters. `potential(x, λ)`
(one parameter) or `potential(x, λ, μ)` (two parameters) is evaluated for every
value in `values1` (and `values2`). The mesh, measures and finite element spaces
are built once and reused.

Returns a matrix `energies[n, i]` (one parameter) or an array `energies[n, i, j]`
(two parameters) with the `nev` eigenvalues closest to `sigma`, in ascending order
at every point.

# Example
```julia
model = CartesianDiscreteModel((-10.0, 10.0), (200,))
V(x, ω) = 0.5 * ω^2 * x[1]^2
E = sweep_eigenvalues(model, V, 0.5:0.1:2.0; nev = 3)    # E[n, i] ≈ ωᵢ (n - 1/2)
```
"""
function sweep_eigenvalues(model, potential, values1; kwargs...)
    energies = sweep_eigenvalues(model, (x, λ, _) -> potential(x, λ), values1, (nothing,); kwargs...)
    return energies[:, :, 1]
end

function sweep_eigenvalues(model, potential, values1, values2; masses = 1.0, order::Integer = 2,
                           quadrature_degree::Integer = 2order + 2, nev::Integer = 10,
                           sigma::Real = 0.0, tol::Real = 1e-10, maxiter::Integer = 1000,
                           guard_states::Integer = 4, boundary_tags = "boundary")
    Ω, dΩ, trial, test = discretization(model; order, quadrature_degree, boundary_tags)
    dimension = num_cell_dims(model)
    energies = zeros(nev, length(values1), length(values2))
    for (j, μ) in enumerate(values2), (i, λ) in enumerate(values1)
        a, b = hamiltonian_forms(masses, x -> potential(x, λ, μ), dΩ, dimension)
        problem = eigen_problem(a, b, trial, test; nev = nev + guard_states, sigma = Float64(sigma),
                                tol = Float64(tol), maxiter = maxiter)
        ε, ϕ = solve(problem)
        energies[:, i, j] = nearest_sorted(ε, ϕ, nev, sigma)[1]
    end
    return energies
end
