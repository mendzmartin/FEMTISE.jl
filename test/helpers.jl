# Shared helpers for the test suite.
#
# Every helper uses only the public low-level API of FEMTISE, so the tests also
# document how that API is meant to be used. Units: Hartree atomic units
# (ħ = mₑ = 1), so the kinetic coefficient of an electron is p = 1/2.

"""
    solve_on_mesh(model; p, q, order=2, quadrature_degree=2order+2, nev=6,
                  sigma=0.0, p2=nothing)

Solve `-∇·(p∇ψ) + qψ = εψ` with homogeneous Dirichlet conditions on the whole
boundary of `model`, using Lagrange elements of degree `order` and a Gauss
quadrature exact for polynomials of degree `quadrature_degree`. If `p2` is
given, the kinetic term is anisotropic, `diag(p, p2)` (two particles of
different masses on a 2D mesh).

Returns `(ε, ϕ, Ω, dΩ, U)` with the eigenvalues as real numbers sorted in
ascending order and the eigenfunctions reordered accordingly.
"""
function solve_on_mesh(model; p, q, order::Int = 2, quadrature_degree::Int = 2order + 2,
                       nev::Int = 6, sigma::Real = 0.0, p2 = nothing)
    Ω, dΩ, _, _ = measures(model, quadrature_degree, "boundary")
    reference_element = ReferenceFE(lagrangian, Float64, order)
    V, U = fe_spaces(model, reference_element; BC_type = "Dirichlet",
                     BC_data = (zero(ComplexF64), "boundary"))
    weight(x) = 1.0
    coefficients = isnothing(p2) ? (p, q, weight) : (p, p2, q, weight)
    solver_params = (nev, 1e-12, 1000, :none, Float64(sigma))   # (nev, tol, maxiter, mode, σ)
    ε, ϕ = eigen_values_and_eigen_vectors(coefficients..., dΩ, U, V; params = solver_params)
    order_of_levels = sortperm(real.(ε))
    return real.(ε[order_of_levels]), ϕ[order_of_levels], Ω, dΩ, U
end

"Kinetic coefficient p = ħ²/2m of an electron in atomic units."
electron_kinetic(x) = 0.5

"Harmonic potential V = x²/2 (1D) or V = (x₁² + x₂²)/2 (2D) with ω = 1."
harmonic_1d(x) = 0.5 * x[1]^2
harmonic_2d(x) = 0.5 * (x[1]^2 + x[2]^2)

"""
    bisection(f, a, b; tol=1e-15)

Root of `f` in `[a, b]` by bisection; `f(a)` and `f(b)` must have opposite signs.
Used to obtain reference eigenvalues from transcendental matching equations.
"""
function bisection(f, a, b; tol = 1e-15)
    fa = f(a)
    @assert fa * f(b) < 0 "f must change sign in [a, b]"
    while b - a > tol * max(1.0, abs(a))
        m = (a + b) / 2
        fm = f(m)
        if fa * fm <= 0
            b = m
        else
            a, fa = m, fm
        end
    end
    return (a + b) / 2
end
