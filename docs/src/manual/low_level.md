# Low-level interface

The low-level interface exposes each step of the calculation. Use it when you
need a different element degree, quadrature, mesh, mass, boundary geometry, or when
FEMTISE is part of a larger program (e.g. potential energy curves inside a loop).

## Step by step

```julia
using FEMTISE, Gridap

# 1. Mesh: 400 cells on (-20, 20). In 2D: CartesianDiscreteModel((x1,x2,y1,y2), (nx,ny)).
model = CartesianDiscreteModel((-20.0, 20.0), (400,))

# 2. Domain, boundary and quadrature measures. The second argument is the quadrature
#    degree k_q; "boundary" is the tag of the whole boundary of a Cartesian model.
Ω, dΩ, Γ, dΓ = measures(model, 6, "boundary")

# 3. Finite element spaces: Lagrange elements of degree 2, ψ = 0 on the boundary.
reff = ReferenceFE(lagrangian, Float64, 2)
V, U = fe_spaces(model, reff; BC_type = "Dirichlet", BC_data = (0.0im, "boundary"))

# 4. Coefficients of -∇·(p∇ψ) + qψ = ε r ψ
p(x) = 0.5                         # ħ²/2m for an electron
q(x) = -1.0 / sqrt(1.0 + x[1]^2)   # soft-Coulomb potential
r(x) = 1.0

# 5. Forms, matrices and eigenpairs
a, b = bilineal_forms(p, q, r, dΩ)
problem = eigen_problem(a, b, U, V; nev = 8, sigma = -1.0, tol = 1e-10, maxiter = 500)
ε, ϕ = FEMTISE.solve(problem)

# Step 5 in a single call:
ε, ϕ = eigen_values_and_eigen_vectors(p, q, r, dΩ, U, V;
                                      params = (8, 1e-10, 500, :none, -1.0))
```

`params` is the tuple `(nev, tol, maxiter, explicittransform, sigma)` passed to
`Arpack.eigs`. Note that `fe_spaces` returns `(V, U)` (test space first).

!!! note "`solve` and Gridap"
    Gridap also exports a function called `solve`. When both packages are loaded,
    write `FEMTISE.solve(problem)`.

## Choosing the discretization

| Choice | Recommendation |
|---|---|
| Element degree ``\ell`` | 2 is a good default; 3 for high accuracy. The cost per unknown is almost independent of ``\ell``. |
| Quadrature degree ``k_q`` | ``2\ell + s`` for a polynomial potential of degree ``s`` (exact integrals, guaranteed upper bounds). For non-polynomial potentials, ``2\ell+2`` or more. |
| Element size ``h`` | Resolve the shortest wavelength of interest. Halving ``h`` reduces the eigenvalue error by ``2^{2\ell}``. |
| Domain size | Large enough that the states of interest decay to negligible values at the walls (check by enlarging it). |
| Shift ``\sigma`` | Below the lowest level wanted, or at the centre of the energy window of interest. |
| Data type | `ComplexF64` spaces (default) are general. For real Hamiltonians, `fe_spaces(...; TypeData = Float64, BC_data = (0.0, "boundary"))` halves the memory and lets ARPACK use its symmetric driver. |

## Position-dependent mass

Pass a non-constant `p`. For a heterostructure with an interface at ``x=\pm w/2``:

```julia
mstar(z) = abs(z) < w/2 ? 0.067 : 0.092          # effective masses (units of mₑ)
p(x) = 0.5 / mstar(x[1])
q(x) = abs(x[1]) < w/2 ? 0.0 : Vb
```

Make sure that **the interfaces coincide with mesh nodes**. Otherwise convergence
drops to first order (see [Accuracy and best practices](@ref)).

## Two particles with different masses

On a 2D mesh whose coordinates are the positions of two 1D particles, the kinetic
term is anisotropic, ``p_1\partial_1^2 + p_2\partial_2^2``:

```julia
p₁(x) = 0.5          # particle 1: electron mass
p₂(x) = 0.5 / 3.0    # particle 2: mass 3
ε, ϕ = eigen_values_and_eigen_vectors(p₁, p₂, q, r, dΩ, U, V; params = (10, 1e-10, 500, :none, 0.0))
```

## Other meshes

[`make_model`](@ref) builds Gmsh meshes: `"simple_line"` (1D segment with element
size Δx, tags `"left_point"`/`"right_point"`), `"simple_rectangle_v1"` (2D
triangles or quadrilaterals) and `"Cartesian2D"`. Any Gridap discrete model can be
used, including meshes read from Gmsh files with `GmshDiscreteModel("mesh.msh")`.
Pass the corresponding Dirichlet tags to `measures` and `fe_spaces`
([`make_boundary_conditions`](@ref) returns them for the built-in mesh types).

## Repeated solves (potential energy curves, sweeps)

Build the model, measures and spaces once, and call
`eigen_values_and_eigen_vectors` for each potential:

```julia
for R in 0.5:0.1:10.0
    qR(x) = potential(x[1], R)      # potential(z, R): your function
    εR, = eigen_values_and_eigen_vectors(p, qR, r, dΩ, U, V; params = (20, 1e-10, 500, :none, σ))
    # store εR ...
end
```

When many eigenpairs are needed (hundreds), split the spectrum into several windows
with different shifts. The cost grows faster than linearly with `nev`.
