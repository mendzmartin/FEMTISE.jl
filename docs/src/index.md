# FEMTISE.jl

*Time-independent Schrödinger equation by the finite element method, in Julia.*

`FEMTISE.jl` computes eigenvalues and eigenfunctions of Hamiltonians of the form

```math
\hat H = -\nabla\cdot\big(\mathbf{p}(\mathbf{x})\,\nabla\big) + V(\mathbf{x}),
\qquad \mathbf{p} = \tfrac{\hbar^2}{2}\,\boldsymbol{\mathcal{M}}^{-1},
```

with conforming Lagrange finite elements built on
[Gridap.jl](https://github.com/gridap/Gridap.jl) and the shift-and-invert Arnoldi
eigensolver of [ARPACK](https://github.com/JuliaLinearAlgebra/Arpack.jl). Because
the method is variational, every computed eigenvalue is an upper bound to the exact
one (with exact quadrature), and the error decreases as ``h^{2\ell}`` with the
element size ``h`` and the polynomial degree ``\ell``.

## Features

- **Arbitrary potentials** ``V(\mathbf{x})`` written as ordinary Julia functions,
  including discontinuous ones (heterostructures, square wells, lattices).
- **Position-dependent and anisotropic masses**: BenDaniel–Duke operator
  ``-\partial_x[\hbar^2/2m^*(x)]\partial_x`` and two particles of different masses.
- **One- and two-dimensional problems** from a plain-text input file, or any mesh,
  element degree and quadrature supported by Gridap through the low-level interface.
- **Parameter sweeps**: eigenvalues as functions of a potential parameter, reusing
  the mesh for every value.
- **Post-processing** on the finite element eigenfunctions: norms, densities,
  reduced densities, moments, Shannon and Rényi entropies, mutual information.
- **Outputs** as compact binary arrays or JLD2 files that can be reloaded to
  compute new observables without solving again.

## Where to start

| You are… | Read |
|---|---|
| New to the package | [Installation](@ref) → [Quick start](@ref) → one of the tutorials |
| Running calculations from input files | [Input file reference](@ref) and [Output files and reloading results](@ref) |
| Writing your own scripts | [Low-level interface](@ref) and [Post-processing](manual/postprocessing.md) |
| Checking the reliability of a result | [Theory](@ref) and [Accuracy and best practices](@ref) |
| Upgrading between versions | [Changelog](@ref) |
| Looking up a function | [API reference](@ref) |

## A minimal example

The ten lowest levels of the one-dimensional harmonic oscillator ``V=x^2/2`` (atomic
units) with cubic elements:

```julia
using FEMTISE, Gridap

model = CartesianDiscreteModel((-10.0, 10.0), (200,))     # 200 cells on (-10, 10)
Ω, dΩ, Γ, dΓ = measures(model, 8, "boundary")             # quadrature of degree 8
reff = ReferenceFE(lagrangian, Float64, 3)                # cubic Lagrange elements
V, U = fe_spaces(model, reff; BC_type = "Dirichlet",
                 BC_data = (0.0im, "boundary"))           # ψ = 0 on the walls

p(x) = 0.5            # ħ²/2m
q(x) = 0.5 * x[1]^2   # potential
r(x) = 1.0            # weight
ε, ϕ = eigen_values_and_eigen_vectors(p, q, r, dΩ, U, V;
                                      params = (10, 1e-10, 500, :none, 0.0))
real.(ε)   # ≈ 0.5, 1.5, 2.5, ...
```

## Versions of this documentation

The documentation is versioned. Use the version selector at the bottom of the left
menu: **stable** is the latest release, **dev** follows the `main` branch, and each
release (e.g. **v0.1.0**) keeps the documentation that matches its code.

## Citing

If you use `FEMTISE.jl`, please cite it as described in [Citing and contact](@ref).

The package was developed during the PhD studies of M. Mendez at the
[Facultad de Matemática, Astronomía, Física y Computación (FAMAF)](https://www.famaf.unc.edu.ar/)
and the Instituto de Física Enrique Gaviola (IFEG, CONICET–UNC), Córdoba, Argentina.
