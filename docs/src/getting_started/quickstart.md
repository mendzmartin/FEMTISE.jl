# Quick start

There are three ways of using `FEMTISE.jl`. They solve the same mathematical problem
and differ only in how much control you need.

| Route | When to use it | Control |
|---|---|---|
| [1. Interactive console](@ref qs-console) | First contact; built-in potentials | Minimal |
| [2. Input file](@ref qs-input) | Reproducible 1D/2D calculations, parameter sweeps | Domain, mesh, potential, outputs |
| [3. Low-level interface](@ref qs-lowlevel) | Custom meshes, element degree, quadrature, masses, observables | Full |

All quantities are in **atomic units** (``\hbar = m_e = e = 1``). The kinetic
coefficient for an electron is therefore ``p = 1/2``.

## [1. Interactive console](@id qs-console)

```julia
using FEMTISE
run_default_eigen_problem(set_type_potential())
```

A menu offers the built-in potentials:

```
   - Unidimensional Quantum Harmonic Oscillator            --> set (1)
   - Unidimensional Symmetric Finite Kronig-Penney         --> set (2)
   - Unidimensional Finite Well Potential                  --> set (3)
   - Bidimensional Isotropic Quantum Harmonic Oscillator   --> set (4)
   - Ad hoc potential from input file                      --> set (5)
```

and then asks, one by one, for the output path, the number of eigenvalues, the
domain, the mesh and the potential parameters (pressing *Enter* accepts the default
shown in parentheses). Results are written as binary files; see
[Output files and reloading results](@ref).

## [2. Input file](@id qs-input)

Create two files in the same folder.

**`potential.jl`** — the potential as a function of the coordinates `x` (a point:
`x[1]`, and `x[2]` in 2D) and a tuple of parameters:

```julia
function my_qho(x, params::Tuple)
    ω, x0 = params
    return 0.5 * ω^2 * (x[1] - x0)^2
end
```

**`input.dat`** — exactly 21 lines in this order (see [Input file reference](@ref)):

```text
full_path_name              = ./results/qho
dom_type                    = s
nev                         = 10
dimension                   = 1D
sigma                       = 0.0
adhoc_file_name             = ./potential
potential_function_name     = my_qho
params_potential_types      = f f
params_potential            = 1.0 0.0
analysis_param              = false
output_format_type          = jld2 eigen
## ONLY FOR 1D EIGENPROBLEMS
L                           = 20.0
Δx                          = 0.05
## ONLY FOR 2D EIGENPROBLEMS
Lx                          =
Ly                          =
nx                          =
ny                          =
different_masses            =
reduced_density             =
```

Run it and reload the results:

```julia
using FEMTISE
mkpath("./results")
run_default_eigen_problem(set_type_potential("./input"))   # no ".dat" extension

id, res = collect_result_data(true, "./input")
real.(res.ϵ)      # eigenvalues
res.ϕ             # eigenfunctions (Gridap finite element functions)
```

To compute the eigenvalues as a function of, e.g., the frequency `ω` (first
parameter) from 0.5 to 2.0 in steps of 0.1, change one line:

```text
analysis_param              = 1 0.5 2.0 0.1
```

The result then contains `res.ϵ_matrix` (`nev × number of values`) and
`res.λvector`.

## [3. Low-level interface](@id qs-lowlevel)

```julia
using FEMTISE, Gridap

model = CartesianDiscreteModel((-10.0, 10.0), (200,))        # mesh
Ω, dΩ, Γ, dΓ = measures(model, 6, "boundary")                # quadrature measures
V, U = fe_spaces(model, ReferenceFE(lagrangian, Float64, 2);
                 BC_type = "Dirichlet", BC_data = (0.0im, "boundary"))

p(x) = 0.5; q(x) = 0.5 * x[1]^2; r(x) = 1.0
ε, ϕ = eigen_values_and_eigen_vectors(p, q, r, dΩ, U, V;
           params = (10, 1e-10, 500, :none, 0.0))   # (nev, tol, maxiter, mode, σ)

norm_l2(ϕ[1], dΩ)          # observables act directly on the FE eigenfunctions
```

Continue with the [Low-level interface](@ref) chapter, or with a complete tutorial
such as [the 1D harmonic oscillator](../examples/isotropic_1d_harmonic_oscillator_potential.md).
