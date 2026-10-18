# Architecture

`FEMTISE.jl` is organized in two layers on top of
[Gridap.jl](https://github.com/gridap/Gridap.jl), which provides meshes, finite
element spaces, quadrature and assembly.

```
 input file (.dat) + potential (.jl)        interactive console
                 \                             /
                  ▼                           ▼
     ┌──────────────────────────────────────────────────┐
     │ DRIVER   run_default_eigen_problem               │   fixed, validated configuration
     │          set_type_potential, input_data          │   (1D/2D, P2 elements, Dirichlet)
     └──────────────────────────────────────────────────┘
                              │
     ┌──────────────────────────────────────────────────┐
     │ LOW-LEVEL INTERFACE (wraps Gridap)               │   full control
     │  make_model ─► measures ─► fe_spaces             │
     │      ─► bilineal_forms ─► eigen_problem ─► solve │
     └──────────────────────────────────────────────────┘
                              │
     ┌──────────────────────────────────────────────────┐
     │ POST-PROCESSING  norms, densities, moments,      │
     │                  entropies, mutual information   │
     │ OUTPUT           write_bin / JLD2 → collect_result_data │
     └──────────────────────────────────────────────────┘
```

## From mathematics to code

Each object of the [Theory](@ref) has one function that builds it:

| Mathematical object | Function | Returns |
|---|---|---|
| Mesh ``\mathcal T_h`` | [`make_model`](@ref) or `Gridap.CartesianDiscreteModel` | discrete model |
| Boundary data (Dirichlet tags/values) | [`make_boundary_conditions`](@ref) | values, tags |
| Domain ``\Omega``, boundary ``\Gamma`` and quadrature measures ``d\Omega``, ``d\Gamma`` | [`measures`](@ref) | `Ω, dΩ, Γ, dΓ` |
| Test and trial spaces ``V_h``, ``U_h`` | [`fe_spaces`](@ref) | `V, U` |
| Forms ``a(u,v)``, ``b(u,v)`` | [`bilineal_forms`](@ref) | two Julia functions |
| Matrices ``\mathbb A``, ``\mathbb B`` and solver settings | [`eigen_problem`](@ref) | `EigenProblem` |
| Eigenpairs ``(\varepsilon_n, u_{n,h})`` | [`FEMTISE.solve`](@ref) | eigenvalues, `Vector{CellField}` |
| All of the above in one call | [`eigen_values_and_eigen_vectors`](@ref) | eigenvalues, eigenfunctions |

The eigenfunctions are Gridap `CellField`s: they can be evaluated at any point
(`ϕ[1](Point(0.3))`), differentiated (`∇(ϕ[1])`) and integrated
(`sum(∫(ϕ[1]' * ϕ[1]) * dΩ)`).

## Source files

| File (`src/functions/`) | Content |
|---|---|
| `mesh_generator_function.jl` | `make_model`: Gmsh 1D lines and 2D rectangles, native Cartesian 2D meshes |
| `boundary_conditions_function.jl` | `make_boundary_conditions`: Dirichlet tags and values per mesh type |
| `miscellaneous_functions.jl` | `measures`, `fe_spaces`, `bilineal_forms`, `space_coord` |
| `eigen_problem_definition_function.jl` | `EigenOperator`, `EigenProblem`, `eigen_problem`, `solve` |
| `eigen_problem_solve_function.jl` | `eigen_values_and_eigen_vectors` |
| `useful_functions_for_FEM_objects.jl` | post-processing (norms, densities, entropies, moments) |
| `default_running_structures.jl` | `Params1D`, `Params2D`, `AnalysisParam`, `InputData*` |
| `default_running_input_data.jl` | `input_data`: parser of the input file |
| `default_running_miscellaneous_functions.jl` | driver solvers, built-in potentials, console prompts |
| `default_eigen_problem.jl` | `run_default_eigen_problem`, `set_type_potential`, writing of outputs |
| `binary_file_io.jl` | `write_bin`, `read_bin` |
| `post_proccesing_data_structures.jl`, `post_proccesing_data_function.jl` | result containers and `collect_result_data` |

## Configuration of the driver

The driver fixes the following choices. To change any of them, use the
low-level interface.

| Setting | Value |
|---|---|
| Elements | Lagrange, degree ``\ell = 2``, complex-valued spaces (`ComplexF64`) |
| Boundary condition | homogeneous Dirichlet on the whole boundary |
| 1D mesh | Gmsh segment with element size `Δx` on ``(-L/2, L/2)`` or ``(0, L)`` |
| 2D mesh | Cartesian ``n_x\times n_y`` quadrilaterals on a rectangle |
| Quadrature degree | ``k_q = 3`` |
| ARPACK | tolerance ``10^{-8}``, at most 500 restarts (1500 for 1D sweeps), `which = :LM`, shift `sigma` from the input |
| Units | atomic units, ``p = 1/2`` (electron mass) |

## Design notes

- **Coefficient functions are plain Julia functions.** They are compiled and inlined
  into the assembly loops, so a user-defined potential is as fast as a built-in one.
- **Potentials loaded by the driver** are `include`d inside the `FEMTISE` module.
  Functions of `SpecialFunctions` (e.g. `erfc`) are available there.
- **Parameter sweeps** build the mesh once; for each value only the matrices and
  the factorization are recomputed.
