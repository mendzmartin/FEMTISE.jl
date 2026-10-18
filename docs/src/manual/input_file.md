# Input file reference

The input-file driver solves 1D and 2D problems described by two files:

- an **input file** `name.dat` with the settings, and
- a **potential file** `name.jl` with a Julia function for the potential.

```julia
using FEMTISE
run_default_eigen_problem(set_type_potential("path/to/name"))   # without ".dat"
```

!!! tip "Calling the driver from inside a function"
    `set_type_potential` defines the potential by `include`-ing the potential file.
    At the top level of a script or the REPL this is transparent. If the driver is
    called from inside a function (or a `do` block), that function was compiled
    before the potential existed and Julia reports *"method too new to be called
    from this world context"*. Call the driver in the latest world:
    ```julia
    function my_workflow(path)
        simulation = set_type_potential(path)
        Base.invokelatest(run_default_eigen_problem, simulation)
    end
    ```

## Format rules

!!! warning "The layout is fixed"
    The parser reads the file **by line number** and takes the **third
    whitespace-separated token** of each line as the value. Therefore:
    - the file must have exactly the 21 lines shown below, in this order,
      including the two `## ONLY FOR ...` lines;
    - the separator must be ` = ` with spaces on both sides;
    - fields that do not apply are left empty, but the line must be present;
    - lists of values (`params_potential`, `analysis_param`, `output_format_type`)
      are separated by spaces.

```text
full_path_name              = ./results/my_run
dom_type                    = s
nev                         = 10
dimension                   = 2D
sigma                       = -10.0
adhoc_file_name             = ./my_potential
potential_function_name     = my_potential
params_potential_types      = f f f
params_potential            = 1.0 0.0 0.0
analysis_param              = false
output_format_type          = jld2 eigen
## ONLY FOR 1D EIGENPROBLEMS
L                           =
Δx                          =
## ONLY FOR 2D EIGENPROBLEMS
Lx                          = 12.0
Ly                          = 12.0
nx                          = 72
ny                          = 72
different_masses            = false
reduced_density             = false
```

## Fields

| Line | Field | Type | Meaning |
|---|---|---|---|
| 1 | `full_path_name` | path | Prefix of all output files, e.g. `./results/run` → `./results/run_eigen_data.jld2`. The folder must exist. |
| 2 | `dom_type` | `s`, `ns` (+ `sns`, `nss` in 2D) | `s`: symmetric domain ``(-L/2, L/2)``; `ns`: ``(0, L)``. In 2D, `sns` is symmetric in ``x`` and ``(0,L_y)`` in ``y``; `nss` the opposite. |
| 3 | `nev` | integer | Number of eigenpairs. |
| 4 | `dimension` | `1D`, `2D` or `3D` | Dimension of the problem (`3D` since v0.2.0, see the optional settings). |
| 5 | `sigma` | real | Shift ``\sigma`` of the eigensolver: the `nev` eigenvalues **closest to `sigma`** are returned. Put it below the ground state to obtain the lowest levels. |
| 6 | `adhoc_file_name` | path | Potential file **without** the `.jl` extension. |
| 7 | `potential_function_name` | name | Name of the potential function defined in that file. |
| 8 | `params_potential_types` | list of `f`, `i`, `c` | Type of each parameter: `f` = `Float64`, `i` = `Int`, `c` = `ComplexF64`. |
| 9 | `params_potential` | list | Values of the parameters, passed to the potential as a `Tuple`. |
| 10 | `analysis_param` | `false` or 4 values | Parameter sweep, see below. |
| 11 | `output_format_type` | 2 words | `bin eigen`, `jld2 eigen` or `jld2 all`, see [Output files and reloading results](@ref). Sweeps accept `bin eigen` and `jld2 eigen`. |
| 13 | `L` | real | 1D: length of the domain. |
| 14 | `Δx` | real | 1D: element size. `L/Δx` should be an integer. |
| 16, 17 | `Lx`, `Ly` | real | 2D: side lengths. |
| 18, 19 | `nx`, `ny` | integer | 2D: number of cells in each direction (element size `Lx/nx`, `Ly/ny`). |
| 20 | `different_masses` | `false` or real | `false`: both coordinates have the electron mass. A real number ``\mu`` sets ``p_2 = \mu/2``, i.e. ``\mu = M_1/M_2`` with ``M_1 = 1``. |
| 21 | `reduced_density` | `true`/`false` | 2D: also compute the reduced densities ``\rho^{(1)}(x)``, ``\rho^{(2)}(y)``. |

## Optional settings (since v0.2.0)

Lines after line 21 are optional and are read **by name** (`key = value`), in any
order. Lines starting with `#` are comments. A misspelled key raises an error instead
of being ignored.

```text
## OPTIONAL SETTINGS
masses                      = 1836.15
element_order               = 3
quadrature_degree           = 8
potential_data_file         = ./potential.dat
interpolation               = cubic
analysis_param_2            = 2 0.0 1.0 0.1
Lz                          = 12.0
nz                          = 16
```

| Key | Values | Meaning |
|---|---|---|
| `masses` | one value, or one per coordinate | masses in units of the electron mass (default 1) |
| `element_order` | integer | degree of the Lagrange elements (default 2) |
| `quadrature_degree` | integer | degree of exactness of the quadrature (default 3); use ``2\ell+2`` or more for guaranteed upper bounds |
| `potential_data_file` | path | tabulated potential (see [Masses, data files and sweeps](new_inputs.md)); leave lines 6–9 empty |
| `interpolation` | `linear`, `cubic` | interpolation of the tabulated potential |
| `analysis_param_2` | `j λᵢ λ_f Δλ` | second swept parameter (needs `analysis_param`) |
| `Lz`, `nz` | real, integer | third side and number of cells when `dimension = 3D` |

Input files with optional settings, and all 3D problems, are solved with the general
solver [`solve_tise`](@ref). They require `output_format_type = jld2 eigen` and return
the eigenvalues in ascending order. Input files without optional settings are solved
exactly as in v0.1.

**Three-dimensional problems.** Set `dimension = 3D`, use the 2D lines for `Lx`, `Ly`,
`nx`, `ny` (lines 20–21 are ignored) and add `Lz` and `nz`. The potential receives
`x[1]`, `x[2]`, `x[3]`.

## The potential file

The potential receives the point `x` (`x[1]`, and `x[2]` in 2D) and the tuple of
parameters, in the order given in the input file:

```julia
# my_potential.jl  —  params = (ω, x0, y0)
function my_potential(x, params::Tuple)
    ω, x0, y0 = params
    return 0.5 * ω^2 * ((x[1] - x0)^2 + (x[2] - y0)^2)
end
```

Guidelines:
- Use only arithmetic and functions available in `FEMTISE` (Base, `SpecialFunctions`);
  the file is loaded inside the `FEMTISE` module.
- Keep the function type-stable and allocation-free: it is evaluated at every
  quadrature point of every cell.
- In 2D, the two coordinates may be two particles (e.g. two electrons on a line).
  Then the potential contains one-body terms ``V(x_1)+V(x_2)`` and the interaction
  ``W(x_1-x_2)``.

## Parameter sweeps

`analysis_param = j λᵢ λ_f Δλ` computes the `nev` eigenvalues for a sweep of the
`j`-th potential parameter. The mesh is built only once. The computed values are

```math
\lambda_k = \lambda_i + k\,\Delta\lambda,\qquad k = 1, 2, \dots, \operatorname{round}\!\big(|\lambda_f-\lambda_i|/\Delta\lambda\big).
```

!!! note "First value of the sweep"
    In version 0.1.0 the first computed value is ``\lambda_i+\Delta\lambda``, not
    ``\lambda_i``. To include ``\lambda_i``, start the sweep one step earlier.

Example: `analysis_param = 1 0.5 2.0 0.1` sweeps the first parameter over
``0.6, 0.7, \dots, 2.0``. Only eigenvalues are stored (`ϵ_matrix`, of size
`nev × number of values`, and `λvector`).

## Interactive console

`run_default_eigen_problem(set_type_potential())` asks for the same information
interactively. Option 5 asks for the path of an input file. Options 1–4 run
built-in potentials:

| Option | Potential | Parameters |
|---|---|---|
| 1 | 1D harmonic oscillator ``\tfrac12\omega^2(x-x_1)^2`` | ``\omega``, ``x_1`` |
| 2 | 1D symmetric finite Kronig–Penney chain | depth ``V_0`` (negative), period, well width, number of wells (odd) |
| 3 | 1D finite square well | depth ``V_0`` (negative), width |
| 4 | 2D isotropic harmonic oscillator | ``\omega``, ``x_1``, ``y_1`` |

Interactive runs also write `name_eigen_problem_attributes.dat`, a summary of
the settings that `collect_result_data(false, name)` uses to reload the results.
