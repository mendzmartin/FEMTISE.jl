# Output files and reloading results

The driver writes its results to files whose names start with `full_path_name`.
The format is selected with `output_format_type` in the input file.

## Formats

| `output_format_type` | Files | Contents | Use when |
|---|---|---|---|
| `bin eigen` | `_eigen_values.bin`, `_eigen_vectors.bin`, `_coordinates.bin` (+ `_reduced_density_DOF1.bin`, `_reduced_density_DOF2.bin` in 2D if requested) | Plain arrays: eigenvalues, eigenfunctions sampled on a uniform grid, grid coordinates | Analysis in any language; smallest files |
| `jld2 eigen` | `_eigen_data.jld2` | Eigenvalues, eigenfunctions as finite element functions, grid | Further analysis in Julia with Gridap |
| `jld2 all` | `_eigen_data.jld2` | As above plus the discrete model (mesh) | Computing new observables later without solving again |

For parameter sweeps (`analysis_param`):

| `output_format_type` | Files | Contents |
|---|---|---|
| `bin eigen` | `_eigen_values.bin`, `_param_values.bin` | `ϵ_matrix` (`nev × number of values`), `λvector` |
| `jld2 eigen` | `_eigen_data.jld2` | keys `"ϵ_matrix"`, `"λvector"` |

The sampling grid of the binary format starts at the lower-left corner of the
domain with spacing `Δx` (1D) or `(Lx/nx, Ly/ny)` (2D). The eigenfunctions are
stored as one column per eigenstate; in 2D the points are ordered with the second
coordinate running fastest.

## Reloading results

```julia
using FEMTISE
id, res = collect_result_data(true, "path/to/input")    # same path given to set_type_potential
```

- `id` is the parsed input (`InputData1D` or `InputData2D`).
- `res` is a container whose fields depend on the format:

| Container | Fields |
|---|---|
| `DefaultBinEigenProblem` | `ϵ`, `ϕ`, `r` — see the known issue below |
| `DefaultBinEigenProblemReducedDensity` | the above + `rhoDOF1`, `rhoDOF2` — see the known issue below |
| `DefaultJLD2EigenProblem` | `ϵ`, `ϕ` (`Vector{CellField}`), `r`, `pts` |
| `DefaultJLD2AllEigenProblem` | the above + `Ω`, `dΩ`, `Γ`, `dΓ`, `USpace`, `VSpace`, `model` |
| `AnalysisParams` (sweeps) | `ϵ_matrix`, `λvector` |

The eigenvalues are stored as complex numbers whose imaginary parts are at the level
of round-off; use `real.(res.ϵ)`. They are ordered by distance to the shift
`sigma`, so sort them if needed.

With `jld2 all`, the measures and spaces are rebuilt, and observables can be computed
directly:

```julia
ρ = density(res.ϕ)
x2 = coord_second_moment(res.ϕ, res.USpace, res.Ω, res.dΩ, 1)
```

!!! warning "Known issue in v0.1.0: reloading the `bin eigen` format"
    For non-sweep runs with `output_format_type = bin eigen`, `collect_result_data`
    fails with a `MethodError`: the result containers expect finite element functions,
    while the binary files contain sampled values. Read the files directly with
    [`read_bin`](@ref), as shown below. JLD2 outputs and parameter sweeps are not
    affected.

Results of the interactive console (no input file) are written in the binary
format; in v0.1.0 read them with `read_bin` as well. (`collect_result_data(false, prefix)`
is the intended reloader for that case and is affected by the same issue.)

## Reading binary files directly

```julia
ϵ = read_bin("run_eigen_values.bin")                                 # vector
ϕ = read_bin("run_eigen_vectors.bin"; matrix_data = true, c_dim = nev)  # one column per state
```

[`write_bin`](@ref) stores the length, the element type and the raw data, so files
can also be read from Python or C by following that layout.
