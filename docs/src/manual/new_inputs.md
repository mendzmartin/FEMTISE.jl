# Masses, potential data files and parameter sweeps

This page describes three features introduced in v0.2.0, available both in the
low-level interface and in the input-file driver:

- particles of any mass;
- potentials given as tables of values (data files);
- eigenvalues as functions of one or two potential parameters, with
  two-dimensional plots.

## The general solver `solve_tise`

[`solve_tise`](@ref) is the recommended entry point of the low-level interface. It
works on 1D, 2D and 3D meshes and returns the eigenvalues sorted by energy with
their eigenfunctions:

```julia
using FEMTISE, Gridap
model = CartesianDiscreteModel((-10.0, 10.0), (200,))
sol = solve_tise(model, x -> 0.5 * x[1]^2; masses = 1.0, order = 2, nev = 10, sigma = 0.0)
sol.energies        # Vector{Float64}, ascending
sol.states          # eigenfunctions (normalized finite element functions)
sol.dΩ, sol.trial   # measure and space, for post-processing
```

By default the quadrature integrates polynomials of degree ``2\ell+2`` exactly, so for
potentials up to quadratic order the eigenvalues are guaranteed upper bounds (see
[Theory](@ref)).

## Particle masses

Masses are given in units of the electron mass (atomic units: ``m_e = 1``). Each
coordinate can have its own mass.

| System | `masses` |
|---|---|
| Electron (default) | `1.0` |
| Muon | `206.77` |
| Proton | `1836.15` |
| H₂ vibration (reduced mass) | `918.08` |
| Electron in GaAs (effective mass) | `0.067` |
| Two particles, masses ``M_1, M_2`` (2D mesh) | `(M₁, M₂)` |
| Nucleus, electron, photon (3D mesh) | `(M_N, 1.0, 1.0)` |

In an input file, add the optional line `masses = 1836.15` (one value for all
coordinates) or one value per coordinate, e.g. `masses = 1.0 4.0`. For
position-dependent effective masses use the low-level interface with a non-constant
kinetic coefficient (see [Low-level interface](@ref)).

## Potentials from data files

[`tabulated_potential`](@ref) reads a potential tabulated on a grid and interpolates
it:

| Format | Content | Interpolation |
|---|---|---|
| 1D | two columns, ``x`` (increasing) and ``V(x)`` | `:linear` or `:cubic` (natural spline) |
| 2D | three columns, ``x``, ``y``, ``V(x,y)``, all points of a regular grid in any order | bilinear |

Lines starting with `#` are comments. Typical sources are *ab initio* potential
energy curves or surfaces, or measured profiles.

```julia
V = tabulated_potential("pec.dat"; interpolation = :cubic)
sol = solve_tise(CartesianDiscreteModel((0.8, 6.0), (400,)), V; masses = 918.08, nev = 10,
                 sigma = minimum(V.values))
```

In an input file, leave `adhoc_file_name`, `potential_function_name`,
`params_potential_types` and `params_potential` empty and add

```text
## OPTIONAL SETTINGS
potential_data_file         = ./pec.dat
interpolation               = cubic
```

!!! warning "Accuracy of tabulated potentials"
    The interpolation error usually dominates over the finite element error. For the
    1D harmonic oscillator tabulated every ``0.1`` a.u., linear interpolation gives
    eigenvalues with relative errors of ``\sim10^{-3}``, cubic splines ``\sim10^{-6}``,
    and the analytic potential ``\sim10^{-7}`` on the same mesh. Prefer cubic splines
    for smooth potentials and a table at least as fine as the mesh. The table must
    cover the whole mesh: evaluating outside it raises an error, unless
    `extrapolation = :flat` is chosen.

## Parameter sweeps and energy maps

[`sweep_eigenvalues`](@ref) computes the eigenvalues as functions of one or two
parameters of the potential, reusing the mesh and the finite element spaces:

```julia
V(x, λ, μ) = λ * (x[1]^2 - 1)^2 + μ * x[1]          # double well: barrier λ, tilt μ
λ = 2.0:0.25:8.0
μ = 0.0:0.05:1.0
E = sweep_eigenvalues(model, V, λ, μ; nev = 4, sigma = -2.0)   # E[n, i, j]
```

In an input file, the first parameter is set with `analysis_param` (line 10) and the
second one with the optional line

```text
analysis_param_2            = 2 0.0 1.0 0.05
```

with the same meaning (index of the parameter, initial and final value, step). The
results are reloaded with `collect_result_data` as an `AnalysisParams2D` with fields
`ϵ_tensor`, `λvector` and `μvector`.

### Plots

With a Makie backend loaded (`using CairoMakie`), FEMTISE provides two plotting
functions:

```julia
using CairoMakie
fig = plot_energy_map(λ, μ, E; gap = (1, 2),
                      labels = ("Barrier height, λ (a.u.)", "Tilt, μ (a.u.)"))
save("tunnel_splitting.pdf", fig)
plot_energy_curves(λ, E[:, :, 1])    # one-parameter cut
```

!!! note "Name clash with Makie"
    Makie also exports a function called `density`. When both packages are loaded
    with `using`, write `FEMTISE.density(ϕ)` for the probability densities.

[`plot_energy_map`](@ref) draws an eigenvalue or a gap ``\varepsilon_n-\varepsilon_m``
over the two parameters, with a perceptually uniform colour map and level lines.
[`plot_energy_curves`](@ref) draws the levels along one parameter. Both are designed
to remain readable in black and white.

## Other optional settings of the input file

| Key | Meaning | Default |
|---|---|---|
| `element_order` | degree of the Lagrange elements | 2 |
| `quadrature_degree` | degree of exactness of the quadrature | 3 |
| `Lz`, `nz` | third side and number of cells of 3D problems | — |

Input files that use any optional setting are solved with the general solver and
must use `output_format_type = jld2 eigen`. Their eigenvalues are sorted in ascending
order. Input files without optional settings behave exactly as in v0.1.
