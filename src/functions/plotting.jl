# =============================================================================
# Plotting functions. They are implemented in the package extension
# `FEMTISEMakieExt` (ext/FEMTISEMakieExt.jl), which is loaded automatically when a
# Makie backend is loaded, e.g. `using CairoMakie`. FEMTISE itself does not
# depend on any plotting package.
# =============================================================================

"""
    plot_energy_curves(λ, energies; states = 1:size(energies, 1),
                       parameter_label = "Parameter, λ", energy_label = "Energy, εₙ (a.u.)",
                       size = (400, 260))

Eigenvalues as functions of one parameter, `energies[n, i]` for the values `λ[i]`
(e.g. from [`sweep_eigenvalues`](@ref) or a one-parameter sweep of the input-file
driver). The ground state is drawn as a thick coloured line and the excited states
as thin gray lines, so the figure can be printed in black and white. `size` is the
figure size in points.

Requires a Makie backend: `using CairoMakie`. Returns a `Makie.Figure`.
"""
function plot_energy_curves end

"""
    plot_energy_map(λ, μ, energies; level = 1, gap = nothing,
                    labels = ("First parameter, λ", "Second parameter, μ"),
                    energy_label = nothing, colormap = :inferno, contours = 8,
                    colorscale = identity, size = (400, 300))

Two-dimensional map of an eigenvalue over two parameters, from a two-parameter
sweep (`energies[n, i, j]` for `λ[i]`, `μ[j]`):
- `gap = nothing`: the eigenvalue `energies[level, :, :]`;
- `gap = (m, n)`: the energy difference `energies[n, :, :] - energies[m, :, :]`
  (e.g. `(1, 2)` for the gap between the two lowest levels).

States are numbered from 0 in the default colour-bar label (``\\varepsilon_0`` is
the ground state); `energy_label` replaces it.

`colorscale = log10` shows quantities that vary over orders of magnitude (e.g.
tunnelling splittings); contour lines are then equally spaced in ``\\log_{10}``.

The map uses a perceptually uniform colour map with monotonic lightness and
`contours` level lines (black over a white halo, visible on any background), so
that it remains readable in black and white.

Requires a Makie backend: `using CairoMakie`. Returns a `Makie.Figure`.
"""
function plot_energy_map end
