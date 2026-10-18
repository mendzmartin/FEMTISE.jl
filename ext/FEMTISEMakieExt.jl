"""
    FEMTISEMakieExt

Plotting functions of FEMTISE, loaded when a Makie backend (CairoMakie, GLMakie,
WGLMakie) is loaded.
"""
module FEMTISEMakieExt

using FEMTISE
using Makie

"""
    label(name, symbol, unit = "a.u.")

LaTeX label "Name, symbol (unit)", e.g. `label("Energy", "\\varepsilon_0")`.
"""
function label(name, symbol, unit = "a.u.")
    text = "\\mathrm{" * replace(name, " " => "\\ ") * "},\\ " * symbol
    isempty(unit) || (text *= "\\ \\mathrm{(" * unit * ")}")
    return Makie.LaTeXString("\$" * text * "\$")
end

"""
    log_ticks(values)

Colour-bar ticks for a logarithmic scale: values 1, 2, 5 × 10ᵏ inside the range of
`values`, refined to 1, 2, 3, 5 × 10ᵏ or 1, 1.5, 2, 3, 5, 7 × 10ᵏ for short ranges.
"""
function log_ticks(values)
    low, high = extrema(values)
    decades = floor(Int, log10(low)):ceil(Int, log10(high))
    candidates(mantissas) = [m * 10.0^k for k in decades for m in mantissas if low <= m * 10.0^k <= high]
    ticks = candidates((1, 2, 5))
    length(ticks) < 4 && (ticks = candidates((1, 2, 3, 5)))
    length(ticks) < 3 && (ticks = candidates((1, 1.5, 2, 3, 5, 7)))
    return ticks, [string(round(t; sigdigits = 2)) for t in ticks]
end

function FEMTISE.plot_energy_curves(λ, energies; states = 1:size(energies, 1),
                                    parameter_label = label("Parameter", "\\lambda", ""),
                                    energy_label = label("Energy", "\\varepsilon_n"),
                                    size = (400, 260))
    Base.size(energies, 2) == length(λ) || throw(ArgumentError(
        "energies must have one column per parameter value"))
    figure = Figure(; size)
    axis = Axis(figure[1, 1]; xlabel = parameter_label, ylabel = energy_label)
    excited = nothing
    for n in states[2:end]
        excited = lines!(axis, λ, energies[n, :]; color = :gray45, linewidth = 0.8)
    end
    ground = lines!(axis, λ, energies[states[1], :]; color = "#0072B2", linewidth = 2.0)
    xlims!(axis, extrema(λ)...)
    entries, names = Any[ground], ["ground state"]
    if !isnothing(excited)
        push!(entries, excited); push!(names, "excited states")
    end
    axislegend(axis, entries, names; position = :lt, framevisible = false)
    return figure
end

function FEMTISE.plot_energy_map(λ, μ, energies; level::Integer = 1, gap = nothing,
                                 labels = (label("First parameter", "\\lambda", ""),
                                           label("Second parameter", "\\mu", "")),
                                 energy_label = nothing, colormap = :inferno,
                                 contours::Integer = 8, colorscale = identity,
                                 size = (400, 300))
    Base.size(energies, 2) == length(λ) && Base.size(energies, 3) == length(μ) ||
        throw(ArgumentError("energies must be an nev × length(λ) × length(μ) array"))
    # states are numbered from 0 in the labels (ε₀ is the ground state)
    values, default_label = if isnothing(gap)
        energies[level, :, :], label("Energy", "\\varepsilon_{$(level - 1)}")
    else
        m, n = gap
        energies[n, :, :] .- energies[m, :, :],
        label("Energy gap", "\\varepsilon_{$(n - 1)}-\\varepsilon_{$(m - 1)}")
    end
    colorscale === log10 && any(<=(0), values) && throw(ArgumentError(
        "colorscale = log10 needs positive values (e.g. an energy gap)"))
    figure = Figure(; size)
    axis = Axis(figure[1, 1]; xlabel = labels[1], ylabel = labels[2])
    map = heatmap!(axis, λ, μ, values; colormap, colorscale, interpolate = true)
    if contours > 0
        # level lines equally spaced in the colour scale; black over a white halo so
        # that they are visible on dark and light regions
        scaled = colorscale.(values)
        contour!(axis, λ, μ, scaled; levels = contours, color = :white, linewidth = 1.6)
        contour!(axis, λ, μ, scaled; levels = contours, color = :black, linewidth = 0.6)
    end
    ticks = colorscale === log10 ? log_ticks(values) : Makie.automatic
    Colorbar(figure[1, 2], map; label = something(energy_label, default_label), ticks)
    return figure
end

end # module
