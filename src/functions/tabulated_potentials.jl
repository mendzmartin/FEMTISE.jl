# =============================================================================
# Potentials given as tables of values (data files).
#
# File formats (plain text, whitespace separated, lines starting with # ignored):
#   1D: two columns    x  V(x)              x strictly increasing
#   2D: three columns  x  y  V(x, y)        all points of a regular grid, any order
#
# The finite element method needs V at arbitrary quadrature points, so the table
# is interpolated. The interpolation error is usually larger than the
# discretization error: tabulate V on a grid at least as fine as the mesh.
# =============================================================================

"""
    TabulatedPotential1D

Potential ``V(x)`` interpolated from a table (see [`tabulated_potential`](@ref)).
Callable as `V(x)` with `x` a mesh point or a number.
"""
struct TabulatedPotential1D <: Function
    x::Vector{Float64}
    values::Vector{Float64}
    second_derivatives::Vector{Float64}    # zeros for linear interpolation
    cubic::Bool
    extrapolation::Symbol
end

"""
    TabulatedPotential2D

Potential ``V(x, y)`` interpolated bilinearly from a table on a regular grid (see
[`tabulated_potential`](@ref)). Callable as `V(x)` with `x` a 2D mesh point.
"""
struct TabulatedPotential2D <: Function
    x::Vector{Float64}
    y::Vector{Float64}
    values::Matrix{Float64}                 # values[i, j] = V(x[i], y[j])
    extrapolation::Symbol
end

"""
    tabulated_potential(path; interpolation = :linear, extrapolation = :error)
    tabulated_potential(x, values; interpolation = :linear, extrapolation = :error)
    tabulated_potential(x, y, values; extrapolation = :error)

Build a potential function from tabulated values, read from the data file `path`
or given as arrays.

# File formats
- **1D**: two columns, ``x`` and ``V(x)``, with ``x`` strictly increasing.
- **2D**: three columns, ``x``, ``y`` and ``V(x,y)``, containing every point of a
  regular (tensor-product) grid in any order.

Lines starting with `#` are comments.

# Keyword arguments
- `interpolation`: `:linear` (1D and 2D, bilinear in 2D) or `:cubic` (natural cubic
  spline, 1D only). Cubic splines are recommended for smooth potentials: their
  error decreases as ``\\Delta^4`` with the spacing ``\\Delta`` of the table,
  against ``\\Delta^2`` for linear interpolation.
- `extrapolation`: `:error` (default) raises an error if the potential is evaluated
  outside the table, which usually means that the mesh is larger than the table;
  `:flat` uses the value at the nearest edge of the table.

# Example
```julia
V = tabulated_potential("potential.dat"; interpolation = :cubic)
sol = solve_tise(CartesianDiscreteModel((-8.0, 8.0), (320,)), V; nev = 5)
```
"""
function tabulated_potential(path::AbstractString; interpolation::Symbol = :linear,
                             extrapolation::Symbol = :error)
    isfile(path) || throw(ArgumentError("potential data file not found: $path"))
    data = readdlm(path, Float64; comments = true, comment_char = '#')
    columns = size(data, 2)
    if columns == 2
        return tabulated_potential(data[:, 1], data[:, 2]; interpolation, extrapolation)
    elseif columns == 3
        interpolation == :linear || throw(ArgumentError(
            "2D tables support only `interpolation = :linear` (bilinear)"))
        x = sort(unique(data[:, 1])); y = sort(unique(data[:, 2]))
        length(x) * length(y) == size(data, 1) || throw(ArgumentError(
            "the 2D table in $path is not a complete regular grid: found $(size(data, 1)) rows " *
            "for $(length(x)) × $(length(y)) distinct (x, y) values"))
        values = fill(NaN, length(x), length(y))
        for row in eachrow(data)
            values[searchsortedfirst(x, row[1]), searchsortedfirst(y, row[2])] = row[3]
        end
        any(isnan, values) && throw(ArgumentError("the 2D table in $path has repeated points"))
        return tabulated_potential(x, y, values; extrapolation)
    else
        throw(ArgumentError("the potential data file must have 2 columns (1D) or 3 columns (2D); " *
                            "$path has $columns"))
    end
end

function tabulated_potential(x::AbstractVector, values::AbstractVector;
                             interpolation::Symbol = :linear, extrapolation::Symbol = :error)
    check_options(interpolation, extrapolation)
    length(x) == length(values) || throw(ArgumentError("x and V must have the same length"))
    length(x) >= 3 || throw(ArgumentError("a tabulated potential needs at least 3 points"))
    all(diff(x) .> 0) || throw(ArgumentError("the x values of the table must be strictly increasing"))
    xs, vs = Float64.(x), Float64.(values)
    cubic = interpolation == :cubic
    second_derivatives = cubic ? natural_spline_second_derivatives(xs, vs) : zeros(length(xs))
    return TabulatedPotential1D(xs, vs, second_derivatives, cubic, extrapolation)
end

function tabulated_potential(x::AbstractVector, y::AbstractVector, values::AbstractMatrix;
                             extrapolation::Symbol = :error)
    check_options(:linear, extrapolation)
    size(values) == (length(x), length(y)) || throw(ArgumentError(
        "values must be a length(x) × length(y) matrix"))
    all(diff(x) .> 0) && all(diff(y) .> 0) || throw(ArgumentError(
        "the x and y values of the table must be strictly increasing"))
    return TabulatedPotential2D(Float64.(x), Float64.(y), Float64.(values), extrapolation)
end

function check_options(interpolation::Symbol, extrapolation::Symbol)
    interpolation in (:linear, :cubic) || throw(ArgumentError(
        "interpolation must be :linear or :cubic, got :$interpolation"))
    extrapolation in (:error, :flat) || throw(ArgumentError(
        "extrapolation must be :error or :flat, got :$extrapolation"))
end

"""
    natural_spline_second_derivatives(x, y)

Second derivatives at the nodes of the natural cubic spline through `(x, y)`
(zero second derivative at both ends), from the standard tridiagonal system.
"""
function natural_spline_second_derivatives(x::Vector{Float64}, y::Vector{Float64})
    n = length(x)
    h = diff(x)
    rhs = [6 * ((y[i+1] - y[i]) / h[i] - (y[i] - y[i-1]) / h[i-1]) for i in 2:n-1]
    system = SymTridiagonal([2 * (h[i-1] + h[i]) for i in 2:n-1], h[2:n-2])
    return vcat(0.0, system \ rhs, 0.0)
end

"Index `i` of the interval [grid[i], grid[i+1]] that contains `z`, applying the extrapolation rule."
function locate(grid::Vector{Float64}, z::Real, extrapolation::Symbol, coordinate::String)
    tolerance = 1e-10 * (grid[end] - grid[1])
    if !(grid[1] - tolerance <= z <= grid[end] + tolerance)
        extrapolation == :error && throw(DomainError(z,
            "the tabulated potential is evaluated at $coordinate = $z, outside the table " *
            "[$(grid[1]), $(grid[end])]. Make the mesh fit inside the table or use " *
            "`extrapolation = :flat`."))
        z = clamp(z, grid[1], grid[end])
    end
    i = clamp(searchsortedlast(grid, z), 1, length(grid) - 1)
    return i, z
end

function (V::TabulatedPotential1D)(point)
    i, z = locate(V.x, point[1], V.extrapolation, "x")
    h = V.x[i+1] - V.x[i]
    t = (z - V.x[i]) / h
    linear = (1 - t) * V.values[i] + t * V.values[i+1]
    V.cubic || return linear
    # cubic spline correction on [xᵢ, xᵢ₊₁] (natural spline in Hermite-free form)
    return linear + h^2 / 6 * (((1 - t)^3 - (1 - t)) * V.second_derivatives[i] +
                               (t^3 - t) * V.second_derivatives[i+1])
end

function (V::TabulatedPotential2D)(point)
    i, x = locate(V.x, point[1], V.extrapolation, "x")
    j, y = locate(V.y, point[2], V.extrapolation, "y")
    s = (x - V.x[i]) / (V.x[i+1] - V.x[i])
    t = (y - V.y[j]) / (V.y[j+1] - V.y[j])
    return (1 - s) * (1 - t) * V.values[i, j] + s * (1 - t) * V.values[i+1, j] +
           (1 - s) * t * V.values[i, j+1] + s * t * V.values[i+1, j+1]
end
