struct Params1D
    dimension::String
    L::Float64
    dom_type::String
    Δx::Float64
    nev::Int
    sigma::Float64
    potential_function_name::String
    params_potential::Tuple
end

struct Params2D
    dimension::String
    Lx::Float64
    Ly::Float64
    dom_type::String
    nx::Int
    ny::Int
    nev::Int
    sigma::Float64
    potential_function_name::String
    params_potential::Tuple
end

"""
    Params3D

Settings of a three-dimensional problem read from an input file (Cartesian mesh of
`nx × ny × nz` cells on a box of sides `Lx`, `Ly`, `Lz`).
"""
struct Params3D
    dimension::String
    Lx::Float64
    Ly::Float64
    Lz::Float64
    dom_type::String
    nx::Int
    ny::Int
    nz::Int
    nev::Int
    sigma::Float64
    potential_function_name::String
    params_potential::Tuple
end

"""
    AnalysisParam

Parameter sweep of an input file: the `λindex`-th potential parameter takes the
values `λi + k Δλ`, `k = 1, …, round(|λf − λi| / Δλ)`.
"""
struct AnalysisParam{T}
    λindex::Int
    λi::T
    λf::T
    Δλ::T
end

"""
    DriverOptions

Optional settings of an input file (lines after the 21 mandatory lines, introduced
in v0.2.0). Default values reproduce the behavior of earlier versions.

# Fields
- `masses`: mass of each coordinate in units of the electron mass, or `nothing`
  (electron mass for every coordinate).
- `element_order`: degree of the Lagrange elements (default 2).
- `quadrature_degree`: degree of exactness of the quadrature (default 3).
- `potential_data_file`: data file with a tabulated potential, or `""`.
- `interpolation`: `:linear` or `:cubic` interpolation of the tabulated potential.
- `analysis_param_2`: second swept parameter (`AnalysisParam`) or `false`.
- `Lz`, `nz`: third side length and number of cells of 3D problems.
"""
struct DriverOptions
    masses::Union{Nothing,Vector{Float64}}
    element_order::Int
    quadrature_degree::Int
    potential_data_file::String
    interpolation::Symbol
    analysis_param_2::Union{Bool,AnalysisParam}
    Lz::Float64
    nz::Int
end

DriverOptions() = DriverOptions(nothing, 2, 3, "", :linear, false, 0.0, 0)

struct InputData{S,T,U}
    full_path_name::String
    adhoc_file_name::String
    params::S
    analysis_param::T
    different_masses::U
end

struct InputData1D{R,S,T}
    full_path_name::String
    adhoc_file_name::String
    params::S
    analysis_param::T
    output_format_type::R
    options::DriverOptions
end
InputData1D(full_path_name, adhoc_file_name, params, analysis_param, output_format_type) =
    InputData1D(full_path_name, adhoc_file_name, params, analysis_param, output_format_type,
                DriverOptions())

struct InputData2D{R,S,T,U}
    full_path_name::String
    adhoc_file_name::String
    params::S
    analysis_param::T
    different_masses::U
    reduced_density::Bool
    output_format_type::R
    options::DriverOptions
end
InputData2D(full_path_name, adhoc_file_name, params, analysis_param, different_masses,
            reduced_density, output_format_type) =
    InputData2D(full_path_name, adhoc_file_name, params, analysis_param, different_masses,
                reduced_density, output_format_type, DriverOptions())

"""
    InputData3D

Parsed input file of a three-dimensional problem (`dimension = 3D`).
"""
struct InputData3D{R,T}
    full_path_name::String
    adhoc_file_name::String
    params::Params3D
    analysis_param::T
    output_format_type::R
    options::DriverOptions
end