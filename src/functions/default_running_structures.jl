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

struct AnalysisParam{T}
    λindex::Int
    λi::T
    λf::T
    Δλ::T
end

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
end

struct InputData2D{R,S,T,U}
    full_path_name::String
    adhoc_file_name::String
    params::S
    analysis_param::T
    different_masses::U
    reduced_density::Bool
    output_format_type::R
end