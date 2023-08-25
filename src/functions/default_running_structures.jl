struct Params1D
    dimension::String
    L::Float64
    dom_type::String
    Δx::Float64
    nev::Int64
    sigma::Float64
    potential_function_name::String
    params_potential::Tuple
end

struct Params2D
    dimension::String
    L::Float64
    dom_type::String
    nx::Int64
    ny::Int64
    nev::Int64
    sigma::Float64
    potential_function_name::String
    params_potential::Tuple
end

struct AnalysisParam
    λindex::Int64
    λi
    λf
    Δλ
end

struct InputData
    full_path_name::String
    adhoc_file_name::String
    params
    analysis_param
    different_masses
end