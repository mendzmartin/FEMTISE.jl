"""
    input_data(data_file_name)

# Aim
- Read the input file `data_file_name * ".dat"` of the input-file driver.

# Format
- Lines 1–21 are mandatory and are read **by position** (see the "Input file
  reference" page of the documentation): the value is the third
  whitespace-separated token of each line.
- Optional lines after line 21 (since v0.2.0) are read **by name**, as
  `key = value(s)`; lines starting with `#` are comments. Keys: `masses`,
  `element_order`, `quadrature_degree`, `potential_data_file`, `interpolation`,
  `analysis_param_2`, `Lz`, `nz` (see [`DriverOptions`](@ref)).

# Arguments
- `data_file_name::String`: name of the input file, without the `.dat` extension.

# Returns
- `InputData1D`, `InputData2D` or `InputData3D`.

# Errors
- Throws an `ArgumentError` describing the problem if the file is inconsistent.
"""
function input_data(data_file_name::String)
    attributes = readdlm("$(data_file_name).dat", String)
    size(attributes, 1) >= 21 || throw(ArgumentError(
        "$(data_file_name).dat has $(size(attributes, 1)) lines; the input file needs the 21 " *
        "mandatory lines (see the input file reference)"))
    value(row) = attributes[row, 3]
    values(row) = attributes[row, 3:end]

    full_path_name::String          = value(1)
    dom_type::String                = value(2)
    nev::Int                        = parse(Int, value(3))
    dimension::String               = value(4)
    sigma::Float64                  = parse(Float64, value(5))
    adhoc_file_name::String         = value(6)
    potential_function_name::String = value(7)
    params_potential = create_tuple(values(8), values(9))
    analysis_param = value(10) == "false" ? false : create_tuple(values(10), params_potential)
    output_format_type::Tuple{String,String} = (value(11), attributes[11, 4])
    check_output_format_type(output_format_type, analysis_param)

    options = read_driver_options(attributes, params_potential)

    if dimension == "1D"
        L::Float64 = parse(Float64, value(13))
        Δx = parse(Float64, value(14))
        params = Params1D("1D", L, dom_type, Δx, nev, sigma, potential_function_name, params_potential)
        return InputData1D(full_path_name, adhoc_file_name, params, analysis_param,
                           output_format_type, options)
    elseif dimension == "2D" || dimension == "3D"
        Lx::Float64 = parse(Float64, value(16))
        Ly::Float64 = parse(Float64, value(17))
        nx = parse(Int, value(18))
        ny = parse(Int, value(19))
        if dimension == "3D"
            options.nz > 0 && options.Lz > 0 || throw(ArgumentError(
                "3D problems need the optional lines `Lz = ...` and `nz = ...`"))
            params = Params3D("3D", Lx, Ly, options.Lz, dom_type, nx, ny, options.nz, nev, sigma,
                              potential_function_name, params_potential)
            return InputData3D(full_path_name, adhoc_file_name, params, analysis_param,
                               output_format_type, options)
        end
        params = Params2D("2D", Lx, Ly, dom_type, nx, ny, nev, sigma, potential_function_name, params_potential)
        different_masses = value(20) == "false" ? (false, nothing) : (true, parse(Float64, value(20)))
        reduced_density = value(21) != "false"
        return InputData2D(full_path_name, adhoc_file_name, params, analysis_param, different_masses,
                           reduced_density, output_format_type, options)
    else
        throw(ArgumentError("dimension must be 1D, 2D or 3D, got \"$dimension\""))
    end
end

"""
    read_driver_options(attributes, params_potential)

Optional settings in the lines after line 21 of an input file, read by key name.
Unknown keys raise an error, so that a misspelled option is never ignored silently.
"""
function read_driver_options(attributes::Matrix{String}, params_potential::Tuple)
    defaults = DriverOptions()
    masses = defaults.masses
    element_order = defaults.element_order
    quadrature_degree = defaults.quadrature_degree
    potential_data_file = defaults.potential_data_file
    interpolation = defaults.interpolation
    analysis_param_2 = defaults.analysis_param_2
    Lz, nz = defaults.Lz, defaults.nz
    for row in 22:size(attributes, 1)
        key = attributes[row, 1]
        (isempty(key) || startswith(key, "#")) && continue
        attributes[row, 2] == "=" || throw(ArgumentError(
            "optional line $row of the input file must have the form `key = value`, got `$key ...`"))
        tokens = filter(!isempty, attributes[row, 3:end])
        isempty(tokens) && continue
        if key == "masses"
            masses = parse.(Float64, tokens)
        elseif key == "element_order"
            element_order = parse(Int, tokens[1])
        elseif key == "quadrature_degree"
            quadrature_degree = parse(Int, tokens[1])
        elseif key == "potential_data_file"
            potential_data_file = tokens[1]
        elseif key == "interpolation"
            interpolation = Symbol(tokens[1])
            interpolation in (:linear, :cubic) || throw(ArgumentError(
                "interpolation must be linear or cubic, got $(tokens[1])"))
        elseif key == "analysis_param_2"
            tokens[1] == "false" || (analysis_param_2 = create_tuple(tokens, params_potential))
        elseif key == "Lz"
            Lz = parse(Float64, tokens[1])
        elseif key == "nz"
            nz = parse(Int, tokens[1])
        else
            throw(ArgumentError("unknown option `$key` in line $row of the input file. Valid options: " *
                                "masses, element_order, quadrature_degree, potential_data_file, " *
                                "interpolation, analysis_param_2, Lz, nz"))
        end
    end
    return DriverOptions(masses, element_order, quadrature_degree, potential_data_file,
                         interpolation, analysis_param_2, Lz, nz)
end

function check_output_format_type(output_format_type::Tuple{String,String}, analysis_param)
    valid = analysis_param == false ? [("bin", "eigen"), ("jld2", "eigen"), ("jld2", "all")] :
                                      [("bin", "eigen"), ("jld2", "eigen")]
    output_format_type in valid || throw(ArgumentError(
        "output_format_type = $(join(output_format_type, " ")) is not valid" *
        (analysis_param == false ? "" : " for a parameter sweep") *
        "; use one of: " * join([join(v, " ") for v in valid], ", ")))
    return true
end

function create_tuple(params_potential_type_array::Vector{String}, params_potential_array::Vector{String})
    n_types = params_quantity(params_potential_type_array)
    n_values = params_quantity(params_potential_array)
    n_types == n_values || throw(ArgumentError(
        "params_potential_types lists $n_types types but params_potential has $n_values values; " *
        "both lines need one entry per parameter"))
    return Tuple(parse(define_type(params_potential_type_array[i]), params_potential_array[i])
                 for i in 1:n_types)
end

function create_tuple(analysis_param_array::Vector{String}, params_potential::Tuple)
    data_index::Int = parse(Int, analysis_param_array[1])
    1 <= data_index <= length(params_potential) || throw(ArgumentError(
        "the swept parameter index $data_index does not exist: the potential has " *
        "$(length(params_potential)) parameters"))
    T::Type = typeof(params_potential[data_index])
    return AnalysisParam(data_index, parse(T, analysis_param_array[2]),
                         parse(T, analysis_param_array[3]), parse(T, analysis_param_array[4]))
end

function define_type(type_data::String)
    type_data == "f" && return Float64
    type_data == "c" && return ComplexF64
    type_data == "i" && return Int
    throw(ArgumentError("unknown parameter type \"$type_data\"; use f (Float64), i (Int) or c (ComplexF64)"))
end

function params_quantity(vector::Vector{String})
    return count(!isempty, vector)
end
