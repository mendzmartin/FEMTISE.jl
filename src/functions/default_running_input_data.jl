"""
    input_data(data_file_name)

# Aim
- Definition of input data form input.dat file using specific type structures.

# Arguments
- `data_file_name::String`: name of input data file
"""
function input_data(data_file_name::String)

#=
    `full_path_name::String`:           Full path name where you want to write problem results
    `L::Float64`:                       Finite element domain length [au]
    `dom_type::String`                  Domain type (symetric or non-symetric domain)
    `nev::Int`:                       Number of eigenvalues
    `dimension::String`:                Dimension of eigen value problem
    `sigma::Float64`:                   Level shift used in inverse iteration [au]
    `adhoc_file_name::String`:          Julia file name with ad hoc potential
    `potential_function_name::String`:  Name of ad hoc potential function
    `params_potential::Tuple`:          Parameters of ad hoc potential function
    only if dimension=="1D"
        `Δx::Float64`:                  Finite element size [au]
    only if dimension=="2D"
        `nx::Int`:                    Number of finite element of x direction
        `ny::Int`:                    Number of finite element of y direction
=#

    attributes = readdlm("$(data_file_name).dat",String);
    init_row::Int = 1; # without count white spaces
    init_column::Int = 3;
    full_path_name::String          = "$(attributes[init_row,init_column])"
    dom_type::String                = "$(attributes[init_row+1,init_column])"
    nev::Int                      = parse(Int,attributes[init_row+2,init_column])
    dimension::String               = "$(attributes[init_row+3,init_column])"
    sigma::Float64                  = parse(Float64,attributes[init_row+4,init_column])
    adhoc_file_name::String         = "$(attributes[init_row+5,init_column])"
    potential_function_name::String = "$(attributes[init_row+6,init_column])"

    params_potential = create_tuple(attributes[init_row+7,init_column:end],attributes[init_row+8,init_column:end])

    if "$(attributes[init_row+9,init_column])" == "false"
        analysis_param = false
    else
        analysis_param = create_tuple(attributes[init_row+9,init_column:end],params_potential)
    end
    
    if dimension == "1D"
        init_row = 12
        L::Float64 = parse(Float64,attributes[init_row,init_column])
        Δx = parse(Float64,attributes[init_row+1,init_column])

        params = Params1D("1D",L,dom_type,Δx,nev,sigma,potential_function_name,params_potential)
        data = InputData1D(full_path_name,adhoc_file_name,params,analysis_param)
    elseif dimension == "2D"
        init_row = 15
        Lx::Float64 = parse(Float64,attributes[init_row,init_column])
        Ly::Float64 = parse(Float64,attributes[init_row+1,init_column])
        nx = parse(Int,attributes[init_row+2,init_column])
        ny = parse(Int,attributes[init_row+3,init_column])
        params = Params2D("2D",Lx,Ly,dom_type,nx,ny,nev,sigma,potential_function_name,params_potential)
        if "$(attributes[init_row+4,init_column])" == "false"
            different_masses = (false,nothing)
        else
            different_masses = tuple(true,parse(Float64,attributes[init_row+4,init_column]))
        end
        if "$(attributes[init_row+5,init_column])" == "false"
            reduced_density = false
        else
            reduced_density = true
        end
        data = InputData2D(full_path_name,adhoc_file_name,params,analysis_param,different_masses,reduced_density)
    end

    return data
end

function create_tuple(params_potential_type_array::Vector{String},params_potential_array::Vector{String})
    params_potential_vec = Vector{Any}(undef,params_quantity(params_potential_type_array))
    for i in eachindex(params_potential_vec)
        params_potential_vec[i] = parse(define_type(params_potential_type_array[i]),params_potential_array[i])
    end
    return Tuple(params_potential_vec)
end

function create_tuple(analysis_param_array::Vector{String},params_potential::Tuple)
    data_index::Int=parse(Int,analysis_param_array[1])
    T::Type=typeof(params_potential[data_index])
    analysis_param = AnalysisParam(data_index,parse(T,analysis_param_array[2]),
        parse(T,analysis_param_array[3]),parse(T,analysis_param_array[4]))
    return analysis_param
end

function define_type(type_data::String)
    if type_data == "f"
        T = Float64
    elseif type_data == "c"
        T = ComplexF64
    elseif type_data == "i"
        T = Int
    end
    return T::Type
end

function params_quantity(vector::Vector{String})
    val::Int=0
    for i in eachindex(vector)
        vector[i]=="" ? continue : nothing
        val+=1
    end
    return val
end