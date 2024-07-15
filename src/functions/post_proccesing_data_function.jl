"""
    charge_input_data(full_path_input_data)

# Aim
- Charge input data from a file

# Arguments
- `full_path_input_data::String`: full path name to input data file
"""
function charge_input_data(full_path_input_data::String)
    return input_data(full_path_input_data)
end

"""
    build_input_data(full_path_name)

# Aim
- Build input data from a file

# Arguments
- `full_path_name::String`: full path name to input data file
"""
function build_input_data(full_path_name::String)
    attributes=readdlm(full_path_name*"_eigen_problem_attributes.dat",String);
    dimension=attributes[2,8];
    if dimension == "1D"
        params=Params1D(dimension,parse(Float64,attributes[4,8]),attributes[5,5],parse(Float64,attributes[7,7]),
            parse(Float64,attributes[3,6]),parse(Float64,attributes[6,10]),"",tuple(nothing))
    elseif dimension == "2D"
        params=Params2D(dimension,parse(Float64,attributes[4,9]),parse(Float64,attributes[4,10]),attributes[5,5],parse(Int64,attributes[7,10]),
            parse(Int64,attributes[8,10]),parse(Float64,attributes[3,6]),parse(Float64,attributes[6,10]),
            "",tuple(nothing))
    end
    adhoc_file_name = ""; analysis_param = false; different_masses = false;
    return InputData(full_path_name,adhoc_file_name,params,analysis_param,different_masses);
end

"""
    charge_results(id)

# Aim
- Charge results from a file

# Arguments
- `id::InputData`: input data
"""
function charge_results(id::InputData)
    if id.params.dimension == "1D"
        eigen_vectors_output= id.full_path_name*"_eigen_vectors.bin"
        eigen_values_output= id.full_path_name*"_eigen_values.bin"
        coordinates_output = id.full_path_name*"_coordinates.bin"
        eigen_energies=read_bin(eigen_values_output;matrix_data=false)
        eigen_states=read_bin(eigen_vectors_output;matrix_data=true,c_dim=id.params.nev)
        coordinates=read_bin(coordinates_output;matrix_data=false,c_dim=1);
    elseif id.params.dimension == "2D"
        eigen_vectors_output= id.full_path_name*"_eigen_vectors.bin"
        eigen_values_output= id.full_path_name*"_eigen_values.bin"
        coordinates_output = id.full_path_name*"_coordinates.bin"
        eigen_energies=read_bin(eigen_values_output;matrix_data=false)
        eigen_states=read_bin(eigen_vectors_output;matrix_data=true,c_dim=id.params.nev)
        coordinates=read_bin(coordinates_output;matrix_data=false,c_dim=1);
    end
    return DefaultBinEigenProblem(eigen_energies,eigen_states,coordinates)
end

"""
    interpolation_eigenstates!(eigen_states,USpace)

# Aim
- Interpolate eigenstates

# Arguments
- `eigen_states::Vector{CellField}`: eigenstates
- `USpace::FESpace`: FE space
"""
function interpolation_eigenstates!(eigen_states::Vector{CellField},USpace::FESpace)
    for i in eachindex(eigen_states)
        eigen_state_i = Interpolable(eigen_states[i])
        eigen_states[i] = interpolate_everywhere(eigen_state_i,USpace)
    end
    return eigen_states
end

"""
    triangulation_repair(model,grid_type)

# Aim
- Repair triangulation

# Arguments
- `model::Model`: FE grid model
- `grid_type::String`: FE grid type (could be "simple_line" or "Cartesian2D")
"""
function triangulation_repair(model,grid_type::String)
    FullDirichlet_values,FullDirichlet_tags=make_boundary_conditions(grid_type,"FullDirichlet",ComplexF64);
    Ω,dΩ,Γ,dΓ=measures(model,3,FullDirichlet_tags);
    reff = ReferenceFE(lagrangian,Float64,2);
    VSpace,USpace=fe_spaces(model,reff;BC_data=(FullDirichlet_values,FullDirichlet_tags),BC_type="Dirichlet")
    return Ω,dΩ,Γ,dΓ,VSpace,USpace
end


"""
    charge_results(id)

# Aim
- Charge results from a file

# Arguments
- `id::InputData1D`: input data
"""
function charge_results(id::InputData1D)
    if id.analysis_param == false
        if id.output_format_type == ("bin","eigen") # use DelimitedFiles package
            eigen_vectors_output= id.full_path_name*"_eigen_vectors.bin"
            eigen_values_output= id.full_path_name*"_eigen_values.bin"
            coordinates_output = id.full_path_name*"_coordinates.bin"
            ϵ=read_bin(eigen_values_output;matrix_data=false)
            ϕ=read_bin(eigen_vectors_output;matrix_data=true,c_dim=id.params.nev)
            r=read_bin(coordinates_output;matrix_data=false,c_dim=1);
            result = DefaultBinEigenProblem(ϵ,ϕ,r)
        elseif id.output_format_type == ("jld2","eigen") # use JLD2 package
            ϵ = load("$(id.full_path_name)_eigen_data.jld2", "ϵ")
            ϕ = load("$(id.full_path_name)_eigen_data.jld2", "ϕ")
            r = load("$(id.full_path_name)_eigen_data.jld2", "r")
            pts = load("$(id.full_path_name)_eigen_data.jld2", "pts")
            result = DefaultJLD2EigenProblem(ϵ,ϕ,r,pts)
        elseif id.output_format_type == ("jld2","all")
            ϵ = load("$(id.full_path_name)_eigen_data.jld2", "ϵ")
            ϕ = load("$(id.full_path_name)_eigen_data.jld2", "ϕ")
            r = load("$(id.full_path_name)_eigen_data.jld2", "r")
            pts = load("$(id.full_path_name)_eigen_data.jld2", "pts")
            model = load("$(id.full_path_name)_eigen_data.jld2", "model")
            Ω,dΩ,Γ,dΓ,VSpace,USpace = triangulation_repair(model,"simple_line")
            ϕ = interpolation_eigenstates!(ϕ,USpace)
            result = DefaultJLD2AllEigenProblem(ϵ,ϕ,r,pts,Ω,dΩ,Γ,dΓ,USpace,VSpace,model)
        end
    elseif id.analysis_param ≠ false
        if id.output_format_type == ("jld2","eigen") # use JLD2 package
            ϵ_matrix = load("$(id.full_path_name)_eigen_data.jld2", "ϵ_matrix")
            λvector = load("$(id.full_path_name)_eigen_data.jld2", "λvector")
        elseif id.output_format_type == ("bin","eigen") # use DelimitedFiles package
            eigen_values_output= id.full_path_name*"_eigen_values.bin";
            param_values_output = id.full_path_name*"_param_values.bin";
            λvector=read_bin(param_values_output;matrix_data=false);
            ϵ_matrix=read_bin(eigen_values_output;matrix_data=true,c_dim=length(λvector));
        end
        result = AnalysisParams(ϵ_matrix,λvector)
    end
    return result
end

"""
    charge_results(id)

# Aim
- Charge results from a file

# Arguments
- `id::InputData2D`: input data
"""
function charge_results(id::InputData2D)
    if id.analysis_param == false
        if id.output_format_type == ("bin","eigen") # use DelimitedFiles package
            eigen_vectors_output= id.full_path_name*"_eigen_vectors.bin"
            eigen_values_output= id.full_path_name*"_eigen_values.bin"
            coordinates_output = id.full_path_name*"_coordinates.bin"
            ϵ=read_bin(eigen_values_output;matrix_data=false)
            ϕ=read_bin(eigen_vectors_output;matrix_data=true,c_dim=id.params.nev)
            r=read_bin(coordinates_output;matrix_data=false,c_dim=1);
            if id.reduced_density
                reduced_density_DOF1_output = id.full_path_name*"_reduced_density_DOF1.bin";
                reduced_density_DOF2_output = id.full_path_name*"_reduced_density_DOF2.bin";
                rhoDOF1=read_bin(reduced_density_DOF1_output;matrix_data=true,c_dim=id.params.nev);
                rhoDOF2=read_bin(reduced_density_DOF2_output;matrix_data=true,c_dim=id.params.nev);
            end
            id.reduced_density ? (result = DefaultBinEigenProblemReducedDensity(ϵ,ϕ,r,rhoDOF1,rhoDOF2)) : 
                (result = DefaultBinEigenProblem(ϵ,ϕ,r))
        elseif id.output_format_type == ("jld2","eigen") # use JLD2 package
            ϵ = load("$(id.full_path_name)_eigen_data.jld2", "ϵ")
            ϕ = load("$(id.full_path_name)_eigen_data.jld2", "ϕ")
            r = load("$(id.full_path_name)_eigen_data.jld2", "r")
            pts = load("$(id.full_path_name)_eigen_data.jld2", "pts")
            if id.reduced_density
                rhoDOF1 = load("$(id.full_path_name)_eigen_data.jld2", "rhoDOF1")
                rhoDOF2 = load("$(id.full_path_name)_eigen_data.jld2", "rhoDOF2")
            end
            id.reduced_density ? (result = DefaultJLD2EigenProblemReducedDensity(ϵ,ϕ,r,pts,rhoDOF1,rhoDOF2)) : 
                (result = DefaultJLD2EigenProblem(ϵ,ϕ,r,pts))
        elseif id.output_format_type == ("jld2","all")
            ϵ = load("$(id.full_path_name)_eigen_data.jld2", "ϵ")
            ϕ = load("$(id.full_path_name)_eigen_data.jld2", "ϕ")
            r = load("$(id.full_path_name)_eigen_data.jld2", "r")
            pts = load("$(id.full_path_name)_eigen_data.jld2", "pts")
            model = load("$(id.full_path_name)_eigen_data.jld2", "model")
            Ω,dΩ,Γ,dΓ,VSpace,USpace = triangulation_repair(model,"Cartesian2D")
            ϕ = interpolation_eigenstates!(ϕ,USpace)
            if id.reduced_density
                rhoDOF1 = load("$(id.full_path_name)_eigen_data.jld2", "rhoDOF1")
                rhoDOF2 = load("$(id.full_path_name)_eigen_data.jld2", "rhoDOF2")
            end
            id.reduced_density ? (result = DefaultJLD2AllEigenProblemReducedDensity(ϵ,ϕ,r,pts,Ω,dΩ,Γ,dΓ,USpace,VSpace,model,rhoDOF1,rhoDOF2)) : 
                (result = DefaultJLD2AllEigenProblem(ϵ,ϕ,r,pts,Ω,dΩ,Γ,dΓ,USpace,VSpace,model))
        end
    else
        if id.output_format_type == ("jld2","eigen") # use JLD2 package
            ϵ_matrix = load("$(id.full_path_name)_eigen_data.jld2", "ϵ_matrix")
            λvector = load("$(id.full_path_name)_eigen_data.jld2", "λvector")
        elseif id.output_format_type == ("bin","eigen") # use DelimitedFiles package
            eigen_values_output= id.full_path_name*"_eigen_values.bin";
            param_values_output = id.full_path_name*"_param_values.bin";
            λvector=read_bin(param_values_output;matrix_data=false);
            ϵ_matrix=read_bin(eigen_values_output;matrix_data=true,c_dim=length(λvector));
        end
        result = AnalysisParams(ϵ_matrix,λvector)
    end
    return result
end

"""
    collect_result_data(switch_input_file,full_path_name)

# Aim
- Collect result data

# Arguments
- `switch_input_file::Bool`: switch input file
- `full_path_name::String`: full path name
"""
function collect_result_data(switch_input_file::Bool,full_path_name::String)
    switch_input_file ? (id = charge_input_data(full_path_name)) : (id = build_input_data(full_path_name))
    return id, charge_results(id)
end