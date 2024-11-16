struct OutputData{R,S}
    data_type::String
    output_format_type::S
    file_name::String
    data::R
end

function remove_and_write_data(output::OutputData)
    println("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    println("Moving result data to trash if data exists ...");
    rm_existing_file(output.file_name)
    println("Saving data ...")

    if output.data_type == "bin_default"
        output.output_format_type == ("bin","eigen") ? write_bin(output.data,output.file_name) : nothing
    elseif output.data_type == "jld2_analysis_param"
        jldsave(output.file_name;ϵ_matrix=output.data[1],λvector=output.data[2])
    elseif output.data_type == "jld2_reduced_density"
        if output.output_format_type == ("jld2","eigen")
            jldsave(output.file_name;ϵ=output.data[1],ϕ=output.data[2],r=output.data[3],pts=output.data[4],
                rhoDOF1=output.data[5],rhoDOF2=output.data[6])
        elseif output.output_format_type == ("jld2","all")
            jldsave(output.file_name;ϵ=output.data[1],ϕ=output.data[2],r=output.data[3],pts=output.data[4],
                model=output.data[5],rhoDOF1=output.data[6],rhoDOF2=output.data[7])
        end
    elseif output.data_type == "jld2_default"
        if output.output_format_type == ("jld2","eigen")
            jldsave(output.file_name;ϵ=output.data[1],ϕ=output.data[2],r=output.data[3],pts=output.data[4])
        elseif output.output_format_type == ("jld2","all")
            jldsave(output.file_name;ϵ=output.data[1],ϕ=output.data[2],r=output.data[3],pts=output.data[4],
                model=output.data[5])
        end
    end
    println("Saved data.")
end

function dispatch_solver(dimension::String,type_simulation::String,parameters)
    if dimension == "1D"
        if type_simulation == "analysis_params"
            params,model = parameters
            eigen_results = solver_eigen_problem_with_analysis_param(params,model)
        elseif type_simulation == "default"
            params = parameters
            eigen_results = default_solver_eigen_problem(params)
        elseif type_simulation == "full_export"
            params = parameters
            eigen_results = default_solver_eigen_problem(params;export_all_data=true)
        end
    elseif dimension == "2D"
        if type_simulation == "analysis_params"
            params,different_masses,model = parameters
            eigen_results = solver_eigen_problem_with_analysis_param(params,different_masses,model)
        elseif type_simulation == "reduced_density"
            params,different_masses,switch_reduced_density = parameters
            eigen_results = default_solver_eigen_problem(params,different_masses;
                switch_reduced_density=switch_reduced_density)
        elseif type_simulation == "default"
            params,different_masses = parameters
            eigen_results = default_solver_eigen_problem(params,different_masses)
        elseif type_simulation == "full_export"
            params,different_masses = parameters
            eigen_results = default_solver_eigen_problem(params,different_masses;export_all_data=true)
        end
    end
    return eigen_results
end

"""
    run_default_eigen_problem

# Aim
- Function to simulate eigen problems from default potential or custom potential
"""

function run_default_eigen_problem(simulation_data::Tuple)
    type_potential,id=simulation_data
    println("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    if type_potential ≠ "5"
        print("Set full path name (e.g: \"./my_directory/my_name\") where you want to write problem results and press Enter = ")
        full_path_name = readline()
        io=open(full_path_name*"_eigen_problem_attributes.dat","w")

        println("Mandatory input data")
        print("Number of eigenvalues: nev::Int = ")
        nev = get_input(Int;default_data=false)

        if type_potential in ["1","2","3"]
            print("Finite element domain length [au] (default L=30.0 press Enter): L::Float64 = ")
            L = get_input(Float64,default_value=30.0)
            print("Set domain type (for symetric domain {-L/2,L/2} set \"s\" and for non-symetric domain {0,L} set \"ns\": dom_type::String = ")
            dom_type = get_input(["s","ns"])
            dimension="1D"
            print("Finite element size [au] (default Δx=0.1 press Enter): Δx::Float64 = ")
            Δx = get_input(Float64;default_value=0.1)
            output_format_type = ("bin","eigen")
            if type_potential=="1"
                potential_function_name="qho_1d"
                print("Harmonic Oscillator frecuency [au]: ω::Float64 = ")
                ω = get_input(Float64;default_data=false)
                print("Harmonic Oscillator center [au]: x₁::Float64 = ")
                x₁ = get_input(Float64;default_data=false)
                print("Level shift used in inverse iteration [au] (default sigma=0.0 press Enter): sigma::Float64 = ")
                sigma = get_input(Float64;default_value=0.0)
                params_potential=(ω,x₁)
                write(io,"Quantum Harmonic Oscillator 1D\n")
                write(io,"Dimension of eigen value problem              dimension::String   = $(dimension)\n")
                write(io,"Number of eigenvalues                         nev::Int          = $(nev)\n")
                write(io,"Finite element domain length [au]             L::Float64          = $(L)\n")
                write(io,"Domain type                                   dom_type::String    = $(dom_type)\n")
                write(io,"Level shift used in inverse iteration [au]    sigma::Float64      = $(sigma)\n")
                write(io,"Finite element size [au]                      Δx::Float64         = $(Δx)\n")
                write(io,"Harmonic Oscillator frecuency [au]            ω::Float64          = $(ω)\n")
                write(io,"Harmonic Oscillator center [au]               x₁::Float64         = $(x₁)\n")
            elseif type_potential=="2"
                potential_function_name="kronig_penney_1d"
                print("Depths (with sign) of potential wells [au]: V₀::Float64 = ")
                V₀ = get_input(Float64;default_data=false)
                print("Distance between well centers [au]: t::Float64 = ")
                t=get_input(Float64;default_data=false)
                print("Widths of potential wells [au]: b::Float64 = ")
                b=get_input(Float64;default_data=false)
                print("Number of ions or potential wells: num_ions::Int = ")
                num_ions=get_input(Int,num_ions_BoolCondition)
                print("Level shift used in inverse iteration [au] (default sigma=V₀ press Enter): sigma::Float64 = ")
                sigma=get_input(Float64;default_value=V₀)
                params_potential=(V₀,t,b,num_ions)
                write(io,"Finite Kronig-Penney 1D\n")
                write(io,"Dimension of eigen value problem              dimension::String   = $(dimension)\n")
                write(io,"Number of eigenvalues                         nev::Int          = $(nev)\n")
                write(io,"Finite element domain length [au]             L::Float64          = $(L)\n")
                write(io,"Domain type                                   dom_type::String    = $(dom_type)\n")
                write(io,"Level shift used in inverse iteration [au]    sigma::Float64      = $(sigma)\n")
                write(io,"Finite element size [au]                      Δx::Float64         = $(Δx)\n")
                write(io,"Depths of potential wells [au]                V₀::Float64         = $(V₀)\n")
                write(io,"Distance between well centers [au]            t::Float64          = $(t)\n")
                write(io,"Widths of potential wells [au]                b::Float64          = $(b)\n")
                write(io,"Number of ions or potential wells             num_ions::Int     = $(num_ions)\n")
            elseif type_potential=="3"
                potential_function_name="finite_well_1d"
                print("Depth (with sign) of potential well [au]: V₀::Float64 = ")
                V₀ = get_input(Float64;default_data=false)
                print("Width of potential well [au]: b::Float64 = ")
                b=get_input(Float64;default_data=false)
                print("Level shift used in inverse iteration [au] (default sigma=V₀ press Enter): sigma::Float64 = ")
                sigma=get_input(Float64;default_value=V₀)
                params_potential=(V₀,b)
                write(io,"Finite Well Potential 1D\n")
                write(io,"Dimension of eigen value problem              dimension::String   = $(dimension)\n")
                write(io,"Number of eigenvalues                         nev::Int          = $(nev)\n")
                write(io,"Finite element domain length [au]             L::Float64          = $(L)\n")
                write(io,"Domain type                                   dom_type::String    = $(dom_type)\n")
                write(io,"Level shift used in inverse iteration [au]    sigma::Float64      = $(sigma)\n")
                write(io,"Finite element size [au]                      Δx::Float64         = $(Δx)\n")
                write(io,"Depth of potential well [au]                  V₀::Float64         = $(V₀)\n")
                write(io,"Width of potential well [au]                  b::Float64          = $(b)\n")
            end
            params=Params1D(dimension,L,dom_type,Δx,nev,sigma,potential_function_name,params_potential)
        elseif type_potential=="4"
            print("Finite element domain length of DOF1 [au] (default L=30.0 press Enter): Lx::Float64 = ")
            Lx = get_input(Float64,default_value=30.0)
            print("Finite element domain length of DOF2 [au] (default L=30.0 press Enter): Ly::Float64 = ")
            Ly = get_input(Float64,default_value=30.0)
            print("Set domain type (for symetric domain {-L/2,L/2} set \"s\" and for non-symetric domain {0,L} set \"ns\": dom_type::String = ")
            dom_type = get_input(["s","ns"])
            potential_function_name="qho_2d"
            dimension="2D"
            print("Number of finite element of x direction (default nx=100 press Enter): nx::Int = ")
            nx=get_input(Int;default_value=100)
            print("Number of finite element of y direction (default ny=100 press Enter): ny::Int = ")
            ny=get_input(Int;default_value=100)
            print("Harmonic Oscillator frecuency [au]: ω::Float64 = ")
            ω = get_input(Float64;default_data=false)
            print("Harmonic Oscillator center of x direction [au]: x₁::Float64 = ")
            x₁=get_input(Float64;default_data=false)
            print("Harmonic Oscillator center of y direction [au]: y₁::Float64 = ")
            y₁=get_input(Float64;default_data=false)
            print("Level shift used in inverse iteration [au] (default sigma=0.0 press Enter): sigma::Float64 = ")
            sigma=get_input(Float64;default_value=0.0)
            params_potential=(ω,x₁,y₁)
            params=Params2D(dimension,Lx,Ly,dom_type,nx,ny,nev,sigma,potential_function_name,params_potential)
            different_masses = tuple(false, nothing)
            switch_reduced_density = false
            output_format_type = ("bin","eigen")
            write(io,"Quantum Harmonic Oscillator 2D\n")
            write(io,"Dimension of eigen value problem                  dimension::String       = $(dimension)\n")
            write(io,"Number of eigenvalues                             nev::Int                = $(nev)\n")
            write(io,"Finite element domain lengths [au]                Lx::Float64 Ly::Float64 = $(Lx) $(Ly)\n")
            write(io,"Domain type                                       dom_type::String        = $(dom_type)\n")
            write(io,"Level shift used in inverse iteration [au]        sigma::Float64          = $(sigma)\n")
            write(io,"Number of finite element of x direction           nx::Int                 = $(nx)\n")
            write(io,"Number of finite element of y direction           ny::Int                 = $(ny)\n")
            write(io,"Harmonic Oscillator frecuency [au]                ω::Float64              = $(ω)\n")
            write(io,"Harmonic Oscillator center of x direction [au]    x₁::Float64             = $(x₁)\n")
            write(io,"Harmonic Oscillator center of y direction [au]    y₁::Float64             = $(y₁)\n")
        end
        write(io,"\nNumber of threads                                                           = $(Threads.nthreads())\n")
        close(io)
    end

    if (type_potential=="5" && id.analysis_param ≠ false)
        λvector=[id.analysis_param.λi+i*id.analysis_param.Δλ 
            for i in 1:round(Int,abs(id.analysis_param.λf-id.analysis_param.λi)/id.analysis_param.Δλ)]
        ϵ_matrix=Matrix{ComplexF64}(undef,id.params.nev,length(λvector))

        model=create_and_remove_model(id.params)

        for i in eachindex(λvector)
            params_potential_array=Vector{Any}(undef, length(id.params.params_potential))
            for j in eachindex(id.params.params_potential)
                if (j == id.analysis_param.λindex)
                    params_potential_array[j] = λvector[i] # (λ=λi+i*Δλ)
                else
                    params_potential_array[j] = id.params.params_potential[j]
                end
            end
            params_potential = Tuple(params_potential_array)
            if id.params.dimension=="1D"
                paramsλ=Params1D(id.params.dimension,id.params.L,id.params.dom_type,
                    id.params.Δx,id.params.nev,id.params.sigma,
                    id.params.potential_function_name,params_potential)

                    println("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
                    println("... Running for λ$(i) of $(length(λvector)) λ's...")
                    parameters = (paramsλ,model)
            elseif id.params.dimension=="2D"
                paramsλ=Params2D(id.params.dimension,id.params.Lx,id.params.Ly,id.params.dom_type,
                    id.params.nx,id.params.ny,id.params.nev,id.params.sigma,
                    id.params.potential_function_name,params_potential)

                    println("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
                    println("... Running for λ$(i) of $(length(λvector)) λ's...")
                    parameters = (paramsλ,id.different_masses,model)
            end
            ϵ,ϕ = dispatch_solver(id.params.dimension,"analysis_params",parameters)
            ϵ_matrix[:,i] = ϵ[:]
        end

        if (id.output_format_type)== ("bin","eigen")
            eigen_values_output = id.full_path_name*"_eigen_values.bin"
            param_values_output = id.full_path_name*"_param_values.bin"
            output=OutputData("bin_default",id.output_format_type,eigen_values_output,ϵ_matrix)
            remove_and_write_data(output)
            output=OutputData("bin_default",id.output_format_type,param_values_output,λvector)
            remove_and_write_data(output)

        elseif (id.output_format_type)== ("jld2","eigen")
            output_file_name = id.full_path_name*"_eigen_data.jld2"
            output=OutputData("jld2_analysis_param",id.output_format_type,output_file_name,(ϵ_matrix,λvector))
            remove_and_write_data(output)
        end
    elseif (type_potential=="5" && id.analysis_param == false)
        if id.params.dimension == "2D"
            switch_reduced_density=id.reduced_density
        end

        println("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
        println("... Running ...")

        if id.params.dimension == "1D"
            parameters = id.params
            switch_reduced_density=false
            if id.params.dom_type=="s"
                dom=(-0.5*id.params.L,0.5*id.params.L)
            elseif id.params.dom_type=="ns"
                dom=(0.0,id.params.L)
            end
            r,pts=space_coord(dom,id.params.Δx,round(Int,abs(id.params.L/id.params.Δx));dimension="1D")
        elseif id.params.dimension == "2D"
            if (id.output_format_type)==("jld2","all")
                parameters = (id.params,id.different_masses,switch_reduced_density)
            else
                parameters = (id.params,id.different_masses,switch_reduced_density)
            end
            if id.params.dom_type=="s"
                dom=(-0.5*id.params.Lx,0.5*id.params.Lx,-0.5*id.params.Ly,0.5*id.params.Ly)
            elseif params.dom_type=="ns"
                dom=(0.0,id.params.Lx,0.0,id.params.Ly)
            end
            r,pts=space_coord(dom,(id.params.Lx/id.params.nx,id.params.Ly/id.params.ny),(id.params.nx,id.params.ny))
        end

        if (id.output_format_type)==("jld2","all")
            ϵ,ϕ,model = dispatch_solver(id.params.dimension,"full_export",parameters)
        else
            if switch_reduced_density
                ϵ,ϕ,model = dispatch_solver(id.params.dimension,"reduced_density",parameters)
            else
                ϵ,ϕ = dispatch_solver(id.params.dimension,"default",parameters)
            end
        end

        if switch_reduced_density && id.params.dimension == "2D"
            rho_DOF1_matrix,rho_DOF2_matrix = reduced_density(ϕ,r,model)
        end

        if (id.output_format_type)== ("bin","eigen")
            eigen_vectors_output = id.full_path_name*"_eigen_vectors.bin"
            eigen_values_output = id.full_path_name*"_eigen_values.bin"
            coordinates_output = id.full_path_name*"_coordinates.bin"

            if switch_reduced_density && id.params.dimension == "2D"
                reduced_density_DOF1_output = id.full_path_name*"_reduced_density_DOF1.bin"
                reduced_density_DOF2_output = id.full_path_name*"_reduced_density_DOF2.bin"
            end
            if id.params.dimension == "1D"
                ϕ_matrix=Matrix{ComplexF64}(undef,length(r),length(ϵ))
            elseif id.params.dimension == "2D"
                ϕ_matrix=Matrix{ComplexF64}(undef,length(r[1])*length(r[2]),length(ϵ))
            end
            ϵ_vector=Vector{ComplexF64}(undef,length(ϵ))

            println("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
            println("Transforming data from FE object to complex value...")
            Threads.@threads for i in eachindex(ϵ)
                ϕ_matrix[:,i] = ϕ[i].(pts)
                ϵ_vector[i] = ϵ[i]
            end

            output=OutputData("bin_default",id.output_format_type,eigen_vectors_output,ϕ_matrix)
            remove_and_write_data(output)
            output=OutputData("bin_default",id.output_format_type,eigen_values_output,ϵ_vector)
            remove_and_write_data(output)

            if switch_reduced_density && id.params.dimension == "2D"
                output=OutputData("bin_default",id.output_format_type,reduced_density_DOF1_output,rho_DOF1_matrix)
                remove_and_write_data(output)
                output=OutputData("bin_default",id.output_format_type,reduced_density_DOF2_output,rho_DOF2_matrix)
                remove_and_write_data(output)
            end

            if id.params.dimension == "1D"
                output=OutputData("bin_default",id.output_format_type,coordinates_output,r)
                remove_and_write_data(output)
            elseif id.params.dimension == "2D"
                r_matrix=Matrix{Float64}(undef,max(length(r[1]),length(r[2])),2)
                Threads.@threads for i in [1,2]
                    r_matrix[1:length(r[i]),i] = r[i]
                end
                output=OutputData("bin_default",id.output_format_type,coordinates_output,r_matrix)
                remove_and_write_data(output)
            end
        elseif (id.output_format_type)==("jld2","eigen")
            if switch_reduced_density && id.params.dimension == "2D"
                output_file_name = id.full_path_name*"_eigen_data.jld2"
                output=OutputData("jld2_reduced_density",id.output_format_type,output_file_name,
                    (ϵ,ϕ,r,pts,rho_DOF1_matrix,rho_DOF2_matrix))
                remove_and_write_data(output)
            else
                output_file_name = id.full_path_name*"_eigen_data.jld2"
                output=OutputData("jld2_default",id.output_format_type,output_file_name,(ϵ,ϕ,r,pts))
                remove_and_write_data(output)
            end
        elseif (id.output_format_type)==("jld2","all")
            if switch_reduced_density && id.params.dimension == "2D"
                output_file_name = id.full_path_name*"_eigen_data.jld2"
                output=OutputData("jld2_reduced_density",id.output_format_type,output_file_name,
                    (ϵ,ϕ,r,pts,model,rho_DOF1_matrix,rho_DOF2_matrix))
                remove_and_write_data(output)
            else
                output_file_name = id.full_path_name*"_eigen_data.jld2"
                output=OutputData("jld2_default",id.output_format_type,output_file_name,(ϵ,ϕ,r,pts,model))
                remove_and_write_data(output)
            end
        end
    elseif type_potential in ["1","2","3","4"]
        println("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
        println("... Running ...")
        if params.dimension == "2D"
            parameters = (params,different_masses)
            if params.dom_type=="s"
                dom=(-0.5*params.Lx,0.5*params.Lx,-0.5*params.Ly,0.5*params.Ly)
            elseif params.dom_type=="ns"
                dom=(0.0,params.Lx,0.0,params.Lx)
            end
            r,pts=space_coord(dom,(params.Lx/params.nx,params.Ly/params.ny),(params.nx,params.ny))
        elseif params.dimension == "1D"
            parameters = params
            if params.dom_type=="s"
                dom=(-0.5*params.L,0.5*params.L)
            elseif params.dom_type=="ns"
                dom=(0.0,params.L)
            end
            r,pts=space_coord(dom,params.Δx,round(Int,abs(params.L/params.Δx));dimension="1D")
        end

        ϵ,ϕ=dispatch_solver(params.dimension,"default",parameters)

        if (output_format_type)== ("bin","eigen")
            eigen_vectors_output = full_path_name*"_eigen_vectors.bin"
            eigen_values_output = full_path_name*"_eigen_values.bin"
            coordinates_output = full_path_name*"_coordinates.bin"
            if params.dimension == "1D"
                ϕ_matrix=Matrix{ComplexF64}(undef,length(r),length(ϵ))
            elseif params.dimension == "2D"
                ϕ_matrix=Matrix{ComplexF64}(undef,length(r[1])*length(r[2]),length(ϵ))
            end
            ϵ_vector=Vector{ComplexF64}(undef,length(ϵ))
            println("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
            println("Transforming data from FE object to complex value...")
            Threads.@threads for i in eachindex(ϵ)
                ϕ_matrix[:,i] = ϕ[i].(pts)
                ϵ_vector[i] = ϵ[i]
            end
            output=OutputData("bin_default",output_format_type,eigen_vectors_output,ϕ_matrix)
            remove_and_write_data(output)
            output=OutputData("bin_default",output_format_type,eigen_values_output,ϵ_vector)
            remove_and_write_data(output)
            if dimension == "1D"
                output=OutputData("bin_default",output_format_type,coordinates_output,r)
                remove_and_write_data(output)
            elseif dimension == "2D"
                r_matrix=Matrix{Float64}(undef,max(length(r[1]),length(r[2])),2)
                for i in [1,2]
                    r_matrix[1:length(r[i]),i] = r[i]
                end
                output=OutputData("bin_default",output_format_type,coordinates_output,r_matrix)
                remove_and_write_data(output)
            end
        elseif (output_format_type)==("jld2","eigen")
            output_file_name = full_path_name*"_eigen_data.jld2"
            output=OutputData("jld2_default",output_format_type,output_file_name,(ϵ,ϕ,r,pts))
            remove_and_write_data(output)
        elseif (output_format_type)==("jld2","all")
            output_file_name = full_path_name*"_eigen_data.jld2"
            output=OutputData("jld2_default",output_format_type,output_file_name,(ϵ,ϕ,r,pts,model))
            remove_and_write_data(output)
        end
    end
end

function set_include()
    println("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    print("Set full path name (e.g: \"./my_directory/my_input_data\") where the data is specified and press Enter = ")
    full_path_input_data_name = readline()
    id=input_data(full_path_input_data_name)
    include((id.adhoc_file_name)*".jl")
    return id
end

function set_type_potential()
    println("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    println("####### ####### #     # #######   ###    #####  #######")
    println("#       #       ##   ##    #       #    #     # #                     #  #")
    println("#       #       # # # #    #       #    #       #                     #  #")
    println("#####   #####   #  #  #    #       #     #####  #####                 #  #")
    println("#       #       #     #    #       #          # #         ###         #  #")
    println("#       #       #     #    #       #    #     # #         ###    #    #  #")
    println("#       ####### #     #    #      ###    #####  #######   ###     ####   ######")
    println("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    println("The type of default potential can be:")
    println("   - Unidimensional Quantum Harmonic Oscillator            --> set (1)")
    println("   - Unidimensional Symmetric Finite Kronig-Penney         --> set (2)")
    println("   - Unidimensional Finite Well Potential                  --> set (3)")
    println("   - Bidimensional Isotropic Quantum Harmonic Oscillator   --> set (4)")
    println("   - Ad hoc potential from input file                      --> set (5)")
    print("Please, set some number to specify the type potential: ")
    type_potential = readline()
    type_potential == "5" ? id=set_include() : id=nothing
    return tuple(type_potential,id)
end

function set_type_potential(full_path_input_data_name::String)
    id=input_data(full_path_input_data_name)
    include((id.adhoc_file_name)*".jl")
    return tuple("5",id)
end