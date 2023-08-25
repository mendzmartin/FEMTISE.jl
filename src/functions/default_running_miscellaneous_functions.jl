function sturm_liouville_formulation(params::Tuple)
    # params=(potential_function_name::String,params_potential::Tuple)
    potential_function_name::Function=eval(Symbol(params[1]))
    ħ=1.0;m=1.0;
    p(x) = 0.5*(ħ*ħ)*(1.0/m);
    q(x) = potential_function_name(x,params[2]);
    r(x) = 1.0;
    return p,q,r;
end

function sturm_liouville_formulation(params::Tuple,different_masses::Tuple)
    # params=(potential_function_name::String,params_potential::Tuple)
    potential_function_name::Function=eval(Symbol(params[1]))
    ħ=1.0;m=1.0;
    if different_masses[1]
        p₁(x) = 0.5*(ħ*ħ)*(1.0/m)
        p₂(x) = p₁(x)*different_masses[2]
        p=tuple(p₁,p₂)
    else
        p₃(x) = 0.5*(ħ*ħ)*(1.0/m)
        p=tuple(p₃)
    end
    q(x) = potential_function_name(x,params[2]);
    r(x) = 1.0;
    return p,q,r;
end

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# auxiliary functions to define potentials
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function heaviside(x)
    return 0.5*(sign(x)+1)==true
 end

"""
    sym_rect_pot_barr(x,b,V₀)

    to compute symetrical rectangular potential barrer
"""
function sym_rect_pot_barr(x,b::Real,V₀::Real)
   return V₀*(heaviside(x+0.5*b)-heaviside(x-0.5*b))
end

function kronig_penney_center(x,b::Real,V₀::Real)
    return sym_rect_pot_barr.(x,b,V₀)
end

function kronig_penney_left(x,num_ions::Integer,a::Real,b::Real,V₀::Real)
    result=0.0
    for i in 1:num_ions
        result = result .+ sym_rect_pot_barr.(x.+i*a,b,V₀)
    end
    return result
end

function kronig_penney_right(x,num_ions::Integer,a::Real,b::Real,V₀::Real)
    return kronig_penney_left(-x,num_ions,a,b,V₀)
end

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Potential functions
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function finite_well_1d(x,params::Tuple)
    V₀::Real,b::Real=params
    return kronig_penney_center(x[1],b,V₀)
end

# symmetric Kronig-Peney potential 1D
function kronig_penney_1d(x,params::Tuple)
    V₀::Real,t::Real,b::Real,num_ions::Integer=params
    a = t+b
    if (mod(num_ions,2) == 0)
        error("num_ions keyword need to be odd")
        stop()
    end
    return kronig_penney_center(x[1],b,V₀) .+ kronig_penney_left(x[1],convert(Int,(num_ions-1)/2),a,b,V₀) .+ kronig_penney_right(x[1],convert(Int,(num_ions-1)/2),a,b,V₀)
end

function qho_1d(x,params::Tuple)
    ω,x₁=params
    return 0.5*(ω*ω)*((x[1]-x₁)*(x[1]-x₁))
end

function qho_2d(x,params::Tuple)
    ω,x₁,y₁=params
    return 0.5*(ω*ω)*((x[1]-x₁)*(x[1]-x₁)+(x[2]-y₁)*(x[2]-y₁))
end

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# function to resolve eigenproblem
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function default_solver_eigen_problem(params::Params1D)
    grid_type="simple_line";
    if params.dom_type=="s"
        dom=(-0.5*params.L,0.5*params.L)
    elseif params.dom_type=="ns"
        dom=(0.0,params.L)
    end
    params_model=("./","model1D",dom,params.Δx);  

    println("Building the grid model ...")

    model=make_model(grid_type,params_model);
    rm(params_model[1]*params_model[2]*".msh")
    FullDirichlet_values,FullDirichlet_tags=make_boundary_conditions(grid_type,"FullDirichlet",ComplexF64);
    # Ω,dΩ,Γ,dΓ=measures(model,3,FullDirichlet_tags);
    dΩ=measures(model,3,FullDirichlet_tags)[2]
    reff = ReferenceFE(lagrangian,Float64,2);
    # VSpace,USpace=fe_spaces(model,reff,grid_type;BC_type="FullDirichlet",TypeData=ComplexF64,conf_type=:H1);
    VSpace,USpace=fe_spaces(model,reff;BC_data=(FullDirichlet_values,FullDirichlet_tags),BC_type="Dirichlet")
    params_sturm_liouville=(params.potential_function_name,params.params_potential)
    p,q,r = sturm_liouville_formulation(params_sturm_liouville)

    println("Solving eigen problem ...")
    ϵ,ϕ = eigen_values_and_eigen_vectors(p,q,r,dΩ,USpace,VSpace;params=(params.nev,10e-9,500,:none,params.sigma));
    
    return tuple(ϵ,ϕ);
end

function default_solver_eigen_problem(params::Params2D,different_masses::Tuple)
    if params.dom_type=="s"
        dom=(-0.5*params.L,0.5*params.L,-0.5*params.L,0.5*params.L)
    elseif params.dom_type=="ns"
        dom=(0.0,params.L,0.0,params.L)
    end
    nxy=(params.nx,params.ny);  # define number of FE to use in each coordinate
    params_model=(dom,nxy);     # organize data in a tuple
    grid_type="Cartesian2D";    # define type of FE grid

    println("Building the grid model ...")

    model=make_model(grid_type,params_model);
    FullDirichlet_values,FullDirichlet_tags=make_boundary_conditions(grid_type,"FullDirichlet",ComplexF64);
    # Ω,dΩ,Γ,dΓ=measures(model,3,FullDirichlet_tags);
    dΩ=measures(model,3,FullDirichlet_tags)[2]
    reff = ReferenceFE(lagrangian,Float64,2);
    # VSpace,USpace=fe_spaces(model,reff,grid_type;BC_type="FullDirichlet",TypeData=ComplexF64,conf_type=:H1);
    VSpace,USpace=fe_spaces(model,reff;BC_data=(FullDirichlet_values,FullDirichlet_tags),BC_type="Dirichlet")
    params_sturm_liouville=(params.potential_function_name,params.params_potential)
    p,q,r = sturm_liouville_formulation(params_sturm_liouville,different_masses)

    println("Solving eigen problem ...")

    if different_masses[1]
        ϵ,ϕ = eigen_values_and_eigen_vectors(p[1],p[2],q,r,dΩ,USpace,VSpace;
            params=(params.nev,10e-9,500,:none,params.sigma));
    else
        ϵ,ϕ = eigen_values_and_eigen_vectors(p[1],q,r,dΩ,USpace,VSpace;
            params=(params.nev,10e-9,500,:none,params.sigma));
    end

    return tuple(ϵ,ϕ);
end

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# function to ask about params of adhoc potential function
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function ask_for_params(length_params::Int64)
    params_vector=Any[]
    for i in 1:length_params
        print("Set parameter Nº $(i)\n")
        print("Specify type of parameter Nº $(i): c (if ComplexF64), f (if Float64), i (if Int64) = ")
        type_λ_str=readline()
        print("Set parameter Nº $(i) λ = ")
        λ_str = readline()
        if type_λ_str=="c"
            λ=parse(ComplexF64,λ_str)
        elseif type_λ_str=="f"
            λ=parse(Float64,λ_str)
        elseif type_λ_str=="i"
            λ=parse(Int64,λ_str)
        end
        push!(params_vector,λ)
    end
    return Tuple(params_vector)
end

function get_input(T::Type;default_data::Bool=true,default_value=nothing)
    while true
        input = readline()
        if  default_data==true && input=="" && default_value≠nothing
            return default_value
        end
        try
            return parse(T, input)
        catch e
            if isa(e, ArgumentError)
                println("Invalid input. Please try again.")
            else
                rethrow()
            end
        end
    end
end

function get_input(possible_data::Array)
    input = readline()
    input in possible_data ? again=false : again=true
    while again
        println("Invalid input. Please try again.")
        input = readline()
        input in possible_data ? again=false : again=true
    end
    return input
end

function get_input(T::Type,BoolCondition::Function)
    input = get_input(T;default_data=false)
    BoolCondition(input;odd=true) ? again=true : again=false
    while again
        println("Invalid input. Please try again.")
        input = get_input(T;default_data=false)
        BoolCondition(input;odd=true) ? again=true : again=false
    end
    return input
end

function num_ions_BoolCondition(num_ions::Int64;odd::Bool=true)
    odd ? bool_data=(mod(num_ions,2)==0) : bool_data=(mod(num_ions,2)==1)
    return bool_data
end

function rm_existing_file(full_path_file::String)
    isfile(full_path_file) ? rm(full_path_file) : nothing
end