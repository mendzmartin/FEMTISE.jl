export default_potential_sturm_liouville
function default_potential_sturm_liouville(type_potential::String,params::Tuple)
    if type_potential=="qho_1d"
        p,q,r=default_qho1d_sturm_liouville(params)
    elseif type_potential=="qho_2d"
        p,q,r=default_qho2d_sturm_liouville(params)
    elseif type_potential=="kronig_penney_1d"
        p,q,r=default_kronig_penney_sturm_liouville(params;fwp=false)
    elseif type_potential=="finite_well_1d"
        p,q,r=default_kronig_penney_sturm_liouville(params;fwp=true)
    end
    return p,q,r;
end

export default_solver_eigen_problem
function default_solver_eigen_problem(type_potential::String,params::Tuple)

    if type_potential=="qho_1d"
        L,Δx,ω,x₁,nev,sigma=params
        params_sturm_liouville=(ω,x₁)
        grid_type="simple_line";
        params_model=("./","model1D",(-0.5*L,0.5*L),Δx);
        dimension="1D"
    elseif type_potential=="qho_2d"
        L,nx,ny,ω,x₁,y₁,nev,sigma=params
        params_sturm_liouville=(ω,x₁,y₁)
        dom=(-0.5*L,0.5*L,-0.5*L,0.5*L);  # define domain 2D symetric concern origin
        nxy=(nx,ny);                        # define number of FE to use in each coordinate
        params_model=(dom,nxy);       # organize data in a tuple
        grid_type="Cartesian2D";        # define type of FE grid
        dimension="2D"
    elseif type_potential=="kronig_penney_1d"
        L,V₀,t,b,num_ions,Δx,nev,sigma=params
        params_sturm_liouville=(num_ions,b+t,b,V₀)
        grid_type="simple_line";
        params_model=("./","model1D",(-0.5*L,0.5*L),Δx);
        dimension="1D"
    elseif type_potential=="finite_well_1d"
        L,V₀,a,Δx,nev,sigma=params
        params_sturm_liouville=(a,V₀)
        grid_type="simple_line";
        params_model=("./","model1D",(-0.5*L,0.5*L),Δx);
        dimension="1D"
    end

    model=make_model(grid_type,params_model);
    dimension=="1D" ? rm(params_model[1]*params_model[2]*".msh") : nothing
    BC_type="FullDirichlet";
    FullDirichlet_values,FullDirichlet_tags=make_boundary_conditions(grid_type,BC_type,ComplexF64);
    Ω,dΩ,Γ,dΓ=measures(model,3,FullDirichlet_tags);
    reff = ReferenceFE(lagrangian,Float64,2);
    VSpace,USpace=fe_spaces(model,reff,grid_type;BC_type=BC_type,TypeData=ComplexF64,conf_type=:H1);
    p,q,r = default_potential_sturm_liouville(type_potential,params_sturm_liouville);
    ϵ,ϕ = eigen_values_and_eigen_vectors(p,q,r,dΩ,USpace,VSpace;params=(nev,10e-9,500,:none,sigma));

    @testset "Check normalization" begin
        for i in eachindex(ϕ)
            @test norm_l2(ϕ[1],dΩ) ≈ 1.0
        end
    end
    
    return ϵ,ϕ;
end