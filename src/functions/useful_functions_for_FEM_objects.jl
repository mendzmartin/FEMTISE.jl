# +++++++++++++++++++++++++++++++++++++++++++++++
# General maths calculations
# +++++++++++++++++++++++++++++++++++++++++++++++
function uniform_trapezoidal_integration_method(x_vec::Vector,fx_vec::Vector)
    dim_x=length(x_vec);
    coef_vec=ones(dim_x);
    coef_vec[2:(dim_x-1)].=2.0;
    @views function_vec=copy(fx_vec);
    Δx=abs(x_vec[2]-x_vec[1]); # válido para cuando Δx es constante
    return 0.5*Δx*dot(coef_vec,function_vec);
end

function aprox_dirac_delta(x,params::Tuple{Float64,Float64,Int,Float64})
    x₀,δnorm,component,Δx=params
    (abs(x[component]-x₀)≤(0.5*Δx)) ? δ=(1.0/Δx)*(1.0/δnorm) : δ=0.0
    return δ
end

function reduced_integration(FE_function::Vector{CellField},r_vector::Tuple{Vector{Float64},Vector{Float64}},
    Ω,dΩ::Gridap.CellData.GenericMeasure)
    N_DOF1=abs(r_vector[1][end]-r_vector[1][1]);
    N_DOF2=abs(r_vector[2][end]-r_vector[2][1]);

    Δr_DOF1=abs(r_vector[1][2]-r_vector[1][1]);
    Δr_DOF2=abs(r_vector[2][2]-r_vector[2][1]);

    reduced_function_DOF1=zeros(Float64,length(r_vector[1]),length(FE_function))
    reduced_function_DOF2=zeros(Float64,length(r_vector[2]),length(FE_function))

    Threads.@threads for i in eachindex(r_vector[1])
        params=(r_vector[1][i],1.0,1,Δr_DOF1)
        gridap_dirac_delta=CellField(x->aprox_dirac_delta(x,params),Ω);

        δnorm=sum(integrate(gridap_dirac_delta,dΩ));

        params=(r_vector[1][i],δnorm/N_DOF2,1,Δr_DOF1);
        gridap_dirac_delta=CellField(x->aprox_dirac_delta(x,params),Ω);

        for j in eachindex(FE_function)
            reduced_function_DOF1[i,j]=sum(integrate(FE_function[j]*gridap_dirac_delta,dΩ))
        end
    end

    Threads.@threads for i in eachindex(r_vector[2])
        params=(r_vector[2][i],1.0,2,Δr_DOF2)
        gridap_dirac_delta=CellField(x->aprox_dirac_delta(x,params),Ω);

        δnorm=sum(integrate(gridap_dirac_delta,dΩ));

        params=(r_vector[2][i],δnorm/N_DOF1,2,Δr_DOF2);
        gridap_dirac_delta=CellField(x->aprox_dirac_delta(x,params),Ω);

        for j in eachindex(FE_function)
            reduced_function_DOF2[i,j]=sum(integrate(FE_function[j]*gridap_dirac_delta,dΩ))
        end
    end

    return reduced_function_DOF1,reduced_function_DOF2
end

# +++++++++++++++++++++++++++++++++++++++++++++++
# Testing mandatory quantum properties
# +++++++++++++++++++++++++++++++++++++++++++++++
"""
    norm_l2(𝜳,dΩ)

# Aim
- Compute de L2 norm for specific FE wave function

# Arguments
- `𝜳::CellField`: specific FE wave function
- `dΩ::Gridap.CellData.GenericMeasure`: integration domain
"""
function norm_l2(𝜳::CellField,dΩ::Gridap.CellData.GenericMeasure)
    return sqrt(real(sum(∫(𝜳'*𝜳)*dΩ)));
end

function orthogonality_check(ϕ::Vector{CellField},dΩ::Gridap.CellData.GenericMeasure)
    nev::Integer=length(ϕ);
    orthogonality_vector=zeros(Float64,round(Int,(nev^2-nev)/2));
    index::Int=1;
    for i in 2:nev
        for j in 1:(i-1)
            orthogonality_vector[index]=abs(sum(∫(ϕ[j]'*ϕ[i])*dΩ));
            index+=1;
        end
    end
    return orthogonality_vector;
end

# multifield option
function orthogonality_check(ϕ::Vector{CellField},dΩ::Gridap.CellData.GenericMeasure,TrialSpace::FESpace)
    nev::Integer=length(ϕ);
    orthogonality_vector=zeros(Float64,round(Int,(nev^2-nev)/2));
    index::Int=1;
    for i in 2:nev
        ϕᵢ=interpolate_everywhere(ϕ[i],TrialSpace);
        ϕ¹ᵢ,ϕ²ᵢ=ϕᵢ;
        for j in 1:(i-1)
            ϕⱼ=interpolate_everywhere(ϕ[j],TrialSpace);
            ϕ¹ⱼ,ϕ²ⱼ=ϕⱼ;
            orthogonality_vector[index]=abs(sum(∫(ϕ¹ⱼ'*ϕ¹ᵢ)*dΩ)+sum(∫(ϕ²ⱼ'*ϕ²ᵢ)*dΩ));
            index+=1;
        end
    end
    return orthogonality_vector;
end

function eigenstates_normalization(ϕ::Vector{CellField},dΩ::Gridap.CellData.GenericMeasure)
    nom_vector=zeros(Float64,length(ϕ));
    Threads.@threads for i in eachindex(ϕ)
        nom_vector[i]=norm_l2(ϕ[i],dΩ);
    end
    return nom_vector;
end

# multifield option
function eigenstates_normalization(ϕ::Vector{CellField},dΩ::Gridap.CellData.GenericMeasure,
    TrialSpace::FESpace)
    nom_vector=zeros(Float64,length(ϕ));
    Threads.@threads for i in eachindex(ϕ)
        ϕᵢ=interpolate_everywhere(ϕ[i],TrialSpace);
        ϕ¹ᵢ,ϕ²ᵢ=ϕᵢ;
        nom_vector[i]=norm_l2(ϕ¹ᵢ,dΩ)+norm_l2(ϕ²ᵢ,dΩ);
    end
    return nom_vector;
end


# +++++++++++++++++++++++++++++++++++++++++++++++
# Quantum probability densities
# +++++++++++++++++++++++++++++++++++++++++++++++

function density(ϕ::Vector{CellField})
    rho=Vector{CellField}(undef,length(ϕ))
    Threads.@threads for i in eachindex(ϕ)
        rho[i] = real(conj(ϕ[i])*ϕ[i])
    end
    return rho
end

function reduced_density(ϕ::Vector{CellField},r::Tuple{Vector{Float64},Vector{Float64}},model::CartesianDiscreteModel)   
    grid_type="Cartesian2D";
    FullDirichlet_tags=make_boundary_conditions(grid_type,"FullDirichlet",ComplexF64)[2];
    Ω,dΩ=measures(model,3,FullDirichlet_tags)[1:2];
    rho = density(ϕ)
    return reduced_integration(rho,r,Ω,dΩ)
end

# +++++++++++++++++++++++++++++++++++++++++++++++
# Quantum Information
# +++++++++++++++++++++++++++++++++++++++++++++++

# Shannon case
"""
    time_indep_entropy(psi,TrialSpace,dOmega)

# Aim
- Compute time independent Shannon entropy for specific FE wave function

# Arguments
- `psi::Vector{CellField}`: specific FE wave function
- `TrialSpace::FESpace`: Trial FE space
- `dOmega::Gridap.CellData.GenericMeasure`: integration domain
"""
function time_indep_entropy(psi::Vector{CellField},TrialSpace::FESpace,dOmega::Gridap.CellData.GenericMeasure)
    S=zeros(Float64,length(psi))
    Threads.@threads for i in eachindex(psi)
        psi_index=interpolate_everywhere(psi[i],TrialSpace);
        psi_index=psi_index*(1.0/norm_l2(psi_index,dOmega));
        rho_index=real(conj(psi_index)*psi_index)
        (sum(integrate(rho_index,dOmega))==0.0) ? (S[i]=0.0) : (S[i]=-sum(integrate(rho_index*(log∘rho_index),dOmega)))
    end
    return S;
end

# Renyi case
"""
    time_indep_entropy(psi,TrialSpace,dOmega,renyi_factor)

# Aim
- Compute time independent Renyi entropy for specific FE wave function and Renyi factor

# Arguments
- `psi::Vector{CellField}`: specific FE wave function
- `TrialSpace::FESpace`: Trial FE space
- `dOmega::Gridap.CellData.GenericMeasure`: integration domain
- `renyi_factor::Real`: Renyi factor (need to be a Real number)
"""
function time_indep_entropy(psi::Vector{CellField},TrialSpace::FESpace,dOmega::Gridap.CellData.GenericMeasure,renyi_factor::Real)
    S=zeros(Float64,length(psi))
    Threads.@threads for i in eachindex(psi)
        psi_index=interpolate_everywhere(psi[i],TrialSpace);
        psi_index=psi_index*(1.0/norm_l2(psi_index,dOmega));
        rho_index=real(conj(psi_index)*psi_index)
        (sum(integrate(rho_index,dOmega))==0.0) ? (S[i]=0.0) : (S[i]=log(sum(integrate(exp∘(renyi_factor*(log∘rho_index)),dOmega))))
    end
    return S .* (1.0/(1.0-renyi_factor));
end

# Shannon case
function integration_argument_diff_entropy(rho::Vector)
    factor=similar(rho);
    Threads.@threads for index in eachindex(rho)
        rho[index]==0.0 ? factor[index]=0.0 : factor[index]=rho[index]*log(rho[index])
    end
    return factor
end

# Renyi case
function integration_argument_diff_entropy(rho::Vector,renyi_factor::Real)
    factor=similar(rho);
    Threads.@threads for index in eachindex(rho)
        rho[index]==0.0 ? factor[index]=0.0 : factor[index]=rho[index]^renyi_factor
    end
    return factor
end

# Shannon case
function reduced_time_indep_entropy(x::Vector,rho::Matrix)
    Sx_vector=similar(rho[1,:]);
    Threads.@threads for i in eachindex(Sx_vector)
        Sx_vector[i]=uniform_trapezoidal_integration_method(x,integration_argument_diff_entropy(rho[:,i]));
    end
    return -1.0 .* Sx_vector;
end

# Renyi case
function reduced_time_indep_entropy(x::Vector,rho::Matrix,renyi_factor::Real)
    Sx_vector=similar(rho[1,:]);
    Threads.@threads for i in eachindex(Sx_vector)
        Sx_vector[i]=uniform_trapezoidal_integration_method(x,integration_argument_diff_entropy(rho[:,i],renyi_factor));
    end
    return Sx_vector .* (1.0/(1.0-renyi_factor));
end

# compute mutual information
function time_indep_diff_mutual_information(Sx::Vector,Sy::Vector,Sxy::Vector)
    return Sx .+ Sy .- Sxy;
end

# +++++++++++++++++++++++++++++++++++++++++++++++
# Quantum Optics
# +++++++++++++++++++++++++++++++++++++++++++++++

# multifield
# to compute the population of specific quantum state
function state_population(psi::Vector{CellField},TrialSpace::FESpace,dOmega::Gridap.CellData.GenericMeasure)
    prob_1=zeros(Float64,length(psi));
    prob_2=copy(prob_1);
    for i in eachindex(psi)
        psi_index=interpolate_everywhere(psi[i],TrialSpace);
        psi_index_1,psi_index_2=psi_index
        prob_1[i]=real(sum(integrate(conj(psi_index_1)*psi_index_1,dOmega)))*(1.0/norm_l2(psi_index_1,dOmega))
        prob_2[i]=real(sum(integrate(conj(psi_index_2)*psi_index_2,dOmega)))*(1.0/norm_l2(psi_index_2,dOmega))
    end
    return prob_1,prob_2;
end

# +++++++++++++++++++++++++++++++++++++++++++++++
# Quantum quantities
# +++++++++++++++++++++++++++++++++++++++++++++++

# to compute first moment (expectated value for coordinate or position)
function coord_first_moment(psi::Vector{CellField},TrialSpace::FESpace,Omega,dOmega::Gridap.CellData.GenericMeasure,x_component::Int)
    Gridap_factor=CellField(x->x[x_component],Omega);
    expval=zeros(Float64,length(psi));
    Threads.@threads for index in eachindex(psi)
        psi_index=interpolate_everywhere(𝛹ₓₜ[index],TrialSpace)
        # ojo! tomamos la parte real porque se trata de la coord. espacial, pero puede ser complejo
        expval[index]=real(sum(integrate(conj(psi_index)*Gridap_factor*psi_index,dOmega)))
    end
    return expval;
end

# to compute second moment (coordinate or position variance)
function coord_second_moment(psi::Vector{CellField},TrialSpace::FESpace,Omega,dOmega::Gridap.CellData.GenericMeasure,x_component::Int)
    Gridap_factor=CellField(x->x[x_component]*x[x_component],Omega);
    var=zeros(Float64,length(psi));
    Threads.@threads for index in eachindex(psi)
        psi_index=interpolate_everywhere(𝛹ₓₜ[index],TrialSpace)
        # ojo! tomamos la parte real porque se trata de la coord. espacial, pero puede ser complejo
        var[index]=real(sum(integrate(conj(psi_index)*Gridap_factor*psi_index,dOmega)))
    end
    return var;
end

# to standar deviation (coordinate or position standar deviation)
function coord_second_moment(expval::Vector{Float64},var::Vector{Float64})
    return sqrt.(var.-(expval.*expval))
end