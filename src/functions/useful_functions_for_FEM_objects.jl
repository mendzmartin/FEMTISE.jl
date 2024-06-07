# +++++++++++++++++++++++++++++++++++++++++++++++
# General maths calculations
# +++++++++++++++++++++++++++++++++++++++++++++++
"""
    uniform_trapezoidal_integration_method(x_vec,fx_vec)

# Aim
- Numeric integration function using trapezoidal method considering uniform discretization.

# Arguments
- `x_vec::Vector`: integration domain (array of values)
- `fx_vec::Vector`: array function evaluated in integration domain

# Output
- `result::Vector`: integral value
"""
function uniform_trapezoidal_integration_method(x_vec::Vector,fx_vec::Vector)
    dim_x=length(x_vec);
    coef_vec=ones(dim_x);
    coef_vec[2:(dim_x-1)].=2.0;
    @views function_vec=copy(fx_vec);
    Δx=abs(x_vec[2]-x_vec[1]); # válido para cuando Δx es constante
    return 0.5*Δx*dot(coef_vec,function_vec);
end

"""
    aprox_dirac_delta(x,params)

# Aim
- Approximation of Dirac delta function as rectangular function

# Arguments
- `x`: independent variable (is a FE coordinate object)
- `params::Tuple{Float64,Float64,Int,Float64}`: parameters
  - `x₀::Float64`: specific point to centre dirac function
  - `δnorm::Float64`: norm value to obtain normalized delta function
  - `component::Int`: component of x coordinate
  - `Δx::Float64`: 	thickness of rectangular function

# Output
- `result::Vector`: integral value
"""
function aprox_dirac_delta(x,params::Tuple{Float64,Float64,Int,Float64})
    x₀,δnorm,component,Δx=params
    (abs(x[component]-x₀)≤(0.5*Δx)) ? δ=(1.0/Δx)*(1.0/δnorm) : δ=0.0
    return δ
end

"""
    reduced_integration(FE_function,r_vector,Ω,dΩ)

# Aim
- Integration over reduced coordinate (partial domain) of specific two dimensional FE function

# Arguments
- `FE_function::Vector{CellField}`: specific 2D FE function to integrate over partial domain.
- `r_vector::Tuple{Vector{Float64},Vector{Float64}}`
  - `x_vector::Vector{Float64}`:: array values of first coordinate (first partial domain)
  - `y_vector::Vector{Float64}`:: array values of second coordinate (second partial domain)
- `Omega`: Finite element domain
- `dΩ::Gridap.CellData.GenericMeasure`: full integration domain

# Output
- `result::Tuple{Vector{Float64},Vector{Float64}}`
  - `reduced_function_DOF1::Vector{Float64}`: integration of 2D FE function over first partial domain.
  - `reduced_function_DOF2::Vector{Float64}`: integration of 2D FE function over second partial domain.
"""
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

# Output
- `orthogonality_vector::Real`: L2 norm value
"""
function norm_l2(𝜳::CellField,dΩ::Gridap.CellData.GenericMeasure)
    return sqrt(real(sum(∫(𝜳'*𝜳)*dΩ)));
end

"""
orthogonality_check(ϕ,dΩ)

# Aim
- Check eigenstates orthogonality

# Arguments
- `ϕ::Vector{CellField}`: array of FE eigenstates
- `dΩ::Gridap.CellData.GenericMeasure`: integration domain

# Output
- `orthogonality_vector::Vector{Float64}`: array of inner product between differents eigenstates
"""
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

"""
orthogonality_check(ϕ,dΩ,TrialSpace)

# Aim
- Check eigenstates orthogonality. Each eigenstate is a multifield FE object.

# Arguments
- `ϕ::Vector{CellField}`: array of FE eigenstates
- `dΩ::Gridap.CellData.GenericMeasure`: integration domain
- `TrialSpace::FESpace`: Trial FE space

# Output
- `orthogonality_vector::Vector{Float64}`: array of inner product between differents eigenstates
"""
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

"""
eigenstates_normalization(ϕ,dΩ)

# Aim
- Check eigenstates normalization.

# Arguments
- `ϕ::Vector{CellField}`: array of FE eigenstates
- `dΩ::Gridap.CellData.GenericMeasure`: integration domain

# Output
- `nom_vector::Vector{Float64}`: array of norm value for specific eigenstate
"""
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
"""
    density(phi)

# Aim
- Compute probability density from set of specific wave functions (FE object)

# Arguments
- `phi::Vector{CellField}`: array of eigenstates (array of FE objects)

# Output
- `rho::Vector{CellField}`: Array of probability densities (array of FE objects)
"""
function density(phi::Vector{CellField})
    rho=Vector{CellField}(undef,length(phi))
    Threads.@threads for i in eachindex(phi)
        rho[i] = real(conj(phi[i])*phi[i])
    end
    return rho
end

"""
    reduced_density(ϕ,r,model)

# Aim
- Compute partial probability density from set of specific wave functions (FE object)

# Arguments
- `phi::Vector{CellField}`: array of eigenstates (array of FE objects)
- `r::Tuple{Vector{Float64},Vector{Float64}}`:
    - `x_vector::Vector{Float64}`:: array values of first coordinate (first partial domain)
    - `y_vector::Vector{Float64}`:: array values of second coordinate (second partial domain)
- `model::CartesianDiscreteModel`:: Cartesian discrete 2D model

# Output
- `result::Tuple{Vector{Float64},Vector{Float64}}`
  - `reduced_function_DOF1::Vector{Float64}`: partial probability density of first coordinate
  - `reduced_function_DOF2::Vector{Float64}`: partial probability density of second coordinate
"""
function reduced_density(phi::Vector{CellField},r::Tuple{Vector{Float64},Vector{Float64}},model::CartesianDiscreteModel)   
    grid_type="Cartesian2D";
    FullDirichlet_tags=make_boundary_conditions(grid_type,"FullDirichlet",ComplexF64)[2];
    Omega,dOmega=measures(model,3,FullDirichlet_tags)[1:2];
    rho = density(phi)
    return reduced_integration(rho,r,Omega,dOmega)
end

# +++++++++++++++++++++++++++++++++++++++++++++++
# Quantum Information
# +++++++++++++++++++++++++++++++++++++++++++++++

"""
    time_indep_entropy(psi,TrialSpace,dOmega)

# Aim
- Compute time independent Shannon entropy for specific FE wave function

# Arguments
- `psi::Vector{CellField}`: specific FE wave function
- `TrialSpace::FESpace`: Trial FE space
- `dOmega::Gridap.CellData.GenericMeasure`: integration domain

# Output
- `S::Vector{Float64}`: Array of entropy values.
"""
function time_indep_entropy(psi::Vector{CellField},TrialSpace::FESpace,dOmega::Gridap.CellData.GenericMeasure)
    S=zeros(Float64,length(psi))
    Threads.@threads for i in eachindex(psi)
        psi_index=interpolate_everywhere(psi[i],TrialSpace);
        psi_index=psi_index*(1.0/norm_l2(psi_index,dOmega));
        rho_index=real(conj(psi_index)*psi_index)
        (sum(integrate(rho_index,dOmega))==0.0) ? (S[i]=0.0) : (S[i]=-sum(integrate(rho_index*(log2∘rho_index),dOmega)))
    end
    return S;
end

"""
    time_indep_entropy(psi,TrialSpace,dOmega,renyi_factor)

# Aim
- Compute time independent Rényi entropy for specific FE wave function and Rényi factor

# Arguments
- `psi::Vector{CellField}`: specific FE wave function
- `TrialSpace::FESpace`: Trial FE space
- `dOmega::Gridap.CellData.GenericMeasure`: integration domain
- `renyi_factor::Real`: Rényi factor (have to be a Real number)

# Output
- `S::Vector{Float64}`: Array of entropy values.
"""
function time_indep_entropy(psi::Vector{CellField},TrialSpace::FESpace,dOmega::Gridap.CellData.GenericMeasure,renyi_factor::Real)
    S=zeros(Float64,length(psi))
    Threads.@threads for i in eachindex(psi)
        psi_index=interpolate_everywhere(psi[i],TrialSpace);
        psi_index=psi_index*(1.0/norm_l2(psi_index,dOmega));
        rho_index=real(conj(psi_index)*psi_index)
        (sum(integrate(rho_index,dOmega))==0.0) ? (S[i]=0.0) : (S[i]=log2(sum(integrate(exp∘(renyi_factor*(log2∘rho_index)),dOmega))))
    end
    return S .* (1.0/(1.0-renyi_factor));
end

# Shannon case
function integration_argument_diff_entropy(rho::Vector)
    factor=similar(rho);
    Threads.@threads for index in eachindex(rho)
        rho[index]==0.0 ? factor[index]=0.0 : factor[index]=rho[index]*log2(rho[index])
    end
    return factor
end

# Rényi case
function integration_argument_diff_entropy(rho::Vector,renyi_factor::Real)
    factor=similar(rho);
    Threads.@threads for index in eachindex(rho)
        rho[index]==0.0 ? factor[index]=0.0 : factor[index]=rho[index]^renyi_factor
    end
    return factor
end

"""
    reduced_time_indep_entropy(r,rho)

# Aim
- Compute time independent Reduced Shannon entropy for specific 2D FE wave function

# Arguments
- `r::Tuple{Vector{Float64},Vector{Float64}}`:
  - `x::Vector{Float64}`: array values of first coordinate
  - `ỳ::Vector{Float64}`: array values of second coordinate
- `rho::Matrix{Float64}`: matrix values of probability density (where rows are rho(x) with fixed y and columns are rho(y) with fixed x)

# Output
- `S::Tuple{Vector{Float64},Vector{Float64}}`:
  - `Sx_vector::Vector{Float64}`: reduced entropy as integration over first coordinate.
  - `Sy_vector::Vector{Float64}`: reduced entropy as integration over second coordinate.
"""
function reduced_time_indep_entropy(r::Tuple{Vector{Float64},Vector{Float64}},rho::Matrix{Float64})
    Sx_vector=similar(r[1]);
    Sy_vector=similar(r[2]);
    Threads.@threads for i in eachindex(Sx_vector)
        Sx_vector[i]=uniform_trapezoidal_integration_method(r[2],integration_argument_diff_entropy(rho[i,:]));
    end
    Threads.@threads for i in eachindex(Sy_vector)
        Sy_vector[i]=uniform_trapezoidal_integration_method(r[1],integration_argument_diff_entropy(rho[:,i]));
    end
    return tuple(-1.0 .* Sx_vector,-1.0 .* Sy_vector);
end

"""
    reduced_time_indep_entropy(r,rho,renyi_factor)

# Aim
- Compute time independent Reduced Rényi entropy for specific 2D FE wave function and Rényi factor

# Arguments
- `r::Tuple{Vector{Float64},Vector{Float64}}`:
  - `x::Vector{Float64}`: array values of first coordinate
  - `ỳ::Vector{Float64}`: array values of second coordinate
- `rho::Matrix{Float64}`: matrix values of probability density (where rows are rho(x) with fixed y and columns are rho(y) with fixed x)
- `renyi_factor::Float64`: Rényi factor (need to be a Real number)

# Output
- `S::Tuple{Vector{Float64},Vector{Float64}}`:
  - `Sx_vector::Vector{Float64}`: reduced entropy as integration over first coordinate.
  - `Sy_vector::Vector{Float64}`: reduced entropy as integration over second coordinate.
"""
function reduced_time_indep_entropy(r::Tuple{Vector{Float64},Vector{Float64}},rho::Matrix{Float64},renyi_factor::Real)
    Sx_vector=similar(r[1]);
    Sy_vector=similar(r[2]);
    Threads.@threads for i in eachindex(Sx_vector)
        Sx_vector[i]=uniform_trapezoidal_integration_method(r[2],integration_argument_diff_entropy(rho[i,:],renyi_factor));
    end
    Threads.@threads for i in eachindex(Sy_vector)
        Sy_vector[i]=uniform_trapezoidal_integration_method(r[1],integration_argument_diff_entropy(rho[:,i],renyi_factor));
    end
    return tuple(Sx_vector,Sy_vector) .* (1.0/(1.0-renyi_factor));
end

"""
    time_indep_diff_mutual_information(Sx,Sy,Sxy)

# Aim
- Compute time independent Differential Mutual Information

# Arguments
- `Sx::Vector{Float64}`: reduced entropy associated with first coordinate
- `Sy::Vector{Float64}`: reduced entropy associated with second coordinate
- `Sxy::Vector{Float64}`: total entropy associated with both coordinate

# Output
- `I::Vector{Float64}`: array of mutual information values
"""
function time_indep_diff_mutual_information(Sx::Vector{Float64},Sy::Vector{Float64},Sxy::Vector{Float64})
    return Sx .+ Sy .- Sxy;
end

# +++++++++++++++++++++++++++++++++++++++++++++++
# Quantum Optics
# +++++++++++++++++++++++++++++++++++++++++++++++

"""
    state_population(psi,TrialSpace,dOmega)

# Aim
- Compute the population of specific quantum state vector (2D)

# Arguments
- `psi::Vector{CellField}`: multifield quantum state (2D)
- `TrialSpace::FESpace: Trial finite element space
- `dOmega::Gridap.CellData.GenericMeasure`: integration domain

# Output
- `P::Tuple{Array{Float64},Array{Float64}}`:
  - `P1::Array{Float64}`: population of first state
  - `P2::Array{Float64}`: population of second state
"""
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

"""
    coord_first_moment(psi,TrialSpace,Omega,dOmega,x_component)

# Aim
- Compute first moment of coordinate (specific component)

# Arguments
- `psi::Vector{CellField}`: array of FE wave functions
- `TrialSpace::FESpace: Trial finite element space
- `Omega`: Finite element domain
- `dOmega::Gridap.CellData.GenericMeasure`: integration domain
- `x_component::Int`: component of specific coordinate variable

# Output
- `expval::Array{Float64}`: array of expectation values of specific coordinate
"""
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

"""
    coord_second_moment(psi,TrialSpace,Omega,dOmega,x_component)

# Aim
- Compute first moment of coordinate (specific component) or position variance

# Arguments
- `psi::Vector{CellField}`: array of FE wave functions
- `TrialSpace::FESpace: Trial finite element space
- `Omega`: Finite element domain
- `dOmega::Gridap.CellData.GenericMeasure`: integration domain
- `x_component::Int`: component of specific coordinate variable

# Output
- `var::Array{Float64}`: array of variance values of specific coordinate
"""
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


"""
    coord_standar_deviation(psi,TrialSpace,Omega,dOmega,x_component)

# Aim
- Compute standar deviation of coordinate (specific component)

# Arguments
- `psi::Vector{CellField}`: array of FE wave functions
- `TrialSpace::FESpace: Trial finite element space
- `Omega`: Finite element domain
- `dOmega::Gridap.CellData.GenericMeasure`: integration domain
- `x_component::Int`: component of specific coordinate variable

# Output
- `sigma::Array{Float64}`: array of standar deviation of specific coordinate
"""
function coord_standar_deviation(psi::Vector{CellField},TrialSpace::FESpace,Omega,dOmega::Gridap.CellData.GenericMeasure,x_component::Int)
    expval=coord_first_moment(psi,TrialSpace,Omega,dOmega,x_component)
    var=coord_second_moment(psi,TrialSpace,Omega,dOmega,x_component)
    return sqrt.(var.-(expval.*expval))
end

"""
    coord_standar_deviation(expval,var)

# Aim
- Compute standar deviation of coordinate (specific component)

# Arguments
- `expval::Array{Float64}`: array of expectation values of specific coordinate
- `var::Array{Float64}`: array of variance values of specific coordinate

# Output
- `sigma::Array{Float64}`: array of standar deviation of specific coordinate
"""
function coord_standar_deviation(expval::Vector{Float64},var::Vector{Float64})
    return sqrt.(var.-(expval.*expval))
end