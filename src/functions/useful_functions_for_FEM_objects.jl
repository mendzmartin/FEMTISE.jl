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