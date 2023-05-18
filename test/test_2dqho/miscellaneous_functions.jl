function eigenvalue_problem_functions(params;switch_potential = "QHO_1D")
    ħ::Float64=1.0; m::Float64=1.0;
    if (switch_potential == "QHO_1D")
        # 1D quantum harmonic oscillator potential
        println("Set quantum harmonic oscillator 1D potential");
        ω,x₁=params;
        p_QHO_1D(x) = 0.5*(ħ*ħ)*(1.0/m);                                      # kinetic energy
        q_QHO_1D(x) = 0.5*m*(ω*ω)*(x[1]-x₁)*(x[1]-x₁);                        # 1D harmonic oscillator centered in x₁
        r_QHO_1D(x) = 1.0;
        return p_QHO_1D,q_QHO_1D,r_QHO_1D
    elseif (switch_potential == "QHO_2D")
        # 2D quantum harmonic oscillator potential
        println("Set quantum harmonic oscillator 2D potential");
        ω,x₁,y₁=params;
        p_QHO_2D(x) = 0.5*(ħ*ħ)*(1.0/m);                                       # kinetic energy
        q_QHO_2D(x) = 0.5*m*(ω*ω)*((x[1]-x₁)*(x[1]-x₁)+(x[2]-y₁)*(x[2]-y₁));   # 2D harmonic osicllator centered in (x₁,y₁)
        r_QHO_2D(x) = 1.0;
        return p_QHO_2D,q_QHO_2D,r_QHO_2D
    end
end

function exactly_eigenvalues_2dqho(num_eigval::Integer)
    ϵ_real_aux=Array{Float64}(undef, num_eigval^2)
    index::Int64=1
    for i in 1:num_eigval
        for j in 1:num_eigval
            ϵ_real_aux[index]=((i-1)+(j-1)+1)
            index+=1
        end
    end
    ϵ_real_aux=sort(ϵ_real_aux)
    ϵ_real = ϵ_real_aux[1:num_eigval];
    return ϵ_real
end

function ortho_check_v2(ϕ::Vector{CellField},TrialSpace::FESpace,dΩ::Gridap.CellData.GenericMeasure)
    nev=length(ϕ)
    OrthoVector=zeros(Float64,round(Int,(nev^2-nev)/2));
    index::Int64=1
    for i in 2:nev
        ϕᵢ=interpolate_everywhere(ϕ[i],TrialSpace);
        for j in 1:(i-1)
            ϕⱼ=interpolate_everywhere(ϕ[j],TrialSpace);
            OrthoVector[index]=abs(sum(∫(ϕⱼ'*ϕᵢ)*dΩ))
            index+=1
        end
    end
    return OrthoVector;
end