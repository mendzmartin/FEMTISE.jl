# 2D quantum harmonic oscillator potential
function eigenvalue_problem_functions(params::Tuple)
    ħ::Float64=1.0; m::Float64=1.0;
    ω,x₁,y₁=params;
    p(x) = 0.5*(ħ*ħ)*(1.0/m);                                       # kinetic energy
    q(x) = 0.5*m*(ω*ω)*((x[1]-x₁)*(x[1]-x₁)+(x[2]-y₁)*(x[2]-y₁));   # 2D harmonic osicllator centered in (x₁,y₁)
    r(x) = 1.0;
    return p,q,r
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