"""
    eigenvalue_problem_functions(params)

# Aim
- Define functions for the 2D quantum harmonic oscillator potential

# Arguments
- `params::Tuple`: tuple with parameters of the 2D quantum harmonic oscillator potential
  - `params[1]::Float64`: frequency of the harmonic oscillator
  - `params[2]::Float64`: x coordinate of the harmonic oscillator
  - `params[3]::Float64`: y coordinate of the harmonic oscillator

# Example
```julia
p,q,r=eigenvalue_problem_functions((1.0,0.0,0.0))
```

# Returns
- `p::Function`: kinetic energy
- `q::Function`: potential energy
- `r::Function`: weight function
"""
function eigenvalue_problem_functions(params::Tuple)
    ħ::Float64=1.0; m::Float64=1.0;
    ω,x₁,y₁=params;
    p(x) = 0.5*(ħ*ħ)*(1.0/m);
    q(x) = 0.5*m*(ω*ω)*((x[1]-x₁)*(x[1]-x₁)+(x[2]-y₁)*(x[2]-y₁));
    r(x) = 1.0;
    return p,q,r
end

"""
    exactly_eigenvalues_2dqho(num_eigval)

# Aim
- Compute exactly the eigenvalues of the 2D quantum harmonic oscillator

# Arguments
- `num_eigval::Integer`: number of eigenvalues to compute

# Example
```julia
ϵ_real=exactly_eigenvalues_2dqho(5)
```
"""
function exactly_eigenvalues_2dqho(num_eigval::Integer)
    ϵ_real_aux=Array{Float64}(undef, num_eigval^2)
    index::Int=1
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