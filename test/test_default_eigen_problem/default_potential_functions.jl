export default_qho2d_sturm_liouville
"""
    default_qho2d_sturm_liouville(params)

# Aim
- Define default potential for 2D quantum harmonic oscillator in Sturm Liouville formulation

# Arguments
- `params::Tuple`: tuple with parameters of the 2D quantum harmonic oscillator potential
  - `params[1]::Float64`: frequency of the harmonic oscillator
  - `params[2]::Float64`: x coordinate of the harmonic oscillator
  - `params[3]::Float64`: y coordinate of the harmonic oscillator

# Example
```julia
p,q,r=default_qho2d_sturm_liouville((1.0,0.0,0.0))
```

# Returns
- `p::Function`: kinetic energy
- `q::Function`: potential energy
- `r::Function`: weight function
"""
function default_qho2d_sturm_liouville(params::Tuple)
    ħ::Float64=1.0; m::Float64=1.0;
    ω,x₁,y₁=params;
    p(x) = 0.5*(ħ*ħ)*(1.0/m);                                       # kinetic energy
    q(x) = 0.5*m*(ω*ω)*((x[1]-x₁)*(x[1]-x₁)+(x[2]-y₁)*(x[2]-y₁));   # 2D harmonic osicllator centered in (x₁,y₁)
    r(x) = 1.0;
    return p,q,r
end

export default_qho1d_sturm_liouville
"""
    default_qho1d_sturm_liouville(params)

# Aim
- Define default potential for 1D quantum harmonic oscillator in Sturm Liouville formulation

# Arguments
- `params::Tuple`: tuple with parameters of the 1D quantum harmonic oscillator potential
  - `params[1]::Float64`: frequency of the harmonic oscillator
  - `params[2]::Float64`: x coordinate of the harmonic oscillator

# Example
```julia
p,q,r=default_qho1d_sturm_liouville((1.0,0.0))
```

# Returns
- `p::Function`: kinetic energy
- `q::Function`: potential energy
- `r::Function`: weight function
"""
function default_qho1d_sturm_liouville(params::Tuple)
    ħ::Float64=1.0; m::Float64=1.0;
    ω,x₁=params;
    p(x) = 0.5*(ħ*ħ)*(1.0/m);                   # kinetic energy
    q(x) = 0.5*m*(ω*ω)*((x[1]-x₁)*(x[1]-x₁));   # 1D harmonic osicllator centered in (x₁,y₁)
    r(x) = 1.0;
    return p,q,r
end

export default_kronig_penney_sturm_liouville
"""
    default_kronig_penney_sturm_liouville(params; <keyword arguments>)

# Aim
- Define default potential for Kronig-Penney potential in Sturm Liouville formulation

# Arguments
- `params::Tuple`: tuple with parameters of the Kronig-Penney potential
  - `params[1]::Integer`: number of ions
  - `params[2]::Real`: lattice parameter
  - `params[3]::Real`: width of the potential barrer
  - `params[4]::Real`: height of the potential barrer
- <keyword arguments>
  - `fwp::Bool=false`: finite well potential as specific case of Kronig-Penney potential

# Example
```julia
p,q,r=default_kronig_penney_sturm_liouville((4,1.0,0.1,1.0))
```

# Returns
- `p::Function`: kinetic energy
- `q::Function`: potential energy
- `r::Function`: weight function
"""
function default_kronig_penney_sturm_liouville(params::Tuple;fwp::Bool=false)
    if fwp # finite well potential as specific case of Kronig-Penney potential
        b,V₀=params
        num_ions=1; # need to be only one
        a=1.0;      # could be any number
    else
        num_ions,a,b,V₀=params
    end
    ħ=1.0;m=1.0;
    p(x) = 0.5*(ħ*ħ)*(1.0/m);
    q(x) = symetric_kronig_penney(x[1],num_ions,a,b,V₀)
    r(x) = 1.0;
    return p,q,r;
end

# auxiliary functions to define potentials

"""
    heaviside(x)

# Aim
    - Heaviside step function

# Arguments
    - `x::Real`: real number

# Example
```julia
heaviside(-1.0)
```
"""
function heaviside(x)
    return 0.5*(sign(x)+1)==true
 end

"""
    sym_rect_pot_barr(x,b,V₀)

# Aim
    - Symetric rectangular potential barrer

# Arguments
    - `x::Real`: real number
    - `b::Real`: width of the potential barrer
    - `V₀::Real`: height of the potential barrer

# Example
```julia
sym_rect_pot_barr(0.0,1.0,1.0)
```
"""
function sym_rect_pot_barr(x,b::Real,V₀::Real)
   return V₀*(heaviside(x+0.5*b)-heaviside(x-0.5*b))
end

"""
    kronig_penney_center(x,b,V₀)

# Aim
    - Kronig-Penney potential centered in the origin

# Arguments
    - `x::Real`: real number
    - `b::Real`: width of the potential barrer
    - `V₀::Real`: height of the potential barrer
"""
function kronig_penney_center(x,b::Real,V₀::Real)
    return sym_rect_pot_barr.(x,b,V₀)
end

"""
    kronig_penney_left(x,num_ions,a,b,V₀)

# Aim
    - Kronig-Penney potential to the left of the origin

# Arguments
    - `x::Real`: real number
    - `num_ions::Integer`: number of ions
    - `a::Real`: lattice parameter
    - `b::Real`: width of the potential barrer
    - `V₀::Real`: height of the potential barrer
"""
function kronig_penney_left(x,num_ions::Integer,a::Real,b::Real,V₀::Real)
    result=0.0
    for i in 1:num_ions
        result = result .+ sym_rect_pot_barr.(x.+i*a,b,V₀)
    end
    return result
end

"""
    kronig_penney_right(x,num_ions,a,b,V₀)

# Aim
    - Kronig-Penney potential to the right of the origin

# Arguments
    - `x::Real`: real number
    - `num_ions::Integer`: number of ions
    - `a::Real`: lattice parameter
    - `b::Real`: width of the potential barrer
    - `V₀::Real`: height of the potential barrer
"""
function kronig_penney_right(x,num_ions::Integer,a::Real,b::Real,V₀::Real)
    return kronig_penney_left(-x,num_ions,a,b,V₀)
end

"""
    symetric_kronig_penney(x,num_ions,a,b,V₀)

# Aim
    - Symetric Kronig-Penney potential

# Arguments
    - `x::Real`: real number
    - `num_ions::Integer`: number of ions
    - `a::Real`: lattice parameter
    - `b::Real`: width of the potential barrer
    - `V₀::Real`: height of the potential barrer

# Example
```julia
symetric_kronig_penney(0.0,4,1.0,0.1,1.0)
```
"""
function symetric_kronig_penney(x,num_ions::Integer,a::Real,b::Real,V₀::Real)
    if (mod(num_ions,2) == 0)
        error("num_ions keyword need to be odd")
        stop()
    elseif (num_ions==1)
        kp = kronig_penney_center(x,b,V₀)
    else
        kp = kronig_penney_center(x,b,V₀) .+ kronig_penney_left(x,convert(Int,(num_ions-1)/2),a,b,V₀) .+ kronig_penney_right(x,convert(Int,(num_ions-1)/2),a,b,V₀)
    end
    return kp
end