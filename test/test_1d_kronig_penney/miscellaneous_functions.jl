"""
    kronig_penney_sturm_liouville(params;<keyword arguments>)

# Aim
    - Compute Kronig-Penney potential as Sturm-Liouville problem

# Arguments
    - `params::Tuple`: tuple with parameters
    - `<keyword arguments>`:
        - `fwp::Bool=false`: finite well potential as specific case of Kronig-Penney potential
"""
function kronig_penney_sturm_liouville(params::Tuple;fwp::Bool=false)
    if fwp # finite well potential as specific case of Kronig-Penney potential
        a,V₀=params
        num_ions=1; # need to be only one
        b=1.0;      # could be any number
    else
        num_ions,a,b,V₀=params
    end
    ħ=1.0;m=1.0;
    p(x) = 0.5*(ħ*ħ)*(1.0/m);
    q(x) = symetric_kronig_penney(x[1],num_ions,a,b,V₀)
    r(x) = 1.0;
    return p,q,r;
end

"""
    heaviside(x)

# Aim
    - Compute Heaviside function

# Arguments
    - `x::Real`: input value

# Example
```julia
heaviside(0.0)
```
"""
function heaviside(x)
    return 0.5*(sign(x)+1)==true
 end

"""
    sym_rect_pot_barr(x,b,V₀)

# Aim
    - Compute symmetric rectangular potential barrier

# Arguments
    - `x::Real`: input value
    - `b::Real`: width of barrier
    - `V₀::Real`: height of barrier

# Example
```julia
sym_rect_pot_barr(0.0,1.0,10.0)
```
"""
function sym_rect_pot_barr(x,b::Real,V₀::Real)
   return V₀*(heaviside(x+0.5*b)-heaviside(x-0.5*b))
end

"""
    kronig_penney_center(x,b,V₀)

# Aim
    - Compute Kronig-Penney potential center

# Arguments
    - `x::Real`: input value
    - `b::Real`: width of barrier
    - `V₀::Real`: height of barrier
"""
function kronig_penney_center(x,b::Real,V₀::Real)
    return sym_rect_pot_barr.(x,b,V₀)
end

"""
    kronig_penney_left(x,num_ions,a,b,V₀)

# Aim
    - Compute Kronig-Penney potential left

# Arguments
    - `x::Real`: input value
    - `num_ions::Integer`: number of ions
    - `a::Real`: distance between ions
    - `b::Real`: width of barrier
    - `V₀::Real`: height of barrier
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
    - Compute Kronig-Penney potential right

# Arguments
    - `x::Real`: input value
    - `num_ions::Integer`: number of ions
    - `a::Real`: distance between ions
    - `b::Real`: width of barrier
    - `V₀::Real`: height of barrier
"""
function kronig_penney_right(x,num_ions::Integer,a::Real,b::Real,V₀::Real)
    return kronig_penney_left(-x,num_ions,a,b,V₀)
end

"""
    symetric_kronig_penney(x,num_ions,a,b,V₀)

# Aim
    - Compute symetric Kronig-Penney potential

# Arguments
    - `x::Real`: input value
    - `num_ions::Integer`: number of ions
    - `a::Real`: distance between ions
    - `b::Real`: width of barrier
    - `V₀::Real`: height of barrier
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