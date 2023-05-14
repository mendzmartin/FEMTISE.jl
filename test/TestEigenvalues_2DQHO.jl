using Pkg
Pkg.activate("../")

using Gridap
using GridapGmsh
using Gmsh
using FileIO
using LinearAlgebra
using SparseArrays
using SuiteSparse
using Arpack

include("../src/SteadyStateSchrodingerEquation.jl")
using .SteadyStateSchrodingerEquation

dom2D=(-25.0,25.0,-25.0,25.0);nxy=(50,50);params_model=((dom2D,nxy));
grid_type="Cartesian2D";

model2D=make_model(grid_type,params_model);

BC_type="FullDirichlet";
FullDirichlet_values,FullDirichlet_tags=make_boundary_conditions(grid_type,BC_type,ComplexF64);

Ω,dΩ,Γ,dΓ = measures(model2D,3,FullDirichlet_tags)

reff = ReferenceFE(lagrangian,Float64,2)

VSpace,USpace = FESpaces(model2D,reff,grid_type;BC_type=BC_type,TypeData=ComplexF64)

const ħ=1.0;
const m=1.0;
function eigenvalue_problem_functions(params;switch_potential::String = "QHO_1D")
    if (switch_potential == "QHO_1D")
        # caso de potencial tipo quantum harmonic oscillator 1D (QHO)
        println("Set quantum harmonic oscillator 1D potential");
        ω,x₁=params;
        p_QHO_1D(x) = 0.5*(ħ*ħ)*(1.0/m);                                      # factor para energía cinética
        q_QHO_1D(x) = 0.5*m*(ω*ω)*(x[1]-x₁)*(x[1]-x₁);                        # oscilador armónico 1D centrado en x₁
        r_QHO_1D(x) = 1.0;
        return p_QHO_1D,q_QHO_1D,r_QHO_1D
    elseif (switch_potential == "QHO_2D")
        # caso de potencial tipo quantum harmonic oscillator 2D (QHO)
        println("Set quantum harmonic oscillator 2D potential");
        ω,x₁,y₁=params;
        p_QHO_2D(x) = 0.5*(ħ*ħ)*(1.0/m);                                       # factor para energía cinética
        q_QHO_2D(x) = 0.5*m*(ω*ω)*((x[1]-x₁)*(x[1]-x₁)+(x[2]-y₁)*(x[2]-y₁));   # oscilador armónico 2D centrado en (x₁,y₁)
        r_QHO_2D(x) = 1.0;
        return p_QHO_2D,q_QHO_2D,r_QHO_2D
    end
end

p,q,r = eigenvalue_problem_functions((1.0,0.0,0.0);switch_potential = "QHO_2D")

ϵ,ϕ = EigenValuesAndEigenVectors(p,q,r,dΩ,USpace,VSpace;params=(4,10e-9,500,:none,0.0))

function exactly_eigenvalues_2DQHO(num_eigval::Integer)
    ϵ_real_aux=Array{Float64}(undef, num_eigval*num_eigval)
    index::Int64=1
    for i in 1:num_eigval
        for j in 1:num_eigval
            ϵ_real_aux[index]=((i-1)+(j-1)+1)
            index+=1
        end
    end
    ϵ_real_aux=sort(ϵ_real_aux);
    ϵ_real = ϵ_real_aux[1:num_eigval];
    return ϵ_real
end

ϵ_real=exactly_eigenvalues_2DQHO(length(ϵ))

using Test
@test real(ϵ) ≈ ϵ_real atol=0.1