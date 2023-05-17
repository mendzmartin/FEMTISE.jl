module QuantumHarmonicOscillator2D

using Pkg
Pkg.activate("../")

using SteadyStateSchrodingerEquation
using Gridap,Test

include("MiscellaneousFunctions.jl")

dom2D=(-25.0,25.0,-25.0,25.0);nxy=(50,50);params_model=((dom2D,nxy));
grid_type="Cartesian2D";

model2D=make_model(grid_type,params_model);

BC_type="FullDirichlet";
FullDirichlet_values,FullDirichlet_tags=make_boundary_conditions(grid_type,BC_type,ComplexF64);

Ω,dΩ,Γ,dΓ = measures(model2D,3,FullDirichlet_tags)

reff = ReferenceFE(lagrangian,Float64,2)

VSpace,USpace = FESpaces(model2D,reff,grid_type;BC_type=BC_type,TypeData=ComplexF64)

p,q,r = eigenvalue_problem_functions((1.0,0.0,0.0);switch_potential = "QHO_2D")

ϵ,ϕ = EigenValuesAndEigenVectors(p,q,r,dΩ,USpace,VSpace;params=(4,10e-9,500,:none,0.0))

ϵ_real=exactly_eigenvalues_2DQHO(length(ϵ))

@test real(ϵ) ≈ ϵ_real atol=0.1

end