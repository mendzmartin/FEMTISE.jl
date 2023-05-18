module QuantumHarmonicOscillator2D

using Pkg
Pkg.activate("../")

using TimeIndependentSchrodingerEquation
using Gridap,Test

include("../useful_functions_to_FEM.jl")
include("miscellaneous_functions.jl")

dom2D=(-25.0,25.0,-25.0,25.0);nxy=(50,50);params_model=((dom2D,nxy));
grid_type="Cartesian2D";

model2D=make_model(grid_type,params_model);

BC_type="FullDirichlet";
FullDirichlet_values,FullDirichlet_tags=make_boundary_conditions(grid_type,BC_type,ComplexF64);

Ω,dΩ,Γ,dΓ = measures(model2D,3,FullDirichlet_tags)

reff = ReferenceFE(lagrangian,Float64,2)

VSpace,USpace = fe_spaces(model2D,reff,grid_type;BC_type=BC_type,TypeData=ComplexF64)

p,q,r = eigenvalue_problem_functions((1.0,0.0,0.0))

ϵ,ϕ = eigen_values_and_eigen_vectors(p,q,r,dΩ,USpace,VSpace;params=(4,10e-9,500,:none,0.0))

ϵ_real=exactly_eigenvalues_2dqho(length(ϵ))

@testset "Check eigenenergies" begin
    @test real(ϵ) ≈ ϵ_real atol=0.1
end

@testset "Check normalization" begin
    for i in eachindex(ϕ)
        @test norm_l2(ϕ[1],dΩ) ≈ 1.0
    end
end

@testset "Check orthogonalization" begin
    for i in eachindex(ϕ)
        @test Float64(ortho_check_v2(ϕ,USpace,dΩ)[i]) ≈ zero(Float64) atol=Float64(1e-1)
    end
end

end