module QuantumHarmonicOscillator2D

using Pkg
Pkg.activate("../"); # activate specific environment from TimeIndependentSchrodingerEquation package
Pkg.instantiate()

using TimeIndependentSchrodingerEquation
using Gridap,Test

# include auxiliary codes
include("../useful_functions_to_FEM.jl")
include("miscellaneous_functions.jl")

dom2D=(-30.0,30.0,-30.0,30.0);  # define domain 2D symetric concern origin
nxy=(100,100);                  # define number of FE to use in each coordinate
params_model=(dom2D,nxy);       # organize data in a tuple
grid_type="Cartesian2D";        # define type of FE grid

# create mesh from specific parameters
model2D=make_model(grid_type,params_model);

BC_type="FullDirichlet"; # define type of boundary conditions (BC)
# create BC with data and labels
FullDirichlet_values,FullDirichlet_tags=make_boundary_conditions(grid_type,BC_type,ComplexF64);

# create measures for interior and boundary domains
Ω,dΩ,Γ,dΓ = measures(model2D,3,FullDirichlet_tags);
# create FE polinomial interpolation
reff = ReferenceFE(lagrangian,Float64,2);
# create trial and test FE spaces
VSpace,USpace = fe_spaces(model2D,reff,grid_type;BC_type=BC_type,TypeData=ComplexF64,conf_type=:H1);

# create function to Sturm-Liouville formulation
p,q,r = eigenvalue_problem_functions((1.0,0.0,0.0));

# resolve eigen problem
ϵ,ϕ = eigen_values_and_eigen_vectors(p,q,r,dΩ,USpace,VSpace;params=(4,10e-9,500,:none,0.0));
# compute exact solution for eigenenergies
ϵ_real=exactly_eigenvalues_2dqho(length(ϵ));

# testing results

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
        @test Float64(orthogonality_check(ϕ,dΩ)[i]) ≈ zero(Float64) atol=Float64(0.1)
    end
end

end