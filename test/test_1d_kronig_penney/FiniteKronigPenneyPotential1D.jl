module FiniteKronigPenneyPotential1D

using Pkg
Pkg.activate("../")

using TimeIndependentSchrodingerEquation
using Gridap,Test

include("../useful_functions_to_FEM.jl")
include("miscellaneous_functions.jl")


# Set atomic unit system
const Bohr_radius_meter=5.29177210903e−11;                              # [m]
const Angstrom_to_meter=1e−10;                                          # [m/Å]
const Angstrom_to_au=Angstrom_to_meter*(1.0/Bohr_radius_meter);         # [au/Å]
const Nanometer_to_au=(1e-9)*(1.0/Angstrom_to_meter)*Angstrom_to_au;    # [au/nm]
const Electronvolt_to_au=0.036749405469679;                             # [au/Ev]

L=10*Nanometer_to_au;
V₀=-30.0*Electronvolt_to_au;
t=0.3*Nanometer_to_au;
b=0.15*Nanometer_to_au;
a=t+b;
num_ions=21;
Δx=0.001*Nanometer_to_au;

grid_type="simple_line";
params_model=("./","model1D",(-0.5*L,0.5*L),Δx);
model1D=make_model(grid_type,params_model);
rm(params_model[1]*params_model[2]*".msh")

BC_type="FullDirichlet";
FullDirichlet_values,FullDirichlet_tags=make_boundary_conditions(grid_type,BC_type,ComplexF64);

Ω,dΩ,Γ,dΓ = measures(model1D,3,FullDirichlet_tags)

reff = ReferenceFE(lagrangian,Float64,2)

VSpace,USpace = fe_spaces(model1D,reff,grid_type;BC_type=BC_type,TypeData=ComplexF64)

p,q,r = kronig_penney_sturm_liouville((num_ions,a,b,V₀))

ϵ,ϕ = eigen_values_and_eigen_vectors(p,q,r,dΩ,USpace,VSpace;params=(50,1e-9,500,:none,V₀))

@testset "Check eigenenergies 1st band" begin
    for i in eachindex(ϵ[1:21])
        @test Float64(-0.828775) ≤ real(ϵ[i]) ≤ Float64(-0.828381)
    end
end

@testset "Check eigenenergies 2nd band" begin
    for i in eachindex(ϵ[22:42])
        @test Float64(-0.164768) ≤ real(ϵ[i+21]) ≤ Float64(-0.138487)
    end
end

@testset "Check normalization" begin
    for i in eachindex(ϕ)
        @test norm_l2(ϕ[1],dΩ) ≈ 1.0
    end
end

end