module FiniteWellPotential1D

using Pkg
Pkg.activate("../")
Pkg.instantiate()

# Load libraries
using FEMTISE
using Gridap,Test

# Include useful functions
include("../useful_functions_to_FEM.jl")
include("../test_1d_kronig_penney/miscellaneous_functions.jl")


# Set atomic unit system
const Bohr_radius_meter=5.29177210903e−11;                              # [m]
const Angstrom_to_meter=1e−10;                                          # [m/Å]
const Angstrom_to_au=Angstrom_to_meter*(1.0/Bohr_radius_meter);         # [au/Å]
const Nanometer_to_au=(1e-9)*(1.0/Angstrom_to_meter)*Angstrom_to_au;    # [au/nm]
const Electronvolt_to_au=0.036749405469679;                             # [au/Ev]

# Set parameters
L=10*Nanometer_to_au;
V₀=-10.0*Electronvolt_to_au;
a=4.2*Nanometer_to_au;
Δx=0.01*Nanometer_to_au;

# Set model
grid_type="simple_line";
params_model=("./","model1D",(-0.5*L,0.5*L),Δx);
model1D=make_model(grid_type,params_model);
rm(params_model[1]*params_model[2]*".msh")

# Set boundary conditions
BC_type="FullDirichlet";
FullDirichlet_values,FullDirichlet_tags=make_boundary_conditions(grid_type,BC_type,ComplexF64);

# Set measures
Ω,dΩ,Γ,dΓ = measures(model1D,3,FullDirichlet_tags)

# Set reference finite element
reff = ReferenceFE(lagrangian,Float64,2)

# Set finite element spaces
VSpace,USpace = fe_spaces(model1D,reff,grid_type;BC_type=BC_type,TypeData=ComplexF64)

# Set bilinear form
p,q,r = kronig_penney_sturm_liouville((a,V₀);fwp=true)

# Set eigenvalues and eigenvectors
ϵ,ϕ = eigen_values_and_eigen_vectors(p,q,r,dΩ,USpace,VSpace;params=(10,1e-9,500,:none,V₀))


@testset "Check normalization" begin
    for i in eachindex(ϕ)
        @test norm_l2(ϕ[1],dΩ) ≈ 1.0
    end
end

end