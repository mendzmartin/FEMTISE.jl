module DefaultEigenProblem

using Pkg
Pkg.activate("../")
Pkg.instantiate()

# Load libraries
using FEMTISE
using Gridap,Test

# Include useful functions
include("../useful_functions_to_FEM.jl")
include("./default_potential_functions.jl")
include("./miscellaneous_functions.jl")

# Define type of potential
type_potential="qho_1d"

# Define parameters of the potential
if type_potential=="qho_1d"
    L=30.0;Δx=0.1;ω=1.0;x₁=0.0;nev=4;sigma=0.0;
    params=(L,Δx,ω,x₁,nev,sigma)
elseif type_potential=="qho_2d"
    L=30.0;n=100;ω=1.0;x₁=1.0;y₁=0.0;nev=4;sigma=0.0;
    params=(L,n,n,ω,x₁,y₁,nev,sigma)
elseif type_potential=="kronig_penney_1d"
    L=30.0;V₀=-10;t=0.1;b=2.0;num_ions=3;Δx=0.1;nev=4;sigma=V₀;
    params=(L,V₀,t,b,num_ions,Δx,nev,sigma)
elseif type_potential=="finite_well_1d"
    L=30.0;V₀=-10;a=2.1;Δx=0.1;nev=4;sigma=V₀;
    params=(L,V₀,a,Δx,nev,sigma)
end

ϵ,ϕ = default_solver_eigen_problem(type_potential,params)

end

#=
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        README
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

This module defines the default eigenvalue problem for the quantum harmonic oscillator 1D, 
quantum harmonic oscillator 2D, Kronig-Penney 1D, and finite well potential 1D.

# Set atomic unit system
const Bohr_radius_meter=5.29177210903e−11;                              # [m]
const Angstrom_to_meter=1e−10;                                          # [m/Å]
const Angstrom_to_au=Angstrom_to_meter*(1.0/Bohr_radius_meter);         # [au/Å]
const Nanometer_to_au=(1e-9)*(1.0/Angstrom_to_meter)*Angstrom_to_au;    # [au/nm]
const Electronvolt_to_au=0.036749405469679;                             # [au/Ev]

set type of default potential could be the next
    - Quantum Harmonic Oscillator 1D
        - type_potential="qho_1d"
        -  default values -> L=30.0,Δx=0.1,sigma=0.0
    - Quantum Harmonic Oscillator 2D
        - type_potential="qho_2d"
        -  default values -> L=30.0,nx=ny=100,sigma=0.0
    - Finit Kronig-Penney 1D
        - type_potential="kronig_penney_1d"
        -  default values -> L=30.0,Δx=0.1,sigma=V₀
    - Finit Well Potential 1D
        - type_potential="finite_well_1d"
        -  default values -> L=30.0,Δx=0.1,sigma=V₀
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
=#