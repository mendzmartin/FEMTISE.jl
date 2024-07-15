using Pkg
Pkg.activate("../")
Pkg.instantiate()

# Load libraries
using FEMTISE
using Test

# Include useful functions
@testset "Test eigenvalue problem for Quantum Harmonic Oscillator 2D" begin
    include("test_2d_qho/QuantumHarmonicOscillator2D.jl")
end
@testset "Test eigenvalue problem for Finite Kronig-Penney potential 1D" begin
    include("test_1d_kronig_penney/FiniteKronigPenneyPotential1D.jl")
end

@testset "Test eigenvalue problem for Finite Well Potential 1D" begin
    include("test_1d_fwp/FiniteWellPotential1D.jl")
end

@testset "Test eigenvalue problem for default potentials" begin
    include("test_default_eigen_problem/DefaultEigenProblem.jl")
end