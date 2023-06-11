using Pkg
Pkg.activate("../")

using TimeIndependentSchrodingerEquation
using Test

@testset "Test eigenvalue problem for Quantum Harmonic Oscillator 2D" begin
    include("test_2d_qho/QuantumHarmonicOscillator2D.jl")
end

@testset "Test eigenvalue problem for Finite Kronig-Penney potential 1D" begin
    include("test_1d_kronig_penney/FiniteKronigPenneyPotential1D.jl")
end

@testset "Test eigenvalue problem for Finite Well Potential 1D" begin
    include("test_1d_fwp/FiniteWellPotential1D.jl")
end