using Pkg
Pkg.activate("../")

using TimeIndependentSchrodingerEquation
using Test

@testset "Test eigenvalue problem of 2Dqho" begin
    include("test_2dqho/QuantumHarmonicOscillator2D.jl")
end

@testset "Test eigenvalue problem of kronig-penney 1D" begin
    include("test_1d_kronig_penney/Finite1DKronigPenneyPotential.jl")
end

