using Pkg
Pkg.activate("../")

using TimeIndependentSchrodingerEquation
using Test

@testset "Test eigenvalue problem of 2Dqho" begin
    include("test_2dqho/QuantumHarmonicOscillator2D.jl")
end

