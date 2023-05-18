using Pkg
Pkg.activate("../")

using TimeIndependentSchrodingerEquation
using Test

@testset "Test Eigenvalues 2DQHO" begin
    include("test_2dqho/QuantumHarmonicOscillator2D.jl")
end

