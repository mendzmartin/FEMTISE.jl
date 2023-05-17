using Pkg
Pkg.activate("../")

using SteadyStateSchrodingerEquation
using Test

@testset "TestEigenvalues_2DQHO" begin
    include("Test_2DQHO/QuantumHarmonicOscillator2D.jl")
end

