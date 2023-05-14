#= using SteadyStateSchrodingerEquation
using Test

@testset "SteadyStateSchrodingerEquation.jl" begin
    # Write your tests here.
end =#

# using SteadyStateSchrodingerEquation
using Test

@testset "TestEigenvalues_2DQHO" begin
    include("TestEigenvalues_2DQHO.jl")
end

