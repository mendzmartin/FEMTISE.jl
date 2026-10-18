# Test suite of FEMTISE.jl
#
# Run with `pkg> test FEMTISE` (or `julia --project -e 'using Pkg; Pkg.test()'`).
#
# Organization:
#   * Validation tests (one file per property, each with its physical reference
#     stated at the top of the file):
#       test_low_level_api.jl            upper bounds, h^(2ℓ) convergence, orthonormality
#       test_position_dependent_mass.jl  BenDaniel–Duke operator vs exact levels
#       test_two_masses.jl               anisotropic kinetic term (different masses)
#       test_input_file_driver.jl        input-file driver end to end (1D, 2D, sweeps)
#       test_post_processing.jl          densities, moments, entropies
#       test_binary_io.jl                write_bin / read_bin round trip
#       test_v020_features.jl            masses, tabulated potentials, sweeps,
#                                        autocorrelation, Pauli–Fierz (3D), plots
#       test_extended_driver.jl          optional settings of the input file
#   * Example tests (test_*/ folders): complete worked examples kept from earlier
#     versions; they run the package as a user would.

using FEMTISE
using Gridap
using LinearAlgebra
import Makie       # loads the plotting extension of FEMTISE (import avoids name clashes)
using Test

include("helpers.jl")

@testset "FEMTISE.jl" begin
    @testset "Validation" begin
        include("test_low_level_api.jl")
        include("test_position_dependent_mass.jl")
        include("test_two_masses.jl")
        include("test_post_processing.jl")
        include("test_binary_io.jl")
        include("test_input_file_driver.jl")
    end

    @testset "v0.2.0 features" begin
        include("test_v020_features.jl")
        include("test_extended_driver.jl")
    end

    @testset "Examples" begin
        @testset "2D quantum harmonic oscillator" begin
            include("test_2d_qho/QuantumHarmonicOscillator2D.jl")
        end
        @testset "Finite Kronig-Penney potential 1D" begin
            include("test_1d_kronig_penney/FiniteKronigPenneyPotential1D.jl")
        end
        @testset "Finite well potential 1D" begin
            include("test_1d_fwp/FiniteWellPotential1D.jl")
        end
        @testset "Default potentials" begin
            include("test_default_eigen_problem/DefaultEigenProblem.jl")
        end
    end
end
