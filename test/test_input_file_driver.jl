# Input-file driver, end to end.
#
# Writes a potential file and an input file in a temporary folder, runs
# `run_default_eigen_problem(set_type_potential(path))`, reloads the results with
# `collect_result_data` and compares them with the analytic spectrum of the
# harmonic oscillator. Covers the 1D and 2D drivers, the JLD2 and binary output
# formats and the parameter sweep (`analysis_param`).
#
# The input file has a fixed layout of 21 lines (see the "Input file reference"
# page of the documentation); `write_input_file` produces it from named fields.

"""
    write_input_file(path; fields...)

Write `path * ".dat"` with the 21-line layout expected by `input_data`.
Missing fields are left empty.
"""
function write_input_file(path; fields...)
    f = Dict{Symbol,Any}(fields)
    value(key) = string(get(f, key, ""))
    lines = [
        "full_path_name              = $(value(:full_path_name))",
        "dom_type                    = $(value(:dom_type))",
        "nev                         = $(value(:nev))",
        "dimension                   = $(value(:dimension))",
        "sigma                       = $(value(:sigma))",
        "adhoc_file_name             = $(value(:adhoc_file_name))",
        "potential_function_name     = $(value(:potential_function_name))",
        "params_potential_types      = $(value(:params_potential_types))",
        "params_potential            = $(value(:params_potential))",
        "analysis_param              = $(value(:analysis_param))",
        "output_format_type          = $(value(:output_format_type))",
        "## ONLY FOR 1D EIGENPROBLEMS",
        "L                           = $(value(:L))",
        "Δx                          = $(value(:Δx))",
        "## ONLY FOR 2D EIGENPROBLEMS",
        "Lx                          = $(value(:Lx))",
        "Ly                          = $(value(:Ly))",
        "nx                          = $(value(:nx))",
        "ny                          = $(value(:ny))",
        "different_masses            = $(value(:different_masses))",
        "reduced_density             = $(value(:reduced_density))",
    ]
    write(path * ".dat", join(lines, "\n") * "\n")
end

const POTENTIAL_SOURCE = """
    # Harmonic potentials used by the driver tests; params = (ω, center...)
    function driver_test_qho_1d(x, params::Tuple)
        ω, x0 = params
        return 0.5 * ω^2 * (x[1] - x0)^2
    end
    function driver_test_qho_2d(x, params::Tuple)
        ω, x0, y0 = params
        return 0.5 * ω^2 * ((x[1] - x0)^2 + (x[2] - y0)^2)
    end
    """

"""
    run_driver(folder, input_path)

Run the input-file driver inside `folder` (the 1D mesher writes a temporary
`.msh` file in the working directory).

`set_type_potential` defines the potential function by `include`-ing the
potential file. This function was compiled before that definition existed, so
the driver must be called in the latest "world age" with `Base.invokelatest`.
User scripts that call the driver at top level do not need this (see the
"Input file reference" page of the documentation).
"""
function run_driver(folder, input_path)
    cd(folder) do
        simulation = set_type_potential(input_path)
        Base.invokelatest(run_default_eigen_problem, simulation)
    end
end

@testset "Input-file driver" begin
    mktempdir() do folder
        write(joinpath(folder, "potentials.jl"), POTENTIAL_SOURCE)
        potential_file = joinpath(folder, "potentials")    # without ".jl", as in the input file

        @testset "1D, JLD2 output" begin
            input = joinpath(folder, "input_1d")
            write_input_file(input; full_path_name = joinpath(folder, "qho1d"), dom_type = "s",
                             nev = 5, dimension = "1D", sigma = 0.0,
                             adhoc_file_name = potential_file,
                             potential_function_name = "driver_test_qho_1d",
                             params_potential_types = "f f", params_potential = "1.0 0.0",
                             analysis_param = false, output_format_type = "jld2 eigen",
                             L = 20.0, Δx = 0.1)
            run_driver(folder, input)
            id, results = collect_result_data(true, input)

            ε = sort(real.(results.ϵ))
            @test id isa InputData1D
            @test ε ≈ [n + 0.5 for n in 0:4] rtol = 1e-4
            @test all(ε .> [n + 0.5 for n in 0:4])           # upper bounds
            @test length(results.ϕ) == 5
        end

        @testset "1D, parameter sweep with binary output" begin
            input = joinpath(folder, "input_sweep")
            write_input_file(input; full_path_name = joinpath(folder, "sweep"), dom_type = "s",
                             nev = 2, dimension = "1D", sigma = 0.0,
                             adhoc_file_name = potential_file,
                             potential_function_name = "driver_test_qho_1d",
                             params_potential_types = "f f", params_potential = "1.0 0.0",
                             analysis_param = "1 0.5 1.0 0.25",    # sweep ω
                             output_format_type = "bin eigen", L = 20.0, Δx = 0.1)
            run_driver(folder, input)
            _, results = collect_result_data(true, input)

            ω = real.(results.λvector)
            # Documented behavior of v0.1.0: the sweep starts at λᵢ + Δλ.
            @test ω ≈ [0.75, 1.0]
            ground = real.(results.ϵ_matrix[1, :])
            @test ground ≈ ω ./ 2 rtol = 1e-4                # ε₀ = ω/2
        end

        @testset "2D, binary output with reduced densities" begin
            input = joinpath(folder, "input_2d")
            write_input_file(input; full_path_name = joinpath(folder, "qho2d"), dom_type = "s",
                             nev = 3, dimension = "2D", sigma = 0.0,
                             adhoc_file_name = potential_file,
                             potential_function_name = "driver_test_qho_2d",
                             params_potential_types = "f f f", params_potential = "1.0 0.0 0.0",
                             analysis_param = false, output_format_type = "bin eigen",
                             Lx = 12.0, Ly = 12.0, nx = 40, ny = 40,
                             different_masses = false, reduced_density = true)
            run_driver(folder, input)

            # Known issue in v0.1.0: the containers of the "bin" format declare the
            # eigenfunctions as Vector{CellField}, so `collect_result_data` cannot
            # rebuild them from the sampled matrix. When fixed, replace
            # @test_broken by @test and use `collect_result_data` below.
            reload_works = try
                collect_result_data(true, input)
                true
            catch
                false
            end
            @test_broken reload_works

            # Documented workaround: read the binary files directly.
            prefix = joinpath(folder, "qho2d")
            nev = 3
            ε = sort(real.(read_bin(prefix * "_eigen_values.bin")))
            @test ε ≈ [1.0, 2.0, 2.0] rtol = 1e-3

            # Reduced density of the ground state, sampled on the output grid
            # (first coordinate column). Exact result: ρ⁽¹⁾(x) = e^{-x²}/√π.
            ε_unsorted = real.(read_bin(prefix * "_eigen_values.bin"))
            ground = argmin(ε_unsorted)
            x = read_bin(prefix * "_coordinates.bin")[1:40]
            Δx = x[2] - x[1]
            ρ1 = read_bin(prefix * "_reduced_density_DOF1.bin"; matrix_data = true, c_dim = nev)[:, ground]
            integral = sum(ρ1) * Δx

            # Shape: after normalization, ρ⁽¹⁾ matches the Gaussian.
            @test maximum(abs.(ρ1 ./ integral .- exp.(-x .^ 2) ./ sqrt(π))) < 1e-2
            # Known issue in v0.1.0: `reduced_integration` normalizes with the span of
            # the sampling points, r[end] - r[1] = L - Δx, instead of the domain length
            # L, so the integral is (L - Δx)/L = 0.975 instead of 1.
            @test integral ≈ 1 - Δx / 12.0 atol = 2e-3
            @test_broken integral ≈ 1.0 atol = 2e-3
        end
    end
end
