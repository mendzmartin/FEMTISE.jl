# Optional settings of the input file (v0.2.0) and the extended driver.
#
# Uses `write_input_file`, `POTENTIAL_SOURCE` and `run_driver` from
# test_input_file_driver.jl. Every case is compared with an exact result:
#   * masses:            oscillator of mass M = 4, ε₀ = 1/(2√M);
#   * potential data:    tabulated oscillator, ε₀ = 1/2;
#   * two-parameter sweep: ε₀(ω, x₀) = ω/2;
#   * 3D:                isotropic 3D oscillator, ε₀ = 3/2.

"Append optional `key = value` lines to an input file."
function append_options(path, options::Pair...)
    open(path * ".dat", "a") do io
        println(io, "## OPTIONAL SETTINGS (v0.2.0)")
        for (key, value) in options
            println(io, rpad(key, 28), "= ", value)
        end
    end
end

@testset "Extended input-file driver" begin
    mktempdir() do folder
        write(joinpath(folder, "potentials.jl"), POTENTIAL_SOURCE * """
            function driver_test_qho_3d(x, params::Tuple)
                ω = params[1]
                return 0.5 * ω^2 * (x[1]^2 + x[2]^2 + x[3]^2)
            end
            """)
        potential_file = joinpath(folder, "potentials")
        common_1d = (dom_type = "s", dimension = "1D", sigma = 0.0, adhoc_file_name = potential_file,
                     potential_function_name = "driver_test_qho_1d", params_potential_types = "f f",
                     params_potential = "1.0 0.0", output_format_type = "jld2 eigen", L = 20.0, Δx = 0.1)

        @testset "Particle mass" begin
            input = joinpath(folder, "mass")
            write_input_file(input; merge(common_1d, (; full_path_name = joinpath(folder, "mass"), nev = 3,
                             analysis_param = false))...)
            append_options(input, "masses" => "4.0", "quadrature_degree" => "6")
            run_driver(folder, input)
            id, results = collect_result_data(true, input)
            @test id.options.masses == [4.0]
            @test real.(results.ϵ) ≈ [(n + 0.5) / 2 for n in 0:2] rtol = 1e-5
        end

        @testset "Potential data file" begin
            data = joinpath(folder, "qho.dat")
            open(data, "w") do io
                for x in -10.0:0.05:10.0
                    println(io, x, " ", 0.5x^2)
                end
            end
            input = joinpath(folder, "tabulated")
            write_input_file(input; merge(common_1d, (; full_path_name = joinpath(folder, "tabulated"), nev = 3,
                             adhoc_file_name = "", potential_function_name = "",
                             params_potential_types = "", params_potential = "", analysis_param = false))...)
            append_options(input, "potential_data_file" => data, "interpolation" => "cubic")
            run_driver(folder, input)
            _, results = collect_result_data(true, input)
            @test real.(results.ϵ) ≈ [0.5, 1.5, 2.5] rtol = 1e-4
        end

        @testset "Two-parameter sweep" begin
            input = joinpath(folder, "sweep2")
            # ω = 1.0, 1.5 (first sweep) and x₀ = 0.0, 1.0 (second sweep)
            write_input_file(input; merge(common_1d, (; full_path_name = joinpath(folder, "sweep2"), nev = 2,
                             analysis_param = "1 0.5 1.5 0.5"))...)
            append_options(input, "analysis_param_2" => "2 -1.0 1.0 1.0")
            run_driver(folder, input)
            _, results = collect_result_data(true, input)
            @test results isa AnalysisParams2D
            @test results.λvector ≈ [1.0, 1.5]
            @test results.μvector ≈ [0.0, 1.0]
            @test size(results.ϵ_tensor) == (2, 2, 2)
            for j in 1:2
                @test results.ϵ_tensor[1, :, j] ≈ results.λvector ./ 2 rtol = 1e-4
            end
        end

        @testset "Three dimensions" begin
            input = joinpath(folder, "qho3d")
            write_input_file(input; full_path_name = joinpath(folder, "qho3d"), dom_type = "s", nev = 4,
                             dimension = "3D", sigma = 0.0, adhoc_file_name = potential_file,
                             potential_function_name = "driver_test_qho_3d",
                             params_potential_types = "f", params_potential = "1.0",
                             analysis_param = false, output_format_type = "jld2 eigen",
                             Lx = 12.0, Ly = 12.0, nx = 14, ny = 14, different_masses = false,
                             reduced_density = false)
            append_options(input, "Lz" => "12.0", "nz" => "14", "quadrature_degree" => "6")
            run_driver(folder, input)
            id, results = collect_result_data(true, input)
            @test id isa InputData3D
            # the triply degenerate level must be complete (guard states of solve_tise)
            @test real.(results.ϵ) ≈ [1.5, 2.5, 2.5, 2.5] rtol = 5e-3
        end

        @testset "Errors are reported, not ignored" begin
            input = joinpath(folder, "bad")
            write_input_file(input; merge(common_1d, (; full_path_name = joinpath(folder, "bad"), nev = 2,
                             analysis_param = false))...)
            append_options(input, "massess" => "4.0")                   # misspelled option
            @test_throws ArgumentError input_data(input)
            write_input_file(input; merge(common_1d, (; full_path_name = joinpath(folder, "bad"), nev = 2,
                             analysis_param = false, dimension = "4D"))...)
            @test_throws ArgumentError input_data(input)
        end
    end
end
