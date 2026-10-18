# =============================================================================
# Extended input-file driver (v0.2.0).
#
# Input files that use any optional setting (masses, element degree, quadrature,
# tabulated potential, two-parameter sweep) or describe a 3D problem are solved
# here, with the general solver `solve_tise` / `sweep_eigenvalues`. Input files
# without optional settings keep using the original driver, so their results do
# not change.
# =============================================================================

"True if the input file needs the extended driver (optional settings or 3D)."
uses_extended_driver(id::Union{InputData1D,InputData2D}) = id.options != DriverOptions()
uses_extended_driver(::InputData3D) = true
uses_extended_driver(::Any) = false

"Domain `(x₁, x₂, y₁, y₂, ...)` of the Cartesian mesh described by the input parameters."
function driver_domain(params)
    sides = params isa Params1D ? (params.L,) :
            params isa Params2D ? (params.Lx, params.Ly) : (params.Lx, params.Ly, params.Lz)
    params.dom_type in ("s", "ns") || throw(ArgumentError(
        "the extended driver supports dom_type = s (symmetric) or ns (starting at 0), got $(params.dom_type)"))
    symmetric = params.dom_type == "s"
    return Tuple(Iterators.flatten(symmetric ? (-L / 2, L / 2) : (0.0, L) for L in sides))
end

"Number of cells in each direction."
driver_cells(params::Params1D) = (round(Int, params.L / params.Δx),)
driver_cells(params::Params2D) = (params.nx, params.ny)
driver_cells(params::Params3D) = (params.nx, params.ny, params.nz)

"Masses of the coordinates, from the optional `masses` line or the legacy `different_masses` field."
function driver_masses(id)
    dimension = num_coordinates(id.params)
    legacy_ratio = id isa InputData2D && id.different_masses[1] ? id.different_masses[2] : nothing
    masses = id.options.masses
    if isnothing(masses)
        # different_masses = μ sets p₂ = μ p₁, i.e. M₂ = 1/μ (electron mass for x₁)
        return isnothing(legacy_ratio) ? ntuple(_ -> 1.0, dimension) : (1.0, 1.0 / legacy_ratio)
    end
    isnothing(legacy_ratio) || throw(ArgumentError(
        "give the masses either with `different_masses` or with the optional `masses` line, not both"))
    masses = length(masses) == 1 ? fill(masses[1], dimension) : masses
    kinetic_coefficients(masses, dimension)          # checks number and sign of the masses
    return Tuple(masses)
end

num_coordinates(::Params1D) = 1
num_coordinates(::Params2D) = 2
num_coordinates(::Params3D) = 3

"""
    driver_potential(id)

Potential family `V(x, λ, μ)` of the input file: the tabulated potential if a data
file is given (it does not depend on λ, μ), otherwise the user function
`potential_function_name(x, params)` with the swept parameters replaced by λ and μ.
"""
function driver_potential(id)
    if !isempty(id.options.potential_data_file)
        id.analysis_param == false && id.options.analysis_param_2 == false || throw(ArgumentError(
            "parameter sweeps need a potential function with parameters, not a potential data file"))
        V = tabulated_potential(id.options.potential_data_file; interpolation = id.options.interpolation)
        return (x, λ, μ) -> V(x)
    end
    name = Symbol(id.params.potential_function_name)
    isdefined(@__MODULE__, name) || throw(ArgumentError(
        "the potential function `$name` is not defined; check potential_function_name and adhoc_file_name"))
    f = getproperty(@__MODULE__, name)
    base = id.params.params_potential
    index1 = id.analysis_param == false ? 0 : id.analysis_param.λindex
    index2 = id.options.analysis_param_2 == false ? 0 : id.options.analysis_param_2.λindex
    index1 != 0 && index1 == index2 && throw(ArgumentError("the two swept parameters must be different"))
    with_parameters(λ, μ) = ntuple(k -> k == index1 ? λ : k == index2 ? μ : base[k], length(base))
    return (x, λ, μ) -> f(x, with_parameters(λ, μ))
end

"Values λᵢ + kΔλ, k = 1, …, round(|λf − λi|/Δλ), of a parameter sweep (same rule as the original driver)."
sweep_values(p::AnalysisParam) = [p.λi + k * p.Δλ for k in 1:round(Int, abs(p.λf - p.λi) / p.Δλ)]

"Regular sampling grid of the outputs (same points as the original driver for 1D and 2D)."
function sampling_grid(params, domain)
    cells = driver_cells(params)
    lower = domain[1:2:end]
    spacing = [(domain[2k] - domain[2k-1]) / cells[k] for k in eachindex(cells)]
    coordinates = Tuple([lower[k] + spacing[k] * (i - 1) for i in 1:cells[k]] for k in eachindex(cells))
    points = length(cells) == 1 ? [Point(x) for x in coordinates[1]] :
             length(cells) == 2 ? [Point(x, y) for x in coordinates[1] for y in coordinates[2]] :
                                  [Point(x, y, z) for x in coordinates[1] for y in coordinates[2]
                                                  for z in coordinates[3]]
    return (length(cells) == 1 ? coordinates[1] : coordinates), points
end

"""
    run_extended_eigen_problem(id)

Solve the problem described by an input file that uses optional settings or is
three-dimensional, and write the results to `id.full_path_name * "_eigen_data.jld2"`.
Called by [`run_default_eigen_problem`](@ref); not meant to be called directly.
"""
function run_extended_eigen_problem(id)
    id.output_format_type == ("jld2", "eigen") || throw(ArgumentError(
        "input files with optional settings or dimension = 3D need output_format_type = jld2 eigen"))
    domain = driver_domain(id.params)
    model = CartesianDiscreteModel(domain, driver_cells(id.params))
    potential = driver_potential(id)
    settings = (masses = driver_masses(id), order = id.options.element_order,
                quadrature_degree = id.options.quadrature_degree, nev = id.params.nev,
                sigma = id.params.sigma)
    output_file = id.full_path_name * "_eigen_data.jld2"
    rm_existing_file(output_file)
    sweep1, sweep2 = id.analysis_param, id.options.analysis_param_2
    sweep2 != false && sweep1 == false && throw(ArgumentError(
        "analysis_param_2 needs analysis_param (the first swept parameter)"))

    println("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    if sweep1 == false
        println("... Running (extended driver, $(id.params.dimension)) ...")
        solution = solve_tise(model, x -> potential(x, nothing, nothing); settings...)
        r, pts = sampling_grid(id.params, domain)
        data = (ϵ = ComplexF64.(solution.energies), ϕ = solution.states, r = r, pts = pts)
        if id isa InputData2D && id.reduced_density
            rho1, rho2 = reduced_density(solution.states, r, model)
            data = merge(data, (rhoDOF1 = rho1, rhoDOF2 = rho2))
        end
        jldsave(output_file; data...)
    elseif sweep2 == false
        λ = sweep_values(sweep1)
        println("... Running parameter sweep: $(length(λ)) values (extended driver) ...")
        energies = sweep_eigenvalues(model, (x, λ) -> potential(x, λ, nothing), λ; settings...)
        jldsave(output_file; ϵ_matrix = ComplexF64.(energies), λvector = ComplexF64.(λ))
    else
        λ, μ = sweep_values(sweep1), sweep_values(sweep2)
        println("... Running two-parameter sweep: $(length(λ)) × $(length(μ)) values ...")
        energies = sweep_eigenvalues(model, potential, λ, μ; settings...)
        jldsave(output_file; ϵ_tensor = energies, λvector = λ, μvector = μ)
    end
    println("Saved data in $output_file.")
    return nothing
end
