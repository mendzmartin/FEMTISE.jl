module FEMTISE

using Gridap, GridapGmsh, Gmsh, FileIO
using Gridap.CellData
using LinearAlgebra, SparseArrays, Arpack
using DelimitedFiles, SpecialFunctions
using JLD2

# FEMTISE adds a method `solve(::EigenProblem)` to Gridap's own `solve`, so that
# `using FEMTISE, Gridap` exposes a single, unambiguous `solve` function.
import Gridap.Algebra: solve

include("functions/default_running_structures.jl")
include("functions/post_proccesing_data_structures.jl")
include("functions/mesh_generator_function.jl")
include("functions/boundary_conditions_function.jl")
include("functions/miscellaneous_functions.jl")
include("functions/eigen_problem_definition_function.jl")
include("functions/eigen_problem_solve_function.jl")
include("functions/useful_functions_for_FEM_objects.jl")
include("functions/default_running_miscellaneous_functions.jl")
include("functions/default_running_input_data.jl")
include("functions/default_eigen_problem.jl")
include("functions/binary_file_io.jl")
include("functions/post_proccesing_data_function.jl")
# since v0.2.0
include("functions/general_solver.jl")
include("functions/tabulated_potentials.jl")
include("functions/autocorrelation.jl")
include("functions/pauli_fierz.jl")
include("functions/extended_driver.jl")
include("functions/plotting.jl")

# from default_running_structures
export Params1D,Params2D,Params3D,AnalysisParam,InputData,InputData1D,InputData2D,InputData3D,
 DriverOptions
# from post_proccesing_data_structures
export DefaultBinEigenProblemReducedDensity,DefaultJLD2EigenProblem,
 DefaultJLD2EigenProblemReducedDensity,DefaultJLD2AllEigenProblem,
 DefaultJLD2AllEigenProblemReducedDensity,AnalysisParams,DefaultBinEigenProblem,AnalysisParams2D
# from mesh_generator_function
export make_model
# from boundary_conditions_function
export make_boundary_conditions
# from miscellaneous_functions
export measures,space_coord,bilineal_forms,fe_spaces
# from eigen_problem_definition_function
export eigen_problem,solve
# from eigen_problem_solve_function
export eigen_values_and_eigen_vectors
# from useful_functions_for_FEM_objects
export uniform_trapezoidal_integration_method,aprox_dirac_delta,reduced_integration,
 norm_l2,orthogonality_check,eigenstates_normalization,
 density,reduced_density,
 time_indep_entropy,reduced_time_indep_entropy,time_indep_diff_mutual_information,
 state_population,
 coord_first_moment,coord_second_moment,coord_standar_deviation
# from default_running_miscellaneous_functions
export default_solver_eigen_problem,ask_for_params,
 get_input,num_ions_BoolCondition,rm_existing_file,
 create_and_remove_model,solver_eigen_problem_with_analysis_param
# from default_running_input_data
export input_data
# from default_eigen_problem
export run_default_eigen_problem,run_default_eigen_problem_jld2,
 set_include,set_type_potential
# from binary_file_io
export write_bin,read_bin
# from post_proccesing_data_function
export collect_result_data
# from general_solver (v0.2.0)
export solve_tise,TISESolution,sweep_eigenvalues,kinetic_coefficients,hamiltonian_forms,discretization
# from tabulated_potentials (v0.2.0)
export tabulated_potential,TabulatedPotential1D,TabulatedPotential2D
# from autocorrelation (v0.2.0)
export expansion_coefficients,autocorrelation,autocorrelation_spectrum,fourier_transform
# from pauli_fierz (v0.2.0)
export pauli_fierz_potential
# from plotting (v0.2.0, implemented in the Makie extension)
export plot_energy_curves,plot_energy_map

end
