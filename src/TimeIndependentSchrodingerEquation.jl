module TimeIndependentSchrodingerEquation

using Gridap,GridapGmsh,Gmsh,FileIO
using LinearAlgebra,SparseArrays,Arpack
using DelimitedFiles,SpecialFunctions,KeywordDispatch

include("functions/mesh_generator_function.jl");
include("functions/boundary_conditions_function.jl")
include("functions/miscellaneous_functions.jl")
include("functions/eigen_problem_definition_function.jl")
include("functions/eigen_problem_solve_function.jl")
include("functions/default_running_miscellaneous_functions.jl")
include("functions/default_running_structures.jl")
include("functions/default_running_input_data.jl")
include("functions/useful_functions_for_FEM_objects.jl")
include("functions/binary_file_io.jl")
include("functions/default_eigen_problem.jl")

export make_model
export make_boundary_conditions
export measures,space_coord,bilineal_forms,fe_spaces
export default_solver_eigen_problem,ask_for_params,
        get_input,num_ions_BoolCondition,rm_existing_file
export eigen_problem,solve
export eigen_values_and_eigen_vectors
export input_data
export norm_l2,orthogonality_check,eigenstates_normalization
export write_bin,read_bin
export run_default_eigen_problem

end
