module SteadyStateSchrodingerEquation

using Gridap
using GridapGmsh
using Gmsh
using Gridap.CellData  # para construir condición inicial interpolando una función conocida
using Gridap.FESpaces  # para crear matrices afines a partir de formas bilineales
using Gridap.Algebra   # para utilizar operaciones algebraicas con Gridap
using FileIO
using LinearAlgebra
using SparseArrays
using SuiteSparse
using Arpack

include("Functions/MeshGeneratorFunction.jl")
include("Functions/BoundaryConditionsFunction.jl")
include("Functions/MiscellaneousFunctions.jl")
include("Functions/EigenProblemDefinitionFunction.jl")
include("Functions/EigenProblemSolveFunction.jl")

end
