module SteadyStateSchrodingerEquation

# Write your package code here.

include("Functions/MeshGeneratorFunction.jl")
include("Functions/BoundaryConditionsFunction.jl")
include("Functions/MiscellaneousFunctions.jl")
include("Functions/EigenProblemDefinitionFunction.jl")
include("Functions/EigenProblemSolveFunction.jl")



end
