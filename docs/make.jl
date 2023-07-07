using Pkg
Pkg.activate("../")

push!(LOAD_PATH,"../src/")

using TimeIndependentSchrodingerEquation
using Documenter

makedocs(
         sitename = "TimeIndependentSchrodingerEquation.jl",
         modules  = [TimeIndependentSchrodingerEquation],
         pages=[
                "Home" => "index.md"
               ])

deploydocs(;
    repo="github.com/mendzmartin/TimeIndependentSchrodingerEquation.jl",
)

# using Pkg
# Pkg.activate("../")

# using Documenter
# using TimeIndependentSchrodingerEquation
# makedocs(sitename="My Documentation")