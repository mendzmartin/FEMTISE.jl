using Pkg
Pkg.activate("../")

push!(LOAD_PATH,"../src/")

# using TimeIndependentSchrodingerEquation
# using Documenter

# # makedocs(
# #          sitename = "TimeIndependentSchrodingerEquation.jl",
# #          modules  = [TimeIndependentSchrodingerEquation],
# #          pages=[
# #                 "Home" => "index.md"
# #                ])

# makedocs(
#     sitename = "TimeIndependentSchrodingerEquation.jl",
#     modules  = [TimeIndependentSchrodingerEquation],
#     pages = ["Home" => "index.md"],
#     format = Documenter.HTML(;
#     prettyurls = get(ENV, "CI", "false") == "true",
#     assets = String[],),
#     )

# deploydocs(;
#     repo="github.com/mendzmartin/TimeIndependentSchrodingerEquation.jl.git",
#     )

using Documenter
using TimeIndependentSchrodingerEquation

makedocs(
    sitename = "TimeIndependentSchrodingerEquation",
    format = Documenter.HTML(),
    modules = [TimeIndependentSchrodingerEquation]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/mendzmartin/TimeIndependentSchrodingerEquation.jl.git"
)