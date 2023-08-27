using Documenter
using TimeIndependentSchrodingerEquation

makedocs(
    sitename = "TimeIndependentSchrodingerEquation.jl",
    modules  = [TimeIndependentSchrodingerEquation],
    pages = [
        "Home" => "index.md",
        "Guide Information"=>"guide_information",
        "Function Information" => "function_information.md",
        "Simulation Example" => "simulation_example.md",
        "Contact Information" => "contact_information.md"
        ],
    format = Documenter.HTML(;
    prettyurls = get(ENV, "CI", "false") == "true",
    assets = String[],),
    )

deploydocs(;
    repo="github.com/mendzmartin/TimeIndependentSchrodingerEquation.jl.git",
    )