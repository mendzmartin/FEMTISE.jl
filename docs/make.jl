using Documenter
using FEMTISE

makedocs(
    sitename = "FEMTISE.jl",
    modules  = [FEMTISE],
    pages = [
        "Home" => "index.md",
        "Guide Information"=>"guide_information.md",
        "Function Information" => "function_information.md",
        "Simulation Example" => "simulation_example.md",
        "Contact Information" => "contact_information.md"
        ],
    format = Documenter.HTML(;
    prettyurls = get(ENV, "CI", "false") == "true",
    assets = String[],),
    )

deploydocs(;
    repo="github.com/mendzmartin/FEMTISE.jl.git",
    )