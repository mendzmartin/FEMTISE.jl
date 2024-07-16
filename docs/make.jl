using Documenter
using FEMTISE

makedocs(
    sitename = "FEMTISE.jl",
    modules  = [FEMTISE],
    pages = [
        "Home" => "index.md",
        "Guide Information"=>"guide_information.md",
        "Function Information" => "function_information.md",
        "Simulation Example" => [
            "General Simulation Examples" => "./examples/simulation_example.md",
            "Symmetric Finite 1D Kronig-Penney Potential" => "./examples/symmetric_finite_1d_kronig_penney_potential .md",
            "Isotropic 1D Quantum Harmonic Oscillator Potential" => "./examples/isotropic_1D_harmonic_oscillator_potential.md",
            "Isotropic 2D Quantum Harmonic Oscillator Potential" => "./examples/isotropic_2D_harmonic_oscillator_potential.md",
            "Two-Electrons Coulomb Interaction Potential" => "./examples/coulomb_interaction_2d_potential.md"
        ],
        "Contact Information" => "contact_information.md"
        ],
    format = Documenter.HTML(;
    prettyurls = get(ENV, "CI", "false") == "true",
    assets = String[],),
    )

deploydocs(;
    repo="github.com/mendzmartin/FEMTISE.jl.git",
    )