# Build the documentation of FEMTISE.jl
#
#   julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
#   julia --project=docs docs/make.jl
#
# The changelog of the repository root is copied into the documentation so that
# both always agree.

using Documenter
using FEMTISE

const DOCS_SRC = joinpath(@__DIR__, "src")
cp(joinpath(@__DIR__, "..", "CHANGELOG.md"), joinpath(DOCS_SRC, "changelog.md"); force = true)

makedocs(
    sitename = "FEMTISE.jl",
    authors  = "Martín Mendez and Federico M. Pont",
    modules  = [FEMTISE],
    repo     = Remotes.GitHub("mendzmartin", "FEMTISE.jl"),
    format   = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical  = "https://mendzmartin.github.io/FEMTISE.jl/stable/",
        edit_link  = "main",
        assets     = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Getting started" => [
            "Installation" => "getting_started/installation.md",
            "Quick start"  => "getting_started/quickstart.md",
        ],
        "Manual" => [
            "Theory"                             => "manual/theory.md",
            "Architecture"                       => "manual/architecture.md",
            "Input file reference"               => "manual/input_file.md",
            "Low-level interface"                => "manual/low_level.md",
            "Output files and reloading results" => "manual/outputs.md",
            "Post-processing"                    => "manual/postprocessing.md",
            "Accuracy and best practices"        => "manual/accuracy.md",
        ],
        "Tutorials" => [
            "Built-in potentials (console)"          => "examples/simulation_example.md",
            "Finite Kronig–Penney chain"             => "examples/symmetric_finite_1d_kronig_penney_potential.md",
            "1D harmonic oscillator"                 => "examples/isotropic_1d_harmonic_oscillator_potential.md",
            "2D harmonic oscillator"                 => "examples/isotropic_2d_harmonic_oscillator_potential.md",
            "Two electrons: helium model"            => "examples/coulomb_interaction_2d_potential.md",
        ],
        "API reference" => "api.md",
        "Changelog"     => "changelog.md",
        "Citing and contact" => "citing.md",
    ],
    # Internal functions without docstrings and cross-references in the legacy
    # tutorials must not stop the build; they are reported as warnings.
    warnonly = [:missing_docs, :cross_references],
)

# `devbranch` publishes `main` as /dev; every tag vX.Y.Z is published as /vX.Y.Z
# and the latest one as /stable.
deploydocs(;
    repo      = "github.com/mendzmartin/FEMTISE.jl.git",
    devbranch = "main",
    push_preview = false,
)
