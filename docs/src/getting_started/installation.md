# Installation

`FEMTISE.jl` requires **Julia 1.10 or newer**. It is not yet registered in the
General registry, so it is installed directly from GitHub.

## Install a released version (recommended)

Pinning a release guarantees that your scripts keep producing the same results when
the package evolves. In the Julia REPL press `]` to enter the package manager, then:

```julia
(@v1.10) pkg> activate .                                          # project of your calculation
(my_project) pkg> add https://github.com/mendzmartin/FEMTISE.jl#v0.1.1
```

The `#v0.1.1` suffix selects the release tag. The exact version is recorded in the
`Manifest.toml` of your project; keep that file next to your scripts. The list of
releases and what changed in each one is in the [Changelog](@ref).

Equivalently, from a script:

```julia
using Pkg
Pkg.activate(".")
Pkg.add(url = "https://github.com/mendzmartin/FEMTISE.jl", rev = "v0.1.1")
```

## Install the development version

To follow the `main` branch (latest changes, possibly unreleased):

```julia
(my_project) pkg> add https://github.com/mendzmartin/FEMTISE.jl
```

`pkg> update` then fetches the newest commit. Read the *Unreleased* section of the
[Changelog](@ref) before updating.

## Work on the package itself

If you want to modify `FEMTISE.jl`, clone it and add it in development mode:

```bash
git clone https://github.com/mendzmartin/FEMTISE.jl.git ~/src/FEMTISE.jl
```

```julia
(my_project) pkg> dev ~/src/FEMTISE.jl
```

Changes to the source are picked up the next time `using FEMTISE` is executed. The
test suite runs with `pkg> test FEMTISE`.

## Running scripts and using threads

Save your calculation in a file, e.g. `my_calculation.jl`:

```julia
using Pkg
Pkg.activate(@__DIR__)      # the project where FEMTISE was added
using FEMTISE
# ... your calculation ...
```

and run it with

```bash
julia my_calculation.jl          # serial
julia -t 4 my_calculation.jl     # 4 threads
```

Threads are used by the post-processing functions and by the conversion of finite
element functions to arrays. The eigensolver itself is serial except for the
multithreaded BLAS.

## Checking the installation

```julia
julia> using FEMTISE
julia> ?eigen_values_and_eigen_vectors     # help mode shows the docstring
```

Then follow the [Quick start](@ref).
