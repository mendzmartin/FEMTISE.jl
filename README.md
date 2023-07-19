# **TimeIndependentSchrodingerEquation**

| **GitHub Actions - Workflows** |
|:------------ |
| [![Documentation](https://github.com/mendzmartin/TimeIndependentSchrodingerEquation.jl/actions/workflows/documentation.yml/badge.svg)](https://github.com/mendzmartin/TimeIndependentSchrodingerEquation.jl/actions/workflows/documentation.yml) |
| [![Build Status](https://github.com/mendzmartin/TimeIndependentSchrodingerEquation.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mendzmartin/TimeIndependentSchrodingerEquation.jl/actions/workflows/CI.yml?query=branch%3Amain) |
| [![CompatHelper](https://github.com/mendzmartin/TimeIndependentSchrodingerEquation.jl/actions/workflows/CompatHelper.yml/badge.svg)](https://github.com/mendzmartin/TimeIndependentSchrodingerEquation.jl/actions/workflows/CompatHelper.yml) |
| [![TagBot](https://github.com/mendzmartin/TimeIndependentSchrodingerEquation.jl/actions/workflows/TagBot.yml/badge.svg)](https://github.com/mendzmartin/TimeIndependentSchrodingerEquation.jl/actions/workflows/TagBot.yml) |

<!-- [![Codecov](https://app.codecov.io/gh/mendzmartin/TimeIndependentSchrodingerEquation.jl/branch/main/graph/badge.svg)](https://app.codecov.io/gh/mendzmartin/TimeIndependentSchrodingerEquation.jl) -->

Julia package repository to resolve Time-Independent Schrödinger equation by (F)init (E)lement (M)ethod. This is an implementation of Gridap package for unidimensional and bidimensional grids.

Project carried out in my PhD studies of Physics at:
* [Faculty of Mathematics, Astronomy, Physics and Computation (FaMAF)](https://www.famaf.unc.edu.ar/)
* The Enrique Gaviola Institute of Physics (IFEG)

[*Check here The GitHub Page*](https://mendzmartin.github.io/TimeIndependentSchrodingerEquation.jl/)

## **Warning!**
***->> The package is under construction <<-***

## **How we can use `TimeIndependentSchrodingerEquation` package ?**
### **1. Clone `TimeIndependentSchrodingerEquation` package**

First we need to clone package from GitHub repository as follow
```bash
    @prompt$: cd ~/my_directory/
    @my_directory$: git clone https://github.com/mendzmartin/TimeIndependentSchrodingerEquation.jl.git
```

This will download a folder called `TimeIndependentSchrodingerEquation.jl`, it is important to keep the `.jl` extension in the repository name. And, in case we have already cloned the repository, we must update it by running `git pull`.

### **2. Build Julia code to use `TimeIndependentSchrodingerEquation` package**
Build Julia code with structure:
```julia
module MyModule

using Pkg
Pkg.activate("./")
Pkg.instantiate()

Pkg.add(path="~/my_directory/TimeIndependentSchrodingerEquation.jl")
using TimeIndependentSchrodingerEquation

#=
... ...
    code block where we use function
    from TimeIndependentSchrodingerEquation package
... ....
=#

end
```
Then we save de Julia code below name `MyModule.jl`. After that, we can run the code doing this
```bash
julia MyModule.jl
```
Also, if we have activated multi-thread configuration we can use the next command to activate parallelism:
```bash
julia -t 4 MyModule.jl
```
where we have specified 4 threads as parallelization.

## **Examples**
See some use examples in [test folder](https://github.com/mendzmartin/TimeIndependentSchrodingerEquation.jl/tree/main/test) where you can find some uses of `TimeIndependentSchrodingerEquation` package's functions.

## **Contact**
Please, contact the project administrator [Méndez Martín](mailto:martinmendez@mi.unc.edu.ar) for any improve suggestion or questions about package use.