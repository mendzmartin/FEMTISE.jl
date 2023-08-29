# How to run default eigen problems defined inside TimeIndependentSchrodingerEquation Package

## Coling TimeIndependentSchrodingerEquation package repository
First of all you need to clone package from GitHub repository as follow
```bash
    @prompt$: cd ~/my_directory/
    @my_directory$: git clone https://github.com/mendzmartin/TimeIndependentSchrodingerEquation.jl.git
```

This will download a folder called `TimeIndependentSchrodingerEquation.jl`, it is important to keep the `.jl` extension in the repository name. And, in case we have already cloned the repository, we must update it by running `git pull`.

## Simulate predefined potential
Create a folder where you want to save simuation data, from Bash terminal write
```bash
    @prompt$: mkdir ~/my_folder
    @prompt$: cd ~/my_folder
    @my_folder$: julia
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.9.0 (2023-05-07)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |
```

```julia
    julia> Ctrl+]
    (@v1.9) pkg>
    (@v1.9) pkg> activate .
    (@v1.9) pkg> dev ~/my_directory/TimeIndependentSchrodingerEquation.jl
    (@v1.9) pkg> add Revise
    julia> exit()
```

Then, open a Julia file `@my_folder$: vi my_script.jl` and write the following comands:
```julia
    begin
        using Pkg
        Pkg.activate("./")
        Pkg.instantiate()
        
        using Revise
        using TimeIndependentSchrodingerEquation;
        
        run_default_eigen_problem(set_type_potential())
    end
```
After save the changes you can run the script and see de following
```bash
    @my_folder$: julia my_script.jl
  Activating project at `~/test_default_eigen_problem_from_input`
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
The types of default potential can be:
   - Unidimensional Quantum Harmonic Oscillator            --> set (1)
   - Unidimensional Symmetric Finit Kronig-Penney          --> set (2)
   - Unidimensional Finit Well Potential                   --> set (3)
   - Bidimensional Isotropic Quantum Harmonic Oscillator   --> set (4)
   - Ad hoc potential                                      --> set (5)
   - Ad hoc potential from input file                      --> set (6)
Please, set some number to specify the type potential:
```

Then you can specify the set number to run a specific potential. Be carefull when you set the numbers 5 or 6 because you need to do a little more steps before to run over this simulations (see the next).


## Simulate custom potential

### Definitons of specific potentials
First you need to create a folder called `adhoc_potential` like this
```bash
    @my_folder$: mkdir ./adhoc_potential
    @my_folder$: cd ./adhoc_potential
```

Inside `adhoc_potentials/` folder we need to create a Julia file with name `my_julia_file.jl` where we can find custom potential functions in specific format.

Inside our `my_julia_file.jl` folder the format to write custom potential function need to be like follows:

**Unidimensional potential**
```
    export my_potential_1d
    function my_potential_1d(x,params::Tuple)
        λ₁,λ₂,λ₃...=params
        return f(x[1],λ₁,λ₂,λ₃...)
    end
```

Here params is a Tuple with potential parameters (can be Integers, Floats, Complex, etc) and `f(x[1],λ₁,λ₂,λ₃...)` is an expresion(function) of `x[1]` (unidmensional DOF) and `λᵢ`'s (parameters).

**Bidimensional potential**
```
    export my_potential_2d
    function my_potential_2d(x,params::Tuple)
        λ₁,λ₂,λ₃...=params
        return f(x[1],x[2],λ₁,λ₂,λ₃...)
    end
```

Here params is a Tuple with potential parameters (can be Integers, Floats, Complex, etc) and `f(x[1],x[2],λ₁,λ₂,λ₃...)` is an expresion(function) of `x[1]` (unidmensional DOF), `x[2]` (unidmensional DOF) and `λᵢ`'s (parameters).

### Building input data for custom potential simulations

In some specific folder we need to create a data folder `@my_folder$: vi my_input.dat` with custom potential input information behing the following format (the data next to equal sign is only for example information)

```dat
full_path_name              = ./qho2D_different_masses/results/qho2D_different_masses
L                           = 30.0
dom_type                    = s
nev                         = 3
dimension                   = 2D
sigma                       = 0.0
adhoc_file_name             = mendez_adhoc_potentials
potential_function_name     = qho2D_different_masses
params_potential_types      = f f f f
params_potential            = 1.0 10.0 0.0 -1.0
analysis_param              = false
Δx                          = 
nx                          = 100
ny                          = 100
different_masses            = 10.0

# #################################################################################################
    full_path_name::String:           Full path name where you want to write problem results
    L::Real:                          Finite element domain length [au]
    dom_type::String                  Domain type (symetric (s) or non-symetric (ns) domain)
    nev::Integer:                     Number of eigenvalues
    dimension::String:                Dimension of eigen value problem
    sigma::Real:                      Level shift used in inverse iteration [au]
    adhoc_file_name::String:          Julia file name with ad hoc potential
    potential_function_name::String:  Name of ad hoc potential function
    params_potential:                 Parameters of ad hoc potential function
    params_potential_types:           Paremters type -> f (Float), i (Integer), c (Complex)
    analysis_param
        only if want to sweap a parameter from params_potential
            analysis_param = λindex::Integer λi λf Δλ
        else
            analysis_param = false::Boolean
    only if dimension == 1D
        Δx::Real:                     Finite element size [au]
    only if dimension == 2D
        nx::Integer:                  Number of finite element of x direction
        ny::Integer:                  Number of finite element of y direction
    different_masses
        only if want to simulate 2D eigenproblem and two particles with different masses
            different_masses = mass2::Real
        only if want to simulate 2D eigenproblem or two particles with equal masses
            different_masses = false::Boolean
# #################################################################################################

# #################################################################################################
    WARNIG! Beware of whitespace between data values and do not use spaces in variables named
    full_path_name, adhoc_file_name and potential_function_name.
# #################################################################################################
```