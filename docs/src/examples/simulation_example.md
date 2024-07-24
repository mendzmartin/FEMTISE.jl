# How to run default eigen problems defined inside FEMTISE Package

## Clone FEMTISE package repository
First of all you need to clone package from GitHub repository as follow
```bash
    @prompt$: cd ~/my_directory/
    @my_directory$: git clone https://github.com/mendzmartin/FEMTISE.jl.git
```

This will download a folder called `FEMTISE.jl`, it is important to keep the `.jl` extension in the repository name. And, in case we have already cloned the repository, we must update it by running `git pull`.

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
    (@v1.9) pkg> dev ~/my_directory/FEMTISE.jl
    julia> exit()
```

Then, open a Julia file `@my_folder$: vi my_script.jl` and write the following comands:
```julia
    begin
        using Pkg
        Pkg.activate("./")
        Pkg.instantiate()
        
        # if you did the previous step from REPL you don't need this
        Pkg.develop(path="~/my_directory/FEMTISE.jl")

        using FEMTISE;
        
        run_default_eigen_problem(set_type_potential())
    end
```
After save the changes you can run the script and see de following
```bash
    @my_folder$: julia my_script.jl
  Activating project at `~/my_folder`
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
####### ####### #     # #######   ###    #####  #######
#       #       ##   ##    #       #    #     # #                     #  #
#       #       # # # #    #       #    #       #                     #  #
#####   #####   #  #  #    #       #     #####  #####                 #  #
#       #       #     #    #       #          # #         ###         #  #
#       #       #     #    #       #    #     # #         ###    #    #  #
#       ####### #     #    #      ###    #####  #######   ###     ####   ######
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
The type of default potential can be:
   - Unidimensional Quantum Harmonic Oscillator            --> set (1)
   - Unidimensional Symmetric Finit Kronig-Penney          --> set (2)
   - Unidimensional Finite Well Potential                  --> set (3)
   - Bidimensional Isotropic Quantum Harmonic Oscillator   --> set (4)
   - Ad hoc potential from input file                      --> set (5)
Please, set some number to specify the type potential:
```

Then you can specify the set number to run a specific potential. Then the running is interactive so in each step you need to configurate the potential and simulation properties.

### 1D Harmonic Oscillator Example

Just to show how we can simulate the Unidimensional Quantum Hamonic Oscillator we need to specify befor the set number 1, then the properties to configurate are like following:

```bash
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Set full path name (e.g: "./my_directory/my_name") where you want to write problem results and press Enter = ./qho1d
Mandatory input data
Number of eigenvalues: nev::Int = 10
Finite element domain length [au] (default L=30.0 press Enter): L::Float64 = 
Set domain type (for symetric domain {-L/2,L/2} set "s" and for non-symetric domain {0,L} set "ns": dom_type::String = s
Finite element size [au] (default Δx=0.1 press Enter): Δx::Float64 = 
Harmonic Oscillator frecuency [au]: ω::Float64 = 1
Harmonic Oscillator center [au]: x₁::Float64 = 0
Level shift used in inverse iteration [au] (default sigma=0.0 press Enter): sigma::Float64 = 0
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
... Running ...
Building the grid model ...
Info    : Meshing 1D...
Info    : Meshing curve 1 (Line)
Info    : Done meshing 1D (Wall 0.000265359s, CPU 0.000267s)
Info    : 301 nodes 302 elements
Info    : Writing './model1D.msh'...
Info    : Done writing './model1D.msh'
Info    : Reading './model1D.msh'...
Info    : 3 entities
Info    : 301 nodes
Info    : 302 elements
Info    : Done reading './model1D.msh'
Solving eigen problem ...
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Transforming data from FE object to complex value...
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Moving result data to trash if data exists ...
Saving data ...
Saved data.
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Moving result data to trash if data exists ...
Saving data ...
Saved data.
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Moving result data to trash if data exists ...
Saving data ...
Saved data.
```

Note that if we don't specify explicit data we choose the default propertie. Finally the output are four files called:
```bash
    @my_folder$: ls *.dat
qho1d_eigen_problem_attributes.dat 
    @my_folder$: ls *.bin
qho1d_coordinates.bin  qho1d_eigen_values.bin  qho1d_eigen_vectors.bin
```
and we note that the results are in binary otuput data file and the attributes of simulation are inside data file (see bellow the file information).

```bash
    @my_folder$: cat qho1d_eigen_problem_attributes.dat
Quantum Harmonic Oscillator 1D
Dimension of eigen value problem              dimension::String   = 1D
Number of eigenvalues                         nev::Int          = 10
Finite element domain length [au]             L::Float64          = 30.0
Domain type                                   dom_type::String    = s
Level shift used in inverse iteration [au]    sigma::Float64      = 0.0
Finite element size [au]                      Δx::Float64         = 0.1
Harmonic Oscillator frecuency [au]            ω::Float64          = 1.0
Harmonic Oscillator center [au]               x₁::Float64         = 0.0

Number of threads                                                           = 4
´´´

Be carefull when you set the numbers 5 or 6 because you need to do a little more steps before to run over this simulations (see the next).

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
full_path_name              = ../my_folder_name/my_file_name
dom_type                    = s
nev                         = 10
dimension                   = 1D
sigma                       = 0.0
adhoc_file_name             = my_adhoc_potentials
potential_function_name     = my_potential
params_potential_types      = f f f
params_potential            = 1.0 0.1 10.0
analysis_param              = 2 0.0 0.1 0.01
## ONLY FOR 1D EIGENPROBLEMS
L                           = 100.0
Δx                          = 0.1
## ONLY FOR 2D EIGENPROBLEMS
Lx                          = 
Ly                          = 
nx                          = 
ny                          = 
different_masses            = 
reduced_density             = 


# #################################################################################################
    full_path_name::String:           Full path name where you want to write problem results
    dom_type::String                  Domain type (symetric (s) or non-symetric (ns) domain)
    nev::Integer:                     Number of eigenvalues
    dimension::String:                Dimension of eigen value problem
    sigma::Real:                      Level shift used in inverse iteration [au]
    adhoc_file_name::String:          Julia file name with ad hoc potential
    potential_function_name::String:  Name of ad hoc potential function
    params_potential_types:           Paremters type -> f (Float), i (Integer), c (Complex)
    params_potential:                 Parameters of ad hoc potential function
    analysis_param
        only if want to sweap a parameter from params_potential
            analysis_param = λindex::Integer λi λf Δλ
        else
            analysis_param = false
    only if dimension == 1D
        L::Real:                      Finite element domain length [au]
        Δx::Real:                     Finite element size [au]
    only if dimension == 2D
        Lx::Real:                     Finite element domain length of x direction [au]
        Ly::Real:                     Finite element domain length of y direction [au]
        nx::Integer:                  Number of finite element of x direction
        ny::Integer:                  Number of finite element of y direction
        different_masses
            only if want to simulate 2 particles with different masses
                different_masses::Real = DOF2mass
            else
                different_masses::Bool = false
        reduced_density
            if want to compute reduced densities
                reduced_density::Bool = true
            else
                reduced_density::Bool = false
# #################################################################################################

# #################################################################################################
    WARNIG! Beware of whitespace between data values and do not use spaces in variables named
    full_path_name, adhoc_file_name and potential_function_name.
```

Then you can run the script and set options `5` or `6` like this:
```bash
    @my_folder$: julia my_script.jl
  Activating project at `~/my_folder`
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
####### ####### #     # #######   ###    #####  #######
#       #       ##   ##    #       #    #     # #                     #  #
#       #       # # # #    #       #    #       #                     #  #
#####   #####   #  #  #    #       #     #####  #####                 #  #
#       #       #     #    #       #          # #         ###         #  #
#       #       #     #    #       #    #     # #         ###    #    #  #
#       ####### #     #    #      ###    #####  #######   ###     ####   ######
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
The type of default potential can be:
   - Unidimensional Quantum Harmonic Oscillator            --> set (1)
   - Unidimensional Symmetric Finit Kronig-Penney          --> set (2)
   - Unidimensional Finite Well Potential                  --> set (3)
   - Bidimensional Isotropic Quantum Harmonic Oscillator   --> set (4)
   - Ad hoc potential from input file                      --> set (5)
Please, set some number to specify the type potential: 6
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Set full path name (e.g: "./my_directory/my_input_data") where the data is specified and press Enter =
```

Where we were specified the full path name of our input data (only for option `6`). Also, we can specify this full path name inside Julia script like following:

First open script file `@my_folder$: vi my_script.jl` and write the following comands:
```julia
    begin
        using Pkg
        Pkg.activate("./")
        Pkg.instantiate()
        
        using FEMTISE;
        
        path_input_data = "./my_input"
        run_default_eigen_problem(set_type_potential(path_input_data))
    end
```
Then you can run the script `@my_folder$: julia my_script.jl` and wait for simulation results.