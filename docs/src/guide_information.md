## **How we can use `TimeIndependentSchrodingerEquation` package ?**
### **1. Clone `TimeIndependentSchrodingerEquation` package**

First we need to clone package from GitHub repository as follow
```bash
    @prompt$: cd ~/my_directory/
    @my_directory$: git clone https://github.com/mendzmartin/TimeIndependentSchrodingerEquation.jl.git
```

This will download a folder called `TimeIndependentSchrodingerEquation.jl`, it is important to keep the `.jl` extension in the repository name. And, in case we have already cloned the repository, we must update it by running `git pull`.

### **2. Build Julia code to use `TimeIndependentSchrodingerEquation` package**

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
    (@my_folder) pkg> instantiate()
    (@my_folder) pkg> dev ~/my_directory/TimeIndependentSchrodingerEquation.jl
    (@my_folder) pkg> add Revise
    julia> exit()
```
then build Julia code with following structure:

```julia
module MyModule

using Pkg
Pkg.activate("./")
Pkg.instantiate()

using TimeIndependentSchrodingerEquation
using Revise

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
    @my_directory$: julia MyModule.jl
```
Also, if we have activated multi-thread configuration we can use the next command to activate parallelism:
```bash
    @my_directory$: julia -t 4 MyModule.jl
```
where we have specified 4 threads as parallelization.

## Run from Julia REPL
We can also run the package directly from Julia REPL by opening the terminal `Ctrl+Alt+T` inside the package folder and typing the following commands inside the terminal:
```bash
    @prompt$: cd my_directory/TimeIndependentSchrodingerEquation.jl/
    @TimeIndependentSchrodingerEquation.jl$: julia
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
    (@my_folder) pkg> instantiate()
    (@my_folder) pkg> dev ~/my_directory/TimeIndependentSchrodingerEquation.jl
    (@my_folder) pkg> add Revise
    julia> using TimeIndependentSchrodingerEquation
    julia> using Revise
```
and then we can, for example, access Julia's help mode to ask for specific package functions such as the following:
```julia
    julia> Shift+?
    help?> space_coord
    search: space_coord

    if dimension=="1D" ⇒ dom=(x₁,x₂); Δr=Δx; n=nx
    if dimension=="2d" ⇒ dom=(x₁,x₂,y₁,y₂); Δr=(Δx,Δy); n=(nx,ny)
```