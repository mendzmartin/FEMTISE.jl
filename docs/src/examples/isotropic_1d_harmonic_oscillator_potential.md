# Tutorial To Simulate Isotropic One Dimensional Harmonic Oscillator

### Create simulation directory
First of all we need to create a specific directory to save this specific simulation results 

```bash
@prompt:~$ mkdir ~/my_directory_path/QHO1D
@prompt:~$ cd ~/my_directory_path/QHO1D
```

### Create function potential

We need to create a aspecific function potential for quantum harmonic oscillator 1D as following:

```bash
@prompt:~/my_directory_path/QHO1D$ vi adhoc_potential_function.jl
```

Inside `adhoc_potential_function` write the following:

```julia
export qho_1d
"""
    qho_1d(x,params)

# Aim: 
    This function is a simple implementation of the 1D quantum harmonic oscillator potential. 
    It is used to test the simulation of the isotropic quantum harmonic oscillator in 1D.

# Arguments
    x::Array{Float64,1} : The position of the particle in 1D.
    params::Tuple : A tuple containing the parameters of the potential. 
        params[1]::Float64 : The frequency of the oscillator.
        params[2]::Float64 : The position of the minimum of the potential.

# Returns
    Float64 : The value of the potential at the position x.
"""
function qho_1d(x,params::Tuple)
    ω,x₁=params
    return 0.5*(ω*ω)*((x[1]-x₁)*(x[1]-x₁))
end
```

### Input file

We need to create an input file to simulate using default solver function inside FEMTISE package

```bash
@prompt:~/my_directory_path/QHO1D$ vi input.dat
```
Inside `input.dat` we need to write the following.

```text
full_path_name              = ~/my_directory_path/QHO1D/name_output_file
dom_type                    = s
nev                         = 10
dimension                   = 1D
sigma                       = 0.0
adhoc_file_name             = ~/my_directory_path/QHO1D/adhoc_potential_function
potential_function_name     = qho_1d
params_potential_types      = f f
params_potential            = 1.0 0.0
analysis_param              = false
output_format_type          = jld2 all
## ONLY FOR 1D EIGENPROBLEMS
L                           = 100
Δx                          = 0.1
## ONLY FOR 2D EIGENPROBLEMS
Lx                          = 
Ly                          = 
nx                          = 
ny                          = 
different_masses            = 
reduced_density             = 
```

### Run script

Create Julia code as
```bash
@prompt:~/my_directory_path/QHO1D$ vi run.jl
```
Inside `run.jl` we need to write the following.

```julia
begin
    using Pkg
    Pkg.activate("~/my_directory_path/QHO1D/")
    develop_package = true; develop_package ? Pkg.develop(path="~/my_path_repo/FEMTISE.jl") : nothing
    Pkg.instantiate()
    using FEMTISE;
    run_default_eigen_problem(set_type_potential("~/my_directory_path/QHO1D/input.dat"))
end
```

After this we can run the simulation using Julia compiler (for example: using multithread running with four threads)

```bash
@prompt:~/my_directory_path/QHO1D$ julia -t 4 run.jl 
```

### Analysis

After running we obtain an output data file in jld2 format called `name_output_file_eigen_data.jld2`.

Then using Jupyter Notebook (by intermediate Visual Studio Code) we can analyse output file so:

```bash
@prompt:~/my_directory_path/QHO1D$ code QHO1D.ipynb
```
Inside `QHO1D.ipynb` we need to write the following:

#### Environment

Activate Julia environment

```julia
using Pkg
Pkg.activate("~/my_directory_path/QHO1D/")
Pkg.instantiate()
```
Is necessary to mark FEMTISE package as developed package using specific path repository:

```julia
develop_package = true; develop_package ? Pkg.develop(path="~/my_path_repo/FEMTISE.jl") : nothing
```

Now we install package (if is nesseary) and use specific packages to analyse output data:

```julia
install_pkg = true
if install_pkg
    Pkg.add("Plots")
    Pkg.add("PlotlyJS")
end
using FEMTISE;
using Plots;
```

#### Read output data

All the information that we need to specify is where we find input file then using specific functions we can collect output data

```julia
path_input_file_name="~/my_directory_path/QHO1D/input.dat"
simulation_info, output_data = collect_result_data(true,path_input_file_name)
```
#### Plotting figures

Defining functions to plot data as following:

```julia
"""
    plot_eigenvalues(id,results;<keyword arguments>)

# Aim
- Plot the eigenvalues of the Hamiltonian operator.
The eigenvalues are obtained from the diagonalization of the Hamiltonian operator.
The keyword arguments are used to set the title, xlabel, ylabel, and legend of the plot.

# Arguments
- `id`: InputData or InputData1D or InputData2D object.
- `results`: Results object.
- `keyword arguments`:
    - `set_title::String`: Title of the plot.
    - `set_xlabel::String`: Label of the x-axis.
    - `set_ylabel::String`: Label of the y-axis.
    - `set_legend::Symbol`: Position of the legend.
"""
function plot_eigenvalues(id,results;
    set_title::String="",
    set_xlabel::String="Energy level (n)",set_ylabel::String="Eigen-energies (ϵn [au])",set_legend::Symbol=:bottomright)
    if id.analysis_param == false
        plotlyjs()
        figure = scatter(real(results.ϵ),title=set_title,xlabel=set_xlabel,ylabel=set_ylabel,legend=set_legend)
    else
        println("PLOT ERROR.")
        println("Check attributes, you are using the wrong function method. Analysis parameter is activated.")
        figure = nothing
    end
    return figure
end

"""
    plot_eigenstates(id,results,range_to_show;<keyword arguments>)

# Aim
- Plot the eigenstates of the Hamiltonian operator.
The eigenstates are obtained from the diagonalization of the Hamiltonian operator.
The eigenstates are plotted for the range of energy levels specified by the range_to_show variable.
The keyword arguments are used to set the title, xlabel, ylabel, and legend of the plot.

# Arguments
- `id`: InputData or InputData1D or InputData2D object.
- `results`: Results object.
- `range_to_show::StepRange{Int, Int}`: Range of energy levels to plot.
- `keyword arguments`:
    - `set_xlabel::String`: Label of the x-axis.
    - `set_ylabel::String`: Label of the y-axis.
"""
function plot_eigenstates(id,results,range_to_show::StepRange{Int, Int};
    set_xlabel::String="Coordinate (x [au])",set_ylabel::String="Probability density (ρ(x))")
    if id.params.dimension == "1D"
        plotlyjs();
        figure = plot()
        rho=zeros(Float64,length(results.r))
        Threads.@threads for i in range_to_show
            if ((typeof(id) <: InputData) || (typeof(id) <: InputData1D && id.output_format_type == ("bin","eigen")))
                rho=real.(conj.((results.ϕ)[:,i]).*((results.ϕ)[:,i]))
            elseif (typeof(id) <: InputData1D && id.output_format_type in [("jld2","eigen"),("jld2","all")])
                rho = real.(conj.((results.ϕ[i]).(results.pts)).*((results.ϕ[i]).(results.pts)))
            end
            figure = plot!(results.r,rho,lw=2,lc=:black,label="")
            figure = scatter!(results.r,rho,label="n=$(i)",lw=0.1)
        end
        figure = plot!(xlabel=set_xlabel,ylabel=set_ylabel,ticks=:native)
    else
        println("PLOT ERROR.")
        println("Check attributes, you are using the wrong function method. 2D dimension problem is activated.")
        figure = nothing
    end
    return figure
end
```

Now we can plot eigenenergies:
```julia
fig1 = plot_eigenvalues(simulation_info, output_data)
display(fig1)
```

and eigenfunctions:
```julia
range_to_show=range(1,step=1,length=3)
fig2 = plot_eigenstates(simulation_info, output_data,range_to_show)
display(fig2)
```