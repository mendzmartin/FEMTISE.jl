# Tutorial To Simulate Isotropic Two Dimensional Harmonic Oscillator

### Create simulation directory
First of all we need to create a specific directory to save this specific simulation results 

```bash
@prompt:~$ mkdir ~/my_directory_path/QHO2D
@prompt:~$ cd ~/my_directory_path/QHO2D
```

### Create function potential

We need to create a aspecific function potential for quantum harmonic oscillator 2D as following:

```bash
@prompt:~/my_directory_path/QHO2D$ vi adhoc_potential_function.jl
```

Inside `adhoc_potential_function` write the following:

```julia
export qho_2d
"""
    qho_2d(x,params)

# Aim:
    This function is a simple implementation of the 2D quantum harmonic oscillator potential. 
    It is used to test the simulation of the isotropic quantum harmonic oscillator in 2D.

# Arguments
    x::Array{Float64,1} : The position of the particle in 2D.
    params::Tuple : A tuple containing the parameters of the potential. 
        params[1]::Float64 : The frequency of the oscillator.
        params[2]::Float64 : The position of the minimum of the potential in the x direction.
        params[3]::Float64 : The position of the minimum of the potential in the y direction.
"""
function qho_2d(x,params::Tuple)
    ω,x₁,y₁=params
    return 0.5*(ω*ω)*((x[1]-x₁)*(x[1]-x₁)+(x[2]-y₁)*(x[2]-y₁))
end
```

### Input file

We need to create an input file to simulate using default solver function inside FEMTISE package

```bash
@prompt:~/my_directory_path/QHO2D$ vi input.dat
```
Inside `input.dat` we need to write the following.

```text
full_path_name              = ~/my_directory_path/QHO2D/name_output_file
dom_type                    = s
nev                         = 10
dimension                   = 2D
sigma                       = 0.0
adhoc_file_name             = ~/my_directory_path/QHO2D/adhoc_potential_function
potential_function_name     = qho_2d
params_potential_types      = f f f
params_potential            = 1.0 0.0 0.0
analysis_param              = false
output_format_type          = jld2 eigen
## ONLY FOR 1D EIGENPROBLEMS
L                           = 
Δx                          = 
## ONLY FOR 2D EIGENPROBLEMS
Lx                          = 10
Ly                          = 10
nx                          = 100
ny                          = 100
different_masses            = false
reduced_density             = false
```
Note that `~/my_path_repo/` is the directory path where we can find FEMTISE.jl package.

### Run script

Create Julia code as
```bash
@prompt:~/my_directory_path/QHO2D$ vi run.jl
```
Inside `run.jl` we need to write the following.

```julia
begin
    using Pkg
    Pkg.activate("~/my_directory_path/QHO2D/")
    develop_package = true; develop_package ? Pkg.develop(path="~/my_path_repo/FEMTISE.jl") : nothing
    Pkg.instantiate()
    using FEMTISE;
    run_default_eigen_problem(set_type_potential("~/my_directory_path/QHO2D/input.dat"))
end
```

After this we can run the simulation using Julia compiler (for example: using multithread running with four threads)

```bash
@prompt:~/my_directory_path/QHO2D$ julia -t 4 run.jl 
```

### Analysis

After running we obtain an output data file in jld2 format called `name_output_file_eigen_data.jld2`.

Then using Jupyter Notebook (by intermediate Visual Studio Code) we can analyse output file so:

```bash
@prompt:~/my_directory_path/QHO2D$ code QHO2D.ipynb
```
Inside `QHO2D.ipynb` we need to write the following:

#### Environment

Activate Julia environment

```julia
using Pkg
Pkg.activate("~/my_directory_path/QHO2D/")
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
path_input_file_name="~/my_directory_path/QHO2D/input.dat"
simulation_info, output_data = collect_result_data(true,path_input_file_name)
```
#### Plotting figures

Defining functions to plot data we have:
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
    density2D(x,y,phi_n)

# Aim
- Calculate the probability density of the eigenstates in 2D.

# Arguments
- `x::Array{Float64,1}`: Array of x-coordinates.
- `y::Array{Float64,1}`: Array of y-coordinates.
- `phi_n::Array{Complex{Float64},1}`: Array of eigenstates.
"""
function density2D(x,y,phi_n)
    density=zeros(length(y),length(x))
    for i in eachindex(y)
        for j in eachindex(x)
            index=(j-1)*length(y)+i
            density[i,j]=real.(conj.(phi_n[index]).*phi_n[index])
        end
    end
    return density
end

"""
    plot_eigenstates(id,results,index_nev;<keyword arguments>)

# Aim
- Plot the eigenstates of the Hamiltonian operator.
The eigenstates are obtained from the diagonalization of the Hamiltonian operator.
The eigenstates are plotted for the energy level specified by the index_nev variable.
The keyword arguments are used to set the color map of the plot (only for 2D plot).

# Arguments
- `id`: InputData or InputData1D or InputData2D object.
- `results`: Results object.
- `index_nev::Int`: Energy level to plot.
- `keyword arguments`:
    - `mapcolor::Symbol=:rainbow1`: Color map of the plot (only for 2D plot).
"""
function plot_eigenstates(id,results,index_nev::Int;mapcolor::Symbol=:rainbow1)
    if id.analysis_param == false
        if id.params.dimension == "1D"
            plotlyjs();
            if ((typeof(id) <: InputData) || (typeof(id) <: InputData1D && id.output_format_type == ("bin","eigen")))
                rho = real.(conj.(results.ϕ[:,index_nev]).*(results.ϕ[:,index_nev]))
            elseif (typeof(id) <: InputData1D && id.output_format_type in [("jld2","eigen"),("jld2","all")])
                rho = real.(conj.(results.ϕ[index_nev].(results.pts)).*(results.ϕ[index_nev].(results.pts)))
            end
            figure = plot(results.r,rho,lw=2,lc=:black,label="")
            figure = scatter!(results.r,rho,label="n=$(index_nev)",lw=0.1)
            figure = plot!(xlabel="Coordinate (x [au])",ylabel="Probability density (ρn(x))",ticks = :native)
        elseif id.params.dimension == "2D"
            gr();
            if ((typeof(id) <: InputData) || (typeof(id) <: InputData2D && id.output_format_type == ("bin","eigen")))
                if id.params.nx==id.params.ny
                    rho = density2D(results.r[:,1],results.r[:,2],results.ϕ[:,index_nev])
                    figure1 = contour(results.r[:,1],results.r[:,2],rho,levels=10, color=mapcolor, fill=true,lw=0)
                else
                    rho = density2D(results.r[1:id.params.nx,1],results.r[1:id.params.ny,2],results.ϕ[:,index_nev])
                    figure1 = contour(results.r[1:id.params.nx,1],results.r[1:id.params.ny,2],rho, levels=10, color=mapcolor, fill=true,lw=0)
                end
            elseif (typeof(id) <: InputData2D && id.output_format_type in [("jld2","eigen"),("jld2","all")])
                rho = density2D(results.r[1],results.r[2],results.ϕ[index_nev].(results.pts))
                figure1 = contour(results.r[1],results.r[2],rho,levels=10, color=mapcolor, fill=true,lw=0)
            end
            figure1 = plot(figure1,title="Probability density (ρ$(index_nev)(x))",xlabel="Coordinate (x [au])", ylabel="Coordinate (y [au])")
    
            if id.reduced_density
                plotlyjs()
                if ((typeof(id) <: InputData) || (typeof(id) <: InputData1D && id.output_format_type == ("bin","eigen")))
                    if id.params.nx==id.params.ny
                        figure2=plot(results.r[:,1],results.rhoDOF1[:,index_nev],label="ρn(x)")
                        figure2=plot!(results.r[:,2],results.rhoDOF2[:,index_nev],label="ρn(y)")
                    else
                        figure2=plot(results.r[1:id.params.nx,1],results.rhoDOF1[:,index_nev],label="ρn(x)")
                        figure2=plot!(results.r[1:id.params.ny,2],results.rhoDOF2[:,index_nev],label="ρn(y)")
                    end
                elseif (typeof(id) <: InputData2D && id.output_format_type in [("jld2","eigen"),("jld2","all")])
                    figure2 = plot(results.r[1],results.rhoDOF1[:,index_nev],label="ρn(x)")
                    figure2 = plot!(results.r[2],results.rhoDOF2[:,index_nev],label="ρn(y)")
                end
                figure2=plot!(title="Reduced probability density for level n=$(index_nev)")
                figure2=plot!(xlabel="Coordinate (x or y [au])",ylabel="",ticks = :native)

                figure = plot(figure1,figure2,layout=2)
            else
                figure = figure1
            end
        end
    else
        println("PLOT ERROR.")
        println("You can not plot eigenstate with activated analysis params.")
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

Also, we can export figures as `*pdf` format using
```julia
savefig(fig1,"./eigen_energies.pdf")
```

and eigenfunctions:
```julia
eigenstate_to_show=2
fig2=plot_eigenstates(simulation_info, output_data,eigenstate_to_show;mapcolor=:turbo)
display(fig2)
```