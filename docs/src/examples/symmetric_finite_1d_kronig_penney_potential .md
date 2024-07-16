# Tutorial To Simulate Symmetric Finite One Dimensional Kronig-Penney Potential

### Create simulation directory
First of all we need to create a specific directory to save this specific simulation results 

```bash
@prompt:~$ mkdir ~/my_directory_path/SFKP1D
@prompt:~$ cd ~/my_directory_path/SFKP1D
```

### Create function potential

We need to create a aspecific function potential for Kronig-Penney potential as following:

```bash
@prompt:~/my_directory_path/SFKP1D$ vi adhoc_potential_function.jl
```

Inside `adhoc_potential_function` write the following:

```julia
"""
    kronig_penney_sturm_liouville(params)

# Aim
    - Compute Kronig-Penney potential as Sturm-Liouville problem

# Arguments
    - `params::Tuple`: tuple with parameters
"""
function kronig_penney_sturm_liouville(params::Tuple)
    num_ions,a,b,V₀=params
    ħ=1.0;m=1.0;
    p(x) = 0.5*(ħ*ħ)*(1.0/m);
    q(x) = symetric_kronig_penney(x[1],num_ions,a,b,V₀)
    r(x) = 1.0;
    return p,q,r;
end

function heaviside(x)
    return 0.5*(sign(x)+1)==true
 end

function sym_rect_pot_barr(x,b::Real,V₀::Real)
   return V₀*(heaviside(x+0.5*b)-heaviside(x-0.5*b))
end

function kronig_penney_center(x,b::Real,V₀::Real)
    return sym_rect_pot_barr.(x,b,V₀)
end

function kronig_penney_left(x,num_ions::Integer,a::Real,b::Real,V₀::Real)
    result=0.0
    for i in 1:num_ions
        result = result .+ sym_rect_pot_barr.(x.+i*a,b,V₀)
    end
    return result
end

function kronig_penney_right(x,num_ions::Integer,a::Real,b::Real,V₀::Real)
    return kronig_penney_left(-x,num_ions,a,b,V₀)
end

export symetric_kronig_penney
"""
    symetric_kronig_penney(x,num_ions,a,b,V₀)

# Aim
    - Compute symetric Kronig-Penney potential

# Arguments
    - `x::Real`: input value
    - `num_ions::Integer`: number of ions
    - `a::Real`: distance between ions
    - `b::Real`: width of barrier
    - `V₀::Real`: height of barrier
"""
function symetric_kronig_penney(x,num_ions::Integer,a::Real,b::Real,V₀::Real)
    if (mod(num_ions,2) == 0)
        error("num_ions keyword need to be odd")
        stop()
    elseif (num_ions==1)
        kp = kronig_penney_center(x,b,V₀)
    else
        kp = kronig_penney_center(x,b,V₀) .+ kronig_penney_left(x,convert(Int,(num_ions-1)/2),a,b,V₀) .+ kronig_penney_right(x,convert(Int,(num_ions-1)/2),a,b,V₀)
    end
    return kp
end
```

Then using Jupyter Notebook (by intermediate Visual Studio Code) we can analyse output file so:

```bash
@prompt:~/my_directory_path/SFKP1D$ code SFKP1D.ipynb
```
Inside `SFKP1D.ipynb` we need to write the following:

### Environment Activation

Activate the specific environment for this simulation, note that within the `activate` function we must place the path where we want to locate this environment.

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

#### Develop Package

In case it is necessary, we must add the FEMTISE package to the environment with the local location.

```julia
develop_package = false
path_repo="~/my_path_repo/"
develop_package ? Pkg.develop(path=path_repo*"FEMTISE.jl") : nothing
```

#### Adding Packages

We install (if necessary) and use specific packages for the simulation.

```julia
install_pkg = false
if install_pkg
    Pkg.add("Gridap")
    Pkg.add("Plots")
end
using FEMTISE;
using Gridap;
using Plots;
```

### Miscellaneous functions

We include the potential function defined in specific Julia file:

```julia
include("~/my_directory_path/SFKP1D/adhoc_potential_function.jl")
```

### Potential parameters

We define the properties of the potential

```julia
grid_size_length=276;
potential_depth=-0.5;
distance_between_wells=11;
well_width=1;
num_ions=23;
space_discretization=0.01;

unit_cell_potential=distance_between_wells+well_width;
```

#### Checking representation

Taking into account the finite size of the FE grid and the dimensions of the potential wells (width and separation), we can check if the number of sites can be represented correctly.

```julia
quantity_check = num_ions*unit_cell_potential;
(grid_size_length ≥ quantity_check) ? println("The value of grid size length is ok (≥ $(quantity_check)).") : println("Increase grid size length, must be grid_size_length ≥ $(quantity_check).")
```

### Grid Building

We create the one-dimensional FE grid.

```julia
grid_type="simple_line";
params_model=("./","model1D",(-0.5*grid_size_length,0.5*grid_size_length),space_discretization);
model1D=make_model(grid_type,params_model);
rm(params_model[1]*params_model[2]*".msh")
```

The last step could be omitted if you want to save the grid in `.msh` format for external visualization.

#### Grid points

We construct the point vectors (grid evaluation points) and coordinate vectors.

```julia
point_number=round(Int,abs(grid_size_length/space_discretization)+1)
space_coordinate,points=space_coord((-0.5*grid_size_length,0.5*grid_size_length),space_discretization,point_number-1;dimension="1D")
```

#### Plotting potential function

```julia
fig = plot(space_coordinate,symetric_kronig_penney(space_coordinate,num_ions,unit_cell_potential,well_width,potential_depth),label="")
fig = plot!(fig,xlabel="space coordinate (x [au])",ylabel="potential depth (v [au])")
display(fig)
```

### Boundary conditions

We define the boundary conditions of the system, in our case we define homogeneous boundary conditions throughout the boundary.

```julia
BC_type="FullDirichlet";
FullDirichlet_values,FullDirichlet_tags=make_boundary_conditions(grid_type,BC_type,ComplexF64);
```

### FE Domains

We construct the FE domains of the grid: interior and boundary. Additionally, we construct the differentials of these domains.

```julia
interior_FE_domain,differential_interior_FE_domain,boundary_FE_domain,differential_boundary_FE_domain = measures(model1D,3,FullDirichlet_tags)
```

### FE Reference

We define the interpolation polynomials to be used and create the Test and Trial spaces associated with the weak formulations of the problem.

```julia
reff = ReferenceFE(lagrangian,Float64,2)
TestSpace,TrialSpace = fe_spaces(model1D,reff,grid_type;BC_type=BC_type,TypeData=ComplexF64)
```

### Sturm-Liouville formulation

We define the functions to use the Sturm-Liouville type formulation.

```julia
p,q,r = kronig_penney_sturm_liouville((num_ions,unit_cell_potential,well_width,potential_depth))
```

### Eigen value problem

We solve the eigenvalue problem

```julia
eigen_energies,eigen_states = eigen_values_and_eigen_vectors(p,q,r,differential_interior_FE_domain,TrialSpace,TestSpace;
params=(10,1e-9,500,:none,potential_depth))
```

### Plotting results

```julia
fig=plot()
probability_densities=FEMTISE.density(eigen_states)
for i in 1:3#eachindex(ϕ)
    fig=plot!(space_coordinate,probability_densities[i].(points),
    label="probability density of energy=$(round(real(eigen_energies[i]),digits=4))",
    legend=:bottomleft)
end

fig=plot!(xlabel="space coordinate (x [au])",ylabel=" ")
fig=plot!(space_coordinate,0.1 .* symetric_kronig_penney(space_coordinate,num_ions,unit_cell_potential,well_width,potential_depth),
label="0.1*Kronig-Penney potential [au]")

display(fig)
```

We can save de figure using `savefig(fig,"example_kronig_penney.pdf")`.