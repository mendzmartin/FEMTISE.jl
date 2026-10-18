# API reference

Functions and types are grouped by the stage of the calculation in which they are
used. The same docstrings are available in the REPL help mode (`?name`).

## Meshes and boundary conditions

```@autodocs
Modules = [FEMTISE]
Pages   = ["functions/mesh_generator_function.jl", "functions/boundary_conditions_function.jl"]
```

## Measures, spaces and forms

```@autodocs
Modules = [FEMTISE]
Pages   = ["functions/miscellaneous_functions.jl"]
```

## Eigenvalue problem

```@autodocs
Modules = [FEMTISE]
Pages   = ["functions/eigen_problem_definition_function.jl", "functions/eigen_problem_solve_function.jl"]
```

## Post-processing

```@autodocs
Modules = [FEMTISE]
Pages   = ["useful_functions_for_FEM_objects.jl"]
```

## Input-file driver

```@autodocs
Modules = [FEMTISE]
Pages   = ["functions/default_eigen_problem.jl", "functions/default_running_input_data.jl",
           "functions/default_running_miscellaneous_functions.jl", "functions/default_running_structures.jl"]
```

## Output and reloading

```@autodocs
Modules = [FEMTISE]
Pages   = ["functions/binary_file_io.jl", "functions/post_proccesing_data_function.jl",
           "functions/post_proccesing_data_structures.jl"]
```

## General solver (v0.2.0)

```@autodocs
Modules = [FEMTISE]
Pages   = ["functions/general_solver.jl"]
```

## Tabulated potentials (v0.2.0)

```@autodocs
Modules = [FEMTISE]
Pages   = ["functions/tabulated_potentials.jl"]
```

## Autocorrelation and spectra (v0.2.0)

```@autodocs
Modules = [FEMTISE]
Pages   = ["functions/autocorrelation.jl"]
```

## Pauli–Fierz Hamiltonian (v0.2.0)

```@autodocs
Modules = [FEMTISE]
Pages   = ["functions/pauli_fierz.jl"]
```

## Extended input-file driver (v0.2.0)

```@autodocs
Modules = [FEMTISE]
Pages   = ["functions/extended_driver.jl"]
```

## Plotting (v0.2.0, requires a Makie backend)

```@autodocs
Modules = [FEMTISE]
Pages   = ["functions/plotting.jl"]
```
