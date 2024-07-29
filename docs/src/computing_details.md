# Comparison between ARPACK, LAPACK and BLAS packages

| **Library** | **Purpose** | **Usage** | **Performance** | **Dependence** |
|-------------|-------------|-----------|-----------------|----------------|
| **ARPACK**  | Solves large-scale eigenvalue problems for sparse matrices. | Suitable for large, sparse matrices. | Optimized for sparse matrices, efficient for large-scale problems. | Depends on LAPACK and BLAS for some operations[1][2]. |
| **LAPACK**  | Provides routines for solving systems of linear equations, least squares problems, and eigenvalue problems for dense and banded matrices. | Suitable for dense matrices or sparse matrices with trivial structure. | Optimized for dense matrices, efficient for high-performance computing. | Often used with BLAS for low-level operations[1]. |
| **BLAS**    | Provides low-level routines for performing basic linear algebra operations such as vector and matrix operations. | Underlies LAPACK and other higher-level linear algebra libraries. | Optimized for low-level operations, highly efficient for matrix multiplications and vector operations. | Often used as a building block for other libraries like LAPACK and ATLAS[3]. |

### Additional Notes:
- **ARPACK**: Designed for solving large-scale eigenvalue problems efficiently, especially for sparse matrices. It depends on LAPACK and BLAS for some operations[1][2].
- **LAPACK**: A comprehensive library for dense and banded matrices, providing routines for matrix factorizations and solving linear systems. It is often used with BLAS for optimized performance[1][3].
- **BLAS**: The foundation for many linear algebra libraries, providing low-level operations like matrix multiplication and vector operations. It is highly optimized for performance and is used as a building block for other libraries[3].

This table highlights the key differences in the purpose, usage, and performance characteristics of these libraries, which are essential for scientific computing and high-performance computing applications.

Citations:

+ [[1] Lehoucq, R.B., Sorensen, D.C. and Yang, C., 1998. ARPACK users' guide: solution of large-scale eigenvalue problems with implicitly restarted Arnoldi methods. Society for Industrial and Applied Mathematics.](http://li.mit.edu/Archive/Activities/Archive/CourseWork/Ju_Li/MITCourses/18.335/Doc/ARPACK/Lehoucq97.pdf)
+ [[2] http://rsuib.cc.rsu.ru/libraries/ARPACK/node1.html](http://rsuib.cc.rsu.ru/libraries/ARPACK/node1.html)
+ [[3] https://stackoverflow.com/questions/17858104/what-is-the-relation-between-blas-lapack-and-atlas](https://stackoverflow.com/questions/17858104/what-is-the-relation-between-blas-lapack-and-atlas)

### Implementation of FEMTISE package

[Here](https://github.com/mendzmartin/FEMTISE.jl/blob/c4c72d603e9e8516f08a37f966d3ee3b91e7f719/src/functions/eigen_problem_definition_function.jl#L67-L79) you can see the specific implementation of solve function to resolve eigen value problems.

# Comparison between FDM, FEM and DVR

In summary, FDM (Finite Difference Method) is the simplest to implement but limited to regular grids, FEM (Finite Element Method) is the most flexible and accurate but requires more sophisticated mathematics, and DVR (Discrete Variable Representation) is tailored for quantum mechanical applications. The choice depends on the specific problem, geometry, and desired accuracy.

| Property | FDM | FEM | DVR |
|----------|-----|-----|-----|
| Type | grid method | grid method | pseudo spectral method |
| Discretization | Divides the domain into a grid of points and approximates derivatives using differences between adjacent points | Divides the domain into small elements of simple geometric shapes and uses basis functions to approximate the solution within each element | Represents the wavefunction on a discrete grid and uses finite difference approximations for derivatives. Use a diagonal potential energy matrix, allowing for efficient computation of quantum systems. |
| Geometry | Works best on regular, rectangular domains | Requires a mesh of elements (triangles, quadrilaterals, tetrahedra, etc.) to discretize the domain. Can handle complex, arbitrarily shaped domains by using unstructured meshes | Can use both regular and irregular grids. Does not require a mesh in the traditional sense; it uses a grid of points to represent the wavefunction and potential energy. |
| Formulation | Directly discretizes the differential equations | Uses a variational formulation based on minimizing a functional | Discretizes the Schrödinger equation using finite differences |
| Computing objects | use structured and sparse matrices | use sparse matrices and need to compute potential matrix which is multidimensional intensive | potential matrix is diagonal and kinetic matrix is a full but has analytical expression |
| Accuracy | Typically lower order of accuracy than FEM | Can achieve higher accuracy by using higher-order polynomial approximations | Accuracy depends on the grid spacing and order of finite difference approximations used |
| Conservation | Does not inherently conserve quantities like mass, momentum, energy | Conserves quantities like mass and momentum if the basis functions satisfy the conservation laws | Inherently conserves probability in quantum mechanics applications |
| Applications | Commonly used for heat transfer, fluid flow, and electromagnetics problems | Widely used for structural analysis, heat transfer, fluid dynamics, and electromagnetics | Primarily used for solving the Schrödinger equation in quantum mechanics and quantum chemistry |