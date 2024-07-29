# `FEMTISE.jl`: A High Performance Julia Package for Time-Independent Schrödinger Equation by Finite Element Method

> Méndez, Martín¹² & Pont, Federico M.¹²

¹Instituto de Física Enrique Gaviola (IFEG-FAMAF-CONICET)

²Facultad de Matemática, Astronomía, Física y Computación (FAMAF), Universidad Nacional de Córdoba (UNC), 5000, Córdoba, Argentina

# Abstract

`FEMTISE.jl` is a package to resolve Time-Independent Schrödinger Equation (TISE) by Finit Element Method (FEM). This is an implementation and extension over [GRIDAP.jl](https://github.com/gridap/Gridap.jl) package using high performance protocol. The package is under constructions and currently can solve one and two dimensional problems for arbitrary potentials.

# Statement of need

`FEMTISE` allows solving the TISE in one-dimensional and two-dimensional systems. However, it is possible to extend the package to solve systems of arbitrary dimensions (by adapting existing functions in the [GRIDAP.jl](https://github.com/gridap/Gridap.jl) package to create multidimensional cartesian grids). By using FEM to solve TISE, and because this is a variational method, we ensure that the numerically calculated ground state energy is bounded from below by its theoretical value. Additionally, it allows for quickly solving the TISE using an input file that defines arbitrary potentials. Another very valuable feature is that it enables the calculation of how eigenenergies of a Hamiltonian change when varying any potential parameter.

# Capability of the package
This package is under construction and is focus to resolve the TISE over high performance protocol using functions from [GRIDAP.jl](https://github.com/gridap/Gridap.jl) package. Main specific features of the package are: possibility of multi-thread parallelization; solver function for usual potential (e.g.: one and two dimensional quantum harmonic oscillator, one dimensional finite well potential and symmetric finite one dimensional Kronig-Penney potential); compute eigenvalue as a function of arbitrary potential parameter; solver function for two particles problems with different masses; and involve [JLD2](https://github.com/JuliaIO/JLD2.jl) output format which allow us keep going the simulation or binary output which allow us write result data as value arrays to post-processing.

# Weak formulation and eigenvalue problem (Galerkin Method)

Considering the TISE we have $\hat{H}\ket{\psi } =\epsilon \ket{\psi }$ in coordinate representation we could write the equation by their Sturm-Liouville form as

$\underbrace{\hat{H}\ket{\psi } =\epsilon \ket{\psi }}_{\mathrm{TISE}} \Rightarrow \overbrace{-\nabla \cdotp ( p\nabla \psi ) +q\psi =\lambda r\psi }^{\mathrm{Sturm-Liouville}} \text{  if} \ \hat{H} =-\frac{\hbar^2}{2m}\nabla^2 +V( r) \Rightarrow \begin{cases} p( x) =\frac{\hbar ^{2}}{2m}(  >0)\\ q( x) =V( r)\\ r( x) =1(  >0) \end{cases}$

For a variational problem first we approximate the wave function $\psi \in \mathcal{H}$ by a function $u\in U_{\mathrm{trial}}$ (where $\mathcal{H}$ is a Hilbert space of infinite dimension in principle and where $U_{\mathrm{trial}}$ is a approximate space of finite dimension for $\mathcal{H}$) and multiplying the PDE by a test function $v\in V_{\mathrm{space}}\left( C^{1}\right)$ (notice that we don't need that $v$ satisfy some boundary condition) and integrating over $\Omega$ using the Green formula which says

$\int _{\Omega }( \partial _{j} h) gdx=-\int _{\Omega } h( \partial _{j} g) dx+\int _{\partial \Omega =\Gamma } hgn_{j} ds$

where $n_{j} =\vec{n} \cdotp e_{j}$ is the $j$-th of $\vec{n}$ on the canonical base of $\mathbb{R}^{d}$. So

$-\nabla \cdotp ( p\nabla \psi ) +q\psi =\lambda r\psi \mathrm{\ }(\mathrm{in} \ \Omega ) \\ \mathrm{with} \ \begin{cases} \psi \Bigl|_{\partial \Omega =\Gamma } =0 & (\mathrm{Dirichlet} \ \mathrm{BC})\\ {\frac{\partial \psi }{\partial n}}{\Bigl|}{_{\partial \Omega =\Gamma }}{=\nabla \psi \cdotp }{\vec{n}}{=c} & (\mathrm{Neumann} \ \mathrm{BC}) \end{cases} \Rightarrow -\int _{\Omega }[ \nabla \cdotp ( p\nabla u)] vd\Omega +\int _{\Omega } quvd\Omega =\int _{\Omega } \lambda ruvd\Omega$

> Note: In some cases this method is equivalent to Rayleigh-Ritz method, that is, the Galerkin method is more general but both came from functional minimization.

if $\begin{cases} g:=v & \Rightarrow \partial g=\nabla v\\ \partial h:=\nabla \cdotp ( p\nabla u) & \Rightarrow h=( p\nabla u) \end{cases}$ then

$\Rightarrow \int _{\Omega } p( \nabla u\cdotp \nabla v) d\Omega {-}{\int _{\Gamma }}{p}{(}{\nabla u\cdotp }{\vec{n}}{)}{vd\Gamma } +\int _{\Omega } quvd\Omega = \lambda \int _{\Omega } ruvd\Omega$

and the problem to resolve would be like $a( u,v) =\lambda b( u,v)$ where we had defined 

$a( u,v) :=\int _{\Omega }[ p( \nabla u\cdotp \nabla v) +quv] d\Omega {-}{\int _{\Gamma }}{p}{(}{\nabla u\cdotp }{\vec{n}}{)}{vd\Gamma } \\ b( u,v) :=\int _{\Omega } ruvd\Omega$

> Note: The integral over $\Gamma$ boundary will avoid if we are considering Dirichlet's boundary conditions or has a fixed value is we are considering Neumann's boundary conditions.

Now, we resolve the eigenvalue problem using the **ARPACK** package building affine operators (each affine operator has associated one matrix and one vector) like

$\begin{cases} a( u,v) =0\mathrm{\xrightarrow[affine\ operator]{associated}}\{A,\vec{\alpha }\}\\ b( u,v) =0\mathrm{\xrightarrow[affine\ operator]{associated}}\{B,\vec{\beta }\} \end{cases} \Rightarrow A\vec{\phi } =\lambda B\vec{\phi }$

here $u$ and $v$ functions are such that $u\in U_{\mathrm{trial}} ;v\in V_{\mathrm{space}}$ where $U_{\mathrm{trial}}$ and $V_{\mathrm{space}}$ are finite spaces which approximate the infinite Hilbert spaces where live $\psi$ and $\psi^{*}$.

If we consider a base $\mathcal{B} =\{\phi _{j}( x)\}_{j=1}^{N}$, which expand all of the space $U_{\mathrm{trial}}$, we have $u$ and $v$ function which could be approximated (expanded in $U,V$) as

$\mathcal{B}_{U,V} =\{\phi _{j}( x)\}_{j=1}^{N} \Rightarrow f( x) =\sum\limits _{j=1}^{N} c_{j} \phi _{j}( x) \\ \Rightarrow \sum _{i,j=1}^{N} c_{i} c_{j}\left[\int\limits _{x_{i}}^{x_{f}} p( x)\frac{d\phi _{i}( x)}{dx}\frac{d\phi _{j}( x)}{dx} dx +\int\limits _{x_{i}}^{x_{f}} q( x) \phi _{i}( x) \phi _{j}( x) dx\right] = \dots \\ \dots = \lambda \sum _{i,j=1}^{N} c_{i} c_{j}\left[\int\limits _{x_{i}}^{x_{f}} r( x) \phi _{i}( x) \phi _{j}( x) dx\right] \Rightarrow \hat{A}\vec{\Phi } =\lambda (\hat{B}\vec{\Phi })$

and defining following affine matrices

$A_{ij} := \int\nolimits _{x_{i}}^{x_{f}} p( x)\frac{d\phi _{i}( x)}{dx}\frac{d\phi _{j}( x)}{dx} dx+\int\nolimits _{x_{i}}^{x_{f}} q( x) \phi _{i}( x) \phi _{j}( x) dx \\ B_{ij} := \int\nolimits _{x_{i}}^{x_{f}} r( x) \phi _{i}( x) \phi _{j}( x) dx$

and noting that we are solving at the end is a generalized eigenvalue problem with the form $\hat{A}\vec{\Phi } =\lambda (\hat{B}\vec{\Phi })$.

The computing problem implementation is resolved by [ARPACK](https://github.com/JuliaLinearAlgebra/Arpack.jl/tree/master) package (using LAPACK and BLAS libraries) which give us the possibility to compute only a subset of ordered pair (eigenvalue,eigenvector) like following (this method is called **Shift-Invert mode**):

If $(\vec{\Phi } ,\lambda )$ is a eigenpair for $(\hat{A} ,\hat{B})$ matrices and $\sigma \neq \lambda$ then

$\underbrace{(\hat{A} -\sigma \hat{B})^{-1}\hat{B} \cdotp \vec{\Phi } =\vec{\nu } \cdotp \vec{\Phi }}_{\mathrm{Shift-Invert\ mode}} ;\ \nu _{j} =\frac{1}{( \lambda -\sigma )}$

this allow us to transform the original eigenvalue problem to another one with different eigenvalues, where those $\lambda \approx \sigma$ will have maximum $\nu _{j}$ values.

# References

+ Sun, J. and Zhou, A., 2016. Finite element methods for eigenvalue problems. Chapman and Hall/CRC.
+ Lehoucq, R.B., Sorensen, D.C. and Yang, C., 1998. ARPACK users' guide: solution of large-scale eigenvalue problems with implicitly restarted Arnoldi methods. Society for Industrial and Applied Mathematics.
+ Verdugo, F. and Badia, S., 2022. The software design of Gridap: a finite element package based on the Julia JIT compiler. Computer Physics Communications, 276, p.108341.

# Implementation of FEMTISE package

[Here](https://github.com/mendzmartin/FEMTISE.jl/blob/c4c72d603e9e8516f08a37f966d3ee3b91e7f719/src/functions/eigen_problem_definition_function.jl#L67-L79) you can see the specific implementation of solve function to resolve eigen value problems.