---
title: 'FEMTISE.jl: A High Performance Julia Package for Time-Independent Schrödinger Equation by Finite Element Method'
tags:
  - Julia
  - Quantum Physics
  - Time-Independent Schrödinger equation
  - Finit Element Method
authors:
  - name: Méndez, Martín
    orcid: 0000-0000-0000-0000
    equal-contrib: true
    affiliation: "1, 2"
  - name: Pont, Federico M.
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: "1, 2"
affiliations:
 - name: [Faculty of Mathematics, Astronomy, Physics and Computation (FaMAF)](https://www.famaf.unc.edu.ar/)
   index: 1
 - name: Insituto de Física Enrique Gaviola (IFEG-FAMAF-CONICET).
   index: 2
date: 18 March 2024
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
aas-journal: Physics Journal <- The name of the AAS journal.
---

# `FEMTISE.jl`: A High Performance Julia Package for Time-Independent Schrödinger Equation by Finite Element Method

# Abstract

`FEMTISE.jl` is a package to resolve Time-Independent Schrödinger Equation (TISE) by Finit Element Method (FEM). This is an implementation and extension over \verb|GRIDAP.jl| package using high performance protocol. The package is under constructions and currently can solve one and two dimensional problems for arbitrary potentials.

# Statement of need

`FEMTISE` permite resolver la TISE en sistemas unidimensionales y bidimensionales, sin embargo, es posible extender el paquete para resolver sistema de dimensión arbitraria (adaptando funciones existentes en el paquete \verb|Gridap| para crear grillas cartesianas multidimensionales). Al utilizar el FEM para resolver la TISE, y al ser este un método variacional, nos aseguramos que la energía del estado fundamental calculada numéricamente esté acotada inferiormente por su valor teórico. Además, permite resolver rápidamente la TISE a partir de un input file, definiendo potenciales arbitrarios. Otra utilidad muy valiosa es que permite calcular cómo se modifican las autoenergías de un hamiltoniano al variar cualquier parámetro del potencial.

# Capability of the package
This package is under construction and is focus to resolve the TISE over high performance protocol using functions from **GRIDAP** package. Main specific features of the package are: possibility of multi-thread parallelization; solver function for usual potential (e.g.: one and two dimensional quantum harmonic oscillator, one dimensional finite well potential and symmetric finite one dimensional Kronig-Penney potential); compute eigenvalue as a function of arbitrary potential parameter; solver function for two particles problems with different masses; and involve JLD2 output format which allow us keep going the simulation or binary output which allow us write result data as value arrays to post-processing.

# Weak formulation and eigenvalue problem (Galerkin Method)

Considering the TISE we have $\hat{H}\ket{\psi } =\epsilon \ket{\psi }$ in coordinate representation we could write the equation by their Sturm-Liouville form as

$$
    \underbrace{\hat{H}\ket{\psi } =\epsilon \ket{\psi }}_{\mathrm{TISE}} \Rightarrow \overbrace{-\nabla \cdotp ( p\nabla \psi ) +q\psi =\lambda r\psi }^{\mathrm{Sturm-Liouville}}
    \text{  if} \ \hat{H} =-\frac{i\hbar }{2m}\vec{\nabla } +V( r) \Rightarrow
        \begin{cases}
            p( x) =\frac{\hbar ^{2}}{2m_{e}}(  >0)\\
            q( x) =V( r)\\
            r( x) =1(  >0)
        \end{cases}
$$

For a variational problem first we approximate the wave function $\psi \in \mathcal{H}$ by a function $u\in U_{\mathrm{trial}}$ (where $\mathcal{H}$ is a Hilbert space of infinite dimension in principle and where $U_{\mathrm{trial}}$ is a approximate space of finite dimension for $\mathcal{H}$) and multiplying the PDE by a test function $v\in V_{\mathrm{space}}\left( C^{1}\right)$ (notice that we don't need that $v$ satisfy some boundary condition) and integrating over $\Omega $ using the Green formula which says

$$
    \int _{\Omega }( \partial _{j} h) gdx=-\int _{\Omega } h( \partial _{j} g) dx+\int _{\partial \Omega =\Gamma } hgn_{j} ds
$$

where $n_{j} =\vec{n} \cdotp e_{j}$ is the $j$-th of $\vec{n}$ on the canonical base of $\mathbb{R}^{d}$. So

$$
    -\nabla \cdotp ( p\nabla \psi ) +q\psi =\lambda r\psi \mathrm{\ }(\mathrm{in} \ \Omega ) \\ 
    \mathrm{with} \ \begin{cases}
    \psi \Bigl|_{\partial \Omega =\Gamma } =0 & (\mathrm{Dirichlet} \ \mathrm{BC})\\
    {\frac{\partial \psi }{\partial n}}{\Bigl|}{_{\partial \Omega =\Gamma }}{=\nabla \psi \cdotp }{\vec{n}}{=c} & (\mathrm{Neumann} \ \mathrm{BC})
    \end{cases}
    \Rightarrow -\int _{\Omega }[ \nabla \cdotp ( p\nabla u)] vd\Omega +\int _{\Omega } quvd\Omega =\int _{\Omega } \lambda ruvd\Omega
$$

${}^{*}$ En algunos casos es equivalente al método de Rayleigh-Ritz. El método de Galerkin es más general y ambos parten por minimizar un funcional.

if $\begin{cases} g:=v & \Rightarrow \partial g=\nabla v\\ \partial h:=\nabla \cdotp ( p\nabla u) & \Rightarrow h=( p\nabla u) \end{cases}$ then

$$
    \Rightarrow \int _{\Omega } p( \nabla u\cdotp \nabla v) d\Omega {-}{\int _{\Gamma }}{p}{(}{\nabla u\cdotp }{\vec{n}}{)}{vd\Gamma } +\int _{\Omega } quvd\Omega = \lambda \int _{\Omega } ruvd\Omega
$$

and the problem to resolve would be like $a( u,v) =\lambda b( u,v)$ where we had defined 

$$
    a( u,v) :=\int _{\Omega }[ p( \nabla u\cdotp \nabla v) +quv] d\Omega {-}{\int _{\Gamma }}{p}{(}{\nabla u\cdotp }{\vec{n}}{)}{vd\Gamma } \\
    b( u,v) :=\int _{\Omega } ruvd\Omega 
$$

${}^{*}$ Notice that integral over $\Gamma$ boundary will avoid if we are considering Dirichlet's boundary conditions or has a fixed value is we are considering Neumann's boundary conditions.

Now, we resolve the eigenvalue problem using the **ARPACK** package building affine operators (each affine operator has associated one matrix and one vector) like

$$
    \begin{cases}
        a( u,v) =0\mathrm{\xrightarrow[affine\ operator]{associated}}\{A,\vec{\alpha }\}\\
        b( u,v) =0\mathrm{\xrightarrow[affine\ operator]{associated}}\{B,\vec{\beta }\}
    \end{cases} \Rightarrow A\vec{\phi } =\lambda B\vec{\phi }
$$

here $u$ and $v$ functions are such that $u\in U_{\mathrm{trial}} ;v\in V_{\mathrm{space}}$ where $U_{\mathrm{trial}}$ and $V_{\mathrm{space}}$ are finite spaces which approximate the infinite Hilbert spaces where live $\psi$ and $\psi^{*}$.

If we consider a base $\mathcal{B} =\{\phi _{j}( x)\}_{j=1}^{N}$, which expand all of the space $U_{\mathrm{trial}}$, we have $u$ and $v$ function which could be approximated (expanded in $U,V$) as

$$
    \mathcal{B}_{U,V} =\{\phi _{j}( x)\}_{j=1}^{N} \Rightarrow f( x) =\sum\limits _{j=1}^{N} c_{j} \phi _{j}( x) \\
    \Rightarrow \sum _{i,j=1}^{N} c_{i} c_{j}\left[\int\limits _{x_{i}}^{x_{f}} p( x)\frac{d\phi _{i}( x)}{dx}\frac{d\phi _{j}( x)}{dx} dx +\int\limits _{x_{i}}^{x_{f}} q( x) \phi _{i}( x) \phi _{j}( x) dx\right] =
    \lambda \sum _{i,j=1}^{N} c_{i} c_{j}\left[\int\limits _{x_{i}}^{x_{f}} r( x) \phi _{i}( x) \phi _{j}( x) dx\right] \Rightarrow \hat{A}\vec{\Phi } =\lambda (\hat{B}\vec{\Phi })
$$

and defining following affine matrices

$$
    A_{ij} := \int\nolimits _{x_{i}}^{x_{f}} p( x)\frac{d\phi _{i}( x)}{dx}\frac{d\phi _{j}( x)}{dx} dx+\int\nolimits _{x_{i}}^{x_{f}} q( x) \phi _{i}( x) \phi _{j}( x) dx \\
    B_{ij} := \int\nolimits _{x_{i}}^{x_{f}} r( x) \phi _{i}( x) \phi _{j}( x) dx
$$

and noting that we are solving at the end is a generalized eigenvalue problem with the form $\hat{A}\vec{\Phi } =\lambda (\hat{B}\vec{\Phi })$.

The computing problem implementation is resolved by ARPACK package (using LAPACK and BLAS libraries) which give us the possibility to compute only a subset of ordered pair (eigenvalue,eigenvector) like following (this method is called **Shift-Invert mode**):

If $(\vec{\Phi } ,\lambda )$ is a eigenpair for $(\hat{A} ,\hat{B})$ matrices and $\sigma \neq \lambda$ then

$$
    \underbrace{(\hat{A} -\sigma \hat{B})^{-1}\hat{B} \cdotp \vec{\Phi } =\vec{\nu } \cdotp \vec{\Phi }}_{\mathrm{Shift-Invert\ mode}} ;\ \nu _{j} =\frac{1}{( \lambda -\sigma )}
$$

this allow us to transform the original eigenvalue problem to another one with different eigenvalues, where those $\lambda \approx \sigma$ will have maximum $\nu _{j}$ values.

<!-- # Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References -->