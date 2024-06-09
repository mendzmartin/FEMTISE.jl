---
title: 'FEMTISE: Julia package to resolve Time-Independent Schrödinger equation by (F)init (E)lement (M)ethod. This is an implementation of Gridap package for unidimensional and bidimensional grids.'
tags:
  - Julia
  - Physics
  - Time-Independent Schrödinger equation
  - Quantum mechanics
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
 - name: Insituto de Física Enrique Gaviola (IFEG-CONICET).
   index: 2
date: 18 March 2024
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
aas-journal: Physics Journal <- The name of the AAS journal.
---

# Summary

The forces on stars, galaxies, and dark matter under external gravitational
fields lead to the dynamical evolution of structures in the universe. The orbits
of these bodies are therefore key to understanding the formation, history, and
future state of galaxies. The field of "galactic dynamics," which aims to model
the gravitating components of galaxies to study their structure and evolution,
is now well-established, commonly taught, and frequently used in astronomy.
Aside from toy problems and demonstrations, the majority of problems require
efficient numerical tools, many of which require the same base code (e.g., for
performing numerical orbit integration).

# Statement of need

`Gala` is an Astropy-affiliated Python package for galactic dynamics. Python
enables wrapping low-level languages (e.g., C) for speed without losing
flexibility or ease-of-use in the user-interface. The API for `Gala` was
designed to provide a class-based and user-friendly interface to fast (C or
Cython-optimized) implementations of common operations such as gravitational
potential and force evaluation, orbit integration, dynamical transformations,
and chaos indicators for nonlinear dynamics. `Gala` also relies heavily on and
interfaces well with the implementations of physical units and astronomical
coordinate systems in the `Astropy` package [@astropy] (`astropy.units` and
`astropy.coordinates`).

`Gala` was designed to be used by both astronomical researchers and by
students in courses on gravitational dynamics or astronomy. It has already been
used in a number of scientific publications [@Pearson:2017] and has also been
used in graduate courses on Galactic dynamics to, e.g., provide interactive
visualizations of textbook material [@Binney:2008]. The combination of speed,
design, and support for Astropy functionality in `Gala` will enable exciting
scientific explorations of forthcoming data releases from the *Gaia* mission
[@gaia] by students and experts alike.

# Mathematics: Formulación débil y problema de autovalores (método de Galerkin)

Partiendo de la ecuación de Schrödinger independiente del tiempo tendremos $\hat{H}\ket{\psi } =\epsilon \ket{\psi }$ en representación coordenada podremos escribir esta ecuación en su formulación de Sturm-Liouville como sigue,

$$
\underbrace{\hat{H}\ket{\psi } =\epsilon \ket{\psi }}_{\mathrm{TISE}} \Rightarrow \overbrace{-\nabla \cdotp ( p\nabla \psi ) +q\psi =\lambda r\psi ;}^{\mathrm{Sturm-Liouville}} \ \mathrm{if} \ \hat{H} =-\frac{i\hbar }{2m}\vec{\nabla } +V( r) \Rightarrow \begin{cases}
p( x) =\frac{\hbar ^{2}}{2m_{e}}(  >0)\\
q( x) =V( r)\\
r( x) =1(  >0)
\end{cases}
$$

Para una formulación variacional primero aproximamos la función de onda $\psi \in \mathcal{H}$ por una función $u\in U_{\mathrm{trial}}$ (donde $\mathcal{H}$ es un espacio de Hilbert, en principio, de dimensión infinita y donde $U_{\mathrm{trial}}$ es un espacio que aproxima a $\mathcal{H}$ de dimensión finita) y multiplicamos la PDE por una función de prueba $v\in V_{\mathrm{space}}\left( C^{1}\right)$ (notemos que no se requiere que $v$ satisfaga alguna condición de contorno), e integramos sobre $\Omega$ usando la formula de Green que nos dice

$$
\overbrace{\int _{\Omega }( \partial _{j} h) gdx=-\int _{\Omega } h( \partial _{j} g) dx+\int _{\partial \Omega =\Gamma } hgn_{j} ds}^{\mathrm{Green\ theorem}}
$$

donde $n_{j} =\vec{n} \cdotp e_{j}$ es la coordenada $j$-ésima de $\vec{n}$ en la base canónica de $\mathbb{R}^{d}$. Entonces,


# Citations

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

# References