# **FEMTISE.jl**

<picture>
<img alt="FEMTISE logo" src="/images/logo_FEMTISE.svg" width="200" height="200" align="right">
</picture>

| **GitHub Actions - Workflows** |
| ------------ |
| [![Documentation](https://github.com/mendzmartin/FEMTISE.jl/actions/workflows/documentation.yml/badge.svg)](https://github.com/mendzmartin/FEMTISE.jl/actions/workflows/documentation.yml) |
| [![Build Status](https://github.com/mendzmartin/FEMTISE.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mendzmartin/FEMTISE.jl/actions/workflows/CI.yml?query=branch%3Amain) |
| [![CompatHelper](https://github.com/mendzmartin/FEMTISE.jl/actions/workflows/CompatHelper.yml/badge.svg)](https://github.com/mendzmartin/FEMTISE.jl/actions/workflows/CompatHelper.yml) |
| [![TagBot](https://github.com/mendzmartin/FEMTISE.jl/actions/workflows/TagBot.yml/badge.svg)](https://github.com/mendzmartin/FEMTISE.jl/actions/workflows/TagBot.yml) |

<!-- [![Codecov](https://app.codecov.io/gh/mendzmartin/FEMTISE.jl/branch/main/graph/badge.svg)](https://app.codecov.io/gh/mendzmartin/FEMTISE.jl) -->

## **WARNING!**
***->> The package is under construction <<-***

## Table of Contents

- [**FEMTISE.jl**](#femtisejl)
  - [**WARNING!**](#warning)
  - [Table of Contents](#table-of-contents)
  - [**Usage**](#usage)
  - [**Installation**](#installation)
  - [**Documentation**](#documentation)
  - [**Tutorials**](#tutorials)
  - [**To cite this work**](#to-cite-this-work)
    - [Publications](#publications)
  - [**Contact**](#contact)

## **Usage**
Julia package repository to resolve Time-Independent Schrödinger equation (TISE) by Finite Element Method (FEM). This is an implementation of [GRIDAP](https://github.com/gridap/Gridap.jl) package for unidimensional and bidimensional grids.

Project carried out in my PhD studies of Physics at:
* [Faculty of Mathematics, Astronomy, Physics and Computation (FAMAF)](https://www.famaf.unc.edu.ar/)
* The Enrique Gaviola Institute of Physics (IFEG)

## **Installation**
- [Clone Package](https://mendzmartin.github.io/FEMTISE.jl/dev/guide_information/#**1.-Clone-FEMTISE-package**)
- [Build Julia code to use FEMTISE package](https://mendzmartin.github.io/FEMTISE.jl/dev/guide_information/#**2.-Build-Julia-code-to-use-FEMTISE-package**)
- [Run from Julia REPL](https://mendzmartin.github.io/FEMTISE.jl/dev/guide_information/#Run-from-Julia-REPL)

## **Documentation**
- Check [here](https://mendzmartin.github.io/FEMTISE.jl/) the GitHub Page to see documentation.
- Check [here](https://deepwiki.com/mendzmartin/FEMTISE.jl) the DeepWiki Page to learn about FEMTISE software.

## **Tutorials**
See some use examples in [test folder](https://github.com/mendzmartin/FEMTISE.jl/tree/main/test) where you can find some tests of specific `FEMTISE` package functions:
- Finite Well Potential 1D
- Symmetric Finite Kronig-Penney Potential 1D
- Isotropic Quantum Harmonic Oscillator 2D
- Default Eigen Value Problem

Also, you can find tutorials in [FEMTISE_TUTORIAL](https://github.com/mendzmartin/FEMTISE_TUTORIAL) repository of specific problems like:
- Tutorial To Simulate Symmetric Finite One Dimensional Kronig-Penney Potential
- Tutorial To Simulate Isotropic One Dimensional Harmonic Oscillator
- Tutorial To Simulate Isotropic Two Dimensional Harmonic Oscillator
- Tutorial To Simulate Coulomb 2D Interaction potential (Helium atom model)

or check FEMTISE Package [GitHub Page](https://mendzmartin.github.io/FEMTISE.jl/) to see documentation.

## **To cite this work**
If you use this package in your research, please cite it using the following BibTeX entry:
```bibtex
  @misc{Mendez2024FEMTISE,
    title = {mendzmartin/{FEMTISE}.jl},
    copyright = {MIT},
    url = {https://github.com/mendzmartin/FEMTISE.jl},
    bstract = {Variational approximation by Gridap package to resolve T.I.S.E. in Julia},
    urldate = {2024-07-19},
    author = {Mendez, Martin},
    year = {2024},
  }
```
### Publications
This package is used in the following publications:
+ [Mendez, M., & Pont, F. M. (2025). Dynamics of correlations and entanglement generation in electron-molecule inelastic scattering. Physical Review A, 111(5), 052411.](https://doi.org/10.1103/PhysRevA.111.052411)


## **Contact**
Please, contact the project administrator [Mendez Martin](mailto:martinmendez@unc.edu.ar) for any improve suggestion or questions about package use.