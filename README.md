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

## **Status**
Current release: **v0.1.0** (see [CHANGELOG.md](CHANGELOG.md)). The package is under
active development: while the version is `0.x`, minor releases may change the API.
Pin a release in your projects to keep results reproducible.

## Table of Contents

- [**FEMTISE.jl**](#femtisejl)
  - [**Status**](#status)
  - [**Usage**](#usage)
  - [**Installation**](#installation)
  - [**Documentation**](#documentation)
  - [**Tutorials**](#tutorials)
  - [**To cite this work**](#to-cite-this-work)
    - [Publications](#publications)
  - [**Contact**](#contact)

## **Usage**
Julia package to solve the Time-Independent Schrödinger Equation (TISE) by the Finite
Element Method (FEM), built on [Gridap](https://github.com/gridap/Gridap.jl) and
[ARPACK](https://github.com/JuliaLinearAlgebra/Arpack.jl). It solves 1D and 2D problems
with arbitrary potentials, position-dependent masses and two particles of different
masses, computes eigenvalues as functions of potential parameters, and provides
post-processing tools (densities, moments, entropies, mutual information).

Project carried out in my PhD studies of Physics at:
* [Faculty of Mathematics, Astronomy, Physics and Computation (FAMAF)](https://www.famaf.unc.edu.ar/)
* The Enrique Gaviola Institute of Physics (IFEG)

## **Installation**
Requires Julia ≥ 1.10. In the package manager (`]`) of your project:
```julia
pkg> add https://github.com/mendzmartin/FEMTISE.jl#v0.1.0
```
See the [installation guide](https://mendzmartin.github.io/FEMTISE.jl/stable/getting_started/installation/)
for development installs and running scripts with threads.

## **Documentation**
- [Stable documentation](https://mendzmartin.github.io/FEMTISE.jl/stable/) (latest release) and
  [development documentation](https://mendzmartin.github.io/FEMTISE.jl/dev/) (`main` branch):
  quick start, theory, input-file reference, low-level interface, post-processing,
  accuracy guidelines and API reference.
- [DeepWiki page](https://deepwiki.com/mendzmartin/FEMTISE.jl) about the FEMTISE software.

## **Tutorials**
See the [test folder](https://github.com/mendzmartin/FEMTISE.jl/tree/main/test) for
validation tests against analytic solutions, and the
[FEMTISE_TUTORIAL](https://github.com/mendzmartin/FEMTISE_TUTORIAL) repository for
complete tutorials:
- Symmetric finite one-dimensional Kronig–Penney potential
- Isotropic one-dimensional harmonic oscillator
- Isotropic two-dimensional harmonic oscillator
- Coulomb 2D interaction potential (helium atom model)

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