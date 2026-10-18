# Changelog

All notable changes to `FEMTISE.jl` are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the
project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
While the version is `0.x`, an increase of the *minor* number (e.g. `0.1 → 0.2`)
may contain breaking changes; *patch* releases (e.g. `0.1.0 → 0.1.1`) only contain
backward-compatible bug fixes.

Each entry is classified as **Added**, **Changed**, **Deprecated**, **Removed**,
**Fixed** or **Security**. Changes that modify numerical results are always marked
explicitly.

## [Unreleased]

Nothing yet.

## [0.1.0] - 2026-10-18

First versioned release. It freezes the state of the package used in the
publications listed in the README, so that their results can be reproduced with
`add https://github.com/mendzmartin/FEMTISE.jl#v0.1.0`.

### Added
- Solution of the time-independent Schrödinger equation in Sturm–Liouville form
  with conforming Lagrange finite elements (`Gridap.jl`) and the shift-and-invert
  Arnoldi eigensolver of `ARPACK` (`Arpack.jl`).
- Low-level interface: `make_model`, `make_boundary_conditions`, `measures`,
  `fe_spaces`, `bilineal_forms`, `eigen_problem`, `solve`,
  `eigen_values_and_eigen_vectors`, `space_coord`.
- Anisotropic kinetic term for two particles of different masses in 2D
  (`bilineal_forms(p₁, p₂, q, r, dΩ)`).
- Input-file driver (`run_default_eigen_problem`, `set_type_potential`,
  `input_data`) for 1D and 2D problems with user-defined potentials, and an
  interactive console with built-in potentials (1D/2D harmonic oscillator, finite
  well, symmetric finite Kronig–Penney chain).
- Parameter sweeps (`analysis_param`) that reuse the mesh for every value of a
  potential parameter.
- Output as binary arrays (`write_bin`/`read_bin`) or JLD2 files with the finite
  element objects; reloading with `collect_result_data`.
- Post-processing: `norm_l2`, `orthogonality_check`, `eigenstates_normalization`,
  `density`, `reduced_density`, `time_indep_entropy` (Shannon and Rényi),
  `reduced_time_indep_entropy`, `time_indep_diff_mutual_information`,
  `state_population`, `coord_first_moment`, `coord_second_moment`,
  `coord_standar_deviation`.
- Documentation reorganized into *Getting started*, *Manual* (theory,
  architecture, input-file reference, low-level API, outputs, post-processing,
  accuracy and best practices) and an API reference grouped by module;
  versioned documentation (`dev`, `stable` and one page per release).
- This changelog.

### Changed
- Test suite reorganized into validation tests with analytic references: upper
  bounds and ``h^{2\ell}`` convergence for ``\ell = 1, 2, 3``, position-dependent
  mass (BenDaniel–Duke), two particles of different masses, post-processing
  observables, binary I/O, and the input-file driver end to end (1D, 2D, sweeps,
  JLD2 and binary outputs). Test files no longer activate environments themselves.
- The orthogonality check of the 2D oscillator example skips pairs of degenerate
  eigenvectors (they are not guaranteed to be orthogonal, see Known issues); the
  previous check could fail randomly.
- CI: updated GitHub Actions (`checkout@v4`, `setup-julia@v2`, `cache@v2`) and
  coverage upload (`julia-processcoverage` + `codecov-action@v4`).
- Removed the stale `test/Manifest.toml` (Julia 1.9).

### Known issues
These issues are present in 0.1.0 and are scheduled for the next releases. They
are documented in the *Accuracy and best practices* page of the manual.
- `time_indep_entropy(ψ, U, dΩ, α)` (Rényi entropy) evaluates
  `exp∘(α*(log2∘ρ))`, i.e. ``ρ^{α/\ln 2}`` instead of ``ρ^α``. The Shannon entropy
  is not affected. *The fix will change numerical results.*
- `make_boundary_conditions("simple_rectangle_v2", ...)` uses `==` instead of `=`
  when setting the boundary tags.
- `FEMTISE.solve` and `Gridap.solve` are both exported; code that loads both
  packages must call `FEMTISE.solve` explicitly.
- The input-file driver integrates with a quadrature of degree 3 for quadratic
  elements, which is not exact for polynomial potentials. Convergence remains
  ``O(h^4)`` but the error is about twice that of exact quadrature, and the
  eigenvalues are no longer guaranteed to be upper bounds.
- In parameter sweeps, the first computed value is ``λ_i + Δλ`` rather than
  ``λ_i``.
- `collect_result_data` cannot reload non-sweep results written in the
  `bin eigen` format (`DefaultBinEigenProblem*` declare the eigenfunctions as
  `Vector{CellField}`); read the files with `read_bin`.
- `reduced_density` (through `reduced_integration`) normalizes with the span of
  the sampling points, ``L-\Delta x``, instead of the domain length ``L``: the
  reduced densities are too small by the factor ``1-\Delta x/L``. Their shape is
  correct.
- Eigenvectors belonging to a degenerate eigenvalue are linearly independent but
  not necessarily orthogonal (complex, non-Hermitian ARPACK driver).
- When the input-file driver is called from inside a function, the potential
  defined by the potential file is "too new" for that function; use
  `Base.invokelatest(run_default_eigen_problem, set_type_potential(path))`.

[Unreleased]: https://github.com/mendzmartin/FEMTISE.jl/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/mendzmartin/FEMTISE.jl/releases/tag/v0.1.0
