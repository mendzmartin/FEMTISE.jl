# Accuracy and best practices

This page explains how to obtain reliable numbers and how to check them. The
figures quoted come from the benchmarks of the accompanying paper (1D harmonic
oscillator, BenDaniel–Duke quantum well, 2D oscillator); the corresponding checks
are part of the test suite.

## Checklist for a converged result

1. **Mesh convergence.** Repeat the calculation with ``h/2``. With exact quadrature
   the eigenvalues decrease and the change should shrink by ``\approx 2^{2\ell}``
   (16 for quadratic elements). The difference between the two runs estimates the
   error of the coarse one.
2. **Domain size.** Enlarge ``\Omega``. States that feel the walls shift upward.
   For bound states the shift decays exponentially with the distance to the wall.
   High-lying states and continuum states are always box states.
3. **Number of eigenpairs and shift.** Check that `sigma` lies below the levels you
   want (or at the centre of your energy window) and that `nev` covers them all.
4. **Normalization and orthogonality.** `eigenstates_normalization` and
   `orthogonality_check` should return 1 and ≈0.

## Element degree and cost

At the same number of unknowns ``N_h`` the cost is nearly independent of ``\ell``,
while the error drops dramatically with ``\ell``. For the ground state of the 1D
oscillator:

| ``\ell`` | ``N_h`` | relative error |
|---|---|---|
| 1 | 2559 | ``3.8\times10^{-6}`` |
| 2 | 319 | ``6.3\times10^{-7}`` |
| 3 | 479 | ``2.5\times10^{-10}`` |

Use ``\ell=2`` (the driver default) or ``\ell=3`` through the low-level interface.
Second-order finite differences behave like ``\ell=1``, but approach the exact
values from below.

## Quadrature degree

The upper-bound property and the optimal rate assume exact integrals. With a
polynomial potential of degree ``s``, the integrals are exact for
``k_q \ge 2\ell + s``.

!!! note "Driver quadrature"
    The input-file driver uses ``k_q = 3`` with quadratic elements. For
    ``V=x^2/2`` this keeps the ``h^4`` rate but roughly **doubles the error**, and
    the eigenvalues are no longer guaranteed to be upper bounds. When certified
    bounds are needed, use the low-level interface with ``k_q = 2\ell+2`` or more.

## Discontinuous potentials and masses: align the mesh

The weak form handles jumps of ``V`` and ``m^*`` exactly **if the discontinuity
lies on a mesh node** (1D) or a cell face (2D). For a GaAs/AlGaAs quantum well
with quadratic elements:

| mesh | ``h`` (bohr) | relative error of ``\varepsilon_0`` |
|---|---|---|
| aligned | 6.25 | ``+9.3\times10^{-7}`` |
| aligned | 1.56 | ``+3.7\times10^{-9}`` (``h^4``) |
| not aligned | 1.56 | ``-2.0\times10^{-3}`` (``h^1``, below the exact value) |

With the uniform meshes of the driver, choose `L` and `Δx` such that every interface
lies at an integer multiple of `Δx` from the edge of the domain.

## Degeneracies on Cartesian meshes

Symmetries of the Hamiltonian that are not symmetries of the mesh are only
recovered as ``h\to 0``. For the 2D isotropic oscillator on a square Cartesian mesh,
the ``(n+1)``-fold degeneracies are split by ``\sim 10^{-5}``–``10^{-4}`` for
``h=0.25``. States related by ``x_1\leftrightarrow x_2`` remain exactly degenerate.

## Performance tips

- The eigensolver dominates the cost. Assembly is 10–20 times cheaper.
- The solve time grows roughly linearly with ``N_h`` in 1D and as
  ``N_h^{1.4}`` in 2D. Typical times on a laptop: ``0.1`` s for ``N_h = 3\times10^4``
  (1D) and ``1.5`` s for ``N_h = 6\times10^4`` (2D), ``n_{\rm ev}=10``.
- The cost grows faster than linearly with `nev`: ``0.6`` s for 10 and ``23`` s for
  160 eigenpairs at ``N_h\approx 2.5\times10^4``. Split large windows into several
  shifts.
- In loops over a parameter, build the mesh, measures and spaces once.
- Start Julia with threads (`julia -t 4`) to accelerate post-processing.
- The first call of each function includes compilation time. Time the second call.

## Known issues

Open issues in the current release (v0.1.1):

| Issue | Workaround |
|---|---|
| Driver quadrature ``k_q=3`` | low-level interface with ``k_q\ge2\ell+2`` |
| Sweeps start at ``\lambda_i+\Delta\lambda`` | start one step earlier |
| Eigenvectors of degenerate levels are not mutually orthogonal (complex ARPACK driver) | orthonormalize within the degenerate subspace (see [Post-processing](postprocessing.md)) |
| Driver called from inside a function: "method too new" error | `Base.invokelatest(run_default_eigen_problem, set_type_potential(path))` (see [Input file reference](@ref)) |

Fixed in v0.1.1 (present in v0.1.0): Rényi entropy, normalization of reduced
densities, reloading of `bin eigen` results, the `solve` name clash with Gridap,
and the boundary tags of `"simple_rectangle_v2"`. See the [Changelog](@ref).
