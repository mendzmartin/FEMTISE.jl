# Autocorrelation function and spectra

The eigenpairs computed by FEMTISE give the exact time evolution of any wave packet
in the same Hamiltonian. This page describes how to obtain from them the
**autocorrelation function** and its **Fourier transform**, the spectrum, which is
the quantity measured in many experiments (absorption, photoelectron and
Franck–Condon spectra). A complete worked example is the tutorial
[Autocorrelation function and spectrum](../examples/autocorrelation_spectrum.md).

## Definitions

In atomic units (``\hbar = 1``) a wave packet ``\psi_0`` evolves as
``\psi(t) = e^{-i\hat H t}\psi_0``. Expanding it in the eigenstates
``\hat H\phi_n = \varepsilon_n\phi_n``,

```math
\psi_0 = \sum_n c_n\,\phi_n,\qquad c_n = \langle\phi_n|\psi_0\rangle ,
```

the time evolution is ``\psi(t) = \sum_n c_n e^{-i\varepsilon_n t}\phi_n``, and

```math
C(t) = \langle\psi_0|\psi(t)\rangle = \sum_n |c_n|^2\, e^{-i\varepsilon_n t},
\qquad
\sigma(\omega) = \frac{1}{2\pi}\int_{-\infty}^{\infty} C(t)\,e^{i\omega t}\,dt
= \sum_n |c_n|^2\,\delta(\omega-\varepsilon_n).
```

- ``|C(t)|^2`` is the **survival probability**: the probability of finding the
  system in its initial state at time ``t``.
- ``\sigma(\omega)`` is the **spectrum** of the wave packet: a line at each
  eigenvalue, weighted by the overlap ``|c_n|^2`` (Franck–Condon factor).
- The two contain the same information. Recurrences of ``|C(t)|`` at times
  ``T = 2\pi/\Delta\varepsilon`` correspond to regularly spaced lines with spacing
  ``\Delta\varepsilon``.

## Finite times and line shapes

A measurement or a simulation lasts a finite time. Multiplying ``C(t)`` by a window
``w(t)`` before the transform turns every delta function into a line of finite
width. Since ``C(-t) = C(t)^*``,
``\sigma(\omega) = \frac{1}{\pi}\mathrm{Re}\int_0^\infty C(t)\,w(t)\,e^{i\omega t}dt``, and

| Window ``w(t)`` | Line shape ``L(\omega-\varepsilon_n)`` | Parameter |
|---|---|---|
| ``e^{-\gamma t}`` | Lorentzian ``\frac{1}{\pi}\frac{\gamma}{(\omega-\varepsilon_n)^2+\gamma^2}`` | half width at half maximum ``\gamma`` (a lifetime ``1/2\gamma``) |
| ``e^{-\gamma^2t^2/2}`` | Gaussian ``\frac{1}{\sqrt{2\pi}\gamma}e^{-(\omega-\varepsilon_n)^2/2\gamma^2}`` | standard deviation ``\gamma`` (resolution) |

Every line has unit area, so ``\int\sigma\,d\omega = \sum_n|c_n|^2``.

## Functions

| Function | Computes |
|---|---|
| [`expansion_coefficients`](@ref)`(ψ₀, states, dΩ)` | ``c_n = \langle\phi_n|\psi_0\rangle`` |
| [`autocorrelation`](@ref)`(c, energies, t)` | ``C(t)`` at one or many times |
| [`autocorrelation_spectrum`](@ref)`(c, energies, ω; width, lineshape)` | broadened ``\sigma(\omega)`` directly from the eigenpairs |
| [`fourier_transform`](@ref)`(t, C, ω; window, width)` | ``\sigma(\omega)`` from a sampled ``C(t)`` |

`autocorrelation_spectrum` and `fourier_transform` use the same conventions. With
matching window and line shape they agree (to the truncation error
``e^{-\gamma T}`` of the exponential window), which the test suite checks.
`fourier_transform` is also useful for autocorrelations computed by other means,
e.g. a direct time propagation or experimental data.

```julia
using FEMTISE, Gridap
model = CartesianDiscreteModel((-10.0, 10.0), (200,))
sol = solve_tise(model, x -> 0.5 * x[1]^2; nev = 30)

ψ₀(x) = exp(-0.5 * (x[1] - 2.0)^2) / π^0.25          # displaced ground state
c = expansion_coefficients(ψ₀, sol.states, sol.dΩ)
sum(abs2, c)                                         # ≈ 1: 30 states are enough

t = 0.0:0.05:60.0
C = autocorrelation(c, sol.energies, t)
ω = 0.0:0.01:6.0
σ = autocorrelation_spectrum(c, sol.energies, ω; width = 0.1)               # Lorentzian
σ_from_C = fourier_transform(collect(t), C, ω; window = :exponential, width = 0.1)
```

## Practical checks

- **Completeness.** Compute enough eigenstates that ``\sum_n|c_n|^2`` is close to the
  norm of ``\psi_0`` (usually 1). Missing weight means missing lines.
- **Degenerate levels.** The eigenvectors of a degenerate level are not necessarily
  orthogonal (see [Post-processing](postprocessing.md)). Orthonormalize them before
  computing ``c_n``. Their total weight ``\sum|c_n|^2`` within the level does not
  depend on the choice of basis.
- **Time step and range** (for `fourier_transform`). Sampling with step ``\Delta t``
  resolves frequencies up to ``\pi/\Delta t``. With a window of width ``\gamma``,
  propagate to ``T \gtrsim 10/\gamma`` so that the window has decayed.
- **Accuracy at long times.** Errors ``\delta\varepsilon_n`` of the eigenvalues give
  phase errors ``\delta\varepsilon_n t`` in ``C(t)``. For long times, use more accurate
  eigenvalues (finer mesh, higher element degree).
