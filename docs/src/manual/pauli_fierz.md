# Three degrees of freedom and cavity QED (Pauli–Fierz)

[`solve_tise`](@ref) works on one-, two- and three-dimensional meshes, and each mesh
coordinate may belong to a different particle with its own mass. This makes it
possible to solve, without any approximation beyond the discretization, problems with
three one-dimensional degrees of freedom, e.g. a **nucleus, an electron and a cavity
photon mode**. A worked example is the tutorial
[Nucleus, electron and photon](../examples/pauli_fierz_cavity.md).

## Masses

The kinetic energy is ``-\sum_k \frac{1}{2M_k}\partial_k^2``, with masses in units of
the electron mass:

```julia
masses = (1836.15, 1.0, 1.0)     # proton (R), electron (z), photon mode (q)
sol = solve_tise(model, V; masses = masses, nev = 10)
```

In the input-file driver, the same is set with the optional line
`masses = 1836.15 1.0 1.0` (see [Input file reference](@ref)).

For a diatomic molecule, the nuclear coordinate ``R`` uses the reduced mass
``\mu = M_AM_B/(M_A+M_B)``.

## The Pauli–Fierz Hamiltonian

In the length gauge and the dipole approximation, matter with dipole moment
``d(\mathbf r)`` coupled to one cavity mode of frequency ``\omega`` is described by
(atomic units)

```math
\hat H = \hat H_{\rm matter}(\mathbf r)
- \frac{1}{2}\frac{\partial^2}{\partial q^2}
+ \frac{1}{2}\omega^2 q^2
+ g\sqrt{2\omega}\; q\, d(\mathbf r)
+ \frac{g^2}{\omega}\, d(\mathbf r)^2 .
```

- ``q`` is the displacement coordinate of the photon mode, with unit mass. Its
  harmonic levels ``\omega(n+\tfrac12)`` are the photon-number states.
- ``g`` is the light–matter coupling strength (it grows as the cavity volume
  decreases).
- The last term, the **dipole self-energy**, keeps the Hamiltonian bounded from below.
  The three cavity terms combine into
  ``\frac12\omega^2\big(q + g\sqrt{2\omega}\,d/\omega^2\big)^2``.

[`pauli_fierz_potential`](@ref) builds the potential of this Hamiltonian from the
matter potential and the dipole function. The **last mesh coordinate is ``q``**:

```julia
# nucleus R and electron z, charge Z = 1: dipole d = R - z
V_matter(r) = my_potential(r[1], r[2])          # r = (R, z)
dipole(r) = r[1] - r[2]
V = pauli_fierz_potential(V_matter, dipole; frequency = 0.1, coupling = 0.01)

model = CartesianDiscreteModel((0.5, 8.0, -10.0, 10.0, -20.0, 20.0), (24, 32, 32))
sol = solve_tise(model, V; masses = (918.0, 1.0, 1.0), nev = 8, sigma = -1.0)
```

The eigenstates are **polaritons**: mixed light–matter states. Since
``g\sqrt{2\omega}\,q = g(\hat a+\hat a^\dagger)``, when a matter excitation of energy
``\omega_{01}`` is resonant with the cavity (``\omega=\omega_{01}``) the two uncoupled
states "excited matter, no photon" and "ground-state matter, one photon" split into a
lower and an upper polariton separated by the Rabi splitting
``\Omega_R \approx 2g\,|d_{01}|``, with ``d_{01}`` the transition dipole (for weak
coupling).

## Validation

The test suite and the tutorial check the 3D solver against a **quadratic model** with
a harmonic nucleus and electron and a linear dipole. Its potential is a quadratic
form, so the exact levels are ``\sum_k\Omega_k(n_k+\tfrac12)``, with ``\Omega_k^2`` the
eigenvalues of the mass-weighted Hessian. With ``16^3`` quadratic elements, the lowest
levels agree to about ``10^{-3}`` and are upper bounds.

## Cost and practical advice

Three-dimensional meshes grow quickly: a mesh of ``n^3`` cells with quadratic elements
has about ``8n^3`` unknowns, and the sparse factorization used by the eigensolver
costs much more than in 1D or 2D. For reference, ``16^3`` cells (``3\times10^4``
unknowns) take about 16 s on a laptop; memory and time grow quickly with ``n``.
Recommendations:

- Choose the box for each coordinate from the extent of the states of interest; the
  photon coordinate needs about ``\pm\sqrt{2n_{\max}+1}/\sqrt{\omega}``.
- Heavy coordinates (nuclei) have much narrower wave functions than light ones and
  need finer meshes in their direction. Adjust ``n_x, n_y, n_z`` independently.
- Check convergence with two meshes, as described in
  [Accuracy and best practices](accuracy.md).
- For realistic molecules, a common strategy is to treat the nucleus adiabatically
  (2D electron–photon problems along a grid of ``R``, a parameter sweep), and to
  solve the full 3D problem only to check the adiabatic approximation.
