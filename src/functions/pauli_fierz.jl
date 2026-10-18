# =============================================================================
# Pauli–Fierz Hamiltonian for matter coupled to one cavity mode.
#
# In the length gauge and the dipole approximation, a single cavity mode of
# frequency ω (photon displacement coordinate q, mass 1) coupled with strength g to
# the dipole moment d(x) of the matter degrees of freedom x adds to the matter
# Hamiltonian (atomic units, ħ = 1)
#
#   H_cav = −½ ∂²/∂q² + ½ ω² q² + g √(2ω) q d(x) + (g²/ω) d(x)²
#         = −½ ∂²/∂q² + ½ ω² (q + g √(2ω) d(x) / ω²)² ,
#
# the last term of the first line being the dipole self-energy. The matter
# coordinates and q are the coordinates of the finite element mesh: e.g. a nucleus
# (R), an electron (z) and the photon mode (q) form a 3D problem.
# =============================================================================

"""
    pauli_fierz_potential(matter_potential, dipole; frequency, coupling,
                          dipole_self_energy = true)

Potential energy of the Pauli–Fierz Hamiltonian (length gauge, dipole
approximation) for matter coupled to one cavity mode, as a function of the mesh
point ``\\mathbf{x} = (x_1, \\dots, x_{d-1}, q)``, whose **last** coordinate is the
photon displacement ``q``:

```math
V(\\mathbf{x}) = V_{\\rm matter}(\\mathbf{r}) + \\frac{\\omega^2 q^2}{2}
+ g\\sqrt{2\\omega}\\, q\\, d(\\mathbf{r}) + \\frac{g^2}{\\omega}\\, d(\\mathbf{r})^2,
\\qquad \\mathbf{r} = (x_1,\\dots,x_{d-1}).
```

# Arguments
- `matter_potential(r)`, `dipole(r)`: functions of the matter coordinates
  `r = (x₁, …)` (a tuple), in atomic units. For an electron at ``z`` and a nucleus of
  charge ``Z`` at ``R``, the dipole is ``d = ZR - z``.
- `frequency`: cavity frequency ``\\omega``.
- `coupling`: light–matter coupling ``g`` (a.u.).
- `dipole_self_energy`: include the term ``(g^2/\\omega)\\,d^2`` (default). Omitting it
  makes the Hamiltonian unbounded from below for large couplings.

The photon coordinate has unit mass: pass masses `(M₁, …, 1.0)` to
[`solve_tise`](@ref).

# Example: nucleus–electron–photon
```julia
V_matter(r) = my_molecular_potential(r[1], r[2])   # user function, r = (R, z)
d(r) = r[1] - r[2]                              # Z = 1
V = pauli_fierz_potential(V_matter, d; frequency = 0.1, coupling = 0.01)
model = CartesianDiscreteModel((0.5, 6.0, -8.0, 8.0, -20.0, 20.0), (30, 40, 40))
sol = solve_tise(model, V; masses = (918.0, 1.0, 1.0), nev = 6, sigma = -1.0)
```
"""
function pauli_fierz_potential(matter_potential, dipole; frequency::Real, coupling::Real,
                               dipole_self_energy::Bool = true)
    frequency > 0 || throw(ArgumentError("the cavity frequency must be positive"))
    ω, g = Float64(frequency), Float64(coupling)
    bilinear_factor = g * sqrt(2ω)
    self_energy_factor = dipole_self_energy ? g^2 / ω : 0.0
    # Split the mesh point into matter coordinates r and photon coordinate q.
    # Dispatching on the dimension D keeps the split type-stable (this function is
    # evaluated at every quadrature point).
    function potential(x::VectorValue{D}) where {D}
        r = ntuple(k -> x[k], Val(D - 1))
        q = x[D]
        d = dipole(r)
        return matter_potential(r) + 0.5 * ω^2 * q^2 + bilinear_factor * q * d + self_energy_factor * d^2
    end
    return potential
end
