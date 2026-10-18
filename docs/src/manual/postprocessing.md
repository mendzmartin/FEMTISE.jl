# Post-processing

The eigenfunctions returned by the solver are Gridap finite element functions
(`CellField`). Every observable that can be written as an integral over the mesh is
computed with the same quadrature used to build the Hamiltonian. The functions
below act on a vector of eigenfunctions `ϕ` and return one value per state.
`U` is the trial space, `Ω` the triangulation and `dΩ` the measure (see
[Low-level interface](@ref)).

## Checks

| Function | Computes |
|---|---|
| [`norm_l2`](@ref)`(ψ, dΩ)` | ``\big(\int|\psi|^2\,d\Omega\big)^{1/2}`` |
| [`eigenstates_normalization`](@ref)`(ϕ, dΩ)` | the norm of every state |
| [`orthogonality_check`](@ref)`(ϕ, dΩ)` | ``|\langle\phi_j|\phi_i\rangle|`` for all pairs ``j<i`` |

The eigenvectors returned by ARPACK are normalized, and eigenvectors of different
eigenvalues are orthogonal. Within a **degenerate** subspace (e.g. the ``(n+1)``-fold
levels of the 2D isotropic oscillator) the complex ARPACK driver returns linearly
independent vectors that are *not necessarily orthogonal*. If an orthonormal basis of
the subspace is needed, orthonormalize it, e.g. with Gram–Schmidt using
`sum(∫(ϕ[j]' * ϕ[i]) * dΩ)` as inner product.

## Densities

| Function | Computes |
|---|---|
| [`density`](@ref)`(ϕ)` | ``\rho_n(\mathbf x) = |\phi_n(\mathbf x)|^2`` as finite element functions |
| [`reduced_density`](@ref)`(ϕ, r, model)` | 2D: ``\rho^{(1)}(x)=\int\rho\,dy`` and ``\rho^{(2)}(y)=\int\rho\,dx`` sampled at the points `r = (x_vector, y_vector)` |

!!! warning "Known issue in v0.1.0: normalization of reduced densities"
    `reduced_density` returns ``(1-\Delta x/L)\,\rho^{(1)}`` instead of
    ``\rho^{(1)}`` (``\Delta x`` is the spacing of the sampling points and ``L`` the
    length of the domain). The shape is correct: divide by the integral of the
    returned values to normalize.

## Moments of the coordinates

| Function | Computes |
|---|---|
| [`coord_first_moment`](@ref)`(ϕ, U, Ω, dΩ, k)` | ``\langle x_k\rangle`` |
| [`coord_second_moment`](@ref)`(ϕ, U, Ω, dΩ, k)` | ``\langle x_k^2\rangle`` |
| [`coord_standar_deviation`](@ref)`(ϕ, U, Ω, dΩ, k)` | ``\sigma_{x_k}=\sqrt{\langle x_k^2\rangle-\langle x_k\rangle^2}`` |

## Information-theoretic quantities

| Function | Computes |
|---|---|
| [`time_indep_entropy`](@ref)`(ϕ, U, dΩ)` | Shannon entropy ``S=-\int\rho\log_2\rho\,d\Omega`` (bits) |
| [`time_indep_entropy`](@ref)`(ϕ, U, dΩ, α)` | Rényi entropy ``S_\alpha=\frac{1}{1-\alpha}\log_2\int\rho^\alpha d\Omega`` — see the known issue below |
| [`reduced_time_indep_entropy`](@ref)`(r, ρ)` | for a density ``\rho(x_i, y_j)`` sampled on a grid: ``-\int\rho(x_i,y)\log_2\rho(x_i,y)\,dy`` for every ``x_i``, and the analogous quantity for every ``y_j`` (trapezoidal rule) |
| [`time_indep_diff_mutual_information`](@ref)`(Sx, Sy, Sxy)` | ``I = S_x + S_y - S_{xy}`` |

To obtain the mutual information between the two coordinates of a 2D state,
compute the entropies of the reduced densities ``\rho^{(1)}``, ``\rho^{(2)}`` (one
number each) and the total entropy ``S``, then
``I = S^{(1)} + S^{(2)} - S``. For two particles, ``I`` measures the correlation
between them; it vanishes for product states.

!!! warning "Known issue in v0.1.0: Rényi entropy"
    `time_indep_entropy(ϕ, U, dΩ, α)` evaluates ``\rho^{\alpha/\ln 2}`` instead of
    ``\rho^\alpha``. Until it is fixed, compute the Rényi entropy directly:
    `log2(sum(∫((z -> z^α) ∘ ρ) * dΩ)) / (1 - α)` with `ρ = density(ϕ)[n]`. The Shannon entropy
    is correct.

## Writing your own observables

Any integral can be written with Gridap syntax. For example the expectation value of
the potential and the kinetic energy of state `n`:

```julia
ψ = ϕ[n]
V_mean = real(sum(∫(ψ' * q * ψ) * dΩ))
T_mean = real(sum(∫(p * (∇(ψ)' ⋅ ∇(ψ))) * dΩ))
```

and the overlap of two states computed on the same mesh,
`sum(∫(ψ1' * ψ2) * dΩ)`.
