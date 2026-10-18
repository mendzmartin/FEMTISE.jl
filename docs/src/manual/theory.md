# Theory

This page summarizes the mathematics implemented by `FEMTISE.jl`: what problem is
solved, how it is discretized, and which properties of the result are guaranteed.
A detailed account, with proofs and references, is given in the accompanying paper.

## The problem

The package solves the Sturm–Liouville eigenvalue problem

```math
-\nabla\cdot\big(\mathbf{p}(\mathbf{x})\,\nabla\psi\big) + q(\mathbf{x})\,\psi
  = \lambda\, r(\mathbf{x})\,\psi \qquad \text{in } \Omega\subset\mathbb{R}^d ,
```

which contains the time-independent Schrödinger equation (TISE) ``\hat H\psi=\varepsilon\psi`` for

```math
\hat H = -\nabla\cdot\big(\mathbf{p}\,\nabla\big) + V,\qquad
\mathbf{p} = \frac{\hbar^2}{2}\,\boldsymbol{\mathcal{M}}^{-1},\quad q = V,\quad r = 1,\quad \lambda=\varepsilon .
```

| Physical situation | Kinetic coefficient ``\mathbf p`` | How to pass it |
|---|---|---|
| One particle of mass ``m`` | ``\hbar^2/2m`` (``=1/2`` for an electron in a.u.) | `p(x) = 0.5` |
| Effective mass ``m^*(x)`` in a heterostructure (BenDaniel–Duke) | ``\hbar^2/2m^*(x)`` | `p(x) = 0.5 / mstar(x[1])` |
| Two 1D particles of masses ``M_1, M_2`` on a 2D mesh | ``\mathrm{diag}(\hbar^2/2M_1,\ \hbar^2/2M_2)`` | `p₁`, `p₂` (see [Low-level interface](@ref)) |

The coefficients may be discontinuous. ``r`` is kept general because it appears
after changes of variables.

## Weak formulation

Physical states have finite norm and finite kinetic energy, i.e. they belong to the
Sobolev space ``H^1(\Omega)``. Multiplying by a test function ``v`` and integrating
by parts gives the weak problem: find ``(\lambda, u)`` with ``u \neq 0`` such that

```math
a(u,v) = \lambda\, b(u,v)\quad \forall v,\qquad
a(u,v) = \int_\Omega \big[(\mathbf p\nabla u)\cdot\nabla v^* + q\,u\,v^*\big]\,d\Omega,\qquad
b(u,v) = \int_\Omega r\,u\,v^*\,d\Omega .
```

For regular ``u``, ``a(u,v)=\langle v|\hat H|u\rangle`` and ``b(u,v)=\langle v|u\rangle``.
Only first derivatives appear. The forms are built by [`bilineal_forms`](@ref).

**Boundary conditions.** The package imposes homogeneous Dirichlet conditions
``\psi|_{\partial\Omega}=0`` (impenetrable walls): they are built into the spaces
(``U=V=H^1_0``) and the boundary term of the integration by parts vanishes.

**Interfaces.** If ``\mathbf p`` or ``V`` jump across a surface, the weak form
automatically enforces continuity of ``\psi`` and of the flux
``\hat n\cdot\mathbf p\nabla\psi``, i.e. the BenDaniel–Duke conditions
``[\psi]=0``, ``[\partial_n\psi/m^*]=0``. Nothing has to be imposed by hand.

## Discretization

The domain is divided into a mesh of cells (segments in 1D, quadrilaterals or
triangles in 2D) of maximum size ``h``. In each cell the wave function is a
polynomial of degree ``\ell`` determined by its values at the nodes (Lagrange
elements). The global basis functions ``\phi_j`` are continuous, piecewise
polynomial and nonzero only on the cells around node ``j``. Writing
``u_h=\sum_j\Phi_j\phi_j`` turns the weak problem into the **generalized
eigenvalue problem**

```math
\mathbb{A}\,\boldsymbol\Phi = \varepsilon\,\mathbb{B}\,\boldsymbol\Phi,\qquad
\mathbb{A}_{ij} = a(\phi_j,\phi_i),\quad \mathbb{B}_{ij} = b(\phi_j,\phi_i).
```

``\mathbb A`` (Hamiltonian = stiffness + potential matrix) and ``\mathbb B``
(overlap or mass matrix) are sparse and Hermitian, and ``\mathbb B`` is positive
definite. The components of an eigenvector are the values of the eigenfunction at
the nodes, ``\Phi_j=u_h(\mathbf x_j)``. The integrals are computed cell by cell with
Gauss quadrature exact for polynomials up to a chosen degree ``k_q``.

## What is guaranteed

**Upper bounds.** The discrete space is a subspace of ``H^1_0``. By the min–max
principle, every computed eigenvalue is an upper bound to the exact one,

```math
\varepsilon_n \le \varepsilon_{n,h},
```

and the bounds decrease monotonically when the mesh is refined. This holds when
the integrals are exact; see [Accuracy and best practices](@ref) for the effect of
an insufficient quadrature degree.

**Convergence.** For eigenfunctions that are smooth on every cell,

```math
0 \le \frac{\varepsilon_{n,h}-\varepsilon_n}{\varepsilon_n} \le C_n\,h^{2\ell},
\qquad \|u_n-u_{n,h}\|_{H^1} = \mathcal O(h^{\ell}),\qquad \|u_n-u_{n,h}\|_{L^2} = \mathcal O(h^{\ell+1}).
```

Eigenvalues converge twice as fast as eigenfunctions. The constant ``C_n`` grows
with ``n`` because excited states oscillate faster. In terms of the number of
unknowns ``N_h \propto h^{-d}``, the eigenvalue error decays as ``N_h^{-2\ell/d}``.

**Discrete spectrum.** Because the domain is a finite box, the spectrum is
discrete. The continuum of an open system appears as a *pseudocontinuum* of box
states whose density increases with the size of ``\Omega``.

## The eigensolver

Only a few eigenpairs of matrices with ``10^3``–``10^5`` rows are needed, so the
package uses the implicitly restarted Arnoldi method of ARPACK in
**shift-and-invert** mode:

```math
(\mathbb A-\sigma\mathbb B)^{-1}\mathbb B\,\boldsymbol\Phi = \nu\,\boldsymbol\Phi,
\qquad \nu = \frac{1}{\varepsilon-\sigma}.
```

The eigenvalues ``\varepsilon`` closest to the shift ``\sigma`` become the largest
``\nu`` and converge first. ``\mathbb A-\sigma\mathbb B`` is factorized once by a
sparse direct method (`SuiteSparse`), and each iteration costs two triangular
solves.

**Choosing ``\sigma``**: a value below the ground state returns the ``n_{\rm ev}``
lowest levels; a value inside the spectrum returns the levels around that energy
(e.g. a band, or states near a threshold).

ARPACK computes few eigenpairs of large sparse matrices; LAPACK (dense) is preferable
only for small matrices or when the full spectrum is needed. Both rely on BLAS.

## Further reading

- J. Sun and A. Zhou, *Finite Element Methods for Eigenvalue Problems* (CRC, 2016).
- I. Babuška and J. Osborn, *Eigenvalue problems*, in Handbook of Numerical Analysis Vol. II (1991).
- G. Strang and G. J. Fix, *An Analysis of the Finite Element Method* (2008).
- R. B. Lehoucq, D. C. Sorensen and C. Yang, *ARPACK Users' Guide* (SIAM, 1998).
- F. Verdugo and S. Badia, *The software design of Gridap*, Comput. Phys. Commun. **276**, 108341 (2022).
