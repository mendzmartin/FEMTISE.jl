"""
    eigen_values_and_eigen_vectors(p,q,r,dΩ,USpace,VSpace; <keyword arguments>)

# Aim
- Compute matrix eigen problem and return eigenvalues and eigenvectors

# Arguments
- `p,q,r::Function`: functions from Sturm Liouville formulation
- `dΩ::Gridap.CellData.GenericMeasure`: measure of FE grid
- `USpace::FESpace`: trial FE Space
- `VSpace::FESpace`: test FE Space
- `params::Tuple=(nev,tol,maxiter,explicittransform,sigma)`: params to Arpack eigs function.
  - `nev::Integer=10`: quantity of eigendata to calculate
  - `tol::Float64=10e-9`: relative tolerance for convergence of Ritz values
  - `maxiter::Integer=100`: maximum number of iterations
  - `explicittransform::Symbol=:none`: shift and invert should be explicitly invoked in julia code
  - `sigma::Float64=1.0`: the level shift used in inverse iteration.

# Example
```julia
using Gridap
using FEMTISE
using LinearAlgebra

# Define functions
p(x)=1.0
q(x)=0.0
r(x)=1.0
# Define measure
Ω=1.0:0.1:10.0
dΩ=Gridap.CellData.CellData(Ω)
# Define FE Spaces
USpace=fe_space(1,1,Ω)
VSpace=fe_space(1,1,Ω)
# Compute eigenvalues and eigenvectors
ϵ,ϕ=eigen_values_and_eigen_vectors(p,q,r,dΩ,USpace,VSpace)
```
"""
function eigen_values_and_eigen_vectors(p::Function,q::Function,r::Function,
  dΩ::Gridap.CellData.GenericMeasure,USpace::FESpace,VSpace::FESpace;
  params=(10::Integer,10e-9::Float64,100,:none,1.0::Float64))
    # Define bilinear forms
    a,b=bilineal_forms(p,q,r,dΩ)
    # solve eigenvalue problem
    prob=eigen_problem(a,b,USpace,VSpace;nev=params[1],tol=params[2],maxiter=params[3],explicittransform=params[4],sigma=params[5]);
    ϵ,ϕ=solve(prob);
    return ϵ,ϕ;
end

"""
    eigen_values_and_eigen_vectors(p₁,p₂,q,r,dΩ,USpace,VSpace; <keyword arguments>)

# Aim
- Compute matrix eigen problem and return eigenvalues and eigenvectors

# Arguments
- `p₁,p₂,q,r::Function`: functions from Sturm Liouville formulation
- `dΩ::Gridap.CellData.GenericMeasure`: measure of FE grid
- `USpace::FESpace`: trial FE Space
- `VSpace::FESpace`: test FE Space
- `params::Tuple=(nev,tol,maxiter,explicittransform,sigma)`: params to Arpack eigs function.
  - `nev::Integer=10`: quantity of eigendata to calculate
  - `tol::Float64=10e-9`: relative tolerance for convergence of Ritz values
  - `maxiter::Integer=100`: maximum number of iterations
  - `explicittransform::Symbol=:none`: shift and invert should be explicitly invoked in julia code
  - `sigma::Float64=1.0`: the level shift used in inverse iteration.

# Example
```julia
using Gridap
using FEMTISE
using LinearAlgebra

# Define functions
p₁(x)=1.0
p₂(x)=0.0
q(x)=0.0
r(x)=1.0
# Define measure
Ω=1.0:0.1:10.0
dΩ=Gridap.CellData.CellData(Ω)
# Define FE Spaces
USpace=fe_space(1,1,Ω)
VSpace=fe_space(1,1,Ω)
# Compute eigenvalues and eigenvectors
ϵ,ϕ=eigen_values_and_eigen_vectors(p₁,p₂,q,r,dΩ,USpace,VSpace)
```
"""
function eigen_values_and_eigen_vectors(p₁::Function,p₂::Function,q::Function,r::Function,
  dΩ::Gridap.CellData.GenericMeasure,USpace::FESpace,VSpace::FESpace;
  params=(10::Integer,10e-9::Float64,100,:none,1.0::Float64))
  # Define bilinear forms
  a,b=bilineal_forms(p₁,p₂,q,r,dΩ)
  # solve eigenvalue problem
  prob=eigen_problem(a,b,USpace,VSpace;nev=params[1],tol=params[2],maxiter=params[3],explicittransform=params[4],sigma=params[5]);
  ϵ,ϕ=solve(prob);
  return ϵ,ϕ;
end