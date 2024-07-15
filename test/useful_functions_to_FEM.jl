export norm_l2
"""
    norm_l2(𝜳,dΩ)

# Aim
- Compute the L2 norm of a field

# Arguments
- `𝜳::CellField`: field to compute the norm
- `dΩ::Gridap.CellData.GenericMeasure`: measure of FE grid

# Example
```julia
using Gridap
using FEMTISE
using LinearAlgebra

# Define measure
Ω=1.0:0.1:10.0
dΩ=Gridap.CellData.CellData(Ω)
# Define FE Spaces
USpace=fe_space(1,1,Ω)
VSpace=fe_space(1,1,Ω)
# Define functions
p(x)=1.0
q(x)=0.0
r(x)=1.0
# Compute eigenvalues and eigenvectors
ϵ,ϕ=eigen_values_and_eigen_vectors(p,q,r,dΩ,USpace,VSpace)
# Compute norm of eigenstates
norm_vector=[norm_l2(ϕ[i],dΩ) for i in eachindex(ϕ)]
```
"""
function norm_l2(𝜳::CellField,dΩ::Gridap.CellData.GenericMeasure)
    return sqrt(real(sum(∫(𝜳'*𝜳)*dΩ)));
end

export orthogonality_check
"""
    orthogonality_check(ϕ,dΩ)

# Aim
- Compute the orthogonality of the eigenstates  

# Arguments
- `ϕ::Vector{CellField}`: eigenstates
- `dΩ::Gridap.CellData.GenericMeasure`: measure of FE grid

# Example
```julia
using Gridap
using FEMTISE
using LinearAlgebra

# Define measure
Ω=1.0:0.1:10.0
dΩ=Gridap.CellData.CellData(Ω)
# Define FE Spaces
USpace=fe_space(1,1,Ω)
VSpace=fe_space(1,1,Ω)
# Define functions
p(x)=1.0
q(x)=0.0
r(x)=1.0
# Compute eigenvalues and eigenvectors
ϵ,ϕ=eigen_values_and_eigen_vectors(p,q,r,dΩ,USpace,VSpace)
# Compute orthogonality of eigenstates
orthogonality_vector=orthogonality_check(ϕ,dΩ)
```
"""
function orthogonality_check(ϕ::Vector{CellField},dΩ::Gridap.CellData.GenericMeasure)
    nev::Integer=length(ϕ);
    orthogonality_vector=zeros(Float64,round(Int,(nev^2-nev)/2));
    index::Int=1;
    for i in 2:nev
        for j in 1:(i-1)
            orthogonality_vector[index]=abs(sum(∫(ϕ[j]'*ϕ[i])*dΩ));
            index+=1;
        end
    end
    return orthogonality_vector;
end

export multifield_orthogonality_check
"""
    multifield_orthogonality_check(ϕ,dΩ,TrialSpace)

# Aim
- Compute the orthogonality of the eigenstates

# Arguments
- `ϕ::Vector{CellField}`: eigenstates
- `dΩ::Gridap.CellData.GenericMeasure`: measure of FE grid
- `TrialSpace::FESpace`: trial FE Space
"""
function multifield_orthogonality_check(ϕ::Vector{CellField},dΩ::Gridap.CellData.GenericMeasure,TrialSpace::FESpace)
    nev::Integer=length(ϕ);
    orthogonality_vector=zeros(Float64,round(Int,(nev^2-nev)/2));
    index::Int=1;
    for i in 2:nev
        ϕᵢ=interpolate_everywhere(ϕ[i],TrialSpace);
        ϕ¹ᵢ,ϕ²ᵢ=ϕᵢ;
        for j in 1:(i-1)
            ϕⱼ=interpolate_everywhere(ϕ[j],TrialSpace);
            ϕ¹ⱼ,ϕ²ⱼ=ϕⱼ;
            orthogonality_vector[index]=abs(sum(∫(ϕ¹ⱼ'*ϕ¹ᵢ)*dΩ)+sum(∫(ϕ²ⱼ'*ϕ²ᵢ)*dΩ));
            index+=1;
        end
    end
    return orthogonality_vector;
end

export eigenstates_normalization
"""
    eigenstates_normalization(ϕ,dΩ)

# Aim
- Compute the normalization of the eigenstates

# Arguments
- `ϕ::Vector{CellField}`: eigenstates
- `dΩ::Gridap.CellData.GenericMeasure`: measure of FE grid

# Example
```julia
using Gridap
using FEMTISE
using LinearAlgebra

# Define measure
Ω=1.0:0.1:10.0
dΩ=Gridap.CellData.CellData(Ω)
# Define FE Spaces
USpace=fe_space(1,1,Ω)
VSpace=fe_space(1,1,Ω)
# Define functions
p(x)=1.0
q(x)=0.0
r(x)=1.0
# Compute eigenvalues and eigenvectors
ϵ,ϕ=eigen_values_and_eigen_vectors(p,q,r,dΩ,USpace,VSpace)
# Compute normalization of eigenstates
nom_vector=eigenstates_normalization(ϕ,dΩ)
```
"""
function eigenstates_normalization(ϕ::Vector{CellField},dΩ::Gridap.CellData.GenericMeasure)
    return [norm_l2(ϕ[i],dΩ) for i in eachindex(ϕ)]
end

export multifield_eigenstates_normalization
"""
    multifield_eigenstates_normalization(ϕ,dΩ,TrialSpace)

# Aim
- Compute the normalization of the eigenstates

# Arguments
- `ϕ::Vector{CellField}`: eigenstates
- `dΩ::Gridap.CellData.GenericMeasure`: measure of FE grid
- `TrialSpace::FESpace`: trial FE Space
"""
function multifield_eigenstates_normalization(ϕ::Vector{CellField},dΩ::Gridap.CellData.GenericMeasure,TrialSpace::FESpace)
    nom_vector=zeros(Float64,length(ϕ));
    for i in eachindex(ϕ)
        ϕᵢ=interpolate_everywhere(ϕ[i],TrialSpace);
        ϕ¹ᵢ,ϕ²ᵢ=ϕᵢ;
        nom_vector[i]=norm_l2(ϕ¹ᵢ,dΩ)+norm_l2(ϕ²ᵢ,dΩ);
    end
    return nom_vector;
end