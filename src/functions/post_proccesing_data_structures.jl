"""
    DefaultBinEigenProblem

Results reloaded from the `bin eigen` output format: eigenvalues `ϵ`, eigenfunctions
`ϕ` sampled on the output grid (one column per state) and grid coordinates `r`.
"""
struct DefaultBinEigenProblem{T}
    ϵ::Vector{ComplexF64}
    ϕ::Matrix{ComplexF64}
    r::T
end

"""
    DefaultBinEigenProblemReducedDensity

As [`DefaultBinEigenProblem`](@ref), plus the reduced densities `rhoDOF1` and
`rhoDOF2` (one column per state) of a two-dimensional problem.
"""
struct DefaultBinEigenProblemReducedDensity{T}
    ϵ::Vector{ComplexF64}
    ϕ::Matrix{ComplexF64}
    r::T
    rhoDOF1::Matrix{Float64}
    rhoDOF2::Matrix{Float64}
end

struct DefaultJLD2EigenProblem{T1,T2}
    ϵ::Vector{ComplexF64}
    ϕ::Vector{CellField}
    r::T1
    pts::T2
end

struct DefaultJLD2EigenProblemReducedDensity{T1,T2}
    ϵ::Vector{ComplexF64}
    ϕ::Vector{CellField}
    r::T1
    pts::T2
    rhoDOF1::Vector{Float64}
    rhoDOF2::Vector{Float64}
end

struct DefaultJLD2AllEigenProblem{T1,T2,T3,T4,T5}
    ϵ::Vector{ComplexF64}
    ϕ::Vector{CellField}
    r::T1
    pts::T2
    Ω::T3
    dΩ::Gridap.CellData.GenericMeasure
    Γ::T4
    dΓ::Gridap.CellData.GenericMeasure
    USpace::FESpace
    VSpace::FESpace
    model::T5
end

struct DefaultJLD2AllEigenProblemReducedDensity{T1,T2,T3,T4,T5}
    ϵ::Vector{ComplexF64}
    ϕ::Vector{CellField}
    r::T1
    pts::T2
    Ω::T3
    dΩ::Gridap.CellData.GenericMeasure
    Γ::T4
    dΓ::Gridap.CellData.GenericMeasure
    USpace::FESpace
    VSpace::FESpace
    model::T5
    rhoDOF1::Vector{Float64}
    rhoDOF2::Vector{Float64}
end

struct AnalysisParams
    ϵ_matrix::Matrix{ComplexF64}
    λvector::Vector{ComplexF64}
end