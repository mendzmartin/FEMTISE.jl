struct DefaultBinEigenProblem{T}
    ϵ::Vector{ComplexF64}
    ϕ::Vector{CellField}
    r::T
end

struct DefaultBinEigenProblemReducedDensity{T}
    ϵ::Vector{ComplexF64}
    ϕ::Vector{CellField}
    r::T
    rhoDOF1::Vector{Float64}
    rhoDOF2::Vector{Float64}
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