function norm_l2(𝜳::CellField,dΩ::Gridap.CellData.GenericMeasure)
    return sqrt(real(sum(∫(𝜳'*𝜳)*dΩ)));
end

function ortho_check_v2(ϕ::Vector{CellField},dΩ::Gridap.CellData.GenericMeasure)
    nev::Integer=length(ϕ)
    OrthoVector=zeros(Float64,round(Int,(nev^2-nev)/2));
    index::Int64=1
    for i in 2:nev
        for j in 1:(i-1)
            OrthoVector[index]=abs(sum(∫(ϕ[j]'*ϕ[i])*dΩ))
            index+=1
        end
    end
    return OrthoVector;
end