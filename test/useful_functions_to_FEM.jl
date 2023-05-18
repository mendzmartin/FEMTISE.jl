function norm_l2(𝜳::CellField,dΩ::Gridap.CellData.GenericMeasure)
    return sqrt(real(sum(∫(𝜳'*𝜳)*dΩ)));
end

function ortho_check_v2(ϕ::Vector{CellField},TrialSpace::FESpace,dΩ::Gridap.CellData.GenericMeasure)
    nev::Integer=length(ϕ)
    OrthoVector=zeros(Float64,round(Int,(nev^2-nev)/2));
    index::Int64=1
    for i in 2:nev
        ϕᵢ=interpolate_everywhere(ϕ[i],TrialSpace);
        for j in 1:(i-1)
            ϕⱼ=interpolate_everywhere(ϕ[j],TrialSpace);
            OrthoVector[index]=abs(sum(∫(ϕⱼ'*ϕᵢ)*dΩ))
            index+=1
        end
    end
    return OrthoVector;
end