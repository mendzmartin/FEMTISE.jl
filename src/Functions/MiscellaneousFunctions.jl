# the triangulation and integration aproximated Lebesgue measure
export measures
function measures(model,degree::Integer,tags_boundary)
    Ω=Triangulation(model); dΩ=Measure(Ω,degree);
    Γ=BoundaryTriangulation(model,tags=tags_boundary); dΓ=Measure(Γ,degree)
    return Ω,dΩ,Γ,dΓ;
end

export space_coord
"""
    if dimension=="1D" ⇒ dom=(x₁,x₂); Δr=Δx; n=nx
    if dimension=="2d" ⇒ dom=(x₁,x₂,y₁,y₂); Δr=(Δx,Δy); n=(nx,ny)
"""
function space_coord(dom::Tuple,Δr,n;dimension::String="2D")
    if (dimension=="1D")
        r=[dom[1]+Δr*(i-1) for i in 1:n];
        pts=[Point(r[i]) for i in 1:n];
    elseif (dimension=="2D")
        r=([dom[1]+Δr[1]*(i-1) for i in 1:n[1]],[dom[3]+Δr[2]*(i-1) for i in 1:n[2]]);
        pts=[Point(r[1][i],r[2][j]) for i in 1:n[1] for j in 1:n[2]];
    end
    return r,pts;
end

# Formas bilineales para problema de autovalores
export bilineal_forms
function bilineal_forms(p::Function,q::Function,r::Function,dΩ::Gridap.CellData.GenericMeasure)
    a(u,v) = ∫(p*(∇(v)⋅∇(u))+q*v*u)*dΩ;
    b(u,v) = ∫(r*u*v)*dΩ;
    return a,b;
end

export FESpaces
"""
    dom=(x₁,x₂,y₁,y₂)
    n=(nx,ny)
"""
function FESpaces(model,reff::Tuple,grid_type::String;BC_type::String="FullDirichlet",TypeData::Type=ComplexF64)
    BC_values,BC_tags = make_boundary_conditions(grid_type,BC_type,TypeData);
    VSpace=TestFESpace(model,reff;vector_type=Vector{TypeData},conformity=:H1,dirichlet_tags=BC_tags);
    USpace=TrialFESpace(VSpace,BC_values);
    return VSpace,USpace;
end