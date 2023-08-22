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

# Bilineals forms for eigenvalues 1D or 2D (equal masses)
export bilineal_forms
function bilineal_forms(p::Function,q::Function,r::Function,dΩ::Gridap.CellData.GenericMeasure)
    a(u,v) = ∫(p*(∇(v)⋅∇(u))+q*v*u)*dΩ
    b(u,v) = ∫(r*u*v)*dΩ
    return a,b;
end

export fe_spaces
"""
    fe_spaces(model,reff,grid_type; <keyword arguments>)

Create finite element (FE) spaces (Trial and Test spaces).

...
# Arguments
- `BC_type::String="FullDirichlet"`: the type of boundary condition.
- `TypeData::Type=ComplexF64`: the type of data to define FE spaces.
- `conf_type::Symbol=:H1`: the regularity of the interpolation at the boundaries of cells in the mesh. (e.g.:L2,:H1,:C0,:Hgrad,)
...
"""
function fe_spaces(model,reff::Tuple,grid_type::String;BC_type::String="FullDirichlet",TypeData::Type=ComplexF64,conf_type::Symbol=:H1)
    BC_values,BC_tags = make_boundary_conditions(grid_type,BC_type,TypeData);
    VSpace=TestFESpace(model,reff;vector_type=Vector{TypeData},conformity=conf_type,dirichlet_tags=BC_tags);
    USpace=TrialFESpace(VSpace,BC_values);
    return VSpace,USpace;
end