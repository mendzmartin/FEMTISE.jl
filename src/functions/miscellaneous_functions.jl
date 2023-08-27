"""
    measures(model,degree,tags_boundary)
    
# Aim
The triangulation and integration aproximated Lebesgue measures

# Arguments
- `model`: FE grid model.
- `degree::Integer`: degree of quadrature rule to use in the cells of triangulation.
- `tags_boundary`: tags of boundary conditions.
"""
function measures(model,degree::Integer,tags_boundary)
    Ω=Triangulation(model); dΩ=Measure(Ω,degree);
    Γ=BoundaryTriangulation(model,tags=tags_boundary); dΓ=Measure(Γ,degree)
    return Ω,dΩ,Γ,dΓ;
end

"""
    space_coord(dom,Δr,n;<keyword arguments>)

# Aim
- Returns coordinate vector (r) and discrete points (pts) for 1D or 2D spaces.
  - if dimension=="1D" ⇒ dom=(x₁,x₂); Δr=Δx; n=nx
  - if dimension=="2d" ⇒ dom=(x₁,x₂,y₁,y₂); Δr=(Δx,Δy); n=(nx,ny)

# Arguments
- `dom::Tuple`: FE cartesian domain.
- `Δr`:: discretization of FE space.
- `n`:: number of FE in each direction.
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

"""
    bilineal_forms(p,q,r,dΩ)
    
# Aim
- Returns bilineals forms (a(u,v) and b(u,v)) for eigenvalues 1D or 2D (equal masses)

# Arguments
- `p::Function`: kinetic energy function from Sturm Liouville_Formulation.
- `q::Function`: potential energy function from Sturm Liouville_Formulation.
- `r::Function`: weight or density function from Sturm Liouville_Formulation.
- `dΩ::Gridap.CellData.GenericMeasure`: differential FE domain
"""
function bilineal_forms(p::Function,q::Function,r::Function,dΩ::Gridap.CellData.GenericMeasure)
    a(u,v) = ∫(p*(∇(v)⋅∇(u))+q*v*u)*dΩ
    b(u,v) = ∫(r*u*v)*dΩ
    return a,b;
end
# Bilineals forms for eigenvalues 2D (different masses)
function bilineal_forms(p₁::Function,p₂::Function,q::Function,r::Function,dΩ::Gridap.CellData.GenericMeasure)
    e₁ = VectorValue(1.0,0.0)
    e₂ = VectorValue(0.0,1.0)
    a(u,v) = ∫(p₁*((∇(v)⋅e₁)⋅(∇(u)⋅e₁))+p₂*((∇(v)⋅e₂)⋅(∇(u)⋅e₂))+q*v*u)*dΩ;
    b(u,v) = ∫(r*u*v)*dΩ;
    return a,b;
end

"""
    fe_spaces(model,reff,grid_type; <keyword arguments>)

# Aim
- Create finite element (FE) spaces (Trial and Test spaces).

# Arguments
- `BC_type::String="FullDirichlet"`: the type of boundary condition.
- `TypeData::Type=ComplexF64`: the type of data to define FE spaces.
- `conf_type::Symbol=:H1`: the regularity of the interpolation at the boundaries of cells in the mesh. (e.g.:L2,:H1,:C0,:Hgrad,)
"""
function fe_spaces(model,reff::Tuple,grid_type::String;
    BC_type::String="FullDirichlet",TypeData::Type=ComplexF64,conf_type::Symbol=:H1)

    BC_values,BC_tags = make_boundary_conditions(grid_type,BC_type,TypeData);
    VSpace=TestFESpace(model,reff;vector_type=Vector{TypeData},conformity=conf_type,dirichlet_tags=BC_tags);
    USpace=TrialFESpace(VSpace,BC_values);
    return VSpace,USpace;
end

function fe_spaces(model,reff::Tuple;
    BC_data::Tuple=(nothing,nothing),BC_type::String="Neumann",
    TypeData::Type=ComplexF64,conf_type::Symbol=:H1)
    
    if (BC_type == "Neumann")
        VSpace=TestFESpace(model,reff;vector_type=Vector{TypeData},conformity=conf_type);
        USpace=TrialFESpace(VSpace);
    elseif (BC_type == "Dirichlet")
        BC_values,BC_tags=BC_data
        VSpace=TestFESpace(model,reff;vector_type=Vector{TypeData},conformity=conf_type,dirichlet_tags=BC_tags);
        USpace=TrialFESpace(VSpace,BC_values);
    end
    return VSpace,USpace;
end

