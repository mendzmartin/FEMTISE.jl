# the triangulation and integration aproximated Lebesgue measure
export measures
function measures(model,degree,tags_boundary)
    Ω=Triangulation(model); dΩ=Measure(Ω,degree);
    Γ=BoundaryTriangulation(model,tags=tags_boundary); dΓ=Measure(Γ,degree)
    return Ω,dΩ,Γ,dΓ;
end

"""
Reference FE space definition
    method  = lagrangian
    type    = Float64,ComplexF64
    order   = 1,2,3
"""
export reference_FEspaces
function reference_FEspaces(method,type,order)
    reff=ReferenceFE(method,type,order);
    return reff;
end

"""
    if grid_type == "simple_line" => params=(path,name,dom,MeshSize)
    if grid_type == "simple_rectangle_v1" => params=(path,name,dom,MeshSize,quad_state)
    if grid_type == "simple_rectangle_v2" => params=(path,name,side_x,side_y,lc,numNodesHE,quad_state,structured_mesh,bumpFactor)
    if grid_type == "Cartesian2D" => params=(dom,n)
"""
export model_FE
function model_FE(grid_type,params)
    model=make_model(grid_type,params)
    return model;
end

"""
    if dimension=="1D" ⇒ dom=(x₁,x₂); Δr=Δx; n=nx
    if dimension=="2d" ⇒ dom=(x₁,x₂,y₁,y₂); Δr=(Δx,Δy); n=(nx,ny)
"""
export space_coord
function space_coord(dom,Δr,n;dimension="2D")
    if (dimension=="1D")
        r=[dom[1]+Δr*(i-1) for i in 1:n[1]];
        pts=[Point(r[i]) for i in 1:n[1]];
    else if (dimension=="2D")
        r=([dom[1]+Δr[1]*(i-1) for i in 1:n[1]],[dom[3]+Δr[2]*(i-1) for i in 1:n[2]]);
        pts=[Point(r[1][i],r[2][j]) for i in 1:n[1] for j in 1:n[2]];
    end
    return r,pts;
end

"""
    dom=(x₁,x₂,y₁,y₂)
    n=(nx,ny)
"""
export FESpaces
function FESpaces(model,reffe,grid_type;BC_type="FullDirichlet",TypeData=ComplexF64)
    BC_values,BC_tags = make_boundary_conditions(grid_type,BC_type,TypeData);
    VSpace=TestFESpace(model,reffe;vector_type=Vector{TypeData},conformity=:H1,dirichlet_tags=BC_tags);
    USpace=TrialFESpace(VSpace,BC_values);
    return VSpace,USpace;
end