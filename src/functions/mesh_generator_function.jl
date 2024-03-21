#=
    References:
        https://gitlab.onelab.info/gmsh/gmsh/blob/master/api/gmsh.jl
=#
function make_model(grid_type::String,params::Tuple)
    GridapGmsh.gmsh.initialize()
    GridapGmsh.gmsh.option.setNumber("General.Terminal",1)

    # 1D line grid
    if (grid_type == "simple_line")
        path,name,dom,MeshSize=params
        GridapGmsh.gmsh.model.add(name)

        GridapGmsh.gmsh.model.geo.addPoint(dom[1],0,0,MeshSize,1)    # 1 punto vértice izquierdo
        GridapGmsh.gmsh.model.geo.addPoint(dom[2],0,0,MeshSize,2)    # 2 punto vértice derecho

        GridapGmsh.gmsh.model.geo.addLine(1,2,1)   # linea que une puntos 1 y 2
        GridapGmsh.gmsh.model.geo.synchronize()
        
        GridapGmsh.gmsh.model.geo.addPhysicalGroup(0,[1],3)        # grupo formado por punto izquierdo
        GridapGmsh.gmsh.model.setPhysicalName(0,3,"left_point")    # le damos nombre al grupo
        GridapGmsh.gmsh.model.geo.synchronize()

        GridapGmsh.gmsh.model.geo.addPhysicalGroup(0,[2],4)        # grupo formado por punto derecho
        GridapGmsh.gmsh.model.setPhysicalName(0,4,"right_point")   # le damos nombre al grupo
        GridapGmsh.gmsh.model.geo.synchronize()
        
        GridapGmsh.gmsh.model.geo.addPhysicalGroup(1,[1],1)    # grupo formado por linea
        GridapGmsh.gmsh.model.setPhysicalName(1,1,"segment")   # le damos nombre a la linea
        GridapGmsh.gmsh.model.geo.synchronize()

        GridapGmsh.gmsh.model.mesh.generate(1)
        GridapGmsh.gmsh.write(path*name*".msh")
        GridapGmsh.gmsh.finalize()
        model=GmshDiscreteModel(path*name*".msh")
    # 2D simple rectangle grid
    elseif (grid_type == "simple_rectangle_v1")

        path,name,dom,MeshSize,quad_state=params

        # GridapGmsh.gmsh.option.setNumber(name, value)
        GridapGmsh.gmsh.option.setNumber("General.Terminal",0)
        GridapGmsh.gmsh.model.add(name)

        # creamos puntos vértice - GridapGmsh.gmsh.model.geo.addPoint(x,y,z,meshSize=0.,tag=-1)
        GridapGmsh.gmsh.model.geo.addPoint(dom[1],dom[3],0,MeshSize,1) # 1 vertice inferior izquierdo
        GridapGmsh.gmsh.model.geo.addPoint(dom[2],dom[3],0,MeshSize,2) # 2 vértice inferior derecho
        GridapGmsh.gmsh.model.geo.addPoint(dom[2],dom[4],0,MeshSize,3) # 3 vértice superior derecho
        GridapGmsh.gmsh.model.geo.addPoint(dom[1],dom[4],0,MeshSize,4) # 4 vértice superior izquierdo
        # creamos lineas de unión entre vértices - GridapGmsh.gmsh.model.geo.addLine(startTag,endTag,tag=-1)
        GridapGmsh.gmsh.model.geo.addLine(1,2,1) # 1 linea inferior
        GridapGmsh.gmsh.model.geo.addLine(2,3,2) # 2 línea lateral derecha
        GridapGmsh.gmsh.model.geo.addLine(3,4,3) # 3 linea superior
        GridapGmsh.gmsh.model.geo.addLine(4,1,4) # 4 linea lateral izquierda
        # creamos curva de unión entre lineas
        GridapGmsh.gmsh.model.geo.addCurveLoop([1,2,3,4],100) # the rectangle
        GridapGmsh.gmsh.model.geo.synchronize()

        # make the surface - GridapGmsh.gmsh.model.geo.addPlaneSurface(wireTags,tag=-1)
        GridapGmsh.gmsh.model.geo.addPlaneSurface([100],101) # the surface
        GridapGmsh.gmsh.model.geo.synchronize()

        #=
            creamos grupos para definir condiciones de bordes
            GridapGmsh.gmsh.model.geo.addPhysicalGroup(dim,tags,tag=-1,name="")
            GridapGmsh.gmsh.model.setPhysicalName(dim,tag,name)
        =#
        # creamos grupo físico de puntos vértices
        GridapGmsh.gmsh.model.geo.addPhysicalGroup(0,[1,2,3,4],300)  # grupo 0D formado por cuatro puntos
        GridapGmsh.gmsh.model.setPhysicalName(0,300,"ext_points")
        GridapGmsh.gmsh.model.geo.synchronize()

        # creamos grupo físico de curva de lineas externa
        # GridapGmsh.gmsh.model.geo.addPhysicalGroup(1,[100],301)      # de esta forma no funciona!!
        GridapGmsh.gmsh.model.geo.addPhysicalGroup(1,[1,2,3,4],301)    # grupo 1D formado por cuatro lineas
        GridapGmsh.gmsh.model.setPhysicalName(1,301,"ext_lines")
        GridapGmsh.gmsh.model.geo.synchronize()

        # creamos grupo físico con superficie interna del rectángulo
        GridapGmsh.gmsh.model.addPhysicalGroup(2,[101],302)
        GridapGmsh.gmsh.model.setPhysicalName(2,302,"surface") # grupo 2D formado por una superficie
        GridapGmsh.gmsh.model.geo.synchronize()

        if quad_state
            # GridapGmsh.gmsh.model.mesh.setRecombine(dim,tag,angle=45.)
            GridapGmsh.gmsh.model.mesh.setRecombine(2,302) # for 2D quadrilaterals
        end

        # GridapGmsh.gmsh.model.mesh.generate(dim=3)
        GridapGmsh.gmsh.model.mesh.generate(2)
        GridapGmsh.gmsh.write(path*name*".msh")
        GridapGmsh.gmsh.finalize()
        model=GmshDiscreteModel(path*name*".msh")
    #= ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        Bloque de código agregado para crear un nuevo test que permita
        utilizar elementos cuadriláteros y mallas estructuradas.
    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ =#
    elseif grid_type == "simple_rectangle_v2"   # simple square Quad and triangle elements
        
        # parámetros de entrada
        path,name,side_x,side_y,lc,numNodesHE,quad_state,structured_mesh,bumpFactor=params
        # Asignamos número de nodos en los bordes horizontal y vertical
        numNodesHE_hor,numNodesHE_ver = numNodesHE

        if (quad_state == true)
            GridapGmsh.gmsh.option.setNumber("Mesh.Algorithm", 5) # delquad
            GridapGmsh.gmsh.option.setNumber("Mesh.RecombineAll", 1)
        end
        GridapGmsh.gmsh.model.add(name)

        # first we build the rectangular boundary
        lc_x = lc
        lc_y = lc*side_y/side_x
        # GridapGmsh.gmsh.model.geo.addPoint(x,y,z,meshSize=0.,tag=-1)
        GridapGmsh.gmsh.model.occ.addPoint(0,0,0,lc_x,1)           # 1 vertice inferior izq
        GridapGmsh.gmsh.model.occ.addPoint(side_x,0,0,lc_x,2)      # 2 vértice inferior der
        GridapGmsh.gmsh.model.occ.addPoint(side_x,side_y,0,lc_y,3) # 3 vértice superior der
        GridapGmsh.gmsh.model.occ.addPoint(0,side_y,0,lc_y,4)      # 4 vértice superior izq

        # make the square boundary
        GridapGmsh.gmsh.model.occ.addLine(1, 2, 1) # 1 linea inferior
        GridapGmsh.gmsh.model.occ.addLine(2, 3, 2) # 2 linea lateral der 
        GridapGmsh.gmsh.model.occ.addLine(3, 4, 3) # 3 linea superior
        GridapGmsh.gmsh.model.occ.addLine(4, 1, 4) # 4 linea lateral izq

        GridapGmsh.gmsh.model.occ.addCurveLoop([1,2,3,4], 10) #the rectangle
        GridapGmsh.gmsh.model.occ.synchronize()

        # creamos puntos internos
        lc_f=lc/4.0;
        x₁=side_x/3;x₂=2*side_x/3;
        y₁=side_y/3;y₂=2*side_y/3;
        GridapGmsh.gmsh.model.occ.addPoint(x₁,y₁,0,lc_f,5)     # 5 vertice interno inferior izq
        GridapGmsh.gmsh.model.occ.addPoint(x₂,y₁,0,lc_f,6)     # 6 vertice interno inferior der
        GridapGmsh.gmsh.model.occ.addPoint(x₂,y₂,0,lc_f,7)     # 7 vertice interno superior der
        GridapGmsh.gmsh.model.occ.addPoint(x₁,y₂,0,lc_f,8)     # 8 vertice interno superior izq
        # make internal rectangle
        GridapGmsh.gmsh.model.occ.addLine(5,6,5) # 1 linea interna inferior
        GridapGmsh.gmsh.model.occ.addLine(6,7,6) # 2 linea interna lateral der
        GridapGmsh.gmsh.model.occ.addLine(7,8,7) # 4 linea interna superior 
        GridapGmsh.gmsh.model.occ.addLine(8,5,8) # 4 linea interna lateral izq    
        # make internal loop
        GridapGmsh.gmsh.model.occ.addCurveLoop([5,6,7,8],11) # the internal rectangle
        GridapGmsh.gmsh.model.occ.synchronize()

        # make the surface
        GridapGmsh.gmsh.model.occ.addPlaneSurface([10,11], 100) #the surface
        # GridapGmsh.gmsh.model.occ.addPlaneSurface([10], 100) #the surface
        GridapGmsh.gmsh.model.occ.synchronize()
        
        type_structured_mesh="AlternateLeft";

        if (structured_mesh == true)
            #= ++++++++++++++++++++++++++++++++++++++++++++++++++++++
            The `setTransfiniteCurve()' meshing constraints explicitly specifies the
            location of the nodes on the curve.
            ++++++++++++++++++++++++++++++++++++++++++++++++++++++ =# 
            # Creamos curvas de interpolación inferior y superior
            # bumpFactor = 0.20 (default)
            GridapGmsh.gmsh.model.mesh.setTransfiniteCurve(1, numNodesHE_hor, "Bump", bumpFactor)
            GridapGmsh.gmsh.model.mesh.setTransfiniteCurve(3, numNodesHE_hor, "Bump", bumpFactor)
            # Creamos curvas de interpolación lateral izquierda y derecha
            GridapGmsh.gmsh.model.mesh.setTransfiniteCurve(2, numNodesHE_ver, "Bump", bumpFactor)
            GridapGmsh.gmsh.model.mesh.setTransfiniteCurve(4, numNodesHE_ver, "Bump", bumpFactor)
            
            #= ++++++++++++++++++++++++++++++++++++++++++++++++++++++
            The `setTransfiniteSurface()' meshing constraint uses a transfinite
            interpolation algorithm in the parametric plane of the surface to connect
            the nodes on the boundary using a structured grid. If the surface has more
            than 4 corner points, the corners of the transfinite interpolation have to
            be specified by hand:
            The way triangles are generated can be controlled by specifying "Left",
            "Right" or "Alternate" in `setTransfiniteSurface()' command.
            ++++++++++++++++++++++++++++++++++++++++++++++++++++++ =#
            # creamos malla estructurada en la cara 2D
            if (type_structured_mesh == "Default")
                GridapGmsh.gmsh.model.mesh.setTransfiniteSurface(100) # for structured mesh
            elseif (type_structured_mesh == "Left" && quad_state == false)
                GridapGmsh.gmsh.model.mesh.setTransfiniteSurface(100,"Left")
            elseif (type_structured_mesh == "Right" && quad_state == false)
                GridapGmsh.gmsh.model.mesh.setTransfiniteSurface(100,"Right")
            elseif (type_structured_mesh == "Alternate" && quad_state == false)
                GridapGmsh.gmsh.model.mesh.setTransfiniteSurface(100,"Alternate")
            elseif (type_structured_mesh == "AlternateLeft" && quad_state == false)
                GridapGmsh.gmsh.model.mesh.setTransfiniteSurface(100,"AlternateLeft")
            elseif (type_structured_mesh == "AlternateRight" && quad_state == false)
                GridapGmsh.gmsh.model.mesh.setTransfiniteSurface(100,"AlternateRight")
            end
        end

        if (quad_state == true)
            println("Choose FE-quadrilaterals");
            GridapGmsh.gmsh.model.mesh.setRecombine(2, 100) # for 2D quadrilaterals
        else
            println("Choose FE-triangles (default)");
        end

        # creamos grupos para definir condiciones de bordes
        # GridapGmsh.gmsh.model.addPhysicalGroup(dimensión,elementos,tag)

        GridapGmsh.gmsh.model.geo.addPhysicalGroup(0,[1,2,3,4],12,"ext_vertices")  # grupo formado por los 4 vértices externos
        # GridapGmsh.gmsh.model.setPhysicalName(0, 12, "ext_vertices")             # le damos nombre al grupo (dim = 0)
        GridapGmsh.gmsh.model.geo.synchronize()                                    # sincronizamos para que sea visible

        GridapGmsh.gmsh.model.addPhysicalGroup(1, [1,2,3,4], 11 )  # grupo formado por las 4 lineas del borde
        GridapGmsh.gmsh.model.setPhysicalName(1, 11, "ext")        # le damos nombre al grupo (dim = 1)
        GridapGmsh.gmsh.model.occ.synchronize()                    # sincronizamos para que sea visible   

        GridapGmsh.gmsh.model.addPhysicalGroup(1,[5,6,7,8],14) # grupo para el loop interno
        GridapGmsh.gmsh.model.setPhysicalName(1,14,"int_loop")
        GridapGmsh.gmsh.model.occ.synchronize()

        GridapGmsh.gmsh.model.addPhysicalGroup(0,[5,6,7,8],13)  # grupo formado por los dos puntos internos
        GridapGmsh.gmsh.model.setPhysicalName(0,13,"int_points")
        GridapGmsh.gmsh.model.occ.synchronize()

        GridapGmsh.gmsh.model.addPhysicalGroup(2, [100], 101)      # grupo formado por cara 2D
        GridapGmsh.gmsh.model.setPhysicalName(2, 101, "surface")   # sincronizamos para que sea visible
        GridapGmsh.gmsh.model.occ.synchronize()                    # sincronizamos para que sea visible

        # generamos mesh 2D
        GridapGmsh.gmsh.model.mesh.generate(2)
        if (path ≠ " ")
            GridapGmsh.gmsh.write(path*name*".msh")
        end

        # # esto abre una consola interactiva con gmsh
        # if !("-nopopup" in ARGS)
        #     GridapGmsh.gmsh.fltk.run()
        # end

        # finalizamos armado de grilla
        GridapGmsh.gmsh.finalize()

        # guardamos mesh en una variable
        model = GmshDiscreteModel(path*name*".msh")

    elseif (grid_type=="Cartesian2D")
        dom,n = params
        model=CartesianDiscreteModel(dom,n);
    end

    return model;
end