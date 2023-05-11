#=
    utils links: https://gitlab.onelab.info/gmsh/gmsh/blob/master/api/gmsh.jl
=#

function make_model(grid_type::String,params::Tuple)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal",1)

    # 1D line grid
    if (grid_type == "simple_line")
        path,name,dom,MeshSize=params
        gmsh.model.add(name)

        gmsh.model.geo.addPoint(dom[1],0,0,MeshSize,1)    # 1 punto vértice izquierdo
        gmsh.model.geo.addPoint(dom[2],0,0,MeshSize,2)    # 2 punto vértice derecho

        gmsh.model.geo.addLine(1,2,1)   # linea que une puntos 1 y 2
        gmsh.model.geo.synchronize()
        
        gmsh.model.geo.addPhysicalGroup(0,[1],3)        # grupo formado por punto izquierdo
        gmsh.model.setPhysicalName(0,3,"left_point")    # le damos nombre al grupo
        gmsh.model.geo.synchronize()

        gmsh.model.geo.addPhysicalGroup(0,[2],4)        # grupo formado por punto derecho
        gmsh.model.setPhysicalName(0,4,"right_point")   # le damos nombre al grupo
        gmsh.model.geo.synchronize()
        
        gmsh.model.geo.addPhysicalGroup(1,[1],1)    # grupo formado por linea
        gmsh.model.setPhysicalName(1,1,"segment")   # le damos nombre a la linea
        gmsh.model.geo.synchronize()

        gmsh.model.mesh.generate(1)
        if (path ≠ " ")
            gmsh.write(path*name*".msh")
        end
        gmsh.finalize()
        model=GmshDiscreteModel(path*name*".msh")
    # 2D simple rectangle grid
    elseif (grid_type == "simple_rectangle_v1")

        path,name,dom,MeshSize,quad_state=params

        # gmsh.option.setNumber(name, value)
        gmsh.option.setNumber("General.Terminal",0)
        gmsh.model.add(name)

        # creamos puntos vértice - gmsh.model.geo.addPoint(x,y,z,meshSize=0.,tag=-1)
        gmsh.model.geo.addPoint(dom[1],dom[3],0,MeshSize,1) # 1 vertice inferior izquierdo
        gmsh.model.geo.addPoint(dom[2],dom[3],0,MeshSize,2) # 2 vértice inferior derecho
        gmsh.model.geo.addPoint(dom[2],dom[4],0,MeshSize,3) # 3 vértice superior derecho
        gmsh.model.geo.addPoint(dom[1],dom[4],0,MeshSize,4) # 4 vértice superior izquierdo
        # creamos lineas de unión entre vértices - gmsh.model.geo.addLine(startTag,endTag,tag=-1)
        gmsh.model.geo.addLine(1,2,1) # 1 linea inferior
        gmsh.model.geo.addLine(2,3,2) # 2 línea lateral derecha
        gmsh.model.geo.addLine(3,4,3) # 3 linea superior
        gmsh.model.geo.addLine(4,1,4) # 4 linea lateral izquierda
        # creamos curva de unión entre lineas
        gmsh.model.geo.addCurveLoop([1,2,3,4],100) # the rectangle
        gmsh.model.geo.synchronize()

        # make the surface - gmsh.model.geo.addPlaneSurface(wireTags,tag=-1)
        gmsh.model.geo.addPlaneSurface([100],101) # the surface
        gmsh.model.geo.synchronize()

        #=
            creamos grupos para definir condiciones de bordes
            gmsh.model.geo.addPhysicalGroup(dim,tags,tag=-1,name="")
            gmsh.model.setPhysicalName(dim,tag,name)
        =#
        # creamos grupo físico de puntos vértices
        gmsh.model.geo.addPhysicalGroup(0,[1,2,3,4],300)  # grupo 0D formado por cuatro puntos
        gmsh.model.setPhysicalName(0,300,"ext_points")
        gmsh.model.geo.synchronize()

        # creamos grupo físico de curva de lineas externa
        # gmsh.model.geo.addPhysicalGroup(1,[100],301)      # de esta forma no funciona!!
        gmsh.model.geo.addPhysicalGroup(1,[1,2,3,4],301)    # grupo 1D formado por cuatro lineas
        gmsh.model.setPhysicalName(1,301,"ext_lines")
        gmsh.model.geo.synchronize()

        # creamos grupo físico con superficie interna del rectángulo
        gmsh.model.addPhysicalGroup(2,[101],302)
        gmsh.model.setPhysicalName(2,302,"surface") # grupo 2D formado por una superficie
        gmsh.model.geo.synchronize()

        if quad_state
            # gmsh.model.mesh.setRecombine(dim,tag,angle=45.)
            gmsh.model.mesh.setRecombine(2,302) # for 2D quadrilaterals
        end

        # gmsh.model.mesh.generate(dim=3)
        gmsh.model.mesh.generate(2)
        if (path ≠ " ")
            gmsh.write(path*name*".msh")
        end
        gmsh.finalize()
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
            gmsh.option.setNumber("Mesh.Algorithm", 5) # delquad
            gmsh.option.setNumber("Mesh.RecombineAll", 1)
        end
        gmsh.model.add(name)

        # first we build the rectangular boundary
        lc_x = lc
        lc_y = lc*side_y/side_x
        # gmsh.model.geo.addPoint(x,y,z,meshSize=0.,tag=-1)
        gmsh.model.occ.addPoint(0,0,0,lc_x,1)           # 1 vertice inferior izq
        gmsh.model.occ.addPoint(side_x,0,0,lc_x,2)      # 2 vértice inferior der
        gmsh.model.occ.addPoint(side_x,side_y,0,lc_y,3) # 3 vértice superior der
        gmsh.model.occ.addPoint(0,side_y,0,lc_y,4)      # 4 vértice superior izq

        # make the square boundary
        gmsh.model.occ.addLine(1, 2, 1) # 1 linea inferior
        gmsh.model.occ.addLine(2, 3, 2) # 2 linea lateral der 
        gmsh.model.occ.addLine(3, 4, 3) # 3 linea superior
        gmsh.model.occ.addLine(4, 1, 4) # 4 linea lateral izq

        gmsh.model.occ.addCurveLoop([1,2,3,4], 10) #the rectangle
        gmsh.model.occ.synchronize()

        # creamos puntos internos
        lc_f=lc/4.0;
        x₁=side_x/3;x₂=2*side_x/3;
        y₁=side_y/3;y₂=2*side_y/3;
        gmsh.model.occ.addPoint(x₁,y₁,0,lc_f,5)     # 5 vertice interno inferior izq
        gmsh.model.occ.addPoint(x₂,y₁,0,lc_f,6)     # 6 vertice interno inferior der
        gmsh.model.occ.addPoint(x₂,y₂,0,lc_f,7)     # 7 vertice interno superior der
        gmsh.model.occ.addPoint(x₁,y₂,0,lc_f,8)     # 8 vertice interno superior izq
        # make internal rectangle
        gmsh.model.occ.addLine(5,6,5) # 1 linea interna inferior
        gmsh.model.occ.addLine(6,7,6) # 2 linea interna lateral der
        gmsh.model.occ.addLine(7,8,7) # 4 linea interna superior 
        gmsh.model.occ.addLine(8,5,8) # 4 linea interna lateral izq    
        # make internal loop
        gmsh.model.occ.addCurveLoop([5,6,7,8],11) # the internal rectangle
        gmsh.model.occ.synchronize()

        # make the surface
        gmsh.model.occ.addPlaneSurface([10,11], 100) #the surface
        # gmsh.model.occ.addPlaneSurface([10], 100) #the surface
        gmsh.model.occ.synchronize()
        
        type_structured_mesh="AlternateLeft";

        if (structured_mesh == true)
            #= ++++++++++++++++++++++++++++++++++++++++++++++++++++++
            The `setTransfiniteCurve()' meshing constraints explicitly specifies the
            location of the nodes on the curve.
            ++++++++++++++++++++++++++++++++++++++++++++++++++++++ =# 
            # Creamos curvas de interpolación inferior y superior
            # bumpFactor = 0.20 (default)
            gmsh.model.mesh.setTransfiniteCurve(1, numNodesHE_hor, "Bump", bumpFactor)
            gmsh.model.mesh.setTransfiniteCurve(3, numNodesHE_hor, "Bump", bumpFactor)
            # Creamos curvas de interpolación lateral izquierda y derecha
            gmsh.model.mesh.setTransfiniteCurve(2, numNodesHE_ver, "Bump", bumpFactor)
            gmsh.model.mesh.setTransfiniteCurve(4, numNodesHE_ver, "Bump", bumpFactor)
            
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
                gmsh.model.mesh.setTransfiniteSurface(100) # for structured mesh
            elseif (type_structured_mesh == "Left" && quad_state == false)
                gmsh.model.mesh.setTransfiniteSurface(100,"Left")
            elseif (type_structured_mesh == "Right" && quad_state == false)
                gmsh.model.mesh.setTransfiniteSurface(100,"Right")
            elseif (type_structured_mesh == "Alternate" && quad_state == false)
                gmsh.model.mesh.setTransfiniteSurface(100,"Alternate")
            elseif (type_structured_mesh == "AlternateLeft" && quad_state == false)
                gmsh.model.mesh.setTransfiniteSurface(100,"AlternateLeft")
            elseif (type_structured_mesh == "AlternateRight" && quad_state == false)
                gmsh.model.mesh.setTransfiniteSurface(100,"AlternateRight")
            end
        end

        if (quad_state == true)
            println("Choose FE-quadrilaterals");
            gmsh.model.mesh.setRecombine(2, 100) # for 2D quadrilaterals
        else
            println("Choose FE-triangles (default)");
        end

        # creamos grupos para definir condiciones de bordes
        # gmsh.model.addPhysicalGroup(dimensión,elementos,tag)

        gmsh.model.geo.addPhysicalGroup(0,[1,2,3,4],12,"ext_vertices")  # grupo formado por los 4 vértices externos
        # gmsh.model.setPhysicalName(0, 12, "ext_vertices")             # le damos nombre al grupo (dim = 0)
        gmsh.model.geo.synchronize()                                    # sincronizamos para que sea visible

        gmsh.model.addPhysicalGroup(1, [1,2,3,4], 11 )  # grupo formado por las 4 lineas del borde
        gmsh.model.setPhysicalName(1, 11, "ext")        # le damos nombre al grupo (dim = 1)
        gmsh.model.occ.synchronize()                    # sincronizamos para que sea visible   

        gmsh.model.addPhysicalGroup(1,[5,6,7,8],14) # grupo para el loop interno
        gmsh.model.setPhysicalName(1,14,"int_loop")
        gmsh.model.occ.synchronize()

        gmsh.model.addPhysicalGroup(0,[5,6,7,8],13)  # grupo formado por los dos puntos internos
        gmsh.model.setPhysicalName(0,13,"int_points")
        gmsh.model.occ.synchronize()

        gmsh.model.addPhysicalGroup(2, [100], 101)      # grupo formado por cara 2D
        gmsh.model.setPhysicalName(2, 101, "surface")   # sincronizamos para que sea visible
        gmsh.model.occ.synchronize()                    # sincronizamos para que sea visible

        # generamos mesh 2D
        gmsh.model.mesh.generate(2)
        if (path ≠ " ")
            gmsh.write(path*name*".msh")
        end

        # # esto abre una consola interactiva con gmsh
        # if !("-nopopup" in ARGS)
        #     gmsh.fltk.run()
        # end

        # finalizamos armado de grilla
        gmsh.finalize()

        # guardamos mesh en una variable
        model = GmshDiscreteModel(path*name*".msh")

    elseif (grid_type=="Cartesian2D")
        dom,n = params
        model=CartesianDiscreteModel(dom,n);
    end

    return model;
end