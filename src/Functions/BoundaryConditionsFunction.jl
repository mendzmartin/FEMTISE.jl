export make_boundary_conditions
function make_boundary_conditions(grid_type::String,BC_type::String,TypeData::Type)
    if (grid_type == "simple_line")
        if (BC_type == "FullDirichlet")
            BC_values = [zero(TypeData),zero(TypeData)];
            BC_tags = ["left_point","right_point"];
        end
    else if grid_type == "simple_rectangle_v1"
        if (BC_type=="FullDirichlet")
            BC_values = [zero(TypeData),zero(TypeData)];
            BC_tags = ["ext_points","ext_lines"]
        end
    else if grid_type == "simple_rectangle_v2"
        if (BC_type=="FullDirichlet")
            BC_values = [zero(TypeData),zero(TypeData)];
            BC_tags == ["ext_vertices","ext"]
        end
    else if (grid_type == "Cartesian2D")
        if (BC_type == "FullDirichlet")
            BC_values = (zero(TypeData));
            BC_tags = "boundary";
        end
    end

    return BC_values,BC_tags
end
