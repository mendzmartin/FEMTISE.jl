export make_boundary_conditions
function make_boundary_conditions(grid_type::String,BC_type::String,TypeData::Type;homogeneous::Bool=true,params::Tuple=(nothing,nothing))
    if (grid_type == "simple_line")
        if (BC_type == "FullDirichlet")
            if homogeneous
                BC_values = [zero(TypeData),zero(TypeData)]
            else
                left_value,right_value=params
                BC_values = [left_value,right_value]
            end
            BC_tags = ["left_point","right_point"];
        end
    elseif grid_type == "simple_rectangle_v1"
        if (BC_type=="FullDirichlet")
            BC_values = [zero(TypeData),zero(TypeData)];
            BC_tags = ["ext_points","ext_lines"]
        end
    elseif grid_type == "simple_rectangle_v2"
        if (BC_type=="FullDirichlet")
            BC_values = [zero(TypeData),zero(TypeData)];
            BC_tags == ["ext_vertices","ext"]
        end
    elseif (grid_type == "Cartesian2D")
        if (BC_type == "FullDirichlet")
            BC_values = (zero(TypeData));
            BC_tags = "boundary";
        end
    end

    return BC_values,BC_tags
end

# we need to change TypeData argument as optional argument