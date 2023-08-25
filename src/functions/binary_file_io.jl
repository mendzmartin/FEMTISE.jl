# Reference: https://gist.github.com/dataPulverizer/3dc0af456a427aeb704a437e31299242
# run command = julia binary_file_io.jl
"""
    write_bin(x,fileName;<keyword arguments>)

Write binary file from array data

...
# Arguments
- `x::Array{T,1}`: array/vector data
- `fileName::String`: name of file data
- `existing_file::Bool=false`: optional boolean keyword to delete or not delete existing data
...
"""
function write_bin(x::Array{T,1},fileName::String;existing_file::Bool=false)::Int where T
    if existing_file
        rm(fileName)
    end
    # Open the file
    io=open(fileName,"w")
    # Cast this number to make sure we know its type
    write(io,Int(size(x)[1]))
    # Get the type as a string
    typ=repr(T)
    # Write the length of the type string
    write(io,Int(length(typ)))
    # Now write the type string
    for i in eachindex(typ)
        write(io,Char(typ[i]))
    end
    # Now write the array
    for i in eachindex(x)
        write(io,x[i])
    end
    # Clean up
    close(io)
    return 0;
end

"""
    write_bin(x,fileName;<keyword arguments>)

Write binary file from matrix data

...
# Arguments
- `x::Matrix{T}`: matrix data
- `fileName::String`: name of file data
- `existing_file::Bool=false`: optional boolean keyword to delete or not delete existing data
...
"""
function write_bin(x::Matrix{T},fileName::String;existing_file::Bool=false)::Int where T
    if existing_file
        rm(fileName)
    end
    # Open the file
    io=open(fileName,"w")
    # Cast this number to make sure we know its type
    write(io,Int(size(x)[1]))
    # Get the type as a string
    typ=repr(T)
    # Write the length of the type string
    write(io,Int(length(typ)))
    # Now write the type string
    for i in eachindex(typ)
        write(io,Char(typ[i]))
    end
    # Now write the array
    for i in eachindex(x)
        write(io,x[i])
    end
    # Clean up
    close(io)
    return 0;
end

"""
    read_bin(io,::Type{T},n,matrix_data,c_dim)

Speeds up the read binary file

...
# Arguments
- `io::IO`: in/output variable
- `::Type{T}`: data type
- `n::Int`: total number of elements in array/matrix data
- `matrix_data::Bool`: boolean keyword to specify matrix or array data
- `c_dim::Int`: column number of matrix data
...
"""
function read_bin(io::IO,::Type{T},n::Int,matrix_data::Bool,c_dim::Int) where T
    (matrix_data==true) ? (x=Array{T,2}(undef,n,c_dim)) : (x=Array{T,1}(undef, n))
    for i in eachindex(x)
        x[i]=read(io,T)
    end
    close(io)
    return x
end

"""
    read_bin(fileName;<keyword arguments>)

Read binary file

...
# Arguments
- `fileName::String`: name of file data
- `matrix_data::Bool`: optional boolean keyword to specify matrix or array data
- `c_dim::Int`: optional column number of matrix data
...
"""
function read_bin(fileName::String;matrix_data=false,c_dim::Int=1)
    # Open the file
    io=open(fileName,"r")
    # Read the total number of elements in the resulting array
    n=read(io,Int)
    # Read the length of the type name
    nt=read(io,Int)
    # Then read the type name
    cName=Array{Char}(undef,nt)
    for i in eachindex(cName)
        cName[i]=read(io,Char)
    end
    # The return type
    T=eval(Symbol(String(cName)))
    # The data
    x=read_bin(io,T,n,matrix_data,c_dim)
    return x
end