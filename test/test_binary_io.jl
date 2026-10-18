# Binary input/output used by the "bin" output format of the driver.
#
# `write_bin` stores the length, the element type and the raw data; `read_bin`
# must return exactly the same array (vectors) or matrix (given its number of
# columns).

@testset "Binary I/O round trip" begin
    mktempdir() do folder
        vector = ComplexF64[1.0 + 2.0im, -0.5, 3.25im]
        vector_file = joinpath(folder, "vector.bin")
        write_bin(vector, vector_file)
        @test read_bin(vector_file; matrix_data = false) == vector

        matrix = reshape(collect(1.0:12.0), 4, 3)
        matrix_file = joinpath(folder, "matrix.bin")
        write_bin(matrix, matrix_file)
        @test read_bin(matrix_file; matrix_data = true, c_dim = 3) == matrix

        # `existing_file = true` replaces a previous file
        write_bin(2 .* vector, vector_file; existing_file = true)
        @test read_bin(vector_file; matrix_data = false) == 2 .* vector
    end
end
