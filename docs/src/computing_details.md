# Comparison between ARPACK, LAPACK and BLAS packages

| **Library** | **Purpose** | **Usage** | **Performance** | **Dependence** |
|-------------|-------------|-----------|-----------------|----------------|
| **ARPACK**  | Solves large-scale eigenvalue problems for sparse matrices. | Suitable for large, sparse matrices. | Optimized for sparse matrices, efficient for large-scale problems. | Depends on LAPACK and BLAS for some operations[1][2]. |
| **LAPACK**  | Provides routines for solving systems of linear equations, least squares problems, and eigenvalue problems for dense and banded matrices. | Suitable for dense matrices or sparse matrices with trivial structure. | Optimized for dense matrices, efficient for high-performance computing. | Often used with BLAS for low-level operations[3][5]. |
| **BLAS**    | Provides low-level routines for performing basic linear algebra operations such as vector and matrix operations. | Underlies LAPACK and other higher-level linear algebra libraries. | Optimized for low-level operations, highly efficient for matrix multiplications and vector operations. | Often used as a building block for other libraries like LAPACK and ATLAS[4][5]. |

### Additional Notes:
- **ARPACK**: Designed for solving large-scale eigenvalue problems efficiently, especially for sparse matrices. It depends on LAPACK and BLAS for some operations[1][2].
- **LAPACK**: A comprehensive library for dense and banded matrices, providing routines for matrix factorizations and solving linear systems. It is often used with BLAS for optimized performance[3][5].
- **BLAS**: The foundation for many linear algebra libraries, providing low-level operations like matrix multiplication and vector operations. It is highly optimized for performance and is used as a building block for other libraries[4][5].

This table highlights the key differences in the purpose, usage, and performance characteristics of these libraries, which are essential for scientific computing and high-performance computing applications.

Citations:

+ [1] [http://li.mit.edu/Archive/Activities/Archive/CourseWork/Ju_Li/MITCourses/18.335/Doc/ARPACK/Lehoucq97.pdf](http://li.mit.edu/Archive/Activities/Archive/CourseWork/Ju_Li/MITCourses/18.335/Doc/ARPACK/Lehoucq97.pdf)
+ [2] [http://rsuib.cc.rsu.ru/libraries/ARPACK/node1.html](http://rsuib.cc.rsu.ru/libraries/ARPACK/node1.html)
+ [3] [https://citeseerx.ist.psu.edu/document?doi=02b94ccb30bbefadc557938aee084c204e0629a1&repid=rep1&type=pdf](https://citeseerx.ist.psu.edu/document?doi=02b94ccb30bbefadc557938aee084c204e0629a1&repid=rep1&type=pdf)
+ [4] [https://stackoverflow.com/questions/33895423/what-is-high-performance-version-of-lapack-and-blas](https://stackoverflow.com/questions/33895423/what-is-high-performance-version-of-lapack-and-blas)
+ [5] [https://stackoverflow.com/questions/17858104/what-is-the-relation-between-blas-lapack-and-atlas](https://stackoverflow.com/questions/17858104/what-is-the-relation-between-blas-lapack-and-atlas)