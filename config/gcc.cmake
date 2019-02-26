# for GCC Compiler
set(CMAKE_C_COMPILER "gcc" CACHE STRING "" FORCE)
set(CMAKE_Fortran_COMPILER "gfortran" CACHE STRING "" FORCE)
set(CMAKE_Fortran_FLAGS_RELEASE " -O2 -fopenmp -g -fbacktrace -ffree-line-length-none")
set(CMAKE_Fortran_FLAGS_DEBUG " -O0 -fopenmp -g -fbacktrace -ffree-line-length-none")
