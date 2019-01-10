# for Intel Compiler
# to detect omp flag
set(CMAKE_C_COMPILER "icc" CACHE STRING "" FORCE)
set(CMAKE_MPI_Fortran_COMPILER "mpiifort" CACHE STRING "" FORCE)
set(CMAKE_MPI_Fortran_FLAGS_RELEASE " -O2 -qopenmp -xHost -mcmodel=large -shared-intel")  
set(CMAKE_Fortran_COMPILER "ifort" CACHE STRING "" FORCE)
set(CMAKE_Fortran_FLAGS_RELEASE " -O2 -g -fbacktrace -DHAVE_SSE2 -DNDEBUG -mcmodel=large -shared-intel")

# for Intel MKL
set(BLA_VENDOR "Intel10_64lp" CACHE STRING "" FORCE)
