# for Intel Compiler
# to detect omp flag
set(CMAKE_C_COMPILER "icc" CACHE STRING "" FORCE)
set(CMAKE_Fortran_COMPILER "ifort" CACHE STRING "" FORCE)
set(CMAKE_Fortran_FLAGS_RELEASE " -O2 -g -fbacktrace -DNDEBUG -xCORE-AVX2 -mcmodel=large -shared-intel")

# for Intel MKL
set(BLA_VENDOR "Intel10_64lp" CACHE STRING "" FORCE)
