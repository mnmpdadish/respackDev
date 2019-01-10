cmake_minimum_required(VERSION 2.8.10)
if(${CMAKE_PROJECT_NAME} STREQUAL "Project")
  message(FATAL_ERROR "cmake should be executed not for 'src' subdirectory, but for the top directory of RESPACK.")
endif(${CMAKE_PROJECT_NAME} STREQUAL "Project")

find_package(MPI REQUIRED)
add_definitions(${MPI_Fortran_COMPILE_FLAGS})
include_directories(${MPI_Fortran_INCLUDE_PATH})
link_directories(${MPI_Fortran_LIBRARY_DIRS})
set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})

set(SOURCES_CHIQW est_NTK.f90 est_nwx2.f90 gen_grid.f90 inv.F90 invZGE.F90
m_fft3d_20150826.f90 m_rdinput.f90 m_tetrahedron_20170325.F chiqw_sub.F
wrt_optical_property.f90)
add_library(chiqw_com STATIC ${SOURCES_CHIQW} )
target_link_libraries(chiqw_com ${LAPACK_LIBRARIES} m)

add_custom_target(chiqw ALL)
add_executable(calc_chiqw chiqw.F)
add_dependencies(chiqw calc_chiqw)
target_link_libraries(calc_chiqw chiqw_com)
target_link_libraries(calc_chiqw ${MPI_Fortran_LIBRARIES})

install(TARGETS calc_chiqw RUNTIME DESTINATION bin)

