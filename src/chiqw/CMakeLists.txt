cmake_minimum_required(VERSION 2.8)
if(${CMAKE_PROJECT_NAME} STREQUAL "Project")
  message(FATAL_ERROR "cmake should be executed not for 'src' subdirectory, but for the top directory of RESPACK.")
endif(${CMAKE_PROJECT_NAME} STREQUAL "Project")

add_definitions(${MPI_Fortran_COMPILE_FLAGS})
link_directories(${MPI_Fortran_LIBRARY_DIRS})
include_directories(${MPI_Fortran_INCLUDE_PATH})

set(SOURCES_CHIQW est_NTK.f90 est_nwx2.f90 gen_grid.f90 inv.F90 invZGE.F90
m_fft3d_20150826.f90 m_rdinput.f90 m_tetrahedron_20170325.F chiqw_sub.F)
add_library(chiqw_com ${SOURCES_CHIQW} )
target_link_libraries(chiqw_com ${LAPACK_LIBRARIES} m)

add_executable(chiqw chiqw.F)

target_link_libraries(chiqw chiqw_com)
target_link_libraries(chiqw ${MPI_Fortran_LIBRARIES})

install(TARGETS chiqw ARCHIVE DESTINATION lib LIBRARY DESTINATION lib RUNTIME DESTINATION bin)
