cmake_minimum_required(VERSION 2.8.10)
if(${CMAKE_PROJECT_NAME} STREQUAL "Project")
  message(FATAL_ERROR "cmake should be executed not for 'src' subdirectory, but for the top directory of RESPACK.")
endif(${CMAKE_PROJECT_NAME} STREQUAL "Project")

find_package(MPI REQUIRED)
add_definitions(${MPI_Fortran_COMPILE_FLAGS})
include_directories(${MPI_Fortran_INCLUDE_PATH})
link_directories(${MPI_Fortran_LIBRARY_DIRS})
set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})

set(SOURCES_GW 
calc_band_disp.f90
calc_gwakw.f90
calc_gwdos.f90
calc_ksdos.f90
gw.f90
m_fft3d_20150826.f90
m_rd_dat_eps.f90
m_rd_dat_wan.f90
m_rd_dat_wfn.f90
m_rd_input.f90
make_kpts.f90
make_sgmw.f90
search_Rmin.f90
sub_det_shift.f90
sub_eigenvalue.f90
sub_eps.f90
sub_gw.f90
sub_gwdos.f90
sub_ksdos.f90
sub_mkidx.f90
sub_wfn.f90
util.f90
m_tetrahedron_20170325.F
)

add_custom_target(gw ALL)
add_executable(calc_gw ${SOURCES_GW})
target_link_libraries(calc_gw ${MPI_Fortran_LIBRARIES})
target_link_libraries(calc_gw ${LAPACK_LIBRARIES} m)
add_dependencies(gw calc_gw)

install(TARGETS calc_gw RUNTIME DESTINATION bin)
