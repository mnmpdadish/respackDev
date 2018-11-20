cmake_minimum_required(VERSION 2.8)
if(${CMAKE_PROJECT_NAME} STREQUAL "Project")
  message(FATAL_ERROR "cmake should be executed not for 'src' subdirectory, but for the top directory of RESPACK.")
endif(${CMAKE_PROJECT_NAME} STREQUAL "Project")

set(SOURCES
calc_dos.f90
calc_hist.f90
est_ef.f90
est_nkbi.f90
m_rd_dat_zvo.f90
m_rd_transdef.f90
make_eig.f90
make_grd.f90
make_kdata.f90
sub_eigenvalue.f90
sub_mkidx.f90
sub_search_Rmin.f90
transfer_analysis.f90
truncation.f90
util.f90
wrt.f90
m_tetrahedron_20170325.F
)

add_executable(transfer_analysis ${SOURCES})
target_link_libraries(transfer_analysis ${LAPACK_LIBRARIES})
install(TARGETS transfer_analysis RUNTIME DESTINATION bin)
