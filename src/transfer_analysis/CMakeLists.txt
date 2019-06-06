cmake_minimum_required(VERSION 2.8)
if(${CMAKE_PROJECT_NAME} STREQUAL "Project")
  message(FATAL_ERROR "cmake should be executed not for 'src' subdirectory, but for the top directory of RESPACK.")
endif(${CMAKE_PROJECT_NAME} STREQUAL "Project")

set(SOURCES
m_band.f90
m_dos.f90
m_eigenstate.f90
m_frmsf.f90
m_hist.f90
m_rd_dat_zvo.f90
m_rd_transdef.f90
m_tetrahedron_20170325.F
m_truncation.f90
transfer_analysis.f90
)

add_custom_target(transfer_analysis ALL)
add_executable(calc_tr ${SOURCES})
target_link_libraries(calc_tr ${LAPACK_LIBRARIES})
add_dependencies(transfer_analysis calc_tr)
install(TARGETS calc_tr RUNTIME DESTINATION bin)
