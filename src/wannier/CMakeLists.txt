cmake_minimum_required(VERSION 2.8)
if(${CMAKE_PROJECT_NAME} STREQUAL "Project")
  message(FATAL_ERROR "cmake should be executed not for 'src' subdirectory, but for the top directory of RESPACK.")
endif(${CMAKE_PROJECT_NAME} STREQUAL "Project")

set(SOURCES_COMMON det.F90 diagV.F90 est_NTK.f90 est_latparam.f90 est_nwx2.f90 frmsf_wan.F90 gencif.f90 inv.F90 
inv_ge_LAPACK.F90 m_bvector_20171208.f90 m_dos.f90 m_eigenstate.f90 m_fft3d_20150826.f90 m_gencif_sub.f90 
m_rdinput.f90 m_tetrahedron_20170325.F makekpts.f90 search_Rmin.f90 wrt_frmsf.f90 wrt_mvmc.f90) 
add_custom_target(wannier ALL)
add_executable(calc_wannier wannier.F wannier_sub.F ${SOURCES_COMMON})
target_link_libraries(calc_wannier ${LAPACK_LIBRARIES} m)
add_dependencies(wannier calc_wannier)

install(TARGETS calc_wannier RUNTIME DESTINATION bin)
