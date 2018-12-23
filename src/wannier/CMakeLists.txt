cmake_minimum_required(VERSION 2.8)
if(${CMAKE_PROJECT_NAME} STREQUAL "Project")
  message(FATAL_ERROR "cmake should be executed not for 'src' subdirectory, but for the top directory of RESPACK.")
endif(${CMAKE_PROJECT_NAME} STREQUAL "Project")

set(SOURCES_COMMON m_rdinput.f90 m_fft3d_20150826.f90 m_bvector_20171208.f90
diagV.F90 det.F90 inv.F90 inv_ge_LAPACK.F90 makekpts.f90 est_NTK.f90 est_nwx2.f90 est_latparam.f90 wrt_mvmc.f90 wrt_frmsf.f90 search_Rmin.f90 frmsf_wan.F90 
gencif.f90 m_gencif_sub.f90 m_dos.f90 m_eigenstate.f90 m_tetrahedron_20170325.F
)

add_executable(calc_wannier wannier.F wannier_sub.F ${SOURCES_COMMON})
target_link_libraries(calc_wannier ${LAPACK_LIBRARIES} m)

install(TARGETS calc_wannier RUNTIME DESTINATION bin)
