if(${CMAKE_PROJECT_NAME} STREQUAL "Project")
  message(FATAL_ERROR "cmake should be executed not for 'src' subdirectory, but for the top directory of RESPACK.")
endif(${CMAKE_PROJECT_NAME} STREQUAL "Project")

add_executable(strconv strconv.F90 subr_fmtconv.F90 subr_readstr.F90)
add_executable(wfn2respack wfn2respack.F90)
install(TARGETS strconv wfn2respack RUNTIME DESTINATION bin)

install(FILES xtapp2respack.sh DESTINATION bin
        PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE
        WORLD_READ WORLD_EXECUTE)
