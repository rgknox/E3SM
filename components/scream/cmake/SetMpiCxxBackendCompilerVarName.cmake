MACRO (SET_MPI_CXX_BACKEND_COMPILER_VAR_NAME)

  execute_process(COMMAND ${CMAKE_CXX_COMPILER} 
                  ${CMAKE_SOURCE_DIR}/cmake/TryCompileOpenMPI.cpp
                  OUTPUT_QUIET ERROR_QUIET
                  RESULT_VARIABLE HAVE_OPENMPI_MPI)

  IF (HAVE_OPENMPI_MPI)
    SET (MPI_CXX_BACKEND_COMPILER_VAR_NAME  OMPI_CXX)
  ELSE ()
    execute_process(COMMAND ${CMAKE_CXX_COMPILER} 
                    ${CMAKE_SOURCE_DIR}/cmake/TryCompileMPICH.cpp
                    OUTPUT_QUIET ERROR_QUIET
                    RESULT_VARIABLE HAVE_MPICH_MPI)

    IF (HAVE_MPICH_MPI)
      SET (MPI_CXX_BACKEND_COMPILER_VAR_NAME MPICH_CXX)
    ELSE ()
      MESSAGE (FATAL_ERROR "Unsupported MPI distribution.")
    ENDIF()
  ENDIF()
ENDMACRO()
