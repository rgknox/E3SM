if (MPILIB STREQUAL mpi-serial AND NOT compile_threaded)
  set(PFUNIT_PATH "/projects/ccsm/pfunit/3.2.9/mpi-serial")
endif()
set(PIO_FILESYSTEM_HINTS "lustre")
