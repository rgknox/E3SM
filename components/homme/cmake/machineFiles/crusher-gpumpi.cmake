#interactive job
#bsub -W 2:00 -nnodes 1 -P cli115 -Is /bin/bash


#cmake -C ~/acme-fork-lb/components/homme/cmake/machineFiles/summit.cmake -DHOMMEXX_MPI_ON_DEVICE=FALSE ~/acme-fork-lb/components/homme/

#SET (HOMMEXX_MPI_ON_DEVICE FALSE CACHE BOOL "")
SET (HOMMEXX_CUDA_MAX_WARP_PER_TEAM "16" CACHE STRING  "")

SET (NETCDF_DIR $ENV{OLCF_NETCDF_FORTRAN_ROOT} CACHE FILEPATH "")
SET (HDF5_DIR $ENV{OLCF_HDF5_ROOT} CACHE FILEPATH "")

SET(BUILD_HOMME_WITHOUT_PIOLIBRARY TRUE CACHE BOOL "")

SET(HOMME_FIND_BLASLAPACK TRUE CACHE BOOL "")

SET(WITH_PNETCDF FALSE CACHE FILEPATH "")

SET(USE_QUEUING FALSE CACHE BOOL "")

SET(ENABLE_CUDA FALSE CACHE BOOL "")

SET(BUILD_HOMME_PREQX_KOKKOS TRUE CACHE BOOL "")
SET(BUILD_HOMME_THETA_KOKKOS TRUE CACHE BOOL "")
SET(HOMME_ENABLE_COMPOSE FALSE CACHE BOOL "")

#SET (HOMMEXX_BFB_TESTING TRUE CACHE BOOL "")

SET(USE_TRILINOS OFF CACHE BOOL "")

SET(HIP TRUE CACHE BOOL "")
#SET(Kokkos_ENABLE_OPENMP OFF CACHE BOOL "")
#SET(Kokkos_ENABLE_CUDA OFF CACHE BOOL "")
#SET(Kokkos_ENABLE_CUDA_LAMBDA OFF CACHE BOOL "")
#SET(Kokkos_ARCH_VEGA908 ON CACHE BOOL "")
#SET(Kokkos_ENABLE_EXPLICIT_INSTANTIATION OFF CACHE BOOL "")

SET(CMAKE_C_COMPILER "cc" CACHE STRING "")
SET(CMAKE_Fortran_COMPILER "ftn" CACHE STRING "")
SET(CMAKE_CXX_COMPILER "hipcc" CACHE STRING "")
#SET(CMAKE_CXX_COMPILER "/ccs/home/onguba/kokkos/bin/nvcc_wrapper" CACHE STRING "")
SET(E3SM_KOKKOS_PATH "/ccs/home/onguba/kokkos-crusher-june2022/bld-hipcc" CACHE STRING "")

#not the proper way!!!
SET(MPICH_DIR "/opt/cray/pe/mpich/8.1.12/ofi/crayclang/10.0" CACHE STRING "")

SET(Extrae_LIBRARY "-I${MPICH_DIR}/include -L${MPICH_DIR}/lib -lmpi -L/opt/cray/pe/mpich/8.1.12/gtl/lib -lmpi_gtl_hsa" CACHE STRING "")

SET(ADD_Fortran_FLAGS "-O3 -DNDEBUG ${Extrae_LIBRARY} -I${E3SM_KOKKOS_PATH}/include -L${E3SM_KOKKOS_PATH}/lib64" CACHE STRING "")
SET(ADD_C_FLAGS "-O3 -DNDEBUG ${Extrae_LIBRARY} -I${E3SM_KOKKOS_PATH}/include -L${E3SM_KOKKOS_PATH}/lib64" CACHE STRING "")
SET(ADD_CXX_FLAGS "-std=c++14 -O3 -DNDEBUG --amdgpu-target=gfx90a -fno-gpu-rdc ${Extrae_LIBRARY} -I${E3SM_KOKKOS_PATH}/include -L${E3SM_KOKKOS_PATH}/lib64" CACHE STRING "")
SET(ADD_LINKER_FLAGS "-O3 -DNDEBUG ${Extrae_LIBRARY} -I${E3SM_KOKKOS_PATH}/include -L${E3SM_KOKKOS_PATH}/lib64" CACHE STRING "")


set (ENABLE_OPENMP OFF CACHE BOOL "")
set (ENABLE_COLUMN_OPENMP OFF CACHE BOOL "")
set (ENABLE_HORIZ_OPENMP OFF CACHE BOOL "")

set (HOMME_TESTING_PROFILE "dev" CACHE STRING "")

set (USE_NUM_PROCS 4 CACHE STRING "")

#set (OPT_FLAGS "-mcpu=power9 -mtune=power9" CACHE STRING "")
SET (USE_MPI_OPTIONS "--bind-to core" CACHE FILEPATH "")
