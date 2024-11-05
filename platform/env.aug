#!/bin/bash -f

PKG=/shares/software/aug-dv/moduledata/Intel/oneapi
IFORT_BIN=${PKG}/compiler/2023.1.0/linux/bin/intel64
export I_MPI_ROOT=${PKG}/mpi/2021.9.0
export MKL_LIBDIR=${PKG}/mkl/2023.1.0/lib/intel64
export COMP_LIBDIR=$INTEL_HOME/compiler/2023.1.0/linux/compiler/lib/intel64_lin
export LD_LIBRARY_PATH=${COMP_LIBDIR}
export PATH=${PATH}:${IFORT_BIN}

export ASTRA_EXT=/shares/departments/AUG/users/git/ASTRA_LIBRARIES_EXT
export AFC=ifort
export ACC="cc -O"
export MPIFC=${I_MPI_ROOT}/bin/mpiifort
