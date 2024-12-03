#!/bin/bash

PKG=/shares/software/aug-dv/moduledata/Intel/oneapi
COMP_BIN=${PKG}/compiler/2023.1.0/linux/bin/intel64
export I_MPI_ROOT=${PKG}/mpi/2021.9.0
export MKL_LIBDIR=${PKG}/mkl/2023.1.0/lib/intel64
export COMP_LIBDIR=${PKG}/compiler/2023.1.0/linux/compiler/lib/intel64_lin
export LD_LIBRARY_PATH=${COMP_LIBDIR}
MPI_BIN=$I_MPI_ROOT/bin
PYTHON_BIN=/shares/software/aug-dv/moduledata/miniforge3/2024.09/bin
if [[ ":$PATH:" != *":$PYTHON_BIN:"* ]]
then
    export PATH=${COMP_BIN}:${I_MPI_ROOT}/bin:${PYTHON_BIN}:${PATH}
fi
if [[ ":$PATH:" != *":$COMP_BIN:"* ]]
then
    export PATH=${COMP_BIN}:${PATH}
fi
if [[ ":$PATH:" != *":$MPI_BIN:"* ]]
then
    export PATH=${MPI_BIN}:${PATH}
fi
export ASTRA_EXT=/shares/departments/AUG/users/git/ASTRA_LIBRARIES_EXT
export FC=ifort
export CC=/usr/bin/x86_64-linux-gnu-gcc-11
export MPIFC=mpiifort
