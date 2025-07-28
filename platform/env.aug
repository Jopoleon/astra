#!/bin/bash

PKG=/shares/departments/AUG/users/git/soft/intel/oneapi
COMP_BIN=${PKG}/compiler/2025.2/bin
export I_MPI_ROOT=${PKG}/mpi/2021.16
export I_MPI_FABRICS=shm:ofi
export I_MPI_DEBUG=5
export I_MPI_LINK=opt
export MKL_LIBDIR=${PKG}/mkl/2025.2/lib
export COMP_LIBDIR=${PKG}/compiler/2025.2/lib
LIBFABRIC_DIR=$I_MPI_ROOT/opt/mpi/libfabric/lib
export LD_LIBRARY_PATH=${COMP_LIBDIR}

export FI_PROVIDER_PATH=$I_MPI_ROOT/opt/mpi/libfabric/lib/prov # necessary

MPI_BIN=$I_MPI_ROOT/bin
PYTHON_BIN=/shares/software/aug-dv/moduledata/miniforge3/2024.09/bin
if [[ ":$PATH:" != *":$PYTHON_BIN:"* ]]
then
    export PATH=${PYTHON_BIN}:${PATH}
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
export FC=ifx
export CC=/usr/bin/x86_64-linux-gnu-gcc-11
#export CC=icx
export MPIFC=mpiifx
