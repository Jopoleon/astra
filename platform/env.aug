#!/bin/bash

PKG=/shares/software/aug-dv/moduledata/Intel/oneapi
export I_MPI_ROOT=${PKG}/mpi/2021.9.0
export I_MPI_FABRICS=shm
export MKL_LIBDIR=${PKG}/mkl/2023.1.0/lib/intel64
export COMP_LIBDIR=${PKG}/compiler/2023.1.0/linux/compiler/lib/intel64_lin
export LD_LIBRARY_PATH=${COMP_LIBDIR}
export FC=ifort
export CC=icc
export MPIFC=mpiifort
export MPICC=mpiicc
export FC_FLAGS="-r8 -no-prec-div -traceback -O0 -w -fPIC -qopenmp"
export LD_FLAGS="-qopenmp -traceback"

COMP_BIN=${PKG}/compiler/2023.1.0/linux/bin/intel64
MPI_BIN=$I_MPI_ROOT/bin
if [[ ":$PATH:" != *":$COMP_BIN:"* ]]
then
    export PATH=${COMP_BIN}:${MPI_BIN}:${PATH}
fi

export PYTHON_BIN=/shares/software/aug-dv/moduledata/miniforge3/2024.09/bin

export ASTRA_EXT=/shares/departments/AUG/users/git/ASTRA_LIBRARIES_EXT
