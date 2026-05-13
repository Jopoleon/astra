#!/bin/bash

rootdir=`dirname $0`      # may be relative path
export AWD=`cd $rootdir && pwd`  # ensure absolute path
if [ $# -ge 1 ]
then
    comp=$1
fi

chmod 744 $AWD/get_platform
platform=`$AWD/get_platform`
if [[ -v comp ]]
then
    platform=${platform}_${comp}
else
    comp=""
fi

echo $platform
source $AWD/platform/env.${platform}

FC_SERIAL=$FC

#------------------
# Packages versions
#------------------

INTEL_VERSION=2025.02
CMAKE_VERSION=cmake-3.30.3-linux-x86_64
JSON_VERSION=9.0.2
NETCDF_VERSION=mar26
RABBIT_VERSION=unstable
TORBEAM_VERSION=unstable
SPIDER_VERSION=unstable
QLK_VERSION=unstable
QLKNN_VERSION=unstable
GA_VERSION=unstable
STRAHL_VERSION=unstable

#--------------------
# ASTRA install paths
#--------------------

SOFT_ROOT=$HOME/soft

if [ "$FC" = 'gfortran' ]
then
    JSON_INSTALL=$ASTRA_EXT/json/gcc_$JSON_VERSION
    NETCDF_INSTALL=$ASTRA_EXT/netcdf/gcc_$NETCDF_VERSION
    export FFLAGS="-fopenmp"
else
    JSON_INSTALL=$ASTRA_EXT/json/$JSON_VERSION
    NETCDF_INSTALL=$ASTRA_EXT/netcdf/$NETCDF_VERSION
    export FFLAGS="-qopenmp"
fi
RABBIT_INSTALL=$ASTRA_EXT/rabbit/$RABBIT_VERSION
TORBEAM_INSTALL=$ASTRA_EXT/torbeam/$TORBEAM_VERSION
SPIDER_INSTALL=$ASTRA_EXT/spider/$SPIDER_VERSION
QLK_INSTALL=$ASTRA_EXT/qualikiz/$QLK_VERSION
QLKNN_INSTALL=$ASTRA_EXT/qlk_nn/$QLKNN_VERSION
TGLF_INSTALL=$ASTRA_EXT/tglf/$GA_VERSION
NEO_INSTALL=$ASTRA_EXT/neo/$GA_VERSION
STRAHL_INSTALL=$ASTRA_EXT/strahl/$STRAHL_VERSION

# To restore at the script end
PATH_OLD=$PATH

export PATH=$SOFT_ROOT/$CMAKE_VERSION/bin:$PATH

CMAKE=$SOFT_ROOT/$CMAKE_VERSION/bin/cmake

mkdir -p $SOFT_ROOT

#------------------
# Intel + MKL + MPI
#------------------

read -p "Install INTEL (y/n) " INTEL_FLAG

if [ "$INTEL_FLAG" = "y" ]
then
    cd $SOFT_ROOT
    INTEL_HOME=$SOFT_ROOT/intel/$INTEL_VERSION/oneapi
    mkdir -p $INTEL_HOME
    cd $INTEL_HOME
    wget https://registrationcenter-download.intel.com/akdlm/IRC_NAS/e974de81-57b7-4ac1-b039-0512f8df974e/intel-oneapi-hpc-toolkit-2025.2.0.575_offline.sh # 2GB
    chmod 744 intel-oneapi-hpc-toolkit-2025.2.0.575_offline.sh
    ./intel-oneapi-hpc-toolkit-2025.2.0.575_offline.sh -a -s --eula accept --install-dir $INTEL_HOME
    echo Installed intel in $INTEL_HOME
fi

#------
# CMAKE
#------

read -p "Install CMAKE (y/n) " CMAKE_FLAG

if [ "$CMAKE_FLAG" = "y" ]
then
    cd $SOFT_ROOT
    wget https://github.com/Kitware/CMake/releases/download/v3.30.3/$CMAKE_VERSION.sh
    sh $CMAKE_VERSION.sh
    rm $CMAKE_VERSION.sh
    chmod u+x $CMAKE
    echo cmake installed in $SOFT_ROOT/$CMAKE_VERSION
fi

#-------------
# json-fortran
#-------------

read -p "Install JSON (y/n) " JSON_FLAG

if [ "$JSON_FLAG" = "y" ]
then
    cd $SOFT_ROOT
    JSON_ROOT=$SOFT_ROOT/json-fortran-${JSON_VERSION}
    rm -rf $JSON_ROOT
    wget https://github.com/jacobwilliams/json-fortran/archive/refs/tags/${JSON_VERSION}.tar.gz
    tar xvf ${JSON_VERSION}.tar.gz
    cd $JSON_ROOT
    mkdir -p build
    $CMAKE -S $JSON_ROOT $JSON_ROOT/build
    cd build
    make

    mkdir -p $JSON_INSTALL/lib
    mkdir -p $JSON_INSTALL/inc
    cp $JSON_ROOT/build/lib/libjsonfortran.a $JSON_INSTALL/lib/
    cp $JSON_ROOT/build/*.mod $JSON_INSTALL/inc/
    cp $AWD/platform/env.$platform $JSON_INSTALL/
    rm $SOFT_ROOT/${JSON_VERSION}.tar.gz
    echo JSON built in $JSON_ROOT installed in $JSON_INSTALL
fi

#-------
# NetCDF
#-------

read -p "Install NetCDF (y/n) " NETCDF_FLAG

if [ "$NETCDF_FLAG" = "y" ]
then
    cd $NETCDF_INSTALL
    FFLAGS_IN=$FFLAGS
    if [[ "$FC" == "ifx" ]]; then
        export CXX=icpx
        export FFLAGS="-O2 $FFLAGS_IN"
        export CFLAGS="-O2"
    fi
    export LD_LIBRARY_PATH="$NETCDF_INSTALL/lib:$LD_LIBRARY_PATH"
    export CPPFLAGS="-I$NETCDF_INSTALL/include"
    export LDFLAGS="-L$NETCDF_INSTALL/lib"

    ZLIB_VER=1.3.2
    HDF5_VER=1.14.3
    NETCDFC_VER=4.9.2
    NETCDFF_VER=4.6.1

    echo "==== Building zlib ===="
    curl -L -o zlib.tar.gz https://zlib.net/zlib-$ZLIB_VER.tar.gz
    tar --same-permissions -xf zlib.tar.gz
    cd zlib-$ZLIB_VER
    chmod +x configure
    ./configure --prefix="$NETCDF_INSTALL"
    make -j$NPROC
    make install
    cd ..

    echo "==== Building HDF5 serial ===="
    curl -L -o hdf5.tar.gz https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.14/hdf5-$HDF5_VER/src/hdf5-$HDF5_VER.tar.gz
    tar --same-permissions -xf hdf5.tar.gz
    cd hdf5-$HDF5_VER
    chmod +x configure
    ./configure --prefix="$NETCDF_INSTALL" --enable-hl --with-pthread=yes --with-zlib="$NETCDF_INSTALL" --enable-shared
    make -j$NPROC
    make install
    cd ..

    echo "==== Building NetCDF-C ===="
    curl -L -o netcdf-c.tar.gz https://github.com/Unidata/netcdf-c/archive/refs/tags/v$NETCDFC_VER.tar.gz
    tar --same-permissions -xf netcdf-c.tar.gz
    cd netcdf-c-$NETCDFC_VER
    chmod +x configure
    ./configure --prefix="$NETCDF_INSTALL" --enable-netcdf-4 --disable-dap
    make -j$NPROC
    make install
    cd ..

    echo "==== Building NetCDF-Fortran ===="
    curl -L -o netcdf-fortran.tar.gz https://github.com/Unidata/netcdf-fortran/archive/refs/tags/v$NETCDFF_VER.tar.gz
    tar -xf netcdf-fortran.tar.gz
    cd netcdf-fortran-$NETCDFF_VER
    chmod +x configure
    ./configure --prefix="$NETCDF_INSTALL"
    make -j$NPROC
    make install
    cd ..

    export FFLAGS=FFLAGS_IN

    echo "========================"
    echo "Build complete!"
    echo "Libraries in $NETCDF_INSTALL/lib"
    echo "========================"
fi

#-------
# RABBIT
#-------

read -p "Install RABBIT (y/n) " RABBIT_FLAG

if [ "$RABBIT_FLAG" = "y" ]
then
    cd $SOFT_ROOT
    rm -rf rabbit

    export LD_LIBRARY_PATH="$NETCDF_INSTALL/lib:$LD_LIBRARY_PATH"
    git clone git@gitlab.mpcdf.mpg.de:markusw/rabbit.git
    RABBIT_HOME=$SOFT_ROOT/rabbit
    cd $RABBIT_HOME
    RABBIT_HASH=`git rev-parse HEAD`
    mkdir build
    cd build
    $CMAKE .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=$FC -DOpenMP_Fortran_FLAGS=$FFLAGS -DNETCDF_HOME=$NETCDF_INSTALL
#    $CMAKE .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=$FC -DOpenMP_Fortran_FLAGS=$FFLAGS -DNETCDF_HOME=selfmade
    make

    mkdir -p $RABBIT_INSTALL/lib
    mkdir -p $RABBIT_INSTALL/inc
    cp $RABBIT_HOME/build/librabbit.so $RABBIT_INSTALL/lib/
    cp $RABBIT_HOME/build/modules/*.mod $RABBIT_INSTALL/inc/
    cp $AWD/platform/env.$platform $RABBIT_INSTALL/
    echo $RABBIT_HASH | cat > $RABBIT_INSTALL/hash
    echo "RABBIT built in $RABBIT_HOME installed in $RABBIT_INSTALL"
fi

#--------
# TORBEAM
#--------

read -p "Install TORBEAM (y/n) " TORBEAM_FLAG

if [ "$TORBEAM_FLAG" = "y" ]
then
    cd $SOFT_ROOT
    rm -rf torbeam
    git clone git@gitlab.mpcdf.mpg.de:ipp-aug/torbeam
    TORBEAM_HOME=$SOFT_ROOT/torbeam
    cd $TORBEAM_HOME
    TORBEAM_HASH=`git rev-parse HEAD`
    which cmake
    make

    mkdir -p $TORBEAM_INSTALL/lib
    cp $TORBEAM_HOME/build-generic/lib/libtorbeamB.so $TORBEAM_INSTALL/lib
    cp $AWD/platform/env.$platform $TORBEAM_INSTALL/
    echo $TORBEAM_HASH | cat > $TORBEAM_INSTALL/hash
    echo TORBEAM built in $TORBEAM_HOME installed in $TORBEAM_INSTALL
fi


#-------
# SPIDER
#-------

read -p "Install SPIDER (y/n) " SPIDER_FLAG

if [ "$SPIDER_FLAG" = "y" ]
then
    cd $SOFT_ROOT
    rm -rf spider
    git clone git@gitlab.mpcdf.mpg.de:git/spider
    SPIDER_HOME=$SOFT_ROOT/spider
    cd $SPIDER_HOME
    SPIDER_HASH=`git rev-parse HEAD`
    make

    mkdir -p $SPIDER_INSTALL/lib
    mkdir -p $SPIDER_INSTALL/inc
    cp $SPIDER_HOME/lib/libspider.a $SPIDER_INSTALL/lib/
    cp $SPIDER_HOME/inc/spider_params.mod $SPIDER_INSTALL/inc/
    cp $AWD/platform/env.$platform $SPIDER_INSTALL/
    echo $SPIDER_HASH | cat > $SPIDER_INSTALL/hash
    echo SPIDER built in $SPIDER_HOME installed in $SPIDER_INSTALL
fi

#---------
# QuaLiKiz
#---------

echo $MPIFC
which $MPIFC
read -p "Install QuaLiKiz (y/n) " QLK_FLAG

if [ "$QLK_FLAG" = "y" ]
then
    cd $SOFT_ROOT
    rm -rf QuaLiKiz
    git clone https://gitlab.com/qualikiz-group/QuaLiKiz.git
    QLK_HOME=$SOFT_ROOT/QuaLiKiz
    cd $QLK_HOME
    QLK_HASH=`git rev-parse HEAD`
    git submodule init
    git submodule update

    if [ "$comp" = "gcc" ]
    then
        export TOOLCHAIN=gcc
    fi
    export FC=$MPIFC
    export LINK=$MPIFC
    export QLK_HAVE_NAG=0
    export TUBSCFG_MPI=0
    export VERBOSE=1
    export BUILD=release

    make

    mkdir -p $QLK_INSTALL/lib
    mkdir -p $QLK_INSTALL/inc
    cp $QLK_HOME/lib/libQLK-intel-release-default.a $QLK_INSTALL/lib/
    cp $QLK_HOME/include/intel-release-default/* $QLK_INSTALL/inc/
    cp $AWD/platform/env.$platform $QLK_INSTALL/
    echo $QLK_HASH | cat > $QLK_INSTALL/hash
    echo QuaLiKiZ built in $QLK_HOME installed in $QLK_INSTALL
fi
#------------
# QuaLiKiz-NN
#------------

read -p "Install QuaLiKiz NN (y/n) " QLKNN_FLAG

if [ "$QLKNN_FLAG" = "y" ]
then
    cd $SOFT_ROOT
    rm -rf QLKNN-fortran
    git clone https://gitlab.com/qualikiz-group/QLKNN-fortran.git
    QLKNN_HOME=$SOFT_ROOT/QLKNN-fortran
    cd $QLKNN_HOME
    QLKNN_HASH=`git rev-parse HEAD`
    git submodule init
    git submodule update

    if [ "$comp" = "gcc" ]
    then
        export TOOLCHAIN=gcc
    fi

    export FC=$FC_SERIAL
    export LINK=$FC_SERIAL
    export QLK_HAVE_NAG=0
    export TUBSCFG_MPI=0
    export VERBOSE=1
    export BUILD=release

    make

    mkdir -p $QLKNN_INSTALL/lib
    mkdir -p $QLKNN_INSTALL/inc
    if [ "$comp" = "gcc" ]
    then
        cp $QLKNN_HOME/lib/libQLKNN-gcc-release-default.a $QLKNN_INSTALL/lib
        cp $QLKNN_HOME/include/gcc-release-default/* $QLKNN_INSTALL/inc/
    else
        cp $QLKNN_HOME/lib/libQLKNN-intel-release-default.a $QLKNN_INSTALL/lib
        cp $QLKNN_HOME/include/intel-release-default/* $QLKNN_INSTALL/inc/
    fi
    cp $AWD/platform/env.$platform $QLKNN_INSTALL/
    cd $QLKNN_INSTALL/
    rm -rf qlknn-hyper-namelists
    rm -rf qlknn-jetexp-namelists
    rm -rf qlknn-hornnet-namelists
    rm -rf qlknn-fullflux-namelists
    git clone https://gitlab.com/qualikiz-group/qlknn-hyper-namelists.git
    git clone https://gitlab.com/qualikiz-group/qlknn-jetexp-namelists.git
    git clone https://gitlab.com/qualikiz-group/qlknn-hornnet-namelists.git
    git clone https://gitlab.com/qualikiz-group/qlknn-fullflux-namelists.git
    echo $QLKNN_HASH | cat > $QLKNN_INSTALL/hash
    echo QuaLiKiZ-NN built in $QLKNN_HOME installed in $QLKNN_INSTALL
fi

#-------
# GACODE
#-------

read -p "Install TGLF/NEO (y/n) " GA_FLAG

if [ "$GA_FLAG" = "y" ]
then
    cd $SOFT_ROOT
    rm -rf gacode
    git clone git@github.com:gafusion/gacode.git
    export GACODE_ROOT=$SOFT_ROOT/gacode
    export GACODE_PLATFORM=MYLOC
    cd $GACODE_ROOT
    GACODE_HASH=`git rev-parse HEAD`

    if [ "$FC" = 'gfortran' ]
    then
        cat << EOT > ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}
IDENTITY="IPP linux cluster"
CORES_PER_NODE=16
NUMAS_PER_NODE=1

FC  = ${MPIFC} -J${GACODE_ROOT}/modules
F77 = ${FC}
FOMP   = ${FFLAGS}
FMATH  =-fdefault-real-8
FOPT   =-Ofast -fallow-argument-mismatch -Wno-error
FDEBUG =-eD -Ktrap=fp -m 1
LMATH = -L${MKLROOT}/lib/intel64 -Wl,-rpath,${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm -ldl
FFTW_INC=${FFTW_INC}
ARCH = ar cr
EOT
    else # Intel
        cat << EOT > ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}
IDENTITY="IPP linux cluster"
CORES_PER_NODE=16
NUMAS_PER_NODE=1

FC  = ${MPIFC} -module ${GACODE_ROOT}/modules
F77 = ${FC}
FOMP   = ${FFLAGS}
FMATH  =-real-size 64
FOPT   =-Ofast
FDEBUG =-eD -Ktrap=fp -m 1
LMATH = -qmkl
FFTW_INC=${FFTW_INC}
ARCH = ar cr
EOT
    fi

    cd $GACODE_ROOT/tglf
    make
    cd $GACODE_ROOT/neo
    make

    mkdir -p $TGLF_INSTALL/lib
    mkdir -p $TGLF_INSTALL/inc
    mkdir -p $NEO_INSTALL/lib
    mkdir -p $NEO_INSTALL/inc
    cp $GACODE_ROOT/tglf/src/tglf_lib.a $TGLF_INSTALL/lib/libtglf.a
    cp $GACODE_ROOT/neo/src/neo_lib.a $NEO_INSTALL/lib/
    cp $GACODE_ROOT/shared/UMFPACK/UMFPACK_lib.a $NEO_INSTALL/lib/
    cp $GACODE_ROOT/shared/math/math_lib.a $NEO_INSTALL/lib/
    cp $GACODE_ROOT/shared/nclass/nclass_lib.a $NEO_INSTALL/lib/
    cp $GACODE_ROOT/f2py/geo/geo_lib.a $NEO_INSTALL/lib/
    cp $GACODE_ROOT/f2py/expro/expro_lib.a $NEO_INSTALL/lib/
    cp $GACODE_ROOT/modules/tglf*.mod $TGLF_INSTALL/inc/
    cp $GACODE_ROOT/modules/neo_interface.mod $NEO_INSTALL/inc/.
    cp $AWD/platform/env.$platform $TGLF_INSTALL/
    cp $AWD/platform/env.$platform $NEO_INSTALL/
    echo $GACODE_HASH | cat > $TGLF_INSTALL/hash
    echo $GACODE_HASH | cat > $NEO_INSTALL/hash
    echo GACODE built in $GACODE_ROOT
    echo TGLF installed in $TGLF_INSTALL
    echo NEO installed  in $NEO_INSTALL
fi

#-------
# STRAHL
#-------

read -p "Install STRAHL (y/n) " STRAHL_FLAG

if [ "$STRAHL_FLAG" = "y" ]
then
    STRAHL_HOME=$SOFT_ROOT/strahl
    cd $SOFT_ROOT
    rm -rf strahl
    git clone git@gitlab.mpcdf.mpg.de:rld/strahl.git
    cd $STRAHL_HOME
    STRAHL_HASH=`git rev-parse HEAD`
    rm source/compile/machine
    cat << EOT > source/compile/machine
F90C=ifort
#the compiler flags
FFLAGS=-u -m64 -O
NCDFLIB=$ASTRA_EXT/netcdf/mar26/libnetcdff.a $ASTRA_EXT/netcdf/mar26/libnetcdf.a
F90FLAGS=-m64 -O -fPIC
#the flags for the shared library
SHAREDFLAGS=-G -fPIC -B symbolic -zdefs
LIB=-lc -lm
BIN=.
EOT
    cd $STRAHL_HOME/source/strahl
    make
    make result_to_astra
    mkdir -p $STRAHL_INSTALL/bin
    cp strahl $STRAHL_INSTALL/bin/
    cp result_to_astra $STRAHL_INSTALL/bin/
    echo $STRAHL_HASH | cat > $STRAHL_INSTALL/hash
    echo STRAHL built in $STRAHL_HOME installed in $STRAHL_INSTALL
fi

#------
# ASTRA
#------

#read -p "Install ASTRA (y/n) " ASTRA_FLAG
ASTRA_FLAG=n
if [ "$ASTRA_FLAG" = "y" ]
then
    cd $AWD
    chmod u+x install.sh
    ./install.sh
fi

chmod -R a+rx $ASTRA_EXT

# Restore initial paths

export PATH=$PATH_OLD
