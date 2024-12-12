#!/bin/bash

rootdir=`dirname $0`      # may be relative path
export AWD=`cd $rootdir && pwd`  # ensure absolute path

chmod 744 $AWD/get_platform
platform=`$AWD/get_platform`

#--------------------
# User dependent part
#--------------------

source $AWD/platform/env.$platform
FC_SERIAL=$FC

#------------------
# Packages versions
#------------------

CMAKE_VERSION=cmake-3.30.3-linux-x86_64
JSON_VERSION=9.0.2
RABBIT_VERSION=unstable
TORBEAM_VERSION=unstable
QLK_VERSION=unstable
QLKNN_VERSION=unstable
GA_VERSION=unstable

#--------------------
# ASTRA install paths
#--------------------

SOFT_ROOT=$HOME/soft
JSON_INSTALL=$ASTRA_EXT/json/$JSON_VERSION
RABBIT_INSTALL=$ASTRA_EXT/rabbit/$RABBIT_VERSION
TORBEAM_INSTALL=$ASTRA_EXT/torbeam/$TORBEAM_VERSION
QLK_INSTALL=$ASTRA_EXT/qualikiz/$QLK_VERSION
QLKNN_INSTALL=$ASTRA_EXT/qlk_nn/$QLKNN_VERSION
TGLF_INSTALL=$ASTRA_EXT/tglf/$GA_VERSION
NEO_INSTALL=$ASTRA_EXT/neo/$GA_VERSION

# To restore at the script end
PATH_OLD=$PATH

export PATH=$SOFT_ROOT/$CMAKE_VERSION/bin:$PATH
export FFLAGS="-qopenmp"
CMAKE=$SOFT_ROOT/$CMAKE_VERSION/bin/cmake
mkdir -p $SOFT_ROOT

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
    echo Installed cmake in $SOFT_ROOT/$CMAKE_VERSION
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
    echo Installed JSON in $JSON_ROOT
fi

#-------
# NetCDF
#-------

read -p "Install NetCDF (y/n) " NETCDF_FLAG
if [ "$NETCDF_FLAG" = "y" ]
then
    cd $SOFT_ROOT
    rm -rf netCDF
    git clone https://github.com/erdc/netCDF.git
    cd netCDF
    chmod u+x configure
    ./configure --with-pic --disable-netcdf-4 --disable-dap
    make
    echo Installed NetCDF in $SOFT_ROOT/netCDF
fi

#-------
# RABBIT
#-------

read -p "Install RABBIT (y/n) " RABBIT_FLAG
if [ "$RABBIT_FLAG" = "y" ]
then
    cd $SOFT_ROOT
    rm -rf rabbit
# git clone https://gitlab.mpcdf.mpg.de/markusw/rabbit
    git clone git@gitlab.mpcdf.mpg.de:markusw/rabbit.git
    RABBIT_HOME=$SOFT_ROOT/rabbit
    cd $RABBIT_HOME
    RABBIT_HASH=`git rev-parse HEAD`
    mkdir build
    cd build
    $CMAKE .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=$FC -DOpenMP_Fortran_FLAGS=-qopenmp -DNETCDF_HOME=$SOFT_ROOT/netCDF
    make

    mkdir -p $RABBIT_INSTALL/lib
    mkdir -p $RABBIT_INSTALL/inc
    cp $RABBIT_HOME/build/librabbit.so $RABBIT_INSTALL/lib/
    cp $RABBIT_HOME/build/modules/*.mod $RABBIT_INSTALL/inc/
    cp $AWD/platform/env.$platform $RABBIT_INSTALL/
    echo $RABBIT_HASH | cat > $RABBIT_INSTALL/hash
    echo Installed RABBIT in $RABBIT_HOME
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
    echo Installed TORBEAM in $TORBEAM_HOME
fi

#---------
# QuaLiKiz
#---------

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
    export FC=$MPIFC
    export LINK=$MPIFC
    export QLK_HAVE_NAG=0
    export TUBSCFG_MPI=0
    make

    mkdir -p $QLK_INSTALL/lib
    mkdir -p $QLK_INSTALL/inc
    cp $QLK_HOME/lib/libQLK-intel-release-default.a $QLK_INSTALL/lib/
    cp $QLK_HOME/include/intel-release-default/* $QLK_INSTALL/inc/
    cp $AWD/platform/env.$platform $QLK_INSTALL/
    echo $QLK_HASH | cat > $QLK_INSTALL/hash
    echo Installed qualikiz in $QLK_HOME
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
    export FC=$FC_SERIAL
    export LINK=$FC_SERIAL
    export QLK_HAVE_NAG=0
    export TUBSCFG_MPI=0
    make

    mkdir -p $QLKNN_INSTALL/lib
    mkdir -p $QLKNN_INSTALL/inc
    cp $QLKNN_HOME/lib/libQLKNN-intel-release-default.a $QLKNN_INSTALL/lib
    cp $QLKNN_HOME/include/intel-release-default/* $QLKNN_INSTALL/inc/
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
    echo Installed qualikiz in $QLKNN_HOME
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
    sed -i "s#zgeev\ =\ .false#zgeev\ =\ .true#g" $GACODE_ROOT/tglf/src/tglf_eigensolver.f90

    cat << EOT > ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}
IDENTITY="IPP linux cluster"
CORES_PER_NODE=16
NUMAS_PER_NODE=1

FC  = ${MPIFC} -module ${GACODE_ROOT}/modules
F77 = ${FC}
FOMP   =-qopenmp
FMATH  =-real-size 64
FOPT   =-Ofast
FDEBUG =-eD -Ktrap=fp -m 1
LMATH = -qmkl -mkl
FFTW_INC=${FFTW_INC}
ARCH = ar cr
EOT

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
    echo Installed GACODE in $GACODE_ROOT
fi

#------
# ASTRA
#------

read -p "Install ASTRA (y/n) " ASTRA_FLAG
if [ "$ASTRA_FLAG" = "y" ]
then
    cd $AWD
    chmod u+x install.sh
    ./install.sh
fi

chmod -R a+rx $ASTRA_EXT

# Restore initial paths

export PATH=$PATH_OLD
