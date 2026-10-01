#!/bin/bash

rootdir=`dirname $0`      # may be relative path
export AWD=`cd $rootdir && pwd`  # ensure absolute path
if [ $# -ge 1 ]
then
    comp=$1
fi

chmod 744 $AWD/get_platform
platform=`$AWD/get_platform`
if [[ -n ${comp+x} ]] && [ -f "$AWD/platform/env.${platform}_${comp}" ]
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

mkdir -p $ASTRA_EXT
JSON_INSTALL=$ASTRA_EXT/json
NETCDF_INSTALL=$ASTRA_EXT/netcdf
NETCDFF_INSTALL=$ASTRA_EXT/netcdf_fortran
RABBIT_INSTALL=$ASTRA_EXT/rabbit
TORBEAM_INSTALL=$ASTRA_EXT/torbeam
QLK_INSTALL=$ASTRA_EXT/qualikiz
QLKNN_INSTALL=$ASTRA_EXT/qlk_nn
TGLF_INSTALL=$ASTRA_EXT/tglf
NEO_INSTALL=$ASTRA_EXT/neo
STRAHL_INSTALL=$ASTRA_EXT/strahl
export FFLAGS="-fopenmp"

# To restore at the script end
PATH_OLD=$PATH

if [[ "$platform" == "linux" ]]; then

    if ! command -v apt-get >/dev/null 2>&1; then
        echo "apt-get is necessary on Linux, interrupting installation script"
        exit 1
    fi

    if dpkg -s libnetcdf-fortran-dev >/dev/null 2>&1; then
        NETCDFF_PKG=libnetcdf-fortran-dev
    else
        NETCDFF_PKG=libnetcdff-dev
    fi

    apt_pkg_prompt() {
        pkg=$1
        label=$2
        required=$3
        if dpkg -s "$pkg" >/dev/null 2>&1; then
            echo "$label already installed"
            return
        fi
        read -p "Install $label with apt (y/n) " APT_FLAG
        if [ "$APT_FLAG" = "y" ]; then
            sudo apt-get install -y "$pkg"
            echo "$label installed"
            return
        fi
        if [ "$required" = "required" ]; then
            echo "$label is necessary on Linux, interrupting installation script"
            exit 1
        fi
    }

    apt_pkg_prompt gcc GCC required
    apt_pkg_prompt cmake CMAKE required
    apt_pkg_prompt libomp-dev libomp optional
    apt_pkg_prompt libopenmpi-dev open-mpi optional
    apt_pkg_prompt libopenblas-dev openblas optional

    CMAKE_MIN_VERSION=3.18.0
    CMAKE_PATH=""
    CMAKE_FOUND_VERSION=""
    if command -v cmake >/dev/null 2>&1; then
        CMAKE_PATH=$(which cmake)
        CMAKE_FOUND_VERSION=$($CMAKE_PATH --version | awk 'NR==1{print $3}')
    fi
    if [ -n "$CMAKE_FOUND_VERSION" ] && [ "$(printf '%s\n%s\n' "$CMAKE_MIN_VERSION" "$CMAKE_FOUND_VERSION" | sort -V | head -n 1)" = "$CMAKE_MIN_VERSION" ]; then
        CMAKE=$CMAKE_PATH
    else
        mkdir -p $SOFT_ROOT
        if [ "$(uname -m)" = "aarch64" ]; then
            CMAKE_VERSION=cmake-3.30.3-linux-aarch64
        else
            CMAKE_VERSION=cmake-3.30.3-linux-x86_64
        fi
        CMAKE_HOME=$SOFT_ROOT/$CMAKE_VERSION
        CMAKE=$CMAKE_HOME/bin/cmake
        if [ ! -x "$CMAKE" ]; then
            cd $SOFT_ROOT
            mkdir -p $CMAKE_HOME
            wget https://github.com/Kitware/CMake/releases/download/v3.30.3/$CMAKE_VERSION.sh
            status=$?
            if [ $status -ne 0 ]; then
                echo "Error downloading CMAKE, interrupting installation script"
                exit 1
            fi
            sh $CMAKE_VERSION.sh --skip-license --prefix=$CMAKE_HOME
            status=$?
            rm -f $CMAKE_VERSION.sh
            if [ $status -ne 0 ] || [ ! -x "$CMAKE" ]; then
                echo "Error installing CMAKE, interrupting installation script"
                exit 1
            fi
        fi
        CMAKE_PATH=$CMAKE
    fi
    echo cmake installed in $CMAKE_PATH

    read -p "Install JSON (y/n) " JSON_FLAG

    if [ "$JSON_FLAG" = "y" ]; then
        cd /tmp
        JSON_ROOT=/tmp/json-fortran-${JSON_VERSION}
        rm -rf json-fortran-${JSON_VERSION}
        wget https://github.com/jacobwilliams/json-fortran/archive/refs/tags/${JSON_VERSION}.tar.gz
        tar xvf ${JSON_VERSION}.tar.gz
        cd $JSON_ROOT
        mkdir -p build
        $CMAKE -S $JSON_ROOT $JSON_ROOT/build
        status=$?
        if [ $status -ne 0 ]; then
            echo "Error configuring JSON, interrupting installation script"
            exit 1
        fi
        cd build
        make
        status=$?
        if [ $status -ne 0 ]; then
            echo "Error building JSON, interrupting installation script"
            exit 1
        fi
        if [ ! -f $JSON_ROOT/build/lib/libjsonfortran.a ] || ! compgen -G "$JSON_ROOT/build/*.mod" >/dev/null; then
            echo "JSON build did not produce expected library/module files, interrupting installation script"
            exit 1
        fi
        echo Installed JSON in $JSON_ROOT

        rm -rf $JSON_INSTALL
        mkdir -p $JSON_INSTALL/lib
        mkdir -p $JSON_INSTALL/inc
        cp $JSON_ROOT/build/lib/libjsonfortran.a $JSON_INSTALL/lib/
        cp $JSON_ROOT/build/*.mod $JSON_INSTALL/inc/
        cp $AWD/platform/env.$platform $JSON_INSTALL/
        rm "/tmp/${JSON_VERSION}.tar.gz"
        echo JSON installed in $JSON_ROOT
    fi

    if dpkg -s libnetcdf-dev >/dev/null 2>&1 && dpkg -s "$NETCDFF_PKG" >/dev/null 2>&1; then
        NETCDF_FLAG=y
    else
        read -p "Install NetCDF (y/n) " NETCDF_FLAG
    fi

    if [ "$NETCDF_FLAG" = "y" ]; then
        if ! dpkg -s libnetcdf-dev >/dev/null 2>&1; then
            sudo apt-get install -y libnetcdf-dev
        fi
        if ! dpkg -s "$NETCDFF_PKG" >/dev/null 2>&1; then
            sudo apt-get install -y "$NETCDFF_PKG"
        fi
        if ! dpkg -s libhdf5-dev >/dev/null 2>&1; then
            sudo apt-get install -y libhdf5-dev
        fi

        ARCH="$(uname -m)"
        if [[ "$ARCH" == "x86_64" ]]; then
            LIB_ARCH="x86_64-linux-gnu"
        elif [[ "$ARCH" == "aarch64" ]]; then
            LIB_ARCH="aarch64-linux-gnu"
        else
            echo "Unsupported architecture: $ARCH"
            exit 1
        fi

        NETCDF_SYS="/usr"
        NETCDFF_SYS="/usr"

        rm -rf $NETCDF_INSTALL
        mkdir -p $NETCDF_INSTALL/include
        mkdir -p $NETCDF_INSTALL/bin
        mkdir -p $NETCDF_INSTALL/lib
        cp -r $NETCDF_SYS/include/netcdf* $NETCDF_INSTALL/include/
        cp -r $NETCDF_SYS/bin/nc-config $NETCDF_INSTALL/bin/ 2>/dev/null || true
        cp -r $NETCDF_SYS/bin/ncdump $NETCDF_INSTALL/bin/ 2>/dev/null || true
        cp -r $NETCDF_SYS/bin/ncgen* $NETCDF_INSTALL/bin/ 2>/dev/null || true
        cp -r $NETCDF_SYS/lib/$LIB_ARCH/libnetcdf* $NETCDF_INSTALL/lib/
        cp $AWD/platform/env.$platform $NETCDF_INSTALL/
        echo "NetCDF library copied to $NETCDF_INSTALL"

        rm -rf $NETCDFF_INSTALL
        mkdir -p $NETCDFF_INSTALL/include
        mkdir -p $NETCDFF_INSTALL/bin
        mkdir -p $NETCDFF_INSTALL/lib
        cp -r $NETCDFF_SYS/include/netcdf* $NETCDFF_INSTALL/include/
        cp -r $NETCDFF_SYS/bin/nf-config $NETCDFF_INSTALL/bin/ 2>/dev/null || true
        cp -r $NETCDFF_SYS/lib/$LIB_ARCH/libnetcdff* $NETCDFF_INSTALL/lib/
        cp $AWD/platform/env.$platform $NETCDFF_INSTALL/
        if ! compgen -G "$NETCDF_INSTALL/lib/libnetcdf*" >/dev/null || ! compgen -G "$NETCDFF_INSTALL/lib/libnetcdff*" >/dev/null; then
            echo "NetCDF copy did not produce expected library files, interrupting installation script"
            exit 1
        fi

        echo "NetCDF-Fortran library copied to $NETCDFF_INSTALL"
        echo "========================"
        echo "Build complete!"
        echo "Libraries in $NETCDF_INSTALL/lib and $NETCDFF_INSTALL/lib"
        echo "========================"
    fi

    read -p "Install RABBIT (y/n) " RABBIT_FLAG

    if [ "$RABBIT_FLAG" = "y" ]; then
        cd /tmp
        rm -rf rabbit
        git clone git@gitlab.mpcdf.mpg.de:markusw/rabbit.git
        RABBIT_HOME=/tmp/rabbit
        cd $RABBIT_HOME
        mkdir build
        cd build
        env -u IMAS_VERSION $CMAKE .. -DCMAKE_BUILD_TYPE=release -DCMAKE_Fortran_COMPILER=$FC -DNETCDFF_HOME=$NETCDFF_INSTALL -DNETCDF_HOME=$NETCDF_INSTALL -DOpenMP_Fortran_FLAGS=-fopenmp
        status=$?
        if [ $status -ne 0 ]; then
            echo "Error configuring RABBIT, interrupting installation script"
            exit 1
        fi
        make
        status=$?
        if [ $status -ne 0 ]; then
            echo "Error building RABBIT, interrupting installation script"
            exit 1
        fi
        if [ ! -f $RABBIT_HOME/build/librabbit.so ] || ! compgen -G "$RABBIT_HOME/build/modules/*.mod" >/dev/null; then
            echo "RABBIT build did not produce expected library/module files, interrupting installation script"
            exit 1
        fi
        echo RABBIT installed in $RABBIT_HOME

        rm -rf $RABBIT_INSTALL
        mkdir -p $RABBIT_INSTALL/lib
        mkdir -p $RABBIT_INSTALL/inc
        cp $RABBIT_HOME/build/librabbit.so $RABBIT_INSTALL/lib/
        cp $RABBIT_HOME/build/modules/*.mod $RABBIT_INSTALL/inc/
        cp $AWD/platform/env.$platform $RABBIT_INSTALL/
        echo RABBIT compiled library copied to $RABBIT_INSTALL
    fi

    read -p "Install TORBEAM (y/n) " TORBEAM_FLAG

    if [ "$TORBEAM_FLAG" = "y" ]; then
        cd /tmp
        rm -rf torbeam
        git clone ssh://gerrit.ipp.mpg.de:29418/ipp/e1/ecrh/libtorbeam torbeam
        TORBEAM_HOME=/tmp/torbeam
        cd $TORBEAM_HOME
        if [[ "$(uname -m)" == "aarch64" ]]; then
            perl -i -0777 -pe 's|SET\(CMAKE_SHARED_LINKER_FLAGS "-m64 -z muldefs"\)|if(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")\n    SET(CMAKE_SHARED_LINKER_FLAGS "-m64 -z muldefs")\nelse()\n    SET(CMAKE_SHARED_LINKER_FLAGS "-z muldefs")\nendif()|' CMakeLists.txt
            perl -i -0777 -pe 's|set\(CMAKE_EXE_COMPILE_FLAGS_LINUX "-m64"\)|if(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")\n  set(CMAKE_EXE_COMPILE_FLAGS_LINUX "-m64")\nelse()\n  set(CMAKE_EXE_COMPILE_FLAGS_LINUX "")\nendif()|' src/exesrc/CMakeLists.txt
            perl -i -0777 -pe 's|    add_custom_command\(TARGET torbeam\$\{DEFAULT_OPT\} POST_BUILD COMMAND \$\{CMAKE_COMMAND\}\n                       -E copy libtorbeam\$\{DEFAULT_OPT\}\.so libtorbeam\.so\)|    add_custom_command(TARGET torbeam\${DEFAULT_OPT} POST_BUILD COMMAND \${CMAKE_COMMAND}\n                   -E copy libtorbeam\${DEFAULT_OPT}.so libtorbeam.so)|' src/libsrc/CMakeLists.txt
            perl -i -0777 -pe 's|set_target_properties\(\$\{EXE\}_\$\{Opt\} PROPERTIES COMPILE_FLAGS \$\{CMAKE_EXE_COMPILE_FLAGS\}\)|if(CMAKE_EXE_COMPILE_FLAGS)\n    set_target_properties(${EXE}_${Opt} PROPERTIES COMPILE_FLAGS ${CMAKE_EXE_COMPILE_FLAGS})\n  endif()|' src/exesrc/CMakeLists.txt
            perl -i -0777 -pe 's|SET\(COMMON_Fortran_FLAGS "-I\$\{SRCLIB\} -m64 -fPIC"\)|if(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")\n    SET(COMMON_Fortran_FLAGS "-I${SRCLIB} -m64 -fPIC")\nelse()\n    SET(COMMON_Fortran_FLAGS "-I${SRCLIB} -fPIC")\nendif()|' CMakeLists.txt
            perl -i -0777 -pe 's|set\(CMAKE_C_FLAGS "\$\{CMAKE_C_FLAGS\} -m64"\)|if(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")\n    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m64")\nendif()|' tools/c_wrapper/CMakeLists.txt
            perl -i -0777 -pe 's|set\(CMAKE_SHARED_LINKER_FLAGS "-m64 -z muldefs"\)|if(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")\n    set(CMAKE_SHARED_LINKER_FLAGS "-m64 -z muldefs")\nelse()\n    set(CMAKE_SHARED_LINKER_FLAGS "-z muldefs")\nendif()|' tools/c_wrapper/CMakeLists.txt
            perl -i -0777 -pe 's|PROPERTIES COMPILE_FLAGS \$\{CMAKE_EXE_COMPILE_FLAGS\} LINK_FLAGS "-m64"|if(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")\n    set_target_properties(wrapper.exe PROPERTIES COMPILE_FLAGS ${CMAKE_EXE_COMPILE_FLAGS} LINK_FLAGS "-m64")\n  else()\n    set_target_properties(wrapper.exe PROPERTIES COMPILE_FLAGS ${CMAKE_EXE_COMPILE_FLAGS})\n  endif()|' src/exesrc/CMakeLists.txt
        fi
        env -u IMAS_VERSION $CMAKE .
        status=$?
        if [ $status -ne 0 ]; then
            echo "Error configuring TORBEAM, interrupting installation script"
            exit 1
        fi
        make
        status=$?
        if [ $status -ne 0 ]; then
            echo "Error building TORBEAM, interrupting installation script"
            exit 1
        fi
        if [ ! -f $TORBEAM_HOME/lib/libtorbeamB.so ]; then
            echo "TORBEAM build did not produce expected library file, interrupting installation script"
            exit 1
        fi
        echo TORBEAM installed in $TORBEAM_HOME

        rm -rf $TORBEAM_INSTALL
        mkdir -p $TORBEAM_INSTALL/lib
        cp $TORBEAM_HOME/lib/libtorbeamB.so $TORBEAM_INSTALL/lib/
        cp $AWD/platform/env.$platform $TORBEAM_INSTALL/
        echo TORBEAM compiled library copied to $TORBEAM_INSTALL
    fi

    read -p "Install QuaLiKiz (y/n) " QLK_FLAG

    if [ "$QLK_FLAG" = "y" ]; then
        cd /tmp
        rm -rf QuaLiKiz
        git clone https://gitlab.com/qualikiz-group/QuaLiKiz.git
        QLK_HOME=/tmp/QuaLiKiz
        cd $QLK_HOME
        QLK_HASH=`git rev-parse HEAD`
        git submodule init
        git submodule update
        if [ -f fruit.make ]; then
            perl -0pi -e 's/\$\(call LOCAL_mod_dep, tests\/test_d01ahf\.f90, fruit_mpi\.mod\)/\$(call LOCAL_mod_dep, tests\/test_d01ahf.f90, fruit.mod)/g; s/\$\(call LOCAL_mod_dep, tests\/test_integration\.f90, fruit_mpi\.mod\)/\$(call LOCAL_mod_dep, tests\/test_integration.f90, fruit.mod)/g' fruit.make
        fi

        FC_IN=$FC
        LINK_IN=$LINK
        export TOOLCHAIN=gcc
        export FC=$MPIFC
        export LINK=$MPIFC
        export QLK_HAVE_NAG=0
        export TUBSCFG_MPI=0
        export VERBOSE=1
        export BUILD=release

        make
        status=$?
        if [ $status -ne 0 ]; then
            echo "Error building QuaLiKiz, interrupting installation script"
            exit 1
        fi

        rm -rf $QLK_INSTALL
        mkdir -p $QLK_INSTALL/lib
        mkdir -p $QLK_INSTALL/inc
        if [ -f $QLK_HOME/lib/libQLK-gcc-release-default.a ]
        then
            cp $QLK_HOME/lib/libQLK-gcc-release-default.a $QLK_INSTALL/lib/
            cp $QLK_HOME/include/gcc-release-default/* $QLK_INSTALL/inc/
        elif [ -f $QLK_HOME/lib/libQLK-intel-release-default.a ]
        then
            cp $QLK_HOME/lib/libQLK-intel-release-default.a $QLK_INSTALL/lib/
            cp $QLK_HOME/include/intel-release-default/* $QLK_INSTALL/inc/
        else
            echo "QuaLiKiZ build did not produce expected library/include files, interrupting installation script"
            exit 1
        fi
        cp $AWD/platform/env.$platform $QLK_INSTALL/
        echo $QLK_HASH | cat > $QLK_INSTALL/hash
        echo QuaLiKiZ built in $QLK_HOME installed in $QLK_INSTALL

        export FC=$FC_IN
        export LINK=$LINK_IN
    fi

    read -p "Install QuaLiKiz NN (y/n) " QLKNN_FLAG

    if [ "$QLKNN_FLAG" = "y" ]; then
        cd /tmp
        rm -rf QLKNN-fortran
        git clone https://gitlab.com/qualikiz-group/QLKNN-fortran.git
        QLKNN_HOME=/tmp/QLKNN-fortran
        cd $QLKNN_HOME
        QLKNN_HASH=`git rev-parse HEAD`
        git submodule init
        git submodule update
        if [ -f fruit.make ]; then
            perl -0pi -e 's/\$\(call LOCAL_mod_dep, tests\/test_d01ahf\.f90, fruit_mpi\.mod\)/\$(call LOCAL_mod_dep, tests\/test_d01ahf.f90, fruit.mod)/g; s/\$\(call LOCAL_mod_dep, tests\/test_integration\.f90, fruit_mpi\.mod\)/\$(call LOCAL_mod_dep, tests\/test_integration.f90, fruit.mod)/g' fruit.make
        fi

        FC_IN=$FC
        LINK_IN=$LINK
        export TOOLCHAIN=gcc
        export FC=$FC_SERIAL
        export LINK=$FC_SERIAL
        export QLK_HAVE_NAG=0
        export TUBSCFG_MPI=0
        export VERBOSE=1
        export BUILD=release

        make
        status=$?
        if [ $status -ne 0 ]; then
            echo "Error building QuaLiKiZ-NN, interrupting installation script"
            exit 1
        fi

        rm -rf $QLKNN_INSTALL
        mkdir -p $QLKNN_INSTALL/lib
        mkdir -p $QLKNN_INSTALL/inc
        if [ -f $QLKNN_HOME/lib/libQLKNN-gcc-release-default.a ]
        then
            cp $QLKNN_HOME/lib/libQLKNN-gcc-release-default.a $QLKNN_INSTALL/lib/
            cp $QLKNN_HOME/include/gcc-release-default/* $QLKNN_INSTALL/inc/
        elif [ -f $QLKNN_HOME/lib/libQLKNN-intel-release-default.a ]
        then
            cp $QLKNN_HOME/lib/libQLKNN-intel-release-default.a $QLKNN_INSTALL/lib/
            cp $QLKNN_HOME/include/intel-release-default/* $QLKNN_INSTALL/inc/
        else
            echo "QuaLiKiZ-NN build did not produce expected library/include files, interrupting installation script"
            exit 1
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

        export FC=$FC_IN
        export LINK=$LINK_IN
    fi

    read -p "Install TGLF/NEO (y/n) " GA_FLAG

    if [ "$GA_FLAG" = "y" ]; then
        cd /tmp
        rm -rf gacode
        git clone git@github.com:gafusion/gacode.git
        FC_IN=$FC
        MPIFC_IN=$MPIFC
        export GACODE_ROOT=/tmp/gacode
        export GACODE_PLATFORM=CONDA_OMPI_GNU
        CONDA_BASE=$(conda info --base 2>/dev/null || echo "")
        if [ -z "$CONDA_BASE" ]; then
            echo "Cannot find conda installation, interrupting"
            exit 1
        fi
        source "$CONDA_BASE/etc/profile.d/conda.sh"
        PREV_CONDA_ENV="${CONDA_DEFAULT_ENV:-}"
        conda deactivate 2>/dev/null || true
        conda env remove -n gacode -y >/dev/null 2>&1 || true
        conda create -n gacode -c conda-forge make python fftw netcdf4
        conda activate gacode
        GACODE_ARGMISMATCHFLAG=$($FC_IN --help=common 2>/dev/null | grep -q -- '-fallow-argument-mismatch' && echo -fallow-argument-mismatch || echo -Wno-argument-mismatch)
        perl -i -pe 's|FOPT\t= -Ofast -m64|FOPT\t= -Ofast \$(shell uname -m \| grep -q x86_64 \&\& echo -m64)|' $GACODE_ROOT/platform/build/make.inc.$GACODE_PLATFORM
        perl -i -pe "s|F77\t= \\\$\\{MF90\\} -g\\s*\$|F77\t= \\\$\\{MF90\\} -g ${GACODE_ARGMISMATCHFLAG}|" "$GACODE_ROOT/platform/build/make.inc.$GACODE_PLATFORM"
        . $GACODE_ROOT/shared/bin/gacode_setup
        . $GACODE_ROOT/platform/env/env.CONDA_OMPI_GNU
        export FC=$FC_IN
        export MPIFC=$MPIFC_IN
        export OMPI_FC=$FC_IN
        cd $GACODE_ROOT
        GACODE_HASH=$(git rev-parse HEAD)
        cd $GACODE_ROOT/tglf
        make
        status=$?
        if [ $status -ne 0 ]; then
            echo "Error building TGLF, interrupting installation script"
            exit 1
        fi
        cd $GACODE_ROOT/neo
        make
        status=$?
        if [ $status -ne 0 ]; then
            echo "Error building NEO, interrupting installation script"
            exit 1
        fi
        if [ ! -f $GACODE_ROOT/tglf/src/tglf_lib.a ] || [ ! -f $GACODE_ROOT/neo/src/neo_lib.a ] || [ ! -f $GACODE_ROOT/modules/neo_interface.mod ]; then
            echo "TGLF/NEO build did not produce expected library/module files, interrupting installation script"
            exit 1
        fi

        rm -rf $TGLF_INSTALL
        mkdir -p $TGLF_INSTALL/lib
        mkdir -p $TGLF_INSTALL/inc
        rm -rf $NEO_INSTALL
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
        echo GACODE installed in $GACODE_ROOT

        conda deactivate
        if [ -n "$PREV_CONDA_ENV" ]; then
            conda activate "$PREV_CONDA_ENV"
        fi
    fi

    read -p "Install STRAHL (y/n) " STRAHL_FLAG

    if [ "$STRAHL_FLAG" = "y" ]; then
        cd /tmp
        rm -rf strahl
        git clone git@gitlab.mpcdf.mpg.de:azito/strahl.git strahl
        STRAHL_HOME=/tmp/strahl
        perl -i -0777 -pe '
        s|FFLAGS \+= -fimplicit-none -m64 |FFLAGS += -fimplicit-none \$(shell uname -m \| grep -q x86_64 \&\& echo -m64) |;
        s|F90FLAGS \+= -m64 |F90FLAGS += \$(shell uname -m \| grep -q x86_64 \&\& echo -m64) |;
        s|LIB = -L/opt/homebrew/lib|LIB = -L\$(shell uname -m \| grep -q x86_64 \&\& echo /usr/lib/x86_64-linux-gnu/hdf5/serial \|\| echo /usr/lib/aarch64-linux-gnu/hdf5/serial)|
        ' "$STRAHL_HOME/source/compile/my_setup"
        cd $STRAHL_HOME/source/compile
        ./make_all my_setup
        status=$?
        if [ $status -ne 0 ]; then
            echo "Error building STRAHL, interrupting installation script"
            exit 1
        fi
        if [ ! -f $STRAHL_HOME/bin/strahl ] || [ ! -f $STRAHL_HOME/bin/result_to_astra ]; then
            echo "STRAHL build did not produce expected binaries, interrupting installation script"
            exit 1
        fi

        rm -rf $STRAHL_INSTALL
        mkdir -p $STRAHL_INSTALL/bin
        cp $STRAHL_HOME/bin/strahl $STRAHL_INSTALL/bin/
        cp $STRAHL_HOME/bin/result_to_astra $STRAHL_INSTALL/bin/
        cp $AWD/platform/env.$platform $STRAHL_INSTALL/
        echo Strahl installed in $STRAHL_HOME
    fi

elif [[ "$platform" == "darwin" ]]; then

    if ! command -v brew >/dev/null 2>&1; then
        echo "Homebrew is necessary on macOS, interrupting installation script"
        exit 1
    fi

    brew_pkg_prompt() {
        formula=$1
        label=$2
        required=$3

        if brew ls --versions $formula >/dev/null 2>&1; then
            echo "$label already installed in $(brew --prefix $formula)"
            return
        fi

        read -p "Install $label with brew (y/n) " BREW_FLAG
        if [ "$BREW_FLAG" = "y" ]; then
            brew install $formula
            echo "$label installed in $(brew --prefix $formula)"
            return
        fi

        if [ "$required" = "required" ]; then
            echo "$label is necessary on macOS, interrupting installation script"
            exit 1
        fi
    }

    brew_pkg_prompt gcc GCC required
    brew_pkg_prompt cmake CMAKE required
    brew_pkg_prompt make "GNU make" required
    brew_pkg_prompt libomp libomp optional
    brew_pkg_prompt open-mpi open-mpi optional
    brew_pkg_prompt openblas openblas optional
    brew_pkg_prompt netcdf NetCDF optional
    brew_pkg_prompt netcdf-fortran NetCDF-Fortran optional

    CMAKE_PATH=$(brew --prefix cmake)
    echo cmake installed in $CMAKE_PATH
    CMAKE=$CMAKE_PATH/bin/cmake
    GMAKE=$(brew --prefix make)/bin/gmake

    read -p "Install JSON (y/n) " JSON_FLAG

    if [ "$JSON_FLAG" = "y" ]; then
        brew install wget
        cd /tmp
        JSON_ROOT=/tmp/json-fortran-${JSON_VERSION}
        rm -rf json-fortran-${JSON_VERSION}
        wget https://github.com/jacobwilliams/json-fortran/archive/refs/tags/${JSON_VERSION}.tar.gz
        tar xvf ${JSON_VERSION}.tar.gz
        cd $JSON_ROOT
        mkdir -p build
        $CMAKE -S $JSON_ROOT $JSON_ROOT/build
        status=$?
        if [ $status -ne 0 ]; then
            echo "Error configuring JSON, interrupting installation script"
            exit 1
        fi
        cd build
        make
        status=$?
        if [ $status -ne 0 ]; then
            echo "Error building JSON, interrupting installation script"
            exit 1
        fi
        if [ ! -f $JSON_ROOT/build/lib/libjsonfortran.a ] || ! compgen -G "$JSON_ROOT/build/*.mod" >/dev/null; then
            echo "JSON build did not produce expected library/module files, interrupting installation script"
            exit 1
        fi
        echo Installed JSON in $JSON_ROOT

        rm -rf $JSON_INSTALL
        mkdir -p $JSON_INSTALL/lib
        mkdir -p $JSON_INSTALL/inc
        cp $JSON_ROOT/build/lib/libjsonfortran.a $JSON_INSTALL/lib/
        cp $JSON_ROOT/build/*.mod $JSON_INSTALL/inc/
        cp $AWD/platform/env.$platform $JSON_INSTALL/
        rm "/tmp/${JSON_VERSION}.tar.gz"
        echo JSON installed in $JSON_ROOT
    fi

    if brew ls --versions netcdf >/dev/null 2>&1 && brew ls --versions netcdf-fortran >/dev/null 2>&1; then
        NETCDF_FLAG=y
    else
        read -p "Install NetCDF (y/n) " NETCDF_FLAG
    fi

    if [ "$NETCDF_FLAG" = "y" ]; then
        if ! brew ls --versions netcdf >/dev/null 2>&1; then
            read -p "Install NetCDF with brew (y/n) " BREW_NETCDF_FLAG
            if [ "$BREW_NETCDF_FLAG" != "y" ]; then
                echo "NetCDF is necessary on macOS for this section, interrupting installation script"
                exit 1
            fi
            brew install netcdf
        fi
        NETCDF_PATH=$(brew --prefix netcdf)
        if ! brew ls --versions netcdf-fortran >/dev/null 2>&1; then
            read -p "Install NetCDF-Fortran with brew (y/n) " BREW_NETCDFF_FLAG
            if [ "$BREW_NETCDFF_FLAG" != "y" ]; then
                echo "NetCDF-Fortran is necessary on macOS for this section, interrupting installation script"
                exit 1
            fi
            brew install netcdf-fortran
        fi
        NETCDFF_PATH=$(brew --prefix netcdf-fortran)
        echo NetCDF installed in $NETCDF_PATH and NetCDF-Fortran installed in $NETCDFF_PATH

        rm -rf $NETCDF_INSTALL
        mkdir -p $NETCDF_INSTALL/include
        mkdir -p $NETCDF_INSTALL/bin
        mkdir -p $NETCDF_INSTALL/lib
        cp -r $NETCDF_PATH/include/* $NETCDF_INSTALL/include/
        cp -r $NETCDF_PATH/bin/* $NETCDF_INSTALL/bin/
        cp -r $NETCDF_PATH/lib/* $NETCDF_INSTALL/lib/
        cp $AWD/platform/env.$platform $NETCDF_INSTALL/
        echo NetCDF compiled library copied to $NETCDF_INSTALL

        rm -rf $NETCDFF_INSTALL
        mkdir -p $NETCDFF_INSTALL/include
        mkdir -p $NETCDFF_INSTALL/bin
        mkdir -p $NETCDFF_INSTALL/lib
        cp -r $NETCDFF_PATH/include/* $NETCDFF_INSTALL/include/
        cp -r $NETCDFF_PATH/bin/* $NETCDFF_INSTALL/bin/
        cp -r $NETCDFF_PATH/lib/* $NETCDFF_INSTALL/lib/
        cp $AWD/platform/env.$platform $NETCDFF_INSTALL/
        if ! compgen -G "$NETCDF_INSTALL/lib/libnetcdf*" >/dev/null || ! compgen -G "$NETCDFF_INSTALL/lib/libnetcdff*" >/dev/null; then
            echo "NetCDF copy did not produce expected library files, interrupting installation script"
            exit 1
        fi
        echo NetCDF-Fortran compiled library copied to $NETCDFF_INSTALL
    fi

    read -p "Install RABBIT (y/n) " RABBIT_FLAG

    if [ "$RABBIT_FLAG" = "y" ]; then
        cd /tmp
        rm -rf rabbit
        git clone git@gitlab.mpcdf.mpg.de:markusw/rabbit.git
        RABBIT_HOME=/tmp/rabbit
        cd $RABBIT_HOME
        mkdir build
        cd build
        env -u IMAS_VERSION $CMAKE .. -DCMAKE_BUILD_TYPE=release -DCMAKE_Fortran_COMPILER=$FC -DNETCDFF_HOME=$NETCDFF_INSTALL -DNETCDF_HOME=$NETCDF_INSTALL -DOpenMP_Fortran_FLAGS=-fopenmp
        status=$?
        if [ $status -ne 0 ]; then
            echo "Error configuring RABBIT, interrupting installation script"
            exit 1
        fi
        make
        status=$?
        if [ $status -ne 0 ]; then
            echo "Error building RABBIT, interrupting installation script"
            exit 1
        fi
        if [ ! -f $RABBIT_HOME/build/librabbit.dylib ] || ! compgen -G "$RABBIT_HOME/build/modules/*.mod" >/dev/null; then
            echo "RABBIT build did not produce expected library/module files, interrupting installation script"
            exit 1
        fi
        echo RABBIT installed in $RABBIT_HOME

        rm -rf $RABBIT_INSTALL
        mkdir -p $RABBIT_INSTALL/lib
        mkdir -p $RABBIT_INSTALL/inc
        cp $RABBIT_HOME/build/librabbit.dylib $RABBIT_INSTALL/lib/
        cp $RABBIT_HOME/build/modules/*.mod $RABBIT_INSTALL/inc/
        cp $AWD/platform/env.$platform $RABBIT_INSTALL/
        echo RABBIT compiled library copied to $RABBIT_INSTALL
    fi

    read -p "Install TORBEAM (y/n) " TORBEAM_FLAG

    if [ "$TORBEAM_FLAG" = "y" ]; then
        cd /tmp
        rm -rf torbeam
        git clone ssh://gerrit.ipp.mpg.de:29418/ipp/e1/ecrh/libtorbeam torbeam
        TORBEAM_HOME=/tmp/torbeam
        cd $TORBEAM_HOME
        perl -i -0777 -pe '
        s|SET\(CMAKE_SHARED_LINKER_FLAGS "-m64 -z muldefs"\)|if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")\n    SET(CMAKE_SHARED_LINKER_FLAGS "-m64")\nelse()\n    SET(CMAKE_SHARED_LINKER_FLAGS "-m64 -z muldefs")\nendif()|
        ' CMakeLists.txt
        perl -i -0777 -pe '
        s|if\(\$\{CMAKE_SYSTEM_NAME\} STREQUAL "Linux"\)\n  set\(CMAKE_EXE_COMPILE_FLAGS \$\{CMAKE_EXE_COMPILE_FLAGS_LINUX\}\)|if(\${CMAKE_SYSTEM_NAME} STREQUAL "Linux" OR \${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")\n  set(CMAKE_EXE_COMPILE_FLAGS \${CMAKE_EXE_COMPILE_FLAGS_LINUX})|
        ' src/exesrc/CMakeLists.txt
        perl -i -0777 -pe '
        s|    add_custom_command\(TARGET torbeam\$\{DEFAULT_OPT\} POST_BUILD COMMAND \$\{CMAKE_COMMAND\}\n                       -E copy libtorbeam\$\{DEFAULT_OPT\}\.so libtorbeam\.so\)|    if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")\n        add_custom_command(TARGET torbeam\${DEFAULT_OPT} POST_BUILD COMMAND \${CMAKE_COMMAND}\n                   -E copy libtorbeam\${DEFAULT_OPT}.dylib libtorbeam.dylib)\n    else()\n        add_custom_command(TARGET torbeam\${DEFAULT_OPT} POST_BUILD COMMAND \${CMAKE_COMMAND}\n                   -E copy libtorbeam\${DEFAULT_OPT}.so libtorbeam.so)\n    endif()|
        ' src/libsrc/CMakeLists.txt
        env -u IMAS_VERSION $CMAKE .
        status=$?
        if [ $status -ne 0 ]; then
            echo "Error configuring TORBEAM, interrupting installation script"
            exit 1
        fi
        make
        status=$?
        if [ $status -ne 0 ]; then
            echo "Error building TORBEAM, interrupting installation script"
            exit 1
        fi
        if [ ! -f $TORBEAM_HOME/lib/libtorbeamB.dylib ]; then
            echo "TORBEAM build did not produce expected library file, interrupting installation script"
            exit 1
        fi
        echo TORBEAM installed in $TORBEAM_HOME

        rm -rf $TORBEAM_INSTALL
        mkdir -p $TORBEAM_INSTALL/lib
        cp $TORBEAM_HOME/lib/libtorbeamB.dylib $TORBEAM_INSTALL/lib/
        cp $AWD/platform/env.$platform $TORBEAM_INSTALL/
        echo TORBEAM compiled library copied to $TORBEAM_INSTALL
    fi

    read -p "Install QuaLiKiz (y/n) " QLK_FLAG

    if [ "$QLK_FLAG" = "y" ]; then
        cd /tmp
        rm -rf QuaLiKiz
        git clone https://gitlab.com/qualikiz-group/QuaLiKiz.git
        QLK_HOME=/tmp/QuaLiKiz
        cd $QLK_HOME
        QLK_HASH=`git rev-parse HEAD`
        git submodule init
        git submodule update
        if [ -f fruit.make ]; then
            perl -0pi -e 's/\$\(call LOCAL_mod_dep, tests\/test_d01ahf\.f90, fruit_mpi\.mod\)/\$(call LOCAL_mod_dep, tests\/test_d01ahf.f90, fruit.mod)/g; s/\$\(call LOCAL_mod_dep, tests\/test_integration\.f90, fruit_mpi\.mod\)/\$(call LOCAL_mod_dep, tests\/test_integration.f90, fruit.mod)/g' fruit.make
        fi

        FC_IN=$FC
        LINK_IN=$LINK
        LNKGRPBEG_IN=$LNKGRPBEG
        LNKGRPEND_IN=$LNKGRPEND
        export TOOLCHAIN=gcc
        export FC=$MPIFC
        export LINK=$MPIFC
        export QLK_HAVE_NAG=0
        export TUBSCFG_MPI=0
        export VERBOSE=1
        export BUILD=release
        export LNKGRPBEG=
        export LNKGRPEND=
        PATH_IN=$PATH
        GCC_AR=$(ls "$(dirname "$CC")"/gcc-ar-* 2>/dev/null | sort -V | tail -n 1)
        GCC_NM=$(ls "$(dirname "$CC")"/gcc-nm-* 2>/dev/null | sort -V | tail -n 1)
        GCC_RANLIB=$(ls "$(dirname "$CC")"/gcc-ranlib-* 2>/dev/null | sort -V | tail -n 1)
        if [ -x "$GCC_AR" ]; then
            mkdir -p /tmp/astra_gcc_bin
            ln -sf "$GCC_AR" /tmp/astra_gcc_bin/gcc-ar
        fi
        if [ -x "$GCC_NM" ]; then
            mkdir -p /tmp/astra_gcc_bin
            ln -sf "$GCC_NM" /tmp/astra_gcc_bin/gcc-nm
        fi
        if [ -x "$GCC_RANLIB" ]; then
            mkdir -p /tmp/astra_gcc_bin
            ln -sf "$GCC_RANLIB" /tmp/astra_gcc_bin/gcc-ranlib
        fi
        if [ -d /tmp/astra_gcc_bin ]; then
            export PATH=/tmp/astra_gcc_bin:$PATH
        fi

        $GMAKE
        status=$?
        if [ $status -ne 0 ]; then
            echo "Error building QuaLiKiz, interrupting installation script"
            exit 1
        fi

        rm -rf $QLK_INSTALL
        mkdir -p $QLK_INSTALL/lib
        mkdir -p $QLK_INSTALL/inc
        if [ -f $QLK_HOME/lib/libQLK-gcc-release-default.a ]
        then
            cp $QLK_HOME/lib/libQLK-gcc-release-default.a $QLK_INSTALL/lib/
            cp $QLK_HOME/include/gcc-release-default/* $QLK_INSTALL/inc/
        elif [ -f $QLK_HOME/lib/libQLK-intel-release-default.a ]
        then
            cp $QLK_HOME/lib/libQLK-intel-release-default.a $QLK_INSTALL/lib/
            cp $QLK_HOME/include/intel-release-default/* $QLK_INSTALL/inc/
        else
            echo "QuaLiKiZ build did not produce expected library/include files, interrupting installation script"
            exit 1
        fi
        cp $AWD/platform/env.$platform $QLK_INSTALL/
        echo $QLK_HASH | cat > $QLK_INSTALL/hash
        echo QuaLiKiZ built in $QLK_HOME installed in $QLK_INSTALL

        export FC=$FC_IN
        export LINK=$LINK_IN
        export LNKGRPBEG=$LNKGRPBEG_IN
        export LNKGRPEND=$LNKGRPEND_IN
        export PATH=$PATH_IN
    fi

    read -p "Install QuaLiKiz NN (y/n) " QLKNN_FLAG

    if [ "$QLKNN_FLAG" = "y" ]; then
        cd /tmp
        rm -rf QLKNN-fortran
        git clone https://gitlab.com/qualikiz-group/QLKNN-fortran.git
        QLKNN_HOME=/tmp/QLKNN-fortran
        cd $QLKNN_HOME
        QLKNN_HASH=`git rev-parse HEAD`
        git submodule init
        git submodule update
        if [ -f fruit.make ]; then
            perl -0pi -e 's/\$\(call LOCAL_mod_dep, tests\/test_d01ahf\.f90, fruit_mpi\.mod\)/\$(call LOCAL_mod_dep, tests\/test_d01ahf.f90, fruit.mod)/g; s/\$\(call LOCAL_mod_dep, tests\/test_integration\.f90, fruit_mpi\.mod\)/\$(call LOCAL_mod_dep, tests\/test_integration.f90, fruit.mod)/g' fruit.make
        fi

        FC_IN=$FC
        LINK_IN=$LINK
        LNKGRPBEG_IN=$LNKGRPBEG
        LNKGRPEND_IN=$LNKGRPEND
        export TOOLCHAIN=gcc
        export FC=$FC_SERIAL
        export LINK=$FC_SERIAL
        export QLK_HAVE_NAG=0
        export TUBSCFG_MPI=0
        export VERBOSE=1
        export BUILD=release
        export LNKGRPBEG=
        export LNKGRPEND=
        PATH_IN=$PATH
        GCC_AR=$(ls "$(dirname "$CC")"/gcc-ar-* 2>/dev/null | sort -V | tail -n 1)
        GCC_NM=$(ls "$(dirname "$CC")"/gcc-nm-* 2>/dev/null | sort -V | tail -n 1)
        GCC_RANLIB=$(ls "$(dirname "$CC")"/gcc-ranlib-* 2>/dev/null | sort -V | tail -n 1)
        if [ -x "$GCC_AR" ]; then
            mkdir -p /tmp/astra_gcc_bin
            ln -sf "$GCC_AR" /tmp/astra_gcc_bin/gcc-ar
        fi
        if [ -x "$GCC_NM" ]; then
            mkdir -p /tmp/astra_gcc_bin
            ln -sf "$GCC_NM" /tmp/astra_gcc_bin/gcc-nm
        fi
        if [ -x "$GCC_RANLIB" ]; then
            mkdir -p /tmp/astra_gcc_bin
            ln -sf "$GCC_RANLIB" /tmp/astra_gcc_bin/gcc-ranlib
        fi
        if [ -d /tmp/astra_gcc_bin ]; then
            export PATH=/tmp/astra_gcc_bin:$PATH
        fi

        $GMAKE
        status=$?
        if [ $status -ne 0 ]; then
            echo "Error building QuaLiKiZ-NN, interrupting installation script"
            exit 1
        fi

        rm -rf $QLKNN_INSTALL
        mkdir -p $QLKNN_INSTALL/lib
        mkdir -p $QLKNN_INSTALL/inc
        if [ -f $QLKNN_HOME/lib/libQLKNN-gcc-release-default.a ]
        then
            cp $QLKNN_HOME/lib/libQLKNN-gcc-release-default.a $QLKNN_INSTALL/lib/
            cp $QLKNN_HOME/include/gcc-release-default/* $QLKNN_INSTALL/inc/
        elif [ -f $QLKNN_HOME/lib/libQLKNN-intel-release-default.a ]
        then
            cp $QLKNN_HOME/lib/libQLKNN-intel-release-default.a $QLKNN_INSTALL/lib/
            cp $QLKNN_HOME/include/intel-release-default/* $QLKNN_INSTALL/inc/
        else
            echo "QuaLiKiZ-NN build did not produce expected library/include files, interrupting installation script"
            exit 1
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
        if [ ! -d qlknn-hyper-namelists ] || [ ! -d qlknn-jetexp-namelists ] || [ ! -d qlknn-hornnet-namelists ] || [ ! -d qlknn-fullflux-namelists ]; then
            echo "QuaLiKiZ-NN install did not produce expected namelist directories, interrupting installation script"
            exit 1
        fi
        echo $QLKNN_HASH | cat > $QLKNN_INSTALL/hash
        echo QuaLiKiZ-NN built in $QLKNN_HOME installed in $QLKNN_INSTALL

        export FC=$FC_IN
        export LINK=$LINK_IN
        export LNKGRPBEG=$LNKGRPBEG_IN
        export LNKGRPEND=$LNKGRPEND_IN
        export PATH=$PATH_IN
    fi

    read -p "Install TGLF/NEO (y/n) " GA_FLAG

    if [ "$GA_FLAG" = "y" ]; then
        brew install fftw
        cd /tmp
        rm -rf gacode
        git clone git@github.com:gafusion/gacode.git
        export GACODE_ROOT=/tmp/gacode
        export GACODE_PLATFORM=GFORTRAN_OSX_BREW
        . $GACODE_ROOT/shared/bin/gacode_setup
        cat <<EOT >$GACODE_ROOT/platform/env/env.GFORTRAN_OSX_BREW
#!/bin/sh

export BREW_LIB=$(brew --prefix)/lib
export LD_LIBRARY_PATH=$(brew --prefix)/lib:$LD_LIBRARY_PATH
export FFTW_INC=$(brew --prefix fftw)/include
export NETCDF_INC=$(brew --prefix netcdf-fortran)/include
EOT
        . $GACODE_ROOT/platform/env/env.GFORTRAN_OSX_BREW
        cd $GACODE_ROOT
        GACODE_HASH=$(git rev-parse HEAD)
        cd $GACODE_ROOT/tglf
        make
        status=$?
        if [ $status -ne 0 ]; then
            echo "Error building TGLF, interrupting installation script"
            exit 1
        fi
        cd $GACODE_ROOT/neo
        make
        status=$?
        if [ $status -ne 0 ]; then
            echo "Error building NEO, interrupting installation script"
            exit 1
        fi
        if [ ! -f $GACODE_ROOT/tglf/src/tglf_lib.a ] || [ ! -f $GACODE_ROOT/neo/src/neo_lib.a ] || [ ! -f $GACODE_ROOT/modules/neo_interface.mod ]; then
            echo "TGLF/NEO build did not produce expected library/module files, interrupting installation script"
            exit 1
        fi

        rm -rf $TGLF_INSTALL
        mkdir -p $TGLF_INSTALL/lib
        mkdir -p $TGLF_INSTALL/inc
        rm -rf $NEO_INSTALL
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
        echo GACODE installed in $GACODE_ROOT
    fi

    read -p "Install STRAHL (y/n) " STRAHL_FLAG

    if [ "$STRAHL_FLAG" = "y" ]; then
        cd /tmp
        rm -rf strahl
        git clone git@gitlab.mpcdf.mpg.de:azito/strahl.git strahl
        STRAHL_HOME=/tmp/strahl
        cd $STRAHL_HOME/source/compile
        ./make_all my_setup
        status=$?
        if [ $status -ne 0 ]; then
            echo "Error building STRAHL, interrupting installation script"
            exit 1
        fi
        if [ ! -f $STRAHL_HOME/bin/strahl ] || [ ! -f $STRAHL_HOME/bin/result_to_astra ]; then
            echo "STRAHL build did not produce expected binaries, interrupting installation script"
            exit 1
        fi

        rm -rf $STRAHL_INSTALL
        mkdir -p $STRAHL_INSTALL/bin
        cp $STRAHL_HOME/bin/strahl $STRAHL_INSTALL/bin/
        cp $STRAHL_HOME/bin/result_to_astra $STRAHL_INSTALL/bin/
        cp $AWD/platform/env.$platform $STRAHL_INSTALL/
        echo Strahl installed in $STRAHL_HOME
    fi

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
