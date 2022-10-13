#!/bin/bash -f

source /etc/profile.d/modules.sh

rootdir=`dirname $0`          # may be relative path
AREP=`cd $rootdir && pwd`  # ensure absolute path

platform=`$AREP/get_platform`
echo Server installation for platform $platform

module purge

if [[ $platform == "gway" ]]
then
    module load cineca intel
elif [[ $platform == "iter" ]]
then
    module load intel
elif [[ $platform == "cz" ]]
then
    module load intelstudio/18
elif [[ $platform == "ga" ]]
then
    module load intel/18
else
    module load intel
fi

mkdir -p ${AREP}/inc

make -f Makefor
make -f Makenbi
make -f Makeequil
