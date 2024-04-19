#!/bin/bash -f

rootdir=`dirname $0`       # may be relative path
AWD=`cd $rootdir && pwd`  # ensure absolute path

chmod 744 $AWD/get_platform
platform=`$AWD/get_platform`

cd $AWD
for DIR in $(ls repo_user_area)
do
    cp -r repo_user_area/$DIR .
done

cp $AWD/repo_user_area/exp/nml/aug34954_${platform} $AWD/exp/nml/aug34954
cp $AWD/repo_user_area/exp/nml/AUG33040_2500_${platform} $AWD/exp/nml/AUG33040_2500

chmod 744 $AWD/exe/Build
chmod 744 $AWD/exe/as_exe
chmod 744 $AWD/exe/nc_concat.py
chmod 744 $AWD/exe/wr_nml
chmod 744 $AWD/exe/CheckObjs
chmod 744 $AWD/pyparse/parser_main.py
chmod 744 $AWD/green/greenMatrices.py
chmod 744 $AWD/clean.sh

make -f exe/Makefile clean
exe/as_exe -m fluxes -v aug34954 -s 4 -e 5

# Generate Green-functions for FBE
green/greenMatrices.py
