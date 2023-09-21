#!/bin/bash -f

rootdir=`dirname $0`       # may be relative path
AWD=`cd $rootdir && pwd`  # ensure absolute path

platform=`$AWD/get_platform`

cd $AWD
for DIR in $(ls repo_user_area)
do
    cp -r repo_user_area/$DIR .
done

cp $AWD/repo_user_area/exp/nml/aug34954_${platform} $AWD/exp/nml/aug34954
chmod 744 $AWD/get_platform
chmod 744 $AWD/exe/*
chmod 744 $AWD/pyparse/*.py

$AWD/exe/as_exe -m fluxes -v aug34954 -s 4 -e 5
