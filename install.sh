#!/bin/bash -f

source /etc/profile.d/modules.sh

rootdir=`dirname $0`       # may be relative path
AWD=`cd $rootdir && pwd`  # ensure absolute path

platform=`$AWD/get_platform`

cd $AWD
for DIR in $(ls repo_user_area)
do
    cp -r repo_user_area/$DIR .
done

cp $AWD/repo_user_area/exe/platforms/astra_rc_${platform} $AWD/exe/astra_rc
cp $AWD/repo_user_area/exp/nml/aug34954_${platform}       $AWD/exp/nml/aug34954
cp $AWD/repo_user_area/tmp/astra_${platform}.nml          $AWD/tmp/astra.nml

$AWD/exe/as_exe -m fluxes -v aug34954 -s 4 -e 5
