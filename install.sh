#!/bin/bash -f

source /etc/profile.d/modules.sh

rootdir=`dirname $0`       # may be relative path
AWD=`cd $rootdir && pwd`  # ensure absolute path

platform=`$AWD/get_platform`

cd $AWD
for DIR in $(ls astra_user)
do
    cp -r astra_user/$DIR .
done

cp $AWD/astra_user/exe/platforms/astra_rc_${platform} $AWD/exe/astra_rc
cp $AWD/astra_user/exp/nml/aug34954_${platform}       $AWD/exp/nml/aug34954
cp $AWD/astra_user/tmp/astra_${platform}.nml          $AWD/tmp/astra.nml

$AWD/exe/as_exe -m fluxes -v aug34954 -s 4 -e 5
