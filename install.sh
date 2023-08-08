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

cp $AWD/repo_user_area/exp/nml/aug34954_${platform}       $AWD/exp/nml/aug34954
cp $AWD/repo_user_area/tmp/astra_${platform}.nml          $AWD/tmp/astra.nml

if [[ $platform == "lx" ]]
then
    module use /afs/.ipp/common/usr/modules/@sys/modulefiles/AUG
elif [[ $platform == "ldaug" ]]
     module use /shares/departments/AUG/users/git/modulefiles
fi
module load astra
status=$?
if [ $status -ne 0 ]
then
    echo module astra was not loaded
    exit 1
fi
$AWD/exe/as_exe -m fluxes -v aug34954 -s 4 -e 5
