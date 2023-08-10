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

cp $AWD/repo_user_area/exp/nml/aug34954_${platform} $AWD/exp/nml/aug34954
cp $AWD/repo_user_area/tmp/astra.nml $AWD/tmp/astra.nml

if [[ $platform == "ldaug" ]]
then
    module use /shares/departments/AUG/users/git/modulefiles
elif [[ $platform == "iter" ]]
then
    module use /home/ITER/tarding/modulefiles
elif [[ $platform == "cz" ]]
then
    module use /compass/home/tardini/modulefiles
elif [[ $platform == "omega" ]]
then
    module use /home/tardinig/modulefiles
elif [[ $platform == "mit" ]]
then
    module use /home/gtardini/modulefiles
fi

if ! command -v module &> /dev/null
then
    echo "module tool not available on this platform, set environment in exe/Build"
else
    module load astra
    RES=$( { module load astra; } 2>&1 )
    echo $RES
    if [[ "$RES" == *"ERROR"* ]];
    then
        echo module astra not loaded, proceed anyway?
        read -p "y/n:" proc
        if [[ $proc == "n" ]]
        then
            exit 1
        fi
    fi
fi

$AWD/exe/as_exe -m fluxes -v aug34954 -s 4 -e 5
