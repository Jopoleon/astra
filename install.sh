#!/bin/bash -f

rootdir=`dirname $0`      # may be relative path
AWD=`cd $rootdir && pwd`  # ensure absolute path

chmod 744 $AWD/get_platform
platform=`$AWD/get_platform`
SAFE=$1

cd $AWD

# Create backup
for DIR in exe fml fnc sbr xpr
do
    if [ ! -d "$DIR" ]
    then
        printf "Subdir '$DIR' exists, creating backup in $AWD/${DIR}_backup\n"
        rm -rf ${DIR}_backup
        cp -r $DIR ${DIR}_backup
    fi
done

# Prompt overwriting option
for DIR in $(ls repo_user_area)
do
    if [ ! -d "$DIR" ] || [ "$SAFE" != "-safe" ] || [ "$DIR" = "equ" ] || [ "$DIR" = "exp" ] || [ "$DIR" = "pyparse" ] || [ "$DIR" = "strahl" ] || [ "$DIR" = "tmp" ] || [ "$DIR" = "udb" ]
    then
        cp -r repo_user_area/$DIR .
    else
        read -p "Overwrite safely '$DIR' from repo_user_area_dir? y/n " OVERWRITE
        if [ "$OVERWRITE" = "y" ]
	then
            cp -r repo_user_area/$DIR .
        fi
        printf "\n"
    fi
done

cp $AWD/repo_user_area/exp/nml/aug34954_${platform} $AWD/exp/nml/aug34954
cp $AWD/repo_user_area/exp/nml/AUG33040_2500_${platform} $AWD/exp/nml/AUG33040_2500

chmod 744 $AWD/exe/Build
chmod 744 $AWD/exe/as_exe
chmod 744 $AWD/exe/json2cdf.py
chmod 744 $AWD/exe/wr_nml
chmod 744 $AWD/exe/CheckObjs
chmod 744 $AWD/pyparse/parser_main.py
chmod 744 $AWD/green/greenMatrices.py
chmod 744 $AWD/clean.sh

make -f exe/Makefile clean
exe/as_exe -m fluxes -v aug34954 -s 4 -e 5
