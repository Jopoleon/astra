#!/bin/bash

rootdir=`dirname $0`      # may be relative path
export AWD=`cd $rootdir && pwd`  # ensure absolute path

chmod 744 $AWD/get_platform
platform=`$AWD/get_platform`
SAFE=$1

cd $AWD

# Prompt overwriting option
for DIR in $(ls repo_user_area)
do
    if [ ! -d "$DIR" ] || [ "$SAFE" != "-safe" ] || [ "$DIR" = "equ" ] || [ "$DIR" = "exp" ] || [ "$DIR" = "pyparse" ] || [ "$DIR" = "strahl" ] || [ "$DIR" = "tmp" ] || [ "$DIR" = "udb" ]
    then
        cp -r repo_user_area/$DIR .
    else
        read -p "'$DIR': backup + copy from repo_user_area_dir? y/n " OVERWRITE
        if [ "$OVERWRITE" = "y" ]
	then
            rm -rf ${DIR}_backup
            cp -r $DIR ${DIR}_backup
            cp -r repo_user_area/$DIR .
            echo "Backup in $AWD/${DIR}_backup"
            printf "Copy $AWD/repo_use_area/$DIR to $AWD/$DIR\n\n"
        fi
    fi
done

if [ "$SAFE" = "-safe" ]
then
    read -p "Do you have a working installation for RABBIT ? y/n " RABBIT
    read -p "Do you have a working installation for TORBEAM? y/n " TORBEAM
    read -p "Do you want QualiKiZ? y/n " QLK
    read -p "Do you want NEO? y/n " NEO
# RABBIT, TORBEAM
    if [ "$RABBIT" = "n" ]
    then
	rm $AWD/sbr/rabbit.f90
	rm $AWD/sbr/a2rabbit.f90
	rm $AWD/sbr/torbeam_rabbit.f90
	sed -i "s#export\ RABBIT_LIB#\#export\ RABBIT_LIB#g" $AWD/exe/astra_rc
	if [ "$TORBEAM" = "n" ]
	then
            rm $AWD/sbr/torba.f90
            rm $AWD/sbr/a2torbeam.f90
            sed -i "s#TORBEAM_RABBIT#\!TORBEAM_RABBIT#g" $AWD/equ/fluxes
            sed -i "s#export\ TORB_LIB#\#export\ TORB_LIB#g" $AWD/exe/astra_rc
	else
            sed -i "s#TORBEAM_RABBIT#TORBA#g" $AWD/equ/fluxes
	fi
    else
	if [ "$TORBEAM" = "n" ]
	then
            rm $AWD/sbr/torba.f90
            rm $AWD/sbr/a2torbeam.f90
            rm $AWD/sbr/torbeam_rabbit.f90
            sed -i "s#TORBEAM_RABBIT#RABBIT#g" $AWD/equ/fluxes
	fi
    fi
# NEO, QUALIKIZ (needing MPI)
    if [ "$QLK" = "n" ]
    then
	if [ "$NEO" = "n" ]
	then
	    sed -i "s#\$(XPR)\/qlki\ \$(XPR)\/neo#\ #g" $AWD/exe/Makexpr
	else
	    sed -i "s#\$(XPR)\/qlki\ \$(XPR)\/neo#\$(XPR)\/neo#g" $AWD/exe/Makexpr
	fi
    else
	if [ "$NEO" = "n" ]
	then
	    sed -i "s#\$(XPR)\/qlki\ \$(XPR)\/neo#\$(XPR)\/qlki#g" $AWD/exe/Makexpr
	fi
    fi    
fi

source $AWD/platform/env.$platform                     # get platform dependent $ASTRA_EXT
echo $ASTRA_EXT

chmod 744 $AWD/exe/Build
chmod 744 $AWD/exe/as_exe
chmod 744 $AWD/exe/wr_nml
chmod 744 $AWD/exe/CheckObjs
chmod 744 $AWD/pyparse/parser_main.py
chmod 744 $AWD/clean.sh
chmod 744 $AWD/compReg.sh

if [ "$SAFE" = "-safe" ]
then
    read -p "Clean dirs for full Make? y/n " CLEAN
    if [ "$CLEAN" = "y" ]
    then
	make -f exe/Makefile clean
    fi
else
    make -f exe/Makefile clean
fi
exe/as_exe -m fluxes -v aug34954 -s 4 -e 5
