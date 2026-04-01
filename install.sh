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
    if [ ! -d "$DIR" ] || [ "$SAFE" != "-safe" ] || [ "$DIR" = "equ" ] || [ "$DIR" = "exp" ] || [ "$DIR" = "udb" ]
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
    read -p "Do you want QualiKiZ NN? y/n " QLKNN
    read -p "Do you want NEO? y/n " NEO
    read -p "Do you want TGLF? y/n " TGLF
# RABBIT, TORBEAM
    if [ "$RABBIT" = "n" ]
    then
	rm $AWD/sbr/rabbit.f90
	cat <<EOF > $AWD/sbr/rabbit.f90
module a2rabbit
contains
subroutine rabbit
end subroutine rabbit
end module a2rabbit
EOF
	sed -i "s#export\ RABBIT_LIB#\#export\ RABBIT_LIB#g" $AWD/exe/astra_rc
	sed -i "s#RABBIT#\!RABBIT#g" $AWD/equ/flux_spider
	sed -i "s#RABBIT#\!RABBIT#g" $AWD/equ/flux_feqis
    fi

    if [ "$TORBEAM" = "n" ]
    then
        rm $AWD/sbr/torba.f90
	cat <<EOF > $AWD/sbr/torbeam.f90
module a2torbeam
contains
subroutine torba
end subroutine torba
end module a2torbeam
EOF
	sed -i "s#export\ TORB_LIB#\#export\ TORB_LIB#g" $AWD/exe/astra_rc
	sed -i "s#TORBA#\!TORBA#g" $AWD/equ/flux_spider
	sed -i "s#TORBA#\!TORBA#g" $AWD/equ/flux_feqis
    fi

    if [ "$QLKNN" = "n" ]
    then
        rm $AWD/sbr/qlknn_serial.f90
	sed -i "s#export\ QLKNN_LIB#\#export\ QLKNN_LIB#g" $AWD/exe/astra_rc
    fi

# NEO, QUALIKIZ, TGLF (needing MPI)
    if [ "$QLK" = "n" ]
    then
	sed -i "s#export\ QLK_LIB#\#export\ QLK_LIB#g" $AWD/exe/astra_rc
	sed -i -e '/all: directories/ s/\$(XPR)\/qlki//g' $AWD/exe/Makexpr
    fi
    if [ "$NEO" = "n" ]
    then
	sed -i "s#export\ NEO_LIB#\#export\ NEO_LIB#g" $AWD/exe/astra_rc
	sed -i -e '/all: directories/ s/\$(XPR)\/neo//g' $AWD/exe/Makexpr
    fi
    if [ "$TGLF" = "n" ]
    then
	sed -i "s#export\ TGLF_LIB#\#export\ TGLF_LIB#g" $AWD/exe/astra_rc
	sed -i -e '/all: directories/ s/\$(XPR)\/tglfi//g' $AWD/exe/Makexpr
	rm $AWD/sbr/tglf_serial.f90
    fi
fi

# get platform dependent $ASTRA_EXT
source $AWD/platform/env.$platform
echo $ASTRA_EXT

mkdir -p $AWD/tmp

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
exe/as_exe -m flux_feqis -v aug34954 -s 4 -e 5
