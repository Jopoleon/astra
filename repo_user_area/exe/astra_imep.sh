#!/bin/bash -f

source /etc/profile.d/modules.sh

rootdir=`dirname $0`        # may be relative path
AWD=`cd $rootdir/.. && pwd` # ensure absolute path

listEqu="imep_pw04 imep_pw08"
exp="30000_3.4"
COUNTER=1
while [[ $COUNTER -ne 42 ]]; do
    for equ in $listEqu
    do
        rm -f "$AWD/ncdf_out/$exp$equ-$COUNTER.cdf"
    done
    COUNTER=$(($COUNTER+1))
done

cd $AWD

for equ in $listEqu
do
    exe/as_exe -m $equ -v $exp -s 0 -e 4 -b
done

module load astra

for equ in $listEqu
do
    fcdf="$AWD/ncdf_out/$exp$equ-41.cdf"
    while [ ! -f $fcdf ]; do
        echo $fcdf does not exist yet, waiting 10s
        sleep 10
    done
    sleep 5
    exe/nc_concat.py -e $exp$equ
    exe/astra2helena.py -m $equ -v $exp
done
