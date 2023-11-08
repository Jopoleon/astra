#!/bin/bash -f

source /etc/profile.d/modules.sh

AWD="/toks/work/git/a82"
COUNTER=1
while [[ $COUNTER -ne 52 ]]; do
    rm "$AWD/ncdf_out/30000_3.4imep-$COUNTER.cdf"
    rm "$AWD/ncdf_out/30000_3.4imep2-$COUNTER.cdf"
    COUNTER=$(($COUNTER+1))
done

cd $AWD

exp="30000_3.4"

for equ in imep imep2
do
    exe/as_exe -m $equ -v $exp -s 4 -e 5 -b
done

module load astra

for equ in imep imep2
do
    fcdf="$AWD/ncdf_out/30000_3.4$equ-51.cdf"
    while [ ! -f $fcdf ]; do
        echo $fcdf does not exist yet, waiting 10s
        sleep 10
    done
    sleep 5
    exe/nc_concat.py -e $exp$equ
    exe/astra2helena.py -m $equ -v $exp
done
