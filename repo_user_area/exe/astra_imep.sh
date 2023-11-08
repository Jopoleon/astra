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

exe/as_exe -m imep  -v 30000_3.4 -s 4 -e 5 -b
exe/as_exe -m imep2 -v 30000_3.4 -s 4 -e 5 -b

file1="$AWD/ncdf_out/30000_3.4imep-51.cdf"
file2="$AWD/ncdf_out/30000_3.4imep2-51.cdf"

module load astra

while [ ! -f $file1 ]; do
    echo $file1 does not exist yet, waiting 10s
    sleep 10
done
sleep 5
exe/nc_concat.py -e 30000_3.4imep

while [ ! -f $file2 ]; do
    echo $file2 does not exist yet, waiting 10s
    sleep 10
done
sleep 5
exe/nc_concat.py -e 30000_3.4imep2
