#!/bin/bash

rootdir=`dirname $0`        # may be relative path
AWD=`cd $rootdir && pwd` # ensure absolute path

platform=`$AWD/get_platform`
source $AWD/platform/env.$platform
if [[ ":$PATH:" != *":$PYTHON_BIN:"* ]]
then
    export PATH=${PYTHON_BIN}:${PATH}
fi

EXP=aug34954
for EQU in fluxes feqis qlk tglf qlknn
do
    $AWD/exe/as_exe -m $EQU -v $EXP -s 4. -e 5.
    python3 $AWD/compareRegressions.py -m $EQU -v $EXP
done

EXP=aug34954_t
for EQU in fluxes
do
    $AWD/exe/as_exe -m $EQU -v $EXP -s 4. -e 5.
    python3 $AWD/compareRegressions.py -m $EQU -v $EXP
done

EQU=fbe
EXP=AUG33040_2500
$AWD/exe/as_exe -m $EQU -v $EXP -s 2.48 -e 2.7
python3 $AWD/compareRegressions.py -m $EQU -v $EXP

EQU=tglf_pid
EXP=AUG36982_3400
$AWD/exe/as_exe -m $EQU -v $EXP -s 4. -e 6.
python3 $AWD/compareRegressions.py -m $EQU -v $EXP

# Slow ones

EXP=aug34954
# for EQU in flux_neo_tglf flux_tlf_serial
for EQU in flux_neo_tglf
do
    $AWD/exe/as_exe -m $EQU -v $EXP -s 4. -e 5.
    python3 $AWD/compareRegressions.py -m $EQU -v $EXP
done
