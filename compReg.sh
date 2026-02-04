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
for EQU in fluxes feqis flux_tglf qlk tglf
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

EXP=30000_3.4
for EQU in imep_pw04 imep_pw08
do
    $AWD/exe/as_exe -m $EQU -v $EXP -s 4. -e 5.
    python3 $AWD/compareRegressions.py -m $EQU -v $EXP
done

# Slow ones

EXP=aug34954
for EQU in flux_tglf_serial flux_neo
do
    $AWD/exe/as_exe -m $EQU -v $EXP -s 4. -e 5.
    python3 $AWD/compareRegressions.py -m $EQU -v $EXP
done
