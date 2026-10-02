#!/bin/bash
# ASTRA 8 -> STRAHL (ASTRA 7 package) wrapper, installed as <ASTRA_EXT>/strahl/sep23/bin/strahl.
# ASTRA 8 runs "strahl a q" in the instance directory; convert the parameter file first.
here=$(cd "$(dirname "$0")" && pwd)
[ -f param_files/sparams.dat ] && python3 "$here/sparams_a8_to_a7.py" param_files/sparams.dat
exec "$here/strahl_a7" "$@"
