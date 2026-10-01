#!/usr/bin/env bash
# Container entrypoint: optional host data folder + astra8.
#
# If a folder is mounted at /data, ASTRA's user folders live there, so models, experiment
# files, U-files and results survive the container:
#   /data/equ  /data/exp  /data/udb  /data/ncdf_out
# Missing files are filled in from the image (existing host files are never overwritten).
#
# Arguments go to astra8 (exe/as_exe). Special first arguments:
#   shell     interactive bash inside the container
#   <cmd>     any other command that is not an option is executed as is

set -euo pipefail

if cp --help | grep -q -- '--update.*none'; then
    cp_keep=(cp -r --update=none)
else
    cp_keep=(cp -rn)
fi

if [[ -d /data ]]; then
    for d in equ exp udb ncdf_out; do
        mkdir -p "/data/$d"
        if [[ -d "$AWD/$d" && ! -L "$AWD/$d" ]]; then
            "${cp_keep[@]}" "$AWD/$d/." "/data/$d/"
            rm -rf "${AWD:?}/$d"
        fi
        ln -sfn "/data/$d" "$AWD/$d"
    done
fi

if [[ $# -gt 0 && "$1" == "shell" ]]; then
    exec bash -l
elif [[ $# -gt 0 && "$1" != -* ]]; then
    exec "$@"
fi
exec astra8 "$@"
