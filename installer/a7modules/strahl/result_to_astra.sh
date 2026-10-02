#!/bin/bash
here=$(cd "$(dirname "$0")" && pwd)
exec python3 "$here/result_to_astra.py" "$@"
