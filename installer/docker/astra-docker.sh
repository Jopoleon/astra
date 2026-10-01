#!/usr/bin/env bash
# Run ASTRA 8 in Docker with the graphic window and a data folder on the host.
#
#   docker/astra-docker.sh -m flux_feqis -v aug34954 -s 4 -e 5          GUI window
#   docker/astra-docker.sh -m flux_feqis -v aug34954 -s 4 -e 5 -batch   no window
#   docker/astra-docker.sh shell                                        bash inside
#
# Environment:
#   ASTRA_IMAGE  image name (default astra8; built automatically if missing)
#   ASTRA_DATA   host folder with equ/ exp/ udb/ ncdf_out/ (default ~/astra-data)

set -euo pipefail

here="$(cd "$(dirname "$0")" && pwd)"
image="${ASTRA_IMAGE:-astra8}"
data="${ASTRA_DATA:-$HOME/astra-data}"

command -v docker >/dev/null || { echo "docker not found" >&2; exit 1; }

if ! docker image inspect "$image" >/dev/null 2>&1; then
    echo "Image '$image' not found, building it (5-15 min, once)..."
    # same UID/GID as the host user, so files in the data folder stay yours
    docker build -f "$here/Dockerfile" -t "$image" \
        --build-arg USER_UID="$(id -u)" --build-arg USER_GID="$(id -g)" "$here/.."
fi

mkdir -p "$data"
args=(--rm -v "$data:/data")
[[ -t 0 && -t 1 ]] && args+=(-it)

# X11 window: WSLg (Windows 11) or a Linux desktop both expose /tmp/.X11-unix
if [[ -n "${DISPLAY:-}" && -d /tmp/.X11-unix ]]; then
    args+=(-e "DISPLAY=$DISPLAY" -v /tmp/.X11-unix:/tmp/.X11-unix)
    # Linux desktop (not WSLg): allow local containers to use the X server
    if [[ ! -d /mnt/wslg ]] && command -v xhost >/dev/null; then
        xhost +local: >/dev/null 2>&1 || true
    fi
fi

exec docker run "${args[@]}" "$image" "$@"
