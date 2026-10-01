#!/usr/bin/env bash
# Installer for ASTRA 8 (https://gitlab.mpcdf.mpg.de/git/astra) on a plain Ubuntu PC.
#
# What it does (each step is idempotent, re-running is safe):
#   1. apt: gfortran, cmake, X11 headers, BLAS/LAPACK, python3 + numpy/scipy/netCDF4
#   2. clones ASTRA (HTTPS, LF line endings) into $PREFIX/a8
#   3. builds json-fortran into $PREFIX/ext (the only mandatory external library)
#   4. adds an "ubuntu" platform (gfortran + apt libraries) and a few portability fixes
#   5. prepares the user area like upstream install.sh, without closed IPP modules
#   6. builds and runs the reference case flux_feqis / aug34954 in batch mode
#   7. installs the `astra8` launcher into ~/.local/bin
#
# Usage: ./install_astra8.sh [options]     (see --help)

set -Eeuo pipefail

# ---------------------------------------------------------------- defaults

PREFIX="${HOME}/astra"
ASTRA_URL="https://gitlab.mpcdf.mpg.de/git/astra.git"
# Last upstream commit this installer was verified against; --ref main for the newest code
ASTRA_REF="e66f6ba13dd0160ed5b3dd37a46d0671cb81e85d"
FULL_HISTORY=0             # 1: full git history (~1 GB) instead of one commit (~90 MB)
BLAS="openblas"             # openblas | mkl (results agree to ~1e-6)
MKL_VERSION="2026.1.0"      # Intel mkl-static wheel from PyPI
JSON_VERSION="9.0.2"
DO_APT=1
DO_TEST=1
ASSUME_YES=0

# ---------------------------------------------------------------- helpers

if [[ -t 1 ]]; then
    C_OK=$'\e[32m'; C_WARN=$'\e[33m'; C_ERR=$'\e[31m'; C_STEP=$'\e[1;36m'; C_OFF=$'\e[0m'
else
    C_OK=""; C_WARN=""; C_ERR=""; C_STEP=""; C_OFF=""
fi
step() { printf '\n%s==> %s%s\n' "$C_STEP" "$*" "$C_OFF"; }
ok()   { printf '%s  ok:%s %s\n' "$C_OK" "$C_OFF" "$*"; }
warn() { printf '%s  warning:%s %s\n' "$C_WARN" "$C_OFF" "$*"; }
die()  { printf '%serror:%s %s\n' "$C_ERR" "$C_OFF" "$*" >&2; exit 1; }
trap 'die "failed at line $LINENO: $BASH_COMMAND (full log: $LOG)"' ERR

usage() {
    cat <<EOF
ASTRA 8 installer for Ubuntu (gfortran build, no closed IPP modules)

Usage: $0 [options]

  --prefix DIR      install root (default: \$HOME/astra)
                      DIR/a8   ASTRA working copy (models, data, results)
                      DIR/ext  external libraries (json-fortran)
  --ref REF         ASTRA branch, tag or full commit hash (default: verified
                    ${ASTRA_REF:0:8}; use "main" for the newest upstream code)
  --full-history    clone the whole ASTRA git history (~1 GB) instead of
                    only the requested commit (~90 MB)
  --blas openblas|mkl
                    LAPACK library (default: openblas from apt; mkl downloads
                    Intel's static MKL $MKL_VERSION from PyPI, ~250 MB)
  --no-apt          do not install system packages (they must already be present)
  --no-test         skip the reference test run
  -y, --yes         do not ask for confirmation
  -h, --help        this help
EOF
}

# ---------------------------------------------------------------- arguments

while [[ $# -gt 0 ]]; do
    case "$1" in
        --prefix)  PREFIX="$2"; shift 2 ;;
        --ref)     ASTRA_REF="$2"; shift 2 ;;
        --full-history) FULL_HISTORY=1; shift ;;
        --blas)    BLAS="$2"; shift 2 ;;
        --no-apt)  DO_APT=0; shift ;;
        --no-test) DO_TEST=0; shift ;;
        -y|--yes)  ASSUME_YES=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) usage; die "unknown option: $1" ;;
    esac
done
[[ "$BLAS" == "mkl" || "$BLAS" == "openblas" ]] || die "--blas must be openblas or mkl"

PREFIX="$(mkdir -p "$PREFIX" && cd "$PREFIX" && pwd)"
AWD="$PREFIX/a8"
EXT="$PREFIX/ext"
BUILD="$PREFIX/build"
LOG="$PREFIX/install.log"
: > "$LOG"

# Long, noisy commands go to the log; the screen shows only the steps.
run() { "$@" >>"$LOG" 2>&1; }

# ---------------------------------------------------------------- 0. checks

step "Checking the system"
[[ "$(uname -s)" == "Linux" ]] || die "Linux only"
[[ "$(uname -m)" == "x86_64" ]] || warn "tested only on x86_64, found $(uname -m)"
if [[ -r /etc/os-release ]]; then
    . /etc/os-release
    ok "$PRETTY_NAME"
    [[ "${ID:-}" == "ubuntu" || "${ID_LIKE:-}" == *debian* ]] || warn "not Ubuntu/Debian: apt step may fail"
fi
if [[ $EUID -eq 0 ]]; then
    warn "running as root: ASTRA will be installed for root; normally run as your own user"
    SUDO=""
else
    SUDO="sudo"
fi

cat <<EOF

  ASTRA source : $ASTRA_URL @ $ASTRA_REF
  Install root : $PREFIX
  BLAS/LAPACK  : $BLAS
  Install log  : $LOG

EOF
if [[ $ASSUME_YES -eq 0 ]]; then
    read -r -p "Continue? [Y/n] " answer
    [[ -z "$answer" || "$answer" =~ ^[YyДд] ]] || exit 0
fi

# ---------------------------------------------------------------- 1. apt

if [[ $DO_APT -eq 1 ]]; then
    step "Installing system packages (sudo password may be asked)"
    packages=(
        git curl ca-certificates make cmake pkg-config
        gcc g++ gfortran
        libx11-dev
        python3 python3-numpy python3-scipy python3-netcdf4 python3-matplotlib
        libnetcdf-dev libnetcdff-dev
        libopenblas-dev liblapack-dev
    )
    run $SUDO apt-get update
    run $SUDO env DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends "${packages[@]}"
    ok "packages installed"
fi

for tool in git curl make cmake gcc gfortran python3; do
    command -v "$tool" >/dev/null || die "'$tool' not found (install it or run without --no-apt)"
done
python3 -c "import numpy, scipy" 2>/dev/null || die "python3 numpy/scipy missing"
ok "$(gfortran --version | head -1)"

# ---------------------------------------------------------------- 2. ASTRA source

step "Getting ASTRA source"
# core.autocrlf=false: with autocrlf=true (common on WSL/Windows users) all scripts
# get CRLF endings and fail with "required file not found"
if [[ -d "$AWD/.git" ]]; then
    ok "existing clone in $AWD, fetching updates"
else
    run git -c core.autocrlf=false init "$AWD"
    run git -C "$AWD" remote add origin "$ASTRA_URL"
fi
git -C "$AWD" config core.autocrlf false
if [[ $FULL_HISTORY -eq 1 ]]; then
    git -C "$AWD" rev-parse -q --is-shallow-repository | grep -q true \
        && run git -C "$AWD" fetch --unshallow --tags origin \
        || run git -C "$AWD" fetch --tags origin
    if git -C "$AWD" rev-parse -q --verify "origin/$ASTRA_REF" >/dev/null; then
        target="origin/$ASTRA_REF"
    else
        target="$ASTRA_REF"
    fi
else
    # Only the requested commit: ~90 MB instead of ~1 GB of history.
    # A commit must be given as the full 40-character hash here.
    run git -C "$AWD" fetch --depth 1 origin "$ASTRA_REF" \
        || die "cannot fetch '$ASTRA_REF' (a commit needs the full hash, or use --full-history)"
    target=FETCH_HEAD
fi
# Our own patches are re-applied below, so reset tracked files to the requested ref
run git -C "$AWD" checkout --force --detach "$target"
ok "ASTRA at $(git -C "$AWD" log -1 --format='%h %cd' --date=short)"

# ---------------------------------------------------------------- 3. json-fortran

JSON_HOME="$EXT/json/$JSON_VERSION"
step "json-fortran $JSON_VERSION"
if [[ -f "$JSON_HOME/lib/libjsonfortran.a" ]]; then
    ok "already built in $JSON_HOME"
else
    mkdir -p "$BUILD"
    run curl -fsSL -o "$BUILD/json-fortran.tar.gz" \
        "https://github.com/jacobwilliams/json-fortran/archive/refs/tags/$JSON_VERSION.tar.gz"
    rm -rf "$BUILD/json-fortran-$JSON_VERSION"
    tar -xzf "$BUILD/json-fortran.tar.gz" -C "$BUILD"
    run cmake -S "$BUILD/json-fortran-$JSON_VERSION" -B "$BUILD/json-fortran-$JSON_VERSION/build" \
        -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_Fortran_FLAGS=-fPIC \
        -DSKIP_DOC_GEN=TRUE -DENABLE_TESTS=FALSE
    run cmake --build "$BUILD/json-fortran-$JSON_VERSION/build" -j "$(nproc)"
    mkdir -p "$JSON_HOME/lib" "$JSON_HOME/inc"
    cp "$BUILD/json-fortran-$JSON_VERSION/build/lib/libjsonfortran.a" "$JSON_HOME/lib/"
    cp "$BUILD/json-fortran-$JSON_VERSION/build/include/"*.mod "$JSON_HOME/inc/"
    ok "installed in $JSON_HOME"
fi

# ---------------------------------------------------------------- 3b. MKL (optional)

MKL_LIBDIR="$EXT/mkl/$MKL_VERSION/lib"
if [[ "$BLAS" == "mkl" ]]; then
    step "Intel MKL $MKL_VERSION (static, from PyPI)"
    if [[ -f "$MKL_LIBDIR/libmkl_core.a" ]]; then
        ok "already in $MKL_LIBDIR"
    else
        # apt's libmkl-dev is ~1 GB and Ubuntu mirrors can be very slow; the PyPI wheel
        # is Intel's official build and holds exactly the three static libraries needed
        mkdir -p "$BUILD" "$MKL_LIBDIR"
        wheel_url="$(curl -fsSL "https://pypi.org/pypi/mkl-static/$MKL_VERSION/json" | python3 -c '
import json, sys
for f in json.load(sys.stdin)["urls"]:
    if "manylinux" in f["filename"] and f["filename"].endswith("x86_64.whl"):
        print(f["url"]); break')"
        [[ -n "$wheel_url" ]] || die "mkl-static $MKL_VERSION wheel not found on PyPI"
        run curl -fSL -o "$BUILD/mkl_static.whl" "$wheel_url"
        python3 - "$BUILD/mkl_static.whl" "$MKL_LIBDIR" <<'PY'
import os, shutil, sys, zipfile
wheel, dest = sys.argv[1:]
with zipfile.ZipFile(wheel) as z:
    for name in z.namelist():
        if os.path.basename(name) in ("libmkl_intel_lp64.a", "libmkl_core.a", "libmkl_sequential.a"):
            with z.open(name) as src, open(os.path.join(dest, os.path.basename(name)), "wb") as dst:
                shutil.copyfileobj(src, dst)
PY
        rm -f "$BUILD/mkl_static.whl"
        [[ -f "$MKL_LIBDIR/libmkl_core.a" ]] || die "MKL libraries not found in the wheel"
        ok "installed in $MKL_LIBDIR"
    fi
fi

# ---------------------------------------------------------------- 4. platform + fixes

step "Configuring the 'ubuntu' platform"

if [[ "$BLAS" == "mkl" ]]; then
    blas_lines="export MKL_LIBDIR=$MKL_LIBDIR"
else
    # The Makefiles link static MKL unconditionally (OpenBLAS only on macOS).
    # Variables in MAKEFLAGS act as command-line overrides for every make call.
    # OPENBLAS_NUM_THREADS=1: ASTRA is OpenMP-parallel itself and upstream links the
    # sequential MKL; threaded OpenBLAS inside OpenMP regions oversubscribes the CPU.
    blas_lines='export MAKEFLAGS="MKL_LIB=-lopenblas"
export OPENBLAS_NUM_THREADS=1'
fi

cat > "$AWD/platform/env.ubuntu" <<EOF
#!/bin/bash
# Generic Ubuntu PC: gfortran and libraries from apt. Written by install_astra8.sh.

export ASTRA_EXT=$EXT
export FC=gfortran
export CC=gcc
export MPIFC=mpif90
export MPICC=mpicc
# -ffree-line-length-none: some sources (e.g. src/astra/a2vmec.f90) exceed 132 columns
export FC_FLAGS="-O0 -w -fPIC -fopenmp -ffree-line-length-none"
# --no-as-needed: Ubuntu links with --as-needed, and the Makefile puts -lX11 before
# libgraph.a, so the GUI build fails with "undefined reference to XOpenDisplay"
export LD_FLAGS="-fopenmp -Wl,--no-as-needed"
$blas_lines
export PYTHON_BINDIR=/usr/bin
EOF
# Unknown hosts are detected as "mypc", for which upstream ships no env file
ln -sf env.ubuntu "$AWD/platform/env.mypc"

# get_platform guesses the site from the hostname (e.g. "node*" -> MIT, "tok*" -> IPP).
# Let ASTRA_PLATFORM override the guess so a lab PC name cannot pick a foreign config.
if ! grep -q ASTRA_PLATFORM "$AWD/get_platform"; then
    sed -i '1a\
\
# Added by install_astra8.sh: explicit platform override\
if [[ -n "${ASTRA_PLATFORM:-}" ]]; then echo "$ASTRA_PLATFORM"; exit 0; fi' "$AWD/get_platform"
fi
ok "platform/env.ubuntu written ($(ASTRA_PLATFORM=ubuntu "$AWD/get_platform"))"

# Parser fix: U-file lines separated by tabs ("TEX<tab>U-file: f  factor: 1.") and comment
# lines mentioning "U-file" crashed pyparse/exp_parser.py with IndexError.
if ! grep -q "U_LINE = re.compile" "$AWD/pyparse/exp_parser.py"; then
    git -C "$AWD" apply --whitespace=nowarn - <<'PATCH'
diff --git a/pyparse/exp_parser.py b/pyparse/exp_parser.py
--- a/pyparse/exp_parser.py
+++ b/pyparse/exp_parser.py
@@ -49,14 +49,18 @@ def append_x(str_in):
         return var_name(str_in) + 'X'


+# A8 form:  "TE   U-file:34954/TE34954.IDA_AVG: 1.e-3"
+# A7 form:  "TEX<tab>U-file: 00000.tes<tab>factor: 1."
+U_LINE = re.compile(r'^\s*(\S+)\s+U-FILE\s*:\s*([^\s:]+)\s*(?::\s*([-+0-9.EeDd]+)|.*?FACTOR\s*:\s*([-+0-9.EeDd]+))?', re.IGNORECASE)
+
+
 def parse_u_line(str_in):
-    var_name, str1 = str_in.split(' ', 1)
-    words = str1.split(':')
-    path = words[1].strip()
-    if len(words) == 2:
-        factor = 1.
-    elif len(words) == 3:
-        factor = float(words[2])
+    m = U_LINE.match(str_in)
+    if m is None:
+        return None, None, None
+    var_name, path, fac_a8, fac_a7 = m.groups()
+    factor = fac_a8 or fac_a7
+    factor = 1. if factor is None else float(factor.upper().replace('D', 'E'))
     return var_name, path, factor


@@ -197,6 +201,8 @@ class EXP_PARSER:

             if 'U-file' in line_strip: # u-file
                 varName, uf_path, factor = parse_u_line(line)
+                if varName is None: # free text mentioning a U-file, e.g. "****  ... U-file ..." comment
+                    continue
                 if varName in uf_keys: # double variable definition: 2x u-file, or u-file+exp_ascii
                     logger.error('>>> Error: var %s defined twice in u-files', varName)
                     sys.exit(4)
PATCH
    ok "pyparse/exp_parser.py patched"
else
    ok "pyparse/exp_parser.py already patched"
fi

# ---------------------------------------------------------------- 5. user area

step "Preparing the user area (equ, exp, exe, sbr, fnc, udb, ...)"
# Same as upstream install.sh: copy templates from repo_user_area. Existing user models
# and data are kept: only directories that do not exist yet are copied, except exe/
# (build scripts) which is always refreshed.
if cp --help | grep -q -- '--update.*none'; then
    cp_keep=(cp -r --update=none)    # coreutils >= 9.3 (Ubuntu 24.04+)
else
    cp_keep=(cp -rn)                 # older coreutils (Ubuntu 22.04)
fi
for dir in "$AWD"/repo_user_area/*/; do
    name="$(basename "$dir")"
    if [[ "$name" == "exe" || ! -d "$AWD/$name" ]]; then
        cp -r "$dir" "$AWD/"
    else
        # copy only files the user does not have yet
        "${cp_keep[@]}" "$dir". "$AWD/$name/"
    fi
done
mkdir -p "$AWD/tmp"
chmod u+x "$AWD"/exe/{Build,as_exe,wr_nml} "$AWD/pyparse/parser_main.py" "$AWD/clean.sh" "$AWD/compReg.sh" "$AWD/get_platform"

# TGLF / QuaLiKiZ / NEO helper programs need libraries that are not installed here
sed -i -e '/^all: directories/ s#\$(XPR)/qlki##; /^all: directories/ s#\$(XPR)/neo##; /^all: directories/ s#\$(XPR)/tglfi##' \
    "$AWD/exe/Makexpr"
# The reference model calls TORBEAM (ECRH) and RABBIT (NBI): closed IPP codes, not available
sed -i -e 's#^TORBA#!TORBA#' -e 's#^RABBIT#!RABBIT#' "$AWD/equ/flux_feqis"
ok "user area ready; TGLF/QuaLiKiZ/NEO, RABBIT, TORBEAM disabled"

# Full rebuild on every install, like upstream install.sh (compiler/BLAS may have changed)
(cd "$AWD" && export AWD && . platform/env.ubuntu && run make -f exe/Makefile clean)

# ---------------------------------------------------------------- 6. launcher

step "Installing the 'astra8' command"
mkdir -p "$HOME/.local/bin"
cat > "$HOME/.local/bin/astra8" <<EOF
#!/usr/bin/env bash
# ASTRA 8 launcher (install_astra8.sh). All arguments go to exe/as_exe, e.g.
#   astra8 -m flux_feqis -v aug34954 -s 4 -e 5          (GUI window)
#   astra8 -m flux_feqis -v aug34954 -s 4 -e 5 -batch   (no window)
#   astra8                                              (repeat the last run)
#   astra8 -h                                           (all options)
export ASTRA_PLATFORM=\${ASTRA_PLATFORM:-ubuntu}
export AWD=$AWD
cd "\$AWD"
exec exe/as_exe "\$@"
EOF
chmod +x "$HOME/.local/bin/astra8"
case ":$PATH:" in
    *":$HOME/.local/bin:"*) ok "astra8 -> $HOME/.local/bin/astra8" ;;
    *)
        # Ubuntu's ~/.profile adds ~/.local/bin only at the next login; add it for new terminals now
        if ! grep -qs 'astra8: ~/.local/bin' "$HOME/.bashrc"; then
            printf '\n# astra8: ~/.local/bin in PATH (install_astra8.sh)\nexport PATH="$HOME/.local/bin:$PATH"\n' >> "$HOME/.bashrc"
        fi
        warn "added ~/.local/bin to PATH in ~/.bashrc: open a new terminal to use 'astra8'" ;;
esac

# ---------------------------------------------------------------- 7. test

if [[ $DO_TEST -eq 1 ]]; then
    step "Building ASTRA and running the reference case (flux_feqis / aug34954, t=4..5 s)"
    warn "first build compiles all of ASTRA, usually 2-5 minutes"
    test_log="$PREFIX/test_flux_feqis.log"
    rm -f "$AWD/ncdf_out/aug34954flux_feqis.CDF"
    # as_exe always exits 0, so success is judged from the output
    ASTRA_PLATFORM=ubuntu AWD="$AWD" bash -c "cd '$AWD' && exe/as_exe -m flux_feqis -v aug34954 -s 4 -e 5 -batch" \
        >"$test_log" 2>&1 || true
    if grep -q "ASTRA normal exit" "$test_log" && [[ -f "$AWD/ncdf_out/aug34954flux_feqis.CDF" ]]; then
        ok "ASTRA normal exit; result: $AWD/ncdf_out/aug34954flux_feqis.CDF"
    else
        tail -n 30 "$test_log" >&2
        die "reference run failed, see $test_log"
    fi
fi

step "Done"
cat <<EOF

  ASTRA 8 is installed in $AWD

  Run with the graphic window:
      astra8 -m flux_feqis -v aug34954 -s 4 -e 5
  Run without window (results in $AWD/ncdf_out):
      astra8 -m flux_feqis -v aug34954 -s 4 -e 5 -batch

  Your models go to $AWD/equ, experiment files to $AWD/exp, U-files to $AWD/udb.
EOF
