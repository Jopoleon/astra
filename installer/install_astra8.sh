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
# >>> a7modules: SPIDER/STRAHL from the user's own ASTRA 7 copy (never downloaded)
A7_SRC=""                  # ASTRA 7 root (contains EQUIL/SPIDER/SRC and MAstraN/for)
SPIDER_SRC=""              # ASTRA 7 EQUIL/SPIDER/SRC
STRAHL_SRC=""              # ASTRA_STRAHL/STRAHL (str/strahl.f, atomdat, ...)
# <<< a7modules
# >>> gacode: TGLF / NEO from GACODE (public, Apache-2.0, github.com/gafusion/gacode)
WITH_TGLF=0
WITH_NEO=0
GACODE_URL="https://github.com/gafusion/gacode.git"
# master of 2026-09-30, verified with flux_tglf / flux_tglf_serial (full hash: shallow fetch)
GACODE_REF="f34a9dd2846afedba7a829d70f6194c2d8d8b304"
# <<< gacode

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
  --a7-src DIR      your local ASTRA 7 copy (e.g. .../ASTRA_7.02): builds SPIDER
                    from DIR/EQUIL/SPIDER/SRC (IPEQL=4); DIR/MAstraN/for is
                    needed by --strahl-src. Closed IPP code: never downloaded
  --spider-src DIR  SPIDER sources (ASTRA 7 EQUIL/SPIDER/SRC), instead of --a7-src
  --strahl-src DIR  STRAHL from the ASTRA 7 package (ASTRA_STRAHL/STRAHL) for
                    A2STRAHL; needs --a7-src for the ASTRA 7 include files
  --with-tglf       build TGLF (GACODE, open source) for the models tglf, flux_tglf,
                    flux_tglf_serial, ...; installs OpenMPI from apt (build only)
  --with-neo        build NEO (GACODE) for flux_neo (experimental: needs a lot of RAM)
  --gacode-ref REF  GACODE branch, tag or full commit hash (default: verified
                    ${GACODE_REF:0:8})
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
        --a7-src)     A7_SRC="$2"; shift 2 ;;       # a7modules
        --spider-src) SPIDER_SRC="$2"; shift 2 ;;   # a7modules
        --strahl-src) STRAHL_SRC="$2"; shift 2 ;;   # a7modules
        --with-tglf)  WITH_TGLF=1; shift ;;         # gacode
        --with-neo)   WITH_NEO=1; shift ;;          # gacode
        --gacode-ref) GACODE_REF="$2"; shift 2 ;;   # gacode
        -y|--yes)  ASSUME_YES=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) usage; die "unknown option: $1" ;;
    esac
done
[[ "$BLAS" == "mkl" || "$BLAS" == "openblas" ]] || die "--blas must be openblas or mkl"
# >>> a7modules: resolve and check the local ASTRA 7 paths before anything is installed
A7MOD_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/a7modules"
A7_FOR=""
if [[ -n "$A7_SRC" ]]; then
    A7_SRC="$(cd "$A7_SRC" && pwd)" || die "--a7-src: no such directory"
    if [[ -z "$SPIDER_SRC" && -d "$A7_SRC/EQUIL/SPIDER/SRC" ]]; then SPIDER_SRC="$A7_SRC/EQUIL/SPIDER/SRC"; fi
    if [[ -f "$A7_SRC/MAstraN/for/parameter.inc" ]]; then A7_FOR="$A7_SRC/MAstraN"; fi
fi
if [[ -n "$SPIDER_SRC" ]]; then
    SPIDER_SRC="$(cd "$SPIDER_SRC" && pwd)" || die "--spider-src: no such directory"
    [[ -f "$SPIDER_SRC/interfac.f" && -f "$SPIDER_SRC/parameters.f90" ]] \
        || die "--spider-src: $SPIDER_SRC has no interfac.f/parameters.f90 (need the full ASTRA 7 EQUIL/SPIDER/SRC)"
fi
if [[ -n "$STRAHL_SRC" ]]; then
    STRAHL_SRC="$(cd "$STRAHL_SRC" && pwd)" || die "--strahl-src: no such directory"
    [[ -f "$STRAHL_SRC/str/strahl.f" ]] || die "--strahl-src: $STRAHL_SRC/str/strahl.f not found (need ASTRA_STRAHL/STRAHL)"
    [[ -n "$A7_FOR" ]] || die "--strahl-src needs --a7-src (STRAHL includes MAstraN/for/*.inc of ASTRA 7)"
fi
if [[ -n "$SPIDER_SRC$STRAHL_SRC" ]]; then
    [[ -d "$A7MOD_DIR" ]] || die "$A7MOD_DIR not found: run install_astra8.sh from the repository checkout"
fi
# <<< a7modules

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

# ASTRA 7 compatibility (exp files, models and fml of A7 user areas; gfortran instead of ifort).
# All of it is checked to leave the A8 test cases unchanged (same pyparse output, same Fortran code):
#  pyparse/exp_parser.py  U-file lines with tabs/"factor:", 1D radial U-files, "POINTS n" blocks,
#                         omitted GRIDTYPE, comments and text inside data blocks, "2.0d0" numbers,
#                         values crossing the fixed columns; never undercounts nr_x_max
#  src/astra/read_input.f90  inline "! comments" (a "points" in a comment opened a group), values
#                         crossing the fixed columns (ifort read "0.45d" cut out of "0.45d0", gfortran
#                         does not), no garbage size added for non-U-file lines
#  src/astra/parse_utils.f90  split2array2 wrote past its 20 words on long A7 data lines (segfault)
#  pyparse/parse_as.py    "((A)*2)/(B)" lost a bracket (ValueError in getInnermostBracket)
#  pyparse/const_text.py  inivar did not "use nclass_mod", so "cc_nc(j)" etc. failed in initial values
#  repo_user_area/sbr/a2eqdsk.f90  EQDSK name in a 120-char buffer: "End of record" (SPIDER runs)
#                         when the install path is longer than ~70 chars
if ! grep -q "fortran_float" "$AWD/pyparse/exp_parser.py"; then
    git -C "$AWD" apply --whitespace=nowarn - <<'PATCH' || die "A7 compatibility patch does not apply to this ASTRA version (--ref); use the default ref"
diff --git a/pyparse/exp_parser.py b/pyparse/exp_parser.py
--- a/pyparse/exp_parser.py
+++ b/pyparse/exp_parser.py
@@ -49,14 +49,18 @@
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
 
 
@@ -66,7 +70,8 @@
     for i, word in enumerate(words[:-1]):   # avoid overflow
         if word in prof_attr_types:
             try:
-                result[word] = prof_attr_types[word](words[i + 1])
+                typ = prof_attr_types[word]
+                result[word] = fortran_float(words[i + 1]) if typ is float else typ(words[i + 1])
             except ValueError:
                 raise ValueError(f"Invalid value for {word}: {words[i+1]}")
     if 'NAMEXP' in result:
@@ -75,10 +80,49 @@
 
 
 def read_float_block(lines, start_line, end_line):
-    block_lines = lines[start_line:end_line]   # end exclusive
-    block_text = " ".join(block_lines)         # normalize spacing
-    values = np.fromstring(block_text, sep=' ')
-    return values
+# Numbers only: skip comment lines ('*', '!'), and stop each line at the first non-number
+# (A7 files carry text after the data, e.g. ".2  .3  these are time[s]")
+    values = []
+    for lin in lines[start_line:end_line]:   # end exclusive
+        if lin.lstrip()[:1] in ('*', '!'):
+            continue
+        for word in lin.split():
+            try:
+                values.append(float(word.upper().replace('D', 'E')))
+            except ValueError:
+                break
+    return np.array(values)
+
+
+def exp_field(line, c1, c2, cmin):
+# Fixed-column field line[c1:c2], widened to the whole word it cuts (not left of cmin);
+# same as exp_field in read_input.f90. A7 files have e.g. "0.45d0" crossing column 22.
+    n = len(line.rstrip())
+    j1, j2 = max(c1, cmin), min(c2, n)
+    if j2 <= j1:
+        return '', c2
+    while j1 > cmin and not line[j1].isspace() and not line[j1-1].isspace():
+        j1 -= 1
+    while j2 < n and not line[j2-1].isspace() and not line[j2].isspace():
+        j2 += 1
+    return line[j1: j2].strip(), max(j2, c2)
+
+
+def fortran_float(str_in):
+# Accepts Fortran exponents, e.g. "2.0d0" (A7 files)
+    return float(str_in.strip().upper().replace('D', 'E'))
+
+
+def uf_scalar_time(f_in):
+# A7 1D U-file: time is the first associated scalar ("; Scalar value, LABEL FOLLOWS:")
+    with open(f_in, 'r') as f:
+        for lin in f:
+            if 'Scalar value' in lin or '-SCALAR, LABEL FOLLOWS' in lin:
+                try:
+                    return float(lin.split()[0])
+                except ValueError:
+                    break
+    return 0.
 
 
 def n_times(attr_d):
@@ -129,9 +173,9 @@
 
         for line in lines:
             varName = var_name(line[:6]) # remove trailing 'X', if any
-            tim   = line[8: 14]
-            data  = line[16: 22]
-            error = line[24: 30]
+            tim, jcol   = exp_field(line,  8, 14, 6)
+            data, jcol  = exp_field(line, 16, 22, jcol)
+            error, jcol = exp_field(line, 24, 30, jcol)
             if varName not in self.variables + ['NA1', 'TSTART', 'TEND']:
                 continue
             if varName in ('NA1', 'TSTART', 'TEND', 'RTOR', 'AWALL', 'AB', 'ELONM', 'TRICH', 'AMJ', 'ZMJ'):
@@ -152,9 +196,9 @@
                 self.scalars[varName]['error'] = len(self.scalars[varName]['data'])*[0.]
                 uf_keys.append(varName)
             else: # ASCII blocks in exp file
-                self.scalars[varName]['time'].append((float(tim) if tim.strip() else 0.))
-                self.scalars[varName]['error'].append((float(error) if error.strip() else 0.))
-                self.scalars[varName]['data'].append(float(data))
+                self.scalars[varName]['time'].append((fortran_float(tim) if tim.strip() else 0.))
+                self.scalars[varName]['error'].append((fortran_float(error) if error.strip() else 0.))
+                self.scalars[varName]['data'].append(fortran_float(data))
                 
 
     def parse_exp2d(self, exp2d):
@@ -197,29 +241,40 @@
 
             if 'U-file' in line_strip: # u-file
                 varName, uf_path, factor = parse_u_line(line)
+                if varName is None: # free text mentioning a U-file, e.g. "****  ... U-file ..." comment
+                    continue
                 if varName in uf_keys: # double variable definition: 2x u-file, or u-file+exp_ascii
                     logger.error('>>> Error: var %s defined twice in u-files', varName)
                     sys.exit(4)
                 f_in = '%s/udb/%s' %(config.awd, uf_path)
                 uf = ufiles.UFILE(fin=f_in)
-                gridtype = grid_d[uf.Y['label'].strip().upper()]
+                if 'Y' in uf.__dict__: # 2D: time x radius
+                    tvec, rad = uf.X['data'], uf.Y
+                else: # A7 1D profile: radius only, time is the associated scalar
+                    tvec, rad = np.atleast_1d(uf_scalar_time(f_in)), uf.X
+                gridtype = grid_d.get(rad['label'].replace(' ', '').upper()[:12], 10)
                 varName = append_x(varName)
                 self.profiles['grid_type'].append(gridtype)
                 self.profiles['label'].append(varName)
-                self.profiles['time'].append(uf.X['data'].tolist())
-                self.profiles['rho'].append(uf.Y['data'].tolist())
+                self.profiles['time'].append(tvec.tolist())
+                self.profiles['rho'].append(rad['data'].tolist())
                 self.profiles['data'].append((factor*uf.f['data']).ravel().tolist())
                 self.profiles['filter'].append(alpha)
                 uf_keys.append(varName)
 
             else: # A line of the kind "NAMEXP ...  GRIDTYPE ...", not u-file
+                if line_strip.startswith('*'): # A7 comment line
+                    continue
                 attr_d = parse_line2d(line)
                 if 'FILTER' in attr_d:
                     alpha = attr_d['FILTER']
                     if line.startswith('FILTER'):
                         alpha_glob = alpha
                 if 'NAMEXP' not in attr_d:
+                    if 'POINTS' in attr_d: # A7 old standard: "POINTS n", then "NAMEX time v1 ... vn" lines
+                        self.add_a7_block(lines[jlin+1: jnext], attr_d['POINTS'])
                     continue
+                attr_d.setdefault('GRIDTYPE', 1) # A7: GRIDTYPE may be omitted
 
                 varName = attr_d['NAMEXP'] # already with trailing 'X'
                 if varName in uf_keys: # double variable definition: 2x u-file, or u-file+exp_ascii
@@ -269,13 +324,42 @@
                         n_grid = nrho
                     else:
                         n_grid = 2*nrho
-                    assert dataStream.size == nt + n_grid + nt*nrho
+                    if dataStream.size != nt + n_grid + nt*nrho:
+                        # A7 forms (GRIDTYPE 18/19 with an extra R/z line, text after numbers, ...):
+                        # the Fortran reader handles them; here only the grid size matters (nr_x_max)
+                        logger.warning('%s: block not parsed by pyparse (A7 form?), left to the Fortran reader', varName)
+                        self.profiles['time'].append(dataStream[:nt].tolist())
+                        self.profiles['rho'].append(n_grid*[0.])
+                        self.profiles['data'].append([])
+                        continue
                     time, rho, data = np.split(dataStream, [nt, nt+n_grid])
                     self.profiles['time'].append(time.tolist())
                     self.profiles['rho'].append(rho.tolist())
                     self.profiles['data'].append(data.tolist())
 
 
+    def add_a7_block(self, block, nrho):
+# A7 old standard: every line "NAMEX  time  v1 ... v_nrho" is one time slice on the default grid
+
+        for blin in block:
+            words = blin.split()
+            if not words or blin[0] in (' ', '\t', '*', '!'):
+                continue
+            varName = append_x(words[0])
+            if varName not in self.profx:
+                continue
+            if varName not in self.profiles['label']:
+                self.profiles['label'].append(varName)
+                self.profiles['grid_type'].append([0])
+                self.profiles['filter'].append([0.])
+                self.profiles['time'].append([])
+                self.profiles['rho'].append([])
+                self.profiles['data'].append([])
+            jvar = self.profiles['label'].index(varName)
+            if len(self.profiles['rho'][jvar]) < nrho:
+                self.profiles['rho'][jvar] = nrho*[0.]
+
+
     def write_json(self):
 
         data = {"coils": self.coils, "scalars": self.scalars, "boundary": self.boundary, "profiles": self.profiles}
diff --git a/pyparse/parse_as.py b/pyparse/parse_as.py
--- a/pyparse/parse_as.py
+++ b/pyparse/parse_as.py
@@ -345,9 +345,13 @@
     if pieces_in[jleft-1] == 'AFX':
         jleft -= 2
         jright += 1
-    pieces_within = pieces_in[jleft-1: jright+1] # function, '(', ..., ')'
+# Take the preceding piece only if it is a function/array name: in "((A)*2)/(B)" it is '(',
+# and swallowing it left an unmatched ')' that broke the next bracket
+    if jleft > 0 and re.match(r'[A-Za-z_]\w*$', pieces_in[jleft-1]):
+        jleft -= 1
+    pieces_within = pieces_in[jleft: jright+1] # [function,] '(', ..., ')'
     str_mid = parse_pieces(pieces_within, parse, idx=idx)
-    pieces_out = pieces_in[:jleft-1] + [str_mid] + pieces_in[jright+1:]
+    pieces_out = pieces_in[:jleft] + [str_mid] + pieces_in[jright+1:]
     return pieces_out
 
 
diff --git a/pyparse/const_text.py b/pyparse/const_text.py
--- a/pyparse/const_text.py
+++ b/pyparse/const_text.py
@@ -157,6 +157,7 @@
 use scalars
 use status
 use pi_const
+use nclass_mod
 use standard_functions
 use debugger, only: markloc
 use json_vars, only: profxNames, n_profx
diff --git a/src/astra/read_input.f90 b/src/astra/read_input.f90
--- a/src/astra/read_input.f90
+++ b/src/astra/read_input.f90
@@ -233,7 +233,9 @@
     double precision, allocatable, dimension(:) :: t_u, x_u, var_u, bnd_rz, bnd_r, bnd_z
     double precision :: XBDRY, YB, YB1, YXB, YXB1, ALFA, ALFA_GLOB, &
         VRDATA, FACTOR, TIMEVR, VRERR
-    character(len=6) :: VNAM, VNAMO, VNAMU, VNAMX, VTIM, VDAT, VERR, keyword
+    character(len=6) :: VNAM, VNAMO, VNAMU, VNAMX, keyword
+    character(len=32) :: VTIM, VDAT, VERR
+    integer :: jcol
     character(len=31) :: rholbl
     character(len=132) :: strarray(20), STRI, lin_upper, &
         err_msg, err_msg_exp, file_exp, ufile_in, uname, uvar
@@ -263,6 +265,8 @@
     eof_1d = .false.
     set_dims_1d: do
         read(n_unit, '(A132)', iostat=ios) STRI
+        jcol = index(STRI, '!')  ! cut inline comments (A7 files)
+        if (jcol > 1) STRI(jcol:) = ' '
 ! An exp that states only scalars has no group header, so end of file is a
 ! normal exit here.  eof_1d tells the 2D pass below that STRI holds no group
 ! header and that nothing is left to read.
@@ -350,6 +354,8 @@
             skip_read = .false.
         else
             read(n_unit, '(A132)', iostat=ios) STRI
+            jcol = index(STRI, '!')  ! cut inline comments (A7 files)
+            if (jcol > 1) STRI(jcol:) = ' '
             if (ios < 0) EXIT set_dims_2d
             if (ios > 0) then
                 write(*, '(A)') err_format
@@ -397,9 +403,9 @@
             call parse_u_line(STRI, uvar, ufile_in, factor)
             if (LEN_TRIM(ufile_in) > 0) then ! string 'U-FILE' found in this line
                 call ufheader(TRIM(ufile_in), nscal_u, ndim_u, nt_u, nx_u, rholbl)
+                len_profs_data = len_profs_data + nt_u*nx_u + nx_u
+                len_profs_time = len_profs_time + nt_u
             endif
-            len_profs_data = len_profs_data + nt_u*nx_u + nx_u
-            len_profs_time = len_profs_time + nt_u
         endif
     enddo set_dims_2d
     close(n_unit)
@@ -453,6 +459,8 @@
     parse_exp_1d: do
 
         read(n_unit, '(A132)', iostat=ios) STRI
+        jcol = index(STRI, '!')  ! cut inline comments (A7 files)
+        if (jcol > 1) STRI(jcol:) = ' '
         if (ios < 0) then
             close(n_unit)
             raw_profiles%n_groups = NGR
@@ -514,9 +522,9 @@
         else
 ! Variable name, time, value, error
             VNAM = VARNAM(lin_upper(1: 6), ier_tab)
-            VTIM = STRI(9 : 14)
-            VDAT = STRI(17: 22)
-            VERR = STRI(25: 30)
+            call exp_field(STRI,  9, 14, 7       , VTIM, jcol)
+            call exp_field(STRI, 17, 22, jcol + 1, VDAT, jcol)
+            call exp_field(STRI, 25, 30, jcol + 1, VERR, jcol)
         endif
 
 ! A special treatment required for time independent NA1, TSTART, TEND
@@ -643,6 +651,8 @@
             skip_read = .false.
         else
             read(n_unit, '(A132)', iostat=ios) STRI
+            jcol = index(STRI, '!')  ! cut inline comments (A7 files)
+            if (jcol > 1) STRI(jcol:) = ' '
             if (ios < 0) EXIT parse_exp_2d
             if (ios > 0) then
                 write(*, '(A)') err_format
@@ -993,6 +1003,38 @@
     end subroutine read_exp
 
 !---------------------------------------------------------------------
+    subroutine exp_field(line, c1, c2, cmin, field, jend)
+! Fixed-column field line(c1:c2) of an old-standard exp line, widened to the whole word it
+! cuts (not left of cmin). A7 files were read by ifort, which accepts "0.45d" cut out of
+! "0.45d0"; gfortran does not. jend: last column used, to start the next field after it.
+
+    character(len=*), intent(in)  :: line
+    integer,          intent(in)  :: c1, c2, cmin
+    character(len=*), intent(out) :: field
+    integer,          intent(out) :: jend
+
+    integer :: j1, j2, n
+
+    n  = LEN_TRIM(line)
+    j1 = max(c1, cmin)
+    j2 = min(c2, n)
+    field = ' '
+    jend  = c2
+    if (j2 < j1) return
+    do while (j1 > cmin)
+        if (line(j1: j1) == ' ' .or. line(j1-1: j1-1) == ' ') EXIT
+        j1 = j1 - 1
+    enddo
+    do while (j2 < n)
+        if (line(j2: j2) == ' ' .or. line(j2+1: j2+1) == ' ') EXIT
+        j2 = j2 + 1
+    enddo
+    field = ADJUSTL(line(j1: j2))
+    jend  = max(j2, c2)
+
+    end subroutine exp_field
+
+
     subroutine str2dbl(line, dbl_out, ierr)
 
     character(len=*), intent(in)  :: line
diff --git a/src/astra/parse_utils.f90 b/src/astra/parse_utils.f90
--- a/src/astra/parse_utils.f90
+++ b/src/astra/parse_utils.f90
@@ -147,6 +147,8 @@
     m = 1
     nout = 0
     do i=1, len(strtmp)
+! Only header words are needed; A7 data lines starting in column 1 can hold more words
+        if (nout >= nwords_max) return
         jpos_null = index(strtmp(m:), ' ')
         jpos_tab  = index(strtmp(m:), tab_ch)
         if (jpos_null > 0) then
@@ -166,7 +168,7 @@
     enddo
 
 ! After the last delimiter
-    if ( LEN_TRIM(strtmp(m: )) > 0 ) then
+    if ( LEN_TRIM(strtmp(m: )) > 0 .and. nout < nwords_max ) then
         nout = nout + 1
         strarray(nout) = TRIM(ADJUSTL( strtmp(m:) )) 
     endif
diff --git a/repo_user_area/sbr/a2eqdsk.f90 b/repo_user_area/sbr/a2eqdsk.f90
--- a/repo_user_area/sbr/a2eqdsk.f90
+++ b/repo_user_area/sbr/a2eqdsk.f90
@@ -28,7 +28,7 @@
      fdia_rect, q_rect, pprime_rect, fprime_rect
 double precision, dimension(nzRect) :: z_rect
 double precision, dimension(nrRect, nzRect) :: psi_rect
-character(len=120) :: f_eqdsk
+character(len=512) :: f_eqdsk
 
 double precision :: rdim, zdim, rcentr, rleft, zmid, rmaxis, zmaxis, &
     simag, sibry, bcentr, current, xdum
PATCH
    ok "A7 compatibility patch applied (pyparse, read_input, parse_utils)"
else
    ok "A7 compatibility patch already applied"
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

# >>> a7modules ---------------------------------------------------- 5b. SPIDER / STRAHL
# SPIDER (equilibrium, IPEQL=4) and STRAHL (impurities, A2STRAHL) are closed IPP codes.
# Nothing is downloaded: the sources come from the user's own ASTRA 7 copy (--a7-src,
# --spider-src, --strahl-src), are copied to $BUILD and compiled there. The glue in
# installer/a7modules is ours (no IPP code). Details: docs/research/spider-strahl-from-a7-2026-10-02.md
HAVE_SPIDER=0
HAVE_STRAHL=0
if [[ -n "$SPIDER_SRC" ]]; then
    step "Building SPIDER from $SPIDER_SRC (ASTRA 7 sources)"
    # exe/astra_rc expects $ASTRA_EXT/spider/<version>/{lib/libspider.a,inc}
    spider_ver="$(sed -n 's#^SPIDER_HOME=.*/spider/\([A-Za-z0-9._-]*\).*#\1#p' "$AWD/exe/astra_rc" | head -n 1)"
    spider_ver="${spider_ver:-latest}"
    sb="$BUILD/spider"
    rm -rf "$sb"
    mkdir -p "$sb"
    cp "$SPIDER_SRC"/*.f "$SPIDER_SRC"/*.f90 "$SPIDER_SRC"/*.inc "$sb"/
    # schemas.f90 is replaced by ours (A8 imas_ids types are a superset of the A7 CPO types)
    cp "$AWD/src/misc/imas_ids.f90" "$A7MOD_DIR"/spider/*.f90 "$sb"/
    # A8 renamed teta2d -> theta2d; A8 calls spider_run with its own parameter type,
    # so the A7 entry is renamed and wrapped by spider_a8_adapter.f90
    sed -i 's/%teta2d/%theta2d/g' "$sb/interfac.f"
    sed -i 's/^\(\s*subroutine\s\+\)spider_run(/\1spider_run_a7(/I' "$sb/interfac.f"
    grep -qi 'subroutine spider_run_a7(' "$sb/interfac.f" || die "SPIDER: no 'subroutine spider_run' in $SPIDER_SRC/interfac.f"
    # object list = esp.a of the A7 SPIDER Makefile + glue; imas_ids is only for its .mod
    spider_objs=(schemas vol2d_eq_face parameters spider_params bnd_modul ppf_modul durs_d_modul
        _bnd _eqa_m_r _eqain _ext_m _extp _fluxt _mesh _metric _rig _sol _wrd B_eqb b_extrp
        B_grid B_MESH B_metric_O B_rig B_sol B_wrd cf_ini com_sub controller Curap curfit_f
        DPCGRCV ENDFL EQ1_m Eq2_m Eq_m Ev1 Ev2 Ev3 interfac interp.cubsplinetau_nonsym Matrix
        sfl out42_dummy exp_promat_e_I_f0 RIG scu solc solver SPAR1 Spider_call sstepon
        sstepon_r Trecur_2 Urs wrc coil2spider spidupdate spider_a8_adapter)
    spider_ff=(-O2 -fPIC -w -std=legacy -fallow-argument-mismatch -I.)
    pushd "$sb" >/dev/null
    run gfortran "${spider_ff[@]}" -c imas_ids.f90 -o imas_ids.o
    for o in "${spider_objs[@]}"; do
        if [[ -f "$o.f90" ]]; then f="$o.f90"; else f="$o.f"; fi
        run gfortran "${spider_ff[@]}" -c "$f" -o "$o.o"
    done
    sp="$EXT/spider/$spider_ver"
    mkdir -p "$sp/lib" "$sp/inc"
    rm -f "$sp/lib/libspider.a" "$sp"/inc/*.mod
    run ar rcs "$sp/lib/libspider.a" "${spider_objs[@]/%/.o}"
    for m in *.mod; do
        if [[ "$m" != imas_ids.mod ]]; then cp "$m" "$sp/inc/"; fi
    done
    popd >/dev/null
    printf 'SPIDER from %s, built by install_astra8.sh\n' "$SPIDER_SRC" > "$sp/source"
    HAVE_SPIDER=1
    ok "$sp/lib/libspider.a (ASTRA builds with -DSPIDER; IPEQL=4 selects SPIDER)"
fi
if [[ -n "$STRAHL_SRC" ]]; then
    step "Building STRAHL from $STRAHL_SRC (ASTRA 7 package)"
    # sbr/a2strahl.f90 runs $ASTRA_EXT/strahl/<version>/bin/{strahl,result_to_astra}
    strahl_ver="$(grep -o "strahl/[A-Za-z0-9._-]*/bin/strahl'" "$AWD/sbr/a2strahl.f90" | head -n 1 | cut -d/ -f2 || true)"
    strahl_ver="${strahl_ver:-sep23}"
    tb="$BUILD/strahl"
    rm -rf "$tb"
    mkdir -p "$tb"
    cp -a "$STRAHL_SRC"/. "$tb"/
    cp "$A7MOD_DIR/strahl/strahl_main.f" "$tb/str/"
    # netcdf.inc must match the installed libnetcdff
    cp "$(nf-config --includedir 2>/dev/null || echo /usr/include)/netcdf.inc" "$tb/str/netcdf.inc"
    # gfortran does not accept .eq. between logicals (an ifort extension)
    sed -i 's/neoclass\.eq\.\.false\./neoclass.eqv..false./' "$tb/str/read_grid.f"
    # STRAHL as a standalone program ("strahl a q", as A8 calls it); the A7 in-process
    # coupling routines are stubbed in strahl_main.f
    strahl_objs=(strahl time_steps read_parameter neutrals sput_yield save_param atomic_data
        plad read_atomdat emissiv dtfus math num_recip strahl_util netcdfutil impden
        neo_trans_atp neo_trans neoarchive saw_mix read_grid radiation_diag radiation_total
        save_dist strahl_main)
    pushd "$tb/str" >/dev/null
    for o in "${strahl_objs[@]}"; do
        run gfortran -O2 -w -std=legacy -fallow-argument-mismatch -I. -I.. -I"$A7_FOR" -c "$o.f" -o "$o.o"
    done
    run gfortran -o strahl_a7 "${strahl_objs[@]/%/.o}" -lnetcdff -lnetcdf
    popd >/dev/null
    stb="$EXT/strahl/$strahl_ver/bin"
    mkdir -p "$stb"
    cp "$tb/str/strahl_a7" "$stb/"
    cp "$A7MOD_DIR/strahl/strahl.sh" "$stb/strahl"
    cp "$A7MOD_DIR/strahl/result_to_astra.sh" "$stb/result_to_astra"
    cp "$A7MOD_DIR/strahl/sparams_a8_to_a7.py" "$A7MOD_DIR/strahl/result_to_astra.py" "$stb/"
    chmod +x "$stb"/*
    printf 'STRAHL from %s, built by install_astra8.sh\n' "$STRAHL_SRC" > "$stb/../source"
    HAVE_STRAHL=1
    ok "$stb/strahl (+ result_to_astra)"
    if [[ $HAVE_SPIDER -eq 0 ]]; then
        warn "equ/astrahl_simple_msp uses IPEQL=4: without SPIDER the equilibrium is not updated"
    fi
fi
# flux_spider calls TORBEAM/RABBIT (closed, not installed); SPIDER writes its
# work files (outp.wr, q.wr, ...) into exp/equ/ and stops if it does not exist
if [[ $HAVE_SPIDER -eq 1 ]]; then
    mkdir -p "$AWD/exp/equ"
    sed -i -e 's#^TORBA#!TORBA#' -e 's#^RABBIT#!RABBIT#' "$AWD/equ/flux_spider"
fi
# <<< a7modules

# >>> gacode -------------------------------------------------------- 5c. TGLF / NEO (GACODE)
# TGLF and NEO are open (github.com/gafusion/gacode, Apache-2.0). Only the libraries ASTRA
# links are built (no cgyro etc.): a sparse shallow checkout of tglf/, neo/, shared/, f2py/
# (~13 MB) in $BUILD/gacode, compiled with GACODE's own makefiles and our platform file.
# Layout expected by exe/astra_rc:  $ASTRA_EXT/tglf/<ver>/{lib/libtglf.a,inc/*.mod}
#   $ASTRA_EXT/neo/<ver>/{lib/{neo,nclass,UMFPACK,math,expro,geo}_lib.a,inc/*.mod}
# Once installed, TGLF/NEO stay enabled on later runs of the installer without the flags.
gacode_ver() {   # module -> version directory from exe/astra_rc (e.g. tglf -> jun25)
    local v
    v="$(sed -n "s#^${1^^}_HOME=.*/$1/\([A-Za-z0-9._-]*\).*#\1#p" "$AWD/exe/astra_rc" | head -n 1)"
    echo "${v:-latest}"
}
TGLF_DIR="$EXT/tglf/$(gacode_ver tglf)"
NEO_DIR="$EXT/neo/$(gacode_ver neo)"
if [[ $WITH_TGLF -eq 1 || $WITH_NEO -eq 1 ]]; then
    step "GACODE (TGLF/NEO) @ $GACODE_REF"
    if [[ $DO_APT -eq 1 ]]; then
        # mpif90 is needed to compile GACODE (tglf_init_mpi.f90 "use mpi") and to link xpr/neo;
        # xpr/tglfi itself links without MPI
        run $SUDO env DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
            openmpi-bin libopenmpi-dev
    fi
    command -v mpif90 >/dev/null || die "mpif90 not found (apt install openmpi-bin libopenmpi-dev, or run without --no-apt)"
    GACODE_ROOT="$BUILD/gacode"
    if [[ ! -d "$GACODE_ROOT/.git" ]]; then
        mkdir -p "$GACODE_ROOT"
        run git -c core.autocrlf=false init "$GACODE_ROOT"
        run git -C "$GACODE_ROOT" remote add origin "$GACODE_URL"
    fi
    git -C "$GACODE_ROOT" config core.autocrlf false
    run git -C "$GACODE_ROOT" fetch --depth 1 --filter=blob:none origin "$GACODE_REF" \
        || die "cannot fetch GACODE '$GACODE_REF' (a commit needs the full hash)"
    run git -C "$GACODE_ROOT" sparse-checkout set tglf neo shared f2py/expro f2py/geo platform/build
    run git -C "$GACODE_ROOT" checkout --force --detach FETCH_HEAD
    gacode_sha="$(git -C "$GACODE_ROOT" rev-parse HEAD)"
    ok "GACODE at ${gacode_sha:0:8} ($(git -C "$GACODE_ROOT" log -1 --format=%cd --date=short))"
    mkdir -p "$GACODE_ROOT/modules"
    # Our platform: gfortran via the OpenMPI wrapper. -fallow-argument-mismatch: gfortran >= 10
    # rejects the old LAPACK-style calls; -fPIC as for json-fortran.
    cat > "$GACODE_ROOT/platform/build/make.inc.ASTRA_UBUNTU" <<'EOF'
# Written by install_astra8.sh: gfortran + OpenMPI from Ubuntu apt
FC     = mpif90 -fPIC -I$(GACODE_ROOT)/modules -J$(GACODE_ROOT)/modules
F77    = mpif90 -fPIC
CC     = mpicc -fPIC
FMATH  = -fdefault-real-8 -fdefault-double-8
FOPT   = -O2 -fallow-argument-mismatch
FDEBUG = -O0 -g -fcheck=all -fbacktrace -fallow-argument-mismatch
FOMP   = -fopenmp
LMATH  = -lopenblas
ARCH   = ar cr
EOF
    # GACODE makefiles have no module dependencies: build serially, in dependency order
    gacode_make() { run env GACODE_ROOT="$GACODE_ROOT" GACODE_PLATFORM=ASTRA_UBUNTU make -C "$GACODE_ROOT/$1" "${@:2}"; }
    if [[ $WITH_TGLF -eq 1 ]]; then
        if [[ -f "$TGLF_DIR/lib/libtglf.a" && "$(cat "$TGLF_DIR/GACODE_COMMIT" 2>/dev/null)" == "$gacode_sha" ]]; then
            ok "TGLF already built in $TGLF_DIR"
        else
            gacode_make tglf/src clean
            gacode_make tglf/src tglf_lib.a
            rm -rf "$TGLF_DIR"
            mkdir -p "$TGLF_DIR/lib" "$TGLF_DIR/inc"
            cp "$GACODE_ROOT/tglf/src/tglf_lib.a" "$TGLF_DIR/lib/libtglf.a"
            cp "$GACODE_ROOT"/modules/tglf*.mod "$TGLF_DIR/inc/"
            echo "$gacode_sha" > "$TGLF_DIR/GACODE_COMMIT"
            ok "TGLF installed in $TGLF_DIR"
        fi
    fi
    if [[ $WITH_NEO -eq 1 ]]; then
        if [[ -f "$NEO_DIR/lib/neo_lib.a" && "$(cat "$NEO_DIR/GACODE_COMMIT" 2>/dev/null)" == "$gacode_sha" ]]; then
            ok "NEO already built in $NEO_DIR"
        else
            for d in shared/math shared/nclass shared/UMFPACK f2py/geo f2py/expro neo/src; do
                gacode_make "$d" clean
            done
            for d in shared/math shared/nclass shared/UMFPACK f2py/geo f2py/expro; do
                gacode_make "$d"
            done
            gacode_make neo/src neo_lib.a
            rm -rf "$NEO_DIR"
            mkdir -p "$NEO_DIR/lib" "$NEO_DIR/inc"
            cp "$GACODE_ROOT"/neo/src/neo_lib.a "$GACODE_ROOT"/shared/{nclass/nclass,UMFPACK/UMFPACK,math/math}_lib.a \
               "$GACODE_ROOT"/f2py/{expro/expro,geo/geo}_lib.a "$NEO_DIR/lib/"
            cp "$GACODE_ROOT"/modules/*.mod "$NEO_DIR/inc/"
            echo "$gacode_sha" > "$NEO_DIR/GACODE_COMMIT"
            ok "NEO installed in $NEO_DIR"
        fi
    fi
fi
HAVE_TGLF=0
HAVE_NEO=0
if [[ -f "$TGLF_DIR/lib/libtglf.a" ]]; then HAVE_TGLF=1; fi
if [[ -f "$NEO_DIR/lib/neo_lib.a" ]]; then HAVE_NEO=1; fi
if [[ $HAVE_TGLF -eq 1 || $HAVE_NEO -eq 1 ]]; then
    # step 5 removed the IPC programs from exe/Makexpr: put back the ones that can be linked
    if [[ $HAVE_TGLF -eq 1 ]]; then
        sed -i '/^all: directories/ { /tglfi/! s#$# $(XPR)/tglfi# }' "$AWD/exe/Makexpr"
    fi
    if [[ $HAVE_NEO -eq 1 ]]; then
        sed -i '/^all: directories/ { /XPR)\/neo/! s#$# $(XPR)/neo# }' "$AWD/exe/Makexpr"
    fi
    # exe/Build overwrites MAX_NWORKERS with the number of physical cores for GUI runs, but the
    # IPC wrappers (sbr/tglf_ipc.f90, neo_ipc.f90) always start 64 workers and stop when
    # MAX_NWORKERS < 64. Let ASTRA_MAX_NWORKERS (set in env.ubuntu) win.
    if ! grep -q ASTRA_MAX_NWORKERS "$AWD/exe/Build"; then
        if grep -q "^    echo 'Max nworkers'" "$AWD/exe/Build"; then
            sed -i "/^    echo 'Max nworkers'/i\\
    export MAX_NWORKERS=\${ASTRA_MAX_NWORKERS:-\$MAX_NWORKERS}   # install_astra8.sh (TGLF/NEO)" "$AWD/exe/Build"
        else
            warn "exe/Build: MAX_NWORKERS line not found; GUI runs of TGLF/NEO models may stop"
        fi
    fi
    {
        echo ""
        echo "# TGLF/NEO from GACODE (install_astra8.sh --with-tglf/--with-neo)."
        echo "# The IPC wrappers always start 64 worker processes (one per radial point) and stop"
        echo "# if MAX_NWORKERS < 64; on a PC they simply share the cores."
        echo 'export ASTRA_MAX_NWORKERS=${ASTRA_MAX_NWORKERS:-64}'
        echo 'export MAX_NWORKERS=$ASTRA_MAX_NWORKERS'
        if [[ $HAVE_TGLF -eq 1 ]]; then
            # exe/Makesub compiles sbr/tglf_serial.F90 with $(TGLF_FLAG) but never defines it,
            # so tglf_serial would always be the "not installed" stub
            echo 'export TGLF_FLAG=-DTGLF_INSTALLED'
        fi
    } >> "$AWD/platform/env.ubuntu"
    enabled=""
    if [[ $HAVE_TGLF -eq 1 ]]; then enabled+=" TGLF (xpr/tglfi, tglf_serial)"; fi
    if [[ $HAVE_NEO -eq 1 ]]; then enabled+=" NEO (xpr/neo)"; fi
    ok "enabled:$enabled"
fi
# <<< gacode

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
# A7 -> A8 helper (models, fml, equ/log): only when run from a clone of the installer repo
A7TO8="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/tools/a7to8.py"
if [[ -f "$A7TO8" ]]; then
    install -m 755 "$A7TO8" "$HOME/.local/bin/astra8-a7to8"
    ok "astra8-a7to8 -> $HOME/.local/bin/astra8-a7to8"
fi
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
# >>> a7modules: SPIDER / STRAHL test runs
a7_test() {   # model exp t0 t1
    local log="$PREFIX/test_$1.log" cdf="$AWD/ncdf_out/$2$1.CDF"
    step "Test: $1 / $2, t=$3..$4 s"
    rm -f "$cdf"
    rm -rf "$AWD/strahl/$2_$1"
    ASTRA_PLATFORM=ubuntu AWD="$AWD" bash -c "cd '$AWD' && exe/as_exe -m $1 -v $2 -s $3 -e $4 -batch" >"$log" 2>&1 || true
    if grep -q "ASTRA normal exit" "$log" && [[ -f "$cdf" ]]; then
        ok "ASTRA normal exit; result: $cdf"
    else
        tail -n 30 "$log" >&2
        die "$1 run failed, see $log"
    fi
}
if [[ $DO_TEST -eq 1 && $HAVE_SPIDER -eq 1 ]]; then
    a7_test flux_spider aug34954 4 5
fi
if [[ $DO_TEST -eq 1 && $HAVE_SPIDER -eq 1 && $HAVE_STRAHL -eq 1 ]]; then
    a7_test astrahl_simple_msp aug34954 4 4.2
fi
# <<< a7modules

# >>> gacode: TGLF test run (64 TGLF worker processes; 1-5 min depending on the CPU)
if [[ $DO_TEST -eq 1 && $WITH_TGLF -eq 1 ]]; then
    step "Test: flux_tglf / aug34954, t=4..4.1 s (TGLF via xpr/tglfi, 2 calls)"
    tglf_log="$PREFIX/test_flux_tglf.log"
    rm -f "$AWD/ncdf_out/aug34954flux_tglf.CDF"
    ASTRA_PLATFORM=ubuntu AWD="$AWD" bash -c "cd '$AWD' && exe/as_exe -m flux_tglf -v aug34954 -s 4 -e 4.1 -batch" \
        >"$tglf_log" 2>&1 || true
    # "XPR wall time" is printed after each TGLF call; without it TGLF did not run
    if grep -q "ASTRA normal exit" "$tglf_log" && grep -q "XPR wall time" "$tglf_log" \
        && [[ -f "$AWD/ncdf_out/aug34954flux_tglf.CDF" ]]; then
        ok "ASTRA normal exit with TGLF; result: $AWD/ncdf_out/aug34954flux_tglf.CDF"
    else
        tail -n 30 "$tglf_log" >&2
        die "flux_tglf run failed, see $tglf_log"
    fi
fi
# <<< gacode

step "Done"
cat <<EOF

  ASTRA 8 is installed in $AWD

  Run with the graphic window:
      astra8 -m flux_feqis -v aug34954 -s 4 -e 5
  Run without window (results in $AWD/ncdf_out):
      astra8 -m flux_feqis -v aug34954 -s 4 -e 5 -batch

  Your models go to $AWD/equ, experiment files to $AWD/exp, U-files to $AWD/udb.
EOF
