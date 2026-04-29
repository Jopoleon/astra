#!/usr/bin/env python3
"""
eqdsk_to_TORIC_gs.py

Purpose
-------
Convert an EFIT EQDSK (gEQDSK) equilibrium into a TORIC "equi" file by:
  1) Reading ψ(R,Z) and 1D profiles (q, fpol, pprime, ffprim) from EQDSK.
  2) Choosing target flux surfaces using a normalized poloidal-flux radius rho_pol.
  3) Extracting for each surface a closed contour ψ(R,Z)=const that encloses the magnetic axis.
  4) Smoothing the contour (periodic Savitzky–Golay).
  5) Parameterizing the contour by a chosen poloidal angle θ (see --poloidal_angle):
       - arclength: θ = 2π s/L (uniform arclength)
       - geometric: θ = atan2(Z-Zax, R-Rax) (shifted so θ=0 at max R)
       - flux:      a "flux-like" proxy rectification of geometric θ using v~(dl/dθ)/R
  6) Computing Fourier modes for R(θ) and Z(θ).
  7) Optionally regularizing coefficients in radius near the axis (even/odd constraints).
  8) Writing TORIC.equi and optionally plotting reconstruction diagnostics.
  9) Optionally writing a smoothed EQDSK with:
       - monotonicity enforcement (arr >= 0),
       - iterative Grad-Shafranov SOR refinement,
       - GS residual diagnostic.

Usage examples
--------------
1) Basic conversion (default: flux theta, ψ smoothing, radial regularization):
     python eqdsk_to_TORIC_gs.py g12345.00001

2) Write a smoothed EQDSK:
     python eqdsk_to_TORIC_gs.py g12345.00001 --write_smooth_eqdsk

3) Write smoothed EQDSK + GS refinement:
     python eqdsk_to_TORIC_gs.py g12345.00001 --write_smooth_eqdsk --gs_refine

4) Full diagnostic: smoothed EQDSK + GS refinement + GS residual check + plot:
     python eqdsk_to_TORIC_gs.py g12345.00001 --write_smooth_eqdsk --gs_refine \\
         --gs_check --gs_check_plot --plot

5) Use geometric poloidal angle:
     python eqdsk_to_TORIC_gs.py g12345.00001 --poloidal_angle geometric

6) Disable ψ(R,Z) smoothing:
     python eqdsk_to_TORIC_gs.py g12345.00001 --no_smooth_psirz

7) Enable radial regularization and output on a grid starting at rho_pol=0:
     python eqdsk_to_TORIC_gs.py g12345.00001 --reg_smooth --nout 101 --rhoPF_max 0.99

8) GS refinement with tuned parameters:
     python eqdsk_to_TORIC_gs.py g12345.00001 --write_smooth_eqdsk \\
         --gs_refine --gs_n_iter 500 --gs_omega 1.0

Scientific note on "flux-like" theta
------------------------------------
The --poloidal_angle flux option is a proxy "rectification" based on a geometric
weight v ~ (dl/dθ)/R. This is not a full straight-field-line coordinate
(Boozer/Hamada) unless additional equilibrium field quantities are used. It is,
however, often smoother than raw geometric θ for Fourier representation.

Scientific note on GS refinement
---------------------------------
The --gs_refine option applies an iterative successive over-relaxation (SOR)
correction to the reconstructed psirz inside the plasma, driving it toward
consistency with the Grad-Shafranov equation:
    Δ*ψ = -μ₀ R² p'(ψ) - FF'(ψ)
where
    Δ*ψ = ∂²ψ/∂R² - (1/R) ∂ψ/∂R + ∂²ψ/∂Z²

The boundary (LCFS + outside) is kept fixed. This typically reduces the
GS residual from ~5-15% (pure interpolation) to <1-2%.

The SOR relaxation parameter ω controls convergence:
    ω < 1 : under-relaxation (safer, slower)
    ω = 1 : Gauss-Seidel
    ω > 1 : over-relaxation (faster convergence, may oscillate)
For nonlinear GS starting from a good initial guess, ω = 0.5–0.8 is safe.

Author
------
R. Bilato (IPP - Garching)
29-Mar-2026

"""

import re
import numpy as np
import argparse
import matplotlib.pyplot as plt

from skimage import measure
from scipy.interpolate import (LinearNDInterpolator, RegularGridInterpolator,
                               CubicSpline, interp1d)
from scipy.signal import savgol_filter
from scipy.ndimage import gaussian_filter, binary_erosion
from datetime import date

TWOPI = 2 * np.pi


# =============================================================================
# EQDSK parsing utilities
# =============================================================================
def parse_floats_fixed_or_regex(line: str, width: int = 16):
    """
    Parse one line of an EQDSK file into floats.

    EQDSK commonly uses fixed-width 16-character fields with Fortran 'D' exponents,
    but some files can be irregular.  We try fixed-width parsing first; if that
    fails, we fall back to a regex-based float extraction.

    Parameters
    ----------
    line : str
        A single line read from the EQDSK file.
    width : int
        Field width for fixed-width parsing (usually 16).

    Returns
    -------
    list[float]
        Floats parsed from the line.
    """
    s = line.rstrip("\n").replace("D", "E")
    out = []
    if len(s) >= width:
        ok = True
        for i in range(0, len(s), width):
            c = s[i:i + width].strip()
            if not c:
                continue
            try:
                out.append(float(c))
            except Exception:
                ok = False
                break
        if ok and out:
            return out

    pat = r"[+-]?(?:\d+\.\d*|\d*\.\d+|\d+)(?:[Ee][+-]?\d+)?"
    return [float(x) for x in re.findall(pat, s)]


def read_n_floats(f, n, width=16):
    """Read *n* floats from file handle *f*, spanning as many lines as needed."""
    v = []
    while len(v) < n:
        line = f.readline()
        if not line:
            raise EOFError("Unexpected EOF while reading floats")
        v.extend(parse_floats_fixed_or_regex(line, width))
    return np.array(v[:n], float)


def read_n_ints(f, n):
    """Read *n* integers from file handle *f*, spanning as many lines as needed."""
    v = []
    while len(v) < n:
        line = f.readline()
        if not line:
            raise EOFError("Unexpected EOF while reading ints")
        v.extend([int(x) for x in re.findall(r"[+-]?\d+", line)])
    return np.array(v[:n], int)


def read_geqdsk_full(path, width=16):
    """
    Read a gEQDSK file and return a dictionary of relevant arrays/metadata.

    Reads all standard EQDSK sections including:
      - 2D poloidal flux ψ(R,Z) on a rectangular grid
      - 1D profiles on a uniform ψ grid: q(ψ), fpol(ψ), pres(ψ), pprime(ψ), ffprim(ψ)
      - Plasma boundary (nbbbs, rbbbs, zbbbs)
      - Surrounding limiter contour (limitr, rlim, zlim)

    The EQDSK format stores some quantities redundantly (e.g. rmaxis appears
    in both line 2 and line 3).  We use the "duplicate" values from lines 3-4,
    which are often the final converged axis/LCFS values.

    Parameters
    ----------
    path : str
        Path to the gEQDSK file.
    width : int
        Fixed-field width for float parsing (typically 16).

    Returns
    -------
    dict
        Keys include nw, nh, R, Z, psirz, raxis, zaxis, simag, sibry,
        bcentr, current, fpol, pres, ffprim, pprime, qpsi, psi1d,
        nbbbs, rbbbs, zbbbs, limitr, rlim, zlim, and geometric scalars
        (rdim, zdim, rcentr, rleft, zmid).
    """
    with open(path, "r") as f:
        header = f.readline().rstrip("\n")
        ints = [int(x) for x in re.findall(r"[+-]?\d+", header)]
        nw, nh = ints[-2], ints[-1]

        # Line 1: grid dimensions and reference geometry
        rdim, zdim, rcentr, rleft, zmid = read_n_floats(f, 5, width)

        # Line 2: axis/boundary/field (primary values)
        rmaxis, zmaxis, simag, sibry, bcentr = read_n_floats(f, 5, width)

        # Line 3: current and duplicate axis/psi values
        line3 = read_n_floats(f, 5, width)
        current = line3[0]
        simag_2 = line3[1]     # duplicate simag (often more accurate)
        rmaxis_2 = line3[3]    # duplicate rmaxis

        # Line 4: duplicate zmaxis and sibry
        line4 = read_n_floats(f, 5, width)
        zmaxis_2 = line4[0]    # duplicate zmaxis
        sibry_2 = line4[2]     # duplicate sibry

        # Use the duplicate values (often the "final" axis/LCFS values)
        rmaxis = float(rmaxis_2)
        zmaxis = float(zmaxis_2)
        simag = float(simag_2)
        sibry = float(sibry_2)

        # 1D profiles on uniform ψ grid (nw points from simag to sibry)
        fpol = read_n_floats(f, nw, width)      # F = R·Bφ [m·T]
        pres = read_n_floats(f, nw, width)      # pressure [Pa]
        ffprim = read_n_floats(f, nw, width)    # FF' [(m·T)²/(Wb/rad)]
        pprime = read_n_floats(f, nw, width)    # p' [Pa/(Wb/rad)]

        # 2D poloidal flux ψ(R,Z) on the rectangular grid
        # Stored as psirz(nw*nh) in Z-major order: psirz[jZ, iR]
        psirz_flat = read_n_floats(f, nw * nh, width)
        psirz = psirz_flat.reshape((nh, nw))

        # Safety factor q(ψ) on the same uniform ψ grid
        qpsi = read_n_floats(f, nw, width)

        # Boundary and limiter point counts
        nbbbs, limitr = read_n_ints(f, 2)

        # Plasma boundary contour (R,Z pairs)
        rbbbs = zbbbs = None
        if nbbbs > 0:
            pairs = read_n_floats(f, 2 * nbbbs, width).reshape((nbbbs, 2))
            rbbbs, zbbbs = pairs[:, 0], pairs[:, 1]

        # Surrounding limiter contour (R,Z pairs)
        rlim = zlim = None
        if limitr > 0:
            pairs = read_n_floats(f, 2 * limitr, width).reshape((limitr, 2))
            rlim, zlim = pairs[:, 0], pairs[:, 1]

    # Build coordinate arrays
    R = np.linspace(rleft, rleft + rdim, nw)
    Z = np.linspace(zmid - zdim / 2, zmid + zdim / 2, nh)
    psi1d = np.linspace(simag, sibry, nw)

    print("\n=== EQDSK parameters ===")
    print('R of magn. axis [m]:                 {0:6.3f}'.format(rmaxis))
    print('Z of magn. axis [m]:                 {0:6.3f}'.format(zmaxis))
    print('Pol. flux at magn. axis [Weber/rad]: {0:6.3f}'.format(simag))
    print('Pol. flux at separatrix [Weber/rad]: {0:6.3f}'.format(sibry))
    print('R of vacuum tor. magn. field [m]:    {0:6.3f}'.format(rcentr))
    print('Vacuum toroidal magnetic field [T]:  {0:6.3f}'.format(bcentr))
    print('Plasma current [MA]:                 {0:6.3f}'.format(current / 1E6))
    print("=======================\n")

    return dict(
        nw=nw, nh=nh,
        R=R, Z=Z, psirz=psirz,
        rdim=rdim, zdim=zdim, rcentr=rcentr, rleft=rleft, zmid=zmid,
        raxis=rmaxis, zaxis=zmaxis,
        simag=simag, sibry=sibry,
        bcentr=bcentr,
        current=current,
        psi1d=psi1d,
        fpol=fpol, pres=pres, ffprim=ffprim, pprime=pprime, qpsi=qpsi,
        nbbbs=nbbbs, limitr=limitr,
        rbbbs=rbbbs, zbbbs=zbbbs,
        rlim=rlim, zlim=zlim
    )


def write_geqdsk(eq, psirz_out, out_path, width=16, ncol=5,
                 title="SMOOTHED_EQDSK"):
    """
    Write a gEQDSK file with fixed-width floats.

    Formatting:
      - floats: %16.9e (Fortran-style) with 'E' exponent
      - ncol floats per line (default 5, standard in EQDSK)
      - header ends with idum nw nh (required by many readers)

    Data layout written:
      scalars (4 lines × 5 floats) +
      fpol, pres, ffprim, pprime (nw each) +
      psirz (nw*nh, Z-major order; i.e. psirz[jZ, iR] flattened row-major) +
      qpsi (nw) +
      nbbbs, limitr +
      rbbbs/zbbbs pairs +
      rlim/zlim pairs (if present in eq dict)

    Notes:
      - EQDSK variants exist; this matches the structure commonly expected.
      - Ensure eq["nw"], eq["nh"], and psirz_out.shape == (nh, nw).

    Parameters
    ----------
    eq : dict
        EQDSK dictionary with all required scalars and arrays.
    psirz_out : 2D ndarray, shape (nh, nw)
        Poloidal flux to write.
    out_path : str
        Output file path.
    width : int
        Fixed-width field width for floats.
    ncol : int
        Number of floats per line.
    title : str
        Header title string (truncated to 48 chars).
    """
    nw = int(eq["nw"])
    nh = int(eq["nh"])

    ps = np.asarray(psirz_out, float)
    if ps.shape != (nh, nw):
        raise ValueError(
            f"psirz_out must have shape (nh,nw)=({nh},{nw}), got {ps.shape}")

    # Fixed-width float formatting (E exponent)
    # Python's %e uses lowercase 'e'; we force uppercase 'E' for Fortran compat.
    fmt = f"%{width}.9e"

    def ffloat(x):
        return (fmt % float(x)).replace("e", "E")

    def write_floats(f, arr):
        arr = np.asarray(arr, float).ravel()
        for i in range(0, arr.size, ncol):
            chunk = arr[i:i + ncol]
            f.write("".join(ffloat(v) for v in chunk) + "\n")

    def write_ints_line(f, a, b):
        # EFIT commonly uses 5-char integers for boundary/limiter counts
        f.write(f"{int(a):5d}{int(b):5d}\n")

    # Extract scalars from the dictionary
    rdim = float(eq["rdim"])
    zdim = float(eq["zdim"])
    rcentr = float(eq["rcentr"])
    rleft = float(eq["rleft"])
    zmid = float(eq["zmid"])
    raxis = float(eq["raxis"])
    zaxis = float(eq["zaxis"])
    simag = float(eq["simag"])
    sibry = float(eq["sibry"])
    bcentr = float(eq["bcentr"])
    current = float(eq["current"])

    # Extract 1D profile arrays
    fpol = np.asarray(eq["fpol"], float)
    pres = np.asarray(eq["pres"], float)
    ffprim = np.asarray(eq["ffprim"], float)
    pprime = np.asarray(eq["pprime"], float)
    qpsi = np.asarray(eq["qpsi"], float)

    if any(a.size != nw for a in [fpol, pres, ffprim, pprime, qpsi]):
        raise ValueError("One or more 1D arrays do not have length nw.")

    # Boundary points (may be None if absent)
    rbbbs = eq.get("rbbbs", None)
    zbbbs = eq.get("zbbbs", None)
    if rbbbs is None or zbbbs is None:
        nbbbs = 0
    else:
        nbbbs = min(len(rbbbs), len(zbbbs))

    # Limiter points (may be None if absent)
    rlim = eq.get("rlim", None)
    zlim = eq.get("zlim", None)
    if rlim is None or zlim is None:
        limitr = 0
    else:
        limitr = min(len(rlim), len(zlim))

    # Header: many tools parse last two ints as nw nh
    idum = 0
    title48 = (title[:48]).ljust(48)
    with open(out_path, "w") as f:
        f.write(f"{title48}{idum:5d}{nw:5d}{nh:5d}\n")

        # Four 5-float lines (EFIT standard layout)
        # Line 1: grid dimensions and reference geometry
        write_floats(f, [rdim, zdim, rcentr, rleft, zmid])
        # Line 2: axis position, flux values, toroidal field
        write_floats(f, [raxis, zaxis, simag, sibry, bcentr])
        # Line 3: current, duplicate simag, dummy, duplicate rmaxis, dummy
        write_floats(f, [current, simag, 0.0, raxis, 0.0])
        # Line 4: duplicate zmaxis, dummy, duplicate sibry, dummy, dummy
        write_floats(f, [zaxis, 0.0, sibry, 0.0, 0.0])

        # 1D profiles (nw values each)
        write_floats(f, fpol)
        write_floats(f, pres)
        write_floats(f, ffprim)
        write_floats(f, pprime)

        # 2D poloidal flux, Z-major order: flatten row-major [jZ, iR]
        write_floats(f, ps.ravel(order="C"))

        # Safety factor profile
        write_floats(f, qpsi)

        # Boundary / limiter counts
        write_ints_line(f, nbbbs, limitr)

        # Boundary points (R,Z pairs interleaved)
        if nbbbs > 0:
            pairs = np.column_stack(
                [np.asarray(rbbbs, float),
                 np.asarray(zbbbs, float)]).ravel(order="C")
            write_floats(f, pairs)

        # Limiter points (R,Z pairs interleaved)
        if limitr > 0:
            pairs = np.column_stack(
                [np.asarray(rlim, float),
                 np.asarray(zlim, float)]).ravel(order="C")
            write_floats(f, pairs)


# =============================================================================
# Axis-preserving psirz smoothing
# =============================================================================
def smooth_psirz_axis_preserving(R, Z, psirz, raxis, zaxis,
                                 sigma_R=1.0, sigma_Z=1.0,
                                 r0=0.10, dr=0.05,
                                 preserve_value=True):
    """
    Smooth ψ(R,Z) with a Gaussian filter, but preserve ψ near the magnetic axis.

    We blend between original and smoothed ψ using a smoothstep function:
      - near axis (rr < r0): mostly original ψ
      - far from axis (rr > r0+dr): fully smoothed ψ
      - transition region: smoothstep blend w = 3t² - 2t³

    Parameters
    ----------
    R, Z : 1D arrays
        Coordinate arrays defining the ψ(R,Z) grid.
    psirz : 2D array (nh, nw)
        Poloidal flux on the grid.
    raxis, zaxis : float
        Magnetic axis position.
    sigma_R, sigma_Z : float
        Gaussian smoothing sigma in R and Z index directions (grid units, not meters).
    r0 : float
        Inner radius (meters) of the axis-preserving blend region.
    dr : float
        Transition width (meters) from original to smoothed.
    preserve_value : bool
        If True, also preserve the exact ψ value at the nearest grid node to the axis.

    Returns
    -------
    2D array (nh, nw)
        Blended (axis-preserving) smoothed ψ field.
    """
    ps = np.asarray(psirz, float)
    ps_s = gaussian_filter(ps, sigma=(sigma_Z, sigma_R), mode="nearest")

    # Distance from axis at each grid cell (in meters)
    ZZ, RR = np.meshgrid(Z, R, indexing="ij")
    rr = np.sqrt((RR - raxis)**2 + (ZZ - zaxis)**2)

    # Smoothstep weight in [0,1]
    t = (rr - r0) / (dr + 1e-300)
    t = np.clip(t, 0.0, 1.0)
    w = t * t * (3.0 - 2.0 * t)

    ps_new = (1.0 - w) * ps + w * ps_s

    if preserve_value:
        # Preserve the closest grid node to the axis location
        iZ = int(np.argmin(np.abs(Z - zaxis)))
        iR = int(np.argmin(np.abs(R - raxis)))
        ps_new[iZ, iR] = ps[iZ, iR]

    return ps_new


# =============================================================================
# Contour extraction
# =============================================================================
def point_in_polygon(x, y, polyx, polyy):
    """
    Ray casting point-in-polygon test.

    Returns True if the point (x, y) lies inside the polygon defined
    by vertices (polyx, polyy).
    """
    inside = False
    j = len(polyx) - 1
    for i in range(len(polyx)):
        xi, yi = polyx[i], polyy[i]
        xj, yj = polyx[j], polyy[j]
        if ((yi > y) != (yj > y)) and \
           (x < (xj - xi) * (y - yi) / (yj - yi + 1e-300) + xi):
            inside = not inside
        j = i
    return inside


def select_axis_enclosing_contour(Rgrid, Zgrid, psi2d, level, raxis, zaxis,
                                  min_points=8):
    """
    Extract all contours ψ=level and choose the one that encloses the
    magnetic axis and has the largest area.

    skimage.measure.find_contours returns curves in index coordinates;
    we map them to (R,Z) physical coordinates.

    Parameters
    ----------
    Rgrid, Zgrid : 1D arrays
        Coordinate arrays for the psi2d grid.
    psi2d : 2D array
        Poloidal flux field.
    level : float
        Contour level ψ = level.
    raxis, zaxis : float
        Magnetic axis position.
    min_points : int
        Minimum number of points for a contour to be considered valid.

    Returns
    -------
    (R, Z) tuple of 1D arrays, or None if no valid contour found.
    """
    cs = measure.find_contours(psi2d, level)
    if not cs:
        return None

    best = None
    best_area = -1.0
    for c in cs:
        if c.shape[0] < min_points:
            continue

        # Map from index coordinates to physical (R,Z)
        z_idx = c[:, 0]
        r_idx = c[:, 1]
        R = np.interp(r_idx, np.arange(Rgrid.size), Rgrid)
        Z = np.interp(z_idx, np.arange(Zgrid.size), Zgrid)

        # Ensure closed contour
        if np.hypot(R[0] - R[-1], Z[0] - Z[-1]) > 1e-6:
            R = np.r_[R, R[0]]
            Z = np.r_[Z, Z[0]]

        # Axis must lie inside the contour
        if not point_in_polygon(raxis, zaxis, R, Z):
            continue

        # Choose the contour with largest area (robust for multiple candidates)
        area = 0.5 * np.abs(
            np.dot(R, np.roll(Z, -1)) - np.dot(Z, np.roll(R, -1)))
        if area > best_area:
            best_area = area
            best = (R, Z)
    return best


def upsample_psi(R, Z, psirz, upsample):
    """
    Upsample ψ(R,Z) onto a finer grid using RegularGridInterpolator.

    Upsampling helps contour extraction accuracy and reduces pixelization
    artifacts that produce jagged contours.

    Parameters
    ----------
    R, Z : 1D arrays
        Original coordinate arrays.
    psirz : 2D array (nh, nw)
        Original poloidal flux.
    upsample : int
        Upsampling factor (e.g. 8 means 8× finer grid in each direction).

    Returns
    -------
    Rf, Zf : 1D arrays
        Fine coordinate arrays.
    psif : 2D array
        Upsampled poloidal flux.
    """
    if upsample <= 1:
        return R, Z, psirz

    nwf = int(R.size * upsample)
    nhf = int(Z.size * upsample)
    Rf = np.linspace(R[0], R[-1], nwf)
    Zf = np.linspace(Z[0], Z[-1], nhf)

    itp = RegularGridInterpolator((Z, R), psirz,
                                  bounds_error=False, fill_value=np.nan)
    ZZ, RR = np.meshgrid(Zf, Rf, indexing="ij")
    pts = np.column_stack([ZZ.ravel(), RR.ravel()])
    psif = itp(pts).reshape((nhf, nwf))

    return Rf, Zf, psif


# =============================================================================
# Curve smoothing + arclength resampling
# =============================================================================
def smooth_closed_curve_periodic_savgol(Rc, Zc, window=21, polyorder=3):
    """
    Periodic Savitzky–Golay smoothing of a closed planar curve.

    Strategy:
      - enforce closure (append first point if not closed),
      - drop the duplicated endpoint,
      - periodic extension by half-window on both sides,
      - apply SavGol filter on the extended data,
      - trim back to original length and re-close.

    Parameters
    ----------
    Rc, Zc : 1D arrays
        Coordinates of the closed curve.
    window : int
        Savitzky–Golay window length (must be odd; adjusted if even).
    polyorder : int
        Polynomial order for SavGol filter (must be < window).

    Returns
    -------
    Rs, Zs : 1D arrays
        Smoothed closed curve coordinates (with closure point appended).
    """
    Rc = np.asarray(Rc, float)
    Zc = np.asarray(Zc, float)

    # Ensure the curve is closed
    if np.hypot(Rc[0] - Rc[-1], Zc[0] - Zc[-1]) > 1e-12:
        Rc = np.r_[Rc, Rc[0]]
        Zc = np.r_[Zc, Zc[0]]

    # Drop the duplicated endpoint for periodic processing
    Rc = Rc[:-1]
    Zc = Zc[:-1]
    n = Rc.size

    # Ensure window is odd
    if window % 2 == 0:
        window += 1
    window = min(window, n - (1 - n % 2))
    if window <= polyorder + 2:
        return np.r_[Rc, Rc[0]], np.r_[Zc, Zc[0]]

    # Periodic extension
    h = window // 2
    Rp = np.r_[Rc[-h:], Rc, Rc[:h]]
    Zp = np.r_[Zc[-h:], Zc, Zc[:h]]

    Rs = savgol_filter(Rp, window_length=window, polyorder=polyorder,
                       mode="interp")
    Zs = savgol_filter(Zp, window_length=window, polyorder=polyorder,
                       mode="interp")

    # Trim back
    Rs = Rs[h:h + n]
    Zs = Zs[h:h + n]

    # Re-close
    return np.r_[Rs, Rs[0]], np.r_[Zs, Zs[0]]


def ensure_closed(Rc, Zc):
    """Ensure closure by appending the first point if needed."""
    Rc = np.asarray(Rc)
    Zc = np.asarray(Zc)
    if np.hypot(Rc[0] - Rc[-1], Zc[0] - Zc[-1]) > 1e-12:
        Rc = np.r_[Rc, Rc[0]]
        Zc = np.r_[Zc, Zc[0]]
    return Rc, Zc


def arclength_param(Rc, Zc):
    """Return cumulative arclength s and total length L for a polyline."""
    dR = np.diff(Rc)
    dZ = np.diff(Zc)
    ds = np.sqrt(dR * dR + dZ * dZ)
    s = np.r_[0.0, np.cumsum(ds)]
    return s, s[-1]


def resample_closed_curve_uniform_arclength(Rc, Zc, npts):
    """
    Resample a closed curve with npts equally spaced in arclength.

    Uses periodic CubicSpline(s) → (R,Z), then evaluates at uniform s in [0, L).

    Parameters
    ----------
    Rc, Zc : 1D arrays
        Closed curve coordinates.
    npts : int
        Number of output sample points.

    Returns
    -------
    su : 1D array
        Uniform arclength samples.
    L : float
        Total arclength.
    Ru, Zu : 1D arrays
        Resampled curve coordinates.
    """
    Rc, Zc = ensure_closed(Rc, Zc)
    s, L = arclength_param(Rc, Zc)

    splR = CubicSpline(s, Rc, bc_type="periodic")
    splZ = CubicSpline(s, Zc, bc_type="periodic")

    su = np.linspace(0.0, L, npts, endpoint=False)
    Ru = splR(su)
    Zu = splZ(su)
    return su, L, Ru, Zu


def phase_lock_theta0_at_maxR(Ru):
    """Choose θ=0 at maximum R to fix an otherwise arbitrary phase rotation."""
    return int(np.argmax(Ru))


def phase_lock_theta0_at_outboard_midplane(Ru, Zu, raxis, zaxis):
    """
    Return index j0 such that theta=0 is at the outboard midplane crossing:
      Z ≈ Zaxis and R > Raxis.
    If no outboard points exist (rare), fall back to max-R.
    """
    Ru = np.asarray(Ru, float)
    Zu = np.asarray(Zu, float)

    mask = Ru > raxis
    if not np.any(mask):
        return int(np.argmax(Ru))

    dz = np.abs(Zu - zaxis)
    dz_masked = np.where(mask, dz, np.inf)
    dzmin = np.min(dz_masked)

    # Candidates closest to midplane on the outboard side
    cand = np.where(dz_masked <= dzmin + 1e-12)[0]
    # Pick the one with largest R (outboard crossing)
    j0 = cand[np.argmax(Ru[cand])]
    return int(j0)


# =============================================================================
# Geometric/flux theta helpers
# =============================================================================
def make_strictly_increasing_theta(th, R, Z, eps=1e-10):
    """
    Make a theta array strictly increasing for spline construction.

    CubicSpline requires strictly increasing x values.  For geometric angles
    obtained from atan2, it is common to have repeated/near-repeated angles
    due to discretization.

    Strategy:
    - sort by theta,
    - merge points whose theta differ by <= eps (average their R and Z),
    - return reduced arrays with strictly increasing theta.

    Parameters
    ----------
    th, R, Z : 1D arrays
        Theta values and corresponding (R,Z) coordinates.
    eps : float
        Merge tolerance in radians.

    Returns
    -------
    th_new, R_new, Z_new : 1D arrays
        Deduplicated, strictly increasing arrays.
    """
    th = np.asarray(th, float)
    R = np.asarray(R, float)
    Z = np.asarray(Z, float)

    idx = np.argsort(th)
    th = th[idx]; R = R[idx]; Z = Z[idx]

    th_new = [th[0]]
    R_new = [R[0]]
    Z_new = [Z[0]]
    cnt = [1]

    for k in range(1, th.size):
        if th[k] - th_new[-1] <= eps:
            R_new[-1] += R[k]
            Z_new[-1] += Z[k]
            cnt[-1] += 1
        else:
            th_new.append(th[k])
            R_new.append(R[k])
            Z_new.append(Z[k])
            cnt.append(1)

    th_new = np.array(th_new)
    R_new = np.array(R_new) / np.array(cnt)
    Z_new = np.array(Z_new) / np.array(cnt)

    if th_new.size < 8:
        raise ValueError("Too few unique theta points after deduplication")

    return th_new, R_new, Z_new


def rectify_to_flux_theta(theta_bar_u, R_u, dR_dth_u, dZ_dth_u):
    """
    Construct a flux-like rectified angle theta_flux from a base angle theta_bar.

    Uses the proxy weight:
        v(θ_bar) ~ (dl/dθ_bar) / R
    with:
        dl/dθ_bar = sqrt( (dR/dθ)² + (dZ/dθ)² )

    The rectified angle is defined by:
        d(θ_flux)/d(θ_bar) = v / <v>
    so that θ_flux runs from 0 to 2π over one poloidal turn.

    Parameters
    ----------
    theta_bar_u : 1D array
        Uniform base angle samples.
    R_u : 1D array
        R coordinate at each base angle sample.
    dR_dth_u, dZ_dth_u : 1D arrays
        Derivatives dR/dθ and dZ/dθ at each sample.

    Returns
    -------
    theta_flux : 1D array
        Rectified flux-like angle in [0, 2π).
    """
    dl_dth = np.sqrt(dR_dth_u**2 + dZ_dth_u**2) + 1e-300
    v = dl_dth / np.maximum(R_u, 1e-300)
    v_avg = np.mean(v)

    dtheta_bar = theta_bar_u[1] - theta_bar_u[0]
    dth_flux = (v / v_avg) * dtheta_bar

    theta_flux = np.cumsum(dth_flux)
    theta_flux -= theta_flux[0]
    theta_flux *= (TWOPI / theta_flux[-1])
    return theta_flux


def sample_surface_uniform_theta(Rc, Zc, *,
                                 raxis, zaxis,
                                 ntheta,
                                 poloidal_angle="arclength",
                                 theta0=-1.0,
                                 theta_eps=1e-10):
    """
    Sample a single flux surface on a uniform theta grid using one of
    three poloidal angles: arclength, geometric, or flux.

    Parameters
    ----------
    Rc, Zc : 1D arrays
        Closed curve coordinates of the flux surface.
    raxis, zaxis : float
        Magnetic axis position.
    ntheta : int
        Number of uniform θ samples.
    poloidal_angle : str
        One of "arclength", "geometric", "flux".
    theta0 : float
        If > 0, θ=0 at outboard midplane; if < 0, θ=0 at max R.
    theta_eps : float
        Deduplication tolerance for geometric/flux θ.

    Returns
    -------
    theta_u : (ntheta,) array
        Uniform theta samples in [0, 2π).
    Ru, Zu : (ntheta,) arrays
        Surface coordinates at those theta samples.

    Notes
    -----
    - arclength: directly resamples curve by arclength (periodic spline in s).
    - geometric: builds periodic spline of R(θ), Z(θ) where θ = atan2(Z-Zax, R-Rax).
    - flux: starts from geometric θ_bar, rectifies to theta_flux using v~(dl/dθ)/R,
            then samples uniformly in theta_flux.
    """
    theta_u = np.linspace(0.0, TWOPI, ntheta, endpoint=False)

    if poloidal_angle == "arclength":
        # Uniform arclength ⇒ uniform theta by construction.
        _, _, Ru, Zu = resample_closed_curve_uniform_arclength(Rc, Zc, ntheta)

        # Phase-lock θ=0
        if theta0 < 0.0:
            j0 = phase_lock_theta0_at_maxR(Ru)
        else:
            j0 = phase_lock_theta0_at_outboard_midplane(Ru, Zu, raxis, zaxis)
        Ru = np.roll(Ru, -j0)
        Zu = np.roll(Zu, -j0)
        return theta_u, Ru, Zu

    # -----------------------------------------------------------------
    # For geometric/flux: build periodic splines in geometric theta first.
    # -----------------------------------------------------------------
    # Compute geometric poloidal angle about the magnetic axis
    th = np.unwrap(np.arctan2(Zc - zaxis, Rc - raxis))

    # Shift so theta=0 at the chosen reference point (phase locking)
    if theta0 < 0.0:
        j0 = int(np.argmax(Rc))
    else:
        j0 = phase_lock_theta0_at_outboard_midplane(Rc, Zc, raxis, zaxis)
    th0 = th[j0]
    th = (th - th0) % TWOPI

    # Ensure strict monotonicity for spline construction
    ths, Rs, Zs = make_strictly_increasing_theta(th, Rc, Zc, eps=theta_eps)

    # Ensure we include endpoints at 0 and 2π for periodic splines
    if ths[0] > 0.0:
        ths = np.r_[0.0, ths]
        Rs = np.r_[Rs[-1], Rs]
        Zs = np.r_[Zs[-1], Zs]
    if ths[-1] < TWOPI:
        ths = np.r_[ths, TWOPI]
        Rs = np.r_[Rs, Rs[0]]
        Zs = np.r_[Zs, Zs[0]]

    # Deduplicate again in case 0 or 2π were already present
    ths, Rs, Zs = make_strictly_increasing_theta(ths, Rs, Zs, eps=theta_eps)

    # Enforce periodic point equality required by CubicSpline(..., bc_type="periodic")
    Rs[-1] = Rs[0]
    Zs[-1] = Zs[0]

    splR = CubicSpline(ths, Rs, bc_type="periodic")
    splZ = CubicSpline(ths, Zs, bc_type="periodic")

    if poloidal_angle == "geometric":
        return theta_u, splR(theta_u), splZ(theta_u)

    if poloidal_angle == "flux":
        # Build mapping theta_flux(theta_bar) on a uniform base grid
        theta_bar_u = theta_u
        Rb = splR(theta_bar_u)
        dR = splR(theta_bar_u, 1)
        dZ = splZ(theta_bar_u, 1)

        theta_flux = rectify_to_flux_theta(theta_bar_u, Rb, dR, dZ)

        # Ensure monotonicity and well-conditioned inversion
        theta_flux = np.maximum.accumulate(theta_flux)
        theta_flux[0] = 0.0
        theta_flux[-1] = TWOPI

        # Invert: theta_flux(theta_bar) → theta_bar(theta_flux)
        inv = CubicSpline(theta_flux, theta_bar_u, bc_type="natural")
        theta_bar_of_flux = inv(theta_u) % TWOPI

        return theta_u, splR(theta_bar_of_flux), splZ(theta_bar_of_flux)

    raise ValueError(
        "Unknown poloidal_angle (must be 'arclength', 'geometric', or 'flux').")


# =============================================================================
# Fourier decomposition
# =============================================================================
def fourier_coeffs(y, mmax):
    """
    Compute Fourier coefficients of y(θ) from uniform θ samples via FFT.

      y(θ) = a0 + Σ_{m≥1} [ a[m] cos(mθ) + b[m] sin(mθ) ]

    Parameters
    ----------
    y : 1D array
        Function values at uniform θ in [0, 2π).
    mmax : int
        Maximum mode number.

    Returns
    -------
    a0 : float
        Mean value.
    a, b : 1D arrays of length mmax+1
        Cosine and sine coefficients.
    """
    n = y.size
    Y = np.fft.rfft(y)
    a0 = Y[0].real / n
    a = np.zeros(mmax + 1)
    b = np.zeros(mmax + 1)
    for m in range(1, mmax + 1):
        if m < Y.size:
            a[m] = 2 * Y[m].real / n
            b[m] = -2 * Y[m].imag / n
    return a0, a, b


def reconstruct(theta, a0, a, b):
    """Reconstruct y(θ) from Fourier coefficients."""
    out = np.full_like(theta, a0, dtype=float)
    for m in range(1, len(a)):
        out += a[m] * np.cos(m * theta) + b[m] * np.sin(m * theta)
    return out


def condensation_mmax(theta_u, R, Z, mmax_max, tol_rms):
    """
    Choose the smallest mmax such that both R(θ) and Z(θ) reconstructions
    achieve RMS geometric error ≤ tol_rms [meters].

    This "condensation" avoids keeping unnecessary high-m modes that fit
    noise rather than shape.

    Parameters
    ----------
    theta_u : 1D array
        Uniform θ samples.
    R, Z : 1D arrays
        Surface coordinates.
    mmax_max : int
        Maximum mode number to try.
    tol_rms : float
        RMS tolerance in meters.

    Returns
    -------
    mmax : int
        Chosen mode number.
    pack : tuple
        (a0R, aR, bR, a0Z, aZ, bZ) full coefficient arrays up to mmax_max.
    """
    a0R, aR_full, bR_full = fourier_coeffs(R, mmax_max)
    a0Z, aZ_full, bZ_full = fourier_coeffs(Z, mmax_max)

    for mmax in range(1, mmax_max + 1):
        Rr = reconstruct(theta_u, a0R, aR_full[:mmax + 1],
                         bR_full[:mmax + 1])
        Zr = reconstruct(theta_u, a0Z, aZ_full[:mmax + 1],
                         bZ_full[:mmax + 1])
        err = np.sqrt(np.mean((R - Rr)**2 + (Z - Zr)**2))
        if err <= tol_rms:
            return mmax, (a0R, aR_full, bR_full, a0Z, aZ_full, bZ_full)

    return mmax_max, (a0R, aR_full, bR_full, a0Z, aZ_full, bZ_full)


# =============================================================================
# Near-axis radial regularization (even/odd parity constraints)
# =============================================================================
def ensureSmoothness_even(x, y, scale, xOut):
    """
    Enforce even near-axis behavior for a radial function y(x):
        y(x) ≈ c0 + c2 x²   for x ≤ scale

    For x > scale, allow independent DOFs at each original knot.

    Implementation:
    - Build a basis matrix A with columns [1, x², e_j for j where x_j > scale].
    - Solve least squares A c ≈ y.
    - Evaluate the same basis on xOut and compute yOut = Aout c.

    This ensures that the Fourier coefficients of even-m modes behave as
    c0 + c2·ρ² near the axis, which is the correct MHD parity.
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    xOut = np.asarray(xOut, float)

    cols = [np.ones_like(x), x**2]
    idx_free = np.where(x > scale)[0]
    for j in idx_free:
        col = np.zeros_like(x)
        col[j] = 1.0
        cols.append(col)

    A = np.vstack(cols).T
    coeffs, *_ = np.linalg.lstsq(A, y, rcond=None)

    cols_out = [np.ones_like(xOut), xOut**2]
    for k in range(2, A.shape[1]):
        f = interp1d(x, A[:, k], fill_value="extrapolate")
        cols_out.append(f(xOut))
    Aout = np.vstack(cols_out).T
    return Aout @ coeffs


def ensureSmoothness_odd(x, y, scale, xOut):
    """
    Enforce odd near-axis behavior for a radial function y(x):
        y(x) ≈ c1 x   for x ≤ scale

    For x > scale, allow independent DOFs at each original knot.

    This ensures that the Fourier coefficients of odd-m modes vanish
    linearly at the axis, which is the correct MHD parity.
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    xOut = np.asarray(xOut, float)

    cols = [x]
    idx_free = np.where(x > scale)[0]
    for j in idx_free:
        col = np.zeros_like(x)
        col[j] = 1.0
        cols.append(col)

    A = np.vstack(cols).T
    coeffs, *_ = np.linalg.lstsq(A, y, rcond=None)

    cols_out = []
    for k in range(A.shape[1]):
        f = interp1d(x, A[:, k], fill_value="extrapolate")
        cols_out.append(f(xOut))
    Aout = np.vstack(cols_out).T
    return Aout @ coeffs


# =============================================================================
# rho_tor interpolator (toroidal flux radius)
# =============================================================================
def build_rho_tor_interpolator(psi1d, qpsi):
    """
    Build an interpolator rho_tor(psi) from q(psi) by computing toroidal flux:
        Φ(ψ) = ∫ q(ψ) dψ   (up to a constant)

    Then normalize:
        Φ_N = (Φ - Φ_axis) / (Φ_lcfs - Φ_axis)

    And define:
        rho_tor = sqrt(Φ_N)

    Uses composite trapezoidal integration on the provided ψ grid.

    Parameters
    ----------
    psi1d : 1D array
        Uniform ψ grid from axis to LCFS.
    qpsi : 1D array
        Safety factor on the same grid.

    Returns
    -------
    interp1d object
        Interpolator rho_tor(psi).
    """
    psi1d = np.asarray(psi1d, float)
    qpsi = np.asarray(qpsi, float)

    dpsi = np.diff(psi1d)
    Phi = np.zeros_like(psi1d)
    Phi[1:] = np.cumsum(0.5 * (qpsi[1:] + qpsi[:-1]) * dpsi)

    PhiN = (Phi - Phi[0]) / (Phi[-1] - Phi[0] + 1e-300)
    rho_tor = np.sqrt(np.clip(PhiN, 0.0, 1.0))

    return interp1d(psi1d, rho_tor, kind="linear",
                    bounds_error=False,
                    fill_value=(rho_tor[0], rho_tor[-1]))


def psi2d_to_psi1d(psi2d, psi_axis_2d, psi_lcfs_2d, simag_1d, sibry_1d):
    """
    Map a ψ value from the 2D ψ(R,Z) field to the 1D ψ grid used by profiles.

    We assume a linear mapping between (axis, LCFS) endpoints:
        t = (ψ_2d - ψ_axis_2d) / (ψ_lcfs_2d - ψ_axis_2d)
        ψ_1d = simag + t · (sibry - simag)

    This is needed because the 2D field may have slightly different axis/LCFS
    values than the 1D profile grid.
    """
    t = (psi2d - psi_axis_2d) / (psi_lcfs_2d - psi_axis_2d + 1e-300)
    return simag_1d + t * (sibry_1d - simag_1d)


# =============================================================================
# COCOS inference (heuristic diagnostic)
# =============================================================================
def _robust_sign(x, frac=0.1):
    """
    Robust sign of an array: ignore values close to 0 by trimming
    to central quantiles.  Returns -1, 0, or +1.
    """
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return 0
    lo = np.quantile(x, frac)
    hi = np.quantile(x, 1.0 - frac)
    y = x[(x >= lo) & (x <= hi)]
    if y.size == 0:
        y = x
    m = np.median(y)
    if np.abs(m) < 1e-30:
        return 0
    return 1 if m > 0 else -1


def infer_cocos_from_eqdsk(eq, psi_axis_2d, psi_lcfs_2d):
    """
    Attempt to infer the COCOS convention from EQDSK content.

    IMPORTANT: This is a *diagnostic/inference* tool.  EQDSK typically does not
    encode all information needed to uniquely determine the COCOS (notably the
    2π normalization and some sign choices).  Therefore, the routine returns a
    *set of candidates* and a detailed diagnostic report.

    Uses:
    - sign(Ip) from eq["current"]
    - sign(B0) from eq["bcentr"]
    - sign of Δψ = ψ_lcfs - ψ_axis
    - sign of q(ψ), F(ψ)=RBφ, p'(ψ)

    Returns
    -------
    report : dict
        Diagnostics (signs, basic consistency checks).
    candidates : list[int]
        Candidate COCOS indices among {1..8, 11..18}.
    """
    Ip = float(eq.get("current", np.nan))
    B0 = float(eq.get("bcentr", np.nan))

    s_Ip = 0 if not np.isfinite(Ip) else (1 if Ip > 0 else -1)
    s_B0 = 0 if not np.isfinite(B0) else (1 if B0 > 0 else -1)

    dpsi = float(psi_lcfs_2d - psi_axis_2d)
    s_dpsi = 0 if abs(dpsi) < 1e-30 else (1 if dpsi > 0 else -1)

    s_q = _robust_sign(eq.get("qpsi", np.array([])))
    s_F = _robust_sign(eq.get("fpol", np.array([])))
    s_pprime = _robust_sign(eq.get("pprime", np.array([])))

    # Heuristic: if ψ increases outward (Δψ>0), typical p'(ψ) < 0
    pprime_expected = -s_dpsi if s_dpsi != 0 else 0
    pprime_ok = (s_pprime == 0 or pprime_expected == 0
                 or s_pprime == pprime_expected)
    q_warn_abs = (s_q < 0)

    report = dict(
        Ip=Ip, B0=B0,
        psi_axis_2d=float(psi_axis_2d),
        psi_lcfs_2d=float(psi_lcfs_2d),
        dpsi=dpsi,
        sign_Ip=s_Ip, sign_B0=s_B0, sign_dpsi=s_dpsi,
        sign_q=s_q, sign_F=s_F, sign_pprime=s_pprime,
        pprime_expected_sign=pprime_expected,
        pprime_sign_consistent=pprime_ok,
        q_negative_warning=q_warn_abs,
        notes=[]
    )

    if not pprime_ok:
        report["notes"].append(
            "pprime sign looks inconsistent with dpsi sign (heuristic check). "
            "Either EQDSK uses different ψ sign convention, or profiles are "
            "not consistent/noisy.")
    if q_warn_abs:
        report["notes"].append(
            "Median(q) < 0. Many producers store |q|, but if q is truly "
            "negative this constrains COCOS.")

    # Candidate selection (see docstring for reasoning)
    if s_q < 0:
        base_family = [5, 6, 7, 8]
    else:
        base_family = ([1, 2, 3, 4] if s_q > 0
                       else [1, 2, 3, 4, 5, 6, 7, 8])
    candidates = base_family + [k + 10 for k in base_family]

    return report, candidates


# =============================================================================
# TORIC writer
# =============================================================================
def writef(file, fmt, n_per_line, array):
    """Write a 1D array with fixed formatting, n_per_line values per line."""
    arr = np.asarray(array).ravel()
    for i in range(0, len(arr), n_per_line):
        chunk = arr[i:i + n_per_line]
        line = "".join([fmt for _ in range(len(chunk))]) % tuple(chunk)
        file.write(line + "\n")


def write_toric_equi(filename,
                     rtorus, raxis, baxis, ip_A,
                     nmod, nmhd,
                     rho_pol, rho_tor,
                     rc, rs, zc, zs,
                     qsf, jphi_avg, G_covBphi):
    """
    Write TORIC equilibrium file in the expected format.

    The file contains:
      - Scalar parameters (Rtorus, Raxis, Baxis, Ip, nmod, nmhd)
      - Radial grids (rho_pol, rho_tor)
      - Fourier mode coefficients for R(θ) and Z(θ) on each surface
      - Safety factor, current profile, and covariant Bφ profiles
    """
    fmtr = "%18.9e"
    fmti = "%5d"
    with open(filename, "w") as f:
        f.write("Rtorus\n" + fmtr % rtorus + "\n")
        f.write("Raxis\n" + fmtr % rc[0, 0] + "\n")
        f.write("Baxis\n" + fmtr % baxis + "\n")
        f.write("Ip\n" + fmtr % ip_A + "\n")
        f.write("nmod\n" + fmti % (nmod - 1) + "\n")
        f.write("nmhd\n" + fmti % nmhd + "\n")

        f.write("rho poloidal \n"); writef(f, fmtr, 4, rho_pol)
        f.write("rho toroidal \n"); writef(f, fmtr, 4, np.abs(rho_tor))
        f.write("rho at last surface \n")
        writef(f, fmtr, 4, [abs(rho_pol[-1]), abs(rho_tor[-1])])

        # Modes: m=0 first, then m=1..nmod-1
        f.write("Modes\n")
        writef(f, fmtr, 4, rc[:, 0])    # R cosine m=0 (= R₀)
        writef(f, fmtr, 4, zc[:, 0])    # Z cosine m=0 (= Z₀)
        for m in range(1, nmod):
            writef(f, fmtr, 4, rc[:, m])  # R cosine m
            writef(f, fmtr, 4, zs[:, m])  # Z sine m
            writef(f, fmtr, 4, rs[:, m])  # R sine m
            writef(f, fmtr, 4, zc[:, m])  # Z cosine m

        f.write("Safety factor\n"); writef(f, fmtr, 4, qsf)
        f.write("Current profile\n"); writef(f, fmtr, 4, jphi_avg)
        f.write("Covariant B phi\n"); writef(f, fmtr, 4, G_covBphi)


# =============================================================================
# Flip q sign consistently
# =============================================================================
def flip_q_sign_consistently(eq, flip_ffprim=False, flip_bcentr=True):
    """
    Return a copy of EQDSK dictionary with q sign flipped consistently
    by reversing the toroidal-field sign convention.

    This keeps ψ(R,Z) unchanged and flips:
      qpsi → -qpsi
      fpol → -fpol

    Optionally also flips:
      ffprim → -ffprim   (set flip_ffprim=True if needed)
      bcentr → -bcentr   (header consistency for toroidal field sign)

    Leaves unchanged:
      psirz, simag, sibry, pres, pprime, current, axis position, boundary
    """
    eq2 = dict(eq)
    eq2["qpsi"] = -np.asarray(eq["qpsi"], float)
    eq2["fpol"] = -np.asarray(eq["fpol"], float)
    eq2["ffprim"] = (-np.asarray(eq["ffprim"], float) if flip_ffprim
                     else np.asarray(eq["ffprim"], float).copy())
    eq2["bcentr"] = (-float(eq["bcentr"]) if flip_bcentr
                     else float(eq["bcentr"]))
    eq2["pres"] = np.asarray(eq["pres"], float).copy()
    eq2["pprime"] = np.asarray(eq["pprime"], float).copy()
    eq2["psirz"] = np.asarray(eq["psirz"], float).copy()
    return eq2


# =============================================================================
# Grad-Shafranov residual diagnostic
# =============================================================================
def compute_gs_residual(eq, verbose=True, plot=False, plot_file=None):
    """
    Compute the Grad-Shafranov residual on the EQDSK rectangular grid.
    
    Profiles p'(ψ) and FF'(ψ) are upsampled onto a fine ψ grid via cubic
    spline before evaluation, to avoid piecewise-linear kinks from the
    (possibly coarse) original EQDSK profile grid.
    """
    mu0 = 4.0e-7 * np.pi

    R = np.asarray(eq["R"], float)
    Z = np.asarray(eq["Z"], float)
    psirz = np.asarray(eq["psirz"], float)
    nw = len(R)
    nh = len(Z)

    simag = float(eq["simag"])
    sibry = float(eq["sibry"])

    # ------------------------------------------------------------------
    # Build high-resolution profile interpolators
    # (same logic as in enforce_gs_consistency)
    # ------------------------------------------------------------------
    nw_orig = len(eq["pprime"])
    psi_orig = np.linspace(simag, sibry, nw_orig)
    pprime_orig = np.asarray(eq["pprime"], float)
    ffprim_orig = np.asarray(eq["ffprim"], float)

    nfine = max(4 * max(nw, nh), 1024, nw_orig)

    if nfine > nw_orig:
        psi_fine = np.linspace(simag, sibry, nfine)
        spl_pp = CubicSpline(psi_orig, pprime_orig)
        spl_ff = CubicSpline(psi_orig, ffprim_orig)
        pprime_fine = spl_pp(psi_fine)
        ffprim_fine = spl_ff(psi_fine)
    else:
        psi_fine = psi_orig
        pprime_fine = pprime_orig
        ffprim_fine = ffprim_orig

    pprime_of_psi = interp1d(psi_fine, pprime_fine, kind="linear",
                             bounds_error=False,
                             fill_value=(pprime_fine[0], pprime_fine[-1]))
    ffprim_of_psi = interp1d(psi_fine, ffprim_fine, kind="linear",
                             bounds_error=False,
                             fill_value=(ffprim_fine[0], ffprim_fine[-1]))

    # Grid spacings
    dR = R[1] - R[0]
    dZ = Z[1] - Z[0]

    # Finite-difference derivatives (unchanged)
    d2psi_dR2 = np.zeros_like(psirz)
    d2psi_dR2[:, 1:-1] = (psirz[:, 2:] - 2.0 * psirz[:, 1:-1]
                           + psirz[:, :-2]) / (dR**2)

    dpsi_dR = np.zeros_like(psirz)
    dpsi_dR[:, 1:-1] = (psirz[:, 2:] - psirz[:, :-2]) / (2.0 * dR)

    d2psi_dZ2 = np.zeros_like(psirz)
    d2psi_dZ2[1:-1, :] = (psirz[2:, :] - 2.0 * psirz[1:-1, :]
                           + psirz[:-2, :]) / (dZ**2)

    ZZ, RR = np.meshgrid(Z, R, indexing="ij")
    lhs = d2psi_dR2 - dpsi_dR / np.maximum(RR, 1e-10) + d2psi_dZ2

    # RHS using upsampled profiles
    pp_2d = pprime_of_psi(psirz)
    ff_2d = ffprim_of_psi(psirz)
    rhs = -mu0 * RR**2 * pp_2d - ff_2d

    residual = lhs - rhs

    # Inside mask (unchanged)
    dpsi_norm = sibry - simag
    if abs(dpsi_norm) < 1e-300:
        dpsi_norm = 1.0
    arr = (psirz - simag) / dpsi_norm

    inside_mask = np.zeros_like(psirz, dtype=bool)
    inside_mask[2:-2, 2:-2] = ((arr[2:-2, 2:-2] >= 0.0)
                                & (arr[2:-2, 2:-2] <= 1.0))

    # Statistics (unchanged)
    res_in = residual[inside_mask]
    rhs_in = rhs[inside_mask]
    lhs_in = lhs[inside_mask]

    rms_res = float(np.sqrt(np.mean(res_in**2)))
    max_res = float(np.max(np.abs(res_in)))
    rms_rhs = float(np.sqrt(np.mean(rhs_in**2)))
    rms_lhs = float(np.sqrt(np.mean(lhs_in**2)))
    rms_norm = rms_res / (rms_rhs + 1e-300)

    if verbose:
        print("\n=== Grad-Shafranov residual diagnostic ===")
        print(f"  Grid:  nw={nw}, nh={nh}, dR={dR:.5e} m, dZ={dZ:.5e} m")
        print(f"  Profile ψ grid: {nw_orig} -> {nfine} points (cubic spline)")
        print(f"  Points inside LCFS: {inside_mask.sum()}")
        print(f"  RMS(Δ*ψ)      inside LCFS: {rms_lhs:.6e}")
        print(f"  RMS(RHS)      inside LCFS: {rms_rhs:.6e}")
        print(f"  RMS(residual) inside LCFS: {rms_res:.6e}")
        print(f"  max|residual| inside LCFS: {max_res:.6e}")
        print(f"  RMS(res)/RMS(RHS) (rel.):  {rms_norm:.6e}")
        if rms_norm < 0.01:
            print("  -> Excellent GS consistency (< 1%)")
        elif rms_norm < 0.05:
            print("  -> Good GS consistency (< 5%)")
        elif rms_norm < 0.10:
            print("  -> Acceptable GS consistency (< 10%)")
        else:
            print("  -> WARNING: significant GS residual (>10%)")
        print("============================================\n")

    if plot:
        fig, axes = plt.subplots(1, 4, figsize=(20, 5.5))

        def masked_outside(field):
            f = field.copy()
            f[~inside_mask] = np.nan
            return f

        im0 = axes[0].contourf(R, Z, masked_outside(lhs), levels=60,
                                cmap="RdBu_r")
        axes[0].set_title(r"$\Delta^*\psi$ (LHS)")
        axes[0].set_xlabel("R [m]"); axes[0].set_ylabel("Z [m]")
        axes[0].set_aspect("equal")
        plt.colorbar(im0, ax=axes[0], shrink=0.8)

        im1 = axes[1].contourf(R, Z, masked_outside(rhs), levels=60,
                                cmap="RdBu_r")
        axes[1].set_title(r"$-\mu_0 R^2 p' - FF'$ (RHS)")
        axes[1].set_xlabel("R [m]"); axes[1].set_aspect("equal")
        plt.colorbar(im1, ax=axes[1], shrink=0.8)

        res_disp = masked_outside(residual)
        vmax = np.nanpercentile(np.abs(res_disp), 98)
        im2 = axes[2].contourf(R, Z, res_disp, levels=60,
                                cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        axes[2].set_title("Residual (LHS−RHS)")
        axes[2].set_xlabel("R [m]"); axes[2].set_aspect("equal")
        plt.colorbar(im2, ax=axes[2], shrink=0.8)

        rel = np.abs(residual) / (np.abs(rhs) + 1e-30)
        rel_disp = masked_outside(rel)
        vmax_rel = np.nanpercentile(rel_disp, 95)
        im3 = axes[3].contourf(R, Z, rel_disp, levels=60,
                                cmap="hot_r", vmin=0,
                                vmax=min(vmax_rel, 1.0))
        axes[3].set_title("|res|/|RHS| (local)")
        axes[3].set_xlabel("R [m]"); axes[3].set_aspect("equal")
        plt.colorbar(im3, ax=axes[3], shrink=0.8)

        plt.tight_layout()
        if plot_file:
            plt.savefig(plot_file, dpi=200, bbox_inches="tight")
            print(f"  Saved GS residual plot: {plot_file}")
        plt.show()

    return dict(
        R=R, Z=Z, residual=residual, lhs=lhs, rhs=rhs,
        arr=arr, inside_mask=inside_mask,
        rms_inside=rms_res, max_inside=max_res,
        rms_rhs=rms_rhs, rms_norm=rms_norm
    )

# =============================================================================
# Grad-Shafranov SOR refinement
# =============================================================================
def enforce_gs_consistency(R, Z, psirz, eq,
                           inside_mask,
                           n_iter=200,
                           omega=0.6,
                           fix_boundary=True,
                           verbose=True):
    """
    Iteratively adjust ψ(R,Z) inside the plasma so that it better satisfies
    the Grad-Shafranov equation, using successive over-relaxation (SOR).

    The GS equation:
        ∂²ψ/∂R² - (1/R) ∂ψ/∂R + ∂²ψ/∂Z² = -μ₀ R² p'(ψ) - FF'(ψ)

    is solved iteratively, keeping boundary values (LCFS and outside) fixed.

    The profiles p'(ψ) and FF'(ψ) are resampled onto a fine ψ grid using
    cubic spline interpolation before building the lookup functions.  This
    avoids kinks from piecewise-linear interpolation of the original
    (possibly coarse) EQDSK profile grid.

    Parameters
    ----------
    R, Z : 1D arrays
        Rectangular grid coordinates.
    psirz : 2D array (nh, nw)
        Initial reconstructed ψ(R,Z). Modified in place on return.
    eq : dict
        EQDSK dictionary containing profiles pprime, ffprim, simag, sibry.
    inside_mask : 2D bool array (nh, nw)
        True for grid points inside the plasma where GS should be enforced.
    n_iter : int
        Number of SOR iterations.
    omega : float
        SOR relaxation parameter (0 < ω < 2).
    fix_boundary : bool
        If True, erode inside_mask by 1 pixel to keep LCFS boundary fixed.
    verbose : bool
        Print convergence info every 50 iterations.

    Returns
    -------
    psirz : 2D array
        The corrected field (same array, modified in place).
    history : list of float
        RMS residual at each iteration.
    """
    mu0 = 4.0e-7 * np.pi

    nw = len(R)
    nh = len(Z)
    dR = R[1] - R[0]
    dZ = Z[1] - Z[0]

    simag = float(eq["simag"])
    sibry = float(eq["sibry"])

    # ------------------------------------------------------------------
    # Build high-resolution profile interpolators.
    #
    # The original EQDSK profiles are on a uniform ψ grid with nw_orig
    # points (e.g. 65).  If the spatial grid is much finer (e.g. 257×257),
    # piecewise-linear interpolation of the profiles introduces kinks at
    # each original ψ knot, which show up as localized GS residual spikes.
    #
    # Solution: first build cubic splines of the original profiles on the
    # coarse ψ grid, then evaluate them on a much finer ψ grid (nfine
    # points).  The final lookup interpolators are piecewise-linear on
    # this fine grid, which is effectively smooth at the spatial grid
    # resolution.
    # ------------------------------------------------------------------
    nw_orig = len(eq["pprime"])
    psi_orig = np.linspace(simag, sibry, nw_orig)
    pprime_orig = np.asarray(eq["pprime"], float)
    ffprim_orig = np.asarray(eq["ffprim"], float)

    # Choose fine ψ grid resolution: at least 4× the spatial grid dimension
    # or 1024, whichever is larger, to ensure the profile is smooth at
    # the scale of the spatial grid.
    nfine = max(4 * max(nw, nh), 1024, nw_orig)

    if nfine > nw_orig:
        # Cubic spline interpolation of profiles onto fine ψ grid
        psi_fine = np.linspace(simag, sibry, nfine)

        # Use "not-a-knot" boundary conditions (robust default)
        spl_pp = CubicSpline(psi_orig, pprime_orig)
        spl_ff = CubicSpline(psi_orig, ffprim_orig)

        pprime_fine = spl_pp(psi_fine)
        ffprim_fine = spl_ff(psi_fine)

        if verbose:
            print(f"  GS-SOR: profiles upsampled from {nw_orig} to "
                  f"{nfine} ψ points (cubic spline)")
    else:
        # Original grid is already fine enough
        psi_fine = psi_orig
        pprime_fine = pprime_orig
        ffprim_fine = ffprim_orig

    # Build piecewise-linear lookup on the fine grid
    # (fast evaluation, effectively smooth at spatial resolution)
    pprime_interp = interp1d(psi_fine, pprime_fine, kind="linear",
                             bounds_error=False,
                             fill_value=(pprime_fine[0], pprime_fine[-1]))
    ffprim_interp = interp1d(psi_fine, ffprim_fine, kind="linear",
                             bounds_error=False,
                             fill_value=(ffprim_fine[0], ffprim_fine[-1]))

    # 2D coordinate mesh
    ZZ, RR = np.meshgrid(Z, R, indexing="ij")

    # Build interior mask: exclude grid edges (FD stencil needs neighbors)
    interior = inside_mask.copy()
    interior[0, :] = False
    interior[-1, :] = False
    interior[:, 0] = False
    interior[:, -1] = False

    if fix_boundary:
        # Erode by 1 pixel: the outermost inside ring stays fixed,
        # preserving the LCFS boundary shape
        interior = binary_erosion(interior, iterations=1)

    n_interior = int(interior.sum())
    if verbose:
        print(f"  GS-SOR: {n_interior} interior points, "
              f"omega={omega:.3f}, n_iter={n_iter}")

    # ------------------------------------------------------------------
    # Precompute FD stencil coefficients (constant across iterations)
    # ------------------------------------------------------------------
    inv_dR2 = 1.0 / (dR**2)
    inv_dZ2 = 1.0 / (dZ**2)
    D = 2.0 * inv_dR2 + 2.0 * inv_dZ2   # diagonal

    Rclip = np.maximum(RR, 1e-10)  # avoid division by zero at R=0
    aR_p = inv_dR2 - 1.0 / (2.0 * Rclip * dR)   # R+dR neighbor
    aR_m = inv_dR2 + 1.0 / (2.0 * Rclip * dR)   # R-dR neighbor
    aZ = inv_dZ2                                    # Z±dZ neighbors

    psi = psirz.copy()
    history = []

    for it in range(n_iter):
        # Compute GS source term from current ψ (nonlinear: depends on ψ)
        pp = pprime_interp(psi)
        ff = ffprim_interp(psi)
        source = -mu0 * RR**2 * pp - ff

        # GS-implied ψ at each interior point (Jacobi-like estimate)
        psi_gs = np.zeros_like(psi)
        psi_gs[1:-1, 1:-1] = (
            aR_p[1:-1, 1:-1] * psi[1:-1, 2:]       # R+dR neighbor
            + aR_m[1:-1, 1:-1] * psi[1:-1, :-2]     # R-dR neighbor
            + aZ * (psi[2:, 1:-1] + psi[:-2, 1:-1]) # Z±dZ neighbors
            - source[1:-1, 1:-1]                      # GS source
        ) / D

        # SOR update: only on interior points, boundary stays fixed
        correction = psi_gs - psi
        psi_new = psi.copy()
        psi_new[interior] = psi[interior] + omega * correction[interior]

        # ------------------------------------------------------------------
        # Compute residual for convergence monitoring
        # ------------------------------------------------------------------
        d2R = (psi_new[:, 2:] - 2 * psi_new[:, 1:-1]
               + psi_new[:, :-2]) * inv_dR2
        dR1 = (psi_new[:, 2:] - psi_new[:, :-2]) / (2.0 * dR)
        d2Z = (psi_new[2:, :] - 2 * psi_new[1:-1, :]
               + psi_new[:-2, :]) * inv_dZ2

        lhs_inner = np.zeros_like(psi_new)
        lhs_inner[1:-1, 1:-1] = (
            d2R[1:-1, :]
            - dR1[1:-1, :] / np.maximum(RR[1:-1, 1:-1], 1e-10)
            + d2Z[:, 1:-1]
        )

        pp_new = pprime_interp(psi_new)
        ff_new = ffprim_interp(psi_new)
        rhs_inner = -mu0 * RR**2 * pp_new - ff_new

        res = lhs_inner - rhs_inner
        res_inside = res[interior]
        rms = float(np.sqrt(np.mean(res_inside**2)))
        history.append(rms)

        psi = psi_new

        if verbose and (it % 50 == 0 or it == n_iter - 1):
            print(f"    iter {it:4d}: RMS residual = {rms:.6e}")

    # Copy corrected field back into the input array
    psirz[:] = psi[:]
    return psirz, history

# =============================================================================
# Build smooth ψ(R,Z) from TORIC Fourier surfaces
# =============================================================================
def build_psirz_from_toric_surfaces(eq, psi_axis_2d, psi_lcfs_2d,
                                    rho_pol, rc, rs, zc, zs,
                                    ntheta=901, nrho_dense=901,
                                    nw_out=None, nh_out=None,
                                    outside_mode="original",
                                    enforce_monotone_arr=True,
                                    arr_floor=0.0,
                                    gs_refine=False,
                                    gs_n_iter=200,
                                    gs_omega=0.6,
                                    gs_verbose=True):
    """
    Smoothly reconstruct ψ(R,Z) from TORIC Fourier surfaces on a rectangular grid.

    Pipeline:
      1) Build output rectangular grid (optionally different resolution from input).
      2) Build dense rho-theta representation of reconstructed surfaces by
         interpolating Fourier coefficients onto a fine radial grid.
      3) Assign ψ = ψ_axis + (ψ_lcfs - ψ_axis)·ρ² on each surface.
      4) Scatter-interpolate the (R,Z,ψ) cloud onto the rectangular grid
         using LinearNDInterpolator.
      5) Determine inside/outside LCFS from the outermost reconstructed surface.
      6) Fill outside region using original EQDSK, constant ψ_lcfs, or NaN.
      7) Enforce monotonicity: arr = (ψ - ψ_axis)/(ψ_lcfs - ψ_axis) >= arr_floor
         inside the LCFS (removes interpolation undershoots near the axis).
      8) Optionally refine ψ with GS-SOR iteration to improve Grad-Shafranov
         consistency while keeping boundary values fixed.
      9) Re-enforce monotonicity after GS refinement.

    Parameters
    ----------
    eq : dict
        Original EQDSK dictionary (used for outside fill and profiles).
    psi_axis_2d, psi_lcfs_2d : float
        Axis and LCFS psi values for the ρ² ↔ ψ mapping.
    rho_pol : 1D array
        Radial grid of reconstructed surfaces (from TORIC output).
    rc, rs, zc, zs : 2D arrays, shape (nmhd, nmod)
        Fourier coefficients for R(θ) and Z(θ) on each surface.
    ntheta : int
        Number of poloidal points for surface reconstruction.
    nrho_dense : int
        Number of dense radial points for the scattered interpolation cloud.
    nw_out, nh_out : int or None
        Output EQDSK grid dimensions. If None, use original nw, nh.
    outside_mode : str
        How to fill outside LCFS: "original", "lcfs", or "nan".
    enforce_monotone_arr : bool
        If True, clip arr >= arr_floor inside the LCFS.
    arr_floor : float
        Lower bound for arr (normalized flux) inside LCFS.
    gs_refine : bool
        If True, apply GS-SOR iterative refinement after interpolation.
    gs_n_iter : int
        Number of SOR iterations for GS refinement.
    gs_omega : float
        SOR relaxation parameter (0 < ω < 2).
    gs_verbose : bool
        Print GS convergence info.

    Returns
    -------
    Rg_new, Zg_new : 1D arrays
        New rectangular grid coordinate arrays.
    psirz_new : 2D array (nh_new, nw_new)
        Reconstructed poloidal flux on the new grid.
    """
    from matplotlib.path import Path

    # ------------------------------------------------------------------
    # Step 1: Build output rectangular grid
    # ------------------------------------------------------------------
    if nw_out is None:
        Rg_new = eq["R"].copy()
    else:
        Rg_new = np.linspace(eq["R"][0], eq["R"][-1], int(nw_out))

    if nh_out is None:
        Zg_new = eq["Z"].copy()
    else:
        Zg_new = np.linspace(eq["Z"][0], eq["Z"][-1], int(nh_out))

    # ------------------------------------------------------------------
    # Step 2: Dense rho-theta surface reconstruction
    # Interpolate Fourier coefficients onto a finer radial grid,
    # then reconstruct (R,Z) at each (rho, theta) point.
    # ------------------------------------------------------------------
    theta = np.linspace(0.0, TWOPI, ntheta, endpoint=False)
    rho_dense = np.linspace(rho_pol[0], rho_pol[-1], nrho_dense)

    def interp_coeff(arr2d):
        """Interpolate 2D Fourier coefficient array onto dense rho grid."""
        out = np.zeros((nrho_dense, arr2d.shape[1]))
        for m in range(arr2d.shape[1]):
            out[:, m] = np.interp(rho_dense, rho_pol, arr2d[:, m])
        return out

    rc_d = interp_coeff(rc)
    rs_d = interp_coeff(rs)
    zc_d = interp_coeff(zc)
    zs_d = interp_coeff(zs)

    RRp = np.zeros((nrho_dense, ntheta))
    ZZp = np.zeros((nrho_dense, ntheta))

    for i in range(nrho_dense):
        Rtmp = np.full(ntheta, rc_d[i, 0], float)
        Ztmp = np.full(ntheta, zc_d[i, 0], float)
        for m in range(1, rc.shape[1]):
            cm = np.cos(m * theta)
            sm = np.sin(m * theta)
            Rtmp += rc_d[i, m] * cm + rs_d[i, m] * sm
            Ztmp += zc_d[i, m] * cm + zs_d[i, m] * sm
        RRp[i, :] = Rtmp
        ZZp[i, :] = Ztmp

    # ------------------------------------------------------------------
    # Step 3: Assign ψ from the ρ² relation
    # ψ(ρ) = ψ_axis + (ψ_lcfs - ψ_axis) · ρ²
    # ------------------------------------------------------------------
    psi_dense = psi_axis_2d + (psi_lcfs_2d - psi_axis_2d) * rho_dense**2
    PP = np.repeat(psi_dense[:, None], ntheta, axis=1)

    # ------------------------------------------------------------------
    # Step 4: Scattered interpolation (R,Z) → ψ inside the LCFS
    # ------------------------------------------------------------------
    pts = np.column_stack([RRp.ravel(), ZZp.ravel()])
    vals = PP.ravel()
    interp_inside = LinearNDInterpolator(pts, vals, fill_value=np.nan)

    ZZg, RRg = np.meshgrid(Zg_new, Rg_new, indexing="ij")
    psi_inside = interp_inside(RRg, ZZg)

    # ------------------------------------------------------------------
    # Step 5: Inside/outside LCFS mask from outermost reconstructed surface
    # ------------------------------------------------------------------
    R_lcfs = RRp[-1, :]
    Z_lcfs = ZZp[-1, :]
    lcfs_path = Path(np.column_stack([R_lcfs, Z_lcfs]))
    inside = lcfs_path.contains_points(
        np.column_stack([RRg.ravel(), ZZg.ravel()])
    ).reshape(ZZg.shape)

    # ------------------------------------------------------------------
    # Step 6: Fill outside region
    # ------------------------------------------------------------------
    if outside_mode == "original":
        # Interpolate original EQDSK ψ onto new grid for outside region
        interp_orig = RegularGridInterpolator(
            (eq["Z"], eq["R"]), eq["psirz"],
            bounds_error=False, fill_value=np.nan)
        psi_orig = interp_orig(
            np.column_stack([ZZg.ravel(), RRg.ravel()])
        ).reshape(ZZg.shape)
        psirz_new = np.array(psi_orig, copy=True)
        m = inside & np.isfinite(psi_inside)
        psirz_new[m] = psi_inside[m]
    elif outside_mode == "lcfs":
        psirz_new = np.full_like(psi_inside, psi_lcfs_2d)
        m = inside & np.isfinite(psi_inside)
        psirz_new[m] = psi_inside[m]
    elif outside_mode == "nan":
        psirz_new = np.full_like(psi_inside, np.nan)
        m = inside & np.isfinite(psi_inside)
        psirz_new[m] = psi_inside[m]
    else:
        raise ValueError("outside_mode must be 'original', 'lcfs', or 'nan'")

    # ------------------------------------------------------------------
    # Step 7: Monotonicity enforcement (pre-GS)
    # Ensure arr = (ψ - ψ_axis)/(ψ_lcfs - ψ_axis) >= arr_floor inside LCFS.
    # This removes interpolation undershoots near the magnetic axis where
    # ψ could be slightly lower than the axis value.
    # ------------------------------------------------------------------
    if enforce_monotone_arr:
        dpsi = psi_lcfs_2d - psi_axis_2d
        if abs(dpsi) > 1e-300:
            arr = (psirz_new - psi_axis_2d) / dpsi
            print("\n ===> max ",arr[inside].max(), arr[inside].min(),"\n")
            arr_in = arr[inside]
            arr_in = np.maximum(arr_in, arr_floor)
            psirz_new[inside] = psi_axis_2d + dpsi * arr_in
            arr = (psirz_new - psi_axis_2d) / dpsi
            print("\n ===> max ",arr[inside].max(), arr[inside].min(),"\n")

    # ------------------------------------------------------------------
    # Step 8: Optional GS-SOR refinement
    # Iteratively adjust ψ at interior points to better satisfy the
    # Grad-Shafranov equation, while keeping boundary values fixed.
    # ------------------------------------------------------------------
    if gs_refine:
        if gs_verbose:
            print("\n  Starting GS-consistent SOR refinement...")

        # Build GS interior mask (exclude grid edges for FD stencil)
        inside_for_gs = inside.copy()
        inside_for_gs[0, :] = False
        inside_for_gs[-1, :] = False
        inside_for_gs[:, 0] = False
        inside_for_gs[:, -1] = False

        psirz_new, gs_history = enforce_gs_consistency(
            Rg_new, Zg_new, psirz_new, eq,
            inside_mask=inside_for_gs,
            n_iter=gs_n_iter,
            omega=gs_omega,
            fix_boundary=True,
            verbose=gs_verbose
        )

        if gs_verbose:
            print(f"  GS refinement done. "
                  f"Final RMS residual: {gs_history[-1]:.6e}\n")

        # ------------------------------------------------------------------
        # Step 9: Re-enforce monotonicity after GS refinement
        # The SOR iteration may introduce new small undershoots.
        # ------------------------------------------------------------------
        if enforce_monotone_arr:
            dpsi = psi_lcfs_2d - psi_axis_2d
            if abs(dpsi) > 1e-300:
                arr = (psirz_new - psi_axis_2d) / dpsi
                arr_in = arr[inside]
                arr_in = np.maximum(arr_in, arr_floor)
                psirz_new[inside] = psi_axis_2d + dpsi * arr_in

    return Rg_new, Zg_new, psirz_new, inside
  
  
def monotonicity(psirz,inside,psi_axis,psi_lcfs,arr_floor):
    # ------------------------------------------------------------------
    # Monotonicity enforcement:
    #     arr = (psi - psi_axis)/(psi_lcfs - psi_axis) must be >= arr_floor
    #     inside the plasma.
    # ------------------------------------------------------------------
    dpsi = psi_lcfs - psi_axis
    if abs(dpsi) > 1e-300:
      arr = (psirz - psi_axis) / dpsi
      
      # Only enforce inside the reconstructed LCFS
      arr_inside = arr[inside]
      arr_inside = np.maximum(arr_inside, arr_floor)
      psirz[inside] = psi_axis + dpsi * arr_inside

    return psirz
  
# =============================================================================
# Build and write consistent smooth EQDSK
# =============================================================================
def build_and_write_consistent_smooth_eqdsk(
        eq, rho_pol, rc, rs, zc, zs,
        out_path,
        ntheta_geom=721,
        nrho_dense=201,
        nw_out=None, nh_out=None,
        outside_mode="original",
        width=16,
        flip_q_sign=False,
        flip_ffprim_too=False,
        title="rbb-SMOOTHED_EQDSK",
        enforce_monotone_arr=True,
        arr_floor=0.0,
        gs_refine=False,
        gs_n_iter=200,
        gs_omega=0.6):
    """
    Build a smoothed EQDSK from the final TORIC surface representation
    and write it self-consistently.

    This routine:
      1) Reconstructs a smooth ψ(R,Z) on a chosen rectangular grid from
         the Fourier surface coefficients (with optional monotonicity
         enforcement and GS-SOR refinement).
      2) Recomputes simag and sibry from the new 2D field by evaluating
         ψ at the magnetic axis and at the original LCFS boundary points.
      3) Resamples all 1D profiles (fpol, pres, ffprim, pprime, qpsi)
         onto the new ψ grid using linear interpolation from the old grid.
      4) Writes a strict EQDSK file (including limiter contour if present).

    Parameters
    ----------
    eq : dict
        Original EQDSK dictionary.
    rho_pol : 1D array
        Radial grid of the final TORIC surfaces.
    rc, rs, zc, zs : 2D arrays
        Fourier coefficients for R(θ), Z(θ) on each surface.
    out_path : str
        Output EQDSK filename.
    ntheta_geom : int
        Number of θ points for surface reconstruction.
    nrho_dense : int
        Number of dense ρ points for scattered interpolation.
    nw_out, nh_out : int or None
        Output grid dimensions. If None, keep original.
    outside_mode : str
        How to fill ψ outside the LCFS: "original", "lcfs", or "nan".
    width : int
        Fixed-width float format width for the EQDSK writer.
    flip_q_sign : bool
        If True, flip q sign consistently (qpsi, fpol, optionally bcentr).
    flip_ffprim_too : bool
        If True, also flip ffprim when flipping q.
    title : str
        Header title for the new EQDSK file.
    enforce_monotone_arr : bool
        If True, clip normalized flux arr >= arr_floor inside LCFS.
    arr_floor : float
        Lower bound for arr inside LCFS.
    gs_refine : bool
        If True, apply GS-SOR refinement.
    gs_n_iter : int
        Number of SOR iterations.
    gs_omega : float
        SOR relaxation parameter.

    Returns
    -------
    eq_new : dict
        New self-consistent EQDSK dictionary.
    psirz_new : 2D ndarray
        Rebuilt smooth ψ(R,Z).
    """
    # ------------------------------------------------------------------
    # Estimate axis and LCFS ψ from the ORIGINAL field.
    # These are used as the initial ρ² ↔ ψ mapping for rebuilding psirz.
    # ------------------------------------------------------------------
    itp_old = RegularGridInterpolator(
        (eq["Z"], eq["R"]), eq["psirz"],
        bounds_error=False, fill_value=np.nan)
    psi_axis_old = float(itp_old([[eq["zaxis"], eq["raxis"]]])[0])

    if eq["rbbbs"] is not None and eq["zbbbs"] is not None:
        bb_old = itp_old(np.column_stack([eq["zbbbs"], eq["rbbbs"]]))
        bb_old = bb_old[np.isfinite(bb_old)]
        psi_lcfs_old = float(np.median(bb_old))
    else:
        # Fallback if no boundary is present
        psi_lcfs_old = float(eq["sibry"])

    # ------------------------------------------------------------------
    # Rebuild smooth ψ(R,Z) from the final Fourier surfaces
    # ------------------------------------------------------------------
    Rg_new, Zg_new, psirz_new, inside_lcfs = build_psirz_from_toric_surfaces(
        eq=eq,
        psi_axis_2d=psi_axis_old,
        psi_lcfs_2d=psi_lcfs_old,
        rho_pol=rho_pol,
        rc=rc, rs=rs, zc=zc, zs=zs,
        ntheta=ntheta_geom,
        nrho_dense=nrho_dense,
        nw_out=nw_out,
        nh_out=nh_out,
        outside_mode=outside_mode,
        gs_refine=gs_refine,
        gs_n_iter=gs_n_iter,
        gs_omega=gs_omega,
        gs_verbose=True
    )

    # ------------------------------------------------------------------
    # Recompute simag and sibry from the NEW 2D field
    # ------------------------------------------------------------------
    itp_new = RegularGridInterpolator(
        (Zg_new, Rg_new), psirz_new,
        bounds_error=False, fill_value=np.nan)
    psi_axis_new = float(itp_new([[eq["zaxis"], eq["raxis"]]])[0])

    if eq["rbbbs"] is not None and eq["zbbbs"] is not None:
        bb_new = itp_new(np.column_stack([eq["zbbbs"], eq["rbbbs"]]))
        bb_new = bb_new[np.isfinite(bb_new)]
        psi_lcfs_new = float(np.median(bb_new))
    else:
        psi_lcfs_new = float(np.nanmax(psirz_new))
        
    # ------------------------------------------------------------------
    # Optional enforce monotonicity
    # ------------------------------------------------------------------
    if enforce_monotone_arr:
      psirz_new = monotonicity(psirz_new,inside_lcfs,psi_axis_new,psi_lcfs_new,arr_floor)
      
    # ------------------------------------------------------------------
    # Build a new EQDSK dictionary with updated grid and ψ
    # ------------------------------------------------------------------
    eq_new = dict(eq)
    eq_new["R"] = np.asarray(Rg_new, float)
    eq_new["Z"] = np.asarray(Zg_new, float)
    eq_new["psirz"] = np.asarray(psirz_new, float)
    eq_new["nw"] = len(Rg_new)
    eq_new["nh"] = len(Zg_new)
    eq_new["rleft"] = float(Rg_new[0])
    eq_new["rdim"] = float(Rg_new[-1] - Rg_new[0])
    eq_new["zmid"] = float(0.5 * (Zg_new[0] + Zg_new[-1]))
    eq_new["zdim"] = float(Zg_new[-1] - Zg_new[0])
    eq_new["simag"] = float(psi_axis_new)
    eq_new["sibry"] = float(psi_lcfs_new)

    # ------------------------------------------------------------------
    # Resample 1D profiles onto the NEW ψ grid
    # IMPORTANT:
    #   old profiles live on psi_old = linspace(old simag, old sibry, old nw)
    #   new profiles must live on psi_new = linspace(new simag, new sibry, new nw)
    # ------------------------------------------------------------------
    psi_old = np.linspace(eq["simag"], eq["sibry"], eq["nw"])
    psi_new = np.linspace(eq_new["simag"], eq_new["sibry"], eq_new["nw"])

    for key in ["fpol", "pres", "ffprim", "pprime", "qpsi"]:
        eq_new[key] = np.interp(psi_new, psi_old,
                                np.asarray(eq[key], float))
    eq_new["psi1d"] = psi_new

    # Optional q sign flip
    if flip_q_sign:
        eq_new = flip_q_sign_consistently(
            eq_new, flip_ffprim=flip_ffprim_too, flip_bcentr=True)

    # ------------------------------------------------------------------
    # Write strict EQDSK (including limiter contour if present)
    # ------------------------------------------------------------------
    write_geqdsk(eq_new, eq_new["psirz"], out_path=out_path,
                 width=width, title=title)

    return eq_new, psirz_new


# =============================================================================
# Main
# =============================================================================
def main():
    ap = argparse.ArgumentParser(
        description=(
            "Extract EQDSK flux surfaces, parameterize by a chosen poloidal "
            "angle, fit Fourier modes, optionally regularize near-axis, "
            "and write TORIC.equi.  Optionally write a smoothed EQDSK with "
            "monotonicity enforcement and GS-consistent SOR refinement."
        )
    )

    # -------------------------------------------------------------------------
    # Required positional argument
    # -------------------------------------------------------------------------
    ap.add_argument("eqdsk",
                    help="Path to input gEQDSK file (EFIT equilibrium).")

    # -------------------------------------------------------------------------
    # Poloidal angle choice
    # -------------------------------------------------------------------------
    ap.add_argument(
        "--poloidal_angle", default="flux",
        choices=["arclength", "geometric", "flux"],
        help=(
            "Choice of poloidal angle θ used to sample each surface before "
            "Fourier decomposition. "
            "'arclength': θ=2π s/L (uniform arclength). "
            "'geometric': θ=atan2(Z-Zax, R-Rax), shifted so θ=0 at max R. "
            "'flux': flux-like proxy rectification of geometric θ using "
            "v~(dl/dθ)/R."
        )
    )
    ap.add_argument(
        "--theta0", type=float, default=1.0,
        help=(
            "(d. 1.0) If > 0, θ=0 at the crossing with the midplane at "
            "Zaxis; otherwise θ=0 at maximum R."
        )
    )
    ap.add_argument(
        "--theta_eps", type=float, default=1e-10,
        help=(
            "Only used for --poloidal_angle geometric/flux. "
            "Deduplication tolerance (radians) to merge near-duplicate θ "
            "values so CubicSpline can be built."
        )
    )

    # -------------------------------------------------------------------------
    # Radial surface selection (rho_pol)
    # -------------------------------------------------------------------------
    ap.add_argument(
        "--rho_min", type=float, default=0.08,
        help=(
            "(d. 0.08) Minimum normalized poloidal-flux radius rho_pol used "
            "for surface extraction (0=axis, 1=LCFS). Avoid too small values "
            "if contours are noisy near the axis."
        )
    )
    ap.add_argument(
        "--rho_max", type=float, default=0.99,
        help=(
            "(d. 0.99) Maximum rho_pol used for surface extraction (<=1). "
            "Use <1 to stay slightly inside the LCFS if boundary is noisy."
        )
    )
    ap.add_argument(
        "--nmhd", type=int, default=301,
        help=(
            "(d. 301) Number of flux surfaces (radial points) to extract "
            "between rho_min and rho_max."
        )
    )

    # -------------------------------------------------------------------------
    # Contour extraction controls
    # -------------------------------------------------------------------------
    ap.add_argument(
        "--upsample", type=int, default=8,
        help=(
            "(d. 8) Upsampling factor applied to the ψ(R,Z) grid before "
            "contour extraction. Higher values give smoother/more accurate "
            "contours but cost more CPU/memory."
        )
    )
    ap.add_argument(
        "--min_points", type=int, default=8,
        help=(
            "(d. 8) Minimum number of points in a candidate contour. Shorter "
            "contours are discarded (helps reject spurious tiny contours)."
        )
    )

    # -------------------------------------------------------------------------
    # Poloidal discretization
    # -------------------------------------------------------------------------
    ap.add_argument(
        "--ntheta", type=int, default=512,
        help=(
            "(d. 512) Number of uniform θ samples per flux surface after "
            "reparameterization. This is also the FFT sample count."
        )
    )

    # -------------------------------------------------------------------------
    # Fourier mode selection
    # -------------------------------------------------------------------------
    ap.add_argument(
        "--mmax_max", type=int, default=12,
        help="(d. 12) Maximum poloidal mode number m allowed during condensation."
    )
    ap.add_argument(
        "--condense_tol", type=float, default=5.e-3,
        help=(
            "(d. 5.e-3) RMS geometric reconstruction error tolerance [m] used "
            "by condensation to pick mmax. Smaller values keep more modes."
        )
    )
    ap.add_argument(
        "--mmax_final", type=int, default=0,
        help=(
            "(d. 0) If >0, force a fixed mmax for all surfaces (disables "
            "condensation). If 0, mmax is chosen per surface via condensation."
        )
    )

    # -------------------------------------------------------------------------
    # Contour pre-smoothing
    # -------------------------------------------------------------------------
    ap.add_argument(
        "--pre_window", type=int, default=21,
        help=(
            "(d. 21) Savitzky–Golay window length (points) for periodic "
            "smoothing of contour coordinates. Must be odd."
        )
    )
    ap.add_argument(
        "--pre_poly", type=int, default=3,
        help=(
            "(d. 3) Savitzky–Golay polynomial order for contour smoothing "
            "(must be < pre_window)."
        )
    )

    # -------------------------------------------------------------------------
    # Optional ψ(R,Z) smoothing (applied to input before contour extraction)
    # -------------------------------------------------------------------------
    ap.add_argument(
        "--no_smooth_psirz", dest="smooth_psirz",
        action="store_false", default=True,
        help=(
            "Disable axis-preserving Gaussian smoothing of ψ(R,Z) before "
            "contour extraction. This smoothing is useful if ψ(R,Z) is noisy "
            "and contours look jagged."
        )
    )
    ap.add_argument(
        "--smooth_sigma_R", type=float, default=2.0,
        help=(
            "(d. 2.0) Gaussian smoothing sigma along the R-index direction "
            "(grid units, not meters)."
        )
    )
    ap.add_argument(
        "--smooth_sigma_Z", type=float, default=2.0,
        help=(
            "(d. 2.0) Gaussian smoothing sigma along the Z-index direction "
            "(grid units, not meters)."
        )
    )
    ap.add_argument(
        "--smooth_r0", type=float, default=0.12,
        help=(
            "(d. 0.12) Inner radius (meters) of the axis-preserving blend "
            "region. For rr < r0: mostly original ψ."
        )
    )
    ap.add_argument(
        "--smooth_dr", type=float, default=0.06,
        help=(
            "(d. 0.06) Transition width (meters) for blending from original "
            "to smoothed ψ. For rr > r0+dr: fully smoothed ψ."
        )
    )

    # -------------------------------------------------------------------------
    # Optional radial regularization of Fourier coefficients
    # -------------------------------------------------------------------------
    ap.add_argument(
        "--no_reg_smooth", dest="reg_smooth",
        action="store_false", default=True,
        help=(
            "Disable near-axis parity regularization of Fourier coefficients "
            "across radius. Even m modes: constant+quadratic near axis; "
            "odd m: linear near axis."
        )
    )
    ap.add_argument(
        "--nout", type=int, default=101,
        help=(
            "(d. 101) Number of radial output points (rhoPF_out) used when "
            "radial regularization is enabled."
        )
    )
    ap.add_argument(
        "--rhoPF_max", type=float, default=-1.0,
        help=(
            "(d. -1.0) Maximum rho_pol value for the regularized output grid. "
            "If <0, uses rho_max."
        )
    )
    ap.add_argument(
        "--rho_reg", type=float, default=-1.0,
        help=(
            "(d. -1.0 → 0.3) Near-axis regularization scale in rho_pol. "
            "For rho <= rho_reg, parity constraints are enforced."
        )
    )

    # -------------------------------------------------------------------------
    # COCOS identifier
    # -------------------------------------------------------------------------
    ap.add_argument(
        "--cocos_diag", action="store_true",
        help="Print COCOS inference diagnostics (heuristic; may be non-unique)."
    )

    # -------------------------------------------------------------------------
    # Plotting options
    # -------------------------------------------------------------------------
    ap.add_argument(
        "--plot", action="store_true",
        help=(
            "Generate a reconstruction plot comparing TORIC reconstruction "
            "vs EQDSK contours."
        )
    )
    ap.add_argument(
        "--plot_n_solid", type=int, default=20,
        help="(d. 20) Number of reconstructed (solid) TORIC surfaces to plot."
    )
    ap.add_argument(
        "--plot_dashed_between", type=int, default=1,
        help=(
            "(d. 1) Number of dashed EQDSK contours plotted between "
            "successive solid surfaces. 0 disables."
        )
    )
    ap.add_argument(
        "--theta_spokes", type=int, default=16,
        help=(
            "(d. 16) Number of constant-θ spokes drawn across radius on "
            "the plot. 0 disables."
        )
    )

    # -------------------------------------------------------------------------
    # EQDSK parsing/writing detail
    # -------------------------------------------------------------------------
    ap.add_argument(
        "--width", type=int, default=16,
        help="(d. 16) Fixed-field width used by EQDSK float parsing."
    )
    ap.add_argument(
        "--write_smooth_eqdsk", action="store_true",
        help="Write a smoothed EQDSK rebuilt from reconstructed TORIC surfaces."
    )
    ap.add_argument(
        "--smooth_eqdsk_nw", type=int, default=129,
        help="(d. 129) Number of R grid points in the smoothed EQDSK."
    )
    ap.add_argument(
        "--smooth_eqdsk_nh", type=int, default=129,
        help="(d. 129) Number of Z grid points in the smoothed EQDSK."
    )
    ap.add_argument(
        "--smooth_eqdsk_outside", default="original",
        choices=["original", "lcfs", "nan"],
        help="(d. 'original') How to fill psi outside the LCFS."
    )

    # -------------------------------------------------------------------------
    # q sign flip options
    # -------------------------------------------------------------------------
    ap.add_argument(
        "--flip_q_sign", action="store_true",
        help=(
            "Flip q sign consistently by reversing toroidal-field sign "
            "convention (qpsi, fpol and optionally bcentr)."
        )
    )
    ap.add_argument(
        "--flip_ffprim_too", action="store_true",
        help="Also flip ffprim when using --flip_q_sign."
    )

    # -------------------------------------------------------------------------
    # Monotonicity enforcement
    # -------------------------------------------------------------------------
    ap.add_argument(
        "--no_enforce_monotone_arr", dest="enforce_monotone_arr",
        action="store_false", default=True,
        help=(
            "Disable clipping of normalized poloidal flux "
            "arr=(ψ-simag)/(sibry-simag) to keep arr >= arr_floor "
            "inside the LCFS."
        )
    )
    ap.add_argument(
        "--arr_floor", type=float, default=0.0,
        help=(
            "(d. 0.0) Lower bound imposed on arr inside the LCFS when "
            "monotonicity enforcement is enabled."
        )
    )

    # -------------------------------------------------------------------------
    # GS refinement options
    # -------------------------------------------------------------------------
    ap.add_argument(
        "--gs_refine", action="store_true",
        help=(
            "Apply iterative GS-consistent SOR refinement to the smoothed "
            "psirz. This drives the reconstructed field toward Grad-Shafranov "
            "consistency while keeping boundary values fixed."
        )
    )
    ap.add_argument(
        "--gs_n_iter", type=int, default=200,
        help=(
            "(d. 200) Number of SOR iterations for GS refinement. "
            "More iterations = better GS consistency but slower."
        )
    )
    ap.add_argument(
        "--gs_omega", type=float, default=0.6,
        help=(
            "(d. 0.6) SOR relaxation parameter (0 < ω < 2). "
            "Use 0.5-0.8 for safety, 1.0-1.5 for speed. "
            "If convergence oscillates, reduce ω."
        )
    )

    # -------------------------------------------------------------------------
    # GS residual check options
    # -------------------------------------------------------------------------
    ap.add_argument(
        "--gs_check", action="store_true",
        help="Compute and print Grad-Shafranov residual diagnostic."
    )
    ap.add_argument(
        "--gs_check_plot", action="store_true",
        help="Also produce a 4-panel plot of the GS residual."
    )

    args = ap.parse_args()

    # =====================================================================
    # Read EQDSK
    # =====================================================================
    eq = read_geqdsk_full(args.eqdsk, args.width)
    R, Z, psirz = eq["R"], eq["Z"], eq["psirz"]
    raxis, zaxis = eq["raxis"], eq["zaxis"]
    rbbbs, zbbbs = eq["rbbbs"], eq["zbbbs"]

    # =====================================================================
    # Optional ψ(R,Z) smoothing (before contour extraction)
    # =====================================================================
    if args.smooth_psirz:
        psirz = smooth_psirz_axis_preserving(
            R, Z, psirz, raxis, zaxis,
            sigma_R=args.smooth_sigma_R,
            sigma_Z=args.smooth_sigma_Z,
            r0=args.smooth_r0,
            dr=args.smooth_dr
        )

    # Build a 2D interpolator for the (possibly smoothed) ψ field
    itpPsi2D = RegularGridInterpolator(
        (Z, R), psirz, bounds_error=False, fill_value=np.nan)

    # Evaluate ψ at axis and at boundary points
    psi_axis_2d = float(itpPsi2D([[zaxis, raxis]])[0])

    bb = itpPsi2D(np.column_stack([zbbbs, rbbbs]))
    bb = bb[np.isfinite(bb)]
    psi_lcfs_2d = float(np.median(bb))

    # =====================================================================
    # COCOS diagnostic / inference (non-unique in general)
    # =====================================================================
    if args.cocos_diag:
        report, candidates = infer_cocos_from_eqdsk(eq, psi_axis_2d, psi_lcfs_2d)
        print("\n=== COCOS inference diagnostics (heuristic) ===")
        for k in ["sign_Ip", "sign_B0", "sign_dpsi", "sign_pprime",
                  "pprime_expected_sign", "pprime_sign_consistent",
                  "sign_q", "sign_F"]:
            print(f"{k:>26s}: {report[k]}")
        print("COCOS candidate set (cannot usually separate k vs k+10 from EQDSK):")
        print("  ", candidates)
        if report["notes"]:
            print("Notes:")
            for n in report["notes"]:
                print("  -", n)
        print("=============================================\n")

    # =====================================================================
    # Choose target surfaces in rho_pol and compute target ψ levels
    # =====================================================================
    rho_in = np.linspace(args.rho_min, args.rho_max, args.nmhd)
    psi_targets_2d = psi_axis_2d + (psi_lcfs_2d - psi_axis_2d) * (rho_in**2)
    psi_targets_1d = psi2d_to_psi1d(
        psi_targets_2d, psi_axis_2d, psi_lcfs_2d, eq["simag"], eq["sibry"])

    # =====================================================================
    # Profile interpolators on ψ1d
    # =====================================================================
    q_of_psi = interp1d(eq["psi1d"], eq["qpsi"], kind="linear",
                        bounds_error=False,
                        fill_value=(eq["qpsi"][0], eq["qpsi"][-1]))
    fpol_of_psi = interp1d(eq["psi1d"], eq["fpol"], kind="linear",
                           bounds_error=False,
                           fill_value=(eq["fpol"][0], eq["fpol"][-1]))
    pprime_of_psi = interp1d(eq["psi1d"], eq["pprime"], kind="linear",
                             bounds_error=False,
                             fill_value=(eq["pprime"][0], eq["pprime"][-1]))
    ffprim_of_psi = interp1d(eq["psi1d"], eq["ffprim"], kind="linear",
                             bounds_error=False,
                             fill_value=(eq["ffprim"][0], eq["ffprim"][-1]))

    q_in = q_of_psi(psi_targets_1d)
    G_in = fpol_of_psi(psi_targets_1d)

    # =====================================================================
    # Compute rho_tor on requested surfaces
    # =====================================================================
    rho_tor_in = build_rho_tor_interpolator(
        eq["psi1d"], eq["qpsi"])(psi_targets_1d)

    # =====================================================================
    # Upsample ψ grid for contour extraction
    # =====================================================================
    Rg, Zg, psig = upsample_psi(R, Z, psirz, args.upsample)

    # Uniform θ samples used for Fourier reconstruction
    theta_u = np.linspace(0.0, TWOPI, args.ntheta, endpoint=False)

    # Storage arrays for Fourier coefficients
    ok = np.ones(args.nmhd, dtype=bool)
    mmax_used = []

    a0R = np.full(args.nmhd, np.nan)
    a0Z = np.full(args.nmhd, np.nan)
    aR = np.zeros((args.nmhd, args.mmax_max + 1))
    bR = np.zeros((args.nmhd, args.mmax_max + 1))
    aZ = np.zeros((args.nmhd, args.mmax_max + 1))
    bZ = np.zeros((args.nmhd, args.mmax_max + 1))
    jphi_avg = np.full(args.nmhd, np.nan)

    # =====================================================================
    # Loop over surfaces: contour extraction → smoothing → θ param → Fourier
    # =====================================================================
    for i in range(args.nmhd):
        # Extract the contour ψ = psi_targets_2d[i] that encloses the axis
        c = select_axis_enclosing_contour(
            Rg, Zg, psig, float(psi_targets_2d[i]),
            raxis, zaxis, min_points=args.min_points)
        if c is None:
            ok[i] = False
            continue

        Rc, Zc = c

        # Smooth contour coordinates with periodic Savitzky–Golay
        Rc, Zc = smooth_closed_curve_periodic_savgol(
            Rc, Zc, window=args.pre_window, polyorder=args.pre_poly)

        # Sample surface on uniform θ using the chosen poloidal angle
        theta_u, Ru, Zu = sample_surface_uniform_theta(
            Rc, Zc, raxis=raxis, zaxis=zaxis,
            ntheta=args.ntheta,
            poloidal_angle=args.poloidal_angle,
            theta0=args.theta0,
            theta_eps=args.theta_eps)

        # Compute Fourier coefficients
        if args.mmax_final and args.mmax_final > 0:
            # Fixed mmax for all surfaces
            mmax_i = args.mmax_final
            a0Ri, aRi, bRi = fourier_coeffs(Ru, mmax_i)
            a0Zi, aZi, bZi = fourier_coeffs(Zu, mmax_i)
            a0R[i] = a0Ri; a0Z[i] = a0Zi
            aR[i, :mmax_i + 1] = aRi; bR[i, :mmax_i + 1] = bRi
            aZ[i, :mmax_i + 1] = aZi; bZ[i, :mmax_i + 1] = bZi
        else:
            # Adaptive mmax via condensation
            mmax_i, pack = condensation_mmax(
                theta_u, Ru, Zu, args.mmax_max, args.condense_tol)
            a0Ri, aR_full, bR_full, a0Zi, aZ_full, bZ_full = pack
            a0R[i] = a0Ri; a0Z[i] = a0Zi
            aR[i, :len(aR_full)] = aR_full; bR[i, :len(bR_full)] = bR_full
            aZ[i, :len(aZ_full)] = aZ_full; bZ[i, :len(bZ_full)] = bZ_full

        mmax_used.append(mmax_i)

        # TORIC current proxy: <J_φ> = <R p' + FF'/R>
        pp = float(pprime_of_psi(psi_targets_1d[i]))
        ff = float(ffprim_of_psi(psi_targets_1d[i]))
        Jtheta = Ru * pp + ff / np.maximum(Ru, 1e-30)
        jphi_avg[i] = float(np.mean(Jtheta))

    if ok.sum() < args.nmhd * 0.8:
        raise RuntimeError("Too many surfaces failed contour extraction.")

    # Global nmod is the maximum mmax encountered (or forced mmax_final)
    mmax_global = (args.mmax_final if (args.mmax_final and args.mmax_final > 0)
                   else int(np.max(mmax_used)))
    nmod = mmax_global + 1
    print("global mmax:", mmax_global)

    # Truncate coefficient arrays to nmod columns
    aR = aR[:, :nmod]; bR = bR[:, :nmod]
    aZ = aZ[:, :nmod]; bZ = bZ[:, :nmod]

    # =====================================================================
    # Optional radial regularization of Fourier coefficients
    # =====================================================================
    if args.reg_smooth:
        rhoPF_max = (args.rhoPF_max if args.rhoPF_max > 0
                     else float(args.rho_max))
        rhoPF_out = np.linspace(0.0, rhoPF_max, args.nout)
        rho_reg = args.rho_reg if args.rho_reg > 0 else 0.3

        # Interpolate scalar profiles onto the new radial grid
        qsf = interp1d(rho_in, q_in, fill_value="extrapolate")(rhoPF_out)
        rho_tor = interp1d(rho_in, rho_tor_in,
                           fill_value="extrapolate")(rhoPF_out)
        G_covBphi = interp1d(rho_in, G_in,
                             fill_value="extrapolate")(rhoPF_out)
        jphi_avg_full = interp1d(rho_in, jphi_avg,
                                 fill_value="extrapolate")(rhoPF_out)

        # Regularize m=0 (even parity)
        a0R_out = ensureSmoothness_even(rho_in, a0R, rho_reg, rhoPF_out)
        a0Z_out = ensureSmoothness_even(rho_in, a0Z, rho_reg, rhoPF_out)

        rc = np.zeros((args.nout, nmod))
        rs = np.zeros((args.nout, nmod))
        zc = np.zeros((args.nout, nmod))
        zs = np.zeros((args.nout, nmod))

        rc[:, 0] = a0R_out
        zc[:, 0] = a0Z_out

        # Regularize m>=1: even m → even parity, odd m → odd parity
        for m in range(1, nmod):
            if m % 2 == 0:
                rc[:, m] = ensureSmoothness_even(rho_in, aR[:, m],
                                                 rho_reg, rhoPF_out)
                rs[:, m] = ensureSmoothness_even(rho_in, bR[:, m],
                                                 rho_reg, rhoPF_out)
                zc[:, m] = ensureSmoothness_even(rho_in, aZ[:, m],
                                                 rho_reg, rhoPF_out)
                zs[:, m] = ensureSmoothness_even(rho_in, bZ[:, m],
                                                 rho_reg, rhoPF_out)
            else:
                rc[:, m] = ensureSmoothness_odd(rho_in, aR[:, m],
                                                rho_reg, rhoPF_out)
                rs[:, m] = ensureSmoothness_odd(rho_in, bR[:, m],
                                                rho_reg, rhoPF_out)
                zc[:, m] = ensureSmoothness_odd(rho_in, aZ[:, m],
                                                rho_reg, rhoPF_out)
                zs[:, m] = ensureSmoothness_odd(rho_in, bZ[:, m],
                                                rho_reg, rhoPF_out)

        rho_pol = rhoPF_out
        print(f"Regularized output grid: rhoPF_out in [0,{rho_pol[-1]}], "
              f"nout={args.nout}, rho_reg={rho_reg}")
    else:
        # No regularization: use the extracted grid directly
        rho_pol = rho_in
        qsf = q_in
        rho_tor = rho_tor_in
        G_covBphi = G_in
        jphi_avg_full = jphi_avg

        rc = np.zeros((args.nmhd, nmod))
        rs = np.zeros((args.nmhd, nmod))
        zc = np.zeros((args.nmhd, nmod))
        zs = np.zeros((args.nmhd, nmod))

        rc[:, 0] = a0R; zc[:, 0] = a0Z
        rc[:, 1:] = aR[:, 1:]; rs[:, 1:] = bR[:, 1:]
        zc[:, 1:] = aZ[:, 1:]; zs[:, 1:] = bZ[:, 1:]

    # Report axis position consistency
    print("rmaxis =", eq["raxis"], " rc[0,0] =", rc[0, 0],
          " diff", np.abs(eq["raxis"] - rc[0, 0]))
    print("zmaxis =", eq["zaxis"], " zc[0,0] =", zc[0, 0],
          " diff", np.abs(eq["zaxis"] - zc[0, 0]))

    # Force the position of the magnetic axis to the extremum of the
    # Fourier-reconstructed innermost surface (more self-consistent)
    eq["raxis"] = rc[0, 0]
    eq["zaxis"] = zc[0, 0]

    # =====================================================================
    # Write TORIC file
    # =====================================================================
    fileout = args.eqdsk + '.equigs'
    write_toric_equi(
        filename=fileout,
        rtorus=eq["rcentr"],
        raxis=eq["raxis"],
        baxis=abs(eq["bcentr"]),
        ip_A=float(eq["current"]),
        nmod=nmod,
        nmhd=len(rho_pol),
        rho_pol=rho_pol,
        rho_tor=rho_tor,
        rc=rc, rs=rs, zc=zc, zs=zs,
        qsf=qsf,
        jphi_avg=jphi_avg_full,
        G_covBphi=G_covBphi
    )
    print(f"Wrote {fileout} (poloidal_angle={args.poloidal_angle})")

    # =====================================================================
    # Optional: write a self-consistent smoothed EQDSK
    # =====================================================================
    eq_neq = None

    if args.write_smooth_eqdsk:
        nw_out = args.smooth_eqdsk_nw if args.smooth_eqdsk_nw > 0 else None
        nh_out = args.smooth_eqdsk_nh if args.smooth_eqdsk_nh > 0 else None

        smooth_eqdsk_out = args.eqdsk + '.smooth.eqdsk'
        eq_neq, psirz_neq = build_and_write_consistent_smooth_eqdsk(
            eq=eq,
            rho_pol=rho_pol,
            rc=rc, rs=rs, zc=zc, zs=zs,
            out_path=smooth_eqdsk_out,
            ntheta_geom=max(args.ntheta, 721),
            nrho_dense=201,
            nw_out=nw_out,
            nh_out=nh_out,
            outside_mode=args.smooth_eqdsk_outside,
            width=args.width,
            flip_q_sign=args.flip_q_sign,
            flip_ffprim_too=args.flip_ffprim_too,
            title='Smooth ' + args.eqdsk + ' ' + str(date.today()),
            enforce_monotone_arr=args.enforce_monotone_arr,
            arr_floor=args.arr_floor,
            gs_refine=args.gs_refine,
            gs_n_iter=args.gs_n_iter,
            gs_omega=args.gs_omega
        )

        print(f"\nWrote self-consistent smoothed EQDSK: {smooth_eqdsk_out}")
        print("  new simag =", eq_neq["simag"])
        print("  new sibry =", eq_neq["sibry"])
        print("  new nw,nh =", eq_neq["nw"], eq_neq["nh"])

        # Monotonicity diagnostic: report the minimum normalized flux
        dpsi_diag = eq_neq["sibry"] - eq_neq["simag"]
        if abs(dpsi_diag) > 1e-300:
            arr_diag = ((eq_neq["psirz"] - eq_neq["simag"]) / dpsi_diag)
            idx = np.unravel_index(np.argmin(arr_diag), arr_diag.shape)
            print(f"  min(arr)  = {arr_diag[idx]:.6e}")
            print(f"  psi(min)  = {eq_neq['psirz'][idx]:.6e}")
            # idx[0] is Z-index, idx[1] is R-index (Z-major storage)
            print(f"  R,Z(min)  = {eq_neq['R'][idx[1]]:.4f}, "
                  f"{eq_neq['Z'][idx[0]]:.4f}")

    # =====================================================================
    # Optional: GS residual check
    # =====================================================================
    if args.gs_check:
        print("\n --- GS check on ORIGINAL EQDSK ---")
        compute_gs_residual(
            eq, verbose=True,
            plot=args.gs_check_plot,
            plot_file=(args.eqdsk + ".gs_orig.png"
                       if args.gs_check_plot else None))

        if eq_neq is not None:
            print("--- GS check on SMOOTHED EQDSK ---")
            compute_gs_residual(
                eq_neq, verbose=True,
                plot=args.gs_check_plot,
                plot_file=(args.eqdsk + ".gs_smooth.png"
                           if args.gs_check_plot else None))

    # =====================================================================
    # Plot reconstruction check (optional)
    # =====================================================================
    if args.plot:
        rho_solid = np.linspace(rho_pol[0], rho_pol[-1], args.plot_n_solid)

        # Build dashed contour rho values between solid surfaces
        dashed_rhos = []
        if args.plot_dashed_between > 0:
            for i in range(len(rho_solid) - 1):
                a = rho_solid[i]
                b = rho_solid[i + 1]
                for k in range(1, args.plot_dashed_between + 1):
                    dashed_rhos.append(
                        a + (b - a) * k / (args.plot_dashed_between + 1))
        dashed_rhos = np.array(dashed_rhos)
        dashed_psis = (psi_axis_2d
                       + (psi_lcfs_2d - psi_axis_2d) * (dashed_rhos**2))

        # ---------------------------------------------------------------------
        # Figure layout:
        #   left column  : geometry (spans both rows)
        #   right column : top = q(rho), bottom = rho_tor(rho_pol)
        # ---------------------------------------------------------------------
        fig = plt.figure(figsize=(13.5, 7.6))
        gs_layout = fig.add_gridspec(
            2, 2, width_ratios=[3.2, 1.4], height_ratios=[1, 1],
            wspace=0.28, hspace=0.28)

        ax = fig.add_subplot(gs_layout[:, 0])    # geometry spans both rows
        axq = fig.add_subplot(gs_layout[0, 1])   # top-right
        axr = fig.add_subplot(gs_layout[1, 1])   # bottom-right

        # =====================================================================
        # Left panel: geometry reconstruction
        # =====================================================================
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("R [m]"); ax.set_ylabel("Z [m]")
        ax.grid(True, alpha=0.25)
        ax.set_title(f"Angle {args.poloidal_angle}: "
                     f"solid=recon, dashed=EQDSK, spokes=equal θ")

        # Plot LCFS and axis
        ax.plot(eq["rbbbs"], eq["zbbbs"], "k-", lw=2, label="LCFS")
        ax.plot([raxis], [zaxis], "rx", ms=8, mew=2, label="axis")

        # Dashed EQDSK contours (intermediate surfaces)
        for psi_lev in dashed_psis:
            c = select_axis_enclosing_contour(
                Rg, Zg, psig, float(psi_lev),
                raxis, zaxis, min_points=args.min_points)
            if c is not None:
                ax.plot(c[0], c[1], "k--", lw=0.9, alpha=0.6)

        # Solid reconstructed TORIC surfaces
        theta_plot = np.linspace(0, TWOPI, args.ntheta, endpoint=False)
        Rsol = []; Zsol = []
        for rp in rho_solid:
            def lincoef(arr):
                return np.interp(rp, rho_pol, arr)

            Rrec = np.full_like(theta_plot, lincoef(rc[:, 0]), float)
            Zrec = np.full_like(theta_plot, lincoef(zc[:, 0]), float)
            for m in range(1, nmod):
                Rrec += (lincoef(rc[:, m]) * np.cos(m * theta_plot)
                         + lincoef(rs[:, m]) * np.sin(m * theta_plot))
                Zrec += (lincoef(zc[:, m]) * np.cos(m * theta_plot)
                         + lincoef(zs[:, m]) * np.sin(m * theta_plot))
            ax.plot(Rrec, Zrec, "C0-", lw=1.2, alpha=0.95)
            Rsol.append(Rrec); Zsol.append(Zrec)

        Rsol = np.array(Rsol); Zsol = np.array(Zsol)

        # Constant-θ spokes across radius
        if args.theta_spokes > 0:
            thetas_spoke = np.linspace(0, TWOPI, args.theta_spokes,
                                       endpoint=False)
            idx = ((thetas_spoke / TWOPI * args.ntheta).astype(int)
                   % args.ntheta)
            for j in idx:
                ax.plot(Rsol[:, j], Zsol[:, j], color="C1", lw=0.9,
                        alpha=0.85)

        ax.legend(loc="best")

        # =====================================================================
        # Top-right panel: q(rho)
        # =====================================================================
        # Original EQDSK q(psi) mapped to rho_pol
        psi1d_eq = np.linspace(eq["simag"], eq["sibry"], len(eq["qpsi"]))
        rho_eq = np.sqrt(np.clip(
            (psi1d_eq - eq["simag"]) / (eq["sibry"] - eq["simag"] + 1e-300),
            0.0, 1.0))
        axq.plot(rho_eq, eq["qpsi"], "k-", lw=1.8, label="qpsi original")

        # Smoothed EQDSK q(psi) if available
        if eq_neq is not None and "qpsi" in eq_neq:
            psi1d_neq = np.linspace(
                eq_neq["simag"], eq_neq["sibry"], len(eq_neq["qpsi"]))
            rho_neq = np.sqrt(np.clip(
                (psi1d_neq - eq_neq["simag"])
                / (eq_neq["sibry"] - eq_neq["simag"] + 1e-300), 0.0, 1.0))
            axq.plot(rho_neq, eq_neq["qpsi"], "C2--", lw=1.8,
                     label="qpsi smoothed")

        # q profile written to TORIC
        axq.plot(rho_pol, qsf, "C3-", lw=2.2, label="qsf TORIC")

        # Mark the same solid surfaces
        for rp in rho_solid:
            axq.plot([rp], [np.interp(rp, rho_pol, qsf)], "ko", ms=3)

        axq.set_xlabel(r"$\rho_{\rm pol}$"); axq.set_ylabel("q")
        axq.set_title("Safety factor profile")
        axq.grid(True, alpha=0.3)
        axq.set_xlim(rho_pol[0], rho_pol[-1])
        axq.legend(loc="best", fontsize=8)

        # =====================================================================
        # Bottom-right panel: rho_tor(rho_pol)
        # =====================================================================
        axr.plot(rho_pol, rho_tor, "C4-", lw=2.2,
                 label=r"$\rho_{\rm tor}(\rho_{\rm pol})$")

        for rp in rho_solid:
            axr.plot([rp], [np.interp(rp, rho_pol, rho_tor)], "ko", ms=3)

        # Reference diagonal
        axr.plot([rho_pol[0], rho_pol[-1]], [rho_pol[0], rho_pol[-1]],
                 color="0.6", ls="--", lw=1.0,
                 label=r"$\rho_{\rm tor}=\rho_{\rm pol}$")

        axr.set_xlabel(r"$\rho_{\rm pol}$")
        axr.set_ylabel(r"$\rho_{\rm tor}$")
        axr.set_title(r"Toroidal vs poloidal radius")
        axr.grid(True, alpha=0.3)
        axr.set_xlim(rho_pol[0], rho_pol[-1])
        axr.legend(loc="best", fontsize=8)

        plt.tight_layout()
        plt.savefig(args.eqdsk + ".png", dpi=300, bbox_inches="tight")
        plt.show()


if __name__ == "__main__":
    main()
