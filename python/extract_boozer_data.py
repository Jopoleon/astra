import sys, os, traceback
import numpy as np
from scipy.interpolate import PchipInterpolator
from parse_fortran_nml import parse_fortran_namelist

def astra_rho(n_rho):
    hrox = 1.0 / (n_rho - 0.5)
    grid = np.arange(n_rho, dtype=np.float32) + 0.5*hrox
    return grid


def pchip_clamped(rho_anchor, b00_anchor, astra_rho_grid):
    order = np.argsort(rho_anchor)
    rho_anchor = rho_anchor[order]
    b00_anchor = b00_anchor[order]
    f = PchipInterpolator(rho_anchor, b00_anchor, extrapolate=False)
    q = np.clip(astra_rho_grid, rho_anchor[0], rho_anchor[-1])
    return f(q)


# ARCHIVE PATH  (text Boozer, all surfaces, dense)

def read_boozer_b00_text(boozer_file, astra_rho_grid):
    """Archive text Boozer, all surfaces. """
    b00_values, s_values, current_s = [], [], None
    if not os.path.exists(boozer_file):
        raise FileNotFoundError(f"Boozer file not found: {boozer_file}")
    with open(boozer_file, 'r') as f:
        while True:
            line = f.readline()
            if not line:
                break
            if "m0b" in line and "n0b" in line:
                break
        for line in f:
            parts = line.split()
            if not parts:
                continue
            try:
                val1 = float(parts[0])
                if "." in parts[0] or "E" in parts[0].upper():
                    current_s = val1
                    continue
            except ValueError:
                continue
            if len(parts) >= 6:
                try:
                    m, n = int(parts[0]), int(parts[1])
                    if m == 0 and n == 0 and current_s is not None:
                        s_values.append(current_s)
                        b00_values.append(float(parts[5]))
                except ValueError:
                    continue

    s = np.array(s_values)
    b = np.array(b00_values)
    if len(s) == 0:
        raise ValueError("No B00 data found in text Boozer file, check format")

    order = np.argsort(s)
    s = s[order]
    b = b[order]

    # axis point
    if s[0] > 1e-5:
        s = np.insert(s, 0, 0.0)
        if len(b) >= 2:
            slope = (b[1] - b[0]) / (s[1] - s[0])
            b = np.insert(b, 0, b[0] - slope * s[0])
        else:
            b = np.insert(b, 0, b[0])

    rho = np.sqrt(s)
    return pchip_clamped(rho, b, astra_rho_grid)


def read_minor_radius_text(boozer_file):
    """Read W7AS minor radius from the text Boozer header (6th value on the
    values line that follows the 'm0b n0b' descriptor line)."""
    def sscan(strng):
        return np.fromstring(strng, float, -1, sep=' ')
    with open(boozer_file, 'r') as f:
        for line in f:
            if "m0b" in line and "n0b" in line:
                header = sscan(f.readline())
                if len(header) > 5:
                    return float(header[5])
                raise ValueError(
                    f"header line after 'm0b n0b' has only {len(header)} "
                    f"values; expected minor radius at index 5")
    raise ValueError("'m0b n0b' descriptor line not found in Boozer file")


# GENERATED PATH  (NetCDF boozmn, 7 surfaces from true_surfaces.txt)
def _read_true_surfaces(path="vmec_io/true_surfaces.txt"):
    if not os.path.exists(path):
        raise FileNotFoundError(f"{path} not found (needed for generated mode)")
    with open(path) as f:
        next(f)                        # line 1 = count
        return [int(x) for x in f.read().split()]


def read_boozer_b00_nc(boozer_nc, wout_nc, astra_rho_grid):
    import xarray as xr
    dsb = xr.open_dataset(boozer_nc)
    dsv = xr.open_dataset(wout_nc)

    m_arr = dsb["ixm_b"].values
    n_arr = dsb["ixn_b"].values
    bmnc  = dsb["bmnc_b"].values            # (pack_rad, mn_modes)
    idx00 = np.where((m_arr == 0) & (n_arr == 0))[0][0]
    b00_rows = bmnc[:, idx00]               # one value per transformed surface

    good = b00_rows > 1.0e-8
    if good.sum() < 2:
        raise ValueError(f"only {good.sum()} non-zero B00 surfaces in boozmn; "
                         f"Boozer transform produced no usable data")
    b00_rows = b00_rows[good]

    true_surf = _read_true_surfaces()        # VMEC indices of computed surfaces
    if len(true_surf) != len(b00_rows):
        raise ValueError(f"true_surfaces ({len(true_surf)}) != nonzero boozer "
                         f"rows ({len(b00_rows)})")

    # half-grid s
    phi = np.abs(dsv["phi"].values).squeeze()
    s_vmec = phi / phi[-1]
    s_booz = s_vmec.copy()
    s_booz[0] = 0.0
    s_booz[1:] = 0.5 * (s_vmec[1:] + s_vmec[:-1])

    idx = np.array([t - 1 for t in true_surf])   # 0-based VMEC index
    rho_anchor = np.sqrt(s_booz[idx])

    return pchip_clamped(rho_anchor, b00_rows, astra_rho_grid)


def minor_radius_modesum(wout_nc):
    import xarray as xr
    ds = xr.open_dataset(wout_nc)
    xm = ds["xm"].values
    R = ds["rmnc"].values[-1, :]
    Z = ds["zmns"].values[-1, :]
    a2 = 0.0
    for i in range(len(xm)):
        if xm[i] == 0:
            continue
        a2 += xm[i] * R[i] * Z[i]
    return float(np.sqrt(a2))


def is_netcdf(path):
    """True if the file is a NetCDF (classic 'CDF' or HDF5-based '.nc').
    Archive Boozer files are plain text, generated boozmn are NetCDF."""
    try:
        with open(path, 'rb') as f:
            magic = f.read(4)
    except OSError:
        print('OSerror')
        return False
    print(magic[:3])
    return magic[:3] == b'CDF' or magic == b'\x89HDF'


if __name__ == "__main__":

    nl_vmec = parse_fortran_namelist('vmec_io/stell_files.nml', 'VMEC_TO_ASTRA_INPUTS')
    nl_dkes = parse_fortran_namelist('vmec_io/stell_files.nml', 'ASTRA_DKES_INTERFACE')

    # NA1: optional command-line arg (a2vmec passes it) overrides the namelist
    # default; the interface reads B00_PHYSICAL_PROFILE(NA1), so the profile
    # length MUST equal NA1.
    if len(sys.argv) > 1:
        npts = int(sys.argv[1])
    else:
        npts = nl_vmec.get('astra_nrad', 91)
    grid = astra_rho(npts)

    boozer_file = nl_vmec['boozer_file']
    wout_file   = nl_vmec['vmec_wout_file']

    if not os.path.exists(boozer_file):
        print('boozer_file %s not found' %boozer_file)
        input('Press <Enter> to continue')
    if not is_netcdf(boozer_file):
        print('boozer_file %s is not NetCDF nor HDF5' %boozer_file)
        input('Press <Enter> to continue')

    b00_out = "dat/b00_profile_boozer.txt"
    rad_out = "dat/minorradiusW7AS.txt"

    try:
        mode = "generated"
        b00 = read_boozer_b00_nc(boozer_file, wout_file, grid)
        rad = minor_radius_modesum(wout_file)
        np.savetxt(b00_out, b00)
        with open(rad_out, 'w') as f:
            f.write(str(rad))
        print(f"[{mode}] wrote {b00_out}  (first 5: {b00[:5]})")
        print(f"[{mode}] wrote {rad_out}  (minorradiusW7AS = {rad})")

    except Exception as e:
        traceback.print_exc()
        sys.stderr.write(
            "\n" + "!" * 70 + "\n"
            f"!! extract_boozer_data FAILED: {e}\n"
            f"!! Boozer file: {boozer_file}\n"
            "!! Writing ZERO B00 and ZERO minor radius as placeholders.\n"
            "!! The DKES interface CANNOT produce a valid Er from these.\n"
            "!! Fix the Boozer file / namelist path before trusting results.\n"
            + "!" * 70 + "\n\n")
        np.savetxt(b00_out, np.zeros(npts))
        with open(rad_out, 'w') as f:
            f.write("0.0")
