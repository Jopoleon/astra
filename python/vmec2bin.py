import sys, os, re, argparse, logging, time
import numpy as np
from scipy.io import netcdf_file
from scipy.interpolate import interp1d
from parse_fortran_nml import parse_fortran_namelist


fmt = logging.Formatter('%(asctime)s | %(name)s | %(levelname)s: %(message)s', '%H:%M:%S')
logger = logging.getLogger('vmec2a')
if len(logger.handlers) == 0:
    hnd = logging.StreamHandler()
    hnd.setFormatter(fmt)
    logger.addHandler(hnd)
#logger.setLevel(logging.DEBUG)
logger.setLevel(logging.INFO)

awd = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

# Interpolation
linterp = lambda x, y, xi: interp1d(x, y, kind='linear'   , fill_value='extrapolate')(xi)
qinterp = lambda x, y, xi: interp1d(x, y, kind='quadratic', fill_value='extrapolate')(xi)


def h2f(var_half):
    """Half to full grid"""
    var_full = np.empty_like(var_half)
    var_full[0]  = 1.5 * var_half[1] - 0.5 * var_half[2]
    var_full[1:-1] = 0.5 * (var_half[1:-1] + var_half[2:])
    var_full[-1] = 2. * var_half[-1] - var_full[-2]
    return var_full


def f2h(var_full):
    """Full to half grid"""
    var_half = np.zeros_like(var_full)
    var_half[1] = 0.5 * (var_full[0] + var_full[1])
    for i in range(2, len(var_full)):
        var_half[i] = 2. * var_full[i-1] - var_half[i-1]
    return var_half


def h2fmn(var_half, mmode, ns):
    """Half to full grid with Fourier interpolation

    This routine takes a 1D field and interpolates it from the half
    to the full grid taking care of even and odd mode parity. For 
    an ns sized array we assumes that the first index [0]=0 and 
    is just a placeholder.

    Parameters
    ----------
    var_half : list
    	Variable on half grid
    mmode : int
    	Poloidal mode number
    Returns
    ----------
    var_full : list
    	Variable on full grid
    """

    temp = var_half.copy()
    if np.mod(mmode, 2) == 1:
        factlo = 0.5*np.sqrt(np.linspace(0, ns-1.0, ns)/np.linspace(-0.5, ns-1.5, ns))
        facthi = 0.5*np.sqrt(np.linspace(0, ns-1.0, ns)/np.linspace( 0.5, ns-0.5, ns))
        temp[1:-1] = factlo[1:-1]*temp[1:-1] + facthi[2:]*temp[2:]
        factlo = 2.*np.sqrt((ns-1)/(ns - 1.5))
        facthi =   -np.sqrt((ns-1)/(ns - 2.0))
        temp[-1] = factlo*temp[-1] + facthi*temp[-2]
        temp[0] = 0.0
    else:
        temp[1:-1] = 0.5*(temp[1:-1] + temp[2:])
        temp[-1] = 2.*temp[-1] - temp[-2]
        temp[0]  = 2.*temp[ 1] - temp[ 2]
    return temp


def c_funct(theta, phi, fmnc, xm, xn):
    (ns, mn) = fmnc.shape
    lt = len(theta)
    lz = len(phi)
    mt = np.outer(xm, theta)
    nz = np.outer(xn, phi)
    cosmt = np.cos(mt)
    sinmt = np.sin(mt)
    cosnz = np.cos(nz)
    sinnz = np.sin(nz)

    f = np.empty((ns, lt, lz))
    for k in range(ns):
        coeff = fmnc[k, :]
        f[k] = (
            (coeff[:, None] * cosmt).T @ cosnz -
            (coeff[:, None] * sinmt).T @ sinnz )

    return f


def s_funct(theta, phi, fmnc, xm, xn):
    (ns, mn) = fmnc.shape
    lt = len(theta)
    lz = len(phi)
    mt = np.outer(xm, theta)
    nz = np.outer(xn, phi)
    cosmt = np.cos(mt)
    sinmt = np.sin(mt)
    cosnz = np.cos(nz)
    sinnz = np.sin(nz)

    f = np.empty((ns, lt, lz))
    for k in range(ns):
        coeff = fmnc[k, :]
        f[k] = (
            (coeff[:, None] * sinmt).T @ cosnz +
            (coeff[:, None] * cosmt).T @ sinnz )

    return f


class VMEC():


    def __init__(self, wout_file):
        self.readWout(wout_file)


    def readWout(self, wout_file):
        with netcdf_file(wout_file, 'r', mmap=False) as f:
            cv = f.variables
        for label in ('raxis_cc', 'zaxis_cs', 'b0', 'rbtor', 'rmnc', 'zmns',
                      'phi', 'phipf', 'iotaf', 'jdotb', 'xm', 'xm_nyq'):
            setattr(self, label, cv[label].data.astype(np.float64))
        for label in ('nfp', 'ns', 'mnmax', 'mnmax_nyq'):
            setattr(self, label, cv[label].data.astype(np.int32))
        for label in ('bvco', 'vp'):
            setattr(self, label, h2f(cv[label].data).astype(np.float64))
        self.lasym = cv['lasym__logical__'].data
        if self.lasym:
            self.rmns = cv['rmns'].data.astype(np.float64)
            self.zmnc = cv['zmnc'].data.astype(np.float64)

        self.rmajor = cv['Rmajor_p'].data.astype(np.float64)
        self.aminor = cv['Aminor_p'].data.astype(np.float64)
        self.volume = cv['volume_p'].data.astype(np.float64)
        self.xn     = -cv['xn'    ].data.astype(np.float64)
        self.xn_nyq = -cv['xn_nyq'].data.astype(np.float64)

        self.lmns = np.zeros_like(cv['lmns'].data.astype(np.float64))
        self.bmnc = np.zeros_like(cv['bmnc'].data.astype(np.float64))
        self.gmnc = np.zeros_like(cv['gmnc'].data.astype(np.float64))
        for mn in range(self.mnmax):
            self.lmns[:, mn] = h2fmn(cv['lmns'][:, mn], int(self.xm[mn]), self.ns)
        for mn in range(self.mnmax_nyq):
            self.bmnc[:, mn] = h2fmn(cv['bmnc'][:, mn], int(self.xm_nyq[mn]), self.ns)
            self.gmnc[:, mn] = h2fmn(cv['gmnc'][:, mn], int(self.xm_nyq[mn]), self.ns)


    def calcMoms(self, nu=64, nv=128):
        t1 = time.time()
        theta = np.linspace(0, 2*np.pi, nu)
        zeta  = np.linspace(0, 2*np.pi, nv)
        self.theta = theta
        self.zeta  = zeta
# Create derivatives
        rumns = - self.rmnc * self.xm
        rvmns = - self.rmnc * self.xn
        zumnc =   self.zmns * self.xm
        zvmnc =   self.zmns * self.xn
        lumnc =   self.lmns * self.xm
        lvmnc =   self.lmns * self.xn
        t2 = time.time()
        self.r  = c_funct(theta, zeta, self.rmnc, self.xm, self.xn)
        self.g  = c_funct(theta, zeta, self.gmnc, self.xm_nyq, self.xn_nyq)
        self.b  = c_funct(theta, zeta, self.bmnc, self.xm_nyq, self.xn_nyq)
        self.ru = s_funct(theta, zeta, rumns, self.xm, self.xn)
        self.rv = s_funct(theta, zeta, rvmns, self.xm, self.xn)
        self.zu = c_funct(theta, zeta, zumnc, self.xm, self.xn)
        self.zv = c_funct(theta, zeta, zvmnc, self.xm, self.xn)
        self.lu = c_funct(theta, zeta, lumnc, self.xm, self.xn)
        self.lv = c_funct(theta, zeta, lvmnc, self.xm, self.xn)
        if self.lasym:
            rumnc =   self.rmns * self.xm
            rvmnc =   self.rmns * self.xn
            zumns = - self.zmnc * self.xm
            zvmns = - self.zmnc * self.xn
            lumns = - self.lmnc * self.xm
            lvmns = - self.lmnc * self.xn
            self.r  += s_funct(theta, zeta, self.rmns, self.xm, self.xn)
            self.g  += s_funct(theta, zeta, self.gmns, self.xm_nyq, self.xn_nyq)
            self.b  += s_funct(theta, zeta, self.bmns, self.xm_nyq, self.xn_nyq)
            self.ru += c_funct(theta, zeta, rumnc, self.xm, self.xn)
            self.rv += c_funct(theta, zeta, rvmnc, self.xm, self.xn)
            self.zu += s_funct(theta, zeta, zumns, self.xm, self.xn)
            self.zv += s_funct(theta, zeta, zvmns, self.xm, self.xn)
            self.lu += s_funct(theta, zeta, lumns, self.xm, self.xn)
            self.lv += s_funct(theta, zeta, lvmns, self.xm, self.xn)
        t3 = time.time()
        self.g_sum = np.sum(self.g, axis=(1, 2))
        t4 = time.time()
        print(self.lasym, t2-t1, t3-t2, t4-t3)
#        input('Time')


    def surfAverage(self, arr_in):
        arr_out = np.sum(arr_in * self.g, axis=(1, 2)) / self.g_sum
        arr_out[0] = 2. * arr_out[1] - arr_out[2]
        return arr_out


    def calcGrad_rho(self):
# Calc metrics
        gsr = -self.zu * self.r
        gsp =  self.zu * self.rv - self.ru * self.zv
        gsz =  self.ru * self.r
        gs  = ( gsr**2 + gsp**2 + gsz**2) / self.g**2
        rho_tor = np.sqrt(self.phi/self.phi[-1])
        gs_rho2 = 0.25 * gs / rho_tor[:, None, None]**2
        gs_sqrt = np.sqrt(gs_rho2)
        self.avg_grad_rho2 = self.surfAverage(gs_rho2)
        self.avg_grad_rho  = self.surfAverage(gs_sqrt)


    def calcRmaj(self):
        """Returns <R>"""
        self.Rmaj = self.surfAverage(self.r)


    def calcVol(self, vp):
        """Returns volume V(s)"""
        dels = (self.phi[1] - self.phi[0])/self.phi[-1]
        nrho = len(self.phi)
        self.Vol = np.append(0., 4.*dels*np.pi**2 * np.cumsum(vp[1: nrho]))


    def calcFtrap(self):
        """Compute trapped fraction f_t according to 
        H. Maassberg; C. D. Beidler; Y. Turkin; Phys. Plasmas 16, 072504 (2009) equation 12

        Returns
        ----------
        ftrap : ndarray
        """

        nlambda = 256
        b_sq = self.b**2
        b2 = self.surfAverage(b_sq)              # <B^2>
        b2_bmax2 = b2/np.max(b_sq, axis=(1, 2))  # <B^2/Bmax^2>
        b_bmax   = self.b/np.max(self.b, axis=(1, 2), keepdims=True)
# Introduce normalised global magnetic moment lambda
        dlambda = 1./(nlambda - 1.)   #stepwidth in lambda
        integrand = np.empty((nlambda, len(self.g_sum)))    #lambda/<sqrt(1-lambda*b_bmax)>
        for mn in range(nlambda):
            integrand[mn, :] = mn*dlambda*self.g_sum/np.sum(np.sqrt(1 - mn*dlambda*b_bmax)*self.g, axis=(1, 2))
        integral = np.sum(integrand, axis=0)*dlambda   #integral over lambda
        self.ftrap = 1. - 0.75*b2_bmax2*integral


    def calcSusceptance(self):

        scale_fact = 1.0 / ( 4 * np.pi**2)
        S11 = ( self.ru**2 + self.zu**2)
        S21 = ( self.ru * self.rv + self.zu * self.zv)
        S12 = ( S21 * (1. + self.lu) - S11 * self.lv )
        S22 = ( (self.rv**2 + self.zv**2 + self.r**2) * (1.0 + self.lu) - S21*self.lv )
        S11 = np.trapezoid(S11/self.g, x=self.zeta, axis=2)
        S12 = np.trapezoid(S12/self.g, x=self.zeta, axis=2)
        S21 = np.trapezoid(S21/self.g, x=self.zeta, axis=2)
        S22 = np.trapezoid(S22/self.g, x=self.zeta, axis=2)
        S11 = np.trapezoid(S11, x=self.theta, axis=1)*scale_fact
        S12 = np.trapezoid(S12, x=self.theta, axis=1)*scale_fact
        S21 = np.trapezoid(S21, x=self.theta, axis=1)*scale_fact
        S22 = np.trapezoid(S22, x=self.theta, axis=1)*scale_fact

        return S11, S12, S21, S22

# ---------------------------------------------------------------------------
# Pointwise VMEC geometry + straight-line pellet chord
# ---------------------------------------------------------------------------
# The flux-surface-averaged transforms above give ASTRA its 1D metric.  The NGS
# pellet model instead needs the LOCAL geometry along the real straight line the
# pellet flies: at each (s, theta, zeta) it wants modB and the contravariant
# basis grad_s / grad_theta / grad_zeta, so a fixed Cartesian velocity can be
# projected onto the flux coordinates (ds/dl = vhat . grad_s, ...).  This is
# exactly what simsopt's vmec_compute_geometry returns, but here it is done in
# pure numpy from the same wout Fourier data, reusing this warm worker instead
# of pulling in simsopt/pySTEL.  See docu/ngs_pellet_todo.md, section 3.
#
# VMEC real-space embedding (stellarator symmetric): with zeta the cylindrical
# toroidal angle,  x = (R cos zeta, R sin zeta, Z), 
#   R = sum rmnc(s) cos(m theta + n zeta),   Z = sum zmns(s) sin(m theta + n zeta), 
#   modB = sum bmnc(s) cos(m_nyq theta + n_nyq zeta)   (Nyquist mode set).
# The covariant basis e_s, e_theta, e_zeta -> contravariant via the reciprocal
# relations; sqrt(g) = e_s . (e_theta x e_zeta).
# ---------------------------------------------------------------------------

class PointGeom(object):
    """Evaluate R, Z, modB and grad_s/theta/zeta at an arbitrary (s, theta, zeta).

    Coefficients are linearly interpolated on the uniform full-mesh s grid
    s = linspace(0, 1, ns); the s-derivatives d(rmnc)/ds, d(zmns)/ds are taken
    once by np.gradient.  modB uses the Nyquist mode set on the same s grid
    (half/full-mesh offset of half a cell is neglected -- consistent with the
    trapped-fraction transform above and small in the ablation region)."""

    def __init__(self, vmc):
        self.xm   = vmc.xm
        self.xn   = vmc.xn
        self.xmn  = vmc.xm_nyq
        self.xnn  = vmc.xn_nyq
        self.rmnc = vmc.rmnc
        self.zmns = vmc.zmns
        self.bmnc = vmc.bmnc
        self.iasym = getattr(vmc, 'iasym', 0)
        if self.iasym == 1:
            self.rmns = vmc.rmns
            self.zmnc = vmc.zmnc
            self.bmns = vmc.bmns
        self.ns = self.rmnc.shape[0]
        self.sgrid = np.linspace(0., 1., self.ns)

    def _bracket(self, s):
        s = min(max(s, 0.), 1.)
        x = s * (self.ns - 1)
        j = int(np.floor(x))
        if j >= self.ns - 1:
            j = self.ns - 2
        return j, x - j

    def _interp_row(self, arr, s):
        """Linear interp of arr[ns, mn] at scalar s on the uniform s grid."""
        j, f = self._bracket(s)
        return arr[j] * (1. - f) + arr[j + 1] * f

    def _deriv_row(self, arr, s):
        """d/ds of the linearly-interpolated arr[ns, mn] (exact per-cell slope, 
        self-consistent with _interp_row so grad_i.e_j = delta_ij holds)."""
        j, _ = self._bracket(s)
        return (arr[j + 1] - arr[j]) * (self.ns - 1)

    def eval(self, s, theta, zeta):
        """Return dict with R, Z, modB and grad_s/grad_t/grad_z (each 3-vectors)."""
        rc  = self._interp_row(self.rmnc, s)
        zc  = self._interp_row(self.zmns, s)
        rc_s = self._deriv_row(self.rmnc, s)
        zc_s = self._deriv_row(self.zmns, s)
        a = self.xm*theta + self.xn*zeta
        ca, sa = np.cos(a), np.sin(a)

        R = np.dot(rc, ca)
        Z = np.dot(zc, sa)
        R_t = np.dot(rc, -self.xm * sa)
        R_z = np.dot(rc, -self.xn * sa)
        Z_t = np.dot(zc,  self.xm * ca)
        Z_z = np.dot(zc,  self.xn * ca)
        R_s = np.dot(rc_s, ca)
        Z_s = np.dot(zc_s, sa)

        bc = self._interp_row(self.bmnc, s)
        an = self.xmn * theta + self.xnn * zeta
        modB = np.dot(bc, np.cos(an))

        if self.iasym == 1:
            rsc  = self._interp_row(self.rmns, s)
            zcc  = self._interp_row(self.zmnc, s)
            rsc_s = self._deriv_row(self.rmns, s)
            zcc_s = self._deriv_row(self.zmnc, s)
            R += np.dot(rsc, sa)
            Z += np.dot(zcc, ca)
            R_t += np.dot(rsc,  self.xm*ca)
            R_z += np.dot(rsc,  self.xn*ca)
            Z_t += np.dot(zcc, -self.xm*sa)
            Z_z += np.dot(zcc, -self.xn*sa)
            R_s += np.dot(rsc_s, sa)
            Z_s += np.dot(zcc_s, ca)
            bsc = self._interp_row(self.bmns, s)
            modB += np.dot(bsc, np.sin(an))

        cz, sz = np.cos(zeta), np.sin(zeta)
        # covariant basis vectors in Cartesian (x = R cos z, R sin z, Z)
        e_s = np.array([R_s*cz,          R_s*sz,          Z_s])
        e_t = np.array([R_t*cz,          R_t*sz,          Z_t])
        e_z = np.array([R_z*cz - R*sz,   R_z *sz + R*cz,  Z_z])
        sqrtg = np.dot(e_s, np.cross(e_t, e_z))
        if abs(sqrtg) < 1e-30:
            sqrtg = 1e-30
        grad_s = np.cross(e_t, e_z) / sqrtg
        grad_t = np.cross(e_z, e_s) / sqrtg
        grad_z = np.cross(e_s, e_t) / sqrtg
        return dict(R=R, Z=Z, modB=abs(modB), grad_s=grad_s, grad_t=grad_t, grad_z=grad_z)


def compute_pellet_chord(vmc, launch_theta, launch_phi, v_xyz, max_len=None, nstep=2000):
    """Trace a straight real-space line from the LCFS launch point along v_xyz
    and tabulate the geometry the Fortran ablation ODE needs.

    The pellet flies with fixed Cartesian velocity v_xyz [m/s]; along the line
    the flux coordinates evolve by  d(s, theta, zeta)/dl = vhat . grad(s, theta, zeta), 
    integrated with RK4 in real arclength l (vhat = v/|v|).  Integration starts
    at s=1, (theta, zeta)=(launch_theta, launch_phi) and stops when the pellet
    exits (s back >= 1 after entering) or after max_len.

    Returns arrays (l, rho, s, modB, R, Z), l ascending from the launch point.
    rho = sqrt(s) is the ASTRA flux label.
    """

    pg = PointGeom(vmc)
    v = np.asarray(v_xyz, dtype=np.float64)
    vmag = np.linalg.norm(v)
    if vmag <= 0.:
        raise ValueError("pellet velocity magnitude is zero")
    vhat = v/vmag

# a generous default path length: a few times the device size

    if max_len is None:
        max_len = 4.*vmc.aminor
    dl = max_len/nstep

    def rhs(y):
        s, th, ze = y
        g = pg.eval(s, th, ze)
        return np.array([np.dot(vhat, g['grad_s']), 
                         np.dot(vhat, g['grad_t']), 
                         np.dot(vhat, g['grad_z'])])

    y = np.array([1., launch_theta, launch_phi])
    ss  = [y[0]]
    ths = [y[1]]
    zes = [y[2]]
    entered = False
    for i in range(nstep):
        k1 = rhs(y)
        k2 = rhs(y + 0.5 * dl * k1)
        k3 = rhs(y + 0.5 * dl * k2)
        k4 = rhs(y + dl * k3)
        y += (dl / 6.) * (k1 + 2*k2 + 2*k3 + k4)
# clamp s to the valid [0, 1]; reflect tiny overshoot at the axis
        if y[0] < 0.:
            y[0] = -y[0]
        ss.append(y[0])
        ths.append(y[1])
        zes.append(y[2])
        if y[0] < 0.98:
            entered = True
        if entered and y[0] >= 1.:
            break

    ls = dl*np.arange(len(ths))
    ss = np.clip(np.array(ss), 0., 1.)
    ths = np.array(ths)
    zes = np.array(zes)
    rho = np.sqrt(ss)
    modB = np.empty_like(ls)
    R = np.empty_like(ls)
    Z = np.empty_like(ls)
    for i in range(len(ls)):
        g = pg.eval(ss[i], ths[i], zes[i])
        modB[i] = g['modB']
        R[i] = g['R']
        Z[i] = g['Z']
    return ls, rho, ss, modB, R, Z


def write_pellet_chord(path, l, rho, s, modB, R, Z, meta):
    """Write the chord table read by ABLATION_NGS (3D mode)."""

    with open(path, 'w') as f:
        f.write("# ABLATION_NGS 3D straight-line chord (VMEC2ASTRA)\n")
        for k, v in meta.items():
            f.write(f"# {k} = {v}\n")
        f.write("# columns: l[m]   rho[-]   s[-]   modB[T]   R[m]   Z[m]\n")
        f.write(f"{len(l)}\n")

        for i in range(len(l)):
            f.write("%14.8e %12.8f %12.8f %12.6f %12.6f %12.6f\n"
                    % (l[i], rho[i], s[i], modB[i], R[i], Z[i]))


def _parse_pellet_group(nl_path):
    """Parse &PELLET_CHORD, allowing comma-separated per-pellet lists for
    LAUNCH_THETA/LAUNCH_PHI/VX/VY/VZ (parse_fortran_namelist only does scalars)."""

    if not os.path.exists(nl_path):
        return {}
    txt = open(nl_path).read()
    # terminate the group on a '/' at the start of a line (namelist convention), 
    # NOT on the '/' inside a path value like 'dat/pellet_chord.dat'
    m = re.search(r'&pellet_chord\b(.*?)^\s*/', txt, re.S | re.M | re.I)
    if not m:
        return {}
    body = m.group(1)

    def grab(name):
        vals = []
        for line in body.splitlines():
            line = line.split('!', 1)[0]
            mm = re.match(r'\s*' + name + r'\s*=\s*(.+)', line, re.I)
            if mm:
                for tok in mm.group(1).split(', '):
                    tok = tok.strip().strip("'\"")
                    if tok:
                        vals.append(tok)
        return vals

    out = {}
    for key in ('launch_theta', 'launch_phi', 'vx', 'vy', 'vz'):
        v = grab(key)
        if v:
            out[key] = [float(x.replace('d', 'e').replace('D', 'e')) for x in v]
    co = grab('chord_out')
    if co:
        out['chord_out'] = co[0]
    npel = grab('npel')
    out['npel'] = int(float(npel[0])) if npel else None
    return out


def maybe_write_pellet_chord(vmc, nl_path='vmec_io/stell_files.nml'):
    """If a &PELLET_CHORD group is present, trace the straight-line chord for
    each configured pellet on the current equilibrium and write the table(s).
    Multi-pellet: LAUNCH_THETA/PHI/VX/VY/VZ may be comma lists of length NPEL;
    pellet ip -> <stem>_<ip>.<ext> (and, for a single pellet, the plain <stem>
    too, so ABLATION_NGS finds it either way).  Failures are logged and
    swallowed so the metric conversion is never broken."""

    try:
        p = _parse_pellet_group(nl_path)
    except Exception:
        return
    if not p or 'vx' not in p:
        return
    try:
        vx = p['vx']
        npel = p.get('npel') or len(vx)
        stem = p.get('chord_out', 'dat/pellet_chord.dat')
        base, ext = os.path.splitext(stem)

        def pick(key, ip, default=0.):
            lst = p.get(key)
            if not lst:
                return default
            return lst[ip] if ip < len(lst) else lst[-1]

        for ip in range(npel):
            theta = pick('launch_theta', ip)
            phi = pick('launch_phi', ip)
            v = [pick('vx', ip), pick('vy', ip), pick('vz', ip)]
            l, rho, s, modB, R, Z = compute_pellet_chord(vmc, theta, phi, v)
            meta = dict(pellet=ip + 1, launch_theta=theta, launch_phi=phi, 
                        v_xyz=f"[{v[0]}, {v[1]}, {v[2]}]", 
                        vmag=np.linalg.norm(v), nfp=vmc.nfp, 
                        penetration_rho_min=np.min(rho))
            paths = [f"{base}_{ip + 1}{ext}"]
            if npel == 1:
                paths.append(stem)
            for out in paths:
                write_pellet_chord(out, l, rho, s, modB, R, Z, meta)
                logger.info(f"Wrote pellet chord {out}: {len(l)} pts, deepest rho={np.min(rho):.4f}, |v|={meta['vmag']:.1f} m/s")

    except Exception as e:
        logger.error(f"pellet chord skipped: {e}")


def write_out_files(wout_file, metric_file, b00_output_file, header_file, radius_out):

    logger.info('Reading VMEC NetCDF output %s', wout_file)

    vmc = VMEC(wout_file)
    logger.info('Calculating moments')
    vmc.calcMoms(nu=64, nv=128)

# NGS pellet: if a launch is configured, trace the straight-line chord for
# this (updated) equilibrium and refresh the table the Fortran model reads.
    maybe_write_pellet_chord(vmc)

    mn_index_b00 = np.where((vmc.xm_nyq == 0) & (vmc.xn_nyq == 0))[0]

    if mn_index_b00.size == 0:
        logger.error("ERROR: Could not find the (m=0, n=0) mode for 'bmnc'.")
        sys.exit(1)

    B00_physical_profile = vmc.bmnc[:, mn_index_b00].squeeze()

# Constants from ASTRA

    NA = NA1 - 1
    GP2 = 2. * np.pi
    HROX = 1./(NA1 - 0.5)

# Define ASTRA arrays
    SXHO = HROX*np.arange(1, NA1+1)
    XRHO = SXHO - 0.5*HROX

#some quantities from VMEC output
    lasym  = vmc.lasym
    RTOR   = vmc.rmajor
    ABC    = vmc.aminor
    volume = vmc.volume
    GVAC = np.abs(vmc.rbtor)  # RB-t edge, coincides with F_b = F0 for tokamak
    if lasym:
        BTOR = GVAC/RTOR      # tokamak
    else:
        BTOR = np.abs(vmc.b0) # stellarator
    phi = np.abs(vmc.phi)     # change from zeta = phi (VMEC) to zeta = -phi (ASTRA)
    phi_a = vmc.phi[-1] / GP2 # Toroidal flux at edge / 2pi

#SHIFT, ELONG, TRIAN, UPDWN
    Fboundary = np.abs(vmc.bvco[-1])  # F at boundary for a tokamak would be RBphi. For stellarator it is more complicated
#calculate RHO and SRHO
    RHOVMECmesh = np.sqrt(phi/(np.pi*BTOR))
    ROC = RHOVMECmesh[-1]
    dphidsb = np.abs(vmc.phipf[-1])
    RHO  = XRHO*ROC
    SRHO = SXHO*ROC
    HRO  = RHO[1] - RHO[0]

# calculate susceptance matrix (on full mesh)
    S11, S12, S21, S22 = vmc.calcSusceptance()
#s31, s32 = vmc.calc_dphidpvb_elements()

# Rewrite into ASTRA coordinate system (minus due to sign convention of jacobian); still VMEC mesh
    fac1 = -phi[-1]/(GP2*BTOR)
    fac2 = fac1/RHOVMECmesh[1:]
    S11[1:] *= fac2
    S12[1:] *= fac2
    S21[1:] *= fac2
    S22 *= fac1

#    Ip_contrib = -phi[-1]*s31[-1]/(GP2*BTOR*RHOVMECmesh)
#    F0_contrib = -phi[-1]*s32[-1]/(GP2*BTOR*RHOVMECmesh)
    Ip_contrib = np.abs(0.*vmc.rbtor)
    F0_contrib = np.abs(0.*vmc.rbtor)

# Force S11[0], S12[0], S21[0] = 0.
    S11[0] = 0.
    S12[0] = 0.
    S21[0] = 0.

# interpolate gp 1 and 2
    RHOVMECmeshsmall = np.delete(RHOVMECmesh, [1, 2])
    S11small = np.delete(S11, [1, 2])
    S12small = np.delete(S12, [1, 2])
    S21small = np.delete(S21, [1, 2])

    S11 = qinterp(RHOVMECmeshsmall, S11small, RHOVMECmesh)
    S12 = qinterp(RHOVMECmeshsmall, S12small, RHOVMECmesh)
    S21 = qinterp(RHOVMECmeshsmall, S21small, RHOVMECmesh)

# calculate V'
    Vp = -vmc.vp * GP2**2 * RHOVMECmesh/fac1   # vp is on full mesh
    Vph = f2h(vmc.vp)              # (dV/ds)/(4*pi*pi) on VMEC half mesh

# calculate <|grad(rho)|>
    vmc.calcGrad_rho()
    GRADROVMEC = vmc.avg_grad_rho * ROC

# calculate <|grad(rho)|^2>
    G1 = vmc.avg_grad_rho2
    G1 = phi[-1]*G1/(np.pi*BTOR)
    g11 = Vp * G1

# calculate <R>, V, aeff
    vmc.calcRmaj()
    vmc.calcVol(Vph)
    Amineff = np.sqrt(2*vmc.Vol/(GP2**2 * vmc.Rmaj))

# calculate ftrapped
    vmc.calcFtrap()

# calculate jpar and iota (not taken to astra automatically)
    jpar   = vmc.jdotb/(1.e6*BTOR)
    muVMEC = vmc.iotaf

# interpolate to ASTRA grid
    SG11   = linterp(RHOVMECmesh, S11, SRHO)
    SG12   = linterp(RHOVMECmesh, S12, SRHO)
    SG21   = linterp(RHOVMECmesh, S21, SRHO)
    SG22   = linterp(RHOVMECmesh, S22, SRHO)
    VR     = linterp(RHOVMECmesh, Vp, RHO)
    VRS    = linterp(RHOVMECmesh, Vp, SRHO)
    GRADRO = linterp(RHOVMECmesh, GRADROVMEC, SRHO)
    G11    = linterp(RHOVMECmesh, g11, SRHO)
    Rmaj   = linterp(RHOVMECmesh, vmc.Rmaj, RHO)
    VOLUM  = linterp(RHOVMECmesh, vmc.Vol, RHO)
    AMETR  = linterp(RHOVMECmesh, Amineff, RHO)
    FTPT   = linterp(RHOVMECmesh, vmc.ftrap, RHO)
    CU     = linterp(RHOVMECmesh, jpar, RHO)
    MU     = linterp(RHOVMECmesh, muVMEC, SRHO)
    SLAT   = VR*GRADRO

# Approximate values at NA1 for shifted grid
    def extrapolate_end(arr):
        return 1.5 * arr[-2] - 0.5 * arr[-3]

    SG11[-1] = extrapolate_end(SG11)
    SG12[-1] = extrapolate_end(SG12)
    SG21[-1] = extrapolate_end(SG21)
    SG22[-1] = extrapolate_end(SG22)
    MV = -SG12/SG11
    VRS[-1]    = extrapolate_end(VRS)
    GRADRO[-1] = extrapolate_end(GRADRO)
    G11[-1]    = extrapolate_end(G11)
    MU[-1]     = extrapolate_end(MU)
    B00_ASTRA_grid = linterp(RHOVMECmesh, B00_physical_profile, RHO)
    xm = vmc.xm.astype(np.int32)   # poloidal mode numbers
    xn = vmc.xn.astype(np.int32)   # toroidal mode numbers

    nn = len(vmc.raxis_cc)

# Writing bin file for ASTRA
    with open(metric_file, 'wb') as f:

        np.array(
            [
                HROX, HRO, ROC, RTOR, ABC, BTOR,
                volume, GVAC, Fboundary,
                vmc.nfp, dphidsb,
                Ip_contrib, F0_contrib
            ],
            dtype=np.float64
        ).tofile(f)

        for arr in [
            RHO, SRHO, SG11, SG12,
            SG21, SG22, MV, VR,
            VRS, GRADRO, G11, Rmaj,
            VOLUM, AMETR, SLAT, FTPT
        ]:
            arr.tofile(f) # They are np.float64, checked

        np.array([vmc.mnmax], dtype=np.int32).tofile(f)

        xm.tofile(f)
        xn.tofile(f)

        vmc.rmnc[-1, :].tofile(f)
        vmc.zmns[-1, :].tofile(f)
        np.array(lasym, dtype=np.int32).tofile(f)
        if lasym:
            vmc.rmns[-1, :].tofile(f)
            vmc.zmnc[-1, :].tofile(f)

        np.array([nn], dtype=np.int32).tofile(f)
        vmc.raxis_cc.tofile(f)
        vmc.zaxis_cs.tofile(f)

        ns_all = vmc.rmnc.shape[0]   # number of flux surfaces
        np.array([ns_all], dtype=np.int32).tofile(f)

        vmc.rmnc.tofile(f)
        vmc.zmns.tofile(f)
        if lasym:
            vmc.rmns.tofile(f)
            vmc.zmnc.tofile(f)
        vmc.phi.tofile(f)

    logger.info(f'Written {metric_file}')

    np.savetxt(b00_output_file, B00_ASTRA_grid)
    logger.info(f'Written B00 profile to {b00_output_file}')

    with open(header_file, 'w') as f:
        f.write(f'{ABC:.10e}\n')
        f.write(f'{phi_a:.10e}\n')

    logger.info(f'Written ABC and phi_a to {header_file}')
    logger.info(f"Major radius (RTOR): {vmc.rmajor:.4f} m")
    logger.info(f"Minor radius (ABC and AB):  {vmc.aminor:.4f} m")
    logger.info(f"Toroidal field on axis (BTOR): {vmc.b0:.4f} T")

    a2 = np.sum(vmc.rmnc[-1, :]*vmc.zmns[-1, :]*xm)
    a_booz = np.sqrt(a2)
    with open(radius_out, 'w') as f:
        f.write(str(a_booz.item()))
    logger.info(f"Written minor radius to: {radius_out}")


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "astra_nrad", type=int, nargs="?", default=91, help="Number of ASTRA radial grid points")

    args = parser.parse_args()

    NA1 = args.astra_nrad

    namelist_path = f'{awd}/vmec_io/stell_files.nml'
    nl_vmec = parse_fortran_namelist(namelist_path, 'VMEC_TO_ASTRA_INPUTS')
    nl_dkes = parse_fortran_namelist(namelist_path, 'ASTRA_DKES_INTERFACE')
    vmec_wd     = f'{awd}/{nl_vmec["vmec_wd"]}'
    wout_file   = f'{vmec_wd}/{nl_vmec["vmec_wout_file"]}'
    metric_file = f'{vmec_wd}/{nl_vmec["vmec2a_metric"]}'
    b00_output_file = f'{vmec_wd}/{nl_dkes["b00_profile_file"]}'
    header_file     = f'{vmec_wd}/{nl_dkes["vmec_header_file"]}'
    radius_out      = f'{vmec_wd}/{nl_dkes["minor_radius_w7as_file"]}'

    write_out_files(wout_file, metric_file, b00_output_file, header_file, radius_out)
