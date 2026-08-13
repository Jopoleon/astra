import os, logging, traceback
import numpy as np
from parse_fortran_nml import parse_fortran_namelist

logger = logging.getLogger('vmec2a.pellet')


def bracket(s, ns):
    s = min(max(s, 0.), 1.)
    x = s * (ns - 1)
    j = int(np.floor(x))
    if j >= ns - 1:
        j = ns - 2
    return j, x - j


def interp_row(arr, s, ns):
    """Linear interp of arr[ns, mn] at scalar s on the uniform s grid."""
    j, f = bracket(s, ns)
    return arr[j] * (1. - f) + arr[j + 1] * f


def deriv_row(arr, s, ns):
    """d/ds of the linearly-interpolated arr[ns, mn] (exact per-cell slope,
       self-consistent with interp_row so grad_i.e_j = delta_ij holds)."""
    j, _ = bracket(s, ns)
    return (arr[j + 1] - arr[j]) * (ns - 1)

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

def PointGeom(vmc, s, theta, zeta):
    """Return dict with R, Z, modB and grad_s/grad_t/grad_z (each 3-vectors)."""

    ns = vmc.rmnc.shape[0]
    rc  = interp_row(vmc.rmnc, s, ns)
    zc  = interp_row(vmc.zmns, s, ns)
    rc_s = deriv_row(vmc.rmnc, s, ns)
    zc_s = deriv_row(vmc.zmns, s, ns)
    a = vmc.xm*theta + vmc.xn*zeta
    ca, sa = np.cos(a), np.sin(a)

    R = np.dot(rc, ca)
    Z = np.dot(zc, sa)
    R_t = np.dot(rc, -vmc.xm * sa)
    R_z = np.dot(rc, -vmc.xn * sa)
    Z_t = np.dot(zc,  vmc.xm * ca)
    Z_z = np.dot(zc,  vmc.xn * ca)
    R_s = np.dot(rc_s, ca)
    Z_s = np.dot(zc_s, sa)

    bc = interp_row(vmc.bmnc, s, ns)
    an = vmc.xm_nyq * theta + vmc.xn_nyq * zeta
    modB = np.dot(bc, np.cos(an))

    if vmc.lasym:
        rsc  = interp_row(vmc.rmns, s, ns)
        zcc  = interp_row(vmc.zmnc, s, ns)
        rsc_s = deriv_row(vmc.rmns, s, ns)
        zcc_s = deriv_row(vmc.zmnc, s, ns)
        R += np.dot(rsc, sa)
        Z += np.dot(zcc, ca)
        R_t += np.dot(rsc,  vmc.xm*ca)
        R_z += np.dot(rsc,  vmc.xn*ca)
        Z_t += np.dot(zcc, -vmc.xm*sa)
        Z_z += np.dot(zcc, -vmc.xn*sa)
        R_s += np.dot(rsc_s, sa)
        Z_s += np.dot(zcc_s, ca)
        bsc = interp_row(vmc.bmns, s, ns)
        modB += np.dot(bsc, np.sin(an))

    cz, sz = np.cos(zeta), np.sin(zeta)
    # covariant basis vectors in Cartesian (x = R cos z, R sin z, Z)
    e_s = np.array([R_s*cz,        R_s*sz,         Z_s])
    e_t = np.array([R_t*cz,        R_t*sz,         Z_t])
    e_z = np.array([R_z*cz - R*sz, R_z *sz + R*cz, Z_z])
    sqrtg = np.dot(e_s, np.cross(e_t, e_z))
    if abs(sqrtg) < 1e-30:
        sqrtg = 1e-30
    grad_s = np.cross(e_t, e_z) / sqrtg
    grad_t = np.cross(e_z, e_s) / sqrtg
    grad_z = np.cross(e_s, e_t) / sqrtg

    return dict(R=R, Z=Z, modB=abs(modB), grad_s=grad_s, grad_t=grad_t, grad_z=grad_z)


class PELLET():


    def __init__(self):
        pass


    def parse_pellet_nml(self, nl_path='vmec_io/stell_files.nml'):
        """Parse the &PELLET_CHORD namelist group."""

        logger.info(f'Parse PELLET_CHORD namelit group in {nl_path}')
        try:
            self.pars_d = parse_fortran_namelist(nl_path, "PELLET_CHORD")
            for key in ('vx', 'vy', 'vz', 'launch_theta', 'launch_phi'):
                self.pars_d[key] = np.atleast_1d(self.pars_d[key])
        except FileNotFoundError:
            return {}


    def chords(self, vmc):

        try:
            npel = self.pars_d.get('npel') or len(self.pars_d['vx'])
            stem = self.pars_d.get('chord_out', 'dat/pellet_chord.dat')
            base, ext = os.path.splitext(stem)
            vel = np.stack((self.pars_d['vx'], self.pars_d['vy'], self.pars_d['vz']))
            for jpel in range(npel):
                self.theta = self.pars_d['launch_theta'][jpel]
                self.phi   = self.pars_d['launch_phi'  ][jpel]
                self.v_xyz = vel[:, jpel]
                self.compute_chord(vmc)
                self.meta_d = dict(pellet=jpel + 1, launch_theta=self.theta, launch_phi=self.phi,
                        v_xyz=self.v_xyz, vmag=np.linalg.norm(self.v_xyz), nfp=vmc.nfp,
                        penetration_rho_min=np.min(self.rho))
                f_out = f"{base}_{jpel + 1}{ext}"
                self.write_chord_file(f_out)
                logger.info(f"Wrote pellet chord {f_out}: {len(self.rho)} pts, deepest rho={np.min(self.rho):.4f}, |v|={self.meta_d['vmag']:.1f} m/s")

        except Exception as e:
            traceback.print_exc()
            logger.error(f"pellet chord skipped: {e}")


    def compute_chord(self, vmc, max_len=None, nstep=2000):
        """Trace a straight real-space line from the LCFS launch point along v_xyz
and tabulate the geometry the Fortran ablation ODE needs.

The pellet flies with fixed Cartesian velocity v_xyz [m/s]; along the line
the flux coordinates evolve by  d(s, theta, zeta)/dl = vhat . grad(s, theta, zeta),
integrated with RK4 in real arclength l (vhat = v/|v|).  Integration starts
at s=1, (theta, zeta)=(launch_theta, launch_phi) and stops when the pellet
exits (s back >= 1 after entering) or after max_len.

Returns arrays (l, rho, s, modB, R, Z), l ascending from the launch point.
rho = sqrt(s) is the ASTRA flux label."""

        vmag = np.linalg.norm(self.v_xyz)
        if vmag <= 0.:
            raise ValueError("pellet velocity magnitude is zero")
        vhat = self.v_xyz/vmag

# A generous default path length: a few times the device size

        if max_len is None:
            max_len = 4.*vmc.aminor
        dl = max_len/nstep

        def vhat_dot(vhat, g):
            return np.array([np.dot(vhat, g['grad_s']), np.dot(vhat, g['grad_t']), np.dot(vhat, g['grad_z'])])

        trajectory = []
        self.modB  = []
        self.R     = []
        self.Z     = []
        entered = False
        y = np.array([1., self.theta, self.phi])

        for i in range(nstep):

            g = PointGeom(vmc, *y)
            trajectory.append(y.copy())
            self.modB.append(g['modB'])
            self.R.append(g['R'])
            self.Z.append(g['Z'])
            k1 = vhat_dot(vhat, g)

            y2 = y + 0.5*dl*k1
            g2 = PointGeom(vmc, *y2)
            k2 = vhat_dot(vhat, g2)

            y3 = y + 0.5*dl*k2
            g3 = PointGeom(vmc, *y3)
            k3 = vhat_dot(vhat, g3)
            
            y4 = y + dl*k3
            g4 = PointGeom(vmc, *y4)
            k4 = vhat_dot(vhat, g4)

            y += (dl / 6.) * (k1 + 2*k2 + 2*k3 + k4)
            if y[0] < 0.:  # protect from tiny negative overshoot at the axis
                y[0] = -y[0]
            if y[0] < 0.98:
                entered = True
            if entered and y[0] >= 1.: # trajectory exiting the plasma after being inside
                break

# Add last point
        g = PointGeom(vmc, *y)
        trajectory.append(y.copy())
        self.modB.append(g['modB'])
        self.R.append(g['R'])
        self.Z.append(g['Z'])

        trajectory = np.array(trajectory)
        self.ss  = np.minimum(trajectory[:, 0], 1.)
        self.ths = trajectory[:, 1]
        self.zes = trajectory[:, 2]
        self.dl = dl
        self.rho = np.sqrt(self.ss)


    def write_chord_file(self, f_out):
        """Write the chord table read by ABLATION_NGS (3D mode)."""

        nl = len(self.rho)

        with open(f_out, 'w') as f:
            f.write("# ABLATION_NGS 3D straight-line chord (VMEC2ASTRA)\n")
            for k, v in self.meta_d.items():
                f.write(f"# {k} = {v}\n")
            f.write("# columns: l[m]   rho[-]   s[-]   modB[T]   R[m]   Z[m]\n")
            f.write(f"{nl}\n")

            for i in range(nl):
                f.write("%14.8e %12.8f %12.8f %12.6f %12.6f %12.6f\n"
                    % (i*self.dl, self.rho[i], self.ss[i], self.modB[i], self.R[i], self.Z[i]))
