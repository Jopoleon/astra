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

    return dict(R=R, Z=Z, modB=abs(modB), grad_s=grad_s, grad_t=grad_t, grad_z=grad_z,
                e_t=e_t, e_z=e_z)


def lcfs_points(vmc, theta, zeta):
    """Cartesian points on the LCFS for flat arrays theta, zeta of equal length."""

    a = np.outer(theta, vmc.xm) + np.outer(zeta, vmc.xn)
    R = np.cos(a) @ vmc.rmnc[-1]
    Z = np.sin(a) @ vmc.zmns[-1]
    if vmc.lasym:
        R = R + np.sin(a) @ vmc.rmns[-1]
        Z = Z + np.cos(a) @ vmc.zmnc[-1]
    return np.stack((R*np.cos(zeta), R*np.sin(zeta), Z), axis=1)


def lcfs_entry(vmc, p0, vdir, tol=1e-10, nseed=12):
    """Where a straight injector line first crosses the LCFS.

Returns (theta, zeta, t) for the crossing nearest the injector in front of it,
t being the distance from p0 along the unit direction, or None when the line
misses the plasma. p0 and vdir are Cartesian [m] in the VMEC frame,
x = R cos(zeta), y = R sin(zeta), z = Z.

Solves X_LCFS(theta, zeta) - p0 - t*vhat = 0 by Newton from the best points of
a coarse surface scan; the Jacobian columns are the covariant e_theta, e_zeta
and -vhat."""

    vhat = np.asarray(vdir, dtype=float)
    vmag = np.linalg.norm(vhat)
    if vmag <= 0.:
        return None
    vhat = vhat/vmag
    p0 = np.asarray(p0, dtype=float)

    nth, nze = 64, 32*max(int(vmc.nfp), 1) + 32
    th, ze = np.meshgrid(np.linspace(0., 2.*np.pi, nth, endpoint=False),
                         np.linspace(0., 2.*np.pi, nze, endpoint=False), indexing='ij')
    th, ze = th.ravel(), ze.ravel()

    rel = lcfs_points(vmc, th, ze) - p0
    tt = rel @ vhat
    dperp = np.linalg.norm(rel - np.outer(tt, vhat), axis=1)

    ahead = np.flatnonzero(tt > 0.)
    if ahead.size == 0:
        return None
    seeds = ahead[np.argsort(dperp[ahead])[:nseed]]

    hits = []
    for k in seeds:
        y = np.array([th[k], ze[k], tt[k]])
        for _ in range(40):
            g = PointGeom(vmc, 1., y[0], y[1])
            res = np.array([g['R']*np.cos(y[1]), g['R']*np.sin(y[1]), g['Z']]) \
                  - p0 - y[2]*vhat
            if np.linalg.norm(res) < tol:
                break
            jac = np.column_stack((g['e_t'], g['e_z'], -vhat))
            try:
                step = np.linalg.solve(jac, -res)
            except np.linalg.LinAlgError:
                break
            # Damped, so a seed on the far side of the surface cannot jump the
            # solver into a different crossing.
            big = np.max(np.abs(step[:2]))
            if big > 0.3:
                step = step*(0.3/big)
            y = y + step
        else:
            continue
        if np.all(np.isfinite(y)) and y[2] > 0. and np.linalg.norm(res) < 1e-8:
            hits.append(y)

    if not hits:
        return None
    theta, zeta, t = min(hits, key=lambda h: h[2])
    return theta % (2.*np.pi), zeta % (2.*np.pi), t


class PELLET():


    def __init__(self):
        self.pars_d = {}


    def parse_pellet_nml(self, f_setting):
        """Read the pellet launch from &pellet_ngs in the exp namelist.

Two ways to aim a pellet. Either pel_theta and pel_phi place it on the LCFS and
pel_ux/pel_uy/pel_uz give the direction, of any length; or pel_x0/y0/z0 and
pel_x1/y1/z1 give two Cartesian points [m] on the injector axis, the start and
the tip, and the line through them is the path. The point pair sets both the
direction and the entry point, so it overrides the other four keys. Velocity is
pel_vp times the normalised direction."""

        logger.info('Parse &pellet_ngs group in %s', f_setting)
        self.pars_d = {}
        try:
            nml_d = parse_fortran_namelist(f_setting, "pellet_ngs")
            npel = int(nml_d.get('npel', 0))
        except FileNotFoundError:
            logger.info('No chord: %s does not exist', f_setting)
            return
        except Exception as err:
            # This parser reads a subset of namelist syntax, and vmec2bin writes
            # the metric after this call: a value it cannot read costs the chord,
            # not the run.
            logger.warning('No chord: cannot read &pellet_ngs in %s: %s', f_setting, err)
            return

        # Without a chord file ABLATION_NGS runs 1D and never reads one.
        stem = nml_d.get('pel_chord_file', '')
        if npel <= 0 or not stem:
            logger.info('No chord: npel = %d, pel_chord_file = %r', npel, stem)
            return

        # A key shorter than npel leaves the rest at zero, as Fortran does, so
        # pad rather than broadcast; chords() skips what stays unset.
        pars_d = {}
        for key, out in (('pel_theta', 'launch_theta'), ('pel_phi', 'launch_phi'),
                         ('pel_vp', 'vp'), ('pel_ux', 'ux'), ('pel_uy', 'uy'),
                         ('pel_uz', 'uz'),
                         ('pel_x0', 'x0'), ('pel_y0', 'y0'), ('pel_z0', 'z0'),
                         ('pel_x1', 'x1'), ('pel_y1', 'y1'), ('pel_z1', 'z1')):
            val = np.atleast_1d(nml_d.get(key, 0.)).astype(float)
            pars_d[out] = np.pad(val, (0, max(0, npel - val.size)))[:npel]

        p0 = np.stack([pars_d.pop(k) for k in ('x0', 'y0', 'z0')])
        seg = np.stack([pars_d.pop(k) for k in ('x1', 'y1', 'z1')]) - p0
        seglen = np.linalg.norm(seg, axis=0)
        has_p0 = seglen > 0.
        with np.errstate(invalid='ignore', divide='ignore'):
            for j, key in enumerate(('ux', 'uy', 'uz')):
                pars_d[key] = np.where(has_p0, seg[j]/np.where(has_p0, seglen, 1.),
                                       pars_d[key])

        unorm = np.sqrt(pars_d['ux']**2 + pars_d['uy']**2 + pars_d['uz']**2)
        if not np.any(unorm > 0.):
            logger.info('No chord: neither pel_ux/uy/uz nor a pel_x0/x1 point pair'
                        ' in &pellet_ngs')
            return
        with np.errstate(invalid='ignore', divide='ignore'):
            for key in ('ux', 'uy', 'uz'):
                pars_d[key] = np.where(unorm > 0., pars_d[key]/unorm, 0.)

        for vkey, ukey in (('vx', 'ux'), ('vy', 'uy'), ('vz', 'uz')):
            pars_d[vkey] = pars_d['vp'] * pars_d[ukey]
        pars_d['p0'] = p0
        pars_d['has_p0'] = has_p0
        pars_d['npel'] = npel
        pars_d['chord_out'] = stem
        self.pars_d = pars_d


    def chords(self, vmc):

        if not self.pars_d:
            return

        try:
            npel = self.pars_d.get('npel') or len(self.pars_d['vx'])
            stem = self.pars_d.get('chord_out', 'dat/pellet_chord.dat')
            base, ext = os.path.splitext(stem)
            vel = np.stack((self.pars_d['vx'], self.pars_d['vy'], self.pars_d['vz']))
            for jpel in range(npel):
                if np.linalg.norm(vel[:, jpel]) <= 0.:
                    logger.info('Pellet %d has no speed or no direction, no chord', jpel + 1)
                    continue
                if self.pars_d['has_p0'][jpel]:
                    hit = lcfs_entry(vmc, self.pars_d['p0'][:, jpel], vel[:, jpel])
                    if hit is None:
                        logger.info('Pellet %d: the injector line misses the plasma,'
                                    ' no chord', jpel + 1)
                        continue
                    self.theta, self.phi, _ = hit
                else:
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
