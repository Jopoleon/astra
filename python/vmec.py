import logging
import numpy as np
from scipy.io import netcdf_file
from scipy.interpolate import interp1d

logger = logging.getLogger('vmec2a.VMEC')

GP2 = 2.*np.pi

# Quadratic interpolation
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
# Bottle neck at runtime
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
        self.g_sum = np.sum(self.g, axis=(1, 2))


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


    def calcAll(self):

#some quantities from VMEC output
        self.phi_abs = np.abs(self.phi) 

        self.Fboundary = np.abs(self.bvco[-1])  # F at boundary for a tokamak would be RBphi. For stellarator it is more complicated

        self.GVAC = np.abs(self.rbtor)  # RB-t edge, coincides with F_b = F0 for tokamak
        if self.lasym:
            self.BTOR = self.GVAC/self.rbtor      # tokamak
        else:
            self.BTOR = np.abs(self.b0) # stellarator

        self.RHOVMECmesh = np.sqrt(self.phi_abs/(np.pi*self.BTOR))
        self.dphidsb = np.abs(self.phipf[-1])

        fac1 = -self.phi_abs[-1]/(GP2*self.BTOR)
        fac2 = fac1/self.RHOVMECmesh[1:]

# calculate susceptance matrix (on full mesh)
        S11, S12, S21, S22 = self.calcSusceptance()
        
        S11[1:] *= fac2
        S12[1:] *= fac2
        S21[1:] *= fac2
        S22 *= fac1

# Force S11[0], S12[0], S21[0] = 0.
        S11[0] = 0.
        S12[0] = 0.
        S21[0] = 0.

# interpolate gp 1 and 2
        RHOVMECmeshsmall = np.delete(self.RHOVMECmesh, [1, 2])
        S11small = np.delete(S11, [1, 2])
        S12small = np.delete(S12, [1, 2])
        S21small = np.delete(S21, [1, 2])

        self.S11a = qinterp(RHOVMECmeshsmall, S11small, self.RHOVMECmesh)
        self.S12a = qinterp(RHOVMECmeshsmall, S12small, self.RHOVMECmesh)
        self.S21a = qinterp(RHOVMECmeshsmall, S21small, self.RHOVMECmesh)
        self.S22a = S22

# calculate V'
        self.Vpf = -self.vp * GP2**2 * self.RHOVMECmesh/fac1   # vp is on full mesh
        Vph = f2h(self.vp)              # (dV/ds)/(4*pi*pi) on VMEC half mesh

# calculate <|grad(rho)|>
        self.calcGrad_rho()

# calculate <|grad(rho)|^2>
        G1 = self.avg_grad_rho2
        G1 = self.phi_abs[-1]*G1/(np.pi*self.BTOR)
        self.g11 = self.Vpf * G1

# calculate <R>, V, aeff
        self.calcRmaj()
        self.calcVol(Vph)
        self.Amineff = np.sqrt(2*self.Vol/(GP2**2 * self.Rmaj))

# calculate ftrapped
        self.calcFtrap()

# calculate jpar and iota (not taken to astra automatically)
        self.jpar = self.jdotb/(1.e6*self.BTOR)


    def write_header(self, header_file):
        phi_a = self.phi[-1] / GP2 # Toroidal flux at edge / 2pi, with sign
        with open(header_file, 'w') as f:
            f.write(f'{self.aminor:.10e}\n')
            f.write(f'{phi_a:.10e}\n')
            logger.info(f'Written ABC and phi_a to {header_file}')


    def write_amin(self, radius_out):
        a2 = np.sum(self.rmnc[-1, :]*self.zmns[-1, :]*self.xm)
        a_booz = np.sqrt(a2)
        with open(radius_out, 'w') as f:
            f.write(str(a_booz.item()))
            logger.info(f"Written minor radius to: {radius_out}")
       
