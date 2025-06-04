#!/usr/bin/env python3

import os, logging, argparse, json
import numpy as np
from scipy.special import ellipk, ellipe
from scipy.constants import mu_0


fmt = logging.Formatter('%(asctime)s | %(name)s | %(levelname)s: %(message)s', '%H:%M:%S')
logger = logging.getLogger('green')
if len(logger.handlers) == 0:
    hnd = logging.StreamHandler()
    hnd.setFormatter(fmt)
    logger.addHandler(hnd)

logger.setLevel(logging.DEBUG)
#logger.setLevel(logging.INFO)

loc_dir = os.path.dirname(os.path.realpath(__file__))
awd = os.path.dirname(loc_dir)
grIOdir = '%s/exp/cnf' %awd

gr_flt = np.float64
gr_int = np.int32

def format_numpy_json(data, indent=0, indent_step=2, fmt="%.6e"):
    pad = ' ' * indent
    next_pad = ' ' * (indent + indent_step)
    if isinstance(data, dict):
        items = []
        for i, (k, v) in enumerate(data.items()):
            formatted = format_numpy_json(v, indent + indent_step, indent_step, fmt)
            items.append(f'{next_pad}"{k}": {formatted}')
        return "{\n" + ",\n".join(items) + f"\n{pad}}}"
    elif isinstance(data, (list, tuple, np.ndarray)):
        arr = np.array(data)
        if arr.ndim == 1:
            items = [fmt % x if isinstance(x, float) else str(x) for x in arr]
            return "[" + ", ".join(items) + "]"
        else:
            blocks = [format_numpy_json(sub, indent + indent_step, indent_step, fmt)
                      for sub in arr]
            return "[\n" + ",\n".join(next_pad + block for block in blocks) + f"\n{pad}]"
    elif isinstance(data, float):
        return fmt % data
    else:
        return str(data)

def greenFunction(Rloc, Zloc, Rcoil, Zcoil):
    sum_sq = (Rloc + Rcoil)**2 + (Zloc - Zcoil)**2
    k_sq = np.clip(4.*Rloc*Rcoil/sum_sq, 0.0, 1.0)
    pre_fac = np.sqrt(sum_sq)*0.25             # Consistent with ASTRA matrices
    K = ellipk(k_sq)
    E = ellipe(k_sq)
    return pre_fac*((2. - k_sq)*K - 2.*E)

def GreensBr_2piR(Rloc, Zloc, Rcoil, Zcoil, eps=1e-10):
    z_zo_sq = (Zcoil - Zloc)**2
    Rsq = Rcoil**2 + Rloc**2 
    Rprod = 2*Rcoil*Rloc
    sum_sq = Rsq + Rprod + z_zo_sq
    rsum4 = np.minimum(2.*Rprod/sum_sq, 1)
# 1/R is omitted because there's "* PFRzw(ki, 1)" in the fortran call
# 2.*pi is omitted because there's a "TWOPI" in the fortran call
    DS = np.maximum(Rsq - Rprod + z_zo_sq, eps)
    return mu_0/np.sqrt(sum_sq)*(Zloc - Zcoil) *(-ellipk(rsum4) + (Rsq + z_zo_sq)/DS * ellipe(rsum4))

def greenBoundary(Rgrid, Zgrid):
    nr2 = len(Rgrid) - 2
    nz2 = len(Zgrid) - 2
    dr = (Rgrid[-2] - Rgrid[0])/nr2
    dz = (Zgrid[-2] - Zgrid[0])/nz2
    llp = 2*(nr2 + nz2)
    der = np.concatenate((
        np.full(nr2, dr, dtype=gr_flt), np.full(nz2, dz, dtype=gr_flt),
        np.full(nr2, dr, dtype=gr_flt), np.full(nz2, dz, dtype=gr_flt)
    ))

    r_up = r_down = Rgrid[1:-1]
    r_left  = np.full(nz2, Rgrid[ 0], dtype=gr_flt)
    r_right = np.full(nz2, Rgrid[-1], dtype=gr_flt)
    z_down  = np.full(nr2, Zgrid[ 0], dtype=gr_flt)
    z_up    = np.full(nr2, Zgrid[-1], dtype=gr_flt)
    z_left = z_right = Zgrid[1:-1]
    rbnd = np.concatenate((r_down, r_right, r_up, r_left))
    zbnd = np.concatenate((z_down, z_right, z_up, z_left))

    gfm = np.zeros((llp, llp), dtype=gr_flt)
    for j in range(llp-1):
        gfm[j, j+1:] = greenFunction(rbnd[j], zbnd[j], rbnd[j+1:], zbnd[j+1:])
    gfm += np.triu(gfm, 1).T # symmetrise matrix

    tt = der/(rbnd*4)
    diag_values = -tt*(np.log(tt) - 2.*np.log(2.) + 1.) * 2 * rbnd**2/der
    np.fill_diagonal(gfm, diag_values)
    return gfm

def identity(R, dr, dz, ntype):
# below replicates what is done in SPIDER (A. A. Ivanov and S. Yu. Medvedev, KIAM Preprint Nr 39, 2009 Moscow)
    (ind1, ) = np.where(ntype == 1)
    (ind2, ) = np.where(ntype == 2)
    gfi = np.zeros_like(R)
    if len(ind2) > 0:
        gfi[ind2] = R[ind2]*(np.log(8.*R[ind2]/(0.2236*(dr[ind2] + dz[ind2]))) - 2.)/2. # From A. Kavin,  used in SPIDER (A. A. Ivanov and S. Yu. Medvedev, KIAM Preprint Nr 39, 2009 Moscow),  self inductance of a rectangular coil in toroidal direction
    if len(ind1) > 0:
        dl  = np.hypot(dr, dz)
        for i in (-1, 0, 1):
            for j in (-1, 0, 1):
                if i == j:
                    gfi[ind1] += (R[ind1] + dr[ind1]*i/3.) * (np.log( 8.*(R[ind1] + dr[ind1]*i/3.)/(dl[ind1]/3.) ) - 0.5)/2.
                else:
                    gfi[ind1] += greenFunction(R[ind1] + dr[ind1]*i/3., dz[ind1]*i/3., R[ind1] + dr[ind1]*j/3., dz[ind1]*j/3.)
                
        gfi[ind1] /= 9.
    return gfi

def nonIdentity(R1loc, Z1loc, R2, Z2, dr1, dz1, dr2, dz2, ntype1, ntype2):
    (ind1, ) = np.where(ntype2 == 1)
    (ind2, ) = np.where(ntype2 == 2)
    ntype2range = {1: (-1, 0, 1), 2: (0,)}
    irange = ntype2range[ntype1]
    gf = np.zeros_like(R2)
    for i in irange:
        if len(ind1) > 0:
            for j in (-1, 0, 1):
                gf[ind1] += greenFunction(R1loc + dr1*i/3., Z1loc + dz1*i/3., R2[ind1] + dr2[ind1]*j/3., Z2[ind1] + dz2[ind1]*j/3.)
        if len(ind2) > 0:
            gf[ind2] += greenFunction(R1loc + dr1*i/3., Z1loc + dz1*i/3., R2[ind2], Z2[ind2])
    if len(ind1) > 0:
        gf[ind1] /= 3.
    return gf/float(len(irange))


class GREEN_MATRICES:


    def __init__(self):
        pass


    def fromMachineInput(self, f_json):

        logger.info('Reading %s', f_json)

# Read machine input file

        with open(f_json) as fjson:
            in_d = json.load(fjson)
        for key, val in in_d.items():
            if isinstance(val, list):
                if isinstance(val[0], int):
                    setattr(self, key, np.array(val, dtype=gr_int))
                else:
                    setattr(self, key, np.array(val, dtype=gr_flt))
            else:
                setattr(self, key, val)

# Spatial grids
        if not hasattr(self, 'Rmin'):
            logger.warning('"Rmin" not found, skipping Green-function calculation for %s', f_json)
            return
        self.Rgrid = np.linspace(self.Rmin, self.Rmax, self.nR, endpoint=True, dtype=gr_flt)
        self.Zgrid = np.linspace(self.Zmin, self.Zmax, self.nZ, endpoint=True, dtype=gr_flt)
        dr = (self.Rmax - self.Rmin)/float(self.nR - 1)
        dz = (self.Zmax - self.Zmin)/float(self.nZ - 1)

# Coils geometry
        self.angh_coil = np.radians(self.angh_coil)
        self.ang_coil  = np.radians(self.ang_coil)
        nCoils = len(self.R_coil)
        self.nActive = np.max(self.m_equiv)
        numeqcump   = np.zeros(self.nActive, dtype=gr_int)
        self.R_cond = np.zeros(self.nActive, dtype=gr_flt)
        self.Z_cond = np.zeros_like(self.R_cond)
        for jequiv1 in range(1, self.nActive+1):
            mask = self.m_equiv == jequiv1
            if np.any(mask):
                self.R_cond[jequiv1 - 1] = np.mean(self.R_coil[mask])
                self.Z_cond[jequiv1 - 1] = np.mean(self.Z_coil[mask])

# Coil resistivity

        resConduc_diag = []

# Limiter geometry

        indR = ((self.Rlim - self.Rmin)/dr + 0.5).astype(int)
        indZ = ((self.Zlim - self.Zmin)/dz + 0.5).astype(int)
        self.Rlim = self.Rgrid[indR]
        self.Zlim = self.Zgrid[indZ]

        ilim_maxR = int((self.lim_maxR - self.Rmin)/dr + 0.5)
        ilim_minR = int((self.lim_minR - self.Rmin)/dr + 0.5)
        ilim_maxZ = int((self.lim_maxZ - self.Zmin)/dz + 0.5)
        ilim_minZ = int((self.lim_minZ - self.Zmin)/dz + 0.5)
        self.limRZ = (self.Rgrid[ilim_maxR], self.Rgrid[ilim_minR], self.Zgrid[ilim_maxZ], self.Zgrid[ilim_minZ])
        self.zLimPotential = np.zeros((self.nR, self.nZ), dtype=gr_int)
        self.zLimPotential[ilim_minR: ilim_maxR, ilim_minZ: ilim_maxZ] = 1
        self.Rg, self.Zg = np.meshgrid(self.Rgrid, self.Zgrid)

# Blanket
        self.nConduc = self.nActive
        if hasattr(self, 'x1'):
            x7 = self.x3 - self.x1
            x8 = self.x4 - self.x2
            x9 = np.hypot(x7, x8)
            self.R_blan = np.ravel([self.x1 + 0.25*x7, self.x1 + 0.75*x7], order='F')
            self.Z_blan = np.ravel([self.x2 + 0.25*x8, self.x2 + 0.75*x8], order='F')
            self.dhoriz = np.ravel([0.5*x7, 0.5*x7], order='F')
            self.dvert  = np.ravel([0.5*x8, 0.5*x8], order='F')
            area_blan   = np.ravel([self.width_blan*x9, self.width_blan*x9], order='F')
            ssfw = np.sum(area_blan/self.R_blan)
            n_blanket = len(self.R_blan)

            self.R_cond = np.append(self.R_cond[:self.nActive], self.R_blan)
            self.Z_cond = np.append(self.Z_cond[:self.nActive], self.Z_blan)
            if self.res_blan > 0:
                resConduc_diag.extend(self.res_blan*self.R_blan/area_blan*ssfw)
            self.nConduc += n_blanket

# Blanketpc

        if hasattr(self, 'R_blan_pc'):
            self.R_cond = np.append(self.R_cond, self.R_blan_pc)
            self.Z_cond = np.append(self.Z_cond, self.Z_blan_pc)
            if self.res_blan < 0:
                resConduc_diag.extend(self.res_blan_pc)
            self.nConduc += len(self.R_blan_pc)

        self.resConduc_diag = np.array(resConduc_diag, dtype=gr_flt)

# Ferromagnet stuff

        if hasattr(self, 'r0ferro'):
            self.angleferro = np.radians(self.angleferro)
        
        logger.debug('nactive, ncoils, nconduc %d %d %d', self.nActive, nCoils, self.nConduc)


    def setCoilProperties(self):

#------
# Coils

        logger.info('Setting coil properties')
        self.nBlocks = 0
        self.nConduc = 0
        nCoils = len(self.R_coil)

        identcoil = np.ones(nCoils, dtype=bool)
        sin_coil  = np.sin(self.ang_coil)
        cos_coil  = np.cos(self.ang_coil)
        sin_coilh = np.sin(self.angh_coil)
        cos_coilh = np.cos(self.angh_coil)

        x12 = self.dR_coil/self.dZ_coil
        indi = np.sqrt(x12   *self.n_elem_coil).astype(gr_int) + 1
        indj = np.sqrt(1./x12*self.n_elem_coil).astype(gr_int) + 1
        dr1 = self.dR_coil/indi
        dz1 = self.dZ_coil/indj
        dr_cs = dr1*cos_coilh
        dz_cs = dz1*cos_coil
        dr_ss = dr1*sin_coilh
        dz_ss = dz1*sin_coil
        r1 = self.R_coil - 0.5*( self.dR_coil*cos_coilh + self.dZ_coil*cos_coil) + 0.5*(dr_cs + dz_cs)
        z1 = self.Z_coil - 0.5*( self.dR_coil*sin_coilh + self.dZ_coil*sin_coil) + 0.5*(dr_ss + dz_ss)

        Rce   = []
        Zce   = []
        dRce  = []
        dZce  = []
        tatmp = []
        equivforce = []
        equivtmp   = []

        for icoil in range(nCoils):
            if identcoil[icoil]:
                self.nConduc += 1
                (ind, ) = np.where(self.m_equiv[icoil+1: nCoils] == self.m_equiv[icoil])
                ind += icoil + 1   # correct array index
                indk = np.append(icoil, ind)
                identcoil[ind] = 0
                for j in indk:
                    indij = indi[j]*indj[j]
                    dRce.extend(indij*[dr_cs[j]]) #list(dr_cs[j] + flt0_ij)
                    dZce.extend(indij*[dz_ss[j]])
                    tatmp.extend(indij*[self.m_turns[j]/float(indij)])
                    equivforce.extend(indij*[j])
                    equivtmp.extend(indij*[self.nConduc-1])
                    r_array = r1[j] + dr_cs[j] * np.arange(indi[j])
                    z_array = z1[j] + dr_ss[j] * np.arange(indi[j])
                    for jj in range(indj[j]):
                        Rce.extend(r_array + jj*dz_cs[j])
                        Zce.extend(z_array + jj*dz_ss[j])
        self.nBlocks += nCoils
        nctype = len(Rce)*[2]

#--------
# Blanket

        if hasattr(self, 'R_blan_pc'):
            n_blanket_pc = len(self.R_blan_pc)
            n_elem_blanket_pc = 9 # sub element of a blanket element, hardwired to 9 for now
            ind_pc = np.round(np.sqrt(n_elem_blanket_pc + 1.e-6)).astype(gr_int)
            dx = np.sqrt(self.area_blan_pc)/ind_pc
            ind_pc2 = ind_pc**2
            dx_area = 1./ind_pc2
            nctype.extend(ind_pc2*n_blanket_pc*[2])
            for j in range(n_blanket_pc):
                self.nConduc += 1
                self.nBlocks += 1
                equivtmp.extend(ind_pc2*[self.nConduc-1])
                equivforce.extend(ind_pc2*[self.nBlocks-1])
                dRce.extend(ind_pc2*[dx[j]])
                dZce.extend(ind_pc2*[dx[j]])
                tatmp.extend(ind_pc2*[dx_area[j]])
                for jjj in range(ind_pc):
                    Zce.extend(ind_pc*[self.Z_blan_pc[j] + jjj*dx[j]])
                    Rce.extend(self.R_blan_pc[j] + np.arange(ind_pc)*dx[j])

        if hasattr(self, 'R_blan'):
            n_blanket = len(self.R_blan)
            Rce.extend(self.R_blan)
            Zce.extend(self.Z_blan)
            dRce.extend(self.dhoriz)
            dZce.extend(self.dvert)
            nctype.extend(n_blanket*[1])
            tatmp.extend(n_blanket*[1.])
            equivtmp.extend(  self.nConduc + np.arange(n_blanket))
            equivforce.extend(self.nBlocks + np.arange(n_blanket))
            self.nConduc += n_blanket
            self.nBlocks += n_blanket

        self.Rce   = np.array(  Rce, dtype=gr_flt)
        self.Zce   = np.array(  Zce, dtype=gr_flt)
        self.dRce  = np.array( dRce, dtype=gr_flt)
        self.dZce  = np.array( dZce, dtype=gr_flt)
        self.tatmp = np.array(tatmp, dtype=gr_flt)
        self.equivforce = np.array(equivforce, dtype=gr_int)
        self.equivtmp   = np.array(  equivtmp, dtype=gr_int)
        self.nctype     = np.array(    nctype, dtype=gr_int)


    def calcGreenMatrices(self):

        self.selfInductances()
        self.gridInductances()
        self.interBlockForces()
        self.greenBnd = greenBoundary(self.Rgrid, self.Zgrid)
        self.ferromagnets()


    def selfInductances(self):

        logger.debug('self induct')
        gf_diag = identity(self.Rce, self.dRce, self.dZce, self.nctype)
        ielem = len(self.tatmp)

        self.indConduc = np.zeros((self.nConduc, self.nConduc))
        for i in range(ielem):
            iii = self.equivtmp[i]
            self.indConduc[iii, iii] += gf_diag[i]*self.tatmp[i]**2
            gf_ta = self.tatmp*nonIdentity(self.Rce[i], self.Zce[i], self.Rce, self.Zce, self.dRce[i], self.dZce[i], self.dRce, self.dZce, self.nctype[i], self.nctype)
            for j in range(ielem):
                if j != i:
                    jjj = self.equivtmp[j]
                    self.indConduc[iii, jjj] += gf_ta[j]*self.tatmp[i]
        self.indConduc *= 0.8*np.pi


    def gridInductances(self):

        logger.info('Grid inductances')
        nR = len(self.Rgrid)
        nZ = len(self.Zgrid)
        self.greeni = np.zeros((nR, nZ, self.nConduc), dtype=gr_flt)
        offsets = np.array([-1, 0, 1], dtype=gr_flt) / 3.
        for i, nctyp in enumerate(self.nctype):
            iii = self.equivtmp[i]
            R0, Z0 = self.Rce[i], self.Zce[i]
            dR, dZ = self.dRce[i], self.dZce[i]
            scale  = self.tatmp[i]
            if nctyp == 1:
                Rlocs = R0 + dR * offsets
                Zlocs = Z0 + dZ * offsets
                gf_off = sum(greenFunction(Rloc, Zloc, self.Rg.T, self.Zg.T) for Rloc, Zloc in zip(Rlocs, Zlocs))
                self.greeni[:, :, iii] += (gf_off / 3.) * scale
            elif nctyp == 2:
                gf_off = greenFunction(self.Rce[i], self.Zce[i], self.Rg.T, self.Zg.T)
                self.greeni[:, :, iii] += gf_off * scale
        self.greeni *= 0.4


    def interBlockForces(self):

        logger.debug('forces %d', self.nBlocks)
        nR = len(self.Rgrid)
        nZ = len(self.Zgrid)
        self.dGreeniRpl = np.zeros((nR, nZ, self.nBlocks), dtype=gr_flt)
        self.dGreeniZpl = np.zeros_like(self.dGreeniRpl)
        self.dGreeniRj  = np.zeros((self.nBlocks, self.nBlocks), dtype=gr_flt)
        self.dGreeniZj  = np.zeros_like(self.dGreeniRj)
        dr = (self.Rgrid[-1] - self.Rgrid[0])/float(nR - 1)
        dz = (self.Zgrid[-1] - self.Zgrid[0])/float(nZ - 1)
        ielem = len(self.tatmp)

        for i in range(ielem):
            iii = self.equivforce[i]
            prefac_r = -0.8*np.pi*self.tatmp[i]/dr
            prefac_z = -0.8*np.pi*self.tatmp[i]/dz
            g_r1 = greenFunction(self.Rce[i] + dr/2., self.Zce[i], self.Rce, self.Zce)
            g_r2 = greenFunction(self.Rce[i] - dr/2., self.Zce[i], self.Rce, self.Zce)
            g_z1 = greenFunction(self.Rce[i], self.Zce[i] + dz/2., self.Rce, self.Zce)
            g_z2 = greenFunction(self.Rce[i], self.Zce[i] - dz/2., self.Rce, self.Zce)
            tar21 = self.tatmp*(g_r1 - g_r2)
            taz21 = self.tatmp*(g_z1 - g_z2)
            for j in range(ielem):
                jjj = self.equivforce[j]
                if jjj != iii:
                    self.dGreeniRj[iii, jjj] += prefac_r*tar21[j]
                    self.dGreeniZj[iii, jjj] += prefac_z*taz21[j]

            gr1 = greenFunction(self.Rce[i] + dr/2., self.Zce[i], self.Rg, self.Zg)
            gr2 = greenFunction(self.Rce[i] - dr/2., self.Zce[i], self.Rg, self.Zg)
            gz1 = greenFunction(self.Rce[i], self.Zce[i] + dz/2., self.Rg, self.Zg)
            gz2 = greenFunction(self.Rce[i], self.Zce[i] - dz/2., self.Rg, self.Zg)
            self.dGreeniRpl[:, :, iii] += prefac_r*(gr1 - gr2)
            self.dGreeniZpl[:, :, iii] += prefac_z*(gz1 - gz2)

                        
    def ferromagnets(self):

        if hasattr(self, 'Lferro'):
            x1 = self.Lferro/self.Rcurvferro
            x6 = r0ferro(i) - self.rcurvferro*np.cos(self.angleferro)
            x7 = z0ferro(i) - self.rcurvferro*np.sin(self.angleferro)
            r1 = z1/self.nferrosub
            n_ferro_mag = self.magValue.shape[1]
            n_sub_mag = np.max(self.nferrosub)
            self.r_ferro   = np.zeros((n_sub_mag, n_ferro_mag), dtype=gr_flt)
            self.z_ferro   = np.zeros_like(self.r_ferro)
            self.ang_ferro = np.zeros_like(self.r_ferro) # tananglferro
            self.len_ferro = np.zeros_like(self.r_ferro)
            self.mferro_ferro = np.zeros((n_sub_mag, n_sub_mag, n_ferro_mag), dtype=gr_flt)
            for i in range(self.n_ferro_mag):
                nfsi = self.nferrosub[i]
                ang = self.angleferro[i] - x1[i]/2. + r1[i]/2. + r1[i]*np.arange(nfsi)
                self.r_ferro  [:nfsi, i] = x6[i] + x3[i]*np.cos(ang)
                self.z_ferro  [:nfsi, i] = x7[i] + x3[i]*np.sin(ang)
                self.tan_ferro[:nfsi, i] = ang
                self.len_ferro[:nfsi, i] = r1[i]*self.rcurvferro[i]

                r_out = np.outer(self.r_ferro[:, i], self.r_ferro[:, i])
                z_out = np.outer(self.z_ferro[:, i], self.z_ferro[:, i])
                y1 = np.hypot(r_out, z_out)
                y2 = np.arctan2(z_out, r_out)
                cos_ferro = np.cos(self.ang_ferro[:, i])
                sin_ferro = np.sin(self.ang_ferro[:, i])
                y3 = np.cos(y2)*cos_ferro[None, :] - np.sin(y2)*sin_ferro[None, :]
                self.mferro_ferro[:, :, i] = y3/y1


    def dumpMachineJson(self, f_out='aug_description_full.json'):

        logger.debug('Dumping %s', f_out)
        nConduc  = self.indConduc.shape[0]
        data = { 'Rmin': self.Rgrid[0], 'Rmax': self.Rgrid[-1],
                 'Zmin': self.Zgrid[0], 'Zmax': self.Zgrid[-1] }
        data['lim_maxR'], data['lim_minR'], data['lim_maxZ'], data['lim_minZ'] = self.limRZ
        for lbl in ( 'alpsep',
            'R_coil', 'Z_coil', 'dR_coil', 'dZ_coil', 'angh_coil', 'ang_coil', 'm_equiv',
            'Rlim', 'Zlim', 'resConduc', 'resConduc_diag',
            'indConduc', 'zLimPotential', 'greenBnd', 'greeni',
            'dGreeniRj', 'dGreeniZj', 'dGreeniRpl', 'dGreeniZpl'):
            data[lbl] = getattr(self, lbl)
        data['R_cond'] = self.R_cond[self.nActive: nConduc]
        data['Z_cond'] = self.Z_cond[self.nActive: nConduc]
        with open(f_out, 'w') as f:
            f.write(format_numpy_json(data))
        logger.info('Stored machine file %s', f_out)


def write_green(f_in, f_out):
    gm = GREEN_MATRICES()
    gm.fromMachineInput(f_in)
    if hasattr(gm, 'Rmin'):
        gm.setCoilProperties()
        gm.calcGreenMatrices()
        gm.dumpMachineJson(f_out=f_out)

    
def main():

    for f_in in os.listdir(grIOdir):
        if os.path.splitext(f_in)[1] == '.json':
            tok = f_in.split('_')[0]
            f_machineIn  = '%s/%s_description_in.json'  %(grIOdir, tok)
            f_machineOut = '%s/%s_description_full.json' %(grIOdir, tok)

            if os.path.isfile(f_machineIn):
                if os.path.isfile(f_machineOut):
                    t_out = os.stat(f_machineOut).st_mtime
                    t_in  = os.stat(f_machineIn).st_mtime
                    t_py1 = os.stat('%s/greenMatrices.py' %loc_dir).st_mtime
                    fsize = os.path.getsize(f_machineOut)
                    if fsize == 0:
                        logger.info('File %s exists, but it has zero size. Removing', f_machineOut)
                        write_green(f_machineIn, f_machineOut)
                    elif t_in > t_out:
                        logger.info('Input file %s newer than output %s', f_machineIn, f_machineOut)
                        write_green(f_machineIn, f_machineOut)
                    elif t_py1 > t_out:
                        logger.info('exe/greenMatrifces.py newer than output file %s', f_machineOut)
                        write_green(f_machineIn, f_machineOut)
                else:
                    write_green(f_machineIn, f_machineOut)


if __name__ == '__main__':
    main()
