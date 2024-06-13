#!/usr/bin/env python

import os, logging, argparse, json
import numpy as np
import green_functions as gf


fmt = logging.Formatter('%(asctime)s | %(name)s | %(levelname)s: %(message)s', '%H:%M:%S')
logger = logging.getLogger('green')
if len(logger.handlers) == 0:
    hnd = logging.StreamHandler()
    hnd.setFormatter(fmt)
    logger.addHandler(hnd)

logger.setLevel(logging.DEBUG)
#logger.setLevel(logging.INFO)

awd = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
grIOdir = '%s/exp/cnf' %awd

gr_flt = np.float64
gr_int = np.int32
empty_fltarr = np.array([], dtype=gr_flt)
empty_intarr = np.array([], dtype=gr_int)
expfmt = '%15.8e'


def truncate(arr, ncols=3):
    nx = np.prod(arr.shape)
    nrows = nx//ncols
    arr_flat = arr.ravel()
    block = arr_flat[:ncols*nrows].reshape((nrows, ncols))
    return block, arr_flat[ncols*nrows:]


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
        for jcoil in range(nCoils):
            jequiv = self.m_equiv[jcoil] - 1
            self.R_cond[jequiv] += self.R_coil[jcoil]
            self.Z_cond[jequiv] += self.Z_coil[jcoil]
            numeqcump[jequiv] += 1
        self.R_cond /= numeqcump
        self.Z_cond /= numeqcump

# Coil resistivity

        self.resConduc_diag = empty_fltarr

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

# Blanket
        nConduc = self.nActive
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
               self.resConduc_diag = np.append(self.resConduc_diag, self.res_blan*self.R_blan/area_blan*ssfw)
            nConduc += n_blanket

# Blanketpc

        if hasattr(self, 'R_blan_pc'):
            self.R_cond = np.append(self.R_cond, self.R_blan_pc)
            self.Z_cond = np.append(self.Z_cond, self.Z_blan_pc)
            if self.res_blan < 0:
                self.resConduc_diag = np.append(self.resConduc_diag, self.res_blan_pc)
            nConduc += len(self.R_blan_pc)

# Ferromagnet stuff

        if hasattr(self, 'r0ferro'):
            self.angleferro = np.radians(self.angleferro)
        
        logger.debug('nactive, ncoils, nconduc, ssfw %d %d %d %12.4e', self.nActive, nCoils, nConduc, ssfw)


    def calcGreenf(self):

#------
# Coils

        nBlocks = 0
        nConduc = 0
        nCoils = len(self.R_coil)

        identcoil = np.ones(nCoils, dtype=bool)
        sin_coil  = np.sin(self.ang_coil)
        cos_coil  = np.cos(self.ang_coil)
        sin_coilh = np.sin(self.angh_coil)
        cos_coilh = np.cos(self.angh_coil)

        x12 = self.dR_coil/self.dZ_coil
        indi = np.sqrt(x12   *self.n_elem_coil).astype(gr_int) + 1
        indj = np.sqrt(1./x12*self.n_elem_coil).astype(gr_int) + 1
        dr1 = self.dR_coil/indi        # x3 ; x4 = 0.
        dz1 = self.dZ_coil/indj        # x6
        dr_cs = dr1*cos_coilh  # x5
        dz_cs = dz1*cos_coil  # x5
        dr_ss = dr1*sin_coilh  # x5
        dz_ss = dz1*sin_coil  # x5
        r1 = self.R_coil - 0.5*( self.dR_coil*cos_coilh + self.dZ_coil*cos_coil) + 0.5*(dr_cs + dz_cs)
        z1 = self.Z_coil - 0.5*( self.dR_coil*sin_coilh + self.dZ_coil*sin_coil) + 0.5*(dr_ss + dz_ss)

        Rce   = empty_fltarr
        Zce   = empty_fltarr
        dRce  = empty_fltarr
        dZce  = empty_fltarr
        tatmp = empty_fltarr
        equivforce = empty_intarr
        equivtmp   = empty_intarr

        for icoil in range(nCoils):
            nBlocks += 1
            if identcoil[icoil]:
                nConduc += 1
                (ind, ) = np.where(self.m_equiv[icoil+1: nCoils] == self.m_equiv[icoil])
                ind += icoil + 1   # correct array index
                indk = np.append(icoil, ind)
                identcoil[ind] = 0
                for j in indk:
                    indij = indi[j]*indj[j]
                    flt0_ij = np.zeros(indij, dtype=gr_flt)
                    int0_ij = np.zeros(indij, dtype=gr_int)
                    dRce  = np.append(dRce , dr_cs[j] + flt0_ij) # lot of redundancy, reduce!
                    dZce  = np.append(dZce , dz_ss[j] + flt0_ij)
                    tatmp = np.append(tatmp, self.m_turns[j]/float(indij) + flt0_ij)
                    equivforce = np.append(equivforce, j + 1   + int0_ij)
                    equivtmp   = np.append(equivtmp  , nConduc + int0_ij)
                    for jj in range(indj[j]):
                        Rce = np.append(Rce, r1[j] + dr_cs[j]*np.arange(indi[j]) + jj*dz_cs[j])
                        Zce = np.append(Zce, z1[j] + dr_ss[j]*np.arange(indi[j]) + jj*dz_ss[j])
        nctype = 2 + np.zeros(len(Rce), dtype=gr_int)

#--------
# Blanket

        if hasattr(self, 'R_blan_pc'):
            n_blanket_pc = len(self.R_blan_pc)
            n_elem_blanket_pc = 9 # sub element of a blanket element, hardwired to 9 for now
            ind_pc = np.round(np.sqrt(n_elem_blanket_pc + 1.e-6)).astype(gr_int)
            dx = 0.5*np.sqrt(self.area_blan_pc/ind_pc)
            dx_area = dx**2/self.area_blan_pc
            r1 = self.R_blan_pc - 0.5*dx
            z1 = self.Z_blan_pc - 0.5*dx
            flt0_iipc = np.zeros(ind_pc**2*n_blanket_pc, dtype=gr_flt)
            int0_iipc = np.zeros(ind_pc**2*n_blanket_pc, dtype=gr_int)
            dRce   = np.append(dRce, dx + flt0_iipc)
            dZce   = np.append(dZce, dx + flt0_iipc)
            nctype = np.append(nctype , 2  + int0_iipc)
            for j in range(n_blanket_pc):
                nConduc += 1
                nBlocks += 1
                flt0_ii = np.zeros(indi_pc**2, dtype=gr_flt)
                int0_ii = np.zeros(indi_pc**2, dtype=gr_int)
                equivtmp   = np.append(equivtmp  , nConduc + int0_ii)
                equivforce = np.append(equivforce, nBlocks + int0_ii)
                dRce  = np.append(dRce , dx[j]      + flt0_ii)
                dZce  = np.append(dZce , dx[j]      + flt0_ii)
                tatmp = np.append(tatmp, dx_area[j] + flt0_ii)
                for jjj in range(ind_pc):
                    Zce = np.append(Zce, self.Z_blan_pc[j] + jjj*dx[j] + np.zeros(ind_pc))
                    Rce = np.append(Rce, self.R_blan_pc[j] + np.arange(ind_pc)*dx[j])

        if hasattr(self, 'R_blan'):
            n_blanket = len(self.R_blan)
            Rce  = np.append(Rce , self.R_blan)
            Zce  = np.append(Zce , self.Z_blan)
            dRce = np.append(dRce, self.dhoriz)
            dZce = np.append(dZce, self.dvert)
            nctype     = np.append(nctype, np.ones(n_blanket, dtype=gr_int))
            tatmp      = np.append(tatmp , np.ones(n_blanket, dtype=gr_flt))
            equivtmp   = np.append(equivtmp  , nConduc + np.arange(n_blanket) + 1)
            equivforce = np.append(equivforce, nBlocks + np.arange(n_blanket) + 1)
            nConduc += n_blanket
            nBlocks += n_blanket

        ielem = len(tatmp)

#-----------------
# Self inductances

        gf_diag = gf.identity(Rce, dRce, dZce, nctype)

        self.indConduc = np.zeros((nConduc, nConduc))
        for i in range(ielem):
            iii = equivtmp[i] - 1
            self.indConduc[iii, iii] += gf_diag[i]*tatmp[i]**2/2.
            gf_ta = tatmp*gf.nonIdentity(Rce[i], Zce[i], Rce, Zce, dRce[i], dZce[i], dRce, dZce, nctype[i], nctype)
            for j in range(ielem):
                if j != i:
                    jjj = equivtmp[j] - 1
                    self.indConduc[iii, jjj] += gf_ta[j]*tatmp[i]
        self.indConduc *= 0.8*np.pi

#-----------------
# Grid inductances

        nR = len(self.Rgrid)
        nZ = len(self.Zgrid)
        dr = (self.Rgrid[-1] - self.Rgrid[0])/float(nR - 1)
        dz = (self.Zgrid[-1] - self.Zgrid[0])/float(nZ - 1)

        ntype2range = {1: (-1, 0, 1), 2: (0,)}
        Rg, Zg = np.meshgrid(self.Rgrid, self.Zgrid)

        self.greeni = np.zeros((nR, nZ, nConduc))
        for i in range(ielem):
            iii = equivtmp[i] - 1
            irange = ntype2range[nctype[i]]
            gf_off = np.zeros_like(Rg)
            for iloc in irange: #nctype2=2
                gf_off += gf.greenFunction(Rce[i] + dRce[i]*iloc/3.,
                                           Zce[i] + dZce[i]*iloc/3.,
                                           Rg.T, Zg.T)
            self.greeni[:, :, iii] += gf_off*tatmp[i]
            self.greeni[:, :, iii] /= len(irange)
        self.greeni *= 0.4

#-------------------
# Inter-block forces

        self.dGreeniRpl = np.zeros((nR, nZ, nBlocks), dtype=gr_flt)
        self.dGreeniZpl = np.zeros((nR, nZ, nBlocks), dtype=gr_flt)
        self.dGreeniRj  = np.zeros((nBlocks, nBlocks) , dtype=gr_flt)
        self.dGreeniZj  = np.zeros((nBlocks, nBlocks) , dtype=gr_flt) 

        for i in range(ielem):
            iii = equivforce[i] - 1
            prefac_r = -0.8*np.pi*tatmp[i]/dr
            prefac_z = -0.8*np.pi*tatmp[i]/dz
            g_r1 = gf.greenFunction(Rce[i] + dr/2., Zce[i], Rce, Zce)
            g_r2 = gf.greenFunction(Rce[i] - dr/2., Zce[i], Rce, Zce)
            g_z1 = gf.greenFunction(Rce[i], Zce[i] + dz/2., Rce, Zce)
            g_z2 = gf.greenFunction(Rce[i], Zce[i] - dz/2., Rce, Zce)
            tar21 = tatmp*(g_r1 - g_r2)
            taz21 = tatmp*(g_z1 - g_z2)
            for j in range(ielem):
                jjj = equivforce[j] - 1
                if jjj != iii:
                    self.dGreeniRj[iii, jjj] += prefac_r*tar21[j]
                    self.dGreeniZj[iii, jjj] += prefac_z*taz21[j]

            gr1 = gf.greenFunction(Rce[i] + dr/2., Zce[i], Rg, Zg)
            gr2 = gf.greenFunction(Rce[i] - dr/2., Zce[i], Rg, Zg)
            gz1 = gf.greenFunction(Rce[i], Zce[i] + dz/2., Rg, Zg)
            gz2 = gf.greenFunction(Rce[i], Zce[i] - dz/2., Rg, Zg)
            self.dGreeniRpl[:, :, iii] += prefac_r*(gr1 - gr2)
            self.dGreeniZpl[:, :, iii] += prefac_z*(gz1 - gz2)

#---------
# Boundary

        self.greenBnd = gf.greenBoundary(self.Rgrid, self.Zgrid)

#--------------
# Ferromagnets

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


    def dumpMachineDescr(self, f_out='machine_description_out.aug'):

        logger.debug('Dumping %s', f_out)
        nR, nZ, nBlocks = self.dGreeniRpl.shape
        nLimiter = len(self.Rlim)
        nCoils   = len(self.R_coil)
        nConduc  = self.indConduc.shape[0]
        nPassive = nConduc - self.nActive

        with open(f_out, 'w') as f:
            f.write('%3d %3d\n' %(nR, nZ))
            f.write('%11.8f\n' %self.Rgrid [0])
            f.write('%11.8f\n' %self.Rgrid[-1])
            f.write('%11.8f\n' %self.Zgrid [0])
            f.write('%11.8f\n' %self.Zgrid[-1])
            f.write('%11.8f\n' %self.alpsep)
            f.write('%d %d\n' %(self.nActive, nPassive))

            f.write('%d\n' %nCoils)
            np.savetxt(f, np.c_[self.R_coil, self.Z_coil, self.dR_coil, self.dZ_coil, self.angh_coil, self.ang_coil, self.m_equiv], fmt='%15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %d')

            f.write('%d\n' %nLimiter)
            np.savetxt(f, np.c_[self.Rlim, self.Zlim], fmt='%11.8f %11.8f')
            f.write(4*'%11.8f\n' %self.limRZ)

            np.savetxt(f, np.c_[self.R_cond[self.nActive: nPassive+self.nActive], self.Z_cond[self.nActive: nPassive+self.nActive]], fmt='%11.8f %11.8f')

            f.write('%d\n' %nConduc)
            for jcon in range(nConduc):
                block, tail = truncate(self.indConduc[jcon, :])
                np.savetxt(f, block, fmt=expfmt)
                np.savetxt(f, tail , fmt=expfmt)

            f.write('%d %s\n' %(self.nActive, nConduc))
            np.savetxt(f, self.resConduc, fmt=expfmt)
            np.savetxt(f, self.resConduc_diag, fmt=expfmt)

            for jcon in range(nConduc):
                for jr in range(nR):
                    block, tail = truncate(self.greeni[jr, :, jcon])
                    np.savetxt(f, block, fmt=expfmt)
                    np.savetxt(f, tail , fmt=expfmt)

            f.write('%d\n' %nBlocks)
            for jb in range(nBlocks):
                np.savetxt(f, np.c_[self.dGreeniRj[:, jb], self.dGreeniZj[:, jb]], fmt=expfmt)
            for jb in range(nBlocks):
                for jz in range(nZ):
                    np.savetxt(f, np.c_[self.dGreeniRpl[:, jz, jb], self.dGreeniZpl[:, jz, jb]], fmt=expfmt)

            np.savetxt(f, self.zLimPotential.ravel(), fmt='%1d')

            grBnd = self.greenBnd.ravel()
            nRZ   = len(grBnd)
            f.write('%d\n' %nRZ)
            np.savetxt(f, grBnd, fmt=expfmt)
            if not hasattr(self, 'Lferro'):
                f.write('-1\n')
            else:
                n_ferro_mag = len(self.Lferro)
                chiValue = self.magValue/self.hValue
                f.write('%d\n' %n_ferro_mag)
                for i in range(n_ferro_mag):
                    f.write('%d %d\n' %(self.nferrosub[i], self.imagValues[i]))
                    np.savetxt(f, np.c_[chiValue[:, i], hValue[:, i]], fmt=expfmt)
                    np.savetxt(f, np.c_[self.r_ferro[:, i], self.z_ferro[:, i], self.ang_ferro[:, i], self.len_ferro[:, i]], fmt=expfmt)
                    np.savetxt(f, self.mferro_ferro[:, :, i], fmt=expfmt)

        logger.info('Stored %s', f_out)


if __name__ == '__main__':

#    parser = argparse.ArgumentParser(description='Write Green matrices for FEQIS')
#    parser.add_argument('-t', '--tok', help='tokamak name', required=False, default='aug')
#    args = parser.parse_args()

    for f_in in os.listdir(grIOdir):
        if os.path.splitext(f_in)[1] == '.json':
            tok = f_in.split('_')[0]
            f_machineIn  = '%s/%s_description_in.json'  %(grIOdir, tok)
            f_machineOut = '%s/machine_description_out.%s' %(grIOdir, tok)

            if os.path.isfile(f_machineOut):
                fsize = os.path.getsize(f_machineOut)
                if fsize == 0:
                    logger.info('File %s exists, but it has zero size. Removing', f_machineOut)
                    os.system('rm %s' %fsize)
            if not os.path.isfile(f_machineOut):
                gm = GREEN_MATRICES()
                gm.fromMachineInput(f_machineIn)
                gm.calcGreenf()
                gm.dumpMachineDescr(f_out=f_machineOut)
