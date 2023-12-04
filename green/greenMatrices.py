#!/usr/bin/env python

import os, logging, argparse
import numpy as np
import matplotlib.pylab as plt
from scipy.interpolate import interp1d
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


def to_float(line):
    return float(line.split('!')[0].replace('d', 'e'))


def truncate(arr, ncols=3):
    nx = np.prod(arr.shape)
    nrows = nx//ncols
    arr_flat = arr.ravel()
    block = arr_flat[:ncols*nrows].reshape((nrows, ncols))
    return block, arr_flat[ncols*nrows:]


class GREEN_MATRICES:


    def __init__(self):
        pass


    def fromMachineInput(self, f_machine):

        logger.info('Reading %s', f_machine)

        with open(f_machine, 'r') as f:
            lines = f.readlines()

# Get blocks
        general_block = []
        active_block  = []
        res_block     = []
        limiter_block = []
        passive_block = []
        blanbpc_block = []

        j_block = 0
        for lin in lines:
            line = lin.split('!')[0].strip()
            if not line:
                j_block += 1
                continue
            if j_block == 1:
                general_block.append(line)
            elif j_block == 2:
                active_block.append(line)
            elif j_block == 3:
                res_block.append(line)
            elif j_block == 4:
                limiter_block.append(line)
            elif j_block == 5:
                passive_block.append(line)
            elif j_block == 6:
                blanbpc_block.append(line)

# General settings

        nR = int(general_block[0]) - 2
        nZ = int(general_block[1]) - 2
        Rmin = to_float(general_block[2])
        Rmax = to_float(general_block[3])
        Zmin = to_float(general_block[4])
        Zmax = to_float(general_block[5])
        self.alpsep = to_float(general_block[6])

        nR2 = nR + 2
        nZ2 = nZ + 2
        nR1 = nR + 1
        nZ1 = nZ + 1

        self.Rgrid = np.linspace(Rmin, Rmax, nR2)
        self.Zgrid = np.linspace(Zmin, Zmax, nZ2)

        dr = self.Rgrid[1] - self.Rgrid[0]
        dz = self.Zgrid[1] - self.Zgrid[0]

# Coils geometry

        self.n_elem_coil = np.array([], dtype=np.int32)
        self.m_turns  = np.array([], dtype=np.int32)
        self.m_equiv  = np.array([], dtype=np.int32)
        self.R_coil   = np.array([], dtype=np.float32)
        self.Z_coil   = np.array([], dtype=np.float32)
        self.dR_coil  = np.array([], dtype=np.float32)
        self.dZ_coil  = np.array([], dtype=np.float32)
        self.angh_coil= np.array([], dtype=np.float32)
        self.ang_coil = np.array([], dtype=np.float32)

#        n_coils = int(active_block[0])
        for line in active_block[1:]:
            try:
                self.n_elem_coil = np.append(self.n_elem_coil, int(line))
            except:
                rc, zc, drc, dzc, angh, ang, mt, me = line.split()
                self.R_coil   = np.append(self.R_coil   , float(rc))
                self.Z_coil   = np.append(self.Z_coil   , float(zc))
                self.dR_coil  = np.append(self.dR_coil  , float(drc))
                self.dZ_coil  = np.append(self.dZ_coil  , float(dzc))
                self.angh_coil= np.append(self.angh_coil, float(ang))
                self.ang_coil = np.append(self.ang_coil , float(ang))
                self.m_turns  = np.append(self.m_turns, int(mt))
                self.m_equiv  = np.append(self.m_equiv, int(me))

        self.angh_coil = np.radians(self.angh_coil)
        self.ang_coil  = np.radians(self.ang_coil)
        nCoils = len(self.R_coil)
        numeqcump = np.zeros(nCoils, dtype=np.int32)
        self.R_cond = np.zeros_like(self.R_coil)
        self.Z_cond = np.zeros_like(self.R_coil)
        for jcoil in range(nCoils):
            jequiv = self.m_equiv[jcoil] - 1
            self.R_cond[jequiv] += self.R_coil[jcoil]
            self.Z_cond[jequiv] += self.Z_coil[jcoil]
            numeqcump[jequiv] += 1

        self.nActive = np.max(self.m_equiv)
        self.R_cond = self.R_cond[:self.nActive] / numeqcump[:self.nActive]
        self.Z_cond = self.Z_cond[:self.nActive] / numeqcump[:self.nActive]

# Coil resistivity

        n_res = res_block[0]
        resConduc = []
        for line in res_block[1:]:
            resConduc.append([float(x) for x in line.split()])
        n_res = len(resConduc)
        n_res_max = 300
        self.resConduc = np.zeros((n_res_max, n_res_max), dtype=np.float32)
        self.resConduc[:n_res, :n_res] = np.array(resConduc, dtype=np.float32)

# Limiter geometry

        n_lim = limiter_block[0]
        Rlim = []
        Zlim = []
        for line in limiter_block[1:]:
            pieces = line.split()
            try:
                a, b = pieces
                Rlim.append(float(a))
                Zlim.append(float(b))
            except:
                lim_maxR, lim_minR, lim_maxZ, lim_minZ = (float(x) for x in pieces)

        Rlim1 = np.array(Rlim)
        Zlim1 = np.array(Zlim)
        indR = ((Rlim1 - Rmin)/dr + 0.5).astype(int)
        indZ = ((Zlim1 - Zmin)/dz + 0.5).astype(int)
        self.Rlim = self.Rgrid[indR]
        self.Zlim = self.Zgrid[indZ]

        self.ilim_maxR = int((lim_maxR - Rmin)/dr + 0.5)
        self.ilim_minR = int((lim_minR - Rmin)/dr + 0.5)
        self.ilim_maxZ = int((lim_maxZ - Zmin)/dz + 0.5)
        self.ilim_minZ = int((lim_minZ - Zmin)/dz + 0.5)
        self.lim_maxR = self.Rgrid[self.ilim_maxR]
        self.lim_minR = self.Rgrid[self.ilim_minR]
        self.lim_maxZ = self.Zgrid[self.ilim_maxZ]
        self.lim_minZ = self.Zgrid[self.ilim_minZ]
        self.zLimPotential = np.zeros((nR2, nZ2))
        self.zLimPotential[self.ilim_minR: self.ilim_maxR, self.ilim_minZ: self.ilim_maxZ] = 1.

# Blanket

        res_blan, width_blan = (to_float(x) for x in passive_block[0].split())
        n_blan = int(passive_block[1].split()[0])
        nConduc = self.nActive

        if n_blan > 0:
            x1 = np.array([])
            x2 = np.array([])
            x3 = np.array([])
            x4 = np.array([])
            x5 = np.array([])
            x6 = np.array([])
            for line in passive_block[2:]:
                val = line.split()
                x1 = np.append(x1, to_float(val[1]))
                x2 = np.append(x2, to_float(val[2]))
                x3 = np.append(x3, to_float(val[3]))
                x4 = np.append(x4, to_float(val[4]))
                x5 = np.append(x5, to_float(val[5]))
                x6 = np.append(x6, to_float(val[6]))

            x7 = x3 - x1
            x8 = x4 - x2
            x9 = np.hypot(x7, x8)
            self.R_blan = np.ravel([x1 + 0.25*x7, x1 + 0.75*x7], 'F')
            self.Z_blan = np.ravel([x2 + 0.25*x8, x2 + 0.75*x8], 'F')
            self.dhoriz = np.ravel([0.5*x7, 0.5*x7], 'F')
            self.dvert  = np.ravel([0.5*x8, 0.5*x8], 'F')
            area_blan   = np.ravel([width_blan*x9, width_blan*x9], 'F')
            ssfw = np.sum(area_blan/self.R_blan)
            n_blanket = len(self.R_blan)

            self.R_cond = np.append(self.R_cond[:self.nActive], self.R_blan)
            self.Z_cond = np.append(self.Z_cond[:self.nActive], self.Z_blan)
            for jblan in range(n_blanket):
                jcond = nConduc + jblan
                self.resConduc[jcond, jcond] = res_blan*self.R_blan[jblan]/area_blan[jblan]*ssfw
            nConduc += n_blanket

#
        n_blanket_pc = int(blanbpc_block[0])
        self.R_blan_pc    = np.array([])
        self.Z_blan_pc    = np.array([])
        self.area_blan_pc = np.array([])
        res_blan_pc  = np.array([])
        cur_blan_pc  = np.array([])
        if n_blanket_pc > 0:
            for line in blanpc_block[2:]:
                val = line.split()
                self.R_blan_pc    = np.append( self.R_blan_pc   , to_float(val[0]))
                self.Z_blan_pc    = np.append( self.Z_blan_pc   , to_float(val[1]))
                self.area_blan_pc = np.append( self.area_blan_pc, to_float(val[3]))
                res_blan_pc  = np.append(res_blan_pc , to_float(val[2]))
                cur_blan_pc  = np.append(cur_blan_pc , to_float(val[4]))
            self.R_cond = np.append(self.R_cond, self.R_blan_pc)
            self.Z_cond = np.append(self.Z_cond, self.Z_blan_pc)
            for jblan in range(n_blanket_pc):
                jcond = nConduc + jblan
                self.resConduc[jcond, jcond] = self.res_blan_pc[jblan]
            nConduc += n_blanket_pc

        self.nPassive = nConduc - self.nActive # n_blanket, n_blanket_c
        logger.debug('nactive, npassive, ncoils, nconduc, ssfw %d %d %d %d %12.4e', self.nActive, self.nPassive, nCoils, nConduc, ssfw)


    def calcGreenf(self):

#------
# Coils
#------

        nBlocks = 0
        nConduc = 0
        nCoils = len(self.R_coil)
        if hasattr(self, 'R_blan_pc'):
            n_blanket_pc = len(self.R_blan_pc)

        identcoil = np.ones(nCoils, dtype=bool)
        sin_coil = np.sin(self.ang_coil)
        cos_coil = np.cos(self.ang_coil)

        r1 = self.R_coil - 0.5*( self.dR_coil + self.dZ_coil*cos_coil/sin_coil)
        z1 = self.Z_coil - 0.5*self.dZ_coil
        x12 = self.dR_coil*sin_coil/self.dZ_coil
        indi = np.sqrt(x12   *self.n_elem_coil).astype(np.int32) + 1
        indj = np.sqrt(1./x12*self.n_elem_coil).astype(np.int32) + 1
        dr1 = self.dR_coil/indi        # x3 ; x4 = 0.
        dz1 = self.dZ_coil/indj        # x6
        dz_cs = dz1*cos_coil/sin_coil  # x5
        r1 += 0.5*(dr1 + dz_cs)
        z1 += 0.5*dz1

        Rce  = np.array([])
        Zce  = np.array([])
        dRce = np.array([])
        dZce = np.array([])
        tatmp   = np.array([])
        equivforce = np.array([], dtype=np.int32)
        equivtmp   = np.array([], dtype=np.int32)

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
                    flt0_ij = np.zeros(indij, dtype=np.float32)
                    int0_ij = np.zeros(indij, dtype=np.int32)
                    dRce = np.append(dRce, dr1[j] + flt0_ij) # lot of redundancy, reduce!
                    dZce = np.append(dZce, dz1[j] + flt0_ij)
                    tatmp    = np.append(tatmp   , self.m_turns[j]/float(indij) + flt0_ij)
                    equivforce = np.append(equivforce, j + 1    + int0_ij)
                    equivtmp   = np.append(equivtmp  , nConduc + int0_ij)
                    for jj in range(indj[j]):
                        Rce = np.append(Rce, r1[j] + dr1[j]*np.arange(indi[j]) + jj*dz_cs[j])
                        Zce = np.append(Zce, z1[j] +        np.zeros (indi[j]) + jj*dz1  [j])
        nctype = 2 + np.zeros(len(Rce), dtype=np.int32)

#--------
# Blanket
#--------

        if hasattr(self, 'R_blan_pc'):
            n_blanket_pc = len(self.R_blan_pc)
            n_elem_blanket_pc = 9 # sub element of a blanket element, hardwired to 9 for now
            ind_pc = np.round(np.sqrt(n_elem_blanket_pc + 1.e-6)).astype(np.int32)
            dx = 0.5*np.sqrt(self.area_blan_pc/ind_pc)
            dx_area = dx**2/self.area_blan_pc
            r1 = self.R_blan_pc - 0.5*dx
            z1 = self.Z_blan_pc - 0.5*dx
            flt0_iipc = np.zeros(ind_pc**2*n_blanket_pc, dtype=np.float32)
            int0_iipc = np.zeros(ind_pc**2*n_blanket_pc, dtype=np.int32)
            dRce  = np.append(dRce, dx + flt0_iipc)
            dZce  = np.append(dZce, dx + flt0_iipc)
            nctype   = np.append(nctype , 2  + int0_iipc)
            for j in range(n_blanket_pc):
                nConduc += 1
                nBlocks += 1
                flt0_ii = np.zeros(indi_pc**2, dtype=np.float32)
                int0_ii = np.zeros(indi_pc**2, dtype=np.int32)
                equivtmp   = np.append(equivtmp  , nConduc + int0_ii)
                equivforce = np.append(equivforce, nBlocks + int0_ii)
                dRce = np.append(dRce, dx[j]      + flt0_ii)
                dZce = np.append(dZce, dx[j]      + flt0_ii)
                tatmp   = np.append(tatmp  , dx_area[j] + flt0_ii)
                for jjj in range(ind_pc):
                    Zce = np.append(Zce, self.Z_blan_pc[j] + jjj*dx[j] + np.zeros(ind_pc))
                    Rce = np.append(Rce, self.R_blan_pc[j] + np.arange(ind_pc)*dx[j])

        if hasattr(self, 'R_blan'):
            n_blanket = len(self.R_blan)
            Rce  = np.append(Rce , self.R_blan)
            Zce  = np.append(Zce , self.Z_blan)
            dRce = np.append(dRce, self.dhoriz)
            dZce = np.append(dZce, self.dvert)
            nctype   = np.append(nctype, np.ones(n_blanket, dtype=np.int32  ))
            tatmp    = np.append(tatmp , np.ones(n_blanket, dtype=np.float32))
            equivtmp   = np.append(equivtmp  , nConduc + np.arange(n_blanket) + 1)
            equivforce = np.append(equivforce, nBlocks + np.arange(n_blanket) + 1)
            nConduc += n_blanket
            nBlocks += n_blanket

        ielem = len(tatmp)

# Self inductances
        self.indConduc = np.zeros((nConduc, nConduc))

        gf_diag = gf.identity(Rce, dRce, dZce, nctype)

        for i in range(ielem):
            iii = equivtmp[i] - 1
            self.indConduc[iii, iii] += gf_diag[i]*tatmp[i]**2/2.
            gf_ta = tatmp*gf.nonIdentity(Rce[i], Zce[i], Rce, Zce, dRce[i], dZce[i], dRce, dZce, nctype[i], nctype)
            for j in range(ielem):
                jjj = equivtmp[j] - 1
                if j != i:
                    self.indConduc[iii, jjj] += gf_ta[j]*tatmp[i]

        self.indConduc *= 0.8*np.pi

# Grid inductances
        nR2 = len(self.Rgrid)
        nZ2 = len(self.Zgrid)

        ntype2range = {1: (-1, 0, 1), 2: (0,)}
        Rg, Zg = np.meshgrid(self.Rgrid, self.Zgrid)

        self.greeni = np.zeros((nR2, nZ2, nConduc))
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

# Calculate inter-block forces (check: are they right?)

        self.dGreeniRpl = np.zeros((nR2, nZ2, nBlocks))
        self.dGreeniZpl = np.zeros((nR2, nZ2, nBlocks))

        self.dGreeniRj = np.zeros((nBlocks, nBlocks))
        self.dGreeniZj = np.zeros((nBlocks, nBlocks))

        dr = self.Rgrid[1] - self.Rgrid[0]
        dz = self.Zgrid[1] - self.Zgrid[0]

        for i in range(ielem):
            iii = equivforce[i] - 1
            prefac_r = -0.8*np.pi*tatmp[i]/dr
            prefac_z = -0.8*np.pi*tatmp[i]/dz
            
            xr1 = gf.greenFunction(Rce[i] + dr/2., Zce[i], Rce, Zce)
            xr2 = gf.greenFunction(Rce[i] - dr/2., Zce[i], Rce, Zce)
            xz1 = gf.greenFunction(Rce[i], Zce[i] + dz/2., Rce, Zce)
            xz2 = gf.greenFunction(Rce[i], Zce[i] - dz/2., Rce, Zce)
            for j in range(ielem):
                jjj = equivforce[j] - 1
                if jjj != iii:
                    self.dGreeniRj[iii, jjj] += prefac_r*tatmp[j]*(xr1[j] - xr2[j])
                    self.dGreeniZj[iii, jjj] += prefac_z*tatmp[j]*(xz1[j] - xz2[j])
            
            xr1 = gf.greenFunction(Rce[i] + dr/2., Zce[i], Rg, Zg)
            xr2 = gf.greenFunction(Rce[i] - dr/2., Zce[i], Rg, Zg)
            xz1 = gf.greenFunction(Rce[i], Zce[i] + dz/2., Rg, Zg)
            xz2 = gf.greenFunction(Rce[i], Zce[i] - dz/2., Rg, Zg)
            self.dGreeniRpl[:, :, iii] += prefac_r*(xr1 - xr2)
            self.dGreeniZpl[:, :, iii] += prefac_z*(xz1 - xz2)

# Boundary

        self.greenBnd = gf.greenBoundary(self.Rgrid, self.Zgrid)


    def dumpMachineDescr(self, f_out='machine_description_out.aug'):

        logger.debug('Dumping %s', f_out)
        nR2 = len(self.Rgrid)
        nZ2 = len(self.Zgrid)
        nR1 = nR2 - 1
        nZ1 = nZ2 - 1
        nR  = nR1 - 1
        nZ  = nZ1 - 1

        nLimiter = len(self.Rlim)
        nCoils   = len(self.R_coil)
        nBlocks  = self.dGreeniRj.shape[0]
        nConduc  = self.indConduc.shape[0]

        with open(f_out, 'w') as f:
            f.write('%3d %3d %3d\n' %(nR, nR2, nR1))
            f.write('%3d %3d %3d\n' %(nZ, nZ2, nZ1))
            f.write('%11.8f\n' %self.Rgrid [0])
            f.write('%11.8f\n' %self.Rgrid[-1])
            f.write('%11.8f\n' %self.Zgrid [0])
            f.write('%11.8f\n' %self.Zgrid[-1])
            f.write('%11.8f\n' %self.alpsep)
            f.write('%d %d\n' %(self.nActive, self.nPassive))

            f.write('%d\n' %nCoils)
            np.savetxt(f, np.c_[self.R_coil, self.Z_coil, self.dR_coil, self.dZ_coil, self.angh_coil, self.ang_coil, self.m_equiv], fmt='%15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %d')

            f.write('%d\n' %nLimiter)
            np.savetxt(f, np.c_[self.Rlim, self.Zlim], fmt='%11.8f %11.8f')

            f.write('%3d %11.8f\n' %(self.ilim_maxR+1, self.lim_maxR))
            f.write('%3d %11.8f\n' %(self.ilim_minR+1, self.lim_minR))
            f.write('%3d %11.8f\n' %(self.ilim_maxZ+1, self.lim_maxZ))
            f.write('%3d %11.8f\n' %(self.ilim_minZ+1, self.lim_minZ))

            np.savetxt(f, np.c_[self.R_cond[self.nActive: self.nPassive], self.Z_cond[self.nActive: self.nPassive]], fmt='%11.8f %11.8f')

            f.write('%d\n' %nConduc)
            for jcon in range(nConduc):
                block, tail = truncate(self.indConduc[jcon, :])
                np.savetxt(f, block, fmt='%15.8e')
                np.savetxt(f, tail , fmt='%15.8e')

            f.write('%d\n' %nConduc)
            for jcon in range(nConduc):
                block, tail = truncate(self.resConduc[jcon, :nConduc])
                np.savetxt(f, block, fmt='%15.8e')
                np.savetxt(f, tail, fmt='%15.8e')

            for jcon in range(nConduc):
                for jr in range(nR2):
                    block, tail = truncate(self.greeni[jr, :, jcon])
                    np.savetxt(f, block, fmt='%15.8e')
                    np.savetxt(f, tail , fmt='%15.8e')

            f.write('%d\n' %nBlocks)
            for jb in range(nBlocks):
                np.savetxt(f, np.c_[self.dGreeniRj[:, jb], self.dGreeniZj[:, jb]], fmt='%15.8e')
            for jb in range(nBlocks):
                for jz in range(nZ2):
                    np.savetxt(f, np.c_[self.dGreeniRpl[:, jz, jb], self.dGreeniZpl[:, jz, jb]], fmt='%15.8e')

            np.savetxt(f, self.zLimPotential.ravel(), fmt='%15.8e')

            grBnd = self.greenBnd.ravel()
            nRZ2 = len(grBnd)
            f.write('%d\n' %nRZ2)
            np.savetxt(f, grBnd, fmt='%15.8e')
        logger.info('Stored %s', f_out)


if __name__ == '__main__':


    parser = argparse.ArgumentParser(description='Write Green matrices for FEQIS')
    parser.add_argument('-t', '--tok', help='tokamak name', required=False, default='aug')
    args = parser.parse_args()

    f_machineIn  = '%s/machine_description_in.%s'  %(grIOdir, args.tok)
    f_machineOut = '%s/green/machine_description_out.%s' %(awd, args.tok)
    
    gm = GREEN_MATRICES()
    gm.fromMachineInput(f_machineIn)
    gm.calcGreenf()
    gm.dumpMachineDescr(f_out=f_machineOut)
