#!/usr/bin/env python

import os, logging, argparse
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


def strip_line(line):
    return line.split('!')[0].strip()


def split_line(line):
    return strip_line(line).split()


def to_int(line):
    return int(strip_line(line))


def to_float(line):
    return float(strip_line(line).replace('d', 'e'))


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

# Get blocks from file

        linesDelim = {}
        labels = ['grid', 'active', 'res', 'lim', 'blan', 'blanpc']
        for key in labels:
            linesDelim[key] = [0, 0]

        line_old = ''
        j_block = 0
        for jline, lin in enumerate(lines):
            line = strip_line(lin)            
            if line and not line_old:
                linesDelim[labels[j_block]][0] = jline
            if not line and line_old:
                linesDelim[labels[j_block]][1] = jline
                j_block += 1
            line_old = line

# Spatial grids

        jBeg, jEnd = linesDelim['grid']
        nR2  = to_int(lines[jBeg]  )
        nZ2  = to_int(lines[jBeg+1])
        Rmin = to_float(lines[jBeg+2])
        Rmax = to_float(lines[jBeg+3])
        Zmin = to_float(lines[jBeg+4])
        Zmax = to_float(lines[jBeg+5])
        self.alpsep = to_float(lines[jBeg+6])

        self.Rgrid = np.linspace(Rmin, Rmax, nR2, endpoint=True, dtype=gr_flt)
        self.Zgrid = np.linspace(Zmin, Zmax, nZ2, endpoint=True, dtype=gr_flt)
        dr = (Rmax - Rmin)/float(nR2 - 1)
        dz = (Zmax - Zmin)/float(nZ2 - 1)

# Coils geometry

        jBeg, jEnd = linesDelim['active']
        n_coils = to_int(lines[jBeg])
        self.R_coil, self.Z_coil, self.dR_coil, self.dZ_coil, self.angh_coil, self.ang_coil, \
            self.m_turns, self.m_equiv, self.n_elem_coil = \
            np.genfromtxt(f_machine, unpack=True, skip_header=jBeg+1, max_rows=n_coils, \
            dtype=6*[np.float32] + 3*[np.int32])

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

        jBeg, jEnd = linesDelim['res']
        n_res = to_int(lines[jBeg])
        self.resConduc = np.loadtxt(f_machine, skiprows=jBeg+1, max_rows=n_res, dtype=gr_flt)
        self.resConduc_diag = empty_fltarr

# Limiter geometry

        jBeg, jEnd = linesDelim['lim']
        n_lim = to_int(lines[jBeg])
        Rlim, Zlim = np.loadtxt(f_machine, unpack=True, skiprows=jBeg+1, max_rows=n_lim, dtype=gr_flt)
        lim_maxR, lim_minR, lim_maxZ, lim_minZ = (float(x) for x in split_line(lines[jEnd-1]))
        indR = ((Rlim - Rmin)/dr + 0.5).astype(int)
        indZ = ((Zlim - Zmin)/dz + 0.5).astype(int)
        self.Rlim = self.Rgrid[indR]
        self.Zlim = self.Zgrid[indZ]

        ilim_maxR = int((lim_maxR - Rmin)/dr + 0.5)
        ilim_minR = int((lim_minR - Rmin)/dr + 0.5)
        ilim_maxZ = int((lim_maxZ - Zmin)/dz + 0.5)
        ilim_minZ = int((lim_minZ - Zmin)/dz + 0.5)
        self.lim_maxR = self.Rgrid[ilim_maxR]
        self.lim_minR = self.Rgrid[ilim_minR]
        self.lim_maxZ = self.Zgrid[ilim_maxZ]
        self.lim_minZ = self.Zgrid[ilim_minZ]
        self.zLimPotential = np.zeros((nR2, nZ2), dtype=gr_int)
        self.zLimPotential[ilim_minR: ilim_maxR, ilim_minZ: ilim_maxZ] = 1

# Blanket

        jBeg, jEnd = linesDelim['blan']
        res_blan, width_blan = (to_float(x) for x in split_line(lines[jBeg]))
        n_blan = int(lines[jBeg+1].split()[0])
        if n_blan > 0:
            x1, x2, x3, x4 = np.loadtxt(f_machine, unpack=True, usecols=(1, 2, 3, 4), skiprows=jBeg+2, max_rows=n_blan, dtype=gr_flt)

        nConduc = self.nActive
        if n_blan > 0:
            x7 = x3 - x1
            x8 = x4 - x2
            x9 = np.hypot(x7, x8)
            self.R_blan = np.ravel([x1 + 0.25*x7, x1 + 0.75*x7], order='F')
            self.Z_blan = np.ravel([x2 + 0.25*x8, x2 + 0.75*x8], order='F')
            self.dhoriz = np.ravel([0.5*x7, 0.5*x7], order='F')
            self.dvert  = np.ravel([0.5*x8, 0.5*x8], order='F')
            area_blan   = np.ravel([width_blan*x9, width_blan*x9], order='F')
            ssfw = np.sum(area_blan/self.R_blan)
            n_blanket = len(self.R_blan)

            self.R_cond = np.append(self.R_cond[:self.nActive], self.R_blan)
            self.Z_cond = np.append(self.Z_cond[:self.nActive], self.Z_blan)
            for jblan in range(n_blanket):
                jcond = nConduc + jblan
                self.resConduc_diag = np.append(self.resConduc_diag, res_blan*self.R_blan[jblan]/area_blan[jblan]*ssfw)
            nConduc += n_blanket

# Blanketpc

        jBeg, jEnd = linesDelim['blanpc']
        n_blanket_pc = to_int(lines[jBeg])
        if n_blanket_pc > 0:
            self.R_blan_pc, self.Z_blan_pc, res_blan_pc, self.area_blan_pc = np.loadtxt(f_machine, unpack=True, usecols=(0, 1, 2, 3), skiprows=jBeg+1, max_rows=n_blanket_pc, dtype=gr_flt)

            self.R_cond = np.append(self.R_cond, self.R_blan_pc)
            self.Z_cond = np.append(self.Z_cond, self.Z_blan_pc)
            for jblan in range(n_blanket_pc):
                jcond = nConduc + jblan
                self.resConduc_diag = np.append(self.resConduc_diag, res_blan_pc[jblan])
            nConduc += n_blanket_pc

        logger.debug('nactive, ncoils, nconduc, ssfw %d %d %d %12.4e', self.nActive, nCoils, nConduc, ssfw)


    def calcGreenf(self):

#------
# Coils

        nBlocks = 0
        nConduc = 0
        nCoils = len(self.R_coil)

        identcoil = np.ones(nCoils, dtype=bool)
        sin_coil = np.sin(self.ang_coil)
        cos_coil = np.cos(self.ang_coil)

        x12 = self.dR_coil*sin_coil/self.dZ_coil
        indi = np.sqrt(x12   *self.n_elem_coil).astype(gr_int) + 1
        indj = np.sqrt(1./x12*self.n_elem_coil).astype(gr_int) + 1
        dr1 = self.dR_coil/indi        # x3 ; x4 = 0.
        dz1 = self.dZ_coil/indj        # x6
        dz_cs = dz1*cos_coil/sin_coil  # x5
        r1 = self.R_coil - 0.5*( self.dR_coil + self.dZ_coil*cos_coil/sin_coil) + 0.5*(dr1 + dz_cs)
        z1 = self.Z_coil - 0.5*self.dZ_coil + 0.5*dz1

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
                    dRce  = np.append(dRce , dr1[j] + flt0_ij) # lot of redundancy, reduce!
                    dZce  = np.append(dZce , dz1[j] + flt0_ij)
                    tatmp = np.append(tatmp, self.m_turns[j]/float(indij) + flt0_ij)
                    equivforce = np.append(equivforce, j + 1    + int0_ij)
                    equivtmp   = np.append(equivtmp  , nConduc + int0_ij)
                    for jj in range(indj[j]):
                        Rce = np.append(Rce, r1[j] + dr1[j]*np.arange(indi[j]) + jj*dz_cs[j])
                        Zce = np.append(Zce, z1[j] +        np.zeros (indi[j]) + jj*dz1  [j])
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

        nR2 = len(self.Rgrid)
        nZ2 = len(self.Zgrid)
        dr = (self.Rgrid[-1] - self.Rgrid[0])/float(nR2 - 1)
        dz = (self.Zgrid[-1] - self.Zgrid[0])/float(nZ2 - 1)

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

#-------------------
# Inter-block forces

        self.dGreeniRpl = np.zeros((nR2, nZ2, nBlocks), dtype=gr_flt)
        self.dGreeniZpl = np.zeros((nR2, nZ2, nBlocks), dtype=gr_flt)
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


    def dumpMachineDescr(self, f_out='machine_description_out.aug'):

        logger.debug('Dumping %s', f_out)
        nR2, nZ2, nBlocks = self.dGreeniRpl.shape
        nLimiter = len(self.Rlim)
        nCoils   = len(self.R_coil)
        nConduc  = self.indConduc.shape[0]
        nPassive = nConduc - self.nActive

        with open(f_out, 'w') as f:
            f.write('%3d %3d\n' %(nR2, nZ2))
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

            f.write('%11.8f\n' %self.lim_maxR)
            f.write('%11.8f\n' %self.lim_minR)
            f.write('%11.8f\n' %self.lim_maxZ)
            f.write('%11.8f\n' %self.lim_minZ)

            np.savetxt(f, np.c_[self.R_cond[self.nActive: nPassive], self.Z_cond[self.nActive: nPassive]], fmt='%11.8f %11.8f')

            f.write('%d\n' %nConduc)
            for jcon in range(nConduc):
                block, tail = truncate(self.indConduc[jcon, :])
                np.savetxt(f, block, fmt='%15.8e')
                np.savetxt(f, tail , fmt='%15.8e')

            f.write('%d %s\n' %(self.nActive, nConduc))
            np.savetxt(f, self.resConduc, fmt='%15.8e')
            np.savetxt(f, self.resConduc_diag, fmt='%15.8e')

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

            np.savetxt(f, self.zLimPotential.ravel(), fmt='%1d')

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
    f_machineOut = '%s/machine_description_out.%s' %(grIOdir, args.tok)
    
    gm = GREEN_MATRICES()
    gm.fromMachineInput(f_machineIn)
    gm.calcGreenf()
    gm.dumpMachineDescr(f_out=f_machineOut)
