import numpy as np
import matplotlib.pylab as plt
from scipy.interpolate import interp1d


def to_float(line):

    return float(line.split('!')[0].replace('d', 'e'))


def to_arr(list_in, typ=np.float32):
    return np.array(list_in, dtype=typ)


class FBE:


    def __init__(self):
        pass


    def from_astra_machine_file(self, f_machine):

        with open(f_machine, 'rb') as f:
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

        nR = int(general_block[0])
        nZ = int(general_block[1])
        Rmin   = to_float(general_block[2])
        Rmax   = to_float(general_block[3])
        Zmin   = to_float(general_block[4])
        Zmax   = to_float(general_block[5])
        alpsep = to_float(general_block[6])

        nR2 = nR + 2
        nZ2 = nZ + 2
        nR1 = nR + 1
        nZ1 = nZ + 1
        
        self.Rgrid = np.linspace(Rmin, Rmax, nR2)
        self.Zgrid = np.linspace(Zmin, Zmax, nZ2)

        dr = self.Rgrid[1] - self.Rgrid[0]
        dz = self.Zgrid[1] - self.Zgrid[0]

# Coils geometry

        n_elem_coil = []
        m_turns = []
        m_equiv = []
        R_coil   = []
        Z_coil   = []
        dR_coil  = []
        dZ_coil  = []
        ang_coil = []

#        n_coils = int(active_block[0])
        for line in active_block[1:]:
            try:
                n_elem_coil.append(int(line))
            except:
                rc, zc, drc, dzc, dummy, ang, mt, me = line.split()
                R_coil.append(float(rc))
                Z_coil.append(float(zc))
                dR_coil.append(float(drc))
                dZ_coil.append(float(dzc))
                ang_coil.append(float(ang))
                m_turns.append(int(mt))
                m_equiv.append(int(me) - 1) # -1 due to array labelling in python

        self.ang_coil = np.radians(to_arr(ang_coil))
        self.R_coil  = to_arr(R_coil)
        self.Z_coil  = to_arr(Z_coil)
        self.dR_coil = to_arr(dR_coil)
        self.dZ_coil = to_arr(dZ_coil)

        n_coils = len(self.R_coil)
        numeqcump = np.zeros(n_coils, dtype=np.int32)
        self.R_cond = np.zeros_like(self.R_coil)
        self.Z_cond = np.zeros_like(self.R_coil)
        for jcoil in range(n_coils):
            jequiv = m_equiv[jcoil]
            self.R_cond[jequiv] += self.R_coil[jcoil]
            self.Z_cond[jequiv] += self.Z_coil[jcoil]
            numeqcump[jequiv] += 1

        n_active = np.max(m_equiv) + 1
        self.R_cond = self.R_cond[:n_active] / numeqcump[:n_active]
        self.Z_cond = self.Z_cond[:n_active] / numeqcump[:n_active]
        self.cur_conduc = np.zeros(n_active, dtype=np.float32)

# Coil resistivity

        n_res = res_block[0]
        res_conduc = []
        for line in res_block[1:]:
            res_conduc.append([float(x) for x in line.split()])
        n_res = len(res_conduc)
        n_res_max = 300
        self.res_conduc = np.zeros((n_res_max,n_res_max), dtype=np.float32)
        self.res_conduc[:n_res, :n_res] = to_arr(res_conduc)

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

        self.Rlim1 = np.array(Rlim)
        self.Zlim1 = np.array(Zlim)
        indR = ((self.Rlim1 - Rmin)/dr + 0.5).astype(int)
        indZ = ((self.Zlim1 - Zmin)/dz + 0.5).astype(int)
        self.Rlim = self.Rgrid[indR]
        self.Zlim = self.Zgrid[indZ]

        ilim_maxR = int((lim_maxR - Rmin)/dr + 0.5)
        ilim_minR = int((lim_minR - Rmin)/dr + 0.5)
        ilim_maxZ = int((lim_maxZ - Zmin)/dz + 0.5)
        ilim_minZ = int((lim_minZ - Zmin)/dz + 0.5)
        lim_maxR = self.Rgrid[ilim_maxR]
        lim_minR = self.Rgrid[ilim_minR]
        lim_maxZ = self.Zgrid[ilim_maxZ]
        lim_minZ = self.Zgrid[ilim_minZ]

# Blanket

        res_blan, width_blan = (to_float(x) for x in passive_block[0].split())
        n_blan = int(passive_block[1].split()[0])
        n_conduc = n_active

        if n_blan > 0:
            x1 = []
            x2 = []
            x3 = []
            x4 = []
            x5 = []
            x6 = []
            for line in passive_block[2:]:
                val = line.split()
                x1.append(to_float(val[1]))
                x2.append(to_float(val[2]))
                x3.append(to_float(val[3]))
                x4.append(to_float(val[4]))
                x5.append(to_float(val[5]))
                x6.append(to_float(val[6]))

            x1 = to_arr(x1)
            x2 = to_arr(x2)
            x3 = to_arr(x3)
            x4 = to_arr(x4)
            x5 = to_arr(x5)
            x6 = to_arr(x6)
            x7 = x3 - x1
            x8 = x4 - x2
            x9 = np.hypot(x7, x8)
            self.R_blan = np.ravel([x1 + 0.25*x7, x1 + 0.75*x7], 'F')
            self.Z_blan = np.ravel([x2 + 0.25*x8, x2 + 0.75*x8], 'F')
            self.dhoriz = np.append(0.5*x7, 0.5*x7)
            self.dvert  = np.append(0.5*x8, 0.5*x8)
            area_blan = np.append(width_blan*x9, width_blan*x9)
            self.ssfw = np.cumsum(area_blan/self.R_blan)
            n_blanket = len(self.R_blan)
            self.R_cond = np.append(self.R_cond[:n_active], self.R_blan)
            self.Z_cond = np.append(self.Z_cond[:n_active], self.Z_blan)
            self.cur_conduc = np.append(slef.cur_conduc, np.zeros(n_blanket)
            for jblan in range(n_blanket):
                jcond = n_conduc + jblan
                self.res_conduc[jcond, jcond] = res_blan*self.R_blan[jblan]/area_blan[jblan]*self.ssfw[jblan]
            n_conduc += n_blanket

#

        n_blanket_pc = int(blanbpc_block[0])
        n_elem_blanket_pc = 9 # sub element of a blanket element, hardwired to 9 for now
        R_blan_pc = []
        Z_blan_pc = []
        res_blan_pc  = []
        area_blan_pc = []
        cur_blan_pc  = []
        if n_blanket_pc > 0:
            for line in blanpc_block[2:]:
                val = line.split()
                R_blan_pc.append(   val[0])
                Z_blan_pc.append(   val[1])
                res_blan_pc.append( val[2])
                area_blan_pc.append(val[3])
                cur_blan_pc.append( val[4])
            self.R_cond = np.append(self.R_cond, to_arr(R_blan_pc))
            self.Z_cond = np.append(self.Z_cond, to_arr(Z_blan_pc))
            self.cur_conduc = np.append(self.cur_conduc, to_arr(cur_blanpc))
            for jblan in range(n_blanket_pc):
                jcond = n_conduc + jblan
                self.res_conduc[jcond, jcond] = res_blan_pc[jblan]
            n_conduc += n_blanket_pc

        n_passive = n_conduc - n_active # n_blanket, n_blanket_c
        n_blocks = 0
        n_conduc = 0
        ielem = 0
        tatmp = 0.


if __name__ == '__main__':

    f_machine = 'exp/cnf/machine_description_in.aug'
    fbe = FBE()
    fbe.from_astra_machine_file(f_machine)

    print(fbe.ang_coil)

    fig = plt.figure('Limiter', (6, 9))
    fig.add_subplot(1, 1, 1, aspect='equal')
    plt.plot(fbe.Rlim , fbe.Zlim , 'b-')
    plt.plot(fbe.Rlim1, fbe.Zlim1, 'g-')
    plt.show()
