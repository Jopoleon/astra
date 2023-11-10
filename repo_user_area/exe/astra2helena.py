#!/usr/bin/env python

import numpy as np
from scipy.interpolate import splev, splrep
import os, argparse
from scipy.interpolate import InterpolatedUnivariateSpline, UnivariateSpline
from scipy.optimize import minimize
from scipy.io import netcdf_file

awd = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))


class astra2helena:


    def __init__(self, f_cdf, shottime=None, na1=451):

        self.shottime = shottime

        cv = netcdf_file(f_cdf, 'r', mmap=False).variables

        self.rs = cv['r2d'][-1, -1, :]
        self.zs = cv['z2d'][-1, -1, :]

        self.IPL  = cv['IPL'][-1]
        self.BTOR = cv['BTOR'][-1]
        RTOR = cv['RTOR' ][-1]
        self.Wi = cv['CDWM1'][-1]

        nibm  = cv['NIBM' ][-1, :]
        rhotN = cv['XRHO' ].data
        rhot  = cv['RHO'  ][-1, :]
        Surf  = cv['SLAT' ][-1, :]
        ne    = cv['NE'   ][-1, :]
        ni    = cv['NI'   ][-1, :]
        Te    = cv['TE'   ][-1, :]
        Ti    = cv['TI'   ][-1, :]
        Zef   = cv['ZEF'  ][-1, :]
        Cbs   = cv['CUBS' ][-1, :]
        FP    = cv['FP'   ][-1, :]
        CU    = cv['CU'   ][-1, :]
        MU    = cv['MU'   ][-1, :]
        Pbper = cv['PBPER'][-1, :]
        Pblon = cv['PBLON'][-1, :]

# fml/bpfit
        BPF = 1.6e-4*4.*np.pi*(ne*Te + ni*Ti + 0.5*Pbper + 0.5*Pblon) * (RTOR/(self.BTOR*rhot*np.abs(MU)))**2

# fml/prest
        Ptot = 1.6*(ne*Te + ni*Ti + 0.5*Pblon + 0.5*Pbper - nibm*Ti)

# Normalised Psi from ASTRA (irregular grid)
        FP_norm = ((FP - FP[0])/(FP[-1] - FP[0]))
        rhop = np.sqrt(FP_norm)

        index = np.argmin(np.abs(rhotN - (1. - float(self.Wi))))
        self.Wipol = 1. - rhop[index]
        self.dens_top = ne[index]
        self.Zeff = np.mean(Zef)

        nrho = len(rhot)
        nr_na1 = nrho/na1
        self.psiN = np.linspace(0, 1, nrho, endpoint=True)

        splp = UnivariateSpline(FP[:-2] - FP[0], -1.e-2*Ptot[:-2], s=0.000001)
        dP_n = splp(FP - FP[0], nu=1)*1.05
        spl = InterpolatedUnivariateSpline(FP_norm, dP_n)
        self.dPtot_psiN = spl(self.psiN)

        def err(c, x, y, t, k, w=None):
            diff = y - splev(x, (t, c, k))
            if w is None:
                diff = np.einsum('...i,...i', diff, diff)
            else:
                diff = np.dot(diff*diff, w)
            return np.abs(diff)

        def spline_neumann(x, y, k=3, s=0, w=None):
            t, c0, k = splrep(x, y, w, k=k, s=s)
            x0 = x[0] # point at which zero slope is required
            con = {'type': 'eq',
                   'fun': lambda c: splev(x0, (t, c, k), der=1),
                   #'jac': lambda c: splev(x0, (t, c, k), der=2) # doesn't help, dunno why
                   }
            opt = minimize(err, c0, (x, y, t, k, w), constraints=con)
            copt = opt.x
            return UnivariateSpline._from_tck((t, copt, k))

        if shottime == '30000_3.4':
            j_ax = min(2.2/2.5*self.BTOR*0.95/1.3, CU[0]*0.95)
        elif shottime[0] == '1':
            print('IT IS C-Mod!')
            j_ax = min(2.2/2.5*self.BTOR*1.65/RTOR*1.2, CU[0]*0.95)
        else:
            j_ax = min(2.2/2.5*self.BTOR*0.95, CU[0]*0.95)

        qind = np.argmin(abs(CU[int(100*nr_na1): int(300*nr_na1)] - j_ax*0.85)) + int(100*nr_na1)

        j_ptop_ind = np.argmin(CU[np.where(rhop < 0.97)])

        qind_2 = CU[qind: j_ptop_ind-141] - CU[qind+1: j_ptop_ind-140] #it was 70 instead of 140 but was too close to top
        if np.any(qind_2 < 0):
            qind_2 = np.max(np.where(qind_2<0)) + qind + 20
            qind = qind_2

        x = np.append(np.linspace(0, rhop[qind-50], 50), rhop[qind:j_ptop_ind])
        rampet = np.linspace(j_ax, j_ax*0.9, 50)
        y = np.append(rampet, CU[qind: j_ptop_ind])
        k = 2
        sp = spline_neumann(x, y, k, s=0.4) #n*std)
        X = rhop
        Y = sp(X)

        cu_n = np.append(Y[int(4*nr_na1): (j_ptop_ind-170)], CU  [(j_ptop_ind-100):])
        RU   = np.append(X[int(4*nr_na1): (j_ptop_ind-170)], rhop[(j_ptop_ind-100):])
        CU_f = UnivariateSpline(RU[:-1], cu_n[:-1], s=0.005*np.max(CU[j_ptop_ind:])**2, k=2)
        cun = CU_f(rhop)

        spu = InterpolatedUnivariateSpline(FP_norm, CU)
        Ip = (spu.integral(np.min(Surf), np.max(Surf)))

        spuf = InterpolatedUnivariateSpline(FP_norm, cun)
        Ipf = (spuf.integral(np.min(Surf), np.max(Surf)))
	
        if Ipf/Ip < 0.9:
            cun[0: qind] = cun[0]
            difff = cun[0] - cun[qind+50]
            for ji in range(50):
                cun[qind+ji] = cun[0] - difff/50.*ji
                cun[qind]    = (cun[qind-1]  + cun[qind+1])/2.
                cun[qind-1]  = (cun[qind-2]  + cun[qind])/2.
                cun[qind+50] = (cun[qind+49] + cun[qind+51])/2.
                cun[qind+49] = (cun[qind+48] + cun[qind+50])/2.

                spuf = InterpolatedUnivariateSpline(FP_norm, cun)
                Ipf = (spuf.integral(np.min(Surf), np.max(Surf)))

        spl = InterpolatedUnivariateSpline(FP_norm, Ptot)
        self.Ptot_psiN = spl(self.psiN)

        spl = InterpolatedUnivariateSpline(FP_norm, cun)
        self.CU_psiN = spl(self.psiN)

        spl = InterpolatedUnivariateSpline(FP_norm, Cbs)
        self.Cbs_psiN = spl(self.psiN)

        spl = InterpolatedUnivariateSpline(FP_norm, Te)
        self.Te_psiN = spl(self.psiN)

        spl = InterpolatedUnivariateSpline(FP_norm, Ti)
        self.Ti_psiN = spl(self.psiN)

        spl = InterpolatedUnivariateSpline(FP_norm, ne)
        self.ne_psiN = spl(self.psiN)

        spl = InterpolatedUnivariateSpline(FP_norm, BPF)
        self.BPF_psiN = spl(self.psiN)


    def dumpHelenaInput(self, astra_equ=None, fHelena=None):

        temp = self.shottime.split('_')
        shot = temp[0]
        time = temp[1]

        if fHelena is None:
            fHelena = 'helena/helena-%s-%s.dat' %(self.shottime, astra_equ)

        f_helena = open(fHelena, 'w')

        f_helena.write("# AUG Shot #%s @ %ss    Wped,tor=%s Wped,pol=%1.7f\n" %(shot, time, self.Wi, self.Wipol))

        np.savetxt(f_helena, (self.psiN**0.5, self.Te_psiN, self.Ti_psiN, self.ne_psiN, self.Ptot_psiN, self.dPtot_psiN, \
                self.BPF_psiN, self.CU_psiN, self.Cbs_psiN), \
                fmt='%1.8e', header="Rho_pol    Te[keV] Ti[keV]    ne[10^19/m^3]    p_tot[MJ/m^3]    dp_tot[Pa/Wb] Beta_p    j[A/m^2]    jBS[A/m^2]    Boundary (r) (z) dens_top[10^19/m^3]    mean(Zeff)  Ip[MA]  Btor[T]\n")

        for r in self.rs:
            f_helena.write("%.6E " %r)
        f_helena.write("\n")
        for z in self.zs:
            f_helena.write("%.6E " %z)
        f_helena.write("\n")
        f_helena.write("%.4E\n" %(self.dens_top))
        f_helena.write("%.4E\n" %(np.mean(self.Zeff)))
        f_helena.write("%.4E\n" %self.IPL)
        f_helena.write("%.4E\n" %self.BTOR)

        f_helena.close()
        print ('Stored file %s' %fHelena)
      

def test():

    exp = '30000_3.4'    
    for astra_equ in ('imep', 'imep2'):
        f_cdf = '/toks/work/git/a82/ncdf_out/%s%s.CDF' %(exp, astra_equ)
        a2h = astra2helena(f_cdf, shottime=exp)
        a2h.dumpHelenaInput(astra_equ=astra_equ)

    import matplotlib.pylab as plt
    plt.subplot(1, 1, 1, aspect='equal')
    plt.plot(a2h.rs, a2h.zs)
    plt.show()


if __name__ == '__main__':

#    test()
    parser = argparse.ArgumentParser(description='Write settings & run ASTRA')
    parser.add_argument('-m', '--equ', help='Model file', required=True)
    parser.add_argument('-v', '--exp', help='Exp file'  , required=True)

    args = parser.parse_args()

    f_cdf = '%s/ncdf_out/%s%s.CDF' %(awd, args.exp, args.equ)
    a2h = astra2helena(f_cdf, shottime=args.exp)
    a2h.dumpHelenaInput(astra_equ=args.equ)
