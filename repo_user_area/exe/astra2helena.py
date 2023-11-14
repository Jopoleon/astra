#!/usr/bin/env python

import numpy as np
from scipy.interpolate import splev, splrep
import os, logging, argparse
from scipy.interpolate import InterpolatedUnivariateSpline, UnivariateSpline
from scipy.optimize import minimize
from scipy.io import netcdf_file

fmt = logging.Formatter('%(asctime)s | %(name)s | %(levelname)s: %(message)s', '%H:%M:%S')
logger = logging.getLogger('a2helena')
logger.setLevel(logging.INFO)
if len(logger.handlers) == 0:
    hnd = logging.StreamHandler()
    hnd.setFormatter(fmt)
    logger.addHandler(hnd)

awd = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

helena_cluster = \
'''#!/bin/sh
#SBATCH -J  HELENA_IPED     #Job name
#SBATCH -D %s            #Initial working directory
#SBATCH --partition=s.tok     #Queue/Partition
#SBATCH --qos=s.tok.short  #Quality of Service
#SBATCH --mem=4000
##
#SBATCH --mail-type=none       #Send mail, e.g. for begin/end/fail/none
#SBATCH --mail-user=%s@ipp.mpg.de  #Mail address
   
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
# For pinning threads correctly:
export OMP_PLACES=cores 
export MKL_NUM_THREADS=1
# Run the program:
srun %s/run_helena.sh'
'''

helena_xml = \
'''<?xml version="1.0" encoding="undecided"?>
version="1.0"?>

<?xml-stylesheet type="text/xsl" href="./input_helena.xsl"
charset="ISO-8859-1"?>


<parameters>

<!-- profile_parameters -->

<profile_parameters>
    <input_type> p' and j_tor </input_type>
 <radial_coordinate> psi </radial_coordinate>
<current_averaging> theta </current_averaging>
<hbt> F </hbt>
</profile_parameters>

<!-- shape_parameters -->

<shape_parameters>
 <isol> 0 </isol>
<ias> 1 </ias>
<mfm> 256 </mfm>
<imesh> 2 </imesh>
<n_acc_points> 1 </n_acc_points>
<s_acc> 1.0  </s_acc>
<sig> 0.02  </sig>
<weights> 1.0  </weights>
<equidistant> 0.2 </equidistant>
</shape_parameters>

<!-- global_parameters -->

<global_parameters>
<match> Ip </match>
<Ip> %s </Ip>
<bvac> %s </bvac>
<cpsurfin> %s  </cpsurfin>
 </global_parameters>

<!-- numerical_parameters -->

<numerical_parameters>
<nr> %d </nr>
<np> %d </np>
<nrmap> %d </nrmap>
<npmap> %d </npmap>
<nchi>  %d </nchi>
<niter> 80 </niter>
<nmesh> 50 </nmesh>
 <nouter> 20 </nouter>
<errcur> 1.0E-5 </errcur>
</numerical_parameters>

<!-- diagnostics_parameters -->

<diagnostics_parameters>
<verbosity> 4 </verbosity>
<output> full </output>
<cpo_output> False </cpo_output>
<diagnostics_on> T </diagnostics_on>
<standard_output> T </standard_output>
</diagnostics_parameters>

</parameters>
'''

# str(Ip), str(Bt), str(psibnd), nr, np, nrmap, npmap, nchi

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
            logger.warning('IT IS C-Mod!')
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
        logger.info('Stored file %s', fHelena)


    def toHelenaBoundary(self):

        rgeo = 0.5*(np.min(self.rs) + np.max(self.rs))
        zgeo = 0.5*(np.min(self.zs) + np.max(self.zs))

        theta = np.arctan2(self.zs - zgeo, self.rs - rgeo)
        whneg = np.where(theta < 0)
        theta[whneg] += 2.*np.pi
        rad = np.sqrt((self.rs - rgeo)**2 + (self.zs - zgeo)**2)
        sorting = np.argsort(theta)
        rad   = rad  [sorting]
        theta = theta[sorting]
        nnew = 256
        newtheta = np.linspace(0., 2.*np.pi, num=nnew, endpoint=False)
        newrs = np.interp(newtheta, theta, rad)
        self.rsHelena = newrs * np.cos(newtheta) + rgeo
        self.zsHelena = newrs * np.sin(newtheta) + zgeo


    def writeHelenaProfiles(self, dir_out='helena'):
        '''Dump HELENA input profiles'''

        f_out = {}
        profs = ('p', 'j_tor', 'te', 'ne', 'dp')
        for lbl in profs:
            f_out[lbl] = '%s/%s.in' %(dir_out, lbl)

# Flatten gradients in core region
        pp = self.dPtot_psiN
        whcore    = np.squeeze(np.where(self.PsiN < 0.1))[::-1]
        whoutcore =  np.where((self.PsiN >= 0.1) & (self.PsiN < 0.2))
        gradcore  = np.mean(np.diff(pp[whoutcore]))       
        for ind in whcore:
            pp[ind] = pp[ind+1] - gradcore
            
        np.savetxt(f_out['p']    , np.transpose(1.e3 *self.Ptot_psiN), fmt='%15.8e')
        np.savetxt(f_out['j_tor'], np.transpose(1.e6 *self.CU_psiN  ), fmt='%15.8e')
        np.savetxt(f_out['te']   , np.transpose(1.e3 *self.Te_psiN  ), fmt='%15.8e')
        np.savetxt(f_out['ne']   , np.transpose(1.e19*self.ne_psiN  ), fmt='%15.8e')
        np.savetxt(f_out['dp']   , np.transpose(1.e5 *pp            ), fmt='%15.8e')
        for lbl in profs:
            logger.info('Stored %s', f_out[lbl])

    def writeHelenaSbatch(self, dir_out='helena'):
        '''Dump HELENA SLURM batch script'''

        user = (os.environ['USER'])

        f_sbatch = dir_out + '/helena_cluster.sh'
        with open(f_sbatch, 'w') as fsb:
            txt = helena_cluster %(dir_out, user, dir_out)
            fsb.write(txt)
        logger.info('Stored %s', f_sbatch)

    def writeHelenaXML(self, dir_out='helena', highres=False):
        '''Dump HELENA xml input'''

        psibnd = -1.3
        if highres:
            nr    = 501
            nrmap = 1001
            nnp   = 257
            npmap = 257
            nchi  = 1024
        else:
            nr    = 301
            nrmap = 1001
            nnp   = 257
            npmap = 257
            nchi  = 1024
        f_xml = dir_out + '/helena.xml'
        with open(f_xml, 'w') as fxml:
            txt = helena_xml %(str(self.IPL), str(self.BTOR), str(psibnd), nr, nnp, nrmap, npmap, nchi)
            fxml.write(txt)
        logger.info('Stored %s', f_xml)

    def writeHelenaBoundary(self, dir_out='helena'):
        '''Dump file with separatrix R,z'''

        if not hasattr(self, 'rsHelena'):
            self.toHelenaBoundary()
        f_bnd = dir_out + '/plasma_boundary.in'
        np.savetxt(f_bnd, np.c_[self.rsHelena, self.zsHelena], fmt='%15.8e')
        logger.info('Stored %s', f_bnd)

    def writeHelenaInput(self, dir_out='helena'):

        self.writeHelenaBoundary(dir_out=dir_out)
        self.writeHelenaProfiles(dir_out=dir_out)
        self.writeHelenaSbatch(dir_out=dir_out)
        self.writeHelenaXML(dir_out=dir_out)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Write settings & run ASTRA')
    parser.add_argument('-m', '--equ', help='Model file', required=True)
    parser.add_argument('-v', '--exp', help='Exp file'  , required=True)

    args = parser.parse_args()

    f_cdf = '%s/ncdf_out/%s%s.CDF' %(awd, args.exp, args.equ)
    a2h = astra2helena(f_cdf, shottime=args.exp)
#    a2h.dumpHelenaInput(astra_equ=args.equ)
    helena_dir = '/toks/work/git/a8/helena'
    a2h.writeHelenaInput(dir_out=helena_dir)
    cmd = 'sbatch %s/helena_cluster.sh' %helena_dir
    logger.info('Executing %s', cmd)
#    os.system(cmd)
