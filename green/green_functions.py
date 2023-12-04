import numpy as np
from scipy.special import ellipk, ellipe
from scipy.constants import mu_0

cPsi  = 2.*np.pi
gs_float = np.float64

def greenFunction(Rloc, Zloc, Rcoil, Zcoil):

    sum_sq = (Rloc + Rcoil)**2 + (Zloc - Zcoil)**2
    rsum4 = 4.*Rloc*Rcoil/sum_sq
#    pre_fac = np.sqrt(sum_sq)/(4.*np.pi) * cPsi *mu_0 # Lackner code def
    pre_fac = np.sqrt(sum_sq)/4              # Consistent with ASTRA matrices
    return pre_fac*((2. - rsum4)*ellipk(rsum4) - 2*ellipe(rsum4))


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
    dr = Rgrid[1] - Rgrid[0]
    dz = Zgrid[1] - Zgrid[0]
    llp = 2*nr2 + 2*nz2

    der = np.array(nr2*[dr] + nz2*[dz] + nr2*[dr] + nz2*[dz], dtype=gs_float)

    rnd1 = Rgrid[1:-1]
    rnd2 = Rgrid[ 0] + np.zeros(nz2, dtype=gs_float)
    rnd3 = Rgrid[-1] + np.zeros(nz2, dtype=gs_float)
    znd1 = Zgrid[ 0] + np.zeros(nr2, dtype=gs_float)
    znd2 = Zgrid[-1] + np.zeros(nr2, dtype=gs_float)
    znd3 = Zgrid[1:-1]
    rnd = np.hstack((rnd1, rnd3, rnd1, rnd2))
    znd = np.hstack((znd1, znd3, znd2, znd3))

    gfm = np.zeros((llp, llp), dtype=gs_float)
    for j in range(llp-1):
        gfm[j, j+1:] = greenFunction(rnd[j], znd[j], rnd[j+1:], znd[j+1:])
    gfm += gfm.T

    tt = der/(rnd*4)
    np.fill_diagonal(gfm, -tt*(np.log(tt**2) - 4.*np.log(2.) + 2.)*rnd**2/der)

    return gfm


def identity(R, dr, dz, ntype):

    (ind1, ) = np.where(ntype == 1)
    (ind2, ) = np.where(ntype == 2)
    gfi = np.zeros_like(R)
    if len(ind2) > 0:
        gfi[ind2] = R[ind2]*(np.log(8.*R[ind2]/(0.2236*(dr[ind2] + dz[ind2]))) - 2.)
    if len(ind1) > 0:
        dl  = np.hypot(dr, dz)
        for i in (-1, 0, 1):
            gfi[ind1] += (R[ind1] + dr[ind1]*i/3.) * (np.log( 8.*(R[ind1] + dr[ind1]*i/3.)/(dl[ind1]/3.) ) - 0.5)
            for j in range(-1, i):
                gfi[ind1] += 4.*greenFunction(R[ind1] + dr[ind1]*i/3., dz[ind1]*i/3., R[ind1] + dr[ind1]*j/3., dz[ind1]*j/3.)
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
