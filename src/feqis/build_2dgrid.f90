subroutine build_2dgrid(nrho, ntheta, Rb, Zb, X0, Y0, lambda2d, lambda2dp, psin_grid, &
    PSI, rtor, pressure, btor, ipol, iplasma, &
! Output
    XX, YY, rmin, G2, G3, areat, perim, volum, G1, GRADRO, &
    BMAXT, BMINT, BDB02, BDB0, B0DB2, FOFB, &
    slat, li3, betapol)

use pi_vars, only: GPI, GPI2
use interp_mod, only: qinterp, extrapolate

implicit none

integer, intent(in) :: ntheta, nrho
double precision, intent(in) :: X0, Y0, btor, iplasma, Rtor
double precision, intent(in), dimension(ntheta) :: Rb, Zb
double precision, intent(in), dimension(nrho) :: psin_grid, ipol, pressure
double precision, intent(in), dimension(nrho, ntheta) :: lambda2d, lambda2dp, PSI

double precision, intent(out) :: li3, betapol
double precision, intent(out), dimension(nrho) :: G1, G2, G3, &
    volum, areat, perim, slat, &
    FOFB, GRADRO, BMAXT, BMINT, BDB02, BDB0, B0DB2
double precision, intent(out), dimension(nrho, ntheta) :: XX, YY, rmin

integer :: jrho, jthe, jthe_l, k
double precision :: drdX, drdY, Mdet, dpsi, dthe, ipol_rmaj, z1, z2, z3, rho_interp
double precision, dimension(nrho) :: rhot, rhoa, &
    dPSIdV, dVa, daa, dum1, AMETR, ONEZ
double precision, dimension(ntheta) :: dXb0, dXb0_i, dl_arc, tar1, tar2
double precision, dimension(ntheta+1) :: thetap, thetap_i
double precision, dimension(nrho, ntheta) :: &
    lambda2dpi, X_i, Y_i, dXdr2, dYdr2, dXdh2, dYdh2, &
    rr2, r_i, Rmaj2, Jcbn2, gradr2, &
    gradPSIa, gradVa, dV2da, dA2da, B_pola, B_ABSa, B_Ta

do jthe=1, ntheta
    dXb0(jthe) = sqrt((rb(jthe) - X0)**2 + (zb(jthe) - Y0)**2)
    thetap(jthe) = ATAN2(zb(jthe) - Y0, rb(jthe) - X0)
    if (thetap(jthe) < 0) thetap(jthe) = thetap(jthe) + GPI2
enddo

if (thetap(1) > thetap(ntheta)) then !reorder
    jthe = 1
    do k=1, ntheta-1
        if (thetap(k+1) < thetap(k)) jthe = k + 1   ! j is the first teta above 0
    enddo 
    if (jthe > 1) thetap(1: jthe-1) = thetap(1: jthe-1) - GPI2
endif 
thetap(ntheta+1) = thetap(1) + GPI2

! Now find the intermediate theta grid
do jthe=1, ntheta
    thetap_i(jthe) = 0.5*(thetap(jthe) + thetap(jthe+1))
enddo
thetap_i(ntheta+1) = thetap_i(1) + GPI2

do jthe=1, ntheta-1
    lambda2dpi(:, jthe) = 0.5*(lambda2dp(:, jthe) + lambda2dp(:, jthe+1))
enddo
lambda2dpi(:, ntheta) = 0.5*(lambda2dp(:, ntheta) + lambda2dp(:, 1))
 
do jthe=1, ntheta-1
    dXb0_i(jthe) = 0.5*(dXb0(jthe) + dXb0(jthe+1))
enddo
dXb0_i(ntheta) = dXb0_i(1)

do jthe=1, ntheta
    do jrho=1, nrho
        rmin(jrho, jthe) = lambda2d  (jrho, jthe)*dXb0  (jthe)
        rr2 (jrho, jthe) = lambda2dp (jrho, jthe)*dXb0  (jthe)
        r_i (jrho, jthe) = lambda2dpi(jrho, jthe)*dXb0_i(jthe)
        XX   (jrho, jthe) = X0 + rmin(jrho, jthe)*cos(thetap  (jthe))
        YY   (jrho, jthe) = Y0 + rmin(jrho, jthe)*sin(thetap  (jthe))
        Rmaj2(jrho, jthe) = X0 + rr2 (jrho, jthe)*cos(thetap  (jthe))
        X_i  (jrho, jthe) = X0 + r_i (jrho, jthe)*cos(thetap_i(jthe))
        Y_i  (jrho, jthe) = Y0 + r_i (jrho, jthe)*sin(thetap_i(jthe))
    enddo
enddo
   
! Now computes jacobian
! J = det( {dX/dr dX/dtheta}  {dY/dr dY/dtheta} ) = dX/dr*dY/dtheta - dX/dtheta*dY/dr= J(r, theta)

do jrho=1, nrho-1
    dpsi = psin_grid(jrho+1) - psin_grid(jrho)
    do jthe=1, ntheta
        if (jthe == 1) then
            jthe_l = ntheta
        else
            jthe_l = jthe - 1
        endif
        dthe = thetap_i(jthe_l+1) - thetap_i(jthe_l) 
        dXdr2(jrho, jthe) = (XX(jrho+1, jthe) - XX (jrho, jthe  ))/dpsi 
        dYdr2(jrho, jthe) = (YY(jrho+1, jthe) - YY (jrho, jthe  ))/dpsi
        dXdh2(jrho, jthe) = (X_i( jrho, jthe) - X_i(jrho, jthe_l))/dthe
        dYdh2(jrho, jthe) = (Y_i( jrho, jthe) - Y_i(jrho, jthe_l))/dthe
    enddo
enddo

Jcbn2 = dXdr2*dYdh2 - dXdh2*dYdr2   ! i+1/2, j 

do jthe=1, ntheta
    do jrho=1, nrho-1
        Mdet =  dXdr2(jrho, jthe)*dYdh2(jrho, jthe) - dXdh2(jrho, jthe)*dYdr2(jrho, jthe)
        drdX =  dYdh2(jrho, jthe)/Mdet
        drdY = -dXdh2(jrho, jthe)/Mdet
        gradr2(jrho, jthe) = drdX**2 + drdY**2 
    enddo
enddo

! build_2dgrid2

do jrho=1, nrho
    rhoa(jrho) = (jrho - 1.)/(nrho - 1.) ! full grid
enddo
rhot = rhoa + 0.5/(nrho - 1.)  ! half grid

da2da = 0.
dv2da = 0.

! half grid points around area and volume (last point doesn exist)
jthe = 1
do jrho=1, nrho-1
    da2da(jrho, jthe) = Jcbn2(jrho, jthe)*(psin_grid(jrho+1) - psin_grid(jrho)) * &
        (thetap_i(ntheta+1) - thetap_i(ntheta)) 
    dv2da(jrho, jthe) = GPI2*Rmaj2(jrho, jthe)*da2da(jrho, jthe)
enddo
do jthe=2, ntheta
    do jrho=1, nrho-1
        da2da(jrho, jthe) = Jcbn2(jrho, jthe)*(psin_grid(jrho+1) - psin_grid(jrho)) * &
            (thetap_i(jthe) - thetap_i(jthe-1)) 
        dv2da(jrho, jthe) = GPI2*Rmaj2(jrho, jthe)*da2da(jrho, jthe)
    enddo
enddo

dVa = sum(dv2da, 2)
daa = sum(da2da, 2)

slat  = 0.
areat = 0.
volum = 0.
do jrho=2, nrho
    volum(jrho) = volum(jrho-1) + dva(jrho-1)
    areat(jrho) = areat(jrho-1) + daa(jrho-1)
enddo

dl_arc = 0.0
do jrho=1, nrho
    dl_arc = 0.0
    do jthe=2, ntheta 
        dl_arc(jthe) = rmin(jrho, jthe)*(thetap_i(jthe) - thetap_i(jthe-1)) !on the full grid
    enddo
    dl_arc(1) = rmin(jrho, 1)*(thetap_i(ntheta+1) - thetap_i(ntheta))
    perim(jrho) = sum(dl_arc)
enddo

!Compute Bpol, gradPSI, gradV

do jrho=1, nrho-1
    dPSIdv(jrho) = (PSI(jrho+1, 1) - PSI(jrho, 1))/(volum(jrho+1) - volum(jrho)) !on the psigp grid
enddo

! there are on half grid
do jrho=1, nrho-1
    onez(jrho) = 0.5*(IPOL(jrho) + IPOL(jrho+1))
enddo
do jrho=1, nrho-1
    z1 = dVa(jrho)/(psin_grid(jrho+1) - psin_grid(jrho))
    z3 = 0.5*(IPOL(jrho) + IPOL(jrho+1))
    do jthe=1, ntheta 
        z2 = (PSI(jrho+1, jthe) - PSI(jrho, jthe))/(psin_grid(jrho+1) - psin_grid(jrho))
        B_pola(jrho, jthe) = z2/(GPI2*Rmaj2(jrho, jthe))*sqrt(gradr2(jrho, jthe))
        gradPSIa(jrho, jthe) = z2*sqrt(gradr2(jrho, jthe))
        gradVa  (jrho, jthe) = z1*sqrt(gradr2(jrho, jthe))
        B_Ta(jrho, jthe) = z3/Rmaj2(jrho, jthe)
    enddo
enddo   
B_absa = sqrt(B_Ta**2 + B_pola**2)

li3 = 2.*sum(B_pola**2 * dV2da)/rtor/(0.4*GPI*iplasma)**2
do jrho=1, nrho-1
    onez(jrho) = 0.5*(pressure(jrho) + pressure(jrho+1))
enddo
betapol = 0.4*GPI2*1.e-6*sum(onez*dva)/sum(B_pola**2 * dV2da)
slat = 0.
do jrho=2, nrho
    do jthe=1, ntheta-1 
        slat(jrho) = slat(jrho) + XX(jrho, jthe)*rmin(jrho, jthe)*(thetap_i(jthe+1) - thetap_i(jthe))
    enddo
    slat(jrho) = slat(jrho) + XX(jrho, ntheta)*rmin(jrho, ntheta)*(thetap_i(1) - thetap_i(ntheta) + GPI2)
enddo
slat = slat*GPI2   !full grid

do jrho=1, nrho
    ametr(jrho) = 0.5*(maxval(XX(jrho, 1: ntheta)) - minval(XX(jrho, 1: ntheta)))
enddo

!Cycle over positions  ! half grid
do jrho=1, nrho-1
    tar2 = dV2da(jrho, 1: ntheta)/dVa(jrho)
    tar1 = 1./(Rmaj2(jrho, 1: ntheta)**2)
    z1 = sum(tar1*tar2)
    G3(jrho) = z1
    z1 = sum(tar2)
    ONEZ(jrho) = z1
    tar1 = (gradVa(jrho, 1: ntheta)/Rmaj2(jrho, 1: ntheta))**2
    z1 = sum(tar1*tar2)
    G2(jrho) = z1

    tar1 = gradVa(jrho, 1: ntheta)**2
    z1 = sum(tar1*tar2)
    G1(jrho) = z1

    tar1 = gradVa(jrho, 1: ntheta)
    z1 = sum(tar1*tar2)
    GRADRO(jrho) = z1

    tar1 = B_ABSa(jrho, 1: ntheta)**2
    z1 = sum(tar1*tar2)
    BDB02(jrho) = z1

    tar1 = B_ABSa(jrho, 1: ntheta)
    z1 = sum(tar1*tar2)
    BDB0(jrho) = z1

    tar1 = 1./(B_ABSa(jrho, 1: ntheta)**2)
    z1 = sum(tar1*tar2)

    B0DB2(jrho) = z1
    BMAXT(jrho) = maxval(B_ABSa(jrho, : ))
    BMINT(jrho) = minval(B_ABSa(jrho, : ))

    tar2 = B_ABSa(jrho, 1: ntheta)/BMAXT(jrho)
    tar1 = ((btor/B_ABSa(jrho, 1: ntheta))**2) * &
        (  1. - (sqrt(1. - tar2)) * (1. + 0.5*tar2)  )
    tar2 = dV2da(jrho, 1: ntheta)/dVa(jrho)
    z1 = sum(tar1*tar2)
    FOFB(jrho) = z1

enddo

rho_interp = rhoa(nrho)

G1(nrho)     = EXTRAPOLATE(rho_interp, nrho-1, nrho-2, nrho-3, nrho, rhot, G1)
G2(nrho)     = EXTRAPOLATE(rho_interp, nrho-1, nrho-2, nrho-3, nrho, rhot, G2)
G3(nrho)     = EXTRAPOLATE(rho_interp, nrho-1, nrho-2, nrho-3, nrho, rhot, G3)
GRADRO(nrho) = EXTRAPOLATE(rho_interp, nrho-1, nrho-2, nrho-3, nrho, rhot, GRADRO)
FOFB(nrho)   = EXTRAPOLATE(rho_interp, nrho-1, nrho-2, nrho-3, nrho, rhot, FOFB)
BMAXT(nrho)  = EXTRAPOLATE(rho_interp, nrho-1, nrho-2, nrho-3, nrho, rhot, BMAXT)
BMINT(nrho)  = EXTRAPOLATE(rho_interp, nrho-1, nrho-2, nrho-3, nrho, rhot, BMINT)
BDB02(nrho)  = EXTRAPOLATE(rho_interp, nrho-1, nrho-2, nrho-3, nrho, rhot, BDB02)
BDB0(nrho)   = EXTRAPOLATE(rho_interp, nrho-1, nrho-2, nrho-3, nrho, rhot, BDB0)
B0DB2(nrho)  = EXTRAPOLATE(rho_interp, nrho-1, nrho-2, nrho-3, nrho, rhot, B0DB2)

call qinterp(rhot(1: nrho-1), G1(1: nrho-1), nrho-1, rhoa(2: nrho-1), dum1(2: nrho-1), nrho-2)
G1(2: nrho-1) = dum1(2: nrho-1)
call qinterp(rhot(1: nrho-1), g2(1: nrho-1), nrho-1, rhoa(2: nrho-1), dum1(2: nrho-1), nrho-2)
g2(2: nrho-1) = dum1(2: nrho-1)
call qinterp(rhot(1: nrho-1), g3(1: nrho-1), nrho-1, rhoa(2: nrho-1), dum1(2: nrho-1), nrho-2)
g3(2: nrho-1) = dum1(2: nrho-1)
call qinterp(rhot(1: nrho-1), gradro(1: nrho-1), nrho-1, rhoa(2: nrho-1), dum1(2: nrho-1), nrho-2)
gradro(2: nrho-1) = dum1(2: nrho-1)
call qinterp(rhot(1: nrho-1), fofb(1: nrho-1), nrho-1, rhoa(2: nrho-1), dum1(2: nrho-1), nrho-2)
fofb(2: nrho-1) = dum1(2: nrho-1)
call qinterp(rhot(1: nrho-1), bmaxt(1: nrho-1), nrho-1, rhoa(2: nrho-1), dum1(2: nrho-1), nrho-2)
bmaxt(2: nrho-1) = dum1(2: nrho-1)
call qinterp(rhot(1: nrho-1), bmint(1: nrho-1), nrho-1, rhoa(2: nrho-1), dum1(2: nrho-1), nrho-2)
bmint(2: nrho-1) = dum1(2: nrho-1)
call qinterp(rhot(1: nrho-1), bdb02(1: nrho-1), nrho-1, rhoa(2: nrho-1), dum1(2: nrho-1), nrho-2)
bdb02(2: nrho-1) = dum1(2: nrho-1)
call qinterp(rhot(1: nrho-1), bdb0(1: nrho-1), nrho-1, rhoa(2: nrho-1), dum1(2: nrho-1), nrho-2)
bdb0(2: nrho-1) = dum1(2: nrho-1)
call qinterp(rhot(1: nrho-1), b0db2(1: nrho-1), nrho-1, rhoa(2: nrho-1), dum1(2: nrho-1), nrho-2)
b0db2(2: nrho-1) = dum1(2: nrho-1)

G1(1) = 0.0
G3(1) = 1./(XX(1, 1)**2)
G2(1) = 0.0
GRADRO(1) = 0.0
ipol_rmaj = IPOL(1)/XX(1, 1)
BDB02(1) = ipol_rmaj**2
BDB0(1)  = ipol_rmaj
B0DB2(1) = 1./ipol_rmaj**2
FOFB(1)  = (btor/ipol_rmaj)**2
BMAXT(1) = ipol_rmaj
BMINT(1) = ipol_rmaj

return
end subroutine build_2dgrid
