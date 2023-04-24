subroutine build_2dgrid2(Nr, Nt, psig, &
    Rmaj, Rmaj2, r, thetap_i, gradr2, Jcbn2, &
    PSI, rtor, pressure, btor, ipol, iplasma, &
! Output
    G2, G3, areat, perim, volum, G1, GRADRO, &
    BMAXT, BMINT, BDB02, BDB0, B0DB2, &
    FOFB, slat, li3, betapol)

implicit none

double precision, parameter :: GPI=3.141592653589793, GPI2=2.*GPI

integer, intent(in) :: Nt, Nr
double precision, intent(in) :: btor, iplasma, Rtor
double precision, intent(in), dimension(Nr) :: psig, ipol, pressure
double precision, intent(in), dimension(Nt+1) :: thetap_i
double precision, intent(in), dimension(Nr, Nt) :: Rmaj, Rmaj2, &
    Jcbn2, r, gradr2, PSI
double precision, intent(out) :: li3, betapol
double precision, intent(out), dimension(Nr) :: G1, G2, G3, &
    volum, areat, perim, slat, &
    FOFB, GRADRO, BMAXT, BMINT, BDB02, BDB0, B0DB2

integer :: jr, jt, i, j
double precision :: ipol_rmaj, z1, z2, rho_interp
double precision, dimension(Nr) :: rhot, rhoa, &
    dPSIdV, dVa, daa, dum1, AMETR, ONEZ
double precision, dimension(Nt) :: dl_arc, tar1, tar2
double precision, dimension(Nr, Nt) :: gradPSIa, gradVa, &
    dV2da, dA2da, B_pola, B_ABSa, B_Ta
double precision, external :: EXTRAPOLATE

do i=1, Nr
    rhoa(i) = (i - 1.)/(Nr - 1.) ! full grid
enddo
rhot = rhoa + 0.5/(Nr - 1.)  ! half grid

da2da = 0.
dv2da = 0.

! half grid points around area and volume (last point doesn exist)
jt = 1
do jr=1, Nr-1
    da2da(jr, jt) = Jcbn2(jr, jt)*(psig(jr+1) - psig(jr)) * &
        (thetap_i(Nt+1) - thetap_i(Nt)) 
    dv2da(jr, jt) = GPI2*Rmaj2(jr, jt)*da2da(jr, jt)
enddo
do jt=2, nt
    do jr=1, Nr-1
        da2da(jr, jt) = Jcbn2(jr, jt)*(psig(jr+1) - psig(jr)) * &
            (thetap_i(jt) - thetap_i(jt-1)) 
        dv2da(jr, jt) = GPI2*Rmaj2(jr, jt)*da2da(jr, jt)
    enddo
enddo

dVa = sum(dv2da, 2)
daa = sum(da2da, 2)

slat  = 0.
areat = 0.
volum = 0.
do i=2, Nr
    volum(i) = volum(i-1) + dva(i-1)
    areat(i) = areat(i-1) + daa(i-1)
enddo

dl_arc = 0.0
do j=1, Nr
    dl_arc = 0.0   
    do i=2, Nt 
        dl_arc(i) = r(j, i)*(thetap_i(i) - thetap_i(i-1)) !on the full grid
    enddo
    dl_arc(1) = r(j, 1)*(thetap_i(nt+1) - thetap_i(nt))
    perim(j) = sum(dl_arc)
enddo

!Compute Bpol, gradPSI, gradV

do i=1, Nr-1
    dPSIdv(i) = (PSI(i+1, 1) - PSI(i, 1))/(volum(i+1) - volum(i)) !on the psigp grid
enddo

! there are on half grid
do i=1, Nr-1
    onez(i) = 0.5*(IPOL(i) + IPOL(i+1))
enddo
do i=1, Nr-1
    do j=1, Nt 
        z1 = dVa(i)/(psig(i+1) - psig(i))
        z2 = (PSI(i+1, j) - PSI(i, j))/(psig(i+1) - psig(i))
        B_pola(i, j) = z2/(GPI2*Rmaj2(i, j))*sqrt(gradr2(i, j))
        gradPSIa(i, j) = z2*sqrt(gradr2(i, j))
        gradVa(  i, j) = z1*sqrt(gradr2(i, j))
        B_Ta(i, j) = onez(i)/Rmaj2(i, j)
    enddo
enddo   
B_absa = sqrt(B_Ta**2 + B_pola**2)

li3 = 2.*sum(B_pola**2 * dV2da)/rtor/(0.4*GPI*iplasma)**2
do i=1, Nr-1
    onez(i) = 0.5*(pressure(i) + pressure(i+1))
enddo
betapol = 0.4*GPI2*1.e-6*sum(onez*dva)/sum(B_pola**2 * dV2da)
slat = 0.
do i=2, Nr
    do j=1, Nt-1 
        slat(i) = slat(i) + Rmaj(i, j)*r(i, j)*(thetap_i(j+1) - thetap_i(j))
    enddo
    slat(i) = slat(i) + Rmaj(i, Nt)*r(i, Nt)*(thetap_i(1) - thetap_i(Nt) + GPI2)
enddo
slat = slat*GPI2   !full grid

do i=1, Nr
    ametr(i) = 0.5*(maxval(Rmaj(i, 1: Nt)) - minval(Rmaj(i, 1: Nt)))      ! full grid
enddo

!Cycle over positions  ! half grid
do i=1, Nr-1
    tar2 = dV2da(i, 1: Nt)/dVa(i)
    tar1 = 1./(Rmaj2(i, 1: Nt)**2)
    z1 = sum(tar1*tar2)
    G3(i) = z1
    z1 = sum(tar2)
    ONEZ(i) = z1
    tar1 = (gradVa(i, 1: Nt)/Rmaj2(i, 1: Nt))**2
    z1 = sum(tar1*tar2)
    G2(i) = z1

    tar1 = gradVa(i, 1: Nt)**2
    z1 = sum(tar1*tar2)
    G1(i) = z1

    tar1 = gradVa(i, 1: Nt)
    z1 = sum(tar1*tar2)
    GRADRO(i) = z1

    tar1 = B_ABSa(i, 1: Nt)**2
    z1 = sum(tar1*tar2)
    BDB02(i) = z1

    tar1 = B_ABSa(i, 1: Nt)
    z1 = sum(tar1*tar2)
    BDB0(i) = z1

    tar1 = 1./(B_ABSa(i, 1: Nt)**2)
    z1 = sum(tar1*tar2)

    B0DB2(i) = z1
    BMAXT(i) = maxval(B_ABSa(i, : ))
    BMINT(i) = minval(B_ABSa(i, : ))

    tar2 = B_ABSa(i, 1: Nt)/BMAXT(i)
    tar1 = ((btor/B_ABSa(i, 1: Nt))**2) * &
        (  1. - (sqrt(1. - tar2)) * (1. + 0.5*tar2)  )
    tar2 = dV2da(i, 1: Nt)/dVa(i)
    z1 = sum(tar1*tar2)
    FOFB(i) = z1

enddo

rho_interp = rhoa(Nr)

G1(Nr)     = EXTRAPOLATE(rho_interp, Nr-1, Nr-2, Nr-3, Nr, rhot, G1)
G2(Nr)     = EXTRAPOLATE(rho_interp, Nr-1, Nr-2, Nr-3, Nr, rhot, G2)
G3(Nr)     = EXTRAPOLATE(rho_interp, Nr-1, Nr-2, Nr-3, Nr, rhot, G3)
GRADRO(Nr) = EXTRAPOLATE(rho_interp, Nr-1, Nr-2, Nr-3, Nr, rhot, GRADRO)
FOFB(Nr)   = EXTRAPOLATE(rho_interp, Nr-1, Nr-2, Nr-3, Nr, rhot, FOFB)
BMAXT(Nr)  = EXTRAPOLATE(rho_interp, Nr-1, Nr-2, Nr-3, Nr, rhot, BMAXT)
BMINT(Nr)  = EXTRAPOLATE(rho_interp, Nr-1, Nr-2, Nr-3, Nr, rhot, BMINT)
BDB02(Nr)  = EXTRAPOLATE(rho_interp, Nr-1, Nr-2, Nr-3, Nr, rhot, BDB02)
BDB0(Nr)   = EXTRAPOLATE(rho_interp, Nr-1, Nr-2, Nr-3, Nr, rhot, BDB0)
B0DB2(Nr)  = EXTRAPOLATE(rho_interp, Nr-1, Nr-2, Nr-3, Nr, rhot, B0DB2)

call qinterp_feqis(rhot(1: Nr-1), G1(1: Nr-1), Nr-1, rhoa(2: Nr-1), dum1(2: Nr-1), Nr-2)
G1(2: Nr-1) = dum1(2: Nr-1)
call qinterp_feqis(rhot(1: Nr-1), g2(1: Nr-1), Nr-1, rhoa(2: Nr-1), dum1(2: Nr-1), Nr-2)
g2(2: Nr-1) = dum1(2: Nr-1)
call qinterp_feqis(rhot(1: Nr-1), g3(1: Nr-1), Nr-1, rhoa(2: Nr-1), dum1(2: Nr-1), Nr-2)
g3(2: Nr-1) = dum1(2: Nr-1)
call qinterp_feqis(rhot(1: Nr-1), gradro(1: Nr-1), Nr-1, rhoa(2: Nr-1), dum1(2: Nr-1), Nr-2)
gradro(2: Nr-1) = dum1(2: Nr-1)
call qinterp_feqis(rhot(1: Nr-1), fofb(1: Nr-1), Nr-1, rhoa(2: Nr-1), dum1(2: Nr-1), Nr-2)
fofb(2: Nr-1) = dum1(2: Nr-1)
call qinterp_feqis(rhot(1: Nr-1), bmaxt(1: Nr-1), Nr-1, rhoa(2: Nr-1), dum1(2: Nr-1), Nr-2)
bmaxt(2: Nr-1) = dum1(2: Nr-1)
call qinterp_feqis(rhot(1: Nr-1), bmint(1: Nr-1), Nr-1, rhoa(2: Nr-1), dum1(2: Nr-1), Nr-2)
bmint(2: Nr-1) = dum1(2: Nr-1)
call qinterp_feqis(rhot(1: Nr-1), bdb02(1: Nr-1), Nr-1, rhoa(2: Nr-1), dum1(2: Nr-1), Nr-2)
bdb02(2: Nr-1) = dum1(2: Nr-1)
call qinterp_feqis(rhot(1: Nr-1), bdb0(1: Nr-1), Nr-1, rhoa(2: Nr-1), dum1(2: Nr-1), Nr-2)
bdb0(2: Nr-1) = dum1(2: Nr-1)
call qinterp_feqis(rhot(1: Nr-1), b0db2(1: Nr-1), Nr-1, rhoa(2: Nr-1), dum1(2: Nr-1), Nr-2)
b0db2(2: Nr-1) = dum1(2: Nr-1)

G1(1) = 0.0
G3(1) = 1./(Rmaj(1, 1)**2)
G2(1) = 0.0
GRADRO(1) = 0.0
ipol_rmaj = IPOL(1)/Rmaj(1, 1)
BDB02(1) = ipol_rmaj**2
BDB0(1)  = ipol_rmaj
B0DB2(1) = 1./ipol_rmaj**2
FOFB(1)  = (btor/ipol_rmaj)**2
BMAXT(1) = ipol_rmaj
BMINT(1) = ipol_rmaj

return
end subroutine build_2dgrid2

!-------------------------------------------------------------------
double precision function EXTRAPOLATE(x_interp, j1, j2, j3, narr, rho, y_in)

implicit none

integer, intent(in) :: j1, j2, j3, narr
double precision, intent(in) :: x_interp
double precision, intent(in), dimension(narr) :: rho, y_in

double precision :: x1, x2, x3, f1, f2, f3, dfdx, d2fdx

f1 = y_in(j1)
f2 = y_in(j2)
f3 = y_in(j3)
x1 = rho(j1)
x2 = rho(j2)
x3 = rho(j3)
dfdx = (f3 - f1)/(x3 - x1)
d2fdx = ((f3 - f2)/(x3 - x2) - (f2 - f1)/(x2 - x1))/(x3 - x1)

EXTRAPOLATE = f2 + (x_interp - x2)*(dfdx + d2fdx*(x_interp - x2))
      
return
end function EXTRAPOLATE
