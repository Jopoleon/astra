subroutine build_2dgrid2(Nr, Nt, psig, &
    Rmaj, Rmaj2, r, thetap_i, gradr2, Jcbn2, &
    PSI, rtor, pressure, btor, ipol, iplasma, G2, G3, &
    areat, perim, volum, G1, GRADRO, &
    BMAXT, BMINT, BDB02, BDB0, B0DB2, &
    FOFB, slat, li3, betapol)

use pi_grec_vars, only: GPI, GPI2
use debugger_ef, only: markloc_ef, debug

implicit none

integer, intent(in) :: Nt, Nr
double precision, intent(in) :: btor, iplasma
double precision, intent(in), dimension(Nr) :: psig
double precision, intent(in), dimension(Nt+1) :: thetap_i
double precision, intent(in), dimension(Nr, Nt) :: Rmaj, Rmaj2, &
    Jcbn2, r, gradr2

integer :: jr, jt, i, j, k, ip0, ip1, ip2, ip3

double precision :: Rtor, li3, betapol, ipol_rmaj, z1, z2, t1, t2, t3, t4
double precision, dimension(Nr) :: rhot, rhoa, areat, perim, slat, FOFB, ipol, &
    dPSIdV, dVa, daa, dum1, pressure, G1, G2, G3, volum, GRADRO, &
    BMAXT, BMINT, BDB02, BDB0, B0DB2, AMETR, ONEZ
double precision, dimension(Nt) :: dl_arc, tar1, tar2
double precision, dimension(Nr, Nt) :: PSI, gradPSI, gradV, gradPSIa, gradVa, &
    dV2da, dA2da, B_pol, B_R, B_Z, B_pola, B_ABSa, B_Ta

call markloc_ef('build_2dgrid2', debug_lev=debug)

do i=1, Nr
    rhoa(i) =(i -1.      )/(Nr-1.) !full grid
    rhot(i) =(i -1. + 0.5)/(Nr-1.) !half grid
enddo

da2da = 0.
dv2da = 0.

! half grid points around area and volume (last point doesn exist)
jt=1
do jr=1, Nr-1
    da2da(jr, jt) = Jcbn2(jr, jt)*(psig(jr+1) - psig(jr)) * &
        (thetap_i(Nt+1) - thetap_i(Nt)) 
    dv2da(jr, jt) = GPI2*Rmaj2(jr, jt)*da2da(jr, jt)
enddo
do jt=2, nt
    do jr=1, Nr-1
        da2da(jr, jt) = Jcbn2(jr, jt)*(psig(jr+1) - psig(jr)) * &
            (thetap_i(jt)-thetap_i(jt-1)) 
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

dl_arc=0.0   
do j=1, Nr
    dl_arc = 0.0   
    do i=2, Nt 
        dl_arc(i) = r(j, i)*(thetap_i(i) - thetap_i(i-1)) !on the full grid
    enddo
    dl_arc(1) = r(j, 1)*(thetap_i(nt+1) - thetap_i(nt))
    perim(j) = sum(dl_arc)
enddo

gradPSI = 0.0
gradV   = 0.0
B_pol   = 0.0
B_R     = 0.0
B_Z     = 0.0

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
B_absa = sqrt(B_Ta**2. + B_pola**2.)

li3 = 2.*sum(B_pola**2.*dV2da)/rtor/(.4*GPI*iplasma)**2.
do i=1, Nr-1
    onez(i) = 0.5*(pressure(i) + pressure(i+1))
enddo
betapol = 2.*0.4*GPI*1.e-6*sum(onez*dva)/sum(B_pola**2.*dV2da)
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
    tar1 = 1./(Rmaj2(i, 1: Nt)**2.0)
    z1 = sum(tar1*tar2)
    G3(i) = z1
    z1 = sum(tar2)
    ONEZ(i) = z1
    tar1 = (gradVa(i, 1: Nt)/Rmaj2(i, 1: Nt))**2.0
    z1 = sum(tar1*tar2)
    G2(i) = z1

    tar1 = gradVa(i, 1: Nt)**2.0
    z1 = sum(tar1*tar2)
    G1(i) = z1

    tar1 = gradVa(i, 1: Nt)
    z1 = sum(tar1*tar2)
    GRADRO(i) = z1

    tar1 = B_ABSa(i, 1: Nt)**2.0
    z1 = sum(tar1*tar2)
    BDB02(i) = z1

    tar1 = B_ABSa(i, 1: Nt)
    z1 = sum(tar1*tar2)
    BDB0(i) = z1

    tar1 = 1./(B_ABSa(i, 1: Nt)**2.0)
    z1 = sum(tar1*tar2)

    B0DB2(i) = z1
    BMAXT(i) = maxval(B_ABSa(i, : ))
    BMINT(i) = minval(B_ABSa(i, : ))

    tar2 = B_ABSa(i, 1: Nt)/BMAXT(i)
    tar1 = ((btor/B_ABSa(i, 1: Nt))**2.0) * &
        (  1. - ((1. - tar2)**0.5) * (1. + 0.5*tar2)  )
    tar2 = dV2da(i, 1: Nt)/dVa(i)
    z1 = sum(tar1*tar2)
    FOFB(i)=z1

enddo

ip0 = Nr
ip1 = Nr-3
ip2 = Nr-2
ip3 = Nr-1

t1 = rhot(ip3)
t2 = rhot(ip2)
t3 = rhot(ip1)
t4 = rhoa(ip0)

call b_extrp_ef(t4, t1, t2, t3, ip3, ip2, ip1, ip0, Nr, G1    , k)
call b_extrp_ef(t4, t1, t2, t3, ip3, ip2, ip1, ip0, Nr, G2    , k)
call b_extrp_ef(t4, t1, t2, t3, ip3, ip2, ip1, ip0, Nr, G3    , k)
call b_extrp_ef(t4, t1, t2, t3, ip3, ip2, ip1, ip0, Nr, GRADRO, k)
call b_extrp_ef(t4, t1, t2, t3, ip3, ip2, ip1, ip0, Nr, FOFB  , k)
call b_extrp_ef(t4, t1, t2, t3, ip3, ip2, ip1, ip0, Nr, BMAXT , k)
call b_extrp_ef(t4, t1, t2, t3, ip3, ip2, ip1, ip0, Nr, BMINT , k)
call b_extrp_ef(t4, t1, t2, t3, ip3, ip2, ip1, ip0, Nr, BDB02 , k)
call b_extrp_ef(t4, t1, t2, t3, ip3, ip2, ip1, ip0, Nr, BDB0  , k)
call b_extrp_ef(t4, t1, t2, t3, ip3, ip2, ip1, ip0, Nr, B0DB2 , k)
      
call qinterp_ef_feqis(rhot(1: Nr-1), G1(1: Nr-1), Nr-1, &
    rhoa(2: Nr-1), dum1(2: Nr-1), Nr-2)
G1(2: Nr-1) = dum1(2: Nr-1)
call qinterp_ef_feqis(rhot(1: Nr-1), g2(1: Nr-1), Nr-1, &
    rhoa(2: Nr-1), dum1(2: Nr-1), Nr-2)
g2(2: Nr-1) = dum1(2: Nr-1)
call qinterp_ef_feqis(rhot(1: Nr-1), g3(1: Nr-1), Nr-1, &
    rhoa(2: Nr-1), dum1(2: Nr-1), Nr-2)
g3(2: Nr-1) = dum1(2: Nr-1)
call qinterp_ef_feqis(rhot(1: Nr-1), gradro(1: Nr-1), Nr-1, &
    rhoa(2: Nr-1), dum1(2: Nr-1), Nr-2)
gradro(2: Nr-1) = dum1(2: Nr-1)
call qinterp_ef_feqis(rhot(1: Nr-1), fofb(1: Nr-1), Nr-1, &
    rhoa(2: Nr-1), dum1(2: Nr-1), Nr-2)
 fofb(2: Nr-1) = dum1(2: Nr-1)
call qinterp_ef_feqis(rhot(1: Nr-1), bmaxt(1: Nr-1), Nr-1, &
    rhoa(2: Nr-1), dum1(2: Nr-1), Nr-2)
bmaxt(2: Nr-1) = dum1(2: Nr-1)
call qinterp_ef_feqis(rhot(1: Nr-1), bmint(1: Nr-1), Nr-1, &
    rhoa(2: Nr-1), dum1(2: Nr-1), Nr-2)
bmint(2: Nr-1) = dum1(2: Nr-1)
call qinterp_ef_feqis(rhot(1: Nr-1), bdb02(1: Nr-1), Nr-1, &
    rhoa(2: Nr-1), dum1(2: Nr-1), Nr-2)
bdb02(2: Nr-1) = dum1(2: Nr-1)
call qinterp_ef_feqis(rhot(1: Nr-1), bdb0(1: Nr-1), Nr-1, &
    rhoa(2: Nr-1), dum1(2: Nr-1), Nr-2)
bdb0(2: Nr-1) = dum1(2: Nr-1)
call qinterp_ef_feqis(rhot(1: Nr-1), b0db2(1: Nr-1), Nr-1, &
    rhoa(2: Nr-1), dum1(2: Nr-1), Nr-2)
b0db2(2: Nr-1) = dum1(2: Nr-1)

G1(1) = 0.0
G3(1) = 1./(Rmaj(1, 1)**2.0)
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

!--------------------------------------------------------------------
subroutine b_extrp_ef(x0, X1, X2, X3, j1, j2, j3, j, narr, yarr, jerr)

use debugger_ef, only: markloc_ef, debug

implicit none

integer, intent(in) :: j1, j2, j3, j, narr
double precision, intent(in) :: X0, X1, X2, X3 
double precision, intent(inout), dimension(narr) :: yarr
integer, intent(out) :: jerr

double precision :: u0, u1, u2, u3

call markloc_ef('b_extrp_ef', debug_lev=debug)

jerr = 0 
u1 = yarr(j1)
u2 = yarr(j2)
u3 = yarr(j3)
call EXTRP2_ef(x0, u0, X1, X2, X3, u1, u2, u3)
if (j > narr .or. j < 1) then
    write(*, *) 'Equil:  in b_extrp j = ', j, ' out of range 1 <', j, ' < ', narr
    jerr = 1
    return
endif 
yarr(j) = u0

return
end subroutine b_extrp_ef

!-------------------------------------------------------------------
subroutine EXTRP2_ef(X, F, X1, X2, X3, F1, F2, F3)

use debugger_ef, only: markloc_ef, debug

implicit none

double precision, intent(in) :: X, X1, X2, X3, F1, F2, F3
double precision, intent(out) :: F

double precision :: dfdx, d2fdx

call markloc_ef('extrp2_ef', debug_lev=debug)

dfdx = (f3 - f1)/(x3 - x1)
d2fdx = ( (f3 - f2)/(x3 - x2) - (f2 - f1)/(x2 - x1) )/(x3 - x1)

f = f2 + (x - x2)*( dfdx + d2fdx*(x - x2) )
      
return
end subroutine EXTRP2_ef
