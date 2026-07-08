subroutine build_2dgrid(nrho, ntheta, psin_grid, &
    PSI, jrho2, darea2, YY2, &
    rtor, pressure, btor, ipol, iplasma, &
    XX, YY, Rmaj2, rmin, Jcbn2, thetap_i, &
    q_new, rhoedge_in, gradr2, &
! Output
    G2, G3, areat, perim, volum, G1, GRADRO, &
    BMAXT, BMINT, BDB02, BDB0, B0DB2, FOFB, &
    slat, li3, betapol, psplex, bpcell, bcell, r_out, r_in, &
    elon, tria_u, tria_l, shif, g41, q_out, shiv, square, li_aug, betapol_iter, dl_dt, &
    wkin, bpkin, dpsidvbez)

use pi_vars, only: GPI, GPI2
use numerical_tools, only: qinterp, extrapolate, polyfitcc
use feqis_tools, only: green_function_includingsamepoint, pol_angle
use metric_coefficients_pbe, only: fsa_kernel

implicit none

integer, intent(in) :: ntheta, nrho
double precision, intent(in) :: btor, iplasma, Rtor, rhoedge_in
double precision, intent(in), dimension(ntheta) :: thetap_i
double precision, intent(in), dimension(nrho) :: psin_grid, ipol, pressure, q_new
double precision, intent(in), dimension(nrho, ntheta) :: PSI, jrho2, darea2, yy2, dl_dt

double precision, intent(out) :: li3, betapol, psplex, li_aug, betapol_iter, wkin, bpkin
double precision, intent(out), dimension(nrho) :: G1, G2, G3, &
    volum, areat, perim, slat, &
    FOFB, GRADRO, BMAXT, BMINT, BDB02, BDB0, B0DB2, dpsidvbez
double precision, intent(in), dimension(nrho, ntheta) :: XX, YY, rmin, Rmaj2, & 
    jcbn2, gradr2
double precision, intent(out), dimension(nrho, ntheta) :: bpcell, bcell
double precision, intent(out), dimension(nrho) :: r_out, r_in, elon, shif, &
    g41, q_out, shiv, square, tria_u, tria_l

integer :: jrho, jthe, jthe_l, k, j, i, ji, i1, i2, ip0, ip1, ip2, ip3
double precision :: drdX, drdY, Mdet, dpsi, dthe, ipol_rmaj, z1, z2, z3, rho_interp, &
    dumba1, dumba2, dumba3, dumba4, qedge, rhoedge, greenf, t4, &
    yrzmin, yrzmax, yzmax, yrmin, yrmax, yrr, yzmin, ya
double precision, dimension(3) :: xxxx1, yyyy1, pppp1
double precision, dimension(nrho) :: rhot, rhoa, dPSIdV, dVa, daa, dum1, AMETR, ONEZ
double precision, dimension(ntheta) :: dl_arc, tar1, tar2, theta_special, dl_arc_special, &
    rmaj2_sq, gradVa_sq, B_ABSa_sq
double precision, dimension(nrho, ntheta) :: gradPSIa, gradVa, dV2da, dA2da, &
    B_pola, B_ABSa, B_Ta, dldt_temp

do jrho=1, nrho
    rhoa(jrho) = (jrho - 1.)/(nrho - 1.) ! full grid
enddo
rhot = rhoa + 0.5/(nrho - 1.)  ! half grid

rhoedge = rhoedge_in/sqrt(btor) !from PHI, rhoedge misses a sqrt(btor) at denominator

da2da = 0.
dv2da = 0.

! half grid points around area and volume (last point doesn exist)
jthe = 1
do jrho=1, nrho-1
    da2da(jrho, jthe) = Jcbn2(jrho, jthe)*(psin_grid(jrho+1) - psin_grid(jrho)) * &
        (thetap_i(1)+GPI2 - thetap_i(ntheta)) 
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
dldt_temp = 0.
do jrho=2, nrho
    volum(jrho) = volum(jrho-1) + dva(jrho-1)
    areat(jrho) = areat(jrho-1) + daa(jrho-1)
    dldt_temp(jrho, :) = 0.5*(dl_dt(jrho, :) + dl_dt(jrho-1, :))
enddo
ip0 = nrho
ip1 = nrho-3
ip2 = nrho-2
ip3 = nrho-1
t4 = rhoa(ip0)
do jthe=1, ntheta
    dldt_temp(nrho, jthe) = extrapolate(t4, ip3, ip2, ip1, ip0, rhot, dl_dt(:, jthe))
enddo

dl_arc = 0.0
do jrho=1, nrho
    dl_arc = 0.0
    dl_arc_special = 0.0
    square(jrho) = 0.0
    dumba1 = 0.5*(maxval(XX(jrho, :)) + minval(XX(jrho, :))) !Rgeo
    dumba2 = 0.5*(maxval(XX(jrho, :)) - minval(XX(jrho, :))) !a
    dumba3 = 0.5*(maxval(YY(jrho, :)) + minval(YY(jrho, :))) !Zgeo
    dumba4 = 0.5*(maxval(YY(jrho, :)) - minval(YY(jrho, :)))/dumba2 !elongation
    theta_special(1) = pol_angle(dumba1, dumba3/dumba4, XX(jrho, 1), YY(jrho, 1)/dumba4)
    do jthe=2, ntheta 
        theta_special(jthe) = pol_angle(dumba1, dumba3/dumba4, XX(jrho, jthe), YY(jrho, jthe)/dumba4)
        if (theta_special(jthe) < theta_special(jthe-1) - GPI2/ntheta) then
            theta_special(jthe) = theta_special(jthe) + GPI2
        endif
        z1 = sin(theta_special(jthe) + theta_special(jthe-1))*sin(0.5*(theta_special(jthe) + theta_special(jthe-1)))
        dl_arc(jthe) = dldt_temp(jrho, jthe)*(thetap_i(jthe) - thetap_i(jthe-1)) !on the full grid
        dl_arc_special(jthe) = dldt_temp(jrho, jthe)*(theta_special(jthe) - theta_special(jthe-1)) !on the full grid
        square(jrho) = square(jrho) + (XX(jrho, jthe) - dumba1)/dumba2*z1*dl_arc_special(jthe)
    enddo
    dl_arc(1) = dldt_temp(jrho, 1)*(thetap_i(1) + GPI2 - thetap_i(ntheta))
    z1 = sin((thetap_i(1) + GPI2 + thetap_i(ntheta)))*sin(0.5*(thetap_i(1) + GPI2 + thetap_i(ntheta)))
    square(jrho) = square(jrho) + (XX(jrho, 1)-dumba1)/dumba2*z1*dl_arc(1)
    perim(jrho)  = sum(dl_arc)
    square(jrho) = square(jrho)/perim(jrho)
enddo
square = 4.*square - 1.
square(1) = square(2)

!Compute Bpol, gradPSI, gradV

do jrho=1, nrho-1
    dPSIdv(jrho) = (PSI(jrho+1, 1) - PSI(jrho, 1))/(volum(jrho+1) - volum(jrho)) !on the psigp grid
enddo

! there are on half grid
do jrho=1, nrho-1
    onez(jrho) = 0.5*(IPOL(jrho) + IPOL(jrho+1))
enddo

ip0 = nrho
ip1 = nrho-3
ip2 = nrho-2
ip3 = nrho-1
t4 = rhoa(ip0)

qedge = extrapolate(t4, ip3, ip2, ip1, ip0, rhot, q_new)

q_out(1:nrho-1) = q_new(1:nrho-1)
q_out(nrho) = qedge

do jthe=1, ntheta 
    do jrho=1, nrho-1
        z1 = dVa(jrho)/(psin_grid(jrho+1) - psin_grid(jrho))
        z2 = (PSI(jrho+1, jthe) - PSI(jrho, jthe))/(psin_grid(jrho+1) - psin_grid(jrho))
        B_pola(jrho, jthe) = z2/(GPI2*Rmaj2(jrho, jthe))*sqrt(gradr2(jrho, jthe))
        bpcell(jrho, jthe) = B_pola(jrho, jthe)
        gradPSIa(jrho, jthe) = z2*sqrt(gradr2(jrho, jthe))
        gradVa  (jrho, jthe) = z1*sqrt(gradr2(jrho, jthe))
        B_Ta(jrho, jthe) = onez(jrho)/Rmaj2(jrho, jthe)
        bcell(jrho, jthe) = sqrt(B_pola(jrho, jthe)**2 + B_Ta(jrho, jthe)**2)
    enddo
    bcell  (nrho, jthe) = extrapolate(t4, ip3, ip2, ip1, ip0, rhot, bcell(:, jthe))   !note that these values could also use the last-1 half grid point value for more accurate extrapolation (todo)
    bpcell  (nrho, jthe) = extrapolate(t4, ip3, ip2, ip1, ip0, rhot, bpcell(:, jthe))
    gradVa  (nrho, jthe) = extrapolate(t4, ip3, ip2, ip1, ip0, rhot, gradVa(:, jthe))
    gradpsia(nrho, jthe) = extrapolate(t4, ip3, ip2, ip1, ip0, rhot, gradpsia(:, jthe))
    B_pola  (nrho, jthe) = extrapolate(t4, ip3, ip2, ip1, ip0, rhot, B_pola  (:, jthe))
enddo   
B_absa = sqrt(B_Ta**2 + B_pola**2)
!write(*,*) 'perim extrap', extrapolate(t4, ip3, ip2, ip1, ip0, rhoa(1:nrho), perim(1:nrho)), &
! perim(nrho),perim(nrho-1)

!calculate psplex, the external inductance of the plasma 
dumba1 = 0.
do j=1, ntheta
    do i=1, ntheta
        greenf = green_function_includingsamepoint(xx(nrho, i), yy(nrho, i), xx(nrho, j), yy(nrho, j), dl_arc(i), xx(nrho, i))
        dumba1 = dumba1 + greenf/xx(nrho, i)*gradpsia(nrho, i)*dl_arc(i)*dl_arc(j)    
    enddo
enddo

psplex = dumba1/sum(dl_arc)
psplex = psplex/(1.*sum(B_pola(nrho, 1:ntheta)*dl_arc(1:ntheta))/0.4) ! for LEXT part, alternative

!li3 = 2.*sum(B_pola**2 * dV2da)/rtor/(0.4*GPI*iplasma)**2
dumba1 = 0.5*(maxval(XX(nrho, 1:ntheta)) + minval(XX(nrho, 1:ntheta))) !Rgeo
li3 = 2.*sum(B_pola**2 * dV2da)/dumba1/(0.4*GPI*iplasma)**2
li_aug = li3 * dumba1 * perim(nrho)**2./(2.*volum(nrho))
!write(*,*) 'liaug',li3,dumba1,perim(nrho),volum(nrho)
do jrho=1, nrho-1
    onez(jrho) = 0.5*(pressure(jrho) + pressure(jrho+1))
enddo
betapol = 0.4*GPI2*1.e-6*sum(onez*dva)/sum(B_pola**2 * dV2da)
betapol_iter = 4.*1.e-6*sum(onez*dva)/(0.4*GPI*dumba1*iplasma**2.)
wkin = sum(onez*dva)  !int(pressure dV)
bpkin = sum(B_pola**2 * dV2da)/(0.4*GPI2*1.e-6) ! int(Bp**2/(2mu0)dV)

slat = 0.
do jrho=2, nrho
    do jthe=1, ntheta-1 
        slat(jrho) = slat(jrho) + XX(jrho, jthe)*dldt_temp(jrho, jthe)*(thetap_i(jthe+1) - thetap_i(jthe))
    enddo
    slat(jrho) = slat(jrho) + XX(jrho, ntheta)*dldt_temp(jrho, ntheta)*(thetap_i(1) - thetap_i(ntheta) + GPI2)
enddo
slat = slat*GPI2   !full grid

do jrho=1, nrho
    ametr(jrho) = 0.5*(maxval(XX(jrho, 1: ntheta)) - minval(XX(jrho, 1: ntheta)))
    shiv(jrho)  = 0.5*(maxval(YY(jrho, 1: ntheta)) + minval(YY(jrho, 1: ntheta)))  !avg between zmax and zmin or mean of Z 
enddo

!flux surface average kernel
do jrho=1, nrho-1
    fsa_kernel(jrho,1: ntheta) = dV2da(jrho, 1: ntheta)/dVa(jrho)
enddo

!Cycle over positions  ! half grid
do jrho=1, nrho-1
    rmaj2_sq  =  Rmaj2(jrho, :)**2
    gradVa_sq = gradVa(jrho, :)**2
    B_ABSa_sq = B_ABSa(jrho, :)**2

    G3(jrho)  = sum(fsa_kernel(jrho, :)/rmaj2_sq) ! this is the definition of flux surface average of <f> of f = tar1
    G41(jrho) = sum(rmaj2_sq*fsa_kernel(jrho, :))
    G2(jrho)  = sum(gradVa_sq/rmaj2_sq*fsa_kernel(jrho, :))
    G1(jrho)  = sum(gradVa_sq*fsa_kernel(jrho,: ))
    GRADRO(jrho) = sum(gradVa(jrho, :)*fsa_kernel(jrho, :))
    BDB02(jrho)  = sum(B_ABSa_sq*fsa_kernel(jrho, :))
    BDB0(jrho)   = sum(B_ABSa(jrho, :)*fsa_kernel(jrho, :))
    B0DB2(jrho)  = sum(fsa_kernel(jrho, :)/B_ABSa_sq)
    BMAXT(jrho) = maxval(B_ABSa(jrho, :))
    BMINT(jrho) = minval(B_ABSa(jrho, :))
    tar1 = btor**2/B_ABSa_sq * &
        ( 1. - (sqrt(1. - (B_ABSa(jrho, :)/BMAXT(jrho)))) * (1. + 0.5*(B_ABSa(jrho, :)/BMAXT(jrho))) )
    FOFB(jrho) = sum(tar1*fsa_kernel(jrho, :))
    dpsidvbez(jrho) = dPSIdv(jrho)
enddo

rho_interp = rhoa(nrho)

G1(nrho)     = EXTRAPOLATE(rho_interp, nrho-1, nrho-2, nrho-3, nrho, rhot, G1)
G2(nrho)     = EXTRAPOLATE(rho_interp, nrho-1, nrho-2, nrho-3, nrho, rhot, G2)
G3(nrho)     = EXTRAPOLATE(rho_interp, nrho-1, nrho-2, nrho-3, nrho, rhot, G3)
G41(nrho)    = EXTRAPOLATE(rho_interp, nrho-1, nrho-2, nrho-3, nrho, rhot, G41)
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
call qinterp(rhot(1: nrho-1), g41(1: nrho-1), nrho-1, rhoa(2: nrho-1), dum1(2: nrho-1), nrho-2)
g41(2: nrho-1) = dum1(2: nrho-1)
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
G3(1) = 1./XX(1, 1)**2
G2(1) = 0.0
G41(1) = XX(1, 1)**2
GRADRO(1) = 0.0
ipol_rmaj = IPOL(1)/XX(1, 1)
BDB02(1) = ipol_rmaj**2
BDB0(1)  = ipol_rmaj
B0DB2(1) = 1./ipol_rmaj**2
FOFB(1)  = (btor/ipol_rmaj)**2
BMAXT(1) = ipol_rmaj
BMINT(1) = ipol_rmaj

do ji=1, Nrho
    i = minloc(yy(ji, 1: Ntheta), 1)  
    i1 = i - 1
    i2 = i + 1
    if (i == 1 ) i1 = Ntheta
    if (i == Ntheta) i2 = 1
    xxxx1(1) = xx(ji, i1)
    xxxx1(2) = xx(ji, i)
    xxxx1(3) = xx(ji, i2)
    yyyy1(1) = yy(ji, i1)
    yyyy1(2) = yy(ji, i)
    yyyy1(3) = yy(ji, i2)
    call polyfitcc(xxxx1, yyyy1, pppp1)
    yrzmin = -pppp1(2)/(2.*pppp1(1))
    yzmin = pppp1(1)*yrzmin**2.0 + pppp1(2)*yrzmin + pppp1(3)

    i = maxloc(yy(ji, 1: Ntheta), 1)  
    i1 = i - 1
    i2 = i + 1
    if (i ==  1) i1 = Ntheta
    if (i == Ntheta) i2 = 1
    xxxx1(1) = xx(ji, i1)
    xxxx1(2) = xx(ji, i)
    xxxx1(3) = xx(ji, i2)
    yyyy1(1) = yy(ji, i1)
    yyyy1(2) = yy(ji, i)
    yyyy1(3) = yy(ji, i2)
    call polyfitcc(xxxx1, yyyy1, pppp1)
    yrzmax = -pppp1(2)/(2.*pppp1(1))
    yzmax = pppp1(1)*yrzmax**2.0 + pppp1(2)*yrzmax + pppp1(3)

    i = minloc(xx(ji, 1:Ntheta), 1)  
    yrmin = xx(ji, i)
    i = maxloc(xx(ji, 1:Ntheta), 1)  
    yrmax = xx(ji, i)

    yrr = 0.5*(yrmax + yrmin)
    ya  = 0.5*(yrmax - yrmin)

    SHIF (ji) = yrr - rtor
    ELON(ji) = (yzmax - yzmin)/(yrmax - yrmin)
    TRIA_U(ji) = (yrr - yrzmax)/ya
    TRIA_L(ji) = (yrr - yrzmin)/ya
    r_out(ji) = yrmax
    r_in(ji) = yrmin
enddo

r_out(1) = xx(1, 1)
r_in (1) = xx(1, 1)
ELON(1) = ELON(2)
TRIA_U(1) = 0.d0
TRIA_L(1) = 0.d0
SHIF(1) = XX(1, 1) - rtor

return
end subroutine build_2dgrid
