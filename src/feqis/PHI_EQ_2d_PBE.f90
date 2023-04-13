subroutine PHI_EQ_2d_PBE(Nr, Nt, psin_grid, iplasma, &
    pressure, ffprimp, pprimp, btor, r0, Rb, Zb, Rax, Zax, &
    PSIb, IPOL, XX, YY, PSI, PSIxx, &
    g2, G3, r_out, r_in, volum, G1, G41, GRADRO, &
    BMAXT, BMINT, BDB02, BDB0, B0DB2, FOFB, &
    areat, perim, shif, elon, slat, tria, thetappp, &
    rmin, jrhoteta, nonegcurr, li3, betapol)

use pi_vars, only: GPI, GPI2, GPI4, MUVAC
use parameters_gsef, only: max_iter, diagnostic_gsef, name_gsefdir, &
    file_eqout, file_fields
use metric_coefficients_pbe, only: lambda2d, lambda2dp

implicit none

integer, intent(in) :: Nr, Nt, nonegcurr
double precision, intent(in) :: iplasma, R0, btor, li3, betapol
double precision, intent(in) , dimension(Nr) :: pressure, ipol, &
    g1, g2, g3, volum, gradro, bmaxt, bmint, bdb02, bdb0, b0db2, &
    fofb, areat, perim, slat
double precision, intent(in) , dimension(Nt) :: Rb, Zb
double precision, intent(out), dimension(Nr) :: PSIxx, r_out, r_in, &
    g41, shif, elon, tria
double precision, intent(out), dimension(Nt) :: thetappp
double precision, intent(out), dimension(Nr, Nt) :: Psi, rmin, jrhoteta, XX, YY
double precision, intent(inout) :: PSIb, rax, zax
double precision, intent(inout), dimension(Nr) :: psin_grid, &
    ffprimp, pprimp

integer :: i, i1, i2, j, jr, jt, i_call_save, ji, j_ok, iax, jax, &
    jiterext, jcall, Ndims, LDAB, nan_count
double precision :: X0, Y0, X0o, Y0o, epssol, epslambda, cnorm, psiax, &
    area, UPDWN, yrr, ya, t1, &
    yrmax, yrmin, yzmax, yzmin, yrzmax, yrzmin
double precision, dimension(3) :: xxxx1, yyyy1, pppp1
double precision, dimension(300, 300) :: psisave
double precision, dimension(Nr) :: PSIn_gridp, effprimp, epprimp, r
double precision, dimension(Nt+1) :: thetap, thetap_i
double precision, dimension(Nr, Nt) :: PSI_imd, dArea, Rmaj, Rmaj2, &
    known_term, dt_i, psio, &
    dArc_rp1, dArc_rm1, dArc_rpt1, dArc_rmt1, &
    dArc_tp1, dArc_tm1, dArc_tpr1, dArc_tmr1, &
    ddr, ddr_i, dtp, dtm, r_a, Rmaji, Rmaji1, rr2, r_i, r_i1, &
    gradr, gradh, gradr2(Nr, nt), gradh2, dArea2, &
    B_R, B_Z, B_T, &
    grt, grt2, dphdr, dphdz, ghht2, X_i1, Y_i1

character(len=120) :: fname

data i_call_save /0/
save i_call_save
save PSIsave

if (nonegcurr /= 0) then
    do j=1, Nr
        pprimp( j) = -max(0., -pprimp( j))
        ffprimp(j) = -max(0., -ffprimp(j))
    enddo
endif 

do j=1, Nr
    epprimp(j)  = -GPI4*muvac*pprimp(j)
    effprimp(j) = -GPI4*ffprimp(j)
enddo

epssol    = 1
epslambda = 1

if (i_call_save == 0) then
    PSI_imd = 0.0
    do jr=1, Nr
        do jt=1, Nt
            lambda2d( jr, jt) = (jr - 1      )/(Nr - 1.)
            lambda2dp(jr, jt) = (jr - 1 + 0.5)/(Nr - 1.)
        enddo
    enddo 
  
    PSI = 0.0
    do j=1, Nt
        PSI(1: Nr, j) = psin_grid(1: Nr)
    enddo
else
     psi(1: Nr, 1: Nt) = psisave(1: Nr, 1: Nt)
endif 

do jr=1, Nr-1
    psin_gridp(jr) = 0.5*(psin_grid(jr+1) + psin_grid(jr))
enddo

psio = psi  
X0  = Rax
Y0  = Zax
X0o = X0
Y0o = Y0

jcall = 0
iax   = 1
jax   = 1
j_ok  = 0

!External iterations
max_iter = 250
 
iter_loop: do jiterext=1, max_iter

! recalculate psin_grid based on ffprime
    if (jiterext >= 2 .and. (iax == 1 .and. jax == 1)) then
        gradh(Nr, 1) = btor*r0
        gradh(Nr-1, 1) = sqrt((btor*r0)**2. - ffprimp(Nr) * (psin_grid(Nr) - psin_grid(Nr-1)) * (psib - psiax))
        do j=Nr-2, 1, -1
            gradh(j, 1) = sqrt(gradh(j+1, 1)**2. - ffprimp(j+1) * &
                (psin_grid(j+2) - psin_grid(j)) * (psib - psiax))
        enddo

        tria(1) = 0.
        do jr=2, Nr ! toroidal flux on full grid
            tria(jr) = tria(jr-1) + gradh(jr-1, 1) * &
                sum(dArea2(jr-1, 1: Nt)/Rmaj2(jr-1, 1: Nt))
        enddo

        do jr=1, Nr-1 ! safety factor at half grid
            elon(jr) = (tria(jr+1) - tria(jr))/(psin_grid(jr+1) - psin_grid(jr)) / &
                (psib - psiax)
        enddo

        gradh2(1, 1) = 0.
        do jr=2, Nr ! new psin grid
            gradh2(jr, 1) = gradh2(jr-1, 1) + (2.*jr - 3.)/elon(jr-1)/(Nr - 1.)**2.
        enddo
        psin_grid = 0.5*psin_grid + 0.5*gradh2(:, 1)/gradh2(Nr, 1)

        do jr=1, Nr-1
            psin_gridp(jr) = 0.5*(psin_grid(jr+1) + psin_grid(jr))
        enddo
    endif

    call build_2dgrid(Nr, Nt, Rb, Zb, X0, Y0, &
        lambda2d(1: Nr, 1: Nt), lambda2dp(1: Nr, 1: Nt), psin_grid, &
        psiax, psib, j_ok, psin_gridp, PSI, &
        dArea, Rmaj, Rmaj2, Rmaji, Rmaji1, dArea2, &
        dArc_rp1, dArc_rm1, dArc_rpt1, dArc_rmt1, &
        dArc_tp1, dArc_tm1, dArc_tpr1, dArc_tmr1, &
        ddr, ddr_i, dtp, dtm, dt_i, Area, XX, YY, X_i1, Y_i1, &
        r_a, thetap, rr2, thetap_i, r_i, r_i1, dphdr, dphdz, &
        gradr, gradh, grt, gradr2, gradh2, grt2, ghht2)

    do jt=1, Nt
        do jr=1, Nr
            known_term(jr, jt) = (effprimp(jr) * dArea(jr, jt)/Rmaj(jr, jt) + &
                Rmaj(jr, jt)*epprimp(jr)*dArea(jr, jt))
        enddo
    enddo

    cnorm = sum(known_term)/iplasma/GPI2/0.4/GPI
    known_term = known_term/cnorm

    do jt=1, Nt
        PSI(Nr, jt) = PSIb
    enddo

    Ndims = 1 + (Nr - 2)*Nt
    LDAB = 6*Nt + 1

    call solver_inversion_matrix_gsef(PSIb, Nr, Nt, &
        known_term, Ndims, LDAB, &
        dArc_rp1, dArc_rm1, dArc_rpt1, dArc_rmt1, &
        dArc_tp1, dArc_tm1, dArc_tpr1, dArc_tmr1, &
        ddr, ddr_i, dtp, dtm, dt_i, PSI)

    call find_new_X0Y0(Nr, Nt, PSI, XX, YY, X0, Y0, PSIax, iax, jax)

    if (iax == 1 .and. jax == 1) then
        j_ok = 1   
    else
        j_ok = 0
    endif

    epssol = sum(abs(psi - psio))/Nr/Nt
    nan_count = 0
    do jr=1, Nr
        do jt=2, Nt
            if (psi(jr, jt) /= psi(jr, jt)) nan_count = nan_count + 1
        enddo
    enddo
    psio = psi

    if ((abs(X0) >= 20.) .and. (abs(y0) >= 20.)) call err_catch_a

!Check convergence
    t1 = abs(x0 - x0o)/x0 + abs(y0 - y0o)/x0
    X0o = X0
    Y0o = Y0
    if (t1 < 1.e-6) then
        EXIT ! Convergence
    endif
enddo iter_loop

psisave(1: Nr, 1: Nt) = psi(1: Nr, 1: Nt)

call build_2dgrid(Nr, Nt,  Rb, Zb,  X0, Y0, &
    lambda2d(1: Nr, 1: Nt), lambda2dp(1: Nr, 1: Nt), &
    psin_grid, psiax, psib, j_ok,  psin_gridp, PSI, &
    dArea, Rmaj, Rmaj2,  Rmaji, Rmaji1, dArea2, &
    dArc_rp1, dArc_rm1, dArc_rpt1, dArc_rmt1, &
    dArc_tp1, dArc_tm1, dArc_tpr1, dArc_tmr1, &
    ddr, ddr_i, dtp, dtm, dt_i, Area, XX, YY, X_i1, Y_i1, &
    r_a, thetap, rr2, thetap_i, r_i, r_i1, dphdr, dphdz, &
    gradr, gradh, grt, gradr2, gradh2, grt2, ghht2)

!regrid
call build_2dgrid2(Nr, Nt, psin_grid, &
    Rmaj, Rmaj2, r_a, thetap_i, gradr2, gradh2, &
    PSI, r0, pressure, btor, ipol, iplasma, &
    G2, G3, areat, perim, volum, G1, GRADRO, &
    BMAXT, BMINT, BDB02, BDB0, B0DB2, FOFB, &
    slat, li3, betapol)

rmin(1: Nr, 1: Nt) = r_a(1: Nr, 1: Nt)

do jt=1, Nt
    do jr=1, Nr
        jrhoteta(jr, jt) = -GPI2*(ffprimp(jr) * 1./Rmaj(jr, jt)/0.4/GPI + &
            Rmaj(jr, jt) * 1.e-6*pprimp(jr))/cnorm
    enddo
enddo

UPDWN = YY(1, 1)

do ji=1, Nr
    i = minloc(yy(ji, 1: Nt), 1)  
    i1 = i - 1
    i2 = i + 1
    if (i == 1 ) i1 = Nt
    if (i == Nt) i2 = 1
    xxxx1(1) = xx(ji, i1)
    xxxx1(2) = xx(ji, i)
    xxxx1(3) = xx(ji, i2)
    yyyy1(1) = yy(ji, i1)
    yyyy1(2) = yy(ji, i)
    yyyy1(3) = yy(ji, i2)
    call polyfitcc_feqis(xxxx1, yyyy1, pppp1)
    yrzmin = -pppp1(2)/(2.*pppp1(1))
    yzmin = pppp1(1)*yrzmin**2.0 + pppp1(2)*yrzmin + pppp1(3)

    i = maxloc(yy(ji, 1: Nt), 1)  
    i1 = i - 1
    i2 = i + 1
    if (i ==  1) i1 = Nt
    if (i == Nt) i2 = 1
    xxxx1(1) = xx(ji, i1)
    xxxx1(2) = xx(ji, i)
    xxxx1(3) = xx(ji, i2)
    yyyy1(1) = yy(ji, i1)
    yyyy1(2) = yy(ji, i)
    yyyy1(3) = yy(ji, i2)
    call polyfitcc_feqis(xxxx1, yyyy1, pppp1)
    yrzmax = -pppp1(2)/(2.*pppp1(1))
    yzmax = pppp1(1)*yrzmax**2.0 + pppp1(2)*yrzmax + pppp1(3)

    yrmin = MINVAL(xx(ji, :))
    yrmax = MAXVAL(xx(ji, :))

    yrr = .5*(yrmax + yrmin)
    ya  = .5*(yrmax - yrmin)

    SHIF (ji) = yrr - R0
    ELON(ji) = (yzmax - yzmin)/(yrmax - yrmin)
    TRIA(ji) = (yrr - 0.5d0*(yrzmin + yrzmax))/ya
    r_out(ji) = yrmax
    r_in(ji) = yrmin
enddo

r_out(1) = xx(1, 1)
r_in(1) = xx(1, 1)
ELON(1) = ELON(2)
TRIA(1) = 0.d0
SHIF (1) = XX(1, 1) - R0

G41 = G1 ! to be fixed

i_call_save = 1

rax = X0
zax = Y0

psixx(1: Nr) = PSI(1: Nr, 1) 

thetappp(1: Nt) = thetap(1: Nt)

return
end subroutine PHI_EQ_2d_PBE
