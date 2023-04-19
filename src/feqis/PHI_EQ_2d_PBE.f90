subroutine PHI_EQ_2d_PBE(Nr, Nt, psin_grid_in, iplasma, &
    pressure, ffprimp, pprimp, btor, r0, Rb, Zb, Rax, Zax, PSIb, IPOL, &
    XX, YY, PSI, &
    psin_grid_out, g2, G3, r_out, r_in, volum, G1, G41, GRADRO, &
    BMAXT, BMINT, BDB02, BDB0, B0DB2, FOFB, &
    areat, perim, shif, elon, slat, tria, thetap_out, &
    rmin, jrhoteta, li3, betapol)

use pi_vars, only: GPI, GPI2

implicit none

integer, parameter :: max_iter=500
double precision, parameter :: GPI4=GPI2**2, muvac=4.e-7*GPI

integer, intent(in) :: Nr, Nt
double precision, intent(in) :: iplasma, R0, btor, psib
double precision, intent(in) , dimension(Nr) :: pressure, ipol, ffprimp, pprimp, psin_grid_in
double precision, intent(in) , dimension(Nt) :: Rb, Zb

double precision, intent(out) :: li3, betapol
double precision, intent(out), dimension(Nr) :: psin_grid_out, g1, g2, g3, g41, &
    volum, gradro, bmaxt, bmint, bdb02, bdb0, b0db2, fofb, &
    areat, perim, slat, r_out, r_in, shif, elon, tria
double precision, intent(out), dimension(Nt) :: thetap_out
double precision, intent(out), dimension(Nr, Nt) :: Psi, rmin, jrhoteta, XX, YY
double precision, intent(inout) :: rax, zax

integer :: i, i1, i2, j, jt, jrho, j_ok, jrho_axis, jthe_axis, &
    jiter, Ndims, LDAB, nan_count, info, jloc, jmin(2)
double precision :: X0, Y0, X0o, Y0o, cnorm, denom, psiax, &
    yrr, ya, axis_change, &
    yrmax, yrmin, yzmax, yzmin, yrzmax, yrzmin
double precision, dimension(3) :: xxxx1, yyyy1, pppp1
double precision, dimension(Nr) :: PSIn_gridp, effprimp, epprimp, psin_grid
double precision, dimension(Nt+1) :: thetap, thetap_i
double precision, dimension(Nr, Nt) :: dArea, Rmaj2, &
    known_term, lambda2d, lambda2dp, &
    dArc_rp1, dArc_rm1, dArc_rpt1, dArc_rmt1, &
    dArc_tp1, dArc_tm1, dArc_tpr1, dArc_tmr1, &
    ddr, ddr_i, dtp, dtm, dt_i, gradh, gradr2, gradh2, dArea2
double precision :: gpsi(2*Nt+1), work(2*(2*Nt+1)*6), matrix(2*Nt+1, 6)
    
Ndims = 1 + (Nr - 2)*Nt
LDAB = 6*Nt + 1
psin_grid = psin_grid_in

do j=1, Nr
    epprimp(j)  = -GPI4*muvac*pprimp(j)
    effprimp(j) = -GPI4*ffprimp(j)
enddo

do jrho=1, Nr
    do jt=1, Nt
        lambda2d(jrho, jt) = (jrho - 1.)/(Nr - 1.)
    enddo
enddo
lambda2dp = lambda2d + 0.5/(Nr - 1.)

do j=1, Nt
    PSI(1: Nr, j) = psin_grid(1: Nr)
enddo

do jrho=1, Nr-1
    psin_gridp(jrho) = 0.5*(psin_grid(jrho+1) + psin_grid(jrho))
enddo

X0  = Rax
Y0  = Zax
X0o = X0
Y0o = Y0

jrho_axis = 1
jthe_axis = 1
j_ok  = 0

psiax = 0.

iter_loop: do jiter=1, max_iter

! recalculate psin_grid based on ffprime
    if (jiter >= 2 .and. (jrho_axis == 1 .and. jthe_axis == 1)) then
        gradh(Nr, 1) = btor*r0
        gradh(Nr-1, 1) = sqrt((btor*r0)**2 - ffprimp(Nr) * (psin_grid(Nr) - psin_grid(Nr-1)) * (psib - psiax))
        do j=Nr-2, 1, -1
            gradh(j, 1) = sqrt(gradh(j+1, 1)**2 - ffprimp(j+1) * &
                (psin_grid(j+2) - psin_grid(j)) * (psib - psiax))
        enddo

        tria(1) = 0.
        do jrho=2, Nr ! toroidal flux on full grid
            tria(jrho) = tria(jrho-1) + gradh(jrho-1, 1) * sum(dArea2(jrho-1, 1: Nt)/Rmaj2(jrho-1, 1: Nt))
        enddo

        do jrho=1, Nr-1 ! safety factor at half grid
            elon(jrho) = (tria(jrho+1) - tria(jrho))/(psin_grid(jrho+1) - psin_grid(jrho)) / &
                (psib - psiax)
        enddo

        gradh2(1, 1) = 0.
        do jrho=2, Nr ! new psin grid
            gradh2(jrho, 1) = gradh2(jrho-1, 1) + (2.*jrho - 3.)/elon(jrho-1)/(Nr - 1.)**2
        enddo
        psin_grid = 0.5*psin_grid + 0.5*gradh2(:, 1)/gradh2(Nr, 1)

        do jrho=1, Nr-1
            psin_gridp(jrho) = 0.5*(psin_grid(jrho+1) + psin_grid(jrho))
        enddo
    endif

    call build_2dgrid(Nr, Nt, Rb, Zb, X0, Y0, &
        lambda2d, lambda2dp, psin_grid, &
        psiax, psib, j_ok, psin_gridp, PSI, &
        dArea, Rmaj2, dArea2, &
        dArc_rp1, dArc_rm1, dArc_rpt1, dArc_rmt1, &
        dArc_tp1, dArc_tm1, dArc_tpr1, dArc_tmr1, &
        ddr, ddr_i, dtp, dtm, dt_i, XX, YY, &
        rmin, thetap, thetap_i, gradh, gradr2, gradh2)

    do jt=1, Nt
        known_term(1: Nr, jt) = (effprimp(1: Nr) * dArea(1: Nr, jt)/XX(1: Nr, jt) + &
            XX(1: Nr, jt)*epprimp(1: Nr)*dArea(1: Nr, jt))
    enddo

    cnorm = sum(known_term)/iplasma/GPI2/0.4/GPI
    known_term = known_term/cnorm

    call solver_inversion_matrix_gsef(PSIb, Nr, Nt, &
        known_term, Ndims, LDAB, &
        dArc_rp1, dArc_rm1, dArc_rpt1, dArc_rmt1, &
        dArc_tp1, dArc_tm1, dArc_tpr1, dArc_tmr1, &
        ddr, ddr_i, dtp, dtm, dt_i, PSI)

!----------------------------------
! Find new magnetic axis: R, z, psi

    jmin = minloc(PSI)
    jrho_axis = jmin(1)
    jthe_axis = jmin(2)

    if (jrho_axis == 1) then
        gpsi(1)      = PSI(1, 1)
        matrix(1, 1) =  XX(1, 1)**2
        matrix(1, 2) =  XX(1, 1)
        matrix(:, 3) = 1.
        matrix(1, 4) =  YY(1, 1)**2
        matrix(1, 5) =  YY(1, 1)
        matrix(1, 6) =  XX(1, 1)*YY(1, 1)
        do jrho=2, 3
            do jt=1, Nt
                jloc = jt + 1 + (jrho - 2)*Nt
                matrix(jloc, 1) =  XX(jrho, jt)**2
                matrix(jloc, 2) =  XX(jrho, jt)
                matrix(jloc, 4) =  YY(jrho, jt)**2
                matrix(jloc, 5) =  YY(jrho, jt)
                matrix(jloc, 6) =  XX(jrho, jt)*YY(jrho, jt)
                gpsi  (jloc)    = PSI(jrho, jt)
            enddo
        enddo

! Lapack DGELS
        call dgels('N', 2*Nt + 1, 6, 1, matrix, 2*Nt + 1, gpsi, 2*Nt + 1, WORK, 2*(2*Nt + 1)*6, INFO)

        denom = 4*gpsi(1)*gpsi(4) - gpsi(6)**2
        x0 = (gpsi(6)*gpsi(5) - 2*gpsi(4)*gpsi(2))/denom
        y0 = (gpsi(6)*gpsi(2) - 2*gpsi(1)*gpsi(5))/denom
        psiax = gpsi(1)*x0**2 + gpsi(2)*x0 + gpsi(3) +  &
                gpsi(4)*y0**2 + gpsi(5)*y0 + gpsi(6)*x0*y0
    else
        x0 = XX(jrho_axis, jthe_axis)
        y0 = YY(jrho_axis, jthe_axis)
        psiax = PSI(jrho_axis, jthe_axis)
    endif

    if (jrho_axis == 1 .and. jthe_axis == 1) then
        j_ok = 1
    else
        j_ok = 0
    endif

    nan_count = 0
    do jrho=1, Nr
        do jt=2, Nt
            if (psi(jrho, jt) /= psi(jrho, jt)) nan_count = nan_count + 1
        enddo
    enddo

    if ((abs(X0) >= 20.) .and. (abs(y0) >= 20.)) then
        write(*, *) 'Error x0, y0 too large'
        stop
    endif

!------------------
! Check convergence

    axis_change = abs(x0 - x0o)/x0 + abs(y0 - y0o)/x0
    X0o = X0
    Y0o = Y0
    if (axis_change < 1.e-6) then
!        write(*, '(A, i3)') 'FEQIS converged, step #', jiter
        EXIT iter_loop
    endif
enddo iter_loop

call build_2dgrid(Nr, Nt,  Rb, Zb,  X0, Y0, &
    lambda2d, lambda2dp, &
    psin_grid, psiax, psib, j_ok, psin_gridp, PSI, &
    dArea, Rmaj2, dArea2, &
    dArc_rp1, dArc_rm1, dArc_rpt1, dArc_rmt1, &
    dArc_tp1, dArc_tm1, dArc_tpr1, dArc_tmr1, &
    ddr, ddr_i, dtp, dtm, dt_i, XX, YY, &
    rmin, thetap, thetap_i, gradh, gradr2, gradh2)

!regrid
call build_2dgrid2(Nr, Nt, psin_grid, &
    XX, Rmaj2, rmin, thetap_i, gradr2, gradh2, &
    PSI, r0, pressure, btor, ipol, iplasma, &
    G2, G3, areat, perim, volum, G1, GRADRO, &
    BMAXT, BMINT, BDB02, BDB0, B0DB2, FOFB, &
    slat, li3, betapol)

do jt=1, Nt
    jrhoteta(1: Nr, jt) = -GPI2*(ffprimp(1: Nr) * 1./XX(1: Nr, jt)/0.4/GPI + &
            XX(1: Nr, jt) * 1.e-6*pprimp(1: Nr))/cnorm
enddo

do jrho=2, Nr

    i = minloc(yy(jrho, :), 1)
    i1 = i - 1
    i2 = i + 1
    if (i == 1 ) i1 = Nt
    if (i == Nt) i2 = 1
    xxxx1(1) = xx(jrho, i1)
    xxxx1(2) = xx(jrho, i)
    xxxx1(3) = xx(jrho, i2)
    yyyy1(1) = yy(jrho, i1)
    yyyy1(2) = yy(jrho, i)
    yyyy1(3) = yy(jrho, i2)
    call polyfitcc_feqis(xxxx1, yyyy1, pppp1)
    yrzmin = -pppp1(2)/(2.*pppp1(1))
    yzmin = pppp1(1)*yrzmin**2 + pppp1(2)*yrzmin + pppp1(3)

    i = maxloc(yy(jrho, :), 1)
    i1 = i - 1
    i2 = i + 1
    if (i ==  1) i1 = Nt
    if (i == Nt) i2 = 1
    xxxx1(1) = xx(jrho, i1)
    xxxx1(2) = xx(jrho, i)
    xxxx1(3) = xx(jrho, i2)
    yyyy1(1) = yy(jrho, i1)
    yyyy1(2) = yy(jrho, i)
    yyyy1(3) = yy(jrho, i2)
    call polyfitcc_feqis(xxxx1, yyyy1, pppp1)
    yrzmax = -pppp1(2)/(2.*pppp1(1))
    yzmax = pppp1(1)*yrzmax**2 + pppp1(2)*yrzmax + pppp1(3)

    yrmin = MINVAL(xx(jrho, :))
    yrmax = MAXVAL(xx(jrho, :))

    yrr = 0.5*(yrmax + yrmin)
    ya  = 0.5*(yrmax - yrmin)

    SHIF (jrho) = yrr - R0
    ELON (jrho) = (yzmax - yzmin)/(yrmax - yrmin)
    TRIA (jrho) = (yrr - 0.5*(yrzmin + yrzmax))/ya
    r_out(jrho) = yrmax
    r_in (jrho) = yrmin
enddo

r_out(1) = xx(1, 1)
r_in(1)  = xx(1, 1)
ELON(1)  = ELON(2)
TRIA(1)  = 0.d0
SHIF(1)  = XX(1, 1) - R0

G41 = G1 ! to be fixed

thetap_out(1: Nt) = thetap(1: Nt)
psin_grid_out = psin_grid
rax = X0
zax = Y0

return
end subroutine PHI_EQ_2d_PBE
