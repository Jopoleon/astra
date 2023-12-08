subroutine PHI_EQ_2d_PBE(nrho, ntheta, psin_grid_in, iplasma, &
    ffprimp, pprimp, rbphi, Rb, Zb, Rax, Zax, psiax_in, PSIb, &
    solve_fix, i_prevc, &
! Output
    XX, YY, PSI, &
    psin_grid, lambda2d, thetap_out, &
    psiax, cnorm, X0, Y0, thetap_i_out, rmaj2, jcbn2, q_new, rhoedge, &
    darea2,epprim_out,efprim_out,r_min, yy2, gradr2, darea, ierr)

use pi_vars, only: GPI, GPI2, GPI4, muvac
implicit none

integer :: max_iter

integer, intent(in) :: nrho, ntheta, solve_fix, i_prevc !solve_fix = 0 uses max_iter =500 ; i_prevc = 0 to not use previous condition
double precision, intent(in) :: iplasma, rbphi, rax, zax, psiax_in, psib
double precision, intent(in) , dimension(nrho) :: ffprimp, pprimp, psin_grid_in
double precision, intent(in) , dimension(ntheta) :: Rb, Zb

integer, intent(out) :: ierr
double precision, intent(out) :: psiax, cnorm, X0, Y0, rhoedge
double precision, intent(out), dimension(nrho) :: psin_grid, q_new,epprim_out,efprim_out
double precision, intent(out), dimension(ntheta) :: thetap_out, thetap_i_out
double precision, intent(out), dimension(nrho, ntheta) :: XX, YY, rmaj2, jcbn2, &
    darea2, r_min, yy2, gradr2, dArea

double precision, intent(inout), dimension(nrho, ntheta) :: Psi, lambda2d

logical :: relambda_flag=.FALSE.
integer :: jthe, jrho, jrho_axis, jthe_axis, &
    jiter, Ndims, LDAB, nan_count, info, jloc, jmin(2)
double precision :: X0o, Y0o, denom, axis_change, dphi, qhalf
double precision, dimension(nrho) :: ddr, ddr_i, PSIn_gridp, effprimp, epprimp, &
    fpol, fpol2, phi_flux
double precision, dimension(ntheta) :: dtp, dtm, dt_i
double precision, dimension(ntheta+1) :: thetap, thetap_i
double precision, dimension(nrho, ntheta) :: lambda2dp, known_term, &
    dArc_rp1, dArc_rm1, dArc_rpt1, dArc_rmt1, &
    dArc_tp1, dArc_tm1, dArc_tpr1, dArc_tmr1
double precision :: gpsi(2*ntheta+1), work(2*(2*ntheta+1)*6), matrix(2*ntheta+1, 6)

Ndims = 1 + (nrho - 2)*ntheta
LDAB = 6*ntheta + 1
psin_grid = psin_grid_in

ierr=0

do jrho=1, nrho
    epprimp(jrho)  = -GPI4*muvac*pprimp(jrho)
    effprimp(jrho) = -GPI4*ffprimp(jrho)
    if (i_prevc == 0) then
        lambda2d(jrho, :) = (jrho - 1.)/(nrho - 1.)
        PSI(jrho, :) = psin_grid(jrho)
    endif
enddo
do jrho=1, nrho-1
    psin_gridp(jrho) = 0.5*(psin_grid(jrho+1) + psin_grid(jrho))
    lambda2dp(jrho, :) = 0.5*(lambda2d(jrho+1, :) + lambda2d(jrho, :))
enddo

X0  = Rax
Y0  = Zax
X0o = X0
Y0o = Y0

jrho_axis = 1
jthe_axis = 1

psiax = psiax_in

if (solve_fix > 0) then
    max_iter = solve_fix
else if (solve_fix == -2) then
    max_iter=250  ! uses fbe
else
    max_iter=250
endif

iter_loop: do jiter=1, max_iter+1

! recalculate psin_grid based on ffprime
    if (jiter >= 2 .and. (jrho_axis == 1 .and. jthe_axis == 1)) then
        fpol(nrho) = rbphi
        fpol(nrho-1) = sqrt(rbphi**2 - ffprimp(nrho) * (psin_grid(nrho) - psin_grid(nrho-1)) * (psib - psiax))
        do jrho=nrho-2, 1, -1
            fpol(jrho) = sqrt(fpol(jrho+1)**2 - ffprimp(jrho+1) * &
                (psin_grid(jrho+2) - psin_grid(jrho)) * (psib - psiax))
        enddo
        fpol2(1) = 0.
        phi_flux(1) = 0.
        do jrho=2, nrho ! toroidal flux on full grid
            dphi = fpol(jrho-1) * sum(dArea2(jrho-1, :)/Rmaj2(jrho-1, :))
            phi_flux(jrho)=phi_flux(jrho-1)+dphi
            qhalf = dphi/(psin_grid(jrho) - psin_grid(jrho-1))
            q_new(jrho-1)=qhalf/(psib-psiax)
            fpol2(jrho) = fpol2(jrho-1) + (2.*jrho - 3.)/qhalf
        enddo

        rhoedge = sqrt(phi_flux(nrho)/GPI)

        psin_grid = 0.5*psin_grid + 0.5*fpol2/fpol2(nrho)
        do jrho=1, nrho-1
            psin_gridp(jrho) = 0.5*(psin_grid(jrho+1) + psin_grid(jrho))
        enddo
    endif

    call jacobians(nrho, ntheta, Rb, Zb, X0, Y0, &
        lambda2d, lambda2dp, psin_grid, &
        psiax, psib, relambda_flag, psin_gridp, PSI, &
        dArea, Rmaj2, dArea2, &
        dArc_rp1, dArc_rm1, dArc_rpt1, dArc_rmt1, &
        dArc_tp1, dArc_tm1, dArc_tpr1, dArc_tmr1, &
        ddr, ddr_i, dtp, dtm, dt_i, XX, YY, &
        thetap, thetap_i, lambda2d, lambda2dp, jcbn2, r_min, yy2, gradr2)

    do jthe=1, ntheta
        known_term(:, jthe) = (effprimp(:) * dArea(:, jthe)/XX(:, jthe) + &
            XX(:, jthe)*epprimp(:)*dArea(:, jthe))
    enddo

    cnorm = sum(known_term)/iplasma/GPI2/0.4/GPI
    known_term = known_term/cnorm

    call solver_inversion_matrix_gsef(PSIb, nrho, ntheta, &
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
            do jthe=1, ntheta
                jloc = jthe + 1 + (jrho - 2)*ntheta
                matrix(jloc, 1) =  XX(jrho, jthe)**2
                matrix(jloc, 2) =  XX(jrho, jthe)
                matrix(jloc, 4) =  YY(jrho, jthe)**2
                matrix(jloc, 5) =  YY(jrho, jthe)
                matrix(jloc, 6) =  XX(jrho, jthe)*YY(jrho, jthe)
                gpsi  (jloc)    = PSI(jrho, jthe)
            enddo
        enddo

! Lapack DGELS
        call dgels('N', 2*ntheta + 1, 6, 1, matrix, 2*ntheta + 1, gpsi, 2*ntheta + 1, WORK, 2*(2*ntheta + 1)*6, INFO)

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

    relambda_flag = (jrho_axis == 1 .and. jthe_axis == 1)

    nan_count = 0
    do jrho=1, nrho
        do jthe=2, ntheta
            if (psi(jrho, jthe) /= psi(jrho, jthe)) nan_count = nan_count + 1
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
    if (axis_change < 1.e-8) then
        write(*, '(A, i3)') 'FEQIS converged, step #', jiter
        EXIT iter_loop
    endif
enddo iter_loop

if (jiter >= max_iter-1 .and. solve_fix == 0) ierr = 1

thetap_out   = thetap  (1: ntheta)
thetap_i_out = thetap_i(1: ntheta)
epprim_out = epprimp
efprim_out = effprimp

return
end subroutine PHI_EQ_2d_PBE
