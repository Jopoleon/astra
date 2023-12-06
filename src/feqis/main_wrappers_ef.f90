subroutine full_system_advance_feqis(j_init)

use errors_params, only: err_epsilon, err_circ_plasma_iter
use feqis_circuit, only: nr2, nz2, nconduc, iplasma, psi_cur_old, &
    psiplasmatoconduc, curconduc, jrz, area_eff, psi_external_calc
use astra2fbe, only: fast_mode, execute_plasma
use parameters_a2equil, only: max_iter
use green_matrix, only: greeni

implicit none

integer, intent(inout) :: j_init

integer :: i, j_iter
double precision :: error_temp
double precision, dimension(300) :: cur_temp

! time stepping
! at iteration 0, dpsidt = 0
if (j_init == 0) then
    write(*, *) 'init full system'
! First do full equilibrium solution at time t=0
    call psi_external_calc
    call solve_gse2d_fbe_full_feqis(0)
    do i=1, nconduc
        psiplasmatoconduc(i) = sum(jrz(1: nr2, 1: nz2) * area_eff(1: nr2, 1: nz2) * greeni(1: nr2, 1: nz2, i))
    enddo
    psi_cur_old(1:nconduc) = psiplasmatoconduc(1:nconduc)

    write(*, *) 'init done'
    call circuit_eq_advance_feqis(0)
    j_init = -1
    return
endif

! Full iterations
write(*, *) 'iter full system'

if (fast_mode == 1 .and. execute_plasma == 1) then
    psi_cur_old(1:nconduc) = psiplasmatoconduc(1:nconduc)
    call psi_external_calc
    call solve_gse2d_fbe_full_feqis_1turn(1, 0, 0.d0, 0.d0)
    do i=1, nconduc
        psiplasmatoconduc(i) = sum(jrz(1: nr2, 1: nz2) * area_eff(1: nr2, 1: nz2) * greeni(1: nr2, 1: nz2, i))
    enddo
endif

do j_iter=1, 2*max_iter

    cur_temp(1:nconduc) = curconduc(1:nconduc)
    call circuit_eq_advance_feqis(1)
    if (fast_mode == 0) then
        call psi_external_calc
        call solve_gse2d_fbe_full_feqis_1turn(1, 0, 0.d0, 0.d0)
        do i=1, nconduc
            psiplasmatoconduc(i) = sum(jrz(1: nr2, 1: nz2) * area_eff(1: nr2, 1: nz2) * greeni(1: nr2, 1: nz2, i))
        enddo
    endif

    error_temp = sum(abs(cur_temp(1:nconduc) - curconduc(1:nconduc)))/(nconduc + err_epsilon)/iplasma

    if (j_iter == 1 .and. fast_mode == 0) then
        error_temp = 100.
    endif

    if (error_temp <= err_circ_plasma_iter) then
        if (execute_plasma == 1 .or. fast_mode == 0) then
            j_init = -1
        endif
        write(*, *) 'circuit eq done'
        return
    endif

    if (j_iter > max_iter) then
        write(*, *) 'circuit equations not converging, max number of iterations override... stopping', error_temp
        stop
    endif

enddo

return
end subroutine full_system_advance_feqis

!--------------------------------------------------------------------
subroutine solve_gse2d_fbe_full_feqis(j_init)

use errors_params, only: err_find_psistab
use feqis_circuit, only: nr, nz, nr2, nz2, psiextrz, redo_bnd, &
    r, z, dr, dz, dr_factor_init, dz_factor_init, rax, zax, raxp, zaxp, &
    trax, tzax, iaxis, jaxis, psistabr, psistabz, &
    restab_axis_with_fourier_wall, restab_boundary_with_fourier_wall, & !doesnt work well
    restab_F_function_full_fonfit
use astra2fbe, only: refit_mode, n_of_newton_iterations
use feqis_tools, only: closest_index

implicit none

integer, intent(in) :: j_init

integer :: j_iter, j_iter2, j_cyclo, jeppa
double precision :: temp_err, raxold, zaxold, temp_err2, raxoldo, zaxoldo, &
    raxtmp, zaxtmp, det, psistab1o, psistab2o, psro, pszo, dist1, dist2, &
    cibapr, cibazr, rleft, rright, zup, zdown, dcrdr, dcrdz, dczdr, dczdz
double precision, dimension(300, 300) :: g000

! First, initialized initial guess coming from prescribed boundary current density: jrhoteta

SELECT CASE(refit_mode)

CASE(-1) ! 1 turn only
    call solve_gse2d_fbe_full_feqis_1turn(j_init, 0, 0.d0, 0.d0)

CASE(0)
! Start iterations to find self-consistent solution
    g000(1:nr2, 1:nz2) = psiextrz(1:nr2, 1:nz2)
    iaxis = closest_index(raxp, r(1), dr)
    jaxis = closest_index(zaxp, z(1), dz)
    rax = r(iaxis)
    zax = z(jaxis)
    raxold = rax
    zaxold = zax
    raxoldo = rax
    zaxoldo = zax
    psistabR = 0.
    psistabZ = 0.
    psistab1o = 0.
    psistab2o = 0.
    temp_err = 100.
    temp_err2 = 100.
    rleft = raxold
    rright = raxold
    zup = zaxold
    zdown = zaxold
! outer cycle, calculate new axis
    j_cyclo = 0
    jeppa = 1
    dist1 = dr*dr_factor_init
    dist2 = dz*dz_factor_init

    do j_iter2=1, 10000000

        SELECT CASE(j_cyclo)
        CASE(1)
            raxold = raxoldo + dist1
            zaxold = zaxoldo
            psistab1o = psistabr
            psistab2o = psistabz

        CASE(2)
            raxold = raxoldo
            zaxold = zaxoldo + dist2
            dcrdr = (psistabr - psistab1o)/dist1
            dczdr = (psistabz - psistab2o)/dist1

        CASE(3) ! inverse of (dcrdr dcrdz  ; dczdr dczdz ) = ( dczdz -dcrdz ; -dczdr dcrdr)
            dcrdz = (psistabr - psistab1o)/dist2
            dczdz = (psistabz - psistab2o)/dist2
            det = (dcrdr*dczdz - dcrdz*dczdr)
            cibapr = ( dczdz*psistab1o - dcrdz*psistab2o)/det
            cibazr = (-dczdr*psistab1o + dcrdr*psistab2o)/det
            raxold = raxoldo - cibapr
            zaxold = zaxoldo - cibazr
            raxoldo = raxold
            zaxoldo = zaxold
        END SELECT

! Inner cycle, calculate psi1 and psi2

        psro = 1000.
        pszo = 1000.
        redo_bnd = 1

        do j_iter=1, 10000
            raxtmp = trax
            zaxtmp = tzax
            call solve_gse2d_fbe_full_feqis_1turn(j_iter - 1 + j_iter2 - 1, 1, raxold, zaxold)
            temp_err = (abs(psro - psistabr) + abs(pszo - psistabz))
            if (temp_err <= err_find_psistab) EXIT
            psro = psistabr
            pszo = psistabz
        enddo

        if (j_iter >= 10000) then
            write(*, *) 'total iterations passed, stopping!'
            stop
        endif

        if (j_cyclo == 3) then
            temp_err2 = abs(psistabr) + abs(psistabZ)
            j_cyclo = 1
            write(*, *) n_of_newton_iterations, nint((0. + j_iter2)/3.)
            if (nint((0. + j_iter2)/3.) >= n_of_newton_iterations) EXIT
        else
            j_cyclo = j_cyclo + 1
        endif

        if (j_iter2 >= 400000) EXIT
        if (temp_err2 <= err_find_psistab) EXIT

    enddo

CASE(101) ! refit_mode=101: only vertical stab
!start iterations to find self-consistent solution
    g000(1:nr2, 1:nz2) = psiextrz(1:nr2, 1:nz2)
    iaxis = closest_index(raxp, r(1), dr)
    jaxis = closest_index(zaxp, z(1), dz)
    rax = r(iaxis)
    zax = z(jaxis)
    raxold    = rax
    zaxold    = zax
    raxoldo   = rax
    zaxoldo   = zax
    psistabR  = 0.
    psistabZ  = 0.
    psistab1o = 0.
    psistab2o = 0.
    temp_err  = 100.
    temp_err2 = 100.
    rleft  = raxold
    rright = raxold
    zup    = zaxold
    zdown  = zaxold
! outer cycle, calculate new axis
    j_cyclo = 0
    jeppa   = 1
    dist1 = dr
    dist2 = dz

    do j_iter2=1, 300
        if (j_iter2 > 250) stop

        if (j_cyclo == 1) then
            zaxold = zaxoldo + dist2
            psistab2o = psistabz
        endif
        if (j_cyclo == 2) then      ! inverse of (dcrdr dcrdz  ; dczdr dczdz ) = ( dczdz - dcrdz ; -dczdr dcrdr)
            dczdz = (psistabz - psistab2o)/dist2
            cibazr = psistab2o/dczdz
            zaxold  = zaxoldo - cibazr
            zaxoldo = zaxold
        endif

! Inner cycle, calculate psi1 and psi2
        psro = 1000.
        pszo = 1000.
        redo_bnd = 1

        do j_iter=1, 10000
            raxtmp = trax
            zaxtmp = tzax
            call solve_gse2d_fbe_full_feqis_1turn(j_iter - 1 + j_iter2 - 1, 1, -1.d6, zaxold)
            temp_err = (abs(pszo - psistabz))
            if (temp_err <= err_find_psistab) EXIT
            psro = psistabr
            pszo = psistabz
        enddo

        if (j_iter >= 10000) then
            write(*, *) 'total iterations passed!'
            stop
        endif

        if (j_cyclo == 2) then
            temp_err2 = abs(psistabZ)
            j_cyclo = 1
            write(*, *) n_of_newton_iterations, nint((0. + j_iter2)/2.)
            if (nint((0. + j_iter2)/2.) >= n_of_newton_iterations) EXIT
        else
            j_cyclo = j_cyclo + 1
        endif

        if (j_iter2 >= 400000) EXIT
        if (temp_err2 <= err_find_psistab) EXIT

    enddo

    write(*, *) 'total iterations passed!'

CASE(1) ! refit_mode=1
    call restab_axis_with_fourier_wall

CASE(2)
    call restab_boundary_with_fourier_wall !doesnt work well

CASE(3)
    call restab_F_function_full_fonfit

END SELECT

return
end subroutine solve_gse2d_fbe_full_feqis

!--------------------------------------------------------------------
subroutine solve_gse2d_fbe_full_feqis_1turn(j_init, j_stab, raxold, zaxold)

use feqis_circuit, only : nr2, nz2, iaxis, jaxis, &
    r, z, dr, dz, rax, zax, raxp, zaxp, &
    iplasma, jrz, psirz, psiextrz, psiplasrz, psistabr, psistabz, &
    nine_point_coeffs_only, boundary, interp_j_fromrhotorz, &
    find_new_axis_part1, find_psi_boundary, new_jrz_feqis

use feqis_tools, only: closest_index

implicit none

integer, intent(in) :: j_init, j_stab
double precision, intent(in) :: raxold, zaxold

integer :: i, j
double precision :: curr, dum1, dum2, zum1, zum2, delr, delz
double precision, dimension(9) :: c
double precision, dimension(300, 300) :: g

!first, initialized initial guess coming from prescribed boundary current density: jrhoteta
if (j_init == 0) then
    write(*, *) 'Re-interp density'
    call interp_j_fromrhotorz
! Rescale current density
    curr = SUM(jrz(1:nr2, 1:nz2)) *dr*dz
    jrz = jrz/curr*iplasma
    rax = raxp
    zax = zaxp
    iaxis = closest_index(rax, r(1), dr)
    jaxis = closest_index(zax, z(1), dz)
    rax = r(iaxis)
    zax = z(jaxis)
    write(*, *) raxp, zaxp, rax, zax
endif

g = 0.
call solve_gs2d(g) !jrz as right hand side
g = boundary(g)   ! gbound = integral (Green*dg/dn) over the boundary
call solve_gs2d(g) ! again jrz as right hand side
psiplasrz(1:nr2, 1:nz2) = g(1:nr2, 1:nz2)

if (j_stab == 1) then
    psistabr = 0.
    psistabz = 0.
    delr = 0.
    delz = 0.
    psirz(1:nr2, 1:nz2) = psiplasrz(1:nr2, 1:nz2) + psiextrz(1:nr2, 1:nz2) !total flux
    call find_new_axis_part1
    write(*, *) 'natural ax', rax, zax

    call nine_point_coeffs_only(raxold, zaxold, c, zum1, zum2)
! dpsidr
    zum1 = (raxold - zum1)/dr
    zum2 = (zaxold - zum2)/dz

    dum1 = 2.*c(2)*(zum1*zum2**2 + 2.*c(2)*zum1*zum2) +  &
              c(3)*zum2**2 + c(4)*zum2 + 2*c(5)*zum1 + c(7)
 
! dpsidz
   dum2 = 2.*c(1)*zum1**2*zum2 + c(2)*zum1**2 +  &
          2.*c(3)*zum1*zum2    + c(4)*zum1 + 2.*c(6)*zum2 + c(8)

    psistabr = -1./(2.*raxold)*dum1/dr
    psistabz = -dum2/dz

    do i=1, nr2
        do j=1, nz2
            psirz(i, j) = psiplasrz(i, j) + psiextrz(i, j) + psistabr*r(i)**2 + psistabz*z(j) ! Total flux
        enddo
    enddo
    call find_new_axis_part1
else
    psirz(1:nr2, 1:nz2) = psiplasrz(1:nr2, 1:nz2) + psiextrz(1:nr2, 1:nz2) ! Total flux
endif

call find_new_axis_part1
call find_psi_boundary
call new_jrz_feqis

return
end subroutine solve_gse2d_fbe_full_feqis_1turn

!--------------------------------------------------------------------
subroutine FEQISUPDATE(coilzzz, nccc)
       
use pi_vars, only: GPI2
use feqis_circuit, only: nconduc, cur_con_old, curconduc, &
    psi_cur_old, psiplasmatoconduc
use astra2fbe, only: fast_mode

implicit none

integer, intent(in) :: nccc
double precision, dimension(nccc) :: coilzzz

coilzzz(1:nccc) = curconduc(1:nccc)*1.e3
cur_con_old(1:nconduc) = curconduc(1:nconduc)

if (fast_mode == 0) then
    psi_cur_old(1:nconduc) = psiplasmatoconduc(1:nconduc)
endif

return
end subroutine FEQISUPDATE

!--------------------------------------------------------------------
subroutine circuit_eq_advance_feqis(j_init)

use pi_vars, only: GPI, GPI2
use feqis_circuit, only: nconduc, i_dim1, tau_new, &
    tau_old, cur_con_old, curconduc, voltage, &
    dpc, psiplasmatoconduc, resconduc, indconduc, psi_cur_old
use feqis_tools, only: solve_circuit_equations
use astra2fbe, only: tau_circuit_feqis, tau_gseq_feqis, activate_coil_feqis, &
    n_equivalence, reconnect_circuits, new_equivalence, &
    use_reduce_circuit, current_limit_feqis, force_coil

implicit none

integer, intent(in) :: j_init

integer :: i, ic, j, k, iii, jjj, invertcommand, i_equivalence, i_cnew, firstcall
integer, dimension(200) :: rem_coils
double precision, dimension(i_dim1) :: dpctemp, vtemp, curotemp, curtemp
double precision, dimension(200, 200) :: restemp, indtemp

data firstcall/0/
save restemp, indtemp, rem_coils, i_cnew, firstcall

tau_new = tau_circuit_feqis
invertcommand = 0

if (tau_new /= tau_old) then
    invertcommand = 1
    tau_old = tau_new
endif
if (j_init == 0 .or. firstcall == 0) then
    invertcommand = 1
    tau_old = tau_new
    rem_coils = 0
endif

ic = nconduc

if (j_init == 0) cur_con_old = curconduc
if (j_init == 0) i_cnew = ic
if (j_init == 0) return

firstcall = 1

dpc(1:ic) = GPI2*(psiplasmatoconduc(1:ic) - psi_cur_old(1:ic))/tau_gseq_feqis !plasma contribution

do j=1, nconduc
    if (activate_coil_feqis(j) == 0) cur_con_old(j) = 0.
    if (activate_coil_feqis(j) == 0) dpc(j) = 0.
enddo

if (reconnect_circuits == 1) then
    invertcommand = 1
    i_cnew = ic
    indtemp(1:ic, 1:ic) = indconduc(1:ic, 1:ic)
    restemp(1:ic, 1:ic) = resconduc(1:ic, 1:ic)
    do k=1, n_equivalence
        i_equivalence = 1
        do i=1, nconduc
            if (new_equivalence(i, k) > 0 .and. i_equivalence == 1) then
                i_equivalence = 0
                do j=i + 1, nconduc
                    if (new_equivalence(j, k) == new_equivalence(i, k)) then
                        indtemp(i, i) = indconduc(i, i) + indconduc(j, j) + indconduc(i, j) + indconduc(j, i)
                        restemp(i, i) = resconduc(i, i) + resconduc(j, j) + resconduc(i, j) + resconduc(j, i)
                        do iii=1, nconduc
                            if (iii /= i .and. iii /= j) then
                                indtemp(i, iii) = indconduc(i, iii) + indconduc(j, iii)
                                restemp(i, iii) = resconduc(i, iii) + resconduc(j, iii)
                                indtemp(iii, i) = indconduc(iii, i) + indconduc(iii, j)
                                restemp(iii, i) = resconduc(iii, i) + resconduc(iii, j)
                            endif
                        enddo
                        rem_coils(j) = new_equivalence(j, k)
                        ic = ic - 1
                    endif
                enddo
    ! Remove obsolete columns and rows
                do jjj=1, nconduc
                    j = nconduc - jjj + 1
                    if (rem_coils(j) > 0) then
                        indtemp(1:nconduc, j:nconduc) = indtemp(1:nconduc, j + 1:nconduc + 1)
                        restemp(1:nconduc, j:nconduc) = restemp(1:nconduc, j + 1:nconduc + 1)
                        indtemp(j:nconduc, 1:nconduc) = indtemp(j + 1:nconduc + 1, 1:nconduc)
                        restemp(j:nconduc, 1:nconduc) = restemp(j + 1:nconduc + 1, 1:nconduc)
                    endif
                enddo
            endif
        enddo
    enddo
    i_cnew = ic
endif

if (use_reduce_circuit == 1) then
    ic = nconduc
    dpctemp(1:ic) = dpc(1:ic)
    vtemp(1:ic) = voltage(1:ic)
    curotemp(1:ic) = cur_con_old(1:ic)
    do k=1, n_equivalence
        i_equivalence = 1
        do i=1, nconduc
            if (new_equivalence(i, k) > 0 .and. i_equivalence == 1) then
                i_equivalence = 0
                do jjj=i+1, nconduc
                    j = nconduc - jjj + i + 1
                    if (new_equivalence(j, k) == new_equivalence(i, k)) then
                        dpctemp(i) = dpctemp(i) + dpctemp(j)
                        ic = ic - 1
                    endif
                enddo
            endif
        enddo
    enddo
    do jjj=1, nconduc
        j = nconduc - jjj + 1
        if (rem_coils(j) > 0) then
            dpctemp (j:nconduc) = dpctemp (j + 1:nconduc + 1)
            vtemp   (j:nconduc) = vtemp   (j + 1:nconduc + 1)
            curotemp(j:nconduc) = curotemp(j + 1:nconduc + 1)
        endif
    enddo
endif

i = i_cnew

if (use_reduce_circuit == 0) then
    curconduc(1:i) = solve_circuit_equations(i, indconduc(1:i, 1:i), resconduc(1:i, 1:i), & 
        cur_con_old(1:i), voltage(1:i), dpc(1:i), tau_new, invertcommand)
else
    curtemp(1:i) = solve_circuit_equations(i, indtemp(1:i, 1:i), restemp(1:i, 1:i), & 
        curotemp(1:i), vtemp(1:i), dpctemp(1:i), tau_new, invertcommand)

! Re-adapt currents
    do i=1, nconduc
        if (rem_coils(i) == 0) then
            curconduc(i) = curtemp(i)
        endif
        if (rem_coils(i) > 0) then
            curconduc(i) = curtemp(rem_coils(i))
            curtemp(i + 1:i_cnew + 1) = curtemp(i:i_cnew)
        endif
    enddo
endif

ic = nconduc
i = ic
do j=1, ic
    if (activate_coil_feqis(j) == 0) curconduc(j) = 0.
    curconduc(j) = max(curconduc(j), current_limit_feqis(j, 2))
    curconduc(j) = min(curconduc(j), current_limit_feqis(j, 1))
    if (sum(force_coil(j, 1:i)) > 0.5) then
        do k=1, i
            if (force_coil(j, k) == 1.) then
                curconduc(j) = curconduc(k)
            endif
        enddo
    endif
enddo

reconnect_circuits = 0

return
end subroutine circuit_eq_advance_feqis

!--------------------------------------------------------------------
subroutine feqis_init(equil_in, params, j_call, ifplasma)

use pi_vars, only: GPI, GPI2, mu0
use errors_params, only: err_circ_plasma_iter, err_find_oxpoints_derivs, &
    err_find_psistab, err_find_delr, err_find_biquad, err_epsilon, &
    err_gaptolez, err_fix_boundary, err_find_oxpoints
use parameters_a2equil, only: type_parameters, max_iter, &
    err_circ_in, err_find_oxpoints_in, err_find_psistab_in, err_find_delr_in, &
    err_find_biquad_in, err_epsilon_in, err_gaptolez_in, err_fix_boundary_in, &
    err_find_oxpoints_derivs_in
use imas_ids, only: type_equilibrium
use feqis_circuit, only: nrho, nrho2d, nteta, use_limiter_yesno, &
    dr_factor_init, dz_factor_init, &
    rexp, zexp, raxp, zaxp, rbnd, zbnd, rbndp, zbndp, &
    teta, dteta, tetaexp, &
    iplasma, Rgeom0, Btor0, voltage, voltage_old, omega_pl, &
    pressure, ipol, pprime, ffprime, &
    psia_2d, ffp_2d, ppp_2d, &
    psistabr, psistabz, psigrid, psigrida, psibnd
use astra2fbe, only: dr_factor_init_astra, dz_factor_init_astra, &
    tau_circuit_feqis, tau_gseq_feqis, activate_coil_feqis, current_limit_feqis, &
    raxis_astra, zaxis_astra, psi0_astra, psib_astra, use_limiter_astra
use numerical_tools, only: linterp
use feqis_tools, only: pol_angle

implicit none

integer, intent(in) :: j_call, ifplasma
type(type_parameters ), intent(in) :: params
type(type_equilibrium), intent(in) :: equil_in

integer :: i, j, k
double precision, dimension(700) :: rdum, zdum, tdum

if (j_call == 0) then
    nteta = equil_in%eqgeometry%boundary%npoints
    nrho = params%neql
    psistabR = 0.
    psistabZ = 0.
! Normalized psi from 0 axis to 1 edge, equispaced
    do i=1, nrho
        psigrid(i) = (i - 1.)/(nrho - 1.)
    enddo

    dr_factor_init = dr_factor_init_astra
    dz_factor_init = dz_factor_init_astra

!errors
    err_circ_plasma_iter     = err_circ_in
    err_find_oxpoints        = err_find_oxpoints_in
    err_find_oxpoints_derivs = err_find_oxpoints_derivs_in
    err_find_psistab         = err_find_psistab_in
    err_find_delr            = err_find_delr_in
    err_find_biquad          = err_find_biquad_in
    err_epsilon              = err_epsilon_in
    err_gaptolez             = err_gaptolez_in
    err_fix_boundary         = err_fix_boundary_in

    Rgeom0 = equil_in%global_param%toroid_field%r0
    tau_circuit_feqis = 0.001 !default value
    tau_gseq_feqis = 0.001  !default value
    activate_coil_feqis = 1. ! when 0., coil is forced to 0 current
    current_limit_feqis(:, 1) = 1e6 ! cant be higher than 1e6 MA
    current_limit_feqis(:, 2) = -1e6 ! cant be lower than -1e6 MA
    max_iter = 10000 !hardwired
! teta for polar grid, goes from 0 to 2*pi-dteta, but point nt + 1 is the periodic one
    omega_pl = 0.
    psi0_astra = equil_in%profiles_1d%psi(1)
    psib_astra = equil_in%profiles_1d%psi(nrho)
    raxp = raxis_astra
    zaxp = zaxis_astra
    psibnd = -1.e6
    use_limiter_yesno = use_limiter_astra
    voltage_old = 0.
    voltage = 0.
endif

if (ifplasma == 1) then
    dteta = GPI2/(nteta + 1.)
    do i=1, nteta + 1
        teta(i) = GPI2*(i - 1.)/(nteta + 0.)
    enddo

    btor0   = equil_in%global_param%toroid_field%b0
    iplasma = equil_in%global_param%i_plasma/1.e6
    pressure(1:nrho) = equil_in%profiles_1d%pressure(1:nrho)
    pprime(1:nrho)   = equil_in%profiles_1d%pprime(1:nrho)
    ffprime(1:nrho)  = equil_in%profiles_1d%ffprime(1:nrho)
    psigrida(1:nrho) = equil_in%profiles_1d%psi(1:nrho) !unnormalized
    psigrida(1:nrho) = (psigrida(1:nrho) - psigrida(1))/(psigrida(nrho) - psigrida(1)) !normalized 0 axis, 1 boundary

    do i=1, nrho2d
        psia_2d(i) = (i - 1.)/(nrho2d - 1.)
    enddo
    call linterp(psigrida(1:nrho), ffprime(1:nrho), nrho, psia_2d, ffp_2d, nrho2d)
    call linterp(psigrida(1:nrho), pprime (1:nrho), nrho, psia_2d, ppp_2d, nrho2d)
    ffp_2d = -GPI2/mu0*ffp_2d
    ppp_2d = -GPI2*1.e-6*ppp_2d
    ipol(1:nrho) = equil_in%profiles_1d%F_dia(1:nrho)

    if (params%k_fixfree == 0) then ! If 1, comes from free boundary
        rexp(1:nteta) = equil_in%eqgeometry%boundary%r(1:nteta)
        zexp(1:nteta) = equil_in%eqgeometry%boundary%z(1:nteta)
! Define angle not based on mag axis, but on geometrical center
        if (j_call == 0) then
            raxp = raxis_astra
            zaxp = zaxis_astra
        endif

        do i=1, nteta
            tetaexp(i) = pol_angle(raxp, zaxp, rexp(i), zexp(i))
        enddo

! Order points
        rdum(1:nteta) = rexp(1:nteta)
        zdum(1:nteta) = zexp(1:nteta)
        tdum(1:nteta) = tetaexp(1:nteta)
        do i=1, nteta
            j = minloc(tdum(1:nteta), 1)
            rexp(i) = rdum(j)
            zexp(i) = zdum(j)
            tetaexp(i) = tdum(j)
            tdum(j) = 1.e6
        enddo

! Reorder
        if (tetaexp(1) > tetaexp(nteta)) then ! Reorder
            j = 1
            do k=1, nteta-1
                if (tetaexp(k + 1) < tetaexp(k)) j = k + 1   ! j is the first teta above 0
            enddo
            if (j > 1) tetaexp(1:j-1) = tetaexp(1:j-1) - GPI2
        endif

        do i=1, nteta
            rexp(nteta + i) = rexp(i)
            zexp(nteta + i) = zexp(i)
            tetaexp(nteta + i) = tetaexp(i) + GPI2
        enddo
        call linterp(tetaexp(1:nteta*2), rexp(1:nteta*2), nteta*2, teta(1:nteta), rbndp(1:nteta), nteta)
        call linterp(tetaexp(1:nteta*2), zexp(1:nteta*2), nteta*2, teta(1:nteta), zbndp(1:nteta), nteta)
        rbnd(1:nteta) = rbndp(1:nteta)
        zbnd(1:nteta) = zbndp(1:nteta)
    endif
endif

return
end subroutine feqis_init

!--------------------------------------------------------------------
subroutine equil_feqis_init_circ

use pi_vars, only: GPI
use fft_mod_eff, only: sintable, costable
use feqis_circuit, only: nr, nr1, nr2, nz, nz1, nz2, &
    nactive, npassive, ncoils, nconduc, nlimiter, nblocks, ngbnd, &
    lim_minr, lim_maxr, lim_minz, lim_maxz, &
    rmin, rmax, zmin, zmax, r, z, dr, dz, rcomp, zcomp, r_cond, z_cond, &
    rcoil, zcoil, drcoil, dzcoil, anglecoil, anglehcoil, mequivalence, &
    limiterr, limiterz, alpsep, curconduc, resconduc, indconduc, &
    zlimpotential, green_bnd_f
use green_matrix, only: greeni, dgreenirj, dgreenizj, dgreenirpl, dgreenizpl
use outcmn_inc, only: machine
use astra2fbe, only: cur_init

implicit none

integer :: i, j, ii, jj
character(len=80) :: fname

fname = 'exp/cnf/machine_description_out.'//trim(machine)
open(32, file=TRIM(fname))
    read(32, *) nr2, nz2
    read(32, *) rmin
    read(32, *) rmax
    read(32, *) zmin
    read(32, *) zmax
    read(32, *) alpsep
    nr1 = nr2 - 1
    nz1 = nz2 - 1
    nr  = nr1 - 1
    nz  = nz1 - 1
    do i=1, nr2
        r(i) = rmin + (i - 1.)*(rmax - rmin)/nr1     ! computational domain is r(2:nr + 1), boundaries are r(1) and r(nr + 2)
    enddo
    do i=1, nz2
        z(i) = zmin + (i - 1.)*(zmax - zmin)/nz1
    enddo
    rcomp(1:nr) = r(2:nr1)
    zcomp(1:nz) = z(2:nz1)
    dr = r(2) - r(1)
    dz = z(2) - z(1)

    do i=1, nz
        do j=1, nz
            sintable(i, j) = sin(i*j*GPI/(nz + 1))
            costable(i, j) = cos(i*j*GPI/(nz + 1))
        enddo
    enddo

! Load everything from file
    read(32, *) nactive, npassive
    read(32, *) ncoils
    do i=1, ncoils
        read(32, *) rcoil(i), zcoil(i), drcoil(i), dzcoil(i), anglehcoil(i), anglecoil(i), mequivalence(i)
    enddo
    read(32, *) nlimiter
    do i=1, nlimiter
        read(32, *) limiterr(i), limiterz(i)
    enddo
    read(32, *) lim_maxR
    read(32, *) lim_minR
    read(32, *) lim_maxZ
    read(32, *) lim_minZ
    do i=nactive + 1, npassive
        read(32, *) r_cond(i), z_cond(i)
    enddo
    read(32, *) nconduc
    do i=1, nconduc
        read(32, *) indconduc(i, 1:nconduc)
    enddo
    read(32, *) nconduc
    do i=1, nconduc
        read(32, *) resconduc(i, 1:nconduc)
    enddo
    do i=1, nconduc
        do j=1, nr2
            read(32, *) greeni(j, 1:nz2, i)
        enddo
    enddo
    read(32, *) nblocks
    do j=1, nblocks
        do i=1, nblocks
            read(32, *) dgreenirj(i, j), dgreenizj(i, j)
        enddo
    enddo
    do ii=1, nblocks
        do j=1, nz2
            do i=1, nr2
                read(32, *) dgreenirpl(i, j, ii), dgreenizpl(i, j, ii)
            enddo
        enddo
    enddo
    do jj=1, nz2
        do ii=1, nr2
            read(32, *) zlimpotential(ii, jj)
        enddo
    enddo
    read(32, *) ngbnd
    read(32, *) green_bnd_f(1:ngbnd)
close(32)

write(*, *) nactive 

! assign initial currents from astra
curconduc(1:nconduc) = cur_init(1:nconduc)

return
end subroutine equil_feqis_init_circ

!--------------------------------------------------------------------
subroutine fix_boundary_feqis(j_init)

use feqis_circuit, only: nr, nrho, nteta, raxp, zaxp, rbndp, zbndp, rho, teta, &
    psiaxisp, psirhoteta, psigrida, psplex, psibndp, &
    psia_1d, ffp_1d, ppp_1d, &
    ffprime, pprime, pressure, ipol, &
    Rgeom0, Btor0, Rpol, Zpol, Rpul, Zpul, jrhoteta, li3, betapol, iplasma
use astra2fbe, only: raxis_astra, zaxis_astra, psi0_astra, psib_astra, &
    solve_fix
use metric_coefficients_pbe, only: lambda2d, R_curr_0d, Z_curr_0d, dator
use transfer_functions, only: rpbez, zpbez, psibez, t2dbez, &
    g1bez, g2bez, gm1bez, gm4bez, gm41bez, gm5bez, ggrhobez, &
    areatbez, surfbez, perimbez, vbez, qbez, phibez, &
    bmaxbez, bminbez, bdb0bez, fofbbez, bcell2dbez, bpcell2dbez, &
    ffprimebez, pprimebez, pressbez, ipolbez, rinbez, routbez, &
    kbez, triaubez, shifbez, rbp2_b2bez, rmin2dbez, dpsidvbez, &
    jrhobez, shivbez, squarebez
use pi_vars, only: GPI, GPI2, GPI4, muvac

implicit none

integer, intent(in) :: j_init
integer :: i, j, jr, jt, ierr
double precision :: psiaxis_new, cnorm, rax_new, zax_new, rhoedge
double precision, dimension(nrho) :: q_new, effprimp, epprimp
double precision, dimension(nteta) :: thetap_i
double precision, dimension(512, 512) :: psisave
double precision, dimension(nrho, nteta) :: rmaj2, jcbn2, darea2, &
    r_min, yy2, jrho2, gradr2, darea

data ierr/0/
save psisave, ierr


! initial guess
if (j_init == 0) then
    raxp = raxis_astra
    zaxp = zaxis_astra
    psiaxisp = psi0_astra
    psibndp = psib_astra
    do jr=1, Nrho
        do jt=1, Nteta
            lambda2d(jr, jt) = (jr - 1)/(Nr - 1.)
        enddo
    enddo 
    do jt=1, Nteta
         psirhoteta(1: Nrho, jt) = psigrida(1: Nrho)
    enddo
else
    psirhoteta(1:nrho, 1:nteta) = psisave(1:nrho, 1:nteta)
    psiaxisp = psirhoteta(   1, 1)
    psibndp  = psirhoteta(nrho, 1)
endif

!boundary from previous time step

call PHI_EQ_2d_PBE(nrho, nteta, psigrida(1:nrho), iplasma, &
    ffprime(1:nrho), pprime(1:nrho), btor0*rgeom0, &
    rbndp(1:nteta), zbndp(1:nteta), Raxp, Zaxp, psiaxisp, psibndp, &
    solve_fix, j_init, rpbez(1:nrho, 1:nteta), zpbez(1:nrho, 1:nteta), & 
    psirhoteta(1:nrho, 1:nteta), psibez(1:nrho), &   ! psinorm new
    lambda2d(1:nrho, 1:nteta), t2dbez(1:nteta), &
    psiaxis_new, cnorm, rax_new, zax_new, thetap_i, rmaj2, & 
    jcbn2, q_new, rhoedge, darea2, epprimp, effprimp, r_min, yy2, gradr2, darea, ierr)
 
raxp = rax_new
zaxp = zax_new
psiaxisp = psiaxis_new
psisave(1:nrho, 1:nteta) = psirhoteta(1:nrho, 1:nteta)

! notice that qedge is the extrapolate of q_new at neql
ffprimebez(1:nrho) = effprimp(1:nrho)/cnorm
pprimebez (1:nrho) = epprimp (1:nrho)/cnorm
pressbez  (1:nrho) = pressure(1:nrho)
ipolbez   (1:nrho) = ipol    (1:nrho)

!additional calculations

do jt=1, nteta
    do jr=1, nrho
        jrhoteta(jr, jt) = (effprimp(jr)/rpbez(jr, jt) + &
            rpbez(jr, jt)*epprimp(jr))/GPI2/0.4/GPI/cnorm
    enddo
enddo

do jt=1, nteta
    do jr=1, nrho-1
        jrho2(jr, jt) = ( 0.5*(effprimp(jr ) + effprimp(jr+1)) * 1./Rmaj2(jr, jt) + &
            Rmaj2(jr, jt)*0.5*(epprimp(jr+1) + epprimp(jr)))/GPI2/0.4/GPI/cnorm
    enddo
enddo

jr = nrho - 1 !nint(1.*nx)
Z_curr_0D = sum(jrho2(1:jr, :)*YY2(1:jr, :)*darea2(1:jr, :))/ &
            sum(jrho2(1:jr, :)*darea2(1:jr, :))
R_curr_0D = sqrt(sum(jrho2(1:jr, :)*Rmaj2(1:jr, :)**2*darea2(1:jr, :))/ &
                 sum(jrho2(1:jr, :)*darea2(1:jr, :)))

!regrid
call build_2dgrid(nrho, nteta, psibez(1:nrho), &
    psirhoteta(1:nrho, 1:nteta), rgeom0, pressure(1:nrho), &
    btor0, ipol(1:nrho), iplasma, rpbez(1:nrho, 1:nteta), zpbez(1:nrho, 1:nteta), &
    rmaj2(1:nrho, 1:nteta), r_min(1:nrho, 1:nteta), jcbn2(1:nrho, 1:nteta), & 
    thetap_i(1:nteta), q_new(1:nrho), rhoedge, gradr2(1:nrho, 1:nteta), &
    g2bez(1:nrho), gm1bez(1:nrho), areatbez(1:nrho), perimbez(1:nrho), vbez(1:nrho), &
    g1bez(1:nrho), ggrhobez(1:nrho), bmaxbez(1:nrho), bminbez(1:nrho), &
    gm4bez(1:nrho), bdb0bez(1:nrho), gm5bez(1:nrho), fofbbez(1:nrho), surfbez(1:nrho), & ! lateral surface
    li3, betapol, psplex, &
    bpcell2dbez(1:nrho, 1:nteta), bcell2dbez(1:nrho, 1:nteta), &
    routbez(1:nrho), rinbez(1:nrho), kbez(1:nrho), triaubez(1:nrho), shifbez(1:nrho), &
    gm41bez(1:nrho), qbez(1:nrho), shivbez(1:nrho), squarebez(1:nrho)) 

phibez(1:nrho) = 0.
rbp2_b2bez(1:nrho) = 0.

rmin2dbez(1:nrho, 1:nteta) = r_min   (1:nrho, 1:nteta)
dator    (1:nrho, 1:nteta) = darea   (1:nrho, 1:nteta)
jrhobez  (1:nrho, 1:nteta) = jrhoteta(1:nrho, 1:nteta)

rpol(1:nrho, 1:nteta) = rpbez(1:nrho, 1:nteta)
zpol(1:nrho, 1:nteta) = zpbez(1:nrho, 1:nteta)
rpul(1:nrho, 1:nteta) = rpbez(1:nrho, 1:nteta)
zpul(1:nrho, 1:nteta) = zpbez(1:nrho, 1:nteta)

gm41bez  (1:nrho) = 0.
dpsidvbez(1:nrho) = 0.

! additional info from rectangular grid
do j=1, nteta
    do i=1, nrho
        rho(i, j) = sqrt((rpol(i, j) - raxp)**2 + (zpol(i, j) - zaxp)**2)
    enddo
enddo
teta(1:nteta) = t2dbez(1:nteta)
rho(1:nrho, nteta + 1) = rho(1:nrho, 1)
teta(nteta + 1) = teta(1) + GPI2

psia_1d(1:nrho) = psigrida(1:nrho)
ffp_1d (1:nrho) = ffprime (1:nrho)
ppp_1d (1:nrho) = pprime  (1:nrho)

return
end subroutine fix_boundary_feqis

!--------------------------------------------------------------------
subroutine equil_assignments(equil_out)

use pi_vars, only: GPI, GPI2
use imas_ids, only: type_equilibrium
use feqis_circuit, only: nr2, nz2, nrho, nteta, &
    psplex, psiaxis, psibnd, psirhoteta, psirz, &
    r, z, betapol, li3, iplasma
use transfer_functions, only: rpbez, zpbez, t2dbez, &
    rinbez, routbez, rmin2dbez, vbez, areatbez, perimbez, surfbez, &
    kbez, shifbez, triaubez, qbez, phibez, &
    g1bez, g2bez, g2ibez, gm1bez, gm4bez, gm41bez, gm5bez, ggrhobez, &
    bcell2dbez, bpcell2dbez, bminbez, bmaxbez, bdb0bez, fofbbez, &
    psibez, dpsidvbez, rbp2_b2bez, &
    ffprimebez, pprimebez, pressbez, ipolbez, jrhobez, shivbez, squarebez
use metric_coefficients_pbe, only: dator

implicit none

type(type_equilibrium), intent(inout) :: equil_out

equil_out%global_param%psplex   = psplex
equil_out%global_param%psibound = -GPI2*psibnd
equil_out%global_param%psiaxis  = -GPI2*psiaxis
equil_out%global_param%li3      = li3
equil_out%global_param%betpol   = betapol
equil_out%global_param%i_plasma = iplasma*1.e6

equil_out%coord_sys%position%teta2d(1:nteta) = t2dbez(1:nteta)
equil_out%coord_sys%position%psirz(1:nrho, 1:nteta) = psirhoteta (1:nrho, 1:nteta)/GPI2
equil_out%coord_sys%position%r    (1:nrho, 1:nteta) = rpbez      (1:nrho, 1:nteta)
equil_out%coord_sys%position%z    (1:nrho, 1:nteta) = zpbez      (1:nrho, 1:nteta)
equil_out%coord_sys%position%rmin (1:nrho, 1:nteta) = rmin2dbez  (1:nrho, 1:nteta)
equil_out%coord_sys%bpcell        (1:nrho, 1:nteta) = bpcell2dbez(1:nrho, 1:nteta)
equil_out%coord_sys%bcell         (1:nrho, 1:nteta) = bcell2dbez (1:nrho, 1:nteta)
equil_out%coord_sys%darea         (1:nrho, 1:nteta) = dator(1:nrho, 1:nteta)
equil_out%coord_sys%jphi          (1:nrho, 1:nteta) = jrhobez(1:nrho, 1:nteta)

equil_out%eqgeometry%rectgrid%npointsr = nr2
equil_out%eqgeometry%rectgrid%npointsz = nz2
equil_out%eqgeometry%rectgrid%r2d(1:nr2) = r(1:nr2)
equil_out%eqgeometry%rectgrid%z2d(1:nz2) = z(1:nz2)
equil_out%eqgeometry%rectgrid%psirz2d(1:nr2, 1:nz2) = psirz(1:nr2, 1:nz2)
equil_out%eqgeometry%rectgrid%psi_axis     = psiaxis
equil_out%eqgeometry%rectgrid%psi_boundary = psibnd

equil_out%profiles_1d%ffprime   (1:nrho) = ffprimebez(1:nrho)
equil_out%profiles_1d%pprime    (1:nrho) = pprimebez(1:nrho)
equil_out%profiles_1d%pressure  (1:nrho) = pressbez (1:nrho)
equil_out%profiles_1d%rho_tor   (1:nrho) = 0.
equil_out%profiles_1d%F_dia     (1:nrho) = ipolbez(1:nrho)
equil_out%profiles_1d%dPSIdV    (1:nrho) = dpsidvbez(1:nrho)
equil_out%profiles_1d%psi       (1:nrho) = psibez(1:nrho) ! psi nnormalized from 0 to 1 (rhopol^2)
equil_out%profiles_1d%g2        (1:nrho) = g2bez(1:nrho)
equil_out%profiles_1d%g2int     (1:nrho) = g2ibez(1:nrho)
equil_out%profiles_1d%gm1       (1:nrho) = gm1bez(1:nrho)
equil_out%profiles_1d%r_outboard(1:nrho) = routbez(1:nrho)
equil_out%profiles_1d%r_inboard (1:nrho) = rinbez(1:nrho)
equil_out%profiles_1d%volume    (1:nrho) = vbez(1:nrho)
equil_out%profiles_1d%g1        (1:nrho) = g1bez(1:nrho)
equil_out%profiles_1d%gm41      (1:nrho) = gm41bez(1:nrho)
equil_out%profiles_1d%ggradro   (1:nrho) = ggrhobez(1:nrho)
equil_out%profiles_1d%bmaxt     (1:nrho) = bmaxbez(1:nrho)
equil_out%profiles_1d%bmint     (1:nrho) = bminbez(1:nrho)
equil_out%profiles_1d%gm4       (1:nrho) = gm4bez(1:nrho)
equil_out%profiles_1d%bdb0      (1:nrho) = bdb0bez(1:nrho)
equil_out%profiles_1d%gm5       (1:nrho) = gm5bez(1:nrho)
equil_out%profiles_1d%fofb      (1:nrho) = fofbbez(1:nrho)
equil_out%profiles_1d%areat     (1:nrho) = areatbez(1:nrho)
equil_out%profiles_1d%perim     (1:nrho) = perimbez(1:nrho)
equil_out%profiles_1d%shif      (1:nrho) = shifbez(1:nrho)
equil_out%profiles_1d%shiv      (1:nrho) = shivbez(1:nrho)
equil_out%profiles_1d%squareness(1:nrho) = squarebez(1:nrho)
equil_out%profiles_1d%elongation(1:nrho) = kbez(1:nrho)
equil_out%profiles_1d%surface   (1:nrho) = surfbez(1:nrho) ! lateral surface
equil_out%profiles_1d%tria_upper(1:nrho) = triaubez(1:nrho)
equil_out%profiles_1d%tria_lower(1:nrho) = triaubez(1:nrho)
equil_out%profiles_1d%phi       (1:nrho) = phibez(1:nrho)
equil_out%profiles_1d%q         (1:nrho) = qbez(1:nrho)
equil_out%profiles_1d%rbp_b2    (1:nrho) = rbp2_b2bez(1:nrho)

return
end subroutine equil_assignments

!--------------------------------------------------------------------
subroutine convert_boundary_to_pbe

use pi_vars, only: GPI2
use feqis_circuit, only: nr, nr2, nz, nteta, nbnd, i_dim5, iaxis, jaxis, &
    teta, dteta, raus, rinner, zbot, ztop, &
    r, z, dr, dz, rax, zax, raxp, zaxp, rbnd, zbnd, rbndp, zbndp, &
    psiaxis, psibnd, psiaxisp, psibndp, psirz
use feqis_tools, only: pol_angle, interp2d_psi

implicit none

integer :: i, j, k, j4
double precision :: x1, x2, t1, t2, t3, z1, z2, z3, x11, dx
double precision, dimension(i_dim5) :: teta_fbe

do i=1, nteta + 1
    teta(i) = GPI2*(i - 1.)/(nteta + 0.)
enddo
dteta = teta(2) - teta(1)

dx = sqrt(dr**2 + dz**2)

! Find boundary
j = jaxis
do i=iaxis, nr2
    if (psirz(i, j) >= psibnd) k = i
    if (psirz(i, j) <= psibnd) EXIT
enddo

if (psirz(k, j) == psibnd) then
    rbnd(1) = r(k)
else
    rbnd(1) = r(k) - (psirz(k, j) - psibnd)/(psirz(k, j) - psirz(k-1, j))*dr
endif
    zbnd(1) = z(j)
teta_fbe(1) = pol_angle(rax, zax, rbnd(1), zbnd(1))

theta_loop: do i=2, nteta
    dx = sqrt(dr**2 + dz**2)
    x1 = sqrt((rbnd(i-1) - rax)**2 + (zbnd(i-1) - zax)**2)
    teta_fbe(i) = teta_fbe(i-1) + dteta

    t1 = rax + x1*cos(teta_fbe(i))
    t2 = zax + x1*sin(teta_fbe(i))

    if (t1 >= raus) then
        x1 = (raus - rax)/cos(teta_fbe(i))
    endif
    if (t1 <= rinner) then
        x1 = (rinner - rax)/cos(teta_fbe(i))
    endif
    if (t2 >= ztop) then
        x1 = (ztop - zax)/sin(teta_fbe(i))
    endif
    if (t2 <= zbot) then
        x1 = (zbot - zax)/sin(teta_fbe(i))
    endif
    t1 = rax + x1*cos(teta_fbe(i))
    t2 = zax + x1*sin(teta_fbe(i))
    t3 = interp2d_psi(t1, t2, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))

    if (t3 == psibnd) then
        rbnd(i) = t1
        zbnd(i) = t2
    elseif (t3 < psibnd) then
        do
            z1 = rax + (x1 - dx)*cos(teta_fbe(i))
            z2 = zax + (x1 - dx)*sin(teta_fbe(i))
            z3 = interp2d_psi(z1, z2, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
            if (z3 < psibnd) then
                dx = 1.1*dx
            else
                EXIT
            endif
        enddo
        x11 = x1 - (t3 - psibnd)/(t3 - z3)*dx
        rbnd(i) = rax + x11*cos(teta_fbe(i))
        zbnd(i) = zax + x11*sin(teta_fbe(i))
    else
        do
            j4 = 0
            x2 = x1 + dx
            z1 = rax + x2*cos(teta_fbe(i))
            z2 = zax + x2*sin(teta_fbe(i))
            if (z1 >= raus) then
                x2 = (raus - rax)/cos(teta_fbe(i))
                j4 = 1
            endif
            if (z1 <= rinner) then
                x2 = (rinner - rax)/cos(teta_fbe(i))
                j4 = 1
            endif
            if (z2 >= ztop) then
                x2 = (ztop - zax)/sin(teta_fbe(i))
                j4 = 1
            endif
            if (z2 <= zbot) then
                x2 = (zbot - zax)/sin(teta_fbe(i))
                j4 = 1
            endif
            z1 = rax + x2*cos(teta_fbe(i))
            z2 = zax + x2*sin(teta_fbe(i))
            z3 = interp2d_psi(z1, z2, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
            if (z3 > psibnd) then
                if (j4 == 1) then
                    rbnd(i) = z1
                    zbnd(i) = z2
                    CYCLE theta_loop
                endif
                dx = dx*1.1
            else
                EXIT
            endif
        enddo
        x11 = x1 + (t3 - psibnd)/(t3 - z3)*dx
        rbnd(i) = rax + x11*cos(teta_fbe(i))
        zbnd(i) = zax + x11*sin(teta_fbe(i))
   endif

enddo theta_loop

nbnd = nteta

teta (1:nteta) = teta_fbe(1:nteta)
rbndp(1:nteta) = rbnd(1:nteta)
zbndp(1:nteta) = zbnd(1:nteta)
teta (nteta + 1) = teta(nteta) + dteta
rbndp(nteta + 1) = rbnd(nteta) !shouldnt this be rbnd(1) ???
zbndp(nteta + 1) = zbnd(nteta)

raxp = rax
zaxp = zax
psibndp  = psibnd
psiaxisp = psiaxis

return
end subroutine convert_boundary_to_pbe

!--------------------------------------------------------------------
subroutine solve_gs2d(g)

use pi_vars, only: mu0
use feqis_circuit, only: nr, nr1, nr2, nz, nz1, nz2, i_dim2, &
    r, dr, dz, rcomp, jrz
use fft_mod_eff, only: costable
use feqis_tools, only: discrete_sine_transform, solve_tridiag_fbe

double precision, intent(inout), dimension(i_dim2, i_dim2) :: g

integer :: i, j, k, j_init
double precision :: x1, x2, r1m_1, r2m_1
double precision, dimension(258) :: A, B, C, z_fourier
double precision, dimension(i_dim2, i_dim2) :: gt, rhs, wrhs

data j_init/0/
save A, B, C, j_init, z_fourier

do j=1, nz2
    do i=1, nr2
        rhs(i, j) = -mu0*r(i)*jrz(i, j)
    enddo
enddo

r1m_1 = r(  2)/dr**2/((r(  1) + r(  2))/2.)
r2m_1 = r(nr1)/dr**2/((r(nr1) + r(nr2))/2.)

rhs(2:nr1,   2) = rhs(2:nr1,   2) - g(2:nr1,   1)/dz**2
rhs(2:nr1, nz1) = rhs(2:nr1, nz1) - g(2:nr1, nz2)/dz**2

rhs(  2, 2:nz1) = rhs(  2, 2:nz1) - g(  1, 2:nz1)*r1m_1
rhs(nr1, 2:nz1) = rhs(nr1, 2:nz1) - g(nr2, 2:nz1)*r2m_1

wrhs = rhs
! CALL CPU_TIME(tin)
do i=2, nr1
    wrhs(i, 2:nz1) = discrete_sine_transform(nz, wrhs(i, 2:nz1))
enddo

! Create inverse matrix for gs2d
if (j_init == 0) then
    A = 0.
    B = 0.
    C = 0.
    do i=2, nr-1
        x1 = 0.5*(rcomp(  i) + rcomp(i-1))
        x2 = 0.5*(rcomp(i+1) + rcomp(i  ))
        B(i) = -rcomp(i)/dr**2*(1./x2 + 1./x1) 
        A(i) =  rcomp(i)/x2/dr**2 
        C(i) =  rcomp(i)/x1/dr**2 
    enddo
    i = 1
    x1 = 0.5*(rcomp(  i) + r(1))
    x2 = 0.5*(rcomp(i+1) + rcomp(i))
    B(i) = -rcomp(i)/dr**2 * (1./x2 + 1./x1) 
    A(i) =  rcomp(i)/x2/dr**2 
    i = nr
    x1 = 0.5*(rcomp(i) + rcomp(i-1))
    x2 = 0.5*(r(nr2)   + rcomp(i)  )
    B(i) = -rcomp(i)/dr**2 * (1./x2 + 1./x1) 
    C(i) =  rcomp(i)/x1/dr**2
    j_init = 1
    do k=2, nz1
        z_fourier(k) = 2./dz**2 * (costable(1, k-1) - 1.)
    enddo
endif

gt = 0.
!solve matrix
do k=2, nz1
    gt(2:nr1, k) = solve_tridiag_fbe(C(1:nr), B(1:nr) + z_fourier(k), A(1:nr), wrhs(2:nr1, k), nr)
enddo
!invert fourier from gt(1:nr, 1:kfourier) to g(2:nr1, 2:nz1)
!  gt(i, k)=sum(invMM_gs2d(i-1, 1:nr, k-1)*wrhs(2:nr1, k))

do i=2, nr1
    gt(i, 2:nz1) = discrete_sine_transform(nz, gt(i, 2:nz1))
enddo
g(2:nr1, 2:nz1) = 2./(nz + 1)*gt(2:nr1, 2:nz1)

return
end subroutine solve_gs2d

!---------------------------------------------------------------------
subroutine coil_forces_feqis(ncoilz, force_R, force_Z, plasma_state)

use feqis_circuit, only: nblocks, npassive, jrz, nr2, nz2, area_eff, &
    curconduc, mequivalence
use green_matrix, only: dgreenirpl, dgreenizpl, dgreenirj, dgreenizj

integer, intent(in) :: ncoilz, plasma_state
double precision, intent(out), dimension(ncoilz) :: force_R, force_Z

integer :: i, j, k, nblock_a
double precision :: x1

force_R = 0.
force_Z = 0.
nblock_a = nblocks - npassive

if (plasma_state == 1) then !not sure about the plasma response...
    do i=1, nblock_a
        x1 =  sum(jrz(1:nr2, 1:nz2) * area_eff(1:nr2, 1:nz2) * dgreeniRpl(1:nr2, 1:nz2, i))
        force_R(i) = force_R(i) + curconduc(mequivalence(i)) * x1
        x1 = -sum(jrz(1:nr2, 1:nz2) * area_eff(1:nr2, 1:nz2) * dgreeniZpl(1:nr2, 1:nz2, i))
        force_Z(i) = force_Z(i) + curconduc(mequivalence(i)) * x1
    enddo
endif

!block-to-block
do i=1, nblock_a
    do j=1, nblock_a
        if (i /= j) then
            force_R(i) = force_R(i) + curconduc(mequivalence(j)) * curconduc(mequivalence(i)) * dgreeniRj(i, j)
            force_Z(i) = force_Z(i) - curconduc(mequivalence(j)) * curconduc(mequivalence(i)) * dgreeniZj(i, j)
        endif
    enddo
enddo

return
end subroutine coil_forces_feqis

