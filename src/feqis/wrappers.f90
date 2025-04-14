subroutine full_system_advance_feqis(j_init, no_circuit_eq)

use errors_params, only: err_epsilon, err_circ_plasma_iter
use feqis_circuit, only: psi_cur_old, &
    psiplasmatoconduc
use fbe_core, only: nr2, nz2, nconduc, curconduc, jrz, area_eff, psi_external_calc
use global_params, only: iplasma
use transport2fbe, only: fast_mode, execute_plasma
use parameters_a2equil, only: max_iter
use green_function, only: greeni

implicit none

integer, intent(in) :: no_circuit_eq
integer, intent(inout) :: j_init

integer :: i, j_iter
double precision :: error_temp
double precision, dimension(300) :: cur_temp

if (j_init == -818) then !!!refit mode 818, special case
    call psi_external_calc
    call solve_gse2d_fbe_full_feqis(0)
    return
endif

! time stepping
! at iteration 0, dpsidt = 0
if (j_init == 0 .or. no_circuit_eq == 1) then
! First do full equilibrium solution at time t=0
    call psi_external_calc
    call solve_gse2d_fbe_full_feqis(0)
    do i=1, nconduc
        psiplasmatoconduc(i) = sum(jrz(1: nr2, 1: nz2) * area_eff(1: nr2, 1: nz2) * greeni(1: nr2, 1: nz2, i))
    enddo
    psi_cur_old(1:nconduc) = psiplasmatoconduc(1:nconduc)

    call circuit_eq_advance_feqis(0)
    j_init = -1
    return
endif

if (fast_mode == 1 .and. execute_plasma == 1) then
    psi_cur_old(1:nconduc) = psiplasmatoconduc(1:nconduc)
    call psi_external_calc
    call solve_gse2d_fbe_full_feqis_1turn(1, 0, 0.d0, 0.d0)
    do i=1, nconduc
        psiplasmatoconduc(i) = sum(jrz(1: nr2, 1: nz2) * area_eff(1: nr2, 1: nz2) * greeni(1: nr2, 1: nz2, i))
    enddo
endif

do j_iter=1, max_iter

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
        return
    endif

enddo

return
end subroutine full_system_advance_feqis

!--------------------------------------------------------------------
subroutine solve_gse2d_fbe_full_feqis(j_init)

use errors_params, only: err_find_psistab
use fbe_core, only: nr, nz, nr2, nz2, psiextrz, redo_bnd, &
    Rrect, Zrect, dr, dz, dr_factor_init, dz_factor_init, rax, zax, &
    trax, tzax, iaxis, jaxis, psistabr, psistabz, jrz, solve_fbe_static_iterations_curgiven
use pbe_core, only: raxp, zaxp
use feqis_circuit, only: restab_axis_with_fourier_wall, restab_boundary_with_fourier_wall, & !doesnt work well
    restab_F_function_full_fonfit, restab_F_function_full_fonfit_xpoints, restab_F_function_full_currents, &
    restab_2_timepoints_evolution, restab_F_function_full_currents_forces, &
    restab_F_function_full_currents_limits,restab_2_timepoints_evolution_limits, &
    restab_j_timepoints_evolution_limits_xpoints_boundariz, & 
    restab_1_timepoint_limits_xpoints_boundariz, interp_j_fromrhotorz
use transport2fbe, only: refit_mode, n_of_newton_iterations, use_isoflux
use feqis_tools, only: closest_index
use global_params, only: iplasma


implicit none

integer, intent(in) :: j_init

integer :: j_iter, j_iter2, j_cyclo, jeppa
double precision :: temp_err, raxold, zaxold, temp_err2, raxoldo, zaxoldo, &
    raxtmp, zaxtmp, det, psistab1o, psistab2o, psro, pszo, dist1, dist2, &
    cibapr, cibazr, rleft, rright, zup, zdown, dcrdr, dcrdz, dczdr, dczdz, curr

! First, initialized initial guess coming from prescribed boundary current density: jrhoteta

SELECT CASE(refit_mode)

CASE(-1) ! 1 turn only, linearized solution
    call solve_gse2d_fbe_full_feqis_1turn(j_init, 0, 0.d0, 0.d0)

CASE(0) ! full static convergent solution with given currents
! Start iterations to find self-consistent solution
! first, initialized initial guess coming from prescribed boundary current density: jrhoteta
    if (j_init == 0) then
        call interp_j_fromrhotorz
! Rescale current density
        curr = SUM(jrz(1:nr2, 1:nz2)) *dr*dz
        jrz = jrz/curr*iplasma
        rax = raxp
        zax = zaxp
        iaxis = closest_index(rax, Rrect(1), dr)
        jaxis = closest_index(zax, Zrect(1), dz)
        rax = Rrect(iaxis)
        zax = Zrect(jaxis)
    endif
    call solve_fbe_static_iterations_curgiven(j_init, raxp, zaxp, n_of_newton_iterations)

CASE(1) ! refit eddy currents using fourier method for axis stability. doesnt respect boundary. fixed given active currents
    call restab_axis_with_fourier_wall

CASE(2) ! doesn't work
    call restab_boundary_with_fourier_wall !doesnt work well

CASE(3) ! refits all currents (active and passive) with F minimization cost function, no constraints. valid also for limited plasmas.
    call restab_F_function_full_fonfit 

CASE(313) ! refits all currents (active and passive) with F minimization cost function, no constraints. Fits also additional X points positions. Use only if n_xpoint_fit > 0
    call restab_F_function_full_fonfit_xpoints 

CASE(4) !finds active currents from scratch, eddy currents zero
    call restab_F_function_full_currents

CASE(41) !refit active currents from scratch for 1 point also isoflux
    call restab_1_timepoint_limits_xpoints_boundariz

CASE(5) !finds active currents from scratch including evolution from time t1 to time t2, with constraint on the consumed flux. eddy currents = 0.
    call restab_2_timepoints_evolution

CASE(6) !finds active currents from scratch, eddy currents zero. minimize magnetic energy and intercoil forces.
    call restab_F_function_full_currents_forces

CASE(7) !finds active currents from scratch, eddy currents zero. minimize magnetic energy and respect current limits.
    call restab_F_function_full_currents_limits

CASE(8) !finds active currents from scratch including evolution from time t1 to time t2, with constraint on the consumed flux. eddy currents = 0. also respect current limits
    call restab_2_timepoints_evolution_limits

CASE(818) !finds active currents from scratch including evolution from time j-1 to time j, with constraint on the consumed flux. eddy currents = 0. also respect current limits, with xpoints. This one does isoflux
    call restab_j_timepoints_evolution_limits_xpoints_boundariz ! uses full boundary

END SELECT

return
end subroutine solve_gse2d_fbe_full_feqis

!--------------------------------------------------------------------
subroutine solve_gse2d_fbe_full_feqis_1turn(j_init, j_stab, raxold, zaxold)

use fbe_core, only : nr2, nz2, iaxis, jaxis, &
    Rrect, Zrect, dr, dz, rax, zax, &
    jrz, psirz, psiextrz, psiplasrz, psistabr, psistabz, &
    solve_fbe_instantaneous
use feqis_circuit, only : interp_j_fromrhotorz
use global_params, only: iplasma
use pbe_core, only : raxp, zaxp

use feqis_tools, only: closest_index

implicit none

integer, intent(in) :: j_init, j_stab
double precision, intent(in) :: raxold, zaxold

integer :: i, j
double precision :: curr, dum1, dum2, zum1, zum2, delr, delz
double precision, dimension(9) :: c
double precision, dimension(nr2, nz2) :: g

!first, initialized initial guess coming from prescribed boundary current density: jrhoteta
if (j_init == 0) then
    call interp_j_fromrhotorz
! Rescale current density
    curr = SUM(jrz(1:nr2, 1:nz2)) *dr*dz
    jrz = jrz/curr*iplasma
    rax = raxp
    zax = zaxp
    iaxis = closest_index(rax, Rrect(1), dr)
    jaxis = closest_index(zax, Zrect(1), dz)
    rax = Rrect(iaxis)
    zax = Zrect(jaxis)
endif

call solve_fbe_instantaneous(j_init, j_stab, raxold, zaxold)

return
end subroutine solve_gse2d_fbe_full_feqis_1turn

!--------------------------------------------------------------------
subroutine FEQISUPDATE(coilzzz, nccc)

use pi_vars, only: GPI2
use fbe_core, only: nconduc, curconduc
use feqis_circuit, only: cur_con_old, &
    psi_cur_old, psiplasmatoconduc
use transport2fbe, only: fast_mode

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
use fbe_core, only: nconduc, curconduc
use feqis_circuit, only: tau_new, tau_old, cur_con_old, voltage, &
    dpc, psiplasmatoconduc, resconduc, indconduc, psi_cur_old
use feqis_tools, only: solve_circuit_equations
use transport2fbe, only: tau_circuit_feqis, tau_gseq_feqis, activate_coil_feqis, &
    n_equivalence, reconnect_circuits, new_equivalence, &
    use_reduce_circuit, current_limit_feqis, force_coil, resistance_change, &
    new_resistance

implicit none

integer, intent(in) :: j_init

integer :: i, ic, j, k, iii, jjj, invertcommand, i_equivalence, i_cnew, firstcall
integer, dimension(500) :: rem_coils
double precision, dimension(500) :: dpctemp, vtemp, curotemp, curtemp
double precision, dimension(500, 500) :: restemp, indtemp

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

if (resistance_change == 1) then
    invertcommand = 1
    do k=1, nconduc
        if (new_resistance(k)>0.) then
            resconduc(k,k) = new_resistance(k)
        endif
    enddo
endif

if (reconnect_circuits == 1.and.use_reduce_circuit == 1) then
    rem_coils=0
    invertcommand = 1
    i_cnew = ic
    indtemp = 0.
    restemp = 0.
    indtemp(1:ic, 1:ic) = indconduc(1:ic, 1:ic)
    restemp(1:ic, 1:ic) = resconduc(1:ic, 1:ic)
    do k=1, n_equivalence
        i_equivalence = 1
        do i=1, nconduc
            if (new_equivalence(i, k) > 0 .and. i_equivalence == 1) then
                i_equivalence = 0
                do j=i + 1, nconduc
                    if (new_equivalence(j, k) == new_equivalence(i, k)) then
                        indtemp(i, i) = indtemp(i, i) + indconduc(j, j) + indconduc(i, j) + indconduc(j, i)
                        restemp(i, i) = restemp(i, i) + resconduc(j, j) + resconduc(i, j) + resconduc(j, i)
                        do iii=1, nconduc
                            if (iii /= i .and. iii /= j) then
                                indtemp(i, iii) = indtemp(i, iii) + indconduc(j, iii)
                                restemp(i, iii) = restemp(i, iii) + resconduc(j, iii)
                                indtemp(iii, i) = indtemp(iii, i) + indconduc(iii, j)
                                restemp(iii, i) = restemp(iii, i) + resconduc(iii, j)
                            endif
                        enddo
                        rem_coils(j) = new_equivalence(j, k)
                        ic = ic - 1
                    endif
                enddo
            endif !found equivalence
            if (activate_coil_feqis(i) == 0) then
                rem_coils(i) = 999999
                ic = ic - 1
            endif
        enddo !conductors
    ! Remove obsolete columns and rows
        do jjj=1, nconduc
            j=nconduc-jjj+1
            if (rem_coils(j) > 0) then
                indtemp(1:nconduc, j:nconduc) = indtemp(1:nconduc, j + 1:nconduc + 1)
                restemp(1:nconduc, j:nconduc) = restemp(1:nconduc, j + 1:nconduc + 1)
                indtemp(j:nconduc, 1:nconduc) = indtemp(j + 1:nconduc + 1, 1:nconduc)
                restemp(j:nconduc, 1:nconduc) = restemp(j + 1:nconduc + 1, 1:nconduc)
            endif
        enddo
    enddo ! n equivalences
    i_cnew = ic
endif

if (use_reduce_circuit == 1) then
    dpctemp  = 0.
    vtemp    = 0.
    curotemp = 0.
    dpctemp (1:nconduc) = dpc        (1:nconduc)
    vtemp   (1:nconduc) = voltage    (1:nconduc)
    curotemp(1:nconduc) = cur_con_old(1:nconduc)
    do k=1, n_equivalence
        i_equivalence = 1
        do i=1, nconduc
            if (new_equivalence(i, k) > 0 .and. i_equivalence == 1) then
                i_equivalence = 0
                do j=i+1, nconduc
                    if (new_equivalence(j, k) == new_equivalence(i, k)) then
                        dpctemp(i) = dpctemp(i) + dpctemp(j)
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
            if (rem_coils(i) == 999999) then !non-activated
                curconduc(i) = 0.
            else
                curconduc(i) = curtemp(rem_coils(i))
            endif
            curtemp(i + 1:nconduc + 1) = curtemp(i:nconduc)
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
use pbe_core, only: nrho, nteta, raxp, zaxp, rbndp, zbndp, &
    teta, dteta, tetaexp, raxp, zaxp, &
    pressure, ipol, pprime, ffprime, &
    psigrid, psigrida, rexp, zexp
use fbe_core, only: nrho2d, use_limiter, &
    dr_factor_init, dz_factor_init, &
    rbnd, zbnd, &
    psia_2d, ffp_2d, ppp_2d, &
    psistabr, psistabz, psibnd
use feqis_circuit, only: ncoils
use global_params, only: iplasma, Rgeom0, Btor0
use transport2fbe, only: dr_factor_init_astra, dz_factor_init_astra, &
    tau_circuit_feqis, tau_gseq_feqis, activate_coil_feqis, current_limit_feqis, &
    raxis_astra, zaxis_astra, psi0_astra, psib_astra, use_limiter_astra, &
    fix_shape_after_fbe_off
use numerical_tools, only: linterp
use feqis_tools, only: pol_angle

implicit none

integer, intent(in) :: j_call, ifplasma
type(type_parameters ), intent(in) :: params
type(type_equilibrium), intent(in) :: equil_in

integer :: it_was_fbe_before, i, j, k
double precision, dimension(700) :: rdum, zdum, tdum

data it_was_fbe_before/0/
save it_was_fbe_before

if (j_call == 0) then
    nteta = equil_in%eqgeometry%boundary%npoints
    nrho = params%neql
    allocate(psigrid(nrho))
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
    max_iter = 1000 !hardwired
! teta for polar grid, goes from 0 to 2*pi-dteta, but point nt + 1 is the periodic one
    psi0_astra = equil_in%profiles_1d%psi(1)
    psib_astra = equil_in%profiles_1d%psi(nrho)
    raxp = raxis_astra
    zaxp = zaxis_astra
    psibnd = -1.e6
    use_limiter = use_limiter_astra
    allocate(teta(nteta+1))
    allocate(pressure(nrho))
    allocate(pprime(nrho))
    allocate(ffprime(nrho))
    allocate(psigrida(nrho))
    allocate(ipol(nrho))
    allocate(rexp(2*nteta))
    allocate(zexp(2*nteta))
    allocate(tetaexp(2*nteta))
    allocate(rbndp(2*nteta))
    allocate(zbndp(2*nteta))
    allocate(rbnd(2*nteta))
    allocate(zbnd(2*nteta))
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

    if (params%k_fixfree == 1.and.fix_shape_after_fbe_off == 1) it_was_fbe_before = 1

    if (params%k_fixfree == 0 .and. it_was_fbe_before == 0) then ! If 1, comes from free boundary
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
use feqis_circuit, only: nactive, npassive, ncoils, nblocks, &
    r_cond, z_cond, &
    rcoil, zcoil, drcoil, dzcoil, anglecoil, anglehcoil, mequivalence, &
    resconduc, indconduc, psiplasmatoconduc, &
    voltage, voltage_old, cur_con_old, &
    psi_cur_old, dpc, nferromag
use fbe_core, only: nr, nr1, nr2, nz, nz1, nz2, &
    nconduc, nlimiter, ngbnd, &
    lim_minr, lim_maxr, lim_minz, lim_maxz, &
    rmin, rmax, zmin, zmax, Rrect, Zrect, dr, dz, rcomp, zcomp, &
    limiterr, limiterz, alpsep, curconduc, &
    zlimpotential, green_bnd_f, &
    jrz, psirz, psiextrz, psiplasrz, u_n, area_eff, &
    psiferro
use feqis_circuit, only: nferromag, psiplasmatoconduc, &
    voltage, voltage_old, cur_con_old, psi_cur_old, dpc
use ferromagstructure, only: type_ferromag
use green_function, only: greeni
use green_matrix, only: dgreenirj, dgreenizj, dgreenirpl, dgreenizpl
use outcmn_inc, only: machine
use transport2fbe, only: cur_init, use_isoflux, n_isoflux, r_isoflux, z_isoflux, which_x_point, &
    voltage_limits_active_coils, sigma_isoflux

implicit none

type(type_ferromag), dimension(:), allocatable :: ferromag
integer :: i, j, ii, jj, nferrosub, imagvalues
integer, dimension(:), allocatable :: n_sames
character(len=80) :: fname, dummy

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
if (allocated(Rrect)) deallocate(Rrect)
if (allocated(Zrect)) deallocate(Zrect)
allocate(Rrect(nr2))
allocate(Zrect(nz2))
allocate(rcomp(nr))
allocate(zcomp(nz))
do i=1, nr2
    Rrect(i) = rmin + (i - 1.)*(rmax - rmin)/nr1     ! computational domain is r(2:nr + 1), boundaries are r(1) and r(nr + 2)
enddo
do i=1, nz2
    Zrect(i) = zmin + (i - 1.)*(zmax - zmin)/nz1
enddo
rcomp(1:nr) = Rrect(2:nr1)
zcomp(1:nz) = Zrect(2:nz1)
dr = Rrect(2) - Rrect(1)
dz = Zrect(2) - Zrect(1)

!some allocate
allocate(sintable(nz, nz))
allocate(costable(nz))

do i=1, nz
    costable(i) = cos(i*GPI/(nz + 1))
    do j=1, nz
        sintable(i, j) = sin(i*j*GPI/(nz + 1))
    enddo
enddo

! Load everything from file
read(32, *) nactive, npassive
read(32, *) ncoils
allocate(rcoil(ncoils))
allocate(zcoil(ncoils))
allocate(drcoil(ncoils))
allocate(dzcoil(ncoils))
allocate(anglehcoil(ncoils))
allocate(anglecoil(ncoils))
allocate(mequivalence(ncoils))
do i=1, ncoils
    read(32, *) rcoil(i), zcoil(i), drcoil(i), dzcoil(i), anglehcoil(i), anglecoil(i), mequivalence(i)
enddo
read(32, *) nlimiter
allocate(limiterr(nlimiter))
allocate(limiterz(nlimiter))
do i=1, nlimiter
    read(32, *) limiterr(i), limiterz(i)
enddo
read(32, *) lim_maxR
read(32, *) lim_minR
read(32, *) lim_maxZ
read(32, *) lim_minZ
allocate(r_cond(nactive+npassive))
allocate(z_cond(nactive+npassive))
allocate(n_sames(nactive))
r_cond = 0.
z_cond = 0.
n_sames = 0
do i=1, ncoils
    r_cond(mequivalence(i)) = r_cond(mequivalence(i)) + rcoil(i)
    z_cond(mequivalence(i)) = z_cond(mequivalence(i)) + zcoil(i)
    n_sames(mequivalence(i)) = n_sames(mequivalence(i)) + 1
enddo
do i=1, nactive
    r_cond(i) = r_cond(i)/(0.+n_sames(i))
    z_cond(i) = z_cond(i)/(0.+n_sames(i))
enddo
deallocate(n_sames)
do i=1, npassive
    read(32, *) r_cond(nactive+i), z_cond(nactive+i)
enddo
read(32, *) nconduc
allocate(curconduc(nconduc))
allocate(voltage(nconduc))
allocate(voltage_old(nconduc))
allocate(cur_con_old(nconduc))
allocate(indconduc(nconduc, nconduc))
allocate(resconduc(nconduc, nconduc))
allocate(psiplasmatoconduc(nconduc))
do i=1, nconduc
    read(32, *) indconduc(i, 1:nconduc)
enddo
read(32, *) nactive, nconduc
resconduc = 0.d0
do i=1, nactive
    read(32, *) (resconduc(i, j), j=1, nactive)
enddo
do i=nactive+1, nconduc
    read(32, *) resconduc(i, i)
enddo
allocate(greeni(nr2, nz2, nconduc))
do i=1, nconduc
    do j=1, nr2
        read(32, *) greeni(j, 1:nz2, i)
    enddo
enddo
read(32, *) nblocks
allocate(dgreenirj(nblocks, nblocks))
allocate(dgreenizj(nblocks, nblocks))
allocate(dgreenirpl(nr2, nz2, nblocks))
allocate(dgreenizpl(nr2, nz2, nblocks))
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
allocate(zlimpotential(nr2, nz2))
do jj=1, nz2
    do ii=1, nr2
        read(32, *) zlimpotential(ii, jj)
    enddo
enddo
read(32, *) ngbnd
allocate(green_bnd_f(ngbnd))
read(32, *) green_bnd_f(1:ngbnd)

read(32, *) nferromag
if (nferromag >= 1) then
    allocate(ferromag(nferromag))
    do i=1, nferromag
        read(32,*) nferrosub, imagvalues, ferromag(i)%position%sigma_surface 
        ferromag(i)%position%npoints   = nferrosub
        ferromag(i)%mhrelation%nvalues = imagvalues
        allocate(ferromag(i)%position%r(nferrosub))
        allocate(ferromag(i)%position%z(nferrosub))
        allocate(ferromag(i)%position%tanangl(nferrosub))
        allocate(ferromag(i)%position%length(nferrosub))
        allocate(ferromag(i)%position%magnetizationchi(nferrosub))
        allocate(ferromag(i)%position%Btangfield(nferrosub))
        allocate(ferromag(i)%position%current(nferrosub))
        allocate(ferromag(i)%mhrelation%chi(imagvalues))
        allocate(ferromag(i)%mhrelation%h(imagvalues))
        allocate(ferromag(i)%mutual_matrix%Mij(nferrosub, nferrosub))
        do j=1, imagvalues
            read(32, *) ferromag(i)%mhrelation%chi(j), ferromag(i)%mhrelation%h(j)
        enddo
        do j=1, nferrosub
           read(32, *) ferromag(i)%position%r(j), ferromag(i)%position%z(j), &
               ferromag(i)%position%tanangl(j), ferromag(i)%position%length(j)
        enddo
        do jj=1, nferrosub
            do ii=1, nferrosub
                read(32, *) ferromag(i)%mutual_matrix%Mij(ii, jj)
            enddo
        enddo
    enddo
endif

close(32)

! assign initial currents from astra
curconduc(1:nconduc) = cur_init(1:nconduc)

allocate(jrz(nr2, nz2))
if (allocated(psirz)) deallocate(psirz)
allocate(psirz(nr2, nz2))
psirz = 0.
allocate(psiextrz(nr2, nz2))
psiextrz = 0.
allocate(psiplasrz(nr2, nz2))
psiplasrz = 0.
allocate(psiferro(nr2, nz2))
allocate(u_n(nr2, nz2))
allocate(area_eff(nr2, nz2))
allocate(psi_cur_old(nconduc))
allocate(dpc(nconduc))
allocate(voltage_limits_active_coils(nactive, 2))

if (n_isoflux > 0 .and. not(allocated(r_isoflux))) then
    allocate(r_isoflux(n_isoflux))
    allocate(z_isoflux(n_isoflux))
    allocate(which_x_point(n_isoflux))
    allocate(sigma_isoflux(n_isoflux))
endif

return
end subroutine equil_feqis_init_circ

!--------------------------------------------------------------------
subroutine fix_boundary_feqis(j_init)

use pbe_core, only: nrho, nteta, raxp, zaxp, rbndp, zbndp, rho, teta, &
    psiaxisp, psirhoteta, psigrida, psibndp, &
    psia_1d, ffp_1d, ppp_1d, &
    ffprime, pprime, pressure, ipol, &
    Rpol, Zpol, Rpul, Zpul, jrhoteta
use global_params, only: Rgeom0, Btor0, iplasma
use scalars, only: li3, li_aug, betapol, betapol_iter, wkin, bpkin, psplex
use transport2fbe, only: raxis_astra, zaxis_astra, psi0_astra, psib_astra, &
    solve_fix
use metric_coefficients_pbe, only: lambda2d, R_curr_0d, Z_curr_0d, dator, fsa_kernel
use transfer_functions, only: rpbez, zpbez, psibez, t2dbez, &
    g1bez, g2bez, gm1bez, gm4bez, gm41bez, gm5bez, ggrhobez, &
    areatbez, surfbez, perimbez, vbez, qbez, phibez, &
    bmaxbez, bminbez, bdb0bez, fofbbez, bcell2dbez, bpcell2dbez, &
    ffprimebez, pprimebez, pressbez, ipolbez, rinbez, routbez, &
    kbez, triaubez, trialbez, shifbez, rbp2_b2bez, rmin2dbez, dpsidvbez, &
    jrhobez, shivbez, squarebez, g2ibez, &
    rminbez, bpcellbez, bcellbez
use pi_vars, only: GPI, GPI2, GPI4, muvac

implicit none

integer, intent(in) :: j_init
integer :: i, j, jr, jt, ierr
double precision :: psiaxis_new, cnorm, rax_new, zax_new, rhoedge
double precision, dimension(nrho) :: q_new, effprimp, epprimp
double precision, dimension(nteta) :: thetap_i
double precision, dimension(512, 512) :: psisave
double precision, dimension(nrho, nteta) :: rmaj2, jcbn2, darea2, &
    r_min, yy2, jrho2, gradr2, darea, dl_dt

data ierr/0/
save psisave, ierr

! initial guess
if (j_init == 0) then
    allocate(rho(nrho, nteta+1))
    allocate(rpol(nrho, nteta))
    allocate(zpol(nrho, nteta))
    allocate(rpul(nrho, nteta))
    allocate(zpul(nrho, nteta))
    allocate(psirhoteta(nrho, nteta))
    allocate(jrhoteta(nrho, nteta+1))
    allocate(psia_1d(nrho))
    allocate(ppp_1d(nrho))
    allocate(ffp_1d(nrho))
    allocate(ffprimebez(nrho))
    allocate(dpsidvbez(nrho))
    allocate(psibez(nrho))
    allocate(g2bez(nrho))
    allocate(g2ibez(nrho))
    allocate(gm1bez(nrho))
    allocate(routbez(nrho))
    allocate(rinbez(nrho))
    allocate(vbez(nrho))
    allocate(g1bez(nrho))
    allocate(gm41bez(nrho))
    allocate(ggrhobez(nrho))
    allocate(bmaxbez(nrho))
    allocate(bminbez(nrho))
    allocate(gm4bez(nrho))
    allocate(bdb0bez(nrho))
    allocate(gm5bez(nrho))
    allocate(fofbbez(nrho))
    allocate(areatbez(nrho))
    allocate(perimbez(nrho))
    allocate(shifbez(nrho))
    allocate(kbez(nrho))
    allocate(surfbez(nrho))
    allocate(triaubez(nrho))
    allocate(trialbez(nrho))
    allocate(phibez(nrho))
    allocate(qbez(nrho))
    allocate(t2dbez(nteta))
    allocate(rbp2_b2bez(nrho))
    allocate(pprimebez(nrho))
    allocate(pressbez(nrho))
    allocate(ipolbez(nrho))
    allocate(shivbez(nrho))
    allocate(squarebez(nrho))
    allocate(rpbez(nrho, nteta))
    allocate(zpbez(nrho, nteta))
    allocate(rminbez(nrho, nteta))
    allocate(bpcellbez(nrho, nteta))
    allocate(bcellbez(nrho, nteta))
    allocate(rmin2dbez(nrho, nteta))
    allocate(jrhobez(nrho, nteta))
    allocate(bpcell2dbez(nrho, nteta))
    allocate(bcell2dbez(nrho, nteta))
    allocate(lambda2d(nrho, nteta))
    allocate(fsa_kernel(nrho, nteta))
    allocate(dator(nrho, nteta))

    raxp = raxis_astra
    zaxp = zaxis_astra
    psiaxisp = psi0_astra
    psibndp = psib_astra
    do jr=1, Nrho
        do jt=1, Nteta
            lambda2d(jr, jt) = (jr - 1)/(Nrho - 1.)
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
    jcbn2, q_new, rhoedge, darea2, epprimp, effprimp, r_min, yy2, gradr2, darea, ierr, dl_dt)

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
R_curr_0D = sqrt(sum(jrho2(1:jr, :)*Rmaj2(1:jr, :)**2 * darea2(1:jr, :))/ &
                 sum(jrho2(1:jr, :)*darea2(1:jr, :)))

!regrid
call build_2dgrid(nrho, nteta, psibez(1:nrho), &
    psirhoteta(1:nrho, 1:nteta),jrho2(1:nrho,1:nteta),darea2(1:nrho,1:nteta),YY2(1:nrho,1:nteta), &
    rgeom0, pressure(1:nrho), &
    btor0, ipol(1:nrho), iplasma, rpbez(1:nrho, 1:nteta), zpbez(1:nrho, 1:nteta), &
    rmaj2(1:nrho, 1:nteta), r_min(1:nrho, 1:nteta), jcbn2(1:nrho, 1:nteta), &
    thetap_i(1:nteta), q_new(1:nrho), rhoedge, gradr2(1:nrho, 1:nteta), &
    g2bez(1:nrho), gm1bez(1:nrho), areatbez(1:nrho), perimbez(1:nrho), vbez(1:nrho), &
    g1bez(1:nrho), ggrhobez(1:nrho), bmaxbez(1:nrho), bminbez(1:nrho), &
    gm4bez(1:nrho), bdb0bez(1:nrho), gm5bez(1:nrho), fofbbez(1:nrho), surfbez(1:nrho), & ! lateral surface
    li3, betapol, psplex, &
    bpcell2dbez(1:nrho, 1:nteta), bcell2dbez(1:nrho, 1:nteta), &
    routbez(1:nrho), rinbez(1:nrho), kbez(1:nrho), triaubez(1:nrho), trialbez(1:nrho), shifbez(1:nrho), &
    gm41bez(1:nrho), qbez(1:nrho), shivbez(1:nrho), squarebez(1:nrho), li_aug, betapol_iter, dl_dt, &
    wkin, bpkin)

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
rho(1:nrho, nteta+1) = rho(1:nrho, 1)
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
use fbe_core, only: nr2, nz2, psirz, &
    Rrect, Zrect, psiaxis, psibnd
use pbe_core, only: nrho, nteta, &
    psiaxisp, psibndp, psirhoteta
use global_params, only: iplasma
use scalars, only: betapol, li3, li_aug, betapol_iter, wkin, bpkin, psplex
use transfer_functions, only: rpbez, zpbez, t2dbez, &
    rinbez, routbez, rmin2dbez, vbez, areatbez, perimbez, surfbez, &
    kbez, shifbez, triaubez, trialbez, qbez, phibez, &
    g1bez, g2bez, g2ibez, gm1bez, gm4bez, gm41bez, gm5bez, ggrhobez, &
    bcell2dbez, bpcell2dbez, bminbez, bmaxbez, bdb0bez, fofbbez, &
    psibez, dpsidvbez, rbp2_b2bez, &
    ffprimebez, pprimebez, pressbez, ipolbez, jrhobez, shivbez, squarebez
use metric_coefficients_pbe, only: dator

implicit none

integer :: i
type(type_equilibrium), intent(inout) :: equil_out

equil_out%global_param%psplex   = psplex
equil_out%global_param%psibound = -GPI2*psibndp
equil_out%global_param%psiaxis  = -GPI2*psiaxisp
equil_out%global_param%li3      = li3
equil_out%global_param%li_aug   = li_aug
equil_out%global_param%betpol   = betapol
equil_out%global_param%betpol_iter   = betapol_iter
equil_out%global_param%i_plasma = iplasma*1.e6
equil_out%global_param%wkin   = wkin ! volume integral of pressure
equil_out%global_param%bpkin   = bpkin ! volume integral of Bp**2/(2mu0)

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
if (allocated(Rrect)) then
    equil_out%eqgeometry%rectgrid%r2d(1:nr2) = Rrect(1:nr2)
    equil_out%eqgeometry%rectgrid%z2d(1:nz2) = Zrect(1:nz2)
    equil_out%eqgeometry%rectgrid%psirz2d(1:nr2, 1:nz2) = psirz(1:nr2, 1:nz2)
    equil_out%eqgeometry%rectgrid%psi_axis     = psiaxis
    equil_out%eqgeometry%rectgrid%psi_boundary = psibnd
endif

equil_out%profiles_1d%ffprime   (1:nrho) = ffprimebez(1:nrho)
equil_out%profiles_1d%pprime    (1:nrho) = pprimebez(1:nrho)
equil_out%profiles_1d%pressure  (1:nrho) = pressbez (1:nrho)
equil_out%profiles_1d%rho_tor   (1:nrho) = (/ ((i-1.)/(nrho-1.), i=1, nrho) /)
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
equil_out%profiles_1d%tria_lower(1:nrho) = trialbez(1:nrho)
equil_out%profiles_1d%phi       (1:nrho) = phibez(1:nrho)
equil_out%profiles_1d%q         (1:nrho) = qbez(1:nrho)
equil_out%profiles_1d%rbp_b2    (1:nrho) = rbp2_b2bez(1:nrho)

return
end subroutine equil_assignments

!--------------------------------------------------------------------
subroutine convert_boundary_to_pbe

use pi_vars, only: GPI2
use fbe_core, only: nr, nr2, nz, nbnd, iaxis, jaxis, &
    raus, rinner, zbot, ztop, &
    Rrect, Zrect, dr, dz, rax, zax, rbnd, zbnd, &
    psiaxis, psibnd, psirz
use pbe_core, only: nteta, teta, dteta, &
    raxp, zaxp, rbndp, zbndp, &
    psiaxisp, psibndp
use feqis_tools, only: pol_angle, interp2d_psi

implicit none

integer :: i, j, k, j4
double precision :: x1, x2, t1, t2, t3, z1, z2, z3, x11, dx
double precision, dimension(500) :: teta_fbe

do i=1, nteta + 1
    teta(i) = GPI2*(i - 1.)/(nteta + 0.)
enddo
dteta = teta(2) - teta(1)

! Find boundary
j = jaxis
do i=iaxis, nr2
    if (psirz(i, j) >= psibnd) k = i
    if (psirz(i, j) <= psibnd) EXIT
enddo

if (psirz(k, j) == psibnd) then
    rbnd(1) = Rrect(k)
else
    rbnd(1) = Rrect(k) - (psirz(k, j) - psibnd)/(psirz(k, j) - psirz(k-1, j))*dr
endif
zbnd(1) = Zrect(j)
teta_fbe(1) = pol_angle(rax, zax, rbnd(1), zbnd(1))

theta_loop: do i=2, nteta
    dx = sqrt((dr*cos(teta(i)))**2 + (dz*sin(teta(i)))**2)
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
    t3 = interp2d_psi(t1, t2, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))

    if (t3 == psibnd) then
        rbnd(i) = t1
        zbnd(i) = t2
    elseif (t3 < psibnd) then
        do
            z1 = rax + (x1 - dx)*cos(teta_fbe(i))
            z2 = zax + (x1 - dx)*sin(teta_fbe(i))
            z3 = interp2d_psi(z1, z2, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
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
            z3 = interp2d_psi(z1, z2, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
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

teta (1:nteta) = teta_fbe(1:nteta)
rbndp(1:nteta) = rbnd(1:nteta)
zbndp(1:nteta) = zbnd(1:nteta)
teta (nteta + 1) = teta(nteta) + dteta
rbndp(nteta + 1) = rbnd(1) 
zbndp(nteta + 1) = zbnd(1)

raxp = rax
zaxp = zax
psibndp  = psibnd
psiaxisp = psiaxis

return
end subroutine convert_boundary_to_pbe

!--------------------------------------------------------------------
subroutine estimate_boundary_to_pbe(rbnd, zbnd, ntetaz)

use pi_vars, only: GPI2
use fbe_core, only: nr, nr2, nz, iaxis, jaxis, &
    raus, rinner, zbot, ztop, &
    Rrect, Zrect, dr, dz, rax, zax,  &
    psiaxis, psibnd, psirz
use pbe_core, only: teta, dteta
use feqis_tools, only: pol_angle, interp2d_psi

implicit none

integer, intent(IN) :: ntetaz
double precision, intent(OUT), dimension(ntetaz) :: rbnd, zbnd

integer :: i, j, k, j4, nteta
double precision :: x1, x2, t1, t2, t3, z1, z2, z3, x11, dx
double precision, dimension(500) :: teta_fbe

nteta = ntetaz
rbnd = 0.
zbnd = 0.

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
    rbnd(1) = Rrect(k)
else
    rbnd(1) = Rrect(k) - (psirz(k, j) - psibnd)/(psirz(k, j) - psirz(k-1, j))*dr
endif
zbnd(1) = Zrect(j)
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
    t3 = interp2d_psi(t1, t2, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))

    if (t3 == psibnd) then
        rbnd(i) = t1
        zbnd(i) = t2
    elseif (t3 < psibnd) then
        do
            z1 = rax + (x1 - dx)*cos(teta_fbe(i))
            z2 = zax + (x1 - dx)*sin(teta_fbe(i))
            z3 = interp2d_psi(z1, z2, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
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
            z3 = interp2d_psi(z1, z2, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
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

return
end subroutine estimate_boundary_to_pbe


!---------------------------------------------------------------------
subroutine coil_forces_feqis(ncoilz, force_R, force_Z, plasma_state)
! this one is only between coils and coils
use feqis_circuit, only: nblocks, npassive, mequivalence
use fbe_core, only: jrz, nr2, nz2, area_eff, &
    curconduc
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
        x1 =  sum(jrz(1:nr2, 1:nz2) * area_eff(1:nr2, 1:nz2) * dgreeniZpl(1:nr2, 1:nz2, i))
        force_Z(i) = force_Z(i) + curconduc(mequivalence(i)) * x1
    enddo
endif

!block-to-block
do i=1, nblock_a
    do j=1, nblock_a
        if (i /= j) then
            force_R(i) = force_R(i) + curconduc(mequivalence(j)) * curconduc(mequivalence(i)) * dgreeniRj(i, j)
            force_Z(i) = force_Z(i) + curconduc(mequivalence(j)) * curconduc(mequivalence(i)) * dgreeniZj(i, j)
        endif
    enddo
enddo

! force_R and force_Z are F_R and F_Z components in [N] for each block , does not include forces from the passive elements or on the passive elements.

return
end subroutine coil_forces_feqis

!---------------------------------------------------------------------
subroutine all_forces_feqis(ncoilz, force_R, force_Z, plasma_state, plasma_force)

!this one is between everything , including plasma. ncoilz = nblocks (active subcoils and passive leemnets)

use feqis_circuit, only: nblocks, npassive, mequivalence
use fbe_core, only: jrz, nr2, nz2, area_eff, &
    curconduc
use green_matrix, only: dgreenirpl, dgreenizpl, dgreenirj, dgreenizj

integer, intent(in) :: ncoilz, plasma_state
double precision, intent(out), dimension(ncoilz) :: force_R, force_Z
double precision, intent(out) :: plasma_force(2) !index 1 is Radial, index 2 is vertical

integer :: i, j, k, nblock_a
double precision :: x1

force_R = 0.
force_Z = 0.
nblock_a = nblocks - npassive

if (plasma_state == 1) then !not sure about the plasma response...
    do i=1, nblocks
        x1 =  sum(jrz(1:nr2, 1:nz2) * area_eff(1:nr2, 1:nz2) * dgreeniRpl(1:nr2, 1:nz2, i))
        force_R(i) = force_R(i) + curconduc(mequivalence(i)) * x1
        x1 =  sum(jrz(1:nr2, 1:nz2) * area_eff(1:nr2, 1:nz2) * dgreeniZpl(1:nr2, 1:nz2, i))
        force_Z(i) = force_Z(i) + curconduc(mequivalence(i)) * x1
    enddo
endif

!full force on plasma
plasma_force(1) = -sum(force_R)
plasma_force(2) = -sum(force_Z)

!block-to-block
do i=1, nblocks
    do j=1, nblocks
        if (i /= j) then
            force_R(i) = force_R(i) + curconduc(mequivalence(j)) * curconduc(mequivalence(i)) * dgreeniRj(i, j)
            force_Z(i) = force_Z(i) + curconduc(mequivalence(j)) * curconduc(mequivalence(i)) * dgreeniZj(i, j)
        endif
    enddo
enddo

return
end subroutine all_forces_feqis

!---------------------------------------------------------------------
subroutine all_forces_feqis_components(ncoilz, force_R, force_Z, force_tot, plasma_state, plasma_contrib_R, plasma_contrib_Z)

!this one is between everything , including plasma. ncoilz = nblocks-npassive (active subcoils only)

use feqis_circuit, only: nblocks, npassive, mequivalence
use fbe_core, only: jrz, nr2, nz2, area_eff, curconduc
use green_matrix, only: dgreenirpl, dgreenizpl, dgreenirj, dgreenizj

integer, intent(in) :: ncoilz, plasma_state
double precision, intent(out), dimension(ncoilz) :: force_R, force_Z, force_tot, plasma_contrib_R, plasma_contrib_Z

integer :: i, j, k, nblock_a
double precision :: x1

force_R = 0.
force_Z = 0.
force_tot = 0.
plasma_contrib = 0.
nblock_a = nblocks - npassive

if (ncoilz /= nblock_a) then
    write(*, *) 'wrong nr of coils. please set the first argument to the call to: ', nblock_a
    stop
endif

if (plasma_state == 1) then !not sure about the plasma response...
    do i=1, nblock_a
        x1 =  sum(jrz(1:nr2, 1:nz2) * area_eff(1:nr2, 1:nz2) * dgreeniRpl(1:nr2, 1:nz2, i))
        force_R(i) = force_R(i) + curconduc(mequivalence(i)) * x1
        x1 =  sum(jrz(1:nr2, 1:nz2) * area_eff(1:nr2, 1:nz2) * dgreeniZpl(1:nr2, 1:nz2, i))
        force_Z(i) = force_Z(i) + curconduc(mequivalence(i)) * x1
        plasma_contrib_R(i) = force_R(i)
        plasma_contrib_Z(i) = force_Z(i)
    enddo
endif

!block-to-block
do i=1, nblock_a
    do j=1, nblock_a
        if (i /= j) then
            force_R(i) = force_R(i) + curconduc(mequivalence(j)) * curconduc(mequivalence(i)) * dgreeniRj(i, j)
            force_Z(i) = force_Z(i) + curconduc(mequivalence(j)) * curconduc(mequivalence(i)) * dgreeniZj(i, j)
        endif
    enddo
enddo

force_tot = sqrt(force_R**2 + force_Z**2)

return
end subroutine all_forces_feqis_components
