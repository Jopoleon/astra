subroutine full_system_advance_ef(j_init)

use ef_circuit, only: psi_cur_old, psiplasmatoconduc, max_iter, &
    nconduc, curconduc, cur_con_old, rax, zax, iplasma
use exchange_with_astra, only: execute_plasma, fast_mode
use debugger_ef, only: markloc_ef, debug

implicit none

integer, intent(inout) :: j_init

integer :: j_iter
double precision :: error_temp
double precision :: cur_temp(300)

call markloc_ef('full_system_advance_ef', debug_lev=debug)

! time stepping
! at iteration 0, dpsidt = 0
if (j_init == 0) then
    write(*, *) 'init full system'
!first do full equilibrium solution at time t=0
    call psi_external_calc_ef
    call solve_gse2d_fbe_full_ef(0)
    call plasma_psi_to_coils_ef
    call circuit_eq_advance_ef(0)
    j_init = -1
    return
endif

! full iterations
write(*, *) 'iter full system'

if (execute_plasma == 1) then
    psi_cur_old = psiplasmatoconduc
endif
    
do j_iter=1, 2*max_iter
    call psi_external_calc_ef
    if (execute_plasma == 1) then
        call solve_gse2d_fbe_full_ef_1turn(1, 0, 0.d0, 0.d0)    
        call plasma_psi_to_coils_ef
    endif

    cur_temp(1: nconduc) = curconduc(1: nconduc)
    call circuit_eq_advance_ef(1)    
    write(*, *) 'iterations current ', j_iter, psiplasmatoconduc(5), &
        psi_cur_old(5), cur_con_old(5), curconduc(5)

    if (j_iter >= 90) pause
    write(8898, *) rax, zax

    if (fast_mode == 1) then
        j_init = -1
        return
    endif

    error_temp = sum(abs(cur_temp(1: nconduc) - curconduc(1: nconduc))) / &
                 (nconduc + 1.e-6)/iplasma
    write(*, *) 'error', error_temp
    if (error_temp <= 1.e-8) then
        if (execute_plasma == 1) then
            j_init = -1
        endif
        return    
    endif

    if (j_iter.gt.max_iter) then
        write(*, *) 'max number of iterations override... stopping', error_temp
        if (execute_plasma == 1) then
            j_init = -1
        endif
        return    
    endif
enddo

return
end subroutine full_system_advance_ef

!--------------------------------------------------------
subroutine solve_gse2d_fbe_full_ef(j_init)

use ef_circuit, only: nr2, nz2, psiextrz, r, z, rax, zax, raxp, zaxp, &
    iaxis, jaxis, psistabr, psistabz, dr, dz, redo_bnd, trax, tzax, &
    psirz, jrz
use exchange_with_astra, only: refit_mode, n_of_newton_iterations
use debugger_ef, only: markloc_ef, debug

implicit none

integer, intent(in) :: j_init
integer :: j_iter, j_iter2, j_cyclo, jeppa
double precision :: temp_err, temp_err2, raxold, zaxold, raxoldo, zaxoldo, &
    raxtmp, zaxtmp, det, psistab1o, psistab2o, psro, pszo, dist1, dist2, &
    cibapr, cibazr, rleft, rright, zup, zdown, dcrdr, dcrdz, dczdr, dczdz
double precision, dimension(300, 300) :: g000

call markloc_ef('solve_gse2d_fbe_full_ef', debug_lev=debug)

SELECT CASE(refit_mode)
CASE(-1) ! 1 turn only
    call solve_gse2d_fbe_full_ef_1turn(j_init, 0, 0.d0, 0.d0)    

CASE(1)
    call restab_axis_with_furier_wall

CASE(0)
!start iterations to find self-consistent solution
    g000(1: nr2, 1: nz2) = psiextrz(1: nr2, 1: nz2)
    call find_actual_index_ef(raxp, zaxp, iaxis, jaxis)
    rax = r(iaxis)
    zax = z(jaxis)
    raxold  = rax
    zaxold  = zax
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
    dist1 = dr
    dist2 = dz

    do j_iter2 = 1, 100000

        if (j_cyclo == 1) then
            raxold = raxoldo + dist1
            zaxold = zaxoldo
            psistab1o = psistabr
            psistab2o = psistabz
        endif
        if (j_cyclo == 2) then
            raxold = raxoldo
            zaxold = zaxoldo + dist2
            dcrdr = (psistabr - psistab1o)/dist1
            dczdr = (psistabz - psistab2o)/dist1
        endif
        if (j_cyclo == 3) then      ! inverse of (dcrdr dcrdz  ; dczdr dczdz ) = ( dczdz -dcrdz ; -dczdr dcrdr)
            dcrdz = (psistabr - psistab1o)/dist2
            dczdz = (psistabz - psistab2o)/dist2
            det = (dcrdr*dczdz - dcrdz*dczdr)
            cibapr = (dczdz*psistab1o - dcrdz*psistab2o)/det
            cibazr = (-dczdr*psistab1o + dcrdr*psistab2o)/det
            raxold  = raxoldo - cibapr
            zaxold  = zaxoldo - cibazr
            raxoldo = raxold
            zaxoldo = zaxold
        endif

!inner cycle, calculate psi1 and psi2

        psro = 1000.
        pszo = 1000.
        redo_bnd = 1

        do j_iter=1, 100000
            raxtmp = trax
            zaxtmp = tzax

            call solve_gse2d_fbe_full_ef_1turn(j_iter - 1 + j_iter2 - 1, 1, &
                raxold, zaxold)    

            temp_err = (abs(psro - psistabr) + abs(pszo - psistabz))
            write(*, *) 'psistab', j_iter, raxold, zaxold, rax, zax, &
                temp_err, psistabR, psistabZ
            write(444, '(8E25.11)') trax, tzax, rax, zax, psistabR, psistabZ

            if (j_iter >= 30) then
                temp_err = 1.e-10
                psistabr = 0.5*(psro + psistabr)
                psistabz = 0.5*(pszo + psistabz)
            endif

            if (temp_err <= 1.e-9) goto 941

            psro = psistabr
            pszo = psistabz
        enddo
        write(*, *) 'total iterations passed!'
        return

941 continue

        if (j_cyclo == 3) then    
            temp_err2 = abs(psistabr) + abs(psistabZ)
            write(3217, *)  rax, zax
            j_cyclo = 1
            write(*, *) n_of_newton_iterations, nint((0. + j_iter2)/3.)
            if (nint((0.+j_iter2)/3.) >= n_of_newton_iterations) EXIT
        else
            j_cyclo = j_cyclo + 1
        endif

        if (j_iter2 >= 400) EXIT
        if (temp_err2 <= 1.e-8) EXIT

    enddo

    open(32, file='fort.4444')
        write(32, *) r(1: nr2), z(1: nz2), jrz(1: nr2, 1: nz2)
    close(32)
    open(32, file='fort.4445')
        write(32, *) psiextrz(1: nr2, 1: nz2)
    close(32)
    open(32, file='fort.4447')
        write(32, *) psirz(1: nr2, 1: nz2)
    close(32)

    write(*, *) 'code converged', j_iter2
    write(*, *) 'total iterations passed!'

END SELECT

return
end subroutine solve_gse2d_fbe_full_ef

!----------------------------------------------------
subroutine restab_axis_with_furier_wall

use ef_circuit, only: nr2, nz2, iaxis, jaxis, jrz, nactive, &
    iplasma, psiplasrz, psistabr, psistabz, psirz, psiextrz, &
    curconduc, nconduc, npassive, &
    dr, dz, rax, zax, raxp, zaxp, r_cond, z_cond, R, Z
use green_matrix, only: greeni
use debugger_ef, only: markloc_ef, debug

implicit none

integer :: i, j, j_iter, iax, jax
double precision :: temp_err, dum1, dum2, zum1, delr, delz, &
    psistab1o, psistab2o, tin, tup, S_00r, C_00r, S_00z, C_00z, angle
double precision, dimension(6) :: ccc
double precision, dimension(8) :: ddipsi
double precision, dimension(9) :: bub, xub, yub
double precision, dimension(npassive) :: cosa, sina, g0_r, g0_z
double precision, dimension(258, 258) :: C_00, S_00
double precision, dimension(300, 300) :: g

call markloc_ef('restab_axis_with_furier_wall', debug_lev=debug)

!first, initialized initial guess coming from prescribed boundary current density:  jrhoteta
write(*, *) 'reinterp curr, restab'
call interp_j_fromrhotorz

!rescale plasma
dum1 = 0.
do j=1, nz2
    do i=1, nr2
        dum1 = dum1 + jrz(i, j)*dr*dz
    enddo
enddo
jrz = jrz/dum1*iplasma

rax = raxp
zax = zaxp
write(*, *) raxp, zaxp, dum1

call find_actual_index_ef(rax, zax, iaxis, jaxis)
iax = iaxis
jax = jaxis

! evaluate coils things
do i=1, npassive
    angle = ATAN2(z_cond(nactive+i) - zax, r_cond(nactive+i) -rax)
    cosa(i) = COS(angle)
    sina(i) = SIN(angle)
!find true axis
    xub(1) = r(iax-1)
    xub(2) = r(iax)
    xub(3) = r(iax+1)
    xub(4) = r(iax)
    xub(5) = r(iax)
    xub(6) = r(iax-1)
    xub(7) = r(iax-1)
    xub(8) = r(iax+1)
    xub(9) = r(iax+1)
    yub(1) = z(jax)
    yub(2) = z(jax)
    yub(3) = z(jax)
    yub(4) = z(jax-1)
    yub(5) = z(jax+1)
    yub(6) = z(jax-1)
    yub(7) = z(jax+1)
    yub(8) = z(jax-1)
    yub(9) = z(jax+1)
    bub(1) = greeni(iax-1, jax  , nactive+i)
    bub(2) = greeni(iax  , jax  , nactive+i)
    bub(3) = greeni(iax+1, jax  , nactive+i)
    bub(4) = greeni(iax  , jax-1, nactive+i)
    bub(5) = greeni(iax  , jax+1, nactive+i)
    bub(6) = greeni(iax-1, jax-1, nactive+i)
    bub(7) = greeni(iax-1, jax+1, nactive+i)
    bub(8) = greeni(iax+1, jax-1, nactive+i)
    bub(9) = greeni(iax+1, jax+1, nactive+i)
    call least_square_biquad_ef(xub, yub, bub, 9, ccc, dum1, dum2, zum1, ddipsi)
    g0_r(i) = ddipsi(1)
    g0_z(i) = ddipsi(2)
enddo

do j = 1, nz2
    do i=1, nr2
        C_00(i, j) = sum(greeni(i, j, nactive+1: nactive+npassive)*cosa)
        S_00(i, j) = sum(greeni(i, j, nactive+1: nactive+npassive)*sina)
    enddo
enddo
C_00r = sum(g0_r*cosa)
S_00r = sum(g0_r*sina)
C_00z = sum(g0_z*cosa)
S_00z = sum(g0_z*sina)

write(*, *) 'coils structure', c_00(iax, jax), s_00(iax, jax), c_00r, s_00r, c_00z, s_00z, iplasma

psistab1o = 1000.
psistab2o = 1000.

do j_iter=1, 30

    CALL CPU_TIME(tin)    
    g = 0.
    call solve_gs2d(g) !jrz as right hand side
    CALL CPU_TIME(tup)
    write(*, *) 'solve g0', tup - tin    
        open(32, file='fort.44491')
        write(32, *) g(1: nr2, 1: nz2)
    close(32)

    write(*, *) 'stop here2'
!find boundary condition using g
    CALL CPU_TIME(tin)    
    call boundary_ef(g)  ! gbound = integral (Green*dg/dn) over the boundary
    write(*, *) 'stop here3'
    CALL CPU_TIME(tup)
    write(*, *) tup - tin    

    CALL CPU_TIME(tin)    
    call solve_gs2d(g) ! again jrz as right hand side
    psiplasrz(1: nr2, 1: nz2) = g(1: nr2, 1: nz2)
    write(*, *) 'stop here4'
    CALL CPU_TIME(tup)
    write(*, *) tup - tin    
        open(32, file='fort.4446')
        write(32, *) g(1: nr2, 1: nz2)
    close(32)

    CALL CPU_TIME(tin)    
    psistabr = 0.
    psistabz = 0.
    delr = 0.
    delz = 0.
    psirz(1: nr2, 1: nz2) = psiplasrz(1: nr2, 1: nz2) + psiextrz(1: nr2, 1: nz2) !total flux
    call find_new_axis_part1    

    dum1 = C_00r*S_00z - C_00z*S_00r    
    call find_fields_interp_ef_psionly(raxp + dr, zaxp, bub(1)) !give back psi, br, bz at r0, z0
    call find_fields_interp_ef_psionly(raxp - dr, zaxp, bub(2)) !give back psi, br, bz at r0, z0
    call find_fields_interp_ef_psionly(raxp, zaxp + dz, bub(3)) !give back psi, br, bz at r0, z0
    call find_fields_interp_ef_psionly(raxp, zaxp - dz, bub(4)) !give back psi, br, bz at r0, z0

    xub(1) = (bub(1) - bub(2))/(2.*dr)
    yub(1) = (bub(3) - bub(4))/(2.*dz)

    delr = (S_00z*(-xub(1)) + (-C_00z)*(-yub(1)))/dum1
    delz = (-S_00r*(-xub(1)) + C_00r*(-yub(1)))/dum1
    psistabr = delr
    psistabz = delz

    write(*, *) j_iter, rax, zax, raxp, zaxp, psistabr, psistabz

    psirz(1: nr2, 1: nz2) = psiplasrz(1: nr2, 1: nz2) + psiextrz(1: nr2, 1: nz2) + &
        psistabr*C_00(1: nr2, 1: nz2) + psistabz*S_00(1: nr2, 1: nz2) !total flux

    CALL CPU_TIME(tup)
    write(*, *) tup - tin    

    write(*, *) 'stop here41'
    CALL CPU_TIME(tin)    
    call find_new_axis_part1    
    CALL CPU_TIME(tup)
    write(*, *) tup - tin    

    write(*, *) j_iter, rax, zax, raxp, zaxp, psistabr, psistabz

    CALL CPU_TIME(tin)    
    write(*, *) 'stop here55'
    call find_psi_boundary
    CALL CPU_TIME(tup)
    write(*, *) tup - tin    

    CALL CPU_TIME(tin)    
    call new_jrz_ef  ! calculate new right hand side
    write(*, *) 'stop here7'
    CALL CPU_TIME(tup)
    write(*, *) 'current', tup - tin    

    temp_err = (abs(psistab1o - psistabr) + abs(psistab2o - psistabz))

    write(*, *) 'temp err', temp_err, psistab1o, psistab2o, psistabr, psistabz

    psistab1o = psistabr
    psistab2o = psistabz

    if (temp_err <= 1.e-9) EXIT

enddo

write(*, *) raxp, zaxp, rax, zax, psistabr, psistabz    

!check
do i=1, npassive
    curconduc(nactive+i) = psistabr*cosa(i) + psistabz*sina(i)
enddo

call psi_external_calc_ef
psirz(1: nr2, 1: nz2) = psiplasrz(1: nr2, 1: nz2) + psiextrz(1: nr2, 1: nz2)
call find_new_axis_part1    
call find_psi_boundary
call new_jrz_ef  ! calculate new right hand side

write(*, *) curconduc(1: nconduc), rax, zax    

return
end subroutine restab_axis_with_furier_wall

!------------------------------------------------------------------------
subroutine solve_gse2d_fbe_full_ef_1turn(j_init, j_stab, raxold, zaxold)    

use ef_circuit, only: nr2, nz2, r, z, rax, zax, raxp, zaxp, trax, tzax, &
    iaxis, jaxis, dr, dz, jrz, iplasma, &
    psiplasrz, psistabr, psistabz, psirz, psiextrz, derivpsi
use debugger_ef, only: markloc_ef, debug

implicit none

integer, intent(in) :: j_init, j_stab
double precision, intent(in) :: raxold, zaxold

integer :: i, j, ii
double precision :: dum1, delr, delz, tin, tup
double precision, dimension(300, 300) :: g

call markloc_ef('solve_gse2d_fbe_full_ef_1turn', debug_lev=debug)

!first, initialized initial guess coming from prescribed boundary current density:  jrhoteta
if (j_init == 0) then
write(*, *) 'reinterp curr'
    call interp_j_fromrhotorz
!rescale plasma
    dum1 = 0.
    do j = 1, nz2
    do i=1, nr2
        dum1 = dum1 + jrz(i, j)*dr*dz
    enddo
    enddo
    jrz = jrz/dum1*iplasma

    rax = raxp
    zax = zaxp
    call find_actual_index_ef(rax, zax, iaxis, jaxis)
    rax = r(iaxis)
    zax = z(jaxis)
    write(*, *) raxp, zaxp, rax, zax
endif

CALL CPU_TIME(tin)    
g=0.
call solve_gs2d(g) !jrz as right hand side
CALL CPU_TIME(tup)
write(*, *) 'solve g0', tup - tin    
write(*, *) 'stop here2'
!find boundary condition using g
CALL CPU_TIME(tin)    
call boundary_ef(g)  ! gbound = integral (Green*dg/dn) over the boundary
write(*, *) 'stop here3'
CALL CPU_TIME(tup)
write(*, *) tup - tin    

CALL CPU_TIME(tin)    
call solve_gs2d(g) ! again jrz as right hand side
psiplasrz(1: nr2, 1: nz2) = g(1: nr2, 1: nz2)
write(*, *) 'stop here4'
CALL CPU_TIME(tup)
write(*, *) tup - tin    

CALL CPU_TIME(tin)    
if (j_stab == 1) then
psistabr = 0.
psistabz = 0.
delr = 0.
delz = 0.
        psirz(1: nr2, 1: nz2) = psiplasrz(1: nr2, 1: nz2) + psiextrz(1: nr2, 1: nz2) !total flux
do ii=1, 20000
    do i=1, nr2
    do j=1, nz2
        psirz(i, j) = psirz(i, j) + delr*r(i)**2. + delz*z(j) !total flux
    enddo
    enddo
call find_new_axis_part1    
delr = -(derivpsi(3)*(raxold - rax) + derivpsi(5)*(zaxold - zax))*0.5/raxold
delz = -(derivpsi(4)*(zaxold - zax) + derivpsi(5)*(raxold - rax))
    psistabr = psistabr + delr
    psistabz = psistabz + delz

if (abs(delr) + abs(delz).lt.1.e-9) goto 2010

enddo

2010 continue

write(*, *) rax, zax, trax, tzax, raxold, zaxold, psistabr, psistabz

else
    psirz(1: nr2, 1: nz2) = psiplasrz(1: nr2, 1: nz2) + psiextrz(1: nr2, 1: nz2) !total flux
endif
CALL CPU_TIME(tup)
write(*, *) tup - tin
write(*, *) 'stop here41'
CALL CPU_TIME(tin)    
call find_new_axis_part1    
CALL CPU_TIME(tup)
write(*, *) tup - tin    

CALL CPU_TIME(tin)    
write(*, *) 'stop here55'
call find_psi_boundary
CALL CPU_TIME(tup)
write(*, *) tup - tin    

CALL CPU_TIME(tin)    
call new_jrz_ef  ! calculate new right hand side
write(*, *) 'stop here7'
CALL CPU_TIME(tup)
write(*, *) 'current', tup - tin    
!    write(*, *) psiaxis, psibnd
open(32, file='fort.444421')
    write(32, *) jrz(1: nr2, 1: nz2)
close(32)

return
end subroutine solve_gse2d_fbe_full_ef_1turn
