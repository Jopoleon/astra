! M. Bergmann 20.01.2026
!
! ER_omp calculates the Er profiles usinf outboard midplane (omp) quantities
! The profiles should be closer to the experimental profiles and should
! be the one used for TGLF/GK simulations
!OUTPUTS:
! ->wexb_tglf MUST BE PASSED TO ER BUT IS NOT THE RADIAL ELECTRIC FIELD
!   IT IS ER adjusted to have local quantities in TGLF. It's in V/m
! To be passed to vper_m = ER/(bpolz*RTOR) in tglf_interf.f

! ->er_lfs is the ER profile in V/m, to use for plots etc.
! ->vdia_lfs is the diamagneic velocity in m/s
! ->bpol_lfs is the local Bpol at the lfs from SPIDER 

! INPUTS:
! ->er_min is the location of the Er minimum in rho_pol
! ->er_sep is the value to impose to Er at rho=1

! GEOM2D extracts directly from SPIDER: GEOM2D(1:NEQUIL,1:MEQUIL,x)
! NEQUIL, MEQUIL -> radial and poloidal grid from SPIDER      
! x=1->PSI, x=2->theta, x=3->R, x=4->Z, x=5->Bpol, x=6->Btot, x=7->r
!

subroutine er_omp(er_min, er_sep, wexb_lfs, er_lfs, vdia_lfs, bp_lfs)

use scalars, only: RTOR, BTOR, NA1, TIME, TSTART, NEQUIL, MEQUIL, AWALL
use status, only: TI, NMAIN, ZMAIN, VTOR, AMETR, MU, rho_pol, VPOL
use parameters_a2equil, only: equil_now
use numerical_tools, only: qinterp

implicit none

double precision, intent(in) :: er_min, er_sep
double precision, dimension(na1), intent(out) :: wexb_lfs, er_lfs, vdia_lfs, bp_lfs

integer :: ispan, jrho, neq, meq
double precision :: er0
double precision, allocatable, dimension(:) :: rmin_sp, rmaj_sp, bp_sp, psi_sp, rpol_sp, zispan
double precision, dimension(na1) :: rmaj_as, rmin_as, bp_as, vdia_as, pi_as, psi_as
double precision, dimension(na1-1) :: rhalf, vdhalf

!TSTART impose 0 otherwise might have NaN in TGLF
if (TIME <= TSTART + 0.0001) then
    write(6,*) 'Wexb TGLF', TIME, TSTART + 0.001
    er_lfs(1:na1) = 0.0001
    wexb_lfs(1:na1) = 0.0001
    return
endif
!NEQUIL and MEQUIL are NOT integers!
neq = int(nequil)
meq = int(mequil)
allocate(rmin_sp(neq))
allocate(rmaj_sp(neq))
allocate(bp_sp(neq))
allocate(psi_sp(neq))
allocate(rpol_sp(neq))
allocate(zispan(neq))

ispan = minloc(abs(equil_now%coord_sys%position%r(neq,1:meq) - (RTOR + AWALL)), 1) ! Z=0, lfs
psi_sp(1:neq) = equil_now%coord_sys%position%psirz(1:neq, ispan) !PSI
bp_sp(1:neq)  = equil_now%coord_sys%bpcell(1:nequil, ispan) !BPOL
bp_sp(neq) = bp_sp(neq-1) !defined up to neq-1
rmaj_sp(1:neq) = equil_now%coord_sys%position%r(1:neq, ispan)
rmin_sp(1:neq) = equil_now%coord_sys%position%r(1:neq, ispan) - (rmaj_sp(1))
zispan(1:neq)  = equil_now%coord_sys%position%z(1:neq, ispan)
rpol_sp(1:neq) = sqrt((psi_sp(1:neq) - psi_sp(1))/(psi_sp(neq) - psi_sp(1))); !rpol spider

!interpolate on astra grid using rpol
 !pi_as(1:na1) = TI(1:na1)*NI(1:na1) !ion pressure
pi_as(1:na1) = TI(1:na1)*NMAIN(1:na1) !ion pressure of main ion species

call qinterp(rpol_sp(1:neq), rmaj_sp(1:neq), neq, rho_pol(1:NA1), rmaj_as(1:NA1), NA1) !r
call qinterp(rpol_sp(1:neq), rmin_sp(1:neq), neq, rho_pol(1:NA1), rmin_as(1:NA1), NA1) !r (defined but not used in original) 
call qinterp(rpol_sp(1:neq), bp_sp(  1:neq), neq, rho_pol(1:NA1), bp_as(1:NA1), NA1)  !Bpol

bp_lfs(1:na1) = bp_as(1:na1)

call qinterp(rpol_sp(1:neq), psi_sp(1:neq), neq, rho_pol(1:NA1), psi_as(1:NA1), NA1)  !Psi

!diamagnetic term = R*B_pol/qi/ni*grad_Psi(Pi) (V/m)
do jrho=1,NA1-1
    rhalf(jrho) = (rho_pol(jrho) + rho_pol(jrho+1))/2.0
    vdhalf(jrho) = 2.0E3/(NMAIN(jrho+1) + NMAIN(jrho))/ZMAIN(jrho) * &
        (pi_as(jrho+1) - pi_as(jrho))/(psi_as(jrho+1) - psi_as(jrho)) * &
        (bp_as(jrho+1) + bp_as(jrho))/2.0 * (rmaj_as(jrho+1) + rmaj_as(jrho))/2.0 
enddo

call qinterp(rhalf(1:NA1-1), vdhalf(1:NA1-1), NA1-1, rho_pol(1:na1), vdia_as(1:na1), na1)

!vdia in m/s
vdia_lfs(1:na1) = vdia_as(1:na1)/(BTOR*RTOR/rmaj_as(1:na1)) !(m/s)

!local Er at the LFS
er_lfs(1:na1) = vdia_as(1:na1) - VPOL(1:na1)*BTOR*RTOR/rmaj_as(1:na1) + bp_as(1:na1)*VTOR(1:na1) !V/m

!Impose Er well (Er=er_sep at the separatrix)
ispan = 0
do jrho=1, na1
    if (rho_pol(jrho) >= er_min) then      
        if (ispan == 0)then
            ispan = 1
            er0 = er_lfs(jrho)
        endif
        er_lfs(jrho) = er0 - (er0 - er_sep)*((rho_pol(jrho) - er_min)/(1. - er_min))**2
    endif
enddo

!WExB in 1/s to pass to TGLF
!RTOR*BPOLZ = RTOR*(BTOR*AMETR(jrho)*MU(jrho)/RTOR)

wexb_lfs(1:na1) = er_lfs(1:na1)/(bp_as(1:na1)*rmaj_as(1:na1))*BTOR*AMETR(1:na1)*MU(1:na1) !ER in V/m to pass to TGLF

ispan = 0
do jrho=NA1, 1,-1 !this is just to avoid zigzags in the inner core
    if (rho_pol(jrho) <= 0.3) then
        if (ispan == 0)then
            ispan = 1
            er0 = wexb_lfs(jrho)
        endif
        wexb_lfs(jrho) = er0 - er0*(1. - rho_pol(jrho)/0.3)**2
    endif
enddo

deallocate(rmin_sp)
deallocate(rmaj_sp)
deallocate(bp_sp)
deallocate(psi_sp)
deallocate(rpol_sp)
deallocate(zispan)

end subroutine er_omp
