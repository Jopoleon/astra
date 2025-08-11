subroutine neo_ipc(nworkers, rho_bnd)

use io_mod, only: equ_file, exp_file
use ipc_mod, only: n_sbp_arr_in, n_sbp_arr_out, mem_neo
use const_inc, only: NA1, DEVAR, n_bouncon
use debugger, only: markloc

implicit none

integer, intent(in) :: nworkers
double precision, intent(in) :: rho_bnd

integer, dimension(nworkers) :: SBP_JBEG, SBP_JEND
character(len=64), dimension(nworkers) :: SBP_NAMES

logical :: first_call=.True.
integer :: j, t_wall1, t_wall2, rate, nrho_step

save first_call

SBP_NAMES = "xpr/neo"//char(0)

nrho_step = int(rho_bnd*NA1/dble(nworkers)) + 1
SBP_JBEG = (/ (nrho_step*j + 1, j=0, nworkers-1) /)
SBP_JEND(1:nworkers-1) = SBP_JBEG(2:nworkers) - 1
SBP_JEND(nworkers) = int(rho_bnd*NA1) + 1
print*, first_call, nworkers
if (first_call) then
    print*, SBP_JBEG
    print*, SBP_JEND
    call markloc("initialise_ipc")
    call initialise_ipc(NA1, n_sbp_arr_in, n_sbp_arr_out, nworkers, equ_file, exp_file)
    call markloc("send_ipc_jobs")
    call send_ipc_jobs(nworkers, 64, SBP_NAMES, SBP_JBEG, SBP_JEND)
    first_call = .False.
endif

! **** Fill shared memory segments
call SYSTEM_CLOCK(t_wall1, rate)
call markloc("setvars")
call setvars(DEVAR, n_bouncon)
call markloc("set_sbp_input_arrays")
call set_sbp_input_arrays

! **** Free each semaphore
do j=1, nworkers
   call unlock_sbp(j)
enddo

! **** Synchronisation point
call wait4all

! **** Collect data from ShMem
do j=1, nworkers
    call sbp2astra(SBP_JBEG(j), SBP_JEND(j), j, mem_neo(1, 1))
enddo

call SYSTEM_CLOCK(t_wall2, rate)
print*, "XPR wall time", dble(t_wall2 - t_wall1)/dble(rate)

return
end subroutine neo_ipc
