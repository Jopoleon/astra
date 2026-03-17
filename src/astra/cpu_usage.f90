module cpu_usage

use io_mod, only: n_sbr, sbr_name
use scalars, only: NSTEPS, TIME, TSTART
use debugger, only: markloc

implicit none

integer :: wall_start, wallTime_equ=0
integer, allocatable :: wallTime_sbr(:)
double precision :: cpu_start, cpuTime_equ=0.
double precision, allocatable :: cpuTime_sbr(:)

contains

!---------------------------------------------------------------------
    subroutine CPU_init

    allocate(wallTime_sbr(n_sbr), cpuTime_sbr(n_sbr))
    cpuTime_sbr = 0.
    wallTime_sbr = 0

    end subroutine CPU_init

!---------------------------------------------------------------------
    subroutine CPU_report(str_in)

    integer, parameter :: nch=6

    character(len=*), intent(in) :: str_in

    integer :: j, j2, wall_now, rate
    double precision :: time_step, cpu_now, cpuTime_tot
    real :: wall_tot

    call markloc('CPU_usage')
    write(nch, '(A)') TRIM(str_in)
    call SYSTEM_CLOCK(wall_now, rate)
    wall_tot = real(wall_now - wall_start) / real(rate)
    call CPU_TIME(cpu_now)
    cpuTime_tot = cpu_now - cpu_start

    call formatTime(nch, 'Total CPU time',   cpuTime_tot, -1.d0)
    write(nch, '(4X, A, I8)') 'Total time steps  ', NSTEPS
    if (NSTEPS == 0) return
    time_step = (TIME - TSTART)/NSTEPS
    if (time_step < 1.d-1) then
        write(nch, '(4X, A, F6.3, A)') 'Average time step   ', 1.d3*time_step, ' msec'
    else
        write(nch, '(4X, A, F6.3, A)') 'Average time step   ', time_step, ' sec'
    endif
    call formatTime(nch, '>>> Astra wall time', wall_tot, -1.d0)
    call formatTime(nch, 'Equilibrium', dble(wallTime_equ)/dble(rate), wall_tot)
    j2 = 1
    do j=1, n_sbr
        call formatTime(nch, 'Subroutine ' // sbr_name(j), dble(wallTime_sbr(j))/dble(rate), wall_tot)
    enddo
    write(nch, *)

    end subroutine CPU_report

!---------------------------------------------------------------------
    subroutine formatTime(nch, string_in, tim, time)

    integer, intent(in) :: nch
    character(len=*), intent(in) :: string_in
    double precision, intent(in) :: tim, time

    integer :: jh, jm, js
    double precision :: ts
    character(len=30) :: padded

    call markloc('formatTime')

    js = tim
    jh = js/3600
    jm = (js - 3600*jh)/60
    ts = tim - 60*jm - 3600*jh
    js = ts

    padded = ADJUSTL(string_in)
    if (time >= 0.) then
        write(nch, '(4X, A30, I4.2, 2(A1, I2.2), F8.1, A1)') padded, &
            jh, ':', jm, ':', js, 100.*tim/time, '%'
    else
        write(nch, '(4X, A30, I4.2, 2(A1, I2.2))') padded, &
            jh, ':', jm, ':', js
    endif

    end subroutine formatTime

end module cpu_usage
