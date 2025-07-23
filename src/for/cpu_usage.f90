module cpu_usage

use parameter_inc, only: NSDELOUT, NSBMX
use outcmn_inc, only: NSBR, DTNAME, IFSBX, tab_ch
use const_inc, only: NSTEPS, TIME, TSTART
use debugger, only: markloc

implicit none

double precision :: cpu_start, cpuTime_equ=0., cpuTime_sbr(NSBMX)=0.
integer :: wall_start

contains

!---------------------------------------------------------------------
  subroutine CPU_report(str_in)

  integer, parameter :: nch=6

  character(len=*), intent(in) :: str_in

  integer :: j, j1, j2, wall_now, rate
  double precision :: Y, cpu_now, cpuTime_tot
  real :: wall_tot

  call markloc('CPU_usage')
  write(nch, '(A)') TRIM(str_in)
  call SYSTEM_CLOCK(wall_now, rate)
  wall_tot = real(wall_now - wall_start) / real(rate)
  call formatTime(nch, '  >>> Astra wall time  ', wall_tot, -1.)

  write(nch, '(A, I8)')    "    Total time steps  ", NSTEPS
  if (NSTEPS == 0) return
  Y = (TIME - TSTART)/NSTEPS
  if (Y < 1.d-1) then
      Y = 1.d3*Y
      write(nch, '(A, F6.3, A)')"    Average time step   ", Y, " msec"
  else
      write(nch, '(A, F6.3, A)')"    Average time step   ", Y, " sec"
  endif
  call CPU_TIME(cpu_now)
  cpuTime_tot = cpu_now - cpu_start
  call formatTime(nch, '    Total CPU time' // char(0), cpuTime_tot, cpuTime_tot)
  call formatTime(nch, '    Equilibrium   ' // char(0), cpuTime_equ, cpuTime_tot)
  j2 = 1
  do j1=1, NSBR
      j = min(6, LEN_TRIM(DTNAME(NSDELOUT+4*j1)))
      if (j1 == IFSBX(j2)) then
          call formatTime(nch, '    Xroutine   "' // &
              DTNAME(NSDELOUT+4*j1)(1: j) // '"', cpuTime_sbr(j1), cpuTime_tot)
          j2 = j2 + 1
      else
          call formatTime(nch, '    Subroutine "' // &
              DTNAME(NSDELOUT+4*j1)(1: j) // '"', cpuTime_sbr(j1), cpuTime_tot)
      endif
  enddo
  write(nch, *)

  return
  end subroutine CPU_report

!---------------------------------------------------------------------
  subroutine formatTime(nch, string, tim, time)

  integer, intent(in) :: nch
  character(len=*), intent(in) :: string
  double precision, intent(in) :: tim, time

  integer :: jh, jm, js
  double precision :: t1

  call markloc('formatTime')

  js = tim
  jh = js/3600
  jm = (js - 3600*jh)/60
  t1 = tim - 60*jm - 3600*jh
  js = t1

  if (time >= 0.) then
      write(nch, '(2A, I4.2, 2(A1, I2.2), F8.1, A1)') TRIM(string), &
          tab_ch, jh, ':', jm, ':', js, 100.*tim/time, '%'
  else
      write(nch, '(2A, I4.2, 2(A1, I2.2))') TRIM(string), &
         tab_ch, jh, ':', jm, ':', js
  endif

  return
  end subroutine formatTime

endmodule cpu_usage
