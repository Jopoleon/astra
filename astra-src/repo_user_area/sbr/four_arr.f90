SUBROUTINE FOUR_ARR(tim, n_cycle, ARRIN, FREQ, n_harm, arrout)
!-----------------------------------------------------------------------|
! Description:  Fourier expansion in time for the function ARRIN(a,t)
!               of 2 arguments
! Input:
!    N_CYCLE  - numbers of period for integration
!    ARRIN(*) - periodic radial dependent function of time to be expanded
!    FREQ     - fundamental harmonic frequency [Hz]
!    N_HARM   - number of harmonics in output is 2*NUMBER+1 <= 21
!
! Output:ARROUT(ro,#) - a0(ro), a1(ro), fi1(ro), a2(ro), fi2(ro), ...
! amplitudes and phases of Fourier harmonics
! as functions of radius and #
!
! Example call from a model:FOUR_ARR(TIME, 3.d0, 1, TE, CF2, 5, CAR6):;
! The Fourier harmonics ao(r), a1(r), fi1(r),..., fi5(r) of TE(r,t)
! will be returned in CAR6, CAR7, CAR8,..., CAR16 respectively.
!
! Important: The subroutine can be call only once, otherwise the 
! second call will spoil the previous. 
! To avoid the interference the same subroutine can be called
! under another name (see sbr/fourex.for)
!-----------------------------------------------------------------------|

use pi_const, only: GP2
use status, only: NRD
use scalars, only: NA1

implicit none

integer, intent(in) :: n_cycle, n_harm
double precision, intent(in) :: tim, freq
double precision, intent(in), dimension(na1) :: ARRIN
double precision, intent(out), dimension(NRD, n_harm*2+1) :: ARROUT

integer :: jharm, jrho, jtime, jtimeold, jn
double precision, dimension(NA1, 11) :: prof_sum
double precision :: arg, dt, tim_old=0., yyaa, phase, ycycle

ARG  = GP2*FREQ*tim    ! Phase of the main frequency
DT   = tim - TIM_OLD
ycycle = n_cycle
jtime    = tim    *FREQ/ycycle ! Integer number of completed periods 
jtimeold = tim_old*FREQ/ycycle

prof_sum = 0.d0

if (jtime <= jtimeold) then
! Time integration step for the 0th Fourier coefficient (time average):
    prof_sum(:, 1) = prof_sum(:, 1) + ARRIN(1: NA1)*DT
! Integration for the Fourier coefficients:
    do jharm = 1, n_harm  !jharm = armonic number
        JN = 2*jharm
        PHASE = jharm*ARG !Phase of the jharm-th harmonic
        do jrho = 1, NA1
            YYAA = (ARRIN(jrho) - ARROUT(jrho, 1))*DT
            prof_sum(jrho, JN)   = prof_sum(jrho, JN)   + YYAA*cos(PHASE)
            prof_sum(jrho, JN+1) = prof_sum(jrho, JN+1) + YYAA*sin(PHASE)
        enddo
    enddo
    TIM_OLD = tim
    return
endif

! If the period has finished, define ARROUT and return it to ASTRA:
! Division by period, return Fourier coefficients to ASTRA:

ARROUT(1: NA1, 1) = FREQ*prof_sum(:, 1)/ycycle !Division by period
prof_sum(:, 1) = ARRIN(1: NA1)*DT ! Reset for the next integration

do jharm = 1, n_harm
    JN = 2*jharm
    PHASE = jharm*ARG
    do jrho = 1, NA1
        ARROUT(jrho, JN)   = 2.*FREQ*sqrt(prof_sum(jrho, JN)**2 + prof_sum(jrho, JN+1)**2)/ycycle
        ARROUT(jrho, JN+1) = ATAN2(prof_sum(jrho, JN+1), prof_sum(jrho,JN))
! Begin for the next integration:
        YYAA = (ARRIN(jrho) - ARROUT(jrho, 1))*DT
        prof_sum(jrho, JN)   = YYAA*cos(PHASE)
        prof_sum(jrho, JN+1) = YYAA*sin(PHASE)
    enddo
enddo

TIM_OLD = tim

end subroutine four_arr
