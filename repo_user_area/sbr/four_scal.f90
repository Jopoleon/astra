SUBROUTINE FOUR_SCAL(tim, n_cycle, VARIN, FREQ, n_harm, VAROUT)
!-----------------------------------------------------------------------|
! Description:  Fourier expansion in time for the function VARIN(a,t)
!               of 2 arguments
! Input:
!    N_CYCLE  - numbers of period for integration
!    VARIN(*) - periodic radial dependent function of time to be expanded
!    FREQ     - fundamental harmonic frequency [Hz]
!    N_HARM   - number of harmonics in output is 2*NUMBER+1 <= 21
!
! Output:VAROUT(ro,#) - a0(ro), a1(ro), fi1(ro), a2(ro), fi2(ro), ...
! amplitudes and phases of Fourier harmonics
! as functions of radius and #
!
! Example call from a model:FOUR_SCAL(TIME, 1, CF1, CF2, 5, CAR1):;
! The Fourier harmonics ao(r), a1(r), fi1(r),..., fi5(r) of TE(r,t)
! will be returned in the first 5 entries of CAR1
!
! Important: The subroutine can be call only once, otherwise the 
! second call will spoil the previous. 
! To avoid the interference the same subroutine can be called
! under another name (see sbr/fourex.for)
!-----------------------------------------------------------------------|

implicit none

double precision, parameter :: gp_2=6.283185d0
integer, intent(in) :: n_cycle, n_harm
double precision, intent(in) :: tim, freq, VARIN
double precision, intent(out) :: VAROUT(*)

integer :: jharm, jtime, jtimeold, jn
double precision YA(21)
double precision yarg, ydt, ytold, yyaa, yjarg, ycycle

save YA

YARG  = gp_2*FREQ*tim    ! Phase of the main frequency
YDT   = tim - YTOLD
ycycle = n_cycle
jtime = tim*FREQ/ycycle ! Integer number of completed periods 

! If the period has finished, define VAROUT and return it to ASTRA:
if(jtime .le. jtimeold) then
! Time integration step for the 0th Fourier coefficient (time average):
   YA(1) = YA(1) + VARIN*YDT
! Integration for the Fourier coefficients:
   do jharm = 1, n_harm  !jharm = armonic number
      JN = 2*jharm
      YJARG = jharm*YARG !Phase of the jharm-th harmonic
      YYAA = (VARIN - VAROUT(1))*YDT
      YA(JN)   = YA(JN)   + YYAA*cos(YJARG)
      YA(JN+1) = YA(JN+1) + YYAA*sin(YJARG)
   enddo
   jtimeold = jtime
   YTOLD = tim
   return
endif

! Division by period, return Fourier coefficients to ASTRA:

VAROUT(1) = FREQ*YA(1)/ycycle !Division by period
YA(1) = VARIN*YDT !Begin for the next integration

do jharm = 1, n_harm
   JN = 2*jharm
   YJARG = jharm*YARG
   VAROUT(JN)   = 2.*FREQ*sqrt(YA(JN)**2 + YA(JN+1)**2)/ycycle
   VAROUT(JN+1) = ATAN2(YA(JN+1), YA(JN))
! Begin for the next integration:
   YYAA = (VARIN - VAROUT(1))*YDT
   YA(JN)   = YYAA*cos(YJARG)
   YA(JN+1) = YYAA*sin(YJARG)
enddo

YTOLD = tim
jtimeold = jtime

return
END
