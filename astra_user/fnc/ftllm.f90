! FTLLMR [%]: 
! Effective trapped particle fraction Lin-Liu and Miller GA-A21820
!                                                        (Oct 94)
! Y.R.Lin-Liu and R.L.Miller,  Phys.Plasmas 2(5),  May 1995,  pp.1666-1668
double precision function FTLLMR(YR)

use const_inc, only: NA, BTOR, HRO
use status_inc, only: BDB0, BDB02, FOFB, BMAXT

implicit none

double precision YR, YYR, YH, YFTUP, YFTLO
integer JK

FTLLMR = 0.1
IF (YR <= 0.) RETURN
YYR = YR
JK = YYR/HRO + 1
IF (JK > NA) THEN
   YYR = HRO*NA
   JK = NA
ENDIF
IF (ABS(BMAXT(JK)) < 0.0001 .OR. ABS(BDB0(JK)) < 0.0001) THEN
   FTLLMR = 0.1
else
   YH = min(.999999d0, BDB0(JK)/BMAXT(JK)*BTOR)
   YFTUP = 1. - BDB02(JK)/BDB0(JK)**2 * (1. - SQRT(1. - YH)*(1. + .5*YH))
   YFTLO = 1. - BDB02(JK)*FOFB(JK)
   FTLLMR = .75*YFTUP + .25*YFTLO
endif

return
end function FTLLMR

