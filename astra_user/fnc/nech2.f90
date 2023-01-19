! NECH2 [10#19/m#3]: Chord average density (r) [m]
! Integral {0, r} (NE) dl / a
!   (Yushmanov 11-MAY-87)

double precision function NECH2R(r_in)

use const_inc, only: HRO, NA, NA1
use status_inc, only: NE

implicit none

double precision, intent(in) :: r_in
integer :: J, JK, jbeg
double precision :: YR, YXO, YX, YDEL

if(r_in >= HRO*NA) then
   NECH2R = 0.
   return
endif

! Look for closest radius <= r_in
do J=1, NA1
   YR = HRO*(J - 0.5)
   if(J == NA1) YR = HRO*NA
   if(YR > r_in) EXIT
enddo

! Initial value of the integral
YXO = 0.
jbeg = J - 1
IF(jbeg == 1) then
   NECH2R = (YX - YXO)*2.*NE(1)
else
   YR = HRO*(jbeg - 0.5)
   YX = SQRT(YR**2 - r_in**2)
   YDEL = (YR - r_in)/HRO
   NECH2R = (YX - YXO)*(YDEL*NE(jbeg-1) + (2. - YDEL)*NE(jbeg))
   YXO = YX
endif

do J=jbeg+1, NA1
   YR = HRO*(J - 0.5)
   YX = SQRT(YR**2 - r_in**2)
   NECH2R = NECH2R + (YX - YXO)*(NE(J-1) + NE(J))
   YXO = YX
enddo

NECH2R = NECH2R*0.5/(HRO*NA1)

return
end function NECH2R
