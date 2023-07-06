! NELA [10#19/m#3]: Density Line Average in the mid-plane (r) [m]
!	Integral {0,r} ( NE ) dl / a
!			(Pereverzev 22-12.97)
! 
double precision FUNCTION NELAR(YD)

use const_inc, only: NA1, HRO, ABC
use status_inc, only: AMETR, NE

implicit none

double precision, intent(in) :: YD
integer :: J, JK
double precision :: YX, YXO, YDEL, YR

if (YD >= ABC) then
    NELAR = 0.
    return
endif

JK = 0
YXO = 0.

do J=1, NA1
    YR = AMETR(j)
    if (J == NA1) YR = ABC
    if (YR <= YD) CYCLE
    YX = SQRT(YR**2 - YD**2)
    if (JK /= 0) then
        NELAR = NELAR + (YX - YXO)*(NE(J-1) + NE(J))
    else
        JK = 1
        if (J == 1) then
            NELAR = (YX - YXO)*2.*NE(1)
        else
            YDEL = (YR - YD)/HRO
            NELAR = (YX - YXO)*(YDEL*NE(J-1) + (2. - YDEL)*NE(J))
        endif
    endif
    YXO = YX
enddo

NELAR = NELAR*0.5/(ABC)

end function NELAR
