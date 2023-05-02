! QTOK [MW]: Integral {0,R} ( PETOT+PITOT ) dV
!            (Yushmanov 11-JAN-89)
double precision FUNCTION QTOKR(YR)

use status_inc
use const_inc

implicit none

double precision, intent(in) :: yr
integer :: JK,J
	
JK = min(na1, nint(ROC/HRO))

QTOKR=0.
do J=1, JK
    QTOKR = QTOKR + (PE(J) + PI(J))*VR(J)
enddo
QTOKR = HRO*(QTOKR - (PE(JK) + PI(JK))*HRO)

return
end function QTOKR


