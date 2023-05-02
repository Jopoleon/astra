! QDT [MW]: Integral {0,R} ( PDT ) dV
!           (Yushmanov 11-JAN-89)
double precision FUNCTION QDTR(YR)

use status_inc
use const_inc

implicit none

double precision, intent(in) :: YR
integer :: j,jk
double precision :: PDT, SVDT

!rho units
JK = nint(YR/HRO)

QDTR = 0.
do J=1, JK
    include 'fml/pdt'
    QDTR = QDTR + PDT*VR(J)
enddo
QDTR = HRO*(QDTR - PDT*HRO)

return
end function QDTR
