! QEIAD [MW]:  Integral {0,R} ( PEIAD ) dV
!              (Yushmanov 11-JAN-89)
double precision FUNCTION QEIADR(YR)

use const_inc, only: HRO, AIM1, AIM2, AIM3
use status_inc, only: NE, TE, TI, VR, ZIM1, ZIM2, ZIM3, &
    F1, F2, F3, F4, F5, F6, F7, F8, F9 

implicit none

double precision, intent(in) :: YR
integer :: JK, J
double precision :: PEIAD, COULG, ZIWOL, T, Z

JK = nint(YR/HRO)
QEIADR = 0.
do J=1, JK
    include 'fml/peiad'
    QEIADR = QEIADR + PEIAD*(TE(J) - TI(J))*VR(J)*HRO
enddo

return
end function QEIADR
