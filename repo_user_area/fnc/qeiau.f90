! QEIAU [MW]:  Integral {0,R} ( PEIAU ) dV
!              (Yushmanov 11-JAN-89)
double precision FUNCTION QEIAUR(YR)

use scalars, only: HRO, AIM1, AIM2, AIM3
use status, only: NE, TE, TI, VR, ZIM1, ZIM2, ZIM3, &
    F1, F2, F3, F4, F5, F6, F7, F8, F9 

implicit none

double precision, intent(in) :: YR
integer :: JK, J
double precision :: PEIAU, COULG, ZINE, ZIKR, ZIAR, T, Z

JK = nint(YR/HRO)
QEIAUR = 0.
do J=1, JK
    include 'fml/peiau'
    QEIAUR = QEIAUR + PEIAU*(TE(J) - TI(J))*VR(J)*HRO
enddo

end function QEIAUR
