! QDT [MW]: Integral {0,R} ( PDT ) dV
!           (Yushmanov 11-JAN-89)
double precision FUNCTION QDTR(YR)

use status, only: VR, TI, NDEUT, NTRIT
use scalars, only: HRO

implicit none

double precision, intent(in) :: YR
integer :: j, jk
double precision :: PDT, SVDT

!rho units
JK = nint(YR/HRO)

QDTR = 0.
do J=1, JK
    SVDT = TI(J) ** (-0.33333333)
    SVDT = 8.972*EXP(-19.9826*SVDT)*SVDT**2 * &
        ((TI(J) + 1.0134)/(1. + 6.386E-3*(TI(J) + 1.0134)**2) + 1.877*EXP(-0.16176*TI(J)*SQRT(TI(J))))
    PDT = 5.632*NDEUT(J)*NTRIT(J)*SVDT
    QDTR = QDTR + PDT*VR(J)
enddo
QDTR = HRO*(QDTR - PDT*HRO)

end function QDTR
