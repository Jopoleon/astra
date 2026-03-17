! QJOUL [MW]:  Integral {0, R} (PJOUL) dV
!   (Yushmanov 11-JAN-89)
double precision FUNCTION QJOULR(YR)

use pi_const, only: GP2
use scalars, only: NA1, RTOR
use status, only: CUTOR, UPL
use standard_functions, only: VINT

implicit none

double precision, intent(in) :: YR
integer :: j 
double precision :: PJOUL(NA1)

do j=1, NA1
    PJOUL(j) = CUTOR(J)*UPL(J)/(GP2*RTOR)
enddo
QJOULR = VINT(PJOUL, YR)

end function QJOULR
