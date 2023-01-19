! QJOUL [MW]:  Integral {0, R} (PJOUL) dV
!   (Yushmanov 11-JAN-89)
double precision FUNCTION QJOULR(YR)

use const_inc, only: NA1, GP2, RTOR
use status_inc, only: CUTOR, UPL

implicit none

double precision, intent(in) :: YR
integer :: J
double precision :: VINT, PJOUL(NA1)

do J=1, NA1
   PJOUL(j) = CUTOR(J)*UPL(J)/(GP2*RTOR)
enddo

QJOULR = VINT(PJOUL, YR)

end function QJOULR
