! IOHM [MA]: Toroidal ohmic current inside  {O, R}
!   (Pereverzev 17-FEB-00)
double precision function IOHMR(YR)

use pi_const, only: GP2
use scalars, only: NA1, RTOR, HRO
use status, only: CC, ULON, RHO
use standard_functions, only: IINT

implicit none

double precision, intent(in) :: YR
double precision :: A(NA1)
integer :: j

do j=1, NA1
   A(j) = CC(j)*ULON(j)/(RTOR*GP2)
   if (YR < RHO(j) + HRO) EXIT
enddo

IOHMR = IINT(A, YR)

end function IOHMR
