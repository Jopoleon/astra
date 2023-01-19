! IOHM [MA]: Toroidal ohmic current inside  {O, R}
!   (Pereverzev 17-FEB-00)
double precision function IOHMR(YR)

use const_inc, only: NA1, RTOR, GP2, HRO
use status_inc, only: CC, ULON, RHO

implicit none

double precision, intent(in) :: YR
double precision :: A(NA1), IINT
integer :: j

do j=1, NA1
   A(j) = CC(j)*ULON(j)/(RTOR*GP2)
   if (YR < RHO(j) + HRO) EXIT
enddo

IOHMR = IINT(A, YR)

return
end function IOHMR
