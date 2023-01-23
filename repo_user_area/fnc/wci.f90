! WCI [MJ]: Ion thermal enegry stored in the plasma core (WI - Pedestal)
!  0.0016 * 3/2* Integral {0, R} (Ni*Ti-Nib*Tib) dV
!   (Pereverzev 26-MAR-00)
double precision function WCIR(YR)

use parameter_inc, only: NRD
use const_inc, only: NA1
use status_inc, only: NI, TI

implicit none

double precision, intent(in) :: YR
integer  J
double precision :: VINT, YY(NRD)

do j=1, NA1
   YY(j) = NI(J)*TI(J) - NI(NA1)*TI(NA1)
enddo
WCIR = 0.0024*VINT(YY, YR)

end function WCIR
