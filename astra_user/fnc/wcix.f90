! WCIX [MJ]: Ion thermal enegry stored in the plasma core (WIX - Pedestal)
!  0.0016 * 3/2* Integral {0, R} (Nix*Tix-Nixb*Tixb) dV
!   (Pereverzev 26-MAR-00)
double precision function WCIXR(YR)

use parameter_inc, only: NRD
use const_inc, only: NA1
use status_inc, only: NIX, TIX

implicit none

double precision, intent(in) :: YR
integer  J
double precision :: VINT, YY(NRD)

do j=1, NA1
   YY(j) = NIX(J)*TIX(J) - NIX(NA1)*TIX(NA1)
enddo
WCIXR = 0.0024*VINT(YY, YR)

end function WCIXR
