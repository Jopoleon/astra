! WCIX [MJ]: Ion thermal enegry stored in the plasma core (WIX - Pedestal)
!  0.0016 * 3/2* Integral {0, R} (Nix*Tix-Nixb*Tixb) dV
!   (Pereverzev 26-MAR-00)
double precision function WCIXR(YR)

use scalars, only: NA1
use status, only: NIX, TIX
use standard_functions, only: VINT

implicit none

double precision, intent(in) :: YR
integer :: j
double precision, dimension(NA1) :: YY

do j=1, NA1
   YY(j) = NIX(J)*TIX(J) - NIX(NA1)*TIX(NA1)
enddo
WCIXR = 0.0024*VINT(YY, YR)

end function WCIXR
