! WCX [MJ]: Exp. thermal enegry stored in the plasma core (WTOTX - Pedestal)
! 0.0016*3/2*Integral {0, R} (Nix*Tix+Nex*Tex-Nixb*Tixb-Nexb*Texb)dV
!   (Pereverzev 26-MAR-00)
double precision function WCXR(YR)

use parameter_inc, only: NRD
use const_inc, only: NA1
use status_inc, only: NIX, TIX, NEX, TEX

implicit none

double precision, intent(in) :: YR
integer j
double precision :: VINT, YY(NRD)

do j = 1, NA1
   YY(j) = NIX(J)*TIX(J) - NIX(NA1)*TIX(NA1) + NEX(J)*TEX(J) - NEX(NA1)*TEX(NA1)
enddo
WCXR = 0.0024*VINT(YY, YR)

end function WCXR
