! TENDN [keV]: Volume-density averaged electron temperature (r)
! <Te*Ne*dV>/<Ne*dV>  
!    (Polevoy 28.09.89)
double precision function TENDNR(YR)

use status, only: NE, TE
use standard_functions, only: VINT

implicit none

double precision, intent(in) :: YR

TENDNR = VINT(TE*NE, YR)/VINT(NE, YR)

end function TENDNR
