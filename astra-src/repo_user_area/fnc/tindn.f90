! TINDN [keV]: Volume-density averaged electron temperature (r)
! <Ti*Ne*dV>/<Ne*dV>  
!    (Polevoy 28.09.89)
double precision function TINDNR(YR)

use status, only: NE, TI
use standard_functions, only: VINT

implicit none

double precision, intent(in) :: YR

TINDNR = VINT(TI*NE, YR)/VINT(NE, YR)

end function TINDNR
