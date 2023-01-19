! TINDN [keV]: Volume-density averaged electron temperature (r)
! <Ti*Ne*dV>/<Ne*dV>  
!    (Polevoy 28.09.89)
double precision function TINDNR(YR)

use status_inc, only: NE, TI
implicit none

double precision, intent(in) :: YR
double precision VINT

TINDNR = VINT(TI*NE, YR)/VINT(NE, YR)

end function TINDNR
