! TENDN [keV]: Volume-density averaged electron temperature (r)
! <Te*Ne*dV>/<Ne*dV>  
!    (Polevoy 28.09.89)
double precision function TENDNR(YR)

use status_inc, only: NE, TE
implicit none

double precision, intent(in) :: YR
double precision VINT

TENDNR = VINT(TE*NE, YR)/VINT(NE, YR)

end function TENDNR
