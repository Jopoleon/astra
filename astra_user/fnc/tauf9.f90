! TAUF9 [s]: Confinement_time@(rho) for the quantity F9
!       TAUF9=Vint(F9)/(Vint(SF9TOT)-dVint(F9)/dt)
!   (Pereverzev 09-OCT-08)
double precision function TAUF9R(YR)

use const_inc, only: TAU
use status_inc, only: SF9TOT, F9, F9O

implicit none

double precision, intent(in) :: YR
double precision :: VINT, VINTO, YQ, YW, YWO

YQ  = VINT(SF9TOT, YR)
YW  = VINT(F9, YR)
YWO = VINTO(F9O, YR)

if (YQ == 0.) then
   TAUF9R = 0.
else
   TAUF9R = TAU*YW/((YQ*TAU + YWO) - YW)
endif

end function TAUF9R
