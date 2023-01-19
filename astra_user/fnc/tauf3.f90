! TAUF3 [s]: Confinement_time@(rho) for the quantity F3
!       TAUF3=Vint(F3)/(Vint(SF3TOT)-dVint(F3)/dt)
!   (Pereverzev 09-OCT-08)
double precision function TAUF3R(YR)

use const_inc, only: TAU
use status_inc, only: SF3TOT, F3, F3O

implicit none

double precision, intent(in) :: YR
double precision :: VINT, VINTO, YQ, YW, YWO

YQ  = VINT(SF3TOT, YR)
YW  = VINT(F3, YR)
YWO = VINTO(F3O, YR)

if (YQ == 0.) then
   TAUF3R = 0.
else
   TAUF3R = TAU*YW/((YQ*TAU + YWO) - YW)
endif

end function TAUF3R
