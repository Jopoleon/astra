! TAUF5 [s]: Confinement_time@(rho) for the quantity F5
!       TAUF5=Vint(F5)/(Vint(SF5TOT)-dVint(F5)/dt)
!   (Pereverzev 09-OCT-08)
double precision function TAUF5R(YR)

use const_inc, only: TAU
use status_inc, only: SF5TOT, F5, F5O

implicit none

double precision, intent(in) :: YR
double precision :: VINT, VINTO, YQ, YW, YWO

YQ  = VINT(SF5TOT, YR)
YW  = VINT(F5, YR)
YWO = VINTO(F5O, YR)

if (YQ == 0.) then
   TAUF5R = 0.
else
   TAUF5R = TAU*YW/((YQ*TAU + YWO) - YW)
endif

end function TAUF5R
