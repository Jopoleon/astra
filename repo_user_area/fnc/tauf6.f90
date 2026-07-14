! TAUF6 [s]: Confinement_time@(rho) for the quantity F6
!       TAUF6=Vint(F6)/(Vint(SF6TOT)-dVint(F6)/dt)
!   (Pereverzev 09-OCT-08)
double precision function TAUF6R(YR)

use scalars, only: TAU
use status, only: SF6TOT, F6, F6O
use standard_functions, only: VINT, VINTO

implicit none

double precision, intent(in) :: YR
double precision :: YQ, YW, YWO

YQ  = VINT(SF6TOT, YR)
YW  = VINT(F6, YR)
YWO = VINTO(F6O, YR)

if (YQ == 0.) then
   TAUF6R = 0.
else
   TAUF6R = TAU*YW/((YQ*TAU + YWO) - YW)
endif

end function TAUF6R
