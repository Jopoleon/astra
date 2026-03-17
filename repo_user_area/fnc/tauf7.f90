! TAUF7 [s]: Confinement_time@(rho) for the quantity F7
!       TAUF7=Vint(F7)/(Vint(SF7TOT)-dVint(F7)/dt)
!   (Pereverzev 09-OCT-08)
double precision function TAUF7R(YR)

use scalars, only: TAU
use status, only: SF7TOT, F7, F7O
use standard_functions, only: VINT, VINTO

implicit none

double precision, intent(in) :: YR
double precision :: YQ, YW, YWO

YQ  = VINT(SF7TOT, YR)
YW  = VINT(F7, YR)
YWO = VINTO(F7O, YR)

if (YQ == 0.) then
   TAUF7R = 0.
else
   TAUF7R = TAU*YW/((YQ*TAU + YWO) - YW)
endif

end function TAUF7R
