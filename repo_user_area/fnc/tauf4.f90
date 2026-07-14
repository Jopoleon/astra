! TAUF4 [s]: Confinement_time@(rho) for the quantity F4
!       TAUF4=Vint(F4)/(Vint(SF4TOT)-dVint(F4)/dt)
!   (Pereverzev 09-OCT-08)
double precision function TAUF4R(YR)

use scalars, only: TAU
use status, only: SF4TOT, F4, F4O
use standard_functions, only: VINT, VINTO

implicit none

double precision, intent(in) :: YR
double precision :: YQ, YW, YWO

YQ  = VINT(SF4TOT, YR)
YW  = VINT(F4, YR)
YWO = VINTO(F4O, YR)

if (YQ == 0.) then
   TAUF4R = 0.
else
   TAUF4R = TAU*YW/((YQ*TAU + YWO) - YW)
endif

end function TAUF4R
