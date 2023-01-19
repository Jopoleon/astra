! TAUF4 [s]: Confinement_time@(rho) for the quantity F4
!       TAUF4=Vint(F4)/(Vint(SF4TOT)-dVint(F4)/dt)
!   (Pereverzev 09-OCT-08)
double precision function TAUF4R(YR)

use const_inc, only: TAU
use status_inc, only: SF4TOT, F4, F4O

implicit none

double precision, intent(in) :: YR
double precision :: VINT, VINTO, YQ, YW, YWO

YQ  = VINT(SF4TOT, YR)
YW  = VINT(F4, YR)
YWO = VINTO(F4O, YR)

if (YQ == 0.) then
   TAUF4R = 0.
else
   TAUF4R = TAU*YW/((YQ*TAU + YWO) - YW)
endif

end function TAUF4R
