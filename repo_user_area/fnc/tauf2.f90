! TAUF2 [s]: Confinement_time@(rho) for the quantity F2
!       TAUF2=Vint(F2)/(Vint(SF2TOT)-dVint(F2)/dt)
!   (Pereverzev 09-OCT-08)
double precision function TAUF2R(YR)

use const_inc, only: TAU
use status_inc, only: SF2TOT, F2, F2O

implicit none

double precision, intent(in) :: YR
double precision :: VINT, VINTO, YQ, YW, YWO

YQ  = VINT(SF2TOT, YR)
YW  = VINT(F2, YR)
YWO = VINTO(F2O, YR)

if (YQ == 0.) then
   TAUF2R = 0.
else
   TAUF2R = TAU*YW/((YQ*TAU + YWO) - YW)
endif

end function TAUF2R
