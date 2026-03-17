! TAUF1 [s]: Confinement_time@(rho) for the quantity F1
!       TAUF1=Vint(F1)/(Vint(SF1TOT)-dVint(F1)/dt)
!   (Pereverzev 09-OCT-08)
double precision function TAUF1R(YR)

use scalars, only: TAU
use status, only: SF1TOT, F1, F1O
use standard_functions, only: VINT, VINTO

implicit none

double precision, intent(in) :: YR
double precision :: YQ, YW, YWO

YQ  = VINT(SF1TOT, YR)
YW  = VINT(F1, YR)
YWO = VINTO(F1O, YR)

if (YQ == 0.) then
   TAUF1R = 0.
else
   TAUF1R = TAU*YW/((YQ*TAU + YWO) - YW)
endif

end function TAUF1R
