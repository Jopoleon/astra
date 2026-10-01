! TAUF8 [s]: Confinement_time@(rho) for the quantity F8
!       TAUF8=Vint(F8)/(Vint(SF8TOT)-dVint(F8)/dt)
!   (Pereverzev 09-OCT-08)
double precision function TAUF8R(YR)

use scalars, only: TAU
use status, only: SF8TOT, F8, F8O
use standard_functions, only: VINT, VINTO

implicit none

double precision, intent(in) :: YR
double precision :: YQ, YW, YWO

YQ  = VINT(SF8TOT, YR)
YW  = VINT(F8, YR)
YWO = VINTO(F8O, YR)

if (YQ == 0.) then
   TAUF8R = 0.
else
   TAUF8R = TAU*YW/((YQ*TAU + YWO) - YW)
endif

end function TAUF8R
