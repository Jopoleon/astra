! IBS [MA]: Toroidal bootstrap current inside {0,R} 
double precision FUNCTION IBSR(YR)

use status_inc, only: CUBS

implicit none

double precision, intent(in) :: YR
double precision, external :: IINT

IBSR = IINT(CUBS, YR)

return
end function IBSR
