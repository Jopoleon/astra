! IBS [MA]: Toroidal bootstrap current inside {0,R} 
double precision FUNCTION IBSR(YR)

use status, only: CUBS
use standard_functions, only: IINT

implicit none

double precision, intent(in) :: YR

IBSR = IINT(CUBS, YR)

end function IBSR
