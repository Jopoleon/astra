! IBS [MA]: Total current inside {0,R}
double precision FUNCTION ITOTR(YR)

use status_inc, only: CU
use standard_functions, only: IINT

implicit none

double precision, intent(in) :: YR

ITOTR = IINT(CU, YR)

return
end function ITOTR
