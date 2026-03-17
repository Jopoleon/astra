! IBS [MA]: Total current inside {0,R}
double precision FUNCTION ITOTR(YR)

use status, only: CU
use standard_functions, only: IINT

implicit none

double precision, intent(in) :: YR

ITOTR = IINT(CU, YR)

end function ITOTR
