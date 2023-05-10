! IBS [MA]: Total current inside {0,R}
double precision FUNCTION ITOTR(YR)

use status_inc, only: CU

implicit none

double precision, intent(in) :: YR
double precision, external :: IINT

ITOTR = IINT(CU, YR)

return
end function ITOTR
