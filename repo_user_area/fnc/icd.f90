! ICD [MA]: Driven current inside {0,R} 
double precision FUNCTION ICDR(YR)

use status_inc, only: CD

implicit none

double precision, intent(in) :: YR
double precision, external :: IINT

ICDR = IINT(CD, YR)

return
end function ICDR
