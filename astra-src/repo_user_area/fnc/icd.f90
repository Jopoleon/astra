! ICD [MA]: Driven current inside {0,R} 
double precision FUNCTION ICDR(YR)

use status, only: CD
use standard_functions, only: IINT

implicit none

double precision, intent(in) :: YR

ICDR = IINT(CD, YR)

end function ICDR
