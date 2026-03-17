! QUTOT 

double precision FUNCTION QUTOTR(YR)

use status, only: TTRQ
use standard_functions, only: VINT

implicit none

double precision, intent(in) :: YR

QUTOTR = VINT(TTRQ, YR)

end function QUTOTR
