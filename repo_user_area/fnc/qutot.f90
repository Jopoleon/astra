! QUTOT 

double precision FUNCTION QUTOTR(YR)

use status_inc, only: TTRQ

implicit none

double precision, intent(in) :: YR
double precision :: VINT

QUTOTR = VINT(TTRQ, YR)

end function QUTOTR
