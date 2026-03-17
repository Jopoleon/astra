! QDN [10^19/s]:  Integral {0, R} ( SDN ) dV
!   (Pereverzev 9-OCT-07)
double precision FUNCTION QDNR(YR)

use status, only: SDN
use standard_functions, only: VINT

implicit none

double precision, intent(in) :: YR

QDNR = VINT(SDN, YR)

end function QDNR
