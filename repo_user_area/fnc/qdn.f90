! QDN [10^19/s]:  Integral {0, R} ( SDN ) dV
!   (Pereverzev 9-OCT-07)
double precision FUNCTION QDNR(YR)

use status_inc, only: SDN

implicit none

double precision, intent(in) :: YR
double precision :: VINT

QDNR = VINT(SDN, YR)

end function QDNR
