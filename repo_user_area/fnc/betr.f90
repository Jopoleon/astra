! BETR  Troyon factor: beta_toroidal[%]/(I[MA]/a[m]B[T])
!  
! Example: beTr_BETRB;
!     (Pereverzev 02-MAY-2006)
double precision function BETRR(YR)

use scalars, only: BTOR
use standard_functions, only: AFR

implicit  none

double precision :: YR
double precision, external :: BETAR, ITOTR

BETRR = BETAR(YR)*AFR(YR)*BTOR/ITOTR(YR)

end function BETRR
