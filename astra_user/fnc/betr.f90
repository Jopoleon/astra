! BETR  Troyon factor: beta_toroidal[%]/(I[MA]/a[m]B[T])
!  
! Example: beTr_BETRB;
!     (Pereverzev 02-MAY-2006)
double precision function BETRR(YR)

use const_inc, only: BTOR

implicit  none

double precision YR, AFR, BETAR, ITOTR
external AFR, BETAR, ITOTR

BETRR = BETAR(YR)*AFR(YR)*BTOR/ITOTR(YR)

return
end function BETRR
