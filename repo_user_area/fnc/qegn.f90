! QEGN [MW]:  Integral {0, R} (PEGN) dV [=5/2*V'(R)*Gn(R)*Te(R)/625]
!   (Yushmanov 20-FEB-89)
double precision FUNCTION QEGNR(YR)

use const_inc, only: NA, NA1, HRO
use status_inc, only: G11, GN, TE

implicit none

double precision, intent(in) :: YR
integer :: JK
double precision :: YDR

JK = YR/HRO
IF(JK >= NA) THEN
   QEGNR = .001*(G11(NA) + G11(NA1))*GN(NA)*(TE(NA+1) + TE(NA))
   return
ELSEIF(JK < 0) THEN
   QEGNR = 0.
   return
ELSEIF(JK == 0) THEN
   QEGNR = YR*(G11(1) + G11(2))*GN(1)*(TE(2) + TE(1))/(HRO*1000.)
   return
ENDIF
YDR = YR/HRO - JK
QEGNR = ((1. - YDR)*(G11(JK) + G11(JK+1))*GN(JK)*(TE(JK+1) + TE(JK)) + &
       YDR*(G11(JK+2) + G11(JK+1))*GN(JK+1)*(TE(JK+2) + TE(JK+1)))*.001

end function QEGNR
