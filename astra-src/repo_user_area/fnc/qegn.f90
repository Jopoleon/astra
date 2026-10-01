! QEGN [MW]:  Integral {0, R} (PEGN) dV [=5/2*V'(R)*Gn(R)*Te(R)/625]
!   (Yushmanov 20-FEB-89)
double precision FUNCTION QEGNR(YR)

use scalars, only: NA, NA1, HRO
use status, only: G11, GN, TE

implicit none

double precision, intent(in) :: YR
integer :: JK
double precision :: YDR

JK = YR/HRO
if (JK >= NA) then 
   QEGNR = .001*(G11(NA) + G11(NA1))*GN(NA)*(TE(NA+1) + TE(NA))
else if (JK < 0) then
   QEGNR = 0.
else if (JK == 0) then 
   QEGNR = YR*(G11(1) + G11(2))*GN(1)*(TE(2) + TE(1))/(HRO*1000.)
else
    YDR = YR/HRO - JK
    QEGNR = ((1. - YDR)*(G11(JK) + G11(JK+1))*GN(JK)*(TE(JK+1) + TE(JK)) + &
       YDR*(G11(JK+2) + G11(JK+1))*GN(JK+1)*(TE(JK+2) + TE(JK+1)))*.001
endif

end function QEGNR
