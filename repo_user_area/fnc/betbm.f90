! BETBMR [ ]: Beta beam poloidal (r)
!
!   2*c*c      2*c*c
!  - - - - - - *{Int(p*dS) - p} = -  - - - - - - * Int[S*dp]
!    J*J       J*J
!
double precision function BETBMR(YR)

use const_inc, only: NA, BTOR, RTOR, HRO, GP2
use status_inc, only: ELON, AMETR, IPOL, PBPER, MU, G22

implicit none

double precision YR, YB
integer  J, JK

YB  = 0.
JK  = YR/HRO + 0.5
IF(JK < 1) JK = 2
IF(JK > NA) JK = NA
DO J = 2, JK
   YB  = YB + ELON(J)*AMETR(J)**2*(PBPER(J - 1) - PBPER(J + 1))
enddo
BETBMR  = 6.4E-4*GP2*YB*0.5*(RTOR/(G22(JK)*IPOL(JK)*BTOR*JK*HRO*MU(JK)))**2

return
end function BETBMR
