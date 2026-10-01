! betpol []: Internal inductance li(r) [4]
!            Pereverzev 03-APR-08
!
! betpol = Wp/Wm
!
! Wp = int(P), Wm = int(Btheta2/(2*mu0))
!----------------------------------------------------------------------|
double precision function BETP3R(YR)

use pi_const, only: GP2
use status, only: NE, TE, NI, TI, VR, MU, RHO, PBLON, PBPER, PFAST
use scalars, only: NA1, ROC, HRO, BTOR, RTOR

implicit none

double precision, intent(in) :: yr
integer :: JK, J
double precision :: Q, V, YRO, YK, YIPL, YWBP

JK = min(NA1, nint(ROC/HRO))
Q = 0.
YWBP = 0.
V = 0.
do J=1, JK
    V = V + VR(J)
    Q = Q + 1602.*(NE(J)*TE(J) + NI(J)*TI(J) + &
       pfast(j) + 0.5*(pblon(j) + pbper(j)))*VR(J)
    YWBP = YWBP + (BTOR*RHO(J)/RTOR*MU(J))**2/(4*GP2*1.e-7)*VR(J)
enddo

BETP3R = Q/YWBP

end function BETP3R
