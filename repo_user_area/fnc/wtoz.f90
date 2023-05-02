! WTOZ [MJ]:	Integral {0:R} ( 3/2*(NE*(TE+Ti)+pfast+0.5*(pblon+pbper) )) dV
!        (Fable Nov 2018) total energy including fast ions
double precision FUNCTION WTOZR(YR)

use const_inc
use status_inc

implicit none

double precision, intent(in) :: yr
integer :: jk,j

jk = nint(yr/hro)

WTOZR = 0.
do J=1,JK
    WTOZR = WTOZR + (NE(J)*TE(J) + NI(J)*TI(J) + &
        0.5*(PBLON(J) + PBPER(J)) + PFAST(J))*VR(J)
enddo
WTOZR = HRO*(WTOZR - (NE(JK)*TE(JK) + NI(JK)*TI(JK) + &
     & 0.5*(PBLON(J) + PBPER(J)) + PFAST(J))*HRO)*.0024

return
end function WTOZR
