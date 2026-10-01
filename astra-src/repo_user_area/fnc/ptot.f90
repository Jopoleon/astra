! PTOT [MW/m#2]: Total plasma heating power (r) [m]
!   (Yushmanov 12-MAY-87)
double precision FUNCTION PTOTR(YR)

use scalars, only: HRO, NA1
use status, only: PETOT, PITOT

implicit none

double precision, intent(in) :: YR
integer :: JK
double precision :: YDH

JK = int(YR/HRO + .5)
if (JK >= NA1) then
   PTOTR = PETOT(NA1) + PITOT(NA1)
else if(JK <= 0) then
   PTOTR = PETOT(1) + PITOT(1)
else
   YDH = YR/HRO + (.5 - JK)
   PTOTR = YDH*(PETOT(JK+1) + PITOT(JK+1)) + (1. - YDH)*(PETOT(JK) + PITOT(JK))
endif

end function PTOTR
