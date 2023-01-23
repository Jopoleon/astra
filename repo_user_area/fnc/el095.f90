double precision function el095r(YD)
!---------------------------------------- Polevoi,  10-SEP-99 --------

use const_inc, only: NA1, BTOR
use status_inc, only: MU, ELON, FP, RHO

double precision FPC0, YM, YYQ, YM1, YQ, YDF
integer JJ1

! el95d=elon(r=r95) & FP(r95)=.95*FP(a)
YQ = 1.d0
YM1 = ELON(NA1)
JJ1 = NA1
FPC0 = FP(1) - BTOR*GP*MU(1)*RHO(1)**2

YDF = 1./(FPC0 - FP(NA1))
do
   JJ1 = JJ1 - 1
   YYQ = YQ
   YM = YM1
   YM1 = ELON(JJ1)
   YQ = 1.d0 - (FP(JJ1) - FP(NA1))*YDF
   if(YQ <= 0.95d0 .or. JJ1 <= 1) EXIT
enddo
el095r = YM1 + (.95d0 - YQ)*(YM - YM1)/(YYQ - YQ)

return
end function el095r
