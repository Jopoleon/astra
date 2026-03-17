! ALIM [ ]: Alpha_MHD at the ballooning limit
!       The function ALIM is created on the basis of G.W.Pacher
!       subroutine ZALF.F
!
! Input: ELON, TRIA, SHEAR, MU, SQEPS
double precision function ALIMR(YR)

use status, only: ELON, SQEPS, SHEAR, TRIA, MU

implicit none

double precision :: YEL, YE2, YEPS, YF2, YEK, YALF, YSHR, YF1, YR
integer :: j, node

j = node(YR)
YEL = ELON(J)
YEPS = SQEPS(J)**2  ! = AMETR(j)/(RTOR+SHIF(j))
! shear for alpha-min=0.5, so limit to this value
! gmsh multiplies shear on S-alfa diagram (is set to 1.)
YSHR = max(0.05, SHEAR(j))
YE2 = YEL**2
YF1 = 4.*YEL**1.5*(3. + YE2)/((1. + YE2)*(1. + 3.*YE2))
YF2 = ((1. + YE2)/(1. + YEL))**3/(YE2*(1. + 3.*YE2))
YEK = (YE2 - 1.)/(YE2 + 1.)
YALF =   (1. - 1.5*YEK*(1. + YEK)/(2. + YEK) + &
   6.*TRIA(J)/YEPS*YEK*(1. - YEK)/(2. + YEK))
YF1 = -0.25*(2.*YSHR*YF1 + YEK**2*YF2)
YF2 = YEPS*YALF - YEPS*MU(J)**2 - 1.5*exp(-1./abs(YSHR))/sqrt(YEL)
YALF = -YF2/(2.*YF1) + sqrt((YF2/(2.*YF1))**2 - YSHR**2/2./YF1)
! include correction for triangularity
YEPS = ((1.17988 + 5.03449*TRIA(J)**2)/1.17988)**2
! normalize to 1.00 at delta=0.155
YEPS  = YEPS/1.216
ALIMR = YEPS*YALF*2. ! Factor 2 due to GWP

end function alimr
