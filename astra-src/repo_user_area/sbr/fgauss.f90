!======================================================================|
! Gaussian deposition function P=P0*exp(-((r-r0*a)/a/width)^2)
!  r    [m] - coordinate in the equatorial cross-section 
!  r0   [d/l] - relative position of the centre of distribution 
! width [d/l] - relative width, normalized to the plasma radius 
!  a [m] - minor radius
! Normalization: INT(PdV)=1
!
! Usage example: FGAUSS(CF1, CF2, CAR1):; PE=...+CF3*CAR1;

subroutine FGAUSS(YCENTR, YWIDTH, YPROF)

use scalars, only: NA1, ROC, HRO
use status, only: RHO, VR

implicit none

double precision :: YCENTR, YWIDTH, YPROF(*), YR, YRO, YWH, YPOW

integer j

YRO = YCENTR*ROC
YWH = YWIDTH*ROC
YPOW = 0.
do j=1, NA1
   YR = (RHO(j) - YRO)/YWH 
   YPROF(j)= exp(-YR*YR)
   if (j < NA1) YPOW = YPOW + YPROF(j)*VR(j)
enddo
YPOW = YPOW*HRO
do j=1, NA1
   YPROF(j)= YPROF(j)/YPOW
enddo

end subroutine fgauss
