!======================================================================|
! Gaussian deposition function P=P0*exp(-((r-r0*a)/a/width)^2)
!  r    [m] - coordinate in the equatorial cross-section 
!  r0   [d/l] - relative position of the centre of distribution 
! width [d/l] - relative width, normalized to the plasma radius 
!  a [m] - minor radius
! Normalization: INT(PdV)=1
!
! Usage example: FGAUSS(CF1, CF2, CAR1):; PE=...+CF3*CAR1;

subroutine FGAUSS_RHOPOL(rhop_center, rhop_width, gauss)

use parameter_inc, only: NRD
use const_inc, only: NA1, ROC, HRO
use status_inc, only: RHO, VR, FP

implicit none

double precision, intent(in) :: rhop_center, rhop_width
double precision, intent(out), dimension(NRD) :: gauss
double precision :: rho_pol, pow, YR, gauss_vol_int
integer :: j

gauss_vol_int = 0.
do j=1, NA1
   rho_pol = SQRT((FP(j) - FP(1))/(FP(NA1) - FP(1)))
   YR = (rho_pol - rhop_center)/rhop_width
   gauss(j)= exp(-YR*YR)
   if (j < NA1) gauss_vol_int = gauss_vol_int + gauss(j)*VR(j)
enddo
gauss_vol_int = gauss_vol_int*HRO
do j=1, NA1
   gauss(j)= gauss(j)/gauss_vol_int
enddo

return
end subroutine fgauss_rhopol
