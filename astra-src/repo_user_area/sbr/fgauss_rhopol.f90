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

use scalars, only: NA1, ROC, HRO
use status, only: NRD, RHO, VR, rho_pol

implicit none

double precision, intent(in) :: rhop_center, rhop_width
double precision, intent(out), dimension(NRD) :: gauss
double precision :: pow, YR, gauss_vol_int
integer :: j

gauss_vol_int = 0.
do j=1, NA1
   YR = (rho_pol(j) - rhop_center)/rhop_width
   gauss(j)= exp(-YR**2)
   if (j < NA1) gauss_vol_int = gauss_vol_int + gauss(j)*VR(j)
enddo
gauss_vol_int = gauss_vol_int*HRO
do j=1, NA1
   gauss(j)= gauss(j)/gauss_vol_int
enddo

end subroutine fgauss_rhopol
