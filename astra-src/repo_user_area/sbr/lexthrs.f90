subroutine lexthrs(lext)

! external inductance LEXT from Hirshman, Neilson, Phys fluids 20 3 1986
! lext in microHenri

use scalars, only: RTOR, ABC, ELONG
use numerical_tools, only: qinterp

implicit none

double precision, intent(out) :: lext
double precision, dimension(9) :: eps, a, b
double precision, dimension(1) :: z1, a0, b0

z1(1) = ABC/RTOR

eps = (/0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9/)
a = (/2.975, 2.285, 1.848, 1.507, 1.217, 0.957, 0.716, 0.487, 0.272/)
b = (/0.228, 0.325, 0.403, 0.465, 0.512, 0.542, 0.553, 0.538, 0.508/)

call qinterp(eps, a, 9, z1, a0, 1)
call qinterp(eps, b, 9, z1, b0, 1)

a0(1) = (1. + 1.81*sqrt(z1(1)) + 2.05*z1(1))*log(8/z1(1)) - &
 (2. + 9.25*sqrt(z1(1)) - 1.21*z1(1))
b0(1) = 0.73*sqrt(z1(1))*(1. + 2.*z1(1)**4 - 6*z1(1)**5 + 3.7*z1(1)**6)

lext = 0.4*3.141592*rtor*a0(1)*(1. - z1(1))/(1. - z1(1) + b0(1)*ELONG)

end subroutine lexthrs
