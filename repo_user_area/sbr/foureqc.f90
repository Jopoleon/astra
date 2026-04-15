subroutine foureqc()

use const_inc, only: NA1, NEQUIl, MEQUIL, IPART, GP2, ABC
use status_inc, only: FP_NORM
use parameters_a2equil, only: equil_now
use numerical_tools, only: qinterp

implicit none

integer, parameter :: four_order_max=7, four_types=8, unit_out=223 ! #Four. coeff: MAX 11 for num stability. 7-9 Optimal. Including mode 0. 
integer :: jrho, jthe, mom_order, mom_type, nrho, ntheta
double precision :: Z0
double precision :: four_coef_as(four_types, four_order_max, NA1)
double precision, allocatable, dimension(:) :: theta, theta_half, dtheta, cos_mthe, sin_mthe, damin, psi_n
double precision, allocatable, dimension(:, :) :: RR, ZZ, drdr, dzdr
double precision, allocatable, dimension(:, :, :) :: four_coef
character(len=120) :: f_four

if (IPART == 1) return

! nequil, mequil are reals
nrho   = int(nequil)
ntheta = int(mequil)

allocate(theta(ntheta), dtheta(ntheta), cos_mthe(ntheta), sin_mthe(ntheta), damin(ntheta))
allocate(theta_half(ntheta+1))
allocate(psi_n(nrho))
allocate(RR(nrho, ntheta), ZZ(nrho, ntheta), drdr(nrho, ntheta), dzdr(nrho, ntheta))
allocate(four_coef(four_types, four_order_max, nrho))

! Set in xpr/tglf_interf.f90: n_fourier_modes = four_order_max; tglf_n_fourier_in = four_order_max - 1 
! (for file reading and tglf_four_coeff setting respectevely)

do jthe=1, ntheta ! R, Z, r, theta, psin values
    theta(jthe) = equil_now%coord_sys%position%theta2d(jthe) ! theta(ntheta)
    do jrho=1, nrho
        RR(jrho, jthe) = equil_now%coord_sys%position%r(jrho, jthe)    ! R(nrho, ntheta)
        ZZ(jrho, jthe) = equil_now%coord_sys%position%z(jrho, jthe)    ! Z(nrho, ntheta)
        psi_n(jrho)    = equil_now%profiles_1d%psi(jrho) ! psin(nrho)
    enddo
enddo

! Z of plasma center
Z0 = ZZ(1, 1)

theta(1:ntheta) = theta(1:ntheta) - theta(1) !Theta < 2PI

! dtheta calculation
do jthe=2, ntheta
    theta_half(jthe) = 0.5*(theta(jthe) + theta(jthe-1))
enddo
theta_half(1) = 0.5*(theta(1) + (theta(ntheta) - GP2))
theta_half(ntheta+1) = theta_half(1) + GP2

do jthe=1, ntheta
    dtheta(jthe) = abs(theta_half(jthe+1) - theta_half(jthe)) !d_teheta
enddo

! Radial derivatives
do jrho=2, nrho-1
    psi_n(jrho) = (psi_n(jrho) - psi_n(1))/(psi_n(nrho) - psi_n(1)) !psiN
    damin(:) = equil_now%coord_sys%position%rmin(jrho+1, 1:ntheta) - equil_now%coord_sys%position%rmin(jrho-1, 1:ntheta)
    drdr(jrho, 1:ntheta) = (RR(jrho+1, :) - RR(jrho-1, :))/damin(:)
    dzdr(jrho, 1:ntheta) = (ZZ(jrho+1, :) - ZZ(jrho-1, :))/damin(:)
enddo
drdr(1, :) = drdr(2, :)
dzdr(1, :) = dzdr(2, :)
drdr(nrho, :) = drdr(nrho-1, :)
dzdr(nrho, :) = dzdr(nrho-1, :)
psi_n(1) = 0.
psi_n(nrho) = 1.
ZZ = ZZ - Z0

! Calculate Fourier coefficients
do mom_order=1, four_order_max
    cos_mthe(:) = cos((mom_order-1.)*theta(:))*dtheta(:)
    sin_mthe(:) = sin((mom_order-1.)*theta(:))*dtheta(:)
    do jrho=1, nrho
        four_coef(1, mom_order, jrho) = sum(RR  (jrho, :)*cos_mthe)  ! <Rcos>
        four_coef(2, mom_order, jrho) = sum(RR  (jrho, :)*sin_mthe)  ! <Rsin>
        four_coef(3, mom_order, jrho) = sum(ZZ  (jrho, :)*cos_mthe)  ! <Zcos>
        four_coef(4, mom_order, jrho) = sum(ZZ  (jrho, :)*sin_mthe)  ! <Rsin>
        four_coef(5, mom_order, jrho) = sum(drdr(jrho, :)*cos_mthe)  ! <Zsin>
        four_coef(6, mom_order, jrho) = sum(drdr(jrho, :)*sin_mthe)  ! <Rsin>
        four_coef(7, mom_order, jrho) = sum(dzdr(jrho, :)*cos_mthe)  ! <Zsin>
        four_coef(8, mom_order, jrho) = sum(dzdr(jrho, :)*sin_mthe)  ! <Rsin>
    enddo
enddo

! ABC normalization as in ASTRA-TGLF
four_coef = four_coef*2./sum(dtheta)
four_coef(1:4, :, :) = four_coef(1:4, :, :)/ABC

! Interpolation to ASTRA grid
do mom_order=1, four_order_max
    do mom_type=1, four_types
        call qinterp(psi_n(:), four_coef(mom_type, mom_order, :), nrho, FP_NORM(1: NA1), four_coef_as(mom_type, mom_order, :), NA1)
    enddo
enddo

! Write file xpr/fort.four_coef for TGLF
f_four = 'xpr/fort.four_coef'
open(unit_out, file=TRIM(f_four))
do jrho=1, NA1
    do mom_order=1, four_order_max
        write(unit_out, '(8e14.6)') four_coef_as(:, mom_order, jrho)
    enddo
enddo
write(*, *) 'Written Fourier coefficients for TGLF into '//TRIM(f_four)
close(unit_out)

deallocate(theta, dtheta, cos_mthe, sin_mthe, damin, theta_half, psi_n, RR, ZZ, drdr, dzdr, four_coef)

return
end subroutine foureqc
