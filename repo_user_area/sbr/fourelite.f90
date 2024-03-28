subroutine fourelite

use const_inc, only: NA1, NEQUIl, MEQUIL, IPART, GP2, ABC, PSIAX
use status_inc, only: FP
use parameters_a2equil, only: equil_now
                                         
implicit none

integer, parameter :: four_order_max=7, four_types=8, unit_out=223 ! #Four. coeff: MAX 11 for num stability. 7-9 Optimal. Including mode 0. 
integer :: jrho, jthe, mom_order, mom_type, nrho, ntheta, nELITE
double precision :: Z0, psi_as(NA1)
double precision :: four_coef_as(four_types, four_order_max, NA1)
double precision, allocatable, dimension(:) :: theta, theta_half, dtheta, cos_mthe, sin_mthe, damin, psi_n, theta_el
double precision, allocatable, dimension(:, :) :: RR, ZZ, drdr, dzdr !Fourier
double precision, allocatable, dimension(:, :, :) :: four_coef
double precision, allocatable, dimension(:, :) :: BBpel,RRel,ZZel,BBpel_as,RRel_as,ZZel_as,BBpel_tg,RRel_tg,ZZel_tg!Elite
character(len=120) :: f_four


if (IPART == 1) return

! nequil, mequil are reals
nrho   = int(nequil)
ntheta = int(mequil)
nELITE = 300

allocate(dtheta(ntheta), cos_mthe(ntheta), sin_mthe(ntheta), damin(ntheta))
allocate(theta(ntheta+1), theta_half(ntheta+1), theta_el(nELITE))
allocate(psi_n(nrho))
allocate(RR(nrho, ntheta), ZZ(nrho, ntheta), drdr(nrho, ntheta), dzdr(nrho, ntheta))
allocate(RRel(nrho,ntheta+1), ZZel(nrho,ntheta+1), BBpel(nrho,ntheta+1)) !ELITE FEQUIS
allocate(RRel_as(NA1,ntheta+1), ZZel_as(NA1,ntheta+1), BBpel_as(NA1,ntheta+1)) !ELITE ASTRA-TGLF
allocate(RRel_tg(NA1,nELITE+1), ZZel_tg(NA1,nELITE+1), BBpel_tg(NA1,nELITE+1)) !ELITE ASTRA-TGLF
allocate(four_coef(four_types, four_order_max, nrho))

! Set in xpr/tglf_interf.f90: n_fourier_modes = four_order_max; tglf_n_fourier_in = four_order_max - 1 
! (for file reading and tglf_four_coeff setting respectevely)

do jthe=1, ntheta ! R, Z, r, theta, psin values
    theta(jthe) = equil_now%coord_sys%position%teta2d(jthe) ! theta(ntheta)
    do jrho=1, nrho
        RR(jrho, jthe) = equil_now%coord_sys%position%r(jrho, jthe)    ! R(nrho, ntheta)
        ZZ(jrho, jthe) = equil_now%coord_sys%position%z(jrho, jthe)    ! Z(nrho, ntheta)
        psi_n(jrho)    = equil_now%profiles_1d%psi(jrho) ! psin(nrho)
        BBpel(jrho, jthe) = equil_now%coord_sys%bpcell(jrho, jthe) !BPOL
    enddo
enddo

BBpel(nrho,1:ntheta) = BBpel(nrho-1,1:ntheta) !Defined up to nequil-1

! Z of plasma center
Z0 = ZZ(1, 1)

theta(ntheta+1) = theta(1)+GP2 !close loop on theta
theta(1:ntheta+1) = theta(1:ntheta+1) - theta(1) !Theta < 2PI

!ELITE geometry
RRel(:,1:ntheta) = RR/ABC !normalize for TGLF
ZZel(:,1:ntheta) = ZZ/ABC
RRel(:,ntheta+1) = RRel(:,1) !close loop R
ZZel(:,ntheta+1) = ZZel(:,1) !close loop Z
BBpel(:,ntheta+1)= BBpel(:,1) !close loop Bpol, TO BE NORMALIZED IN TGLF_INTERF

psi_as(1:NA1) = (FP(1:NA1)-PSIAX)/(FP(NA1)-PSIAX) !ASTRA PSIN

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
    cos_mthe(:) = cos((mom_order-1.)*theta(1:ntheta))*dtheta(:)
    sin_mthe(:) = sin((mom_order-1.)*theta(1:ntheta))*dtheta(:)
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

! Interpolation FOURIER to ASTRA grid
do mom_order=1, four_order_max
    do mom_type=1, four_types
        !call qinterp(psi_n(:), four_coef(mom_type, mom_order, :), nrho, FP_NORM(1: NA1), four_coef_as(mom_type, mom_order, :), NA1)
        call qinterp(psi_n(:), four_coef(mom_type, mom_order, :), nrho, psi_as(1: NA1), four_coef_as(mom_type, mom_order, :), NA1)
    enddo
enddo

! Interpolation ELITE to ASTRA-TGLF grids
!PSI interpolation for ASTRA
do jthe=0,nELITE+1
  theta_el(jthe+1) = jthe*GP2/(nELITE+0.d0)
enddo
!Interpolation on ASTRA grid
!do jthe=1,ntheta+1
!call qinterp(psi_n(:), RRel(:,jthe), nrho, FP_NORM(1:NA1), RRel_as(:,jthe), NA1)
!call qinterp(psi_n(:), ZZel(:,jthe), nrho, FP_NORM(1:NA1), ZZel_as(:,jthe), NA1)
!call qinterp(psi_n(:), BBpel(:,jthe), nrho, FP_NORM(1:NA1), BBpel_as(:,jthe), NA1)
!enddo

!Interpolation on ASTRA grid
do jthe=1,ntheta+1
call qinterp(psi_n(:), RRel(:,jthe), nrho, psi_as(1:NA1), RRel_as(:,jthe), NA1)
call qinterp(psi_n(:), ZZel(:,jthe), nrho, psi_as(1:NA1), ZZel_as(:,jthe), NA1)
call qinterp(psi_n(:), BBpel(:,jthe), nrho, psi_as(1:NA1), BBpel_as(:,jthe), NA1)
enddo


!Interpolation on equispaced theta
do jrho=1,NA1
call qinterp(theta(:), RRel_as(jrho,:), ntheta+1, theta_el, RRel_tg(jrho,:), nELITE+1)
call qinterp(theta(:), ZZel_as(jrho,:), ntheta+1, theta_el, ZZel_tg(jrho,:), nELITE+1)
call qinterp(theta(:), BBpel_as(jrho,:), ntheta+1, theta_el, BBpel_tg(jrho,:), nELITE+1)
enddo


! Write file xpr/fort.four_coef for TGLF
f_four = 'xpr/fort.four_coef'
open(unit_out, file=TRIM(f_four))
do jrho=1, NA1
    do mom_order=1, four_order_max
        write(unit_out, '(8F)') four_coef_as(:, mom_order, jrho)
    enddo
enddo
!write(*, *) 'Written Fourier coefficients for TGLF into '//TRIM(f_four)
close(unit_out)

f_four = 'xpr/fort.elite_geom'
open(unit_out, file=TRIM(f_four))
do jrho=1, NA1
    do jthe=1,nELITE+1
        write(unit_out, '(3F)') RRel_tg(jrho,jthe), ZZel_tg(jrho,jthe), BBpel_tg(jrho,jthe)
    enddo
enddo
!write(*, *) 'Written Elite geometry for TGLF into '//TRIM(f_four)
close(unit_out)

f_four = 'xpr/fort.elite_geom_flip'
open(unit_out, file=TRIM(f_four))
do jrho=1, NA1
    do jthe=nELITE+1,1,-1
        write(unit_out, '(3F)') RRel_tg(jrho,jthe), ZZel_tg(jrho,jthe), BBpel_tg(jrho,jthe)
    enddo
enddo
!write(*, *) 'Written Elite geometry for TGLF into '//TRIM(f_four)
close(unit_out)


f_four = 'xpr/fort.elite_fequis'
open(unit_out, file=TRIM(f_four))
    do jthe=1,ntheta+1
        write(unit_out, '(4F)') theta(jthe), RRel(20,jthe), ZZel(20,jthe), BBpel(nrho,jthe)
    enddo
!write(*, *) 'Written Elite geometry for TGLF into '//TRIM(f_four)
close(unit_out)

f_four = 'xpr/fort.elite_tglf'
open(unit_out, file=TRIM(f_four))
    do jthe=1,nELITE+1
        write(unit_out, '(4F)') theta_el(jthe), RRel_tg(36,jthe), ZZel_tg(36,jthe), BBpel_tg(NA1,jthe)
    enddo
!write(*, *) 'Written Elite geometry for TGLF into '//TRIM(f_four)
close(unit_out)


deallocate(theta, dtheta, cos_mthe, sin_mthe, damin, theta_half, psi_n, RR, ZZ, drdr, dzdr, four_coef)
deallocate(BBpel, RRel, ZZel, BBpel_as, RRel_as, ZZel_as, BBpel_tg, RRel_tg, ZZel_tg, theta_el) 

return
end subroutine fourelite
