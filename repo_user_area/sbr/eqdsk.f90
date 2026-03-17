subroutine EQDSK(coco_number)

use pi_const, only: GP2
use scalars, only: RTOR, BTOR, IPL, TIME, TSTART, SGNBT, SGNIP, NA1
use parameters_a2equil, only: equil_now
use read_input, only: awd, exp_file, equ_file
use status, only: MU, FP_NORM
use numerical_tools, only: qinterp
use surface_contours, only: ctr2rz_fun3

implicit none

integer, parameter :: nrRect=257, nzRect=257, eqdsk_unit=11
integer, dimension(8), parameter :: coco_dpsi_sign=(/1, 1, -1, -1, 1, 1, -1, -1/)

integer, intent(in) :: coco_number

integer :: i, j, nrho_surf, nthe_surf
double precision :: dpsin_rect, Rmin, Rmax, zmin, zmax, dr, dz, dpsi_sgn, psi_2pi
double precision, allocatable, dimension(:) :: psin_eq
double precision, dimension(nrRect) :: r_rect, psin_rect, pres_rect, &
     fdia_rect, q_rect, pprime_rect, fprime_rect
double precision, dimension(nzRect) :: z_rect
double precision, dimension(nrRect, nzRect) :: psi_rect
character(len=120) :: f_eqdsk

if (TIME <= TSTART) return

dpsi_sgn = coco_dpsi_sign(coco_number)
if (coco_number < 10) then
    psi_2pi = 1./GP2
else
    psi_2pi = 1.
endif
nrho_surf = SIZE(equil_now%coord_sys%position%r, 1)
nthe_surf = SIZE(equil_now%coord_sys%position%r, 2)
allocate(psin_eq(nrho_surf))

dpsin_rect = 1./(nrRect - 1.d0)
psin_rect = (/ (dpsin_rect*(i - 1.d0), i=1, nrRect) /)

! Getting 1d profiles from equilibrium
psin_eq = (equil_now%profiles_1d%psi - equil_now%profiles_1d%psi(1))/(equil_now%profiles_1d%psi(nrho_surf) - equil_now%profiles_1d%psi(1))

Rmin = MINVAL(equil_now%coord_sys%position%r(nrho_surf, :)) - 0.03
Rmax = MAXVAL(equil_now%coord_sys%position%r(nrho_surf, :)) + 0.03
zmin = MINVAL(equil_now%coord_sys%position%z(nrho_surf, :)) - 0.03
zmax = MAXVAL(equil_now%coord_sys%position%z(nrho_surf, :)) + 0.03
dr = (Rmax - Rmin)/(nrRect - 1.d0)
dz = (zmax - zmin)/(nzRect - 1.d0)
r_rect = (/ (Rmin + dr*(i - 1.d0), i=1, nrRect) /)
z_rect = (/ (zmin + dz*(i - 1.d0), i=1, nzRect) /)

! Biquadratic interpolation from psi(rho, theta) to psi(R, Z)
call ctr2rz_fun3(nrho_surf, nthe_surf, equil_now%profiles_1d%psi, &
    equil_now%coord_sys%position%r, equil_now%coord_sys%position%z, nrRect, nzRect, r_rect, z_rect, psi_rect)

! Quadratic interpolation ro nrRect grid (eqdsk requires that)
call qinterp(psin_eq, equil_now%profiles_1d%pressure, nrho_surf, psin_rect,   pres_rect, nrRect)
call qinterp(psin_eq, equil_now%profiles_1d%F_dia   , nrho_surf, psin_rect,   fdia_rect, nrRect)
call qinterp(psin_eq, equil_now%profiles_1d%pprime  , nrho_surf, psin_rect, pprime_rect, nrRect)
call qinterp(psin_eq, equil_now%profiles_1d%ffprime , nrho_surf, psin_rect, fprime_rect, nrRect)
call qinterp(FP_NORM(1:na1), 1./MU(1:NA1), NA1, psin_rect, q_rect, nrRect)

! EQDSk file output
if (TIME < 10.) then
    write(f_eqdsk, '(5A, f5.3, A)') TRIM(awd), '/ncdf_out/', TRIM(exp_file), TRIM(equ_file), '0', TIME, '.eqdsk'
else
    write(f_eqdsk, '(4A, f6.3, A)') TRIM(awd), '/ncdf_out/', TRIM(exp_file), TRIM(equ_file), TIME, '.eqdsk'
endif
 
write(*, '(2A, i4)') 'Storing ' // TRIM(f_eqdsk), '   nR =', nrRect

open(eqdsk_unit, file=TRIM(f_eqdsk))
write(eqdsk_unit, '(A48, 3i4)') 'ASTRA', 3, nrRect, nzRect
! Boxdim(R, m), BOxdim(Z, m), R0(vacuum), Rmin(box, m), Zmid(box, m)
write(eqdsk_unit, '(5E16.9)') R_rect(nrRect) - R_rect(1), Z_rect(nzRect) - Z_rect(1), RTOR, &
    R_rect(1), 0.5*(Z_rect(1) + Z_rect(nzRect))
! Rmagnaxis(m), Zmagnaxis(m)
write(eqdsk_unit, '(5E16.9)') equil_now%coord_sys%position%r(1, 1), equil_now%coord_sys%position%z(1, 1), &
    dpsi_sgn*SGNIP*equil_now%profiles_1d%psi(1)*psi_2pi, SGNIP*equil_now%profiles_1d%psi(nrho_surf)*psi_2pi, SGNBT*BTOR
write(eqdsk_unit, '(5E16.9)') SGNIP*IPL*1.d6, dpsi_sgn*SGNIP*equil_now%profiles_1d%psi(1)*psi_2pi, 0., equil_now%coord_sys%position%r(1, 1), 0.
write(eqdsk_unit, '(5E16.9)') equil_now%coord_sys%position%z(1, 1), 0., dpsi_sgn*SGNIP*equil_now%profiles_1d%psi(nrho_surf)*psi_2pi, 0., 0.
write(eqdsk_unit, '(5E16.9)') (SGNBT*fdia_rect(i), i=1, nrRect)
write(eqdsk_unit, '(5E16.9)') (pres_rect(i), i=1, nrRect)
write(eqdsk_unit, '(5E16.9)') (dpsi_sgn*SGNIP*fprime_rect(i)/psi_2pi, i=1, nrRect)
write(eqdsk_unit, '(5E16.9)') (dpsi_sgn*SGNIP*pprime_rect(i)/psi_2pi, i=1, nrRect)
write(eqdsk_unit, '(5E16.9)') ((dpsi_sgn*SGNIP*psi_rect(i, j)*psi_2pi, i=1, nrRect), j=1, nzRect)
write(eqdsk_unit, '(5E16.9)') (SGNBT*SGNIP*q_rect(i), i=1, nrRect)
write(eqdsk_unit, '(2i5)') nthe_surf, nthe_surf
write(eqdsk_unit, '(5E16.9)') (equil_now%coord_sys%position%r(nrho_surf, i), equil_now%coord_sys%position%z(nrho_surf, i), i=1, nthe_surf)
write(eqdsk_unit, '(5E16.9)') (equil_now%coord_sys%position%r(nrho_surf, i), equil_now%coord_sys%position%z(nrho_surf, i), i=1, nthe_surf)
close(eqdsk_unit)

deallocate(psin_eq)

end subroutine EQDSK
