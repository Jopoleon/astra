subroutine EQDSK(coco_number)

use const_inc, only: RTOR, BTOR, IPL, TIME, TSTART, SGNBT, SGNIP
use parameters_a2equil, only: equil_now, GP2
use outcmn_inc, only: awd, exp_file, equ_file

implicit none

integer, parameter :: Nrrect=257, Nzrect=257, eqdsk_unit=11
integer, dimension(8), parameter :: coco_psi_sign=(/1, 1, -1, -1, 1, 1, -1, -1/)

integer, intent(in) :: coco_number

integer :: i, j, nrho_surf, nthe_surf, n_Rrect, n_Zrect
double precision :: dpsin_rect, Rmin, Rmax, zmin, zmax, dr, dz, psi_sgn
double precision, allocatable, dimension(:) :: r_rect, z_rect, &
    psin_eq, psin_rect, pres_rect, fdia_rect, q_rect, pprime_rect, fprime_rect
double precision, allocatable, dimension(:, :) :: psi_rect
character(len=120) :: f_eqdsk

if (TIME <= TSTART) return

psi_sgn = coco_psi_sign(coco_number)
n_Rrect = Nrrect
n_Zrect = Nzrect

nrho_surf = SIZE(equil_now%coord_sys%position%r, 1)
nthe_surf = SIZE(equil_now%coord_sys%position%r, 2)
allocate(psin_eq(nrho_surf))
allocate(psin_rect(n_Rrect), pres_rect(n_Rrect), fdia_rect(n_Rrect), q_rect(n_Rrect), &
    pprime_rect(n_Rrect), fprime_rect(n_Rrect))
allocate(r_rect(n_Rrect), z_rect(n_Zrect), psi_rect(n_Rrect, n_Zrect))

dpsin_rect = 1./(n_Rrect - 1.d0)
psin_rect = (/ (dpsin_rect*(i - 1.d0), i=1, n_Rrect) /)

! Getting 1d profiles from equilibrium
psin_eq = (equil_now%profiles_1d%psi - equil_now%profiles_1d%psi(1))/(equil_now%profiles_1d%psi(nrho_surf) - equil_now%profiles_1d%psi(1))

Rmin = MINVAL(equil_now%coord_sys%position%r(nrho_surf, :)) - 0.03 ! 1.08
Rmax = MAXVAL(equil_now%coord_sys%position%r(nrho_surf, :)) + 0.03 ! 2.26
zmin = MINVAL(equil_now%coord_sys%position%z(nrho_surf, :)) - 0.03 ! -1.0
zmax = MAXVAL(equil_now%coord_sys%position%z(nrho_surf, :)) + 0.03 ! 1.0
dr = (Rmax - Rmin)/(n_Rrect - 1.d0)
dz = (zmax - zmin)/(n_Zrect - 1.d0)
r_rect = (/ (Rmin + dr*(i - 1.d0), i=1, n_Rrect) /)
z_rect = (/ (zmin + dz*(i - 1.d0), i=1, n_Zrect) /)

! Biquadratic interpolation from psi(rho,theta) to psi(R, Z)
call ctr2rz_fun3(nrho_surf, nthe_surf, equil_now%profiles_1d%psi, &
    equil_now%coord_sys%position%r, equil_now%coord_sys%position%z, n_Rrect, n_Zrect, r_rect, z_rect, psi_rect)

! Quadratic interpolation ro n_Rrect grid (eqdsk requires that)
call qinterp(psin_eq, equil_now%profiles_1d%pressure, nrho_surf, psin_rect,   pres_rect, n_Rrect)
call qinterp(psin_eq, equil_now%profiles_1d%F_dia   , nrho_surf, psin_rect,   fdia_rect, n_Rrect)
call qinterp(psin_eq, equil_now%profiles_1d%pprime  , nrho_surf, psin_rect, pprime_rect, n_Rrect)
call qinterp(psin_eq, equil_now%profiles_1d%ffprime , nrho_surf, psin_rect, fprime_rect, n_Rrect)
call qinterp(psin_eq, equil_now%profiles_1d%q       , nrho_surf, psin_rect,      q_rect, n_Rrect)

! EQDSk file output
if (TIME < 10.) then
    write(f_eqdsk, '(5A, f5.3, A)') TRIM(awd), '/ncdf_out/', TRIM(exp_file), TRIM(equ_file), '0', TIME, '.eqdsk'
else
    write(f_eqdsk, '(4A, f6.3, A)') TRIM(awd), '/ncdf_out/', TRIM(exp_file), TRIM(equ_file), TIME, '.eqdsk'
endif
 
write(*, '(2A, i4)') 'Storing ' // TRIM(f_eqdsk), '   nR =', n_Rrect

open(eqdsk_unit, file=TRIM(f_eqdsk))
write(eqdsk_unit, '(A48, 3i4)') 'ASTRA', 3, n_Rrect, n_Zrect
! Boxdim(R, m), BOxdim(Z, m), R0(vacuum), Rmin(box, m), Zmid(box, m)
write(eqdsk_unit, '(5E16.9)') R_rect(n_Rrect) - R_rect(1), Z_rect(n_Zrect) - Z_rect(1), RTOR, &
    R_rect(1), 0.5*(Z_rect(1) + Z_rect(n_Zrect))
! Rmagnaxis(m), Zmagnaxis(m)
write(eqdsk_unit, '(5E16.9)') equil_now%coord_sys%position%r(1, 1), equil_now%coord_sys%position%z(1, 1), &
    psi_sgn*SGNIP*equil_now%profiles_1d%psi(1)/GP2, SGNIP*equil_now%profiles_1d%psi(nrho_surf)/GP2, SGNBT*BTOR
write(eqdsk_unit, '(5E16.9)') SGNIP*IPL*1.d6, psi_sgn*SGNIP*equil_now%profiles_1d%psi(1)/GP2, 0., equil_now%coord_sys%position%r(1, 1), 0.
write(eqdsk_unit, '(5E16.9)') equil_now%coord_sys%position%z(1, 1), 0., psi_sgn*SGNIP*equil_now%profiles_1d%psi(nrho_surf)/GP2, 0., 0.
write(eqdsk_unit, '(5E16.9)') (SGNBT*fdia_rect(i), i=1, n_Rrect)
write(eqdsk_unit, '(5E16.9)') (pres_rect(i), i=1, n_Rrect)
write(eqdsk_unit, '(5E16.9)') (psi_sgn*SGNIP*fprime_rect(i)*GP2, i=1, n_Rrect)
write(eqdsk_unit, '(5E16.9)') (psi_sgn*SGNIP*pprime_rect(i)*GP2, i=1, n_Rrect)
write(eqdsk_unit, '(5E16.9)') ((psi_sgn*SGNIP*psi_rect(i, j)/GP2, i=1, n_Rrect), j=1, n_Zrect)
write(eqdsk_unit, '(5E16.9)') (SGNBT*SGNIP*q_rect(i), i=1, n_Rrect)
write(eqdsk_unit, '(2i5)') nthe_surf, nthe_surf
write(eqdsk_unit, '(5E16.9)') (equil_now%coord_sys%position%r(nrho_surf, i), equil_now%coord_sys%position%z(nrho_surf, i), i=1, nthe_surf)
write(eqdsk_unit, '(5E16.9)') (equil_now%coord_sys%position%r(nrho_surf, i), equil_now%coord_sys%position%z(nrho_surf, i), i=1, nthe_surf)
close(eqdsk_unit)

deallocate(psin_rect, pres_rect, fdia_rect, q_rect, pprime_rect, fprime_rect)
deallocate(r_rect, z_rect, psi_rect)

return
end subroutine EQDSK
