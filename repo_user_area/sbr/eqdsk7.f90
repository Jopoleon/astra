subroutine EQDSK7

use const_inc, only: NA1, RTOR, BTOR, IPL, UPDWN, TIME, TSTART, MEQUIL, NEQUIL, &
    PSIAX, PSIBO, IPEQL, SGNBT, SGNIP
use status_inc, only: SHIF, FP, XRHO
use parameters_a2equil, only: equil_now, GP, GP2
use outcmn_inc, only: AWD, exp_file, equ_file

implicit none

integer, parameter :: Nrrect=256, Nzrect=256, eqdsk_unit=11
integer :: i, j, nrho_surf, nthe_surf, n_Rrect, n_Zrect
double precision :: raxis, zaxis, drhot_eq, dpsin_rect, psib, Rmin, Rmax, zmin, zmax, dr, dz
double precision, allocatable, dimension(:) :: r_rect, z_rect, rhotor1d, pf1d, &
    rhot_eq, pf_eq, psin_eq, pres_eq, fdia_eq, q_eq, pprime_eq, fprime_eq, &
    psin_rect, pres_rect, fdia_rect, q_rect, pprime_rect, fprime_rect, psi_g_norm
double precision, allocatable, dimension(:, :) :: psi_rect
character(len=120) :: f_eqdsk, checks_nico

if (TIME <= TSTART) return

allocate(rhotor1d(NA1+1), pf1d(NA1+1))
rhotor1d(1) = 0.
rhotor1d(2: ) = XRHO(1: NA1)
pf1d(1) = PSIAX
pf1d(2: ) = FP(1:NA1)

n_Rrect = Nrrect
n_Zrect = Nzrect

raxis = RTOR + SHIF(1)
zaxis = UPDWN

nrho_surf = SIZE(equil_now%coord_sys%position%r, 1)
nthe_surf = SIZE(equil_now%coord_sys%position%r, 2)
allocate(pf_eq(nrho_surf), rhot_eq(nrho_surf), psin_eq(nrho_surf), pres_eq(nrho_surf), &
    fdia_eq(nrho_surf), q_eq(nrho_surf), pprime_eq(nrho_surf), fprime_eq(nrho_surf), &
    psi_g_norm(nrho_surf))
allocate(psin_rect(n_Rrect), pres_rect(n_Rrect), fdia_rect(n_Rrect), q_rect(n_Rrect), &
    pprime_rect(n_Rrect), fprime_rect(n_Rrect))
allocate(r_rect(n_Rrect), z_rect(n_Zrect), psi_rect(n_Rrect, n_Zrect))

! rho_tor grids

drhot_eq = 1./(nrho_surf - 1.d0)
rhot_eq = (/ (drhot_eq*(i - 1.d0), i=1, nrho_surf) /)

dpsin_rect = 1./(n_Rrect - 1.d0)
psin_rect = (/ (dpsin_rect*(i - 1.d0), i=1, n_Rrect) /)

! Getting 1d profiles from equilibrium
psin_eq = equil_now%profiles_1d%psi
do i=1, nrho_surf
    psi_g_norm(i) = (psin_eq(i) - psin_eq(1))/(psin_eq(nrho_surf) - psin_eq(1)) 
enddo
psin_eq   = psi_g_norm
pres_eq   = equil_now%profiles_1d%pressure
q_eq      = equil_now%profiles_1d%q
fdia_eq   = equil_now%profiles_1d%F_dia
pprime_eq = equil_now%profiles_1d%pprime
fprime_eq = equil_now%profiles_1d%ffprime
if (IPEQL == 4) then ! SPIDER
    fdia_eq = -fdia_eq
endif

call qinterp(rhotor1d, pf1d, NA1+1, rhot_eq, pf_eq, nrho_surf)

Rmin = MINVAL(equil_now%coord_sys%position%r(nrho_surf, :)) - 0.03 ! 1.08
Rmax = MAXVAL(equil_now%coord_sys%position%r(nrho_surf, :)) + 0.03 ! 2.26
zmin = MINVAL(equil_now%coord_sys%position%z(nrho_surf, :)) - 0.03 ! -1.0
zmax = MAXVAL(equil_now%coord_sys%position%z(nrho_surf, :)) + 0.03 ! 1.0
dr = (Rmax - Rmin)/(n_Rrect - 1.d0)
dz = (zmax - zmin)/(n_Zrect - 1.d0)
r_rect = (/ (Rmin + dr*(i - 1.d0), i=1, n_Rrect) /)
z_rect = (/ (zmin + dz*(i - 1.d0), i=1, n_Zrect) /)

! Biquadratic interpolation from psi(rho,theta) to psi(R, Z)
call ctr2rz_fun(nrho_surf, nthe_surf, pf_eq, &
    equil_now%coord_sys%position%r, equil_now%coord_sys%position%z, n_Rrect, n_Zrect, r_rect, z_rect, psi_rect)

! Quadratic interpolation ro n_Rrect grid (eqdsk requires that)
call qinterp(psin_eq,   pres_eq, nrho_surf, psin_rect,   pres_rect, n_Rrect)
call qinterp(psin_eq,   fdia_eq, nrho_surf, psin_rect,   fdia_rect, n_Rrect)
call qinterp(psin_eq, pprime_eq, nrho_surf, psin_rect, pprime_rect, n_Rrect)
call qinterp(psin_eq, fprime_eq, nrho_surf, psin_rect, fprime_rect, n_Rrect)
call qinterp(psin_eq,      q_eq, nrho_surf, psin_rect,      q_rect, n_Rrect)

! EQDSk file output
write(f_eqdsk, '(4A, f5.3, A)') TRIM(awd), '/ncdf_out/', TRIM(exp_file), TRIM(equ_file), TIME, '.eqdsk'
write(*, *) 'Storing '//TRIM(f_eqdsk), '   nR =', n_Rrect

open(eqdsk_unit, file=TRIM(f_eqdsk))
write(eqdsk_unit, '(A48, 3i4)') 'ASTRA', 3, n_Rrect, n_Zrect
! Boxdim(R, m), BOxdim(Z, m), R0(vacuum), Rmin(box, m), Zmid(box, m)
write(eqdsk_unit, '(5E16.9)') R_rect(n_Rrect) - R_rect(1), Z_rect(n_Zrect) - Z_rect(1), RTOR, &
    R_rect(1), 0.5*(Z_rect(1) + Z_rect(n_Zrect))
! Rmagnaxis(m), Zmagnaxis(m)
write(eqdsk_unit, '(5E16.9)') equil_now%coord_sys%position%r(1, 1), equil_now%coord_sys%position%z(1, 1), -SGNIP*PSIAX/GP2, -SGNIP*PSIBO/GP2, SGNBT*BTOR
write(eqdsk_unit, '(5E16.9)') SGNIP*IPL*1.d6, -SGNIP*PSIAX/GP2, 0., equil_now%coord_sys%position%r(1, 1), 0.
write(eqdsk_unit, '(5E16.9)') equil_now%coord_sys%position%z(1, 1), 0., SGNIP*PSIBO/GP2, 0., 0.
write(eqdsk_unit, '(5E16.9)') (SGNBT*fdia_rect(i), i=1, n_Rrect)
write(eqdsk_unit, '(5E16.9)') (pres_rect(i), i=1, n_Rrect)
write(eqdsk_unit, '(5E16.9)') (-SGNIP*fprime_rect(i)*GP2, i=1, n_Rrect)
write(eqdsk_unit, '(5E16.9)') (-SGNIP*pprime_rect(i)*GP2, i=1, n_Rrect)
write(eqdsk_unit, '(5E16.9)') ((-SGNIP*psi_rect(i, j)/GP2, i=1, n_Rrect), j=1, n_Zrect)
write(eqdsk_unit, '(5E16.9)') (SGNBT*SGNIP*q_rect(i), i=1, n_Rrect)
write(eqdsk_unit, '(2i5)') nthe_surf, nthe_surf
write(eqdsk_unit, '(5E16.9)') (equil_now%coord_sys%position%r(nrho_surf, i), equil_now%coord_sys%position%z(nrho_surf, i), i=1, nthe_surf)
write(eqdsk_unit, '(5E16.9)') (equil_now%coord_sys%position%r(nrho_surf, i), equil_now%coord_sys%position%z(nrho_surf, i), i=1, nthe_surf)
close(eqdsk_unit)

deallocate(pf_eq, rhot_eq, psin_eq, pres_eq, fdia_eq, q_eq, pprime_eq, fprime_eq)
deallocate(psin_rect, pres_rect, fdia_rect, q_rect, pprime_rect, fprime_rect)
deallocate(r_rect, z_rect, psi_rect, psi_g_norm)


return
end subroutine EQDSK7
