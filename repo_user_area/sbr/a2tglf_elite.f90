subroutine a2tglf_elite()

use scalars, only: IPART, NA1
use status, only: RHO, FP_NORM
use parameters_a2equil, only: equil_now
use numerical_tools, only: qinterp

implicit none

integer, parameter :: unit_out=21

integer :: jrho, jthe, nrho_equ, nthe_equ
double precision, allocatable, dimension(:) :: pfn_equ
double precision, allocatable, dimension(:, :) :: RR_as, ZZ_as, Bp_as
character(len=120) :: f_elite

if (IPART == 1) return ! Avoid "before equilibrium"

nrho_equ = SIZE(equil_now%coord_sys%position%r, dim=1)
nthe_equ = SIZE(equil_now%coord_sys%position%r, dim=2)

if (.not. allocated(pfn_equ)) allocate(pfn_equ(nrho_equ))
if (.not. allocated(RR_as)) then
    allocate(RR_as(NA1, nthe_equ))
    allocate(ZZ_as(NA1, nthe_equ))
    allocate(Bp_as(NA1, nthe_equ))
endif
pfn_equ = (equil_now%profiles_1d%psi - equil_now%profiles_1d%psi(1))/(equil_now%profiles_1d%psi(nrho_equ) - equil_now%profiles_1d%psi(1))

! Interpolation on ASTRA rho-grid
do jthe=1, nthe_equ
    call qinterp(pfn_equ, equil_now%coord_sys%position%r(:, jthe), nrho_equ, FP_NORM(1:NA1), RR_as(:, jthe), NA1)
    call qinterp(pfn_equ, equil_now%coord_sys%position%z(:, jthe), nrho_equ, FP_NORM(1:NA1), ZZ_as(:, jthe), NA1)
    call qinterp(pfn_equ, equil_now%coord_sys%bpcell    (:, jthe), nrho_equ, FP_NORM(1:NA1), Bp_as(:, jthe), NA1)
enddo

f_elite = 'tglf/tglf4elite.dat'
open(unit_out, file=TRIM(f_elite))
write(unit_out, '(I0)') nthe_equ
do jthe=1, nthe_equ
    write(unit_out, '(e14.6)') equil_now%coord_sys%position%theta2d(jthe)
enddo
write(unit_out, '(2I0)') NA1, nthe_equ
do jrho=1, NA1
    do jthe=1, nthe_equ
        write(unit_out, '(3e14.6)') RR_as(jrho, jthe), ZZ_as(jrho, jthe), Bp_as(jrho, jthe)
    enddo
enddo
close(unit_out) 

write(*, '(2A)') 'Written ', TRIM(f_elite)

end subroutine a2tglf_elite
