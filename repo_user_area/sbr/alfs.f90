SUBROUTINE alfs(a_lfs, dt_tetop, avdte)

use parameter_inc, only: NRD
use const_inc, only: NA1, ROC, CDWM2
use status_inc, only: TE, AMETR
use parameters_a2equil, only: equil_now

implicit none

double precision, intent(out) :: a_lfs(NRD), dt_tetop(NRD), avdte
integer :: j, i, Nx, Nt, n1, n2
double precision :: te_top, Rmag, Zmag
double precision, dimension(556) :: a_lfs_eq, ametr_eq
double precision, external :: RADIAL, GRAD

nx = SIZE(equil_now%coord_sys%position%r, 1)
nt = SIZE(equil_now%coord_sys%position%r, 2)

Rmag = equil_now%coord_sys%position%r(1, 1)
Zmag = equil_now%coord_sys%position%z(1, 1)
i = MINLOC(abs(equil_now%coord_sys%position%z(nx, 1:nt-1) - Zmag) + &
    abs(equil_now%coord_sys%position%r(nx, 1:nt-1) - MAXVAL(equil_now%coord_sys%position%r(nx, 1:nt-1), 1)), 1)
j = MINLOC(abs(equil_now%coord_sys%position%z(nx, 1:nt-1) - Zmag) + &
    abs(equil_now%coord_sys%position%r(nx, 1:nt-1) - MINVAL(equil_now%coord_sys%position%r(nx, 1:nt-1), 1)), 1)

a_lfs_eq(1:nx) = equil_now%coord_sys%position%r(1:nx, i) - Rmag
AMETR_eq(1:nx) = 0.5*(equil_now%coord_sys%position%r(1:nx, i) - equil_now%coord_sys%position%r(1:nx, j))

call qinterp(ametr_eq(1:nx), a_lfs_eq(1:nx), nx, AMETR(1:NA1), a_lfs(1:NA1), NA1)

te_top = 100.*RADIAL(TE, ROC*CDWM2)
do j=1, NA1
   dt_tetop(J) = -GRAD(TE, J)/GRAD(a_lfs, J)/te_top
enddo

n1 = nint(CDWM2*NA1) + 1
n2 = nint(0.999*NA1) - 1
avdte = SUM(dt_tetop(n1: n2))/dble(n2 + 1 - n1)

RETURN
END SUBROUTINE alfs
