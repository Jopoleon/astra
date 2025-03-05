subroutine alfs(rhot_ped_top, a_lfs, dt_tetop, avdte)

use parameter_inc, only: NRD
use numerical_tools, only: qinterp
use const_inc, only: NA1, RTOR, ROC
use status_inc, only: AMETR, TE
use parameters_a2equil, only: equil_now

implicit none

double precision, intent(in) :: rhot_ped_top
double precision, intent(out) :: avdte
double precision, intent(out), dimension(NRD) :: a_lfs, dt_tetop

integer :: jrho, jthe_eq, neq, meq, jlfs, n1, n2
integer, dimension(2) :: dims
double precision :: te_top, min_z_zmag, zmag
double precision, dimension(556) :: a_lfs_eq, ametr_eq
double precision, external :: GRAD, RADIAL, RFA, AFX

dims = SHAPE(equil_now%coord_sys%position%r)
neq = dims(1)
meq = dims(2)

min_z_zmag = 1000.
zmag = equil_now%coord_sys%position%z(1, 1)
do jthe_eq=1, meq
    if (equil_now%coord_sys%position%r(neq, jthe_eq) > RTOR) then
        if (ABS(equil_now%coord_sys%position%z(neq, jthe_eq) - zmag) < min_z_zmag) then
            jlfs = jthe_eq
            min_z_zmag = ABS(equil_now%coord_sys%position%z(neq, jthe_eq) - zmag)
        endif
    endif
enddo

ametr_eq(1:neq) = 0.5*(equil_now%profiles_1d%r_outboard(1:neq) - equil_now%profiles_1d%r_inboard(1:neq))
call qinterp(ametr_eq(1:neq), a_lfs_eq(1:neq), neq, AMETR(1:NA1), a_lfs(1:NA1), NA1)

! ratio of average pedestal electron temperature gradient to electron temperature at pedestal top
te_top = 10.d1*RADIAL(TE, ROC*rhot_ped_top)
do jrho=1, NA1
    dt_tetop(jrho) = -GRAD(TE, jrho)/GRAD(a_lfs, jrho)/te_top
enddo

n1 = nint(rhot_ped_top*NA1) + 1
n2 = nint(0.999*NA1) - 1
avdte = SUM(dt_tetop(n1:n2))/float(n2 -n1 + 1)

return
end subroutine alfs
