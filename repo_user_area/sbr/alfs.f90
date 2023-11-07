subroutine alfs(rhot_ped_top, aminlfsa, dt_tetop, avdte)

use parameter_inc, only: NRD
use numerical_tools, only: qinterp
use const_inc, only: NA1, AWALL, RTOR
use status_inc, only: AMETR, TE
use parameters_a2equil, only: equil_now

implicit none

double precision, intent(in) :: rhot_ped_top
double precision, intent(out) :: avdte
double precision, intent(out), dimension(NRD) :: aminlfsa, dt_tetop

integer :: j, KK, neq, meq, ispan
integer, dimension(2) :: dims
double precision :: tmp
double precision, dimension(556) :: aminlfs, aminm
double precision, external :: GRAD, RADIAL, RFA, AFX

dims = SHAPE(equil_now%coord_sys%position%r)
neq = dims(1)
meq = dims(2)
ispan = minloc(abs(equil_now%coord_sys%position%r(neq, 1:meq) - (RTOR+AWALL)), 1) ! Z=0, lfs
aminlfs(1:neq) = equil_now%coord_sys%position%rmin(1:neq, ispan)
aminm(1:neq) = 0.5*(equil_now%profiles_1d%r_outboard(1:neq) - equil_now%profiles_1d%r_inboard(1:neq)) 

call qinterp(aminm(1:neq), aminlfs(1:neq), neq, AMETR(1:NA1), aminlfsa(1:NA1), NA1)

! ratio of average pedestal electron temperature gradient to electron temperature at pedestal top 
do j=1, NA1
    dt_tetop(J) = (-GRAD(TE, J)/GRAD(aminlfsa, J)/10.d1/RADIAL(TE, RFA(AFX(rhot_ped_top))))
enddo

KK  = 0
tmp = 0.

do j=1, NA1
    if ((j > nint(rhot_ped_top*NA1)) .and. (j < nint(0.999*NA1))) then
        tmp = tmp + dt_tetop(J)
        KK = KK + 1
    else
        continue
    endif
enddo

avdte = tmp/KK


return
end subroutine alfs
