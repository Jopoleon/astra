!----------------------------------------------------------------------|
! Subroutine minimizes the value of functional
! INTEGRAL(alfa*(dU/dx)**2+(U-F)**2)*dx with respect to U(x).
! f_in(1:NA1) is a given array on the grid X(1:NA1)
! The result is a smoothed array f_out(1:NA1) given on the same grid
!    ALFA ~ 0.01*X(NA1)**2  is a regularizator
!    The target function f_out obeys the additional conditions:
!       df_out/dx(x=0)=0 - cylindrical case
!       f_out(XN(NA1))=f_in(XO(NA1))
!----------------------------------------------------------------------|
subroutine SMEARR(ALFA, f_in, f_out)

use const_inc, only: NA1
use status_inc, only: RHO

implicit none

double precision, intent(in) :: ALFA, f_in(*)
double precision, intent(out) :: f_out(*)

call SGLAZH(ALFA, NA1, f_in, RHO, NA1, f_out, RHO)

return
end subroutine SMEARR

!----------------------------------------------------------------------|
subroutine SMEARRX(ALFA, f_in, f_out)

use parameter_inc, only: NRD
use const_inc, only: NA1
use status_inc, only: XRHO

implicit none

double precision, intent(in) :: ALFA, f_in(NRD)
double precision, intent(out) :: f_out(NRD)

f_out = 0.
call SGLAZH(ALFA, NA1, f_in(1:NA1), XRHO(1:NA1), NA1, f_out(1:NA1), XRHO(1:NA1))

write(*, '(6e16.8)') f_in(1:NA1)
write(*, '(6e16.8)') f_out(1:NA1)
write(*, *) ''

return
end subroutine SMEARRX

!----------------------------------------------------------------------|
subroutine SMEARR2(ALFA, f_in, f_out)

use parameter_inc, only: NRD
use const_inc, only: NA1, NA1N, NA1E, NA1I, NA1U, NA10, NA11, NA12, NA13, &
    NA14, NA15, NA16, NA17, NA18, NA19
use status_inc, only: RHO

implicit none

double precision, intent(in) :: ALFA, f_in(NRD)
double precision, intent(out) :: f_out(NRD)
integer :: nrho_max

nrho_max = maxval((/ NA1N, NA1E, NA1I, NA1U, NA10, NA11, NA12, NA13, NA14, NA15, NA16, NA17, NA18, NA19 /))
nrho_max = MIN(nrho_max, NA1)

call SGLAZH(ALFA, nrho_max, f_in(1: nrho_max), RHO(1: nrho_max), nrho_max, f_out, RHO(1:nrho_max))

f_out(nrho_max+1: NA1) = f_in(nrho_max+1: NA1)

return
end subroutine SMEARR2

!----------------------------------------------------------------------|
subroutine SMEARR3(ALFA, f_in, f_out)

use parameter_inc, only: NRD
use const_inc, only: NA1, NA1N, NA1E, NA1I, NA1U, NA10, NA11, NA12, NA13, &
    NA14, NA15, NA16, NA17, NA18, NA19
use status_inc, only: RHO

implicit none

double precision, intent(in) :: ALFA, f_in(NRD)
double precision, intent(out) :: f_out(NRD)
integer :: nrho_max, j

nrho_max = maxval((/ NA1N, NA1E, NA1I, NA1U, NA10, NA11, NA12, NA13, NA14, NA15, NA16, NA17, NA18, NA19 /))
nrho_max = MIN(nrho_max, NA1)

f_out(nrho_max+1: NA1) = f_in(nrho_max)

call SGLAZH(ALFA, nrho_max, f_in(1: nrho_max), RHO(1: nrho_max), nrho_max, f_out, RHO(1:nrho_max))

return
end subroutine SMEARR3

!----------------------------------------------------------------------|
subroutine SGLAZH(ALFA, n_in, f_in, x_in, n_out, f_out, x_out) ! same as SMOOTH

use parameter_inc, only: NRD

implicit none

integer, intent(in) :: n_in, n_out
double precision, intent(in) :: ALFA, x_in(*), f_in(*), x_out(*)
double precision, intent(out) :: f_out(*)
integer :: i, j
double precision :: YF, YX, YP, YQ, YD, FJ, P(NRD)

if (n_out > NRD .or. n_in .le. 0) then
    write(*, *)' >>> SMEARR: array is out of limits'
    stop
endif
if (n_in == 1) then
    do j=1, n_out
        f_out(j) = f_in(1)
    enddo
    return
endif
if (n_in == 2) then
    do j=1, n_out
        f_out(j) = (f_in(2)*(x_out(j) - x_in(1)) - f_in(1)*(x_out(j) - x_in(2)))/(x_in(2) - x_in(1))
    enddo
    return
endif
if (n_out < 2) then
    write(*, *)' >>> SMEARR: no output grid is provided'
    stop
endif
if (abs(x_in(n_in) - x_out(n_out)) > x_out(n_out)/n_out) then
    write(*, *)'>>> SMEARR: grids are not aligned'
    write(*, '(1A23, i4, F8.4)')'     Old grid size/edge', n_in, x_in(n_in)
    write(*, '(1A23, i4, F8.4)')'     New grid size/edge', n_out, x_out(n_out)
    stop
endif
do j=2, n_out
    P(j) = ALFA/(x_out(j) - x_out(j-1))/x_in(n_in)**2
enddo
P(1)  = 0.
f_out(1) = f_in(1) ! git 0.
i = 1
YF = (f_in(2) - f_in(1))/(x_in(2) - x_in(1))
YX = 2./(x_out(2) + x_out(1))
YP = 0.
YQ = 0.
do j=1, n_out-1
    if (x_in(i) <= x_out(j)) then
        do
            i = i + 1
            i = min(i, n_in)
            if (i == n_in .or. x_in(i) >= x_out(j)) EXIT
        enddo
        YF = (f_in(i) - f_in(i-1))/(x_in(i) - x_in(i-1))
    endif
    FJ = f_in(i) + YF*(x_out(j) - x_in(i))
    YD = 1. + YX*(YP + P(j+1))
    P(j) = YX*P(j+1)/YD
    f_out(j) = (FJ + YX*YQ)/YD
    YX = 2./(x_out(j+2) - x_out(j))
    YP = (1. - P(j))*P(j+1)
    YQ = f_out(j)*P(j+1)
enddo

f_out(n_out) = f_in(n_in)

do j=n_out-1, 1, -1
    f_out(j) = P(j)*f_out(j+1) + f_out(j)
enddo

return
end subroutine SGLAZH
