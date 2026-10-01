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

use scalars, only: NA1
use status, only: RHO
use numerical_tools, only: smooth

implicit none

double precision, intent(in) :: ALFA, f_in(*)
double precision, intent(out) :: f_out(*)

call SMOOTH(ALFA, RHO, f_in, NA1, RHO, f_out, NA1)

end subroutine SMEARR

!----------------------------------------------------------------------|
subroutine SMEARRX(ALFA, f_in, f_out)

use scalars, only: NA1
use status, only: NRD, XRHO
use numerical_tools, only: smooth

implicit none

double precision, intent(in) :: ALFA, f_in(NRD)
double precision, intent(out) :: f_out(NRD)

f_out = 0.
call SMOOTH(ALFA, XRHO(1:NA1), f_in(1:NA1), NA1, XRHO(1:NA1), f_out(1:NA1), NA1)

write(*, '(6e16.8)') f_in(1:NA1)
write(*, '(6e16.8)') f_out(1:NA1)
write(*, *) ''

end subroutine SMEARRX

!----------------------------------------------------------------------|
subroutine SMEARR2(ALFA, f_in, f_out)

use scalars, only: NA1, NA1N, NA1E, NA1I, NA1U, NA10, NA11, NA12, NA13, &
    NA14, NA15, NA16, NA17, NA18, NA19
use status, only: NRD, RHO
use numerical_tools, only: smooth

implicit none

double precision, intent(in) :: ALFA, f_in(NRD)
double precision, intent(out) :: f_out(NRD)
integer :: nrho_max

nrho_max = maxval((/ NA1N, NA1E, NA1I, NA1U, NA10, NA11, NA12, NA13, NA14, NA15, NA16, NA17, NA18, NA19 /))
nrho_max = MIN(nrho_max, NA1)
if (nrho_max == 0) then
    write(*, *) 'SMEARR2', nrho_max
    write(*, *) 'No smoothing applied'
else
    call SMOOTH(ALFA, RHO(1: nrho_max), f_in(1: nrho_max), nrho_max, RHO(1:nrho_max), f_out, nrho_max)
endif
f_out(nrho_max+1: NA1) = f_in(nrho_max+1: NA1)

end subroutine SMEARR2

!----------------------------------------------------------------------|
subroutine SMEARR3(ALFA, f_in, f_out)

use scalars, only: NA1, NA1N, NA1E, NA1I, NA1U, NA10, NA11, NA12, NA13, &
    NA14, NA15, NA16, NA17, NA18, NA19
use status, only: NRD, RHO
use numerical_tools, only: smooth

implicit none

double precision, intent(in) :: ALFA, f_in(NRD)
double precision, intent(out) :: f_out(NRD)
integer :: nrho_max, j

nrho_max = maxval((/ NA1N, NA1E, NA1I, NA1U, NA10, NA11, NA12, NA13, NA14, NA15, NA16, NA17, NA18, NA19 /))
nrho_max = MIN(nrho_max, NA1)

f_out(nrho_max+1: NA1) = f_in(nrho_max)

call SMOOTH(ALFA, RHO(1: nrho_max), f_in(1: nrho_max), nrho_max, RHO(1:nrho_max), f_out, nrho_max)

end subroutine SMEARR3
