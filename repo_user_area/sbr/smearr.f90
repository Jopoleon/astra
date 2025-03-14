subroutine SMEARR(ALFA, f_in, f_out)

use parameter_inc, only: NRD
use const_inc, only: NA1
use status_inc, only: RHO

implicit none

double precision, intent(in) :: ALFA, f_in(NRD)
double precision, intent(out) :: f_out(NRD)

f_out = 0.
call SMOOTH(ALFA, NA1, f_in(1:NA1), RHO(1:NA1), NA1, f_out(1:NA1), RHO(1:NA1))

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
call SMOOTH(ALFA, NA1, f_in(1:NA1), XRHO(1:NA1), NA1, f_out(1:NA1), XRHO(1:NA1))

write(*, *) 'SMEARRX'
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

call SMOOTH(ALFA, nrho_max, f_in(1: nrho_max), RHO(1: nrho_max), nrho_max, f_out(1:nrho_max), RHO(1:nrho_max))

f_out(nrho_max+1: NA1) = f_in(nrho_max+1: NA1)

return
end subroutine SMEARR2
