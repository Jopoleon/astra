subroutine err_catch_a

use const_inc, only: TIME,  NA1
use status_inc, only: TE, FP, NE, G11
use debugger, only: astra_stop

implicit none

double precision :: time_ext
integer :: j

data time_ext / 0. /
save time_ext   ! counter to use psi as bc stuff

write(*, *) 'TE    Fp    NE    G11 '
write(*, *) te(1)  , fp(1)  , ne(1)  , g11(1)
write(*, *) te(na1), fp(na1), ne(na1), g11(na1)
write(*,*) 'somethings not right, quit run'
stop

return
end subroutine err_catch_a
