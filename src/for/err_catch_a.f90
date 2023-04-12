subroutine err_catch_a

use const_inc, only: TIME, CDHJ7, CV6, NA1
use status_inc, only: TE, FP, NE, G11
use debugger, only: astra_stop, flightsim

implicit none

double precision :: time_ext
integer :: j 

data time_ext / 0. /
save time_ext   ! counter to use psi as bc stuff

if (flightsim == 1) then
    do j=1, 1
        if (j == 1) time_ext = TIME
        write(*, *) ' STOP THE SIMULINK RUN!'
        write(*, *) 'fporse  ' 
        CDHJ7 = 1. ! astra has crashed, so simulink should stop with this error variable
        TIME = time_ext
        call astra_stop
    enddo
else
    write(*, *) 'TE    Fp    NE    G11 '
    write(*, *) te(1)  , fp(1)  , ne(1)  , g11(1)
    write(*, *) te(na1), fp(na1), ne(na1), g11(na1)
    write(*,*) 'somethings not right, pausing'
    stop
endif

return
end subroutine err_catch_a
