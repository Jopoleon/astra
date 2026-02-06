subroutine INTVAR
!-----------------------------------------------------------------------
! Time evolution of the scalar input data
! For the current time, a value is stored in the array
! varxValues(NTVAR) - (description in the file astra_variables.json)
!
! Output:
!    varxValues, varValues
!
! IFDFVX(N)  - type of variable
! IFDFVX:
!    = 0 - determined by the data file (independent on time),
!    = 1 - determined by the data file (dependent on time),
!    = 2 - determined by the MODEL,
!    = 3 - Keyboard
!    = 4 - any change of the variable is forbidden 
!          eg. (AB, RTOR, ELONM, TRICH or set interactively)
!-----------------------------------------------------------------------

use parameter_inc, only: NTVAR
use io_mod, only: IFDFVX
use const_inc, only: varxValues, varValues, TIME
use exp_data, only: raw_scalars
use debugger, only: markloc

implicit none

integer :: jtvar, N1, N2
double precision :: ydt, ydtr, ydtl

call markloc('INTVAR')

N1 = 0
N2 = 0

do jtvar=1, NTVAR
    if (raw_scalars%var_index(jtvar) == 0) EXIT
    N2 = N1
    N1 = raw_scalars%var_index(jtvar)
    if (IFDFVX(N1) >= 0) then
        if (IFDFVX(N1) == 0 .or. N1 /= N2) varxValues(N1) = raw_scalars%data(jtvar)
        if (N1 == N2) then
            if (TIME >= raw_scalars%time(jtvar-1)) then
                if (TIME < raw_scalars%time(jtvar)) then
                    YDT = raw_scalars%time(jtvar) - raw_scalars%time(jtvar-1)
                    YDTR = (raw_scalars%time(jtvar) - TIME)/YDT
                    YDTL = (TIME - raw_scalars%time(jtvar-1))/YDT
                    varxValues(N1) = raw_scalars%data(jtvar)*YDTL + raw_scalars%data(jtvar-1)*YDTR
                else
                    varxValues(N1) = raw_scalars%data(jtvar)
                endif
            endif
        endif
        if (IFDFVX(N1) <= 1) varValues(N1) = varxValues(N1)
    endif
enddo

return
end subroutine INTVAR
