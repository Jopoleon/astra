subroutine INTVAR
!-----------------------------------------------------------------------
! Time evolution of the scalar input data
! For the current time, a value is stored in the array
! DEVARX(NCONST) - (description in the file src/for/const.f90)
!
! Input:
!    IVAR, raw_scalar
! Output:
!    DEVARX, DEVAR
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
use outcmn_inc, only: IFDFVX
use const_inc, only: DEVARX, DEVAR, TIME
use expdat, only: raw_scalar
use debugger, only: markloc

implicit none

integer :: jtvar, N1, N2
double precision :: ydt, ydtr, ydtl

call markloc('INTVAR')

N1 = 0
N2 = 0

do jtvar=1, NTVAR
    if (raw_scalar%var_index(jtvar) == 0) EXIT
    N2 = N1
    N1 = raw_scalar%var_index(jtvar)
    if (IFDFVX(N1) >= 0) then
        if (IFDFVX(N1) == 0 .or. N1 /= N2) DEVARX(N1) = raw_scalar%value(jtvar)
        if (N1 == N2) then
            if (TIME >= raw_scalar%time(jtvar-1)) then
                if (TIME < raw_scalar%time(jtvar)) then
                    YDT = raw_scalar%time(jtvar) - raw_scalar%time(jtvar-1)
                    YDTR = (raw_scalar%time(jtvar) - TIME)/YDT
                    YDTL = (TIME - raw_scalar%time(jtvar-1))/YDT
                    DEVARX(N1) = raw_scalar%value(jtvar)*YDTL + raw_scalar%value(jtvar-1)*YDTR
                else
                    DEVARX(N1) = raw_scalar%value(jtvar)
                endif
            endif
        endif
        if (IFDFVX(N1) <= 1) DEVAR(N1) = DEVARX(N1)
    endif
enddo

return
end subroutine INTVAR
