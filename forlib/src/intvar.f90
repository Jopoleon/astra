subroutine INTVAR
!-----------------------------------------------------------------------
! The time evolution of the input data is taken from
! VARDAT(1, NTVAR) - time   |
! VARDAT(2, NTVAR) - value   |  for variables
! VARDAT(3, NTVAR) - error (not used) |
!
! For the current time, a value is stored in the array
! DEVARX(NCONST) - (description in the file for/variables_x.txt)
! NTVAR total number of the time slices for all variables in a data file
!
! Input:
!    IVAR, INDVAR, IFDFVX, VARDAT, TIME
! Output:
!    DEVARX, DEVAR
!
! N=INDVAR(jvar) - ordinal number of the quantity in DEVAR(N) & DEVARX(N)
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
use expdat
use debugger, only: markloc

implicit none

integer :: jvar, N1, N2
double precision :: ydt, ydtr, ydtl

call markloc('INTVAR')

N1 = 0
N2 = 0

do jvar=1, IVAR
   N2 = N1
   N1 = INDVAR(jvar)
   if (IFDFVX(N1) >= 0) then
      if (IFDFVX(N1) == 0 .or. N1 /= N2) DEVARX(N1) = VARDAT(2, jvar)
      if (N1 == N2) then
         if (TIME >= VARDAT(1, jvar-1)) then
            if (TIME < VARDAT(1, jvar)) then
               YDT = VARDAT(1, jvar) - VARDAT(1, jvar-1)
               YDTR = (VARDAT(1, jvar) - TIME)/YDT
               YDTL = (TIME - VARDAT(1, jvar-1))/YDT
               DEVARX(N1) = VARDAT(2, jvar)*YDTL + VARDAT(2, jvar-1)*YDTR
            else
               DEVARX(N1) = VARDAT(2, jvar)
            endif
         endif
      endif
      if (IFDFVX(N1) <= 1) DEVAR(N1) = DEVARX(N1)
   endif
enddo

return
end subroutine INTVAR
