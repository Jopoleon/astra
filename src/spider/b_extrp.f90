subroutine b_extrp(x0, X1, X2, X3, j1, j2, j3, j, yarr, jerr)

use sp_parameters, only: nrp

implicit none

integer, intent(in) :: j, j1, j2, j3
real*8, intent(in) :: x0, X1, X2, X3
integer, intent(out) :: jerr
real*8, intent(out) :: yarr(*)

real*8 :: u0, u1, u2, u3

jerr = 0 
u1 = yarr(j1)
u2 = yarr(j2)
u3 = yarr(j3)     
call EXTRP2(x0, u0, X1, X2, X3, u1, u2, u3)
if(j > nrp .or. j < 1) then
   write(*, *) 'Equil: in b_extrp j = ',  j, ' out of range 1 <', j, ' < ', nrp
   jerr = 1
else
   yarr(j) = u0
endif

return
end subroutine b_extrp
