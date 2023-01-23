!----------------------------------------------------------------------|
!  Subroutine minimizes the value of functional
!  INTEGRAL(alfa*(dU/dx)**2+(U-F)**2)*dx with respect to U(x).
!  f_in(1:NA1) is a given array on the grid X(1:NA1)
!  The result is a smoothed array f_out(1:NA1) given on the same grid
!	ALFA ~ 0.01*X(NA1)**2  is a regularizator
!	The target function f_out obeys the additional conditions:
!	     df_out/dx(x=0)=0 - cylindrical case
!	     f_out(XN(NA1))=f_in(XO(NA1))
!----------------------------------------------------------------------|

subroutine SMEARR(ALFA, f_in, f_out)

use const_inc, only: NA1
use status_inc, only: RHO

implicit none
double precision, intent(in) :: ALFA, f_in(*)
double precision, intent(out) :: f_out(*)

call SGLAZH(ALFA, NA1, f_in, RHO, NA1, f_out, RHO)
end subroutine SMEARR

!----------------------------------------------------------------------|
subroutine SGLAZH(ALFA, NO, f_in, XO, N, f_out, XN) ! same as SMOOTH

use parameter_inc, only: NRD

implicit none

integer :: NO, N, J, I
double precision, intent(in) :: ALFA, XO(*), f_in(*), XN(*)
double precision, intent(out) :: f_out(*)
double precision :: YF, YX, YP, YQ, YD, FJ, P(NRD)

if (N > NRD .or. NO .le. 0) then
   write(*, *)' >>> SMEARR: array is out of limits'
   stop
endif
if (NO == 1) then
   do j=1, N
      f_out(j) = f_in(1)
   enddo
   return
endif
if (NO == 2) then
   do j=1, N
      f_out(j) = (f_in(2)*(XN(j) - XO(1)) - f_in(1)*(XN(j) - XO(2)))/(XO(2) - XO(1))
   enddo
   return
endif
if (N < 2) then
   write(*, *)' >>> SMEARR: no output grid is provided'
   stop
endif
if (abs(XO(NO)-XN(N)) > XN(N)/N) then
   write(*, *)'>>> SMEARR: grids are not aligned'
   write(*, '(1A23, I4, F8.4)')'     Old grid size/edge', NO, XO(NO)
   write(*, '(1A23, I4, F8.4)')'     New grid size/edge', N, XN(N)
   stop
endif
do j=2, N
   P(j) = ALFA/(XN(j) - XN(j-1))/XO(NO)**2
enddo
P(1)  = 0.
f_out(1) = f_in(1) ! git 0.
I = 1
YF = (f_in(2) - f_in(1))/(XO(2) - XO(1))
YX = 2./(XN(2) + XN(1))
YP = 0.
YQ = 0.
do j=1, N-1
   if (XO(I) <= XN(j)) then
      do
         I = I + 1
         I = min(I, NO)
         if (I == NO .or. XO(I) >= XN(j)) EXIT
      enddo
      YF = (f_in(I) - f_in(I-1))/(XO(I) - XO(I-1))
   endif
   FJ = f_in(I) + YF*(XN(j) - XO(I))
   YD = 1. + YX*(YP + P(j+1))
   P(j) = YX*P(j+1)/YD
   f_out(j) = (FJ + YX*YQ)/YD
   YX = 2./(XN(j+2) - XN(j))
   YP = (1. - P(j))*P(j+1)
   YQ = f_out(j)*P(j+1)
enddo

f_out(N) = f_in(NO)

do j=N-1, 1, -1
   f_out(j) = P(j)*f_out(j+1) + f_out(j)
enddo

return
end subroutine SGLAZH
