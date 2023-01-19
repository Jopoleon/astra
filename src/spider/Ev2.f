!----------------------------------------------------------------
! NVAR=1 - SQUARED FORM OF SURFACES IS USED:
! A1*X**2+A2*X*Y+A3*Y**2+A4*X+A5*Y+A6=0
!
! NVAR =  2 - BILINEAR FORM OF SURFACES IS USED
! NVAR = 11 - TEST FOR NVAR=1
!----------------------------------------------------------------
      subroutine FUN(X, F)

      implicit none

      real*8, intent(in) :: X
      real*8, intent(out) :: F

      integer :: NVAR, IFAILF
      real*8 :: A1, A2, A3, A4, A5, A6, ALF, B, DET, Y2, F2
      real*8 :: RCP, ZCP, UU(5), DF, Y
      COMMON/CPCOOR/ NVAR, IFAILF, RCP, ZCP, UU, DF
      COMMON/YSCR/ Y

      if(NVAR .EQ. 1) then

         A1 = 0.5D0*UU(3)
         A2 =       UU(4)
         A3 = 0.5D0*UU(5)
         A4 =       UU(1)
         A5 =       UU(2)
         A6 = DF
         ALF = 1.D0
         IFAILF = 0

         B = (A2*X+A5) / (2.D0*A3)
         DET = B**2 - (A1*X**2+A4*X+A6)/A3
         if(DET .LT. 0) then
            DET = - DET
            ALF = 1.D+3
            IFAILF = 1
         endif
         DET = dsqrt( DET ) * ALF
         Y  =  -B + DET
         Y2 =  -B - DET
         F  = (X-RCP)**2 + (Y -ZCP)**2
         F2 = (X-RCP)**2 + (Y2-ZCP)**2
         if (F2 .LT. F ) then
            F = F2
            Y = Y2
         endif
      endif

      if(NVAR .EQ. 11) then
         A1 = 1.D0
         A2 = 0.D0
         A3 = 1.D0
         A4 = 0.D0
         A5 = 0.D0
         A6 = -1.D0
         B = (A2*X+A5) / (2.D0*A3)
         DET = dsqrt( B**2 - (A1*X**2+A4*X+A6)/A3 )
         Y  =  -B + DET
         Y2 =  -B - DET
         F  = (X-RCP)**2 + (Y -ZCP)**2
         F2 = (X-RCP)**2 + (Y2-ZCP)**2
         F = min(F, F2)
      endif

      return
      end subroutine FUN
