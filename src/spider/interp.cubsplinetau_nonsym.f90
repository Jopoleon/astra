subroutine INTRPTAU(PXIN,PYIN,PYINNEW,PY2,KNIN,PXOUT,PYOUT,PYOUTP, &
    PYOUTPP,KNOUT,KOPT,PTAUS,PWORK,PAMAT,MDAMAT,PBCLFT,PBCRGT)
!=================================================================
!
!  NOTE: THIS ROUTINE INCLUDES THE STANDARD CUBIC SPLINE IF PTAUS=0:
!        THEN PYINNEW IS NOT USED AND MDAMAT=3 IS SUFFICIENT
!        (=> PYINNEW(1) OR PYINNEW=PYIN IS OK)
!
!  Interpolate (pxin,pyin) on (pxout,pyout) using
!  Hirshman fitted cubic spline with ptaus value or
!  standard cubic spline if PTAUS=0
!
!  KOPT = 0: ONLY INTERPOLATE FUNCTION INTO PYOUT
!  KOPT = 1: INTERPOLATE FUNCTION INTO PYOUT AND 1ST DER. INTO PYOUTP
!  KOPT = 2: AS KOPT=1 PLUS 2ND DER. INTO PYOUTPP
!
!  BOUNDARY CONDITIONS FOR CUBIC SPLINE: PBCLFT AND PBCRGT
!  THEN IF
!    = 2  : B.C.: 2ND DERIVATIVE (NATURAL)
!    = 1  : B.C.: 1ST DERIVATIVE FROM CUBIC LAGRANGE INTERPOLATION
!    = 0  : B.C.: 1ST DERIVATIVE = 0.0
!    = 3  : B.C.: 1ST DERIVATIVE FROM LINEAR INTERPOLATION
!  OTHERWISE: B.C: 1ST DERIVATIVE = PBCLFT OR PBCRGT

implicit none

integer, intent(in) :: knin, knout, kopt, mdamat
real*8, intent(in) :: ptaus, pbclft, pbcrgt
real*8, intent(in), dimension(knin) :: pxin, pyin, py2, pwork, pyinnew
real*8, intent(in), dimension(knout) :: pxout
real*8, intent(out), dimension(knout) :: pyout, pyoutp, pyoutpp
real*8, intent(out), dimension(mdamat, knin) :: pamat

integer :: i
real*8 :: zeps, zy, zyp, zyp1, zypn, zypp

ZEPS = 1.0D-12

!L  1. PREPARE SPLINE (DEFAULT <=> PKBCOPT.EQ.0)
!   LEFT

ZYP1 = PBCLFT
if (ABS(PBCLFT - 1.D0) .LE. ZEPS) ZYP1 = -1.D+32
if (ABS(PBCLFT - 2.D0) .LE. ZEPS) ZYP1 =  1.D+32
if (ABS(PBCLFT - 3.D0) .LE. ZEPS) ZYP1 = -1.D+34

!   RIGHT

ZYPN = PBCRGT
if (ABS(PBCRGT - 1.D0) .LE. ZEPS) ZYPN = -1.D+32
if (ABS(PBCRGT - 2.D0) .LE. ZEPS) ZYPN=  1.D+32
if (ABS(PBCRGT - 3.D0) .LE. ZEPS) ZYPN = -1.D+34

call CUBSPLFIT(PXIN,PYIN,PYINNEW,KNIN,ZYP1,ZYPN,PY2,PTAUS,PWORK,PAMAT,MDAMAT)

!L    2. COMPUTE INTERPOLATED VALUE AT EACH PXOUT

do I=1,KNOUT
  if (PTAUS .EQ. 0.D0) then
    call SPLINTEND(PXIN,PYIN,PY2,KNIN,PXOUT(I),ZY,ZYP,ZYPP,1)
  else
    call SPLINTEND(PXIN,PYINNEW,PY2,KNIN,PXOUT(I),ZY,ZYP,ZYPP,1)
  endif
  PYOUT(I) = ZY
  if (KOPT .GE. 1) PYOUTP(I) = ZYP
  if (KOPT .EQ. 2) PYOUTPP(I) = ZYPP
enddo

return
end subroutine INTRPTAU

!----------------------------------------------------------------
real*8 function FA3(A1, A2, A3, A4, B1, B2, B3, B4)

implicit none

real*8, intent(in) :: A1, A2, A3, A4, B1, B2, B3, B4

FA3 = (A1-A2) / ((B1-B2)*(B2-B4)*(B2-B3)) + &
      (A1-A3) / ((B4-B3)*(B3-B1)*(B3-B2)) + &
      (A1-A4) / ((B1-B4)*(B2-B4)*(B3-B4))

return
end function FA3

!----------------------------------------------------------------
real*8 function FA2(A1, A2, A3, A4, B1, B2, B3, B4)

implicit none

real*8, intent(in) :: A1, A2, A3, A4, B1, B2, B3, B4
real*8, external :: FA3

FA2 = (A1-A2) / ((B2-B1)*(B3-B2)) + &
      (A3-A1) / ((B3-B1)*(B3-B2)) - &
      (B1+B2+B3) * FA3(A1, A2, A3, A4, B1, B2, B3, B4)

return
end function FA2

!----------------------------------------------------------------
real*8 function FA1(A1, A2, A3, A4, B1, B2, B3, B4)

implicit none

real*8, intent(in) :: A1, A2, A3, A4, B1, B2, B3, B4
real*8, external :: FA2, FA3

FA1 = (A1-A2) / (B1-B2) - &
      (B1+B2) * FA2(A1,A2,A3,A4,B1,B2,B3,B4) - &
      (B1*B1+B1*B2+B2*B2) * FA3(A1, A2, A3, A4, B1, B2, B3, B4)

return
end function FA1

!----------------------------------------------------------------
real*8 function FCCCC1(A1, A2, A3, A4, B1, B2, B3, B4, PX)

implicit none

real*8, intent(in) :: A1, A2, A3, A4, B1, B2, B3, B4, PX
real*8, external :: FA1, FA2, FA3

FCCCC1 = FA1(A1, A2, A3, A4, B1, B2, B3, B4) + &
         PX * (2.D0 * FA2(A1, A2, A3, A4, B1, B2, B3, B4) + &
               3.D0 * PX * FA3(A1, A2, A3, A4, B1, B2, B3, B4))

return
end function FCCCC1

!----------------------------------------------------------------
subroutine CUBSPLFIT(X, Y, YNEW, N, YP1, YPN, Y2, TAUS, WORK, AMAT, MDAMAT)

!     PREPARE SECOND DERIVATIVE OF CUBIC SPLINE INTERPOLATION AND NEW
!     VALUES OF Y AT NODES YNEW FITTED SUCH THAT CHI**2 + TAUS*F''**2
!     IS MINIMIZED ACCORDING TO HIRSHMAN ET AL, PHYS. PLASMAS 1 (1994) 2280.
!     TAUS = TAU*SIGMA_K OF PAPER ASSUMING SIGMA_K CONSTANT.
!
!     SETTING TAUS=0., ONE FINDS THE USUAL CUBIC SPLINE INT. WITH CHI**2=0
!     TAUS LARGE => FIT CLOSER TO STRAIGHT LINE (SECOND DERIV.=0)
!
!     IF LAPACK ROUTINES NOT AVAILABLE, USE NONSYM.F AND REMOVE "c%nonsym"
!
!     IF TAUS=0, YNEW NOT USED => YNEW(1) OR YNEW=Y IS OK

implicit none

integer, parameter :: MNONSYM=501

integer, intent(in) :: n, mdamat
real*8, intent(in), dimension(n) :: x, y
real*8, intent(out), dimension(n) :: y2, work, ynew
real*8, intent(out), dimension(mdamat, n) :: amat

integer :: i, j, k, iband, idiag, idima, idimrhs, info2, irhs, iup
real*8 :: taus, ztaueff, yp1, ypn, zyp1, zypn
integer, dimension(mnonsym) :: ipivot
real*8, dimension(mnonsym*7) :: anonsym
real*8, external :: FCCCC1

! ----------------------------------------------------------------------
! --     STATEMENT FUNCTION FOR CUBIC INTERPOLATION                   --
! --                         23.04.88            AR        CRPP       --
! --                                                                  --
! -- CUBIC INTERPOLATION OF A FUNCTION F(X)                           --
! -- THE EIGHT ARGUMENTS A1,A2,A3,A4,B1,B2,B3,B4 ARE DEFINED BY:      --
! -- F(B1) = A1 , F(B2) = A2 , F(B3) = A3 , F(B4) = A4                --
! ----------------------------------------------------------------------
! -- FCCCC0 GIVES THE VALUE OF THE FUNCTION AT POINT PX:              --
! -- FCCCC0(......,PX) = F(PX)                                        --
! ----------------------------------------------------------------------
!        FCCCC0(A1,A2,A3,A4,B1,B2,B3,B4,PX) =
!     F              FA0(A1,A2,A3,A4,B1,B2,B3,B4) +
!     F              PX * (FA1(A1,A2,A3,A4,B1,B2,B3,B4) +
!     F                    PX * (FA2(A1,A2,A3,A4,B1,B2,B3,B4) +
!     F                          PX * FA3(A1,A2,A3,A4,B1,B2,B3,B4)))
! ----------------------------------------------------------------------
! -- FCCCC1 GIVES THE VALUE OF THE FIRST DERIVATIVE OF F(X) AT PX:    --
! -- FCCCC1(......,PX) = DF/DX (PX)                                   --
! ----------------------------------------------------------------------
! -- FCCCC2 GIVES THE VALUE OF THE SECOND DERIVATIVE OF F(X) AT PX:   --
! -- FCCCC2(......,PX) = D2F/DX2 (PX)                                 --
! ----------------------------------------------------------------------
!         FCCCC2(A1,A2,A3,A4,B1,B2,B3,B4,PX) =
!     F             2.D0 * FA2(A1,A2,A3,A4,B1,B2,B3,B4) +
!     F             6.D0 * FA3(A1,A2,A3,A4,B1,B2,B3,B4) * PX
! ----------------------------------------------------------------------
! -- FCCCC3 GIVES THE VALUE OF THE THIRD DERIVATIVE OF F(X) AT PX:     -
! -- FCCCC3(......,PX) = D3F/DX3 (PX)                                  -
! ----------------------------------------------------------------------
!         FCCCC3(A1,A2,A3,A4,B1,B2,B3,B4,PX) =
!     F                      6.D0* FA3(A1,A2,A3,A4,B1,B2,B3,B4)
!-----------------------------------------------------------------------

do I=1,N
  Y2(I) = 0.D0
enddo
do I=1,MDAMAT
  do J=1,N
    AMAT(I,J) = 0.D0
  enddo
enddo

!     PREPARE 1 / H_K

do K=1,N-1
  WORK(K) = 1.D0/ (X(K+1) - X(K))
enddo

ZTAUEFF = TAUS

!     PREPARE BAND WIDTH

IUP = 2
if (ZTAUEFF .EQ. 0.D0) IUP = 1
IBAND = 2*IUP + 1
!%nonsym
IDIAG = IBAND
if (MNONSYM*7 .LT. (3*IUP+1)*N) then
  PRINT *,' MMNONSYM*7= ',MNONSYM*7,' < (3*IUP+1)*N= ',(3*IUP+1)*N
  STOP
endif
!%nonsym
!
!     CONSTRUCT MATRIX AND R.H.S
!
!.......................................................................
!     AS MATRIX SYMMETRIC, COMPUTE ONLY UPPER PART
!
!     K=2,N-2 (BULK PART)

do K=2,N-2
!     A(K,K)
  AMAT(IDIAG,K) = (1.D0/WORK(K)+1.D0/WORK(K-1)) / 3.D0 &
+ 2.D0*ZTAUEFF*(WORK(K)**2 + WORK(K)*WORK(K-1)+WORK(K-1)**2)
!     A(K,K+1)
  AMAT(IDIAG-1,K+1) = 1.D0/WORK(K)/6.D0 &
- ZTAUEFF*WORK(K)*(WORK(K+1)+2.D0*WORK(K)+WORK(K-1))
!     A(K,K+2)
  if (IUP .EQ. 2) AMAT(IDIAG-2,K+2) = ZTAUEFF*WORK(K+1)*WORK(K)
!     B(K)
  Y2(K) = (Y(K+1)-Y(K))*WORK(K) - (Y(K)-Y(K-1))*WORK(K-1)

enddo

K = 1
AMAT(IDIAG,K) = 1.D0/WORK(K) / 3.D0 &
    + 2.D0*ZTAUEFF*WORK(K)**2
AMAT(IDIAG-1,K+1) = 1.D0/WORK(K)/6.D0 &
    - ZTAUEFF*WORK(K)*(WORK(K+1)+2.D0*WORK(K))
if (IUP .EQ. 2) AMAT(IDIAG-2,K+2) = ZTAUEFF*WORK(K+1)*WORK(K)
Y2(K) = (Y(K+1)-Y(K))*WORK(K)

K = N-1
AMAT(IDIAG,K) = (1.D0/WORK(K)+1.D0/WORK(K-1)) / 3.D0 &
    + 2.D0*ZTAUEFF*(WORK(K)**2 + WORK(K)*WORK(K-1)+WORK(K-1)**2)
AMAT(IDIAG-1,K+1) = 1.D0/WORK(K)/6.D0 &
    - ZTAUEFF*WORK(K)*(2.D0*WORK(K)+WORK(K-1))
Y2(K) = (Y(K+1)-Y(K))*WORK(K) - (Y(K)-Y(K-1))*WORK(K-1)

K = N
AMAT(IDIAG,K) = 1.D0/WORK(K-1) / 3.D0 &
+ 2.D0*ZTAUEFF*WORK(K-1)**2
Y2(K) = - (Y(K)-Y(K-1))*WORK(K-1)

!.......................................................................
!     BOUNDARY CONDITIONS
!
if (YP1 .GT. .99D30) then
!     SECOND DERIVATIVE = 0 (NATURAL B.C.) AND SYMMETRIZE
  Y2(1) = 0.D0
  AMAT(IDIAG,1) = 1.D0
  AMAT(IDIAG-1,2) = 0.D0
  if (IUP .EQ. 2) AMAT(IDIAG-2,3) = 0.D0
  ZYP1 = 0.D0
else if (-YP1 .GT. .99D+34) then
  ZYP1 = (Y(2)-Y(1))/(X(2)-X(1))
else if (-YP1 .GT. .99D+30) then
  ZYP1 = FCCCC1(Y(1),Y(2),Y(3),Y(4),X(1),X(2),X(3),X(4),X(1))
else
  ZYP1 = YP1
endif
if (YPN .GT. .99D30) then
!     SECOND DERIVATIVE = 0 (NATURAL B.C.) AND SYMMETRIZE
  Y2(N) = 0.D0
  AMAT(IDIAG,N) = 1.D0
  AMAT(IDIAG-1,N) = 0.D0
  if (IUP .EQ. 2) AMAT(IDIAG-2,N) = 0.D0
  ZYPN = 0.D0
else if (-YPN .GT. .99D+34) then
  ZYPN = (Y(N)-Y(N-1))/(X(N)-X(N-1))
else if (-YPN .GT. .99D+30) then
  ZYPN = FCCCC1(Y(N-3),Y(N-2),Y(N-1),Y(N), &
X(N-3),X(N-2),X(N-1),X(N),X(N))
else
  ZYPN = YPN
endif

Y2(1) = Y2(1) - ZYP1
Y2(N) = Y2(N) + ZYPN

!     SOLVE SYSTEM

IDIMA = MDAMAT
IDIMRHS = N
IRHS = 1

!%nonsym
do I=1,IBAND*N
  ANONSYM(I) = 0.0
enddo
do I=1,N
!     UPPER PART
  do J=max(1,i),min(i+IUP,n)
    ANONSYM((I-1)*IBAND+J-I+IUP+1) = AMAT(IDIAG+I-J,J)
  enddo
!     LOWER PART
  do J=max(1,i-IUP),min(i-1,n-1)
    ANONSYM((I-1)*IBAND+J-I+IUP+1) = AMAT(IDIAG+J-I,I)
  enddo
enddo
CALL NONSYM(ANONSYM,AMAT,Y2,N,IUP,IUP,1.0D-06,INFO2)
!%nonsym
!
!     ADAPT Y AT NODES
!
if (ZTAUEFF .NE. 0.D0) then
!
  do K=2,N-1
    YNEW(K) = Y(K) - ZTAUEFF* ((Y2(K+1)-Y2(K))*WORK(K) &
  - (Y2(K)-Y2(K-1))*WORK(K-1))
  enddo
  YNEW(1) = Y(1) - ZTAUEFF * (Y2(2)-Y2(1))*WORK(1)
  YNEW(N) = Y(N) + ZTAUEFF * (Y2(N)-Y2(N-1))*WORK(N-1)

endif

if (INFO2 .LT. 0) then
  PRINT *,' ERROR IN SPBTRS: INFO2 = ',INFO2
  STOP 'INFO2'
endif

return
end subroutine CUBSPLFIT

!----------------------------------------------------------------
real*8 function FC3(X1, F1, P1, X2, F2, P2)

implicit none

real*8, intent(in) :: X1, F1, P1, X2, F2, P2

FC3 = (2.D0* (F2 - F1) / (X1 - X2) + (P1 + P2)) / &
      ((X1 - X2) * (X1 - X2))

return
end function FC3

!----------------------------------------------------------------
real*8 function FC2(X1, F1, P1, X2, F2, P2)

implicit none

real*8, intent(in) :: X1, F1, P1, X2, F2, P2

FC2 = (3.D0* (X1 + X2) * (F1 - F2) / (X1 - X2) - &
      P1 * (X1 + 2.D0* X2) - P2 * (X2 + 2.D0* X1)) / &
      ((X1 - X2) * (X1 - X2))

return
end function FC2

!----------------------------------------------------------------
real*8 function FC1(X1, F1, P1, X2, F2, P2)

implicit none

real*8, intent(in) :: X1, F1, P1, X2, F2, P2

FC1 = (6.D0* X1 * X2 * (F2 - F1) / (X1 - X2) + &
      X2 * P1 * (2 * X1 + X2) + X1 * P2 * (X1 + 2.D0* X2)) / &
      ((X1 - X2) * (X1 - X2))

return
end function FC1

!----------------------------------------------------------------
real*8 function FC0(X1, F1, P1, X2, F2, P2)

implicit none

real*8, intent(in) :: X1, F1, P1, X2, F2, P2

FC0 = (F1 * X2**2 + F2 * X1**2 - &
      X1 * X2 * (X2 * P1 + X1 * P2) + 2.D0* X1 * X2 * &
      (F1 * X2 - F2 * X1) / (X1 - X2)) / ((X1 - X2) * (X1 - X2))

return
end function FC0

!----------------------------------------------------------------
real*8 function FCDCD0(X1, F1, P1, X2, F2, P2, PX)

implicit none

real*8, intent(in) :: X1, F1, P1, X2, F2, P2, PX
real*8, external :: FC0, FC1, FC2, FC3

FCDCD0 = FC0(X1, F1, P1, X2, F2, P2)   + &
   PX * (FC1(X1, F1, P1, X2, F2, P2)   + &
   PX * (FC2(X1, F1, P1, X2, F2, P2)   + &
   PX *  FC3(X1, F1, P1, X2, F2, P2)))

return
end function FCDCD0

!----------------------------------------------------------------
real*8 function FCDCD1(X1, F1, P1, X2, F2, P2, PX)

implicit none

real*8, intent(in) :: X1, F1, P1, X2, F2, P2, PX
real*8, external :: FC1, FC2, FC3

FCDCD1 =    FC1(X1,F1,P1,X2,F2,P2)  + &
   PX*(2.D0*FC2(X1,F1,P1,X2,F2,P2)  + &
   3.D0 *PX*FC3(X1,F1,P1,X2,F2,P2))

return
end function FCDCD1

!----------------------------------------------------------------
real*8 function FCDCD2(X1, F1, P1, X2, F2, P2, PX)

implicit none

real*8, intent(in) :: X1, F1, P1, X2, F2, P2, PX
real*8, external :: FC2, FC3

FCDCD2 = 2.D0* FC2(X1, F1, P1, X2, F2, P2)      + &
         6.D0* FC3(X1, F1, P1, X2, F2, P2) * PX

return
end function FCDCD2

!----------------------------------------------------------------
real*8 function FD2(X1, F1, P1, X2, F2)

implicit none

real*8, intent(in) :: X1, F1, P1, X2, F2

FD2 = ((F2 - F1)/(X2 - X1) - P1) / (X2 - X1)

return
end function FD2

!----------------------------------------------------------------
real*8 function FD1(X1, F1, P1, X2, F2)

implicit none

real*8, intent(in) :: X1, F1, P1, X2, F2
real*8, external :: FD2

FD1 = P1 - 2.*X1*FD2(X1, F1, P1, X2, F2)

return
end function FD1

!----------------------------------------------------------------
real*8 function FD0(X1, F1, P1, X2, F2)

implicit none

real*8, intent(in) :: X1, F1, P1, X2, F2
real*8, external :: FD1, FD2

FD0 = F1 - X1*(X1*FD2(X1,F1,P1,X2,F2) + FD1(X1,F1,P1,X2,F2))

return
end function FD0

!----------------------------------------------------------------
real*8 function FQDQ0(X1, F1, P1, X2, F2, PX)

implicit none

real*8, intent(in) :: X1, F1, P1, X2, F2, PX
real*8, external :: FD0, FD1, FD2

FQDQ0 = FD0(X1, F1, P1, X2, F2) + &
    PX*(FD1(X1, F1, P1, X2, F2) + &
    PX* FD2(X1, F1, P1, X2, F2))

return
end function FQDQ0

!----------------------------------------------------------------
real*8 function FQDQ1(X1, F1, P1, X2, F2, PX)

implicit none

real*8, intent(in) :: X1, F1, P1, X2, F2, PX
real*8, external :: FD1, FD2

FQDQ1 = FD1(X1, F1, P1, X2, F2) + &
  2.*PX*FD2(X1, F1, P1, X2, F2)

return
end function FQDQ1

!----------------------------------------------------------------
real*8 function FQDQ2(X1, F1, P1, X2, F2, PX)

implicit none

real*8, intent(in) :: X1, F1, P1, X2, F2, PX
real*8, external :: FD2

FQDQ2 = 2.D0* FD2(X1, F1, P1, X2, F2)

return
end function FQDQ2

!----------------------------------------------------------------
subroutine SPLINTEND(XA, YA, Y2A, N, X, Y, YP, YPP, KOPT)

implicit none

integer, intent(in) :: n, kopt
real*8, intent(in) :: X
real*8, intent(in), dimension(n) :: xa, ya, y2a
real*8, intent(out) :: Y, YP, YPP

integer :: k, khi, klo
real*8 :: h, a, b, zepsilon, ypakhi, ypaklo
real*8, external :: FCDCD0, FCDCD1, FCDCD2, FQDQ0, FQDQ1, FQDQ2

!----------------------------------------------------------------
!     KOPT = 0: NOTHING SPECIAL IF X OUT OF BOUND
!     KOPT = 1: USE QUADRATIC INTERPOLATION IF X OUT OF BOUND
!     KOPT = 2: USE CUBIC INTERPOLATION IF X OUT OF BOUND
!----------------------------------------------------------------
!*COMDECK CUCDCD
!-----------------------------------------------------------------------
! --     STATEMENT FUNCTION FOR CUBIC INTERPOLATION                   --
! --                         19.01.87            AR        CRPP       --
! --                                                                  --
! -- CUBIC INTERPOLATION OF A FUNCTION F(X)                           --
! -- THE SIX ARGUMENTS X1,F1,P1,X2,F2,P2 ARE DEFINED AS FOLLOWS:      --
! -- F(X1) = F1 , F(X2) = F2 , DF/DX(X1) = P1 , DF/DX(X2) = P2        --
! ----------------------------------------------------------------------
! -- FCDCD0 GIVES THE VALUE OF THE FUNCTION AT POINT PX               --
! -- FCDCD0(......,PX) = F(PX)                                        --
! ----------------------------------------------------------------------
! -- FCDCD1 GIVES THE VALUE OF THE FIRST DERIVATIVE OF F(X) AT PX:    --
! -- FCDCD1(......,PX) = DF/DX (PX)                                   --
! ----------------------------------------------------------------------
! -- FCDCD2 GIVES THE VALUE OF THE SECOND DERIVATIVE OF F(X) AT PX:   --
! -- FCDCD2(......,PX) = D2F/DX2 (PX)                                 --
! ----------------------------------------------------------------------
! -- FCDCD3 GIVES THE VALUE OF THE THIRD DERIVATIVE OF F(X) AT PX:    --
! -- FCDCD3(......,PX) = D3F/DX3 (PX)                                 --
! ----------------------------------------------------------------------
!         FCDCD3(X1,F1,P1,X2,F2,P2,PX) =
!     =                      6.D0* FC3(X1,F1,P1,X2,F2,P2)
!-----------------------------------------------------------------------
!*COMDECK QUAQDQ
! ----------------------------------------------------------------------
! --     STATEMENT FUNCTION FOR QUADRATIC INTERPOLATION               --
! --                         19.01.87            AR        CRPP       --
! --                                                                  --
! -- QUADRATIC INTERPOLATION OF A FUNCTION F(X)                       --
! -- THE FIVE PARAMETERS X1,F1,P1,X2,F2    ARE DEFINED AS FOLLOWS:    --
! -- F(X1) = F1 , DF/DX(X1) = P1 , F(X2) = F2                         --
! ----------------------------------------------------------------------
! -- FQDQ0 GIVES THE VALUE OF THE FUNCTION AT POINT PX                --
! -- FQDQ0(......,PX) = F(PX)                                         --
! ----------------------------------------------------------------------
! -- FQDQ1 GIVES THE VALUE OF THE FIRST DERIVATIVE OF F(X) AT PX      --
! -- FQDQ1(......,PX) = DF/DX (PX)                                    --
! ----------------------------------------------------------------------
! -- FQDQ2 GIVES THE VALUE OF THE SECOND DERIVATIVE OF F(X) AT PX     --
! -- FQDQ2(......,PX) = D2F/DX2 (PX)                                  --
! ----------------------------------------------------------------------

KLO=1
KHI=N
1     if (KHI-KLO.GT.1) then
  K=(KHI+KLO)/2
  if(XA(K).GT.X)then
    KHI=K
  else
    KLO=K
  endif
GOTO 1
endif
H=XA(KHI)-XA(KLO)
if (H.EQ.0.) STOP 'BAD XA INPUT.'
A=(XA(KHI)-X)/H
B=(X-XA(KLO))/H
Y=A*YA(KLO)+B*YA(KHI)+ &
           ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.
YP=(YA(KHI)-YA(KLO))/H - &
     ( (3.*A*A-1.)*Y2A(KLO) - (3.*B*B-1.)*Y2A(KHI) )*H/6.
YPP=A*Y2A(KLO)+B*Y2A(KHI)

if (KOPT .EQ. 0) return

! ZEPSILON: DISTANCE OUTSIDE INTERVAL RELATIVE TO H

ZEPSILON = 1.0D-05

! x outside interval: use quadratic with y and yp at edge and y-1

if (KOPT .EQ. 1) then
!     LEFT END
  if (B.lt. -ZEPSILON) then
    print *,' warning points outside interval at left.', &
  ' Use quadratic interpolation'
    A=1.D0
    B=0.D0
    YPAKLO=(YA(KHI)-YA(KLO))/H - &
 ( (3.D0*A*A-1.D0)*Y2A(KLO) - (3.D0*B*B-1.D0)*Y2A(KHI) )*H/6.D0
    A=0.D0
    B=1.D0
    YPAKHI=(YA(KHI)-YA(KLO))/H - &
 ( (3.D0*A*A-1.D0)*Y2A(KLO) - (3.D0*B*B-1.D0)*Y2A(KHI) )*H/6.D0
    Y=FQDQ0(XA(KLO),YA(KLO),YPAKLO,XA(KHI),YA(KHI),X)
    YP=FQDQ1(XA(KLO),YA(KLO),YPAKLO,XA(KHI),YA(KHI),X)
    YPP=FQDQ2(XA(KLO),YA(KLO),YPAKLO,XA(KHI),YA(KHI))
  endif
!     RIGHT END
  if (A .lt. -ZEPSILON) then
    print *,' warning points outside interval at right.', &
  ' Use quadratic interpolation'
    A=1.D0
    B=0.D0
    YPAKLO=(YA(KHI)-YA(KLO))/H - &
 ( (3.D0*A*A-1.D0)*Y2A(KLO) - (3.D0*B*B-1.D0)*Y2A(KHI) )*H/6.D0
    A=0.D0
    B=1.D0
    YPAKHI=(YA(KHI)-YA(KLO))/H - &
 ( (3.D0*A*A-1.D0)*Y2A(KLO) - (3.D0*B*B-1.D0)*Y2A(KHI) )*H/6.D0
    Y=FQDQ0(XA(KHI),YA(KHI),YPAKHI,XA(KLO),YA(KLO),X)
    YP=FQDQ1(XA(KHI),YA(KHI),YPAKHI,XA(KLO),YA(KLO),X)
    YPP=FQDQ2(XA(KHI),YA(KHI),YPAKHI,XA(KLO),YA(KLO))
  endif

  return
endif

! x outside interval: use CUBIC with y and yp

if (KOPT .EQ. 2) then
! LEFT OR RIGHT END
  if (a.lt.-ZEPSILON .OR. B.LT.-ZEPSILON) then
    print *,' warning points outside interval.', &
  ' Use cubic interpolation'
    A=1.D0
    B=0.D0
    YPAKLO=(YA(KHI)-YA(KLO))/H - &
 ( (3.D0*A*A-1.D0)*Y2A(KLO) - (3.D0*B*B-1.D0)*Y2A(KHI) )*H/6.D0
    A=0.D0
    B=1.D0
    YPAKHI=(YA(KHI)-YA(KLO))/H - &
 ( (3.D0*A*A-1.D0)*Y2A(KLO) - (3.D0*B*B-1.D0)*Y2A(KHI) )*H/6.D0
    Y=FCDCD0(XA(KLO),YA(KLO),YPAKLO,XA(KHI),YA(KHI),YPAKHI,X)
    YP=FCDCD1(XA(KLO),YA(KLO),YPAKLO,XA(KHI),YA(KHI),YPAKHI,X)
    YPP=FCDCD2(XA(KLO),YA(KLO),YPAKLO,XA(KHI),YA(KHI),YPAKHI,X)
  endif
  return
endif

return
end subroutine SPLINTEND

!--------------------------------------------------------
subroutine NONSYM(A, X, C, N, MLEFT, MRIGHT, EPS, NCOND)
!
!     5.3  SOLVES A real VALUED NONSYMMETRIC LINEAR SYSTEM  A . X  =  C
!
!     VERSION 1          APRIL 1988     KA      LAUSANNE
!     VERSION 2          MAI 1992       AJ      LAUSANNE
!---------------------------------------------------------

implicit none

integer, intent(in) :: n, mleft, mright
real*8, intent(in) :: eps
real*8, intent(out), dimension(n) :: x, c
integer, intent(out) :: ncond
real*8, intent(out), dimension(N*(MLEFT+MRIGHT+1)) :: a

integer :: j, imess, icorn, iconst, idown, ielem, ihor1, ihor2, &
   iknown, incmpl, inm1, ipiv1, ipiv2, isolut, &
   ja, jdown, jhoriz, jp, jpabs, jpivot, mband, moffdi
real*8 :: ratio, ZPIVOT, ZSUM, ZTOP, zmax

data IMESS / 0 /

!-----------------------------------------------------------------------
!
!     A IS A BANDMATRIX OF RANK N WITH MLEFT/MRIGHT OFF-DIAGONAL
!     ELEMENTS TO THE LEFT/RIGHT.
!     :A: AND :C: ARE DESTROYED BY :NONCYM: 
!     AT THE return :C: CONTAINS THE SOLUTION VECTOR, :X: THE SQUARE
!     ROOT OF THE DIAGONAL ELEMENTS
!
!     THE ELEMENTS ARE HORIZONTALLY NUMBERED ROW AFTER ROW
!
!     FOR EASY NUMBERING, ZERO ELEMENTS HAVE BEEN INTRODUCED IN THE
!     UPPER LEFT-HAND CORNER AND IN THE LOWER RIGHT-HAND CORNER.
!     THESE ELEMENTS
!     M U S T  B E  S E T  T O  Z E R O
!     BEFORE CALLING :NONCYM: 
!
!     EXAMPLE FOR NUMBERING   MLEFT=2,  MRIGHT=3
!     -------
!
!     A(1)   A(2) I A(3)   A(4)   A(5)   A(6)
!     A(7) I A(8)   A(9)   A(10)  A(11)  A(12)
!     I                 .
!     I                        .
!
!     A PIVOT IS CONSIDERED TO BE BAD IF IT IS BY A FACTOR :EPS: SMALLER
!     THAN ITS OFF-DIAGONAL ELEMENTS. IF A BAD PIVOT IS ENCOUNTERED
!     A MESSAGE IS ISSUED (AT MOST TEN TIMES IN A RUN).
!     THE FLAG :NCOND: IS SET TO -1, OTHERWISE 0. 
!
!-----------------------------------------------------------------------
!L    0.        INITIALIZATION

MOFFDI=MLEFT+MRIGHT
MBAND=MOFFDI+1
RATIO=0.D0
NCOND=0
! PRELIMINARY CHECK OF PIVOTS
IPIV1=MLEFT+1
IPIV2=IPIV1+(N-1)*MBAND
do 10 J=IPIV1,IPIV2,MBAND
  if( ABS(A(J)) .EQ. 0.D0) GO TO 510
 10   continue

!-----------------------------------------------------------------------
!L    1.        PRECONDITIONNING: GET ONES ON THE DIAGONAL

do 130 JP=1,N
  JPABS=(MLEFT+1)+(JP-1)*MBAND
  X(JP)=SQRT(ABS(A(JPABS)))

! DIVIDE THE EQUATIONS BY THE SQUARE ROOT OF THE PIVOT
  do 110 J=-MLEFT,MRIGHT
    A(JPABS+J)=A(JPABS+J)/X(JP)
 110    continue
  C(JP)=C(JP)/X(JP)

! CHANGE VARIABLES, DIVIDE COLUMN BY SQUARE ROOT OF THE PIVOT
  do 120 JA=JPABS-MRIGHT*MOFFDI,JPABS+MLEFT*MOFFDI,MOFFDI
    if (JA.GE.IPIV1 .AND. JA.LE.IPIV2) A(JA)=A(JA)/X(JP)
 120    continue
 130  continue

!-----------------------------------------------------------------------
!L    2.        CONSTRUCT UPPER TRIANGULAR MATRIX 

! INITIALIZATION OF VERTICAL COUNTER
IDOWN=MLEFT
! PIVOT INDEX FOR FIRST INCOMPLETE COLUMN
INCMPL=(N-MLEFT)*MBAND+MLEFT+1 
! FIRST PIVOT AND LAST BUT ONE
IPIV1=MLEFT+1
IPIV2=IPIV1+(N-2)*MBAND

do 220 JPIVOT=IPIV1,IPIV2,MBAND
! COLUMN LENGTH 
  if(JPIVOT.GE.INCMPL) IDOWN=IDOWN-1
  ZPIVOT=A(JPIVOT)
! FIRST AND LAST HORIZONTAL ELEMENT INDEX
  IHOR1=JPIVOT+1
  IHOR2=JPIVOT+MRIGHT
! CHECK THE PIVOT
  ZMAX=0.D0
  do 205 JHORIZ=IHOR1,IHOR2
    ZMAX=DMAX1(ZMAX, ABS(A(JHORIZ)))
 205    continue
  if( ABS(ZPIVOT).EQ.0.D0) GO TO 510
  RATIO=DMAX1(RATIO,ZMAX/ ABS(ZPIVOT))

  do 210 JHORIZ=IHOR1,IHOR2
    ZTOP=-A(JHORIZ)/ZPIVOT
    A(JHORIZ)=ZTOP
! INITIALIZATION OF ELEMENT INDEX
    IELEM=JHORIZ
! INITIALIZATION OF RECTANGULAR RULE CORNER ELEMENT
    ICORN=JPIVOT
! LOOP DOWN THE COLUMN
    do 211 JDOWN=1,IDOWN 
      IELEM=IELEM+MOFFDI
      ICORN=ICORN+MOFFDI
      A(IELEM)=A(IELEM)+ZTOP*A(ICORN)
 211      continue
 210    continue
! TREAT THE CONSTANTS
! INDEX OF THE FIRST ONE
  ICONST=JPIVOT/MBAND+1
  ZTOP=-C(ICONST)/ZPIVOT
  C(ICONST)=ZTOP
! INITIALIZATION OF RECTANGULAR RULE CORNER ELEMENT
  ICORN=JPIVOT
! LOOP DOWN THE CONSTANTS 
  do 221 JDOWN=1,IDOWN 
    ICONST=ICONST+1
    ICORN=ICORN+MOFFDI
    C(ICONST)=C(ICONST)+ZTOP*A(ICORN)
 221    continue
 220  continue

!-----------------------------------------------------------------------
!L    3.        BACKSUBSTITUTION

 300  continue
! INITIALIZATION OF SOLUTION INDEX
ISOLUT=N
! LAST PIVOT
JPIVOT=IPIV2+MBAND
! CHECK LAST PIVOT
if( ABS(A(JPIVOT)).EQ.0.D0) go to 510
! LAST UNKNOWN
C(ISOLUT)=C(ISOLUT)/A(JPIVOT)
! INITIALIZATION OF HORIZONTAL RANGE
IHOR1=JPIVOT+1
IHOR2=JPIVOT+MRIGHT

INM1=N-1
do 320 J=1,INM1
  ISOLUT=ISOLUT-1
  IHOR1=IHOR1-MBAND
  IHOR2=IHOR2-MBAND
  ZSUM=-C(ISOLUT)
  IKNOWN=ISOLUT
  do 310 JHORIZ=IHOR1,IHOR2
    IKNOWN=IKNOWN+1
    if(IKNOWN.GT.N) GO TO 310
    ZSUM=ZSUM+A(JHORIZ)*C(IKNOWN)
 310    continue
  C(ISOLUT)=ZSUM
 320  continue

!-----------------------------------------------------------------------
!L    4.        PRECONDITIONNING: BACK TO ORIGINAL VARIABLES

do 410 J=1,N
  C(J)=C(J)/X(J)
 410  continue

!-----------------------------------------------------------------------
!L    5.        HAVE WE ENCOUNTERED A BAD PIVOT ? 

if(RATIO .LT. 1.D0/EPS) return
NCOND=-1
! COUNT THE NUMBER OF MESSAGES ISSUED
IMESS=IMESS+1
if(IMESS.LE.10) PRINT 9500, RATIO
return

! ZERO PIVOT ENCOUNTERED
 510  continue
PRINT 9510 
STOP 'ZERO PIVOT IN NONCYM'

 9500 FORMAT(/////20X,10(1H*),21H  MESSAGE FROM NONCYM,2X,10(1H*)/// 20X, &
       'BAD PIVOT ENCOUNTERED, RATIO: ',1PE12.3/// 20X, 43(1H*)///)
 9510 FORMAT(/////20X,10(1H*),'   ZERO PIVOT IN NONCYM'///)

end subroutine NONSYM
