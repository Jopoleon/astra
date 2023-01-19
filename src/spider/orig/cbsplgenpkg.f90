!
! This package contains the routines required to use cubic spline with tension
! It contains routines:
! cbsplgen0.f90
! cbsplgen.f90
! cbfitbnd.f90
! splibnd.f90
! It also calls Lapack routines
!
! It is based on paper HIRSHMAN ET AL, PHYS. PLASMAS 1 (1994) 2280.
! Routines created by olivier.sauter@epfl.ch
!
SUBROUTINE cbsplgen0(XIN,YIN,nin,xout,yout,youtp,youtpp, &
     &  nout,ioptder,iextrapo,psig,nbc,ybc,IFLAG)
  !
  USE prec_rkind
  implicit none
  integer LENREAL
  parameter(LENREAL=8)
  integer nin,nout, nbc(2)
  REAL(RKIND) ::  psig(nin)
  REAL(RKIND) ::  xin(nin), yin(nin), ybc(6)
  REAL(RKIND) ::  xout(nout),yout(nout), youtp(nout), youtpp(nout)
  INTEGER ioptder, NBCLFT,NBCRGT, iextrapo
  REAL(RKIND) :: XBCLFT,XBCRGT, YBCLFT,YBCRGT
  REAL(RKIND) :: PXEXP0,PXEXPDL
  !      pointer(iptr_pynew,pynew)
  REAL(RKIND), DIMENSION(:), ALLOCATABLE :: pynew
  !      REAL(RKIND) :: pynew(1)
  !      pointer(iptr_pyinpp,pyinpp)
  REAL(RKIND), DIMENSION(:), ALLOCATABLE :: pyinpp
  !      REAL(RKIND) :: pyinpp(1)
  !      pointer(iptr_pamat,pamat)
  REAL(RKIND), DIMENSION(:), ALLOCATABLE :: pamat
  !      REAL(RKIND) :: pamat(1)
  !      pointer(iptr_pwork,pwork)
  REAL(RKIND), DIMENSION(:), ALLOCATABLE :: pwork
  !      REAL(RKIND) :: pwork(1)
  !      pointer(iptr_kwork,kwork)
  INTEGER, DIMENSION(:), ALLOCATABLE :: kwork
  !      integer kwork(1)
  !
  integer iflag, idamat, mdamat, nbytes, infomalloc
  !OS      integer*4 malloc_f
  !
  !
  !%OS      NBCLFT = 1
  !%OS      NBCRGT = 1
  !%OS      YBCLFT = 1.E32_RKIND
  !%OS      YBCRGT = 1.E32_RKIND
  NBCLFT = nbc(1)
  NBCRGT = nbc(2)
  YBCLFT = ybc(1)
  YBCRGT = ybc(2)
  if (NBCLFT .ge. 10) XBCLFT = ybc(3)
  if (NBCRGT .ge. 10) XBCRGT = ybc(4)
  !
  IDAMAT = 3
  IF (PSIG(1) .EQ. 0._RKIND) IDAMAT = 2
  IF (NBCLFT.GE.10 .OR. NBCRGT.GE.10 .OR. NBCLFT.EQ.2 &
       &  .OR. NBCRGT.EQ.2) IDAMAT = 3*(IDAMAT-1)+1
  IF (NBCLFT .GE. 10) IDAMAT = IDAMAT + 1
  IF (NBCRGT .GE. 10) IDAMAT = IDAMAT + 2
  mdamat = IDAMAT*nin
  IF( .NOT. ALLOCATED(pynew) ) ALLOCATE(pynew(nin))
  IF( .NOT. ALLOCATED(pyinpp) ) ALLOCATE(pyinpp(nin))
  IF( .NOT. ALLOCATED(pamat) ) ALLOCATE(pamat(mdamat))
  IF( .NOT. ALLOCATED(pwork) ) ALLOCATE(pwork(2*nin))
  IF( .NOT. ALLOCATED(kwork) ) ALLOCATE(kwork(5*nin))
  PXEXP0 = ybc(5)
  PXEXPDL= ybc(6)
  CALL CBSPLGEN(XIN,YIN,PYNEW,PYINPP,Nin,XOUT,YOUT,YOUTP,YOUTPP, &
       &    nout,ioptder,PSIG,KWORK,PWORK,PAMAT,mdamat,NBCLFT, &
       &    NBCRGT,XBCLFT,XBCRGT,YBCLFT,YBCRGT,iextrapo,PXEXP0,PXEXPDL &
       &    ,IFLAG)
  !
END SUBROUTINE cbsplgen0
!
!.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.     
!
SUBROUTINE CBSPLGEN(PXIN,PYIN,PYINNEW,PYINPP,KNIN,PXOUT,PYOUT, &
     &  PYOUTP,PYOUTPP,KNOUT,KOPT,PSIG,KWORK,PWORK,PAMAT,MDMATOT,NBCLFT &
     &  ,NBCRGT,XBCLFT,XBCRGT,YBCLFT,YBCRGT,KEXTRPOL,PXEXP0,PXEXPDL &
     &  ,KFLAG)
  !     =================================================================
  !
  !     NOTE: THIS ROUTINE INCLUDES THE STANDARD CUBIC SPLINE IF PTAUS=0 (i.e. psig=0):
  !           THEN PYINNEW IS NOT USED AND MDAMAT=3 IS SUFFICIENT
  !           (=> PYINNEW(1) OR PYINNEW=PYIN IS OK)
  !
  !     Interpolate (pxin,pyin) on (pxout,pyout) using
  !     Hirshman fitted cubic spline with ptaus value or
  !     standard cubic spline if PTAUS=0 (psig=ptaus*sigma/sigmamin)
  !
  !     KOPT = 0: ONLY INTERPOLATE FUNCTION INTO PYOUT
  !     KOPT = 1: INTERPOLATE FUNCTION INTO PYOUT AND 1ST DER. INTO PYOUTP
  !     KOPT = 2: AS KOPT=1 PLUS 2ND DER. INTO PYOUTPP
  !
  !     MDMATOT = TOTAL DIMENSION OF PAMAT. THE REQUIRED SPACE DEPENDS IF
  !     .         PTAUS IS ZERO OR NOT AND ON THE B.C. (SYMMETRIC OR NOT)
  !     THUS, MDMATOT CAN VARY BETWEEN 2*KNIN AND 10*KNIN (MAX. VALUE NEEDED)
  !
  !     SEE COMMENTS FOR ROUTINE CBFITBND FOR MORE INFORMATION
  !
  !     IF LAPACK ROUTINES NOT AVAILABLE, USE spgbtrf_s.f
  !     (LAPACK SOURCE COPIED FROM NETLIB.ORG)
  !
  !     PXIN    : INPUT ABSCISSA (GIVEN DATA)
  !     PYIN    : INPUT VALUES OF FUNCTION AT PXIN(I),I=1,KNIN
  !     PYINNEW : IF PTAUS.NE.0, THEN PYNEW CONTAINS ON OUTPUT THE NEW VALUES
  !     .         OF THE FUNCTION AT PXIN(I) FOR THE CUBIC SPLINE FIT
  !     PYINPP  : SECOND DER. OF THE CUBIC SPLINE FIT FOR (PXIN,PYIN) IF PTAUS=0 OR
  !     .         ON (PXIN,PYNEW) OTHERWISE
  !     KNIN    : NUMBER OF INPUT POINTS
  !     PXOUT   : X VALUES AT WHICH THE FUNCTION HAS TO BE INTERPOLATED (INPUT)
  !     PYOUT   : INTERPOLATED VALUES AT PXOUT(I),I=1,KNOUT (OUTPUT)
  !     PYOUTP  : INTERPOLATED VALUES OF 1ST DER. OF FUNCTIONS AT PXOUT(I) (OUTPUT, IF KOPT.GE.1)
  !     PYOUTPP : INTERPOLATED VALUES OF 2ND DER. OF FUNCTIONS AT PXOUT(I) (OUTPUT, IF KOPT.EQ.2)
  !     KNOUT   : NUMBER OF POINTS FOR OUTPUT
  !     KOPT    : SEE ABOVE
  !     PSIG    : SIGMA at each point, normalized by minimum sigma, times PTAUS
  !     PTAUS   : WEIGHT OF SECOND DERIVATIVE IN THE CHI**2 TO BE MINIMIZED. PTAUS=0 GIVES THE
  !     .         STANDARD CUBIC SPLINE. LARGER VALUES OF PTAUS WILL SMOOTH MORE THE 2ND DER.
  !     KWORK   : INTEGER WORK SPACE OF DIMENSION 5*KNIN
  !     PWORK   : REAL WORK SPACE OF DIMENSION 2*KNIN
  !     PAMAT   : REAL WORK SPACE FOR THE MATRIX OF DIMENSION MDMATOT
  !     MDMATOT : DIMENSION OF PAMAT, WHICH VALUE SHOULD BE IN [2*KNIN,10*KNIN], DEPENDING ON
  !     .         PTAUS BEING 0 OR NOT AND ON THE BOUNDARY CONDITION (SEE PARAGRAPH 1. BELOW)
  !     NBCLFT  : FOR LEFT B.C., VALUE SHOULD BE 0,1,2,10,11 OR 12. (SEE ROUTINE CBFITBND BELOW)
  !     NBCRGT  : FOR RIGHT B.C. (SEE ROUTINE CBFITBND BELOW)
  !     XBCLFT  : FOR LEFT B.C., USED ONLY IF NBCLFT.GE.10 (SEE ROUTINE CBFITBND BELOW)
  !     XBCRGT  : FOR RIGHT B.C., USED ONLY IF NBCRGT.GE.10 (SEE ROUTINE CBFITBND BELOW)
  !     YBCLFT  : VALUE OF LEFT B.C.
  !     YBCRGT  : VALUE OF RIGHT B.C.
  !     .         STANDARD B.C. (SECOND DER. = 0) IS OBTAINED WITH:
  !     .         NBCLFT = NBCRGT = 0 AND YBCLFT = YBCRGT = 0.
  !     KEXTRPOL: OPTION ON HOW TO EXTRAPOLATE THE FUNCTION IF PXOUT(I) IS OUTSIDE [PXIN(1),PXIN(KNIN)]
  !     .       = 0: STOP WITH ERROR MESSAGE IF OUT OF BOUND
  !     .       = 1: LINEAR EXTRAPOLATION
  !     .       = 2: USE QUADRATIC INTERPOLATION IF X OUT OF BOUND
  !     .       = 3: USE CUBIC INTERPOLATION IF X OUT OF BOUND
  !     .       = 21: USE QUADRATIC WITHIN ALFA*DELTA_X AND LINEAR FURTHER
  !     .       = 31: USE CUBIC WITHIN ALFA*DELTA_X AND LINEAR    FURTHER
  !     .       = 32: USE CUBIC WITHIN ALFA*DELTA_X AND QUADRATIC FURTHER
  !     PXEXP0  : PTAUS IS WEIGHTED BY AN EXP(-((X-PXEXP0)/PXEXPDL)**2)
  !     PXEXPDL : IF PXEXP0 NOT IN [PXIN(1),PXIN(KNIN)], EXP() IGNORED AND PTAUS=CST
  !     .         (SEE ROUTINE CBFITBND BELOW)
  !     .         GIVE PXEXP0=PXIN(1)-1._RKIND AND PXEXPDL=1. TO GET CONSTANT PTAUS
  !     KFLAG   : ERROR FLAG: IF NOT 0, THERE IS A PROBLEM
  !
  !-----------------------------------------------------------------------
  !c      implicit REAL(RKIND) :: (p)
  USE prec_rkind
  implicit none
  REAL(RKIND) :: EPSILON
  PARAMETER(EPSILON = 1.0E-10_RKIND)
  !
  integer knin, knout, kopt, mdmatot, nbclft, nbcrgt, kextrpol, &
       &  kflag, i, idamat
  REAL(RKIND) :: PXIN(KNIN), PYIN(KNIN), PYINNEW(KNIN), PYINPP(KNIN), &
       &  PXOUT(KNOUT), PYOUT(KNOUT), PYOUTP(KNOUT), PYOUTPP(KNOUT), &
       &  PSIG(KNIN)
  !   
  REAL(RKIND) :: PAMAT(MDMATOT), PWORK(2*KNIN)
  INTEGER :: KWORK(5*KNIN)
  !
  REAL(RKIND) :: xbclft, xbcrgt, ybclft, ybcrgt, pxexp0, pxexpdl, zy, zyp, &
       &  zypp
  !
  !-----------------------------------------------------------------------
  !     0. CHECK INPUT CONSISTENCY
  !
  KFLAG = 0
  IF (PSIG(1) .EQ. 0._RKIND) THEN
     IF (NBCLFT.GE.10 .OR. NBCRGT.GE.10 .OR. NBCLFT.EQ.2 &
          &  .OR. NBCRGT.EQ.2) THEN
        PRINT *,' PTAUS=0, BUT NEED SMOOTHING, WHEN'
        PRINT *,'     NBCLFT.GE.10 .OR. NBCRGT.GE.10 .OR.', &
             &      ' NBCLFT.EQ.2 .OR. NBCRGT.EQ.2'
        PRINT *,' NBCLFT = ',NBCLFT
        PRINT *,' NBCRGT = ',NBCRGT
        !%OS          STOP 'TAU=0'
        KFLAG = 1
        return
     ENDIF
  ENDIF
  !
  !   PXIN in ASCENDING ORDER
  !
  DO i=1,KNIN-1
     if (PXIN(i) .GE. PXIN(i+1)) then
        print *,' xin not in ascending order:'
        print *,' xin(',i,')= ',PXIN(i),'   >=   xin(',i+1,')= ', &
             &      PXIN(i+1)
        KFLAG = 2
        RETURN
     endif
  END DO
  !
  !     1. PREPARE MATRIX DIMENSION.
  !     MINIMUM REQUIRED:
  !     IUP = 1 = IDOWN; IF PTAUS.NE.0 => IUP = IDOWN = 2
  !     IF SYMMETRIC, USE ONLY UPPER BAND AND DIAGONAL =>IDAMAT=IUP+1
  !     IF ASYMMETRIC => IDAMAT = 2*IDOWN + IUP + 1
  !     IF B.C. NOT AT END OF INTERVAL => IUP = IUP + 1, AND/OR IDOWN=IDOWN+1
  !
  !     => ALTOGETHER, MINIMUM VALUE: IDOWN=1, IUP=1, SYM. =>IDAMAT_MAX = 2
  !     => ALTOGETHER, MAXIMUM VALUE: IDOWN=3, IUP=3 =>IDAMAT_MAX = 10
  !
  IDAMAT = 3
  IF (PSIG(1) .EQ. 0._RKIND) IDAMAT = 2
  IF (NBCLFT.GE.10 .OR. NBCRGT.GE.10 .OR. NBCLFT.EQ.2 &
       &  .OR. NBCRGT.EQ.2) IDAMAT = 3*(IDAMAT-1)+1
  IF (NBCLFT .GE. 10) IDAMAT = IDAMAT + 1
  IF (NBCRGT .GE. 10) IDAMAT = IDAMAT + 2
  !
  IF (MDMATOT .LT. IDAMAT*KNIN) THEN
     PRINT *,' '
     PRINT *,' DIMENSION MDMATOT= ',MDMATOT,' IS TOO SMALL,', &
          &    ' NEED IDAMAT*KNIN= ',IDAMAT,'*',KNIN,' = ',IDAMAT*KNIN
     !%OS        STOP
     RETURN
  ENDIF
  !
  CALL CBFITBND(PXIN,PYIN,PYINNEW,KNIN,PYINPP,PSIG,KWORK, &
       &  PWORK(1),PWORK(KNIN+1),PAMAT,IDAMAT,NBCLFT,NBCRGT, &
       &  XBCLFT,XBCRGT,YBCLFT,YBCRGT,PXEXP0,PXEXPDL)
  !
  !L    2. COMPUTE INTERPOLATED VALUE AT EACH PXOUT
  !
  DO I=1,KNOUT
     IF (PSIG(1) .EQ. 0.0_RKIND) THEN
        CALL SPLIBND(PXIN,PYIN   ,PYINPP,KNIN,PXOUT(I),ZY,ZYP,ZYPP, &
             &      KEXTRPOL)
     ELSE
        CALL SPLIBND(PXIN,PYINNEW,PYINPP,KNIN,PXOUT(I),ZY,ZYP,ZYPP, &
             &      KEXTRPOL)
     ENDIF
     PYOUT(I) = ZY
     IF (KOPT .GE. 1) PYOUTP(I) = ZYP
     IF (KOPT .EQ. 2) PYOUTPP(I) = ZYPP
  END DO
  !
  RETURN
END SUBROUTINE CBSPLGEN
!-----------------------------------------------------------------------
SUBROUTINE CBFITBND(PXIN,PYIN,PYINNEW,KNIN,PYINPP,PSIG, &
     &  KPM2,WHK,WOHK,PAMAT,MDAMAT,NBCLFT,NBCRGT,XBCLFT,XBCRGT,YBCLFT, &
     &  YBCRGT,PXEXP0,PXEXPDL)
  !
  !     NON-PERIODIC B.C
  !
  !     PREPARE SECOND DERIVATIVE OF CUBIC SPLINE INTERPOLATION AND NEW
  !     VALUES OF Y AT NODES YNEW FITTED SUCH THAT CHI**2 + TAUS*F''**2
  !     IS MINIMIZED ACCORDING TO HIRSHMAN ET AL, PHYS. PLASMAS 1 (1994) 2280.
  !     SIG = TAU*SIGMA_K/min(SIGMA_K) OF PAPER.
  !
  !     SETTING TAUS=0., ONE FINDS THE USUAL CUBIC SPLINE INT. WITH CHI**2=0
  !     TAUS LARGE => FIT CLOSER TO STRAIGHT LINE (SECOND DERIV.=0)
  !
  !     IF LAPACK ROUTINES NOT AVAILABLE, USE spgbtrf_s.f
  !     (LAPACK SOURCE COPIED FROM NETLIB.ORG)
  !
  !     IF TAUS=0, PYINNEW NOT USED => PYINNEW(1) OR PYINNEW=PYIN IS OK
  !
  !     BOUNDARY CONDITIONS, 3 TYPES DETERMINED BY THE VALUE OF (NBCLFT/RGT):
  !
  !     0) VALUE OF SECOND DERIVATIVE AT XBCLFT OR RGT IS GIVEN (0 OR 10)
  !     1) VALUE OF 1ST        "       "   "     "  "   "   "   (1 OR 11)
  !     2) VALUE OF FUNCTION AT XBCLFT OR RGT IS GIVEN          (2 OR 12)
  !
  !     THE VALUE IS GIVEN BY YBCLFT OR YBCRGT RESPECTIVELY.
  !
  !     FOR TYPE 1: IF (YBCLFT OR YBCRGT > 1E31 THEN DER. FROM LAGRANGIAN INTERP.
  !     FOR TYPE 1: IF (YBCLFT OR YBCRGT <-1E31 THEN DER. FROM LINEAR     INTERP.
  !
  !     IF NBCLFT OR NBCRGT IS < 10, PXIN(1) OR PXIN(KNIN) IS USED INSTEAD
  !     OF XBCLFT OR XBCRGT, RESPECTIVELY => XBCLFT OR XBCRGT NOT USED
  !
  !     IF END POINTS ARE USED FOR THE B.C. AND TYPE 0 OR 1, THEN USE SYMMETRY
  !     OF MATRIX
  !
  !     IF XBCLFT OR XBCRGT ARE USED, IMPOSE B.C. ON NODE CLOSEST TO XBCLFT OR XBCRGT
  !
  !     TENSION TAUS(K) IS GIVEN WITH AN EXPONENTIAL FORM TO BE ABLE TO LOCALIZE
  !     IT:
  !     .     TAU_K = PTAUS * EXP( -COF * ((X-X0)/DX)**2)
  !
  !     WHERE X0 = PXEXP0 AND DX = PXEXPDL, AND:
  !     COF = 1. IF PXEXP0 IN [PXIN(1),PXIN(KNIN)], 0 OTHERWISE
  !     THUS SETTING PXEXP0 OUTSIDE DOMAIN GIVES A CST TAU_K VALUE
  !
  !-----------------------------------------------------------------------
  !
  USE prec_rkind
  implicit none
  integer LENREAL
  parameter(LENREAL=8)
  integer knin, mdamat, nbclft, nbcrgt, iik, &
       &  itauval, isym, n, i, j, k, iup, idown, idiag, ishift, ieff, ikp2 &
       &  ,ikp1, ikm1, ikm2, jk, jkp1, jkp2, jeff, iii, iupsofar, idwnsofa &
       &  , jbc, ik, idiamik, idiapik, iklft, idima, idimrhs, irhs, info, &
       &  info2, jkm1, jkm2, infomalloc
  REAL(RKIND) :: PXIN(KNIN), PYIN(KNIN), PYINNEW(KNIN), &
       &  PYINPP(KNIN), WHK(KNIN), WOHK(KNIN), &
       &  PAMAT(MDAMAT,KNIN), ZYBC(2), PSIG(KNIN)
  INTEGER :: KPM2(KNIN,-2:+2), IXBC(2), IBCTYP(2)
  !OS      pointer(iptr_ftauk,ftauk)
  REAL(RKIND), DIMENSION(KNIN) :: ftauk
  !OS      dimension ftauk(1)
  !
  REAL(RKIND) :: xbclft, xbcrgt, ybclft, ybcrgt, pxexp0, pxexpdl
  REAL(RKIND) :: xtkm1, xohkm1, xohkm2, xtk, xhkm1, &
       &  xohk, xtkp1, xhk, xohkp1, xykp1, xyk, &
       &  xykm1, ztaueff, zcofexp, zxexp0, zxexpdl, a1, a2, a3, a4, b1, &
       &  b2, b3, b4, px, fakk, fakkp1, fakkp2, frhs, zdelx, zero, zvalue, &
       &  zypeff, ztohkk1, zsign, fa2, fa3, fa0, fa1, &
       &  fakkm1, fakkm2, fun_ftauk, fcccc0, fcccc1, fcccc2, fcccc3
  !      integer*4 malloc_f
  !
  !
  !     FUNCTIONS FOR MATRIX COEFFICIENTS
  !
  REAL(RKIND) :: zsix, zthree, ztwo, zone
  PARAMETER(zsix=6._RKIND, zthree=3._RKIND, ztwo=2._RKIND, zone=1._RKIND)
  FAKKM2(XTKM1,XOHKM1,XOHKM2) = XTKM1*XOHKM1*XOHKM2
  FAKKM1(XTK,XTKM1,XHKM1,XOHK,XOHKM1,XOHKM2) = XHKM1/zsix &
       &  - XOHKM1*(XTK*XOHK+(XTK+XTKM1)*XOHKM1 + XTKM1*XOHKM2)
  FAKK(XTKP1,XTK,XTKM1,XHK,XHKM1,XOHK,XOHKM1) = (XHK+XHKM1)/ZTHREE &
       &  + XOHK*XOHK*(XTKP1+XTK) &
       &  + XOHKM1*(ZTWO*XTK*XOHK+(XTK+XTKM1)*XOHKM1)
  FAKKP1(XTKP1,XTK,XHK,XOHKP1,XOHK,XOHKM1) = XHK/zsix &
       &  - XOHK*(XTKP1*XOHKP1+(XTK+XTKP1)*XOHK + XTK*XOHKM1)
  FAKKP2(XTKP1,XOHKP1,XOHK) = XTKP1*XOHKP1*XOHK
  !
  FRHS(XYKP1,XYK,XYKM1,XOHK,XOHKM1) = (XYKP1-XYK)*XOHK &
       &  - (XYK-XYKM1)*XOHKM1
  !
  !     WHEN ONE WANTS AN ARRAY FOR TAU*SIGMA_K**2, THEN ONE SHOULD REPLACE
  !     THE FUNCTION FTAU BY AN ARRAY
  !
  fun_FTAUK(IIK)= ZTAUEFF*EXP(-ZCOFEXP*(PXIN(IIK)-ZXEXP0)**2 &
       &  /ZXEXPDL**2)
  !%OS      FTAUK(IIK) = ZTAUEFF
  !
  !.......................................................................
  !*COMDECK CUCCCC
  ! ----------------------------------------------------------------------
  ! --     STATEMENT FUNCTION FOR CUBIC INTERPOLATION                   --
  ! --                         23.04.88            AR        CRPP       --
  ! --                                                                  --
  ! -- CUBIC INTERPOLATION OF A FUNCTION F(X)                           --
  ! -- THE EIGHT ARGUMENTS A1,A2,A3,A4,B1,B2,B3,B4 ARE DEFINED BY:      --
  ! -- F(B1) = A1 , F(B2) = A2 , F(B3) = A3 , F(B4) = A4                --
  ! ----------------------------------------------------------------------
  !
  FA3(A1,A2,A3,A4,B1,B2,B3,B4) = &
       &        (A1-A2) / ((B1-B2)*(B2-B4)*(B2-B3)) + &
       &        (A1-A3) / ((B4-B3)*(B3-B1)*(B3-B2)) + &
       &        (A1-A4) / ((B1-B4)*(B2-B4)*(B3-B4))
  FA2(A1,A2,A3,A4,B1,B2,B3,B4) = &
       &        (A1-A2) / ((B2-B1)*(B3-B2)) + &
       &        (A3-A1) / ((B3-B1)*(B3-B2)) - &
       &        (B1+B2+B3) * FA3(A1,A2,A3,A4,B1,B2,B3,B4)
  FA1(A1,A2,A3,A4,B1,B2,B3,B4) = &
       &        (A1-A2) / (B1-B2) - &
       &        (B1+B2) * FA2(A1,A2,A3,A4,B1,B2,B3,B4) - &
       &        (B1*B1+B1*B2+B2*B2) * FA3(A1,A2,A3,A4,B1,B2,B3,B4)
  FA0(A1,A2,A3,A4,B1,B2,B3,B4) = &
       &        A1 - &
       &        B1 * (FA1(A1,A2,A3,A4,B1,B2,B3,B4) + &
       &              B1 * (FA2(A1,A2,A3,A4,B1,B2,B3,B4) + &
       &                    B1 * FA3(A1,A2,A3,A4,B1,B2,B3,B4)))
  ! ----------------------------------------------------------------------
  ! -- FCCCC0 GIVES THE VALUE OF THE FUNCTION AT POINT PX:              --
  ! -- FCCCC0(......,PX) = F(PX)                                        --
  ! ----------------------------------------------------------------------
  FCCCC0(A1,A2,A3,A4,B1,B2,B3,B4,PX) = &
       &              FA0(A1,A2,A3,A4,B1,B2,B3,B4) + &
       &              PX * (FA1(A1,A2,A3,A4,B1,B2,B3,B4) + &
       &                    PX * (FA2(A1,A2,A3,A4,B1,B2,B3,B4) + &
       &                          PX * FA3(A1,A2,A3,A4,B1,B2,B3,B4)))
  ! ----------------------------------------------------------------------
  ! -- FCCCC1 GIVES THE VALUE OF THE FIRST DERIVATIVE OF F(X) AT PX:    --
  ! -- FCCCC1(......,PX) = DF/DX (PX)                                   --
  ! ----------------------------------------------------------------------
  FCCCC1(A1,A2,A3,A4,B1,B2,B3,B4,PX) = &
       &              FA1(A1,A2,A3,A4,B1,B2,B3,B4) + &
       &              PX * (ZTWO * FA2(A1,A2,A3,A4,B1,B2,B3,B4) + &
       &                    ZTHREE * PX * FA3(A1,A2,A3,A4,B1,B2,B3,B4))
  ! ----------------------------------------------------------------------
  ! -- FCCCC2 GIVES THE VALUE OF THE SECOND DERIVATIVE OF F(X) AT PX:   --
  ! -- FCCCC2(......,PX) = D2F/DX2 (PX)                                 --
  ! ----------------------------------------------------------------------
  FCCCC2(A1,A2,A3,A4,B1,B2,B3,B4,PX) = &
       &             ZTWO * FA2(A1,A2,A3,A4,B1,B2,B3,B4) + &
       &             zsix * FA3(A1,A2,A3,A4,B1,B2,B3,B4) * PX
  ! ----------------------------------------------------------------------
  ! -- FCCCC3 GIVES THE VALUE OF THE THIRD DERIVATIVE OF F(X) AT PX:     -
  ! -- FCCCC3(......,PX) = D3F/DX3 (PX)                                  -
  ! ----------------------------------------------------------------------
  FCCCC3(A1,A2,A3,A4,B1,B2,B3,B4,PX) = &
       &                      zsix * FA3(A1,A2,A3,A4,B1,B2,B3,B4)
  !.......................................................................
  !
  !-----------------------------------------------------------------------
  !
  !     0. INITIALIZATION
  !
  ITAUVAL = 0
  IF (PSIG(1) .NE. 0._RKIND) ITAUVAL = 1
  ZTAUEFF = abs(PSIG(1))
  ZXEXP0 = PXEXP0
  ZXEXPDL = PXEXPDL
  ZCOFEXP = 1.0_RKIND
  IF (ZXEXP0.LT.PXIN(1) .OR. ZXEXP0.GT.PXIN(KNIN)) ZCOFEXP=0.0_RKIND
  !
  ISYM = 1
  IF (NBCLFT.GE.10 .OR. NBCRGT.GE.10 .OR. NBCLFT.EQ.2 &
       &  .OR. NBCRGT.EQ.2) ISYM = 0
  !
  N = KNIN
  DO I=1,N
     PYINPP(I) = 0.0_RKIND
  END DO
  DO I=1,MDAMAT
     DO J=1,N
        PAMAT(I,J) = 0.0_RKIND
     END DO
  END DO
  !
  !     0.2 PRE-COMPUTE H_K AND 1./H_K, and zftauk
  !
  DO K=1,N-1
     WHK(K)  = (PXIN(K+1) - PXIN(K))
  enddo
  DO K=1,N-1
     WHK(K)  = (PXIN(K+1) - PXIN(K))
     WOHK(K) = zone / WHK(K)
     ftauk(k)=psig(k)
     !%OS        ftauk(k)=fun_ftauk(k)
  END DO
  WHK(N) = 0.0_RKIND
  WOHK(N) = 0.0_RKIND
  ftauk(n)=psig(N)
  !%OS      ftauk(n)=fun_ftauk(N)
  if (PSIG(1).lt.0._RKIND) ftauk(1)=-10._RKIND*ftauk(1)
  !
  !     0.3 PREPARE BAND WIDTH
  !
  IUP   = 2
  IDOWN = 2
  IF (ITAUVAL .EQ. 0) IUP   = 1
  IF (ITAUVAL .EQ. 0) IDOWN = 1
  IDIAG = IUP + 1
  IF (ISYM .EQ. 0) THEN
     IF (NBCLFT .GE. 10) IUP = IUP + 1
     IF (NBCRGT .GE. 10) IDOWN = IDOWN + 1
     IDIAG = IUP + IDOWN + 1
  ENDIF
  IF (MDAMAT .LT. IUP+1+2*(ISYM-1)*IDOWN) THEN
     PRINT *,' MDAMAT= ',MDAMAT,' < IUP+1+2*(ISYM-1)*IDOWN= ', &
          &    IUP+1+2*(ISYM-1)*IDOWN
     !%OS        STOP 'PAMAT'
     print *,'PAMAT'
     RETURN
  ENDIF
  !
  !     0.4 DETERMINE NEIGHBOURS: K-2, K-1, .., K+2
  !     WHEN OUT OF BOUNDS, POINT TO INDEX N, AS WHK, WOHK(N)=0.0
  !
  DO ISHIFT=-2,+2
     DO K=1,N
        KPM2(K,ISHIFT) = K + ISHIFT
     END DO
  END DO
  !     OUT OF INTERVAL: SEND TO N
  KPM2(1,-2)   = N
  KPM2(1,-1)   = N
  KPM2(2,-2)   = N
  KPM2(N-1,+2) = N
  KPM2(N  ,+1) = N
  KPM2(N  ,+2) = N
  !
  !     1. CONSTRUCT MATRIX AND R.H.S
  !     LAPACK SET-UP OF MATRIX:    A(I,J) -> PAMAT(I-J+IDIAG, J)
  !
  !.......................................................................
  !     AS MATRIX SYMMETRIC, COMPUTE ONLY UPPER PART, THAT IS IF J.GE.I
  !
  DO K=1,N-1
     IEFF = K + IDIAG
     IKP2 = KPM2(K,+2)
     IKP1 = KPM2(K,+1)
     IKM1 = KPM2(K,-1)
     IKM2 = KPM2(K,-2)
     !     A(K,K)
     JK = K
     PAMAT(IEFF-JK,JK) = FAKK(FTAUK(IKP1),FTAUK(K),FTAUK(IKM1),WHK(K) &
          &    ,WHK(IKM1),WOHK(K),WOHK(IKM1))
     !
     !     A(K,K+1)
     JKP1 = K + 1
     PAMAT(IEFF-JKP1,JKP1) = FAKKP1(FTAUK(IKP1),FTAUK(K),WHK(K), &
          &    WOHK(IKP1),WOHK(K),WOHK(IKM1))
     !     A(K,K-1)
     !%OS        JKM1 = K - 1
     !%OS        PAMAT(IEFF-JKM1,JKM1) = FAKKM1(FTAUK(K),FTAUK(IKM1),WHK(IKM1),
     !%OS     +    WOHK(K),WOHK(IKM1),WOHK(IKM2))
     !
     IF (ITAUVAL .EQ. 1) THEN
        !     A(K,K+2)
        JKP2 = K + 2
        IF (JKP2 .LE. N) &
             &      PAMAT(IEFF-JKP2,JKP2) = FAKKP2(FTAUK(IKP1),WOHK(IKP1), &
             &      WOHK(K))
        !     A(K,K-2)
        !%OS          JKM2 = K - 2
        !%OS          PAMAT(IEFF-JKM2,JKM2) = FAKKM2(FTAUK(IKM1),WOHK(IKM1),
        !%OS     +      WOHK(IKM2))
     ENDIF
     !     B(K)
     PYINPP(K) = FRHS(PYIN(IKP1),PYIN(K),PYIN(IKM1),WOHK(K), &
          &    WOHK(IKM1))
     !
  END DO
  !
  !     2. BOUNDARY CONDITIONS
  !
  !     2.1 IF NON-SYMMETRIC, COPY TOP PART TO BOTTOM BEFORE APPLYING
  !     B.C.
  !
  IF (ISYM .EQ. 0) THEN
     DO I=1,N
        IEFF = I + IDIAG
        DO J=I+1,MIN(I+MIN(IUP,IDOWN),N)
           JEFF = J + IDIAG
           !     A(J,I) = A(I,J)
           PAMAT(JEFF-I,I) = PAMAT(IEFF-J,J)
        END DO
     END DO
  ENDIF
  !%OS
  !     debug, print matrix and rhs
  !%OS      write(6,'(3a4,a)') 'i ','j1 ','j2 ',' i,j1  i,j1+1,..... i,j2'
  !%OS      do i=1,n
  !%OS        ieff = idiag + i
  !%OS        j1 = i-idown
  !%OS        j2 = i+iup
  !%OSc%OS        j1 = max(i-idown,1)
  !%OSc%OS        j2 = min(i+iup,n)
  !%OS        write(6,'(3i4,1p10e13.4)') i,j1,j2,(pamat(ieff-j,j),j=j1,j2)
  !%OS      end do
  !%OS      write(6,'(a4,a12)') 'i','RHS'
  !%OS      write(6,'(i4,1pe13.4)') (i,pyinpp(i),i=1,n)
  !
  !%OS
  !
  !     2.2 B.C. AT TWO LOCATIONS PXIN(IXBC(JBC)), JBC=1,2
  !     IBCTYP(JBC) = 0, 1 OR 2 (TYPE OF B.C, SEE ABOVE).
  !     SO FAR USES NODE CLOSEST TO XBCLFT/RGT FOR LOCATION
  !     OF B.C., INSTEAD OF ACTUAL VALUE OF XBCLFT/RGT
  !
  IXBC(1) = 1
  IXBC(2) = N
  IF (NBCLFT .GE. 10) THEN
     DO I=1,KNIN
        IF (PXIN(I) .GE. XBCLFT) GO TO 220
     END DO
220  CONTINUE
     ZDELX = ABS(PXIN(I)-XBCLFT)
     IXBC(1) = I
     IF (I .GE. N) THEN
        IXBC(1) = N
        PRINT *,' WARNING: LEFT B.C. AT I=N: XBCLFT=',XBCLFT, &
             &      '  PXIN(N)= ',PXIN(N)
     ELSE IF (ABS(PXIN(I-1)-XBCLFT).LE.ZDELX .AND. I.NE.1) THEN
        IXBC(1) = I-1
     ENDIF
  ENDIF
  !
  IF (NBCRGT .GE. 10) THEN
     DO I=1,KNIN
        IF (PXIN(I) .GE. XBCRGT) GO TO 221
     END DO
221  CONTINUE
     ZDELX = ABS(PXIN(I)-XBCRGT)
     IXBC(2) = I
     IF (I .LE. 1) THEN
        IXBC(2) = 1
        PRINT *,' WARNING: RIGHT B.C. AT I=1: XBCRGT=',XBCRGT, &
             &      '  PXIN(1)= ',PXIN(1)
     ELSE IF (I .GT. N) THEN
        IXBC(2) = N
     ELSE IF (ABS(PXIN(I-1)-XBCRGT) .LE. ZDELX) THEN
        IXBC(2) = I-1
     ENDIF
  ENDIF
  !
  ZYBC(1) = YBCLFT
  ZYBC(2) = YBCRGT
  IBCTYP(1) = MOD(NBCLFT,10)
  IBCTYP(2) = MOD(NBCRGT,10)
  IF (IXBC(1) .EQ. IXBC(2)) THEN
     PRINT *,' ERROR, B.C. AT SAME LOCATIONS: IXBC(1)=IXBC(2)= ', &
          &    IXBC(1)
     !%OS        STOP '1=2'
     RETURN
  ELSE IF (IXBC(1) .GT. IXBC(2)) THEN
     PRINT *,' WARNING, NEEDED TO SWITCH B.C. POINTS AS IXBC(1)= ', &
          &    IXBC(1),' > IXBC(2)= ',IXBC(2)
     III = IXBC(1)
     IXBC(1) = IXBC(2)
     IXBC(2) = III
     ZYBC(1) = YBCRGT
     ZYBC(2) = YBCLFT
     IBCTYP(1) = MOD(NBCRGT,10)
     IBCTYP(2) = MOD(NBCLFT,10)
  ENDIF
  !
  !     2.3 MOVE EQUATIONS UP OR DOWN IF B.C. IS NOT AN END POINT
  !
  IF (IXBC(1) .NE. 1) THEN
     !
     !     MOVE ROW EQ. K=2,..,IXBC(1) UP BY ONE
     IUPSOFAR = IUP - 1
     DO K=2,IXBC(1)
        DO J=MAX(1,K-IDOWN),MIN(N,K+IUPSOFAR)
           PAMAT(IDIAG+(K-1)-J,J) = PAMAT(IDIAG+K-J,J)
        END DO
        PYINPP(K-1) = PYINPP(K)
        !     ZERO A((K-1),(K-1)-IDOWN)
        IF (K-1-IDOWN .GE. 1) PAMAT(IDIAG+IDOWN,K-1-IDOWN) = 0.0_RKIND
     END DO
     !     ZERO ROW IXBC(1) AND RHS
     K = IXBC(1)
     DO J=MAX(1,K-IDOWN),MIN(N,K+IUP)
        PAMAT(IDIAG+K-J,J) = 0.0_RKIND
     END DO
     PYINPP(K) = 0.0_RKIND
  ENDIF
  !
  IF (IXBC(2) .NE. N) THEN
     !     
     !     MOVE EQ. K=IXBC(2),..,N-1 DOWN BY ONE
     IDWNSOFA = IDOWN - 1
     DO K=N-1,IXBC(2),-1
        DO J=MAX(1,K-IDWNSOFA),MIN(N,K+IUP)
           PAMAT(IDIAG+(K+1)-J,J) = PAMAT(IDIAG+K-J,J)
        END DO
        PYINPP(K+1) = PYINPP(K)
        !     ZERO A((K+1),(K+1)+IUP)
        IF (K+1+IUP .LE. N) PAMAT(IDIAG-IUP,K+1+IUP) = 0.0_RKIND
     END DO
     !     ZERO ROW IXBC(2) AND RHS
     K = IXBC(2)
     DO J=MAX(1,K-IDOWN),MIN(N,K+IUP)
        PAMAT(IDIAG+K-J,J) = 0.0_RKIND
     END DO
     PYINPP(K) = 0.0_RKIND
  ENDIF
  !
  !     2.4 FOR ROW=IXBC(), MODIFY MATRIX AND RHS ACCORDING TO B.C. TYPE
  !
  ZERO = 0.0_RKIND
  DO JBC=1,2
     IK = IXBC(JBC)
     ZVALUE = ZYBC(JBC)
     IEFF = IK + IDIAG
     IKP2 = KPM2(IK,+2)
     IKP1 = KPM2(IK,+1)
     IKM1 = KPM2(IK,-1)
     IKM2 = KPM2(IK,-2)
     IF (IBCTYP(JBC) .EQ. 0) THEN
        !
        !     SYMMETRIZE => COL IK GOES TO RIGHT-HAND SIDE AND THEN ZEROED
        !
        IF (ISYM .EQ. 1) THEN
           IDIAMIK = IDIAG - IK
           DO I=MAX(1,IK-IUP),IK-1
              PYINPP(I) = PYINPP(I) - ZVALUE * PAMAT(I+IDIAMIK,IK)
              PAMAT(I+IDIAMIK,IK) = 0.0_RKIND
           END DO
           IDIAPIK = IDIAG + IK
           DO I=IK+1,MIN(N,IK+IUP)
              PYINPP(I) = PYINPP(I) - ZVALUE * PAMAT(IDIAPIK-I,I)
              PAMAT(IDIAPIK-I,I) = 0.0_RKIND
           END DO
        ELSE
           IDIAMIK = IDIAG - IK
           DO I=MAX(1,IK-IUP),MIN(N,IK+IDOWN)
              PYINPP(I) = PYINPP(I) - ZVALUE * PAMAT(I+IDIAMIK,IK)
              PAMAT(I+IDIAMIK,IK) = 0.0_RKIND
           END DO
           !     ZERO ROW IK
           DO J=MAX(1,IK-IDOWN),MIN(N,IK+IUP)
              PAMAT(IEFF-J,J) = 0.0_RKIND
           END DO
        ENDIF
        !
        !     REPLACE ROW IK BY EQUATION: G_K = ZVALUE
        PAMAT(IDIAG,IK) = 1.0_RKIND
        PYINPP(IK) = ZVALUE
        !
     ELSE IF (IBCTYP(JBC) .EQ. 1) THEN
        !
        !     1ST DERIVATIVE GIVEN
        !
        ZYPEFF = ZVALUE
        IF (ZVALUE .GT. 1.E+31_RKIND) THEN
           !     FROM LGRANGIAN INTERPOLATION
           IKLFT = IK - 1
           IF (IK .EQ. 1) IKLFT = IK
           IF (IKLFT+3 .GT. N) IKLFT = N - 3
           ZYPEFF = FCCCC1(PYIN(IKLFT),PYIN(IKLFT+1),PYIN(IKLFT+2), &
                &        PYIN(IKLFT+3),PXIN(IKLFT),PXIN(IKLFT+1),PXIN(IKLFT+2), &
                &        PXIN(IKLFT+3),PXIN(IK))
        ELSE IF (ZVALUE .LT. -1.E+31_RKIND) THEN
           IKLFT = IK
           IF (IK .EQ. N) IKLFT = IK - 1
           ZYPEFF = (PYIN(IKLFT+1)-PYIN(IKLFT)) &
                &        / (PXIN(IKLFT+1)-PXIN(IKLFT))
        ENDIF
        ZTOHKK1 = FTAUK(IK)*WOHK(IK)*WOHK(IKM1)
        !     A(IK,IK)
        IF (IK .NE. N) PAMAT(IEFF-IK,IK) = FAKK(FTAUK(IKP1),FTAUK(IK), &
             &      ZERO,WHK(IK),ZERO,WOHK(IK),ZERO) + ZTOHKK1
        IF (IK .EQ. N) PAMAT(IEFF-IK,IK) = FAKK(ZERO,FTAUK(IK), &
             &      FTAUK(IKM1),ZERO,WHK(IKM1),ZERO,WOHK(IKM1)) + ZTOHKK1
        !     A(IK,IK-1)
        JKM1 = IK - 1
        IF (ISYM.EQ.0 .AND. JKM1.GE.1) THEN
           IF (IK .NE. N) PAMAT(IEFF-JKM1,JKM1) = - ZTOHKK1
           IF (IK .EQ. N) PAMAT(IEFF-JKM1,JKM1) = FAKKM1(FTAUK(IK), &
                &        FTAUK(IKM1),WHK(IKM1),ZERO,WOHK(IKM1),WOHK(IKM2))
        ENDIF
        !     A(IK,IK+1)
        JKP1 = IK + 1
        IF (JKP1 .LE. N) &
             &      PAMAT(IEFF-JKP1,JKP1) = FAKKP1(FTAUK(IKP1), &
             &      FTAUK(IK),WHK(IK),WOHK(IKP1),WOHK(IK),ZERO)
        !
        IF (ITAUVAL .EQ. 1) THEN
           !     A(IK,IK+2)
           JKP2 = IK + 2
           IF (JKP2 .LE. N) &
                &        PAMAT(IEFF-JKP2,JKP2) = FAKKP2(FTAUK(IKP1),WOHK(IKP1), &
                &        WOHK(IK))
           !     A(IK,IK-2)
           JKM2 = IK - 2
           IF (ISYM.EQ.0 .AND. JKM2.GE.1) THEN
              IF (IK .NE. N) PAMAT(IEFF-JKM2,JKM2) = 0.0_RKIND
              IF (IK .EQ. N) PAMAT(IEFF-JKM2,JKM2) = FAKKM2(FTAUK(IKM1), &
                   &          WOHK(IKM1),WOHK(IKM2))
           ENDIF
        ENDIF
        !     RHS
        ZSIGN = -1._RKIND
        IF (IK .EQ. N) ZSIGN = +1._RKIND
        IF (IK .NE. N) PYINPP(IK) = FRHS(PYIN(IKP1),PYIN(IK),ZERO, &
             &      WOHK(IK),ZERO) - ZYPEFF
        IF (IK .EQ. N) PYINPP(IK) = FRHS(ZERO,PYIN(IK),PYIN(IKM1), &
             &      ZERO,WOHK(IKM1)) + ZYPEFF
        !
     ELSE IF (IBCTYP(JBC) .EQ. 2) THEN
        !
        !     FUNCTION IS GIVEN
        !
        !     A(IK,IK)
        PAMAT(IEFF-IK,IK) = - FTAUK(IK) * (WOHK(IK) + WOHK(IKM1))
        !     A(IK,IK+1)
        JKP1 = IK + 1
        IF (JKP1 .LE. N) &
             &      PAMAT(IEFF-JKP1,JKP1) = FTAUK(IK) * WOHK(IK)
        !     A(IK,IK-1)
        JKM1 = IK - 1
        IF (ISYM.EQ.0 .AND. JKM1.GE.1) &
             &      PAMAT(IEFF-JKM1,JKM1) = FTAUK(IK) * WOHK(IKM1)
        !
        IF (ITAUVAL .EQ. 1) THEN
           !     A(IK,IK+2)
           JKP2 = IK + 2
           IF (JKP2 .LE. N) PAMAT(IEFF-JKP2,JKP2) = 0.0_RKIND
           !     A(IK,IK-2)
           JKM2 = IK - 2
           IF (ISYM.EQ.0 .AND. JKM2.GE.1) PAMAT(IEFF-JKM2,JKM2) = 0.0_RKIND
        ENDIF
        !     RHS
        PYINPP(IK) = PYIN(IK) - ZVALUE
        !
     ENDIF
     !
  END DO
  !
  !     3. SOLVE SYSTEM
  !
  !     USE INTEGER WORK SPACE FROM KPM2(0) ARRAY FOR IPIVOT,
  !     AS KPM2(K,0) NOT NEEDED NOR USED
  !
  IDIMA = MDAMAT
  IDIMRHS = N
  IRHS = 1
  !%OS
  !     debug, print matrix and rhs
  !%OS      write(6,'(3a4,a)') 'i ','j1 ','j2 ',' i,j1  i,j1+1,..... i,j2'
  !%OS      do i=1,n
  !%OS        ieff = idiag + i
  !%OS        j1 = i-idown
  !%OS        j2 = i+iup
  !%OSc%OS        j1 = max(i-idown,1)
  !%OSc%OS        j2 = min(i+iup,n)
  !%OS        write(6,'(3i4,1p10e13.4)') i,j1,j2,(pamat(ieff-j,j),j=j1,j2)
  !%OS      end do
  !%OS      write(6,'(a4,a12)') 'i','RHS'
  !%OS      write(6,'(i4,1pe13.4)') (i,pyinpp(i),i=1,n)
  !
  !%OS
  IF (ISYM .EQ. 1) THEN
     CALL DPBTRF('U',N,IUP,PAMAT,IDIMA,INFO)
  ELSE
     CALL DGBTRF(N,N,IDOWN,IUP,PAMAT,IDIMA,KPM2(1,0),INFO)
  ENDIF
  IF (INFO .EQ. 0) THEN
     IF (ISYM .EQ. 1) THEN
        CALL DPBTRS('U',N,IUP,IRHS,PAMAT,IDIMA,PYINPP,IDIMRHS,INFO2)
     ELSE
        CALL DGBTRS('N',N,IDOWN,IUP,IRHS,PAMAT,IDIMA,KPM2(1,0),PYINPP, &
             &      IDIMRHS,INFO2)
     ENDIF
  ELSE
     PRINT *,' ERROR IN SP/GBTRF: INFO = ',INFO
     STOP 'INFO'
     RETURN
  ENDIF
  !
  !     4. COMPUTE NEW VALUES OF Y_K (NON-STANDARD CUBIC SPLINE ONLY)
  !
  IF (ITAUVAL .EQ. 1) THEN
     DO K=1,N
        IKP1 = KPM2(K,+1)
        IKM1 = KPM2(K,-1)
        PYINNEW(K) = PYIN(K) - FTAUK(K) * &
             &      ((PYINPP(IKP1)-PYINPP(K))*WOHK(K) &
             &      - (PYINPP(K)-PYINPP(IKM1))*WOHK(IKM1))
     END DO
     !
  ENDIF
  
  IF (INFO2 .LT. 0) THEN
     PRINT *,' ERROR IN SP/GBTRS: INFO2 = ',INFO2
     !%OS        STOP 'INFO2'
     RETURN
  ENDIF
  !
  RETURN
END SUBROUTINE CBFITBND
!-----------------------------------------------------------------------
SUBROUTINE SPLIBND(XA,YA,Y2A,N,X,Y,YP,YPP,KOPTIN)
  USE prec_rkind
  implicit none
  REAL(RKIND) :: zsix, zthree, ztwo, zone
  PARAMETER(zsix=6._RKIND, zthree=3._RKIND, ztwo=2._RKIND, zone=1._RKIND)
  REAL(RKIND) :: X, Y, YP, YPP
  integer N, KOPTIN
  REAL(RKIND) :: XA(N),YA(N),Y2A(N)
  !
  !   ABS(KOPTIN):
  !     KOPTIN = 0: STOP WITH ERROR MESSAGE IF OUT OF BOUND
  !     KOPTIN = 1: LINEAR EXTRAPOLATION
  !     KOPTIN = 2: USE QUADRATIC INTERPOLATION IF X OUT OF BOUND
  !     KOPTIN = 3: USE CUBIC INTERPOLATION IF X OUT OF BOUND
  !     KOPTIN = 21: USE QUADRATIC WITHIN ALFA*DELTA_X AND LINEAR FURTHER
  !     KOPTIN = 31: USE CUBIC WITHIN ALFA*DELTA_X AND LINEAR    FURTHER
  !     KOPTIN = 32: USE CUBIC WITHIN ALFA*DELTA_X AND QUADRATIC FURTHER
  !
  !   KOPTIN >=0 => INCONTDER = 1
  !   KOPTIN < 0 => INCONTDER = 0
  !     ICONTDER = 1: VALUE AND 1ST DER. CONTINUOUS AT END OF INTERVAL, THUS
  !     .             USES CUBIC SPLINE OF LAST INTERVAL TO CONTINUE
  !     ICONTDER = 0: ONLY VALUE CONTINUOUS AND USES VALUES AT LAST BUT ONE,
  !     .             TWO, THREE POINTS TO EXTRAPOLATE (BETTER IF DER. AT EDGE
  !     .             IS WILD)
  !
  !-----------------------------------------------------------------------
  !
  REAL(RKIND) :: ALFA
  PARAMETER(ALFA = 1._RKIND)
  !
  !.......................................................................
  !*COMDECK CUCCCC
  ! ----------------------------------------------------------------------
  ! --     STATEMENT FUNCTION FOR CUBIC INTERPOLATION                   --
  ! --                         23.04.88            AR        CRPP       --
  ! --                                                                  --
  ! -- CUBIC INTERPOLATION OF A FUNCTION F(X)                           --
  ! -- THE EIGHT ARGUMENTS A1,A2,A3,A4,B1,B2,B3,B4 ARE DEFINED BY:      --
  ! -- F(B1) = A1 , F(B2) = A2 , F(B3) = A3 , F(B4) = A4                --
  ! ----------------------------------------------------------------------
  !
  REAL(RKIND) :: h, a, b, dellim, zxedg, zyedg, zypedg, zxdellim, zyxdl, &
       &  zypxdl, zxedgm1, zyedgm1, zyxdlp, zypedgm1
  integer icontder, kopt, klo, khi, k, ikopt, k1, k2, k3, k4, klohi
  REAL(RKIND) :: fc3, x1, f1, p1, x2, &
       &  f2, p2, fc2, fc1, fc0, fqqq0, fqqq1, fqqq2, &
       &  flinear, flinearp, fcccc0, fcccc1, fcccc2, fcccc3, fqdq0, fqdq1, &
       &  fqdq2, fcdcd0, fcdcd1, fcdcd2, fcdcd3, fb1, &
       &  fb2, fa2, fa3, fd2, fd1
  REAL(RKIND) :: a1, a2, a3, a4, b1, b2, b3, b4, px
  REAL(RKIND) :: fb0, fd0, fa0, fa1
  FA3(A1,A2,A3,A4,B1,B2,B3,B4) = &
       &        (A1-A2) / ((B1-B2)*(B2-B4)*(B2-B3)) + &
       &        (A1-A3) / ((B4-B3)*(B3-B1)*(B3-B2)) + &
       &        (A1-A4) / ((B1-B4)*(B2-B4)*(B3-B4))
  FA2(A1,A2,A3,A4,B1,B2,B3,B4) = &
       &        (A1-A2) / ((B2-B1)*(B3-B2)) + &
       &        (A3-A1) / ((B3-B1)*(B3-B2)) - &
       &        (B1+B2+B3) * FA3(A1,A2,A3,A4,B1,B2,B3,B4)
  FA1(A1,A2,A3,A4,B1,B2,B3,B4) = &
       &        (A1-A2) / (B1-B2) - &
       &        (B1+B2) * FA2(A1,A2,A3,A4,B1,B2,B3,B4) - &
       &        (B1*B1+B1*B2+B2*B2) * FA3(A1,A2,A3,A4,B1,B2,B3,B4)
  FA0(A1,A2,A3,A4,B1,B2,B3,B4) = &
       &        A1 - &
       &        B1 * (FA1(A1,A2,A3,A4,B1,B2,B3,B4) + &
       &              B1 * (FA2(A1,A2,A3,A4,B1,B2,B3,B4) + &
       &                    B1 * FA3(A1,A2,A3,A4,B1,B2,B3,B4)))
  ! ----------------------------------------------------------------------
  ! -- FCCCC0 GIVES THE VALUE OF THE FUNCTION AT POINT PX:              --
  ! -- FCCCC0(......,PX) = F(PX)                                        --
  ! ----------------------------------------------------------------------
  FCCCC0(A1,A2,A3,A4,B1,B2,B3,B4,PX) = &
       &              FA0(A1,A2,A3,A4,B1,B2,B3,B4) + &
       &              PX * (FA1(A1,A2,A3,A4,B1,B2,B3,B4) + &
       &                    PX * (FA2(A1,A2,A3,A4,B1,B2,B3,B4) + &
       &                          PX * FA3(A1,A2,A3,A4,B1,B2,B3,B4)))
  ! ----------------------------------------------------------------------
  ! -- FCCCC1 GIVES THE VALUE OF THE FIRST DERIVATIVE OF F(X) AT PX:    --
  ! -- FCCCC1(......,PX) = DF/DX (PX)                                   --
  ! ----------------------------------------------------------------------
  FCCCC1(A1,A2,A3,A4,B1,B2,B3,B4,PX) = &
       &              FA1(A1,A2,A3,A4,B1,B2,B3,B4) + &
       &              PX * (ZTWO * FA2(A1,A2,A3,A4,B1,B2,B3,B4) + &
       &                    ZTHREE * PX * FA3(A1,A2,A3,A4,B1,B2,B3,B4))
  ! ----------------------------------------------------------------------
  ! -- FCCCC2 GIVES THE VALUE OF THE SECOND DERIVATIVE OF F(X) AT PX:   --
  ! -- FCCCC2(......,PX) = D2F/DX2 (PX)                                 --
  ! ----------------------------------------------------------------------
  FCCCC2(A1,A2,A3,A4,B1,B2,B3,B4,PX) = &
       &             ZTWO * FA2(A1,A2,A3,A4,B1,B2,B3,B4) + &
       &             zsix * FA3(A1,A2,A3,A4,B1,B2,B3,B4) * PX
  ! ----------------------------------------------------------------------
  ! -- FCCCC3 GIVES THE VALUE OF THE THIRD DERIVATIVE OF F(X) AT PX:     -
  ! -- FCCCC3(......,PX) = D3F/DX3 (PX)                                  -
  ! ----------------------------------------------------------------------
  FCCCC3(A1,A2,A3,A4,B1,B2,B3,B4,PX) = &
       &                      zsix * FA3(A1,A2,A3,A4,B1,B2,B3,B4)
  !-----------------------------------------------------------------------
  !.......................................................................
  !*COMDECK CUCDCD
  ! ----------------------------------------------------------------------
  ! --     STATEMENT FUNCTION FOR CUBIC INTERPOLATION                   --
  ! --                         19.01.87            AR        CRPP       --
  ! --                                                                  --
  ! -- CUBIC INTERPOLATION OF A FUNCTION F(X)                           --
  ! -- THE SIX ARGUMENTS X1,F1,P1,X2,F2,P2 ARE DEFINED AS FOLLOWS:      --
  ! -- F(X1) = F1 , F(X2) = F2 , DF/DX(X1) = P1 , DF/DX(X2) = P2        --
  ! ----------------------------------------------------------------------
  !
  FC3(X1,F1,P1,X2,F2,P2) = &
       &      (ZTWO * (F2 - F1) / (X1 - X2) + (P1 + P2)) / &
       &      ((X1 - X2) * (X1 - X2))
  FC2(X1,F1,P1,X2,F2,P2) = &
       &      (ZTHREE * (X1 + X2) * (F1 - F2) / (X1 - X2) - &
       &       P1 * (X1 + ZTWO * X2) - P2 * (X2 + ZTWO * X1)) / &
       &      ((X1 - X2) * (X1 - X2))
  FC1(X1,F1,P1,X2,F2,P2) = &
       &      (zsix * X1 * X2 * (F2 - F1) / (X1 - X2) + &
       &       X2 * P1 * (2 * X1 + X2) + X1 * P2 * (X1 + ZTWO * X2)) / &
       &      ((X1 - X2) * (X1 - X2))
  FC0(X1,F1,P1,X2,F2,P2) = &
       &      (F1 * X2**2 + F2 * X1**2 - X1 * X2 * (X2 * P1 + X1 * P2) + &
       &       ZTWO * X1 * X2 * (F1 * X2 - F2 * X1) / (X1 - X2)) / &
       &      ((X1 - X2) * (X1 - X2))
  ! ----------------------------------------------------------------------
  ! -- FCDCD0 GIVES THE VALUE OF THE FUNCTION AT POINT PX               --
  ! -- FCDCD0(......,PX) = F(PX)                                        --
  ! ----------------------------------------------------------------------
  FCDCD0(X1,F1,P1,X2,F2,P2,PX) = &
       &              FC0(X1,F1,P1,X2,F2,P2) + &
       &              PX * (FC1(X1,F1,P1,X2,F2,P2) + &
       &                    PX * (FC2(X1,F1,P1,X2,F2,P2) + &
       &                          PX * FC3(X1,F1,P1,X2,F2,P2)))
  ! ----------------------------------------------------------------------
  ! -- FCDCD1 GIVES THE VALUE OF THE FIRST DERIVATIVE OF F(X) AT PX:    --
  ! -- FCDCD1(......,PX) = DF/DX (PX)                                   --
  ! ----------------------------------------------------------------------
  FCDCD1(X1,F1,P1,X2,F2,P2,PX) = &
       &              FC1(X1,F1,P1,X2,F2,P2) + &
       &              PX * (ZTWO * FC2(X1,F1,P1,X2,F2,P2) + &
       &                    ZTHREE * PX * FC3(X1,F1,P1,X2,F2,P2))
  ! ----------------------------------------------------------------------
  ! -- FCDCD2 GIVES THE VALUE OF THE SECOND DERIVATIVE OF F(X) AT PX:   --
  ! -- FCDCD2(......,PX) = D2F/DX2 (PX)                                 --
  ! ----------------------------------------------------------------------
  FCDCD2(X1,F1,P1,X2,F2,P2,PX) = &
       &             ZTWO * FC2(X1,F1,P1,X2,F2,P2) + &
       &             zsix * FC3(X1,F1,P1,X2,F2,P2) * PX
  ! ----------------------------------------------------------------------
  ! -- FCDCD3 GIVES THE VALUE OF THE THIRD DERIVATIVE OF F(X) AT PX:    --
  ! -- FCDCD3(......,PX) = D3F/DX3 (PX)                                 --
  ! ----------------------------------------------------------------------
  FCDCD3(X1,F1,P1,X2,F2,P2,PX) = &
       &                      zsix * FC3(X1,F1,P1,X2,F2,P2)
  !
  !.......................................................................
  !*COMDECK QUAQQQ
  ! ----------------------------------------------------------------------
  ! --     STATEMENT FUNCTION FOR QUADRATIC INTERPOLATION               --
  ! --                         19.01.87            AR        CRPP       --
  ! --                                                                  --
  ! -- QUADRATIC INTERPOLATION OF A FUNCTION F(X)                       --
  ! -- THE SIX PARAMETERS A1,A2,A3,B1,B2,B3 ARE DEFINED AS FOLLOWS:     --
  ! -- F(B1) = A1 , F(B2) = A2 , F(B3) = A3                             --
  ! ----------------------------------------------------------------------
  !
  FB2(A1,A2,A3,B1,B2,B3) = &
       &               ((A1-A2)/(B1-B2)-(A1-A3)/(B1-B3))/(B2-B3)
  FB1(A1,A2,A3,B1,B2,B3) = ((A1-A2)/(B1-B2))- &
       &         FB2(A1,A2,A3,B1,B2,B3)*(B1+B2)
  FB0(A1,A2,A3,B1,B2,B3) = A1-FB1(A1,A2,A3,B1,B2,B3)*B1 &
       &         -FB2(A1,A2,A3,B1,B2,B3)*B1*B1
  ! ----------------------------------------------------------------------
  ! -- FQQQ0 GIVES THE VALUE OF THE FUNCTION AT THE POINT PX            --
  ! -- FQQQ0(......,PX) = F(PX)                                         --
  ! ----------------------------------------------------------------------
  FQQQ0(A1,A2,A3,B1,B2,B3,PX) = FB0(A1,A2,A3,B1,B2,B3) + &
       &                                 PX * (FB1(A1,A2,A3,B1,B2,B3) + &
       &                                 PX * FB2(A1,A2,A3,B1,B2,B3))
  ! ----------------------------------------------------------------------
  ! -- FQQQ1 GIVES THE VALUE OF THE FIRST DERIVATIVE OF F(X) AT PX      --
  ! -- FQQQ1(......,PX) = DF/DX (PX)                                    --
  ! ----------------------------------------------------------------------
  FQQQ1(A1,A2,A3,B1,B2,B3,PX) = FB1(A1,A2,A3,B1,B2,B3) + &
       &     ZTWO * PX * FB2(A1,A2,A3,B1,B2,B3)
  ! ----------------------------------------------------------------------
  ! -- FQQQ2 GIVES THE VALUE OF THE SECOND DERIVATIVE OF F(X) AT PX     --
  ! -- FQQQ2(......,PX) = D2F/DX2 (PX)                                  --
  ! ----------------------------------------------------------------------
  FQQQ2(A1,A2,A3,B1,B2,B3) = ZTWO * FB2(A1,A2,A3,B1,B2,B3)
  !.......................................................................
  !*COMDECK QUAQDQ
  ! ----------------------------------------------------------------------
  ! --     STATEMENT FUNCTION FOR QUADRATIC INTERPOLATION               --
  ! --                         19.01.87            AR        CRPP       --
  ! --                                                                  --
  ! -- QUADRATIC INTERPOLATION OF A FUNCTION F(X)                       --
  ! -- THE FIVE PARAMETERS X1,F1,P1,X2,F2    ARE DEFINED AS FOLLOWS:    --
  ! -- F(X1) = F1 , DF/DX(X1) = P1 , F(X2) = F2                         --
  ! ----------------------------------------------------------------------
  !
  FD2(X1,F1,P1,X2,F2) = ((F2-F1)/(X2-X1) - P1) / (X2-X1)
  FD1(X1,F1,P1,X2,F2) = P1 - ZTWO*X1*FD2(X1,F1,P1,X2,F2)
  FD0(X1,F1,P1,X2,F2) = F1 - X1*(X1*FD2(X1,F1,P1,X2,F2) + &
       &                                     FD1(X1,F1,P1,X2,F2))
  ! ----------------------------------------------------------------------
  ! -- FQDQ0 GIVES THE VALUE OF THE FUNCTION AT POINT PX                --
  ! -- FQDQ0(......,PX) = F(PX)                                         --
  ! ----------------------------------------------------------------------
  FQDQ0(X1,F1,P1,X2,F2,PX) = FD0(X1,F1,P1,X2,F2) + &
       &                              PX * (FD1(X1,F1,P1,X2,F2) + &
       &                                    PX * FD2(X1,F1,P1,X2,F2))
  ! ----------------------------------------------------------------------
  ! -- FQDQ1 GIVES THE VALUE OF THE FIRST DERIVATIVE OF F(X) AT PX      --
  ! -- FQDQ1(......,PX) = DF/DX (PX)                                    --
  ! ----------------------------------------------------------------------
  FQDQ1(X1,F1,P1,X2,F2,PX) = FD1(X1,F1,P1,X2,F2) + &
       &                              ZTWO* PX * FD2(X1,F1,P1,X2,F2)
  ! ----------------------------------------------------------------------
  ! -- FQDQ2 GIVES THE VALUE OF THE SECOND DERIVATIVE OF F(X) AT PX     --
  ! -- FQDQ2(......,PX) = D2F/DX2 (PX)                                  --
  ! ----------------------------------------------------------------------
  FQDQ2(X1,F1,P1,X2,F2) = ZTWO * FD2(X1,F1,P1,X2,F2)
  !-----------------------------------------------------------------------
  !.......................................................................
  !     LINEAR
  !
  FLINEAR(X1,F1,X2,F2,PX) = F2 + (PX-X2)/(X2-X1) * (F2-F1)
  FLINEARP(X1,F1,X2,F2) = (F2-F1) / (X2-X1)
  !-----------------------------------------------------------------------
  ICONTDER = 1
  IF (KOPTIN .LT. 0) ICONTDER = 0
  KOPT=ABS(KOPTIN)
  !
  !     1. POINT INSIDE INTERVAL
  !
  IF (X.LT.XA(1) .OR. X.GT.XA(N)) GO TO 200
  !
  KLO=1
  KHI=N
1 IF (KHI-KLO.GT.1) THEN
     K=(KHI+KLO)/2
     IF(XA(K).GT.X)THEN
        KHI=K
     ELSE
        KLO=K
     ENDIF
     GOTO 1
  ENDIF
  H=XA(KHI)-XA(KLO)
  !%OS      IF (H.EQ.0._RKIND) STOP 'BAD XA INPUT.'
  IF (H .EQ. 0._RKIND) THEN
     PRINT *,'BAD XA INPUT.'
     RETURN
  ENDIF
  A=(XA(KHI)-X)/H
  B=(X-XA(KLO))/H
  Y=A*YA(KLO)+B*YA(KHI)+ &
       &      ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/zsix
  YP=(YA(KHI)-YA(KLO))/H - &
       & ((ZTHREE*A*A-zone)*Y2A(KLO)-(ZTHREE*B*B-ZONE)*Y2A(KHI) )*H/zsix
  YPP=A*Y2A(KLO)+B*Y2A(KHI)
  !
  RETURN
  !
  !     2. POINT OUTSIDE INTERVAL
  !
200 CONTINUE
  !
  IKOPT = KOPT
  IF (KOPT .EQ. 0) THEN
     PRINT *,' POINT X=',X,' IS OUTSIDE INTERVAL [',XA(1),',' &
          &    ,XA(N),']'
     IKOPT = 1
     !C        STOP 'KOPT=0'
  ENDIF
  !
  KLO = 1
  KHI = 2
  IF (X .GT. XA(N)) THEN
     KLO = N-1
     KHI = N
  ENDIF
  H=XA(KHI)-XA(KLO)
  DELLIM = ALFA * H
  !%OS      IF (H.EQ.0._RKIND) STOP 'BAD XA INPUT.'
  IF (H .EQ. 0._RKIND) THEN
     PRINT *,'BAD XA INPUT.'
     RETURN
  ENDIF
  !
  !.......................................................................
  !     2.1 LINEAR, IKOPT=1
  !
  IF (IKOPT .EQ. 1) THEN
     !
     !     LINEAR EXTRAPOLATION
     IF (ICONTDER .EQ. 0) THEN
        Y = YA(KHI)
        YP = 0.0_RKIND
        YPP = 0.0_RKIND
        !%OS          Y = FLINEAR(XA(KLO),YA(KLO),XA(KHI),YA(KHI),X)
        !%OS          YP = FLINEARP(XA(KLO),YA(KLO),XA(KHI),YA(KHI))
        !%OS          YPP = 0.0_RKIND
     ELSE
        !     COMPUTE VALUE AND 1ST DER. AT EDGE
        IF (KLO .EQ. 1) THEN
           ZXEDG = XA(KLO)
           ZYEDG = YA(KLO)
           ZYPEDG= (YA(KHI)-YA(KLO))/H-(ZTWO*Y2A(KLO)+Y2A(KHI))*H/zsix
        ELSE
           ZXEDG = XA(KHI)
           ZYEDG = YA(KHI)
           ZYPEDG= (YA(KHI)-YA(KLO))/H+(Y2A(KLO)+ZTWO*Y2A(KHI))*H/zsix
        ENDIF
        Y = FLINEAR(ZXEDG,ZYEDG,ZXEDG+ZONE,ZYEDG+ZYPEDG,X)
        YP = ZYPEDG
        YPP = 0.0_RKIND
     ENDIF
     !
     !.......................................................................
     !     2.2 LINEAR FAR END OF IKOPT=21 AND IKOPT=31
     !
  ELSE IF (MAX(XA(KLO)-X,X-XA(KHI)).GT.DELLIM .AND. &
       &    (IKOPT.EQ.21 .OR. IKOPT.EQ.31)) THEN
     !
     !     LINEAR EXTRAPOLATION OUTSIDE DELLIM
     !     COMPUTE STARTING POINT AND DERIVATIVE FROM END OF QUADR. OR CUBIC
     !     INTERPOLATION IN ALFA*DELTA_X INTERVAL
     !
     IF (IKOPT .EQ. 21) THEN
        !     QUADRATIC
        IF (ICONTDER .EQ. 0) THEN
           K1 = 1
           IF (KHI .EQ. N) K1 = N-2
           K2 = K1 + 1
           K3 = K1 + 2
           ZXDELLIM = XA(1) - DELLIM
           IF (KHI .EQ. N) ZXDELLIM = XA(N) + DELLIM
           ZYXDL = FQQQ0(YA(K1),YA(K2),YA(K3),XA(K1),XA(K2),XA(K3), &
                &        ZXDELLIM)
           ZYPXDL= FQQQ1(YA(K1),YA(K2),YA(K3),XA(K1),XA(K2),XA(K3), &
                &        ZXDELLIM)
        ELSE
           !     COMPUTE VALUE AND 1ST DER. AT EDGE
           IF (KLO .EQ. 1) THEN
              ZXEDGM1 = XA(KHI)
              ZYEDGM1 = YA(KHI)
              ZXEDG = XA(KLO)
              ZYEDG = YA(KLO)
              ZYPEDG=(YA(KHI)-YA(KLO))/H-(ZTWO*Y2A(KLO)+Y2A(KHI))*H/zsix
           ELSE
              ZXEDGM1 = XA(KLO)
              ZYEDGM1 = YA(KLO)
              ZXEDG = XA(KHI)
              ZYEDG = YA(KHI)
              ZYPEDG=(YA(KHI)-YA(KLO))/H+(Y2A(KLO)+ZTWO*Y2A(KHI))*H/zsix
           ENDIF
           ZXDELLIM = XA(1) - DELLIM
           IF (KHI .EQ. N) ZXDELLIM = XA(N) + DELLIM
           ZYXDL = FQDQ0(ZXEDG,ZYEDG,ZYPEDG,ZXEDGM1,ZYEDGM1,ZXDELLIM)
           ZYPXDL= FQDQ1(ZXEDG,ZYEDG,ZYPEDG,ZXEDGM1,ZYEDGM1,ZXDELLIM)
        ENDIF
     ELSE IF (IKOPT .EQ. 31) THEN
        !     CUBIC
        IF (ICONTDER .EQ. 0) THEN
           K1 = 1
           IF (KHI .EQ. N) K1 = N-3
           K2 = K1 + 1
           K3 = K1 + 2
           K4 = K1 + 3
           ZXDELLIM = XA(1) - DELLIM
           IF (KHI .EQ. N) ZXDELLIM = XA(N) + DELLIM
           IF (K1.LT.1 .OR. K4.GT.N) THEN
              IF (KHI .EQ. N) K1 = N-1
              K2 = K1 + 1
              ZYXDL=FLINEAR(XA(K1),YA(K1),XA(K2),YA(K2),ZXDELLIM)
              ZYXDLP=FLINEARP(XA(K1),YA(K1),XA(K2),YA(K2))
           ELSE
              ZYXDL = FCCCC0(YA(K1),YA(K2),YA(K3),YA(K4),XA(K1),XA(K2), &
                   &          XA(K3),XA(K4),ZXDELLIM)
              ZYPXDL= FCCCC1(YA(K1),YA(K2),YA(K3),YA(K4),XA(K1),XA(K2), &
                   &          XA(K3),XA(K4),ZXDELLIM)
           ENDIF
        ELSE
           !     COMPUTE VALUE AND 1ST DER. AT EDGE AND EDGE-1
           IF (KLO .EQ. 1) THEN
              ZXEDGM1 = XA(KHI)
              ZYEDGM1 = YA(KHI)
              ZYPEDGM1 = (YA(KHI)-YA(KLO))/H + (Y2A(KLO)+ZTWO*Y2A(KHI)) &
                   &          *H/zsix
              ZXEDG = XA(KLO)
              ZYEDG = YA(KLO)
              ZYPEDG=(YA(KHI)-YA(KLO))/H-(ZTWO*Y2A(KLO)+Y2A(KHI))*H/zsix
           ELSE
              !%OS              IF (H.LE.1.E-02_RKIND) THEN
              !%OS                KLO=KLO-1
              !%OS                H=XA(KHI)-XA(KLO)
              !%OS              endif
              ZXEDGM1 = XA(KLO)
              ZYEDGM1 = YA(KLO)
              ZYPEDGM1 =(YA(KHI)-YA(KLO))/H-(ZTWO*Y2A(KLO)+Y2A(KHI)) &
                   &          *H/zsix
              ZXEDG = XA(KHI)
              ZYEDG = YA(KHI)
              ZYPEDG=(YA(KHI)-YA(KLO))/H+(Y2A(KLO)+ZTWO*Y2A(KHI))*H/zsix
           ENDIF
           ZXDELLIM = XA(1) - DELLIM
           IF (KHI .EQ. N) ZXDELLIM = XA(N) + DELLIM
           ZYXDL = FCDCD0(ZXEDGM1,ZYEDGM1,ZYPEDGM1,ZXEDG,ZYEDG,ZYPEDG, &
                &        ZXDELLIM)
           ZYPXDL= FCDCD1(ZXEDGM1,ZYEDGM1,ZYPEDGM1,ZXEDG,ZYEDG,ZYPEDG, &
                &        ZXDELLIM)
        ENDIF
     ENDIF
     !
     Y = FLINEAR(ZXDELLIM,ZYXDL,ZXDELLIM+ZONE,ZYXDL+ZYPXDL,X)
     YP = ZYPXDL
     YPP = 0.0_RKIND
     !
     !.......................................................................
     !     2.3 QUADRATIC, IKOPT=2 OR FIRST PART OF IKOPT=21
     !
  ELSE IF (IKOPT.EQ.2 .OR. &
       &    (MAX(XA(KLO)-X,X-XA(KHI)).LE.DELLIM .AND. IKOPT.EQ.21)) THEN
     !
     !     QUADRATIC
     IF (ICONTDER .EQ. 0) THEN
        K1 = 1
        IF (KHI .EQ. N) K1 = N-2
        K2 = K1 + 1
        K3 = K1 + 2
        Y =  FQQQ0(YA(K1),YA(K2),YA(K3),XA(K1),XA(K2),XA(K3),X)
        YP = FQQQ1(YA(K1),YA(K2),YA(K3),XA(K1),XA(K2),XA(K3),X)
        YPP =FQQQ2(YA(K1),YA(K2),YA(K3),XA(K1),XA(K2),XA(K3))
     ELSE
        !     COMPUTE VALUE AND 1ST DER. AT EDGE
        IF (KLO .EQ. 1) THEN
           ZXEDGM1 = XA(KHI)
           ZYEDGM1 = YA(KHI)
           ZXEDG = XA(KLO)
           ZYEDG = YA(KLO)
           ZYPEDG= (YA(KHI)-YA(KLO))/H-(ZTWO*Y2A(KLO)+Y2A(KHI))*H/zsix
        ELSE
           ZXEDGM1 = XA(KLO)
           ZYEDGM1 = YA(KLO)
           ZXEDG = XA(KHI)
           ZYEDG = YA(KHI)
           ZYPEDG= (YA(KHI)-YA(KLO))/H+(Y2A(KLO)+ZTWO*Y2A(KHI))*H/zsix
        ENDIF
        Y  = FQDQ0(ZXEDG,ZYEDG,ZYPEDG,ZXEDGM1,ZYEDGM1,X)
        YP = FQDQ1(ZXEDG,ZYEDG,ZYPEDG,ZXEDGM1,ZYEDGM1,X)
        YPP= FQDQ2(ZXEDG,ZYEDG,ZYPEDG,ZXEDGM1,ZYEDGM1)
     ENDIF
     !
     !.......................................................................
     !     2.4 QUADRATIC, FAR END PART OF IKOPT=32
     !
  ELSE IF (MAX(XA(KLO)-X,X-XA(KHI)).GT.DELLIM .AND. IKOPT.EQ.32)THEN
     !
     !     QUADRATIC FROM X+ALFA*DELTA_X
     !     COMPUTE STARTING POINT AND DERIVATIVE FROM END OF CUBIC
     !     WARNING: MUST BE COMPATIBLE WITH ALFA=0
     IF (ICONTDER .EQ. 0) THEN
        K1 = 1
        IF (KHI .EQ. N) K1 = N-3
        K2 = K1 + 1
        K3 = K1 + 2
        K4 = K1 + 3
        ZXDELLIM = XA(1) - DELLIM
        IF (KHI .EQ. N) ZXDELLIM = XA(N) + DELLIM
        IF (K1.LT.1 .OR. K4.GT.N) THEN
           IF (KHI .EQ. N) K1 = N-1
           K2 = K1 + 1
           ZYXDL=FLINEAR(XA(K1),YA(K1),XA(K2),YA(K2),ZXDELLIM)
           ZYXDLP=FLINEARP(XA(K1),YA(K1),XA(K2),YA(K2))
        ELSE
           ZYXDL = FCCCC0(YA(K1),YA(K2),YA(K3),YA(K4),XA(K1),XA(K2), &
                &        XA(K3),XA(K4),ZXDELLIM)
           ZYPXDL= FCCCC1(YA(K1),YA(K2),YA(K3),YA(K4),XA(K1),XA(K2), &
                &        XA(K3),XA(K4),ZXDELLIM)
        ENDIF
     ELSE
        !     COMPUTE VALUE AND 1ST DER. AT EDGE AND EDGE-1
        IF (KLO .EQ. 1) THEN
           ZXEDGM1 = XA(KHI)
           ZYEDGM1 = YA(KHI)
           ZYPEDGM1 = (YA(KHI)-YA(KLO))/H + (Y2A(KLO)+ZTWO*Y2A(KHI)) &
                &        *H/zsix
           ZXEDG = XA(KLO)
           ZYEDG = YA(KLO)
           ZYPEDG= (YA(KHI)-YA(KLO))/H-(ZTWO*Y2A(KLO)+Y2A(KHI))*H/zsix
        ELSE
           ZXEDGM1 = XA(KLO)
           ZYEDGM1 = YA(KLO)
           ZYPEDGM1 = (YA(KHI)-YA(KLO))/H - (ZTWO*Y2A(KLO)+Y2A(KHI)) &
                &        *H/zsix
           ZXEDG = XA(KHI)
           ZYEDG = YA(KHI)
           ZYPEDG= (YA(KHI)-YA(KLO))/H+(Y2A(KLO)+ZTWO*Y2A(KHI))*H/zsix
        ENDIF
        ZXDELLIM = XA(1) - DELLIM
        IF (KHI .EQ. N) ZXDELLIM = XA(N) + DELLIM
        ZYXDL = FCDCD0(ZXEDGM1,ZYEDGM1,ZYPEDGM1,ZXEDG,ZYEDG,ZYPEDG, &
             &      ZXDELLIM)
        ZYPXDL= FCDCD1(ZXEDGM1,ZYEDGM1,ZYPEDGM1,ZXEDG,ZYEDG,ZYPEDG, &
             &      ZXDELLIM)
     ENDIF
     !     QUADRATIC FROM END OF INTERVAL, END+DELLIM AND 1ST DER. AT END+DELLIM
     KLOHI = 1
     IF (KHI .EQ. N) KLOHI = N
     IF (ABS(ALFA/H).GT.1E-04_RKIND) THEN
        Y  = FQDQ0(ZXDELLIM,ZYXDL,ZYPXDL,XA(KLOHI),YA(KLOHI),X)
        YP = FQDQ1(ZXDELLIM,ZYXDL,ZYPXDL,XA(KLOHI),YA(KLOHI),X)
        YPP =FQDQ2(ZXDELLIM,ZYXDL,ZYPXDL,XA(KLOHI),YA(KLOHI))
     ELSE
        Y  = FQDQ0(ZXDELLIM,ZYXDL,ZYPXDL,XA(KLO),YA(KLO),X)
        YP = FQDQ1(ZXDELLIM,ZYXDL,ZYPXDL,XA(KLO),YA(KLO),X)
        YPP =FQDQ2(ZXDELLIM,ZYXDL,ZYPXDL,XA(KLO),YA(KLO))
     ENDIF
     !
     !.......................................................................
     !     2.5 CUBIC, IKOPT=3 OR FIRST PART OF IKOPT=31 AND IKOPT=32
     !
  ELSE IF (IKOPT.EQ.3 .OR. &
       &    (MAX(XA(KLO)-X,X-XA(KHI)).LE.DELLIM .AND. &
       &    (IKOPT.EQ.32 .OR. IKOPT.EQ.31)) ) THEN
     !
     !     CUBIC
     IF (ICONTDER .EQ. 0) THEN
        K1 = 1
        IF (KHI .EQ. N) K1 = N-3
        K2 = K1 + 1
        K3 = K1 + 2
        K4 = K1 + 3
        IF (K1.LT.1 .OR. K4.GT.N) THEN
           IF (KHI .EQ. N) K1 = N-1
           K2 = K1 + 1
           Y =FLINEAR(XA(K1),YA(K1),XA(K2),YA(K2),X)
           YP=FLINEARP(XA(K1),YA(K1),XA(K2),YA(K2))
           YPP=0._RKIND
        ELSE
           Y   =FCCCC0(YA(K1),YA(K2),YA(K3),YA(K4),XA(K1),XA(K2),XA(K3) &
                &        ,XA(K4),X)
           YP  =FCCCC1(YA(K1),YA(K2),YA(K3),YA(K4),XA(K1),XA(K2),XA(K3) &
                &        ,XA(K4),X)
           YPP =FCCCC2(YA(K1),YA(K2),YA(K3),YA(K4),XA(K1),XA(K2),XA(K3) &
                &        ,XA(K4),X)
        ENDIF
     ELSE
        !     COMPUTE VALUE AND 1ST DER. AT EDGE AND EDGE-1
        IF (KLO .EQ. 1) THEN
           ZXEDGM1 = XA(KHI)
           ZYEDGM1 = YA(KHI)
           ZYPEDGM1=(YA(KHI)-YA(KLO))/H+(Y2A(KLO)+ZTWO*Y2A(KHI))*H/zsix
           ZXEDG = XA(KLO)
           ZYEDG = YA(KLO)
           ZYPEDG= (YA(KHI)-YA(KLO))/H-(ZTWO*Y2A(KLO)+Y2A(KHI))*H/zsix
        ELSE
           !%OS            IF (H.LE.1.E-02_RKIND) THEN
           !%OS              KLO=KLO-1
           !%OS              H=XA(KHI)-XA(KLO)
           !%OS            endif
           ZXEDGM1 = XA(KLO)
           ZYEDGM1 = YA(KLO)
           ZYPEDGM1=(YA(KHI)-YA(KLO))/H-(ZTWO*Y2A(KLO)+Y2A(KHI))*H/zsix
           ZXEDG = XA(KHI)
           ZYEDG = YA(KHI)
           ZYPEDG= (YA(KHI)-YA(KLO))/H+(Y2A(KLO)+ZTWO*Y2A(KHI))*H/zsix
        ENDIF
        Y  = FCDCD0(ZXEDGM1,ZYEDGM1,ZYPEDGM1,ZXEDG,ZYEDG,ZYPEDG,X)
        YP = FCDCD1(ZXEDGM1,ZYEDGM1,ZYPEDGM1,ZXEDG,ZYEDG,ZYPEDG,X)
        YPP= FCDCD2(ZXEDGM1,ZYEDGM1,ZYPEDGM1,ZXEDG,ZYEDG,ZYPEDG,X)
     ENDIF
  ENDIF
  !
  RETURN
END SUBROUTINE SPLIBND
