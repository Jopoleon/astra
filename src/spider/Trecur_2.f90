!***********************************************************************
!    SUBROUTINE  FOR  READING  AND  TREATMENT
!    PFC  PARAMETERS  ( THEIR POSITIONS, CURRENTS, RESISTANCE,
!                       NUMBER OF TURNS, NUMBER OF DIVIDING,
!                       EQUIVALENT SIGNS )
!***********************************************************************
!  OUTPUT DATE:
!  ----------
!   NPFC      - NUMBER OF PF-COILS
!
!   PFCUR1(L) - VALUES OF TOTAL CURRENT OF PFC CROSS-SECTION (in MA)
!               L = 1,2,...,NPFC
!   NEPFC(L)  - EQUIVALENT "SIGNS"  OF PF-COILS
!               L = 1,2,...,NPFC
!   WEPFC(L)  - EQUIVALENT "WEIGHTS" OF PF-COILS
!               L = 1,2,...,NPFC
!   NTURN(L)  - NUMBER OF REAL TURNS  OF PF-COILS
!               L = 1,2,...,NPFC
!   NDIV(L)   - NUMBER OF DIVIDING TURNS OF PF-COILS
!               L = 1,2,...,NPFC
!   NLOC(L)   - LOCATION ARRAY OF DIVIDING TURNS OF PF-COILS
!
!               !!! L = 0,1,2,...,NPFC !!! :
!               NLOC(0) = 0, NLOC(1) = NDIV(1),
!               NLOC(2) = NDIV(1) + NDIV(2),...,
!               NLOC(NPFC) = NDIV(1) + NDIV(2) +...+ NDIV(NPFC) = NCPFC
!  ---------------------------------------------------------------------
!   NEQUI     - NUMBER OF EQUIVALENT PF-COIL GROUPS ( NEQUI .LE. NPFC )
!
!   PFRES(L)  - RESISTANCE OF EQUIVALENT PF-COIL GROUPS; in [micro*Ohm]
!               L = 1,2,...,NEQUI
!   PFVOL1(L) - VOLTAGE    OF EQUIVALENT PF-COIL GROUPS; in [V]
!               L = 1,2,...,NEQUI
!   PFCEQW(L) - CURRENT    OF EQUIVALENT PF-COIL GROUPS; in [MA]
!               L = 1,2,...,NEQUI
!  ---------------------------------------------------------------------
!   NCPFC    - NUMBER OF PFC-CURRENTS (AFTER DIVIDING)
!
!   PC(L)    - VALUES (IN MA) OF PFC-CURRENTS (AFTER DIVIDING)
!              L=1,...,NCPFC.
!   RI(L),ZI(L) - CILINDER COORDINATES (IN M) OF PFC-CURRENTS
!                 (AFTER DIVIDING)  L = 1,2,...,NCPFC
!   NTYPE(L) - NUMBER OF TYPE OF PFC-CURRENT CROSS-SECTIONS
!              (AFTER DIVIDING)  L = 1,2,...,NCPFC
!   VERS(L)  - VERTICAL   SIZE OF PFC-CURRENT CROSS-SECTIONS
!              (AFTER DIVIDING)  L = 1,2,...,NCPFC
!   HORS(L)  - HORIZONTAL SIZE OF PFC-CURRENT CROSS-SECTIONS
!              (AFTER DIVIDING)  L = 1,2,...,NCPFC
!   NECON(L) - EQUIVALENT "SIGNS"  OF PF-CURRENTS
!              (AFTER DIVIDING)  L = 1,2,...,NCPFC
!   WECON(L) - EQUIVALENT "WEIGHTS" OF PF-CURRENTS
!              (AFTER DIVIDING)  L = 1,2,...,NCPFC
!
!***********************************************************************

subroutine TRECUR(NCPFC, RI, ZI, PC, NTYPE, NECON, WECON, &
                  HORS, VERS)

use sp_parameters, only: nloopp, nprobp, njlim, nplim, npfc0
use iopath, only: path
use comevl, only: nequi, npfc, nloc, nturn, nepfc, ndiv, &
    pfcd1, pfceqw, pfcur1, pfcur2, pfcw1, pfres, pfvol1, wepfc

implicit none

integer, intent(out) :: NCPFC
integer, intent(out), dimension(*) :: ntype, necon
real*8, intent(out), dimension(*) :: PC, HORS, VERS, RI, ZI, WECON

integer :: i, j, l, NDIVW, NDIVH, NTIPE, NDI, NDIVRE, NINP, KEYCUR
integer, dimension(300) :: NDW, NDH
real*8 :: horc, verc, PS, awc, ahc, wc, hc, rc, zc
real*8, dimension(700) :: RS, ZS
real*8, dimension(300) :: DELZ, DELR
character(len=80) :: fname

!***********************************************************************

NINP  = 1

write(fname,'(a,a)') TRIM(path), '/coilres.dat'
open(1,file=fname,form='formatted')
   read(1,*)  NEQUI
   do I=1,NEQUI
      read(1,*) L
      read(1,*)(pfres(L,j),j=1,NEQUI)
   enddo
close(1)		
write(*,*) 'res read'

write(fname,'(a,a)') TRIM(path), '/coil.dat'
open(1,file=fname,form='formatted')
   read(1,*) NPFC, KEYCUR
   if(NPFC.gt.NPFC0) then
      WRITE(*,*) '**************************************'
      WRITE(*,*) '  NUMBER OF PF-COILS            NPFC =', NPFC
      WRITE(*,*) '  ABOVE MAXIMAL STOP            NPFC0 =', NPFC0
      stop
   endif

   print *,' npfc ==',npfc

   NCPFC   = 0
   NLOC(0) = 0

   do I=1,NPFC
      read(1,*) NDI
      read(1,*) RC,ZC,WC,HC,AWC,AHC,PFCUR1(I),NTURN(I),NEPFC(I)
      NTIPE = 2
!------------------------------------------------------
! EXPERIMENT FOR TCV PF FAST COILS (numbers 24,...,29)
!+++   IF( I.GE.24 )  NTIPE = 3
!------------------------------------------------------
      write(*,*) RC,ZC,WC,HC,AWC,AHC,PFCUR1(I),NTURN(I),NEPFC(I)
      IF( KEYCUR .EQ. 1 )  PFCUR1(I) = PFCUR1(I) * NTURN(I)
      CALL DIVPAR( RC, ZC, WC, HC, AWC, AHC, PFCUR1(I), NDI, &
                   NDIVRE, NDIVW, NDIVH, RS, ZS, PS, VERC, HORC )
      NDIV(I)   = NDIVRE
      WEPFC(I)  = NTURN(I)  / NDIV(I)
      PFCUR2(I) = PFCUR1(I) 
      PFCW1(I)  = PFCUR1(I) / NTURN(I)
      PFCD1(I)  = PFCUR1(I) / NDIV(I)

      do J=1,NDIVRE
         L  = NCPFC + J
         RI(L) = RS(J)
         ZI(L) = ZS(J)
         PC(L) = PFCD1(I)
         HORS(L) = HORC
         VERS(L) = VERC
         NTYPE(L) = NTIPE
         NECON(L) = NEPFC(I)
         WECON(L) = WEPFC(I)
      enddo
      NCPFC     = NCPFC + NDIVRE
      NLOC(I)   = NCPFC
      NDW(I)    = NDIVW
      NDH(I)    = NDIVH
      DELZ(I)   = VERC
      DELR(I)   = HORC
   enddo
close(1)

do I=1,NEQUI
   do J=1,NPFC
      if( NEPFC(J) .EQ. I ) then
         PFCEQW(I) = PFCW1(J)
         CYCLE
      endif
   enddo
   PFVOL1(I) = PFRES(I,i)*PFCEQW(I) ! [Volt] = [micro*Ohm]*[mega*Amper]
enddo

write(*,*) 'coil read'

return
end subroutine TRECUR

!**************************************************************
!        SUBROUTINE FOR DIVIDING OF PARALLELOGRAM
!                ( FOR PFC CROSS-SECTION )
!**************************************************************
!  INPUT DATE:
!  ----------
!   RC,ZC - CILINDER COORDINATES OF CENTER OF PARALLELOGRAM
!   WC    - PROJECTION OF THE FIRST SIDE OF PARALLELOGRAM
!           ON AXIS "R" (IN METER)
!   HC    - PROJECTION OF THE SECOND SIDE OF PARALLELOGRAM
!           ON AXIS "Z" (IN METER)
!   AWC   - ANGLE BETWEEN THE FIRST SIDE OF PARAL. AND
!           AXIS "R" (IN DEGREES, IT MUST NOT BE EQUAL 90 ,
!                                 IT MUST BE :  -90 < AWC < 90 )
!   AHC   - ANGLE BETWEEN THE SECOND SIDE OF PARAL. AND
!           AXIS "R" (IN DEGREES, IT MUST NOT BE EQUAL 0 ,
!                                 IT MUST BE :  0 < AHC < 180 )
!   CURC  - CURRENT OF PARALLELOGRAM (IN MA)
!   NDIV  - APPROXIMATE NUMBER OF CELLS OF DIVIDING:
!           NDIV=0 - A SPECIAL CASE: AUTOMATICALLY NDIVRE=1,
!                    RS(1)=RC, ZS(1)=ZC, PS=CURC
!           IF NDIV > 0 THEN WE HAVE THE MOST TOTAL ALGORITHM OF
!                            DIVIDING
!           IF NDIV < 0 THEN NDIVW AND NDIVH ARE CUT OFF BY ABS(NDIV)
!
!  OUTPUT DATE:
!  ----------
!   NDIVRE      - REAL NUMBER OF CELLS OF DIVIDING ( = NDIVW*NDIVH )
!   NDIVW       - NUMBER OF DIVIDING OF THE FIRST  SIDE OF PARAL.
!   NDIVH       - NUMBER OF DIVIDING OF THE SECOND SIDE OF PARAL.
!   RS(L),ZS(L) - CILINDER COORDINATES OF CENTERS OF CELLS OF DIVIDING
!                 L = 1,2,...,NDIVRE ! (IN METER)
!   PS          - CURRENT OF EVERY CELL OF DIVIDING  (IN MA)
!   VERS - VERTICAL (OR LINEAR, OR RADIUS) SIZE OF CELL CROSS-S.
!   HORS - HORIZONTAL SIZE OF CELL CROSS-SECTION
!
!**************************************************************

subroutine DIVPAR( RC, ZC, WC, HC, AWC, AHC, CURC, NDIV, NDIVRE, &
                   NDIVW, NDIVH, RS, ZS, PS, VERS, HORS )

use sp_parameters, only: pi

implicit none

integer, intent(in) :: NDIV
real*8, intent(in) :: Rc, ZC, WC, HC, AWC, CURC

integer, intent(out) :: NDIVH, NDIVW, NDIVRE
real*8, intent(out), dimension(*) ::  RS, ZS
real*8, intent(out) :: PS, VERS, HORS
real*8, intent(inout) :: AHC

integer :: NDIVA, i, j, L, IDINT
real*8 :: WR, WZ, HR, HZ, SH, SW, AHCR, AWCR, R0, Z0, HSIZE, WSIZE

if(NDIV.EQ.0) then
   NDIVW  = 1
   NDIVH  = 1
   NDIVRE = 1
   RS(1)  = RC
   ZS(1)  = ZC
   PS     = CURC
   VERS   = HC
   HORS   = WC
   return
endif

if(AHC.LT.0) AHC = AHC + 180.

if(NDIV.GE.0) then
   NDIVA =  NDIV
else
   NDIVA = -NDIV
endif

!***************************************
!   PARAMETERS OF PARALLELOGRAM

AWCR  = AWC * pi/180.
AHCR  = AHC * pi/180.
R0    = RC - 0.5*( WC + HC * DCOS(AHCR)/DSIN(AHCR) )
Z0    = ZC - 0.5*( HC + WC * DSIN(AWCR)/DCOS(AWCR) )
WSIZE = WC / DCOS(AWCR)
HSIZE = HC / DSIN(AHCR)
WR    = WC
WZ    = WC * DSIN(AWCR) / DCOS(AWCR)
HR    = HC * DCOS(AHCR) / DSIN(AHCR)
HZ    = HC

!***************************************
!   CALCULATION  NDIVW, NDIVH, NDIVRE, PS

SW    = DSQRT( NDIVA*WSIZE/HSIZE )
SH    = DSQRT( NDIVA*HSIZE/WSIZE )
SW    = SW + 0.5
SH    = SH + 0.5
NDIVW = IDINT(SW)
NDIVH = IDINT(SH)

IF(NDIVW.EQ.0)  NDIVW = 1
IF(NDIVH.EQ.0)  NDIVH = 1
IF((NDIV.LT.0).AND.(NDIVW.GT.NDIVA)) NDIVW = NDIVA
IF((NDIV.LT.0).AND.(NDIVH.GT.NDIVA)) NDIVH = NDIVA

NDIVRE = NDIVW * NDIVH
PS     = CURC / NDIVRE

!***************************************
!   CALCULATION  RS(L), ZS(L) : L = 1,2,...,NDIVRE

WR  = WR / NDIVW
WZ  = WZ / NDIVW
HR  = HR / NDIVH
HZ  = HZ / NDIVH
HORS = WR
VERS = HZ
RS(1) = R0 + 0.5*(WR + HR)
ZS(1) = Z0 + 0.5*(WZ + HZ)

do I=1,NDIVW
   do J=1,NDIVH
      L     = (I-1)*NDIVH + J
      RS(L) = RS(1) + (I-1)*WR + (J-1)*HR
      ZS(L) = ZS(1) + (I-1)*WZ + (J-1)*HZ
   enddo
enddo

return
end subroutine DIVPAR
