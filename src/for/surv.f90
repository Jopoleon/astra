!----------------------------------------------------------------
subroutine SCAL(NOUT, SN, SO, OUT, NP, NDIM)
!----------------------------------------------------------------
!  Input:
!     NOUT Total number of channels (NROUT or NTOUT)
!     SO(NOUT) Current scales
!     NP Upper index boundary for scale definition (= NDIM-3)
!     NDIM Upper index boundary according to the array description
!     OUT(NDIM, NOUT) Array of curves for scale definition (ROUT or TOUT)
!  Output:
!     SN(NOUT) Returned scales

implicit none

integer, intent(in) :: NOUT, NP, NDIM
double precision, intent(in)  :: SO(*), OUT(NDIM, *)
double precision, intent(out) :: SN(*)

integer :: J, JJ
double precision :: SC, SCALA
external SCALA

! Exclude 3 edge/central points
var_loop: do J=1, NOUT
   if (SO(J) > 0) then
      SN(J) = SO(J)
   else if (SO(J) .EQ. 0) then
      SN(J) = SCALA(OUT(1, J), NP)
   else
      do JJ = 1, J-1
         if(SO(J) .EQ. SO(JJ)) then
            SN(J) = SN(JJ)
            CYCLE var_loop
         endif
      enddo
      SC = 0.
      do JJ=J, NOUT
         if(SO(J) .EQ. SO(JJ)) then
            SC = MAX(SC, SCALA(OUT(1, JJ), NP))
            if(ABS(SO(JJ) + JJ) < .01) then
               SN(J) = SCALA(OUT(1, JJ), NP)
               CYCLE var_loop
            endif
         endif
      enddo
      SN(J) = SC
   endif
enddo var_loop

return
end subroutine SCAL

!-------------------------------------------
double precision function SCALA(Y, NJ)

implicit none

integer, intent(in) :: NJ
double precision, intent(in) :: Y(*)

integer :: j
double precision :: YS(6), YMAX

YMAX = MAXVAL(ABS(Y(1: NJ)))

YS(1) = 1.0d-9
YS(2) = 1.5d-9
YS(3) = 2.0d-9
YS(4) = 3.0d-9
YS(5) = 5.0d-9
YS(6) = 8.0d-9

do while(1.05*YMAX > YS(6))
   YS = 10.*YS
enddo
do j=1, 6
   if (1.05*YMAX <= YS(j)) EXIT
enddo

SCALA = YS(J)

end function SCALA

!-------------------------------------------
subroutine CMARK(NL, JPOS, SC, OS, NAME, STYL)
! Mark variable/scale in 1 & 2 modes

use outcmn_inc, only: frame_hei, DXLET, DYLET
use dbl2char, only: fmt4

implicit none

integer, intent(in) :: STYL, JPOS, NL
double precision, intent(in) :: SC, OS
character(len=4), intent(in) :: NAME

integer :: JY, POINT(2)
character(len=10) :: ST
character(len=4 ) :: F4
character(len=5 ) :: F5

ST(1: 4) = fmt4(SC)
ST(5: 5)  = ' '
ST(6: 10) = NAME
if (OS /= 0) then
   F4 = fmt4(abs(OS))
   if (OS < 0) F5 = '-' // F4
   if (OS > 0) F5 = '+' // F4
   if (NL < frame_hei/2) then
      JY = DYLET + 2
   else
      JY = -DYLET - 2
   endif
   call textvm(JPOS + 4*DXLET, NL + JY, F5, 5)
endif
call textvm(JPOS, NL, ST, 10)
if (STYL > 0) then
   POINT(1) = JPOS + 4*DXLET + 4
   POINT(2) = NL - 4
   call NMARK(POINT, STYL)
endif

return
end subroutine CMARK

!-------------------------------------------
subroutine CMARKT(NL, JPOS, SC, OS, NAME, STYL)
! Mark variable/scale in 6th (time) mode

use outcmn_inc, only: DXLET, DYLET
use dbl2char, only: fmt4

implicit none

integer, intent(in) :: NL, JPOS, STYL
double precision, intent(in) :: SC, OS
character(len=4), intent(in) :: NAME

integer :: J, JP, plot_arr(2)
character(len=10) :: ST
character(len=4 ) :: F4

ST(1: 10) = '          '
ST(2: 5)  = NAME
if (OS == 0) then
   JP = JPOS
   J  = 5
else
   F4 = fmt4(abs(OS))
   if (OS < 0) ST(6: 10) = '-' // F4
   if (OS > 0) ST(6: 10) = '+' // F4
   JP = JPOS - DXLET
   if (ST(1: 1) == ' ') JP = JP - DXLET
   j = LEN_TRIM(ST(2: 10))
   if (j <= 6) JP = JP + DXLET
endif

call textvm(JP, NL, ST, J)   ! type name+yshift
plot_arr(1) = JP + 3
plot_arr(2) = NL - 5
F4 = fmt4(SC)
J = NL + DYLET
JP = JPOS + DXLET
call textvm(JP, J, F4, 4)   ! type scale
if (STYL > 0) call NMARK(plot_arr, STYL)

return
end subroutine CMARKT

!-------------------------------------------
subroutine CMARKP(NL, JPOS, NAME, STYL)
! Mark variable/scale in 4 & 5 (time-radial) modes

implicit none

character(len=4), intent(in) :: NAME
integer, intent(in) :: NL, JPOS, STYL

integer :: NBIT(6)
character(len=5) :: ST

! diamond(1), o(111), +(43), *(42), x(120), #(35), $(36), 
data NBIT / 1, 111, 43, 42, 120, 36 /
ST(2: 5) = NAME
ST(1: 1) = ' '
if (STYL >= 7 .and. STYL <= 12) ST(1:1) = char(NBIT(STYL - 6))
call textvm(JPOS, NL, ST, 5)

end subroutine CMARKP

!-------------------------------------------
subroutine PLOTXY(YARR, NP, JX, IX, IXO, IY, IYO, DMET, STYL, plot_arr)

! The subroutine displays NP points
! of integer array IY vs IX to screen with the style=STYL
! and puts marks with time interval equal to DMET(sec)
! Entry: YARR - time array (TTOUT)
! NP, JX, IY, IX, DMET, STYL

use outcmn_inc, only: EraseColor

implicit none

integer, intent(in) :: NP, STYL, JX, IX(*), IY(*)
integer, intent(out), dimension(*) :: IXO, IYO, plot_arr
double precision, intent(in) :: DMET, YARR(NP)

integer :: J, J0, JMET, JJ, JPOINT, J1
double precision :: DMETO

save DMETO
data DMETO/999999./

call colovm(EraseColor)
JPOINT = 0
do J0=0, 1
   do J=1, NP
      plot_arr(2*J - 1) = JX  + IXO(J)
      plot_arr(2*J)     = 350 - IYO(J)
   enddo
   JMET = 0
   do JJ=1, NP - 1 + J0
      if(YARR(JJ) - YARR(1) >= JMET*DMETO .or. JJ == 1) then
         if(JJ > 1) call curvvm(0, JPOINT + 1, plot_arr(J1))
         J1 = 2*JJ - 1
         call NMARK(plot_arr(J1), STYL)
         JMET = JMET + 1
         JPOINT = 1
      else
         JPOINT = JPOINT+1
      endif
   enddo
   if (JPOINT /= 0) call curvvm(0, JPOINT, plot_arr(J1))
   if (J0 == 1) return   ! J0=0 <- erasing
   do J=1, NP
      IXO(J) = IX(J)
      IYO(J) = IY(J)
   enddo
   DMETO = DMET
   call colovm(1)
enddo

return
end subroutine PLOTXY

!-------------------------------------------
subroutine PLOTCR(NP, NPO, IX, IXOLD, IY, IYOLD, ICOLOR, STYL, plot_arr)

! The subroutine displays NP points of the integer array IY
! NP  is a number of points to plot
! NPO is a number of points to erase, in addition, 
! NPO is a control parameter:
! NPO > 0  the old curve is erased, the drawn one is stored in IYOLD 
! NPO <= 0 a new curve IY(1:NP) is drawn, (IXOLD, IYOLD) are NOT used
! NPO = 0  no erasure, the drawn curve is stored in IYOLD
! The points of the array IYOLD are used for erasing
! curve of the previous call and are determined inside PLOTG1
! STYL
! plot_arr is a working array 2*NB1
! Input: NP, IY, IYOLD, ICOLOR, STYL
! Output: IXOLD, IYOLD

use outcmn_inc, only: EraseColor

implicit none

integer, intent(in) :: STYL, NP, NPO, ICOLOR
integer, intent(in) :: IX(*), IY(*)
integer, intent(out) :: plot_arr(*)
integer, intent(inout) :: IXOLD(*), IYOLD(*)

integer :: J

if (NPO > 0) then
! erase the old curve
   do J=1, NPO
      plot_arr(2*J - 1) = IXOLD(J)
      plot_arr(2*J)     = IYOLD(J)
   enddo
   call colovm(EraseColor)
   call CURV1(NPO, plot_arr, STYL)
endif

! draw a new curve

do J=1, NP
   plot_arr(2*J - 1) = IX(J)
   plot_arr(2*J)     = IY(J)
enddo
! Colors: 1(Red) 2(Blue) 3(MeduimSeeGreen) 4(VioletRed) 5(Brown) 6(LightBlue)
! 7(Turquoise)
! STYL 1 2 3 4 5 6 7  8 9 10 11 12 13 14 15 16 17 18
!olor: 1 2 3 4 5 6 7   1 2 3  4  5  6  7   1  2  3  4  5  6 7
call colovm(ICOLOR)
call CURV1(NP, plot_arr, STYL)
if (NPO < 0) return
do J=1, NP
   IYOLD(J) = IY(J)
   IXOLD(J) = IX(J)
enddo

return
end subroutine PLOTCR

!-------------------------------------------
subroutine CURV1(NP, plot_arr, STYL)
! The subroutine has replaced older subroutine CURV

implicit none

integer, intent(in) :: plot_arr(*), STYL, NP

integer :: LE, NF, J, j1, JJ, NM, PT1(2)

if (STYL < 0) then  ! Draw dashed curves
   jj = -STYL
   if (jj >= 7) jj = jj + 1 - jj/7*7
   if (jj > 1) then
      LE = 8
      do j=1, jj
         j1 = j + 1
         LE = LE + j1
      enddo
      j1 = jj
      if (jj == 2) LE = min(LE, 16)
      if (jj == 3) LE = min(LE, 8)
      if (jj == 4) LE = min(LE, 4)
      NF = max(1, LE/4)
      do j=1, NP, LE
         j1 = min(NP - j + 1, LE - NF)
         call drcurv(0, j1, plot_arr(2*j-1))
      enddo
      return
   endif

else if (STYL > 0) then

   NM = NP/5
   NM = max(10, NP/5)
   LE = NM/5*STYL    ! 1st marker position
   if (LE >= NM+2) LE = LE - NM
   LE = max(1, LE)
! call getcolor(j1)
   do jj=LE, NP, NM
      J  = 2*jj
      PT1(1) = plot_arr(j-1)
      PT1(2) = plot_arr(j)/10
      call NMARK(PT1, STYL)
!    call puto(jx, jy, j1, STYL)
   enddo
endif

call drcurv(0, NP, plot_arr(1))

end subroutine CURV1

!-------------------------------------------
subroutine NMARK_(POINT, STYL)
! used when NMARK_ is called from C

implicit none

integer, intent(in) :: POINT(2), STYL

call NMARK(POINT, STYL)

end subroutine NMARK_

!-------------------------------------------
subroutine NMARK(POINT, STYL)

implicit none

integer, parameter, dimension(7) :: N=(/16, 13, 5, 9, 14, 9, 10/)

integer, intent(in) :: POINT(2), STYL

integer :: plot_arr(32), DX(16, 7), DY(16, 7), J, JJ, IST

! IST definition shoud coincide with NBIT() in CMARK
! IST = 1-filled diamond, 2-o, 3-+, 4-*, 5-x(#), 6-<, 7-filled square
! STYL  <=7,              8    9    10   11      12   >=13

save DX, DY ! Actually parameters
data DX/ &
   0, 3, 0, -3, 0, 0, 2, 0, -2, 0, 0, 1, 0, -1, 0, 0, &
   3, 3, 2, 1, -1, -2, -3, -3, -2, -1, 1, 2, 3, 0, 0, 0, &
   3, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
  -3, -1, 2, -2, 1, -2, 2, -1, 3, 0, 0, 0, 0, 0, 0, 0, &
  -2, -1, 0, 1, 2, 1, 0, 1, 2, 1, 0, -1, -2, 0, 0, 0, &
   0, 0, -3, 0, -2, 0, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, &
  -2, 2, 2, -2, -2, 1, 1, -1, -1, 0, 0, 0, 0, 0, 0, 0/ 
data DY/ &
   3, 0, -3, 0, 3, 2, 0, -2, 0, 2, 1, 0, -1, 0, 1, 0, &
   1, -1, -2, -3, -3, -2, -1, 1, 2, 3, 3, 2, 1, 0, 0, 0, &
   0, 0, 0, 3, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
   0, 0, 3, -3, 0, 3, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
  -2, -1, 0, 1, 2, 1, 0, -1, -2, -1, 0, 1, 2, 0, 0, 0, &
   3, -3, 0, 3, -1, -1, 1, 1, -3, 0, 0, 0, 0, 0, 0, 0, &
  -2, -2, 2, 2, -1, -1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0/  

IST = max(1, min(STYL, 7))
do JJ=1, N(IST)
   J = 2*JJ
   plot_arr(J-1) = POINT(1) + DX(JJ, IST)
   plot_arr(J)   = POINT(2) + DY(JJ, IST)
enddo
call curvvm(0, N(IST), plot_arr(1))

end subroutine NMARK

!-------------------------------------------
subroutine READF6(LINE, F6, IERR)
! Read a number in LINE to 6-positinal field 

implicit none

character(len=6), intent(in) :: LINE
integer, intent(out) :: IERR
double precision, intent(out) :: F6

integer :: N, JE, JP, JM, J
character(len=6) :: STRI

IERR = 0
JE = 0
JP = 0
JM = 0
do J=1, 6
   if(LINE(J: J) == '-') JM = J
   if(LINE(J: J) == 'e' .or. LINE(J: J) == 'E') JE = J
   if(LINE(J: J) == '.') JP = J
enddo

if(JP > 0) then
   READ(LINE, '(F6.3)', ERR=2) F6
   return
endif

if(JE <= 0) then
   if(JM <= 1) then
      READ(LINE, '(I6)', ERR=2) N
      F6 = N
      return
   endif
   STRI = LINE(1: JM-1)
   READ(STRI, '(I6)', ERR=2) N
   F6 = N
   STRI = LINE(JM: 6)
   READ(STRI, '(I6)', ERR=2) N
   F6 = F6 * 10.**N
   return
else
   STRI = LINE(1: JE-1)
   READ(STRI, '(I6)', ERR=2) N
   F6 = N
   STRI = LINE(JE+1: 6)
   READ(STRI, '(I6)', ERR=2) N
   F6 = F6 * 10.**N
endif

2 continue

IERR = 1
write(*, *) '>>> READF6: found ERROR in "', LINE, '"'

return
end subroutine READF6

!---------------------------------------------------------------------=|
double precision function GETNUM(FIELD, ERCODE)
!----------------------------------------------------------------------|

use parameter_inc
use outcmn_inc, only: PRNAME, CFNAME
use const_inc , only: CONSTF, DEVARX
use char_manip, only: str_in_list
use debugger, only: debug

implicit none

integer, intent(out) :: ERCODE
character(len=*), intent(in) :: FIELD

integer :: l, j, jnam, jpos, jpos1, j1, ISHIFT, ios
character(len=6) :: ZNUM

save ISHIFT
data ISHIFT/0/

if (ISHIFT == 0) then
   j = str_in_list('ZRD1  ', PRNAME)
   ISHIFT = max(j-1, 0)
endif

l  = len(FIELD)
jpos = index(TRIM(FIELD), 'ZRD')

if (jpos == 0) then

   jpos1 = index(TRIM(FIELD), 'C')
   ERCODE = 1 ! Only initial, before CFNAME-search
   if (jpos1 >  0) then
      j1 = min(LEN_TRIM(FIELD(jpos1: l)), l - jpos1 + 1)
      ZNUM = FIELD(jpos1: jpos1+j1-1)
      jnam = str_in_list(ZNUM, CFNAME)
      if (jnam > 0) then
         ERCODE = 0
         GETNUM = CONSTF(jnam)
      endif
   endif

else

   j1 = index(FIELD(jpos+3: l), 'X')
   if (j1 /= 0) then
      if (j1 > 3) then
         ERCODE = 1
         return
      endif
      ZNUM = FIELD(jpos+3: jpos+j1+1)
   else
      ZNUM = FIELD(jpos+3:)
   endif

   read(ZNUM(1:), *, iostat=ios) j1

   if (ios /= 0 .or. j1 < 1 .or. j1 > 96) then
      ERCODE = 1
   else
      GETNUM = DEVARX(ISHIFT + j1)
      ERCODE = 0
   endif

endif

return
end function GETNUM

!-------------------------------------------
double precision function ROUNDN(R1, N)

implicit none

integer, intent(in) :: N
double precision, intent(in) :: R1
   
integer :: J
double precision :: R
   
if (R >= 1.e13 .or. R < 1.e-9) then
   ROUNDN = R1
else
   R = R1*1.e9
   do J=1, 25
      if (R < 10.) EXIT
      R = R/10.
   enddo
   ROUNDN = (R + 50./10.**N)*10.**(J - 10)
endif

end function ROUNDN
   
!---------------------------------------------------------------------=|
subroutine STREAD(NCH, NFIELD, ARRAY, ERCODE)
!----------------------------------------------------------------------|
! Reads one record group of the NBINP ("*.nbi") file 
! and fills ARRAY(1:NFIELD) with data.
! Numbers and references to ZRD*, ZRD*X and CONSTF_list
! are allowed as records in the input file.
! ERCODE values:
!   0 - Normal exit
!   1 - Unrecognized variable name
!   2 - Read error (never occurring, protected by ISNUM)
!   3 - Wrong format
!   4 - Array out of limits
!   5 - Missing records
!----------------------------------------------------------------------|

use char_manip, only: to_upper
use dbl2char, only: isnum
use debugger, only: debug

implicit none

integer, intent(in) :: NCH, NFIELD
integer, intent(out) :: ERCODE
double precision, intent(out) :: ARRAY(NFIELD)

integer :: j, ios
double precision :: GETNUM
character(len=12 ) :: SFIELD(20)
character(len=132) :: str_line
!----------------------------------------------------------------------|

if (NFIELD > 20) then
   ERCODE = 4
   return
endif

! Skip all lines beginning with '!'
j = 1
do while (j == 1)
   read(NCH, '(A)', iostat=ios) str_line
   if (ios < 0) then ! EOF encountered
      ERCODE = 5
      return
   else if (ios > 0) then
      ERCODE = 1
      return
   endif
   j = index(str_line, '!')
enddo

if (j /= 0) then
   ERCODE = 3
   write(*, *) 'STREAD error: exclamation marks allowed only at line beginning'
else  ! Read numbers and/or variable names
   ERCODE = 5  ! Missing entries
   read(NCH, '(5A)', iostat=ios) (SFIELD(j), j=1, NFIELD)
   if (ios < 0) then ! EOF encoutnered
      ERCODE = 5
   else if (ios > 0) then
      ERCODE = 1     ! Error reading, probably never occurring
   else
      ERCODE = 0
      do j=1, NFIELD
         SFIELD(J) = to_upper(SFIELD(j))
         if ( ISNUM(SFIELD(j), 12) ) then
            read(SFIELD(j), *) ARRAY(j) ! read err never occurs, protected by ISNUM
         else ! In case it is a variable name, like ZRD*, pick its value
            ARRAY(j) = GETNUM(SFIELD(j), ERCODE)
            if (ERCODE /= 0) then ! Unrecognised variable name
               ERCODE = 1
               EXIT
            endif
         endif
      enddo
   endif
endif

return
end subroutine STREAD

!---------------------------------------------------------------------=|
subroutine ASKINT(NV, NVAR, NAME)

use parameter_inc, only: NRW

implicit none

integer, intent(in) :: NV
integer, intent(inout) :: NVAR(*)
character(len=4), intent(in) :: NAME(*)

integer :: J
double precision :: VAR(NRW)

do J=1, NV
   VAR(J) = NVAR(J)
enddo
call ASKLIS(NV, VAR, NAME, 4)
do J=1, NV
   NVAR(J) = VAR(J)
enddo

return
end subroutine ASKINT

!---------------------------------------------------------------------=|
subroutine ASTWIN(NB, IBOX, NAME, SCALE, SHIFT, MOD10, YMODE)

use parameter_inc, only: NRW
use outcmn_inc, only: IP1, IP2, IP30, IP31, null_ch
use char_manip, only: to_upper

implicit none

integer, parameter :: NRW16=NRW+16, NRW96=NRW-128

integer, intent(in) :: NB, YMODE, MOD10
integer, intent(out) :: IBOX(*)
double precision, intent(in) :: SCALE(*), SHIFT(*)
character(len=4), intent(out) :: NAME(*)

integer :: JMODE, JGR, j, j1, j2, jj, jn, jm, js, jb, jw, jsep, ASKTAB
integer, dimension(NRW) :: IB
character(len=80) :: rows(NRW16), STR, TITLE
character(len=1) :: KEY

TITLE = "Presentation" // null_ch
!              ----5----0----5----0----5----0----5----0----5----0----5
STR = "Name|Box| Scale|Offset||Name|Box| Scale|Offset" // null_ch

jsep = index(STR, '||') + 1
JMODE = MOD10
if (JMODE == 1) then
   JGR = 8
elseif (JMODE == 2 .or. JMODE == 3) then
   JGR = 4
elseif (JMODE == 6) then
    if (YMODE ==  1) JGR = 4
    if (YMODE ==  0) JGR = 4
    if (YMODE == -1) JGR = 2
else
    return
endif

KEY = 'Y'

ask_key: do while(KEY == 'Y')

   do j = 1, NRW
      write(rows(j)(1:80), '(79X, 1A1)') null_ch
   enddo

! j  - ordinal box No.
! jb - box No. in the Astra nominations
! jw - position in the table
! jj - horizontal row
! js - position in the current row

   jn = 0
   jm = 0
   do j = 1, NB
      jb = IBOX(j)
      if (jb <= 0 ) then
         jm = jm + 1
         CYCLE
      endif
      if (JMODE == 1) jw = IP1(jb)
      if (JMODE == 2 .or. JMODE == 3) jw = IP2(jb)
      if (JMODE == 6) then
         if (YMODE == 1) jw = jb
         if (YMODE == 0) jw = IP30(jb)
         if (YMODE  == -1) jw = IP31(jb)
      endif
      js = jsep*(1 - jw + jw/2*2)
      jj = 1 + (jw - 1)/2
      write(rows(jj)(js+1:js+4), '(1A4)') NAME(j)
      write(rows(jj)(js+6:js+8), '(1I3)') jb
      call num2str(SCALE(j), rows(jj)(js+10:js+15), 6)
      call num2str(SHIFT(j), rows(jj)(js+17:js+22), 6)
      jn = max(jn, jj)
   enddo

   if (jm /= 0 ) then
      j2 = jn
      j1 = 2*jn + 1
      do j = 1, NB
         if (IBOX(j) <= 0 ) then
            jw = j1
            js = jsep*(1 - jw + jw/2*2)
            jj = 1 + (jw - 1)/2
            write(rows(jj)(js+1: js+4), '(1A4)') NAME(j)
            write(rows(jj)(js+6: js+8), '(1I3)') IBOX(j)
            call num2str(SCALE(j), rows(jj)(js+10: js+15), 6)
            call num2str(SHIFT(j), rows(jj)(js+17: js+22), 6)
            jn = max(jn, jj)
            j1 = j1 + 1
         endif
      enddo
   endif

   j = 1
   do while(j > 0)
      j1 = (jm + 1)/2
      j = ASKTAB(TITLE, STR, rows, 80, jn, JGR, j1)
   enddo
   
   jn = 0
   do j=1, NB
      if (IBOX(j) >  0 ) then
         if (JMODE == 1) jw = IP1(IBOX(j))
         if (JMODE == 2 .or. JMODE == 3) jw = IP2(IBOX(j))
         if (JMODE == 6) then
            if (YMODE  == 1) jw = IBOX(j)
            if (YMODE  == 0) jw = IP30(IBOX(j))
            if (YMODE == -1) jw = IP31(IBOX(j))
         endif
         js = jsep*(1 - jw + jw/2*2)
         jj = 1 + (jw - 1)/2
         write(NAME(j), '(1A4)', ERR=77) rows(jj)(js+1: js+4)
         read(rows(jj)(js+6 : js+8 ), *, ERR=77) IB(j)
         read(rows(jj)(js+10: js+15), *, ERR=77) SCALE(j)
         read(rows(jj)(js+17: js+22), *, ERR=77) SHIFT(j)
         jn = max(jn, jj)
      endif
   enddo

   if (jm /= 0 ) then
      j1 = 2*jn + 1
      do j=1, NB
         if (IBOX(j) <= 0 ) then
            jw = j1
            js = jsep*(1 - jw + jw/2*2)
            jj = 1 + (jw - 1)/2
            write(NAME(j), '(1A4)', ERR=77) rows(jj)(js+1: js+4)
            read(rows(jj)(js+6: js+8), *, ERR=77) IB(j)
            read(rows(jj)(js+10:js+15), *, ERR=77) SCALE(j)
            read(rows(jj)(js+17:js+22), *, ERR=77) SHIFT(j)
            jn = max(jn, jj)
            j1 = j1 + 1
         endif
      enddo
   endif

   do j = 1, NB
      IBOX(j) = IB(j)
   enddo
   
   return

77 continue

   write(*, *)
   write(*, '(A)') '>>> INPUT ERROR encountered in the dialog window "Presentation"'
   write(*, '(A23, A52, A1)') '                Line: "', rows(j)(1: 52), '"'
   write(*, '(A52, $)') '     Enter "Y" to return, any other key to ignore > '
   KEY = 'X'
   read(*, '(:, A1)') KEY
   KEY = to_upper(KEY)

enddo ask_key

return
end subroutine ASTWIN

!---------------------------------------------------------------------=|
subroutine ASXWIN(NB, IBOX, NAME, SCALE, YSHIFT, XL, XR, MOD10, YMODE)

use parameter_inc, only: NRW
use const_inc, only: AB
use outcmn_inc, only: IP1, IP2, IP30, IP31, null_ch

implicit none

integer, parameter :: NRW16=NRW+16, NRW96=NRW-128

integer, intent(in) :: NB, YMODE, MOD10
integer, intent(out) :: IBOX(*)
double precision, intent(in), dimension(*) :: SCALE, YSHIFT, XL, XR
character(len=4), intent(out) :: NAME(*)

integer :: JMODE, JGR, j, j1, j2, jj, jn, jm, js, jb, jw, jsep, ASKTAB
integer, dimension(NRW) :: IB
double precision :: YY
character(len=80) :: rows(NRW16), STR, TITLE
character(len=1) :: YKEY
character(len=6) :: ABNUM

TITLE = "Presentation" // null_ch
!              ----5----0----5----0----5----0----5----0----5----0----5
STR = "Name|Box| Scale|Offset| [a?, ]| [, a?]||" // &
      "Name|Box| Scale|Offset| [a?, ]| [, a?]"   // null_ch
jsep = index(STR, '||')+1

JMODE = MOD10
if (JMODE == 1) then
   JGR = 8
elseif (JMODE == 2 .or. JMODE == 3) then
   JGR = 4
elseif (JMODE == 6) then
   if (YMODE ==  1) JGR = 4
   if (YMODE ==  0) JGR = 4
   if (YMODE == -1) JGR = 2
else
    return
endif

10 continue

do j = 1, NRW
   write(rows(j)(1:80), '(79X, 1A1)') null_ch
enddo

! j  - ordinal box No.
! jb - box No. in the Astra nominations
! jw - position in the table
! jj - horizontal row
! js - position in the current row

jn = 0
jm = 0
call num2str(AB, ABNUM, 6)
do j = 1, NB
   jb = IBOX(j)
   if (jb <= 0 ) then
      jm = jm + 1
      CYCLE
   endif
   if (JMODE == 1) jw = IP1(jb)
   if (JMODE == 2 .or. JMODE == 3) jw = IP2(jb)
   if (JMODE == 6) then
      if (YMODE == 1) jw = jb
      if (YMODE == 0) jw = IP30(jb)
      if (YMODE  == -1) jw = IP31(jb)
   endif
   YY = XR(j)
   js = jsep*(1 - jw + jw/2*2)
   jj = 1 + (jw - 1)/2
   write(rows(jj)(js+1: js+4), '(1A4)') NAME(j)
   write(rows(jj)(js+6: js+8), '(1I3)') jb
   call num2str(SCALE(j) , rows(jj)(js+10: js+15), 6)
   call num2str(YSHIFT(j), rows(jj)(js+17: js+22), 6)
   call num2str(XL(j)    , rows(jj)(js+24: js+29), 6)
   call num2str(XR(j)    , rows(jj)(js+31: js+36), 6)
   if (YY > AB) rows(jj)(js+31:js+36) = ABNUM
   jn = max(jn, jj)
enddo

if (jm /= 0 ) then
   j2 = jn
   j1 = 2*jn+1
   do j = 1, NB
      if (IBOX(j) <= 0 ) then
         YY = XR(j)
         jw = j1
         js = jsep*(1 - jw + jw/2*2)
         jj = 1 + (jw - 1)/2
         write(rows(jj)(js+1: js+4), '(1A4)') NAME(j)
         write(rows(jj)(js+6: js+8), '(1I3)') IBOX(j)
         call num2str(SCALE(j) , rows(jj)(js+10: js+15), 6)
         call num2str(YSHIFT(j), rows(jj)(js+17: js+22), 6)
         call num2str(XL(j)    , rows(jj)(js+24: js+29), 6)
         call num2str(XR(j)    , rows(jj)(js+31: js+36), 6)
         if (YY > AB) rows(jj)(js+31: js+36) = ABNUM
         jn = max(jn, jj)
         j1 = j1 + 1
      endif
   enddo
endif

j = 1
do while(j > 0)
   j1 = (jm + 1)/2
   j = ASKTAB(TITLE, STR, rows, 80, jn, JGR, j1)
enddo

jn = 0
do j = 1, NB
   if (IBOX(j) > 0) then
      if (JMODE == 1) jw = IP1(IBOX(j))
      if (JMODE == 2 .or. JMODE == 3) jw = IP2(IBOX(j))
      if (JMODE == 6) then
         if (YMODE ==  1) jw = IBOX(j)
         if (YMODE ==  0) jw = IP30(IBOX(j))
         if (YMODE == -1) jw = IP31(IBOX(j))
      endif
      js = jsep*(1 - jw + jw/2*2)
      jj = 1 + (jw - 1)/2
      write(NAME(j), '(1A4)', ERR=77) rows(jj)(js+1: js+4)
      read(rows(jj)(js+6 : js+8), *, ERR=77) IB(j)
      read(rows(jj)(js+10: js+15), *, ERR=77) SCALE(j)
      read(rows(jj)(js+17: js+22), *, ERR=77) YSHIFT(j)
      read(rows(jj)(js+24: js+29), *, ERR=77) XL(j)
      read(rows(jj)(js+31: js+36), *, ERR=77) YY
      if (rows(jj)(js+31: js+36) /= ABNUM) then
         read(rows(jj)(js+31: js+36), *, ERR=77) XR(j)
      endif
      jn = max(jn, jj)
   endif
enddo

if (jm /= 0 ) then
   j1 = 2*jn + 1
   do j=1, NB
      if (IBOX(j) <= 0) then
         jw = j1
         js = jsep*(1 - jw + jw/2*2)
         jj = 1 + (jw - 1)/2
         write(NAME(j), '(1A4)', ERR=77) rows(jj)(js+1: js+4)
         read(rows(jj)(js+6 : js+8), *, ERR=77) IB(j)
         read(rows(jj)(js+10: js+15), *, ERR=77) SCALE(j)
         read(rows(jj)(js+17: js+22), *, ERR=77) YSHIFT(j)
         read(rows(jj)(js+24: js+29), *, ERR=77) XL(j)
         if  (rows(jj)(js+31: js+36) /= ABNUM) then
            read(rows(jj)(js+31:js+36), *, ERR=77) XR(j)
         endif
         jn = max(jn, jj)
         j1 = j1 + 1
      endif
   enddo
endif

do j = 1, NB
   IBOX(j) = IB(j)
enddo

return

77 continue
write(*, *)
write(*, '(A)') '>>> INPUT ERROR encountered in the dialog window "Presentation"'
write(*, '(A, A52, A1)') '                Line: "', rows(j)(1:52), '"'
write(*, '(A52, $)') '     Enter "Y" to return, any other key to ignore > '
YKEY = 'X'
read(*, '(:, A1)')YKEY
if (YKEY == 'Y' .or. YKEY == 'y') goto 10

return
end subroutine ASXWIN

!---------------------------------------------------------------------=|
subroutine ASKXGR(JCHAN, IBOX, NAME, YMODE, JNB, JXMODE, JGR, OUTFIG, OUTNAME )

use parameter_inc, only: NRW
use outcmn_inc, only: MOD10, IP1, IP2, IP30, IP31, null_ch
use const_inc, only: XOUT

implicit none

integer, parameter :: NRW16=NRW+16, NRW96=NRW-128

integer, intent(in) :: JCHAN, YMODE, IBOX(*)
integer, intent(out) :: OUTFIG(*), JNB, JGR, JXMODE
character(len=4), intent(in)  :: NAME(*)
character(len=8), intent(out) :: OUTNAME(*)

integer :: JDONE, ASKGRF, JMODE, j, j1, jj, jn, jm, js, jb, jc, jw, jsep
character(len=80) :: rows(NRW16), STR, TITLE
character(len=1) :: YKEY

jdone = 1
TITLE = "Select profiles for plotting" // null_ch
!              ----5----0----5----0----5----0----5----0----5----0----5
STR = "Box #|  Name  |Marker||" // "Box #|  Name  |Marker"  // null_ch
jsep  = index(STR, '||') + 1
JMODE = MOD10
JXMODE = XOUT + 0.49

SELECT CASE(JMODE)
CASE(1)
   JGR = 8
CASE(2, 3)
   JGR = 4
CASE(6)
   if (YMODE == 1) JGR = 4
   if (YMODE == 0) JGR = 4
   if (YMODE  == -1) JGR = 2
   JXMODE = -1
CASE DEFAULT
   return
END SELECT

 10 continue

do j = 1, NRW
   OUTNAME(j) = '        '
   write(rows(j)(1: 80), '(79X, 1A1)') null_ch
enddo
JNB = JCHAN
if (JMODE  /=  1) JNB = min(jchan, 96)

! j  - ordinal box No.
! jb - box No. in the Astra nominations
! jm - number of empty boxes
! jw - position in the table
! jj - horizontal row
! js - position in the current row
jn = 0
jm = 0
do j=1, JNB
   jb = IBOX(j)
   if (jb <= 0 ) then
      jm = jm + 1
      CYCLE
   endif
   if (JMODE == 1) jw = IP1(jb)
   if (JMODE == 2 .or. JMODE == 3) jw = IP2(jb)
   if (JMODE == 6) then
      if (YMODE ==  1) jw = jb
      if (YMODE ==  0) jw = IP30(jb)
      if (YMODE == -1) jw = IP31(jb)
   endif
   js = jsep*(1 - jw + jw/2*2)
   jj = 1 + (jw - 1)/2
   write(rows(jj)(js+1 : js+5), '(1I4, 1X)') jb
   write(rows(jj)(js+7 : js+14), '(2X, 1A4, 2X)') NAME(j)
   write(rows(jj)(js+16: js+21), '(1I4, 2X)') 0
   jn = max(jn, jj)
enddo

if (jm /= 0 ) then
   jw = 2*jn + 1
   do j = 1, JNB
      if (IBOX(j) <= 0 ) then
         js = jsep*(1 - jw + jw/2*2)
         jj = 1 + (jw - 1)/2
         write(rows(jj)(js+1 : js+5 ), '(1I4, 1X)') IBOX(j)
         write(rows(jj)(js+7 : js+14), '(2X, 1A4, 2X)') NAME(j)
         write(rows(jj)(js+16: js+21), '(1I4, 2X)') 0
         jn = max(jn, jj)
         jw = jw + 1
      endif
   enddo
endif

j = 1
do while(j > 0)
   j1 = (jm + 1)/2  ! jm number of switched off windows
   j = ASKGRF(TITLE, STR, rows, 80, jn, JGR, j1, JXMODE)
enddo

jn = 0
do j=1, JNB
   if (IBOX(j) > 0) then
      if (JMODE == 1) jw = IP1(IBOX(j))
      if (JMODE == 2 .or. JMODE == 3) jw = IP2(IBOX(j))
      if (JMODE == 6) then
         if (YMODE ==  1) jw = IBOX(j)
         if (YMODE ==  0) jw = IP30(IBOX(j))
         if (YMODE == -1) jw = IP31(IBOX(j))
      endif
      js = jsep*(1 - jw + jw/2*2)
      jj = 1 + (jw - 1)/2
      write(OUTNAME(j), '(1A8)', ERR=77) rows(jj)(js+7: js+14)
      if (OUTNAME(j) /= '        ') then
         do while (OUTNAME(j)(1: 1) == ' ')
            OUTNAME(j)(1:) = OUTNAME(j)(2:) // '       '
         enddo
      endif
      j1 = index(rows(jj)(js+16: js+21), '.')
      if (j1 >= 2) then
         read(rows(jj)(js+16: js+14+j1), *, ERR=77) OUTFIG(j)
         read(rows(jj)(js+16+j1: js+21), *, ERR=77) jc
      elseif (j1 == 0) then
         read(rows(jj)(js+16: js+21), *, ERR=77) OUTFIG(j)
         jc = 0
      else
         goto 77
      endif
      if (OUTFIG(j) > 999) goto 78
      OUTFIG(j) = 1000*OUTFIG(j) + jc
      jn = max(jn, jj)
    endif
enddo

if (jm /= 0 ) then
jw = 2*jn + 1
   do j=1, JNB
      if (IBOX(j) <= 0) then
         js = jsep*(1 - jw + jw/2*2)
         jj = 1 + (jw - 1)/2
         write(OUTNAME(j), '(1A8)', ERR=77) rows(jj)(js+7: js+14)
         if (OUTNAME(j) /= '        ') then
            do while (OUTNAME(j)(1:1) == ' ')
               OUTNAME(j)(1:) = OUTNAME(j)(2:) // '       '
            enddo
         endif
         j1 = index(rows(jj)(js+16: js+21), '.')
         if (j1 >= 2) then
            read(rows(jj)(js+16: js+14+j1), *, ERR=77) OUTFIG(j)
            read(rows(jj)(js+16+j1: js+21), *, ERR=77) jc
         elseif (j1 == 0) then
            read(rows(jj)(js+16: js+21), *, ERR=77) OUTFIG(j)
            jc = 0
         else
            goto 77
         endif
         if (OUTFIG(j) > 999) goto 78
         OUTFIG(j) = 1000*OUTFIG(j) + jc
         jn = max(jn, jj)
         jw = jw + 1
      endif
   enddo
endif

jgr = 0
do j=1, JNB
   jb = OUTFIG(j)/1000
   jgr = max(jgr, jb)
enddo
if (jgr == 0) return

if (jxmode < -1 .or. jxmode > 5) then
   write(*, *)"Unknown data type"
   return
endif
return

 77 continue
write(*, *)
write(*, '(A)') '>>> INPUT ERROR encountered in the dialog window "Output control"'
write(*, '(A23, A52, A1)') '                Line: "', rows(j)(1:52), '"'
write(*, '(A52, $)') '     Enter "Y" to return, any other key to ignore > '
YKEY = 'X'
read(*, '(:, A1)') YKEY
if (YKEY == 'Y' .or. YKEY == 'y') goto 10
return

 78 continue
write(*, *)
write(*, '(A)') '>>> INPUT ERROR encountered in the dialog window "Output control"'
write(*, '(A)') '    Fig.# cannot exceed 999'
return

end subroutine ASKXGR

!---------------------------------------------------------------------=|
subroutine WRFIGS(JNB, JXMODE, JGR, OUTFIG, OUTNAME, IBOX, ITIMES, TTOUT, TOUT)
!----------------------------------------------------------------------|
! The subroutine writes file in the directories AWD/out/ and AWD/xmg/ 
!     with curves selected in the dialog window ASKXGR (hot key "O")
!----------------------------------------------------------------------|

use parameter_inc, only: NRW
use outcmn_inc, only: exp_file, equ_file, AWD, RUNID, &
   LTOUT, ROUT, XWH, WarningColor, null_ch
use const_inc, only: NAB, NA1, ABC, ROC, VOLUME
use status_inc, only: AMETR, RHO, FP, VOLUM

implicit none

integer, parameter :: JN0=0

integer, intent(in) :: ITIMES, JNB, JXMODE, JGR, OUTFIG(*), IBOX(*)
double precision, intent(in) :: TTOUT(ITIMES), TOUT(ITIMES, NRW)
character(len=8), intent(in) :: OUTNAME(*)

logical :: EXI
integer :: JLR, j, jj, jr, jn, jm, js, jb, jc, jw, jl, &
   NP1, JNUM(NRW), j0, j1, j2, j3, j4, j5, ios
double precision :: YY
character(len=1) :: CVE
character(len=4) :: ext
character(len=132) :: FNAME, xmg_dir, out_dir, fig_name, STRI

! Create output file name and open logical channel 7:
j = LEN_TRIM(exp_file)

jr = index(RUNID, exp_file(1:j)) + j

xmg_dir = TRIM(AWD) // 'xmg/' // TRIM(equ_file) // '.' // TRIM(exp_file)
out_dir = TRIM(AWD) // 'out/'
call system('mkdir -p ' // TRIM(xmg_dir))
call system('mkdir -p ' // TRIM(out_dir))

FNAME = TRIM(out_dir) // TRIM(exp_file) // '.' // TRIM(equ_file)
call set_filename(FNAME)
call colovm(WarningColor)

STRI = ' >>> Dataset is written in the file: ' // TRIM(FNAME)
write(*, '(/, A)') TRIM(STRI)
JLR = XWH - 125
call textvm(JN0, JLR, STRI, LEN_TRIM(STRI))
open(7, file=TRIM(FNAME), iostat=ios)
if(ios /= 0) then
   write(*, *) '>>> WRFIGS: Data file error'
   stop
endif
fig_name = TRIM(xmg_dir) // '/Fig'

j2 = LEN_TRIM(fig_name)

SELECT CASE(jxmode)
CASE(0)
   STRI(1: 10) = " a, m     "
   NP1 = NAB
CASE(1)
   STRI(1: 10) = " a_N, d/l "
   NP1 = NA1
CASE(2)
   STRI(1: 10) = "rho_N, d/l"
   NP1 = NA1
CASE(3)
   STRI(1: 10) = "  Psi, Vs "
   NP1 = NA1
CASE(4)
   STRI(1: 10) = " rho_V, m "
   NP1 = NA1
CASE(5)
   STRI(1: 10) = "rho_pol, m"
   NP1 = NA1
CASE(-1)
   STRI(1: 10) = " time, s  "
   NP1 = LTOUT-1
END SELECT

! The next block selects all the curves with the same number
! and submits them to a Figure with this number.
! The figures are put in ascending order. The curve numbers are ignored.
jl = 0     ! jl - total amount of curves
j3 = 1
do jj=1, jgr   ! 999 - max Fig #
   jc = 0    ! Amount of curves in one Fig
   do j1=1, JNB   ! NRW - max channel #
      jb = OUTFIG(j1)/1000  ! Fig #
      if (jb ==  0) CYCLE
      if (jb /= jj) CYCLE
      if (jc == 0) then
         FNAME = fig_name ! FNAME = .../Fig
         j5 = 0    ! Check existance
         EXI = .True.
         do while(EXI)
            j5 = j5+1
            write(ext, '(A1, 1I3)') '.', j5
            FNAME = TRIM(fig_name) // TRIM(ext)
            inquire(FILE=TRIM(FNAME) // '.dat', EXIST=EXI)
         enddo
         write(*, '(/3A)')' >>> New Figure: "', FNAME(1:jw), '"'
         open(8, file=TRIM(FNAME)//'.dat', iostat=ios)
         if (ios /= 0) goto 97
         open(9, file=TRIM(FNAME)//'.par', iostat=ios)
         if (ios /= 0) goto 98
         if (index(TRIM(AWD), 'efda-itm') == 0) write(9, 107)
         write(9, 106) FNAME(j2-2:jw), RUNID(1:jr), '"' ! Run & Fig ID
         SELECT CASE(jxmode)
         CASE(0)
            write(9, 101) "a [m]"
         CASE(1)
            write(9, 101) "a_N"
         CASE(2)
            write(9, 101) "rho_N"
         CASE(3)
            write(9, 101) "Psi [Vs]"
         CASE(4)
            write(9, 101) "rho_V [m]"
         CASE(5)
            write(9, 101) "rho_pol [m]"
         CASE(-1)
            write(9, 101) "time [s]"
         END SELECT
         write(7, '(2A)') 'Abscissa: ', STRI(1:10)
      endif
      STRI(11+10*jc:) = OUTNAME(j1) // '  '
      if (jc < 10) then
         write(9, 102) jc, OUTNAME(j1), jc
      else 
         write(9, 103) jc, OUTNAME(j1), jc
      endif
      jc = jc + 1
      CVE = char(96+jc)
      write(*, '(3A, I3, A, I3, 4A, I3)') &
        'Name  "', OUTNAME(j1) &
       , '"     Box', IBOX(j1) &
       , '"     Fig', jj, '(', CVE, ')' &
       , ',     Curve', jc
      write(7, '(2A, 3(A, I3))') &
        'Name  ', OUTNAME(j1) &
       , ' Box', IBOX(j1) &
       , ' Fig', jj &
       , ' Curve', jc
      jl = jl + 1
      JNUM(jl) = IBOX(j1)
   enddo

   if (jc /= 0) then
      j4 = jc + 1
      do j0=1, NP1
         SELECT CASE(jxmode)
         CASE(0)
            YY = AMETR(j0)
         CASE(1)
            YY = AMETR(j0)/ABC
         CASE(2)
            YY = RHO(j0)/ROC
         CASE(3)
            YY = (FP(j0) - FP(1))/(FP(NA1) - FP(1))
         CASE(4)
            YY = sqrt(VOLUM(j0)/VOLUME)
         CASE(5)
            YY = sqrt((FP(j0) - FP(1))/(FP(NA1) - FP(1)))
         CASE(-1)
            YY = TTOUT(j0)
         END SELECT
         j3 = jl - jc + 1
         if (jxmode >= 0) then
            write(8, '(1P, 11E12.4)') YY, (ROUT(j0, ibox(jnum(j))), j=jl-jc+1, jl)
         else
            write(8, '(1P, 11E15.7)') YY, (TOUT(j0, ibox(jnum(j))), j=jl-jc+1, jl)
         endif
      enddo
      close(8)
      close(9)
   endif
enddo

! return  ! If enabled suppress writing file to out/ (unit 7)
goto 13

! This block submits the curves to a Figure according to their numbering.
! Both the figures and curves within a figure are put in ascending order
! in exact correspondence with the pre-defined numbers.
! The curves with repeated numbers are ignored.
j1 = 0
jl = 0
do jj=1, jgr            ! 999 - max Fig #
   jw = -1
   do js=0, 19          ! 19  - max curve #
      do j=1, JNB       ! NRW - max channel #
         jb = OUTFIG(j)/1000   ! Fig #
         if (jb == jj .and. jj /= j1) j1 = jj
         CVE = '`'
         do jm=1, JNB
            jn = OUTFIG(jm)/1000      ! Retrieve Fig #
            if (jn /= jj) CYCLE
            CVE = char(ichar(CVE) + 1)
            jc = OUTFIG(jm) - jn*1000  ! Curve #
            if (js /= jc) CYCLE        ! Ignore curve #
            if (js == jw) CYCLE        ! Skip repeated curves
            CVE = char(97+jc)
            write(*, *) &
               'Name  "', OUTNAME(jm) &
             , '",    Box', IBOX(jm) &
             , ' ,    Curve', jc+1 &
             , '      Fig', jj, '(', CVE, ')'
            write(7, '(2A, 3(A, I3))') &
               'Name  ', OUTNAME(jm) &
             , ' Box', IBOX(jm) &
             , ' Fig', jj &
             , ' Curve', jc+1
            jw = js
            jl = jl+1
            JNUM(jl) = IBOX(jm)
         enddo    ! jm
      enddo    ! j
   enddo    ! js
enddo    ! jj

13 continue

jj = 0
js = 1

do j1=1, NP1
   SELECT CASE(jxmode)
   CASE(0)
      YY = AMETR(j1)
   CASE(1)
      YY = AMETR(j1)/ABC
   CASE(2)
      YY = RHO(j1)/ROC
   CASE(3)
      YY = FP(j1)
   CASE(4)
      YY = sqrt(VOLUM(j1)/VOLUME)
   CASE(5)
      YY = sqrt((FP(j1) - FP(1))/(FP(NA1) - FP(1)))
   CASE(-1)
      YY = TTOUT(j1)
   END SELECT
   if (jxmode >= 0) then
      write(7, 100) YY, (ROUT(j1, ibox(jnum(j))), j=js, jl)
   else
      write(7, 100) YY, (TOUT(j1, ibox(jnum(j))), j=js, jl)
   endif
enddo

close(7)
return

 100 format(1P, 6(E11.3))
 101 format( '    xaxis  label "', A, '"'/ &
     '    xaxis  label char size 1.75'/ &
     '    xaxis  ticklabel char size 1.5'/ &
     '    yaxis  label char size 1.75'/ &
     '    yaxis  ticklabel char size 1.5'/ &
     '    legend char size 1.75')
 102 format( '    s', I1, ' legend  "', A, '"'/ &
     '    s', I1, ' line linewidth 2.')
 103 format( '    s', I2, ' legend  "', A, '"'/ &
     '    s', I2, ' line linewidth 2.')
 104 format('    title "View file:   ', A, '"')
 105 format('    title size 1.000000')
 106 format('    subtitle "', 3A)
 107 format('page size 421, 298')
!     . /'with string'/'    string on'
!     . /'string loctype view'/'string 0.35, 0.5'/'string color 14'
!     . /'string rot 45'/'string font 0'/'string just 0'
!     . /'string char size 1.0'/'string def "Any text"'

 97 write(*, *) 'Xmgrace file "', FNAME(1: j1), '.dat" creation error'
return
 98 write(*, *) 'Xmgrace file "', FNAME(1: j1), '.par" creation error'
return

end subroutine WRFIGS

!---------------------------------------------------------------------==
!   The subroutine re-arranges real array A(1:N)
!-----------------------------------------------------------------------
function SORT(arr_in, nlen)

implicit none

integer, intent(in) :: nlen
double precision, dimension(nlen), intent(in) :: arr_in

logical, dimension(nlen) :: mask
integer :: j
double precision, dimension(nlen) :: SORT

if (nlen <= 1) return
mask = .True.

do j=1, nlen
   SORT(j) = MINVAL(arr_in, mask)
   mask(MINLOC(arr_in, mask)) = .False.
enddo

end function SORT

!---------------------------------------------------------------------==
!   The subroutine re-arranges real arrays A(1:N) and B(1:N)
!-----------------------------------------------------------------------
subroutine SORTAB(A, B, nlen)

implicit none

integer, intent(in) :: nlen
double precision, intent(inout), dimension(nlen) :: A, B

integer :: j, jj, jl, jr
double precision :: YA, YB

if (nlen  .LE.  1) return

do jl=1, nlen
   jr = jl
   YA = A(jl)
   YB = B(jl)
   do jj=jl+1, nlen
      if (A(jj) < YA) then
         YA = A(jj)
         YB = B(jj)
         jr = jj
      endif
   enddo
   if (jr /= jl) then
      do j=jr-1, jl, -1
         A(j+1) = A(j)
         B(j+1) = B(j)
      enddo
      A(jl) = YA
      B(jl) = YB
   endif
enddo

return
end subroutine SORTAB

!---------------------------------------------------------------------==
subroutine curvvm(id, npnts, array)
! The same as drcurv but without the factor 10 
! The chain: PLOTGR(obsolete) -> CURV(obsolete) -> CURVVM
!   -> drawvm(PSADrawLine) is not used any more
! dimension array(2*npnts)

implicit none

integer, intent(in) :: npnts, array(*), id

if (npnts == 1) call drawvm(id, array(1), array(2), array(1), array(2))
if (npnts <= 1) return

call drawline(id, array, npnts)

return
end subroutine curvvm

!---------------------------------------------------------------------==
subroutine drcurv(id, npnts, array)
! The same as curvvm but the supplied integer array is multiplied 
!     by the factor 10 in order to enhance PS resolution
!     This factor is then removed in C function d1line
!       PLOTCR -> CURV1 -> DRCURV -> d1line
! dimension array(2*npnts)

implicit none

integer, intent(in) :: npnts, array(*), id

if (npnts == 1) call d1line(id, array(1), array(2), array(1), array(2))
if (npnts <= 1) return

call d1polyline(id, array, npnts)

return
end subroutine drcurv

!---------------------------------------------------------------------==
character*6 function VARNAM(str_in, ierr)
!-----------------------------------------------------------------------
! If 1st character is tab or space, a blank string is returned
! If tabs are present anywhere else, they are replaced by blanks 
!    and ierr is set to 1
! Eventual trailing 'X' is removed
! Otherwise, VARNAM is returned.
!-----------------------------------------------------------------------

use outcmn_inc, only: esc_ch, tab_ch

implicit none

integer, intent(out) :: ierr
character(len=*), intent(in) :: str_in

integer :: j, nlen
character(len=1) :: symb

ierr = 0
VARNAM = str_in(1:6)
nlen = LEN_TRIM(VARNAM)

symb = str_in(1: 1)
if (symb == ' ' .or. symb == tab_ch) then
!   Ignore names starting with spaces and tabs
   VARNAM = '      '
   return
endif

do j=2, 6
   symb = str_in(j: j)
   if(symb == tab_ch .or. symb == esc_ch) then
      ierr = 1
      VARNAM(j: j) = ' '
   endif
enddo

if (VARNAM(nlen: nlen) == 'X') VARNAM(nlen: nlen) = ' '

return
end function VARNAM

!---------------------------------------------------------------------==
character*6 function ARRNAM(str_in)
!-----------------------------------------------------------------------
! The subroutine analizes a character*6 "string"
! If the 1st position is tab or space the string 6*' ' is returned
! If tabs are encountered on the end of the "string", 
!  they are removed the "string" is appended with spaces
! Trailing "X" is added when not present in string*6 
! Finally ARRNAM in the Astra standard is created, 
!-----------------------------------------------------------------------

use char_manip, only: clean_string, to_upper

implicit none

character(len=*), intent(in) :: str_in

integer :: nlen
character(len=len(str_in)) :: strtmp

strtmp = to_upper(str_in)
strtmp = clean_string(strtmp)
!call clean_string(strtmp, strtmp)
nlen = LEN_TRIM(strtmp)

ARRNAM = strtmp(1: 6)
! Append "X" if absent

if (nlen < 6 .and. strtmp(nlen: nlen) /= 'X') then
   ARRNAM(nlen+1: nlen+1) = 'X'
endif

return
end function ARRNAM

!---------------------------------------------------------------------==
subroutine CHECKU(INTYPE, ABC, AB, XBDRY, YX, jrad, jbdry, STRING, FILENA)
!----------------------------------------------------------------------|
! Consistency check for grid array YX(1:jrad) and plasma boundary AB/ABC
! Input:
! ABC - 
! AB - 
! YX(1:jrad) - array for a "radial" coordinate 
! jrad - YX array dimensionality
! STRING - U-file "Independent variable" description
! FILENA - U-file name
! Analyse array YX and returns proper values for XBDRY and jbdry
! Output:
! INTYPE - 
! XBDRY = ABC or AB depending on INTYPE selected
! jbdry - is determined from {YX(jbdry) <= ABC} or {YX(jbdry) <= AB}
!    for {INTYPE = 10} or {INTYPE = 11}, respectively
!  jbdry = jrad if  if YX(jrad) < ABC <= AB
!----------------------------------------------------------------------|

use char_manip, only: to_upper, clean_string

implicit none

integer, intent(in)  ::  jrad
integer, intent(out) :: jbdry, INTYPE
real*4 , intent(in)  :: YX(jrad)
double precision, intent(in)  :: ABC, AB
double precision, intent(out) :: XBDRY
character(len=*), intent(in)   :: FILENA
character(len=*), intent(inout) :: STRING

integer :: j
character(len=12) :: STRAD

jbdry = 0
XBDRY = YX(jrad)
STRAD = STRING(21:31)
STRING = to_upper(clean_string(STRING))
if (STRING(1:8) == 'MINORRAD') then
   INTYPE = 10
elseif (STRING(1:8) == 'MAJORRAD') then
   INTYPE = 19
elseif (STRING(1:3) == 'RHO') then
   INTYPE = 12
elseif (STRING(1:12) == 'POLOIDALFLUX') then
   INTYPE = 13
else
   write(*, *) '>>> U-file "', TRIM(FILENA), '"', &
         '    Unrecognized "radial" variable. Input ignored.' &
        // '    Allowed options are:' &
        // '  Minor Radius        m' &
        // '  Major Radius        m' &
        // '  Rho Toroidal, normalized' &
        // '  Poloidal Flux, normalized'
   INTYPE = -1
   return
endif
STRAD = to_upper(clean_string(STRAD))
if ((INTYPE == 10 .or. INTYPE == 19) .and. STRAD(1: 1) /= 'M') then
   write(*, *) ">>> Warning: Inconsistency in the input data"
   write(*, *) '>>> U-file "', TRIM(FILENA), '": ', &
         '    radial grid is expected to be given in "m"'
endif

if (INTYPE == 10 .and. abs(XBDRY - AB) > 0.3*AB/jrad)  then
   if(XBDRY < AB) then
      jbdry = jrad
   else
      do j=jrad, 1, -1
         if (YX(j) > AB) jbdry = j
      enddo
      XBDRY = AB
   endif
elseif (INTYPE == 11 .and. abs(XBDRY - ABC) > 0.3*ABC/jrad) then
   if(XBDRY < ABC) then
      jbdry = jrad
   else
      do j=jrad, 1, -1
        if (YX(j) > ABC) jbdry = j
      enddo
      XBDRY = ABC
   endif
endif
if (INTYPE == 18) then
   write(*, *)'>>> U-file "', TRIM(FILENA), '"', &
        " Don't know a distance to the major axis.", &
            '           Set to RTOR'
endif
if (jbdry == 0) then
   jbdry = jrad
endif

return
end subroutine CHECKU

!---------------------------------------------------------------------=|
!  Subroutine minimizes the value of functional
!  INTEGRAL(alfa*P(x)*(dU/dx)**2+(U-F)**2)*dx, 
!  where FO(NO) is a function, given on the grid XOld(NOld)
! P(x) is equal to unit now
! ALFA=alfa<<0.01*XO(NO)**2 is regularizator
! NO - number of old grid points
! N=<NRD - number of new grid points
! 0<=XO(NO) - old grid |     both grids are arbitrary
! 0<=XN(N)  - new grid |     but XO(NO)=XN(N)
! FO(NO) - origin function, given on the grid XO(NO)
! FN(N) - smoothed function on grid XN(N)
!  The result is function FN(XN), given on the new grid
! with additional conditions:
! dFN/dx(x=0)=0 - cylindrical case and
! FN(XN(N))=FO(XO(NO))
!----------------------------------------------------------------------|
subroutine SMOOTH(ALFA, NO, FO, XO, N, FN, XN)
!----------------------------------------------------------------------|

use parameter_inc, only: NRD
use debugger, only: astra_stop

implicit none

integer, intent(in) :: NO, N
double precision, intent(in) :: ALFA, XO(*), FO(*), XN(*)
double precision, intent(out) :: FN(*)

integer :: I, J
double precision :: YF, YX, YP, YQ, YD, FJ, P(NRD)
character(len=132) :: err_msg

if (N > NRD .or. NO <= 0) then
   write(err_msg, '(A, 1X, i, 1X, i)') ' >>> SMOOTH: array is out of limits', NRD, NO
   call err_catch_a
   call astra_stop(err_msg)
endif
if (NO == 1) then
   do j=1, N
      FN(j) = FO(1)
   enddo
   return
endif
if (NO == 2) then
   do j=1, N
      FN(j) = (FO(2)*(XN(j) - XO(1)) - FO(1)*(XN(j) - XO(2)))/(XO(2) - XO(1))
   enddo
   return
endif
if (N < 2) then
   call err_catch_a
   call astra_stop(' >>> SMOOTH: no output grid is provided')
endif
if (abs(XO(NO) - XN(N)) > XN(N)/N) then
   write(*, *) '>>> SMOOTH: grids are not aligned'
   write(*, '(1A23, I4, F8.4)')'     Old grid size/edge', NO, XO(NO)
   write(*, '(1A23, I4, F8.4)')'     New grid size/edge', N, XN(N)
   call err_catch_a
   call astra_stop
endif
do j=2, N
   YP = XN(j) - XN(j-1)
   if (YP <= 0.d0) then
      write(*, *)'>>> SMOOTH: new grid is not increasing monotonically'
      write(*, '(A, I4, A, F8.4)')'Node ', j-1, '   Value', XN(j-1)
      write(*, '(A, I4, A, F8.4)')'Node ', j,  '   Value', XN(j)
      call err_catch_a
      call astra_stop
   endif
   P(j) = ALFA/YP/XO(NO)**2
enddo
P(1) = 0.
FN(1) = 0.
I = 1
YF = (FO(2) - FO(1))/(XO(2) - XO(1))
YX = 2./(XN(2) + XN(1))
YP = 0.
YQ = 0.
do j=1, N-1
   if (XO(I) <= XN(j)) then
      do
         I = I + 1
         if (I > NO) I = NO
         if (I == NO .or. XO(I) >= XN(j)) EXIT
      enddo
      YF = (FO(I) - FO(I-1))/(XO(I) - XO(I-1))
   endif
   FJ = FO(I) + YF*(XN(j) - XO(I))
   YD = 1. + YX*(YP + P(j+1))
   P(j) = YX*P(j+1)/YD
   FN(j) = (FJ + YX*YQ)/YD
   if (j /= N-1) then
      YX = 2./(XN(j+2) - XN(j))
      YP = (1. - P(j))*P(j+1)
      YQ = FN(j)*P(j+1)
   endif
enddo

FN(N) = FO(NO)
do j=N-1, 1, -1
   FN(j) = P(j)*FN(j+1) + FN(j)
enddo
end subroutine SMOOTH

!---------------------------------------------------------------------=|
subroutine SMAP(ALFA, NO, XO, N, XN, F)
!----------------------------------------------------------------------|
! Smooth mapping from grid 
! Similar to SMOOTH but the same array, F, is used for input and output
!----------------------------------------------------------------------|

use parameter_inc, only: NRD

implicit none

integer, intent(in) :: NO, N
double precision, intent(in) :: ALFA, XO(*), XN(*)
double precision, intent(inout) :: F(*)

integer :: J
double precision :: P(NRD)

call SMOOTH(ALFA, NO, F, XO, N, P, XN)
do j=1, N
   F(j) = P(j)
enddo
end subroutine SMAP

!---------------------------------------------------------------------=|
subroutine TAMP(ALFA, X1, X2, NO, FO, XO, N, FN, XN)
!----------------------------------------------------------------------|
! Same as SMOOTH but P=P(ALFA)
!----------------------------------------------------------------------|

use parameter_inc, only: NRD
use debugger, only: astra_stop

implicit none

integer, intent(in) :: NO, N
double precision, intent(in) :: ALFA, X1, X2, XO(*), FO(*), XN(*)
double precision, intent(out) :: FN(*)

integer :: I, J
double precision :: P(NRD), YF, YX, YP, YQ, YD, FJ

if (N > NRD .or. NO <= 0) then
   call astra_stop(' >>> TAMP: array is out of limits')
endif
if (NO == 1) then
   do j=1, N
      FN(j) = FO(1)
   enddo
   return
endif
if (NO == 2) then
   do j=1, N
      FN(j) = (FO(2)*(XN(j) - XO(1)) - FO(1)*(XN(j) - XO(2)))/(XO(2) - XO(1))
   enddo
   return
endif
if (N < 2) then
   call astra_stop(' >>> TAMP: no output grid is provided')
endif
if (abs(XO(NO) - XN(N)) > XN(N)/N) then
   write(*, *) ' >>> TAMP: input grids are not aligned'
   write(*, '(1A23, I4, F8.4)') '     Old grid size/edge', NO, XO(NO)
   write(*, '(1A23, I4, F8.4)') '     New grid size/edge', N, XN(N)
   call astra_stop
endif
do j=2, N
   P(j) = ALFA/(XN(j) - XN(j-1))/XO(NO)**2
   P(j) = P(j)*exp(-((XN(j) - X1)/X2)**2)
enddo
P(1)  = 0.
FN(1) = 0.
I =1
YF = (FO(2) - FO(1))/(XO(2) - XO(1))
YX =2./(XN(2) + XN(1))
YP = 0.
YQ = 0.
do j=1, N-1
   if(XO(I) <= XN(j)) then
      do
         I = I + 1
         if (I > NO) I = NO
         if (I == NO .or. XO(I) >= XN(j)) EXIT
      enddo
      YF = (FO(I) - FO(I-1))/(XO(I) - XO(I-1))
   endif
   FJ = FO(I) + YF*(XN(j) - XO(I))
   YD = 1. + YX*(YP + P(j+1))
   P(j) = YX*P(j+1)/YD
   FN(j) = (FJ + YX*YQ)/YD
   if (j /= N - 1) then
      YX = 2./(XN(j+2) - XN(j))
      YP = (1. - P(j))*P(j+1)
      YQ = FN(j)*P(j+1)
   endif
enddo

FN(N) = FO(NO)
do j=N-1, 1, -1
   FN(j) = P(j)*FN(j+1) + FN(j)
enddo

return
end subroutine TAMP

!---------------------------------------------------------------------=|
double precision function RA2Z(RTOR, YR, YA, YS, YE, YT)
!---------------------------------------------------------09-FEB-97
!  Z [m] as a function of major and minor radii Z(R, a)
!   above the plasma midplane 
! YS = SHIF(A), YE=ELON(A), YT=TRIA(A), YV=SHIV(A)
! NB !  Up-down shift not allowed (G.P.)
!---------------------------------------------------------==Polevoi

implicit none

double precision, intent(in) :: RTOR, YR, YA, YS, YE, YT

double precision :: YF, YF1, YT2, YSIN2

YF = (YR - RTOR - YS)/YA
YF1 = 2.*YF*YT
YT2 = (1. + YF1)
if(YF >= 1. .or. YF < -1. .or. YT2 <= 0.) then
   RA2Z = 0.
else 
   YF1 = 2.*YF*YT
   YSIN2 = (1. - YF)*(1. + YF)/YT2
   if(YSIN2 > .001 .and. YT > .01) then
      YT2 = 2.*YT*YT
      YSIN2 = (SQRT(1. + 2.*(YF1 + YT2)) - (1. + YF1))/YT2
   endif
   RA2Z = YA*YE*SQRT(YSIN2)
endif

return
end function RA2Z

!---------------------------------------------------------------------=|
double precision function RZ2A(YR, YZ, N)
!-----------------------------------------------------------------------
! The function returns A(r, z) where
!
! r = R0 + del(A) + A*[cos(theta)-tri(A)*sin^2(theta)]
! z = UPD + A*elo(A)*sin(theta)
!
! and del(j), elo(j), tri(j) are given as arrays[1:N] 
!      on the grid A=AMETR(j)
!-----------------------------------------------------------------------

use const_inc, only: RTOR, AB
use status_inc, only: AMETR, SHIF, SHIV, ELON, TRIA

implicit none

integer, intent(in) :: N
double precision, intent(in) :: YR, YZ

integer :: j, j1
double precision :: Y1, YAS, YAO, YHOR, YVER, YELO, YTRI, YA, QUADIN

YAS = ((YZ - SHIV(N))/ELON(N))**2
YA = sqrt(YAS + (YR - RTOR - SHIF(N))**2)
do j=1, 20
   YAO = YA
   YHOR = QUADIN(N, AMETR, SHIF, YA, YAS, j1)
   YELO = QUADIN(N, AMETR, ELON, YA, YAS, j1)
   YTRI = QUADIN(N, AMETR, TRIA, YA, YAS, j1)
   YVER = QUADIN(N, AMETR, SHIV, YA, YAS, j1)
   YAS = ((YZ - YVER)/YELO)**2
   Y1 = YAS
   if (YA > 1.E-4) Y1 = YAS/YA
   YA = sqrt(YAS + (YR - RTOR - YHOR + YTRI*Y1)**2)
   j1 = j
   if (abs(YA-YAO) < 1.E-10) EXIT
enddo
RZ2A = min(YA, AB)

end function RZ2A

!---------------------------------------------------------------------==
double precision function QUADIN(NG, XGRID, FUNC, X, PQD, jj)
!-----------------------------------------------------------------------
!   Quadratic interpolation of a function FUNC(NG) 
!     given on a grid XGRID(NG) to a position X
! Input  NG, XGRID(NG), FUNC(NG), X
! Output QUADIN, PQD, jj
!-----------------------------------------------------------------------

implicit none

integer, intent(in) :: NG
integer, intent(out) :: jj
double precision, intent(in) :: XGRID(NG), FUNC(NG), X
double precision, intent(out) :: PQD

integer :: j
double precision :: YF1, YF2, YF3, Y, YY, YX, &
    YD21, YD23, YD31, YDX1, YDX2, YDX3

! No extrapolation
YX = max(X, 0.d0)
YX = min(YX, XGRID(NG))
YY = 0.
jj = 1
do j=1, NG
   Y = XGRID(j)-YX
   jj = j
   if (Y < 0) then
      YY = Y
   else if (Y == 0) then
      EXIT
   else
      if (Y > -YY) jj = jj - 1
      EXIT
   endif
enddo

jj = max(jj, 2)
jj = min(jj, NG-1)

YD21 = XGRID(jj)   - XGRID(jj-1)
YD23 = XGRID(jj)   - XGRID(jj+1)
YD31 = XGRID(jj+1) - XGRID(jj-1)
YDX1 = YX - XGRID(jj-1)
YDX2 = YX - XGRID(jj)
YDX3 = YX - XGRID(jj+1)
if (YD21 <= 0. .or. YD23 >= 0.) then
   QUADIN = FUNC(jj)
   PQD    = 0.
else
   YF1 = FUNC(jj-1)/YD21/YD31
   YF2 = FUNC(jj)/YD21/YD23
   YF3 = FUNC(jj+1)/YD31/YD23
   QUADIN = YF1* YDX2 * YDX3  + YF2* YDX1 * YDX3  - YF3* YDX1 * YDX2
   PQD    = YF1*(YDX2 + YDX3) + YF2*(YDX1 + YDX3) - YF3*(YDX1 + YDX2)
endif

return
end function QUADIN

!---------------------------------------------------------------------==
double precision function RECTAN(NG, XGRID, FUNC, X)
!-----------------------------------------------------------------------
! Step-function FUNC(X) is defined
! Input  NG      grid size
!  XGRID(1:NG) grid
!  FUNC(1:NG)  grid function
!  X     function argument
! Output RECTAN      function value
!  jj the grid point nearest to X (not used)
!-----------------------------------------------------------------------

implicit none

integer, intent(in) :: NG
double precision, intent(in) :: XGRID(NG), FUNC(NG), X

integer :: j, jj
double precision :: Y, YY, YX

! No extrapolation
YX = max(X, 0.d0)
YX = min(YX, XGRID(NG))

YY = 0.
jj = 1
do j=1, NG
   Y = XGRID(j)-YX
   jj = j
   if (Y < 0) then
      YY = Y
   else if (Y == 0) then
      EXIT
   else
      if (Y > -YY) jj = jj - 1
      EXIT
   endif
enddo

jj = max(jj, 1)
jj = min(jj, NG)

RECTAN = FUNC(jj)

return
end function RECTAN

!---------------------------------------------------------------------=|
subroutine get_runid
!-----------------------------------------------------------------------
! The subroutine forms string RUNID and additionally returns 
! date and time when those are not defined (calling from INIT)
!-----------------------------------------------------------------------

use outcmn_inc, only: RUNID, equ_file, exp_file, VERSION

implicit none

integer :: time_arr(8), YEAR, MONTH, DAY, HOUR, MINUTE
integer :: j
character(len=3) :: vers
character(len=15) :: datetime

call date_and_time(VALUES=time_arr)

YEAR   = time_arr(1)
MONTH  = time_arr(2)
DAY    = time_arr(3)
HOUR   = time_arr(5)
MINUTE = time_arr(6)
write(datetime, "(1I2, 2('-', 1I2.2), 1I3, ':', 1I2.2)") &
    DAY, MONTH, YEAR-2000, HOUR, MINUTE

j = index(VERSION, 'Version')
vers = version(j+8: j+10)

RUNID = "ASTRA " // vers // " -- " // datetime // ' -- Model: ' // TRIM(equ_file) &
   // ' -- Data: ' // TRIM(exp_file)

return
end subroutine get_runid

!---------------------------------------------------------------------=|
subroutine makemovie(PNMNAME)
!----------------------------------------------------------------------|

use outcmn_inc, only: XWH, Magenta, null_ch
use debugger, only: markloc

implicit none

character, intent(inout) :: PNMNAME*(*)

call markloc('makemovie')

if (PNMNAME(1:1) /= "*") then
   call colorb(Magenta)
   call textbf(20, XWH - 124, "Saving movie in tmp/movie.mpeg  ", 32)
   call redraw(0)
   call system(".srv/esc2movie " // TRIM(PNMNAME) )
   call system("rm tmp/*.ppm") ! Done in esc2movie
   call textbf(20, XWH - 124, "                                ", 32)
   PNMNAME = '*'
endif

end subroutine makemovie

!---------------------------------------------------------------------=|
subroutine dump_movie_frame(jframe, PNMNAME)

use outcmn_inc, only: XWH, Magenta, null_ch
use debugger, only: markloc

implicit none

integer, intent(inout) :: jframe
character, intent(inout) :: PNMNAME*(*)

character(len=32) :: STRI
integer :: jj

call markloc('dump_movie_frame')

if (PNMNAME(1:1) /= "*") then
   jj = 1
   if (jframe > 9)   jj = 2
   if (jframe > 99)  jj = 3
   if (jframe > 999) jj = 4
   write(STRI, '(A, i0)') "Making Movie:  Frame # ", jframe
   call colorb(Magenta)
   call textbf(20, XWH - 124, TRIM(STRI), 23+jj)
! xwd -silent -id XW_ID | xwdtopnm > frame_name.pnm) >& /dev/null
   call redraw(0)
   call px2pnm(PNMNAME) ! dump an image of X window (man "xwd") 
   jframe = jframe + 1
endif

return
end subroutine dump_movie_frame

!---------------------------------------------------------------------=|
subroutine getnames_(YEQUNAME, YEXPNAME)

character(len=*), intent(out) :: YEQUNAME, YEXPNAME

call getnames(YEQUNAME, YEXPNAME)

return
end subroutine getnames_

!---------------------------------------------------------------------=|
subroutine getnames(YEQUNAME, YEXPNAME)
!----------------------------------------------------------------------|
! The subroutine can be called from C function, returns equ_file, exp_file
! G.V.Pereverzev 16.02.2004
!----------------------------------------------------------------------|

use outcmn_inc, only: equ_file, exp_file, null_ch

implicit none

character(len=*), intent(out) :: YEQUNAME, YEXPNAME

YEQUNAME = TRIM(equ_file) // null_ch
YEXPNAME = TRIM(exp_file) // null_ch

return
end subroutine getnames

!---------------------------------------------------------------------=|
!Efable GRID2GRID computes quantity on shifted grid from main grid
!
!  grid_type: 1 - main to shift, 2 shift to main
!  x_input: x_variable
!  y_input: y_variable
!  y_output
!  nagrid: number of grid points
! iextrap: 1 if yes interpolate last grid point, 0 do not interpolate last grid point
!---------------------------------------------------------------------=|
subroutine GRID2GRID(grid_type, x_input, y_input, y_output, nagrid, iextrap)

implicit none

integer, intent(in) :: nagrid, iextrap, grid_type
double precision, intent(in) , dimension(nagrid) :: x_input, y_input
double precision, intent(out), dimension(nagrid) :: y_output

if (grid_type == 1) then
   call  MAIN2SHIFT(y_input, y_output, nagrid)
endif

if (grid_type == 2) then
   call SHIFT2MAIN(x_input, y_input, y_output, nagrid)
endif

return
end subroutine GRID2GRID

!---------------------------------------------------------------------=|
!Efable MAIN2SHIFT computes quantity on shifted grid from main grid
!
!  x_input: x_variable
!  y_input: y_variable
!  y_output
!  nagrid: number of grid points
!---------------------------------------------------------------------=|
subroutine MAIN2SHIFT(y_input, y_output, nagrid)

implicit none

integer, intent(in) :: nagrid
double precision, intent(in) , dimension(nagrid) :: y_input
double precision, intent(out), dimension(nagrid) :: y_output

integer :: j

! Normalized grid , GRP style
do j=1, nagrid-1
   y_output(j) = 0.5*(y_input(j+1) + y_input(j))
enddo

y_output(nagrid) = 0.5*(3.0*y_input(nagrid) - y_input(nagrid-1))

return
end subroutine MAIN2SHIFT

!---------------------------------------------------------------------=|
!Efable SHIFT2MAIN computes quantity on main grid from shifted grid
!
!  x_input: x_variable
!  y_input: y_variable
!  y_output
!  nagrid: number of grid points

!---------------------------------------------------------------------=|
subroutine SHIFT2MAIN(x_input, y_input, y_output, nagrid)

use numerical_tools, only: polyfitcc

implicit none

integer, intent(in) :: nagrid
double precision, intent(in) , dimension(nagrid) :: x_input, y_input
double precision, intent(out), dimension(nagrid) :: y_output

integer :: j
double precision :: y1tmp, P(3)
  
! Normalized grid , GRP style
do j=2, nagrid
   y_output(j) = 0.5*(y_input(j) + y_input(j-1))
enddo

call polyfitcc(x_input(1: 3), y_input(1: 3), P)
y1tmp = P(3)
y_output(1) = 0.5*(y1tmp + y_input(1))

return
end subroutine SHIFT2MAIN

