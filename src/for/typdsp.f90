subroutine ERASXY()

use outcmn_inc, only: Xwin_height, Xwin_width, DXLET, White

implicit none

integer :: JN0, JLR
character(len=118) :: STRI

JLR = Xwin_height - 125
JN0 = 0
STRI = repeat(' ', 80)

call colovm(White)
call textvm(JN0, JLR, STRI(1: 15), 15)
call textvm(Xwin_width - 83*DXLET + 2, JLR, STRI(1: 80), 80)

end subroutine ERASXY

!---------------------------------------------------------------------
subroutine TYPDSP(YN, ITIMES, TTOUT, TOUT)
! NCH= 5 - terminal, 0 - file (old format), 1 - file (new format)

use parameter_inc, only: NRW
use const_inc, only: XOUT, NAB, NA1, NA, XLINE1, IPL, BTOR, RTOR, ABC, &
    TIME, TAU, CONSTF
use status_inc, only: MU, AMETR, RHO, FP
use outcmn_inc, only: NCFNAM, LTOUT, NTOUT, NROUT, MOD10, RUNID, AWD, &
    equ_file, exp_file, NAMER, NAMET, WarningColor, ROUT, Xwin_height
use dbl2char, only: fmt4, fmt_xf

implicit none

integer, parameter :: NLINSC=50
character(len=40), parameter :: STRMN=' R=     a=     B=     I=     q=     <n>='
character(len=6), dimension(6), parameter :: &
    CONN = (/ ' CF   ', ' CV   ', ' CH   ', ' CCD  ', ' CBND ', ' CRAD ' /)
integer, intent(in) :: ITIMES
double precision, intent(in) :: YN, TTOUT(ITIMES), TOUT(ITIMES, NRW)

integer :: NP1, ITBE, ITEND, ITEN, MODEX, &
    JBE, JEND, J, JEN, JJ, J1, JLR, ios
double precision :: YQ
character(len=6) :: CH6
character(len=118) :: STRI
character(len=132) :: FNAME, dat_dir

! NLINSC - maximum line number 
MODEX = XOUT + 0.49
! MODEX = 0 [0, AB]  against "a"
! MODEX = 1 [0, ABC] against "a"
! MODEX = 2 [0, ROC] against "rho"
! MODEX = 3 [FP(1), FP(NA1)] against "psi"
! otherwise Unknown option => MODEX=0
NP1 = NAB
if (MODEX >= 1 .and. MODEX <= 3 .or. MOD10 == 3) NP1 = NA1

if (MOD10 <= 5) then
    JBE  = 1
    JEND = 16
    do
        JEN = MIN0(NROUT, JEND)
        write(STRI, '(16(1X, 1A4))') (NAMER(J), J=JBE, JEN)
        write(*, '(1X, A)') TRIM(STRI)
        do J=1, NP1
            STRI = ' '
            do JJ=JBE, JEN
                J1 = 5*(JJ - JBE + 1) - 4
                STRI(J1: J1+4) = fmt_xf(ROUT(J, JJ), 4)
            enddo
            write(*, '(1X, A)') TRIM(STRI)
        enddo
        if (JEN == NROUT) return
        JBE  = JEN + 1
        JEND = JEN + 16
    enddo 

endif

if (MOD10 <= 7) then
    JBE  = 1
    JEND = 15
    do
        JEN = MIN(NTOUT, JEND)
        ITBE  = 1
        ITEND = NLINSC
        do
            ITEN = MIN(LTOUT-1, ITEND)
            write(STRI, 308) (NAMET(J), J=JBE, JEN)
            write(*, '(1X, A)') TRIM(STRI)

            do J1=ITBE, ITEN
                STRI = ' '
                STRI(1: 5) = fmt_xf(TTOUT(J1), 4)
                do J=JBE, JEN
                    JJ = 5*(J - JBE) + 6
                    STRI(JJ: JJ+4) = fmt_xf(TOUT(J1, J), 4)
                enddo
                write(*, '(1X, A)') TRIM(STRI)
            enddo
            if (ITEN == LTOUT - 1) EXIT
            ITBE  = ITEN
            ITEND = ITEN + NLINSC - 1
        enddo
        if (JEN == NTOUT) EXIT
        JBE  = JEN + 1
        JEND = JEN + 15
    enddo
endif

101 format(1X, 1A6, 1A111)
102 format('   Time', 16(3X, 1A4))
104 format(1X, 1A120)
308 format(1X, 'Time', 15(1X, 1A4))
408 format(1PE12.3, 64(1PE12.3))

return
end subroutine TYPDSP

!---------------------------------------------------------------------
integer function GETIME(TIME, TIMES, NNOUT)
! The function returns
!  if NNOUT=1  then GETIME=1
!  otherwise
!     GETIME = an index of the array TIMES(1:NNOUT) element >= TIME

implicit none

integer, intent(in) :: NNOUT
double precision, intent(in) :: TIMES(*)
double precision, intent(inout) :: TIME

integer :: j

if (NNOUT <= 1) then
    GETIME = 1
    return
endif

if (TIME <= TIMES(1)) then
    GETIME = 1
    TIME = TIMES(1)
    return
endif

do j=2, NNOUT
    if (TIMES(j) >= TIME) then
        GETIME = j
        return
    endif
enddo

GETIME = NNOUT
TIME = TIMES(NNOUT)

return
end function GETIME

!---------------------------------------------------------------------
subroutine PUTXY(IX, IY, ITIMES, TTOUT, TOUT)
!---------------------------------------------------------------------
! Input: IX
!  IY
!  MODEY
!---------------------------------------------------------------------

use parameter_inc, only: NRW
use outcmn_inc, only: Xwin_height, Xwin_width, MOD10, IY0, IYM, &
    scale_bnd, canv_hei, canv_wid, resizeGraph, &
    IDT, IDX, MODEY, DXLET, DYLET, frame_hei, LTOUT, NTOUT, active_tab, &
    NWIND3, NAMET, White, Red, Blue, nx_canvas, ny_canvas
use status_inc, only: AMETR, SHIF, ELON, TRIA, FP, RHO
use const_inc, only: TIME, TINIT, TSCALE, NA, NA1, NAB, XOUT, AB, ABC, ROC, HRO
use dbl2char, only: fmt5
use numerical_tools, only: QUADIN

implicit none

integer, intent(in) :: ITIMES, IX, IY
double precision, intent(in) :: TTOUT(ITIMES), TOUT(ITIMES, NRW)

integer :: IM, JX, JY, JLR, j, j1, JC, JL, IX0, IXM, MODEX, &
    GETIME, JW, JN0, JN2
double precision :: DX, DY, YX, YX1, YY, YY1, YA, YA1, YD, YE, YT, &
    YRHO, YFP, YFPC, RZ2A
character(len=80) :: STRI

JN0 = 0
JLR = Xwin_height - int(125*resizeGraph)
! DLINER_ theGCA, 0, Xwin_height-128, Xwin_width-1, Xwin_height-128);
! DLINER_ theGCA, 0, Xwin_height-110, Xwin_width-1, Xwin_height-110);
! DLINER_ theGCA, 0, Xwin_height-109, Xwin_width-1, Xwin_height-109);
if (MOD10 <= 0) return
if (MOD10 == 7) call NEGA(IM, ITIMES, TOUT)
call set_frame(IM, IX0, IXM)
STRI = repeat(' ', 80)
STRI(7:25) = '(x, y)=(     ,     )'
JX = IX - 10
JY = IY - 10
if (IX0 > JX .or. JX > IXM .or. IY0 > JY .or. JY > IYM) then
    STRI = repeat(' ', 80)
    call colovm(White)
    call textvm(JN0, JLR, "               ", 15)
    call textvm(Xwin_width - 83*DXLET + 2, JLR, STRI(1: 80), 80)
    return
endif

if (MOD10 == 7) nx_canvas = 2
DX = 1./nx_canvas
DY = 1./ny_canvas
YX1 =      (JX - IX0 + 0.)/(IXM - IX0)
YY1 = 1. - (JY - IY0 + 0.)/(IYM - IY0)

if (MOD10 == 6) then
! (window_width)/(step=IDX=23)/(n_labels)=592/23/25=1.0295652
    YX = TINIT + 1.029565*YX1*abs(TSCALE)
    do j=1, ny_canvas
        YY1 = YY1 - DY
        if (YY1 < 0) EXIT
    enddo
    YY = (YY1 + DY)/DY
    JN2 = frame_hei + 15 + 3*DYLET
    JC = 0
    JL = 0

    YY1 = max(TIME, TTOUT(LTOUT-1), TTOUT(LTOUT))
    YY1 = min(YX, YY1)
    YY1 = max(TTOUT(1), YY1)
    j = GETIME(YX, TTOUT, LTOUT)
    do J1=1, NTOUT
        JW = NWIND3(J1) - 8*active_tab(MOD10)
        if (NAMET(J1) == '    ') JW = 0
        if (JW > 0 .and. JW <= 8) call down_label(j, ITIMES, TOUT) ! for the all modes
    enddo

    call colovm(Red)
    STRI(1 :  5) = 'Time='
    STRI(6 : 10) = fmt5(YY1)
    STRI(11: 11) = 's'
    call textvm(DXLET, JN2 - 3*DYLET + DYLET/2, STRI, 11)
    return
else if (MOD10 == 8) then
!    YX = 5.*YX1*scale_bnd
    YX = YX1*scale_bnd*canv_wid/IDT/IDX
    YY = (YY1 - 0.5)*scale_bnd*canv_hei/IDT/IDX
    YA1= RZ2A(YX, YY, NAB)
    YD = QUADIN(NAB, AMETR, SHIF, YA1)
    YE = QUADIN(NAB, AMETR, ELON, YA1)
    YT = QUADIN(NAB, AMETR, TRIA, YA1)
    STRI( 7: 12) = "(r, z)="
    STRI(32: 37) = "(a, S)="
    STRI(38: 50) = '(     ,     )'
    STRI(39: 43) = fmt5(YA1)
    STRI(45: 49) = fmt5(YD)
    STRI(57: 62) = "(E, T)="
    STRI(63: 75) = '(     ,     )'
    STRI(64: 68) = fmt5(YE)
    STRI(70: 74) = fmt5(YT)
else
    do j=1, nx_canvas
        YX1 = YX1 - DX
        if (YX1 < 0) EXIT
    enddo
    YX = (YX1 + DX)/DX
    if (MOD10 < 5) then
        YA = 0.
        YFP = 0.
        YRHO = 0.
        YFPC = 1.125*FP(1) - 0.125*FP(2)
        MODEX = XOUT + 0.49
        if (MOD10 == 3) MODEX = 3
        if (MOD10 == 4) MODEX = 0

        SELECT CASE(modex)
        CASE(0)
            YA = YX*AB
            YRHO = QUADIN(NA1, AMETR, RHO, YA)
            YFP  = QUADIN(NA1, AMETR, FP , YA)
        CASE(1)
            YA = YX*ABC
            YRHO = QUADIN(NA1, AMETR, RHO, YA)
            YFP  = QUADIN(NA1, AMETR, FP , YA)
        CASE(2)
            YRHO = YX*ROC
            YA   = QUADIN(NA1, RHO, AMETR, YRHO)
            YFP  = QUADIN(NA1, RHO, FP   , YRHO)
        CASE(3)
            YFP  = YFPC + (FP(NA1) - YFPC)*YX
            YA   = QUADIN(NA1, FP, AMETR, YFP)
            YRHO = QUADIN(NA1, FP, RHO  , YFP)
        END SELECT

        j = YRHO/HRO + 1
        if (YRHO > 0.5*(RHO(NA) + ROC)) j = NA1

        STRI(32: 33) = "a="
        STRI(34: 38) = fmt5(YA)
        YRHO = YRHO/ROC
        STRI(39: 49) = 'm,   rho_t='
        STRI(50: 54) = fmt5(YRHO)
        if (YFP > YFPC) then
            YFP = sqrt((YFP - YFPC)/(FP(NA1) - YFPC))
        else
            YFP = 0.
        endif
        STRI(55: 64) = ",   rho_p="
        STRI(65: 69) = fmt5(YFP)
        STRI(70: 77) = ",  Node:"
        write(STRI(78: 80), '(1I3)')j
    endif

    do j=1, ny_canvas
        YY1 = YY1 - DY
        if (YY1 < 0) EXIT
    enddo
    YY = (YY1 + DY)/DY
    if (MOD10 > 1 .and. MOD10 < 6 .and. MODEY == -1) then
        if (JY - IY0 > (IYM - IY0)/ny_canvas) YY = YY - 1.
    endif

endif

STRI(14: 18) = fmt5(YX)
STRI(20: 24) = fmt5(YY)
call colovm(Blue)

call textvm(JN0, JLR, "               ", 15)
call textvm(Xwin_width - 83*DXLET + 2, JLR, STRI(1: 80), 80)

return
end subroutine putxy

!---------------------------------------------------------------------
subroutine set_frame(plot_mode, canv_x_left, canv_x_right)
!----------------------------------------------------------------------|
! Input: MODEY
! Output: canv_x_left - x_left  of the graphic area
!  canv_x_right - x_right of the graphic area
!  plot_mode - for using in set_plot
!  NST - for using in set_plot
!----------------------------------------------------------------------|

use outcmn_inc, only: MOD10, IY0, IYM, canv_wid, canv_hei, MODEY, NST, &
    DXLET, DYLET, frame_wid, frame_hei, nx_canvas, ny_canvas

implicit none

integer, intent(inout) :: plot_mode
integer, intent(out) :: canv_x_left, canv_x_right

integer :: STRUP, STRDN

! IY0, IYM     - upper & lower plot boundaries
! canv_wid, canv_hei     - X & Y window dimension
! STRUP, STRDN - number of text strings up & down
! nx_canvas, ny_canvas       - number of windows in horisontal & vertical axises
! canv_x_left, canv_x_right    - left & right plot boundaries

if (NST >= 1 .or. plot_mode == 11 .or. plot_mode == 12 .or. plot_mode == 21 .or. plot_mode == 22) then

    STRUP = 2
    STRDN = 4
    NST = NST + 1
    canv_x_left = 0
    canv_x_right = frame_wid
    nx_canvas = 1
    ny_canvas = 1
    if (NST == 1) THEN
        canv_x_right = frame_wid/2.
        if (plot_mode == 21 .or. plot_mode == 22) THEN
            nx_canvas = 2
            ny_canvas = 2
        endif
    endif
    if (NST == 2) THEN
        canv_x_left = frame_wid/2.
        if (plot_mode == 12 .or. plot_mode == 22)THEN
            nx_canvas = 2
            ny_canvas = 2
        endif
    endif
    IYM = 320

else

    if (MOD10 <= 1 .or. MOD10 >= 7) then
        plot_mode = MOD10
    elseif (MOD10 >= 2 .and. MOD10 <= 5) then
        if (MODEY == 1)  plot_mode = 2
        if (MODEY == -1) plot_mode = 3
    else
        if (MODEY == 1)  plot_mode = 5
        if (MODEY == 0)  plot_mode = 9
        if (MODEY == -1) plot_mode = 6
    endif

    if (plot_mode <= 0 .or. plot_mode >= 10) then
        plot_mode = 0
        return
    endif

    SELECT CASE(plot_mode)

    CASE(1)
        STRUP = 2
        STRDN = 5
        canv_x_left = 0
        canv_x_right = frame_wid
        nx_canvas = 4
        ny_canvas = 2

! modes 2, 3, 4, 5 at y-mode = +1, dummy mode (plot_mode=4 - not used)
    CASE(2, 4, 7)
        STRUP = 2
        STRDN = 5
        canv_x_left = 0
        canv_x_right = frame_wid
        nx_canvas = 2
        ny_canvas = 1

! modes 2, 3, 4, 5 at y-mode = -1
    CASE(3)
        STRUP = 2
        STRDN = 5
        canv_x_left = 0
        canv_x_right = frame_wid
        nx_canvas = 2
        ny_canvas = 2

! mode # 6 (time) at y-mode=1 (2 windows)
    CASE(5)
        STRUP = 1
        STRDN = 1
        canv_x_left = 6*DXLET
        canv_x_right = frame_wid
        nx_canvas = 1
        ny_canvas = 2

! mode # 6 (time) at y-mode=-1 (4 windows)
    CASE(6)
        STRUP = 1
        STRDN = 1
        canv_x_left = 6*DXLET
        canv_x_right = frame_wid
        nx_canvas = 1
        ny_canvas = 4

! mode 8 (equilibrium)
    CASE(8)
        STRUP = 1
        STRDN = -2
        canv_x_right = 0.7*frame_wid
        canv_x_left = 0
        nx_canvas = 1
        ny_canvas = 1

! mode # 6 (time) at y-mode=0 (1 window), mode 9 (user's plot)
    CASE(9)
        STRUP = 1
        STRDN = 1
        canv_x_left = 6*DXLET
        canv_x_right = frame_wid
        nx_canvas = 1
        ny_canvas = 1

    END SELECT

endif

IY0 = STRUP*DYLET + 1
if (MOD10 == 6) IY0 = IY0 + 1
IYM = frame_hei - STRDN*DYLET - 1
canv_wid = (canv_x_right - canv_x_left)/nx_canvas
canv_hei = (IYM - IY0)/ny_canvas

return
end subroutine set_frame

!---------------------------------------------------------------------
subroutine set_plot(plot_mode, canv_x_left, canv_x_right)
! Subroutine draw frame for different modes

use const_inc, only: TSCALE, TINIT, ABC
use outcmn_inc, only: Black, MOD10, KPRI, &
    Xwin_width, Xwin_height, IY0, IYM, IDX, IDT, DXLET, DYLET, scale_bnd, &
    frame_hei, frame_wid, canv_hei, canv_wid, pixel_ymid, meter2pixel
use dbl2char, only: fmt_xf
use char_manip, only: len_trim_tab

implicit none

integer, intent(in) :: plot_mode, canv_x_left, canv_x_right

integer :: JJ, J, JN0, TIMWIN, JX, JY, &
    LENG, XP, XM, YP, YM, JXSCM, LYM
double precision :: DY, YY, TIND, scale_fac
character(len=5) :: XF4
character(len=6) :: CH6
character(len=80) :: COMMENT

data LENG/3/ JN0/0/

! IY0, IYM     - upper & lower grafic boundary
! canv_wid, canv_hei     - X & Y window dimension
! IDX, DY      - X & Y distance (in points) between X & Y axis labels
! IDT       - distance (in labels) between longer labels in modes 6&8
! LENG        - label length

if (plot_mode <= 0) return

TIMWIN = 0
if (plot_mode == 5 .or. plot_mode == 6 .or. plot_mode == 8 .or. plot_mode == 9)  TIMWIN = 1

call colovm(Black)
call rectvm(0, JN0, JN0, Xwin_width-1, Xwin_height-1)

if (MOD10 == 6)  then
    call colovm(Black)
    j = Xwin_width - 20*DXLET + 1
    jj = frame_hei + DYLET
    call textvm(j, jj, 'time, s', 7)
endif

! Vertical lines & Y-labels
if (KPRI >= 1 .and. KPRI <= 2) then
    write(COMMENT, '(A)') "Vertical lines"
    j = len_trim_tab(COMMENT)
    call pscom(COMMENT, j)
endif
call colovm(Black)

! Skipping from a subplot to the next along x-axis
do JJ=canv_x_left, canv_x_right, canv_wid
    JX = MIN0(canv_x_right, JJ)
    XP = MIN(canv_x_right , JX + LENG)
    XM = MAX(canv_x_left  , JX - LENG)
    call drawvm(0, JX, IY0, JX, IYM)
! Y-line labels
    if (KPRI >= 1 .and. KPRI <= 2) then
        write(COMMENT, '(A)') "Y-line labels"
        j = len_trim_tab(COMMENT)
        call pscom(COMMENT, j)
    endif
    DY = (IYM - IY0)/20.
    YY = dble(IYM)
    do
        JY = YY
        if (plot_mode /= 8) call drawvm(0, XM, JY, XP, JY)
        YY = YY - DY
        if (sign(1.d0, DY)*(YY - dble(IY0)) < 0.d0) EXIT
    enddo
enddo

! Horizontal lines
if (KPRI >= 1 .and. KPRI <= 2) then
    write(COMMENT, '(A)') "Horizontal lines"
    j = len_trim_tab(COMMENT)
    call pscom(COMMENT, j)
endif

do JY=IYM, IY0, -canv_hei
    YM = MAX(IY0, JY - LENG)
    if (plot_mode == 4 .or. plot_mode == 5 .or. plot_mode == 6) then
        YP = JY
    else
        YP = MIN0(IYM, JY + LENG)
    endif
    call drawvm(0, canv_x_left, JY, canv_x_right, JY)
! X-line labels
    if (KPRI >= 1 .and. KPRI <= 2) then
        write(COMMENT, '(A)')"X-line labels"
        j = len_trim_tab(COMMENT)
        call pscom(COMMENT, j)
    endif
    JX = canv_x_left
    IDX = 16
    JXSCM = canv_x_right - IDX
    if (TIMWIN == 1) then
        IDX = 23
        JXSCM = canv_x_right
    endif
    do J=1, 100
        JX = JX + IDX
        if (JX > JXSCM) EXIT
        if (TIMWIN == 1 .and. J/IDT*IDT == J) then
            LYM = YM - 2
        else
            LYM = YM
        endif
        if (LYM > IY0 + LENG)  call drawvm(0, JX, LYM, JX, YP)
    enddo
enddo

if (plot_mode == 8) then
    JY = (IY0 + IYM)/2
    do j=0, 10
        jj = JY + IDX*j
        if (jj < IYM) call drawvm(0, XM, JJ, XP, JJ)
        jj = JY - IDX*j
        if (jj > IY0) call drawvm(0, XM, JJ, XP, JJ)
    enddo
endif

if (MOD10 == 6) then

! time-axis legend:
    call colovm(Black)
    JJ = frame_hei - DYLET + 12
    do J=0, frame_wid, IDX
        JX = (J - IDT)*IDT + canv_x_left
        if (JX > frame_wid) CYCLE
! (right_label_pos)/(n_labels)=575/IDT=115
        YY = abs(TSCALE)
        TIND = TINIT + J*YY/115
        if ( TINIT + YY > 10.0 .or. (TINIT + YY > 1.0 .and. YY < 0.1) .or. YY < 0.01) then
            CH6 = fmt_xf(TIND, 5)
            call textvm(JX - 2, JJ, CH6, 6)
        else
            XF4 = fmt_xf(TIND, 4)
            call textvm(JX, JJ, XF4, 5)
        endif
  enddo

endif

if (MOD10 == 8) then
! horizontal axis labels

    scale_bnd = nint(20.*ABC)/10.
    scale_fac = dble(frame_hei)/350.
    IDX = IDX*scale_fac
    JJ = IYM + DYLET + 2
    call colovm(Black)
    do J=1, 5
        JX = IDX*IDT*J - 24
        if (JX > JXSCM) CYCLE
        YY = J*scale_bnd
        XF4 = fmt_xf(YY, 4)
        call textvm(JX, JJ, XF4, 5)
    enddo
! vertical axis labels
    do J=-1, 1
        JX = (IYM + IY0 + DYLET)/2 + IDT*IDX*J - 0.5*DYLET
        YY = -J*scale_bnd
        XF4 = fmt_xf(YY, 4)
        call textvm(canv_x_right + 2, JX, XF4, 5)
    enddo
endif

if (KPRI >= 1 .and. KPRI <= 2) then
     write(COMMENT, '(A)') "Frame done"
     j = len_trim_tab(COMMENT)
     call pscom(COMMENT, j)
endif

pixel_ymid  = 0.5*(IY0 + IYM)
meter2pixel = dble(IDX*IDT)/scale_bnd

return
end subroutine set_plot

!---------------------------------------------------------------------
subroutine NEGA(J1, ITIMES, TOUT)

use parameter_inc, only: NRW
use ac_neg1, only: NKL1, NKL2, JMIN, JMAX, NUM, MODK

implicit none

integer, intent(in) :: ITIMES
integer, intent(out) :: J1
double precision, intent(in) :: TOUT(ITIMES, NRW)

integer :: IBEG, JJ, J

save IBEG
data IBEG/0/

IBEG = IBEG + 1
J1 = 11
if (IBEG < 3) return
do JJ=NKL1, NKL2
    if (JJ <= 0) CYCLE
    MODK(JJ) = 0
    do  J=JMIN, JMAX
        if (TOUT(J, NUM(2*JJ-1)) < 0.0 .or. TOUT(J, NUM(2*JJ)) < 0.0) then
            MODK(JJ) = 1
            EXIT
        endif
    enddo
enddo
! (0, 0) - 11, (-1, 0) - 21, (0, -1) - 12, (-1, -1) - 22
if (MODK(1) == 0) then
    if (MODK(2) == 0) then
        J1 = 11
    else
        J1 = 12
    endif
else
    if (MODK(2) == 0) then
        J1 = 21
    else
        J1 = 22
    endif
endif

return
end subroutine NEGA

!---------------------------------------------------------------------
subroutine down_label(jt, ITIMES, TOUT)
!----------------------------------------------------------------------|
! Time dependences for radial output
! Curve to digit conversion
!  Input:  jt defines the current time
!
!    if mod10 != 6 or call from run then jt = LTOUT  
!
!  Data used from outcmn.inc
! double precision TOUT(1, 1)
! integer LTOUT, frame_hei, DXLET, DYLET, MOD10, active_tab(*), NTOUT
! integer NWIND3(*)
! character*4 NAMET(*)
!----------------------------------------------------------------------|

use parameter_inc, only: NRW
use outcmn_inc, only: LTOUT, MOD10, DXLET, DYLET, frame_hei, NTOUT, &
    NWIND3, active_tab, NAMET, Black, Blue, curves_per_frame
use dbl2char, only: fmt_xf

implicit none

integer, parameter :: fshift=10

integer, intent(in) :: ITIMES
integer, intent(inout) :: jt
double precision, intent(in) :: TOUT(ITIMES, *)

integer :: JN2, JN0, JEND, JB, JL, JC, JW, JJ, J
character(len=5) :: XF4
character(len=7) :: XF7
character(len=80) :: STRI, STRIN

if (jt == 0) then
    jt = LTOUT
    call colovm(Black)
else
    call colovm(Blue)
endif

if (MOD10 == 6) then

    JN0 = 6*DXLET
    JN2 = frame_hei + FSHIFT + 3*DYLET + 5
    JC = 0
    JL = 0
    write(STRIN, '(79X, 1A1)') ' '
    write(STRI , '(79X, 1A1)') ' '
    do j=1, NTOUT
        JW = NWIND3(j) - curves_per_frame(MOD10)*active_tab(MOD10)
        if (NAMET(j) == '    ') JW = 0
        if (JW <= 0 .or. JW > curves_per_frame(MOD10)) CYCLE
        JC = JC + 1   ! Actual curve number in the mindow
        jj = 8*JC - 6
        if (jj > 74) CYCLE
        if (jj >= 66) JN0 = 5*DXLET
        if (jj == 74) JN0 = -DXLET
! curve #, win #, chan #, mode 6, screen #
        XF7 = fmt_xf(TOUT(jt, j), 6)
        STRI (jj: jj+6) = XF7
        STRIN(jj: jj+6) = '  ' // NAMET(J) // ' '
        JL = max(JL, jj + 6)
    enddo
    call textvm(JN0, JN2, STRI, JL)
    JN2 = JN2 - DYLET + 1
    call textvm(JN0, JN2, STRIN, JL)
else
    JB = 1
    JN0 = 0
    JN2 = frame_hei - 5*DYLET + FSHIFT + 2
    do
        JEND = MIN0(JB + 15, NTOUT)
        do J=JB, JEND
            XF4 = fmt_xf(TOUT(jt, J), 4)
            if (NAMET(J) == ' ') XF4 = '    '
            JJ = 5*(J - JB + 1) - 4
            STRI(JJ: JJ+4) = XF4
        enddo
        JN2 = JN2 + 2*DYLET + 2
        JL = 5*(JEND - JB + 1)
        call textvm(JN0, JN2, STRI, JL)
        write(STRI, '(16(1X, 1A4))') (NAMET(J), J=JB, JEND)
        JN2 = JN2 - DYLET + 1
        call textvm(JN0, JN2, STRI, JL)
        JN2 = JN2 + DYLET - 1
        if (JEND == NTOUT .or. JEND == NRW) EXIT
        JB = JB + 16
    enddo
endif

return
end subroutine down_label

!---------------------------------------------------------------------
! Upper string of a picture
subroutine up_label(YN, YQ)

use outcmn_inc, only: Xwin_height, Xwin_width, DXLET, DYLET, active_tab, MOD10, Black, Blue
use const_inc, only: RTOR, BTOR, IPL, ABC, XLINE1
use dbl2char, only: fmt40

implicit none

integer, parameter :: JN0=0, fshift=2
double precision, intent(in) :: YN, YQ

character(len=2) :: CHR
character(len=62) :: STRMN

STRMN = ' R=     a=     B=     I=     q=     n=     '

STRMN( 4:  7) = fmt40(RTOR)
STRMN(11: 14) = fmt40(ABC)
STRMN(18: 21) = fmt40(BTOR)
STRMN(25: 28) = fmt40(IPL)
STRMN(32: 35) = fmt40(YQ)
STRMN(39: 42) = fmt40(YN)

call colovm(Black)
call textvm(JN0, FSHIFT, XLINE1(1: 19) // STRMN, 61)
call textvm(JN0 + 55*DXLET, FSHIFT - DYLET + 4, '_', 1)
call colovm(Blue)
write(CHR, '(1I2)') active_tab(MOD10) + 1
call textvm(JN0 + 79*DXLET, FSHIFT + DYLET - 2, CHR, 2) ! Screen No.
call colovm(Black)
call rectvm(0, JN0, JN0, Xwin_width - 1, Xwin_height - 1)

return
end subroutine up_label

!---------------------------------------------------------------------
subroutine TIMEDT(TIME, DT)

use outcmn_inc, only: Black, DXLET
use dbl2char, only: fmt50

implicit none

integer, parameter :: FSHIFT=2
double precision, intent(in) :: TIME, DT

character(len=19) :: STRI

STRI(1 : 5 ) = 'Time='
STRI(11: 16) = ' dt='
STRI( 6: 10) = fmt50(TIME)
STRI(15: 19) = fmt50(DT)
call colovm(Black)
call textvm(62*DXLET, FSHIFT, STRI, 19)

return
end subroutine TIMEDT

!---------------------------------------------------------------------
subroutine const2ps
! Appending the list of constants to a PS file

use const_inc, only: CONSTF, DEVAR
use outcmn_inc, only: NCFNAM, NPRNAM, PRNAME, null_ch
use dbl2char, only: fmt_xf

implicit none

character(len=6), dimension(22), parameter :: CONN = (/ &
    'CF1-> ', 'CF5-> ', 'CF9-> ', 'CF13->', &
    'CV1-> ', 'CV5-> ', 'CV9-> ', 'CV13->', &
    'CHE   ', 'CHI   ', 'CNB   ', 'CNBI  ', &
    'CCD   ', 'CRF   ', 'CNEUT ', 'CPEL  ', &
    'CBND  ', 'CFUS  ', 'CIMP  ', 'CMHD  ', &
    'CRAD  ', 'CSOL  ' /)

integer :: J, J1, J2, JJ, JDUM, JNY, JNX
character(len=6 ) :: CH6
character(len=80) :: STRI

! Writing constants

JNY = 540
JNX = 10

ps_loop: do J2=1, 11

    JDUM = 2*(J2 - 1) + 1
    STRI = CONN(JDUM)

    do J=1, 4
        J1 = J1 + 1
        if (J1 > NCFNAM) EXIT ps_loop
        CH6 = fmt_xf(CONSTF(J1), 5)
        JJ = 7*(J - 1) + 8
    enddo

    STRI(JJ: JJ+5) = CH6
    STRI(JJ+12: JJ+18) = CONN(JDUM+1)

    do J=5, 8
        J1 = J1 + 1
        if (J1 > NCFNAM) EXIT ps_loop
        CH6 = fmt_xf(CONSTF(J1), 5)
        JJ = 7*(J - 1) + 20
    enddo

    STRI(JJ: JJ+5) = CH6
    JNY = JNY + 17

    call textvm(JNX, JNY, STRI, 75)

enddo ps_loop

! Writing variables
j1 = 1
JNY = JNY + 20
JDUM = NPRNAM - 48

do j=1, JDUM
    CH6 = fmt_xf(DEVAR(j), 5)
    STRI(j1: j1+19) = PRNAME(j) // '=' // CH6 // '     '
    j1 = j1+20
    if (j1 > 70 .or. j == JDUM) then
        STRI(j1-5:) = null_ch
        JNY = JNY + 17
        call textvm(JNX, JNY, STRI, j1-5)
        j1 = 1
    endif
enddo

return
end subroutine const2ps

!===================================================================
! DEVID - device ID e.g. "AUGD"
! ARUNID - Astra run ID (the lowest string of the graphic output)
! UFNAME - U-file name to write
! TIME - current time value
! SIGNAM - signal name
! NPOINT - number of points
! NASC - number of associated scalar quantities
! PROCOD - processing code
! ARGUM - argument values
! SIGNAL - function values
! RTOR, AB, BTOR, IPL - standard Astra notations
! YN - average density
!-----------------------------------------------------------------------
subroutine UF1DWA(DEVID, ARUNID, UFNAME, TIME, SIGNAM, NPOINT, NASC, &
    PROCOD, ARGUM, SIGNAL, RTOR, AB, BTOR, IPL, YN, MODEX)

implicit none

integer, intent(in) :: NPOINT, NASC, PROCOD, MODEX
double precision, intent(in) :: TIME, YN, IPL, RTOR, AB, BTOR
double precision, intent(in), dimension(NPOINT) :: ARGUM, SIGNAL
character(len=*), intent(in) :: UFNAME
character(len=4 ), intent(in) :: DEVID, SIGNAM
character(len=80), intent(in) :: ARUNID

integer :: j, NSHOT, NCHU

NSHOT = 0
NCHU = 1

open(NCHU, file='udb/' // TRIM(UFNAME), status='UNKNOWN')

write(NCHU, 1001) NSHOT, DEVID, 1, 0
write(NCHU, 1002)
write(NCHU, 1003) NASC
if (NASC == 0) write(NCHU, 1005)
if (NASC == 1) then
    write(NCHU, 1004) TIME
    write(NCHU, 1006)
    if (MODEX == 0 .or. MODEX == 1) write(NCHU, 1007)
    if (MODEX == 2) write(NCHU, 1010)
    if (MODEX == 3) write(NCHU, 1013)
endif
if (NASC > 1 .or. MODEX > 3) then
    write(NCHU, *) "Warning: don't know how to write 1D U-file"
    write(*, *)    "Warning: don't know how to write 1D U-file"
endif
write(NCHU, 1008) SIGNAM
write(NCHU, 1009) PROCOD
if (NASC == 0) write(NCHU, 1011) NPOINT
if (NASC == 1) write(NCHU, 1012) NPOINT
write(NCHU, '(1X, 1P, 6E13.5)') (ARGUM(J) , J=1, NPOINT)
write(NCHU, '(1X, 1P, 6E13.5)') (SIGNAL(J), J=1, NPOINT)
write(NCHU, *) ' ;----END-OF-DATA---------------COMMENTS:-----------'
write(NCHU, '(1A80)') ARUNID
write(NCHU, '(5(A, F6.2), A)') &
    ' R =', RTOR, 'm,   a =', AB, 'm,   B =', BTOR, &
    'T,   I =', IPL, 'MA,   <n_e> =', 0.1*YN, 'E20m^-3'
write(*, *) '>>>  Dataset "', SIGNAM, '" is written in the 1D ' // &
     'U-file "udb/' // TRIM(UFNAME), '"'
close(NCHU)

return

 1001 format(1I7.5, 1A4, 2I2, 16X, '; Shot #  ')
 1002 format(1X, 9HDD-MMM-YY, 21X, '; Shot date -  Ufiles ASCII file system')
 1003 format(1I4, 27X, '; # of associated scalars')
 1004 format(1P, 1E13.6, 18X, '; Scalar value, LABEL FOLLOWS:')
 1005 format(' Time:', 15X, 'sec', 7X, '; Independent variable: time')
 1006 format(' Time:', 15X, 'sec', 7X, ';')
 1007 format(' MINOR RAD', 11X, 'm', 9X, '; Independent variable: a')
 1010 format(' Rho toroidal, normalized', 6X, '; Independent variable: Rho')
 1013 format(' Poloidal flux, normalized', 5X, '; Independent variable: Psi')
 1008 format( ' Function            ??        ; Function name:   "', 1A4, '"')
 1009 format(1I2, 29X, '; Processing code: ..., 4 - Astra output')
 1011 format(1I11, 20X, ';-# of  time  pts  t, F(t) (data follow)')
 1012 format(1I11, 20X, ';-# of radial pts  X, F(X) (data follow)')

end subroutine UF1DWA

!---------------------------------------------------------------------
subroutine UF2DWA(DEVID, UFNAME, SIGNAM, JN, NRP, NASC, PROCOD, &
   PRMARK, TIMOD4, RTOR, AB, BTOR, IPL, YN, YWA, YWB, YWC)
!---------------------------------------------------------------------
! Note: 2D U-file is written with the radial number of points NRP
! as mapped from the full grid size NB1 because NA1 can vary in
! time and is not suitable for the U-file fixed grid
!  YWA(*), YWB(*), YWC(*) working arrays 
!---------------------------------------------------------------------

use parameter_inc, only: NRD
use outcmn_inc, only: IPOUT, rev_file, NXOUT, NGR, NROUT, RUNID

implicit none

integer, intent(in) :: NRP, NASC, PROCOD, JN
double precision, intent(in) :: PRMARK(*), TIMOD4(*), YN, IPL, RTOR, AB, BTOR
double precision, intent(out) :: YWA(*), YWB(*), YWC(*)
  
character(len=4 ), intent(in) :: DEVID, SIGNAM
character(len=40), intent(in) :: UFNAME

integer*2 :: JNT2(NRD)
integer :: NSHOT, NCHU, NCHR, JNT, j, j1, j2, NTP, SKIPM, JAB, I, jj, ios
double precision :: TEMPR, YTIME, SCL, DOWN, ALFA, SIGNAL(NRD)
character(len=4) :: CHAR4

! NTP - total number of time slices: 
NTP = 0
do j = 1, IPOUT-1
    if (PRMARK(j) >= 0) NTP = NTP + 1
enddo
NSHOT = 0
NCHU = 1
NCHR = 2

open(NCHU, file='udb/' // TRIM(UFNAME), status='UNKNOWN')
write(NCHU, 1001) NSHOT, DEVID, 2, 0
write(NCHU, 1002)
write(NCHU, 1003) NASC
if (NASC /= 0) then
    write(NCHU, *) "Warning: don't know how to write 2D U-file"
    write(*, *)    "Warning: don't know how to write 2D U-file"
    write(*, *)    "         associated scalar is not defined"
endif
! NOTE! a radially contiguous u-file is written
write(NCHU, 1007)
write(NCHU, 1005)
write(NCHU, 1008) SIGNAM
write(NCHU, 1009) PROCOD
write(NCHU, 1012) NRP
write(NCHU, 1011) NTP
do j=1, NRP
    YWC(j) = (j - 1.)/(NRP - 1.)
enddo
write(NCHU, '(1X, 1P, 6E13.5)') (YWC(J)*AB, J=1, NRP)

j1 = 0
j2 = 0
do j=1, NTP
    if (PRMARK(j) < 0) CYCLE
    j1 = j1+1
    j2 = j2+1
    if (j2 == 1) then
        write(NCHU, '(1X, 1P, E13.5, $)') TIMOD4(j1)
    elseif (j2 == 6) then
        write(NCHU, '(1P, E13.5)') TIMOD4(j1)
        j2 = 0
    else
        write(NCHU, '(1P, E13.5, $)') TIMOD4(j1)
    endif
enddo
if (j2 /= 0) write(NCHU, *)

open(NCHR, FILE=TRIM(rev_file), FORM='UNFORMATTED', iostat=ios)
if (ios /= 0) then
    write(*, *) '>>> UF2DWA: Output profile error'
    stop
endif
j = 0
j = SKIPM(NCHR, j)  ! Returned value is not used
read(NCHR, ERR=38, END=39) CHAR4
if (NXOUT > 0 .and. NGR > 0) read(NCHR, ERR=38, END=39) JNT

do J2 = 1, IPOUT-1
    read(NCHR, ERR=38, END=39) JNT
    if (JNT /= 0) read(NCHR, ERR=38) TEMPR
    read(NCHR, END=39) YTIME
!  Skip  CONSTF, DEVAR, LINEAV, ROC
    read(NCHR, END=39) TEMPR

    read(NCHR, END=39) JAB, JAB, (JNT, I=1, 10), (TEMPR, I=1, 10)
!  Retrieve  AMETR
    read(NCHR)SCL, DOWN, (JNT2(jj), jj=1, JAB)
    do j=1, JAB
        YWA(J) = (DOWN + SCL*(JNT2(J) + 32768)/65535.)/AB
    enddo
    YWA(JAB) = 1.
!  Skip  SHIF, ELON, TRIA, RHS1, RHS2, FP
    do j=1, 6
        read(NCHR) TEMPR, TEMPR, (JNT2(jj), jj=1, JAB)
    enddo
    do J1=1, NROUT
        read(NCHR, ERR=38) SCL, DOWN, (JNT2(jj), jj=1, JAB)
        if (j1 == JN .and. PRMARK(j2) >= 0) then
            do j=1, JAB
                SIGNAL(J) = DOWN + SCL*(JNT2(J) + 32768)/65535.
            enddo
            ALFA = .0001
            CALL SMOOTH(ALFA, JAB, SIGNAL, YWA, NRP, YWB, YWC)
            write(NCHU, '(1X, 1P, 6E13.5)') (YWB(J), J=1, NRP)
        endif
    enddo
enddo

 39 close(NCHR)

write(NCHU, *) ' ;----END-OF-DATA---------------COMMENTS:-----------'
write(NCHU, '(1A80)') RUNID
write(NCHU, '(1A4, 1F6.2, 3(1A8, 1F6.2), 1A13, 1F6.2, 1A7)') &
    ' R =', RTOR, 'm,   a =', AB, 'm,   B =', BTOR, &
    'T,   I =', IPL, 'MA,   <n_e> =', .1*YN, 'E20m^-3'
write(*, *) '>>>  Dataset "', SIGNAM, &
    '" is written in the 2D U-file "udb/', TRIM(UFNAME), '"'
close(NCHU)

return

 1001 format(1I7.5, 1A4, 2I2, 16X, '; Shot #  ')
 1002 format(1X, 9HDD-MMM-YY, 21X, '; Shot date -  Ufiles ASCII file system')
 1003 format(1I4, 27X, '; # of associated scalars')
 1005 format(' Time:', 15X, 'sec', 7X, '; Independent variable: time')
 1007 format(' MINOR RAD', 11X, 'm', 9X, '; Independent variable: Rho')
 1008 format( ' Function            ??        ; Function name:   "', 1A4, '"')
 1009 format(1I2, 29X, '; Processing code: 4 - Astra output')
 1011 format(1I11, 20X, ';-# of  time  pts  t, F(t) (data follow)')
 1012 format(1I11, 20X, ';-# of radial pts  X, F(X) (data follow)')

38 write(*, *) 'Read file PROFIL.DAT error'
stop

end subroutine UF2DWA

!---------------------------------------------------------------------
integer function SKIPM(NCHR, NCHL)
! Skip model & model.log records
! If NCHL =/= 0 then SCRATCH file containing model.log is writen
! Returns (a returned value is used in the postviewer only)
!   0 for the old format (version before 5.1)
!   1 for the new format (version 5.2 and later)

implicit none

integer, intent(in) :: NCHR
integer, intent(inout) :: NCHL

character(len=1) :: CH1
character(len=32) :: STR
character(len=132) :: STRI

integer :: j, n

read(NCHR) STR
if (STR /= "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^") then
    rewind(NCHR)
    SKIPM = 0
    return ! -> Old format file
endif
SKIPM = 1 ! -> New format file
j = 0
n = 1
STR = "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"

do
    read(NCHR) CH1, STRI(1: ichar(CH1))
    if (ichar(CH1) == 32 .and. STR == STRI(1:32)) then
        j = j + 1 ! No. of adjacent strings 32*"^"
        n = n + 1 ! Total No. of strings 32*"^" already encountered
        if (j == 2) then
            if (n <= 2) write(*, *) "Error in binary output file: wrong format"
            if (n <= 3) then
                if (NCHL /= 0) then
                    close(NCHL)
                    NCHL = 0
                endif
                return
            endif
            if (n <= 4) return
        endif
        if (n == 2 .and. NCHL /= 0) open(NCHL, STATUS='SCRATCH')
    else
        j = 0
    endif
    if (j == 0 .and. n == 2 .and. NCHL /= 0) write(NCHL, '(A)') STRI(1:ichar(CH1))
enddo

return
end function SKIPM

!---------------------------------------------------------------------
subroutine set_filename(FNAME)
! FNAME - input name (without blanks) is appended with an extension.
!   The extension is the ordinal number of the file

implicit none

character(len=*), intent(inout) :: FNAME

integer :: jext
logical :: EXI
character(len=4) :: ext
character(len=140) :: filename

EXI = .True.
jext = 0

do while(EXI)
    jext = jext + 1
    write(ext, '(A, i0)') '.', jext
    filename = TRIM(FNAME) // TRIM(ext)
    inquire(FILE=TRIM(filename), EXIST=EXI)
enddo
FNAME = TRIM(filename)

return
end subroutine set_filename
