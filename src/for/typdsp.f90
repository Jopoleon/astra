subroutine ERASXY()

use outcmn_inc, only: astra_gui, White

implicit none

integer :: JN0, JLR
character(len=118) :: STRI

JLR = astra_gui%Height - 125
JN0 = 0
STRI = repeat(' ', 80)

call setColor(White)
call textvm(JN0, JLR, STRI(1: 15), 15)
call textvm(astra_gui%Width - 83*astra_gui%dxlet + 2, JLR, STRI(1: 80), 80)

end subroutine ERASXY

!---------------------------------------------------------------------
subroutine writeData(CHORDN)

use const_inc, only: XOUT, NAB, NA, NA1, exp_header, RTOR, ABC, BTOR, IPL, TIME, TAU, constValues
use status_inc, only: MU, AMETR, RHO, FP
use io_mod, only: AWD, equ_file, exp_file
use outcmn_inc, only: LTOUT, NTOUT, NROUT, MOD10, NAMER, NAMET, &
    WarningColor, ROUT, RUNID, NRW, NTIMES, TTOUT, TOUT
use dbl2char, only: fmt4, fmt_xf
use json_vars, only: n_const

implicit none

integer, parameter :: NLINSC=50
character(len=40), parameter :: STRMN=' R=     a=     B=     I=     q=     <n>='

double precision, intent(in) :: CHORDN

integer :: NCH=0, NP1, MODEX, JBE, JEND, J, JEN, JJ, J1, ios
double precision :: YQ
character(len=6) :: CH6
character(len=118) :: STRI
character(len=132) :: FNAME, dat_dir, file_out
character(len=:), external :: set_filename

! NLINSC - maximum line number
MODEX = XOUT + 0.49
! MODEX = 0 [0, AB]  against "a"
! MODEX = 1 [0, ABC] against "a"
! MODEX = 2 [0, ROC] against "rho"
! MODEX = 3 [FP(1), FP(NA1)] against "psi"
! otherwise Unknown option => MODEX=0
NP1 = NAB
if (MODEX >= 1 .and. MODEX <= 3 .or. MOD10 == 3) NP1 = NA1

dat_dir = TRIM(AWD) // '/dat/'
call system('mkdir -p ' // TRIM(dat_dir))
FNAME = TRIM(dat_dir) // TRIM(exp_file) // '.' // TRIM(equ_file)
file_out = set_filename(FNAME)

call setColor(WarningColor)
write(*, *) '>>>  Data are written into the file: ' // TRIM(file_out)

open(7, file=TRIM(file_out), iostat=ios)

if (ios /= 0) then
    write(*, *) '>>> TYPDSP: Output file error'
    stop
endif

! Creating UPSTRI
STRI = exp_header(1:16)
STRI(17:) = STRMN
STRI(20: 23) = fmt4(RTOR)
STRI(27: 30) = fmt4(ABC)
STRI(34: 37) = fmt4(BTOR)
STRI(41: 44) = fmt4(IPL)
! Triangularity corrected MHD q (accoding to ITER guidelines)
! YQ =ELON(NA)**2
! YD =TRIA(NA)
! YQ=(1.+YQ*(1.+YD**2*(2.-1.2*YD)))/(MU(NA)*(1.+YQ))
YQ = 1./MU(NA)
STRI(48: 51) = fmt4(YQ)
STRI(57: 60) = fmt4(CHORDN)
write(STRI(62: 76), '(A, 1F6.3, A)') 'Time=', TIME, ' dt='
STRI(77: 80) = fmt4(1000.*TAU)
write(7, 104) STRI

if (NCH == 0) then

    if (MOD10 <= 5) then   ! Writing radial data
        JBE = 1
        JEND = 16
        do
            JEN = MIN0(NTOUT, JEND)
            write(7, 102) (NAMET(J), J=JBE, JEN)
            STRI = ' '
            STRI(1:5) = fmt_xf(TIME, 4)
            do J=JBE, JEN
                JJ = 7*(J - JBE) + 8
                STRI(JJ: JJ+5) = fmt_xf(TOUT(LTOUT, J), 5)
            enddo
            write(7, 104) STRI
            if (JEN == NTOUT) EXIT
            JBE  = JEN + 1
            JEND = JEN + 16
        enddo

        JBE  = 1
        JEND = 16

        do
            JEN = MIN0(NROUT, JEND)
            if (MODEX == 0 .or. MODEX == 1) then
                write(7, '("     a  ", 16(3X, 1A4))') (NAMER(J), J=JBE, JEN)
            elseif (MODEX == 2) then
                write(7, '("     rho", 16(3X, 1A4))') (NAMER(J), J=JBE, JEN)
            elseif (MODEX == 3 .or. MOD10 == 3) then
                write(7, '("     psi", 16(3X, 1A4))') (NAMER(J), J=JBE, JEN)
            else
                write(7, '("     ???", 16(3X, 1A4))') (NAMER(J), J=JBE, JEN)
            endif
            do J=1, NP1
                STRI = ' '
                do JJ=JBE, JEN
                    J1 = 7*(JJ - JBE + 1) + 1
                    STRI(J1: J1+5) = fmt_xf(ROUT(J, JJ), 5)
                enddo
! Different options for a radial variable 
                if (MODEX == 0) then
                    STRI(1: 5) = fmt_xf(AMETR(j), 4)
                elseif (MODEX == 1) then
                    STRI(1: 5) = fmt_xf(AMETR(j), 4)
                elseif (MODEX == 2) then
                    STRI(1: 5) = fmt_xf(RHO(j), 4)
                elseif (MODEX == 3 .or. MOD10 == 3) then
                    STRI(1: 5) = fmt_xf(FP(j), 4)
                else
                    STRI(1: 5) = fmt_xf(AMETR(j), 4)
                endif
                write(7, 104)STRI
            enddo
            if (JEN == NROUT) EXIT
            JBE  = JEN + 1
            JEND = JEN + 16
        enddo

    else if (MOD10 == 6) then ! Writing time data

        JBE  = 1
        JEND = 16
        do
            JEN = MIN0(NTOUT, JEND)
            write(7, 102) (NAMET(J), J=JBE, JEN)
            do J1=1, LTOUT - 1
                STRI = ' '
                STRI(1: 5) = fmt_xf(TTOUT(J1), 4)
                do J=JBE, JEN
                    JJ = 7*(J - JBE) + 8
                    STRI(JJ: JJ+5) = fmt_xf(TOUT(J1, J), 5)
                enddo
                write(7, 104) STRI
            enddo
            if (JEN == NTOUT) EXIT
            JBE  = JEN + 1
            JEND = JEN + 16
        enddo
    endif

endif !NCH=0

! Writing constants

write(7, '(10X, 1A80)') RUNID
write(7, '(A)') 'Constants'
J1 = 0
do JEN=1, 11
    STRI = ' '
    do J=1, 16
        J1 = J1 + 1
        if (J1 > n_const) EXIT
        CH6 = fmt_xf(constValues(J1), 5)
        JJ = 7*(J - 1) + 1
        STRI(JJ: JJ+5) = CH6
    enddo
    write(7, '(A)') TRIM(STRI)
enddo

if (NCH == 1) then

    if (MOD10 <= 5) then   ! Writing radial data
        JBE  = 1
        JEND = 16
        do
            JEN = MIN0(NTOUT, JEND)
            write(7, '(3X, "Time", 16(3X, 1A4))') (NAMET(J), J=JBE, JEN)
            STRI = ' '
            STRI(1: 5) = fmt_xf(TIME, 4)
            do J=JBE, JEN
                JJ = 7*(J - JBE) + 8
                STRI(JJ: JJ+5) = fmt_xf(TOUT(LTOUT, J), 5)
            enddo
            write(7, 104) STRI
            if (JEN == NTOUT) EXIT
            JBE  = JEN + 1
            JEND = JEN + 16
        enddo

        JBE  = 1
        JEND = NRW
        JEN  = MIN0(NROUT, JEND)
! Different options for radial variable 
        if (MODEX == 2) then
            write(7, '(8X, "rho ", 64(8X, 1A4))') (NAMER(J), J=JBE, JEN)
            do j=1, NP1
                write(7, 408) RHO(j), (ROUT(J, JJ), JJ=JBE, JEN)
            enddo
        elseif (MODEX == 3 .or. MOD10 == 3) then
            write(7, '(8X, "psi ", 64(8X, 1A4))') (NAMER(J), J=JBE, JEN)
            do j=1, NP1
                write(7, 408) FP(j), (ROUT(J, JJ), JJ=JBE, JEN)
            enddo
        else
            write(7, '(8X, "a   ", 64(8X, 1A4))') (NAMER(J), J=JBE, JEN)
            do j=1, NP1
                write(7, 408) AMETR(j), (ROUT(J, JJ), JJ=JBE, JEN)
            enddo
        endif

    else if (MOD10 == 6) then  ! Writing time data

        STRI=' '
        write(7, 104) STRI
        write(7, 104) STRI
        write(7, 104) STRI
        write(7, 104) STRI
        JBE  = 1
        JEND = MIN(NTOUT, NRW)
        do while(JBE < JEND)
            JEN = JBE + 7
            write(7, '(8X, "Time", 64(8X, 1A4))') (NAMET(J), J=JBE, JEN)
            do J1=1, LTOUT-1
                STRI = ' '
                STRI(1: 5) = fmt_xf(TTOUT(J1), 4)
                write(7, 408) TTOUT(J1), (TOUT(J1, J), J=JBE, min(JEN, JEND))
            enddo
            JBE = JBE + 8
        enddo
    endif
endif

close(7)

101 format(1X, 1A6, 1A111)
102 format('   Time', 16(3X, 1A4))
104 format(1X, 1A120)
408 format(1PE12.3, 64(1PE12.3))

return
end subroutine writeData

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
subroutine PUTXY(IX, IY)

! Prints x, y coordinates on GUI in "Step" mode

use outcmn_inc, only: astra_gui, plot_area, MOD10, &
    scale_bnd, resizeGraph,  NTIMES, TOUT, TTOUT, &
    IDT, IDX, MODEY, LTOUT, NRW, NTOUT, active_tab, &
    NWIND3, NAMET, White, Red, Blue
use status_inc, only: AMETR, SHIF, ELON, TRIA, FP, RHO
use const_inc, only: TIME, TINIT, TSCALE, NA, NA1, NAB, XOUT, AB, ABC, ROC, HRO
use dbl2char, only: fmt5
use numerical_tools, only: QUADIN

implicit none

integer, parameter :: JN0=0
integer, intent(in) :: IX, IY

integer :: JX, JY, JLR, j, j1, JC, JL, MODEX, JW, JN2
double precision :: DX, DY, YX, YX1, YY, YY1, YA, YA1, YD, YE, YT, &
    YRHO, YFP, YFPC, RZ2A
character(len=80) :: STRI
integer, external :: GETIME

JLR = astra_gui%Height - int(125*resizeGraph)
if (MOD10 <= 0) return

STRI = repeat(' ', 80)
STRI(7:25) = '(x, y)=(     ,     )'
JX = IX - 10
JY = IY - 10
if (plot_area%xmin > JX .or. JX > plot_area%xmax .or. plot_area%ymax > JY .or. JY > plot_area%ymin) then
    STRI = repeat(' ', 80)
    call setColor(White)
    call textvm(JN0, JLR, "               ", 15)
    call textvm(astra_gui%Width - 83*astra_gui%dxlet + 2, JLR, STRI(1: 80), 80)
    return
endif

if (MOD10 == 7) plot_area%nx_canvas = 2
DX = 1./plot_area%nx_canvas
DY = 1./plot_area%ny_canvas

YX1 =      (JX - plot_area%xmin + 0.)/(plot_area%xmax - plot_area%xmin)
YY1 = 1. - (JY - plot_area%ymax + 0.)/(plot_area%ymin - plot_area%ymax)

if (MOD10 == 6) then
   ! (window_width)/(step=IDX=23)/(n_labels)=592/23/25=1.0295652
    YX = TINIT + 1.029565*YX1*abs(TSCALE)
    do j=1, plot_area%ny_canvas
        YY1 = YY1 - DY
        if (YY1 < 0) EXIT
    enddo
    YY = (YY1 + DY)/DY
    JN2 = plot_area%height + 15 + 3*astra_gui%dylet
    JC = 0
    JL = 0
    YY1 = max(TIME, TTOUT(LTOUT-1), TTOUT(LTOUT))
    YY1 = min(YX, YY1)
    YY1 = max(TTOUT(1), YY1)
    j = GETIME(YX, TTOUT, LTOUT)
    do J1=1, NTOUT
        JW = NWIND3(J1) - 8*active_tab(MOD10)
        if (NAMET(J1) == '    ') JW = 0
        if (JW > 0 .and. JW <= 8) call down_label(j, TOUT) ! for all plotting modes
    enddo
    call setColor(Red)
    STRI(1 :  5) = 'Time='
    STRI(6 : 10) = fmt5(YY1)
    STRI(11: 11) = 's'
    call textvm(astra_gui%dxlet, JN2 - int(2.5*astra_gui%dylet), STRI, 11)
    return
else if (MOD10 == 8) then
    YX = YX1*scale_bnd*plot_area%canvas_width/IDT/IDX
    YY = (YY1 - 0.5)*scale_bnd*plot_area%canvas_height/IDT/IDX
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
    do j=1, plot_area%nx_canvas
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

    do j=1, plot_area%ny_canvas
        YY1 = YY1 - DY
        if (YY1 < 0) EXIT
    enddo
    YY = (YY1 + DY)/DY
    if (MOD10 > 1 .and. MOD10 < 6 .and. MODEY == -1) then
        if (JY - plot_area%ymax > (plot_area%ymin - plot_area%ymax)/plot_area%ny_canvas) YY = YY - 1.
    endif

endif

STRI(14: 18) = fmt5(YX)
STRI(20: 24) = fmt5(YY)

call setColor(Blue)
call textvm(JN0, JLR, "               ", 15)
call textvm(astra_gui%Width - 83*astra_gui%dxlet, JLR, STRI(1: 80), 80)

return
end subroutine putxy

!---------------------------------------------------------------------
integer function plotMode(mod_10, mode_y)

implicit none
  
integer, intent(in) :: mod_10, mode_y

if (mod_10 <= 1 .or. mod_10 >= 7) then
    plotMode = mod_10
elseif (mod_10 >= 2 .and. mod_10 <= 5) then
    if (mode_y == 1) then
        plotMode = 2
    else
        plotMode = 3
    endif
else ! mod_10 = 6
    if (mode_y == 1) then
        plotMode = 5
    else if (mode_y == 0) then
        plotMode = 9
    else
        plotMode = 6
    endif
endif

return
end function plotMode

!---------------------------------------------------------------------
subroutine set_plot_area(plot_mode)
!----------------------------------------------------------------------|
! Input: MODEY
! Output: plot_area%xmin - x_left  of the graphic area
!  plot_area%xmax - x_right of the graphic area
!  plot_mode - for using in set_plot
!  NST - for using in set_plot
!----------------------------------------------------------------------|

use outcmn_inc, only: astra_gui, plot_area, MOD10, MODEY, NST

implicit none

integer, intent(in) :: plot_mode

integer :: n_str_up, n_str_down

! n_str_up, n_str_down - number of text strings up & down
SELECT CASE(plot_mode)
CASE(1)
    n_str_up   = 2
    n_str_down = 5
    plot_area%xmin = 0
    plot_area%xmax = plot_area%width
    plot_area%nx_canvas = 4
    plot_area%ny_canvas = 2
! modes 2, 3, 4, 5 at y-mode = +1, dummy mode (plot_mode=4 - not used)
CASE(2, 4, 7)
    n_str_up   = 2
    n_str_down = 5
    plot_area%xmin = 0
    plot_area%xmax = plot_area%width
    plot_area%nx_canvas = 2
    plot_area%ny_canvas = 1
! modes 2, 3, 4, 5 at y-mode = -1
CASE(3)
    n_str_up   = 2
    n_str_down = 5
    plot_area%xmin = 0
    plot_area%xmax = plot_area%width
    plot_area%nx_canvas = 2
    plot_area%ny_canvas = 2
! mode # 6 (time) at y-mode=1 (2 windows)
CASE(5)
    n_str_up   = 1
    n_str_down = 1
    plot_area%xmin = 6*astra_gui%dxlet
    plot_area%xmax = plot_area%width
    plot_area%nx_canvas = 1
    plot_area%ny_canvas = 2
! mode # 6 (time) at y-mode=-1 (4 windows)
CASE(6)
    n_str_up   = 1
    n_str_down = 1
    plot_area%xmin = 6*astra_gui%dxlet
    plot_area%xmax = plot_area%width
    plot_area%nx_canvas = 1
    plot_area%ny_canvas = 4
! mode 8 (equilibrium)
CASE(8)
    n_str_up   = 1
    n_str_down = -2
    plot_area%xmin = 0
    plot_area%xmax = 0.7*plot_area%width
    plot_area%nx_canvas = 1
    plot_area%ny_canvas = 1
! mode # 6 (time) at y-mode=0 (1 window), mode 9 (user's plot)
CASE(9)
    n_str_up   = 1
    n_str_down = 1
    plot_area%xmin = 6*astra_gui%dxlet
    plot_area%xmax = plot_area%width
    plot_area%nx_canvas = 1
    plot_area%ny_canvas = 1
END SELECT

plot_area%ymax = n_str_up*astra_gui%dylet + 1
plot_area%ymin = plot_area%height - n_str_down*astra_gui%dylet - 1
plot_area%canvas_width  = (plot_area%xmax - plot_area%xmin)/plot_area%nx_canvas
plot_area%canvas_height = (plot_area%ymin - plot_area%ymax)/plot_area%ny_canvas

return
end subroutine set_plot_area

!---------------------------------------------------------------------
subroutine set_plot(plot_mode)
! Subroutine draw frame for different modes

use const_inc, only: TSCALE, TINIT, AWALL
use outcmn_inc, only: astra_gui, plot_area, Black, MOD10, KPRI, &
    IDX, IDT, scale_bnd, pixel_ymid, meter2pixel
use dbl2char, only: fmt_xf
use char_manip, only: len_trim_tab

implicit none

integer, parameter :: LENG=3, JN0=0
integer, intent(in) :: plot_mode

integer :: JJ, J, TIMWIN, JX, JY, &
    XP, XM, YP, YM, JXSCM, LYM
double precision :: DY, YY, TIND, scale_fac
character(len=5) :: XF4
character(len=6) :: CH6
character(len=80) :: COMMENT

! IDX, DY  - X & Y distance (in points) between X & Y axis labels
! IDT      - distance (in labels) between longer labels in modes 6&8
! LENG     - label length

if (plot_mode <= 0) return

TIMWIN = 0
if (plot_mode == 5 .or. plot_mode == 6 .or. plot_mode == 8 .or. plot_mode == 9) TIMWIN = 1

call setColor(Black)
call rectvm(0, JN0, JN0, astra_gui%Width-1, astra_gui%Height-1)

if (MOD10 == 6)  then
    call setColor(Black)
    j = astra_gui%Width - 20*astra_gui%dxlet + 1
    jj = plot_area%height + astra_gui%dylet
    call textvm(j, jj, 'time, s', 7)
endif

! Vertical lines & Y-labels
if (KPRI >= 1 .and. KPRI <= 2) then
    write(COMMENT, '(A)') "Vertical lines"
    j = len_trim_tab(COMMENT)
    call pscom(COMMENT, j)
endif
call setColor(Black)

! Skipping from a subplot to the next along x-axis
do JJ=plot_area%xmin, plot_area%xmax, plot_area%canvas_width
    JX = MIN(plot_area%xmax, JJ)
    XP = MIN(plot_area%xmax, JX + LENG)
    XM = MAX(plot_area%xmin, JX - LENG)
    call drawvm(0, JX, plot_area%ymax, JX, plot_area%ymin)
! Y-line labels
    if (KPRI >= 1 .and. KPRI <= 2) then
        write(COMMENT, '(A)') "Y-line labels"
        j = len_trim_tab(COMMENT)
        call pscom(COMMENT, j)
    endif
    DY = (plot_area%ymin - plot_area%ymax)/20.
    YY = dble(plot_area%ymin)
    do
        JY = YY
        if (plot_mode /= 8) call drawvm(0, XM, JY, XP, JY)
        YY = YY - DY
        if (sign(1.d0, DY)*(YY - dble(plot_area%ymax)) < 0.d0) EXIT
    enddo
enddo

! Horizontal lines
if (KPRI >= 1 .and. KPRI <= 2) then
    write(COMMENT, '(A)') "Horizontal lines"
    j = len_trim_tab(COMMENT)
    call pscom(COMMENT, j)
endif

do JY=plot_area%ymin, plot_area%ymax, -plot_area%canvas_height
    YM = MAX(plot_area%ymax, JY - LENG)
    if (plot_mode == 4 .or. plot_mode == 5 .or. plot_mode == 6) then
        YP = JY
    else
        YP = MIN(plot_area%ymin, JY + LENG)
    endif
    call drawvm(0, plot_area%xmin, JY, plot_area%xmax, JY)
! X-line labels
    if (KPRI >= 1 .and. KPRI <= 2) then
        write(COMMENT, '(A)') "X-line labels"
        j = len_trim_tab(COMMENT)
        call pscom(COMMENT, j)
    endif
    JX = plot_area%xmin
    IDX = 16
    JXSCM = plot_area%xmax - IDX
    if (TIMWIN == 1) then
        IDX = 23
        JXSCM = plot_area%xmax
    endif
    do J=1, 100
        JX = JX + IDX
        if (JX > JXSCM) EXIT
        if (TIMWIN == 1 .and. J/IDT*IDT == J) then
            LYM = YM - 2
        else
            LYM = YM
        endif
        if (LYM > plot_area%ymax + LENG)  call drawvm(0, JX, LYM, JX, YP)
    enddo
enddo

if (plot_mode == 8) then
    JY = (plot_area%ymax + plot_area%ymin)/2
    do j=0, 10
        jj = JY + IDX*j
        if (jj < plot_area%ymin) call drawvm(0, XM, JJ, XP, JJ)
        jj = JY - IDX*j
        if (jj > plot_area%ymax) call drawvm(0, XM, JJ, XP, JJ)
    enddo
endif

if (MOD10 == 6) then
! time-axis legend:
    call setColor(Black)
    JJ = plot_area%height - astra_gui%dylet + 12
    do J=0, plot_area%width, IDX
        JX = (J - IDT)*IDT + plot_area%xmin
        if (JX > plot_area%width) CYCLE
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
    scale_bnd = 1.3*AWALL
    scale_fac = dble(plot_area%height)/350.
    IDX = IDX*scale_fac
    JJ = plot_area%ymin + astra_gui%dylet + 2
    call setColor(Black)
    do J=1, 5
        JX = IDX*IDT*J - 24
        if (JX > JXSCM) CYCLE
        YY = J*scale_bnd
        XF4 = fmt_xf(YY, 4)
        call textvm(JX, JJ, XF4, 5)
    enddo
! vertical axis labels
    do J=-1, 1
        JX = (plot_area%ymin + plot_area%ymax + astra_gui%dylet)/2 + IDT*IDX*J - 0.5*astra_gui%dylet
        YY = -J*scale_bnd
        XF4 = fmt_xf(YY, 4)
        call textvm(plot_area%xmax + 2, JX, XF4, 5)
    enddo
endif

if (KPRI >= 1 .and. KPRI <= 2) then
     write(COMMENT, '(A)') "Frame done"
     j = len_trim_tab(COMMENT)
     call pscom(COMMENT, j)
endif

pixel_ymid  = 0.5*(plot_area%ymax + plot_area%ymin)
meter2pixel = dble(IDX*IDT)/scale_bnd

return
end subroutine set_plot

!---------------------------------------------------------------------
subroutine down_label(jt_in, TOUT)
!---------------------------------------------------------------------
! Time dependences for radial output
! Curve to digit conversion
!  Input:  jt defines the current time
!
! if mod10 != 6 or call from run then jt = LTOUT
!---------------------------------------------------------------------

use outcmn_inc, only: astra_gui, plot_area, LTOUT, MOD10, NTOUT, NRW, NTIMES, &
    NWIND3, active_tab, NAMET, Black, Blue, curves_per_frame, resizeGraph
use dbl2char, only: fmt_xf

implicit none

integer, intent(in) :: jt_in
double precision, intent(in) :: TOUT(NTIMES, *)

integer :: jt, JN2, JN0, JEND, JB, JL, JC, JW, JJ, J, fshift
character(len=5) :: XF4
character(len=7) :: XF7
character(len=80) :: STRI, STRIN

if (jt_in == 0) then
    jt = LTOUT
    call setColor(Black)
else
    jt = jt_in
    call setColor(Blue)
endif

fshift = int(10*resizeGraph)

if (MOD10 == 6) then

    JN0 = 6*astra_gui%dxlet
    JN2 = plot_area%height + FSHIFT + 3*astra_gui%dylet + 5
    JC = 0
    JL = 0
    write(STRIN, '(79X, 1A1)') ' '
    write(STRI , '(79X, 1A1)') ' '
    do j=1, NTOUT
        JW = NWIND3(j) - curves_per_frame(MOD10)*active_tab(MOD10)
        if (NAMET(j) == '    ') JW = 0
        if (JW <= 0 .or. JW > curves_per_frame(MOD10)) CYCLE
        JC = JC + 1   ! Actual curve number in the window
        jj = 8*JC - 6
        if (jj > 74) CYCLE
        if (jj >= 66) JN0 = 5*astra_gui%dxlet
        if (jj == 74) JN0 = -astra_gui%dxlet
! curve #, win #, chan #, mode 6, screen #
        XF7 = fmt_xf(TOUT(jt, j), 6)
        STRI (jj: jj+6) = XF7
        STRIN(jj: jj+6) = '  ' // NAMET(J) // ' '
        JL = max(JL, jj + 6)
    enddo
    call textvm(JN0, JN2, STRI, JL)
    JN2 = JN2 - astra_gui%dylet + 1
    call textvm(JN0, JN2, STRIN, JL)
else
    JB = 1
    JN0 = 0
    JN2 = plot_area%height - 5*astra_gui%dylet + FSHIFT + 2
    do
        JEND = MIN(JB + 15, NTOUT)
        do J=JB, JEND
            XF4 = fmt_xf(TOUT(jt, J), 4)
            if (NAMET(J) == ' ') XF4 = '    '
            JJ = 5*(J - JB + 1) - 4
            STRI(JJ: JJ+4) = XF4
        enddo
        JN2 = JN2 + 2*astra_gui%dylet + 2
        JL = 5*(JEND - JB + 1)
        call textvm(JN0, JN2, STRI, JL)
        write(STRI, '(16(1X, 1A4))') (NAMET(J), J=JB, JEND)
        JN2 = JN2 - astra_gui%dylet + 1
        call textvm(JN0, JN2, STRI, JL)
        JN2 = JN2 + astra_gui%dylet - 1
        if (JEND == NTOUT .or. JEND == NRW) EXIT
        JB = JB + 16
    enddo
endif

return
end subroutine down_label

!---------------------------------------------------------------------
! Upper string of the Astra graphic window
subroutine up_label(YN, YQ)

use outcmn_inc, only: astra_gui, active_tab, MOD10, Black, Blue
use const_inc, only: RTOR, BTOR, IPL, ABC, exp_header
use dbl2char, only: fmt40

implicit none

double precision, intent(in) :: YN, YQ

character(len=2) :: CHR
character(len=42) :: STRMN

STRMN(1: 42) = ' R=     a=     B=     I=     q=     n=    '

STRMN( 4:  7) = fmt40(RTOR)
STRMN(11: 14) = fmt40(ABC)
STRMN(18: 21) = fmt40(BTOR)
STRMN(25: 28) = fmt40(IPL)
STRMN(32: 35) = fmt40(YQ)
STRMN(39: 42) = fmt40(YN)

call setColor(Black)
call textvm(0, 2, exp_header(1: 15) // STRMN(1: 42), 56)
call setColor(Blue)
write(CHR, '(1I2)') active_tab(MOD10) + 1
call textvm(astra_gui%width - 2*astra_gui%dxlet, astra_gui%dylet + 1, CHR, 2) ! Screen No.
call setColor(Black)
call rectvm(0, 0, 0, astra_gui%Width - 1, astra_gui%Height - 1) ! Outer frame

return
end subroutine up_label

!---------------------------------------------------------------------
subroutine TIMEDT(TIME, DT)

use outcmn_inc, only: astra_gui, astra_gui_ref, Black
use dbl2char, only: fmt50

implicit none

integer, parameter :: FSHIFT=2, str_len=19
double precision, intent(in) :: TIME, DT

character(len=str_len) :: STRI

STRI(1 : 5 ) = 'Time='
STRI(11: 16) = ' dt='
STRI( 6: 10) = fmt50(TIME)
STRI(15: 19) = fmt50(DT)
call setColor(Black)
call textvm(astra_gui%width - (str_len+3)*astra_gui_ref%dxlet, FSHIFT, STRI, str_len)

return
end subroutine TIMEDT

!---------------------------------------------------------------------
subroutine const2ps
! Appending the list of constants to a PS file

use const_inc, only: constValues, varValues
use char_manip, only: null_ch
use outcmn_inc, only: resizeGraph
use dbl2char, only: fmt_xf
use json_vars, only: n_const, n_var, varNames

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

JNY = 540*resizeGraph
JNX = 10*resizeGraph

ps_loop: do J2=1, 11

    JDUM = 2*(J2 - 1) + 1
    STRI = CONN(JDUM)

    do J=1, 4
        J1 = J1 + 1
        if (J1 > n_const) EXIT ps_loop
        CH6 = fmt_xf(constValues(J1), 5)
        JJ = 7*(J - 1) + 8
    enddo

    STRI(JJ: JJ+5) = CH6
    STRI(JJ+12: JJ+18) = CONN(JDUM+1)

    do J=5, 8
        J1 = J1 + 1
        if (J1 > n_const) EXIT ps_loop
        CH6 = fmt_xf(constValues(J1), 5)
        JJ = 7*(J - 1) + 20
    enddo

    STRI(JJ: JJ+5) = CH6
    JNY = JNY + 17

    call textvm(JNX, JNY, STRI, 75)

enddo ps_loop

! Writing variables
j1 = 1
JNY = JNY + 20
JDUM = n_var - 48

do j=1, JDUM
    CH6 = fmt_xf(varValues(j), 5)
    STRI(j1: j1+19) = varNames(j) // '=' // CH6 // '     '
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

!---------------------------------------------------------------------
function set_filename(fname_in) result(fname_out)
    implicit none

    character(len=*), intent(in) :: fname_in
    character(len=:), allocatable :: fname_out

    integer :: jext
    logical :: fileExists
    character(len=:), allocatable :: filename
    character(len=16) :: ext   ! long enough for big integers

    fileExists = .true.
    jext = 0

    do while (fileExists)
        jext = jext + 1
        write(ext, '(".", i0)') jext
        filename = trim(fname_in) // trim(ext)
        inquire(file=filename, exist=fileExists)
    enddo

    fname_out = filename
end function set_filename
