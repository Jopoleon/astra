module graph_utils

use char_manip, only: null_ch
use status, only: NRD

implicit none

integer, parameter :: White=0, Black=1, Red=2, Blue=3, Green=5, &
     WarningColor=30, EraseColor=31, Magenta=14, Pink=13, &
     nplots_max=32, NRW=128, NTIMES=1024, plot_modes=9

integer, dimension(NRW), parameter :: IP1 = (/ &
      1,   3,   5,   7,   9,  11,  13,  15,   2,   4,   6,  8,   10,  12,  14,  16, &
     17,  19,  21,  23,  25,  27,  29,  31,  18,  20,  22, 24,   26,  28,  30,  32, &
     33,  35,  37,  39,  41,  43,  45,  47,  34,  36,  38, 40,   42,  44,  46,  48, &
     49,  51,  53,  55,  57,  59,  61,  63,  50,  52,  54, 56,   58,  60,  62,  64, &
     65,  67,  69,  71,  73,  75,  77,  79,  66,  68,  70, 72,   74,  76,  78,  80, &
     81,  83,  85,  87,  89,  91,  93,  95,  82,  84,  86, 88,   90,  92,  94,  96, &
     97,  99, 101, 103, 105, 107, 109, 111,  98, 100, 102, 104, 106, 108, 110, 112, &
    113, 115, 117, 119, 121, 123, 125, 127, 114, 116, 118, 120, 122, 124, 126, 128 /)
integer, dimension(NRW), parameter :: IP2 = (/ &
      1,   2,   9,  10,   3,   4,  11,  12,   5,   6,  13,  14,   7,   8,  15,  16, &
     17,  18,  25,  26,  19,  20,  27,  28,  21,  22,  29,  30,  23,  24,  31,  32, &
     33,  34,  41,  42,  35,  36,  43,  44,  37,  38,  45,  46,  39,  40,  47,  48, &
     49,  50,  57,  58,  51,  52,  59,  60,  53,  54,  61,  62,  55,  56,  63,  64, &
     65,  66,  73,  74,  67,  68,  75,  76,  69,  70,  77,  78,  71,  72,  79,  80, &
     81,  82,  89,  90,  83,  84,  91,  92,  85,  86,  93,  94,  87,  88,  95,  96, &
     97,  98, 105, 106,  99, 100, 107, 108, 101, 102, 109, 110, 103, 104, 111, 112, &
    113, 114, 121, 122, 115, 116, 123, 124, 117, 118, 125, 126, 119, 120, 127, 128/)
integer, dimension(NRW), parameter :: IP30 = (/ &
      1,   3,   5,   7,   2,   4,   6,   8,   9,  11,  13,  15,  10,  12,  14,  16, &
     17,  18,  21,  22,  25,  26,  29,  30,  19,  20,  23,  24,  27,  28,  31,  32, &
     33,  34,  37,  38,  41,  42,  45,  46,  35,  36,  39,  40,  43,  44,  47,  48, &
     49,  50,  53,  54,  57,  58,  61,  62,  51,  52,  55,  56,  59,  60,  63,  64, &
     65,  66,  69,  70,  73,  74,  77,  78,  67,  68,  71,  72,  75,  76,  79,  80, &
     81,  82,  85,  86,  89,  90,  93,  94,  83,  84,  87,  88,  91,  92,  95,  96, &
     97,  98, 101, 102, 105, 106, 109, 110,  99, 100, 103, 104, 107, 108, 111, 112, &
    113, 114, 117, 118, 115, 116, 119, 120, 121, 122, 125, 126, 123, 124, 127, 128 /)
integer, dimension(NRW), parameter :: IP31 = (/ &
      1,   2,   5,   6,   3,   4,   7,   8,   9,  10,  13,  14,  11,  12,  15,  16, &
     17,  18,  21,  22,  19,  20,  23,  24,  25,  26,  29,  30,  27,  28,  31,  32, &
     33,  34,  37,  38,  35,  36,  39,  40,  41,  42,  45,  46,  43,  44,  47,  48, &
     49,  50,  53,  54,  51,  52,  55,  56,  57,  58,  61,  62,  59,  60,  63,  64, &
     65,  66,  69,  70,  67,  68,  71,  72,  73,  74,  77,  78,  75,  76,  79,  80, &
     81,  82,  85,  86,  83,  84,  87,  88,  89,  90,  93,  94,  91,  92,  95,  96, &
     97,  98, 101, 102, 105, 106, 109, 110,  99, 100, 103, 104, 107, 108, 111, 112, &
    113, 114, 117, 118, 115, 116, 119, 120, 121, 122, 125, 126, 123, 124, 127, 128 /)

type astra_X11_window
    integer :: Width, Height, Xpos, Ypos, dxlet, dylet, LineWidth, yMessage
    double precision :: resizeGraph
    character(len=128) :: title='Per aspera ad ASTRA' // null_ch
endtype astra_X11_window

type plot_frame
    integer :: width, height, xmin, xmax, ymin, ymax, nx_canvas, ny_canvas, canvas_height, canvas_width
endtype plot_frame

! Colors, array AstraColorNum in Astra2XW.c
integer, dimension(NRW) :: NWIND1, NWIND3, NWIND4, NWIND7, NWINDX
integer :: NTOUT, NROUT, LTOUT, IPOUT, MOD10, NXOUT
integer :: MODEY, IDX, IDT, KPRI, NST, AVERS, ARLEAS, AEDIT
integer, dimension(plot_modes) :: active_tab, curves_per_frame
double precision, dimension(NRW)   :: GRAL, GRAP, OSHIFT, OSHIFR, SCALET, SCALER
double precision, dimension(NRD, NRW) :: ROUT
double precision :: TIM7(4), scale_bnd, pixel_ymid, meter2pixel
double precision :: TTOUT(NTIMES), TOUT(NTIMES, NRW)

character(len=4), dimension(NRW) :: NAMET, NAMER
character(len=6), dimension(NRW) :: NAMEX
character(len=6), allocatable :: DTNAME(:)
character(len=6), dimension(4) :: NAM7
character(132) :: VERSION, RUNID
type(astra_X11_window) :: astra_gui_ref, astra_gui
type(plot_frame) :: plot_area_ref, plot_area

contains

!---------------------------------------------------------------------
    subroutine gui_init()

    use read_input, only: resize, n_sbr
    use scalars, only: AB, TINIT, TSCALE, XOUT
    use json_vars, only: n_control, controlNames

    integer :: i, j, ios, j0, j1, jgrid, jj, plot_mode
    character(len=132) :: STRI

    TTOUT(1) = -1.d10
    TIM7(1) = TINIT
    TIM7(3) = abs(TSCALE)/8.

! Constants

    pixel_ymid  = 0.
    meter2pixel = 0.

    VERSION = repeat(' ', 32)

    IDT = 5

! Screen parameters: default (-resize 1)

    astra_gui_ref%Width  = 660
    astra_gui_ref%Height = 550
    astra_gui_ref%Xpos   = 470
    astra_gui_ref%Ypos   = 10
    astra_gui_ref%dxlet  = 8
    astra_gui_ref%dylet  = 13
    astra_gui_ref%LineWidth = 1
    astra_gui_ref%yMessage  = 426
    astra_gui_ref%resizeGraph = 1.d0

    plot_area_ref%width  = 640
    plot_area_ref%height = 350
    plot_area_ref%nx_canvas = 0
    plot_area_ref%ny_canvas = 0
    plot_area_ref%xmin = 0
    plot_area_ref%xmax = 0
    plot_area_ref%ymin = 0
    plot_area_ref%ymax = 0

! Astra colors: 
!   #0 - background, ##1-7 - plots 1-7
!   #7 - also color for REPORT (user's output)
!   #8 - not used
!   #9 - delete curve, #10 - standard text color (black)
!   #11 - Black
!   #12 - warning messages, axis upper marks in View
!   #13 -> #15 -reserved (black)

    MODEY  = 1
    LTOUT  = 1
    IPOUT  = 1
    active_tab = 0
    OSHIFT = 0.
    OSHIFR = 0.
    GRAL   = 0.
    GRAP   = AB

    NAM7 = (/ 'Tmin  ', 'Tmax  ', 'Tmark ', 'Style ' /)
    TIM7 = (/ 0, 9999, 9999, 1 /)

! Output windows

    NWIND1 = (/ (j, j=1, NRW) /)
    NWIND3 = (/ (j, j=1, NRW) /)
    NWIND7 = (/ (j, j=1, NRW) /)
    do j=1, 4
        NWIND4(j) = 4*j-3
        NWIND4(j+4)  = NWIND4(j) + 2
        NWIND4(j+8)  = NWIND4(j) + 1
        NWIND4(j+12) = NWIND4(j) + 3
    enddo
    do j=17, NRW
        NWIND4(j) = NWIND4(j-16) + 16
    enddo
    curves_per_frame = (/ 16, 8, 8, 2, 2, 8, 4, 0, 0 /)

    DTNAME(1: n_control) = controlNames
    do j=1, n_sbr
        i = (j-1)*4 + n_control
        write(DTNAME(i+1), '(A, i0)') 'DTeq', j
        write(DTNAME(i+2), '(A, i0)') 'BEeq', j
        write(DTNAME(i+3), '(A, i0)') 'ENeq', j
        write(DTNAME(i+4), '(A, i0)') ' Keq', j
    enddo

!-------------------------
! Parse file "exe/version"

    open(131, FILE='exe/version', iostat=ios)
    if (ios /= 0) then
        write(*, *) '>>> Warning: Unknown version'
    else
        do j=1,5
            read(131, '(A)') STRI
        enddo
        j = index(STRI, 'Version')
        VERSION = STRI(j: j+30) // null_ch
        close(131)
        j0 = index(VERSION, '.')
        if (j0 == 0) then
            write(*, *) '>>> Warning: Unknown version'
        else
            read(VERSION(j0-1: j0-1), *) AVERS 
            read(VERSION(j0+1: j0+1), *) ARLEAS 
            j1 = INDEX(VERSION(j0+1:), '.')
            if (j1 == 0) then
                AEDIT = 0
            else
                read(VERSION(j0+j1+1: j0+j1+1), *) AEDIT
            endif
        endif
    endif

    RUNID = runidLabel()
    jj = max(0, (15 + NTOUT - 64)/16)

!-----------------
! Main GUI window

    astra_gui%LineWidth = int(0.85*resize) + astra_gui_ref%LineWidth
    astra_gui%dxlet    = resize*astra_gui_ref%dxlet
    astra_gui%dylet    = resize*astra_gui_ref%dylet
    astra_gui%yMessage = resize*astra_gui_ref%yMessage + 135
    astra_gui%Width    = resize*astra_gui_ref%width
    astra_gui%Height   = resize*(astra_gui_ref%Height + 2*jj*resize*(astra_gui_ref%dylet + 2))
    astra_gui%Xpos  = astra_gui_ref%Xpos
    astra_gui%Ypos  = astra_gui_ref%Ypos
    astra_gui%title = astra_gui_ref%title
    astra_gui%resizeGraph = resize

    plot_area%width  = resize*plot_area_ref%width
    plot_area%height = resize*plot_area_ref%height

    call initvm(astra_gui%xpos, astra_gui%ypos, astra_gui%Width, astra_gui%Height, &
        astra_gui%LineWidth, astra_gui%title, LEN(astra_gui%title)) ! Initialise graphic window
    plot_mode = 1
    NST = 0
    MOD10 = 1
    plot_mode = plotMode(MOD10, MODEY)
    call set_plot_area(plot_mode)
    call set_plot(plot_mode)

    jgrid = XOUT + 0.49

    call taskmenu(jgrid) ! Task menu
    call textbf(0, astra_gui%Height - int(104*astra_gui%resizeGraph), RUNID, 80) ! Task ID

    end subroutine gui_init

!---------------------------------------------------------------------
    function runidLabel() result(runid_label)

    use read_input, only: equ_file, exp_file

    character(len=:), allocatable :: runid_label

    integer :: time_arr(8)
    integer :: year, month, day, hour, minute, j
    character(len=3)  :: vers
    character(len=32) :: datetime
    character(len=256) :: tmp

    call date_and_time(values=time_arr)

    year   = time_arr(1)
    month  = time_arr(2)
    day    = time_arr(3)
    hour   = time_arr(5)
    minute = time_arr(6)

    write(datetime, "(I2.2,'-',I2.2,'-',I2.2,' ',I2.2,':',I2.2)") &
        day, month, mod(year,100), hour, minute

    j = index(VERSION, 'Version')
    vers = VERSION(j+8:j+10)

    tmp = "ASTRA " // vers // " -- " // trim(datetime) // &
          " -- Model: " // trim(equ_file) // &
          " -- Data: "  // trim(exp_file)

    runid_label = trim(tmp)

    end function runidLabel

!---------------------------------------------------------------------
    subroutine SCAL(NOUT, SN, SO, OUT, NP, NDIM)
!---------------------------------------------------------------------
!  Input:
!     NOUT Total number of channels (NROUT or NTOUT)
!     SO(NOUT) Current scales
!     NP Upper index boundary for scale definition (= NDIM-3)
!     NDIM Upper index boundary according to the array description
!     OUT(NDIM, NOUT) Array of curves for scale definition (ROUT or TOUT)
!  Output:
!     SN(NOUT) Returned scales

    integer, intent(in) :: NOUT, NP, NDIM
    double precision, intent(in)  :: SO(*), OUT(NDIM, *)
    double precision, intent(out) :: SN(*)

    integer :: J, JJ
    double precision :: SC

! Exclude 3 edge/central points
    var_loop: do J=1, NOUT
        if (SO(J) > 0) then
            SN(J) = SO(J)
        else if (SO(J) .EQ. 0) then
            SN(J) = SCALA(OUT(1, J), NP)
        else
            do JJ = 1, J-1
                if (SO(J) .EQ. SO(JJ)) then
                    SN(J) = SN(JJ)
                    CYCLE var_loop
                endif
            enddo
            SC = 0.
            do JJ=J, NOUT
                if (SO(J) .EQ. SO(JJ)) then
                    SC = MAX(SC, SCALA(OUT(1, JJ), NP))
                    if (ABS(SO(JJ) + JJ) < .01) then
                        SN(J) = SCALA(OUT(1, JJ), NP)
                        CYCLE var_loop
                    endif
                endif
            enddo
            SN(J) = SC
        endif
    enddo var_loop

    end subroutine SCAL

!---------------------------------------------------------------------
    double precision function SCALA(Y, NJ)

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

!---------------------------------------------------------------------
    subroutine CMARK(xpos_in, ypos_in, prof_yscale, yshift, prof_name, STYL, jplot_parity)
! Mark variable/scale in 1 & 2 modes

    use dbl2char, only: fmt_smart

    integer, intent(in) :: STYL, xpos_in, ypos_in, jplot_parity
    double precision, intent(in) :: prof_yscale, yshift
    character(len=4), intent(in) :: prof_name

    integer :: POINT(2), xpos, ypos, str_len, str_pixels
    character(len=15) :: name_scale_label

    ypos = ypos_in

    name_scale_label(:)  = ' '
    name_scale_label(1: 4) = fmt_smart(prof_yscale, 4)
    name_scale_label(6: 10) = prof_name
    if (yshift /= 0) then
        if (yshift < 0) then
            name_scale_label(11: 11) = '-'
        else
            name_scale_label(11: 11) = '+'
        endif
        name_scale_label(12: 15) = fmt_smart(abs(yshift), 4)
    endif
    str_len = LEN_TRIM(name_scale_label)
    str_pixels = 8*str_len
    xpos = xpos_in + 2 + jplot_parity*(plot_area%canvas_width - str_pixels - 4)
    call textvm(xpos, ypos, TRIM(name_scale_label), str_len)
    if (STYL > 0) then ! If clicking 'Style' in ASTRA graphic window
        POINT(1) = xpos_in + str_pixels + 8 + jplot_parity*(plot_area%canvas_width - 2*str_pixels - 16)
        POINT(2) = ypos - 5
        call NMARK(POINT, STYL)
    endif

    end subroutine CMARK

!---------------------------------------------------------------------
    subroutine CMARKT(xpos_in, ypos_in, sig_yscale, yshift, sig_name, STYL)
! Mark variable/scale in 6th (time) mode

    use dbl2char, only: fmt_smart

    integer, intent(in) :: xpos_in, ypos_in, STYL
    double precision, intent(in) :: sig_yscale, yshift
    character(len=4), intent(in) :: sig_name

    integer :: str_len, xpos, ypos, plot_arr(2)
    character(len=10) :: name_shift_label
    character(len=4 ) :: F4

    xpos = xpos_in
    ypos = ypos_in
    name_shift_label(1: 10) = '          '
    name_shift_label(2: 5)  = sig_name
    if (yshift /= 0) then
        F4 = fmt_smart(abs(yshift), 4)
        if (yshift < 0) name_shift_label(6: 10) = '-' // F4
        if (yshift > 0) name_shift_label(6: 10) = '+' // F4
        xpos = xpos_in - astra_gui%dxlet
        if (name_shift_label(1: 1) == ' ') xpos = xpos - astra_gui%dxlet
        str_len = LEN_TRIM(name_shift_label(2: 10))
        if (str_len <= 6) xpos = xpos + astra_gui%dxlet
    endif

    call textvm(xpos, ypos, TRIM(name_shift_label), LEN_TRIM(name_shift_label))   ! type name+yshift
    plot_arr(1) = xpos + 3
    plot_arr(2) = ypos - 5
    F4 = fmt_smart(sig_yscale, 4)
    ypos = ypos + 15
    xpos = xpos + astra_gui%dxlet
    call textvm(xpos, ypos, F4, 4)   ! type scale
    if (STYL > 0) call NMARK(plot_arr, STYL)

    end subroutine CMARKT

!---------------------------------------------------------------------
    subroutine update_curve(NP, np_old, ICOLOR, STYL, xold, yold, xnew, ynew)

! The subroutine displays NP points of the float array YNEW
! NP  is a number of points to plot
! NPO is a number of points to erase
! The points of the array YOLD are used for erasing

    integer, intent(in) :: STYL, NP, np_old, ICOLOR
    double precision, intent(in), dimension(np) :: xold, yold, xnew, ynew

    if (np_old > 0) then
! erase the old curve
        call setColor(EraseColor)
        call plot_curve(np_old, STYL, xold, yold)
    endif

! draw a new curve

    call setColor(ICOLOR)
    call plot_curve(NP, STYL, xnew, ynew)

    end subroutine update_curve

!---------------------------------------------------------------------
    subroutine plot_curve(np, STYL, xplot, yplot)
      
    integer, intent(in) :: STYL, np
    double precision, intent(in), dimension(*) :: xplot, yplot

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
                call drawcurve(0, j1, xplot(j: j+j1-1), yplot(j: j+j1-1))
            enddo
            return
        endif
    else if (STYL > 0) then
        NM = max(10, NP/5)
        LE = NM/5*STYL    ! 1st marker position
        if (LE >= NM+2) LE = LE - NM
        LE = max(1, LE)
        do jj=LE, NP, NM
            PT1(1) = xplot(jj)
            PT1(2) = yplot(jj)
            call NMARK(PT1, STYL)
        enddo
    endif

    call drawcurve(0, NP, xplot(1:NP), yplot(1:np))

    end subroutine plot_curve

!---------------------------------------------------------------------
    subroutine NMARK(POINT, STYL)

    use read_input, only: resize

    integer, parameter :: n_symbols=7, sym_points=16
    integer, parameter, dimension(n_symbols) :: sym_size = (/16, 13, 5, 9, 14, 9, 10/)
    integer, dimension(sym_points, n_symbols) :: dx, dy

    integer, intent(in) :: POINT(2), STYL

    integer :: J, IST
    double precision, dimension(sym_points) :: xsym, ysym

! IST definition shoud coincide with NBIT() in CMARK
! IST = 1-filled diamond, 2-o, 3-+, 4-*, 5-x(#), 6-<, 7-filled square
! STYL  <=7,              8    9    10   11      12   >=13

    save DX, DY ! Actually parameters
    data DX/ &
        0,  3,  0, -3,  0,  0,  2,  0, -2,  0,  0,  1,  0, -1,  0,  0, &
        3,  3,  2,  1, -1, -2, -3, -3, -2, -1,  1,  2,  3,  0,  0,  0, &
        3, -3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
       -3, -1,  2, -2,  1, -2,  2, -1,  3,  0,  0,  0,  0,  0,  0,  0, &
       -2, -1,  0,  1,  2,  1,  0,  1,  2,  1,  0, -1, -2,  0,  0,  0, &
        0,  0, -3,  0, -2,  0,  0, -2,  0,  0,  0,  0,  0,  0,  0,  0, &
       -2,  2,  2, -2, -2,  1,  1, -1, -1,  0,  0,  0,  0,  0,  0,  0/ 
    data DY/ &
        3,  0, -3,  0,  3,  2,  0, -2,  0,  2,  1,  0, -1,  0,  1,  0, &
        1, -1, -2, -3, -3, -2, -1,  1,  2,  3,  3,  2,  1,  0,  0,  0, &
        0,  0,  0,  3, -3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
        0,  0,  3, -3,  0,  3, -3,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
       -2, -1,  0,  1,  2,  1,  0, -1, -2, -1,  0,  1,  2,  0,  0,  0, &
        3, -3,  0,  3, -1, -1,  1,  1, -3,  0,  0,  0,  0,  0,  0,  0, &
       -2, -2,  2,  2, -1, -1,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0/  

    IST = max(1, min(STYL, 7))
    do J=1, sym_size(IST)
        xsym(J) = POINT(1) + NINT(resize*dx(J, IST))
        ysym(J) = POINT(2) + NINT(resize*dy(J, IST))
    enddo
    call drawcurve(0, sym_size(IST), xsym, ysym)

    end subroutine NMARK

!---------------------------------------------------------------------
    subroutine menutable(arr_size, array_in, var_names, id)

    use char_manip, only: null_ch

    integer, intent(in) :: arr_size, id
    double precision, intent(in), dimension(arr_size) :: array_in
    character(len=6), intent(in), dimension(arr_size) :: var_names

    integer :: editable=1

    call menubox(arr_size, array_in, var_names, id, editable)

    end subroutine menutable

!---------------------------------------------------------------------
    subroutine ASKINT(NV, NVAR, NAME)

    integer, intent(in) :: NV
    integer, intent(inout) :: NVAR(NV)
    character(len=4), intent(in) :: NAME(NV)

    integer :: J
    double precision :: VAR(NV)

    do J=1, NV
        VAR(J) = NVAR(J)
    enddo
    call MENUTABLE(NV, VAR, NAME, 4)
    do J=1, NV
        NVAR(J) = VAR(J)
    enddo

    end subroutine ASKINT

!---------------------------------------------------------------------
    subroutine ASTWIN(NB, IBOX, NAME, yscale, yshift, MOD10, YMODE)

    use char_manip, only: to_upper

    integer, parameter :: NRW16=NRW+16

    integer, intent(in) :: NB, YMODE, MOD10
    integer, intent(out) :: IBOX(*)
    double precision, intent(inout), dimension(NB) :: yscale, yshift
    character(len=4), intent(out) :: NAME(*)

    integer :: JMODE, JGR, j, j1, j2, jj, jn, jm, js, jb, jw, jsep, layoutBox
    integer, dimension(NRW) :: IB
    character(len=80) :: rows(NRW16), STR, TITLE
    character(len=1) :: KEY

    TITLE = "Presentation" // null_ch
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

        do j=1, NRW
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
            write(rows(jj)(js+1:js+4), '(A4)') NAME(j)
            write(rows(jj)(js+6:js+8), '(I3)') jb
            write(rows(jj)(js+10:js+15), '(f6.1)') yscale(j)
            write(rows(jj)(js+17:js+22), '(f6.1)') yshift(j)
            jn = max(jn, jj)
        enddo

        if (jm /= 0 ) then
            j2 = jn
            j1 = 2*jn + 1
            do j=1, NB
                if (IBOX(j) <= 0 ) then
                    jw = j1
                    js = jsep*(1 - jw + jw/2*2)
                    jj = 1 + (jw - 1)/2
                    write(rows(jj)(js+1: js+4), '(A4)') NAME(j)
                    write(rows(jj)(js+6: js+8), '(I3)') IBOX(j)
                    write(rows(jj)(js+10:js+15), '(f6.1)') yscale(j)
                    write(rows(jj)(js+17:js+22), '(f6.1)') yshift(j)
                    jn = max(jn, jj)
                    j1 = j1 + 1
                endif
            enddo
        endif

        j = 1
        do while(j > 0)
            j1 = (jm + 1)/2
            j = layoutBox(TITLE, STR, rows, 80, jn, JGR, j1)
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
                write(NAME(j), '(A4)', ERR=77) rows(jj)(js+1: js+4)
                read(rows(jj)(js+6 : js+8 ), *, ERR=77) IB(j)
                read(rows(jj)(js+10: js+15), *, ERR=77) yscale(j)
                read(rows(jj)(js+17: js+22), *, ERR=77) yshift(j)
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
                    write(NAME(j), '(A4)', ERR=77) rows(jj)(js+1: js+4)
                    read(rows(jj)(js+6: js+8), *, ERR=77) IB(j)
                    read(rows(jj)(js+10:js+15), *, ERR=77) yscale(j)
                    read(rows(jj)(js+17:js+22), *, ERR=77) yshift(j)
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

    end subroutine ASTWIN

!---------------------------------------------------------------------
    subroutine ASXWIN(AB, NB, IBOX, NAME, yscale, yshift, r_min, r_max, MOD10, YMODE)

    integer, parameter :: NRW16=NRW+16

    integer, intent(in) :: NB, YMODE, MOD10
    double precision, intent(in) :: AB
    integer, intent(out) :: IBOX(*)
    double precision, intent(inout), dimension(NB) :: yscale, yshift, r_min, r_max
    character(len=4), intent(out) :: NAME(*)

    integer :: JMODE, JGR, j, j1, j2, jj, jn, jm, js, jb, jw, jsep, layoutBox
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

    YKEY = 'y'

    do while (YKEY == 'Y' .or. YKEY == 'y')

        do j=1, NRW
            write(rows(j)(1:80), '(79X, 1A1)') null_ch
        enddo

! j  - ordinal box No.
! jb - box No. in the Astra nominations
! jw - position in the table
! jj - horizontal row
! js - position in the current row

        jn = 0
        jm = 0
        write(ABNUM, '(f6.3)') AB
        do j=1, NB
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
            YY = r_max(j)
            js = jsep*(1 - jw + jw/2*2)
            jj = 1 + (jw - 1)/2
            write(rows(jj)(js+1: js+4), '(A4)') NAME(j)
            write(rows(jj)(js+6: js+8), '(I3)') jb
            write(rows(jj)(js+10: js+15), '(f6.1)') yscale(j)
            write(rows(jj)(js+17: js+22), '(f6.1)') yshift(j)
            write(rows(jj)(js+24: js+29), '(f6.3)') r_min(j)
            write(rows(jj)(js+31: js+36), '(f6.3)') r_max(j)
            if (YY > AB) rows(jj)(js+31:js+36) = ABNUM
            jn = max(jn, jj)
        enddo

        if (jm /= 0 ) then
            j2 = jn
            j1 = 2*jn+1
            do j = 1, NB
                if (IBOX(j) <= 0 ) then
                    YY = r_max(j)
                    jw = j1
                    js = jsep*(1 - jw + jw/2*2)
                    jj = 1 + (jw - 1)/2
                    write(rows(jj)(js+1: js+4), '(A4)') NAME(j)
                    write(rows(jj)(js+6: js+8), '(I3)') IBOX(j)
                    write(rows(jj)(js+10: js+15), '(f6.1)') yscale(j)
                    write(rows(jj)(js+17: js+22), '(f6.1)') yshift(j)
                    write(rows(jj)(js+24: js+29), '(f6.3)') r_min(j)
                    write(rows(jj)(js+31: js+36), '(f6.3)') r_max(j)
                    if (YY > AB) rows(jj)(js+31: js+36) = ABNUM
                    jn = max(jn, jj)
                    j1 = j1 + 1
                endif
            enddo
        endif

        j = 1
        do while(j > 0)
            j1 = (jm + 1)/2
            j = layoutBox(TITLE, STR, rows, 80, jn, JGR, j1)
        enddo

        jn = 0
        do j=1, NB
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
                read(rows(jj)(js+10: js+15), *, ERR=77) yscale(j)
                read(rows(jj)(js+17: js+22), *, ERR=77) yshift(j)
                read(rows(jj)(js+24: js+29), *, ERR=77) r_min(j)
                read(rows(jj)(js+31: js+36), *, ERR=77) YY
                if (rows(jj)(js+31: js+36) /= ABNUM) then
                    read(rows(jj)(js+31: js+36), *, ERR=77) r_max(j)
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
                    read(rows(jj)(js+10: js+15), *, ERR=77) yscale(j)
                    read(rows(jj)(js+17: js+22), *, ERR=77) yshift(j)
                    read(rows(jj)(js+24: js+29), *, ERR=77) r_min(j)
                    if  (rows(jj)(js+31: js+36) /= ABNUM) then
                        read(rows(jj)(js+31:js+36), *, ERR=77) r_max(j)
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
        read(*, '(:, A1)') YKEY
    enddo

    end subroutine ASXWIN

!---------------------------------------------------------------------
    subroutine show_plots(MARK, JIFNEW, IYO, TT_out, t_out)
!---------------------------------------------------------------------
! Drawing options:
!   X - axis
!      0 <= rho <= ROC=RHO(NA1)
!      0 <= a <= ABC
!      0 <= a <= AB
!
! MARK = 1     Put marks
! MARK = 0     Solid lines
! MARK =-1     Dashed lines
! JIFNEW =  0 Re-draw (erase) the previous curves
! JIFNEW = 1 New curves only
!---------------------------------------------------------------------

    use status, only: AMETR, MU, SHIF, ELON, TRIA
    use scalars, only: XOUT, NAB, NA1, NA1E, ABC, TINIT, TSCALE, RTOR, &
        MEQUIL, IPEQL, TIME
    use read_input, only: raw_profiles, equ_file, nr_x_max, &
        IFDFAX, NPTM, XAXES, DATAX, TOUTX
    use dbl2char, only: fmt_smart
    use char_manip, only: len_trim_tab, str_in_list
    use debugger, only: markloc, debug
    use json_vars, only: profxNames, n_profx
    use standard_functions, only: AFVAL

    integer, parameter :: jzero=0, NRDX=500
    integer, intent(in) :: MARK, JIFNEW
    integer, intent(inout) :: IYO(NTIMES,*)
    double precision, intent(in) :: TT_out(NTIMES)
    double precision, intent(inout) :: t_out(NTIMES, NRW)

    integer :: PTM(2), PTMO(2, NRDX, NRW), &
        IWN(16), fshift, &
        text_posx, text_posy, JS, MODEX, &
        IYM0, LTOUT1, LTOUT2, STYL, x_shift, y_shift, jx_canv, jy_canv, JY, jxout, &
        jplot_in_tab, j_curve, j_canv, &
        IYMN, IYMX, JDSP, NPTMO(NRW), jlx(8), &
        NP1, j, half_wid, &
        j1, jj, jsco, jn, jpnt, jsc, jarr, jtyp, n_canvas, &
        jplot_in_canv, jcol, jsym, jprof, jtrace
    double precision :: SC(NRW), YX, r_out, YA, YL, YR, &
         YZ, ymin, ymax, px_rmag, yq1, xq1, xte, te_bc
    double precision ,dimension(2) :: xbar, xbar_old, ybar, x8bar, y8bar
    double precision, dimension(16) :: xq1_old, xte_old
    double precision, dimension(NRD) :: xplot, yplot
    double precision, dimension(NTIMES) :: xtrace, ytrace, xtrace_old
    double precision, dimension(NRD, nplots_max) :: xold, yold
    double precision, dimension(NTIMES, nplots_max) :: ytrace_old
    character(len=80) :: STRI
    character(len=5 ) :: XF4
    character(len=6 ) :: CHAR6

    save PTMO, NPTMO, IWN, xq1_old, xte_old, xold, yold, xtrace_old, ytrace_old

!---------------------------------------------------------------------

    call markloc('show_plots')

    fshift = 12
    half_wid = plot_area%width/2
    YA = 0.
    if (JIFNEW /= 0) then
        do J=1, 16
            IWN(J) = 0
        enddo
    endif

! Different options for abscissa. 
! Implemented for modes 0, 1, 2, 3
    MODEX = XOUT + 0.49
    SELECT CASE(MODEX)
        CASE(0)        ! Draw up to AB against "a"
            NP1 = NAB
        CASE(1:3)      ! Draw up to ABC against "a", "rho", "Psi"
            NP1 = NA1
        CASE DEFAULT ! Unknown option
            return
    END SELECT

    if (MOD10 == 3) then
        NP1 = NA1
    endif

    IYMN = plot_area%height - plot_area%ymin
    IYM0 = plot_area%ymin - plot_area%canvas_height
    IYMX = plot_area%height - plot_area%ymax
    JY = 10*plot_area%height

    ymin = dble(plot_area%height - plot_area%ymin)
    ymax = dble(plot_area%height - plot_area%ymax)

    n_canvas = plot_area%nx_canvas*plot_area%ny_canvas

    SELECT CASE(MOD10)

!-----------------
    CASE(1: 3)  ! Profiles

        call markloc('Drawing mode 1-3', debug_lev=debug*2)
        do JPROF=1, NROUT
            do J=1, NP1
                ROUT(j, jprof) = ROUT(j, jprof) + OSHIFR(jprof)
            enddo
        enddo
        call SCAL(NROUT, SC, SCALER, ROUT(4, 1), NP1 - 3, NRD)

! Plot curves

        j_curve = 0
        plot_prof: do jprof=1, NROUT
            jplot_in_tab = NWIND1(jprof) - curves_per_frame(MOD10)*active_tab(MOD10) ! 1-16 for mode '1'
            if (NAMER(jprof) == '    ') jplot_in_tab = 0
            if (jplot_in_tab <= 0 .or. jplot_in_tab > curves_per_frame(MOD10)) CYCLE plot_prof
            j_canv = MOD(jplot_in_tab - 1, n_canvas) + 1        ! 1-8 for mode '1'
            jx_canv = MOD(j_canv - 1, plot_area%nx_canvas)      ! 0-3 for mode '1'
            jy_canv = (j_canv - 1)/plot_area%nx_canvas          ! 0-1
            jplot_in_canv = (jplot_in_tab - 1)/n_canvas
            jcol = jplot_in_canv + 2
            y_shift = plot_area%canvas_height * (plot_area%ny_canvas - 1 - jy_canv)
            x_shift = plot_area%canvas_width*jx_canv

            YL = max(0.d0, abscissa(GRAL(jprof)))
            YR = abscissa(GRAP(jprof))
            if (YL >= YR) YL = 0.d0

! Translate curve into pixel
            jxout = 0
            do j=1, NP1
                YX = abscissa(AMETR(j))
                if (YX >= YL .and. YX <= YR) then
                    jxout = jxout + 1
                    if (jxout == 1 .and. j > 1) then ! left edge interpolation
                        YA = ROUT(J, jprof) + (ROUT(J-1, jprof) - ROUT(J, jprof))*(YL - YX)/(YA - YX)
                        r_out = min(max(YA/SC(jprof), -7.d0), 7.d0)
                        JDSP  = 10*(plot_area%canvas_height*r_out + IYMN + y_shift)
                        xplot(jxout) = dble(x_shift)
                        yplot(jxout) = plot_area%height - min(max(dble(plot_area%canvas_height)*r_out + ymin + dble(y_shift), ymin), ymax)
                        jxout = jxout + 1
                    endif
                    r_out = min(max(ROUT(J, jprof)/SC(jprof), -7.d0), 7.d0)
                    JDSP  = 10*(plot_area%canvas_height*r_out + IYMN + y_shift)
                    xplot(jxout) = dble(x_shift) + dble(plot_area%width)/dble(plot_area%nx_canvas)*(YX - YL)/(YR - YL)
                    yplot(jxout) = plot_area%height - min(max(dble(plot_area%canvas_height)*r_out + ymin + dble(y_shift), ymin), ymax)
                endif
                if (YA <= YR .and. YX > YR) then ! right edge interpolation
                    jxout = jxout + 1
                    YA = ROUT(J, jprof) + (ROUT(J-1, jprof) - ROUT(J, jprof))*(YR - YX)/(YA - YX)
                    r_out = min(max(YA/SC(jprof), -7.d0), 7.d0)
                    JDSP  = 10*(plot_area%canvas_height*r_out + IYMN + y_shift)
                    xplot(jxout) = dble(x_shift) + dble(plot_area%width)/dble(plot_area%nx_canvas)
                    yplot(jxout) = dble(plot_area%height) - min(max(dble(plot_area%canvas_height)*r_out + ymin + dble(y_shift), ymin), ymax)
                endif
                YA = YX
            enddo
            if (KPRI >= 1 .and. KPRI <= 2) then
                write(STRI, '(1A6, 1A4, 1A1)') 'Plot "', NAMER(jprof), '"'
                j = len_trim_tab(STRI)
                call pscom(STRI, j)
            endif
            STYL = (jcol - 1)*MARK
            j_curve = j_curve + 1
            if (j_curve <= nplots_max) then
                call update_curve(jxout, IWN(jplot_in_tab), jcol, STYL, xold(1:, j_curve), yold(1:, j_curve), xplot, yplot)
                xold(1: jxout, j_curve) = xplot(1: jxout)
                yold(1: jxout, j_curve) = yplot(1: jxout)
            endif
            IWN(jplot_in_tab) = jxout
            do J=1, NP1
                ROUT(j, jprof) = ROUT(j, jprof) - OSHIFR(jprof)
            enddo
! Variable labels
            call setColor(jcol)
            text_posx = x_shift
            text_posy = astra_gui%dylet + FSHIFT + (plot_area%ymin - plot_area%ymax + 3 + astra_gui%dylet)*jy_canv
            call CMARK(text_posx, text_posy, SC(jprof), OSHIFR(jprof), NAMER(jprof), STYL, jplot_in_canv)

        enddo plot_prof

! Plot dots
        do j=1, n_canvas
            jlx(j) = 0
        enddo
        jsco = 0
        jn = 0

        plot_profx: do jxout=1, NXOUT
            if (NWINDX(jxout) == 0) CYCLE plot_profx ! NWINDX set in ininam.f90
            CHAR6 = NAMEX(jxout)
            if (CHAR6(1: 1) == ' ') CYCLE plot_profx
            do jprof=1, n_profx
                if (profxNames(jprof) == CHAR6) jn = jprof 
            enddo
            if (jn == 0) then
                write(*, '(/3A)') '>>> Error in the model "', TRIM(equ_file), '"'
                write(*, '(A)') '>>> Unknown data set "' // CHAR6 // '" is requested'
                ERROR STOP
            endif

            jpnt = NPTM(jn)
            if (jpnt <= 0) CYCLE plot_profx

            jsc = NWINDX(jxout)
            if (abs(SC(jsc)) < 1.1E-7) call SCAL(1, SC(jsc), SCALER(jsc), DATAX(1, jn), jpnt, nr_x_max)
            jplot_in_tab = NWIND1(jsc) - curves_per_frame(MOD10)*active_tab(MOD10)
            if (jplot_in_tab <= 0 .or. jplot_in_tab > curves_per_frame(MOD10)) CYCLE plot_profx
            j_canv = MOD(jplot_in_tab - 1, n_canvas) + 1        ! 1-8 for mode '1'
            jx_canv = MOD(j_canv - 1, plot_area%nx_canvas)      ! 0-3 for mode '1'
            jy_canv = (j_canv - 1)/plot_area%nx_canvas          ! 0-1
            jcol = (jplot_in_tab - 1)/n_canvas + 2
            y_shift = plot_area%canvas_height * (plot_area%ny_canvas - 1 - jy_canv)
            x_shift = plot_area%canvas_width*jx_canv

            call setColor(jcol)
            if (jsc /= jsco) then
                jsym = 5  ! 13
            else
                jsym = jsym + 1
            endif

            if (JIFNEW == 0) then
                call setColor(EraseColor)
                do j=1, NPTMO(jxout)
                    call NMARK(PTMO(1, j, jxout), jsym)
                enddo
            endif

            if (KPRI >= 1 .and. KPRI <= 2) then
                write(STRI, '(1A6, 1A6, 1A1)')'Dots "', profxNames(jn), '"'
                j = len_trim_tab(STRI)
                call pscom(STRI, j)
            endif
            YL = max(0.d0, abscissa(GRAL(jsc)))
            YR = abscissa(GRAP(jsc)) 
            if (YL >= YR) YL = 0.d0
            j1 = 0
            do j=1, jpnt
                YX = XAXES(j, jn)
                if (MODEX >= 1 .and. YX > ABC) CYCLE
                YA = abscissa(YX)
                if (YA > 1. .or. YA < YL .or. YA > YR) CYCLE
                j1 = j1 + 1
                PTM(1) = x_shift + plot_area%width/plot_area%nx_canvas*(YA-YL)/(YR-YL)
                r_out= max((DATAX(j, jn) + OSHIFR(jsc))/SC(jsc), -7.d0)
                r_out= min(r_out, 7.d0)
                JDSP = plot_area%canvas_height*r_out + IYMN + y_shift
                PTM(2) = plot_area%height - min(max(JDSP, IYMN), IYMX)
                call setColor(jcol)
                call NMARK(PTM, jsym)
                PTMO(1, j1, jxout) = PTM(1)
                PTMO(2, j1, jxout) = PTM(2)
            enddo
            NPTMO(jxout) = j1

            if (KPRI >= 1 .and. KPRI <= 2) then
                write(STRI, '(1A)')'Mark & time for dots'
                j = len_trim_tab(STRI)
                call pscom(STRI, j)
            endif
            jlx(j_canv) = jlx(j_canv) + 1
            text_posx = x_shift + plot_area%canvas_width - astra_gui%dxlet - 45
            text_posy = (1 + jlx(j_canv))*astra_gui%dylet + FSHIFT + (plot_area%ymin - plot_area%ymax - plot_area%canvas_height)*(jy_canv) + 3
            XF4 = fmt_smart(TOUTX(jn), 4)
            call textvm(text_posx, text_posy, XF4, 5) ! Text (time) -> plot legend
            PTM(1) = text_posx + astra_gui%dxlet + 37 ! 12 is fixed, as the font size does not scale
            PTM(2) = text_posy - 0.3*astra_gui%dylet
            jcol = (jplot_in_tab - 1)/n_canvas + 2
            call NMARK(PTM, jsym) ! Marker symbol -> plot legend

            jsco = jsc
            if (jsym == 7) jsym = 1
        enddo plot_profx

! Erase/put q=1 radius, BC for Te
        yq1   = abscissa(AFVAL(MU, 1.d0))
        te_bc = abscissa(AMETR(max(NA1E, 1)))
        ymax = dble(IYM0) - 0.8*plot_area%canvas_height
        do j_canv=1, plot_area%nx_canvas
            xq1 = plot_area%canvas_width*(j_canv -1 + YQ1)
            xte = plot_area%canvas_width*(j_canv -1 + te_bc)
            do jy=1, plot_area%ny_canvas
                ybar = (/ dble(IYM0) - (jy-2)*plot_area%canvas_height, ymax - (jy-2)*plot_area%canvas_height/)
                if (yq1 > 1.d-3 .and. yq1 < 0.999) then
                    xbar_old = xq1_old(j_canv)
                    xbar = xq1
                    call update_curve(2, 2,   Red, 0, xbar_old, ybar, xbar, ybar)
                endif
                if (te_bc > 1.d-2 .and. te_bc < 0.999) then
                    xbar_old = xte_old(j_canv)
                    xbar = xte
                    call update_curve(2, 2, Green, 0, xbar_old, ybar, xbar, ybar)
                endif
            enddo
            xq1_old(j_canv) = xq1
            xte_old(j_canv) = xte
        enddo

!-------------------
    CASE(4:5)
        call markloc('Drawing mode 4/5', debug_lev=2*debug)

!-------------------
    CASE(6)  ! Time traces

        call markloc('Drawing mode 6', debug_lev=2*debug)
        LTOUT1 = 1
        LTOUT2 = LTOUT
        if (LTOUT < 2) return
        ! right_label_position=JDX*JDMX=23*5*5=575 (see typdsp.f)
        do J=1, LTOUT-1
            r_out = (TT_out(J) - TINIT)*575/abs(TSCALE)
            IYO(J, nplots_max+1) = 6*astra_gui%dxlet + r_out
            xtrace(J) = 6*astra_gui%dxlet + r_out
            if (r_out < 0)  LTOUT1 = J + 1
            if (IYO(J, nplots_max+1) <= plot_area%width - 1) LTOUT2 = J - 1
        enddo
        LTOUT2 = LTOUT2 - LTOUT1 + 1
        if (LTOUT2 < 3) then ! No plot for small time
            return
        endif
        do jj=1, min(NTOUT, NRW)
            do j=LTOUT1, LTOUT1 + LTOUT2
                t_out(j, jj) = t_out(j, jj) + OSHIFT(jj)
            enddo
        enddo
        call SCAL(NTOUT, SC, SCALET, t_out(LTOUT1, 1), LTOUT2, NTIMES)

        j_curve = 0
        plot_traces: do jtrace=1, min(NRW, NTOUT)
            jplot_in_tab = NWIND3(jtrace) - curves_per_frame(MOD10)*active_tab(MOD10)
            if (NAMET(jtrace) == '    ') CYCLE plot_traces
            if (jplot_in_tab <= 0 .or. jplot_in_tab > curves_per_frame(MOD10)) CYCLE plot_traces
            jplot_in_canv = (jplot_in_tab - 1)/n_canvas       ! <-> color
            j_canv = MOD(jplot_in_tab - 1, n_canvas) + 1      ! 1-8 for mode '1'
            do J=1, LTOUT
                r_out = max(t_out(J, jtrace)/SC(jtrace), -7.d0)
                r_out = min(r_out, 7.d0)
                JDSP  = 10*(plot_area%canvas_height*r_out + IYMN + (n_canvas - j_canv)*plot_area%canvas_height)
                JDSP  = max(JDSP, 10*IYMN)
                IYO(J, nplots_max+2) = JY - min(JDSP, 10*IYMX)
                ytrace(J) = dble(plot_area%height) - min(max(plot_area%canvas_height*r_out + ymin + (n_canvas - j_canv)*plot_area%canvas_height, ymin), ymax)
            enddo

            if (KPRI >= 1 .and. KPRI <= 2) then
                write(STRI, '(1A6, 1A4, 1A1)') 'Plot "', NAMET(jj), '"'
                j = len_trim_tab(STRI)
                call pscom(STRI, j)
            endif
            jcol = jplot_in_canv + 1 
            STYL = jcol*MARK
            jcol  = jcol  + 1
            j_curve = j_curve + 1
            if (j_curve <= nplots_max) then
                call setColor(White)
                call plot_curve(LTOUT2, STYL, xtrace_old(1:LTOUT2), ytrace_old(1:LTOUT2, j_curve))
                call setColor(jcol)
                call plot_curve(LTOUT2+1, STYL, xtrace(1:LTOUT2+1), ytrace(1:LTOUT2+1) )
                xtrace_old(1: LTOUT2+1) = xtrace(1: LTOUT2+1)
                ytrace_old(1: LTOUT2+1, j_curve) = ytrace(1: LTOUT2+1)
            endif
            text_posx = 0
            text_posy = (j_canv-1)*plot_area%canvas_height + astra_gui%dylet*(jplot_in_canv*2 + 2) + 2*jplot_in_canv
            call CMARKT(text_posx, text_posy, SC(jtrace), OSHIFT(jtrace), NAMET(jtrace), STYL)

            do j=LTOUT1, LTOUT1 + LTOUT2
                t_out(j, jtrace) = t_out(j, jtrace) - OSHIFT(jtrace)
            enddo
        enddo plot_traces

!-------------------
    CASE(7)
        call markloc('Drawing mode 7', debug_lev=2*debug)

!-------------------
    CASE(8)

        call markloc('Drawing mode 8', debug_lev=2*debug)
        px_rmag = (RTOR + SHIF(1))*meter2pixel
        call setColor(Black)
        x8bar = (/0.d0, px_rmag/)
        y8bar = (/pixel_ymid, pixel_ymid/)
        call drawcurve(0, 2, x8bar, y8bar)
! Plot the complete wall structure (Pixmap # 1)
        call plot_wall()
! if (data file includes NAMEXP BND) then (n_bnd > 0);
! or (n_bnd == 8) after calling equil with no boundary points provided;
!     n_bnd == 0 otherwise

        SELECT CASE(IPEQL)
        CASE(3)
            call plot_lcfs(JIFNEW, IYO, TIME)
        CASE(4: 5)
            if (MEQUIL == 0) then
                SHIF = 0.0
                ELON = 1.0
                TRIA = 0.0 
            else
                call plot_lcfs(JIFNEW, IYO, TIME)
            endif
        END SELECT
        call plot_flux_surfaces()
! Plot dots of R, z if an exp profile is input as function of R, z (GRIDTYPE 18-20)
        jcol = -1
        loop8: do jxout=1, NXOUT
            if (NWINDX(jxout) == 0) CYCLE loop8
            CHAR6 = NAMEX(jxout)
            if (CHAR6(1:1) == ' ') CYCLE loop8
            if (jxout > 1) then
                do jj=jxout-1, 1, -1
                    if (NAMEX(jj) == CHAR6) CYCLE loop8 ! double name?
                enddo
            endif

            jn = str_in_list(CHAR6, profxNames)
        
            jarr = IFDFAX(jn)
            if (jarr <= 0)  CYCLE loop8
            jtyp = raw_profiles%grid_type(jarr)
            if (jtyp < 18)  CYCLE loop8
            jpnt = raw_profiles%nrho(jarr)
            if (jpnt <= 0)  CYCLE loop8
            jcol = jcol+1
            js = raw_profiles%jbeg_grid(jarr)

            if (JIFNEW == 0) then
                call setColor(EraseColor)
                do j=1, NPTMO(jxout)
                    call NMARK(PTMO(1, j, jxout), 7)
                enddo
            endif
            call setColor(jcol+2)
            do j=1, jpnt
                if (jtyp == 18)  then
                    YR = raw_profiles%data(js)
                    YZ = raw_profiles%data(js+j)
                elseif (jtyp == 19)  then
                    YR = raw_profiles%data(js+j)
                    YZ = raw_profiles%data(js)
                elseif (jtyp == 20) then
                    YR = raw_profiles%data(js-1+j)
                    YZ = raw_profiles%data(jpnt+js-1+j)
                else
                    write(*, *) 'Unknown input-grid type'
                endif
                PTM(1) = YR*meter2pixel
                PTM(2) = pixel_ymid - YZ*meter2pixel
                call NMARK(PTM, 7)
                PTMO(1, j, jxout) = PTM(1)
                PTMO(2, j, jxout) = PTM(2)
            enddo
! git hardcoded 510
            call textvm(510, 100 + 2*astra_gui%dylet*jcol, NAMEX(jxout), 6)
            NPTMO(jxout) = jpnt
        enddo loop8

    END SELECT

    end subroutine show_plots

!---------------------------------------------------------------------
    subroutine plot_lcfs(ifnew, IYO, time_in)

! Scatter plot of the LCFS
! IFNEW  =  0 Re-draw (erase) the previous curves
! IFNEW =/= 0 New curves only
! IFNEW < 0 Don't mark resonances q=m/n
! IFNEW > 10 Call from Review. (JIFNEW-10) is used to control erasing

    use read_input, only: raw_boundary

    integer, intent(in) :: ifnew
    integer, intent(inout) :: IYO(2, *)
    double precision, intent(in) :: time_in

    integer :: j, j1, j2, jj, PTM(2)
    double precision :: YS, YX, YXL, YXR, YZ

    j2 = 1 + raw_boundary%n_theta/32

    if (IFNEW == 0) then
        call setColor(EraseColor)
        do j=1, raw_boundary%n_theta, j2
            PTM(1) = IYO(1, j)
            PTM(2) = IYO(2, j)
            call NMARK(PTM, 4)
        enddo
    endif

! Boundary points
    call setColor(Red)

    if (raw_boundary%nt <= 1) then
        jj = 1
        do j=1, raw_boundary%n_theta, j2
            j1 = max(1, raw_boundary%nt + (j - 1)*jj)
            PTM(1) = raw_boundary%R(j1)*meter2pixel
            PTM(2) = pixel_ymid - raw_boundary%Z(j1)*meter2pixel
            call NMARK(PTM, 4)   !Use (PTM, 4) for *
            IYO(1, j) = PTM(1)
            IYO(2, j) = PTM(2)
        enddo
    else if (time_in <= raw_boundary%time(1) .or. time_in >= raw_boundary%time(raw_boundary%nt)) then ! extrapolate flat
        jj = raw_boundary%nt
        do j=1, raw_boundary%n_theta, j2
            j1 = raw_boundary%nt + (j - 1)*jj
            PTM(1) = raw_boundary%R(j1)*meter2pixel
            PTM(2) = pixel_ymid - raw_boundary%Z(j1+jj)*meter2pixel
            call NMARK(PTM, 4)   !Use (PTM, 4) for *
            IYO(1, j) = PTM(1)
            IYO(2, j) = PTM(2)
        enddo
    else               ! interpolate linearly
        do j=1, raw_boundary%nt   ! Find current time
            if (time_in > raw_boundary%time(j)) jj = j
        enddo
        YS  = raw_boundary%time(jj+1) - raw_boundary%time(jj)
        YXL = (time_in - raw_boundary%time(jj  ))/YS
        YXR = (time_in - raw_boundary%time(jj+1))/YS
        do j=1, raw_boundary%n_theta, j2 ! Time differentiation
            j1 = jj + (j - 1)*raw_boundary%nt
            YX = YXL*raw_boundary%R(j1+1) - YXR*raw_boundary%R(j1)
            YZ = YXL*raw_boundary%Z(j1+1) - YXR*raw_boundary%Z(j1)
            PTM(1) = YX*meter2pixel
            PTM(2) = pixel_ymid - YZ*meter2pixel
            call NMARK(PTM, 4)
            IYO(1, j) = PTM(1)
            IYO(2, j) = PTM(2)
        enddo
    endif

    end subroutine plot_lcfs

!---------------------------------------------------------------------
    subroutine plot_wall()
! Plot vessel components reading them from json machine file

    use pi_const, only: GP2
    use scalars, only: AB, ELONM, RTOR, TRICH
    use machine_config, only: config, json_cfg, cfg_exists

    integer, parameter :: n_theta=64, ngc_max=750
    integer :: j, j1, jgc, jbeg, NGC, ndim_gc
    integer, allocatable, dimension(:) :: contour_len, contour_color
    double precision :: pol_ang, Rwall, Zwall
    double precision, dimension(n_theta) :: xwall, ywall
    double precision, dimension(ngc_max) :: xGC, yGC
    double precision, allocatable, dimension(:) :: rGC, zGC

    call setColor(Blue)

    if (cfg_exists) then
        call config%get('contour_len', contour_len)
        call config%get('contour_color', contour_color)
        call config%get('Rvessel', rGC)
        call config%get('Zvessel', zGC)
        NGC = SIZE(contour_len)
        jbeg = 0
        do j1=1, NGC
            if (contour_color(j1) /= White) then
                ndim_gc = contour_len(j1)
                do j=1, ndim_gc
                    xgc(j) = meter2pixel*rGC(jbeg+j)
                    ygc(j) = pixel_ymid - meter2pixel*zGC(jbeg+j)
                enddo
                call plot_curve(ndim_gc, 0, xgc(1:ndim_gc), ygc(1:ndim_gc))
            endif
            jbeg = jbeg + contour_len(j1) 
        enddo
    else
        do j=1, n_theta
            pol_ang = GP2*(j - 1)/float(n_theta)
            Zwall = AB*ELONM*SIN(pol_ang)
            Rwall = RTOR + AB*(COS(pol_ang) + 0.5*TRICH*(COS(2.*pol_ang) - 1.))
            xwall(j) = Rwall*meter2pixel
            ywall(j) = pixel_ymid - Zwall*meter2pixel
        enddo
        call plot_curve(jgc, 0, xwall, ywall)
        write(*, *) '>>> plot_wall: problems opening file ' // TRIM(json_cfg)
    endif

    end subroutine plot_wall

!---------------------------------------------------------------------
    double precision function abscissa(YIN)
! Input: MODEX, YIN, FP
! Output: Value a=YIN mapped to the current abscissa

    use status, only: AMETR, FP_NORM
    use scalars, only: XOUT, AB, ABC, ROC, NA1
    use numerical_tools, only: QUADIN
    use standard_functions, only: RFA

    double precision, intent(in) :: YIN

    integer :: MODEX
    double precision :: YAB

    MODEX = XOUT + 0.49

    SELECT CASE(MODEX)
    CASE(0) ! Draw up to AB against "a"
        YAB = AB
    CASE(1) ! Draw up to ABC against "a"
        YAB = ABC
    CASE(2) ! Draw up to ROC against "rho"
        YAB = ROC
    CASE(3) ! Draw up to ROC against "Psi"
        YAB = 1.
    CASE DEFAULT ! Unknown option
        return
    END SELECT

    if (MOD10 == 3 .or. MODEX == 3) then
        abscissa = QUADIN(NA1, AMETR, FP_NORM, YIN)
        return
    endif

    if (MODEX == 0 .or. MODEX == 1) then
        abscissa = min(1.d0, max(0.d0, YIN/YAB))
    elseif (MODEX == 2) then
        abscissa = RFA(YIN)/YAB
    else
        write(*, *) "MODEX is neither 0, nor 1, nor 2, it cannot be"
    endif

    end function abscissa

!---------------------------------------------------------------------
    subroutine plot_flux_surfaces()
! Redraw magnetic surfaces:

    use scalars, only: NEQUIL, MEQUIL
    use parameters_a2equil, only: equil_now

    integer, parameter :: n_surf=556, nrho_plot=12
    integer :: jrho, jr, nskip, n_theta, n_theta1, n_rho_surf
    double precision, dimension(n_surf) :: xplot, yplot
    double precision, dimension(nrho_plot+1, n_surf) :: xplot_old, yplot_old

    save xplot_old, yplot_old

    if (SIZE(equil_now%coord_sys%position%r) == 0) return

    n_rho_surf = NEQUIL
    n_theta    = MEQUIL
    n_theta1 = n_theta + 1
    nskip = 1 + n_rho_surf/nrho_plot

    jr = 1
    do jrho=1, n_rho_surf + nskip - 1, nskip
        if (jrho > n_rho_surf) EXIT
        xplot(1: n_theta) = meter2pixel*equil_now%coord_sys%position%r(jrho, 1: n_theta)
        xplot(n_theta1)   = meter2pixel*equil_now%coord_sys%position%r(jrho, 1)  ! Close polygon
        yplot(1: n_theta) = pixel_ymid - meter2pixel*equil_now%coord_sys%position%z(jrho, 1:n_theta)
        yplot(n_theta1)   = pixel_ymid - meter2pixel*equil_now%coord_sys%position%z(jrho, 1)
        call update_curve(n_theta1, n_theta1, Magenta, 0, xplot_old(jr, 1:n_theta1), &
             yplot_old(jr, 1:n_theta1), xplot(1:n_theta1), yplot(1:n_theta1))
        xplot_old(jr, 1:n_theta1) = xplot(1:n_theta1)
        yplot_old(jr, 1:n_theta1) = yplot(1:n_theta1)
        jr = jr + 1
    enddo

    end subroutine plot_flux_surfaces

!---------------------------------------------------------------------
    integer function plotMode(mod_10, mode_y)

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

    end subroutine set_plot_area

!---------------------------------------------------------------------
    subroutine set_plot(plot_mode)
    ! Subroutine draw frame for different modes

    use scalars, only: TSCALE, TINIT, AWALL
    use dbl2char, only: fmt_smart
    use char_manip, only: len_trim_tab

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
        j = astra_gui%Width - 20*astra_gui%dxlet + 1
        jj = plot_area%height + astra_gui%dylet
        call setColor(Black)
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
                CH6 = fmt_smart(TIND, 5)
                call textvm(JX - 2, JJ, CH6, 6)
            else
                XF4 = fmt_smart(TIND, 4)
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
            XF4 = fmt_smart(YY, 4)
            call textvm(JX, JJ, XF4, 5)
        enddo
    ! vertical axis labels
        do J=-1, 1
            JX = (plot_area%ymin + plot_area%ymax + astra_gui%dylet)/2 + IDT*IDX*J - 0.5*astra_gui%dylet
            YY = -J*scale_bnd
            XF4 = fmt_smart(YY, 4)
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

    end subroutine set_plot

end module graph_utils
