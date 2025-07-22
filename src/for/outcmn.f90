module outcmn_inc

use parameter_inc, only: NRD, NRDX, NRW, NTARR, NARRX, NSBMX, NCNBM, NCNBTM, &
    NCONST, NSDELOUT, plot_modes

implicit none

type astra_xwindow
    integer :: Width, Height, Xpos, Ypos, dxlet, dylet, LineWidth, yMessage
    double precision :: resizeGraph
    character(len=128) :: title='Per aspera ad ASTRA'//char(0)
endtype astra_xwindow

type plot_frame
    integer :: width, height, xmin, xmax, ymin, ymax, nx_canvas, ny_canvas, canvas_height, canvas_width
endtype plot_frame

! Colors, array AstraColorNum in Astra2XW.c
integer, parameter :: White=0, Black=1, Red=2, Blue=3, Green=5, &
     WarningColor=30, EraseColor=31, Magenta=14, Pink=13, ICVMX=32
character(len=1), parameter :: null_ch=char(0), tab_ch=char(9), esc_ch=char(13), backslash=char(92)

integer, dimension(NRW)   :: MARKT, MARKR, NWIND1, NWIND3, NWIND4, NWIND7, NWINDX, IP1, IP2, IP30, IP31
integer, dimension(NARRX) :: IFDFAX, jbeg_arrx, NPTM
integer, dimension(NSBMX) :: IFSBX, IFSBP
integer :: &
    NDTNAM, NTOUT, NROUT, NSBR, NSBP, &
    NBNT, NCNBT, LTOUT, IPOUT, MOD10, NGR, NXOUT, IFDFVX(NCONST)
integer :: MODEY, IDX, IDT, KPRI, NST, AVERS, ARLEAS, AEDIT
integer, dimension(plot_modes) :: active_tab, curves_per_frame
double precision, dimension(NRW)   :: GRAL, GRAP, OSHIFT, OSHIFR, SCALET, SCALER
double precision, dimension(NARRX) :: TOUTX
double precision, dimension(NCNBM) :: CCOIL, VCOIL
double precision, dimension((NCNBM+1)*NCNBTM) :: CCOILX, VCOILX
double precision, dimension(NRDX, NARRX) :: XAXES, DATAX
double precision, dimension(NRD, NRW) :: ROUT
double precision :: TIM7(4), scale_bnd, pixel_ymid, meter2pixel, resizeGraph
double precision :: cpu_start, cpuTime_equ, cpuTime_tra, cpuTime_sbr(NSBMX)
integer :: wall_start

character(len=4), dimension(NRW) :: NAMET, NAMER
character(len=4) :: TASK, machine
character(len=6), dimension(NRW) :: NAMEX
character(len=6) :: DTNAME(NSDELOUT+4*NSBMX), NAM7(4)
character(132) :: nml_file, exp_file, equ_file, rev_file, TASKID, NBFILE, VERSION, RUNID, AWD
type(astra_xwindow) :: astra_gui_ref, astra_gui
type(plot_frame) :: plot_area_ref, plot_area

contains

subroutine outcmn_init

integer :: i, j

! Constants

cpuTime_equ = 0.
cpuTime_sbr = 0.
pixel_ymid  = 0.
meter2pixel = 0.
resizeGraph = 1.

VERSION = repeat(' ', 32)
NBFILE = '***'

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

NWINDX = 0
MARKT  = 0
MARKR  = 0
MODEY  = 1
LTOUT  = 1
IPOUT  = 1

! Boundary

NBNT = 0

active_tab = 0

! Coils
NCNBT = 0

OSHIFT = 0.
OSHIFR = 0.
GRAL   = 0.

NAMEX(:) = '      '
NAM7 = (/ 'Tmin', 'Tmax', 'Tmark', 'Style' /)

IFDFVX = -1
IFDFAX = -1
TIM7 = (/ 0, 9999, 9999, 1 /)

! Coils

VCOILX = 0.
CCOILX = 0.
VCOIL  = 0.
CCOIL  = 0.

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

IP1 = (/ &
  1,   3,   5,   7,   9,  11,  13,  15,   2,   4,   6,  8,   10,  12,  14,  16, &
 17,  19,  21,  23,  25,  27,  29,  31,  18,  20,  22, 24,   26,  28,  30,  32, &
 33,  35,  37,  39,  41,  43,  45,  47,  34,  36,  38, 40,   42,  44,  46,  48, &
 49,  51,  53,  55,  57,  59,  61,  63,  50,  52,  54, 56,   58,  60,  62,  64, &
 65,  67,  69,  71,  73,  75,  77,  79,  66,  68,  70, 72,   74,  76,  78,  80, &
 81,  83,  85,  87,  89,  91,  93,  95,  82,  84,  86, 88,   90,  92,  94,  96, &
 97,  99, 101, 103, 105, 107, 109, 111,  98, 100, 102, 104, 106, 108, 110, 112, &
113, 115, 117, 119, 121, 123, 125, 127, 114, 116, 118, 120, 122, 124, 126, 128 /)

IP2 = (/ &
  1,   2,   9,  10,   3,   4,  11,  12,   5,   6,  13,  14,   7,   8,  15,  16, &
 17,  18,  25,  26,  19,  20,  27,  28,  21,  22,  29,  30,  23,  24,  31,  32, &
 33,  34,  41,  42,  35,  36,  43,  44,  37,  38,  45,  46,  39,  40,  47,  48, &
 49,  50,  57,  58,  51,  52,  59,  60,  53,  54,  61,  62,  55,  56,  63,  64, &
 65,  66,  73,  74,  67,  68,  75,  76,  69,  70,  77,  78,  71,  72,  79,  80, &
 81,  82,  89,  90,  83,  84,  91,  92,  85,  86,  93,  94,  87,  88,  95,  96, &
 97,  98, 105, 106,  99, 100, 107, 108, 101, 102, 109, 110, 103, 104, 111, 112, &
113, 114, 121, 122, 115, 116, 123, 124, 117, 118, 125, 126, 119, 120, 127, 128/)

IP30 = (/ &
  1,   3,   5,   7,   2,   4,   6,   8,   9,  11,  13,  15,  10,  12,  14,  16, &
 17,  18,  21,  22,  25,  26,  29,  30,  19,  20,  23,  24,  27,  28,  31,  32, &
 33,  34,  37,  38,  41,  42,  45,  46,  35,  36,  39,  40,  43,  44,  47,  48, &
 49,  50,  53,  54,  57,  58,  61,  62,  51,  52,  55,  56,  59,  60,  63,  64, &
 65,  66,  69,  70,  73,  74,  77,  78,  67,  68,  71,  72,  75,  76,  79,  80, &
 81,  82,  85,  86,  89,  90,  93,  94,  83,  84,  87,  88,  91,  92,  95,  96, &
 97,  98, 101, 102, 105, 106, 109, 110,  99, 100, 103, 104, 107, 108, 111, 112, &
113, 114, 117, 118, 115, 116, 119, 120, 121, 122, 125, 126, 123, 124, 127, 128 /)

IP31 = (/ &
  1,   2,   5,   6,   3,   4,   7,   8,   9,  10,  13,  14,  11,  12,  15,  16, &
 17,  18,  21,  22,  19,  20,  23,  24,  25,  26,  29,  30,  27,  28,  31,  32, &
 33,  34,  37,  38,  35,  36,  39,  40,  41,  42,  45,  46,  43,  44,  47,  48, &
 49,  50,  53,  54,  51,  52,  55,  56,  57,  58,  61,  62,  59,  60,  63,  64, &
 65,  66,  69,  70,  67,  68,  71,  72,  73,  74,  77,  78,  75,  76,  79,  80, &
 81,  82,  85,  86,  83,  84,  87,  88,  89,  90,  93,  94,  91,  92,  95,  96, &
 97,  98, 101, 102, 105, 106, 109, 110,  99, 100, 103, 104, 107, 108, 111, 112, &
113, 114, 117, 118, 115, 116, 119, 120, 121, 122, 125, 126, 123, 124, 127, 128 /)

DTNAME(1: NSDELOUT) = (/ &
    'dRout ', 'dTout ', 'dPout ', 'Time  ', 'TAUmin', 'TAUmax', &
    'TAUinc', 'DELvar', 'Iterex', 'NiTrEq', 'Tinit ', 'Tscale', &
    'NA1   ', 'NUF   ', 'Xaxis ', 'Xdeflt', 'NB2EQL', 'NEQUIL', &
    'NBND  ', 'Xflag ', 'DTeql ', 'MEQUIL', 'Tpause', 'Tend  ', &
    'Inume1', 'Inume2', 'Inume3', 'Inume4', 'Iprot ', 'Itfbe ', &
    'Itfbp ', 'Icircq', 'Ipctrl', 'Adcmpf', 'Flxdr ', 'Sgnip ', &
    'Sgnbt ', 'Ifbeg ', 'Ipeql ' /)
do j=1, 30
    i = (j-1)*4 + NSDELOUT
    write(DTNAME(i+1), '(A, i0)') 'DTeq', j
    write(DTNAME(i+2), '(A, i0)') 'BEeq', j
    write(DTNAME(i+3), '(A, i0)') 'ENeq', j
    write(DTNAME(i+4), '(A, i0)') ' Keq', j
enddo

end subroutine outcmn_init

end module outcmn_inc
