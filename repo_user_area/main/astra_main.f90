program astra

use parameter_inc, only: NSBMX, NRD
use outcmn_inc, only: VCOIL, CCOIL, CCOILO, outcmn_init, &
    XWH, XWW, XWX, XWY, DXLET, DYLET, LRJJ, frame_wid, frame_hei, &
    TASK, machine, exp_file, equ_file, rev_file, AWD, &
    COLTAB, RUNID, NST, MOD10, NTOUT
use const_inc, only: IPART, const_init, XOUT, NA
use status_inc, only: status_init, MU
use debugger, only: debug, flightsim

implicit none

!-------------------------------------------
! Find self-consistent initial configuration
!-------------------------------------------

integer :: j, jj, IM, ios, XSC0, XSC
double precision :: CHORDN, resize, tbeg_nml, tend_nml, tpause_nml
character(len=64) :: LISTSB(NSBMX)
character(len=132) :: file_in, STRI, win_title
double precision, external :: LINEAV

namelist / astra_log / AWD, exp_file, equ_file, rev_file, TASK, machine, &
debug, tbeg_nml, tend_nml, tpause_nml, flightsim, resize

call fenvex()   !  Enable floating exception handling

!-------------------- Initial settings --------------------------------|

call outcmn_init
call const_init
call status_init

debug = 0 ! Initialise to: no debugging

call ininam(LISTSB)
IPART = 1   ! Mark initial iteration section

call read_input

resize = 1.

!--------------------
! ASTRA graphic frame
!--------------------

call GETENV('expfile', exp_file)
call GETENV('equfile', equ_file)

file_in = 'tmp/' // TRIM(exp_file) // TRIM(equ_file) // '.nml'

OPEN(161, FILE=TRIM(file_in), delim='apostrophe')
READ(161, nml=astra_log, iostat=ios)
CLOSE(161)

if (TASK(1: 3) == 'BGD') then
    STRI = 'BGD'//char(0)
    call initvm(XWX, XWY, XWW, XWH, COLTAB, STRI(1: 3), 3)
else
    call get_runid()
! Resize
    frame_wid = resize*frame_wid
    frame_hei = resize*frame_hei
    XWW   = resize*XWW
    XWH   = resize*XWH
    XWX   = resize*XWX
    XWY   = resize*XWY
    DXLET = resize*DXLET
    DYLET = resize*DYLET
    LRJJ  = resize*LRJJ
    jj = max(0, (15 + NTOUT - 64)/16)
    XWH = XWH + 2*jj*(DYLET + 2)
    win_title = 'Per aspera ad ASTRA'
    call initvm(XWX, XWY, XWW, XWH, COLTAB, TRIM(win_title), LEN_TRIM(win_title)) ! Initialise graphic window

    IM = 1
    NST = 0
    MOD10 = 1
    call set_frame(IM, XSC0, XSC)
    call set_plot(IM, XSC0, XSC)

    j = XOUT + 0.49

    call ASRUMN(j) ! Task menu
    call textbf(0, XWH-104, RUNID, 80) ! Task ID

    CHORDN = LINEAV()
    call UPSTR(CHORDN, 1./MU(NA))
endif

call SETARX(1)
call INIVAR
call SETVAR
call DETVAR_INIT
call EQGUESS
call INIVAR

call CONVERGE_INIT(LISTSB)

! write output file for simulink or whatever control system
if (flightsim == 1) then
    call write_output_diag_file
endif

!---------------
! Time step loop
!---------------

do
    call STEPUP 
enddo

end program astra
