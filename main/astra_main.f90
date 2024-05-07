program astra

use parameter_inc, only: NSBMX, NRD
use outcmn_inc, only: VCOIL, CCOIL, CCOILO, outcmn_init, &
    Xwin_height, Xwin_width, Xwin_xpos, Xwin_ypos, &
    DXLET, DYLET, LRJJ, frame_wid, frame_hei, TASK, &
    RUNID, NST, MOD10, NTOUT, LineWidth, resizeGraph, null_ch
use const_inc, only: IPART, const_init, XOUT, NA, &
    TIME, TEND, DPOUT, TAU, ATREQ, IFBEY, NITOT
use status_inc, only: status_init, MU, defarr
use debugger, only: debug, astra_stop, markloc
use ext_bnd, only: use_ext_bnd
use transport2fbe, only: transport2fbe_init

implicit none

!-------------------------------------------
! Find self-consistent initial configuration
!-------------------------------------------

integer :: j, jj, IM, ios, XSC0, XSC, jt1, jt2, jt3, jt_req, jkey
double precision :: CHORDN, Y, timeb
character(len=64) :: LISTSB(NSBMX)
character(len=132) :: file_in, STRI, Xwin_title
double precision, external :: LINEAV, SWATCH
integer, external :: IFKEY, IFTREQ

call fenvex()   !  Enable floating exception handling

!-------------------- Initial settings --------------------------------|

call outcmn_init
call const_init
call status_init
use_ext_bnd = 0

debug = 0 ! Initialise to: no debugging

call ininam(LISTSB)
IPART = 1   ! Mark initial iteration section

call read_input

!--------------------
! ASTRA graphic frame
!--------------------

if (TASK(1: 3) /= 'BGD') then
    call get_runid()
! Resize
    frame_wid = resizeGraph*frame_wid
    frame_hei = resizeGraph*frame_hei
    DXLET = resizeGraph*DXLET
    DYLET = resizeGraph*DYLET
    LRJJ  = resizeGraph*LRJJ
    jj = max(0, (15 + NTOUT - 64)/16)
    Xwin_width  = resizeGraph*Xwin_width
    Xwin_height = resizeGraph*Xwin_height + 2*jj*(DYLET + 2)
    LineWidth = int(0.8*resizeGraph) + 1
    Xwin_title = 'Per aspera ad ASTRA' // null_ch
    call initvm(Xwin_xpos, Xwin_ypos, Xwin_width, Xwin_height, LineWidth, TRIM(Xwin_title), LEN_TRIM(Xwin_title)) ! Initialise graphic window

    IM = 1
    NST = 0
    MOD10 = 1
    call set_frame(IM, XSC0, XSC)
    call set_plot(IM, XSC0, XSC)

    j = XOUT + 0.49

    call taskmenu(j) ! Task menu
    call textbf(0, Xwin_Width-104, RUNID, 80) ! Task ID

    CHORDN = LINEAV()
    call up_label(CHORDN, 1./MU(NA))
endif

call SETARX(1)
call INIVAR
call SETVAR
call DETVAR_INIT
call EQGUESS
call INIVAR

call transport2fbe_init

jt_req = 0
do while (jt_req == 0) ! Till convergence (jt_req /= 0). Max #iterations is set in IFTREQ (for/defarr.f90)
    if (TASK(1:3) /= 'BGD') jkey = IFKEY(256) 
    call INTVAR      ! Set exp scalars
    call DETVAR_INIT
    call DEFARR
    call SETARX(1)   ! Set X-data w/o time interpolation
    call INIVAR
    call markloc("init")
    NITOT = NITOT + 1

    call INIT_CONVERGE_STEP(LISTSB)
    call markloc("init done")
    
    IFBEY = 0. ! no fbe possible here
    call METRIC
    jt_req = IFTREQ(ATREQ)     ! ++ITREQ; Convergence check; 1 - yes
enddo

!---------------
! Time step loop
!---------------

do while (TIME - TEND + 1.E-8 < DPOUT + TAU)
    call STEPUP 
enddo

timeb = swatch(Y)
jt1 = timeb
jt2 = jt1/3600
jt3 = (jt1 - 3600*jt2)/60
jt1 = timeb - 60*jt3 - 3600*jt2
write(6, '(A, I4.2, 2(A1, I2.2))') '>>> ASTRA normal exit >>>  Run time', jt2, ':', jt3, ':', jt1
call CPUSE(6)
call astra_stop

end program astra
