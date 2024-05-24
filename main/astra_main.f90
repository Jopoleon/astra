program astra

use parameter_inc, only: NSBMX, NRD
use outcmn_inc, only: astra_gui, astra_gui_ref, VCOIL, CCOIL, CCOILO, outcmn_init, TASK
use const_inc, only: IPART, const_init, &
    TIME, TEND, DPOUT, TAU, ATREQ, IFBEY, NITOT
use status_inc, only: status_init, defarr
use debugger, only: debug, astra_stop, markloc
use ext_bnd, only: use_ext_bnd
use transport2fbe, only: transport2fbe_init

implicit none

!-------------------------------------------
! Find self-consistent initial configuration
!-------------------------------------------

integer :: j, jj, IM, ios, XSC0, XSC, jt1, jt2, jt3, jt_req, jkey
double precision :: Y, timeb
character(len=64) :: LISTSB(NSBMX)
character(len=132) :: STRI
double precision, external :: SWATCH
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
    call initMainWindow
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
if (TASK(1:3) /= 'BGD') then
    STRI(1:16) = ' ' ! Erase iteration number, iterations label top right
    call textvm(astra_gui%width-18*astra_gui_ref%dxlet, 2, STRI(1:16), 16)
    call textvm(astra_gui%width-17*astra_gui_ref%dxlet, astra_gui_ref%dylet + 1, STRI(1:14), 14)
endif

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
