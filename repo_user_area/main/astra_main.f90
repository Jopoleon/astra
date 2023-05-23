program astra

use parameter_inc, only: NSBMX, NRD
use outcmn_inc, only: VCOIL, CCOIL, CCOILO, DUMCT, DUMCTP, CTRLM, outcmn_init
use const_inc, only: IPART, const_init
use status_inc, only: status_init
use debugger, only: debug, flightsim

implicit none

!-------------------------------------------
! Find self-consistent initial configuration
!-------------------------------------------

character(len=64) :: LISTSB(NSBMX)

call fenvex()   !  Enable floating exception handling

!-------------------- Initial settings --------------------------------|

call outcmn_init
call const_init
call status_init

debug = 0 ! Initialise to: no debugging

call ininam(LISTSB)
IPART = 1   ! Mark initial iteration section

call read_input

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
