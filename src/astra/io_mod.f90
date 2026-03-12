module io_mod

use parameter_inc, only: NRDX, NARRX, n_coils_max, n_sbr_max
use debugger, only: debug, flightsim
use const_inc, only: TSTART, TEND, TPAUSE

implicit none

integer, parameter :: NCONST=256

integer, dimension(n_sbr_max) :: IFSBX
integer, dimension(NARRX) :: IFDFAX=-1, jbeg_arrx, NPTM
integer, dimension(NCONST) :: IFDFVX=-1
integer :: n_sbr, NGR, restart
double precision :: resize
double precision, dimension(NARRX) :: TOUTX
double precision, dimension(n_coils_max) :: CCOIL=0., VCOIL=0.
double precision, dimension(NRDX, NARRX) :: XAXES, DATAX

character(len=4) :: machine, TASK
character(len=20), dimension(n_sbr_max) :: sbr_name
character(132) :: AWD, astra_ext, nml_file, equ_file, exp_file, NBFILE='***'

contains

    subroutine io_init

    logical :: nml_exists
    integer :: ios
    character(len=132) :: log_file
    double precision :: tbeg_nml, tend_nml, tpause_nml

    namelist / astra_log / equ_file, exp_file, task, machine, &
        debug, tbeg_nml, tend_nml, tpause_nml, resize, restart, flightsim

    tbeg_nml   = -1.
    tend_nml   = -1.
    tpause_nml = -1.
    resize = 1.

    log_file = 'tmp/astra.nml'
    OPEN(161, FILE=TRIM(log_file), delim='apostrophe')
    READ(161, nml=astra_log, iostat=ios)
    CLOSE(161)

! Override TPAUSE, TSTART & TEND: command line (astra.nml) has priority
    if (tbeg_nml   /= -1.) TSTART = tbeg_nml
    if (tend_nml   /= -1.) TEND   = tend_nml
    if (tpause_nml /= -1.) TPAUSE = tpause_nml

! Define namelist file nml_file
    nml_file = 'exp/nml/' // trim(exp_file)
    INQUIRE(FILE=trim(nml_file), EXIST=nml_exists)
    if (.not. nml_exists) then
        nml_file = 'exp/nml/' // trim(machine)
    endif

    CALL getenv('ASTRA_EXT', astra_ext)

    end subroutine io_init

end module io_mod
