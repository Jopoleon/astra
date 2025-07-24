module io_mod

use parameter_inc, only: NRDX, NARRX, NCNBM, NCNBTM, NCONST, NSBMX
use debugger, only: debug, flightsim
use const_inc, only: TSTART, TEND, TPAUSE

implicit none

integer, dimension(NSBMX) :: IFSBX, IFSBP
integer, dimension(NARRX) :: IFDFAX=-1, jbeg_arrx, NPTM
integer, dimension(NCONST) :: IFDFVX=-1
integer :: NSBR, NBNT=0, NCNBT=0, NGR
double precision :: resize
double precision, dimension(NARRX) :: TOUTX
double precision, dimension(NCNBM) :: CCOIL=0., VCOIL=0.
double precision, dimension((NCNBM+1)*NCNBTM) :: CCOILX=0., VCOILX=0.
double precision, dimension(NRDX, NARRX) :: XAXES, DATAX

character(len=4) :: machine, TASK
character(len=20), dimension(NSBMX) :: sbr_name
character(132) :: AWD, nml_file, equ_file, exp_file, NBFILE='***'

contains

    subroutine io_init

    logical :: nml_exists
    integer :: ios
    character(len=132) :: log_file
    double precision :: tbeg_nml, tend_nml, tpause_nml

    namelist / astra_log / equ_file, exp_file, task, machine, &
        debug, tbeg_nml, tend_nml, tpause_nml, resize, flightsim

    print*, 'exp file', exp_file
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

    return
    end subroutine io_init

end module io_mod
