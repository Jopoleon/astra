module io_mod

implicit none

integer, allocatable, dimension(:) :: IFDFAX, jbeg_arrx, NPTM, IFDFVX
integer :: n_sbr, restart, nr_x_max
double precision :: resize
double precision, allocatable :: TOUTX(:)
double precision, allocatable, dimension(:, :) :: XAXES, DATAX

character(len=4) :: machine, TASK
character(len=20), allocatable :: sbr_name(:)
character(132) :: AWD, astra_ext, nml_file, equ_file, exp_file, NBFILE='***'

contains

    subroutine io_init

    use json_vars, only: n_profx, n_var
    use debugger, only: debug, flightsim
    use scalars, only: TSTART, TEND, TPAUSE

    logical :: nml_exists
    integer :: ios
    character(len=132) :: log_file
    double precision :: tbeg_nml, tend_nml, tpause_nml

    allocate(TOUTX(n_profx))
    allocate(XAXES(nr_x_max, n_profx), DATAX(nr_x_max, n_profx))
    allocate(IFDFAX(n_profx), jbeg_arrx(n_profx), NPTM(n_profx))
    allocate(IFDFVX(n_var))
    IFDFAX = -1
    IFDFVX = -1

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
