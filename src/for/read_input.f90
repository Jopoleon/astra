subroutine read_input
!----------------------------------------------------------------------|
!  NCONST   amount of simple variables readable (initiated)
!  NTVAR    maximal number of time slices for all variables
!  IVAR     number of actually defined variables
!  NTARR    maximal number of time slices for all arrays 
!  NGR     number of actually defined groups
!  NARRX    maximal number of arrays recognizable from a data file 
!----------------------------------------------------------------------|
! The subroutine is called once at the start-up, it reads the "exp" file
! and stores the time evolution of all input data in the arrays 
! VARDAT (simple variables) and 
! DATARR (sequential data in groups [grid,quantity]), TIMEX (group time)
! NGRIDX (grid size), NTYPEX(grid type), KTO (internal group name)
! FILTER (smoothing/transfer parameter)
! GDEX   (relative position of grid in DATARR)
! GDEY   (relative position of profile in DATARR)
! KOGDA  (relative position of time in TIMEX)
!        
! VARDAT(1,NTVAR) - time
! VARDAT(2,NTVAR) - value
! VARDAT(3,NTVAR) - error (not used)
! DATARR(NRDX*NTARR) - data array
!       Let   1 <= j <= NTARR is an ordinal number of array in DATARR
! TIMEX(j)  - time for this array
! NGRIDX(j) - number of grid points
! NTYPEX(j) - type of grid
! KTO(j)   - pointer to a position in the array EXARNM
!      so that EXARNM(KTO(j)) gives the name of the quantity
! GDEX(j)   - pointer to grid in DATARR
! GDEY(j)   - pointer to data in DATARR
! XAXES(NRDX,j) - not used here
!       Let   1 <= jx <= NARRX is an ordinal number of q-ty EXARNM(jx)
! KOGDA(jx)  - pointer to a position in the array TIMEX
!       KOGDA(KTO(j)) gives pointer to 1st time for KTO(j)
!----------------------------------------------------------------------|

use parameter_inc, only: NTVAR
use const_inc
use status_inc
use outcmn_inc, only: AWD, exp_file, nml_file, equ_file, rev_file, &
    TASK, machine, CPT, &
    TASKID, VERSION, AVERS, ARLEAS, AEDIT, COLTAB, IFDFVX, IFDFAX, KOGDA, KTO, &
    PRNAME, CFNAME, SRNAME, EXARNM, NBFILE, MSFILE, wall_gc_file, &
    NPRNAM, NCFNAM, NSRNAM, NEXNAM, FILTER, &
    NGR, NBNT, NCNBT, NBDMAX, NBDTMAX, NRDX, NTARR, NGRIDX, NTYPEX, NRW, &
    CCOILX, VCOILX, BNDR, BNDZ, BNDTIM, DATARR, TIMEX, GDEX, GDEY, GRAP, TIM7
    
use expdat
use char_manip, only: to_upper, str_in_list, clean_string
use debugger, only: markloc, debug, astra_stop, flightsim
use parse_utils
use timeoutput_inc, only: NTIMES, TTOUT

use numerical_tools, only: EXTRAP, INTEGR

use plasma_state, only: plasma_up

implicit none

integer, parameter :: MPEX=101, MSIGEX=1, MTEX=50, MSIG=1, MEXT=MPEX*MTEX

logical :: exilog, file_existence

integer :: jarr, INTYPE, jtype, SYSTEM, jbdry, ntim, ntim1
integer :: jj, j, j0, j1, IERR, ier_tab, jexar, jex1, jpos
integer :: KAB, KABC, KAWALL, KRTOR, KELONM, KTRICH
integer :: n_var, n_color, n_words, i_filter_glob
integer :: nt_u, nx_u, ios, ndim_u, jvar, jrt, jt, jthe

double precision :: resize
double precision :: tbeg_nml, tend_nml, tpause_nml
double precision, allocatable :: t_u(:), x_u(:), var_u(:), bnd_rz(:)
double precision :: XBDRY, YB, YB1, YXB, YXB1, ALFA, ALFA_GLOB, &
    VRDATA, FACTOR, TIMEVR, VRERR, ROC3A, YTP=-1.d9
character(len=6) :: VNAM, VNAMO, VNAMU, VNAMX, VTIM, VDAT, VERR, VARNAM, ARRNAM, keyword
character(len=30) :: rholbl
character(len=132) :: strarray(10), STRI, lin_upper, dir_path, fname, &
    err_msg, err_format, err_msg_exp, file_in, uname, uvar

namelist / astra_log / AWD, exp_file, equ_file, rev_file, TASK, machine, &
debug, tbeg_nml, tend_nml, tpause_nml, flightsim, resize

!----------------------------------------------------------------------|
! Fortran tests

call markloc('read_input')

call ADDTIME(CPT)  ! Initialize timer

call getarg(0, STRI)
j = min(len(TASKID), LEN_TRIM(STRI))
TASKID = STRI(1: j)

!----------------------------------------------------------------------|
! Initialisation with default values

tbeg_nml   = -1.
tend_nml   = -1.
tpause_nml = -1.
NITREQ = 1. ! Initialization: g95 does not like it in blockdata
i_filter_glob = 0 ! if i_filter_glob = 1, a global filter is set

plasma_up = 1  ! plasma is up by default, can be set to 0 for breakdown by the user in a user-defined sbr called with "<"

!----------------------------------------------------------------------|
! Parse file ".exe/version"
!----------------------------------------------------------------------|

open(131, FILE='exe/version', iostat=ios)

if (ios /= 0) then
    write(*,*)'>>> Warning: Unknown version'
else
    do j=1,5
        read(131,'(A)') STRI
    enddo
    j = index(STRI, 'Version')
    VERSION = STRI(j: j+30)//char(0)
    close(131)
    j0 = index(VERSION, '.')
    if (j0 == 0) then
        write(*,*)'>>> Warning: Unknown version'
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

!----------------------------------------------------------------------|
! Read tables to set global arrays PRNAME, CFNAME, SRNAME
!----------------------------------------------------------------------|

call set_vars('main/variables.txt', PRNAME , NPRNAM)
call set_vars('main/constants.txt', CFNAME , NCFNAM)
call set_vars('main/internal.txt' , SRNAME , NSRNAM)

do j=1, NPRNAM
    SELECT CASE (PRNAME(j))
    CASE('AB    ')
        KAB    = j
    CASE('ABC   ')
        KABC   = j
    CASE('AWALL ')
        KAWALL = j
    CASE('RTOR  ')
        KRTOR  = j
    CASE('ELONM ')
        KELONM = j
    CASE('TRICH ')
        KTRICH = j
    END SELECT
enddo

TIME = TSTART  ! Here TSTART=0

! Read file status.inc

call set_vars('main/profiles_x.txt', EXARNM, NEXNAM)

!----------------------------------------------------------------------|
! Read file 'tmp/<exp><equ>.nml'
!----------------------------------------------------------------------|

call GETENV('expfile', exp_file)
call GETENV('equfile', equ_file)

file_in = 'tmp/' // TRIM(exp_file) // TRIM(equ_file) // '.nml'

OPEN(161, FILE=TRIM(file_in), delim='apostrophe')
READ(161, nml=astra_log, iostat=ios)
CLOSE(161)

!define namelist file nml_file
nml_file = 'exp/nml/' // trim(exp_file)
INQUIRE(FILE=trim(nml_file), EXIST=file_existence)
if (.not. file_existence) then
   nml_file = 'exp/nml/' // trim(machine)
endif

call path_split(rev_file, dir_path, fname, jpos)
if (LEN_TRIM(fname) == 0) fname = 'profil.dat'
! Disallowing user-defined subdirs for Review file
rev_file = '.res/' // TRIM(fname)

!----------------------------------------------------------------------|
! Read file equ/log/<model>, checking existence of obsolete equ/<model>.log 
!----------------------------------------------------------------------|

jj = LEN_TRIM(equ_file)
if (jj == 0) call astra_stop('>>> read_input: Error, empty model file name')

call path_split(equ_file, dir_path, fname, jpos)

if (jpos == 0) then ! no subdir
    file_in = 'equ/log/'//TRIM(fname)
else ! equ/<subdir>/log/<model>
    file_in = 'equ/'//TRIM(dir_path)//'log/'//TRIM(fname)
endif

inquire(file=TRIM(file_in), exist=EXILOG)

if (.not. EXILOG)  then ! Missing log file

    inquire(file=TRIM(equ_file)//'.log', exist=EXILOG)
    if ( EXILOG )  then
        write(*, *) '>>> Warning: File "' // TRIM(equ_file) // '.log," has been found'
        write(*, *) '             Most probably it should be moved to "log/' // TRIM(file_in) // '"'
    else
        write(*, *) '>>> Warning: file "' // TRIM(file_in) // '" missing'
        write(*, *) 'Press any key to continue at your own risk'
        pause
    endif

else  ! Read log file

    n_var = 37
    call assign_val(file_in, n_var , PRNAME(1: n_var) , DEVAR (1: n_var) , n_color)
    call assign_val(file_in, NCFNAM, CFNAME(1: NCFNAM), CONSTF(1: NCFNAM), n_color)
    call assign_val(file_in, NSRNAM, SRNAME(1: NSRNAM), DELOUT(1: NSRNAM), n_color)

    NA1   = DELOUT(13)
    NUF   = DELOUT(14)
    NBND  = DELOUT(19)
    XFLAG = DELOUT(20)
! Replace equivalence

    read(171, *, iostat=ios) (COLTAB(j), j=1, n_color)
    close(171)

endif

!----------------------------------------------------------------------|
! Skipping EX-file reading for now
!----------------------------------------------------------------------|

!----------------------------------------------------------------------|
! Read experimental file
!----------------------------------------------------------------------|

file_in='exp/' // TRIM(exp_file)

err_msg_exp = '>>> Data file "' // TRIM(exp_file) // '" error:\n    '

open(201, FILE=TRIM(file_in), iostat=ios)
if (ios /= 0) call astra_stop('>>> read_input: No such experimental variant "' // TRIM(exp_file) // '"')

read(201, '(A132)', iostat=ios) XLINE1
if (ios /= 0) call astra_stop(err_msg_exp // 'in header')

read(201, '(A132)') XLINE2

!----------------------------------------------------------------------|
! Read simple variable loop (between the labels "5" and "10"):

VNAMO = ' '
VNAMU = ' '

parse_exp_1d: do

    read(201, '(A132)', end=39) STRI
    err_format  =  err_msg_exp // '".  Format error in group ' // TRIM(STRI)
    lin_upper = to_upper(STRI)
    if (lin_upper(1: 1) == '!') CYCLE parse_exp_1d
    if (lin_upper(1: 3) == 'END') goto 39

    call split2array2(lin_upper, strarray, n_words)

    do j=1, n_words
        keyword = strarray(j)
        jpos = str_in_list(keyword, (/'POINTS', 'GRIDTY', 'FILTER', 'PROFIL'/))
        if (jpos > 0) EXIT parse_exp_1d
    enddo

    if (strarray(1) == 'NAMEXP') then
        VNAM = strarray(2)
        VTIM = strarray(4)
        VERR = '0.'
        jvar = str_in_list(VNAM(1:6), PRNAME) 
        if (jvar == 0) CYCLE parse_exp_1d  ! var doesnt exist
        IVAR = IVAR + 1
        ntim = 0
        read(vtim, *) ntim !from NTIMES
        factor = 1.
        if (IVAR + ntim > NTVAR) then
            write(err_msg, '(A, i)') '>>> read_input: Time dependent variable strings >', NTVAR
            call astra_stop(err_msg)
        endif
        if (ntim <= 1) then
            ntim = 1
            IFDFVX(jvar) = 0
        else
            IFDFVX(jvar) = 1
        endif
! Read "time" array & function array
        read(201, *, iostat=ios) (VARDAT(1, IVAR+jj), jj=1, ntim)
        if (ios /= 0) call astra_stop(err_format)
        read(201, *, iostat=ios) (VARDAT(2, IVAR+jj), jj=1, ntim)
        if (ios /= 0) call astra_stop(err_format)
        DEVAR(jvar) = factor*VARDAT(2, IVAR+1)
        do jj=1, ntim
            IVAR = IVAR + 1
            INDVAR(IVAR) = jvar
            VARDAT(2, IVAR) = factor*VARDAT(2, IVAR)
            VARDAT(3, IVAR) = 0.
        enddo
        VNAMO = VNAM
        CYCLE parse_exp_1d
    else
! Variable name, time, value, error
        VNAM = VARNAM(lin_upper(1: 6), ier_tab)
        VTIM = STRI(9 : 14)
        VDAT = STRI(17: 22)
        VERR = STRI(25: 30)
    endif

! A special treatment required for time independent NA1, TSTART, TEND
    if (VNAM == 'NA1   ' .or. VNAM == 'TSTART' .or. VNAM == 'TEND  ') then

        err_msg = err_msg_exp // '"' // TRIM(VNAM) // '"'
        if (ier_tab /= 0 ) then
            call astra_stop(TRIM(err_msg) // ': tabulation not allowed in this type of input')
        endif
        if (VNAMO == VNAM) then
            call astra_stop(TRIM(err_msg) // ' cannot vary in time')
        endif

        call READF6(VDAT, VRDATA, IERR)
        VNAMO = VNAM
        if (VNAM == 'NA1   ') NA1 = VRDATA
        if (VNAM == 'TSTART') then
            TSTART = VRDATA
            TIME  = TSTART
        endif
        if (VNAM == 'TEND  ') TEND = VRDATA
        CYCLE parse_exp_1d
    endif

    ntim = 0
    factor = 1.

    jvar = str_in_list(VNAM, PRNAME)
    if (jvar == 0) CYCLE parse_exp_1d

! U-file name duplicated:
    if (VNAM == VNAMU) call astra_stop(err_msg_exp // 'ambiguous ' // TRIM(VNAM) // ' definition')

    j1 = index(STRI,':')

    if (j1 == 0) then  ! time-dependent values in exp-file

! Tabulation encountered in the old-standard line:
        if (ier_tab /= 0 ) then
            err_msg = err_msg_exp // '"' // TRIM(VNAM) // &
                '": tabulation not allowed in this type of input'
            call astra_stop(err_msg)
        endif

        err_msg = '>>> read_input: File "' // TRIM(file_in) // '" reading error'
        call READF6(VTIM, TIMEVR, IERR)
        if (IERR /= 0) call astra_stop(err_msg)
        call READF6(VDAT, VRDATA, IERR)
        if (IERR /= 0) call astra_stop(err_msg)
        call READF6(VERR, VRERR, IERR)
        if (IERR /= 0) call astra_stop(err_msg)

        IFDFVX(jvar) = 0
        DEVAR(jvar) = factor*VRDATA
        if (VNAM == VNAMO) IFDFVX(jvar) = 1
! Repeated name
        IVAR = IVAR+1
        if (IVAR > NTVAR) then
            write(err_msg, '(A, i)') '>>> read_input: Time dependent variable strings >', NTVAR
            call astra_stop(err_msg)
        endif

        INDVAR(IVAR) = jvar
        VARDAT(1, IVAR) = TIMEVR
        VARDAT(2, IVAR) = factor*VRDATA
        VARDAT(3, IVAR) = VRERR

    else  ! ":" found in the input string "STRI", pointer to U-file

        if (VNAM == VNAMO) call astra_stop(err_msg_exp // 'ambiguous ' // TRIM(VNAM) // ' definition')

        call parse_u_line(STRI, uvar, uname, factor)
        VNAMU = VNAM
        call ufheader(TRIM(uname), ndim_u, nt_u, nx_u, rholbl)
        allocate(t_u(nt_u))
        allocate(x_u(nx_u))
        allocate(var_u(nt_u*nx_u))
        call ufrd(TRIM(uname), ndim_u, nt_u, nx_u, t_u, x_u, var_u)

        DEVAR(jvar) = factor*var_u(1)
        if (nt_u == 1) then
            IFDFVX(jvar) = 0
        else
            IFDFVX(jvar) = 1
        endif
        do jj=1, nt_u
            IVAR = IVAR + 1
            INDVAR(IVAR) = jvar
            VARDAT(1, IVAR) = t_u(jj) ! git
            VARDAT(2, IVAR) = factor*var_u(jj)
            VARDAT(3, IVAR) = 0.
        enddo
        deallocate(t_u, x_u, var_u)
    endif

    VNAMO = VNAM

enddo parse_exp_1d

close(201)

!----------------------------------------------------------------------|
! Start Astra standard assignments
! GIT allocate all profiles here!!

if (NA1 > NRD) then
    write(err_msg, '(2A, i)') '>>> FATAL ERROR: The radial grid size out of range.\n', &
        '                 Parameter "NA1" cannot exceed', NRD
    call astra_stop(err_msg)
endif

! Override TPAUSE, TSTART & TEND: command line (astra.nml) has priority
if (tbeg_nml   /= -1.) TSTART = tbeg_nml
if (tend_nml   /= -1.) TEND   = tend_nml
if (tpause_nml /= -1.) YTP    = tpause_nml

!----------------------------------------------------------------------|
! Uncomment the next line to set the "WAIT" mode at the start
! TASK='WAIT'

TIME = TSTART
if (YTP > -1.d8) TPAUSE = YTP

call INTVAR

if (.not. EXILOG) then
! Determine ABC & AB if not defined by "exp" file
    if (IFDFVX(KABC) < 0 .and. IFDFVX(KAB) < 0) then
        write(*, *) '>>> The minor radius AB is not defined in "exp/' // TRIM(exp_file) // '"'
        write(*, *) '  This can cause equilibrium convergence problem'
    endif
    if (IFDFVX(KAB) > 0) write(*, *) '>>> Warning: AB cannot depend on time'

    if (IFDFVX(KABC) >= 0 .and. IFDFVX(KAB) <  0) AB = ABC
    if (IFDFVX(KABC) <  0 .and. IFDFVX(KAB) >= 0) then
        ABC = AB
        IFDFVX(KABC) = 0
    endif
! Determine AWALL if not defined by an "exp" file
    if (IFDFVX(KAWALL) < 0) AWALL = AB
endif

if (AWALL < AB) then
     if (IFDFVX(KAWALL) >= 0) write(*, *) '>>> Warning: AWALL < AB.  Setting AWALL = AB'
     AWALL = AB
endif
if (IFDFVX(KAWALL) < 0 .and. (AWALL < AB .or. AWALL > 1.2*AB)) AWALL = AB

if (ABC > AB) then
    write(*, *) '>>> Warning: ABC cannot exceed AB. Setting ABC = AB'
    if (IFDFVX(KABC) < 1) ABC = AB
endif

if (AWALL > 1.2*AB .or. AWALL > RTOR) then
    write(*, *) '>>> Warning: AWALL is set unreasonably large'
    write(*, *) '    Check settings in data and log files'
endif
IFDFVX(KAB)    = 4
IFDFVX(KELONM) = 4
IFDFVX(KRTOR)  = 4
IFDFVX(KTRICH) = 4
IFDFVX(KAWALL) = 4


!----------------------------------------------------------------------|
! Skipping for now: JAMS through Ex-files
!----------------------------------------------------------------------|
! Start reading radial profiles:  Radial grid: jbdry 
!    Option for different input grids can be added

VNAMO = ' '

open(201, FILE=TRIM(file_in), iostat=ios)
read(201, '(A132)', ERR=906, END=39) STRI
read(201, '(A132)', ERR=906, END=39) STRI

NGR = 0
jarr = 0
NBNT = 0
ALFA_GLOB = 0.001

parse_exp_2d: do

    INTYPE = -1
    TIMEVR = .0
    ntim = 0
    jbdry = 0
    factor = 1.
    ALFA = ALFA_GLOB

    read(201, '(A132)', iostat=ios) STRI
    if (ios < 0) EXIT parse_exp_2d
    if (ios > 0) call astra_stop(err_format)
    if (STRI(1:6) == 'FILTER') i_filter_glob = 1
    lin_upper = to_upper(STRI)
    if (LEN_TRIM(lin_upper) == 0) CYCLE parse_exp_2d
    if (lin_upper(1: 1) == '!') CYCLE parse_exp_2d

    VNAM = VARNAM(lin_upper(1: 6), ier_tab)
    VNAMX = ARRNAM(VNAM)
    if (ier_tab /= 0) CYCLE parse_exp_2d    ! Ignore lines starting with a blank
    jex1 = str_in_list(VNAMX, EXARNM) ! Checks if VNAM is in array list

    if (jex1 == 0) then ! 1d, or no u-file
        jpos = str_in_list(VNAM, (/'POINTS', 'NAMEXP', 'GRIDTY', 'NTIMES', 'FILTER', 'FACTOR'/) )
        if (jpos == 0) then ! A 1D variable
            VNAM = ' '
            CYCLE parse_exp_2d
        else

            VNAM = ' '
            call split2array2(lin_upper, strarray, n_words)

            do j=1, n_words, 2
                keyword = strarray(j)(1: 6)
                SELECT CASE(keyword)
                CASE('POINTS')
                    read(strarray(j+1), *, iostat=ios) jbdry
                    if (ios /= 0) call astra_stop(err_format)
                CASE('NAMEXP')      ! New variable, exp-block
                    VNAM = VARNAM(strarray(j+1), ier_tab)
                CASE('GRIDTY')
                    read(strarray(j+1), *, iostat=ios) INTYPE
                    if (ios /= 0) call astra_stop(err_format)
                CASE('NTIMES')
                    read(strarray(j+1), *, iostat=ios) ntim
                    if (ios /= 0) call astra_stop(err_format)
                CASE('FILTER')
                    read(strarray(j+1), *, iostat=ios) ALFA
                    if (ios /= 0) call astra_stop(err_format)
                    if (i_filter_glob == 1) then
	         	alfa_glob = alfa
			i_filter_glob = 0
		    endif
                CASE('FACTOR')
                    read(strarray(j+1), *, iostat=ios) factor
                    if (ios /= 0) call astra_stop(err_format)
                CASE('PROFIL')
                    EXIT
                CASE DEFAULT
                    write(*, *) '>>> read_input error unknown key word in string:'
                    write(*, *) TRIM(STRI)
                    EXIT
                END SELECT
            enddo

! If VNAM not there, ignore this line, such as line 'FILTER 0.001'
! Needed for backward compatibility
            if (LEN_TRIM(VNAM) == 0) then
                CYCLE parse_exp_2d
            else
                VNAMO = VNAM
            endif
        endif
    endif

! VNAM is updated via NAMEXP keyword, in case of old non-u-file array
    VNAMX = ARRNAM(VNAM) ! Appends 'X'

    jexar = str_in_list(VNAMX, EXARNM) ! Checks if VNAMX-string is in array EXARNM

    if (VNAM /= VNAMO .and. IFDFAX(jexar) < 0) then
        KOGDA(jexar) = NGR + 1
    endif

! Special arrays

    SELECT CASE(VNAMX)

    CASE('CCOILX')
        NCNBT = 0
! read only if NCNBT==0, i.e. CCOILX was not defined before
        call read_arrx(201, NCNBT, ntim, NCNB, STRI, CCOILX)
        VNAMO = VNAM

    CASE('VCOILX') !note that both CCOIL and VCOIL need to appear in the exp file with the same number of points and times
        NCNBT = 0
        call read_arrx(201, NCNBT, ntim, NCNB, STRI, VCOILX)
        VNAMO = VNAM

    CASE ('BNDX  ')
        if (NBNT /= 0) then
            call astra_stop(err_msg_exp // 'Boundary must be defined in a single group')
        endif
        j = INDEX(lin_upper, 'POINTS')
        if (j /= 0) read(STRI(j+6:),*) NBND
        if (j == 0) then
            call astra_stop(err_msg_exp // 'Number of boundary points must be defined')
        endif

        NBNT = max(ntim, 1)
        write(*, *) 'Reading BND, dims:', nbnd, ntim

        if (NBND > NBDMAX) then
            write(err_msg, '(2A, i)') err_msg_exp, 'Boundary data #theta must not exceed ', NBDMAX
            call astra_stop(err_msg)
        endif
         if (NBNT > NBDTMAX) then
            write(err_msg, '(2A, i)') err_msg_exp, 'Boundary data #times must not exceed ', NBDTMAX
            call astra_stop(err_msg)
        endif

! Input order:
! t_1  t_2  t_3
! r_1(t_1) r_1(t_2) r_1(t_3) 
! z_1(t_1) z_1(t_2) z_1(t_3) 
! r_2(t_1) r_2(t_2) r_2(t_3) 
! z_2(t_1) z_2(t_2) z_2(t_3)

        read(201, *, iostat=ios) (BNDTIM(j), j=1, NBNT)
        allocate(bnd_rz(2*NBNT*NBND))
        read(201, fmt=*, iostat=ios) (bnd_rz(j), j=1, 2*NBNT*NBND)
        jrt = 1
        do jthe=1, NBND
            do jt=1, NBNT
                BNDR((jthe-1)*NBNT + jt) = bnd_rz(jrt)
                jrt = jrt + 1
            enddo
            do jt=1, NBNT
                BNDZ((jthe-1)*NBNT + jt) = bnd_rz(jrt)
                jrt = jrt + 1
            enddo
        enddo
        deallocate(bnd_rz)

        if (ios /= 0) then
            call astra_stop(err_msg_exp // 'More data items than data values for BND group')
        endif
        VNAMO = VNAM
        CYCLE parse_exp_2d

    CASE('BNDUX ')
        write(*, *) 'Reading BND from u file'
        read(201, '(A)', iostat=ios) STRI ! u-file name in exp-file
        if (ios < 0) EXIT parse_exp_2d
        call ufheader('udb/'//trim(STRI)//'_r', ndim_u, nt_u, nx_u, rholbl)
        NBND = nx_u
        NBNT = nt_u
        if (NBND > NBDMAX) then
            write(err_msg, '(2A, i)') TRIM(err_msg_exp), 'Boundary data #theta must not exceed ', NBDMAX
            call astra_stop(err_msg)
        endif
        if (NBNT > NBDTMAX) then
            write(err_msg, '(2A, i)') TRIM(err_msg_exp), 'Boundary data #times must not exceed ', NBNT
            call astra_stop(err_msg)
        endif

        allocate(x_u(nx_u))
        call ufrd('udb/' // trim(STRI) // '_r', ndim_u, nt_u, nx_u, BNDTIM(1:nt_u), x_u, BNDR(1:nx_u))
        call ufrd('udb/' // trim(STRI) // '_z', ndim_u, nt_u, nx_u, BNDTIM(1:nt_u), x_u, BNDZ(1:nx_u))
        deallocate(x_u)

        VNAMO = VNAM

        CYCLE parse_exp_2d

    CASE('      ') 
        CYCLE parse_exp_2d

    CASE('ENDX  ') 
        EXIT parse_exp_2d

    END SELECT

    if (jexar == 0) then
        write(*, *) 'Array name ' // TRIM(VNAMX) // ' not found, skipping line'
        CYCLE parse_exp_2d
    endif

! Check for U-file command string
    call parse_u_line(STRI, uvar, file_in, factor)

    if (LEN_TRIM(file_in) > 0) then ! string 'U-FILE' found in this line

        call ufheader(TRIM(file_in), ndim_u, nt_u, nx_u, rholbl)
        allocate(t_u(nt_u))
        allocate(x_u(nx_u))
        allocate(var_u(nt_u*nx_u))
        call ufrd(TRIM(file_in), ndim_u, nt_u, nx_u, t_u, x_u, var_u)

        if (nx_u > NRDX) then
            write(err_msg, '(3A, i, A, i)') '>>> U-file "', TRIM(file_in), & 
                '" error: radial grid size ', nx_u, ' is larger than', NRDX
            call astra_stop(err_msg)
        endif

        if (NGR + nt_u + 1 > NTARR) then
            write(err_msg, '(A, i)') &
                '>>> read_input: Number of time dependent arrays cannot exceed', NTARR
            call astra_stop(err_msg)
        endif

        if (jarr + nx_u  > NRDX*NTARR) call astra_stop('>>> read_input: Buffer size exceeded')
        call CHECKU(INTYPE, ABC, AB, XBDRY, DATARR(jarr+1), nx_u, jbdry, rholbl, file_in)

        do j=1, jbdry
            DATARR(jarr + j) = x_u(j)
        enddo
        do j=1, nt_u
            TIMEX(NGR + j) = t_u(j)
        enddo

        if (INTYPE == 18 .or. INTYPE == 19) then
            DATARR(jarr + 2: jarr + jbdry) = DATARR(jarr + 1)
            if (INTYPE == 18) DATARR(jarr + 1) = RTOR
            if (INTYPE == 19) DATARR(jarr + 1) = 0.
        endif

        jarr = jarr + jbdry

        if (jarr + nx_u + (nt_u - 1)*jbdry > NRDX*NTARR) then
            call astra_stop('>>> read_input: Buffer size exceeded')
        endif

        jrt = jarr
        do jj=1, nt_u
            do j=1, jbdry
                jrt = jrt + 1
                DATARR(jrt) = var_u(jj + (j - 1)*nt_u)
            enddo
        enddo

        deallocate(t_u, x_u, var_u)

        if (INTYPE > 13 .and. INTYPE /= 19) write(*, *) 'Unknown U-file type'
        YXB = DATARR(jarr)
        YXB1 = DATARR(jarr-1)
        do j=1, nt_u
            NGR = NGR + 1
            if (j == 1) then
                GDEX(NGR) = jarr - jbdry + 1
                GDEY(NGR) = jarr + 1
            else
                GDEX(NGR) = GDEX(NGR - 1)
                GDEY(NGR) = jarr + 1
            endif
            do j0=1, jbdry
                DATARR(jarr + j0) = DATARR(jarr + j0)*factor
            enddo
            jarr = jarr + jbdry
            if (jbdry /= nx_u) then
                YB = DATARR(jarr)
                YB1 = DATARR(jarr-1)
                YB = (YB*(XBDRY - YXB1) - YB1*(XBDRY - YXB))/(YXB - YXB1)
                DATARR(jarr) = YB
            endif
            KTO(NGR) = jexar
            NGRIDX(NGR) = jbdry
            NTYPEX(NGR) = INTYPE
            FILTER(NGR) = ALFA
        enddo

    else ! read exp-block data
        if (jbdry <= 1) then
            write(err_msg, '(3A, 8X, A)') err_msg, 'Input quantity: ', TRIM(VNAM), &
                'Number of grid points must be > 1'
            call astra_stop(err_msg)
        endif
        if (jbdry > NRDX) then
            write(err_msg, '(2A, i)') err_msg_exp, 'Number of radial points > ', NRDX
            call astra_stop(err_msg)
        endif

        jtype = 1 + INTYPE/10
        ntim1 = max(ntim, 1)
        if (NGR + ntim1 > NTARR) then  
            write(err_msg, '(2A, i)') err_msg_exp, &
                'Number of time dependent arrays cannot exceed', NTARR
            call astra_stop(err_msg)
        endif
        if (jarr + jtype*jbdry + 1 > NRDX*NTARR) then
            call astra_stop('>>> read_input: Buffer size exceeded')
        endif

! INTYPE unknown
        if (INTYPE < 0 .or. INTYPE > 20) then 
            write(err_msg, '(A, i0, A)') '>>> ERROR: Unknown input type =', INTYPE, ',  ignored'
            call astra_stop(err_msg)
        endif

        if (ntim > 0) read(201, *, ERR=906) (TIMEX(NGR+j), j=1, ntim)

        do j=1, ntim1
            NGR = NGR + 1
            KTO(NGR) = jexar
            NGRIDX(NGR) = jbdry
            NTYPEX(NGR) = INTYPE
            FILTER(NGR) = ALFA
            if (j == 1) then
                KOGDA(jexar) = NGR
                GDEX(NGR) = jarr + 1
                if (INTYPE == 18 .or. INTYPE == 19) then
                    jarr = jarr + 1
                    read(201, *, ERR=906) DATARR(jarr)
                endif
                do j1=1, jtype
                    read(201, *, iostat=ios) (DATARR(jarr + jj), jj=1, jbdry)
                    if (ios /= 0) then
                        write(err_msg, '(3A, /, A, 1p, 6e12.4)') err_msg_exp, &
                            '".  Format error in group ', TRIM(STRI), 'Last data read: ', &
                            (DATARR(jarr+jj), jj=1, jbdry)
                            call astra_stop(err_msg)
                    endif
                    jarr = jarr + jbdry
                    if (jtype == 2 .and. INTYPE <= 17 .and. j1 == 1) XBDRY = DATARR(jarr)
                enddo
                GDEY(NGR) = jarr - jbdry + 1
                do jj=1, jbdry
                    DATARR(jarr - jbdry + jj) = factor*DATARR(jarr - jbdry + jj)
                enddo
            else
                GDEX(NGR) = GDEX(NGR-1)
                GDEY(NGR) = jarr + 1
                read(201, *, iostat=ios) (DATARR(jarr + jj), jj=1, jbdry)
                if (ios /= 0) then
                    write(err_msg, '(3A, /, A, 1p, 6e12.4)') err_msg_exp, &
                            '".  Format error in group ', TRIM(STRI), 'Last data read: ', &
                        (DATARR(jarr+jj), jj=1, jbdry)
                    call astra_stop(err_msg)
                endif
                do jj=1, jbdry
                    DATARR(jarr + jj) = factor*DATARR(jarr + jj)
                enddo
                jarr = jarr + jbdry
            endif
        enddo
        if (INTYPE == 17) write(*, *) 'X-axis analysis is not implemented GRIDTYPE =', INTYPE
    endif

    VNAMO = VNAM
    IFDFAX(jexar) = 0

enddo parse_exp_2d

39 continue

close(201)

!-----------------------
! End reading "exp" file
!-----------------------

!if boundary is given, calculates initial geometry from that
if (NBNT > 0) then
!find time index of most proximum boundary
    j=1
    do jt=1,NBNT
        if (BNDTIM(jt) <= TSTART) j = jt
    enddo
    jt=j
		
    allocate(bnd_rz(2*NBND))
    do jthe=1, NBND
        bnd_rz(jthe)      = BNDR((jthe-1)*NBNT + jt)
        bnd_rz(NBND+jthe) = BNDZ((jthe-1)*NBNT + jt) 
    enddo
!calculate ABC
    ABC = (maxval(bnd_rz(1: nbnd)) - minval(bnd_rz(1: nbnd)))/2.
!calculate elong
    YB  = (maxval(bnd_rz(1: nbnd))        + minval(bnd_rz(1: nbnd)       ))/2. !Rgeo
    YB1 = (maxval(bnd_rz(NBND+1: 2*nbnd)) + minval(bnd_rz(nbnd+1: 2*nbnd)))/2. !Zgeo
    ELONG = (maxval(bnd_rz(NBND+1: 2*nbnd)) - minval(bnd_rz(nbnd+1: 2*nbnd)))/(2.*ABC)
    ELONG = max(ELONG, 1.d0)
    deallocate(bnd_rz)
endif

!assign variables here for initialization:

VOLUME = GP2*GP*RTOR*AB*AB*ELONG

ROC  = ROC3A(RTOR, SHIFT, ABC, ELONG, TRIAN)
ROCO = ROC

! Initialization of magnetic quantities
BTN = BTOR
FTO = GP*BTOR*ROC**2
FTN = GP*BTN*ROC**2
IPLN = IPL

if (AWALL > RTOR+SHIFT) then
    ROWALL = AWALL*SQRT(max(ELONM, ELONG))
elseif (ELONM > ELONG) then
    ROWALL = ROC3A(RTOR, SHIFT, AWALL, ELONM, TRICH)
else
    ROWALL = ROC3A(RTOR, 0.d0, AWALL, ELONG, TRIAN)
endif

!Efable
! basic space steps, these do not change
! Two versions: HROX = 1/NA1, XRHO(NA1)<1, HROX=1/(NA1-1/2) gives XRHO(NA1) = 1
! Case with XRHO(NA1)=1-HROX/2, XSHO(NA1) = 1   (fluxes, MU, etc are on LCFS at at last grid point NA1 / NE, TE, and quantities are a bit inside at last grid point NA1)

if (int(FLXDR) == 1) then
    HROX  = 1.0/(NA1)
else if (int(FLXDR) == 0) then
    HROX  = 1.0/(NA1-0.5)
endif

do j=1, NRD
    XRHO(j) = (j - 0.5)*HROX
    SXHO(j) = j*HROX
! real space grids, these are function of ROC 
    RHO (j) = XRHO(j)*ROC
    SRHO(j) = SXHO(j)*ROC
enddo

! real space grids, these are function of ROC 
HRO  = HROX*ROC

! Compute NB1
NB1 = NA1
NA  = NA1 - 1

call SETGEO(0)
call NEW_GRID

do J=1, NB1
    G22(J) = RHO(J)
    VR(J)  = (GP2*(RTOR + SHIFT))**2*RHO(J)/RTOR
    VRS(J) = (GP2*(RTOR + SHIFT))**2*J*HRO/RTOR
    G11(J) = VRS(J)
enddo

call INTEGR(RHO, 1, VR, VOLUM, NA1)

n_bouncon = NA1

PSIBO = FP(NA1)
call EXTRAP(XRHO(1: NA1), FP(1: NA1), 0.0, 1, PSIAX, 2, NA1)

VOLUME = VOLUM(NA1)

NAB = NA1
if (NA1 < NB1 .and. AB > ABC) then
    do j=NA1+1, NB1
        if (AMETR(j) < AB) NAB = j
    enddo
    if (NAB < NB1) NAB = NAB + 1
endif
ROB = RHO(NAB)

AMETR(NAB) = AB
AMETR(NA1) = ABC

NEO   = NE
TEO   = TE
FPO   = FP
UPARO = UPAR
VRO   = VR    
MRHO  = AMAIN*NE
UPS0  = MRHO*RTOR 
UPS0O = UPS0

do j=1, NEXNAM
    if (ARXUSE(j) /= 0) then
        if (IFDFAX(ARXUSE(j)) < 0) then
            write(*, *) '>>> Warning >>> X-array used but not defined: "', &
                TRIM(EXARNM(ARXUSE(j))), '"', ARXUSE(j), IFDFAX(ARXUSE(j))
        endif
    endif
enddo

if (.not. IFDEFX('TEX   ') ) then
    j = system('grep HEXP= ./tmp/*.tmp | grep TEX > /dev/null')
    if (j == 0) write(*,*) '>>> Warning >>> X-array "TEX" is used but not defined'
endif
if (.not. IFDEFX('TIX   ') ) then
    j = system('grep XEXP= ./tmp/*.tmp | grep TIX > /dev/null')
    if (j == 0) write(*,*) '>>> Warning >>> X-array "TIX" is used but XEXP not defined'
    j = system('grep SVCXX ./tmp/*.tmp | grep TIX > /dev/null')
    if (j == 0) write(*,*) '>>> Warning >>> X-array "TIX" is used but not defined'
endif

DELOUT(13) = NA1
DELOUT(14) = NUF
DELOUT(19) = NBND
DELOUT(20) = XFLAG
TIMEQL = TIME - DTEQL - 1.d-7

! Define TAU, TAUMIN, TAUMAX, TSCALE, DROUT, DTOUT, DPOUT
TTOUT(1) = -1.d10
call set_timescale(NTIMES)  ! Set time scale (mode 6)

GRAP(1:NRW) = AB

TIM7(1) = TINIT
if (TIME > TINIT + 1.025*abs(TSCALE)) TINIT = TSTART
TIM7(3) = abs(TSCALE)/8.

call inquire_fname('cnf', TRIM(exp_file), TRIM(machine), wall_gc_file)
call inquire_fname('nbi', TRIM(exp_file), TRIM(machine), NBFILE)
call inquire_fname('mse', TRIM(exp_file), TRIM(machine), MSFILE)

return

906 continue
call astra_stop(err_format)

end subroutine read_input

!---------------------------------------------------------------------
subroutine set_timescale(n_times)
! Define time scales (former SETTSC)

use outcmn_inc, only: equ_file
use const_inc, only: TAUPRP, TAUMIN, TAUMAX, TAU, TSCALE, &
        VOLUME, DTOUT, DROUT, DPOUT

implicit none

integer, intent(in) :: n_times

logical :: EXILOG
integer :: j
double precision, dimension(8) :: YS

inquire(file='equ/log/' // TRIM(equ_file), exist=EXILOG)

if ( EXILOG ) then ! Always the case (equ/log exists)
    TAUPRP = TAUMIN
else
    TAUMAX = 0.01*VOLUME
    YS(1) = 0.00000010
    YS(2) = 0.00000015
    YS(3) = 0.00000020
    YS(4) = 0.00000025
    YS(5) = 0.00000030
    YS(6) = 0.00000040
    YS(7) = 0.00000050
    YS(8) = 0.00000075
    do while (1.1*TAUMAX > YS(8))
        YS(1: 8) = 10.*YS(1: 8)
    enddo
    do J=1, 8
        if (1.1*TAUMAX <= YS(J)) EXIT
    enddo

    TSCALE = 0.1*n_times*YS(J)

! 115=(right_label_position)/IDT=575/5
    do while (TSCALE*115./n_times > YS(8))
        YS(1: 8) = 10.*YS(1: 8)
    enddo
    do J = 1, 8
        if (TSCALE*115./n_times <= YS(J)) EXIT
    enddo

    TSCALE = YS(J)
    TAUMAX = 0.01*VOLUME
    YS(1) = 0.00000010
    YS(2) = 0.00000015
    YS(3) = 0.00000020
    YS(4) = 0.00000025
    YS(5) = 0.00000050
    YS(6) = 0.00000075

    do while (1.1*TAUMAX > YS(6))
        YS(1: 6) = 10.*YS(1: 6)
    enddo
    do J=1, 6
        if (1.1*TAUMAX <= YS(J)) EXIT
    enddo

    TAUMAX = YS(J)
    TAUMIN = 0.001*TAUMAX
    TAU    = TAUMIN
    TAUPRP = TAUMIN
    DTOUT  = 0.01*TAUMAX
    DROUT  = 0.02*TAUMAX
    DPOUT  = 0.2*TAUMAX
endif

end subroutine set_timescale

!---------------------------------------------------------------------
subroutine READF6(LINE, F6, IERR)
! Read a number in LINE to 6-positinal field 

implicit none

character(len=6), intent(in) :: LINE
integer, intent(out) :: IERR
double precision, intent(out) :: F6

integer :: N, JE, JP, JM, J
character(len=6) :: STRI

IERR = 0
JE = 0
JP = 0
JM = 0
do J=1, 6
    if(LINE(J: J) == '-') JM = J
    if(LINE(J: J) == 'e' .or. LINE(J: J) == 'E') JE = J
    if(LINE(J: J) == '.') JP = J
enddo

if(JP > 0) then
    READ(LINE, '(F6.3)', ERR=2) F6
    return
endif

if(JE <= 0) then
    if(JM <= 1) then
        READ(LINE, '(I6)', ERR=2) N
        F6 = N
        return
    endif
    STRI = LINE(1: JM-1)
    READ(STRI, '(I6)', ERR=2) N
    F6 = N
    STRI = LINE(JM: 6)
    READ(STRI, '(I6)', ERR=2) N
    F6 = F6 * 10.**N
    return
else
    STRI = LINE(1: JE-1)
    READ(STRI, '(I6)', ERR=2) N
    F6 = N
    STRI = LINE(JE+1: 6)
    READ(STRI, '(I6)', ERR=2) N
    F6 = F6 * 10.**N
endif

2 continue

IERR = 1
write(*, *) '>>> READF6: found ERROR in "', LINE, '"'

return
end subroutine READF6

!---------------------------------------------------------------------
character(len=6) function VARNAM(str_in, ierr)
!---------------------------------------------------------------------
! If 1st character is tab or space, a blank string is returned
! If tabs are present anywhere else, they are replaced by blanks 
!    and ierr is set to 1
! Eventual trailing 'X' is removed
! Otherwise, VARNAM is returned.

use outcmn_inc, only: esc_ch, tab_ch

implicit none

integer, intent(out) :: ierr
character(len=*), intent(in) :: str_in

integer :: j, nlen
character(len=1) :: symb

ierr = 0
VARNAM = str_in(1:6)
nlen = LEN_TRIM(VARNAM)

symb = str_in(1: 1)
if (symb == ' ' .or. symb == tab_ch) then
!   Ignore names starting with spaces and tabs
    VARNAM = '      '
    return
endif

do j=2, 6
    symb = str_in(j: j)
    if(symb == tab_ch .or. symb == esc_ch) then
        ierr = 1
        VARNAM(j: j) = ' '
    endif
enddo

if (VARNAM(nlen: nlen) == 'X') VARNAM(nlen: nlen) = ' '

return
end function VARNAM

!---------------------------------------------------------------------
character(len=6) function ARRNAM(str_in)
!---------------------------------------------------------------------
! The subroutine analizes a character*6 "string"
! If the 1st position is tab or space the string 6*' ' is returned
! If tabs are encountered on the end of the "string", 
!  they are removed the "string" is appended with spaces
! Trailing "X" is added when not present in string*6 
! Finally ARRNAM in the Astra standard is created, 
!-----------------------------------------------------------------------

use char_manip, only: clean_string, to_upper

implicit none

character(len=*), intent(in) :: str_in

integer :: nlen
character(len=len(str_in)) :: strtmp

strtmp = to_upper(str_in)
strtmp = clean_string(strtmp)
!call clean_string(strtmp, strtmp)
nlen = LEN_TRIM(strtmp)

ARRNAM = strtmp(1: 6)
! Append "X" if absent

if (nlen < 6 .and. strtmp(nlen: nlen) /= 'X') then
    ARRNAM(nlen+1: nlen+1) = 'X'
endif

return
end function ARRNAM

!---------------------------------------------------------------------
subroutine CHECKU(INTYPE, ABC, AB, XBDRY, YX, jrad, jbdry, STRING, FILENA)
!---------------------------------------------------------------------
! Consistency check for grid array YX(1:jrad) and plasma boundary AB/ABC
! Input:
! ABC - 
! AB - 
! YX(1:jrad) - array for a "radial" coordinate 
! jrad - YX array dimensionality
! STRING - U-file "Independent variable" description
! FILENA - U-file name
! Analyse array YX and returns proper values for XBDRY and jbdry
! Output:
! INTYPE - 
! XBDRY = ABC or AB depending on INTYPE selected
! jbdry - is determined from {YX(jbdry) <= ABC} or {YX(jbdry) <= AB}
!    for {INTYPE = 10} or {INTYPE = 11}, respectively
!  jbdry = jrad if  if YX(jrad) < ABC <= AB
!----------------------------------------------------------------------|

use char_manip, only: to_upper, clean_string

implicit none

integer, intent(in)  ::  jrad
integer, intent(out) :: jbdry, INTYPE
real*4 , intent(in)  :: YX(jrad)
double precision, intent(in)  :: ABC, AB
double precision, intent(out) :: XBDRY
character(len=*), intent(in)   :: FILENA
character(len=*), intent(inout) :: STRING

integer :: j
character(len=12) :: STRAD

jbdry = 0
XBDRY = YX(jrad)
STRAD = STRING(21:31)
STRING = to_upper(clean_string(STRING))
if (STRING(1:8) == 'MINORRAD') then
    INTYPE = 10
elseif (STRING(1:8) == 'MAJORRAD') then
    INTYPE = 19
elseif (STRING(1:3) == 'RHO') then
    INTYPE = 12
elseif (STRING(1:12) == 'POLOIDALFLUX') then
    INTYPE = 13
else
    write(*, *) '>>> U-file "', TRIM(FILENA), '"', &
         '    Unrecognized "radial" variable. Input ignored.' &
        // '    Allowed options are:' &
        // '  Minor Radius        m' &
        // '  Major Radius        m' &
        // '  Rho Toroidal, normalized' &
        // '  Poloidal Flux, normalized'
    INTYPE = -1
    return
endif
STRAD = to_upper(clean_string(STRAD))
if ((INTYPE == 10 .or. INTYPE == 19) .and. STRAD(1: 1) /= 'M') then
    write(*, *) ">>> Warning: Inconsistency in the input data"
    write(*, *) '>>> U-file "', TRIM(FILENA), '":     radial grid is expected to be given in "m"'
endif

if (INTYPE == 10 .and. abs(XBDRY - AB) > 0.3*AB/jrad)  then
    if (XBDRY < AB) then
        jbdry = jrad
    else
        do j=jrad, 1, -1
           if (YX(j) > AB) jbdry = j
        enddo
        XBDRY = AB
    endif
elseif (INTYPE == 11 .and. abs(XBDRY - ABC) > 0.3*ABC/jrad) then
    if (XBDRY < ABC) then
        jbdry = jrad
    else
        do j=jrad, 1, -1
            if (YX(j) > ABC) jbdry = j
        enddo
        XBDRY = ABC
    endif
endif
if (INTYPE == 18) then
    write(*, *)'>>> U-file "', TRIM(FILENA), '"', &
        " Don't know a distance to the major axis.           Set to RTOR"
endif
if (jbdry == 0) then
    jbdry = jrad
endif

return
end subroutine CHECKU
