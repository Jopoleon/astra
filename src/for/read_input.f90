!======================================================================|
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
use outcmn_inc
use expdat
use char_manip, only: to_upper, str_in_list, clean_string
use debugger, only: markloc, debug, astra_stop, flightsim
use parse_utils

use numerical_tools, only: EXTRAP, INTEGR

implicit none

integer, parameter :: MPEX=101, MSIGEX=1, MTEX=50, MSIG=1, MEXT=MPEX*MTEX

logical :: exilog

integer :: IFKEY, jarr, INTYPE, jtype, IM, SYSTEM, jbdry, ntim, ntim1
integer :: jj, j, j0, j1, IERR, ier_tab, jexar, jex1, jpos
integer :: KAB, KABC, KAWALL, KRTOR, KELONM, KTRICH
integer :: XSC0, XSC, n_var, n_color, n_words
integer :: nt_u, nx_u, ios, ndim_u, jvar, jkey, jrt, jt, jthe

double precision :: CHORDN, LINEAV, resize
double precision :: tbeg_nml, tend_nml, tpause_nml
double precision, allocatable :: t_u(:), x_u(:), var_u(:), bnd_rz(:)
double precision :: XBDRY, YB, YB1, YXB, YXB1, ALFA, &
    VRDATA, FACTOR, TIMEVR, VRERR, ROC3A, YTP=-1.d9
character(len=6) :: VNAM, VNAMO, VNAMU, VNAMX, VTIM, VDAT, VERR, VARNAM, ARRNAM, keyword
character(len=30) :: rholbl
character(len=132) :: strarray(10), STRI, lin_upper, dir_path, fname, &
    err_msg, err_format, err_msg_exp, file_in, uname, uvar, win_title

namelist / astra_log / AEXT, AWD, exp_file, equ_file, rev_file, TASK, machine, &
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

resize = 1.
tbeg_nml   = -1.
tend_nml   = -1.
tpause_nml = -1.
NITREQ = 1. ! Initialization: g95 does not like it in blockdata

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
! Read file 'tmp/astra.nml'
!----------------------------------------------------------------------|

call GETENV('expfile', exp_file)
call GETENV('equfile', equ_file)
file_in = 'tmp/' // TRIM(exp_file) // TRIM(equ_file) // '.nml'

OPEN(161, FILE=TRIM(file_in), delim='apostrophe')
READ(161, nml=astra_log, iostat=ios)
CLOSE(161)

! Resize ASTRA frame

frame_wid = resize*frame_wid
frame_hei = resize*frame_hei
XWW = resize*XWW
XWH = resize*XWH
XWX = resize*XWX
XWY = resize*XWY
DXLET = resize*DXLET
DYLET = resize*DYLET
LRJJ = resize*LRJJ

call path_split(rev_file, dir_path, fname, jpos)
if (LEN_TRIM(fname) == 0) fname = 'profil.dat'
! Disallowing user-defined subdirs for Review file
rev_file = '.res/' // TRIM(fname)

call get_runid()
! initvm: initialises graphic window
if (TASK(1: 3) == 'BGD') then
    STRI = 'BGD'//char(0)
    call initvm(XWX, XWY, XWW, XWH, COLTAB, STRI(1: 3), 3)
else
    jj = max(0, (15 + NTOUT - 64)/16)
    XWH = XWH + 2*jj*(DYLET + 2)
    win_title = 'Per aspera ad ASTRA'
    call initvm(XWX, XWY, XWW, XWH, COLTAB, TRIM(win_title), LEN_TRIM(win_title))

    IM = 1
    NST = 0
    MOD10 = 1
    call set_frame(IM, XSC0, XSC)
    call set_plot(IM, XSC0, XSC)

    j = XOUT + .49

    call ASRUMN(j) ! Task menu
    call textbf(0, XWH-104, RUNID, 80) ! Task ID
endif


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

parse_exp_2d: do

    INTYPE = -1
    TIMEVR = .0
    ntim = 0
    jbdry = 0
    factor = 1.
    ALFA = 0.001

    read(201, '(A132)', iostat=ios) STRI
    if (ios < 0) EXIT parse_exp_2d
    if (ios > 0) call astra_stop(err_format)
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
! read only if NCNBT==0, i.e. CCOILX was not defined before
        call read_arrx(201, NCNBT, ntim, NCNB, STRI, CCOILX)
        VNAMO = VNAM

    CASE('VCOILX')
        call read_arrx(201, NCNBT, ntim, NCNB, STRI, VCOILX)
        VNAMO = VNAM

    CASE('DUMCTX')
        call read_arrx(201, NCTPT, ntim, NCTP, STRI, DUMCTX)
        VNAMO = VNAM

    CASE('CTRLMX')
        call read_arrx(201, NCRMT, ntim, NCRM, STRI, CTRLMX)
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
! GIT debug
!open(31, file='out')
!write(31, *) 'DATARR', DATARR(1: jarr)
!write(31, *) 'TIMEX', TIMEX(1:NGR)
!write(31, *) 'NGIRDX', NGRIDX(1:NGR)
!write(31, *) 'NTYPEX', NTYPEX(1:NGR)
!write(31, *) 'KTO', KTO(1:NGR)
!write(31, *) 'KOGDA', KOGDA(KTO(1:NGR))
!write(31, *) 'GDEX', GDEX(1:NGR)
!write(31, *) 'GDEY', GDEY(1:NGR)
!close(31)

39 continue

close(201)

!-----------------------
! End reading "exp" file
!-----------------------

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
jkey = ifkey(258)

GRAP(1:NRW) = AB

TIM7(1) = TINIT
if (TIME > TINIT + 1.025*abs(TSCALE)) TINIT = TSTART
TIM7(3) = abs(TSCALE)/8.

if (TASK(1:3) /= 'BGD') then
    CHORDN = LINEAV()
    call UPSTR(CHORDN, 1./MU(NA))
endif

call inquire_fname('cnf', TRIM(exp_file), TRIM(machine), wall_gc_file)
call inquire_fname('nbi', TRIM(exp_file), TRIM(machine), NBFILE)
call inquire_fname('mse', TRIM(exp_file), TRIM(machine), MSFILE)

return

906 continue
call astra_stop(err_format)

end subroutine read_input
