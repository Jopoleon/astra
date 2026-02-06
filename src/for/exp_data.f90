module exp_data

use parameter_inc, only: NTARR, NRDX

implicit none

integer, parameter :: NTVAR=250000, n_bnd_max=256, nt_bnd_max=1500

type rawScalars
    integer :: nt_all
    integer, dimension(NTVAR) :: var_index=0
    double precision, dimension(NTVAR) :: time=0., data=0., error=0.
    character(len=6), dimension(NTVAR) :: label
endtype rawScalars
type rawProfiles
    integer, dimension(NTARR) :: arr_index=0, jbeg_grid=0, jbeg_data=0, grid_type=0, nrho=0
    double precision, dimension(NTARR) :: time=0., filter=0.001
    character(len=6), dimension(NTARR) :: label
    real*4, dimension(NRDX*NTARR) :: data
endtype rawProfiles
type rawBoundary
    integer :: nt, n_theta
    double precision, dimension(nt_bnd_max) :: time=0.
    double precision, dimension(nt_bnd_max*n_bnd_max) :: R=0., Z=0.
endtype rawBoundary

type(rawScalars)  :: raw_scalars
type(rawProfiles) :: raw_profiles
type(rawBoundary) :: raw_boundary

contains

!----------------------------------------------------------------------
    subroutine read_input
!----------------------------------------------------------------------|
!  NTVAR    maximal number of time slices for all variables
!  IVAR     number of actually defined variables
!  NTARR    maximal number of time slices for all arrays
!  NGR     number of actually defined groups
!  NARRX    maximal number of arrays recognizable from a data file
!----------------------------------------------------------------------|
! The subroutine is called once at the start-up, it reads the "exp" file
! and stores the time evolution of all input data in the arrays
! raw_profiles%data(NRDX*NTARR) - data array
!       Let   1 <= j <= NTARR is an ordinal number of array in raw_profiles%data
! jbeg_arrx(jx)  - pointer to a position in the array raw_profiles%time
!----------------------------------------------------------------------|

    use parameter_inc, only: NRD, NRDX, NTARR
    use const_inc, only: NITREQ, NA, NA1, NB1, NAB, &
        TIME, TSTART, TEND, TPAUSE, TAUMIN, TAUPRP, TINIT, TSCALE, TIMEQL, DTEQL, &
        varValues, constValues, internValues, exp_header, ARXUSE, &
        AB, ABC, AWAll, ROC, ROCO, ROWALL,  HRO, HROX, RTOR, &
        ELONG, ELONM, TRIAN, TRICH, SHIFT, VOLUME, &
        GP, GP2, BTOR, BTN, FTO, FTN, IPL, IPLN, PSIAX, PSIBO
    use status_inc, only: XRHO, SXHO, RHO, SRHO, AMETR, &
        G11, G22, VR, VRO, VRS, VOLUM, &
        FP, FPO, FP_NORM, rho_pol, NE, NEO, TE, TEO, UPAR, UPARO, MRHO, &
        AMAIN, UPS0, UPS0O
    use io_mod, only: exp_file, equ_file, machine, NBfile, CCOILX, VCOILX, &
        IFDFVX, IFDFAX, jbeg_arrx, NGR, n_coils, nt_coils

    use char_manip, only: to_upper, str_in_list, clean_string
    use debugger, only: markloc, debug, astra_stop
    use parse_utils, only: path_split, split2array2, &
        ufheader, ufrd, parse_u_line, inquire_fname, assign_val, read_coilx
    use numerical_tools, only: EXTRAP, INTEGR
    use plasma_state, only: plasma_up
    use json_vars, only: read_metadata, internNames, constNames, varNames, profxNames, &
        n_intern, n_const, n_var, n_profx
    use machine_config, only: config_read

    logical :: log_exists
    integer :: jarr, INTYPE, jtype, jbdry, ntim, ntim1, IVAR
    integer, allocatable, dimension(:) :: int_json
    integer :: jj, j, j0, j1, IERR, ier_tab, jexar, jex1, jpos
    integer :: KAB, KAWALL, KRTOR, KELONM, KTRICH
    integer :: nvar, n_color, n_words, i_filter_glob
    integer :: nt_u, nx_u, ios, ndim_u, jvar, jrt, jt, jthe, nbnd

    double precision, allocatable :: t_u(:), x_u(:), var_u(:), bnd_rz(:), bnd_r(:), bnd_z(:)
    double precision :: XBDRY, YB, YB1, YXB, YXB1, ALFA, ALFA_GLOB, &
        VRDATA, FACTOR, TIMEVR, VRERR, ROC3A, YTP=-1.d9
    character(len=6) :: VNAM, VNAMO, VNAMU, VNAMX, VTIM, VDAT, VERR, keyword
    character(len=31) :: rholbl
    character(len=132) :: strarray(20), STRI, lin_upper, dir_path, fname, &
        err_msg, err_format, err_msg_exp, file_in, uname, uvar, workflow

!----------------------------------------------------------------------|
! Fortran tests

    call markloc('read_input')

!----------------------------------------------------------------------|
! Initialisation with default values

    NITREQ = 1. ! Initialization: g95 does not like it in blockdata
    i_filter_glob = 0 ! if i_filter_glob = 1, a global filter is set

    plasma_up = 1  ! plasma is up by default, can be set to 0 for breakdown by the user in a user-defined sbr called with "<"

!----------------------------------------------------------------------|
! Read tables to set global scalars varNames, constNames, internNames
!----------------------------------------------------------------------|

    call read_metadata

    do j=1, n_var
        SELECT CASE (varNames(j))
        CASE('AB    ')
            KAB    = j
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


!----------------------------------------------------------------------
! Read machine configuration, if available (need "machine" variable defined)
    call config_read()

! Look for NBfile
    call inquire_fname('nbi', TRIM(exp_file), TRIM(machine), NBFILE)

!----------------------------------------------------------------------
! Read file equ/log/<model>
!----------------------------------------------------------------------

    jj = LEN_TRIM(equ_file)
    if (jj == 0) call astra_stop('>>> read_input: Error, empty model file name')

    call path_split(equ_file, dir_path, fname, jpos)

    if (jpos == 0) then ! no subdir
        file_in = 'equ/log/'//TRIM(fname)
    else ! equ/<subdir>/log/<model>
        file_in = 'equ/'//TRIM(dir_path)//'log/'//TRIM(fname)
    endif

    inquire(file=TRIM(file_in), exist=LOG_EXISTS)

    if (.not. LOG_EXISTS)  then ! Missing log file
        call astra_stop('>>> Error: file "' // TRIM(file_in) // '" missing')
    endif
! Read log file
    nvar = 37
    call assign_val(file_in, nvar    ,    varNames(1: nvar)    ,    varValues(1: nvar)    , n_color)
    call assign_val(file_in, n_const ,  constNames(1: n_const) ,  constValues(1: n_const) , n_color)
    call assign_val(file_in, n_intern, internNames(1: n_intern), internValues(1: n_intern), n_color)

    close(171)

!----------------------------------------------------------------------
! Read experimental file
!----------------------------------------------------------------------

    file_in='exp/' // TRIM(exp_file)

    err_msg_exp = '>>> Data file "' // TRIM(exp_file) // '" error:\n    '

    open(201, FILE=TRIM(file_in), iostat=ios)
    if (ios /= 0) call astra_stop('>>> read_input: No such experimental variant "' // TRIM(exp_file) // '"')

    read(201, '(A132/)', iostat=ios) exp_header
    if (ios /= 0) call astra_stop(err_msg_exp // 'in header')

!----------------------------------------------------------------------
! Read simple variable loop (between the labels "5" and "10"):

    VNAMO = ' '
    VNAMU = ' '
    IVAR = 0

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
            jvar = str_in_list(VNAM(1:6), varNames)
            if (jvar == 0) CYCLE parse_exp_1d  ! var doesnt exist
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
            read(201, *, iostat=ios) (raw_scalars%time(IVAR+jj), jj=1, ntim)
            if (ios /= 0) call astra_stop(err_format)
            read(201, *, iostat=ios) (raw_scalars%data(IVAR+jj), jj=1, ntim)
            if (ios /= 0) call astra_stop(err_format)
            varValues(jvar) = factor*raw_scalars%data(IVAR+1)
            do jj=1, ntim
                IVAR = IVAR + 1
                raw_scalars%var_index(IVAR) = jvar
                raw_scalars%data(IVAR) = factor*raw_scalars%data(IVAR)
                raw_scalars%error(IVAR) = 0.
                raw_scalars%label(IVAR) = VNAM
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

        jvar = str_in_list(VNAM, varNames)
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
            varValues(jvar) = factor*VRDATA
            if (VNAM == VNAMO) IFDFVX(jvar) = 1
! Repeated name
            IVAR = IVAR+1
            if (IVAR > NTVAR) then
                write(err_msg, '(A, i)') '>>> read_input: Time dependent variable strings >', NTVAR
                call astra_stop(err_msg)
            endif

            raw_scalars%var_index(IVAR) = jvar
            raw_scalars%time(IVAR) = TIMEVR
            raw_scalars%data(IVAR) = factor*VRDATA
            raw_scalars%error(IVAR) = VRERR
            raw_scalars%label(IVAR) = VNAM

        else  ! ":" found in the input string "STRI", pointer to U-file

            if (VNAM == VNAMO) call astra_stop(err_msg_exp // 'ambiguous ' // TRIM(VNAM) // ' definition')

            call parse_u_line(STRI, uvar, uname, factor)
            VNAMU = VNAM
            call ufheader(TRIM(uname), ndim_u, nt_u, nx_u, rholbl)
            allocate(t_u(nt_u))
            allocate(x_u(nx_u))
            allocate(var_u(nt_u*nx_u))
            call ufrd(TRIM(uname), ndim_u, nt_u, nx_u, t_u, x_u, var_u)

            varValues(jvar) = factor*var_u(1)
            if (nt_u == 1) then
                IFDFVX(jvar) = 0
            else
                IFDFVX(jvar) = 1
            endif
            do jj=1, nt_u
                IVAR = IVAR + 1
                raw_scalars%var_index(IVAR) = jvar
                raw_scalars%time(IVAR) = t_u(jj)
                raw_scalars%data(IVAR) = factor*var_u(jj)
                raw_scalars%error(IVAR) = 0.
                raw_scalars%label(IVAR) = VNAM
            enddo
            deallocate(t_u, x_u, var_u)
        endif

        VNAMO = VNAM

    enddo parse_exp_1d

    raw_scalars%nt_all = IVAR

    close(201)

!----------------------------------------------------------------------
! Start Astra standard assignments

    if (NA1 > NRD) then
        write(err_msg, '(2A, i)') '>>> FATAL ERROR: The radial grid size out of range.\n', &
            '                 Parameter "NA1" cannot exceed', NRD
        call astra_stop(err_msg)
    endif

    TIME = TSTART
    if (YTP > -1.d8) TPAUSE = YTP

    call INTVAR

    if (AWALL < AB) then
         if (IFDFVX(KAWALL) >= 0) write(*, *) '>>> Warning: AWALL < AB.  Setting AWALL = AB'
         AWALL = AB
    endif
    if (IFDFVX(KAWALL) < 0 .and. (AWALL < AB .or. AWALL > 1.2*AB)) AWALL = AB

    if (ABC > AB) then
        err_msg = '>>> Error: ABC cannot exceed AB. Check your file exp/'//TRIM(exp_file)
        call astra_stop(TRIM(err_msg))
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

!----------------------------------------------------------------------
! Start reading radial profiles:  Radial grid: jbdry
!    Option for different input grids can be added

    VNAMO = ' '

    open(201, FILE=TRIM(file_in), iostat=ios)
    read(201, '(A132)', ERR=906, END=39) STRI
    read(201, '(A132)', ERR=906, END=39) STRI

    NGR = 0
    jarr = 0
    raw_boundary%nt = 0
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
        if (lin_upper(1: 3) == 'END') EXIT parse_exp_2d

        VNAM = VARNAM(lin_upper(1: 6), ier_tab)
        if (ier_tab /= 0 .or. vnam == '') CYCLE parse_exp_2d    ! Ignore lines starting with a blank
        VNAMX = ARRNAM(VNAM)
        jex1 = str_in_list(VNAMX, profxNames) ! Checks if VNAM is in array list

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

        jexar = str_in_list(VNAMX, profxNames) ! Checks if VNAMX-string is in array profxNames

        if (jexar > 0) then
            if (VNAM /= VNAMO .and. IFDFAX(jexar) < 0) then
               jbeg_arrx(jexar) = NGR + 1
            endif
        endif

! Special arrays

        SELECT CASE(VNAMX)

        CASE('CCOILX')
            nt_coils = 0
! read only if nt_coils==0, i.e. CCOILX was not defined before
            call read_coilx(201, nt_coils, ntim, n_coils, STRI, CCOILX)
            VNAMO = VNAM

        CASE('VCOILX') !note that both CCOIL and VCOIL need to appear in the exp file with the same number of points and times
            nt_coils = 0
            call read_coilx(201, nt_coils, ntim, n_coils, STRI, VCOILX)
            VNAMO = VNAM

        CASE ('BNDX  ')
            if (raw_boundary%nt /= 0) then
                call astra_stop(err_msg_exp // 'Boundary must be defined in a single group')
            endif
            j = INDEX(lin_upper, 'POINTS')
            if (j /= 0) read(STRI(j+6:),*) raw_boundary%n_theta
            if (j == 0) then
                call astra_stop(err_msg_exp // 'Number of boundary points must be defined')
            endif

            raw_boundary%nt = max(ntim, 1)
            write(*, *) 'Reading BND, dims:', raw_boundary%n_theta, raw_boundary%nt

            if (raw_boundary%n_theta > n_bnd_max) then
                write(err_msg, '(2A, i)') err_msg_exp, 'Boundary data #theta must not exceed ', n_bnd_max
                call astra_stop(err_msg)
            endif
             if (raw_boundary%nt > nt_bnd_max) then
                write(err_msg, '(2A, i)') err_msg_exp, 'Boundary data #times must not exceed ', nt_bnd_max
                call astra_stop(err_msg)
            endif

! Input order:
! t_1  t_2  t_3
! r_1(t_1) r_1(t_2) r_1(t_3)
! z_1(t_1) z_1(t_2) z_1(t_3)
! r_2(t_1) r_2(t_2) r_2(t_3)
! z_2(t_1) z_2(t_2) z_2(t_3)

            read(201, *, iostat=ios) (raw_boundary%time(j), j=1, raw_boundary%nt)
            nbnd = 2*raw_boundary%nt*raw_boundary%n_theta
            allocate(bnd_rz(nbnd))
            read(201, fmt=*, iostat=ios) (bnd_rz(j), j=1, nbnd)
            jrt = 1
            do jthe=1, raw_boundary%n_theta
                do jt=1, raw_boundary%nt
                    raw_boundary%R((jthe-1)*raw_boundary%nt + jt) = bnd_rz(jrt)
                    raw_boundary%Z((jthe-1)*raw_boundary%nt + jt) = bnd_rz(jrt+raw_boundary%nt)
                    jrt = jrt + 1
                enddo
                jrt = jrt + raw_boundary%nt
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
            raw_boundary%n_theta = nx_u
            raw_boundary%nt = nt_u
            if (raw_boundary%n_theta > n_bnd_max) then
                write(err_msg, '(2A, i)') TRIM(err_msg_exp), 'Boundary data #theta must not exceed ', n_bnd_max
                call astra_stop(err_msg)
            endif
            if (raw_boundary%nt > nt_bnd_max) then
                write(err_msg, '(2A, i)') TRIM(err_msg_exp), 'Boundary data #times must not exceed ', nt_bnd_max
                call astra_stop(err_msg)
            endif

            allocate(x_u(nx_u))
            call ufrd('udb/' // trim(STRI) // '_r', ndim_u, nt_u, nx_u, raw_boundary%time(1:nt_u), x_u, raw_boundary%R(1:nx_u))
            call ufrd('udb/' // trim(STRI) // '_z', ndim_u, nt_u, nx_u, raw_boundary%time(1:nt_u), x_u, raw_boundary%Z(1:nx_u))
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
            call CHECKU(INTYPE, ABC, AB, XBDRY, raw_profiles%data(jarr+1), nx_u, jbdry, rholbl, file_in)

            do j=1, jbdry
                raw_profiles%data(jarr + j) = x_u(j)
            enddo
            do j=1, nt_u
                raw_profiles%time(NGR+j) = t_u(j)
            enddo

            if (INTYPE == 18 .or. INTYPE == 19) then
                raw_profiles%data(jarr + 2: jarr + jbdry) = raw_profiles%data(jarr + 1)
                if (INTYPE == 18) raw_profiles%data(jarr + 1) = RTOR
                if (INTYPE == 19) raw_profiles%data(jarr + 1) = 0.
            endif

            jarr = jarr + jbdry

            if (jarr + nx_u + (nt_u - 1)*jbdry > NRDX*NTARR) then
                call astra_stop('>>> read_input: Buffer size exceeded')
            endif

            jrt = jarr
            do jj=1, nt_u
                do j=1, jbdry
                    jrt = jrt + 1
                    raw_profiles%data(jrt) = var_u(jj + (j - 1)*nt_u)
                enddo
            enddo

            deallocate(t_u, x_u, var_u)

            if (INTYPE > 13 .and. INTYPE /= 19) write(*, *) 'Unknown U-file type'
            YXB = raw_profiles%data(jarr)
            YXB1 = raw_profiles%data(jarr-1)
            do j=1, nt_u
                NGR = NGR + 1
                if (j == 1) then
                    raw_profiles%jbeg_grid(NGR) = jarr - jbdry + 1
                    raw_profiles%jbeg_data(NGR) = jarr + 1
                else
                    raw_profiles%jbeg_grid(NGR) = raw_profiles%jbeg_grid(NGR-1)
                    raw_profiles%jbeg_data(NGR) = jarr + 1
                endif
                do j0=1, jbdry
                    raw_profiles%data(jarr + j0) = raw_profiles%data(jarr + j0)*factor
                enddo
                jarr = jarr + jbdry
                if (jbdry /= nx_u) then
                    YB = raw_profiles%data(jarr)
                    YB1 = raw_profiles%data(jarr-1)
                    YB = (YB*(XBDRY - YXB1) - YB1*(XBDRY - YXB))/(YXB - YXB1)
                    raw_profiles%data(jarr) = YB
                endif
                raw_profiles%arr_index(NGR) = jexar
                raw_profiles%label    (NGR) = VNAMX
                raw_profiles%nrho     (NGR) = jbdry
                raw_profiles%grid_type(NGR) = INTYPE
                raw_profiles%filter   (NGR) = ALFA
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

            if (ntim > 0) read(201, *, ERR=906) (raw_profiles%time(NGR+j), j=1, ntim)

            do j=1, ntim1
                NGR = NGR + 1
                raw_profiles%arr_index(NGR) = jexar
                raw_profiles%label    (NGR) = VNAMX
                raw_profiles%nrho     (NGR) = jbdry
                raw_profiles%grid_type(NGR) = INTYPE
                raw_profiles%filter   (NGR) = ALFA
                if (j == 1) then
                    jbeg_arrx(jexar) = NGR
                    raw_profiles%jbeg_grid(NGR) = jarr + 1
                    if (INTYPE == 18 .or. INTYPE == 19) then
                        jarr = jarr + 1
                        read(201, *, ERR=906) raw_profiles%data(jarr)
                    endif
                    do j1=1, jtype
                        read(201, *, iostat=ios) (raw_profiles%data(jarr + jj), jj=1, jbdry)
                        if (ios /= 0) then
                            write(err_msg, '(3A, /, A, 1p, 6e12.4)') err_msg_exp, &
                                '".  Format error in group ', TRIM(STRI), 'Last data read: ', &
                                (raw_profiles%data(jarr+jj), jj=1, jbdry)
                                call astra_stop(err_msg)
                        endif
                        jarr = jarr + jbdry
                        if (jtype == 2 .and. INTYPE <= 17 .and. j1 == 1) XBDRY = raw_profiles%data(jarr)
                    enddo
                    raw_profiles%jbeg_data(NGR) = jarr - jbdry + 1
                    do jj=1, jbdry
                        raw_profiles%data(jarr - jbdry + jj) = factor*raw_profiles%data(jarr - jbdry + jj)
                    enddo
                else
                    raw_profiles%jbeg_grid(NGR) = raw_profiles%jbeg_grid(NGR-1)
                    raw_profiles%jbeg_data(NGR) = jarr + 1
                    read(201, *, iostat=ios) (raw_profiles%data(jarr + jj), jj=1, jbdry)
                    if (ios /= 0) then
                        write(err_msg, '(3A, /, A, 1p, 6e12.4)') err_msg_exp, &
                                '".  Format error in group ', TRIM(STRI), 'Last data read: ', &
                            (raw_profiles%data(jarr+jj), jj=1, jbdry)
                        call astra_stop(err_msg)
                    endif
                    do jj=1, jbdry
                        raw_profiles%data(jarr + jj) = factor*raw_profiles%data(jarr + jj)
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
    if (raw_boundary%nt > 0) then
!find time index of most proximum boundary
        j=1
        do jt=1, raw_boundary%nt
            if (raw_boundary%time(jt) <= TSTART) j = jt
        enddo
        jt = j
        allocate(bnd_r(raw_boundary%n_theta), bnd_z(raw_boundary%n_theta))
        do jthe=1, raw_boundary%n_theta
            bnd_r(jthe) = raw_boundary%R((jthe-1)*raw_boundary%nt + jt)
            bnd_z(jthe) = raw_boundary%Z((jthe-1)*raw_boundary%nt + jt)
        enddo
!calculate ABC
        ABC = (maxval(bnd_r) - minval(bnd_r))/2.
!calculate elong
        ELONG = (maxval(bnd_z) - minval(bnd_z))/(2.*ABC)
        ELONG = max(ELONG, 1.d0)
        deallocate(bnd_r)
        deallocate(bnd_z)
    endif

!assign variables here for initialization:

    VOLUME = GP2*GP*RTOR*AB**2 * ELONG

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

    HROX  = 1.0/(NA1 - 0.5)

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

    call SETGEO
    call NEW_GRID

    do J=1, NB1
        G22(J) = RHO(J)
        VR(J)  = (GP2*(RTOR + SHIFT))**2*RHO(J)/RTOR
        VRS(J) = (GP2*(RTOR + SHIFT))**2*J*HRO/RTOR
        G11(J) = VRS(J)
    enddo

    call INTEGR(RHO, 1, VR, VOLUM, NA1)

    PSIBO = FP(NA1)
    PSIAX = EXTRAP(XRHO(1: NA1), FP(1: NA1), 0.0, NA1, 2, .true.)

    FP_NORM = (FP - PSIAX)/(PSIBO - PSIAX)
    rho_pol = SQRT(FP_NORM)

    VOLUME = VOLUM(NA1)

    NAB = NA1
    if (NA1 < NB1 .and. AB > ABC) then
        do j=NA1+1, NB1
            if (AMETR(j) < AB) NAB = j
        enddo
        if (NAB < NB1) NAB = NAB + 1
    endif

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

    do j=1, n_profx
        if (ARXUSE(j) /= 0) then
            if (IFDFAX(ARXUSE(j)) < 0) then
                write(*, *) '>>> Warning >>> X-array used but not defined: "', &
                    TRIM(profxNames(ARXUSE(j))), '"', ARXUSE(j), IFDFAX(ARXUSE(j))
            endif
        endif
    enddo

    TIMEQL = TIME - DTEQL - 1.d-7
    TAUPRP = TAUMIN
    if (TIME > TINIT + 1.025*abs(TSCALE)) TINIT = TSTART

    return

    906 continue
    call astra_stop(err_format)

    end subroutine read_input

!---------------------------------------------------------------------
    subroutine READF6(LINE, F6, IERR)
! Read a number in LINE to 6-positinal field

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

    use char_manip, only: esc_ch, tab_ch

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

end module exp_data
