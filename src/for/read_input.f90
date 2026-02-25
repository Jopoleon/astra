module read_input

implicit none

type rawScalars
    integer :: nt_all
    integer, dimension(:), allocatable :: var_index
    double precision, dimension(:), allocatable :: time, data, error
    character(len=6), dimension(:), allocatable :: label
endtype rawScalars
type rawProfiles
    integer :: nt_arr, n_groups
    integer, dimension(:), allocatable :: arr_index, jbeg_grid, jbeg_data, grid_type, nrho
    double precision, dimension(:), allocatable :: time, filter
    character(len=6), dimension(:), allocatable :: label
    double precision, dimension(:), allocatable :: data
endtype rawProfiles
type rawBoundary
    integer :: nt, n_theta
    double precision, dimension(:), allocatable :: time, R, Z
endtype rawBoundary
type rawCoils
    integer :: nt=0, ncoils=0
    double precision, dimension(:), allocatable :: time
    double precision, dimension(:), allocatable :: current
endtype rawCoils

type(rawScalars)  :: raw_scalars
type(rawProfiles) :: raw_profiles
type(rawBoundary) :: raw_boundary
type(rawCoils) :: raw_cCoil, raw_vCoil

contains

!----------------------------------------------------------------------
    subroutine readInput

    use machine_config, only: config_read
    use io_mod, only: exp_file, equ_file, machine, NBfile
    use parse_utils, only: path_split, inquire_fname, assign_val
    use json_vars, only: internNames, constNames, varNames, n_intern, n_const
    use const_inc, only: varValues, constValues, internValues
    use debugger, only: astra_stop

    logical :: log_exists
    integer :: jj, jpos, nvar, n_color
    character(len=132) :: file_in, dir_path, fname
    
! Read machine configuration, if available (need "machine" variable defined)
    call config_read()

! Look for NBfile
    call inquire_fname('nbi', TRIM(exp_file), TRIM(machine), NBFILE)

! Read file equ/log/<model>
    jj = LEN_TRIM(equ_file)
    if (jj == 0) call astra_stop('>>> read_equ_log: Error, empty model file name')

    call path_split(equ_file, dir_path, fname, jpos)

    if (jpos == 0) then ! no subdir
        file_in = 'equ/log/' // TRIM(fname)
    else ! equ/<subdir>/log/<model>
        file_in = 'equ/' // TRIM(dir_path) // 'log/' // TRIM(fname)
    endif

    inquire(file=TRIM(file_in), exist=LOG_EXISTS)
    if (.not. LOG_EXISTS)  then ! Missing log file
        call astra_stop('>>> Error: file "' // TRIM(file_in) // '" missing')
    endif

    nvar = 37
    call assign_val(file_in, nvar    ,    varNames(1: nvar)    ,    varValues(1: nvar)    , n_color)
    call assign_val(file_in, n_const ,  constNames(1: n_const) ,  constValues(1: n_const) , n_color)
    call assign_val(file_in, n_intern, internNames(1: n_intern), internValues(1: n_intern), n_color)

! Read exp file
    call read_exp

! ASTRA default assignments
    call astra_assignments

    return
    end subroutine readInput

!------------------------------------------------------------
    subroutine read_coilx(nunit, stri_in, coilx_out)

    use io_mod, only: exp_file, n_coils_max
    use debugger, only: markloc, astra_stop

    integer, parameter :: nt_coils_max=25000

    integer, intent(in) :: nunit
    character(len=*), intent(in) :: stri_in
    type(rawCoils), intent(out) :: coilx_out

    integer :: j, ios, nt, n_coils
    character(132) :: err_msg

    call markloc('read_coilx')

    err_msg =  '>>> Data file "' // TRIM(exp_file) // '" error:\n'

    j = INDEX(stri_in, 'NTIMES')
    if (j == 0) then
        err_msg = TRIM(err_msg) // '    Number of COILSX must be defined'
        call astra_stop(err_msg)
    endif
    read(stri_in(j+6:), *) nt

    j = INDEX(stri_in, 'POINTS')
    if (j == 0) then
        err_msg = TRIM(err_msg) // '    Number of COILSX must be defined'
        call astra_stop(err_msg)
    endif
    read(stri_in(j+6:), *) n_coils
    if (n_coils > n_coils_max) then
        write(err_msg, '(2A, i)') TRIM(err_msg), &
           '    Number of coils must be <', n_coils_max
        call astra_stop(err_msg)
    endif

    coilx_out%nt = nt
    coilx_out%ncoils = n_coils
    if (.not. allocated(coilx_out%time)) then
        allocate(coilx_out%time(nt))
        allocate(coilx_out%current(nt*n_coils))
    endif

    if (n_coils*nt > n_coils_max*nt_coils_max) then
        write(err_msg, '(2A)') TRIM(err_msg), &
           '    COILSX data length must be n_coils_max*nt_coils_max <'
        call astra_stop(err_msg)
    endif

    read(nunit, *, iostat=ios)(coilx_out%time(j), j=1, nt)
    read(nunit, *, iostat=ios)(coilx_out%current(j), j=1, n_coils*nt)
    if (ios /= 0) then
        err_msg = TRIM(err_msg) // '    Size mismatch in COILSX group'
        call astra_stop(err_msg)
    endif

    return
    end subroutine read_coilx

!----------------------------------------------------------------------
    subroutine read_exp
!----------------------------------------------------------------------|
! len_data_max   maximal number of time slices for all arrays
! len_scalars    number of actually defined variables
! len_profs_time number of actually defined groups
!----------------------------------------------------------------------|
! The subroutine is called once at the start-up, it reads the "exp" file
! and stores the time evolution of all input data in the arrays raw_*%*
! jbeg_arrx  - pointer to a position in the array raw_profiles%time
!----------------------------------------------------------------------|

    use parameter_inc, only: NRDX
    use const_inc, only: NA1, AB, ABC, RTOR, varValues, exp_header, TSTART, TEND
    use io_mod, only: exp_file, IFDFVX, IFDFAX, jbeg_arrx, NGR
    use char_manip, only: to_upper, str_in_list
    use debugger, only: markloc, debug, astra_stop
    use parse_utils, only: split2array2
    use json_vars, only: varNames, profxNames

    integer, parameter :: len_data_max=250000, n_unit=201, nbnd_max=400000

    logical :: skip_read=.false.
    integer :: jarr, INTYPE, jtype, nr_exp, ntim, ntim1, n_coils, IVAR
    integer, allocatable, dimension(:) :: int_json
    integer :: jj, j, j0, j1, IERR, ier_tab, jexar, jex1, jpos
    integer :: n_words, i_filter_glob, len_profs_data, len_profs_time, len_scalars
    integer :: nt_u, nx_u, ios, ndim_u, nscal_u, jvar, jrt, jt, jthe, nbnd

    double precision, allocatable :: t_u(:), x_u(:), var_u(:), bnd_rz(:)
    double precision :: XBDRY, YB, YB1, YXB, YXB1, ALFA, ALFA_GLOB, &
        VRDATA, FACTOR, TIMEVR, VRERR
    character(len=6) :: VNAM, VNAMO, VNAMU, VNAMX, VTIM, VDAT, VERR, keyword
    character(len=31) :: rholbl
    character(len=132) :: strarray(20), STRI, lin_upper, &
        err_msg, err_format, err_msg_exp, file_exp, ufile_in, uname, uvar

!----------------------------------------------------------------------|
    call markloc('read_exp')

    file_exp = 'exp/' // TRIM(exp_file)
    err_msg_exp = '>>> Data file "' // TRIM(exp_file) // '" error:\n    '

!-----------------------------------------------
! Setting dimension for raw_scalars%* allocation

    len_scalars = 0
    open(n_unit, FILE=TRIM(file_exp), iostat=ios)
    if (ios /= 0) call astra_stop('>>> read_exp: No such experimental variant "' // TRIM(exp_file) // '"')
    read(n_unit, '(A132/)', iostat=ios) exp_header
    if (ios > 0) call astra_stop(err_format)

    set_dims_1d: do
        read(n_unit, '(A132)', iostat=ios) STRI
        if (ios > 0) EXIT set_dims_1d
        err_format  =  err_msg_exp // '".  Format error in group ' // TRIM(STRI)
        lin_upper = to_upper(STRI)
        if (lin_upper(1: 1) == '!') CYCLE set_dims_1d
        if (lin_upper(1: 3) == 'END') EXIT set_dims_1d

        call split2array2(lin_upper, strarray, n_words)

        do j=1, n_words
            keyword = strarray(j)
            jpos = str_in_list(keyword, (/'POINTS', 'GRIDTY', 'FILTER', 'PROFIL'/))
            if (jpos > 0) EXIT set_dims_1d
        enddo

        if (strarray(1) == 'NAMEXP') then
            VNAM = strarray(2)
            VTIM = strarray(4)
            jvar = str_in_list(VNAM(1:6), varNames)
            if (jvar == 0) CYCLE set_dims_1d  ! var doesnt exist
            read(vtim, *) ntim !from NTIMES
            len_scalars = len_scalars + max(ntim, 1)
            CYCLE set_dims_1d
        else
! Variable name, time, value, error
            VNAM = VARNAM(lin_upper(1: 6), ier_tab)
        endif

! A special treatment required for time independent NA1, TSTART, TEND
        if (VNAM == 'NA1   ' .or. VNAM == 'TSTART' .or. VNAM == 'TEND  ') then
            CYCLE set_dims_1d
        endif

        jvar = str_in_list(VNAM, varNames)
        if (jvar == 0) CYCLE set_dims_1d

        j1 = index(STRI, ':')

        if (j1 == 0) then  ! time-dependent values in exp-file
            len_scalars = len_scalars + 1
        else  ! ":" U-file
            call parse_u_line(STRI, uvar, uname, factor)
            call ufheader(TRIM(uname), nscal_u, ndim_u, nt_u, nx_u, rholbl)
            len_scalars = len_scalars + nt_u
        endif
    enddo set_dims_1d

    if (len_scalars > len_data_max) then
        write(err_msg, '(A, i)') &
            '>>> read_exp: Size of time dependent scalars data stream cannot exceed', len_data_max
        call astra_stop(err_msg)
    endif

    if (.not. allocated(raw_scalars%var_index)) then
        allocate(raw_scalars%var_index(len_scalars))
        allocate(raw_scalars%time(len_scalars))
        allocate(raw_scalars%data(len_scalars))
        allocate(raw_scalars%error(len_scalars))
        allocate(raw_scalars%label(len_scalars))
    endif
    raw_scalars%var_index = 0
    raw_scalars%time = 0.
    raw_scalars%data = 0.
    raw_scalars%error = 0.
    raw_scalars%label = ' '

!------------------------------------------------
! Setting dimension for raw_profiles%* allocation
! No rewinding

    len_profs_data = 0
    len_profs_time = 0
    skip_read = .true.
    set_dims_2d: do
        if (skip_read) then
            skip_read = .false.
        else
            read(n_unit, '(A132)', iostat=ios) STRI
            if (ios < 0) EXIT set_dims_2d
            if (ios > 0) call astra_stop(err_format)
        endif
        lin_upper = to_upper(STRI)
        if (LEN_TRIM(lin_upper) == 0) CYCLE set_dims_2d
        if (lin_upper(1: 1) == '!') CYCLE set_dims_2d
        if (lin_upper(1: 3) == 'END') EXIT set_dims_2d
        VNAM = VARNAM(lin_upper(1: 6), ier_tab)
        if (ier_tab /= 0 .or. VNAM == '') CYCLE set_dims_2d    ! Ignore lines starting with a blank
        VNAMX = ARRNAM(VNAM)
        jex1 = str_in_list(VNAMX, profxNames) ! Checks if VNAM is in array list
        if (jex1 == 0) then ! ASCII-exp
                VNAM = ' '
                call split2array2(lin_upper, strarray, n_words)
                ntim = 1
                do j=1, n_words
                    keyword = strarray(j)(1: 6)
                    SELECT CASE(keyword)
                    CASE('POINTS')
                        read(strarray(j+1), *, iostat=ios) nr_exp
                        if (ios /= 0) call astra_stop(err_format)
                    CASE('NAMEXP')      ! New variable, exp-block
                        VNAM = VARNAM(strarray(j+1), ier_tab)
                    CASE('NTIMES')
                        read(strarray(j+1), *, iostat=ios) ntim
                        if (ios /= 0) call astra_stop(err_format)
                    END SELECT
                enddo
                VNAMX = ARRNAM(VNAM)
                jexar = str_in_list(VNAMX, profxNames)
                if (jexar > 0) then ! only profiles_x
                    len_profs_data = len_profs_data + ntim*nr_exp + 2*nr_exp ! Some margin in case GRIDTYPE=18, 19, 20
                    len_profs_time = len_profs_time + ntim
                endif
        else       ! ufile
            call parse_u_line(STRI, uvar, ufile_in, factor)
            if (LEN_TRIM(ufile_in) > 0) then ! string 'U-FILE' found in this line
                call ufheader(TRIM(ufile_in), nscal_u, ndim_u, nt_u, nx_u, rholbl)
            endif
            len_profs_data = len_profs_data + nt_u*nx_u + nx_u
            len_profs_time = len_profs_time + nt_u
        endif
    enddo set_dims_2d
    close(n_unit)

    if (len_profs_data > len_data_max) then
        write(err_msg, '(A, i)') &
            '>>> read_exp: Size of time dependent profiles data stream cannot exceed', len_data_max
        call astra_stop(err_msg)
    endif

    if (.not. allocated(raw_profiles%data)) then
        allocate(raw_profiles%data(len_profs_data))
        allocate(raw_profiles%time(len_profs_time))
        allocate(raw_profiles%filter(len_profs_time))
        allocate(raw_profiles%label(len_profs_time))
        allocate(raw_profiles%arr_index(len_profs_time))
        allocate(raw_profiles%jbeg_grid(len_profs_time))
        allocate(raw_profiles%jbeg_data(len_profs_time))
        allocate(raw_profiles%grid_type(len_profs_time))
        allocate(raw_profiles%nrho(len_profs_time))
    endif
    raw_profiles%time = 0.
    raw_profiles%filter = 0.001
    raw_profiles%arr_index = 0
    raw_profiles%jbeg_grid = 0
    raw_profiles%jbeg_data = 0
    raw_profiles%grid_type = 0
    raw_profiles%nrho = 0

!    print*, 'exp file array lengths: ', len_scalars, len_profs_time, len_profs_data
!    pause

!----------------
! Read 1D data
! Rewind exp file
!----------------

    open(n_unit, FILE=TRIM(file_exp), iostat=ios)
    read(n_unit, '(A132/)', iostat=ios) exp_header
    if (ios /= 0) call astra_stop(err_msg_exp // 'in header')

    VNAMO = ' '
    VNAMU = ' '
    IVAR = 0
    i_filter_glob = 0 ! if i_filter_glob = 1, a global filter is set

! Parse scalar block
    parse_exp_1d: do

        read(n_unit, '(A132)', iostat=ios) STRI
        if (ios < 0) then
            close(n_unit)
            raw_profiles%n_groups = NGR
            return
        endif
        err_format  =  err_msg_exp // '".  Format error in group ' // TRIM(STRI)
        lin_upper = to_upper(STRI)
        if (lin_upper(1: 1) == '!') CYCLE parse_exp_1d
        if (lin_upper(1: 3) == 'END') then
            close(n_unit)
            raw_profiles%n_groups = NGR
            return
        endif           

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
            if (ntim <= 1) then
                ntim = 1
                IFDFVX(jvar) = 0
            else
                IFDFVX(jvar) = 1
            endif
! Read "time" array & function array
            read(n_unit, *, iostat=ios) (raw_scalars%time(IVAR+jj), jj=1, ntim)
            if (ios /= 0) call astra_stop(err_format)
            read(n_unit, *, iostat=ios) (raw_scalars%data(IVAR+jj), jj=1, ntim)
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
            call str2dbl(VDAT, VRDATA, IERR)
            VNAMO = VNAM
            if (VNAM == 'NA1   ') NA1 = VRDATA
            if (VNAM == 'TSTART') then
                TSTART = VRDATA
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

        j1 = index(STRI, ':')

        if (j1 == 0) then  ! time-dependent values in exp-file

! Tabulation encountered in the old-standard line:
            if (ier_tab /= 0 ) then
                err_msg = err_msg_exp // '"' // TRIM(VNAM) // &
                    '": tabulation not allowed in this type of input'
                call astra_stop(err_msg)
            endif

            err_msg = '>>> read_exp: File "' // TRIM(file_exp) // '" reading error'
            call str2dbl(VTIM, TIMEVR, IERR)
            if (IERR /= 0) call astra_stop(err_msg)
            call str2dbl(VDAT, VRDATA, IERR)
            if (IERR /= 0) call astra_stop(err_msg)
            call str2dbl(VERR, VRERR, IERR)
            if (IERR /= 0) call astra_stop(err_msg)

            IFDFVX(jvar) = 0
            varValues(jvar) = factor*VRDATA
            if (VNAM == VNAMO) IFDFVX(jvar) = 1
! Repeated name
            IVAR = IVAR + 1
            raw_scalars%var_index(IVAR) = jvar
            raw_scalars%time(IVAR) = TIMEVR
            raw_scalars%data(IVAR) = factor*VRDATA
            raw_scalars%error(IVAR) = VRERR
            raw_scalars%label(IVAR) = VNAM

        else  ! ":" found in the input string "STRI", pointer to U-file

            if (VNAM == VNAMO) call astra_stop(err_msg_exp // 'ambiguous ' // TRIM(VNAM) // ' definition')

            call parse_u_line(STRI, uvar, uname, factor)
            VNAMU = VNAM
            call ufheader(TRIM(uname), nscal_u, ndim_u, nt_u, nx_u, rholbl)
            call ufrd(TRIM(uname), nscal_u, ndim_u, nt_u, nx_u, t_u, x_u, var_u)

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
        endif

        VNAMO = VNAM

    enddo parse_exp_1d

    raw_scalars%nt_all = IVAR

!-------------
! Read 2D data
! No rewinding    
!-------------

    ios = 0
    VNAMO = ' '
    NGR = 0
    jarr = 0
    raw_boundary%nt = 0
    ALFA_GLOB = 0.001

    skip_read = .true.
    parse_exp_2d: do
        if (skip_read) then
            skip_read = .false.
        else
            read(n_unit, '(A132)', iostat=ios) STRI
            if (ios < 0) EXIT parse_exp_2d
            if (ios > 0) call astra_stop(err_format)
        endif

        INTYPE = -1
        TIMEVR = .0
        ntim = 0
        nr_exp = 0
        factor = 1.
        ALFA = ALFA_GLOB

        if (STRI(1:6) == 'FILTER') i_filter_glob = 1
        lin_upper = to_upper(STRI)
        if (LEN_TRIM(lin_upper) == 0) CYCLE parse_exp_2d
        if (lin_upper(1: 1) == '!')   CYCLE parse_exp_2d
        if (lin_upper(1: 3) == 'END')  EXIT parse_exp_2d

        VNAM = VARNAM(lin_upper(1: 6), ier_tab)
        if (ier_tab /= 0 .or. vnam == '') CYCLE parse_exp_2d    ! Ignore lines starting with a blank
        VNAMX = ARRNAM(VNAM)
        jex1 = str_in_list(VNAMX, profxNames) ! Checks if VNAM is in array list

        if (jex1 == 0) then ! exp ASCII
            VNAM = ' '
            call split2array2(lin_upper, strarray, n_words)

            do j=1, n_words, 2
                keyword = strarray(j)(1: 6)
                SELECT CASE(keyword)
                CASE('POINTS')
                    read(strarray(j+1), *, iostat=ios) nr_exp
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
                    write(*, *) '>>> read_exp error unknown key word in string:'
                    write(*, *) TRIM(STRI)
                    EXIT
                END SELECT
            enddo

            if (LEN_TRIM(VNAM) == 0) then
                CYCLE parse_exp_2d
            else
                VNAMO = VNAM
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
! read only if nt_coils==0, i.e. CCOILX was not defined before
            if (raw_cCoil%nt == 0) then
                call read_coilx(n_unit, STRI, raw_cCoil)
            endif
            VNAMO = VNAM

        CASE('VCOILX') !note that both CCOIL and VCOIL need to appear in the exp file with the same number of points and times
            if (raw_vCoil%nt == 0) then
                call read_coilx(n_unit, STRI, raw_vCoil)
            endif
            VNAMO = VNAM

        CASE ('BNDX  ')
            if (raw_boundary%nt /= 0) then
                call astra_stop(err_msg_exp // 'Boundary must be defined in a single group')
            endif
            j = INDEX(lin_upper, 'POINTS')
            if (j /= 0) read(STRI(j+6:), *) raw_boundary%n_theta
            if (j == 0) then
                call astra_stop(err_msg_exp // 'Number of boundary points must be defined')
            endif

            raw_boundary%nt = max(ntim, 1)
            write(*, *) 'Reading BND, dims:', raw_boundary%n_theta, raw_boundary%nt

            nbnd = raw_boundary%nt*raw_boundary%n_theta
            if (nbnd > nbnd_max) then
                write(err_msg, '(2A, i)') err_msg_exp, 'Boundary data must not exceed ', nbnd_max
                call astra_stop(err_msg)
            endif

! Input order:
! t_1  t_2  t_3
! r_1(t_1) r_1(t_2) r_1(t_3)
! z_1(t_1) z_1(t_2) z_1(t_3)
! r_2(t_1) r_2(t_2) r_2(t_3)
! z_2(t_1) z_2(t_2) z_2(t_3)

            allocate(raw_boundary%time(raw_boundary%nt))
            read(n_unit, *, iostat=ios) (raw_boundary%time(j), j=1, raw_boundary%nt)
            allocate(raw_boundary%R(nbnd), raw_boundary%Z(nbnd))
            allocate(bnd_rz(2*nbnd))
            read(n_unit, fmt=*, iostat=ios) (bnd_rz(j), j=1, 2*nbnd)
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

        CASE('BNDUX ')
            write(*, *) 'Reading BND from u file'
            read(n_unit, '(A)', iostat=ios) STRI ! u-file name in exp-file
            if (ios < 0) EXIT parse_exp_2d
            call ufheader('udb/'//trim(STRI)//'_r', nscal_u, ndim_u, nt_u, nx_u, rholbl)
            raw_boundary%n_theta = nx_u
            raw_boundary%nt = nt_u
            nbnd = nx_u*nt_u
            if (nbnd > nbnd_max) then
                write(err_msg, '(2A, i)') TRIM(err_msg_exp), 'Boundary data #theta must not exceed ', nbnd_max
                call astra_stop(err_msg)
            endif
            call ufrd('udb/' // trim(STRI) // '_r', nscal_u, ndim_u, nt_u, nx_u, raw_boundary%time, x_u, raw_boundary%R)
            call ufrd('udb/' // trim(STRI) // '_z', nscal_u, ndim_u, nt_u, nx_u, raw_boundary%time, x_u, raw_boundary%Z)
            VNAMO = VNAM

        CASE('ENDX  ')
            EXIT parse_exp_2d

        CASE DEFAULT
            if (jexar == 0) then
                write(*, *) 'Array name ' // TRIM(VNAMX) // ' not found, skipping line'
            endif
            
        END SELECT
 
        if (jexar == 0) then ! Variable name not in profx list
            CYCLE parse_exp_2d
        endif

! Check for U-file command string
        call parse_u_line(STRI, uvar, ufile_in, factor)

        if (LEN_TRIM(ufile_in) > 0) then ! string 'U-FILE' found in this line

            call ufheader(TRIM(ufile_in), nscal_u, ndim_u, nt_u, nx_u, rholbl)
            call ufrd(TRIM(ufile_in), nscal_u, ndim_u, nt_u, nx_u, t_u, x_u, var_u)

            if (nx_u > NRDX) then
                write(err_msg, '(3A, i, A, i)') '>>> U-file "', TRIM(ufile_in), &
                    '" error: radial grid size ', nx_u, ' is larger than', NRDX
                call astra_stop(err_msg)
            endif

            call CHECKU(INTYPE, ABC, AB, XBDRY, raw_profiles%data(jarr+1), nx_u, nr_exp, rholbl, ufile_in)

            do j=1, nr_exp
                raw_profiles%data(jarr + j) = x_u(j)
            enddo
            do j=1, nt_u
                raw_profiles%time(NGR+j) = t_u(j)
            enddo

            if (INTYPE == 18 .or. INTYPE == 19) then
                raw_profiles%data(jarr + 2: jarr + nr_exp) = raw_profiles%data(jarr + 1)
                if (INTYPE == 18) raw_profiles%data(jarr + 1) = RTOR
                if (INTYPE == 19) raw_profiles%data(jarr + 1) = 0.
            endif

            jarr = jarr + nr_exp
            jrt = jarr
            do jj=1, nt_u
                do j=1, nr_exp
                    jrt = jrt + 1
                    raw_profiles%data(jrt) = var_u(jj + (j - 1)*nt_u)
                enddo
            enddo

            if (INTYPE > 13 .and. INTYPE /= 19) write(*, *) 'Unknown U-file type'
            YXB = raw_profiles%data(jarr)
            YXB1 = raw_profiles%data(jarr-1)
            do j=1, nt_u
                NGR = NGR + 1
                if (j == 1) then
                    raw_profiles%jbeg_grid(NGR) = jarr - nr_exp + 1
                    raw_profiles%jbeg_data(NGR) = jarr + 1
                else
                    raw_profiles%jbeg_grid(NGR) = raw_profiles%jbeg_grid(NGR-1)
                    raw_profiles%jbeg_data(NGR) = jarr + 1
                endif
                do j0=1, nr_exp
                    raw_profiles%data(jarr + j0) = raw_profiles%data(jarr + j0)*factor
                enddo
                jarr = jarr + nr_exp
                if (nr_exp /= nx_u) then
                    YB = raw_profiles%data(jarr)
                    YB1 = raw_profiles%data(jarr-1)
                    YB = (YB*(XBDRY - YXB1) - YB1*(XBDRY - YXB))/(YXB - YXB1)
                    raw_profiles%data(jarr) = YB
                endif
                raw_profiles%arr_index(NGR) = jexar
                raw_profiles%label    (NGR) = VNAMX
                raw_profiles%nrho     (NGR) = nr_exp
                raw_profiles%grid_type(NGR) = INTYPE
                raw_profiles%filter   (NGR) = ALFA
            enddo

        else ! read exp-block data
            if (nr_exp <= 1) then
                write(err_msg, '(3A, 8X, A)') err_msg, 'Input quantity: ', TRIM(VNAM), &
                    'Number of grid points must be > 1'
                call astra_stop(err_msg)
            endif
            if (nr_exp > NRDX) then
                write(err_msg, '(2A, i)') err_msg_exp, 'Number of radial points > ', NRDX
                call astra_stop(err_msg)
            endif

            jtype = 1 + INTYPE/10
            ntim1 = max(ntim, 1)

! INTYPE unknown
            if (INTYPE < 0 .or. INTYPE > 20) then
                write(err_msg, '(A, i0, A)') '>>> ERROR: Unknown input type =', INTYPE, ',  ignored'
                call astra_stop(err_msg)
            endif

            if (ntim > 0) read(n_unit, *, iostat=ios) (raw_profiles%time(NGR+j), j=1, ntim)
            if (ios > 0) call astra_stop(err_format)

            do j=1, ntim1
                NGR = NGR + 1
                raw_profiles%arr_index(NGR) = jexar
                raw_profiles%label    (NGR) = VNAMX
                raw_profiles%nrho     (NGR) = nr_exp
                raw_profiles%grid_type(NGR) = INTYPE
                raw_profiles%filter   (NGR) = ALFA
                if (j == 1) then
                    jbeg_arrx(jexar) = NGR
                    raw_profiles%jbeg_grid(NGR) = jarr + 1
                    if (INTYPE == 18 .or. INTYPE == 19) then
                        jarr = jarr + 1
                        read(n_unit, *, iostat=ios) raw_profiles%data(jarr)
                        if (ios > 0) call astra_stop(err_format)
                    endif
                    do j1=1, jtype
                        read(n_unit, *, iostat=ios) (raw_profiles%data(jarr + jj), jj=1, nr_exp)
                        if (ios /= 0) then
                            write(err_msg, '(3A, /, A, 1p, 6e12.4)') err_msg_exp, &
                                '".  Format error in group ', TRIM(STRI), 'Last data read: ', &
                                (raw_profiles%data(jarr+jj), jj=1, nr_exp)
                                call astra_stop(err_msg)
                        endif
                        jarr = jarr + nr_exp
                        if (jtype == 2 .and. INTYPE <= 17 .and. j1 == 1) XBDRY = raw_profiles%data(jarr)
                    enddo
                    raw_profiles%jbeg_data(NGR) = jarr - nr_exp + 1
                    do jj=1, nr_exp
                        raw_profiles%data(jarr - nr_exp + jj) = factor*raw_profiles%data(jarr - nr_exp + jj)
                    enddo
                else
                    raw_profiles%jbeg_grid(NGR) = raw_profiles%jbeg_grid(NGR-1)
                    raw_profiles%jbeg_data(NGR) = jarr + 1
                    read(n_unit, *, iostat=ios) (raw_profiles%data(jarr + jj), jj=1, nr_exp)
                    if (ios /= 0) then
                        write(err_msg, '(3A, /, A, 1p, 6e12.4)') err_msg_exp, &
                                '".  Format error in group ', TRIM(STRI), 'Last data read: ', &
                            (raw_profiles%data(jarr+jj), jj=1, nr_exp)
                        call astra_stop(err_msg)
                    endif
                    do jj=1, nr_exp
                        raw_profiles%data(jarr + jj) = factor*raw_profiles%data(jarr + jj)
                    enddo
                    jarr = jarr + nr_exp
                endif
            enddo
            if (INTYPE == 17) write(*, *) 'X-axis analysis is not implemented GRIDTYPE =', INTYPE
        endif

        VNAMO = VNAM
        IFDFAX(jexar) = 0

    enddo parse_exp_2d

    close(n_unit)

    raw_profiles%n_groups = NGR

    return
    end subroutine read_exp

!----------------------------------------------------------------------
    subroutine astra_assignments

    use parameter_inc, only: NRD
    use const_inc, only: NA1, NA, NB1, NAB, AB, ABC, AWALL, TIME, TSTART, TPAUSE, &
         TAUMIN, TAUPRP, VOLUME, IPL, IPLN, HRO, HROX, ROC, ROCO, BTN, FTO, FTN, &
         GP, GP2, PSIAX, PSIBO, RTOR, BTOR, SHIFT, ROWALL, ELONM, ELONG, TRICH, TRIAN, &
         ARXUSE, TINIT, TSCALE, TIMEQL, DTEQL
    use status_inc, only: XRHO, SXHO, RHO, SRHO, AMETR, &
        G11, G22, VR, VRO, VRS, VOLUM, &
        FP, FPO, FP_NORM, rho_pol, NE, NEO, TE, TEO, UPAR, UPARO, MRHO, &
        AMAIN, UPS0, UPS0O
    use json_vars, only: varNames, n_profx, profxNames, n_var
    use io_mod, only: IFDFVX, IFDFAX, exp_file
    use numerical_tools, only: EXTRAP, INTEGR
    use debugger, only: astra_stop
     
    integer :: KAB, KAWALL, KRTOR, KELONM, KTRICH
    integer :: j, jt, jthe
    double precision :: YTP=-1.d9
    double precision, dimension(:), allocatable :: bnd_r, bnd_z
    character(len=132) :: err_msg
    double precision, external :: ROC3A

    if (NA1 > NRD) then
        write(err_msg, '(2A, i)') '>>> FATAL ERROR: The radial grid size out of range.\n', &
            '                 Parameter "NA1" cannot exceed', NRD
        call astra_stop(err_msg)
    endif

    TIME = TSTART
    if (YTP > -1.d8) TPAUSE = YTP

    call INTVAR

    do j=1, n_var
        SELECT CASE(varNames(j))
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
 
    if (AWALL < AB) then
         if (IFDFVX(KAWALL) >= 0) write(*, *) '>>> Warning: AWALL < AB.  Setting AWALL = AB'
         AWALL = AB
    endif
    if (IFDFVX(KAWALL) < 0 .and. (AWALL < AB .or. AWALL > 1.2*AB)) AWALL = AB

    if (ABC > AB) then
        err_msg = '>>> Error: ABC cannot exceed AB. Check your file exp/' // TRIM(exp_file)
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

! If boundary is given, calculates initial geometry from that
    if (raw_boundary%nt > 0) then
! Find time index of most proximum boundary
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
! Calculate ABC
        ABC = (maxval(bnd_r) - minval(bnd_r))/2.
! Calculate elong
        ELONG = (maxval(bnd_z) - minval(bnd_z))/(2.*ABC)
        ELONG = max(ELONG, 1.d0)
        deallocate(bnd_r)
        deallocate(bnd_z)
    endif

! Assign variables here for initialization:

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
    end subroutine astra_assignments

!---------------------------------------------------------------------
    subroutine str2dbl(line, dbl_out, ierr)

    character(len=*), intent(in)  :: line
    double precision, intent(out) :: dbl_out
    integer, intent(out)          :: ierr

    if (len_trim(line) == 0) then
        dbl_out = 0.0d0
        ierr    = 0
        return
    endif

    read(line, *, iostat=ierr) dbl_out
    end subroutine str2dbl

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
        if (symb == tab_ch .or. symb == esc_ch) then
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
    subroutine CHECKU(INTYPE, ABC, AB, XBDRY, YX, jrad, nr_exp, STRING, FILENA)
!---------------------------------------------------------------------
! Consistency check for grid array YX(1:jrad) and plasma boundary AB/ABC
! Input:
! ABC -
! AB -
! YX(1:jrad) - array for a "radial" coordinate
! jrad - YX array dimensionality
! STRING - U-file "Independent variable" description
! FILENA - U-file name
! Analyse array YX and returns proper values for XBDRY and nr_exp
! Output:
! INTYPE -
! XBDRY = ABC or AB depending on INTYPE selected
! nr_exp - is determined from {YX(nr_exp) <= ABC} or {YX(nr_exp) <= AB}
!    for {INTYPE = 10} or {INTYPE = 11}, respectively
!  nr_exp = jrad if  if YX(jrad) < ABC <= AB
!----------------------------------------------------------------------|

    use char_manip, only: to_upper, clean_string

    integer, intent(in)  ::  jrad
    integer, intent(out) :: nr_exp, INTYPE
    double precision, intent(in)  :: YX(jrad)
    double precision, intent(in)  :: ABC, AB
    double precision, intent(out) :: XBDRY
    character(len=*), intent(in)   :: FILENA
    character(len=*), intent(inout) :: STRING

    integer :: j
    character(len=12) :: STRAD

    nr_exp = 0
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
            nr_exp = jrad
        else
            do j=jrad, 1, -1
               if (YX(j) > AB) nr_exp = j
            enddo
            XBDRY = AB
        endif
    elseif (INTYPE == 11 .and. abs(XBDRY - ABC) > 0.3*ABC/jrad) then
        if (XBDRY < ABC) then
            nr_exp = jrad
        else
            do j=jrad, 1, -1
                if (YX(j) > ABC) nr_exp = j
            enddo
            XBDRY = ABC
        endif
    endif
    if (INTYPE == 18) then
        write(*, *)'>>> U-file "', TRIM(FILENA), '"', &
            " Don't know a distance to the major axis.           Set to RTOR"
    endif
    if (nr_exp == 0) then
        nr_exp = jrad
    endif

    return
    end subroutine CHECKU

!------------------------------------------------------------
    subroutine parse_u_line(str_in, var_name, uname, factor)

    use debugger, only: markloc, astra_stop
    use char_manip, only: to_upper, split_string
    use parse_utils, only: split2array 

    character(len=*), intent(in) :: str_in
    character(len=len(str_in)), intent(out) :: var_name, uname
    double precision, intent(out) :: factor

    integer :: n_words
    character(len=132) :: str1, strarray(20), err_msg

    call markloc('parse_u_line')

    uname = repeat(' ', 40)
    err_msg = 'Error in exp-file line ' // TRIM(str_in)

    str1 = repeat(' ', 132)

    call split_string(str_in, ' ', var_name, str1)
    if (LEN_TRIM(var_name) == 0) call astra_stop(err_msg)

    call split2array(str1, ':', strarray, n_words)

    if (to_upper(TRIM(strarray(1))) /= 'U-FILE' ) return

    if (LEN_TRIM(strarray(2)) > 0) then
        uname = 'udb/' // TRIM(strarray(2))
    else
        call astra_stop(err_msg // ' no u-file name found')
    endif

    factor = 1.d0
    if (n_words == 3) then
        if (LEN_TRIM(strarray(3)) > 0) then
            read(strarray(3), *) factor
        endif
    endif

    return
    end subroutine parse_u_line

!------------------------------------------------------------
    subroutine ufheader(uname, n_scal, n_dim, nt, nx, lbl2)

    use debugger, only: markloc, astra_stop
    use char_manip, only: to_upper, split_string

    integer, parameter :: n_unit=11

    integer, intent(out) :: nt, nx, n_scal, n_dim
    character(len=*), intent(in) :: uname
    character(len=30), intent(out) :: lbl2

    integer :: ios, j, n_shot
    character(132) :: err_msg
    character(32) :: STRI
    character(30) :: lbl1, lbl3, var1_lbl, unit1
    character(4) :: sdev

    var1_lbl = repeat(' ', 30)
    unit1    = repeat(' ', 30)

    open(n_unit, FILE=TRIM(uname), iostat=ios)
    err_msg = '>>> U-file "' // TRIM(uname) // '" opening error'
    if (ios /= 0) call astra_stop(err_msg)

! # shot, device, #dimensions
    read(n_unit, '(i7, A4, 1X, i1)', ERR=925) n_shot, sdev, n_dim
    err_msg = '>>> U-file "' // TRIM(uname) // '" wrong dims'
    if (n_dim <= 0 .or. n_dim > 2) call astra_stop(err_msg)
 
    read(n_unit, *) ! Shot date

! Scalar parameters
    read(n_unit, *, ERR=925) n_scal
    if (n_scal > 0) then
        do j=1, 2*n_scal
            read(n_unit, *)
        enddo
    endif

! 1st independent variable label: X-
    read(n_unit, '(A32)', ERR=925) STRI
    STRI = ADJUSTL(STRI)
    lbl1 = to_upper(STRI(1:30))

    call split_string(TRIM(lbl1), ' ', var1_lbl, unit1)
    if (var1_lbl(1:4) /= 'TIME') then
        write(*, *) '>>> U-file "', TRIM(uname), '" 1st independent variable should be time'
        close(n_unit)
        return
    endif
    if (TRIM(unit1) /= 'SECONDS') then
        write(*, *) '>>> U-file "', TRIM(uname), '" time unit should be second, not ', TRIM(unit1)
    endif

! 2nd independent variable label: Y-
    if (n_dim == 2) then
        read(n_unit, '(A32)', ERR=925) STRI
        STRI = ADJUSTL(STRI)
        lbl2 = to_upper(STRI(1:30))
    endif

! Function label
    read(n_unit, '(A32/)', ERR=925) STRI
    STRI = ADJUSTL(STRI)
    lbl3 = to_upper(STRI(1:30))

! Dimensions
    read(n_unit, *, ERR=925) nt
    if (n_dim == 2) then
        read(n_unit, *, ERR=925) nx
    else
        nx = 1
    endif
    close(n_unit)

    return

925 call astra_stop(err_msg)

    end subroutine ufheader

!------------------------------------------------------------
    subroutine ufrd(uname, n_scal, n_dim, nt, nx, t_out, x_out, arr_out)

    integer, parameter :: n_unit=11

    character(len=*), intent(in) :: uname
    integer, intent(in) :: nt, nx, n_scal, n_dim
    double precision, intent(out), dimension(:), allocatable :: t_out, x_out, arr_out

    integer :: ios, jx, jt, jlin, n_header
    character(132) :: err_msg

    write(*, *) 'Reading u-file ' // TRIM(uname)
    allocate(t_out(nt), x_out(nx), arr_out(nt*nx))

    open(n_unit, FILE=TRIM(uname), iostat=ios)

! Skip header
    n_header = 2*n_scal + 5 + 2*n_dim
    do jlin=1, n_header
        read(n_unit, *)
    enddo

! Read data
    read(n_unit, *) (t_out(jt), jt=1, nt)
    if (n_dim == 2) read(n_unit, *) (x_out(jx) , jx=1, nx)
    read(n_unit, '(1X, 6E13.6)') ((arr_out(jt + (jx - 1)*nt), jt=1, nt), jx=1, nx)
    close(n_unit)

    return
    end subroutine ufrd

end module read_input
