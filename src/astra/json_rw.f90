module json_rw

use json_module, only: json_file, json_core, json_value, json_ck

implicit none

integer, parameter :: nunit=25
integer :: jid

contains

!---------------------------------------------------------------------
    subroutine json_load(json_in, fjson_out)

    character(len=*), intent(in) :: json_in
    type(json_file), intent(out) :: fjson_out

    logical :: status_ok
    character(len=:), allocatable :: error_msg

    call fjson_out%initialize()
    call fjson_out%load(filename=trim(json_in))

    if (fjson_out%failed()) then
        print*, 'Error: '
        call fjson_out%check_for_errors(status_ok, error_msg)    
        write(*, *) 'Error: ', error_msg
        return
    endif

    end subroutine json_load

!---------------------------------------------------------------------
    subroutine read_scalar_block(fjson_in, label, values)

    type(json_file), intent(inout) :: fjson_in
    character(len=*), intent(in) :: label
    double precision, intent(out), allocatable, dimension(:) :: values

    logical :: found
    integer :: j, nvars
    double precision :: val
    type(json_value), pointer :: jsonOut, dictPointer
    type(json_core) :: jCore    

    write(*, '(2A)') 'Reading ', TRIM(label)

    call fjson_in%info(TRIM(label), n_children=nvars)
    allocate(values(nvars))

    call fjson_in%get(TRIM(label), jsonOut, found)
    do j=1, nvars
        call jCore%get_child(jsonOut, j, dictPointer, found)
        call jCore%get(dictPointer, val)
        values(j) = val
    enddo

    return
    end subroutine read_scalar_block

!---------------------------------------------------------------------
    subroutine read_array_block(fjson_in, label, profs)

    type(json_file), intent(inout) :: fjson_in
    character(len=*), intent(in) :: label
    double precision, intent(out), allocatable :: profs(:, :)

    logical :: found
    integer :: j, k, nvars, npts
    double precision :: val
    character(len=:), allocatable :: error_msg

    type(json_value), pointer :: jsonOut
    type(json_value), pointer :: arrPointer
    type(json_value), pointer :: valPointer
    type(json_core) :: jCore

    write(*, '(2A)') 'Reading ', trim(label)

! get "profiles"
    call fjson_in%get(trim(label), jsonOut, found)

! number of arrays (arr1, arr2, ...)
    call jCore%info(jsonOut, n_children=nvars)

! get first array to get length
    call jCore%get_child(jsonOut, 1, arrPointer, found)
    call jCore%info(arrPointer, n_children=npts)

    print*, 'DIMS', npts, nvars
    allocate(profs(npts, nvars))

! read values
    do j=1, nvars
        call jCore%get_child(jsonOut, j, arrPointer, found)
        do k=1, npts
            call jCore%get_child(arrPointer, k, valPointer, found)
            call jCore%get(valPointer, val)
            profs(k, j) = val
        enddo
    enddo

    return
    end subroutine read_array_block

!---------------------------------------------------------------------
    subroutine write_json

    use parameters_a2equil, only: equil_now
    use const_inc, only: NA1, varValues, varxValues, constValues, internValues, intern2Values
    use status_inc, only: profiles, profiles_x
    use io_mod, only: awd, exp_file, equ_file
    use debugger, only: debug
    use json_vars, only: equil_sigPtr, equil_profPtr, equil_rectPtr, equil_coordPtr, &
        profPtr, profxPtr, constPtr, internPtr, intern2Ptr, varPtr, varxPtr, n_prof, n_profx

    integer :: j, jrho, ios, j_call=1, nrho_surf, nthe_surf, nR, nZ
    character(len=180) :: json_out
    double precision, dimension(20) :: equil_traces
    character(KIND=JSON_CK, len=:), allocatable :: sunit, sdesc, sname
    type(json_core) :: jCore

    save j_call

    nrho_surf = SIZE(equil_now%profiles_1d%rho_tor_norm)
    nthe_surf = SIZE(equil_now%coord_sys%position%teta2d)
    nR = SIZE(equil_now%eqgeometry%rectgrid%r2d)
    nZ = SIZE(equil_now%eqgeometry%rectgrid%z2d)

    if (debug > 0) write(*, '(A, 3i)') 'Starting a2json', j_call, nrho_surf, nthe_surf

    write(json_out, '(5A, i0, A)') TRIM(awd), '/ncdf_out/', TRIM(exp_file), TRIM(equ_file), '-', j_call, '.json'

    open(nunit, file=TRIM(json_out), iostat=ios)
    write(nunit, '(A/)') '{'

!------
! ASTRA
!------

! Scalars
    call write_scalar_block(varPtr, varValues, label="variables")
    call write_scalar_block(varxPtr, varxValues, label="variables_x")
    call write_scalar_block(constPtr, constValues, label="constants")
    call write_scalar_block(internPtr, internValues, label="internal")
    call write_scalar_block(intern2Ptr, intern2Values, label="intern2")

!-----------
! Profiles
!-----------

! Exp profiles
    jid = 0    
! Get sub-dictionaries dimensions
    write(nunit, '(A/)') '"profiles_x": {'
    do j=1, n_profx-1
        call write_array((/NA1/), profiles_x(1:NA1, j), profxPtr)
    enddo
    call write_array((/NA1/), profiles_x(1:NA1, n_profx), profxPtr, dict_end=.true.)
    write(nunit, '(A/)') '},' ! End of "profiles_x" dictionary

! CAR profiles
    jid = 0    
    write(nunit, '(A/)') '"profiles": {'
    do j=1, n_prof-1
        call write_array((/NA1/), profiles(1:NA1, j), profPtr)
    enddo
    call write_array((/NA1/), profiles(1:NA1, n_prof), profPtr, dict_end=.true.) ! no comma
    write(nunit, '(A/)') '},' ! End of "profiles" dictionary

!------------
! Equilibrium
!------------

! Scalars

    equil_traces(1) = equil_now%global_param%toroid_field%b0
    equil_traces(2) = equil_now%global_param%betpol
    equil_traces(3) = equil_now%global_param%i_plasma
    equil_traces(4) = equil_now%global_param%li3
    equil_traces(5) = equil_now%global_param%psibound
    equil_traces(6) = equil_now%global_param%psiaxis
    equil_traces(7) = equil_now%global_param%toroid_field%r0
    equil_traces(8) = equil_now%global_param%Vloop
    call write_scalar_block(equil_sigPtr, equil_traces(1: 8), label='equil_signals')

! 1 d profiles
    jid = 0
    if (debug > 0) write(*, *) 'Writing equilibrium entries to json file'
    write(nunit, '(A/)') '"equil_profiles": {'
    call write_array((/nrho_surf/), equil_now%profiles_1d%areat  , equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%bdb0   , equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%bmaxt  , equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%bmint  , equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%dpsidv , equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%elongation, equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%f_dia  , equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%ffprime, equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%fofb   , equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%g1     , equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%g2     , equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%ggradro, equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%gm1    , equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%gm4    , equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%gm41   , equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%gm5    , equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%perim  , equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%phi    , equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%pprime , equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%pressure, equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%psi    , equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%q      , equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%r_inboard , equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%r_outboard, equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%rho_tor_norm, equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%shif   , equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%surface, equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%volume , equil_profPtr, dict_end=.true.)
    write(nunit, '(A/)') '},' ! End of "equil_profiles" dictionary

    if (debug > 0) write(*, *) 'Writing mag.surf quantities'
    jid = 0
    write(nunit, '(A/)') '"equil_coord": {'
    call write_array((/nrho_surf, nthe_surf/), equil_now%coord_sys%position%r, equil_coordPtr)
    call write_array((/nrho_surf, nthe_surf/), equil_now%coord_sys%position%rmin, equil_coordPtr)
    call write_array((/nrho_surf, nthe_surf/), equil_now%coord_sys%position%psirz, equil_coordPtr)
    call write_array((/nthe_surf/), equil_now%coord_sys%position%teta2d, equil_coordPtr)
    call write_array((/nrho_surf, nthe_surf/), equil_now%coord_sys%position%z, equil_coordPtr, dict_end=.true.)
    write(nunit, '(A/)') '},' ! End of "equil_coord" dictionary

    jid = 0
    if (debug > 0) write(*, *) 'Writing equil cartesian quantities'
    write(nunit, '(A/)') '"equil_rect": {'
    call write_array((/nR, nZ/), equil_now%eqgeometry%rectgrid%psirz2d, equil_rectPtr)
    call write_array((/nR, nZ/), equil_now%eqgeometry%rectgrid%fdia2d , equil_rectPtr)
    call write_array((/nR/), equil_now%eqgeometry%rectgrid%r2d, equil_rectPtr)
    call write_array((/nZ/), equil_now%eqgeometry%rectgrid%z2d, equil_rectPtr, dict_end=.true.) ! No comma
    write(nunit, '(A)') '}' ! End of "equil_rect" dictionary

!-----------
! Close json
!-----------

    write(nunit, '(A)') '}' ! End of "equil" dictionary
    close(nunit)

    j_call = j_call + 1

    write(*, '(A)') '   Written file ' // TRIM(json_out)

    return
    end subroutine write_json

!---------------------------------------------------------------
    subroutine write_scalar_block(json_in, scalar_list, label, indent)

    double precision, intent(in), dimension(*) :: scalar_list
    character(len=*), intent(in), optional :: label
    type(json_value), intent(in), pointer :: json_in
    integer, intent(in), optional :: indent

    logical :: found
    integer :: nvars, j, n_indent
    character(KIND=JSON_CK, len=:), allocatable :: sname
    type(json_value), pointer :: dictPointer
    type(json_core) :: jCore
    character(len=40) :: space

101 format(A, '"', A, '": ', es16.8e3, ',')
102 format(A, '"', A, '": ', es16.8e3, '')
201 format(   '"', A, '": {')
202 format(A, '"', A, '": {')

    space = ' '
    if (present(indent)) then
        n_indent = indent
    else
        n_indent = 0
    endif
    call jCore%info(json_in, n_children=nvars)

    if (present(label)) then
        if (n_indent == 0) then
            write(nunit, 201) TRIM(label)
        else
            write(nunit, 202) space(1:n_indent), TRIM(label)
        endif
    endif

    do j=1, nvars
        call jCore%get_child(json_in, j, dictPointer, found)
        call jCore%info(dictPointer, name=sname)
        if (j == nvars) then
            write(nunit, 102) space(1:n_indent+4), sname, scalar_list(j)
        else
            write(nunit, 101) space(1:n_indent+4), sname, scalar_list(j)
        endif
    enddo
 
    if (present(label)) then
        if (n_indent == 0) then
            write(nunit, '(A)') '},'
        else
            write(nunit, '(A, A)') space(1: n_indent), '},'
        endif
    endif

    return
    end subroutine write_scalar_block

!---------------------------------------------------------------------
    subroutine prettyFloat(n_u, arr_in, fmt, n_columns)

    integer, intent(in) :: n_u, n_columns
    character(len=*), intent(in) :: fmt
    double precision, intent(in), dimension(:) :: arr_in

    integer :: ndim, i
    character(len=120) :: fmt1, fmt2

    ndim = SIZE(arr_in)
    fmt1 = '(' // TRIM(fmt) // ')'
    fmt2 = '(' // TRIM(fmt) // ', ", ")'
    do i=1, ndim-1
        write(n_u, TRIM(fmt2), advance='no') arr_in(i)
        if (MODULO(i, n_columns) == 0) write(nunit, *)
    enddo
    write(n_u, TRIM(fmt1), advance='no') arr_in(ndim)

    return
    end subroutine prettyFloat

!---------------------------------------------------------------------
    subroutine write_array(dims, arr_in, json_in, dict_end)

    integer, intent(in), dimension(:) :: dims
    double precision, intent(in), dimension(*) :: arr_in
    type(json_value), intent(in), pointer :: json_in
    logical, intent(in), optional :: dict_end

    logical :: found
    integer :: i, j, ij, n_columns, n_remain
    double precision :: abs_val
    double precision, dimension(:), allocatable :: array
    double precision, dimension(:, :), allocatable :: array2
    character(KIND=JSON_CK, len=:), allocatable :: sname
    character(len=120) :: fmt
    type(json_value), pointer :: dictPointer
    type(json_core) :: jCore

102 format('    "', A, '": [')

    jid = jid + 1
    call jCore%get_child(json_in, jid, dictPointer, found)
    call jCore%info(dictPointer, name=sname)
    write(nunit, 102, advance="no") sname

    if (SIZE(dims) == 1) then
        allocate(array(dims(1)))
        do i=1, dims(1)
            abs_val = ABS(arr_in(i))
            if (abs_val < 1e-20 .or. abs_val > 1e20) then
                array(i) = 0.
            else
                array(i) = arr_in(i)
            endif
        enddo
        write(nunit, *)
        if (MAXVAL(ABS(array)) > 0.) then
            n_columns = 6
            fmt = 'es16.8e3'
        else
            n_columns = 12
            fmt = 'f3.1'
        endif
        call prettyFloat(nunit, array, fmt, n_columns)

    else if (SIZE(dims) == 2) then
        allocate(array2(dims(1), dims(2)))
        do i=1, dims(1)
            do j=1, dims(2)
                ij = dims(1)*(j-1) + i
                abs_val = ABS(arr_in(ij))
                if (abs_val < 1e-20 .or. abs_val > 1e20) then
                    array2(i, j) = 0.
                else
                    array2(i, j) = arr_in(ij)
                endif
            enddo
        enddo
        if (MAXVAL(ABS(array2)) > 0.) then
            n_columns = 6
            fmt = 'es16.8e3'
        else
            n_columns = 12
            fmt = 'f3.1'
        endif
        do i=1, dims(1)
            write(nunit, '(A)') '['
            call prettyFloat(nunit, array2(i, :), fmt, n_columns)
            if (i == dims(1)) then
                write(nunit, '(A)') ']'
            else
                write(nunit, '(A)') '],'
            endif
        enddo
    endif

    if (present(dict_end)) then
        write(nunit, '(A)') ']' ! No comma after last array entry
    else
        write(nunit, '(A)') '],'
    endif

    return
    end subroutine write_array

end module json_rw
