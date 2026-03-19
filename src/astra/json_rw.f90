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
    write(*, '(2A)') 'Opening json file ', TRIM(json_in)

    if (fjson_out%failed()) then
        print*, 'Error: '
        call fjson_out%check_for_errors(status_ok, error_msg)    
        write(*, *) 'Error: ', error_msg
    endif

    end subroutine json_load

!---------------------------------------------------------------------
    subroutine read_scalars_int(fjson, label, values)

    type(json_file), intent(inout) :: fjson
    character(len=*), intent(in) :: label
    integer, intent(out), allocatable, dimension(:) :: values

    logical :: found
    integer :: j, nvars
    integer :: val
    type(json_value), pointer :: jsonBlock, dictPointer
    type(json_core) :: jCore    

    write(*, '(2A)') 'Reading json block ', TRIM(label)

    call fjson%info(TRIM(label), n_children=nvars)
    allocate(values(nvars))

    call fjson%get(TRIM(label), jsonBlock, found)
    do j=1, nvars
        call jCore%get_child(jsonBlock, j, dictPointer, found)
        call jCore%get(dictPointer, val)
        values(j) = val
    enddo

    end subroutine read_scalars_int

!---------------------------------------------------------------------
    subroutine read_scalars_float(fjson, label, values)

    type(json_file), intent(inout) :: fjson
    character(len=*), intent(in) :: label
    double precision, intent(out), allocatable, dimension(:) :: values

    logical :: found
    integer :: j, nvars
    double precision :: val
    type(json_value), pointer :: jsonBlock, dictPointer
    type(json_core) :: jCore    

    write(*, '(2A)') 'Reading json block ', TRIM(label)

    call fjson%info(TRIM(label), n_children=nvars)
    allocate(values(nvars))

    call fjson%get(TRIM(label), jsonBlock, found)
    do j=1, nvars
        call jCore%get_child(jsonBlock, j, dictPointer, found)
        call jCore%get(dictPointer, val)
        values(j) = val
    enddo

    end subroutine read_scalars_float

!---------------------------------------------------------------------
    subroutine read_array_block(fjson, label, profs)

    type(json_file), intent(inout) :: fjson
    character(len=*), intent(in) :: label
    double precision, intent(out), allocatable :: profs(:, :)

    logical :: found
    integer :: j, k, nvars, npts
    double precision :: val

    type(json_value), pointer :: jsonBlock, arrPointer, valPointer
    type(json_core) :: jCore

    write(*, '(2A)') 'Reading json block ', trim(label)

! get "profiles"
    call fjson%get(trim(label), jsonBlock, found)

! number of arrays (arr1, arr2, ...)
    call jCore%info(jsonBlock, n_children=nvars)

! get first array to get length
    call jCore%get_child(jsonBlock, 1, arrPointer, found)
    call jCore%info(arrPointer, n_children=npts)

    allocate(profs(npts, nvars))

! read values
    do j=1, nvars
        call jCore%get_child(jsonBlock, j, arrPointer, found)
        do k=1, npts
            call jCore%get_child(arrPointer, k, valPointer, found)
            call jCore%get(valPointer, val)
            profs(k, j) = val
        enddo
    enddo

    end subroutine read_array_block

!---------------------------------------------------------------------
    subroutine read_array_1d(fjson, parent, name, arr1d)

    type(json_file), intent(inout) :: fjson

    character(len=*), intent(in) :: parent
    character(len=*), intent(in) :: name
    double precision, allocatable, intent(out) :: arr1d(:)

    logical :: found
    integer :: n, i
    double precision :: val
    type(json_value), pointer :: parentPtr, arrPtr, valPtr
    type(json_core) :: jCore

! get parent object
    call fjson%get(trim(parent), parentPtr, found)
    if (.not. found) then
        print *, "Parent not found:", parent
        stop
    endif

! get array inside parent
    call jCore%get(parentPtr, trim(name), arrPtr, found)
    if (.not. found) then
        print *, "Array not found:", name
        stop
    endif

! get length
    call jCore%info(arrPtr, n_children=n)
    allocate(arr1d(n))

! read values
    do i=1, n
        call jCore%get_child(arrPtr, i, valPtr, found)
        call jCore%get(valPtr, val)
        arr1d(i) = val
    enddo

    end subroutine read_array_1d
    
!---------------------------------------------------------------------
    subroutine read_array_2d(fjson, parent, name, arr2d)

    type(json_file), intent(inout) :: fjson
    character(len=*), intent(in) :: parent
    character(len=*), intent(in) :: name

    double precision, allocatable, intent(out) :: arr2d(:, :)

    type(json_core) :: jCore
    type(json_value), pointer :: parentPtr, arrPtr, rowPtr, valPtr

    logical :: found
    integer :: nrows, ncols, i, j
    double precision :: val

! -------------------------
! get parent object
! -------------------------
    call fjson%get(trim(parent), parentPtr, found)
    if (.not. found) then
        print *, "Parent not found:", parent
        stop
    endif

! -------------------------
! get array inside parent
! -------------------------
    call jCore%get(parentPtr, trim(name), arrPtr, found)
    if (.not. found) then
        print *, "Array not found:", name
        stop
    endif

! -------------------------
! number of rows
! -------------------------
    call jCore%info(arrPtr, n_children=nrows)
    if (nrows <= 0) then
        print *, "Empty array:", name
        stop
    endif

! -------------------------
! get first row to get ncols
! -------------------------
    call jCore%get_child(arrPtr, 1, rowPtr, found)
    call jCore%info(rowPtr, n_children=ncols)

    allocate(arr2d(nrows, ncols))

! -------------------------
! read values
! -------------------------
    do i=1, nrows
        call jCore%get_child(arrPtr, i, rowPtr, found)
        do j=1, ncols
            call jCore%get_child(rowPtr, j, valPtr, found)
            call jCore%get(valPtr, val)
            arr2d(i, j) = val
        enddo
    enddo

end subroutine read_array_2d

!---------------------------------------------------------------------
    subroutine read_equil(fjson)

    use parameters_a2equil, only: equil_now

    type(json_file), intent(inout) :: fjson

    integer :: nrho_eq, nthe_eq, nr_eq, nz_eq
    double precision, allocatable, dimension(:) :: equil_traces, &
        r2d, z2d, theta2d
    double precision, allocatable, dimension(:, :) :: equil_profiles, &
        psirz2d, fdia2d, r, z, rmin, psirz

! equil_scalars
    call read_scalars_float(fjson, "equil_signals", equil_traces)

    equil_now%global_param%toroid_field%b0 = equil_traces(1)
    equil_now%global_param%betpol          = equil_traces(2)
    equil_now%global_param%i_plasma        = equil_traces(3)
    equil_now%global_param%li3             = equil_traces(4)
    equil_now%global_param%psibound        = equil_traces(5)
    equil_now%global_param%psiaxis         = equil_traces(6)
    equil_now%global_param%toroid_field%r0 = equil_traces(7)
    equil_now%global_param%Vloop           = equil_traces(8)

! equil_profiles
    call read_array_block(fjson, "equil_profiles", equil_profiles)
    nrho_eq = SIZE(equil_profiles, 1)

    allocate(equil_now%profiles_1d%areat     (nrho_eq))
    allocate(equil_now%profiles_1d%bdb0      (nrho_eq))
    allocate(equil_now%profiles_1d%bmaxt     (nrho_eq))
    allocate(equil_now%profiles_1d%bmint     (nrho_eq))
    allocate(equil_now%profiles_1d%dpsidv    (nrho_eq))
    allocate(equil_now%profiles_1d%elongation(nrho_eq))
    allocate(equil_now%profiles_1d%f_dia     (nrho_eq))
    allocate(equil_now%profiles_1d%ffprime   (nrho_eq))
    allocate(equil_now%profiles_1d%fofb      (nrho_eq))
    allocate(equil_now%profiles_1d%g1        (nrho_eq))
    allocate(equil_now%profiles_1d%g2        (nrho_eq))
    allocate(equil_now%profiles_1d%ggradro   (nrho_eq))
    allocate(equil_now%profiles_1d%gm1       (nrho_eq))
    allocate(equil_now%profiles_1d%gm4       (nrho_eq))
    allocate(equil_now%profiles_1d%gm41      (nrho_eq))
    allocate(equil_now%profiles_1d%gm5       (nrho_eq))
    allocate(equil_now%profiles_1d%perim     (nrho_eq))
    allocate(equil_now%profiles_1d%phi       (nrho_eq))
    allocate(equil_now%profiles_1d%pprime    (nrho_eq))
    allocate(equil_now%profiles_1d%pressure  (nrho_eq))
    allocate(equil_now%profiles_1d%psi       (nrho_eq))
    allocate(equil_now%profiles_1d%q         (nrho_eq))
    allocate(equil_now%profiles_1d%r_inboard (nrho_eq))
    allocate(equil_now%profiles_1d%r_outboard(nrho_eq))
    allocate(equil_now%profiles_1d%rho_tor_norm(nrho_eq))
    allocate(equil_now%profiles_1d%shif      (nrho_eq))
    allocate(equil_now%profiles_1d%surface   (nrho_eq))
    allocate(equil_now%profiles_1d%volume    (nrho_eq))

    equil_now%profiles_1d%areat      = equil_profiles(:, 1)
    equil_now%profiles_1d%bdb0       = equil_profiles(:, 2)
    equil_now%profiles_1d%bmaxt      = equil_profiles(:, 3)
    equil_now%profiles_1d%bmint      = equil_profiles(:, 4)
    equil_now%profiles_1d%dpsidv     = equil_profiles(:, 5)
    equil_now%profiles_1d%elongation = equil_profiles(:, 6)
    equil_now%profiles_1d%f_dia      = equil_profiles(:, 7)
    equil_now%profiles_1d%ffprime    = equil_profiles(:, 8)
    equil_now%profiles_1d%fofb       = equil_profiles(:, 9)
    equil_now%profiles_1d%g1         = equil_profiles(:, 10)
    equil_now%profiles_1d%g2         = equil_profiles(:, 11)
    equil_now%profiles_1d%ggradro    = equil_profiles(:, 12)
    equil_now%profiles_1d%gm1        = equil_profiles(:, 13)
    equil_now%profiles_1d%gm4        = equil_profiles(:, 14)
    equil_now%profiles_1d%gm41       = equil_profiles(:, 15)
    equil_now%profiles_1d%gm5        = equil_profiles(:, 16)
    equil_now%profiles_1d%perim      = equil_profiles(:, 17)
    equil_now%profiles_1d%phi        = equil_profiles(:, 18)
    equil_now%profiles_1d%pprime     = equil_profiles(:, 19)
    equil_now%profiles_1d%pressure   = equil_profiles(:, 20)
    equil_now%profiles_1d%psi        = equil_profiles(:, 21)
    equil_now%profiles_1d%q          = equil_profiles(:, 22)
    equil_now%profiles_1d%r_inboard  = equil_profiles(:, 23)
    equil_now%profiles_1d%r_outboard = equil_profiles(:, 24)
    equil_now%profiles_1d%rho_tor_norm = equil_profiles(:, 25)
    equil_now%profiles_1d%shif       = equil_profiles(:, 26)
    equil_now%profiles_1d%surface    = equil_profiles(:, 27)
    equil_now%profiles_1d%volume     = equil_profiles(:, 28)

! equil_rect
    call read_array_1d(fjson, "equil_rect", "r2d", r2d)
    call read_array_1d(fjson, "equil_rect", "z2d", z2d)
    call read_array_2d(fjson, "equil_rect", "psirz2d", psirz2d)
    call read_array_2d(fjson, "equil_rect", "fdia2d", fdia2d)

    nr_eq = SIZE(r2d)
    nz_eq = SIZE(z2d)

    allocate(equil_now%eqgeometry%rectgrid%r2d(nr_eq))
    allocate(equil_now%eqgeometry%rectgrid%z2d(nz_eq))
    allocate(equil_now%eqgeometry%rectgrid%psirz2d(nr_eq, nz_eq))
    allocate(equil_now%eqgeometry%rectgrid%fdia2d(nr_eq, nz_eq))

    equil_now%eqgeometry%rectgrid%r2d = r2d
    equil_now%eqgeometry%rectgrid%z2d = z2d
    equil_now%eqgeometry%rectgrid%psirz2d = psirz2d
    equil_now%eqgeometry%rectgrid%fdia2d  = fdia2d

! equil_coord
    call read_array_1d(fjson, "equil_coord", "theta2d", theta2d)
    call read_array_2d(fjson, "equil_coord", "r", r)
    call read_array_2d(fjson, "equil_coord", "z", z)
    call read_array_2d(fjson, "equil_coord", "rmin", rmin)
    call read_array_2d(fjson, "equil_coord", "psirz", psirz)

    nthe_eq = SIZE(theta2d)

    allocate(equil_now%coord_sys%position%theta2d(nthe_eq))
    allocate(equil_now%coord_sys%position%r(nrho_eq, nthe_eq))
    allocate(equil_now%coord_sys%position%z(nrho_eq, nthe_eq))
    allocate(equil_now%coord_sys%position%rmin(nrho_eq, nthe_eq))
    allocate(equil_now%coord_sys%position%psirz(nrho_eq, nthe_eq))

    equil_now%coord_sys%position%theta2d = theta2d
    equil_now%coord_sys%position%r = r
    equil_now%coord_sys%position%z = z
    equil_now%coord_sys%position%rmin  = rmin
    equil_now%coord_sys%position%psirz = psirz

    end subroutine read_equil

!---------------------------------------------------------------------
    subroutine write_ajson

    use parameters_a2equil, only: equil_now
    use scalars, only: NA1, varValues, varxValues, constValues, &
        internValues, internIntValues, intern2Values
    use status, only: profiles, profiles_x
    use read_input, only: awd, exp_file, equ_file, restart
    use debugger, only: debug
    use json_vars, only: equil_sigPtr, equil_profPtr, equil_rectPtr, equil_coordPtr, &
        profPtr, profxPtr, constPtr, internPtr, internIntPtr, intern2Ptr, &
        varPtr, varxPtr, n_prof, n_profx

    integer :: j, ios, j_call=1, j_out, nrho_surf, nthe_surf, nR, nZ
    character(len=180) :: json_out
    double precision, dimension(20) :: equil_traces

    save j_call

    nrho_surf = SIZE(equil_now%profiles_1d%rho_tor_norm)
    nthe_surf = SIZE(equil_now%coord_sys%position%theta2d)
    nR = SIZE(equil_now%eqgeometry%rectgrid%r2d)
    nZ = SIZE(equil_now%eqgeometry%rectgrid%z2d)

    j_out = j_call + restart
    if (debug > 0) write(*, '(A, 3i)') 'Starting a2json', j_out, nrho_surf, nthe_surf

    write(json_out, '(5A, i0, A)') TRIM(awd), '/ncdf_out/', TRIM(exp_file), TRIM(equ_file), '-', j_out, '.json'

    open(nunit, file=TRIM(json_out), iostat=ios)
    write(nunit, '(A/)') '{'

!------
! ASTRA
!------

! Scalars
    call write_scalars_float(varPtr, varValues, label="variables")
    call write_scalars_float(varxPtr, varxValues, label="variables_x")
    call write_scalars_float(constPtr, constValues, label="constants")
    call write_scalars_float(internPtr, internValues, label="internal")
    call write_scalars_int(internIntPtr, internIntValues, label="internInt")
    call write_scalars_float(intern2Ptr, intern2Values, label="intern2")

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
    call write_array((/NA1/), profiles_x(1:NA1, n_profx), profxPtr, last_array=.true.)
    write(nunit, '(A/)') '},' ! End of "profiles_x" dictionary

! CAR profiles
    jid = 0    
    write(nunit, '(A/)') '"profiles": {'
    do j=1, n_prof-1
        call write_array((/NA1/), profiles(1:NA1, j), profPtr)
    enddo
    call write_array((/NA1/), profiles(1:NA1, n_prof), profPtr, last_array=.true.) ! no comma
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
    call write_scalars_float(equil_sigPtr, equil_traces(1: 8), label='equil_signals')

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
    call write_array((/nrho_surf/), equil_now%profiles_1d%volume , equil_profPtr, last_array=.true.)
    write(nunit, '(A/)') '},' ! End of "equil_profiles" dictionary

    if (debug > 0) write(*, *) 'Writing mag.surf quantities'
    jid = 0
    write(nunit, '(A/)') '"equil_coord": {'
    call write_array((/nrho_surf, nthe_surf/), equil_now%coord_sys%position%r, equil_coordPtr)
    call write_array((/nrho_surf, nthe_surf/), equil_now%coord_sys%position%rmin, equil_coordPtr)
    call write_array((/nrho_surf, nthe_surf/), equil_now%coord_sys%position%psirz, equil_coordPtr)
    call write_array((/nthe_surf/), equil_now%coord_sys%position%theta2d, equil_coordPtr)
    call write_array((/nrho_surf, nthe_surf/), equil_now%coord_sys%position%z, equil_coordPtr, last_array=.true.)
    write(nunit, '(A/)') '},' ! End of "equil_coord" dictionary

    jid = 0
    if (debug > 0) write(*, *) 'Writing equil cartesian quantities'
    write(nunit, '(A/)') '"equil_rect": {'
    call write_array((/nR, nZ/), equil_now%eqgeometry%rectgrid%psirz2d, equil_rectPtr)
    call write_array((/nR, nZ/), equil_now%eqgeometry%rectgrid%fdia2d , equil_rectPtr)
    call write_array((/nR/), equil_now%eqgeometry%rectgrid%r2d, equil_rectPtr)
    call write_array((/nZ/), equil_now%eqgeometry%rectgrid%z2d, equil_rectPtr, last_array=.true.) ! No comma
    write(nunit, '(A)') '}' ! End of "equil_rect" dictionary

!-----------
! Close json
!-----------

    write(nunit, '(A)') '}' ! End of "equil" dictionary
    close(nunit)

    j_call = j_call + 1

    write(*, '(A)') '   Written file ' // TRIM(json_out)

    end subroutine write_ajson

!---------------------------------------------------------------
    subroutine write_scalars_int(json_in, scalar_list, label)

    integer, intent(in), dimension(*) :: scalar_list
    character(len=*), intent(in), optional :: label
    type(json_value), intent(in), pointer :: json_in

    logical :: found
    integer :: nvars, j
    character(KIND=JSON_CK, len=:), allocatable :: sname
    type(json_value), pointer :: dictPointer
    type(json_core) :: jCore

101 format('    "', A, '": ', i0, ',')
102 format('    "', A, '": ', i0, '')
201 format(   '"', A, '": {')

    call jCore%info(json_in, n_children=nvars)

    if (present(label)) write(nunit, 201) TRIM(label)

    do j=1, nvars
        call jCore%get_child(json_in, j, dictPointer, found)
        call jCore%info(dictPointer, name=sname)
        if (j < nvars) then
            write(nunit, 101) sname, scalar_list(j)
        else
            write(nunit, 102) sname, scalar_list(j)
        endif
    enddo
 
    if (present(label)) write(nunit, '(A)') '},'

    end subroutine write_scalars_int

!---------------------------------------------------------------
    subroutine write_scalars_float(json_in, scalar_list, label)

    double precision, intent(in), dimension(*) :: scalar_list
    character(len=*), intent(in), optional :: label
    type(json_value), intent(in), pointer :: json_in

    logical :: found
    integer :: nvars, j
    character(KIND=JSON_CK, len=:), allocatable :: sname
    type(json_value), pointer :: dictPointer
    type(json_core) :: jCore

101 format('    "', A, '": ', es16.8e3, ',')
102 format('    "', A, '": ', es16.8e3, '')
201 format(   '"', A, '": {')

    call jCore%info(json_in, n_children=nvars)

    if (present(label)) write(nunit, 201) TRIM(label)

    do j=1, nvars
        call jCore%get_child(json_in, j, dictPointer, found)
        call jCore%info(dictPointer, name=sname)
        if (j < nvars) then
            write(nunit, 101) sname, scalar_list(j)
        else
            write(nunit, 102) sname, scalar_list(j)
        endif
    enddo
 
    if (present(label)) write(nunit, '(A)') '},'

    end subroutine write_scalars_float

!---------------------------------------------------------------------
    subroutine prettyFloatArray(n_u, arr_in, fmt, n_columns)

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

    end subroutine prettyFloatArray

!---------------------------------------------------------------------
    subroutine write_array(dims, arr_in, json_in, last_array)

    integer, intent(in), dimension(:) :: dims
    double precision, intent(in), dimension(*) :: arr_in
    type(json_value), intent(in), pointer :: json_in
    logical, intent(in), optional :: last_array

    logical :: found
    integer :: i, j, ij, n_columns
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
        call prettyFloatArray(nunit, array, fmt, n_columns)

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
            call prettyFloatArray(nunit, array2(i, :), fmt, n_columns)
            if (i == dims(1)) then
                write(nunit, '(A)') ']'
            else
                write(nunit, '(A)') '],'
            endif
        enddo
    endif

    if (present(last_array)) then
        write(nunit, '(A)') ']' ! No comma after last array of the block
    else
        write(nunit, '(A)') '],'
    endif

    end subroutine write_array

!---------------------------------------------------------------------
    subroutine read_ajson(n_restart)

    use read_input, only: awd, exp_file, equ_file
    use scalars, only: NA1, constValues, varValues, varxValues, &
        internValues, internIntValues, intern2Values
    use status, only: profiles, profiles_x
    use json_vars, only: n_intern

    integer, intent(in) :: n_restart

    double precision, dimension(:), allocatable :: internVal
    double precision, dimension(:, :), allocatable :: profs, profs_x
    character(len=132) :: f_json
    type(json_file) :: fjson

    write(f_json, '(5A, i0, A)') TRIM(awd), '/ncdf_out/', TRIM(exp_file), &
        TRIM(equ_file), '-', n_restart, '.json'
    call json_load(f_json, fjson)
    call read_scalars_float(fjson, "constants", constValues)
    call read_scalars_float(fjson, "variables", varValues)
    call read_scalars_float(fjson, "variables_x", varxValues)
    call read_scalars_float(fjson, "internal", internVal)
    call read_scalars_int(fjson, "internInt", internIntValues)
    call read_scalars_float(fjson, "intern2", intern2Values)
    internValues(1: n_intern) = internVal(1: n_intern)
    call read_array_block(fjson, "profiles", profs)
    call read_array_block(fjson, "profiles_x", profs_x)
    profiles(  1:NA1, :) = profs
    profiles_x(1:NA1, :) = profs_x
    call read_equil(fjson)
    call fjson%destroy()
      end subroutine read_ajson
end module json_rw
