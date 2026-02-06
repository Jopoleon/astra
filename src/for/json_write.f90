module json_write

use json_module, only: json_core, json_value, json_ck

implicit none

integer, parameter :: nunit=25, nunit_x=35
integer :: jid

contains

!---------------------------------------------------------------------
    subroutine write_jsonx

    use io_mod, only: awd, exp_file
    use exp_data, only: raw_scalars, raw_profiles, raw_boundary

    integer :: ios, i, ndim
    character(len=180) :: jsonx_out
    type(json_core) :: jCore

    write(jsonx_out, '(4A)') TRIM(awd), '/ncdf_out/', TRIM(exp_file), '_x.json'

    open(nunit_x, file=TRIM(jsonx_out), iostat=ios)

    write(nunit_x, '(A/)') '{'

! Scalars

    write(nunit_x, '(A/)') '   "scalars": {'

    write(nunit_x, '(A)') '        "timeStream": {"unit": "s", "data": ['
    call prettyFloatArray(nunit_x, raw_scalars%nt_all, raw_scalars%time)

    write(nunit_x, '(A)') '        "dataStream": {"data": ['
    call prettyFloatArray(nunit_x, raw_scalars%nt_all, raw_scalars%data)

    write(nunit_x, '(A/)', advance='no') '        "labels": {"data": ['
    do i=1, raw_scalars%nt_all-1
        write(nunit_x, '(3A)', advance='no') '"', TRIM(raw_scalars%label(i)), '", '
        if (MODULO(i, 10) == 0) write(nunit_x, '(A)') '        '
    enddo
    write(nunit_x, '(3A)') '"', TRIM(raw_scalars%label(raw_scalars%nt_all)), '"]'
    write(nunit_x, '(A/)') '        }'

    write(nunit_x, '(A/)') '    },'

! Profiles

    write(nunit_x, '(A)') '   "profiles": {'
    write(nunit_x, '(A/)') '    },'

! Boundary

    write(nunit_x, '(A, i3, A, i3, A/)') '   "boundary": { "nt": ', raw_boundary%nt, &
         ' "n_theta": ', raw_boundary%n_theta, ','

    write(nunit_x, '(A)') '        "time": {"unit": "s", "data": ['
    call prettyFloatArray(nunit_x, raw_boundary%nt, raw_boundary%time)

    ndim = raw_boundary%nt * raw_boundary%n_theta

    write(nunit_x, '(A)') '        "R": {"unit": "m", "data": ['
    call prettyFloatArray(nunit_x, ndim, raw_boundary%R)

    write(nunit_x, '(A)') '        "Z": {"unit": "m", "data": ['
    call prettyFloatArray(nunit_x, ndim, raw_boundary%Z, dict_end=.true.)

    write(nunit_x, '(A/)') '    }'

! Closing

    write(nunit_x, '(A)') '}'
    close(nunit_x)

    write(*, '(A)') '   Written file ' // TRIM(jsonx_out)

11  format(5('"', A, '",'))

    return
    end subroutine write_jsonx

!---------------------------------------------------------------------
    subroutine prettyFloatArray(n_u, ndim, arr, dict_end)

    integer, intent(in) :: n_u, ndim
    double precision, intent(in) :: arr(*)
    logical, intent(in), optional :: dict_end

    integer :: i

    do i=1, ndim-1
        write(n_u, '(es16.8e3, A)', advance='no') arr(i), ','
        if (MODULO(i, 6) == 0) write(n_u, '(A)') ''
    enddo
    write(n_u, '(es16.8e3, A)') arr(ndim), ']'
    if (present(dict_end)) then
        write(n_u, '(A/)') '        }'
    else
        write(n_u, '(A/)') '        },'
    endif
 
    return
    end subroutine prettyFloatArray

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
    
    write(nunit, '(A/)') '"astra": {'

! Scalars
    call write_scalar_block(varPtr, varValues)
    call write_scalar_block(varxPtr, varxValues)
    call write_scalar_block(constPtr, constValues)
    call write_scalar_block(internPtr, internValues)
    call write_scalar_block(intern2Ptr, intern2Values)

!-----------
! Profiles
!-----------

! Exp profiles
    jid = 0    
! Get sub-dictionaries dimensions
    do j=1, n_profx
        call write_array((/NA1/), profiles_x(1:NA1, j), profxPtr)
    enddo

! CAR profiles
    jid = 0    
    do j=1, n_prof-1
        call write_array((/NA1/), profiles(1:NA1, j), profPtr)
    enddo
    call write_array((/NA1/), profiles(1:NA1, n_prof), profPtr, dict_end=.true.) ! no comma

    write(nunit, '(A/)') '},' ! End of "astra" dictionary

!------------
! Equilibrium
!------------

    write(nunit, '(A/)') '"equil": {'

! Scalars

    equil_traces(1) = equil_now%global_param%toroid_field%b0
    equil_traces(2) = equil_now%global_param%betpol
    equil_traces(3) = equil_now%global_param%i_plasma
    equil_traces(4) = equil_now%global_param%li3
    equil_traces(5) = equil_now%global_param%psibound
    equil_traces(6) = equil_now%global_param%psiaxis
    equil_traces(7) = equil_now%global_param%toroid_field%r0
    equil_traces(8) = equil_now%global_param%Vloop
    call write_scalar_block(equil_sigPtr, equil_traces(1: 8))

! 1 d profiles
    jid = 0
    if (debug > 0) write(*, *) 'Writing equilibrium entries to json file'
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
    call write_array((/nrho_surf/), equil_now%profiles_1d%volume , equil_profPtr)

    if (debug > 0) write(*, *) 'Writing mag.surf quantities'
    jid = 0
    call write_array((/nrho_surf, nthe_surf/), equil_now%coord_sys%position%r, equil_coordPtr)
    call write_array((/nrho_surf, nthe_surf/), equil_now%coord_sys%position%rmin, equil_coordPtr)
    call write_array((/nrho_surf, nthe_surf/), equil_now%coord_sys%position%psirz, equil_coordPtr)
    call write_array((/nthe_surf/), equil_now%coord_sys%position%teta2d, equil_coordPtr)
    call write_array((/nrho_surf, nthe_surf/), equil_now%coord_sys%position%z, equil_coordPtr)

    jid = 0
    if (debug > 0) write(*, *) 'Writing equil cartesian quantities'
    call write_array((/nR, nZ/), equil_now%eqgeometry%rectgrid%psirz2d, equil_rectPtr)
    call write_array((/nR, nZ/), equil_now%eqgeometry%rectgrid%fdia2d , equil_rectPtr)
    call write_array((/nR/), equil_now%eqgeometry%rectgrid%r2d, equil_rectPtr)
    call write_array((/nZ/), equil_now%eqgeometry%rectgrid%z2d, equil_rectPtr, dict_end=.true.) ! No comma
    write(nunit, '(A)') '}' ! End of "equil" dictionary

!-----------
! Close json
!-----------

    write(nunit, '(A)') '}'
    close(nunit)

    j_call = j_call + 1

    write(*, '(A)') '   Written file ' // TRIM(json_out)

    return
    end subroutine write_json

!---------------------------------------------------------------
    subroutine write_scalar_block(json_in, scalar_list)

    double precision, intent(in), dimension(*) :: scalar_list
    type(json_value), intent(in), pointer :: json_in

    logical :: found
    integer :: nvars, j
    character(KIND=JSON_CK, len=:), allocatable :: sunit, sdesc, sname
    type(json_value), pointer :: dictPointer
    type(json_core) :: jCore

101 format('    "', A, '": {"units": "', A, '", "long_name": "', A, '", "data": ', es16.8e3, '},')

    call jCore%info(json_in, n_children=nvars)
    do j=1, nvars
        call jCore%get_child(json_in, j, dictPointer, found)
        call jCore%info(dictPointer, name=sname)
        call jCore%get(dictPointer, 'units', sunit, found)
        call jCore%get(dictPointer, 'desc' , sdesc, found)
        write(nunit, 101) sname, sunit, sdesc, scalar_list(j)
    enddo
    
    return
    end subroutine write_scalar_block

!---------------------------------------------------------------
    subroutine ndim_string(dims, ndim, sdim)

    integer, intent(in), dimension(:) :: dims
    integer, intent(out) :: ndim  
    character(len=120), intent(out) :: sdim

    if (SIZE(dims) == 1) then
        write(sdim, '(A,i0,A)') '[', dims(1), ']'
        ndim = dims(1)
    else if (SIZE(dims) == 2) then
        write(sdim, '(A,i0,A,i0,A)') '[', dims(1), ', ', dims(2), ']'
        ndim = dims(1)*dims(2)
    endif
     
    return
    end subroutine ndim_string

!---------------------------------------------------------------
    subroutine write_array(dims, arr_in, json_in, dict_end)

    integer, intent(in), dimension(:) :: dims
    double precision, intent(in), dimension(*) :: arr_in
    type(json_value), intent(in), pointer :: json_in
    logical, intent(in), optional :: dict_end

    logical :: found
    integer :: i, j, ndim, ij
    double precision :: abs_val
    double precision, dimension(:), allocatable :: array
    double precision, dimension(:, :), allocatable :: array2
    character(KIND=JSON_CK, len=:), allocatable :: sunit, sdesc, sname
    character(len=120) :: sdim
    type(json_value), pointer :: dictPointer
    type(json_core) :: jCore

102 format('    "', A, '": {"dims": ', A, ', "units": "', A, '", "long_name": "', A, '", "data": [')

    call ndim_string(dims, ndim, sdim)
    jid = jid + 1
    call jCore%get_child(json_in, jid, dictPointer, found)
    call jCore%get(dictPointer, 'units', sunit, found)
    call jCore%get(dictPointer, 'desc', sdesc, found)
    call jCore%info(dictPointer, name=sname)
    write(nunit, 102, advance="no") sname, TRIM(sdim), sunit, sdesc

    if (SIZE(dims) == 1) then
        allocate(array(ndim))
        do i=1, ndim
            abs_val = ABS(arr_in(i))
            if (abs_val < 1e-20 .or. abs_val > 1e20) then
                array(i) = 0.
            else
                array(i) = arr_in(i)
            endif
        enddo
        if (MAXVAL(ABS(array)) > 0.) then
            write(nunit, *)
            write(nunit, '(5(es16.8e3, ","))') (array(i), i=1, ndim-1)
            write(nunit, '(es16.8e3)', advance='no') array(ndim) ! No comma after last array entry
        endif
 
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
            do i=1, dims(1)
                write(nunit, '(A)') '['
                write(nunit, '(5(es16.8e3, ","))') (array2(i, j), j=1, dims(2)-1)
                write(nunit, '(es16.8e3)', advance='no') array2(i, dims(2)) ! No comma after last array entry
                if (i == dims(1)) then
                    write(nunit, '(A)') ']'
                else
                    write(nunit, '(A)') '],'
                endif
            enddo
        endif
    endif

    write(nunit, '(A)') ']}' ! No comma after last array entry
    if (.not. present(dict_end)) write(nunit, '(A)') ','

    return
    end subroutine write_array

end module json_write
