module json_vars

use json_module, only : json_file, json_core, json_value, json_ck

implicit none

integer :: n_var, n_varx, n_const, n_intern, n_internInt, n_intern2, &
    n_prof, n_profx, n_equil_sig, n_equil_prof, n_equil_rect, n_equil_coord
character(len=6), allocatable, dimension(:) :: varNames, &
    varxNames, constNames, internNames, internIntNames, intern2Names, &
    profNames, profxNames, equil_sigNames, equil_profNames, &
    equil_rectNames, equil_coordNames
type(json_file) :: astra_vars
type(json_value), pointer :: varPtr, varxPtr, constPtr, &
     internPtr, internIntPtr, intern2Ptr, profPtr, profxPtr, &
     equil_sigPtr, equil_profPtr, equil_rectPtr, equil_coordPtr

contains

!---------------------------------------------------------------------
    subroutine get_subdict(fjson_in, label, nvars, names_out, jsonOut)

    type(json_file), intent(inout) :: fjson_in
    character(len=*), intent(in) :: label
    integer, intent(out) :: nvars  
    character(len=6), intent(out), allocatable, dimension(:) :: names_out
    type(json_value), pointer, intent(out) :: jsonOut

    logical :: found
    integer :: j
    character(KIND=JSON_CK, len=:), allocatable :: sname
    type(json_value), pointer :: dictPointer
    type(json_core) :: jCore

    call fjson_in%info(label, n_children=nvars)
    allocate(names_out(nvars))
    call fjson_in%get(label, jsonOut, found)
    do j=1, nvars
        call jCore%get_child(jsonOut, j, dictPointer, found)
        call jCore%info(dictPointer, name=sname)
        names_out(j) = sname
    enddo

    end subroutine get_subdict
  
!---------------------------------------------------------------------
    subroutine read_metadata

    logical :: status_ok
    character(len=:), allocatable :: error_msg
    character(len=240) :: file_in
    file_in = 'astra_variables.json'
    call astra_vars%initialize()
    call astra_vars%load(filename=trim(file_in))
    if (astra_vars%failed()) then
        write(*, *) 'Error: '
        call astra_vars%check_for_errors(status_ok, error_msg)    
        write(*, *) 'Error: ', error_msg
    endif

! Get sub-dictionaries objects & lists
    call get_subdict(astra_vars, 'variables'     , n_var        ,         varNames,         varPtr)
    call get_subdict(astra_vars, 'variables_x'   , n_varx       ,        varxNames,        varxPtr)
    call get_subdict(astra_vars, 'constants'     , n_const      ,       constNames,       constPtr)
    call get_subdict(astra_vars, 'internal'      , n_intern     ,      internNames,      internPtr)
    call get_subdict(astra_vars, 'internInt'     , n_internInt  ,   internIntNames,   internIntPtr)
    call get_subdict(astra_vars, 'intern2'       , n_intern2    ,     intern2Names,     intern2Ptr)
    call get_subdict(astra_vars, 'profiles'      , n_prof       ,        profNames,        profPtr)
    call get_subdict(astra_vars, 'profiles_x'    , n_profx      ,       profxNames,       profxPtr)
    call get_subdict(astra_vars, 'equil_signals' , n_equil_sig  ,   equil_sigNames,   equil_sigPtr)
    call get_subdict(astra_vars, 'equil_profiles', n_equil_prof ,  equil_profNames,  equil_profPtr)
    call get_subdict(astra_vars, 'equil_rect'    , n_equil_rect ,  equil_rectNames,  equil_rectPtr)
    call get_subdict(astra_vars, 'equil_coord'   , n_equil_coord, equil_coordNames, equil_coordPtr)

    end subroutine read_metadata

end module json_vars
