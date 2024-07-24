module time_functions

implicit none

integer, parameter :: nloc=2200
integer :: f_id = 0
  
contains

!---------------------------------------------------------------------
    subroutine function_id(str_in)

    use const_inc, only: TIME, TSTART

    character(len=*), intent(in) :: str_in
    double precision :: time_loc

    save time_loc

    if (time_loc /= TIME .or. TIME == TSTART) then
        f_id = 1
    else
        f_id = f_id + 1
    endif
    if (f_id > NLOC) then
        write(*, *) 'Too many time functions calls', f_id
        write(*, *) 'Calling from' // trim(str_in)
        return
    endif

    time_loc = TIME

!    write(*, *) TRIM(str_in), TIME, f_id

    return
    end subroutine function_id

end module time_functions
