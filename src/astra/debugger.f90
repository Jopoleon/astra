module debugger

implicit none

integer :: debug=0
character(len=132) :: last_mark, sec_last_mark

contains

    subroutine markloc(str_in, debug_lev)

    character(len=*), intent(in) :: str_in
    integer, optional, intent(in) :: debug_lev

    integer :: verbose

    if (PRESENT(debug_lev)) then
        verbose = debug_lev
    else
        verbose = debug
    endif

    SELECT CASE(verbose)
    CASE(1)
        write(*, '(/A)') 'Trackback:'
        write(*, '(4X, A)') TRIM(str_in)
    CASE(2)
        write(*, '(/A)') 'Trackback:'
        write(*, '(4X, A)') TRIM(last_mark)
        write(*, '(8X, A)') TRIM(str_in)
    CASE(3)
        write(*, '(/A)') 'Trackback:'
        write(*, '( 4X, A)') TRIM(sec_last_mark)
        write(*, '( 8X, A)') TRIM(last_mark)
        write(*, '(12X, A)') TRIM(str_in)
    END SELECT

    sec_last_mark = TRIM(last_mark)
    last_mark = TRIM(str_in)

    end subroutine markloc

end module debugger
