module char_manip

implicit none

character(len=1), parameter :: null_ch=char(0), beep_ch=char(7), tab_ch=char(9), esc_ch=char(13), backslash=char(92)

contains

!----------------------------
    function clean_string(str_in)
! Replaces tabs and <cr> with blanks, adjusts string to left

    character(len=*), intent(in) :: str_in
    character(len=len(str_in)) :: clean_string, str_out

    str_out = str_in
    str_out = replace_char(str_out, tab_ch, ' ')  ! Replace tabs with blanks
    str_out = replace_char(str_out, char(11), ' ')  ! Replace tabs with blanks
    str_out = replace_char(str_out, esc_ch, ' ')  ! Replace <cr> with blanks
    clean_string = ADJUSTL(str_out)

    end function clean_string

!----------------------------
    function replace_char(str_in, char_in, char_out)
! splitstring splits a string to an array of
! substrings based on a selected delimiter
! note any facing space/blank in substrings will be removed

    character(len=*), intent(in) :: str_in
    character(len=1), intent(in) :: char_in, char_out
    character(len=len(str_in)) :: replace_char

    integer :: j
    character(len=1) :: ch

    do j=1, LEN(str_in)
        ch = str_in(j: j)
        if (ch == char_in) ch = char_out
        replace_char(j: j) = ch
    enddo

    end function replace_char

!----------------------------
    function to_upper(str_in)

    integer, parameter :: ch_shift = ichar('A') - ichar('a')

    character(len=*), intent(in) :: str_in
    character(len=len(str_in)) :: to_upper
    integer :: i
    character :: ch

    do i=1, LEN(str_in)
        ch = str_in(i: i)
        if (ch >= 'a' .and. ch <= 'z' ) then
            ch = char(ichar(ch) + ch_shift)
        endif
        to_upper(i: i) = ch
    enddo
    end function to_upper

!----------------------------
    function to_lower(str_in)

    integer, parameter :: ch_shift = ichar('A') - ichar('a')

    character(len=*), intent(in) :: str_in
    character(len=len(str_in)) :: to_lower
    integer :: i
    character :: ch

    do i = 1, LEN(str_in)
        ch = str_in(i: i)
        if (ch >= 'A' .and. ch <= 'Z' ) then
            ch = char(ichar(ch) - ch_shift)
        endif
        to_lower(i: i) = ch
    enddo
    end function to_lower

!----------------------------
    integer function first_non_blank(str_in)

    character(len=*), intent(in) :: str_in

    integer :: j

    first_non_blank = 0
    do j=1, LEN_TRIM(str_in)
        if (str_in(j: j) == null_ch .or. str_in(j: j) == esc_ch) return
        if (str_in(j: j) /= ' '    .and. str_in(j: j) /= tab_ch) then
            first_non_blank = j !9<->tab
            return
        endif
    enddo
    end function first_non_blank

!----------------------------
    integer function len_trim_tab(str_in)

    character(len=*), intent(in) :: str_in

    integer :: j

    len_trim_tab = 0
    do j=1, LEN_TRIM(str_in)
        if (str_in(j: j) == null_ch .or. str_in(j: j) == esc_ch) return
        if (str_in(j: j) /= ' '    .and. str_in(j: j) /= tab_ch) len_trim_tab = j !9<->tab
    enddo
    end function len_trim_tab

!----------------------------
    integer function str_in_list(str_in, str_arr, len_str)

    integer, intent(in), optional :: len_str
    character(len=*), intent(in) :: str_in, str_arr(:)

    integer :: j, n_arr, len

    if (PRESENT(len_str)) then
        len = len_str
    else
        len = 6
    endif
    n_arr = SIZE(str_arr, 1)
    str_in_list = 0
    do j=1, n_arr
        if (str_in(1: len) == str_arr(j)(1: len)) then
            str_in_list = j
            EXIT
        endif
    enddo

    end function str_in_list

!-----------------------------
! split a string into 2 either side of a delimiter (1st occurrence)
    subroutine split_string(str_in, delim, string1, string2)

    character(len=*), intent(in)  :: str_in
    character(1), intent(in) :: delim
    character(len=len(str_in)), intent(out) :: string1, string2

    integer :: jpos
    character(len=len(str_in)) :: strtmp

    strtmp = clean_string(str_in)
    string1 = repeat(' ', len(str_in))
    string2 = repeat(' ', len(str_in))
    jpos = index(strtmp, delim)
    if (jpos == 0) then
        return
    endif

    string1 = ADJUSTL(strtmp(1: jpos-1))
    string2 = ADJUSTL(strtmp(jpos+1:))

    end subroutine split_string


end module char_manip
