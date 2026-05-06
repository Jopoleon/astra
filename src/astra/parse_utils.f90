module parse_utils

use char_manip, only: clean_string, split_string, str_in_list, to_upper, to_lower

implicit none

contains

!--------------------------------------------------
! Set variable list from file parsing
    SUBROUTINE assign_val(file_in, narr, arr_in, arr_out, n_dim_next)

    use debugger, only: markloc

    integer, intent(in) :: narr
    character(len=6), dimension(narr), intent(in) :: arr_in
    character(len=40), intent(in) :: file_in
    integer, intent(out) :: n_dim_next
    double precision, dimension(narr), intent(out) :: arr_out

    integer :: ios, jpos, j
    character(len=132) :: STRI, str_nam, str_val

    call markloc('assign_val')

    do j=1, narr
        open(171, FILE=TRIM(file_in), iostat=ios, status='old')
        line_loop: do
            read(171, '(A132)', iostat=ios) STRI
            if (ios < 0) EXIT line_loop  ! End of file encountered
            if (ios > 0) STOP '>>> READAT: File "' // TRIM(file_in) // '" reading error'
            call split_string(STRI, '=', str_nam, str_val)
            if (LEN_TRIM(str_val) == 0) CYCLE line_loop
            if ( TRIM(arr_in(j)) == str_nam(1: LEN_TRIM(str_nam)) ) then
                read(str_val, *) arr_out(j)
                EXIT
            endif
        enddo line_loop
        close(171)
    enddo

    n_dim_next = 0
    jpos = index(STRI, ':')
    if (jpos > 0) then
        if (LEN_TRIM(STRI(jpos+1: )) > 0) then
            read(STRI(jpos+1: ), *) n_dim_next
        endif
    endif

    end subroutine assign_val

!------------------------------------------------------------
    subroutine path_split(str_path_in, dir_path, fname, jpos)

    character(len=*), intent(in)  :: str_path_in
    integer, intent(out) :: jpos
!character(len=len(str_path_in)), intent(out) :: dir_path, fname
    character(len=*), intent(out) :: dir_path, fname

    integer :: len_str, jstr
    character(len=len(str_path_in)) :: str_path

    jpos = 0
    str_path = clean_string(str_path_in)
    len_str = LEN_TRIM(str_path)
! If delimiter is not there
    jpos = 0
    do jstr=1, len_str
        if (str_path(jstr: jstr) == '/') then
            jpos = jstr
        endif
    enddo
    if (jpos > 0) then
        dir_path = TRIM(str_path(1: jpos-1))
        fname = ADJUSTL(TRIM(str_path(jpos+1: ))) ! whole string if delimiter not found
    else
        dir_path = ''
        fname = TRIM(str_path)
    endif

    end subroutine path_split

!------------------------------------------------------------
    subroutine split2array(str_in, delim, strarray, nout)
! splitstring splits a string to an array of
! substrings based on a selected delimiter
! note any facing space/blank in substrings will be removed

    character(len=*), intent(in) :: str_in
    character, intent(in) :: delim
    integer, intent(out) :: nout
    character(len=len(str_in)), intent(out) :: strarray(20)

    integer :: m, i, jpos
    character(len=len(str_in)):: strtmp
    
    strtmp = clean_string(str_in)

    do i=1, 10
        strarray(i) = repeat(' ', len(str_in))
    enddo

    m = 1
    nout = 0
    do i=2, 10
        jpos = index(strtmp(m:), delim)
        if ( LEN_TRIM(strtmp(m: m+jpos-1)) > 0 ) then
            nout = nout + 1
            strarray(nout) = TRIM(ADJUSTL( strtmp(m: m+jpos-2) )) 
        endif
        m = m + jpos
    enddo

! After the last delimiter
    if ( LEN_TRIM(strtmp(m: )) > 0 ) then
        nout = nout + 1
        strarray(nout) = TRIM(ADJUSTL( strtmp(m:) )) 
    endif

    end subroutine split2array

!------------------------------------------------------------
    subroutine split2array2(str_in, strarray, nout)
! splitstring splits a string to an array of
! substrings using blank and/or tab as delimiter
! note any facing space/blank in substrings will be removed

    use char_manip, only: null_ch, tab_ch

    integer, parameter :: nwords_max=20
    character(len=*), intent(in) :: str_in
    integer, intent(out) :: nout
    character(len=len(str_in)), intent(out) :: strarray(nwords_max)

    integer :: m, i, jpos_null, jpos_tab, jpos
    character(len=len(str_in)):: strtmp

    strtmp = clean_string(str_in)

    do i=1, nwords_max
        strarray(i) = repeat(' ', len(str_in))
    enddo

    m = 1
    nout = 0
    do i=1, len(strtmp)
        jpos_null = index(strtmp(m:), ' ')
        jpos_tab  = index(strtmp(m:), tab_ch)
        if (jpos_null > 0) then
            if (jpos_tab > 0) then
                jpos = min(jpos_null, jpos_tab)
            else
                jpos = jpos_null
            endif
        else
            jpos = jpos_tab
        endif
        if ( LEN_TRIM(strtmp(m: m+jpos-1)) > 0 ) then
            nout = nout + 1
            strarray(nout) = TRIM(ADJUSTL( strtmp(m: m+jpos-2) )) 
        endif
        m = m + jpos
    enddo

! After the last delimiter
    if ( LEN_TRIM(strtmp(m: )) > 0 ) then
        nout = nout + 1
        strarray(nout) = TRIM(ADJUSTL( strtmp(m:) )) 
    endif

    end subroutine split2array2

!------------------------------------------------------------
    subroutine inquire_fname(ftype, fdefault, dev_name, fname)

    character(len=3), intent(in) :: ftype
    character(len=*), intent(in) :: dev_name, fdefault
    character(len=80), intent(out) :: fname

    logical :: EXI

    fname = 'exp/' // TRIM(ftype) // '/' // TRIM(fdefault) ! Check default name
    inquire(FILE=TRIM(fname), EXIST=EXI)

    if (.not. EXI ) then
! Check generic names: AUG, FEAT, ITER, JET
        fname = 'exp/' // TRIM(ftype) // '/' // TRIM(to_lower(dev_name))
        inquire(FILE=TRIM(fname), EXIST=EXI)
        if ( .not. EXI ) fname = '***'
    endif

    end subroutine inquire_fname

end module parse_utils
