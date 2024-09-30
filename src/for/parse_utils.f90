module parse_utils
    use machine_config, only: config
  
implicit none

contains

!------------------------------------------------------------
    logical function IFDEFX(XARNAM)
! Name exists in profxNames, and the array is defined
    
    use parameter_inc, only: NARRX
    use outcmn_inc, only: IFDFAX
    use char_manip, only: str_in_list
    use json_vars, only: profxNames

    character(len=6), intent(in) :: XARNAM

    integer :: j
    
    IFDEFX = .false.
    j = str_in_list(XARNAM, profxNames)
    if (j > 0) then
        if (IFDFAX(j) /= -1) IFDEFX = .true. ! True (X-array is defined)
    endif

    return
    end function IFDEFX

!----------------------------------------------------------
    logical function IFDEFX2(XARNAM)
! Name exists in profxNames, and the array is defined

    use parameter_inc, only: NARRX
    use outcmn_inc, only: IFDFAX
    use char_manip, only: str_in_list
    use json_vars, only: profxNames

    character(len=6), intent(in) :: XARNAM

    integer :: j

    j = str_in_list(XARNAM, profxNames)
    if (j > 0) then
        if (IFDFAX(j) > 0) IFDEFX2 = .true. ! True (X-array is defined)
    endif

    return
    end function IFDEFX2

!--------------------------------------------------
! Set variable list from file parsing
    SUBROUTINE assign_val(file_in, narr, arr_in, arr_out, n_dim_next)

    use debugger, only: markloc, astra_stop

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
            if (ios > 0) then
                call astra_stop('>>> READAT: File "' // TRIM(file_in) // '" reading error')
            endif
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

    return
    END SUBROUTINE assign_val

!------------------------------------------------------------
    subroutine path_split(str_path_in, dir_path, fname, jpos)

    use char_manip, only: clean_string
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

    return
    end subroutine path_split

!------------------------------------------------------------
    subroutine split2array(str_in, delim, strarray, nout)
! splitstring splits a string to an array of
! substrings based on a selected delimiter
! note any facing space/blank in substrings will be removed

    use char_manip, only: clean_string

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

    return
    end subroutine split2array

!------------------------------------------------------------
    subroutine split2array2(str_in, strarray, nout)
! splitstring splits a string to an array of
! substrings based on a selected delimiter
! note any facing space/blank in substrings will be removed

    use char_manip, only: clean_string
    use outcmn_inc, only: null_ch, tab_ch

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
       jpos_tab = index(strtmp(m:), tab_ch)
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

    return
    end subroutine split2array2

!------------------------------------------------------------
    subroutine read_arrx(nunit, nt_io, ntim, nrho, stri_in, var_out)

    use outcmn_inc, only: exp_file, NCNBM, NCNBTM
    use debugger, only: markloc, astra_stop

    integer, intent(in) :: nunit, ntim
    character(len=*), intent(in) :: stri_in
    integer, intent(inout) :: nt_io
    integer, intent(out) :: nrho

    integer :: j, ios
    character(132) :: err_msg
    double precision, dimension((NCNBM+1)*NCNBTM), intent(out) :: var_out

    call markloc('read_arrx')

    err_msg =  '>>> Data file "' // TRIM(exp_file) // '" error:\n'

    if (nt_io /= 0) then
        err_msg = TRIM(err_msg) // '    Boundary must be defined in a single group'
        call astra_stop(err_msg)
    endif

    nt_io = max(ntim, 1)

    j = INDEX(stri_in, 'POINTS')
    if (j == 0) then
        err_msg = TRIM(err_msg) // '    Number of boundary points must be defined'
        call astra_stop(err_msg)
    endif
    if (nrho > NCNBM) then
        write(err_msg, '(2A, i)') TRIM(err_msg), &
           '    Number of boundary points must be NBND <', NCNBM
        call astra_stop(err_msg)
    endif

    read(stri_in(j+6:), *) nrho

    if ((nrho + 1)*nt_io > (NCNBM + 1)*(NCNBTM - 10)) then
        write(err_msg, '(2A)') TRIM(err_msg), &
           '    Boundary data length must be NBND*NBNT <'
        call astra_stop(err_msg)
    endif

    read(nunit, *, iostat=ios)(var_out(j), j=1, (nrho+1)*nt_io)
    if (ios /= 0) then
        err_msg = TRIM(err_msg) // '    More data items than data values for BND group'
        call astra_stop(err_msg)
    endif

    return
    end subroutine read_arrx

!------------------------------------------------------------
    subroutine parse_u_line(str_in, var_name, uname, factor)

    use char_manip, only: to_upper
    use debugger, only: markloc, astra_stop

    character(len=*), intent(in) :: str_in
    character(len=len(str_in)), intent(out) :: var_name, uname
    double precision, intent(out) :: factor

    integer :: n_words
    character(len=132) :: str1, strarray(20), err_msg

    call markloc('parse_u_line')

    uname = repeat(' ', 40)
    err_msg = 'Error in exp-file line ' // TRIM(str_in)

    str1   = repeat(' ', 132)

    call split_string(str_in, ' ', var_name, str1)
    if (LEN_TRIM(var_name) == 0) then
       call astra_stop(err_msg)
    endif

    call split2array(str1, ':', strarray, n_words)

    if (to_upper(TRIM(strarray(1))) /= 'U-FILE' ) then
       return
    endif

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
    subroutine ufheader(uname, n_dim, nt, nx, lbl2)

    use char_manip, only: to_upper
    use debugger, only: markloc, astra_stop

    integer, intent(out) :: nt, nx, n_dim
    character(len=*), intent(in) :: uname
    character(len=30), intent(out) :: lbl2

    integer :: ios, j, n_scal, ISHOT
    character(132) :: err_msg
    character(32) :: STRI
    character(30) :: lbl1, lbl3, var1_lbl, unit1
    character(4) :: sdev

    call markloc('ufheader')

    var1_lbl = repeat(' ', 30)
    unit1    = repeat(' ', 30)

    open(11, FILE=TRIM(uname), iostat=ios)

    if (ios /= 0) then
       err_msg = '>>> READAT: U-file "' // TRIM(uname) // '" reading error'
       call astra_stop(err_msg)
    endif

! # shot, device, #dimensions

    read(11,'(A32)', ERR=925) STRI
    read(STRI(3:7)  , '(1I5)', ERR=925) ISHOT
    read(STRI(8:11) , '(1A4)', ERR=925) sdev
    read(STRI(13:13), '(1I1)', ERR=925) n_dim

    write(*,*) ishot, sdev, n_dim
    if (n_dim <= 0 .or. n_dim > 2) then
        err_msg = '>>> U-file "' // TRIM(uname) // '" error: wrong dimensionality'
        call astra_stop(err_msg)
    endif

    read(11, '(A32)', ERR=925) STRI ! Dummy line

! Scalar quantities

    read(11, *, ERR=925) n_scal
    if (n_scal > 0) then
        do j=1, n_scal
            read(11, '(A32)', ERR=925) STRI
            read(11, '(A32)', ERR=925) STRI
        enddo
    endif

! Continue reading 2D U-file
! 1st independent variable label: X-
    read(11, '(A32)', ERR=925) STRI
    STRI = ADJUSTL(STRI)
    lbl1 = to_upper(STRI(1:30))

    call split_string(TRIM(lbl1), ' ', var1_lbl, unit1)
    if (var1_lbl(1:4) /= 'TIME') then
        write(*,*) '>>> U-file "', TRIM(uname), '" 1st independent variable should be time'
        close(11)
        return
    endif
    if (TRIM(unit1) /= 'SECONDS') then
        write(*, *) '>>> U-file "', TRIM(uname),'" time unit should be second, not ', TRIM(unit1)
    endif

! 2nd independent variable label: Y-
    if (n_dim == 2) then
        read(11, '(A32)', ERR=925) STRI
        STRI = ADJUSTL(STRI)
        lbl2 = to_upper(STRI(1:30))
    endif

! Dependent variable label
    read(11,'(A32)', ERR=925) STRI
    STRI = ADJUSTL(STRI)
    lbl3 = to_upper(STRI(1:30))

! Dummy, "PROC CODE"

    read(11,'(A32)', ERR=925) STRI
! Dimensions
    read(11, *, ERR=925) nt
    if (n_dim == 2) then
        read(11, *, ERR=925) nx
    else
        nx = 1
    endif

    close(11)

    return

925 call astra_stop('>>> U-file "' // TRIM(uname) // '" read error')

    end subroutine ufheader

!------------------------------------------------------------
    subroutine ufrd(uname, n_dim, nt, nx, t_out, x_out, arr_out)

    use debugger, only: markloc, astra_stop

    character(len=*), intent(in) :: uname
    integer, intent(in) :: nt, nx, n_dim
    double precision, intent(out) :: t_out(nt), x_out(nx), arr_out(nt*nx)

    integer :: ios, j, jj, n_scal
    character(132) :: err_msg
    character(32) :: STRI

    call markloc('ufrd')

    write(*, *) 'Reading u-file ' // TRIM(uname)

    open(11, FILE=TRIM(uname), iostat=ios)

    if (ios /= 0) then
        err_msg = '>>> READAT: U-file "' // TRIM(uname) // '" reading error'
        call astra_stop(err_msg)
    endif

!-------
! Header
!-------

! # shot, device, #dimensions
    read(11,'(A32)') STRI

    if (n_dim <= 0 .or. n_dim > 2) then
       err_msg = '>>> U-file "' // TRIM(uname) // '" error: wrong dimensionality'
       call astra_stop(err_msg)
    endif

    read(11, '(A32)') STRI ! Dummy line

! Scalar quantities

    read(11, *) n_scal
    if (n_scal > 0) then
        do j=1, n_scal
            read(11, '(A32)') STRI
            read(11, '(A32)') STRI
        enddo
    endif

! Continue reading 2D U-file
! 1st independent variable label: X-
    read(11, '(A32)') STRI

! 2nd independent variable label: Y-
    if (n_dim == 2) then
        read(11, '(A32)') STRI
    endif

! Dependent variable label
    read(11,'(A32)') STRI

! Dummy, "PROC CODE"
    read(11,'(A32)') STRI

! Dimensions
    read(11, *) STRI
    if (n_dim == 2) then
        read(11, *) STRI
    endif

! Read grid and data arrays

    read(11, *) (t_out(j), j=1, nt)
    if (n_dim == 2) then
        read(11, *) (x_out(j) , j=1, nx)
    endif
    read(11, '(1X, 6E13.6)') ((arr_out(jj + (j - 1)*nt), jj=1, nt), j=1, nx)

    close(11)

    return
    end subroutine ufrd

!------------------------------------------------------------
    subroutine inquire_fname(ftype, fdefault, dev_name, fname)

    use char_manip, only: to_lower

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

    return
    end subroutine inquire_fname

end module parse_utils
