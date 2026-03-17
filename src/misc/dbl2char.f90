module dbl2char

implicit none

contains

!-------------------------------------------
   logical function ISNUM(str_in, N)
   ! Checks whether a string is actually a float

   integer, intent(in) :: N
   character(len=*), intent(in) :: str_in

   integer :: ios
   double precision :: R

   ISNUM = .True.
   if (N < 1) return
   read(str_in(1: N), *, iostat=ios) R

   ISNUM = (ios == 0)

   end function ISNUM

!---------------------------------------------------------------------
    function fmt_smart(x, width) result(out)

    double precision, intent(in) :: x
    integer, intent(in) :: width
    character(len=width) :: out
    character(len=32) :: tmp, s
    integer :: prec, expo
    logical :: fitted
    double precision :: r, mant

    if (width < 4) stop "fmt_smart: width too small"

! Handle very small numbers
    if (abs(x) < 1.0d-12) then
        write(out, '(A)') adjustl('0.'//repeat('0', width-2))
        return
    endif

! Try fixed-point first
    fitted = .false.
    do prec=width-1, 1, -1
       write(tmp, '(F0.'//trim(adjustl(itoa(prec)))//')') x
! If rounding produced zero but x is not zero → reject fixed format
       read(tmp, *) r
       if (r == 0.0d0 .and. x /= 0.0d0) cycle
! Remove leading zero for numbers between -1 and 1
        if (abs(x) < 1.0d0 .and. x /= 0.0d0) then
            if (tmp(1:1) == '0') tmp = tmp(2:)   ! remove leading 0
            if (tmp(1:1) == '-') tmp = '-'//tmp(2:) ! keep minus
        endif
        if (len_trim(tmp) <= width) then
            write(out, '(A)') adjustl(TRIM(tmp))
            fitted = .true.
            exit
        endif
    enddo

! Leave as it is (F) if width=4 and 0.001 < abs(x) < 1000.
    if (.not. ((width == 4 .and. abs(x) >= 1.0d-3 .and. abs(x) < 1.0d3) .or. &
               (width == 5 .and. abs(x) >= 1.0d-3 .and. abs(x) < 1.0d3) .or. &
               (width == 6 .and. abs(x) >= 1.0d-4 .and. abs(x) < 1.0d4))  ) then
        do prec=width-2, 0, -1
            expo = floor(log10(abs(x)))
            mant = x / 10.0d0**expo
            write(s, '(F0.' // trim(adjustl(itoa(prec))) // ')') mant

! remove leading zero
            if (s(1:1) == '0') s = s(2:)
            if (s(1:1) == '-' .and. s(2:2) == '0') s = '-'//s(3:)

! remove trailing decimal point
            if (index(s, '.') > 0) then
                if (s(len_trim(s): len_trim(s)) == '.') then
                    s = s(:len_trim(s)-1)
                endif
            endif
 
            if (expo == 0) then
                tmp = trim(s)
            else
                tmp = trim(s) // 'e' // trim(adjustl(itoa(expo)))
            endif
 
            if (len_trim(tmp) <= width) then
                out = adjustl(tmp)
                fitted = .true.
                exit
            endif
        enddo
    endif 
 
! Fallback: truncate if nothing fits
    if (.not. fitted) then
        write(out, '(A)') tmp(1:width)
    endif

    end function fmt_smart

!---------------------------------------------------------------------
    pure function itoa(i) result(str)

    integer, intent(in) :: i
    character(len=3) :: str
    write(str, '(I0)') i

    end function itoa

end module dbl2char
