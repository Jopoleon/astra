module dbl2char

implicit none

contains

!-------------------------------------------
    function fmt4(R1) result(out)
! Convert real positive R1>0 to a string dbl_4*4

    double precision, intent(in) :: R1
    character(len=4) :: out

    integer :: J
    character(len=6) :: F6

    if (abs(R1) < 1.e-9) then
        out = ' 0. '
        return
    endif
    if(R1 < 0) then
        out = ' <0 '
        return
    endif

    F6 = fmt456(R1, 4)

    if ( index(F6(3: 6), '.') > 0 ) then
        do J=6, 1, -1
            if (F6(J: J) == '.') then
                out = F6(3: J-1)
                return
            else if(F6(J: J) /= '0' ) then
                out = F6(3: J)
                return
            endif
        enddo

    else ! No '.' in F6(3:6)
        out = F6(3:6)
        if (out(2: 3) == '0-') then
            read(out(4: 4), *) J
            write(out, '(A1, A2, I0)') out(1: 1), 'e-', J - 1
       endif
       if (out(1: 1) == ' ' .and. out(3: 4) == '-9') write(out, '(A1, A3)') out(2: 2), 'e-9'
       if (out(2: 4) == '0-4') out ='.00' // out(1: 1)
       if (out(2: 3) == '0e') then
           read(out(4: 4), *) J
           write(out, '(2A1, I0)') out(1: 1), 'e', J+1
       endif
    endif

    return
    end function fmt4

!-------------------------------------------
    function fmt40(R1) result(out)
! The same as FMTF4 but does not suppress leading "0." if possible

    double precision, intent(in) :: R1
    character(len=4) :: out

    out = fmt4(R1)
    if (out(1: 1) == '.' .and. out(4:4) == ' ') out = '0.' // out(2:3)

    return
    end function fmt40

!-------------------------------------------
    function fmt5(R1) result(out)
! Convert real R1 to a string F5*5

    double precision, intent(in) :: R1
    character(len=5) :: out

    character(len=4) :: F4
    character(len=6) :: F6

    if (abs(R1) < 1.e-9) then
        out = ' 0.  '
    else if (R1 < 0) then
        F4 = fmt4(abs(R1))
        out = '-' // F4
    else
        F6 = fmt456(R1, 5)
        out = F6(2: 6)
        if (out(3: 5) == "000") out(3: 5) = '   '
        if (out(4: 5) == "00" ) out(4: 5) = '  '
        if (out(5: 5) == "0"  ) out(5: 5) = ' '
    endif

    return
    end function fmt5

!-------------------------------------------
    function fmt50(R1) result(out)
! Convert real R1 to a string F5*5
! The same as FMTF5 but without suppressing end zeros

    double precision, intent(in) :: R1
    character(len=5) :: out

    character(len=4) :: F4
    character(len=6) :: F6

    if (abs(R1) < 1.e-9) then
        out = ' 0.  '
    else if (R1 < 0) then
        F4 = fmt4(abs(R1))
        out = '-' // F4
    else
        F6 = fmt456(R1, 5)
        out = F6(2: 6)
    endif

    return
    end function fmt50

!-------------------------------------------
    function fmt_xf(R1, nlen) result(out)
! nlen=4, 5 or 6

    integer, intent(in) :: nlen
    double precision, intent(in) :: R1
    character(len=nlen+1) :: out

    character(len=6) :: F6

    out = ' '
    if (R1 < 0.) out = '-'

    F6 = fmt456(abs(R1), nlen)
    if (nlen == 4) then
        if (F6(4: 6) == '0-4') F6(3: 6) = '.00' // F6(3: 3)
    endif
    out(2: nlen+1) = F6(7-nlen: 6)

    return
    end function fmt_xf

!-------------------------------------------
    double precision function ROUNDN(R1, N)

    integer, intent(in) :: N
    double precision, intent(in) :: R1
   
    integer :: J
    double precision :: R
   
    if (R1 >= 1.e13 .or. R1 < 1.e-9) then
        ROUNDN = R1
    else
        R = R1*1.e9
        do J=1, 25
            if (R < 10.) EXIT
            R = R/10.
        enddo
        ROUNDN = (R + 50./10.**N)*10.**(J - 10)
    endif

    return
    end function ROUNDN
   
!-------------------------------------------
    function fmt456(R2, N) result(out)

    use debugger, only: astra_stop

! R2 > 0 only

    integer, intent(in) :: N
    double precision, intent(in) :: R2
    character(len=6) :: out

    integer :: J, JM
    double precision :: R
    character(len=25) :: T ! T*25 in case of crash

    if (R2 < 0) then
        call astra_stop('Improper use of routine fmt456, input must be >= 0')
    endif

    JM = 6 - N
    out = ' '
    if (R2 < 1.e-9) then
        out = ' 0.000'
        return
    endif

    R = ROUNDN(R2, 5 - JM)

    if (R < 1.e-4 * 10.**JM) then
        write(T, '(1F25.11)') R
        do J=15, 20+JM
            if(T(J: J) /= '0') EXIT
        enddo
        J = min(J, 20 + JM) ! At the end of the loop it is J=21+JM
        R = ROUNDN(R2, 4 - JM)
        write(T, '(1F25.11)') R
        out(1+JM: 5) = T(J: J+2-JM) // 'e-'
        write(out(6: 6), '(I1)') J - JM - 12 ! GIT J - 11 - JM
!        write(6, *) 'GIT1 fmt456', R2, out
       return
    endif
   
    R = ROUNDN(R2, 6 - JM)
    if (R < 1.e6/10.**JM) then
        if (R >= 1.e5/10.**JM) R = ROUNDN(R2, 7 - JM)
        write(T, '(1F25.11)') R
        do J = 8 + JM, 17 - JM
            if(T(J:J) /= ' ' .and. T(J:J) /= '0') EXIT
        enddo
        out(1+JM: 6) = T(J: J+5-JM)
        return
    endif

    R = ROUNDN(R2, 5 - JM)
    if (R < 1.e13/10.**JM) then
        write(T, '(1F25.11)') R
        do J=1, 11
            if(T(J:J) /= ' ' .and. T(J:J) /= '0') EXIT
        enddo
        out(1+JM: 5) = T(J: J+3-JM) // 'e'
        write(out(6: 6), '(I1)') 10 - J + JM
        return
    endif

   out = '  1e1 '
   write(out(6: 6), '(I1)') 3 - JM

   return
   end function fmt456

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

end module dbl2char
