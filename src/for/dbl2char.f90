module dbl2char

   implicit none

   contains

!-------------------------------------------
   function fmt4(R1)
! Convert real positive R1>0 to a string dbl_4*4

   double precision, intent(in) :: R1
   character(len=4) :: fmt4

   integer :: J
   character(len=6) :: F6

   if (abs(R1) < 1.e-9) then
      fmt4 = ' 0. '
      return
   endif
   if (R1 < 0) then
      fmt4 = ' <0 '
      return
   endif

   F6 = fmt456(R1, 4)

   if ( index(F6(3: 6), '.') > 0 ) then
      do J=6, 1, -1
         if (F6(J: J) == '.') then
            fmt4 = F6(3: J-1)
            return
         else if (F6(J: J) /= '0' ) then
            fmt4 = F6(3: J)
            return
         endif
      enddo

   else ! No '.' in F6(3:6)

      fmt4 = F6(3:6)
      if (fmt4(2: 3) == '0-') then
         read(fmt4(4: 4), *) J
         write(fmt4, '(A1, A2, I0)') fmt4(1: 1), 'e-', J - 1
      endif
      if (fmt4(1: 1) == ' ' .and. fmt4(3: 4) == '-9') write(fmt4, '(A1, A3)') fmt4(2: 2), 'e-9'
      if (fmt4(2: 4) == '0-4') fmt4 ='.00' // fmt4(1: 1)
      if (fmt4(2: 3) == '0e') then
         read(fmt4(4: 4), *) J
         write(fmt4, '(2A1, I0)') fmt4(1: 1), 'e', J+1
      endif
   endif

   return
   end function fmt4

!-------------------------------------------
   function fmt40(R1)
! The same as FMTF4 but does not suppress leading "0." if possible

   double precision, intent(in) :: R1
   character(len=4) :: fmt40

   fmt40 = fmt4(R1)
   if (fmt40(1: 1) == '.' .and. fmt40(4:4) == ' ') fmt40 = '0.' // fmt40(2:3)

   return
   end function fmt40

!-------------------------------------------
   function fmt5(R1)
! Convert real R1 to a string F5*5

   double precision, intent(in) :: R1
   character(len=5) :: fmt5

   character(len=4) :: F4
   character(len=6) :: F6

   if (abs(R1) < 1.e-9) then
      fmt5 = ' 0.  '
   else if (R1 < 0) then
      F4 = fmt4(abs(R1))
      fmt5 = '-' // F4
      return
   else
      F6 = fmt456(R1, 5)
      fmt5 = F6(2: 6)
      if (fmt5(3: 5) == "000") fmt5(3: 5) = '   '
      if (fmt5(4: 5) == "00" ) fmt5(4: 5) = '  '
      if (fmt5(5: 5) == "0"  ) fmt5(5: 5) = ' '
   endif

   return
   end function fmt5

!-------------------------------------------
   function fmt50(R1)
! Convert real R1 to a string F5*5
! The same as FMTF5 but without suppressing end zeros

   double precision, intent(in) :: R1
   character(len=5) :: fmt50

   character(len=4) :: F4
   character(len=6) :: F6

   if (abs(R1) < 1.e-9) then
      fmt50 = ' 0.  '
   else if (R1 < 0) then
      F4 = fmt4(abs(R1))
      fmt50 = '-' // F4
   else
      F6 = fmt456(R1, 5)
      fmt50 = F6(2: 6)
   endif

   return
   end function fmt50

!-------------------------------------------
   function fmt6(R1)
! Convert real R1 to a string F6*6

   double precision, intent(in) :: R1
   character(len=6) :: fmt6

   character(len=5) :: F5

   if (abs(R1) < 1.e-9) then
      fmt6 = '    0.'
      return
   endif
   if (R1 < 0) then
      F5 = fmt5(abs(R1))
      fmt6 = '-' // F5
   else
      fmt6 = fmt456(R1, 6)
   endif
   ! Right adjust
   if (fmt6(3: 6) == "0000") fmt6 = '    ' // fmt6(1:2)
   if (fmt6(4: 6) == "000" ) fmt6 = '   '  // fmt6(1:3)
   if (fmt6(5: 6) == "00"  ) fmt6 = '  '   // fmt6(1:4)
   if (fmt6(6: 6) == "0"   ) fmt6 = ' '    // fmt6(1:5)
   if (fmt6(3: 6) == "    ") fmt6 = '    ' // fmt6(1:2)
   if (fmt6(4: 6) == "   " ) fmt6 = '   '  // fmt6(1:3)
   if (fmt6(5: 6) == "  "  ) fmt6 = '  '   // fmt6(1:4)
   if (fmt6(6: 6) == " "   ) fmt6 = ' '    // fmt6(1:5)

   return
   end function fmt6

!-------------------------------------------
   function fmt_xf(R1, nlen)
! nlen=4, 5 or 6

   integer, intent(in) :: nlen
   double precision, intent(in) :: R1
   character(len=nlen+1) :: fmt_xf

   character(len=6) :: F6

   fmt_xf = ' '
   if (R1 < 0.) fmt_xf = '-'

   F6 = fmt456(abs(R1), nlen)
   if (nlen == 4) then
      if (F6(4: 6) == '0-4') F6(3: 6) = '.00' // F6(3: 3)
   endif
   fmt_xf(2: nlen+1) = F6(7-nlen: 6)

   return
   end function fmt_xf

!-------------------------------------------
   function fmx_5(R1)
! the same as FMTXF5 but suppresses right "." or "0" when possible

   double precision, intent(in) :: R1
   character(len=6) :: fmx_5

   integer :: J, JM, J1
   character(len=6) :: F6

   if (R1 < 0.) then
      JM = 1
      fmx_5 = '-'
   else
      fmx_5 = ' '
      JM = 0
   endif
   F6 = fmt456(abs(R1), 6 - JM)

   if ( index(F6(1:6), '.') == 0) then
      fmx_5(1+JM: 6) = F6(1+JM: 6)
   else
      do J=6, 1, -1
         if (F6(J:J) == '.') then
            J1 = J - 1
            EXIT
         endif
         if (F6(J:J) /= '0') then
            J1 = J
            EXIT
         endif
      enddo
   endif

   fmx_5(7+JM-J1: 6) = F6(1+JM: J1)
   if (R1 < 0.) then
      fmx_5(1: 1) = ' '
      fmx_5(6+JM-J1: 6+JM-J1) = '-'
   endif

   return
   end function fmx_5

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
   function fmt456(R2, N)

   use debugger, only: astra_stop

! R2>0 only

   integer, intent(in) :: N
   double precision, intent(in) :: R2
   character(len=6) :: fmt456

   integer :: J, JM
   double precision :: R
   character(len=25) :: T ! T*25 in case of crash

   if (R2 < 0) then
      call astra_stop('Improper use of routine fmt456, input must be >= 0')
   endif

   JM = 6 - N
   fmt456 = ' '
   if (R2 < 1.e-9) then
      fmt456 = ' 0.000'
      return
   endif

   R = ROUNDN(R2, 5 - JM)

   if (R < 1.e-4 * 10.**JM) then
      write(T, '(1F25.11)') R
      do J=15, 20+JM
         if (T(J: J) /= '0') EXIT
      enddo
      J = min(J, 20 + JM) ! At the end of the loop it is J=21+JM
      R = ROUNDN(R2, 4 - JM)
      write(T, '(1F25.11)') R
      fmt456(1+JM: 5) = T(J: J+2-JM) // 'e-'
      write(fmt456(6: 6), '(I1)') J - JM - 12 ! GIT J - 11 - JM

!      write(6, *) 'GIT1 fmt456', R2, fmt456
      return
   endif
   
   R = ROUNDN(R2, 6 - JM)
   if (R < 1.e6/10.**JM) then
      if (R >= 1.e5/10.**JM) R = ROUNDN(R2, 7 - JM)
      write(T, '(1F25.11)') R
      do J = 8 + JM, 17 - JM
         if (T(J:J) /= ' ' .and. T(J:J) /= '0') EXIT
      enddo
      fmt456(1+JM: 6) = T(J: J+5-JM)
      return
   endif

   R = ROUNDN(R2, 5 - JM)
   if (R < 1.e13/10.**JM) then
      write(T, '(1F25.11)') R
      do J=1, 11
         if (T(J:J) /= ' ' .and. T(J:J) /= '0') EXIT
      enddo
      fmt456(1+JM: 5) = T(J: J+3-JM) // 'e'
      write(fmt456(6: 6), '(I1)') 10 - J + JM
      return
   endif

   fmt456 = '  1e1 '
   write(fmt456(6: 6), '(I1)') 3 - JM

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
