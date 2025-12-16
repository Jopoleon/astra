!---------------------------------------------------------------------
double precision function GETNUM(FIELD, ERCODE)

use const_inc , only: constValues, varxValues
use char_manip, only: str_in_list
use json_vars, only: constNames, varNames

implicit none

integer, intent(out) :: ERCODE
character(len=*), intent(in) :: FIELD

integer :: l, j, jnam, jpos, jpos1, j1, ISHIFT, ios
character(len=6) :: ZNUM

save ISHIFT
data ISHIFT/0/

if (ISHIFT == 0) then
    j = str_in_list('ZRD1  ', varNames)
    ISHIFT = max(j-1, 0)
endif

l  = len(FIELD)
jpos = index(TRIM(FIELD), 'ZRD')

if (jpos == 0) then
    jpos1 = index(TRIM(FIELD), 'C')
    ERCODE = 1 ! Only initial, before constNames-search
    if (jpos1 >  0) then
        j1 = min(LEN_TRIM(FIELD(jpos1: l)), l - jpos1 + 1)
        ZNUM = FIELD(jpos1: jpos1+j1-1)
        jnam = str_in_list(ZNUM, constNames)
        if (jnam > 0) then
            ERCODE = 0
            GETNUM = constValues(jnam)
        endif
    endif
else
    j1 = index(FIELD(jpos+3: l), 'X')
    if (j1 /= 0) then
        if (j1 > 3) then
            ERCODE = 1
            return
        endif
        ZNUM = FIELD(jpos+3: jpos+j1+1)
    else
        ZNUM = FIELD(jpos+3:)
    endif

    read(ZNUM(1:), *, iostat=ios) j1

    if (ios /= 0 .or. j1 < 1 .or. j1 > 96) then
        ERCODE = 1
    else
        GETNUM = varxValues(ISHIFT + j1)
        ERCODE = 0
    endif
endif

return
end function GETNUM

!---------------------------------------------------------------------
subroutine STREAD(NCH, NFIELD, ARRAY, ERCODE)
!---------------------------------------------------------------------
! Reads one record group of the NBINP ("*.nbi") file 
! and fills ARRAY(1:NFIELD) with data.
! Numbers and references to ZRD*, ZRD*X and constValues_list
! are allowed as records in the input file.
! ERCODE values:
!   0 - Normal exit
!   1 - Unrecognized variable name
!   2 - Read error (never occurring, protected by ISNUM)
!   3 - Wrong format
!   4 - Array out of limits
!   5 - Missing records
!---------------------------------------------------------------------

use char_manip, only: to_upper
use dbl2char, only: isnum

implicit none

integer, intent(in) :: NCH, NFIELD
integer, intent(out) :: ERCODE
double precision, intent(out) :: ARRAY(NFIELD)

integer :: j, ios
double precision, external :: GETNUM
character(len=12 ) :: SFIELD(20)
character(len=132) :: str_line
!---------------------------------------------------------------------

if (NFIELD > 20) then
    ERCODE = 4
    return
endif

! Skip all lines beginning with '!'
j = 1
do while (j == 1)
    read(NCH, '(A)', iostat=ios) str_line
    if (ios < 0) then ! EOF encountered
        ERCODE = 5
        return
    else if (ios > 0) then
        ERCODE = 1
        return
    endif
    j = index(str_line, '!')
enddo

if (j /= 0) then
    ERCODE = 3
    write(*, *) 'STREAD error: exclamation marks allowed only at line beginning'
else  ! Read numbers and/or variable names
    ERCODE = 5  ! Missing entries
    read(NCH, '(5A)', iostat=ios) (SFIELD(j), j=1, NFIELD)
    if (ios < 0) then ! EOF encoutnered
        ERCODE = 5
    else if (ios > 0) then
        ERCODE = 1      ! Error reading, probably never occurring
    else
        ERCODE = 0
        do j=1, NFIELD
            SFIELD(J) = to_upper(SFIELD(j))
            if ( ISNUM(SFIELD(j), 12) ) then
                read(SFIELD(j), *) ARRAY(j) ! read err never occurs, protected by ISNUM
            else ! In case it is a variable name, like ZRD*, pick its value
                ARRAY(j) = GETNUM(SFIELD(j), ERCODE)
                if (ERCODE /= 0) then ! Unrecognised variable name
                    ERCODE = 1
                    EXIT
                endif
            endif
        enddo
    endif
endif

return
end subroutine STREAD
