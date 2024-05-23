subroutine SCAL(NOUT, SN, SO, OUT, NP, NDIM)
!---------------------------------------------------------------------
!  Input:
!     NOUT Total number of channels (NROUT or NTOUT)
!     SO(NOUT) Current scales
!     NP Upper index boundary for scale definition (= NDIM-3)
!     NDIM Upper index boundary according to the array description
!     OUT(NDIM, NOUT) Array of curves for scale definition (ROUT or TOUT)
!  Output:
!     SN(NOUT) Returned scales

implicit none

integer, intent(in) :: NOUT, NP, NDIM
double precision, intent(in)  :: SO(*), OUT(NDIM, *)
double precision, intent(out) :: SN(*)

integer :: J, JJ
double precision :: SC, SCALA
external SCALA

! Exclude 3 edge/central points
var_loop: do J=1, NOUT
    if (SO(J) > 0) then
        SN(J) = SO(J)
    else if (SO(J) .EQ. 0) then
        SN(J) = SCALA(OUT(1, J), NP)
    else
        do JJ = 1, J-1
            if (SO(J) .EQ. SO(JJ)) then
                SN(J) = SN(JJ)
                CYCLE var_loop
            endif
        enddo
        SC = 0.
        do JJ=J, NOUT
            if (SO(J) .EQ. SO(JJ)) then
                SC = MAX(SC, SCALA(OUT(1, JJ), NP))
                if (ABS(SO(JJ) + JJ) < .01) then
                    SN(J) = SCALA(OUT(1, JJ), NP)
                    CYCLE var_loop
                endif
            endif
        enddo
        SN(J) = SC
    endif
enddo var_loop

return
end subroutine SCAL

!---------------------------------------------------------------------
double precision function SCALA(Y, NJ)

implicit none

integer, intent(in) :: NJ
double precision, intent(in) :: Y(*)

integer :: j
double precision :: YS(6), YMAX

YMAX = MAXVAL(ABS(Y(1: NJ)))

YS(1) = 1.0d-9
YS(2) = 1.5d-9
YS(3) = 2.0d-9
YS(4) = 3.0d-9
YS(5) = 5.0d-9
YS(6) = 8.0d-9

do while(1.05*YMAX > YS(6))
    YS = 10.*YS
enddo
do j=1, 6
    if (1.05*YMAX <= YS(j)) EXIT
enddo

SCALA = YS(J)

end function SCALA

!---------------------------------------------------------------------
subroutine CMARK(xpos_in, ypos_in, prof_yscale, yshift, prof_name, STYL, jplot_parity)
! Mark variable/scale in 1 & 2 modes

use outcmn_inc, only: canv_wid
use dbl2char, only: fmt4

implicit none

integer, intent(in) :: STYL, xpos_in, ypos_in, jplot_parity
double precision, intent(in) :: prof_yscale, yshift
character(len=4), intent(in) :: prof_name

integer :: POINT(2), xpos, ypos, str_len, str_pixels
character(len=15) :: name_scale_label

ypos = ypos_in

name_scale_label(:)  = ' '
name_scale_label(1: 4) = fmt4(prof_yscale)
name_scale_label(6: 10) = prof_name
if (yshift /= 0) then
    if (yshift < 0) then
        name_scale_label(11: 11) = '-'
    else
        name_scale_label(11: 11) = '+'
    endif
    name_scale_label(12: 15) = fmt4(abs(yshift))  
endif
str_len = LEN_TRIM(name_scale_label)
str_pixels = 8*str_len
xpos = xpos_in + 2 + jplot_parity*(canv_wid - str_pixels - 4)
call textvm(xpos, ypos, TRIM(name_scale_label), str_len)
if (STYL > 0) then ! If clicking 'Style' in ASTRA graphic window
    POINT(1) = xpos_in + str_pixels + 8 + jplot_parity*(canv_wid - 2*str_pixels - 16)
    POINT(2) = ypos - 5
    call NMARK(POINT, STYL)
endif

return
end subroutine CMARK

!---------------------------------------------------------------------
subroutine CMARKT(xpos_in, ypos_in, sig_yscale, yshift, sig_name, STYL)
! Mark variable/scale in 6th (time) mode

use outcmn_inc, only: astra_gui
use dbl2char, only: fmt4

implicit none

integer, intent(in) :: xpos_in, ypos_in, STYL
double precision, intent(in) :: sig_yscale, yshift
character(len=4), intent(in) :: sig_name

integer :: str_len, xpos, ypos, plot_arr(2)
character(len=10) :: name_shift_label
character(len=4 ) :: F4

xpos = xpos_in
ypos = ypos_in
name_shift_label(1: 10) = '          '
name_shift_label(2: 5)  = sig_name
if (yshift /= 0) then
    F4 = fmt4(abs(yshift))
    if (yshift < 0) name_shift_label(6: 10) = '-' // F4
    if (yshift > 0) name_shift_label(6: 10) = '+' // F4
    xpos = xpos_in - astra_gui%dxlet
    if (name_shift_label(1: 1) == ' ') xpos = xpos - astra_gui%dxlet
    str_len = LEN_TRIM(name_shift_label(2: 10))
    if (str_len <= 6) xpos = xpos + astra_gui%dxlet
endif

call textvm(xpos, ypos, TRIM(name_shift_label), LEN_TRIM(name_shift_label))   ! type name+yshift
plot_arr(1) = xpos + 3
plot_arr(2) = ypos - 5
F4 = fmt4(sig_yscale)
ypos = ypos + 15
xpos = xpos + astra_gui%dxlet
call textvm(xpos, ypos, F4, 4)   ! type scale
if (STYL > 0) call NMARK(plot_arr, STYL)

return
end subroutine CMARKT

!---------------------------------------------------------------------
subroutine update_curve(NP, np_old, ICOLOR, STYL, xold, yold, xnew, ynew)

! The subroutine displays NP points of the float array YNEW
! NP  is a number of points to plot
! NPO is a number of points to erase
! The points of the array YOLD are used for erasing

use outcmn_inc, only: EraseColor

implicit none

integer, intent(in) :: STYL, NP, np_old, ICOLOR
double precision, intent(in), dimension(np) :: xold, yold, xnew, ynew

if (np_old > 0) then
! erase the old curve
    call colovm(EraseColor)
    call plot_curve(np_old, STYL, xold, yold)
endif

! draw a new curve

call colovm(ICOLOR)
call plot_curve(NP, STYL, xnew, ynew)

return
end subroutine update_curve

!---------------------------------------------------------------------
subroutine plot_curve(np, STYL, xplot, yplot)
  
implicit none

integer, intent(in) :: STYL, np
double precision, intent(in), dimension(*) :: xplot, yplot

integer :: LE, NF, J, j1, JJ, NM, PT1(2)

if (STYL < 0) then  ! Draw dashed curves
    jj = -STYL
    if (jj >= 7) jj = jj + 1 - jj/7*7
    if (jj > 1) then
        LE = 8
        do j=1, jj
            j1 = j + 1
            LE = LE + j1
        enddo
        j1 = jj
        if (jj == 2) LE = min(LE, 16)
        if (jj == 3) LE = min(LE, 8)
        if (jj == 4) LE = min(LE, 4)
        NF = max(1, LE/4)
        do j=1, NP, LE
            j1 = min(NP - j + 1, LE - NF)
            call drawcurve(0, j1, xplot(j: j+j1-1), yplot(j: j+j1-1))
        enddo
        return
    endif
else if (STYL > 0) then
    NM = max(10, NP/5)
    LE = NM/5*STYL    ! 1st marker position
    if (LE >= NM+2) LE = LE - NM
    LE = max(1, LE)
    do jj=LE, NP, NM
        PT1(1) = xplot(jj)
        PT1(2) = yplot(jj)
        call NMARK(PT1, STYL)
    enddo
endif

call drawcurve(0, NP, xplot(1:NP), yplot(1:np))

return
end subroutine plot_curve

!---------------------------------------------------------------------
subroutine NMARK(POINT, STYL)

use outcmn_inc, only: resizeGraph

implicit none

integer, parameter :: n_symbols=7, sym_points=16
integer, parameter, dimension(n_symbols) :: sym_size = (/16, 13, 5, 9, 14, 9, 10/)
integer, dimension(sym_points, n_symbols) :: dx, dy

integer, intent(in) :: POINT(2), STYL

integer :: J, IST
double precision, dimension(sym_points) :: xsym, ysym

! IST definition shoud coincide with NBIT() in CMARK
! IST = 1-filled diamond, 2-o, 3-+, 4-*, 5-x(#), 6-<, 7-filled square
! STYL  <=7,              8    9    10   11      12   >=13

save DX, DY ! Actually parameters
data DX/ &
    0,  3,  0, -3,  0,  0,  2,  0, -2,  0,  0,  1,  0, -1,  0,  0, &
    3,  3,  2,  1, -1, -2, -3, -3, -2, -1,  1,  2,  3,  0,  0,  0, &
    3, -3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
   -3, -1,  2, -2,  1, -2,  2, -1,  3,  0,  0,  0,  0,  0,  0,  0, &
   -2, -1,  0,  1,  2,  1,  0,  1,  2,  1,  0, -1, -2,  0,  0,  0, &
    0,  0, -3,  0, -2,  0,  0, -2,  0,  0,  0,  0,  0,  0,  0,  0, &
   -2,  2,  2, -2, -2,  1,  1, -1, -1,  0,  0,  0,  0,  0,  0,  0/ 
data DY/ &
    3,  0, -3,  0,  3,  2,  0, -2,  0,  2,  1,  0, -1,  0,  1,  0, &
    1, -1, -2, -3, -3, -2, -1,  1,  2,  3,  3,  2,  1,  0,  0,  0, &
    0,  0,  0,  3, -3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
    0,  0,  3, -3,  0,  3, -3,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
   -2, -1,  0,  1,  2,  1,  0, -1, -2, -1,  0,  1,  2,  0,  0,  0, &
    3, -3,  0,  3, -1, -1,  1,  1, -3,  0,  0,  0,  0,  0,  0,  0, &
   -2, -2,  2,  2, -1, -1,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0/  

IST = max(1, min(STYL, 7))
do J=1, sym_size(IST)
    xsym(J) = POINT(1) + NINT(resizeGraph*dx(J, IST))
    ysym(J) = POINT(2) + NINT(resizeGraph*dy(J, IST))
enddo
call drawcurve(0, sym_size(IST), xsym, ysym)

end subroutine NMARK

!---------------------------------------------------------------------
double precision function GETNUM(FIELD, ERCODE)

use outcmn_inc, only: PRNAME, CFNAME
use const_inc , only: CONSTF, DEVARX
use char_manip, only: str_in_list

implicit none

integer, intent(out) :: ERCODE
character(len=*), intent(in) :: FIELD

integer :: l, j, jnam, jpos, jpos1, j1, ISHIFT, ios
character(len=6) :: ZNUM

save ISHIFT
data ISHIFT/0/

if (ISHIFT == 0) then
    j = str_in_list('ZRD1  ', PRNAME)
    ISHIFT = max(j-1, 0)
endif

l  = len(FIELD)
jpos = index(TRIM(FIELD), 'ZRD')

if (jpos == 0) then
    jpos1 = index(TRIM(FIELD), 'C')
    ERCODE = 1 ! Only initial, before CFNAME-search
    if (jpos1 >  0) then
        j1 = min(LEN_TRIM(FIELD(jpos1: l)), l - jpos1 + 1)
        ZNUM = FIELD(jpos1: jpos1+j1-1)
        jnam = str_in_list(ZNUM, CFNAME)
        if (jnam > 0) then
            ERCODE = 0
            GETNUM = CONSTF(jnam)
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
        GETNUM = DEVARX(ISHIFT + j1)
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
! Numbers and references to ZRD*, ZRD*X and CONSTF_list
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
double precision :: GETNUM
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

!---------------------------------------------------------------------
subroutine ASKINT(NV, NVAR, NAME)

use parameter_inc, only: NRW

implicit none

integer, intent(in) :: NV
integer, intent(inout) :: NVAR(*)
character(len=4), intent(in) :: NAME(*)

integer :: J
double precision :: VAR(NRW)

do J=1, NV
    VAR(J) = NVAR(J)
enddo
call MENUTABLE(NV, VAR, NAME, 4)
do J=1, NV
    NVAR(J) = VAR(J)
enddo

return
end subroutine ASKINT

!---------------------------------------------------------------------
subroutine ASTWIN(NB, IBOX, NAME, yscale, yshift, MOD10, YMODE)

use parameter_inc, only: NRW
use outcmn_inc, only: IP1, IP2, IP30, IP31, null_ch
use char_manip, only: to_upper

implicit none

integer, parameter :: NRW16=NRW+16, NRW96=NRW-128

integer, intent(in) :: NB, YMODE, MOD10
integer, intent(out) :: IBOX(*)
double precision, intent(inout), dimension(NB) :: yscale, yshift
character(len=4), intent(out) :: NAME(*)

integer :: JMODE, JGR, j, j1, j2, jj, jn, jm, js, jb, jw, jsep, layoutBox
integer, dimension(NRW) :: IB
character(len=80) :: rows(NRW16), STR, TITLE
character(len=1) :: KEY

TITLE = "Presentation" // null_ch
!              ----5----0----5----0----5----0----5----0----5----0----5
STR = "Name|Box| Scale|Offset||Name|Box| Scale|Offset" // null_ch

jsep = index(STR, '||') + 1
JMODE = MOD10
if (JMODE == 1) then
    JGR = 8
elseif (JMODE == 2 .or. JMODE == 3) then
    JGR = 4
elseif (JMODE == 6) then
    if (YMODE ==  1) JGR = 4
    if (YMODE ==  0) JGR = 4
    if (YMODE == -1) JGR = 2
else
    return
endif

KEY = 'Y'

ask_key: do while(KEY == 'Y')

    do j = 1, NRW
        write(rows(j)(1:80), '(79X, 1A1)') null_ch
    enddo

! j  - ordinal box No.
! jb - box No. in the Astra nominations
! jw - position in the table
! jj - horizontal row
! js - position in the current row

    jn = 0
    jm = 0
    do j = 1, NB
        jb = IBOX(j)
        if (jb <= 0 ) then
            jm = jm + 1
            CYCLE
        endif
        if (JMODE == 1) jw = IP1(jb)
        if (JMODE == 2 .or. JMODE == 3) jw = IP2(jb)
        if (JMODE == 6) then
            if (YMODE == 1) jw = jb
            if (YMODE == 0) jw = IP30(jb)
            if (YMODE  == -1) jw = IP31(jb)
        endif
        js = jsep*(1 - jw + jw/2*2)
        jj = 1 + (jw - 1)/2
        write(rows(jj)(js+1:js+4), '(A4)') NAME(j)
        write(rows(jj)(js+6:js+8), '(I3)') jb
        write(rows(jj)(js+10:js+15), '(f6.1)') yscale(j)
        write(rows(jj)(js+17:js+22), '(f6.1)') yshift(j)
        jn = max(jn, jj)
    enddo

    if (jm /= 0 ) then
        j2 = jn
        j1 = 2*jn + 1
        do j=1, NB
            if (IBOX(j) <= 0 ) then
                jw = j1
                js = jsep*(1 - jw + jw/2*2)
                jj = 1 + (jw - 1)/2
                write(rows(jj)(js+1: js+4), '(A4)') NAME(j)
                write(rows(jj)(js+6: js+8), '(I3)') IBOX(j)
                write(rows(jj)(js+10:js+15), '(f6.1)') yscale(j)
                write(rows(jj)(js+17:js+22), '(f6.1)') yshift(j)
                jn = max(jn, jj)
                j1 = j1 + 1
            endif
        enddo
    endif

    j = 1
    do while(j > 0)
        j1 = (jm + 1)/2
        j = layoutBox(TITLE, STR, rows, 80, jn, JGR, j1)
    enddo

    jn = 0
    do j=1, NB
        if (IBOX(j) >  0 ) then
            if (JMODE == 1) jw = IP1(IBOX(j))
            if (JMODE == 2 .or. JMODE == 3) jw = IP2(IBOX(j))
            if (JMODE == 6) then
                if (YMODE  == 1) jw = IBOX(j)
                if (YMODE  == 0) jw = IP30(IBOX(j))
                if (YMODE == -1) jw = IP31(IBOX(j))
            endif
            js = jsep*(1 - jw + jw/2*2)
            jj = 1 + (jw - 1)/2
            write(NAME(j), '(A4)', ERR=77) rows(jj)(js+1: js+4)
            read(rows(jj)(js+6 : js+8 ), *, ERR=77) IB(j)
            read(rows(jj)(js+10: js+15), *, ERR=77) yscale(j)
            read(rows(jj)(js+17: js+22), *, ERR=77) yshift(j)
            jn = max(jn, jj)
        endif
    enddo

    if (jm /= 0 ) then
        j1 = 2*jn + 1
        do j=1, NB
            if (IBOX(j) <= 0 ) then
                jw = j1
                js = jsep*(1 - jw + jw/2*2)
                jj = 1 + (jw - 1)/2
                write(NAME(j), '(A4)', ERR=77) rows(jj)(js+1: js+4)
                read(rows(jj)(js+6: js+8), *, ERR=77) IB(j)
                read(rows(jj)(js+10:js+15), *, ERR=77) yscale(j)
                read(rows(jj)(js+17:js+22), *, ERR=77) yshift(j)
                jn = max(jn, jj)
                j1 = j1 + 1
            endif
        enddo
    endif

    do j = 1, NB
        IBOX(j) = IB(j)
    enddo

    return

77 continue

    write(*, *)
    write(*, '(A)') '>>> INPUT ERROR encountered in the dialog window "Presentation"'
    write(*, '(A23, A52, A1)') '                Line: "', rows(j)(1: 52), '"'
    write(*, '(A52, $)') '     Enter "Y" to return, any other key to ignore > '
    KEY = 'X'
    read(*, '(:, A1)') KEY
    KEY = to_upper(KEY)

enddo ask_key

return
end subroutine ASTWIN

!---------------------------------------------------------------------
subroutine ASXWIN(NB, IBOX, NAME, yscale, yshift, r_min, r_max, MOD10, YMODE)

use parameter_inc, only: NRW
use const_inc, only: AB
use outcmn_inc, only: IP1, IP2, IP30, IP31, null_ch

implicit none

integer, parameter :: NRW16=NRW+16, NRW96=NRW-128

integer, intent(in) :: NB, YMODE, MOD10
integer, intent(out) :: IBOX(*)
double precision, intent(inout), dimension(NB) :: yscale, yshift, r_min, r_max
character(len=4), intent(out) :: NAME(*)

integer :: JMODE, JGR, j, j1, j2, jj, jn, jm, js, jb, jw, jsep, layoutBox
integer, dimension(NRW) :: IB
double precision :: YY
character(len=80) :: rows(NRW16), STR, TITLE
character(len=1) :: YKEY
character(len=6) :: ABNUM

TITLE = "Presentation" // null_ch
!              ----5----0----5----0----5----0----5----0----5----0----5
STR = "Name|Box| Scale|Offset| [a?, ]| [, a?]||" // &
      "Name|Box| Scale|Offset| [a?, ]| [, a?]"   // null_ch
jsep = index(STR, '||')+1

JMODE = MOD10
if (JMODE == 1) then
    JGR = 8
elseif (JMODE == 2 .or. JMODE == 3) then
    JGR = 4
elseif (JMODE == 6) then
    if (YMODE ==  1) JGR = 4
    if (YMODE ==  0) JGR = 4
    if (YMODE == -1) JGR = 2
else
    return
endif

10 continue

do j=1, NRW
    write(rows(j)(1:80), '(79X, 1A1)') null_ch
enddo

! j  - ordinal box No.
! jb - box No. in the Astra nominations
! jw - position in the table
! jj - horizontal row
! js - position in the current row

jn = 0
jm = 0
write(ABNUM, '(f6.3)') ABNUM
do j=1, NB
    jb = IBOX(j)
    if (jb <= 0 ) then
        jm = jm + 1
        CYCLE
    endif
    if (JMODE == 1) jw = IP1(jb)
    if (JMODE == 2 .or. JMODE == 3) jw = IP2(jb)
    if (JMODE == 6) then
        if (YMODE == 1) jw = jb
        if (YMODE == 0) jw = IP30(jb)
        if (YMODE  == -1) jw = IP31(jb)
    endif
    YY = r_max(j)
    js = jsep*(1 - jw + jw/2*2)
    jj = 1 + (jw - 1)/2
    write(rows(jj)(js+1: js+4), '(A4)') NAME(j)
    write(rows(jj)(js+6: js+8), '(I3)') jb
    write(rows(jj)(js+10: js+15), '(f6.1)') yscale(j)
    write(rows(jj)(js+17: js+22), '(f6.1)') yshift(j)
    write(rows(jj)(js+24: js+29), '(f6.3)') r_min(j)
    write(rows(jj)(js+31: js+36), '(f6.3)') r_max(j)
    if (YY > AB) rows(jj)(js+31:js+36) = ABNUM
    jn = max(jn, jj)
enddo

if (jm /= 0 ) then
    j2 = jn
    j1 = 2*jn+1
    do j = 1, NB
        if (IBOX(j) <= 0 ) then
            YY = r_max(j)
            jw = j1
            js = jsep*(1 - jw + jw/2*2)
            jj = 1 + (jw - 1)/2
            write(rows(jj)(js+1: js+4), '(A4)') NAME(j)
            write(rows(jj)(js+6: js+8), '(I3)') IBOX(j)
            write(rows(jj)(js+10: js+15), '(f6.1)') yscale(j)
            write(rows(jj)(js+17: js+22), '(f6.1)') yshift(j)
            write(rows(jj)(js+24: js+29), '(f6.3)') r_min(j)
            write(rows(jj)(js+31: js+36), '(f6.3)') r_max(j)
            if (YY > AB) rows(jj)(js+31: js+36) = ABNUM
            jn = max(jn, jj)
            j1 = j1 + 1
        endif
    enddo
endif

j = 1
do while(j > 0)
    j1 = (jm + 1)/2
    j = layoutBox(TITLE, STR, rows, 80, jn, JGR, j1)
enddo

jn = 0
do j = 1, NB
    if (IBOX(j) > 0) then
        if (JMODE == 1) jw = IP1(IBOX(j))
        if (JMODE == 2 .or. JMODE == 3) jw = IP2(IBOX(j))
        if (JMODE == 6) then
            if (YMODE ==  1) jw = IBOX(j)
            if (YMODE ==  0) jw = IP30(IBOX(j))
            if (YMODE == -1) jw = IP31(IBOX(j))
        endif
        js = jsep*(1 - jw + jw/2*2)
        jj = 1 + (jw - 1)/2
        write(NAME(j), '(1A4)', ERR=77) rows(jj)(js+1: js+4)
        read(rows(jj)(js+6 : js+8), *, ERR=77) IB(j)
        read(rows(jj)(js+10: js+15), *, ERR=77) yscale(j)
        read(rows(jj)(js+17: js+22), *, ERR=77) yshift(j)
        read(rows(jj)(js+24: js+29), *, ERR=77) r_min(j)
        read(rows(jj)(js+31: js+36), *, ERR=77) YY
        if (rows(jj)(js+31: js+36) /= ABNUM) then
            read(rows(jj)(js+31: js+36), *, ERR=77) r_max(j)
        endif
        jn = max(jn, jj)
    endif
enddo

if (jm /= 0 ) then
    j1 = 2*jn + 1
    do j=1, NB
        if (IBOX(j) <= 0) then
            jw = j1
            js = jsep*(1 - jw + jw/2*2)
            jj = 1 + (jw - 1)/2
            write(NAME(j), '(1A4)', ERR=77) rows(jj)(js+1: js+4)
            read(rows(jj)(js+6 : js+8), *, ERR=77) IB(j)
            read(rows(jj)(js+10: js+15), *, ERR=77) yscale(j)
            read(rows(jj)(js+17: js+22), *, ERR=77) yshift(j)
            read(rows(jj)(js+24: js+29), *, ERR=77) r_min(j)
            if  (rows(jj)(js+31: js+36) /= ABNUM) then
                read(rows(jj)(js+31:js+36), *, ERR=77) r_max(j)
            endif
            jn = max(jn, jj)
            j1 = j1 + 1
        endif
    enddo
endif

do j = 1, NB
    IBOX(j) = IB(j)
enddo

return

77 continue
write(*, *)
write(*, '(A)') '>>> INPUT ERROR encountered in the dialog window "Presentation"'
write(*, '(A, A52, A1)') '                Line: "', rows(j)(1:52), '"'
write(*, '(A52, $)') '     Enter "Y" to return, any other key to ignore > '
YKEY = 'X'
read(*, '(:, A1)')YKEY
if (YKEY == 'Y' .or. YKEY == 'y') goto 10

return
end subroutine ASXWIN

!---------------------------------------------------------------------
subroutine get_runid
!---------------------------------------------------------------------
! The subroutine forms string RUNID and additionally returns 
! date and time when those are not defined (calling from INIT)
!---------------------------------------------------------------------

use outcmn_inc, only: RUNID, equ_file, exp_file, VERSION

implicit none

integer :: time_arr(8), YEAR, MONTH, DAY, HOUR, MINUTE
integer :: j
character(len=3) :: vers
character(len=15) :: datetime

call date_and_time(VALUES=time_arr)

YEAR   = time_arr(1)
MONTH  = time_arr(2)
DAY    = time_arr(3)
HOUR   = time_arr(5)
MINUTE = time_arr(6)
write(datetime, "(1I2, 2('-', 1I2.2), 1I3, ':', 1I2.2)") &
    DAY, MONTH, YEAR-2000, HOUR, MINUTE

j = index(VERSION, 'Version')
vers = version(j+8: j+10)

RUNID = "ASTRA " // vers // " -- " // datetime // ' -- Model: ' // &
    TRIM(equ_file) // ' -- Data: ' // TRIM(exp_file)

return
end subroutine get_runid
