module nb_injection

implicit none
contains

!---------------------------------------------------------------------
    subroutine NBINJ(file_nbi, BTOR, RTOR, ABC, AB, ROC, SHIFT, UPDWN, &
        HRO, TAU, NA1, NB1, AIM1, AIM2, AIM3, AMJ, ZMJ, NNCL, NNWM, QNBI, &
        CBM1, CBM2, CBMI3, CBMI1, CBM4, CBMI2, CBM3, CBMI4)

    use nbstatus, only: n_rho, n_fields, n_theta, ISPE, ISPEND, YRIPLR, &
        AMETR, PBEAM, SCUBM, SNNBM, SNEBM, SNIBM1, SNIBM2, SNIBM3, &
        NNBM1, NNBM2, NNBM3, NIBM, PEBM, PIBM, CUBM, CUFI, PBPER, PBLON, &
        NE, NI, NHYDR, NDEUT, NTRIT, NHE3, NALF, NIZ1, NIZ2, NIZ3, &
        TI, ZIM1, ZIM2, ZIM3
    use nbibce, only: nbionr
    use nbicom, only: nbspec, nbsrsr, nbion0, stnbdp, sdnbtp, sdnbdp1, sdnbdp2, RMB, ZB

    integer, parameter :: n_unit=36

    integer, intent(in) :: NA1, NB1
    double precision, intent(in) :: BTOR, RTOR, ABC, AB, ROC, SHIFT, UPDWN, HRO, TAU, &
        AIM1, AIM2, AIM3, AMJ, ZMJ, NNCL, NNWM, CBM2, CBM3, CBM4, CBMI1, CBMI2, CBMI4
    character(len=*), intent(in) :: file_nbi
    double precision, intent(out) :: QNBI
    double precision, intent(inout) :: CBM1, CBMI3

    logical :: file_exists
    integer :: N, j_nbi, J, JN, JN1, JNAX, ERCODE, JWARN, &
        JSRREC, IFLAG, n_nbi, jElevAll
    double precision :: CBMH1, CBMH2, CBMR1, CBMR2, CBMS1, CBMS2, CBMS3, CBMS4, EBEAM, &
        DBM1=1.d0, DBM2=0.d0, DBM3=0.d0, ABEAM=1.d0, HBEAM, CONTR, RBMAX=0.5d0, RBMIN=1.d0, QBEAM
    double precision :: YQBEAM, ZBEAM, YABEAM, YEBEAM, YUD, YHM
    double precision :: ARRAY(n_fields), yEXTARR(n_rho, 9)

!---- Ripple normalized radius of the ripple boundary
!---- banana with RTCRIT > YRIPLR is lost
    double precision :: Y, Y1, Y2, YMAX, YQBM
    real*8, allocatable :: YQBMJE(:, :)

    jElevAll = 0 ! for write to SSFPQL
    ZBEAM = 1.d0 ! Hydrogen isotopes only

    inquire(file=TRIM(file_nbi), exist=file_exists)

    if (.not. file_exists) write(*, *) 'File "', TRIM(file_nbi), '" not found'
    n_nbi = abs(int(CBM1 + 0.5d0))
    if (CBM1 < 0.0) n_nbi = n_nbi + 1
    if (n_nbi == 0) then
        write(*, *) '>>> NBI >>> Zero number of sources'
        write(*, *)"            Don't know what to do"
        stop
    endif

    if (CBMI1 == 0.d0) then ! print for SSFPQL
        if (.not. allocated(YQBMJE)) allocate (YQBMJE(n_nbi, 3))
        YQBMJE = 0.d0
    endif ! print for SSFPQL

! CBMI3 to make NBI internal mesh interval > max Larmor raduis
    Y = 1.d-3
    EBEAM = 1.d-3

    open(2, file=TRIM(file_nbi), status='OLD')
    do j=1, n_nbi
        call STREAD(2, 20, ARRAY, ERCODE)
        if (j == 1) then
            ABEAM = ARRAY(3)
            EBEAM = ARRAY(5)
        else
            if (EBEAM < ARRAY(5)) EBEAM = ARRAY(5)
            if (CBMI1 /= 1.d0 .and. ABEAM /= ARRAY(3)) write(*, *) 'ABEAM must be the same for all NBIs if CBM4  /=  1'
        endif
        if (Y < (ARRAY(3)*ARRAY(5))) Y = ARRAY(3)*ARRAY(5)
    enddo
    close(2)

    if (CBMI1 > 0.d0) then
        CBMI3 = 1.d0
        J = 5.d-3*dsqrt(Y)*(NA1/ABC)/(BTOR*RTOR/(RTOR + SHIFT))
        if (J > 1 .and. J < (NA1-1)) CBMI3 = J
        if (5*J > (NA1-1) .and. NA1 > 6) CBMI3 = (NA1 - 1)/5
    else  ! w/o FP solver FI source on the transport grid
        CBMI3 = 1.d0
    endif

    QNBI = 0.d0
    CBM1 = n_nbi

!---- Ripple

    YUD = UPDWN
    Y = (AMETR(NA1) + AMETR(NA1 - 1))/2.d0
    N = (NA1 - 1)/CBMI3
    Y1 = 0.
    do JN1=N+1, 2, -1
        JN = JN1 - 1
        JNAX = 1 + CBMI3*(JN - 1)
        Y2 = RIPRAD(YUD, JNAX)/Y
        if (Y2 >= Y1) YMAX = Y2
        Y1 = Y2
        YRIPLR(JN) = YMAX
    enddo
    YRIPLR(N+1) = YRIPLR(N)

! Hot ion' source
    do JN=1, NB1
        PBEAM(JN)   = 0.
        SCUBM(JN)   = 0.
        SNNBM(JN)   = 0.
        SNEBM(JN)   = 0.
        NNBM1(JN)   = 0.
        NNBM2(JN)   = 0.
        NNBM3(JN)   = 0.
        stnbdp(JN)  = 0.
        sdnbtp(JN)  = 0.
        sdnbdp1(JN) = 0.
        sdnbdp2(JN) = 0.
        SNIBM1(jn)  = 0.
        SNIBM2(jn)  = 0.
        SNIBM3(jn)  = 0.
        PEBM(JN)    = 0.
        PIBM(JN)    = 0.
        CUBM(JN)    = 0.
        CUFI(JN)    = 0.
        PBPER(JN)   = 0.
        NIBM(JN)    = 0.
        PBLON(JN)   = 0.
    enddo
    ISPE = 0
    yEXTARR = 0.d0

! Set plasma composition
    call NBSPEC(AIM1, AIM2, AIM3, AMJ, ZMJ, NA1, &
        yEXTARR, NE, NI, NHYDR, NDEUT, NTRIT, NHE3, NALF, NIZ1, NIZ2, NIZ3, &
        ZIM1, ZIM2, ZIM3, RMB, ZB, ISPEND, ISPE, IFLAG)

    if (IFLAG /= 0) return

    if (n_nbi == 0) then
        write(*, *) ' -1 < CBM1 < 1 = no NBI sources '
        return
    endif

    JWARN = 0
    jElevAll = 0

    YEBEAM = EBEAM
    YABEAM = ABEAM
    YQBM = 0.d0

    open(2, file=TRIM(file_nbi), status='OLD')

    do j_nbi=1, n_nbi
        call STREAD(2, 20, ARRAY, ERCODE)
        QBEAM = ARRAY(1)
        CONTR = ARRAY(2)
        ABEAM = ARRAY(3)
        EBEAM = ARRAY(5)
        DBM1  = ARRAY(6)
        DBM2  = ARRAY(7)
        DBM3  = ARRAY(8)
        CBMS1 = ARRAY(9)
        CBMS2 = ARRAY(10)
        HBEAM = ARRAY(11)
        RBMAX = ARRAY(12)
        RBMIN = ARRAY(13)
        CBMS3 = ARRAY(14)
        CBMS4 = ARRAY(15)
        CBMH1 = ARRAY(16)
        CBMH2 = ARRAY(17)
        CBMR1 = ARRAY(18)
        CBMR2 = ARRAY(19)
        if (CBMI1 == 0.d0) then ! print for SSFPQL
            if (QBEAM*EBEAM > 0d0) then ! NBI on
                if (DBM1 > 0.d0) then
                    YQBMJE(j_nbi, 1) = QBEAM*DBM1
                    jElevAll = jElevAll + 1
                endif
                if (DBM2 > 0.d0) then
                    YQBMJE(j_nbi, 2) = QBEAM*DBM2
                    jElevAll = jElevAll + 1
                endif
                if (DBM3 > 0.d0) then
                    YQBMJE(j_nbi, 3) = QBEAM*DBM3
                    jElevAll = jElevAll + 1
                endif
            endif
        endif

        if (YABEAM /= ABEAM) JWARN = 1
        QNBI = QNBI + QBEAM
        YQBEAM = QBEAM

! Switch off the ion source
        YQBM = YQBM + YQBEAM
        YHM = (HBEAM - YUD)*100.
        if (CONTR == 0.d0 .or. CONTR == 1.d0) then
            call NBSRSR(j_nbi, 2.d0*CONTR - 1.d0, NA1, RTOR, SHIFT, AB, BTOR, &
                HRO, YHM, CBMH1, CBMH2, CBMS1, CBMS2, CBMS3, CBMS4, &
                CBMR1, CBMR2, CBMI3, CBMI1, EBEAM, DBM1, DBM2, DBM3, ABEAM, &
                QBEAM, RBMAX, RBMIN, JSRREC, YEXTARR, CBMI4)
            if (CBMI1 == 1.d0) call NBION0(NA1, ABEAM, EBEAM, RTOR, CBMI3, yEXTARR)
        else ! balanced injection
! coinj. part
            QBEAM = YQBEAM*(1. - CONTR)
            call NBSRSR(j_nbi, -1.d0, NA1, RTOR, SHIFT, AB, BTOR, HRO, YHM, &
                CBMH1, CBMH2, CBMS1, CBMS2, CBMS3, CBMS4, CBMR1, CBMR2, CBMI3, CBMI1, &
                EBEAM, DBM1, DBM2, DBM3, ABEAM, QBEAM, RBMAX, RBMIN, JSRREC, &
                YEXTARR, CBMI4)
            if (CBMI1 == 1.d0) call NBION0(NA1, ABEAM, EBEAM, RTOR, CBMI3, yEXTARR)

!: counter injection part...
            QBEAM = YQBEAM*CONTR
            call NBSRSR(j_nbi, 1.d0, NA1, RTOR, SHIFT, AB, BTOR, HRO, YHM, &
                CBMH1, CBMH2, CBMS1, CBMS2, CBMS3, CBMS4, CBMR1, CBMR2, CBMI3, CBMI1, &
                EBEAM, DBM1, DBM2, DBM3, ABEAM, QBEAM, RBMAX, RBMIN, JSRREC, &
                YEXTARR, CBMI4)
            if (CBMI1 == 1.d0) call NBION0(NA1, ABEAM, EBEAM, RTOR, CBMI3, yEXTARR)
            QBEAM = YQBEAM
        endif
    enddo

    close(2)
    close(39)

    EBEAM = YEBEAM
    QBEAM = YQBM

! Fast ion' distribution (3D + t)
    if (CBMI1 >= 2.d0)  then
        if (JWARN == 1) then
            write(*, *) 'ABEAM must be the same for CBMI1#1'
            ABEAM = YABEAM
        endif

        call NBIONR(EBEAM, ABEAM, RTOR, NA1, TAU, NNCL, NNWM, CBM1, CBM3, CBM4, &
            CBMI2, CBMI3, JSRREC, YEXTARR)
    endif

! Conversion to rough mesh keeping the intagrals
    if (CBMI3 /= 1.d0)    then
        call smooth_int(NIBM, CBMI3, ROC, NA1)
        call smooth_int(PIBM, CBMI3, ROC, NA1)
        call smooth_int(PEBM, CBMI3, ROC, NA1)
        call smooth_int(PBLON, CBMI3, ROC, NA1)
        call smooth_int(PBPER, CBMI3, ROC, NA1)
        call smooth_int(PBEAM, CBMI3, ROC, NA1)
        call smooth_int(SNNBM, CBMI3, ROC, NA1)
        call smooth_int(SNEBM, CBMI3, ROC, NA1)
        if (ABEAM < 1.5d0) then
            call smooth_int(SNIBM1, CBMI3, ROC, NA1) !H ion source
        else if (ABEAM > 2.5d0) then
            call smooth_int(SNIBM3, CBMI3, ROC, NA1) !T ion source
        else
            call smooth_int(SNIBM2, CBMI3, ROC, NA1) !D ion source
        endif
        call smooth_int(SCUBM  , CBMI3, ROC, NA1)
        call smooth_int(stnbdp , CBMI3, ROC, NA1)
        call smooth_int(sdnbdp1, CBMI3, ROC, NA1)
        call smooth_int(sdnbdp2, CBMI3, ROC, NA1)
        call smooth_int(sdnbtp , CBMI3, ROC, NA1)
        call smooth_int(CUFI, CBMI3, ROC, NA1, surf=.true.)
        call smooth_int(CUBM, CBMI3, ROC, NA1, surf=.true.)
    endif

    do j=1, NA1
        PEBM(J) = PEBM(J) - 2.08E-5*SNEBM(J)
        if (CBM2 > 0.)  then
            PIBM(J) = PIBM(J) - 0.0024*SNNBM(J)*TI(J)
        endif
    enddo

    do j=1, na1
        NNBM2(j)  = SNIBM2(j)
        NNBM3(j)  = stnbdp(j)
        SNIBM1(j) = sdnbdp1(j)
        SNIBM2(j) = sdnbdp2(j)
        SNIBM3(j) = sdnbtp(j)
    enddo

    return

990  write(*, *) '>>> NBI >>> read error in dat/srsfi.dat'
    stop

    end subroutine NBINJ

!---------------------------------------------------------------------
    subroutine smooth_int(YFO, YCI3, YROC, JNA1, surf)
!---------------------------------------------------------------------
!  The subroutine transmits a histogram like function YFO(1:N1)
!  from an N1 grid X(1:N1) to a smooth functon YFO(1:JNA1)
!  keeping the same total volume integral.
! Input: YFO(1:N1)
! Output: YFO(1:JNA1)
!---------------------------------------------------------------------

    use numerical_tools, only: smooth
    use standard_functions, only: VINT, IINT
    use nbicom, only: N1, X, XJ, DRI

    integer, intent(in) :: JNA1
    double precision, intent(in) :: YCI3, YROC
    double precision, intent(inout) :: YFO(*)
    logical, intent(in), optional :: surf

    logical :: surf_int
    integer :: j, JNA, JNAC, JSIGN
    double precision :: Y, YOLD, YNEW, ALFA

    if (present(surf)) then
        surf_int = surf
    else
        surf_int = .false.
    endif
 
    JSIGN = 0
    do j=1, n1-1
        JNA = 1 + YCI3*(j - 1)
        JNAC = JNA - 1 + YCI3
        X(j) = XJ(JNAC)
        DRI(j) = YFO(JNA)
        if (JSIGN < 1) then
            if (DRI(1)*DRI(J) < 0.d0) JSIGN = 1
        endif
    enddo

    DRI(N1)  = 0.d0
    X(N1)    = 1.d0
    XJ(JNA1) = 1.d0

! Total power normalization

    if (surf_int) then
        YOLD = IINT(YFO, YROC)
    else
        YOLD = VINT(YFO, YROC)
    endif

    ALFA = 0.001d0
    call SMOOTH(ALFA, x, dri, n1, xj, yfo, jna1)

! Cut of artificial negatives/positive after smoothing
    if (JSIGN == 0) then !no real change of sign
        do J=1, JNA1
            if (YOLD > 0.d0) then
                YFO(J) = max(YFO(J), 0.d0)
            else
                YFO(J) = min(YFO(J), 0.d0)
            endif
        enddo
    endif

    if (surf_int) then
        YNEW = IINT(YFO, YROC)
    else
        YNEW = VINT(YFO, YROC)
    endif

    if (dabs(YNEW) > 1.d-19) then
        Y = YOLD/YNEW
        do J=1, JNA1
            YFO(J) = Y*YFO(J)
        enddo
    endif

    end subroutine smooth_int

!---------------------------------------------------------------------
    double precision function RIPRAD(YUPDWN, J)
! Ripple loss boundary R[m] for each magnetic surface
! YUPDWN [m] plasma midplane shift in respect to the plane
!            of the ripple losses simmetry
! J - magnetic surface number

    integer, intent(in) :: J
    double precision, intent(in) :: YUPDWN

! Default (No ripple losses)
    RIPRAD = 99999.

    end function RIPRAD

!---------------------------------------------------------------------
    double precision function GETNUM(FIELD, ERCODE)

    use scalars, only: constValues, varxValues
    use char_manip, only: str_in_list
    use json_vars, only: constNames, varNames

    character(len=*), intent(in) :: FIELD
    integer, intent(out) :: ERCODE

    integer :: l, j, jnam, jpos, jpos1, j1, ISHIFT, ios
    character(len=6) :: ZNUM

    save ISHIFT
    data ISHIFT/0/

    if (ISHIFT == 0) then
        j = str_in_list('ZRD1  ', varNames)
        ISHIFT = max(j-1, 0)
    endif

    l = len(FIELD)
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
            else
                ZNUM = FIELD(jpos+3: jpos+j1+1)
            endif
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

    integer, intent(in) :: NCH, NFIELD
    integer, intent(out) :: ERCODE
    double precision, intent(out) :: ARRAY(NFIELD)

    integer :: j, ios
    character(len=12 ) :: SFIELD(20)
    character(len=132) :: str_line
!---------------------------------------------------------------------

    if (NFIELD > 20) then
        ERCODE = 4
    else
! Skip all lines beginning with '!'
        j = 1
        do while (j == 1)
            read(NCH, '(A)', iostat=ios) str_line
            if (ios < 0) then ! EOF encountered
                ERCODE = 5
                EXIT
            else if (ios > 0) then
                ERCODE = 1
                EXIT
            endif
            j = index(str_line, '!')
        enddo

        if (j /= 0 .and. ERCODE == 0) then
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
    endif 

    if (ERCODE > 0) then
        SELECT CASE(ERCODE)
        CASE(1)
            write(*, *) '>>> NBI >>> Error in NBI file: unrecognized variable name'
        CASE(2)
            write(*, *) '>>> NBI >>> NBI File read error'
        CASE(3)
            write(*, *) '>>> NBI >>> Wrong NBI configuration file format. '
            write(*, *) '            More records expected than available.'
        CASE(4)
            write(*, *) '>>> NBI calling STREAD: array out of limits'
        CASE(5)
            write(*, *)'>>> NBI >>> Wrong NBI configuration file format: '
        END SELECT
        stop
    endif

    end subroutine STREAD

end module nb_injection
