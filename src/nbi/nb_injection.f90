module nb_injection

implicit none
contains

!---------------------------------------------------------------------
    subroutine NBINJ(file_nbi, BTOR, RTOR, ABC, AB, ROC, SHIFT, UPDWN, &
        HRO, TAU, NA1, NB1, AIM1, AIM2, AIM3, AMJ, ZMJ, NNCL, NNWM, QNBI, &
        n_nbi, cx_flag, dn_rho, fp_flag, cx_cold, CBMI2, CBM3, calc_fus)

    use nbstatus, only: n_rho, n_fields, n_theta, n_energy, ISPE, ISPEND, YRIPLR, &
        AMETR, PBEAM, SCUBM, SNNBM, SNEBM, SNIBM1, SNIBM2, SNIBM3, &
        NNBM1, NNBM2, NNBM3, NIBM, PEBM, PIBM, CUBM, CUFI, PBPER, PBLON, &
        NE, NI, NHYDR, NDEUT, NTRIT, NHE3, NALF, NIZ1, NIZ2, NIZ3, &
        TI, ZIM1, ZIM2, ZIM3
    use nbibce, only: nbionr
    use nbicom, only: nbsrsr, nbion0, stnbdp, sdnbtp, sdnbdp1, sdnbdp2, RMB, ZB

    integer, parameter :: n_unit=36

    integer, intent(in) :: NA1, NB1, cx_flag, fp_flag, cx_cold, calc_fus
    double precision, intent(in) :: BTOR, RTOR, ABC, AB, ROC, SHIFT, UPDWN, HRO, TAU, &
        AIM1, AIM2, AIM3, AMJ, ZMJ, NNCL, NNWM, CBM3, CBMI2
    character(len=*), intent(in) :: file_nbi
    double precision, intent(out) :: QNBI
    integer, intent(inout) :: n_nbi, dn_rho

    logical :: file_exists
    integer :: N, j_nbi, J, JN, JN1, JNAX, ERCODE, JWARN, j_rho, JSRREC, IFLAG, jE_min
    double precision :: CBMH1, CBMH2, CBMR1, CBMR2, CBMS1, CBMS2, CBMS3, CBMS4, EBEAM, &
        power_frac(n_energy), ABEAM=1.d0, HBEAM, CONTR, RBMAX=0.5d0, RBMIN=1.d0, QBEAM
    double precision :: YQBEAM, YABEAM, YEBEAM, YUD, YHM
    double precision :: ARRAY(n_fields), yEXTARR(n_rho, 9)

!---- Ripple normalized radius of the ripple boundary
!---- banana with RTCRIT > YRIPLR is lost
    double precision :: Y, Y1, Y2, YMAX, YQBM

    power_frac(1) = 1.d0
    inquire(file=TRIM(file_nbi), exist=file_exists)

    if (.not. file_exists) write(*, *) 'File "', TRIM(file_nbi), '" not found'
    if (n_nbi < 0) n_nbi = n_nbi + 1
    if (n_nbi == 0) then
        write(*, *) '>>> NBI >>> Zero sources, stopping'
        stop
    endif

! dn_rho to make NBI internal mesh interval > max Larmor raduis
    Y = 1.d-3
    EBEAM = 1.d-3

    open(2, file=TRIM(file_nbi), status='OLD')
    do j_nbi=1, n_nbi
        call STREAD(2, 20, ARRAY, ERCODE)
        if (j_nbi == 1) then
            ABEAM = ARRAY(3)
            EBEAM = ARRAY(5)
        else
            if (EBEAM < ARRAY(5)) EBEAM = ARRAY(5)
            if (fp_flag /= 1 .and. ABEAM /= ARRAY(3)) write(*, *) 'ABEAM must be the same for all NBIs if cx_cold  /=  1'
        endif
        Y = max(Y, ARRAY(3)*ARRAY(5))
    enddo
    close(2)

    dn_rho = 1
    if (fp_flag > 0) then
        J = 5.d-3*sqrt(Y)*(NA1/ABC)/(BTOR*RTOR/(RTOR + SHIFT))
        if (J > 1 .and. J < (NA1-1)) dn_rho = J
        if (5*J > (NA1-1) .and. NA1 > 6) dn_rho = (NA1 - 1)/5
    endif

    QNBI = 0.d0

!---- Ripple
    YUD = UPDWN
    Y = (AMETR(NA1) + AMETR(NA1 - 1))/2.d0
    N = (NA1 - 1)/dn_rho
    Y1 = 0.
    do JN1=N+1, 2, -1
        JN = JN1 - 1
        JNAX = 1 + dn_rho*(JN - 1)
        Y2 = RIPRAD(YUD, JNAX)/Y
        if (Y2 >= Y1) YMAX = Y2
        Y1 = Y2
        YRIPLR(JN) = YMAX
    enddo
    YRIPLR(N+1) = YRIPLR(N)

! Hot ion' source
    do j_rho=1, NB1
        PBEAM(j_rho)   = 0.
        SCUBM(j_rho)   = 0.
        SNNBM(j_rho)   = 0.
        SNEBM(j_rho)   = 0.
        NNBM1(j_rho)   = 0.
        NNBM2(j_rho)   = 0.
        NNBM3(j_rho)   = 0.
        stnbdp(j_rho)  = 0.
        sdnbtp(j_rho)  = 0.
        sdnbdp1(j_rho) = 0.
        sdnbdp2(j_rho) = 0.
        SNIBM1(j_rho)  = 0.
        SNIBM2(j_rho)  = 0.
        SNIBM3(j_rho)  = 0.
        PEBM(j_rho)    = 0.
        PIBM(j_rho)    = 0.
        CUBM(j_rho)    = 0.
        CUFI(j_rho)    = 0.
        PBPER(j_rho)   = 0.
        NIBM(j_rho)    = 0.
        PBLON(j_rho)   = 0.
    enddo
    ISPE = 0
    yEXTARR = 0.d0

! Set plasma composition
    call NBSPEC(AIM1, AIM2, AIM3, AMJ, ZMJ, NA1, &
        yEXTARR, NE, NI, NHYDR, NDEUT, NTRIT, NHE3, NALF, NIZ1, NIZ2, NIZ3, &
        ZIM1, ZIM2, ZIM3, RMB, ZB, ISPEND, ISPE, IFLAG)

    if (IFLAG /= 0) return

    if (n_nbi == 0) then
        write(*, *) 'No NBI sources '
        return
    endif

    JWARN = 0

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
        power_frac = ARRAY(6: 8)
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

        if (YABEAM /= ABEAM) JWARN = 1
        QNBI = QNBI + QBEAM
        YQBEAM = QBEAM

! Switch off the ion source
        YQBM = YQBM + YQBEAM
        YHM = (HBEAM - YUD)*1.d2

        if (CONTR /= 1.d0) then ! coinjection part
            QBEAM = YQBEAM*(1. - CONTR)
            call NBSRSR(j_nbi, -1.d0, NA1, RTOR, SHIFT, AB, BTOR, HRO, YHM, &
                CBMH1, CBMH2, CBMS1, CBMS2, CBMS3, CBMS4, CBMR1, CBMR2, dn_rho, fp_flag, &
                EBEAM, power_frac, ABEAM, QBEAM, RBMAX, RBMIN, JSRREC, YEXTARR, calc_fus, jE_min)
        endif

        if (CONTR /= 0.d0) then ! counter injection part...
            QBEAM = YQBEAM*CONTR
            call NBSRSR(j_nbi, 1.d0, NA1, RTOR, SHIFT, AB, BTOR, HRO, YHM, &
                CBMH1, CBMH2, CBMS1, CBMS2, CBMS3, CBMS4, CBMR1, CBMR2, dn_rho, fp_flag, &
                EBEAM, power_frac, ABEAM, QBEAM, RBMAX, RBMIN, JSRREC, YEXTARR, calc_fus, jE_min)
        endif

        if (fp_flag == 1) call NBION0(jE_min, NA1, ABEAM, EBEAM, RTOR, dn_rho, yEXTARR)

    enddo

    close(2)
    close(39)

    EBEAM = YEBEAM
    QBEAM = YQBM

! Fast ion' distribution (3D + t)
    if (fp_flag > 1)  then
        if (JWARN == 1) then
            write(*, *) 'ABEAM must be the same for fp_flag=2'
            ABEAM = YABEAM
        endif

        call NBIONR(EBEAM, ABEAM, RTOR, NA1, TAU, NNCL, NNWM, n_nbi, CBM3, cx_cold, &
            CBMI2, dn_rho, JSRREC, YEXTARR)
    endif

! Conversion to rough mesh keeping the intagrals
    if (dn_rho /= 1)    then
        call smooth_int(NIBM , dn_rho, ROC, NA1)
        call smooth_int(PIBM , dn_rho, ROC, NA1)
        call smooth_int(PEBM , dn_rho, ROC, NA1)
        call smooth_int(PBLON, dn_rho, ROC, NA1)
        call smooth_int(PBPER, dn_rho, ROC, NA1)
        call smooth_int(PBEAM, dn_rho, ROC, NA1)
        call smooth_int(SNNBM, dn_rho, ROC, NA1)
        call smooth_int(SNEBM, dn_rho, ROC, NA1)
        if (ABEAM < 1.5d0) then
            call smooth_int(SNIBM1, dn_rho, ROC, NA1) ! H ion source
        else if (ABEAM > 2.5d0) then
            call smooth_int(SNIBM3, dn_rho, ROC, NA1) ! T ion source
        else
            call smooth_int(SNIBM2, dn_rho, ROC, NA1) ! D ion source
        endif
        call smooth_int(SCUBM  , dn_rho, ROC, NA1)
        call smooth_int(stnbdp , dn_rho, ROC, NA1)
        call smooth_int(sdnbdp1, dn_rho, ROC, NA1)
        call smooth_int(sdnbdp2, dn_rho, ROC, NA1)
        call smooth_int(sdnbtp , dn_rho, ROC, NA1)
        call smooth_int(CUFI, dn_rho, ROC, NA1, surf=.true.)
        call smooth_int(CUBM, dn_rho, ROC, NA1, surf=.true.)
    endif

    do j_rho=1, NA1
        PEBM(j_rho) = PEBM(j_rho) - 2.08E-5*SNEBM(j_rho)
        if (cx_flag > 0)  then
            PIBM(j_rho) = PIBM(j_rho) - 0.0024*SNNBM(j_rho)*TI(j_rho)
        endif
        NNBM2(j_rho)  = SNIBM2(j_rho)
        NNBM3(j_rho)  = stnbdp(j_rho)
        SNIBM1(j_rho) = sdnbdp1(j_rho)
        SNIBM2(j_rho) = sdnbdp2(j_rho)
        SNIBM3(j_rho) = sdnbtp(j_rho)
    enddo

    end subroutine NBINJ

!---------------------------------------------------------------------
    subroutine nbspec(AIM1, AIM2, AIM3, AMJ, ZMJ, NA1, &
        yEXTARR, NE, NI, NHYDR, NDEUT, NTRIT, NHE3, NALF, NIZ1, NIZ2, NIZ3,  &
        ZIM1, ZIM2, ZIM3, m_mp, Zspec, n_spec, spc_lbl, IFLAG)
!---------------------------------------------------------------------
! Input (arrays only)
!   NHYDR,  NDEUT,  NTRIT,  NHE3,  NALF,  NIZ1,  NIZ2,  NIZ3,  NI
!   EXTARR(, )
! Output:
!   n_spec    = 1 + total number of ion species
! Arrays: 1(e),  p,  d,  t,  He3,  He4,  Imp1,  Imp2,  Imp3
!   m_mp(9)     m/m_p
!   Zspec(9)    Charge number
!   spc_lbl(9)  number of ion species in the EXTARR
!   EXTARR ! One of arrays for (p, d, t, He3) is spoiled !
!---------------------------------------------------------------------

    use nbstatus, only: n_rho

    double precision, intent(in) :: AMJ, ZMJ, AIM1, AIM2, AIM3
    double precision, intent(in), dimension(*) :: NE, NHYDR, NDEUT, NTRIT, NHE3, NALF, NI, &
        ZIM1, ZIM2, ZIM3, NIZ1, NIZ2, NIZ3
    double precision, intent(out), dimension(*) :: Zspec, m_mp
    double precision, intent(out) :: yEXTARR(n_rho, 9)

    integer :: n_spec, spc_lbl(*), JIHYDR, JIDEUT, JITRIT, JIHE3, JIALF, JN, &
        JIZ1, JIZ2, JIZ3, NA1, JENE, JINI, IFLAG

! Identification of plasma species for NBI
! Polevoy A.R. 18.03.91
    JIHYDR = 0
    JIDEUT = 0
    JITRIT = 0
    JIHE3  = 0
    JIALF  = 0
    JIZ1   = 0
    JIZ2   = 0
    JIZ3   = 0
    JENE   = 0
    JINI   = 0
    n_spec = 0
    IFLAG  = 0

    do JN=1, NA1 ! check for non zero species

        if (NE(JN) > 0.d0) then
            if (JENE < 1) then
                JENE = 1
                n_spec = n_spec + 1
                m_mp(n_spec)  = 1./1836.d0
                Zspec(n_spec) = 1.d0
                spc_lbl(n_spec) = 1
            endif
            yEXTARR(JN, 1) = NE(JN)
        endif

        if (NHYDR(JN) > 0.d0) then
            if (JIHYDR < 1) then
                JIHYDR = 1
                n_spec = n_spec + 1
                m_mp(n_spec)  = 1.d0
                Zspec(n_spec) = 1.d0
                spc_lbl(n_spec) = 2
            endif
            yEXTARR(JN, 2) = NHYDR(JN)
        endif
 
        if (NDEUT(JN) > 0.d0) then
            if (JIDEUT < 1) then
                JIDEUT = 1
                n_spec = n_spec + 1
                m_mp(n_spec)  = 2.d0
                Zspec(n_spec) = 1.d0
                spc_lbl(n_spec) = 3
            endif
            yEXTARR(JN, 3) = NDEUT(JN)
        endif

        if (NTRIT(JN) > 0.d0) then
            if (JITRIT < 1) then
                JITRIT = 1
                n_spec = n_spec + 1
                m_mp(n_spec)  = 3.d0
                Zspec(n_spec) = 1.d0
                spc_lbl(n_spec) = 4
            endif
            yEXTARR(JN, 4) = NTRIT(JN)
        endif

        if (NHE3(JN) > 0.d0) then
            if (JIHE3 < 1) then
                JIHE3 = 1
                n_spec = n_spec + 1
                m_mp(n_spec)  = 3.d0
                Zspec(n_spec) = 2.d0
                spc_lbl(n_spec) = 5
            endif
            yEXTARR(JN, 5) = NHE3(JN)
        endif

        if (NALF(JN) > 0.d0) then
            if (JIALF < 1) then
                JIALF = 1
                n_spec = n_spec + 1
                m_mp(n_spec)  = 4.d0
                Zspec(n_spec) = 2.d0
                spc_lbl(n_spec) = 6
            endif
            yEXTARR(JN, 6)=NALF(JN)
        endif

        if (NIZ1(JN) > 0.d0) then
            if (JIZ1 < 1) then
                JIZ1 = 1
                n_spec = n_spec + 1
                m_mp(n_spec)  = AIM1
                Zspec(n_spec) = ZIM1(1)
                spc_lbl(n_spec) = 7
            endif
            yEXTARR(JN, 7) = NIZ1(JN)
        endif

        if (NIZ2(JN) > 0.d0) then
            if (JIZ2 < 1) then
                JIZ2 = 1
                n_spec = n_spec + 1
                m_mp(n_spec)  = AIM2
                Zspec(n_spec) = ZIM2(1)
                spc_lbl(n_spec) = 8
            endif
            yEXTARR(JN, 8) = NIZ2(JN)
        endif

        if (NIZ3(JN) > 0.d0) then
            if (JIZ3 < 1) then
                JIZ3 = 1
                n_spec = n_spec + 1
                m_mp(n_spec)  = AIM3
                Zspec(n_spec) = ZIM3(1)
                spc_lbl(n_spec)  = 9
            endif
            yEXTARR(JN, 9) = NIZ3(JN)
        endif

        if (NI(JN) > 0.d0) then
            if (JINI < 1) JINI = 1
        endif
    enddo

!  No plasma
    if (spc_lbl(1) < 1) then
        write(*, *) 'quit NBI: no ionised plasma, Ne is not specified'
        IFLAG = 1
        return
    endif

! Main ion specia identification
    if (spc_lbl(2) < 1) then ! no ion species
        n_spec = 2
        spc_lbl(n_spec) = n_spec
        if ((AMJ*ZMJ*JINI) > 0.d0) then
            do jn=1, NA1
                yEXTARR(jn, n_spec) = NI(jn)
                m_mp(n_spec)  = AMJ
                Zspec(n_spec) = ZMJ
            enddo
        else
            do jn=1, NA1
                yEXTARR(jn, n_spec) = NI(jn)
                m_mp(n_spec)  = 1.d0
                Zspec(n_spec) = 1.d0
            enddo
            write(*, *) 'NBI: no ion specie is specified. Continue with hydrogen'
            IFLAG = 0
        endif
    endif

    end subroutine nbspec

!---------------------------------------------------------------------
    subroutine smooth_int(YFO, dn_rho, YROC, JNA1, surf)
!---------------------------------------------------------------------
!  The subroutine transmits a histogram like function YFO(1:N1)
!  from an N1 grid X(1:N1) to a smooth functon YFO(1:JNA1)
!  keeping the same total volume integral.
! Input: YFO(1:N1)
! Output: YFO(1:JNA1)
!---------------------------------------------------------------------

    use numerical_tools, only: smooth
    use standard_functions, only: VINT, IINT
    use nbicom, only: N1, a_norm, a_norm_astra

    integer, intent(in) :: dn_rho, JNA1
    double precision, intent(in) :: YROC
    double precision, intent(inout) :: YFO(*)
    logical, intent(in), optional :: surf

    logical :: surf_int
    integer :: j, JNA, JNAC, JSIGN
    double precision :: Y, YOLD, YNEW, ALFA
    double precision, allocatable, dimension(:) :: dri

    if (present(surf)) then
        surf_int = surf
    else
        surf_int = .false.
    endif

    allocate(dri(n1))

    JSIGN = 0
    do j=1, n1-1
        JNA = 1 + dn_rho*(j - 1)
        JNAC = JNA - 1 + dn_rho
        a_norm(j) = a_norm_astra(JNAC)
        DRI(j) = YFO(JNA)
        if (JSIGN < 1) then
            if (DRI(1)*DRI(J) < 0.d0) JSIGN = 1
        endif
    enddo

    DRI(N1)    = 0.d0
    a_norm(N1) = 1.d0
    a_norm_astra(JNA1) = 1.d0

! Total power normalization

    if (surf_int) then
        YOLD = IINT(YFO, YROC)
    else
        YOLD = VINT(YFO, YROC)
    endif

    ALFA = 0.001d0
    call SMOOTH(ALFA, a_norm, dri, n1, a_norm_astra, yfo, jna1)

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

    deallocate(dri)

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
