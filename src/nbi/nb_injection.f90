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
        TI, RHO, VR, ZIM1, ZIM2, ZIM3
    use nbibce, only: nbionr
    use nbicom, only: nbspec, nbsrsr, nbion0, stnbdp, sdnbtp, sdnbdp1, sdnbdp2, RMB, ZB, YASBA

    integer, intent(in) :: NA1, NB1
    double precision, intent(in) :: BTOR, RTOR, ABC, AB, ROC, SHIFT, UPDWN, HRO, TAU, &
        AIM1, AIM2, AIM3, AMJ, ZMJ, NNCL, NNWM, CBM2, CBM3, CBM4, CBMI1, CBMI2
    character(len=*), intent(in) :: file_nbi
    double precision, intent(out) :: QNBI, CBM1, CBMI3, CBMI4

    logical :: file_exists
    integer :: N, JSRNUM, J, J2FRST, JN, JN1, JNAX, ERCODE, JWARN, INBMS, &
        JSRREC, jABEAM, jZBEAM, IFLAG
    double precision :: CBMH1, CBMH2, CBMR1, CBMR2, CBMS1, CBMS2, CBMS3, CBMS4, &
        EBEAM, DBM1, DBM2, DBM3, ABEAM, HBEAM, CONTR, RBMAX, RBMIN, QBEAM
    double precision :: YQBEAM, ZBEAM, YCBMI3, YCBMI4, YABEAM, YEBEAM, YUD, YHM
    double precision :: ARRAY(n_fields), YCOS(n_theta), yEXTARR(n_rho, 9)
    character(len=25) :: STRI
    character(len=40) :: YVARNAME

!---- Ripple normalized radius of the ripple boundary
!---- banana with RTCRIT > YRIPLR is lost
    double precision :: Y, Y1, Y2, YMAX, YQBM
    real*8, allocatable :: Y4TORIC(:, :, :), YQBMJE(:, :), yetmp(:), ypwtmp(:), ypartmp(:)
    integer jElevAll, je, jiounit, JT

    save J2FRST, INBMS
    data J2FRST, INBMS, YCBMI3, YCBMI4 /0, 0, 1.d0, 1.d0/
    data ABEAM/1.d0/ DBM1/1.d0/ DBM2/0.d0/ DBM3/0.d0/
    data RBMIN/1.d0/ RBMAX/.5d0/

    jElevAll = 0 ! for write to SSFPQL
    ZBEAM = 1.d0 ! Hydrogen isotopes only

    inquire(file=TRIM(file_nbi), exist=file_exists)

    if (.not. file_exists) write(*, *) 'File "', TRIM(file_nbi), '" not found'
    J = abs(int(CBM1 + 0.5d0))
    if (CBM1 < .0) J = J + 1
    if (J == 0) then
        write(*, *) '>>> NBI >>> Zero number of sources'
        write(*, *)"            Don't know what to do"
        stop
    endif

    if (INBMS == 0) INBMS = J ! 1st call

    if (INBMS /= J) then  ! No. of sources changed
        INBMS = J
        call NQUERY(2, TRIM(file_nbi), INBMS, ERCODE)
        if (ERCODE /= 0) then
            write(*, *) "NQUERY Errcode =", ERCODE
            write(*, *) '>>> NBI >>> Error in file "', TRIM(file_nbi), '": unrecognized variable name'
            stop
        endif
    else if (.not. file_exists .or. CBM1 <= 0.d0 ) then
! Interactive viewer/editor of the beam source parameter list
        call NQUERY(2, TRIM(file_nbi), INBMS, ERCODE)
        if (ERCODE /= 0) then
            write(*, *) "NQUERY Errcode =", ERCODE
            write(*, *) '>>> NBI >>> Error in file "', TRIM(file_nbi), '": unrecognized variable name'
            stop
        endif
    endif

    if (CBMI1 == 0.d0) then ! print for SSFPQL
        if (.not. allocated(YQBMJE)) allocate (YQBMJE(INBMS, 3))
        do jn=1, INBMS
            do je=1, 3
                YQBMJE(jn, je) = 0.d0
            enddo
        enddo
    endif ! print for SSFPQL

! CBMI3 to make NBI internal mesh interval > max Larmor raduis
    Y = 1.d-3
    EBEAM = 1.d-3

    open(2, file=TRIM(file_nbi), status='OLD')
    do j=1, INBMS
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
        if (CBMI1 > 1.d0 .and. J2FRST == 0) then
            J2FRST = 1
            YCBMI3 = CBMI3
            YCBMI4 = CBMI4
        endif
        if (CBMI1 > 1.d0 .and. J2FRST /= 0) then
            if (CBMI3 /= YCBMI3 .or. CBMI4 /= YCBMI4) write(*, *) 'CBMI3, CBMI4 can`t be changed after CNB4=2'
            CBMI3 = YCBMI3
            CBMI4 = YCBMI4
        endif
    else  ! w/o FP solver FI source on the transport grid
        CBMI3 = 1.d0 !!!!
    endif

    QNBI = 0.d0
    CBM1 = INBMS

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

!---- Ripple
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

    if (INBMS == 0) then
        write(*, *) ' -1 < CBM1 < 1 = no NBI sources '
        return
    endif

    JWARN = 0
    jElevAll = 0

    YEBEAM = EBEAM
    YABEAM = ABEAM
    YQBM = 0.d0

    open(2, file=TRIM(file_nbi), status='OLD')

    do JN=1, INBMS
        JSRNUM = JN
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
                    YQBMJE(jn, 1) = QBEAM*DBM1
                    jElevAll = jElevAll + 1
                endif
                if (DBM2 > 0.d0) then
                    YQBMJE(jn, 2) = QBEAM*DBM2
                    jElevAll = jElevAll + 1
                endif
                if (DBM3 > 0.d0) then
                    YQBMJE(jn, 3) = QBEAM*DBM3
                    jElevAll = jElevAll + 1
                endif
            endif
        endif

! Beam footprint drawing:
        write(STRI, '(a, i3, a)')" NBINJ: Beam #", JN, "    "//char(0)

        if (YABEAM /= ABEAM) JWARN = 1
        QNBI = QNBI + QBEAM
        YQBEAM = QBEAM

! Switch off the ion source
        YQBM = YQBM + YQBEAM
        YHM = (HBEAM - YUD)*100.
        if (CONTR /= 0.d0 .and. CONTR /= 1.d0) then
! balanced injection
! coinj. part
            QBEAM = YQBEAM*(1. - CONTR)
            call NBSRSR(JSRNUM, -1.d0, NA1, RTOR, SHIFT, AB, BTOR, HRO, YHM, &
                CBMH1, CBMH2, CBMS1, CBMS2, CBMS3, CBMS4, CBMR1, CBMR2, CBMI3, CBMI1, &
                EBEAM, DBM1, DBM2, DBM3, ABEAM, QBEAM, RBMAX, RBMIN, JSRREC, &
                YEXTARR, CBMI4)

            if (CBMI1 == 1.d0) call NBION0(NA1, ABEAM, EBEAM, RTOR, CBMI3, yEXTARR)

!: counter injection part...
            QBEAM = YQBEAM*CONTR
            call NBSRSR(JSRNUM, 1.d0, NA1, RTOR, SHIFT, AB, BTOR, HRO, YHM, &
                CBMH1, CBMH2, CBMS1, CBMS2, CBMS3, CBMS4, CBMR1, CBMR2, CBMI3, CBMI1, &
                EBEAM, DBM1, DBM2, DBM3, ABEAM, QBEAM, RBMAX, RBMIN, JSRREC, &
                YEXTARR, CBMI4)
            if (CBMI1 == 1.d0) call NBION0(NA1, ABEAM, EBEAM, RTOR, CBMI3, yEXTARR)
            QBEAM = YQBEAM
        else
            if (CONTR == 1.d0) call NBSRSR(JSRNUM, 1.d0, NA1, RTOR, SHIFT, AB, BTOR, &
                HRO, YHM, CBMH1, CBMH2, CBMS1, CBMS2, CBMS3, CBMS4, &
                CBMR1, CBMR2, CBMI3, CBMI1, EBEAM, DBM1, DBM2, DBM3, ABEAM, &
                QBEAM, RBMAX, RBMIN, JSRREC, YEXTARR, CBMI4)
            if (CONTR == 0.d0) call NBSRSR(JSRNUM, -1.d0, NA1, RTOR, SHIFT, AB, BTOR, &
                HRO, YHM, CBMH1, CBMH2, CBMS1, CBMS2, CBMS3, CBMS4, &
                CBMR1, CBMR2, CBMI3, CBMI1, EBEAM, DBM1, DBM2, DBM3, ABEAM, &
                QBEAM, RBMAX, RBMIN, JSRREC, YEXTARR, CBMI4)
            if (CBMI1 == 1.d0) call NBION0(NA1, ABEAM, EBEAM, RTOR, CBMI3, yEXTARR)
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

        call NBIONR(EBEAM, ABEAM, RTOR, NA1, TAU, NNCL, NNWM, &
            CBM1, CBM3, CBM4, CBMI2, CBMI3, JSRREC, YEXTARR)
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
        if (ABEAM < 1.5d0) call smooth_int(SNIBM1, CBMI3, ROC, NA1) !H ion source
        if (ABEAM < 2.5d0 .and. ABEAM > 1.5d0) call smooth_int(SNIBM2, CBMI3, ROC, NA1) !D ion source
        if (ABEAM > 2.5d0) call smooth_int(SNIBM3, CBMI3, ROC, NA1) !T ion source
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

    if (CBMI1 == 0.d0) then ! write to SSFPQL
        do j=1, n_theta
            YCOS(j) = -1.d0 + 2.d0*((J - 1.0)/(N_THETA - 1.0))
        enddo
        if (jElevAll == 0) then !no power in NBI
            jElevAll = 1
            if (.not. allocated(Y4TORIC)) allocate(Y4TORIC(jElevAll, n_theta, NA1))
            if (.not. allocated(yetmp)) allocate(yetmp(jElevAll))
            if (.not. allocated(ypwtmp)) allocate(ypwtmp(jElevAll))
            if (.not. allocated(ypartmp)) allocate(ypartmp(jElevAll))
            do je=1, jElevAll
                yetmp(je)   = EBEAM
                ypwtmp(jE)  = 0.d0
                ypartmp(jE) = 0.d0
                do j=1, n_theta
                    do jn=1, NA1-1
                        Y4TORIC(jE, j, jn) = 0.d0
                    enddo
                enddo
            enddo
        else ! power in NBI
            if (.not. allocated(Y4TORIC)) allocate(Y4TORIC(jElevAll, n_theta, NA1))
            if (.not. allocated(yetmp)) allocate(yetmp(jElevAll))
            if (.not. allocated(ypwtmp)) allocate(ypwtmp(jElevAll))
            if (.not. allocated(ypartmp)) allocate(ypartmp(jElevAll))
            do je=1, jElevAll
                yetmp(je)   = 0.d0
                ypwtmp(je)  = 0.d0
                ypartmp(je) = 0.d0
            enddo

            open(35, FILE='dat/srsfi.dat', FORM='UNFORMATTED', STATUS='UNKNOWN', ACCESS='DIRECT', RECL=JSRREC)
            J = 0  ! for JELEVEL
            do JSRNUM=1, INBMS ! cycle for NBI sources
                read(35, REC=JSRNUM, ERR=990) EBEAM, &
                    (((YASBA(JE, JN, JT), JE=1, 3), JN=1, NA1-1), JT=1, n_theta)
! JE=3 for full energy EBEAM
                if (YQBMJE(JSRNUM, 1) > 0.d0) then
                    j = j + 1
                    yetmp(j) = EBEAM
                    do jt=1, n_theta
                        do jn=1, NA1
                            if (jn < na1) then
                                Y = YASBA(3, JN, JT)
                            else
                                Y = 0.d0
                            endif
                            Y4TORIC(j, jt, jn) = Y*1.d19
                            ypwtmp(j) = ypwtmp(j) + Y*1.6d-3*yetmp(j)
                            ypartmp(j) = ypartmp(j) + Y4TORIC(j, jt, jn)
                        enddo
                    enddo
                endif
                if (YQBMJE(JSRNUM, 2) > 0.d0) then
                    j = j + 1
                    yetmp(j) = EBEAM/2.d0
                    do jt=1, n_theta
                        do jn=1, NA1
                            if (jn < na1) then
                                Y = YASBA(3, JN, JT)
                            else
                                Y = 0.d0
                            endif

                            Y4TORIC(j, jt, jn) = Y*1.d19
                            ypwtmp(j) = ypwtmp(j) + Y*1.6d-3*yetmp(j)
                            ypartmp(j) = ypartmp(j) + Y4TORIC(j, jt, jn)
                        enddo
                    enddo
                endif
                if (YQBMJE(JSRNUM, 3) > 0.d0) then
                    j = j + 1
                    yetmp(j) = EBEAM/3.d0
                    do jt=1, n_theta
                        do jn=1, NA1
                            if (jn < na1) then
                                Y = YASBA(3, JN, JT)
                            else
                                Y = 0.d0
                            endif

                            Y4TORIC(j, jt, jn) = Y*1.d19
                            ypwtmp(j) = ypwtmp(j) + Y*1.6d-3*yetmp(j)
                            ypartmp(j) = ypartmp(j) + Y4TORIC(j, jt, jn)
                        enddo
                    enddo
                endif
            enddo
            close(35)
        endif ! end power in NBI

        jsrnum = 1 ! number of species with differet mass
        jiounit = 36
        jABEAM = ABEAM
        jZBEAM = ZBEAM

        open(jiounit, file='dat/toric.nbi')
        YVARNAME = 'ASTRA'
        write(jiounit, '(A40)') YVARNAME
        write(jiounit, '(3I5)') jsrnum, NA1, n_theta !Energy, Space, Pitch angle
        write(jiounit, '(A)')   'Radial mesh step'
        write(jiounit, '(E17.9)') HRO/RHO(NA1)
        write(jiounit, '(A)') 'Radial mesh SQRT(NormTorFlux)'
        write(jiounit, '(6E17.9)') (RHO(j)/RHO(NA1), j=1, NA1)
        write(jiounit, '(A)')   'Specific volumes (m^3)'
        write(jiounit, '(6E17.9)') (VR(j)*HRO, j=1, NA1)
        write(jiounit, '(A)')   'Cos(PitchAngle) Mesh'
        write(jiounit, '(6E17.9)') (YCOS(j), j=1, n_theta)

        do j=1, jsrnum  ! count of NBI species (=1)
            write(jiounit, '(A11, I3)') 'NBI Species', j
            write(jiounit, '(3I5)') jABEAM, jZBEAM, jElevAll
            write(jiounit, '(A)')   'Energy levels, keV'
            write(jiounit, '(6E17.9)') (yetmp(je), je=1, jElevAll)
            write(jiounit, '(A)')  'Powers, MW'
            write(jiounit, '(6E17.9)') (ypwtmp(je), je=1, jElevAll)
            write(jiounit, '(A)')  'Ionization rates, prtcl/sec'
            write(jiounit, '(6E17.9)') (ypartmp(je), je=1, jElevAll)
            write(jiounit, *) ' Particle source(jE, jCos, jRad), 1/sec'
            do je=1, jElevAll
                do jn=1, NA1
                    write(jiounit, '(6E17.9)') (Y4TORIC(je, jt, jn), jt=1, n_theta)
                enddo
            enddo
        enddo

        close(jiounit)
        if (allocated(Y4TORIC)) deallocate(Y4TORIC)
        if (allocated(YQBMJE)) deallocate (YQBMJE)
        if (allocated(yetmp)) deallocate (yetmp)
        if (allocated(ypwtmp)) deallocate (ypwtmp)
        if (allocated(ypartmp)) deallocate (ypartmp)

    endif !write to SSFPQL

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
    subroutine NQUERY(NCH, file_name, NBS, ERCODE)
!---------------------------------------------------------------------
! Input:
!   NBS    - Requested No. of NBI sources
!   NCH    - Logical unit (must be connected to the file NBINP)
! Output:   Data are written into file NBINP
!   ERCODE - Error code (0 for normal exit)
!---------------------------------------------------------------------

    use char_manip, only: to_upper, first_non_blank
    use dbl2char, only: isnum, fmt_smart
    use nbstatus, only: n_fields, n_nbi_max

    integer, intent(in) :: NCH, NBS
    character(len=*), intent(in) :: file_name
    integer, intent(out) :: ERCODE

    integer :: j, j1, jj, i, ios
    double precision :: YCB
    character(len=6) :: val_str
    character(len=12), dimension(n_fields) :: SFIELD
    character(len=80) :: STR, STRI, ADATA(n_nbi_max), CDATA(n_nbi_max)

    if (NBS > n_nbi_max) then
        write(*, '(2(A, I3, A))') &
           '>>> NQUERY: Too many NB sources NBS =', NBS, ' requested', &
           '            Call ignored,  NBSmax =', n_nbi_max, ' is allowed'
        return
    endif

    if (NBS <= 0) then
        write(*, '(2(A, I3))') '>>> NQUERY: Number of sources must be positive, NBS =', NBS
        return
    endif

! Default definition:
    do jj=1, NBS
        write(ADATA(jj)(1:2), '(1I2)') jj
        write(CDATA(jj)(1:2), '(1I2)') jj
        ADATA(jj)(3:) = "   zrd1     0.     2.     1.   100.     .7     .2     .1     1.    10."
        CDATA(jj)(3:) = "     .3    1.6     1.     0.     1.     9.     2.     9.     2.     0."
    enddo

    open(NCH, file=TRIM(file_name), status='OLD', iostat=ios)
    if (ios /= 0) then
        open(NCH, file=TRIM(file_name), status='NEW')
    endif
    do jj=1, NBS
        read(NCH, '(A)', end=10) STR
        j = index(STR, '!')
        if (j == 0) then ! Don't skip if not starting with "!"
            read(NCH, '(5A)', err=99) SFIELD
            i = -3
            do j=1, n_fields
                i = i + 7
                SFIELD(j)(1:12) = to_upper(SFIELD(j)(1:12))
                write(*, *) 'nbi.sfield', SFIELD(j)(1:12)
                if (isnum(SFIELD(j), 12)) then
                    read(SFIELD(j), *) YCB
                    val_str = fmt_smart(YCB, 6)
                    write(*, *) '"nbi.value', val_str, '"', YCB
                else
                    j1 = first_non_blank(SFIELD(j))
                    val_str = SFIELD(j)(j1:)
                endif
                if (j <= 10) then
                    write(ADATA(jj)(i:i+5), '(1A6)') val_str
                else
                    write(CDATA(jj)(i:i+5), '(1A6)') val_str
                endif
                if (j == 10) i = -3
            enddo
        endif
    enddo

 10 continue

! Template for the 1st string of a table
    STRI = "# |QBeam |Contr |ABeam |ZBeam |EBeam |DBeam1|DBeam2|DBeam3|Orb_av|Penc.#" // char(0)
    STR = "NBI configuration file: " // TRIM(file_name) // char(0)
    j = len(ADATA(1))
    call NBIBOX(STR, STRI, ADATA, j, NBS, 0, 0)

    STRI = "# |HBeam |RBmax |RBmin |tg(A) |Aspect|Cver1 |Cver2 |Chor1 |Chor2 |Unused" // char(0)
    STR = "NBI configuration file (cnt.): " // TRIM(file_name) // char(0)
    j = len(CDATA(1))
    call NBIBOX(STR, STRI, CDATA, j, NBS, 0, 0)

    rewind(NCH)  ! Modify input file

    if (ios /= 0) then
        write(NCH, '(1I3)') 1
    else
        j = 1
        do while(j /= 0) ! Skip lines starting with "!"
            read(NCH, '(A)', end=99) STR
            j = index(STR, '!')
        enddo
    endif

    do jj=1, NBS
        i = -3
        do j=1, n_fields
            i = i + 7
            if (j <= 10) then
                 write(val_str, '(1A6)') ADATA(jj)(i:i+5)
            else
                 write(val_str, '(1A6)') CDATA(jj)(i:i+5)
            endif
            val_str(1:6) = to_upper(val_str(1:6))
            if (ISNUM(val_str, 6)) then
                read(val_str, *)YCB
                write(SFIELD(j), '(1P, E12.4)')YCB
            else
                j1 = first_non_blank(val_str)
                val_str = val_str(j1:) // '     '
                SFIELD(j)(1:) = val_str // '      '
            endif
            if (j == 10) i = -3
        enddo
        if (jj /= 1) write(NCH, '(1I3)') jj
        write(NCH, '(3(5A12/), 5A12)')SFIELD
    enddo

    close(NCH)
    ERCODE = 0
    return

99  continue
    close(NCH)
    ERCODE = 1

    end subroutine NQUERY

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

    use scalars , only: constValues, varxValues
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
