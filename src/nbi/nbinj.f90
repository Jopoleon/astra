! Neutral Beam driven current, momentum and power deposition
!    including Ripple Losses
!-----------------------------------------------------------------------
!	entry:	AMETR, SHIF, NA1, RTOR, AB, BTOR, NI, HBEAM,
!     		RBMIN, RBMAX, CBMS2, CBMI3, EBEAM, QBEAM, ABEAM, DBM1, DBM2, DBM3,
!	      MU, TE, TI, NE, NN, NNB, ZEF, AMAIN, PBEAM, SCUBM, CONTR, ELON, TRIAN
! NBI interactive control:
! <VARIABLES>:
! ABEAM=1, 2, 3 Hydrogen, Deuterium, Tritium NBI only  (Def=1)
! CONTR=0-1 fraction of contr-inj. power   (Def=1)
! DBM1=0-1 ~ full energy component EBEAM power fraction  (Def=1)
! DBM2=0-1 ~ half energy component EBEAM/2 power fraction (Def=0)
! DBM3=0-1 ~third energy component EBEAM/2 power fraction (Def=0)
! WARNING: DBM1 > 0 must be used; SUM(DBMi) /= 1. could be used
! Power is renormalized withinNBI: Power(i)~DBMi/SUM(DBMi)
! EBEAM [keV] the main component energy   (Def=1)
! HBEAM [m] NBI footprint center height in respect to Z=0 (Def=0)
! RBMAX [m] max. major radius of the NBI footprint (Def=0)
! RBMIN [m] min. major radius of the NBI footprint (Def=0)
! QBEAM [MW] NBI power      (Def=0)
!
! Comments: the footprint position is interpreted as:
!
! for tangential NBI: (RBMAX+RBMIN)/2 > RTOR - AB
! as NBI cross-section by the meridianal plane perpendicular
! to the vertical plane of the central NBI pencil beam
! with tangential radius of (RBMAX+RBMIN)/2
!
! for perpendicular NBI: (RBMAX+RBMIN)/2 < RTOR - AB
! as NBI cross-section by the verticall plane perpendicular to
! the central NBI pencil beam (RBMAX+RBMIN)/2 vertical plane,
! placed at R = RTOR (the corrections, connected with
! Shafranov's shift are calculated into the NBI block
!
! <CONSTANTS>
! CBM1=1, 2, .. number of sorces with different geometry (Def=1)
!  if CBM1 > 1 source parameters are controlled by
!  special NBMENU rather than VARIABLE, CONSTANTS
!  menues for each source separately. NBMENU is called
!  automatically when you increase CBM1. User can
!  call NBMENU anytime by changing a sign of CBM1
!  to negative (CBM1 = 4 to CBM1 = -4, etc.)
!
! CX heat losses due to NBI control
!
! CBM2=1  heat losses (=0 No losses) from bulk ions due
!  to NBI charge exchange    (Def=1)
!   CBM2=1 corresponds to the explicit calculations
!  of the heat sink from the bulk ion component at
!  the moment when NBI is called.
!  Implicit form is recommended:
!  CBM2=0; PIT=...-.0024*SNNBM
!
! CBM3, 4 are active for CBMI1=2 only
!
! CBM3=1  losses (=0 No losses) from fast ions due to CX
!  with the NBI neutrals itself   (Def=1)
!  The same CX cross-section as for cold neutrals
!  is considered.
! CBM4=1  losses (=0 No losses) from fast ions due to CX
!  with the cold neutrals    (Def=1)
!
! NBI Source control:
!
! CBMH1, 2 reserved for the footprint power profile control(Def=0)
! CBMR1-4 reserved for the footprint power profile control(Def=0)
!
! CBMS1=0 No averaging over drift orbit
! CBMS1=1 Averaging over drift oribt
! CBMS1>1 Zero-width drift orbit
! CBMS2=1, ..., 82 number of 'pencils'   (Def=1)
! CBMS3=tg(angle between central pencil and midplane) (Def=1)
! Comments: CBMS3 is used for perpendicular NBI only
! CBMS4=(HBmax-HBmin)/(RBMAX-RBMIN): HBEAM=(HBmax+HBmin)/2(Def=1)
!
! Fast ion solver control
!
! CBMI1=1(2) steady (time dependent) Fokker-Planck Solver (Def=1)
!  number of angle mesh points NTET1=50/CBMI1+1
! CBMI1=0 no FP solver, only ionisation source Y4TORIC is calculated
! CBMI2-4 are active for CBMI1=2 only
!
! CBMI2=TAUNBI/TAU ratio of the NBI time-step to TAU (Def=1)
! WARNING: TAUNBI does not coinside automatically with the
!  interval between sequental calls of NBI (that
!  enables the consumption of calculations at the
!  NBI-steady state phase)
! CBMI3=1, 2... NBI X-mesh points number N1=(NA1-1)/CBMI3+1(Def=1)
!       3(EB, EB/2, EB/3), 2(EB, EB/2), 1(EB)
! CBMI4=0, 1, 2...control of beam-plasma fusion              (Def=1)
!  0 (no fusion), 1 (finite Te, Ti), 2 (finite Te, Ti=0)
!       NBI V-mesh points number IV1=161
!
! USER's Functions (NBUSER.f) shold be treated from user's sbr/
!
! NBFHZ(NBI_No, Z, CBMH1, CBMH2), NBFRY(NBI_No, Y, CBMR1, CBMR2)
!  power distribution in the footprint profile
!  (see comments in the subroutines)
!
! RIPRAD(ZUPDWN, J) [m] - Ripple loss cone boundary (major
!  radius) for each magnetic surface: ZUPDWN [m]-
!  shift in respect to the midplane of ripple
!  simmetry, J - surface index
! RIPRAD -depends on tokamak mag. field coils' and
!  plasma configurations  (Def: No Ripple losses) (Def=999)
! exit: PBEAM, PEBM, PIBM, NIBM, CUFI, CUBM, PBLON, PBPER,
!  SCUBM, SNEBM, SNNBM, NNBM1, 2, 3 for MAIN
!=============================================================
!  stnbdp [1/s]- intensity of burn out of T from T NBI
!   on D bulk component
!  Tritium sink/n14.1MeV+3.52MeVHe4 source: S= NDEUT(j)*stnbdp(j)
!
!  sdnbtp [1/s]- intensity of burn out of D from D NBI
!   on T bulk component
! Deuterium sink/n14.1MeV+3.52MeVHe4 source: S= NTRIT(j)*sdnbtp(j)
!
! sdnbdp1 [1/s]- intensity of burn out of D from D NBI
!   on D bulk component
! Deuterium sink/t1.008MeV+3.025MeVp source: S= NDEUT(j)*sdnbdp1(j)
!
! sdnbdp2 [1/s]- intensity of burn out of D from D NBI
!   on D bulk component
! Deuterium sink/n2.45MeV+0.87MeVHe3 source: S= NDEUT(j)*sdnbdp2(j)
!======================================================= 01-AUG-2012 Polevoi

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

implicit none

integer, intent(in) :: NA1, NB1
double precision, intent(in) :: BTOR, RTOR, ABC, AB, ROC, SHIFT, UPDWN, HRO, TAU, &
    AIM1, AIM2, AIM3, AMJ, ZMJ, NNCL, NNWM, CBM2, CBM3, CBM4, CBMI1, CBMI2
character(len=*), intent(in) :: file_nbi
double precision, intent(out) :: QNBI, CBM1, CBMI3, CBMI4

logical :: EXI
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
double precision :: RIPRAD, Y, Y1, Y2, YMAX, YQBM
real*8, allocatable :: Y4TORIC(:, :, :), YQBMJE(:, :), yetmp(:), ypwtmp(:), ypartmp(:)
integer jElevAll, je, jiounit, JT

save J2FRST, INBMS
data J2FRST, INBMS, YCBMI3, YCBMI4 /0, 0, 1.d0, 1.d0/
data ABEAM/1.d0/ DBM1/1.d0/ DBM2/0.d0/ DBM3/0.d0/
data RBMIN/1.d0/ RBMAX/.5d0/

jElevAll = 0 ! for write to SSFPQL
ZBEAM = 1.d0 ! Hydrogen isotopes only

inquire(file=TRIM(file_nbi), exist=EXI)

if (.not. EXI) write(*, *) 'File "', TRIM(file_nbi), '" not found'
J = abs(int(CBM1+0.5d0))
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
else if (.not. EXI .or. CBM1 <= 0.d0 ) then
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
Y1 =0.
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
do JN1=1, 9
    ISPE(JN1) = 0
    do JN=1, NB1
        yEXTARR(JN, JN1) = 0.d0
    enddo
enddo

! Set plasma composition
call NBSPEC(AIM1, AIM2, AIM3, AMJ, ZMJ, NA1, NB1, &
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
        call NBSRSR(JSRNUM, -1.d0, NA1, RTOR, SHIFT, AB, BTOR, NNCL, NNWM, HRO, YHM, &
            CBMH1, CBMH2, CBMS1, CBMS2, CBMS3, CBMS4, CBMR1, CBMR2, CBMI3, CBMI1, &
            EBEAM, DBM1, DBM2, DBM3, ABEAM, CONTR, QBEAM, RBMAX, RBMIN, JSRREC, &
            YEXTARR, CBMI4)

        if (CBMI1 == 1.d0) call NBION0(NA1, NNCL, NNWM, ABEAM, EBEAM, RTOR, CBMI3, yEXTARR)

!: counter injection part...
        QBEAM = YQBEAM*CONTR
        call NBSRSR(JSRNUM, 1.d0, NA1, RTOR, SHIFT, AB, BTOR, NNCL, NNWM, HRO, YHM, &
            CBMH1, CBMH2, CBMS1, CBMS2, CBMS3, CBMS4, CBMR1, CBMR2, CBMI3, CBMI1, &
            EBEAM, DBM1, DBM2, DBM3, ABEAM, CONTR, QBEAM, RBMAX, RBMIN, JSRREC, &
            YEXTARR, CBMI4)
        if (CBMI1 == 1.d0) call NBION0(NA1, NNCL, NNWM, ABEAM, EBEAM, RTOR, CBMI3, yEXTARR)
        QBEAM =YQBEAM
    else
        if (CONTR == 1.d0) call NBSRSR(JSRNUM, 1.d0, NA1, RTOR, SHIFT, AB, BTOR, &
            NNCL, NNWM, HRO, YHM, CBMH1, CBMH2, CBMS1, CBMS2, CBMS3, CBMS4, &
            CBMR1, CBMR2, CBMI3, CBMI1, EBEAM, DBM1, DBM2, DBM3, ABEAM, CONTR, &
            QBEAM, RBMAX, RBMIN, JSRREC, YEXTARR, CBMI4)
        if (CONTR == 0.d0) call NBSRSR(JSRNUM, -1.d0, NA1, RTOR, SHIFT, AB, BTOR, &
            NNCL, NNWM, HRO, YHM, CBMH1, CBMH2, CBMS1, CBMS2, CBMS3, CBMS4, &
            CBMR1, CBMR2, CBMI3, CBMI1, EBEAM, DBM1, DBM2, DBM3, ABEAM, CONTR, &
            QBEAM, RBMAX, RBMIN, JSRREC, YEXTARR, CBMI4)
        if (CBMI1 == 1.d0) call NBION0(NA1, NNCL, NNWM, ABEAM, EBEAM, RTOR, CBMI3, yEXTARR)
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

    call NBIONR(EBEAM, ABEAM, AMJ, RTOR, NA1, TAU, NNCL, NNWM, CBM1, CBM2, &
        CBM3, CBM4, CBMI1, CBMI2, CBMI3, CBMI4, JSRREC, YEXTARR)
endif

! Conversion to rough mesh keeping the intagrals
if (CBMI3 /= 1.d0)    then
    call NBSMTV(NIBM, CBMI3, ROC, NA1)
    call NBSMTV(PIBM, CBMI3, ROC, NA1)
    call NBSMTV(PEBM, CBMI3, ROC, NA1)
    call NBSMTV(PBLON, CBMI3, ROC, NA1)
    call NBSMTV(PBPER, CBMI3, ROC, NA1)
    call NBSMTV(PBEAM, CBMI3, ROC, NA1)
    call NBSMTV(SNNBM, CBMI3, ROC, NA1)
    call NBSMTV(SNEBM, CBMI3, ROC, NA1)
    if (ABEAM < 1.5d0) call NBSMTV(SNIBM1, CBMI3, ROC, NA1) !H ion source
    if (ABEAM < 2.5d0 .and. ABEAM > 1.5d0) call NBSMTV(SNIBM2, CBMI3, ROC, NA1) !D ion source
    if (ABEAM > 2.5d0) call NBSMTV(SNIBM3, CBMI3, ROC, NA1) !T ion source
    call NBSMTV(SCUBM, CBMI3, ROC, NA1)
    call NBSMTS(CUFI, CBMI3, ROC, NA1)
    call NBSMTS(CUBM, CBMI3, ROC, NA1)
    call NBSMTV(stnbdp, CBMI3, ROC, NA1)
    call NBSMTV(sdnbdp1, CBMI3, ROC, NA1)
    call NBSMTV(sdnbdp2, CBMI3, ROC, NA1)
    call NBSMTV(sdnbtp, CBMI3, ROC, NA1)
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
            yetmp(je) = EBEAM
            ypwtmp(jE) = 0.d0
            ypartmp(jE) = 0.d0
            do j=1, n_theta
                do jn=1, NA1-1
                    Y4TORIC(jE, j, jn)=0.d0
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

    do j = 1, jsrnum  ! count of NBI species (=1)
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

93 continue
write(*, *)'>>> NBI >>> Wrong NBI configuration file format: '
write(*, '(6X, I2, A, I2, A)') JN-1, ' beam records available, ', &
    int(CBM1), ' records required.'
stop

95 continue
write(*, *) '>>> NBI calling STREAD: array out of limits'
stop

96 continue
write(*, *) '>>> NBI >>> Error in file "', TRIM(file_nbi), '": unrecognized variable name'
stop

97 continue
write(*, *) '>>> NBI >>> File "', TRIM(file_nbi), '" read error'
stop

98 continue
write(*, *) '>>> NBI >>> Wrong NBI configuration file format. '
write(*, *) '            More records expected than available.'
stop

990 write(*, *) '>>> NBI >>> read error in dat/srsfi.dat'
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

implicit none

integer, intent(in) :: NCH, NBS
character(len=*), intent(in) :: file_name
integer, intent(out) :: ERCODE

integer :: j, j1, jj, i, ios


double precision :: YCB
character(len=6) :: val_str
character(len=12), dimension(n_fields) :: SFIELD
character(len=80) :: STR, STRI, ADATA(n_nbi_max), CDATA(n_nbi_max)

!----------------------------------------------------------------------|
! Old NBI file: 25xNBS fields to read
!  CBMH1  CBMH2 CBMH3  CBMH4
!  CBMR1  CBMR2
!  CBMS1  CBMS2 CBMS3  CBMS4
!  EBEAM DBM1 DBM2 DBM3
!  ABEAM HBEAM CONTR QBEAMy
!  RBMAX RBMIN YCBEAM YJBEAM
!  SRSNo
! New NBI file:
!   Parameters
!  QBEAM CONTR ABEAM ZBEAM EBEAM
!  DBM1 DBM2 DBM3 CBMS1 CBMS2
!     Geometry
!  HBEAM RBMAX RBMIN CBMS3 CBMS4
!  CBMH1 CBMH2 CBMR1 CBMR2   Not_used
!----------------------------------------------------------------------|

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

99 continue
close(NCH)
ERCODE = 1

end subroutine NQUERY

!---------------------------------------------------------------------
subroutine NBSMTV(YFO, YCI3, YROC, JNA1)
!---------------------------------------------------------------------
!  The subroutine transmits a histogram like function YFO(1:N1)
!  from an N1 grid X(1:N1) to a smooth functon YFO(1:JNA1)
!  keeping the same total volume integral.
! Input: YFO(1:N1)
! Output: YFO(1:JNA1)
!---------------------------------------------------------------------

use numerical_tools, only: smooth
use standard_functions, only: VINT
use nbicom, only: N1, X, XJ, DRI

implicit none

integer, intent(in) :: JNA1
double precision, intent(in) :: YCI3, YROC
double precision, intent(inout) :: YFO(*)

integer :: j, JNA, JNAC, JSIGN
double precision :: Y, YOLD, YNEW, ALFA

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

YOLD = VINT(YFO, YROC)
ALFA = 0.001d0
call SMOOTH(ALFA, x, dri, n1, xj, yfo, jna1)

! Cut of artificial negatives/positive after smoothing
if (JSIGN == 0) then !no real change of sign
    do J=1, JNA1
        if (YOLD > 0d0) then
            if (YFO(J) < 0.d0) YFO(J) =0.d0 !cut negative
        else
            if (YFO(J) > 0.d0) YFO(J) =0.d0 !cut positive
        endif
    enddo
endif

YNEW = VINT(YFO, YROC)
if (dabs(YNEW) > 1.d-19) then
    Y = YOLD/YNEW
    do J=1, JNA1
        YFO(J) = Y*YFO(J)
    enddo
endif

end subroutine NBSMTV

!---------------------------------------------------------------------
subroutine NBSMTS(YFO, YCI3, YROC, JNA1)
!---------------------------------------------------------------------
! Same as NBSMTV, but for surface integrals
!---------------------------------------------------------------------

use numerical_tools, only: smooth
use standard_functions, only: IINT
use nbicom, only: X, XJ, DRI, N1

implicit none

integer, intent(in) :: JNA1
double precision, intent(in) :: YCI3, YROC
double precision, intent(inout) :: YFO(*)

integer :: j, JNA, JNAC, JSIGN
double precision :: Y, YOLD, YNEW, ALFA

JSIGN = 0         ! change of sign? (1= yes, 0= no)
do j=1, n1-1
    JNA = 1 + YCI3*(j - 1)
    JNAC = JNA - 1 + YCI3
    X(j) = XJ(JNAC)
    DRI(J) = YFO(JNA)
    if (JSIGN < 1) then
        if (DRI(1)*DRI(J) < 0.d0) JSIGN = 1
    endif
enddo

dri(N1)  = 0.d0
X(N1)    = 1.d0
XJ(JNA1) = 1.d0

! Total current normalization
YOLD = IINT(YFO, YROC)
ALFA = 0.001d0
call SMOOTH(ALFA, x, dri, n1, xj, yfo, jna1)

! Cut of artificial negatives/positive after smoothing
if (JSIGN == 0) then !no real change of sign
    do J =1, JNA1
        if (YOLD > 0.d0) then
            if (YFO(J) < 0.d0) YFO(J) = 0.d0 !cut negative
        else
            if (YFO(J) > 0.d0) YFO(J) = 0.d0 !cut positive
        endif
    enddo
endif

YNEW = IINT(YFO, YROC)
if (abs(YNEW) > 1.e-19) then
    Y = YOLD/YNEW
    do J=1, JNA1
        YFO(J) = Y*YFO(J)
    enddo
endif

end subroutine NBSMTS

!---------------------------------------------------------------------
double precision function RIPRAD(YUPDWN, J)
! Ripple loss boundary R[m] for each magnetic surface
! YUPDWN [m] plasma midplane shift in respect to the plane
!            of the ripple losses simmetry
! J - magnetic surface number

implicit none

integer, intent(in) :: J
double precision, intent(in) :: YUPDWN

! Default (No ripple losses)
RIPRAD	= 99999.

end function RIPRAD
