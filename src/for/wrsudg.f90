subroutine wrsudg

use parameter_inc, only: NRD
use const_inc, only: ROC, NA1, IPL, TIME, TAU, TSTART, IPEQL, NEQUIL, MEQUIL, IPLFBE,&
    CRAD4, CMHD2, CHE3, CSCL4, CDWM1, CDWM2, CDWM7, &
    CDJM1, CDJM2, CDJM3, CDJM4, CDJM6, CDJM7, &
    CDMJ5, CDMJ6, CDMJ7, CDHJ7, CNEUT1, CNEUT2, &
    ZRD15, ZRD77, ZRD78, ZRD93, CRAD3
use status_inc, only: NE, TE, TI, ZEF, HE, XI, PE, PI, NIBM, PRAD, &
    MU, CU, UPL, F4, CAR34, CAR54, SHIF
use outcmn_inc, only: machine, ccoil
use fenix_params, only: ipl_bf_bkdw
use flight_sim_geometrics, only: geom1d
use plasma_state, only: plasma_up
use parameters_a2equil, only: equil_now

implicit none

integer :: mmequi, j, ntetap, neqlp, i_error
integer, dimension(6) :: gapnum
double precision :: qalp1, qalp2, qrad1, qrad2, QRADR, QTOKR, QDTR, QEDWTR, QIDWTR, &
    WTOZR, LI3R, BETP3R
double precision, dimension(NA1) :: ne1, ne2
double precision, dimension(NRD) :: te_0, te_1, te_now, ne_0, ne_1, ne_now
double precision, dimension(10) :: ccoil_0, ccoil_1, ccoil_now
double precision, dimension(32) :: geom1d_0, geom1d_1, geom1d_now
double precision, dimension(500) :: magnetics
double precision, dimension(660) :: yroutfull
double precision, dimension(256, 256) :: yrout, yzout
double precision, dimension(10) :: ccoil_scramble

double precision coil_forces(na1,2)
double precision betp3r_0, LI3R_0, IPL_0, upl_0, &
    QRADR_0, QTOKR_0, &
    CRAD4_0, CMHD2_0, CSCL4_0, CDWM1_0, &
    CDWM2_0, CDJM6_0, CDJM1_0, &
    CDJM2_0, CDJM4_0, &
    CDJM7_0, CHE3_0, CDJM3_0, ZRD77_0, &
    CDMJ5_0, CDMJ6_0, CDMJ7_0, &
    CDHJ7_0, &
    time_0
double precision betp3r_1, LI3R_1, IPL_1, upl_1, &
    QRADR_1, QTOKR_1, &
    CRAD4_1, CMHD2_1, CSCL4_1, CDWM1_1, &
    CDWM2_1, CDJM6_1, CDJM1_1, &
    CDJM2_1, CDJM4_1, &
    CDJM7_1, CHE3_1, CDJM3_1, ZRD77_1, &
    CDMJ5_1, CDMJ6_1, CDMJ7_1, &
    CDHJ7_1, &
    time_1, dum1, dum2
double precision betp3r_now, LI3R_now, IPL_now, upl_now, &
    QRADR_now, QTOKR_now, &
    CRAD4_now, CMHD2_now, CSCL4_now, CDWM1_now, &
    CDWM2_now, CDJM6_now, CDJM1_now, &
    CDJM2_now, CDJM4_now, &
    CDJM7_now, CHE3_now, CDJM3_now, ZRD77_now, &
    CDMJ5_now, CDMJ6_now, CDMJ7_now, &
    CDHJ7_now, &
    time_now, time_ext, dt_smlk

save betp3r_0, LI3R_0, IPL_0, upl_0, &
    QRADR_0, QTOKR_0, &
    CRAD4_0, CMHD2_0, CSCL4_0, CDWM1_0, &
    CDWM2_0, CDJM6_0, CDJM1_0, &
    CDJM2_0, CDJM4_0, &
    CDJM7_0, CHE3_0, CDJM3_0, ZRD77_0, &
    CDMJ5_0, CDMJ6_0, CDMJ7_0, &
    CDHJ7_0, &
    te_0, ne_0, &
    ccoil_0, &
    geom1d_0,time_0,betp3r_1, LI3R_1, IPL_1, upl_1, &
    QRADR_1, QTOKR_1, &
    CRAD4_1, CMHD2_1, CSCL4_1, CDWM1_1, &
    CDWM2_1, CDJM6_1, CDJM1_1, &
    CDJM2_1, CDJM4_1, &
    CDJM7_1, CHE3_1, CDJM3_1, ZRD77_1, &
    CDMJ5_1, CDMJ6_1, CDMJ7_1, &
    CDHJ7_1, &
    te_1, ne_1, &
    ccoil_1, &
    geom1d_1,time_1,betp3r_now, LI3R_now, IPL_now, upl_now, &
    QRADR_now, QTOKR_now, &
    CRAD4_now, CMHD2_now, CSCL4_now, CDWM1_now, &
    CDWM2_now, CDJM6_now, CDJM1_now, &
    CDJM2_now, CDJM4_now, &
    CDJM7_now, CHE3_now, CDJM3_now, ZRD77_now, &
    CDMJ5_now, CDMJ6_now, CDMJ7_now, &
    CDHJ7_now, &
    te_now, ne_now, &
    ccoil_now, &
    geom1d_now, time_now


if (MACHINE(1:3) == 'dem') then
    i_error = 0 ! whether add noise latencies errors to diagnostics
    qalp1 = QDTR(ROC)
    qrad1 = QRADR(ROC)
    NE1 = NE(1:NA1)
    if (i_error /= 0) then
        call add_errors(qalp1, qrad1, ne1, qalp2, qrad2, ne2)
! QDTR: fusion power
! QRAD: radiated power
! NE: electron density profile
        qalp1 = qalp2
        qrad1 = qrad2
        ne1 = ne2
    endif

    gapnum(1) = 16
    gapnum(2) = 23
    gapnum(3) = 30
    gapnum(4) = 36
    gapnum(5) = 46
    gapnum(6) = 47
    do j=1,6
        magnetics(j) = geom1d(94 - 52 +1 +gapnum(j) - 1)
    enddo
    magnetics(7) = geom1d(97) !Rcurr
    magnetics(8) = geom1d(98) !Zcurr
    magnetics(9) = ipl*1.e6   !Ipl A

    call shmw( &
        TIME, &
        qalp1, qrad1, &
        IPL, QTOKR(ROC) - qedwtr(ROC) - qidwtr(ROC), &
        CRAD4, CMHD2, CSCL4, CDWM1, &
        CDWM2, CDWM7, CDJM1, CDJM2, CDJM3, CDJM4, ZRD77, &
        CDMJ5, CDMJ6, CDMJ7, &
        CNEUT1, CNEUT2, TE(1:NA1), NE1(1:NA1), &
        1.d0/MU(1:NA1), &
        CU(1:NA1), &
        CCOIL(1:15)*1.e3, magnetics(1:439), & 
        CRAD3, F4(1:NA1), ZEF(1:NA1), TI(1:NA1), SHIF(1:NA1))

elseif (MACHINE(1:4) == 'iter') then
    i_error = 0 ! whether add noise latencies errors to diagnostics
    qalp1 = QDTR(ROC)
    qrad1 = QRADR(ROC)
    NE1 = NE(1:NA1)

! QDTR: fusion power
! QRAD: radiated power
! NE: electron density profile
    qalp1 = qalp2
    qrad1 = qrad2
    ne1 = ne2

!    gapnum(1) = 16
!    gapnum(2) = 23
!    gapnum(3) = 30
!    gapnum(4) = 36
!    gapnum(5) = 46
!    gapnum(6) = 47
!    do j=1,6
!        magnetics(j) = geom1d(94 - 52 +1 +gapnum(j) - 1)
!    enddo
!    magnetics(7) = geom1d(97) !Rcurr
!    magnetics(8) = geom1d(98) !Zcurr!
!    magnetics(9) = ipl*1.e6   !Ipl A

    call shmw( &
        TIME, &
        qalp1, qrad1, &
        IPL, QTOKR(ROC) - qedwtr(ROC) - qidwtr(ROC), &
        CRAD4, CMHD2, CSCL4, CDWM1, &
        CDWM2, CDWM7, CDJM1, CDJM2, CDJM3, CDJM4, ZRD77, &
        CDMJ5, CDMJ6, CDMJ7, &
        CNEUT1, CNEUT2, TE(1:NA1), NE1(1:NA1), &
        1.d0/MU(1:NA1), &
        CU(1:NA1), &
        CCOIL(1:15)*1.e3, magnetics(1:439))

else if (MACHINE(1:3) == 'aug') then

    ccoil_scramble(1: 10) = ccoil(1: 10)
    ccoil(2) = ccoil_scramble(2) - ccoil_scramble(1)
    ccoil(3) = ccoil_scramble(3) - ccoil_scramble(2)

!scramble coil currents

    geom1d(64) = WTOZR(ROC)*1e6 !total Wmhd including fast ions  in MJ
    time_1 = time - tstart
    time_0 = time - tstart - TAU

! get coil forces
    if (nint(IPEQL) == 4) then
        call coil_force2(coil_forces(1:100, 1:2), plasma_up)
    else
        call coil_forces_feqis(100, coil_forces(1:100, 1), coil_forces(1:100, 2), plasma_up)
        coil_forces = -coil_forces
    endif

    car54( 1:21) = coil_forces(1:21, 1)
    car54(22:42) = coil_forces(1:21, 2)

    write(8998, '(333E25.11)') TIME, coil_forces(1:21, 1), coil_forces(1:21, 2)

    neqlp  = abs(nint(NEQUIL))
    ntetap = abs(nint(MEQUIL)) + 1
    yrout(1:neqlp, 1:ntetap-1) = equil_now%coord_sys%position%r(1:neqlp, 1:ntetap)
    yzout(1:neqlp, 1:ntetap-1) = equil_now%coord_sys%position%z(1:neqlp, 1:ntetap)
    yrout(1:neqlp, ntetap  ) = yrout(1:neqlp, 1)
    yzout(1:neqlp, ntetap  ) = yzout(1:neqlp, 1)
    yrout(1:neqlp, ntetap+1) = yrout(1:neqlp, 2)
    yzout(1:neqlp, ntetap+1) = yzout(1:neqlp, 2)
    mmequi = 2*ntetap
    do j=1,ntetap
        yroutfull(2*j-1) = yrout(neqlp, j)
	yroutfull(2*j  ) = yzout(neqlp, j)
    enddo
    if (TIME - TSTART <= ZRD93 +1.e-8) then
        if (plasma_up == 0) then
            call shmw(TIME - TSTART, &
                betp3r(roc), LI3R(ROC), ipl_bf_bkdw, upl(na1), &
                QRADR(ROC), QTOKR(ROC), &
                CRAD4, CMHD2, CSCL4, CDWM1, &
                ZRD15, CDJM6, CDJM1, &
                CDJM4, CDJM7, CHE3, CDMJ5, CDMJ6, CDMJ7, &
                te(1:na1), ne(1:na1), ZEF(1:NA1), &
                F4(1:NA1), &
                coil_forces(1:na1, 1), &
                coil_forces(1:na1, 2), &
                TI(1:na1), XI(1:na1), HE(1:na1), &
                PE(1:na1), PI(1:na1), &
                NIBM(1:na1), &
                PRAD(1:na1), &
                ccoil(1:10), &
                geom1d(51:82), &
                car34(1:12), &
                yroutfull(1:mmequi), &
                CDHJ7)
        else
            call shmw(TIME - TSTART, &
                betp3r(roc), LI3R(ROC), IPLFBE, upl(na1), &
                QRADR(ROC), QTOKR(ROC), &
                CRAD4, CMHD2, CSCL4, CDWM1, &
                ZRD15, CDJM6, CDJM1, &
                CDJM4, &
                CDJM7, CHE3, & 
                CDMJ5, CDMJ6, CDMJ7, &
                te(1:na1), ne(1:na1), &
                ZEF(1:NA1), &
                F4(1:NA1), &
                coil_forces(1:na1, 1), &
                coil_forces(1:na1, 2), &
                TI(1:na1), &
                XI(1:na1), &
                HE(1:na1), &
                PE(1:na1), &
                PI(1:na1), &
                NIBM(1:na1), &
                PRAD(1:na1), &
                ccoil(1:10), &
                geom1d(51:82), &
                car34(1:12), &
                yroutfull(1:mmequi), &
                cdhj7)
        endif

        betp3r_0 = betp3r(roc)
        LI3R_0 = LI3R(ROC)
        if (plasma_up == 0) then
            IPL_0 = ipl_bf_bkdw
        else
            IPL_0 = IPL
        endif
        upl_0 = upl(na1)
        QRADR_0 = QRADR(ROC)
        QTOKR_0 = QTOKR(ROC)
        CRAD4_0 = CRAD4
        CMHD2_0 = CMHD2
        CSCL4_0 = CSCL4
        CDWM1_0 = CDWM1
        CDWM2_0 = CDWM2
        CDJM6_0 = CDJM6
        CDJM1_0 = CDJM1
        CDJM2_0 = CDJM2
        CDJM4_0 = CDJM4
        CDJM7_0 = CDJM7
        CHE3_0 = CHE3
        CDJM3_0 = CDJM3
        ZRD77_0 = ZRD77
        CDMJ5_0 = CDMJ5
        CDMJ6_0 = CDMJ6
        CDMJ7_0 = CDMJ7
        CDHJ7_0 = CDHJ7
        te_0(1:na1) = te(1:na1)
        ne_0(1:na1) = ne(1:na1)
        ccoil_0 = ccoil(1:10)
        geom1d_0 = geom1d(51:82)

    else

200	continue

        betp3r_1 = betp3r(roc)
        LI3R_1 = LI3R(ROC)
        if (plasma_up == 0) then
            IPL_1 = ipl_bf_bkdw
        else
            IPL_1 = IPL
        endif
        upl_1 = upl(na1)
        QRADR_1 = QRADR(ROC)
        QTOKR_1 = QTOKR(ROC)
        CRAD4_1 = CRAD4
        CMHD2_1 = CMHD2
        CSCL4_1 = CSCL4
        CDWM1_1 = CDWM1
        CDWM2_1 = CDWM2
        CDJM6_1 = CDJM6
        CDJM1_1 = CDJM1
        CDJM2_1 = CDJM2
        CDJM4_1 = CDJM4
        CDJM7_1 = CDJM7
        CHE3_1 = CHE3
        CDJM3_1 = CDJM3
        ZRD77_1 = ZRD77
        CDMJ5_1 = CDMJ5
        CDMJ6_1 = CDMJ6
        CDMJ7_1 = CDMJ7
        CDHJ7_1 = CDHJ7
        te_1(1:na1) = te(1:na1)
        ne_1(1:na1) = ne(1:na1)
        ccoil_1 = ccoil(1:10)
        geom1d_1 = geom1d(51:82)

        time_now = ZRD93
        dum1 = (time_1   - time_0)**0.5
        dum2 = (time_now - time_0)**0.5
        betp3r_now = betp3R_0 + (betp3R_1 - betp3R_0)/dum1*dum2
        LI3R_now   = LI3R_0 + (LI3R_1 - LI3R_0)/dum1*dum2
        IPL_now    = IPL_0  + ( IPL_1 -  IPL_0)/dum1*dum2
        upl_now    = upl_0  + ( upl_1 -  upl_0)/dum1*dum2
        QRADR_now = QRADR_0 + (qradr_1 - qradr_0)/dum1*dum2
        QTOKR_now = QTOKR_0 + (QTOKR_1 - qtokr_0)/dum1*dum2
        CRAD4_now = CRAD4_0 + (CRAD4_1 - crad4_0)/dum1*dum2
        CMHD2_now = CMHD2_0 + (CMHD2_1 - CMHD2_0)/dum1*dum2
        CSCL4_now = CSCL4_0 + (CSCL4_1 - CSCL4_0)/dum1*dum2
        CDWM1_now = CDWM1_0 + (CDWM1_1 - CDWM1_0)/dum1*dum2
        CDWM2_now = CDWM2_0 + (CDWM2_1 - CDWM2_0)/dum1*dum2
        CDJM6_now = CDJM6_0 + (CDJM6_1 - CDJM6_0)/dum1*dum2
        CDJM1_now = CDJM1_0 + (CDJM1_1 - CDJM1_0)/dum1*dum2
        CDJM2_now = CDJM2_0 + (CDJM2_1 - CDJM2_0)/dum1*dum2
        CDJM4_now = CDJM4_0 + (CDJM4_1 - CDJM4_0)/dum1*dum2
        CDJM7_now = CDJM7_0 + (CDJM7_1 - CDJM7_0)/dum1*dum2
        CHE3_now  = CHE3_0  + ( CHE3_1 -  CHE3_0)/dum1*dum2
        CDJM3_now = CDJM3_0 + (CDJM3_1 - CDJM3_0)/dum1*dum2
        ZRD77_now = ZRD77_0 + (ZRD77_1 - ZRD77_0)/dum1*dum2
        CDMJ5_now = CDMJ5_0 + (CDMJ5_1 - CDMJ5_0)/dum1*dum2
        CDMJ6_now = CDMJ6_0 + (CDMJ6_1 - CDMJ6_0)/dum1*dum2
        CDMJ7_now = CDMJ7_0 + (CDMJ7_1 - CDMJ7_0)/dum1*dum2
        CDHJ7_now = CDHJ7_0 + (CDHJ7_1 - CDHJ7_0)/dum1*dum2
        te_now(1:na1) = te_0(1:na1) + (te_1(1:na1) - te_0(1:na1))/dum1*dum2
        ne_now(1:na1) = ne_0(1:na1) + (ne_1(1:na1) - ne_0(1:na1))/dum1*dum2
        ccoil_now  = ccoil_0  + ( ccoil_1 -  ccoil_0)/dum1*dum2
        geom1d_now = geom1d_0 + (geom1d_1 - geom1d_0)/dum1*dum2

        call shmw(time_now, &
            betp3r_now, LI3R_now, IPL_now, upl_now, &
            QRADR_now, QTOKR_now, &
            CRAD4_now, CMHD2_now, CSCL4_now, CDWM1_now, &
            ZRD15, CDJM6_now, CDJM1_now, &
            CDJM4_now, &
            CDJM7_now, CHE3_now, &
            CDMJ5_now, CDMJ6_now, CDMJ7_now, &
            te_now(1:na1), ne_now(1:na1), &
            ZEF(1:NA1), &
            F4(1:NA1), &
            coil_forces(1:na1, 1), &
            coil_forces(1:na1, 2), &
            TI(1:na1), &
            XI(1:na1), &
            HE(1:na1), &
            PE(1:na1), &
            PI(1:na1), &
            NIBM(1:na1), &
            PRAD(1:na1), &
            ccoil(1:10), &
            geom1d(51:82), &
            car34(1:12), &
            yroutfull(1:mmequi), &
            cdhj7) !na1+1.d-16  !& geom1d(1:10),

        if (TIME - TSTART.le. ZRD93 + 1.e-8) goto 300
            call read_input_constant_file(time_ext, dt_smlk)
            ZRD93 = time_ext + dt_smlk
            goto 200
        endif

300 continue
    betp3r_0 = betp3r(roc)
    LI3R_0 = LI3R(ROC)
    if (plasma_up == 0) then
        IPL_0 = ipl_bf_bkdw
    else
        IPL_0 = IPL
    endif
    upl_0 = upl(na1)
    QRADR_0 = QRADR(ROC)
    QTOKR_0 = QTOKR(ROC)
    CRAD4_0 = CRAD4
    CMHD2_0 = CMHD2
    CSCL4_0 = CSCL4
    CDWM1_0 = CDWM1
    CDWM2_0 = CDWM2
    CDJM6_0 = CDJM6
    CDJM1_0 = CDJM1
    CDJM2_0 = CDJM2
    CDJM4_0 = CDJM4
    CDJM7_0 = CDJM7
    CHE3_0 = CHE3
    CDJM3_0 = CDJM3
    ZRD77_0 = ZRD77
    CDMJ5_0 = CDMJ5
    CDMJ6_0 = CDMJ6
    CDMJ7_0 = CDMJ7
    CDHJ7_0 = CDHJ7
    te_0(1:na1) = te(1:na1)
    ne_0(1:na1) = ne(1:na1)
    ccoil_0 = ccoil(1:10)
    geom1d_0 = geom1d(51:82)


	ccoil(1:10)=ccoil_scramble(1:10)

endif

return
end subroutine wrsudg

!-------------------------------------------------------------------------
subroutine add_errors(qalp1, qrad1, ne1, qalp2, qrad2, ne2)  !for demo

use const_inc, only: NA1, TIME
use status_inc, only: NE

implicit none

double precision, intent(in) :: qalp1, qrad1
double precision, intent(in), dimension(na1) :: ne1
double precision, intent(out) :: qalp2, qrad2
double precision, intent(out), dimension(na1) :: ne2

integer :: jrho
double precision :: dum1, dum2, dum3, dum4, dum5, dum11, dum12, &
    dum21, dum22, dum31, dum41, dum51, dum52
double precision, dimension(NA1) :: dum32, dum42

save dum4, dum5

qalp2 = qalp1
qrad2 = qrad1
ne2 = ne1

!add error on alpha power
dum1 = 0.01 !error on alpha power reconstruction
dum2 = 0.4 ! integration time in s
if (TIME - dum5 >= dum2) then
    dum5 = TIME
    dum4 = qalp1
endif
qalp2 = dum4
call random_number(dum3)
qalp2 = qalp2*(1. + dum1*(1. - 2.*dum3))
!add error on radiated power
dum11 = 0.1 !error on alpha power reconstruction
dum21 = 0.5 ! integration time in s
if (TIME - dum51 >= dum21) then
    dum51 = TIME
    dum41 = qrad1
endif
    qrad2 = dum41
call random_number(dum31)
qrad2 = qrad2*(1. + dum11*(1. - 2.*dum31))
!add error on density power
dum12 = 0.02 !error on alpha power reconstruction
dum22 = 0.001 ! integration time in s
if (TIME - dum52 >= dum22) then
    dum52 = TIME
    dum42 = ne(1:na1)
endif
    ne2 = dum42
do jrho=1, na1
    call random_number(dum32(jrho))
enddo
ne2 = ne2*(1. + dum12*(1. - 2.*dum32))

return
end subroutine add_errors
