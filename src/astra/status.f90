module status

implicit none

integer, parameter :: NRD=801

double precision, dimension(:), pointer :: &
    TEO, TIO, NEO, UPAR, UPARO, &
    FPO, ULON, QU, VP, GN, SQEPS, &
    PRES, PELH, PETOT, PELON, CULH, PITOT, &
    SNTOT, PEECR, PEFW, PEICR, PEIQI, PIICR, PIFW, PEPER, &
    QE, QI, QN, QF0, QF1, QF2, QF3, QF4, QF5, QF6, QF7, &
    QF8, QF9, CUECR, CUFW, CUICR, CUTOR, TTRQ, VR, VRO, NIO, &
    SHIV, CUBS, &
    VPFP, G22, G33, G41, G42, G43, G44, G45, &
    VOLUM, CV, DRODA, PDE, PDI, SDN, &
    SD0, SD1, SD2, SD3, SD4, SD5, SD6, SD7, SD8, SD9, &
    UPS0, UPS0O, UPS1, UPS1O, DLNEO, SGNEO, UPS2, &
    UPS2O, MRHO, DDNEO, SGNEOD, DLNEOD, SQUARN, &
    B0DB2, BDB02, BDB0, BMAXT, BMINT, FOFB, GRADRO, &
    EQFF, EQPF, SLAT, FV, MV, XRHO, rho_pol, FP_NORM, &
    SXHO, SRHO, G22E, G33E, AREAT, PERIM, &
    AMAIN, AMETR, CU, ELON, ER, FP, G11, IPOL, MU, &
    NALF, NDEUT, NE, NHE3, NHYDR, NI, NIBM, NIZ1, NIZ2, NIZ3, NTRIT, &
    PBLON, PBPER, PFAST, RHO, SHEAR, SHIF, TE, TI, TRIA, UPL, VPOL, VRS, VTOR, &
    ZEF, ZIM1, ZIM2, ZIM3, ZMAIN, ZIMPT, NIMPT, AIMPT, &
    NN, TN, PRAD, PBOL1, PBOL2, PBOL3, PSXR1, &
    PSXR2, PSXR3, ZEF1, ZEF2, ZEF3, &
    DIMP1, DIMP2, DIMP3, VIMP1, VIMP2, VIMP3, NMAIN, &
    PBEAM, SNNBM, SNEBM, NNBM1, NNBM2, NNBM3, &
    SCUBM, CUBM, PEBM, PIBM, &
    CUFI, SNIBM1, SNIBM2, SNIBM3, NRATE, &
    CC,  CD,  CE,  CI,  CN, &
    DC,  DE,  DI,  DN,  HC, &
    HE,  HI,  HN,  PE,  PET, &
    PI,  PIT, SN,  SNN, XC, &
    XE,  XI,  XN,  DVN, DSN, &
    DVE, DSE, DVI, DSI, CNPAR, &
    CNPAP, XUPAR, RUPAR, RUPYR, RUPFR, &
    CNPAD, XUPAP, TTRQI, XUPAD, DDNEOD

double precision, dimension(:), pointer :: &
    CAR1,  CAR2,  CAR3,  CAR4,  CAR5 , CAR6,  CAR7,  CAR8,  CAR9,  CAR10, &
    CAR11, CAR12, CAR13, CAR14, CAR15, CAR16, CAR17, CAR18, CAR19, CAR20, &
    CAR21, CAR22, CAR23, CAR24, CAR25, CAR26, CAR27, CAR28, CAR29, CAR30, &
    CAR31, CAR32, CAR33, CAR34, CAR35, CAR36, CAR37, CAR38, CAR39, CAR40, &
    CAR41, CAR42, CAR43, CAR44, CAR45, CAR46, CAR47, CAR48, CAR49, CAR50, &
    CAR51, CAR52, CAR53, CAR54, CAR55, CAR56, CAR57, CAR58, CAR59, CAR60, &
    CAR61, CAR62, CAR63, CAR64, CAR65, CAR66, CAR67, CAR68, CAR69, CAR70, &
    CAR71, CAR72, CAR73, CAR74, CAR75, CAR76, CAR77, CAR78, CAR79, CAR80, &
    CAR81, CAR82, CAR83, CAR84, CAR85, CAR86, CAR87, CAR88, CAR89, CAR90, &
    CAR91, CAR92, CAR93, CAR94, CAR95, CAR96, CAR97, CAR98, CAR99, CAR100, &
    CAR101, CAR102, CAR103, CAR104, CAR105, CAR106, CAR107, CAR108, CAR109, CAR110, &
    CAR111, CAR112, CAR113, CAR114, CAR115, CAR116, CAR117, CAR118, CAR119, CAR120, &
    CAR121, CAR122, CAR123, CAR124, CAR125, CAR126, CAR127, CAR128

double precision, dimension(:), pointer :: &
    CAR1X,  CAR2X,  CAR3X,  CAR4X,  CAR5X,  CAR6X,  CAR7X,  CAR8X,  CAR9X,  CAR10X, &
    CAR11X, CAR12X, CAR13X, CAR14X, CAR15X, CAR16X, CAR17X, CAR18X, CAR19X, CAR20X, &
    CAR21X, CAR22X, CAR23X, CAR24X, CAR25X, CAR26X, CAR27X, CAR28X, CAR29X, CAR30X, &
    CAR31X, CAR32X, CAR33X, CAR34X, CAR35X, CAR36X, CAR37X, CAR38X, CAR39X, CAR40X, &
    CAR41X, CAR42X, CAR43X, CAR44X, CAR45X, CAR46X, CAR47X, CAR48X, CAR49X, CAR50X, &
    CAR51X, CAR52X, CAR53X, CAR54X, CAR55X, CAR56X, CAR57X, CAR58X, CAR59X, CAR60X, &
    CAR61X, CAR62X, CAR63X, CAR64X, &
    F0X, F1X, F2X, F3X, F4X, F5X, F6X, F7X, F8X, F9X, &
    MUX, MVX, GNX, SNX, PEX, PIX, PRADX, TEX, TIX, NEX, CUX, &
    ZEFX, VRX, SHX, ELX, TRX, G11X, G22X, G33X, DRODAX, IPOLX, &
    NIX, VPOLX, VTORX, SLATX, SHIVX, SQUAX

double precision, allocatable, dimension(:, :), target :: profiles_x, profiles

double precision, dimension(:), pointer :: &
    F0, F1, F2, F3, F4, F5, F6, F7, F8, F9, &
    F0O, F1O, F2O, F3O, F4O, F5O, F6O, F7O, F8O, F9O, &
    DF0, DF1, DF2, DF3, DF4, DF5, DF6, DF7, DF8, DF9, &
    VF0, VF1, VF2, VF3, VF4, VF5, VF6, VF7, VF8, VF9, &
    GF0, GF1, GF2, GF3, GF4, GF5, GF6, GF7, GF8, GF9, &
    SF0, SF1, SF2, SF3, SF4, SF5, SF6, SF7, SF8, SF9, &
    DVF0, DVF1, DVF2, DVF3, DVF4, DVF5, DVF6, DVF7, DVF8, DVF9, &
    DSF0, DSF1, DSF2, DSF3, DSF4, DSF5, DSF6, DSF7, DSF8, DSF9, &
    SFF0, SFF1, SFF2, SFF3, SFF4, SFF5, SFF6, SFF7, SFF8, SFF9, &
    SF0TOT, SF1TOT, SF2TOT, SF3TOT, SF4TOT, SF5TOT, SF6TOT, SF7TOT, SF8TOT, SF9TOT

contains

!---------------------------------------------------------------------
    subroutine status_init()

    integer :: j

! Defaults, rather fall-back than initial values

    profiles_x = 0.d0
    profiles   = 0.d0

! Geometry/equilibrium

    ELON = 1.
    G33 = 1.
    G41 = 1.
    G42 = 1.
    G43 = 0.2
    G44 = 1.
    G45 = 1.
    G22E = 1.
    G33E = 1.
    SQEPS = 0.2 
    AREAT = 1.
    PERIM = 1.
    VOLUM = 1.

! Species

    ZEF = 1.
    AMAIN = 1.
    ZIMPT = 1.
    AIMPT = 1.
    ZMAIN = 1.
    ZIM1  = 1.
    ZIM2  = 1.
    ZIM3  = 1.
    TE = 0.01
    TI = 0.01
    TN = 0.01
    NE = 0.1
    NI = 0.1
    NN = 1.d-5

! Current

    CU = 0.1
    MU = 0.3
    IPOL = 1.

! Sources

    PETOT = 1.d-4
    PITOT = 1.d-4

    MRHO = 1.d0
    F0 = 1.d0
    F1 = 1.d0
    F2 = 1.d0
    F3 = 1.d0
    F4 = 1.d0
    F5 = 1.d0
    F6 = 1.d0
    F7 = 1.d0
    F8 = 1.d0
    F9 = 1.d0

    do j=1, NRD
       FP(j) = 0.1*j
       QE(j) = FP(j)
       QU(j) = FP(j)
       QI(j) = FP(j)
       QN(j) = FP(j)
       GN(j) = FP(j)
    enddo

    end subroutine status_init

!---------------------------------------------------------------------
    subroutine DEFARR()
!---------------------------------------------------------------------
! 1) Check positiveness of Z_eff, n_e, n_i, T_e, T_i.  Stop if negative.
! 2) Extend definition of all standard arrays beyond NA1.
!Efable
! 3) Define upsilons for momentum transport equation
!---------------------------------------------------------------------

    use pi_const, only: GP2
    use read_input, only: exp_file
    use scalars, only: RTOR, BTOR, HRO, ROC, ABC, NA1, NB1, &
        TIME, TAU, TSTART, WTE, WTI, WNE
    use debugger, only: markloc

    integer :: j
    double precision :: YV, YF, YMU, YN, YNE, YNI, YTE, YTI, YZF
    character(len=2) :: str
    character(len=132) :: err_msg

    call markloc('DEFARR')

    YZF = 0.
    YNE = 0.
    YNI = 0.
    YTE = 0.
    YTI = 0.
    YV = GP2*BTOR*HRO
    do j=1, NA1
        SQEPS(j) = SQRT(AMETR(j)/(RTOR + SHIF(j)))
! VP = c<E_||/h>/B_p in Hinton & Hazeltine Rev.Mod.Phys.48(1976) p.297
        if (MU(j) <= 0.d0) then
            write(*,'(/A,F10.6,A)')'>>> ERROR >>>  Time =', TIME, ' sec'
            write(*,'(A/A,F10.6)') &
                '               The rotational transform is less or equal zero ', &
                '               at rho_N =', RHO(j)/ROC
            if (TIME <= TSTART + TAU/2.) write(*, '(2A,1H"/)') &
                '               Check if it is defined in the data file "', &
                TRIM(exp_file)
            call error_catch()
            STOP
        endif
        VP(j) = ULON(j)/(YV*j*MU(j))
        YZF = max(YZF, ZEF(j))
        YNE = max(YNE, NE(j))
        YNI = max(YNI, NI(j))
        YTE = max(YTE, TE(j))
        YTI = max(YTI, TI(j))
    enddo

    str = '3M'
    if (min(YZF, YNE, YNI, YTE, YTI) <= 0.) then

        if (YZF <= 0.) then
            str = "ZF"
        elseif (YNE <= 0.) then
            str = "NE"
        elseif (YNI <= 0.) then
            str = "NI"
        elseif (YTE <= 0.) then
            str = "TE"
        elseif (YTI <= 0.) then
            str = "TI"
        endif

        write(*, '(/A, F6.3, A)') '>>> ERROR >>>  Time =',TIME,'   DEFVAR'

        SELECT CASE(str)
        CASE('3M')
            err_msg = '               The 3M equilibrium solver does not converge'
        CASE('ZF')
            err_msg = '               The variable "ZEF" is less or equal zero'
        CASE DEFAULT
            err_msg = '               The variable "' // str // '" is less or equal zero'
        END SELECT
        STOP err_msg

    endif

    if (NA1 /= NB1) then

! Determine arrays beyond ABC
        YF = GP2*BTOR
        YMU = (MU(NA1) - MV(NA1))*ROC**2

        do j=NA1+1, NB1
            FPO(j)   = FPO(NA1)
            MU(j)    = MU(NA1)
            ULON(j)  = ULON(NA1)
            TTRQ(j)  = TTRQ(NA1)
            CU(j)    = CU(NA1)
            VP(j)    = VP(NA1)
            GN(j)    = GN(NA1)
            SQEPS(j) = SQEPS(NA1)
            PELH(j)  = PELH(NA1)
            PETOT(j) = PETOT(NA1)
            PELON(j) = PELON(NA1)
            CULH(j)  = CULH(NA1)
            PITOT(j) = PITOT(NA1)
            SNTOT(j) = SNTOT(NA1)
            ZEF(j)   = ZEF(NA1)
            PEECR(j) = PEECR(NA1)
            PEFW(j)  = PEFW(NA1)
            PEICR(j) = PEICR(NA1)
            PIICR(j) = PIICR(NA1)
            PIFW(j)  = PIFW(NA1)
            PEPER(j) = PEPER(NA1)

            QE(j)    = QE(NA1)
            QI(j)    = QI(NA1)
            QN(j)    = QN(NA1)
            QF0(j)   = QF0(NA1)
            QF1(j)   = QF1(NA1)
            QF2(j)   = QF2(NA1)
            QF3(j)   = QF3(NA1)
            QF4(j)   = QF4(NA1)
            QF5(j)   = QF5(NA1)
            QF6(j)   = QF6(NA1)
            QF7(j)   = QF7(NA1)
            QF8(j)   = QF8(NA1)
            QF9(j)   = QF9(NA1)

            CUECR(j) = CUECR(NA1)
            CUFW(j)  = CUFW(NA1)
            CUICR(j) = CUICR(NA1)
            CUTOR(j) = CUTOR(NA1)
            UPL(j)   = UPL(NA1)
            ER(j)    = ER(NA1)
            VRS(j)   = VRS(NA1)
            NI(j)    = NI(NA1)
            NIO(j)   = NIO(NA1)
            ZMAIN(j) = ZMAIN(NA1)
            AMAIN(j) = AMAIN(NA1)

            VRO(j)   = VRO(NA1)
            VR(j)    = VR(NA1)
            SHIF(j)  = SHIF(NA1)
            SHIV(j)  = SHIV(NA1)
            SQUARN(j)= SQUARN(NA1)
            ELON(j)  = ELON(NA1)
            TRIA(j)  = TRIA(NA1)
            G11(j)   = G11(NA1)
            G22(j)   = G22(NA1)
            G33(j)   = G33(NA1)

            G41(j)   = G41(NA1)
            G42(j)   = G42(NA1)
            G43(j)   = G43(NA1)
            G44(j)   = G44(NA1)
            G45(j)   = G45(NA1)

            VPFP(j)  = VPFP(NA1)
            IPOL(j)  = IPOL(NA1)
            CUBS(j)  = CUBS(NA1)
            VOLUM(j) = VOLUM(NA1)
            CV(j)    = CV(NA1)
            DRODA(j) = DRODA(NA1)
            PRES(j)  = PRES(NA1)
            PDE(j)   = PDE(NA1)
            PDI(j)   = PDI(NA1)

            SDN(j)   = SDN(NA1)
            SD0(j)   = SD0(NA1)
            SD1(j)   = SD1(NA1)
            SD2(j)   = SD2(NA1)
            SD3(j)   = SD3(NA1)
            SD4(j)   = SD4(NA1)
            SD5(j)   = SD5(NA1)
            SD6(j)   = SD6(NA1)
            SD7(j)   = SD7(NA1)
            SD8(j)   = SD8(NA1)
            SD9(j)   = SD9(NA1)

            UPS0(j)  = UPS0(NA1)
            UPS1(j)  = UPS1(NA1)
            UPS0O(j) = UPS0O(NA1)
            UPS1O(j) = UPS1O(NA1)

            NN(j) = NN(NA1)
            TN(j) = TN(NA1)

            NIZ1(j)  = NIZ1(NA1)
            NIZ2(j)  = NIZ2(NA1)
            NIZ3(j)  = NIZ3(NA1)
            NALF(j)  = NALF(NA1)
            PRAD(j)  = PRAD(NA1)
            VTOR(j)  = VTOR(NA1)
            NHYDR(j) = NHYDR(NA1)
            NDEUT(j) = NDEUT(NA1)
            NTRIT(j) = NTRIT(NA1)
            NHE3(j)  = NHE3(NA1)
            ZIM1(j)  = ZIM1(NA1)
            ZIM2(j)  = ZIM2(NA1)
            ZIM3(j)  = ZIM3(NA1)
            PBOL1(j) = PBOL1(NA1)
            PBOL2(j) = PBOL2(NA1)
            PBOL3(j) = PBOL3(NA1)
            PSXR1(j) = PSXR1(NA1)
            PSXR2(j) = PSXR2(NA1)
            PSXR3(j) = PSXR3(NA1)
            ZEF1(j)  = ZEF1(NA1)
            ZEF2(j)  = ZEF2(NA1)
            ZEF3(j)  = ZEF3(NA1)
            DIMP1(j) = DIMP1(NA1)
            DIMP2(j) = DIMP2(NA1)
            DIMP3(j) = DIMP3(NA1)
            VIMP1(j) = VIMP1(NA1)
            VIMP2(j) = VIMP2(NA1)
            VIMP3(j) = VIMP3(NA1)
            VPOL(j)  = VPOL(NA1)
            NMAIN(j) = NMAIN(NA1)

            YN = exp((ABC - AMETR(j))/WNE)
            NE(j)    = NE(NA1)*YN
            NI(j)    = NI(NA1)*YN
            NALF(j)  = NALF(NA1)*YN
            NHE3(j)  = NHE3(NA1)*YN
            NHYDR(j) = NHYDR(NA1)*YN
            NDEUT(j) = NDEUT(NA1)*YN
            NTRIT(j) = NTRIT(NA1)*YN
            NMAIN(j) = NMAIN(NA1)*YN
            TE(j)     = TE(NA1)*exp((ABC - AMETR(j))/WTE)
            TI(j)     = TI(NA1)*exp((ABC - AMETR(j))/WTI)
            UPAR(j) = UPAR(NA1)*exp((ABC - AMETR(j))/WTI)
! Efable: should be changed when a SOL equilibrium is implemented
            MV(j) = MV(NA1)
            FV(j) = FV(j-1) + YF*RHO(j)*MV(j)*(RHO(j) - RHO(j-1))
            MU(j) = MV(j) + YMU/RHO(j)**2
            FP(j) = FP(j-1) + YF*RHO(j)*MU(j)*(RHO(j) - RHO(j-1))
            CU(j) = 0.
            CUBS(j) = 0.
            CV(j) = 0.
            CD(j) = 0.
            CUTOR(j) = 0.
        enddo
    endif

    return
    end subroutine DEFARR

!---------------------------------------------------------------------
    subroutine SETVAR()

    use scalars, only: NA1, AMJ, ZMJ, AB, ABC, SHIFT, IPEQL

    integer :: j

    do j=1, NA1
        ZEF(J)= max(1.d0, ZEFX(J))
        ZMAIN(J) = ZMJ
        AMAIN(J) = AMJ
        NI(J) = NE(J)/ZMJ
    enddo
    if (ABC+abs(SHIFT) > AB) then
        write(*, *) char(7), ">>> Warning >>> Inconsistent boundary setting."
        if (IPEQL == 3) then
            write(*, *) "    Plasma beyond the vacuum vessel has been cut off"
        else
            write(*, *) "    Plasma boundary intersects the vacuum vessel"
        endif
    endif

    end subroutine setvar

!---------------------------------------------------------------------
    subroutine error_catch()

    use scalars, only: NA1

    write(*, *) 'TE    Fp    NE    G11 '
    write(*, *) te(1)  , fp(1)  , ne(1)  , g11(1)
    write(*, *) te(na1), fp(na1), ne(na1), g11(na1)
    write(*, *) 'somethings not right, quit run'
    stop

    return
    end subroutine error_catch

end module status
