module scalars

implicit none

integer, parameter :: NEQNS=19

double precision, pointer :: &
    AB,    ABC,   AIM1,  AIM2,  AIM3,  AMJ,   AWALL, BTOR, &
    ELONG, ELONM, ENCL,  ENWM,  FECR,  FFW,   FICR,  FLH, &
    GN2E,  GN2I,  IPL,   LEXT,  NNCL,  NNWM,  QECR,  QFW, &
    QICR,  QLH,   QNBI,  RTOR,  SHIFT, TRIAN, TRICH, UEXT, &
    UPDWN, WNE,   WTE,   WTI,   ZMJ, &
    ZRD1,  ZRD2,  ZRD3,  ZRD4,  ZRD5,  ZRD6,  ZRD7,  ZRD8,  ZRD9,  ZRD10, &
    ZRD11, ZRD12, ZRD13, ZRD14, ZRD15, ZRD16, ZRD17, ZRD18, ZRD19, ZRD20, &
    ZRD21, ZRD22, ZRD23, ZRD24, ZRD25, ZRD26, ZRD27, ZRD28, ZRD29, ZRD30, &
    ZRD31, ZRD32, ZRD33, ZRD34, ZRD35, ZRD36, ZRD37, ZRD38, ZRD39, ZRD40, &
    ZRD41, ZRD42, ZRD43, ZRD44, ZRD45, ZRD46, ZRD47, ZRD48, ZRD49, ZRD50, &
    ZRD51, ZRD52, ZRD53, ZRD54, ZRD55, ZRD56, ZRD57, ZRD58, ZRD59, ZRD60, &
    ZRD61, ZRD62, ZRD63, ZRD64, ZRD65, ZRD66, ZRD67, ZRD68, ZRD69, ZRD70, &
    ZRD71, ZRD72, ZRD73, ZRD74, ZRD75, ZRD76, ZRD77, ZRD78, ZRD79, ZRD80, &
    ZRD81, ZRD82, ZRD83, ZRD84, ZRD85, ZRD86, ZRD87, ZRD88, ZRD89, ZRD90, &
    ZRD91, ZRD92, ZRD93, ZRD94, ZRD95, ZRD96

double precision, pointer :: &
    ABX,    ABCX,   AIM1X, AIM2X, AIM3X,  AMJX,   AWALLX, BTORX, &
    ELONGX, ELONMX, ENCLX, ENWMX, FECRX,  FFWX,   FICRX,  FLHX, &
    GN2EX,  GN2IX,  IPLX,  LEXTX, NNCLX,  NNWMX,  QECRX,  QFWX, &
    QICRX,  QLHX,   QNBIX, RTORX, SHIFTX, TRIANX, TRICHX, UEXTX, &
    UPDWNX, WNEX,   WTEX,  WTIX,  ZMJX, &
    ZRD1X,  ZRD2X,  ZRD3X,  ZRD4X,  ZRD5X,  ZRD6X,  ZRD7X,  ZRD8X,  ZRD9X,  ZRD10X, &
    ZRD11X, ZRD12X, ZRD13X, ZRD14X, ZRD15X, ZRD16X, ZRD17X, ZRD18X, ZRD19X, ZRD20X, &
    ZRD21X, ZRD22X, ZRD23X, ZRD24X, ZRD25X, ZRD26X, ZRD27X, ZRD28X, ZRD29X, ZRD30X, &
    ZRD31X, ZRD32X, ZRD33X, ZRD34X, ZRD35X, ZRD36X, ZRD37X, ZRD38X, ZRD39X, ZRD40X, &
    ZRD41X, ZRD42X, ZRD43X, ZRD44X, ZRD45X, ZRD46X, ZRD47X, ZRD48X, ZRD49X, ZRD50X, &
    ZRD51X, ZRD52X, ZRD53X, ZRD54X, ZRD55X, ZRD56X, ZRD57X, ZRD58X, ZRD59X, ZRD60X, &
    ZRD61X, ZRD62X, ZRD63X, ZRD64X, ZRD65X, ZRD66X, ZRD67X, ZRD68X, ZRD69X, ZRD70X, &
    ZRD71X, ZRD72X, ZRD73X, ZRD74X, ZRD75X, ZRD76X, ZRD77X, ZRD78X, ZRD79X, ZRD80X, &
    ZRD81X, ZRD82X, ZRD83X, ZRD84X, ZRD85X, ZRD86X, ZRD87X, ZRD88X, ZRD89X, ZRD90X, &
    ZRD91X, ZRD92X, ZRD93X, ZRD94X, ZRD95X, ZRD96X

double precision, pointer :: &
    CF1,   CF2,   CF3,   CF4,   CF5,   CF6,   CF7,   CF8, &
    CF9,   CF10,  CF11,  CF12,  CF13,  CF14,  CF15,  CF16, &
    CV1,   CV2,   CV3,   CV4,   CV5,   CV6,   CV7,   CV8, &
    CV9,   CV10,  CV11,  CV12,  CV13,  CV14,  CV15,  CV16, &
    CHE1,  CHE2,  CHE3,  CHE4,  CHI1,  CHI2,  CHI3,  CHI4, &
    CNB1,  CNB2,  CNB3,  CNB4,  CNBI1, CNBI2, CNBI3, CNBI4, &
    CCD1,  CCD2,  CCD3,  CCD4,  CRF1,  CRF2,  CRF3,  CRF4, &
    CNEUT1,CNEUT2,CNEUT3,CNEUT4,CPEL1, CPEL2, CPEL3, CPEL4, &
    CBND1, CBND2, CBND3, CBND4, CFUS1, CFUS2, CFUS3, CFUS4, &
    CIMP1, CIMP2, CIMP3, CIMP4, CMHD1, CMHD2, CMHD3, CMHD4, &
    CRAD1, CRAD2, CRAD3, CRAD4, CSOL1, CSOL2, CSOL3, CSOL4, &
    CSCL1, CSCL2, CSCL3, CSCL4, CSCL5, CSCL6, CSCL7, CSCL8, &
    CDWM0, CDWM1, CDWM2, CDWM3, CDWM4, CDWM5, CDWM6, CDWM7, CDWM8, CDWM9, &
    CDYM0, CDYM1, CDYM2, CDYM3, CDYM4, CDYM5, CDYM6, CDYM7, CDYM8, CDYM9, &
    CDVM0, CDVM1, CDVM2, CDVM3, CDVM4, CDVM5, CDVM6, CDVM7, CDVM8, CDVM9, &
    CDBC1, CDBC2, CDBC3, CDBC4, CDBC5, CDBC6, CDBC7, CDBC8, CDBC9, &
    CDJM1, CDJM2, CDJM3, CDJM4, CDJM5, CDJM6, CDJM7, CDJM8, CDJM9, &
    CDmj1, CDmj2, CDmj3, CDmj4, CDmj5, CDmj6, CDmj7, CDmj8, CDmj9, &
    CDhj1, CDhj2, CDhj3, CDhj4, CDhj5, CDhj6, CDhj7, CDhj8, CDhj9

double precision, target, allocatable, dimension(:) :: constValues, &
    varValues, varxValues, internValues, intern2Values
integer, target, allocatable, dimension(:) :: internIntValues

double precision, allocatable :: TEQ(:)

double precision, pointer :: &
    DPOUT, TIME, TAUMIN, TAUMAX, TAUINC, DELVAR, &
    TINIT, TSCALE, XOUT, XINPUT, ITFBE, &
    NB2EQL, DTEQL, TPAUSE, TEND, DTEQ(:, :)

integer, pointer :: &
    NEQUIL, MEQUIL, INUME1, INUME2, INUME3, INUME4, &
    IPROT, ITFBP, ICIRCQ, IPCTRL, &
    SGNIP, SGNBT, IPEQL, IFBEY, IBCPSI

double precision :: IBKDW ! IBKDW=-1 for breakdown yes

double precision, pointer :: &
    HRO, HROX, VOLUME, ROC, ROWALL, FTO, &
    FTN, BTN, PSIFB, RBDOT, ALBPL, &
    TSTART, TAU, TAU_NEW, TAU_OLD, TIMEQL, QBEAM, IPLN, &
    ROCO, RON, ROE, ROI, ROU, &
    RO0, RO1, RO2, RO3, RO4, RO5, RO6, RO7, RO8, RO9, &
    PSIAX, PSIBO, PSIFBO, PSIEXT, PSPLEX, IPLFBE, &
    PSIEXO, PSPLXO, ATREQ, PTREQ, BBDOT, TAUPRP, &
    TEB, TIB, NEB, UPARB, QEB, QIB, QNB, MUB, TTRQB, QETB, QITB, QNNB, &
    F0B, F1B, F2B, F3B, F4B, F5B, F6B, F7B, F8B, F9B, &
    QF0B, QF1B, QF2B, QF3B, QF4B, QF5B, QF6B, QF7B, QF8B, QF9B, &
    QFF0B, QFF1B, QFF2B, QFF3B, QFF4B, QFF5B, QFF6B, QFF7B, QFF8B, QFF9B

integer :: NA, NA1, NAB, NB1, NNCX, KEY, &
    NSTEPS, ITREQ, IPART, LEQ(NEQNS), NITOT

integer :: NA1N, NA1E, NA1I, NA1U, &
    NA10, NA11, NA12, NA13, NA14, NA15, NA16, NA17, NA18, NA19
character(len=132) :: exp_header

double precision :: tbeg_eq, tend_eq

contains

!------------------------------------------
    subroutine scalars_init

    constValues   = 0.
    varValues     = 0.
    varxValues    = 0.
    internValues  = 0.
    intern2Values = 0.

    NA1N = 0
    NA1E = 0
    NA1I = 0
    NA1U = 0
    NA10 = 0
    NA11 = 0
    NA12 = 0
    NA13 = 0
    NA14 = 0
    NA15 = 0
    NA16 = 0 
    NA17 = 0
    NA18 = 0
    NA19 = 0

    DPOUT  = 0.01
    TIME   = 0.
    TAUMIN = 1.e-6
    TAUMAX = 0.05
    TAUINC = 1.1
    DELVAR = 0.1
    TINIT  = 0.
    TSCALE = 1.
    XOUT   = 1.
    XINPUT = 1.
    NB2EQL = 1.
    NEQUIL = 0
    DTEQL  = 0.
    MEQUIL = 0
    TPAUSE = 100.
    TEND   = 1000.
    INUME1 = 22
    INUME2 = 22
    INUME3 = 22
    INUME4 = 22
    IPROT  = 0
    ITFBE  = 1.e6
    ITFBP  = 0
    ICIRCQ = 0
    IPCTRL = 0
    SGNIP  = 1
    SGNBT  = 1
    IPEQL  = 4 ! 4- SPIDER, 5- FEQIS
    DTEQ(1, :) = 0.
    DTEQ(2, :) = -99999.
    DTEQ(3, :) =  99999.
    DTEQ(4, :) = -1.

    TAU    = 0.000001
    TAUPRP = 0.000001

! Global variables:
    AB     = 0.3
    ABC    = 0.3
    AMJ    = 2.
    AWALL  = 0.4
    BTOR   = 3.
    ELONG  = 1.
    ELONM  = 1.
    ENCL   = 0.002
    ENWM   = 0.02
    IPL    = 0.3
    RTOR   = 1.5
    ROC    = 0.3
    ROCO   = 0.3
    SHIFT  = 0.
    UPDWN  = 0.
    NNCL   = 0.001
    NNWM   = 0.0001
    GN2E   = 0.
    GN2I   = 0.
    UEXT   = 0.
    TRIAN  = 0.0
    TRICH  = 0.0
    ZMJ    = 1.
    WNE    = 0.03
    WTE    = 0.03
    WTI    = 0.03
    TSTART = 0.
    ITREQ  = 0
    PSIFB  = 0.
    PSIFBO = 0.
    RBDOT  = 0.
    PSIEXT = 0.
    PSPLEX = 0.
    IPLFBE = 0.
    PSIEXO = 0.
    PSPLXO = 0.
    ATREQ  = 1.E-04
    PTREQ  = 1.E-05
    BBDOT  = 0.0
    IFBEY  = 0
    IBCPSI = 0
    NB1    = 41
    NA1    = 41
    NNCX   = 200
    NAB    = 41
    NA     = 40
    NITOT  = 0
    NSTEPS = 0

    TEQ = -1.e3
    LEQ = -1
    exp_header(:) = ' '

    end subroutine scalars_init

end module scalars
