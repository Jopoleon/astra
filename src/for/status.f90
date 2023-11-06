module status_inc
 
use parameter_inc, only: NRD, NARRX

implicit none

double precision, dimension(NRD) :: &
    TEO, TIO, NEO, UPAR, UPARO, &
    FPO, ULON, QU, VP, GN, SQEPS, &
    PRES, PELH, PETOT, PELON, CULH, PITOT, &
    SNTOT, PEECR, PEFW, PEICR, PEIQI, PIICR, PIFW, PEPER, &
    QE, QI, QN, QF0, QF1, QF2, QF3, QF4, QF5, QF6, QF7, &
    QF8, QF9, CUECR, CUFW, CUICR, CUTOR, TTRQ, VR, VRO, &
    NIO

double precision, dimension(NRD) :: &
    SHIV, CUBS, &
    VPFP, G22, G33, G41, G42, G43, G44, G45, &
    VOLUM, CV, DRODA, PDE, PDI, SDN, &
    SD0, SD1, SD2, SD3, SD4, SD5, SD6, SD7, SD8, SD9, &
    UPS0, UPS0O, UPS1, UPS1O, DLNEO, SGNEO, UPS2, &
    UPS2O, MRHO, DDNEO, SGNEOD, DLNEOD, SQUARN  !47

double precision, dimension(NRD) :: &
    B0DB2, BDB02, BDB0, BMAXT, BMINT, FOFB, GRADRO, &
    EQFF, EQPF, SLAT, FV, MV, XRHO, &
    SXHO, SRHO, G22E, G33E, AREAT, PERIM

double precision, dimension(:), pointer :: &
    AMAIN, AMETR, CU, ELON, ER, FP, G11, IPOL, MU, &
    NALF, NDEUT, NE, NHE3, NHYDR, NI, NIBM, NIZ1, NIZ2, NIZ3, NTRIT, &
    PBLON, PBPER, PFAST, RHO, SHEAR, SHIF, TE, TI, TRIA, UPL, VPOL, VRS, VTOR, &
    ZEF, ZIM1, ZIM2, ZIM3, ZMAIN, ZIMPT, NIMPT, AIMPT
!zimpt, nimpt, aimpt --> average charge, total density, average mass of impurities

double precision, dimension(NRD*41), target :: plasma_profs

double precision, dimension(NRD) :: &
    NN, TN, PRAD, PBOL1, PBOL2, PBOL3, PSXR1, &
    PSXR2, PSXR3, ZEF1, ZEF2, ZEF3, &
    DIMP1, DIMP2, DIMP3, VIMP1, VIMP2, &
    VIMP3, NMAIN

double precision, dimension(:), pointer :: &
    CAR1,  CAR2,  CAR3,  CAR4,  CAR5 , CAR6,  CAR7,  CAR8,  CAR9,  CAR10, &
    CAR11, CAR12, CAR13, CAR14, CAR15, CAR16, CAR17, CAR18, CAR19, CAR20, &
    CAR21, CAR22, CAR23, CAR24, CAR25, CAR26, CAR27, CAR28, CAR29, CAR30, &
    CAR31, CAR32, CAR33, CAR34, CAR35, CAR36, CAR37, CAR38, CAR39, CAR40, &
    CAR41, CAR42, CAR43, CAR44, CAR45, CAR46, CAR47, CAR48, CAR49, CAR50, &
    CAR51, CAR52, CAR53, CAR54, CAR55, CAR56, CAR57, CAR58, CAR59, CAR60, &
    CAR61, CAR62, CAR63, CAR64 
double precision, dimension(NRD, 64), target :: CAR

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

double precision, dimension(NRD, NARRX), target :: EXT

double precision, dimension(NRD) :: &
    PBEAM, SNNBM, SNEBM, NNBM1, NNBM2, NNBM3, &
    SCUBM, CUBM, PEBM, PIBM, &
    CUFI, PBFUS, SNIBM1, SNIBM2, SNIBM3

double precision :: NNBM(NRD, 3), IOBMN(NRD, 17)

double precision, dimension(NRD) :: &
    CC,  CD,  CE,  CI,  CN, &
    DC,  DE,  DI,  DN,  HC, &
    HE,  HI,  HN,  PE,  PET, &
    PI,  PIT, SN,  SNN, XC, &
    XE,  XI,  XN,  DVN, DSN, &
    DVE, DSE, DVI, DSI, CNPAR, &
    CNPAP, XUPAR, RUPAR, RUPYR, RUPFR, &
    CNPAD, XUPAP, TTRQI, XUPAD, DDNEOD

double precision :: WORK(NRD, 100)

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

double precision, dimension(NRD, 0:9), target :: &
    FJ, DFJ, VFJ, GFJ, SFJ, DVFJ, DSFJ, SFFJ, FJO, SFJTOT

double precision, dimension(:), pointer :: &
    VOLDR, VNEWR, CUBTS, R2MTR, G2GR2, G22G, SHIH

contains

!---------------------
subroutine status_init

integer :: j

! Defaults, rather fall-back than initial values

AMAIN => plasma_profs(         1:    NRD)
AMETR => plasma_profs(   NRD + 1:  2*NRD)
CU    => plasma_profs( 2*NRD + 1:  3*NRD)
ELON  => plasma_profs( 3*NRD + 1:  4*NRD)
ER    => plasma_profs( 4*NRD + 1:  5*NRD)
FP    => plasma_profs( 5*NRD + 1:  6*NRD)
G11   => plasma_profs( 6*NRD + 1:  7*NRD)
IPOL  => plasma_profs( 7*NRD + 1:  8*NRD)
MU    => plasma_profs( 8*NRD + 1:  9*NRD)
NALF  => plasma_profs( 9*NRD + 1: 10*NRD)
NDEUT => plasma_profs(10*NRD + 1: 11*NRD)
NE    => plasma_profs(11*NRD + 1: 12*NRD)
NHE3  => plasma_profs(12*NRD + 1: 13*NRD)
NHYDR => plasma_profs(13*NRD + 1: 14*NRD)
NI    => plasma_profs(14*NRD + 1: 15*NRD)
NIBM  => plasma_profs(15*NRD + 1: 16*NRD)
NIZ1  => plasma_profs(16*NRD + 1: 17*NRD)
NIZ2  => plasma_profs(17*NRD + 1: 18*NRD)
NIZ3  => plasma_profs(18*NRD + 1: 19*NRD)
NTRIT => plasma_profs(19*NRD + 1: 20*NRD)
PBLON => plasma_profs(20*NRD + 1: 21*NRD)
PBPER => plasma_profs(21*NRD + 1: 22*NRD)
PFAST => plasma_profs(22*NRD + 1: 23*NRD)
RHO   => plasma_profs(23*NRD + 1: 24*NRD)
SHEAR => plasma_profs(24*NRD + 1: 25*NRD)
SHIF  => plasma_profs(25*NRD + 1: 26*NRD)
TE    => plasma_profs(26*NRD + 1: 27*NRD)
TI    => plasma_profs(27*NRD + 1: 28*NRD)
TRIA  => plasma_profs(28*NRD + 1: 29*NRD)
UPL   => plasma_profs(29*NRD + 1: 30*NRD)
VPOL  => plasma_profs(30*NRD + 1: 31*NRD)
VRS   => plasma_profs(31*NRD + 1: 32*NRD)
VTOR  => plasma_profs(32*NRD + 1: 33*NRD)
ZEF   => plasma_profs(33*NRD + 1: 34*NRD)
ZIM1  => plasma_profs(34*NRD + 1: 35*NRD)
ZIM2  => plasma_profs(35*NRD + 1: 36*NRD)
ZIM3  => plasma_profs(36*NRD + 1: 37*NRD)
ZMAIN => plasma_profs(37*NRD + 1: 38*NRD)
ZIMPT => plasma_profs(38*NRD + 1: 39*NRD)
NIMPT => plasma_profs(39*NRD + 1: 40*NRD)
AIMPT => plasma_profs(40*NRD + 1: 41*NRD)

! Geometry/equilibrium

ELON = 1.
TRIA = 0.
SHIF = 0.
SHIV = 0.
SQUARN = 0.
G33 = 1.
G41 = 1.
G42 = 1.
G43 = .2
G44 = 1.
G45 = 1.
G22E = 1.
G33E = 1.
SQEPS = 0.2 
AREAT = 1.
PERIM = 1.
VOLUM = 1.
ER = 0.

UPS0 = 0.
UPS1 = 0.
UPS2 = 0.
UPS0O = 0.
UPS1O = 0.
UPS2O = 0.

! Species

VTOR = 0.
VPOL = 0.
UPAR = 0.
ZEF = 1.

AMAIN = 1.

ZIMPT = 1.
NIMPT = 0.
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
NHYDR = 0.
NDEUT = 0.
NTRIT = 0.
NHE3  = 0.
NALF  = 0.
NIZ1  = 0.
NIZ2  = 0.
NIZ3  = 0.
NMAIN = 0.
PEIQI = 0.

! Current

CU = 0.1
MU = 0.3
IPOL = 1.

! Transport

PBLON = 0.
PBPER = 0.
MRHO = 1.
DLNEO = 0.
SGNEO = 0.
DDNEO = 0.
DLNEOD = 0.
SGNEOD = 0.
DDNEOD = 0.
CNPAP = 0.
CNPAD = 0.
XUPAP = 0.
XUPAD = 0.
XUPAR = 0.
TTRQ  = 0.
TTRQI = 0.
MV = 0.
CV = 0.
FV = 0.
SNN = 0.

DVN  = 0.d0
DVE  = 0.d0
DVI  = 0.d0

DSN  = 0.d0
DSE  = 0.d0
DSI  = 0.d0

SDN = 0.d0
PDE = 0.d0
PDI = 0.d0

SD0 = 0.d0
SD1 = 0.d0
SD2 = 0.d0
SD3 = 0.d0
SD4 = 0.d0
SD5 = 0.d0
SD6 = 0.d0
SD7 = 0.d0
SD8 = 0.d0
SD9 = 0.d0

! Sources

PET = 0.
PIT = 0.
PETOT = 1.d-4
PITOT = 1.d-4

! Exp input arrays

CAR1  => CAR(:, 1)
CAR2  => CAR(:, 2)
CAR3  => CAR(:, 3)
CAR4  => CAR(:, 4)
CAR5  => CAR(:, 5)
CAR6  => CAR(:, 6)
CAR7  => CAR(:, 7)
CAR8  => CAR(:, 8)
CAR9  => CAR(:, 9)
CAR10 => CAR(:, 10)
CAR11 => CAR(:, 11)
CAR12 => CAR(:, 12)
CAR13 => CAR(:, 13)
CAR14 => CAR(:, 14)
CAR15 => CAR(:, 15)
CAR16 => CAR(:, 16)
CAR17 => CAR(:, 17)
CAR18 => CAR(:, 18)
CAR19 => CAR(:, 19)
CAR20 => CAR(:, 20)
CAR21 => CAR(:, 21)
CAR22 => CAR(:, 22)
CAR23 => CAR(:, 23)
CAR24 => CAR(:, 24)
CAR25 => CAR(:, 25)
CAR26 => CAR(:, 26)
CAR27 => CAR(:, 27)
CAR28 => CAR(:, 28)
CAR29 => CAR(:, 29)
CAR30 => CAR(:, 30)
CAR31 => CAR(:, 31)
CAR32 => CAR(:, 32)
CAR33 => CAR(:, 33)
CAR34 => CAR(:, 34)
CAR35 => CAR(:, 35)
CAR36 => CAR(:, 36)
CAR37 => CAR(:, 37)
CAR38 => CAR(:, 38)
CAR39 => CAR(:, 39)
CAR40 => CAR(:, 40)
CAR41 => CAR(:, 41)
CAR42 => CAR(:, 42)
CAR43 => CAR(:, 43)
CAR44 => CAR(:, 44)
CAR45 => CAR(:, 45)
CAR46 => CAR(:, 46)
CAR47 => CAR(:, 47)
CAR48 => CAR(:, 48)
CAR49 => CAR(:, 49)
CAR50 => CAR(:, 50)
CAR51 => CAR(:, 51)
CAR52 => CAR(:, 52)
CAR53 => CAR(:, 53)
CAR54 => CAR(:, 54)
CAR55 => CAR(:, 55)
CAR56 => CAR(:, 56)
CAR57 => CAR(:, 57)
CAR58 => CAR(:, 58)
CAR59 => CAR(:, 59)
CAR60 => CAR(:, 60)
CAR61 => CAR(:, 61)
CAR62 => CAR(:, 62)
CAR63 => CAR(:, 63)
CAR64 => CAR(:, 64)

! Exp input arrays

CAR1X  => EXT(:, 1)
CAR2X  => EXT(:, 2)
CAR3X  => EXT(:, 3)
CAR4X  => EXT(:, 4)
CAR5X  => EXT(:, 5)
CAR6X  => EXT(:, 6)
CAR7X  => EXT(:, 7)
CAR8X  => EXT(:, 8)
CAR9X  => EXT(:, 9)
CAR10X => EXT(:, 10)
CAR11X => EXT(:, 11)
CAR12X => EXT(:, 12)
CAR13X => EXT(:, 13)
CAR14X => EXT(:, 14)
CAR15X => EXT(:, 15)
CAR16X => EXT(:, 16)
CAR17X => EXT(:, 17)
CAR18X => EXT(:, 18)
CAR19X => EXT(:, 19)
CAR20X => EXT(:, 20)
CAR21X => EXT(:, 21)
CAR22X => EXT(:, 22)
CAR23X => EXT(:, 23)
CAR24X => EXT(:, 24)
CAR25X => EXT(:, 25)
CAR26X => EXT(:, 26)
CAR27X => EXT(:, 27)
CAR28X => EXT(:, 28)
CAR29X => EXT(:, 29)
CAR30X => EXT(:, 30)
CAR31X => EXT(:, 31)
CAR32X => EXT(:, 32)
CAR33X => EXT(:, 33)
CAR34X => EXT(:, 34)
CAR35X => EXT(:, 35)
CAR36X => EXT(:, 36)
CAR37X => EXT(:, 37)
CAR38X => EXT(:, 38)
CAR39X => EXT(:, 39)
CAR40X => EXT(:, 40)
CAR41X => EXT(:, 41)
CAR42X => EXT(:, 42)
CAR43X => EXT(:, 43)
CAR44X => EXT(:, 44)
CAR45X => EXT(:, 45)
CAR46X => EXT(:, 46)
CAR47X => EXT(:, 47)
CAR48X => EXT(:, 48)
CAR49X => EXT(:, 49)
CAR50X => EXT(:, 50)
CAR51X => EXT(:, 51)
CAR52X => EXT(:, 52)
CAR53X => EXT(:, 53)
CAR54X => EXT(:, 54)
CAR55X => EXT(:, 55)
CAR56X => EXT(:, 56)
CAR57X => EXT(:, 57)
CAR58X => EXT(:, 58)
CAR59X => EXT(:, 59)
CAR60X => EXT(:, 60)
CAR61X => EXT(:, 61)
CAR62X => EXT(:, 62)
CAR63X => EXT(:, 63)
CAR64X => EXT(:, 64)

F0X => EXT(:, 65)
F1X => EXT(:, 66)
F2X => EXT(:, 67)
F3X => EXT(:, 68)
F4X => EXT(:, 69)
F5X => EXT(:, 70)
F6X => EXT(:, 71)
F7X => EXT(:, 72)
F8X => EXT(:, 73)
F9X => EXT(:, 74)

MUX    => EXT(:, 75)
MVX    => EXT(:, 76)
GNX    => EXT(:, 77)
SNX    => EXT(:, 78)
PEX    => EXT(:, 79)
PIX    => EXT(:, 80)
PRADX  => EXT(:, 81)
TEX    => EXT(:, 82)
TIX    => EXT(:, 83)
NEX    => EXT(:, 84)
CUX    => EXT(:, 85)
ZEFX   => EXT(:, 86)
VRX    => EXT(:, 87)
SHX    => EXT(:, 88)
ELX    => EXT(:, 89)
TRX    => EXT(:, 90)
G11X   => EXT(:, 91)
G22X   => EXT(:, 92)
G33X   => EXT(:, 93)
DRODAX => EXT(:, 94)
IPOLX  => EXT(:, 95)
NIX    => EXT(:, 96)
VPOLX  => EXT(:, 97)
VTORX  => EXT(:, 98)
SLATX  => EXT(:, 99)
SHIVX  => EXT(:, 100)
SQUAX  => EXT(:, 101)

! F0-F9

F0   => FJ(:, 0)
F1   => FJ(:, 1)
F2   => FJ(:, 2)
F3   => FJ(:, 3)
F4   => FJ(:, 4)
F5   => FJ(:, 5)
F6   => FJ(:, 6)
F7   => FJ(:, 7)
F8   => FJ(:, 8)
F9   => FJ(:, 9)

F0O  => FJO(:, 0)
F1O  => FJO(:, 1)
F2O  => FJO(:, 2)
F3O  => FJO(:, 3)
F4O  => FJO(:, 4)
F5O  => FJO(:, 5)
F6O  => FJO(:, 6)
F7O  => FJO(:, 7)
F8O  => FJO(:, 8)
F9O  => FJO(:, 9)

DF0  => DFJ(:, 0)
DF1  => DFJ(:, 1)
DF2  => DFJ(:, 2)
DF3  => DFJ(:, 3)
DF4  => DFJ(:, 4)
DF5  => DFJ(:, 5)
DF6  => DFJ(:, 6)
DF7  => DFJ(:, 7)
DF8  => DFJ(:, 8)
DF9  => DFJ(:, 9)

SF0  => SFJ(:, 0)
SF1  => SFJ(:, 1)
SF2  => SFJ(:, 2)
SF3  => SFJ(:, 3)
SF4  => SFJ(:, 4)
SF5  => SFJ(:, 5)
SF6  => SFJ(:, 6)
SF7  => SFJ(:, 7)
SF8  => SFJ(:, 8)
SF9  => SFJ(:, 9)

VF0  => VFJ(:, 0)
VF1  => VFJ(:, 1)
VF2  => VFJ(:, 2)
VF3  => VFJ(:, 3)
VF4  => VFJ(:, 4)
VF5  => VFJ(:, 5)
VF6  => VFJ(:, 6)
VF7  => VFJ(:, 7)
VF8  => VFJ(:, 8)
VF9  => VFJ(:, 9)

GF0  => GFJ(:, 0)
GF1  => GFJ(:, 1)
GF2  => GFJ(:, 2)
GF3  => GFJ(:, 3)
GF4  => GFJ(:, 4)
GF5  => GFJ(:, 5)
GF6  => GFJ(:, 6)
GF7  => GFJ(:, 7)
GF8  => GFJ(:, 8)
GF9  => GFJ(:, 9)

DVF0 => DVFJ(:, 0)
DVF1 => DVFJ(:, 1)
DVF2 => DVFJ(:, 2)
DVF3 => DVFJ(:, 3)
DVF4 => DVFJ(:, 4)
DVF5 => DVFJ(:, 5)
DVF6 => DVFJ(:, 6)
DVF7 => DVFJ(:, 7)
DVF8 => DVFJ(:, 8)
DVF9 => DVFJ(:, 9)

DSF0 => DSFJ(:, 0)
DSF1 => DSFJ(:, 1)
DSF2 => DSFJ(:, 2)
DSF3 => DSFJ(:, 3)
DSF4 => DSFJ(:, 4)
DSF5 => DSFJ(:, 5)
DSF6 => DSFJ(:, 6)
DSF7 => DSFJ(:, 7)
DSF8 => DSFJ(:, 8)
DSF9 => DSFJ(:, 9)

SFF0 => SFFJ(:, 0)
SFF1 => SFFJ(:, 1)
SFF2 => SFFJ(:, 2)
SFF3 => SFFJ(:, 3)
SFF4 => SFFJ(:, 4)
SFF5 => SFFJ(:, 5)
SFF6 => SFFJ(:, 6)
SFF7 => SFFJ(:, 7)
SFF8 => SFFJ(:, 8)
SFF9 => SFFJ(:, 9)

SF0TOT => SFJTOT(:, 0)
SF1TOT => SFJTOT(:, 1)
SF2TOT => SFJTOT(:, 2)
SF3TOT => SFJTOT(:, 3)
SF4TOT => SFJTOT(:, 4)
SF5TOT => SFJTOT(:, 5)
SF6TOT => SFJTOT(:, 6)
SF7TOT => SFJTOT(:, 7)
SF8TOT => SFJTOT(:, 8)
SF9TOT => SFJTOT(:, 9)

! Multi-dimensional

WORK  = 0.d0
!CAR   = 0.d0
EXT   = 0.d0

FJ = 1.d0
DVFJ = 0.d0
DSFJ = 0.d0
SFFJ = 0.d0

do j=1, NRD
   FP(j) = 0.1*j
   QE(j) = FP(j)
   QU(j) = FP(j)
   QI(j) = FP(j)
   QN(j) = FP(j)
   GN(j) = FP(j)
enddo

end subroutine status_init

end module status_inc
