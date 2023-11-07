subroutine DoutTGLF

use const_inc, only: NA1, ROC, CFUS1, NNCL, CF14, CF15, CFUS4, CHI2, CIMP4, &
    CRAD4, CHI1, CF3, CSOL3, CMHD1, CMHD3, CDWM1, TIME, GP2, RTOR, HRO, GP, BTOR
use status_inc, only: RHO, AMETR, G11, VOLUM, SLAT, NE, NEX, TE, TEX, &
    NI, NDEUT, TI, TIX, VTOR, VTORX, HE, XI, PE, PI, PEBM, PIBM, &
    PEECR, PEICR, CC, ULON, IPOL, G33, ZMAIN, AMAIN, SNTOT, &
    PIICR, PRAD, SCUBM, SNEBM, SN, SNN, MU, UPL, ZEF, CUBS, PETOT, PITOT, &
    PBPER, PBLON, FP, CU, CD, CN, CE, CI, CAR17, CAR18, NN, EQPF, DN, ER, &
    CAR41, CAR42, CAR52, CAR56, CAR57, CAR58, CAR59, CAR60, CAR44, NIBM, &
    SNNBM, CAR25, CAR26, VR, QE, QI, QN

use outcmn_inc, only: equ_file, exp_file, awd, machine
  
implicit none

integer :: l, j
double precision :: COULG, POH, PEICL, QNTOT, QETOT, QITOT, BPFIT, BETPT, ZF1, ZF2, PREST, DTEQL
double precision, dimension(NA1) :: ZPOH, ZPEICL, ZXIMY, ZXINEO, ZQNTOT, ZQETOT, ZQITOT, &
    ZHEXP, ZXEXP, ZBPFIT, ZBETPT, ZPREST
character(len=80) :: filemod
character(len=6)  :: str_cdwm1, str_time

DTEQL = 0.1

do j=1, NA1
    include 'fml/coulg'
    include 'fml/poh'
    include 'fml/peicl'
    include 'fml/bpfit'
    include 'fml/betpt'
    include 'fml/prest'
    ZPOH(j)   = POH
    ZPEICL(j) = PEICL
    ZQNTOT(j) = QNTOT
    ZQETOT(j) = QETOT
    ZQITOT(j) = QITOT
    ZXEXP(j)  = CAR41(j)
    ZHEXP(j)  = CAR42(j)
    ZXINEO(j) = CAR52(j)
    ZBPFIT(j) = BPFIT
    ZBETPT(j) = BETPT
    ZPREST(j) = PREST
enddo

write(str_cdwm1, '(F6.4)') CDWM1
write(str_time, '(F6.4)') TIME

filemod = trim(awd) // 'out/CPout_' // trim(exp_file) // '_' // &
    trim(machine) // '_W_' // str_cdwm1 // '_at_' // str_time

open(31, file=trim(filemod))
write(31, *)
write(31, *) '%RHO,AMETR,DUMMYS'
write(31, '(83e25.11)') &
    (RHO(l)/ROC, AMETR(L), G11(l), VOLUM(l), &
    SLAT(l), NE(l), NEX(l), TE(l), TEX(l), NI(l), NDEUT(l), &
    TI(l), TIX(l), VTOR(l), VTORX(l), HE(l), XI(l), &
    ZXEXP(l), ZHEXP(l), ZXINEO(l), PE(l), PI(l), &
    ZPOH(l), ZPEICL(l), PEBM(l), PIBM(l), PEECR(l), &
    PEICR(l), PIICR(l), PRAD(l), SCUBM(l), SNEBM(l), &
    SN(l), SNN(l), ZQETOT(l), ZQITOT(l), ZQNTOT(l), MU(l), &
    UPL(l), ZEF(l), CUBS(l), ZBPFIT(l), FP(l), ZBETPT(l), CU(l), CD(l), &
    ZPREST(l), CN(l), CE(l), CI(l), CAR17(l), CAR18(l), CFUS1, NNCL, &
    CF14, CF15, CFUS4, CHI2, CIMP4, CRAD4, NN(l), EQPF(l), DN(l), CHI1, &
    CF3, ER(l), CAR56(l), CAR57(l), CAR58(l), CAR59(l), CAR60(l), &
    CAR44(l), CSOL3, NIBM(l), SNNBM(l), CMHD1, CMHD3, &
    CAR25(l), CAR26(l), VR(l), QE(l), QI(l), QN(l), l=1, NA1)
close(31)

print *, 'Output stored in ', filemod

return
end subroutine DoutTGLF
