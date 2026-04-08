module nbstatus

use status, only: NRD

implicit none

integer, parameter :: n_rho=NRD, n_nbi_max=16, n_fields=20, n_theta=51, n_energy=3
!   n_nbi_max - Max. No. of NBI sources (max No. of groupss in *.nbi file)
!   n_fields - No. of fields in one group of *.nbi file
!   n_theta - #pitch angle
!   n_energy: full, half, 1/3

integer :: ISPEND, ISPE(9)
double precision, dimension(n_rho) :: NE, NHYDR, NDEUT, &
    NTRIT, NHE3, NALF, NI, NIZ1, NIZ2, NIZ3, ZIM1, ZIM2, ZIM3, &
    TE, TI, VR, SHIF,  SHIV,  ELON,  TRIA, AMETR, RHO, &
    FP, MU, AMAIN, NIBM, PIBM, PEBM, PBLON, PBPER, PBEAM, &
    SNEBM, SNNBM, CUFI, CUBM, SCUBM, SNIBM1, SNIBM2, SNIBM3, &
    NNBM1, NNBM2, NNBM3, NN, TN, ZEF, G33, IPOL
double precision :: YRIPLR(n_rho)

contains

!---------------------------------------------------------------------
    subroutine set_input(yNE, yNHYDR, yNDEUT, yNTRIT, &
        yNHE3, yNALF, yNI, yNIZ1, yNIZ2, yNIZ3, yZIM1, yZIM2, yZIM3, yTE, yTI, &
        yVR, ySHIF, ySHIV, yELON, yTRIA, yAMETR, yRHO, yFP, yMU, yAMAIN, &
        yNN, yTN, yZEF, yG33, yIPOL, yNIBM, yPIBM, yPEBM, yPBLON, yPBPER, &
        yPBEAM, ySNEBM, ySNNBM, yCUFI, yCUBM, ySCUBM, &
        ySNIBM1, ySNIBM2, ySNIBM3, yNNBM1, yNNBM2, yNNBM3)
! Called from sbr/nbi.f90

    double precision, intent(in), dimension(*) :: YNE, YNHYDR, YNDEUT, YNTRIT, &
        YNHE3, YNALF, YNI, YNIZ1, YNIZ2, YNIZ3, YZIM1, YZIM2, YZIM3, YTE, YTI, &
        YVR, YSHIF, YSHIV, YELON, YTRIA, YAMETR, YRHO, YFP, YMU, YAMAIN, &
        YG33, YIPOL, YNN, YTN, YZEF, YNIBM, YPIBM, YPEBM, &
        YPBLON, YPBPER, YPBEAM, YSNEBM, YSNNBM, YCUFI, YCUBM, YSCUBM, &
        YSNIBM1, YSNIBM2, YSNIBM3, YNNBM1, YNNBM2, YNNBM3

    integer :: j

    do j=1, n_rho
        NE(j)    = YNE(j)
        NHYDR(j) = YNHYDR(j)
        NDEUT(j) = YNDEUT(j)
        NTRIT(j) = YNTRIT(j)
        NHE3(j)  = YNHE3(j)
        NALF(j)  = YNALF(j)
        NI(j)    = YNI(j)
        NIZ1(j)  = YNIZ1(j)
        NIZ2(j)  = YNIZ2(j)
        NIZ3(j)  = YNIZ3(j)
        ZIM1(j)  = YZIM1(j)
        ZIM2(j)  = YZIM2(j)
        ZIM3(j)  = YZIM3(j)
        TE(j)    = YTE(j)
        TI(j)    = YTI(j)
        VR(j)    = YVR(j)
        SHIF(j)  = YSHIF(j)
        SHIV(j)  = YSHIV(j)
        ELON(j)  = YELON(j)
        TRIA(j)  = YTRIA(j)
        AMETR(j) = YAMETR(j)
        RHO(j)   = YRHO(j)
        FP(j)    = YFP(j)
        MU(j)    = YMU(j)
        AMAIN(j) = YAMAIN(j)
        NIBM(j)  = YNIBM(j)
        PIBM(j)  = YPIBM(j)
        PEBM(j)  = YPEBM(j)
        PBLON(j) = YPBLON(j)
        PBPER(j) = YPBPER(j)
        PBEAM(j) = YPBEAM(j)
        SNEBM(j) = YSNEBM(j)
        SNNBM(j) = YSNNBM(j)
        CUFI(j)  = YCUFI(j)
        CUBM(j)  = YCUBM(j)
        SCUBM(j) = YSCUBM(j)
        SNIBM1(j)= YSNIBM1(j)
        SNIBM2(j)= YSNIBM2(j)
        SNIBM3(j)= YSNIBM3(j)
        NNBM1(j) = YNNBM1(j)
        NNBM2(j) = YNNBM2(j)
        NNBM3(j) = YNNBM3(j)
        NN(j)    = YNN(j)
        TN(j)    = YTN(j)
        ZEF(j)   = YZEF(j)
        G33(j)   = YG33(j)
        IPOL(j)  = YIPOL(j)
    enddo

    end subroutine set_input

!---------------------------------------------------------------------
    subroutine get_output(yNIBM, yPIBM, yPEBM, yPBLON, yPBPER, &
        yPBEAM, ySNEBM, ySNNBM, yCUFI, yCUBM, ySCUBM, &
        ySNIBM1, ySNIBM2, ySNIBM3, yNNBM1, yNNBM2, yNNBM3)
! Called from sbr/nbi.f90

    double precision, intent(out), dimension(*) :: YNIBM, YPIBM, YPEBM, &
        YPBLON, YPBPER, YPBEAM, YSNEBM, YSNNBM, YCUFI, YCUBM, YSCUBM, &
        YSNIBM1, YSNIBM2, YSNIBM3, YNNBM1, YNNBM2, YNNBM3

    integer :: j

    do j=1, n_rho
        YNIBM(j)   = NIBM(j)
        YPIBM(j)   = PIBM(j)
        YPEBM(j)   = PEBM(j)
        YPBLON(j)  = PBLON(j)
        YPBPER(j)  = PBPER(j)
        YPBEAM(j)  = PBEAM(j)
        YSNEBM(j)  = SNEBM(j)
        YSNNBM(j)  = SNNBM(j)
        YCUFI(j)   = CUFI(j)
        YCUBM(j)   = CUBM(j)
        YSCUBM(j)  = SCUBM(j)
        YSNIBM1(j) = SNIBM1(j)
        YSNIBM2(j) = SNIBM2(j)
        YSNIBM3(j) = SNIBM3(j)
        YNNBM1(j)  = NNBM1(j)
        YNNBM2(j)  = NNBM2(j)
        YNNBM3(j)  = NNBM3(j)
    enddo

    end subroutine get_output

end module nbstatus
