!----------------
module nbstatus_inc
 
use status, only: NRD

implicit none

double precision, dimension(NRD) :: NE, NHYDR, NDEUT, &
    NTRIT, NHE3, NALF, NI, NIZ1, NIZ2, NIZ3, ZIM1, ZIM2, ZIM3, &
    TE, TI, VR, SHIF,  SHIV,  ELON,  TRIA, AMETR, RHO, &
    FP, MU, AMAIN, NIBM, PIBM, PEBM, PBLON, PBPER, PBEAM, &
    SNEBM, SNNBM, CUFI, CUBM, SCUBM, SNIBM1, SNIBM2, SNIBM3, &
    NNBM1, NNBM2, NNBM3, NN, TN, ZEF, G33, IPOL

end module nbstatus_inc
