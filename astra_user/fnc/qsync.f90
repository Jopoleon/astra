! QSYNC [MW]:  Integral {0, R} ( PSYNC ) dV
!   (Pereverzev 9-AUG-02)
double precision function QSYNCR(YR)

use const_inc, only: ROC, BTOR, RTOR, AB
use status_inc, only: VOLUM

implicit none

double precision YR, PSYNC, RADIAL, TEAVR, NEAVR

PSYNC = 1.32E-7*(TEAVR(ROC)*BTOR)**2.5  *SQRT(NEAVR(ROC)/AB*(1.+18.*AB/(RTOR*SQRT(TEAVR(ROC)))))
QSYNCR = PSYNC*RADIAL(VOLUM, YR)

end function QSYNCR
