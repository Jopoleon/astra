!ITER 98(y,2): TAU_scal*((Psep+Prad)**0.69) --> tau_scal only uses Pabs, not Psep
double precision function H98ITERR(YR)

use const_inc, only: BTOR, IPL, ABC, NA1, ELONG, AMJ, RTOR
use status_inc, only: NE

implicit none

double precision, intent(in) :: YR

H98ITR = 0.0562 * IPL**0.93 * BTOR**0.15 * &
    sum(NE(1:NA1)/DBLE(NA1))**0.41 * AMJ**0.19 * RTOR**1.97 * &
    (ABC/RTOR)**0.58 * ELONG**0.78

return
end function H98ITERR
