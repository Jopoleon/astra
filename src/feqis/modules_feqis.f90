module dimensions_ef_parameters

implicit none

integer, parameter :: i_dim1=300, i_dim2=300 ! #coils, #rad.grid

end module dimensions_ef_parameters

!-------------------------------------------
module ef_circuit       ! declaration of minimal CPOs

use dimensions_ef_parameters, only: i_dim1, i_dim2

implicit none 

integer :: nr2, nz2, nrho, nteta

double precision :: psibndp, psiaxisp, raxp, zaxp, &
    iplasma, btor0, rgeom0, psplex, li3, betapol
double precision, dimension(i_dim2) :: teta, rbndp, zbndp, rexp, zexp, &
    psia_2d, ffp_2d, ppp_2d, ipol_2d, pres_2d, pprime, ffprime, pressure, psigrida, ipol

end module ef_circuit
