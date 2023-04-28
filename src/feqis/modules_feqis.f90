module pi_vars

implicit none

double precision, parameter :: GPI=3.141592653589793, GPI2=2.*GPI, &
    GPI4=GPI2**2, muvac=4.e-7*GPI

end module pi_vars

!---------------------------------------------------------------------
module feqis_geom

implicit none 

double precision :: raxp, zaxp
double precision, allocatable, dimension(:) :: theta

end module feqis_geom
