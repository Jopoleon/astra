module strahl_mod

use parameter_inc, only: NRD

implicit none

double precision, dimension(NRD) :: prad_tot, ne_source, nneut_imp
double precision, dimension(NRD, 11) :: prad_strahl, nimp_strahl, nesrc_strahl

end module strahl_mod
