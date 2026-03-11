module parameter_inc
 
implicit none

integer, parameter :: NRD=801, NARRX=101, n_sbr_max=60, NRDX=500, n_coils_max=60

end module parameter_inc

!--------------------------------
module plasma_state ! for plasma yes/no (no will not solve the transport equations)

implicit none

integer :: plasma_up, plasma_trig

end module plasma_state

!--------------------------------
module ext_bnd

implicit none

double precision, dimension(:, :), allocatable :: ext_bnd_in ! 50 , 2 boundary values R,Z
integer :: use_ext_bnd

end module ext_bnd
