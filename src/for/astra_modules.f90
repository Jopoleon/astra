module parameter_inc
 
implicit none
save

integer, parameter :: NRD=801, NCONST=256, NARRX=101, n_sbr_max=60, &
    n_coils_max=60, nt_coils_max=25000, NRDX=500, NTVAR=250000, &
    NTARR=250000, NEQNS=19, n_bnd_max=256, nt_bnd_max=1500, plot_modes=9

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
