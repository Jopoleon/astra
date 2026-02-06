module parameter_inc
 
implicit none
save

integer, parameter :: NRD=801, NCONST=256, NARRX=101, n_sbr_max=60, &
    n_coils_max=60, nt_coils_max=25000, NRDX=500, NTVAR=250000, &
    NTARR=250000, NEQNS=19, n_bnd_max=256, nt_bnd_max=1500, plot_modes=9

end module parameter_inc

!--------------------------------
module expdat

use parameter_inc, only: NTVAR, NTARR, NRDX, n_bnd_max, nt_bnd_max

implicit none

type rawScalars
    integer, dimension(NTVAR) :: var_index=0
    double precision, dimension(NTVAR) :: time=0., data=0., error=0.
    character(len=6), dimension(NTVAR) :: label
endtype rawScalars
type rawProfiles
    integer, dimension(NTARR) :: arr_index=0, jbeg_grid=0, jbeg_data=0, grid_type=0, nrho=0
    double precision, dimension(NTARR) :: time=0., filter=0.001
    character(len=6), dimension(NTARR) :: label
    real*4, dimension(NRDX*NTARR) :: data
endtype rawProfiles
type rawBoundary
    double precision, dimension(nt_bnd_max) :: time=0.
    double precision, dimension(nt_bnd_max*n_bnd_max) :: R, Z
endtype rawBoundary

type(rawScalars)  :: raw_scalars
type(rawProfiles) :: raw_profiles
type(rawBoundary) :: raw_boundary

end module expdat

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
