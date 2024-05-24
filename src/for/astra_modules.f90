module parameter_inc
 
implicit none
save

integer, parameter :: NRD=801, NRW=128, NCONST=256, NARRX=101, NSBMX=60, &
    NSDELOUT=39, NCNBM=60, NCNBTM=25000, NRDX=500, NTVAR=250000, &
    NTARR=250000, NEQNS=19, NBDMAX=256, NBDTMAX=1500, plot_modes=9

end module parameter_inc

!--------------------------------
module timeoutput_inc

use parameter_inc, only: NRW

implicit none

integer, parameter :: NTIMES=1024

! TOUT   - Time variables output array
! TTOUT  - time-coordinate array for time output [s] TTOUT(1:LTOUT<=NTIMES)
double precision :: TTOUT(NTIMES), TOUT(NTIMES, NRW), TPOUT

end module timeoutput_inc

!--------------------------------
module expdat

use parameter_inc, only: NTVAR, NTARR, NRDX, NBDMAX, NBDTMAX

implicit none

type rawScalar
    integer, dimension(NTVAR) :: var_index=0
    double precision, dimension(NTVAR) :: time=0., value=0., error=0.
    character(len=6), dimension(NTVAR) :: label
endtype rawScalar
type rawProfileMap
    integer, dimension(NTARR) :: arr_index=0, jbeg_grid=0, jbeg_data=0, grid_type=0, nrho=0
    double precision, dimension(NTARR) :: time=0., filter=0.001
    character(len=6), dimension(NTARR) :: label
endtype rawProfileMap

type(rawScalar) :: raw_scalar
type(rawProfileMap) :: raw_profile_map

real*4 :: DATARR(NRDX*NTARR)
double precision, dimension(NBDTMAX) :: BNDTIM
double precision, dimension(NBDTMAX*NBDMAX) :: BNDR, BNDZ

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
