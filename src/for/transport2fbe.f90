module transport2fbe  !these are coupling variables with the equilibrium solver and astra

use const_inc, only: TAU, TSTART, RTOR, UPDWN, SHIFT, NA1, NCNB, PSIAX, PSIBO
use outcmn_inc, only: MACHINE, CCOILX, NCNBT
use debugger, only: debug

implicit none

integer, parameter :: ncoil_dim=300

integer :: use_limiter_astra  ! 1-uses limiter, 0-ignore limiter
integer :: refit_mode         ! if -1 - 1 pass only , 0 - self-consistent solution, if 1 - stab axis using passive wall currents fourier modes cos and sin, if 2 - same as 1 but uses boundary points using 5 fourier modes, 3-uses full currents fit using analytic F function and fit file efonfit.dat
integer :: solve_fix          ! if 0 - solve full fix boundary problem, if 1 - 1 iteration only , 2 - only contouring
integer :: execute_plasma     ! if 0 - only circuit equations, if 1 - solve plasma gseq too 
integer :: nonegcurr          ! nonegcurr = 0 --> no negative current allowed in plasma
integer :: fast_mode          ! 0 - normale, 1 - domnt do iterationsin fbe gse
integer :: reconnect_circuits ! 0-nothing, 1-recompute matrix with new circuits
integer :: new_equivalence(ncoil_dim, 10)   ! if reconnect, says what is the new_equivalence, for example (1,1,1,0,0,0,0,..) means coil 1,2,3 become 1,1,1
integer :: resistance_change  ! if  = 1, changes resistances of coils using new_resistance
integer :: n_equivalence      ! number of equivalences
integer :: use_reduce_circuit ! 0-all coils solved. 1 - some coils not solved
integer :: n_of_newton_iterations    ! to find actual mag axis. recommended between 5 - 10 
integer :: n_fourier_restab_boundary ! nr of fourier modes for boundary restab, default = 5
integer :: psplex_from_fbe    ! put 1 to get psplex fromfree boundary
integer :: plasma_config      ! 0 if limiter, 1 if xpoint
integer :: simple_plasma_model_breakdown ! if 0, no plasma feedback to coils in vacuum, if 1 yes

integer :: use_isoflux, n_isoflux         ! 0 does nothing, 1 when mode 818 is used to reconstruct coil currents, the boundary is obtained on the isoflux points r_isoflux and z_isoflux , of length n_isoflux

double precision :: tau_circuit_feqis, tau_gseq_feqis, time_astra
double precision :: dr_factor_init_astra, dz_factor_init_astra ! factors of dr and dz for initial iterations
double precision :: raxis_astra, zaxis_astra, psi0_astra, psib_astra, sigma_B, sigma_axis, & 
   sigma_xpoint, r_xpoint_fit(5), z_xpoint_fit(5), sigma_energy, sigma_forces, sigma_limits   ! sigma_B multiplies the boundary, sigma_axis the axis, sigma_energy the block (sum sigma_coil coil_cur**2 induc), sigma_forces multiplies the force block: sum_ij force_ij I_i I_j. sigma_xpoint can be up to 5 x points to fit.
integer :: n_xpoint_fit

double precision :: vloop_avg, L_ext, dIp_dt   ! use tau_gseq_feqis here for refit mode 818
 
double precision :: x_point_save(20, 2) ! R, Z of xpoints, max 20 x points
double precision, dimension(ncoil_dim) :: activate_coil_feqis, cur_init, sigma_coils, sigma_coils_psiext ! initial currents from astra exp, not from coil.dat, in MA/turn
double precision, dimension(ncoil_dim) :: new_resistance ! whichever is > 0, it is used as new resistance.
double precision, dimension(ncoil_dim, 2) :: current_limit_feqis ! 1 is upper, 2 is lower
double precision, dimension(ncoil_dim, ncoil_dim) :: force_coil ! where it is 1, forces coil i,i to current of i,j
character(len=80) :: machine_description ! name of device, in astra it's called MACHINE

!below, allocatable
integer, dimension(:), allocatable :: which_x_point !same length of n_isoflux. where 0 --> continous point, where 1 --> x point
double precision, dimension(:), allocatable :: r_isoflux, z_isoflux
double precision, dimension(:,:), allocatable :: voltage_limits_active_coils

contains

    subroutine transport2fbe_init

    use_limiter_astra = 1  ! 1-uses limiter, 0-ignore limiter
    refit_mode = 0
    solve_fix = 0
    execute_plasma = 1     ! if 0 - only circuit equations, if 1 - solve plasma gseq too 
    nonegcurr = 0          ! nonegcurr = 0 --> no negative current allowed in plasma
    fast_mode = 0          ! 0 - normale, 1 - domnt do iterationsin fbe gse
    reconnect_circuits = 0 ! 0-nothing, 1-recompute matrix with new circuits
    new_equivalence = 0
    resistance_change = 0  ! if  = 1, changes resistances of coils using new_resistance
    n_equivalence = 0      ! number of equivalences
    use_reduce_circuit = 0 ! 0-all coils solved. 1 - some coils not solved
    n_of_newton_iterations = 150    ! to find actual mag axis. recommended between 5 - 10 
    n_fourier_restab_boundary = 8 ! nr of fourier modes for boundary restab, default = 5
    psplex_from_fbe = 0    ! put 1 to get psplex fromfree boundary
    plasma_config = 1      ! 0 if limiter, 1 if xpoint
    simple_plasma_model_breakdown = 0 ! if 0, no plasma feedback to coils in vacuum, if 1 yes

    use_isoflux = 0
    n_isoflux = 0

    tau_circuit_feqis = TAU
    tau_gseq_feqis = TAU
    time_astra = TSTART

    dr_factor_init_astra = 1.
    dz_factor_init_astra = 1.
    raxis_astra = RTOR + SHIFT
    zaxis_astra = UPDWN
    psi0_astra = PSIAX
    psib_astra = PSIBO
    sigma_B = 1.
    sigma_axis = 1.
    sigma_xpoint = 1.
    r_xpoint_fit = 0.
    z_xpoint_fit = 0.
    sigma_energy = 0.001
    sigma_forces = 0.
    sigma_limits = 0.  
    n_xpoint_fit = 0
    vloop_avg = 0.
    L_ext = 0.
    dIp_dt = 0.
    x_point_save = 0. ! R, Z of xpoints, max 20 x points

    activate_coil_feqis = 1

    cur_init = 0.
    cur_init(1: NCNB) = CCOILX(NCNBT+1: NCNBT+NCNB)/1.e3
    sigma_coils = 1.
    sigma_coils_psiext = 0. ! initial currents from astra exp, not from coil.dat, in MA/turn
    new_resistance = 0. ! whichever is > 0, it is used as new resistance.
    current_limit_feqis(:, 1) =  1.e6 ! 1 is upper, 2 is lower
    current_limit_feqis(:, 2) = -1.e6 ! 1 is upper, 2 is lower
    force_coil = 0 ! where it is 1, forces coil i,i to current of i,j
    machine_description = trim(MACHINE(1:4))

    return
    end subroutine transport2fbe_init

end module transport2fbe
    
