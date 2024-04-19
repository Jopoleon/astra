module astra2fbe  !these are coupling variables with the equilibrium solver and astra

use const_inc, only: tau, time, RTOR, shift, psiax, psibo, iplx, taumin, NA1
use status_inc, only: SHIV
use outcmn_inc, only: machine, ccoil

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

integer :: use_isoflux, n_isoflux         ! 0 does nothing, 1 when mode 818 is used to reconstruct coil currents, the boundary is obtained on the isoflux points r_isoflux and z_isoflux , of length n_isoflux
integer, dimension(:), allocatable :: which_x_point !same length of n_isoflux. where 0 --> continous point, where 1 --> x point
double precision, dimension(:), allocatable :: r_isoflux, z_isoflux
double precision, dimension(:,:), allocatable :: voltage_limits_active_coils

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

contains
    subroutine astra2fbe_init

    tau_circuit_feqis = tau
    tau_gseq_feqis    = tau
    time_astra     = time
    activate_coil_feqis = 1 ! if 0, coil is disconnected if use_reduce_circuit and reconnect_circuits is used, otherwise just sets the current to zero (you will get a different result!)
    current_limit_feqis(:, 1) =  1.e6 ! 1 is upper, 2 is lower
    current_limit_feqis(:, 2) = -1.e6 ! 1 is upper, 2 is lower
    force_coil = 0. ! where it is 1, forces coil i,i to current of i,j --> better use reconnect circuits.
    
    use_reduce_circuit = 0 ! run with reconnected circuits, does not reset the matrix
    reconnect_circuits = 0 ! this resets the circuit matrix, to change connections
    n_equivalence = 0
    new_equivalence = 0
    !new_equivalence(1:14,1) = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    resistance_change = 0 ! this resets the circuit matrix, to change resistances (diagonals)
    new_resistance = 0.
    !	new_resistance(1:12) = (/0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./) !diagonal resistance in microOhm
    
    raxis_astra = RTOR + SHIFT
    zaxis_astra = SHIV(1)
    psi0_astra = PSIAX
    psib_astra = PSIBO
    n_fourier_restab_boundary = 5
    use_limiter_astra = 1   ! do not use limiter for DEMO
    refit_mode = 0   ! if -1 - 1 turn only, 0 - stab method, if 1 - restab with prescribed axis , 3 - full fit like spider but only for eddy currents, 101 - only Z stab
    sigma_B = 1.
    sigma_axis = 50000.
    sigma_coils = 1.e5
    sigma_energy = 1
    
    solve_fix = 0   ! if 0 - solve full fix boundary problem, >0 - N pass only, -2 - uses fbe solution 
    execute_plasma = 1   ! if 0 - only circuit equations, if 1 - solve plasma gseq too
    n_of_newton_iterations = 150
    
    !factors of dr and dz for initial iterations
    dr_factor_init_astra = 1.
    dz_factor_init_astra = 1.
    fast_mode = 0
    psplex_from_fbe = 0
    
    execute_plasma = 1
    tau_gseq_feqis = taumin ! spider GS solver time step
    
    if (MACHINE(1:3) == 'aug') then
        cur_init( 1:12) = CCOIL(1:12)/1.e3
        cur_init(13:52) = 0.
    endif
    
    if (TIME > 2.52) fast_mode = 1
    
    write(*, *) 'eqtime', time, fast_mode, execute_plasma, tau_gseq_feqis, cur_init(1:12)
    
    return
    end subroutine astra2fbe_init

end module astra2fbe
    
