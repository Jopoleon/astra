module transport2fbe  !these are coupling variables with the equilibrium solver and astra

implicit none

integer, parameter :: ncoil_dim=300

integer :: use_limiter        ! 1-uses limiter, 0-ignore limiter
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
integer, parameter :: n_x_point=20    ! number of x points to be saved
integer :: n_coils

double precision :: tau_circuit, tau_gseq, time_astra
double precision :: dr_factor_init, dz_factor_init ! factors of dr and dz for initial iterations
double precision :: raxis_astra, zaxis_astra, psi0_astra, psib_astra, sigma_B, sigma_axis, & 
   sigma_xpoint, r_xpoint_fit(n_x_point), z_xpoint_fit(n_x_point), sigma_energy, sigma_forces, sigma_limits   ! sigma_B multiplies the boundary, sigma_axis the axis, sigma_energy the block (sum sigma_coil coil_cur**2 induc), sigma_forces multiplies the force block: sum_ij force_ij I_i I_j. sigma_xpoint can be up to 5 x points to fit.
integer :: n_xpoint_fit
integer :: fix_shape_after_fbe_off

double precision :: vloop_avg, L_ext, dIp_dt   ! use tau_gseq here for refit mode 818
 
double precision :: x_point_save(n_x_point, 2) ! R, Z of xpoints, max n_x_point x points
double precision, dimension(ncoil_dim) :: activate_coil, cur_init, sigma_coils, sigma_coils_ref ! initial currents from astra exp, not from coil.dat, in MA/turn
double precision, dimension(ncoil_dim) :: new_resistance ! whichever is > 0, it is used as new resistance.
double precision, dimension(ncoil_dim, 2) :: current_limit ! 1 is upper, 2 is lower
double precision, dimension(ncoil_dim, ncoil_dim) :: force_coil ! where it is 1, forces coil i,i to current of i,j
character(len=80) :: machine_description ! name of device, in astra it's called MACHINE

!below, allocatable
integer, dimension(:), allocatable :: which_x_point !same length of n_isoflux. where 0 --> continous point, where 1 --> x point
double precision, dimension(:), allocatable :: r_isoflux, z_isoflux, sigma_isoflux !sigma_isoflux is weight of isoflux point to boundary
double precision, dimension(:,:), allocatable :: voltage_limits_active_coils

contains

    subroutine transport2fbe_init(tau_in, tstart_in, R_in, updown_in, shift_in, &
        psi0_in, psib_in, machine_name, ncoils, coil_currents)

    integer, intent(in) :: ncoils
    double precision, intent(in) :: tau_in, tstart_in, R_in, updown_in, shift_in, psi0_in, psib_in
    double precision, intent(in), dimension(ncoils) :: coil_currents
    character(len=4), intent(in) :: machine_name

    use_limiter = 1
    refit_mode = 0
    solve_fix = 0
    execute_plasma = 1
    nonegcurr = 0
    fast_mode = 0
    reconnect_circuits = 0
    new_equivalence = 0
    resistance_change = 0
    n_equivalence = 0
    use_reduce_circuit = 0
    n_of_newton_iterations = 150
    n_fourier_restab_boundary = 8
    psplex_from_fbe = 0
    plasma_config = 0
    simple_plasma_model_breakdown = 0

    use_isoflux = 0
    n_isoflux = 0

    tau_circuit = tau_in
    tau_gseq = tau_in
    time_astra = tstart_in

    dr_factor_init = 1.
    dz_factor_init = 1.
    raxis_astra = R_in + shift_in
    zaxis_astra = updown_in
    psi0_astra = psi0_in
    psib_astra = psib_in
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
    x_point_save = 0.

    activate_coil = 1
    fix_shape_after_fbe_off = 1

    cur_init = 0.
    n_coils = ncoils
    cur_init(1: n_coils) = coil_currents/1.e3
    sigma_coils = 1.
    sigma_coils_ref = 1.
    new_resistance = 0.
    current_limit(:, 1) =  1.e6
    current_limit(:, 2) = -1.e6
    force_coil = 0
    machine_description = trim(machine_name(1:4))

    end subroutine transport2fbe_init

end module transport2fbe
