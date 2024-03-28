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
endtype
type rawProfileMap
    integer, dimension(NTARR) :: arr_index=0, jbeg_grid=0, jbeg_data=0, grid_type=0, nrho=0
    double precision, dimension(NTARR) :: time=0., filter=0.001
    character(len=6), dimension(NTARR) :: label
endtype

type(rawScalar) :: raw_scalar
type(rawProfileMap) :: raw_profile_map

real*4 :: DATARR(NRDX*NTARR)
double precision, dimension(NBDTMAX) :: BNDTIM
double precision, dimension(NBDTMAX*NBDMAX) :: BNDR, BNDZ

end module expdat

!--------------------------------
module ac_neg1

implicit none

integer :: NUM(4), JMIN, JMAX, NKL1, NKL2, MODK(2)

end module ac_neg1

!--------------------------------
module plasma_state ! for flight simulator plasma yes/no

implicit none

integer :: plasma_up, plasma_trig

end module plasma_state

!--------------------------------
module flight_sim_geometrics ! for flight simulator diagnostics

implicit none

double precision, dimension(1000) :: geom1d

end module flight_sim_geometrics

!--------------------------------
module fenix_params !for flight simulator parameters

implicit none

double precision :: wallpos, d_j_m, dw_j_m, dw2_j_m, &
    vsoldiv, solwidth, &
    recycl_wall, boron_wall, predep_w, &
    transp_variance, hmodetransp, lmodetransp, &
    solmod_dt,  &
    saves_dt, savep_dt, neocl_dt, trmod_ty, nbieqmix_dt, &
    torba_dt, simdtmultip, simdtmultip_min, ped_width, chie_chii, &
    D_chie, D_ped_mult, gs2d_tmin_multip, timecirc, timepsi, &
    diohdt, diohdtthreshold, dteqz2, zibkdw, zifbey, &
    dt_adapt, VV, dt_fazt, ipl_bf_bkdw, &
    alp0, alpnew0, rx00, zx00, w_wall_source, pedtop_multiplier, &
    pfast_pressure_coef, denssource_fast, &
    tesep_multip, pellet_pos, pellet_width,tavg_psep, &
    Dped_H_multi, sh_crit_saw, ntm_m, ntm_n, ntm_seed, &
    elm_area, elm_res, elm_fbs, dioh_thresh, iprd_thresh, &
    eq_err_thres, pel_speed, nbi_model, ec_model, ic_model, pel_model, &
    coil_start_time(15)

integer :: lhmodel, btipdirec, pr_clamp, reinitcirc, reinitpsi, &
    ipsmk2, eq_cmd, use_zlim_pot, cmnd_dioh2s, cmnd_dioh2u, nequiz, &
    s_adapt, yesfitcc, s_fazt, kastr2,  &
    isafazt, ispid_contour, res_trigts06, resres_oh6, &
    max_max_iteri, max_max_iterb, max_max_iterj

double precision :: radial_error, vertical_error
double precision, dimension(500, 2) :: psitok

end module fenix_params

!--------------------------------
module fs_coupling_variables

!F2A
double precision, dimension(:), allocatable :: fs_voltages ! max 100 circuits
double precision, dimension(:), allocatable :: fs_pow_IC ! max 10 antennas
double precision, dimension(:), allocatable :: fs_pol_EC, fs_pow_EC, fs_pow_NB ! max 8 independent gyros
double precision, dimension(:), allocatable :: fs_pellet !max 10 injectors - nr of particles in 1e19 per injected pellet. First 2 positions are reserved for D and T.
double precision, dimension(:), allocatable :: fs_valves, fs_pumped_gases !valves for gas puff
double precision, dimension(:, :), allocatable :: fs_geometric_profiles  !500 x 13 max
double precision, dimension(:, :), allocatable :: fs_bnd_in ! 50 , 2 boundary values R,Z
integer :: fs_bnd_yes
double precision :: fs_next_sim_time, fs_pump, &
    fs_NTM_trig, fs_NTM_M, fs_NTM_N, fs_NTM_seed, fs_stop_time, fs_ipl_in, &
    fs_oh_switch

!A2F
double precision, dimension(:), allocatable :: fs_currents ! max 500 magnetics signals
double precision, dimension(:), allocatable :: fs_magnetics ! max 500 magnetics signals
double precision, dimension(:, :), allocatable :: fs_cforces ! 100 , 2 coil forces max
double precision :: fs_paux, fs_prad, fs_psep, fs_pintrinsic, fs_pfus, fs_ipl_out

!interpolants
double precision, dimension(:,:), allocatable :: fs_interpolants

end module fs_coupling_variables

!-------------------------------------
module astra2fbe  !these are coupling variables with the equilibrium solver and astra

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

double precision :: tau_circuit_feqis, tau_gseq_feqis, time_astra
double precision :: dr_factor_init_astra, dz_factor_init_astra ! factors of dr and dz for initial iterations
double precision :: raxis_astra, zaxis_astra, psi0_astra, psib_astra, sigma_B, sigma_axis, & 
 sigma_xpoint, r_xpoint_fit(5), z_xpoint_fit(5), sigma_energy, sigma_forces   ! sigma_B multiplies the boundary, sigma_axis the axis, sigma_energy the block (sum sigma_coil coil_cur**2 induc), sigma_forces multiplies the force block: sum_ij force_ij I_i I_j. sigma_xpoint can be up to 5 x points to fit.
integer :: n_xpoint_fit

double precision :: vloop_avg, L_ext, dIp_dt   ! use tau_gseq_feqis here for refit mode 818
 
double precision :: x_point_save(20, 2) ! R, Z of xpoints, max 20 x points
double precision, dimension(ncoil_dim) :: activate_coil_feqis, cur_init, sigma_coils, sigma_coils_psiext ! initial currents from astra exp, not from coil.dat, in MA/turn
double precision, dimension(ncoil_dim) :: new_resistance ! whichever is > 0, it is used as new resistance.
double precision, dimension(ncoil_dim, 2) :: current_limit_feqis ! 1 is upper, 2 is lower
double precision, dimension(ncoil_dim, ncoil_dim) :: force_coil ! where it is 1, forces coil i,i to current of i,j
character(len=80) :: machine_description ! name of device, in astra it's called MACHINE

end module astra2fbe
