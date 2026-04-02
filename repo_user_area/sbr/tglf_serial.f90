subroutine tglf_serial(chi_i, chi_e, e_pflux, vimp1, vimp2, i_mflux_as, exchi_as, gamma_as, omega_as)

use tglf_interface, only: nsm, tglf_zs_in, tglf_ns_in, tglf_mass_in, &
    tglf_find_width_in, tglf_iflux_in, tglf_use_bper_in, tglf_use_mhd_rule_in, &
    tglf_use_bisection_in, tglf_use_inboard_detrapped_in, tglf_new_eikonal_in, &
    tglf_adiabatic_elec_in, tglf_ibranch_in, tglf_use_bpar_in, tglf_nmodes_in, &
    tglf_nbasis_max_in, tglf_nbasis_min_in, tglf_nxgrid_in, tglf_nky_in, &
    tglf_units_in, tglf_path_in, tglf_use_transport_model_in, &
    tglf_use_ave_ion_grid_in, tglf_sign_Bt_in, tglf_ky_in, tglf_width_in, &
    tglf_width_min_in, tglf_nwidth_in, tglf_geometry_flag_in, tglf_dump_flag_in, &
    tglf_test_flag_in, tglf_nn_max_error_in, tglf_write_wavefunction_flag_in, &
    tglf_theta_trapped_in, tglf_wdia_trapped_in, tglf_park_in, tglf_ghat_in, &
    tglf_gchat_in, tglf_sign_It_in, tglf_wd_zero_in, tglf_linsker_factor_in, &
    tglf_gradB_factor_in, tglf_filter_in, tglf_damp_psi_in, tglf_damp_sig_in, &
    tglf_kx0_loc_in, tglf_alpha_e_in, tglf_alpha_p_in, tglf_alpha_quench_in, &
    tglf_alpha_zf_in, tglf_xnu_factor_in, tglf_debye_factor_in, &
    tglf_etg_factor_in, tglf_sat_rule_in, tglf_kygrid_model_in, &
    tglf_xnu_model_in, tglf_vpar_model_in, tglf_vpar_shear_model_in, &
    tglf_b_model_sa_in, tglf_ft_model_sa_in, tglf_as_in, tglf_taus_in, &
    tglf_rlns_in, tglf_rlts_in, tglf_vpar_in, tglf_vpar_shear_in, &
    tglf_alpha_mach_in, tglf_vexb_shear_in, tglf_vexb_in, tglf_betae_in, &
    tglf_xnue_in, tglf_zeff_in, tglf_debye_in, tglf_rmin_loc_in, &
    tglf_rmaj_loc_in, tglf_zmaj_loc_in, tglf_drmindx_loc_in, tglf_drmajdx_loc_in, &
    tglf_dzmajdx_loc_in, tglf_kappa_loc_in, tglf_s_kappa_loc_in, &
    tglf_delta_loc_in, tglf_s_delta_loc_in, tglf_zeta_loc_in, tglf_s_zeta_loc_in, &
    tglf_q_loc_in, tglf_q_prime_loc_in, tglf_p_prime_loc_in, tglf_rmin_sa_in, &
    tglf_rmaj_sa_in, tglf_q_sa_in, tglf_shat_sa_in, tglf_alpha_sa_in, &
    tglf_xwell_sa_in, tglf_theta0_sa_in, file_dump_local, &
    tglf_elec_eflux_out, tglf_ion_eflux_out, tglf_ion_mflux_out, &
    tglf_elec_pflux_out, tglf_ion_pflux_out, tglf_elec_expwd_out, &
    tglf_q_elite_in, tglf_q_prime_elite_in, tglf_p_prime_elite_in, &
    tglf_n_ELITE_in, tglf_R_ELITE_in, tglf_Z_ELITE_in, tglf_Bp_ELITE_in
use tglf_pkg, only: get_eigenvalue_spectrum_out, get_ky_spectrum_out, &
     get_flux_spectrum_out

use parameter_inc, only: NRD
use const_inc, only: BTOR, RTOR, GP2, AMJ, AIM1, AIM2, AIM3, ZMJ, NA1
use status_inc, only: NE, TE, NI, TI, ZEF, PBLON, PBPER, PFAST, &
    ZIM1, ZIM2, ZIM3, NIZ1, NIZ2, NIZ3, ER, MU, FP_NORM, &
    RHO, AMETR, SHIF, ELON, NDEUT, TRIA, VTOR, G11, VPOL, VRS
use parameters_a2equil, only: equil_now
use numerical_tools, only: qinterp

implicit none

logical, parameter :: debug_elite=.false.
integer, parameter :: nrho_m=64, nthe_elite=400, mpol=6
double precision, parameter :: &
   k0   = 1.6022E-12, &       ! erg/ev
   e0   = 4.8032E-10, &       ! elementary charge (statcoulombs)
   e00  = 1.6020e-19, &       ! elementary charge (C)
   c0   = 2.9979E+10, &       ! speed of light (cm/sec)
   mp   = 1.6726E-24, &       ! proton mass (g)
   mpp  = 1.6726E-27, &       ! proton mass (kg)
   pi   = 3.141592653589793
double precision, parameter :: c_vpol=1.d0

! Parameters used:
! In A_vars & A_arrs:

double precision, dimension(NRD), intent(out) :: chi_i, chi_e, e_pflux, vimp2, vimp1, &
    i_mflux_as, exchi_as, gamma_as, omega_as

!----------------------------------------------------------------------

integer :: jr, jrho, jr_r, jr_l, jgamma_max, jion, kyloop, mom_order
integer :: sat_rule           ! Saturation rule
integer :: geom_flag          ! 1: Miller; 2: Fourier; 3: ELITE
integer :: nmodes_tg          ! number of unstable modes to use in computing fluxes (max=4)
integer :: kygrid_model_tg    ! select version of ky-grid to use 1
integer :: xnu_model_tg       ! select version of trapped-passing 2
integer :: jthe, jthe_rev, nrho_equ, nthe_equ ! for ELITE geometry
integer :: t_wall1, t_wall2, rate

double precision :: bmod, bpolz, alpha_zf_in, ion_eflux, ion_mflux, xstep, rho_min, rho_max
double precision :: dtheta_elite, drmin, drmaj, drho, dti, dte, dne, dq, dptot, &
        delong, dtrian, dvper, drhodr, dstep, dr, dv_r
double precision :: Bunit_gauss, Bunit_T, cs0, cs00, rhos0, omega0, rhostar2, lnlamda, taue, cexb
double precision :: a0_cm, a0_m, T0, N0, m0, rmin_m, drho_cs, drho_nt, nt_cs, gradrhosq_inv
double precision :: wdia_trap_tg          ! parameter for trapped fraction model

double precision, dimension(NRD) :: rmaj_as, q_as, ni_main_as, &
    vexb_as, vpar_as, vper_as, mtori_as, chie_as, chii_as, e_pflux_as, ptot_as
double precision, dimension(nrho_m) :: mtori, chie, chii, exchi, epflux, rho_m, &
    gamma_max, omega_max, kymax, ti_m, te_m, ne_m, vpar_m, vper_m, vexb_m, &
    ametr_m, elon_m, tria_m, rmaj_m, ptot_m, q_m, zef_m, pfn_m
double precision, dimension(nthe_elite) :: theta_elite, RR_elite, ZZ_elite, Bp_elite
double precision, allocatable, dimension(:) :: gamma, omega, kyspectrum, efluxspectrum
double precision, dimension(nsm-1) :: dni
double precision, dimension(nsm-1, nrho_m) :: ni_m, i_pflux
double precision, dimension(nsm-2, nrho_m) :: zimp_m
double precision, dimension(nsm-1, NRD) :: ni_as, i_pflux_as
! ELITE
double precision, allocatable, dimension(:) :: theta_equ, pfn_equ
double precision, allocatable, dimension(:, :) :: RR_tg, ZZ_tg, Bp_tg
character(len=120) :: f_elite

call SYSTEM_CLOCK(t_wall1, rate)

! Interpolate from ASTRA grid to TGLF grid
rho_min = RHO(1)
rho_max = RHO(NA1)
xstep = (rho_max - rho_min)/(nrho_m - 1.)
rho_m = (/ (rho_min + (jr - 1.)*xstep, jr=1, nrho_m) /)

call qinterp(RHO(1:NA1), ZIM1(1:NA1), NA1, rho_m, zimp_m(1, :), nrho_m)
call qinterp(RHO(1:NA1), ZIM2(1:NA1), NA1, rho_m, zimp_m(2, :), nrho_m)
call qinterp(RHO(1:NA1), ZIM3(1:NA1), NA1, rho_m, zimp_m(3, :), nrho_m)
call qinterp(RHO(1:NA1), NIZ1(1:NA1), NA1, rho_m,   ni_m(2, :), nrho_m)
call qinterp(RHO(1:NA1), NIZ2(1:NA1), NA1, rho_m,   ni_m(3, :), nrho_m)
call qinterp(RHO(1:NA1), NIZ3(1:NA1), NA1, rho_m,   ni_m(4, :), nrho_m)
call qinterp(RHO(1:NA1),      TI(1:NA1), NA1, rho_m,    ti_m, nrho_m)
call qinterp(RHO(1:NA1),      TE(1:NA1), NA1, rho_m,    te_m, nrho_m)
call qinterp(RHO(1:NA1),      NE(1:NA1), NA1, rho_m,    ne_m, nrho_m)
call qinterp(RHO(1:NA1),     ZEF(1:NA1), NA1, rho_m,   zef_m, nrho_m)
call qinterp(RHO(1:NA1),   AMETR(1:NA1), NA1, rho_m, ametr_m, nrho_m)
call qinterp(RHO(1:NA1),    ELON(1:NA1), NA1, rho_m,  elon_m, nrho_m)
call qinterp(RHO(1:NA1),    TRIA(1:NA1), NA1, rho_m,  tria_m, nrho_m)
call qinterp(RHO(1:NA1), FP_NORM(1:NA1), NA1, rho_m,   pfn_m, nrho_m)

do jrho=1, NA1
    if (NDEUT(jrho) >= 0.01*NE(jrho)) then
        ni_main_as(jrho) = NDEUT(jrho)
    else ! likely: NDEUT not defined in equ file, hence zero
        ni_main_as(jrho) = NI(jrho)
    endif
    rmaj_as(jrho) = RTOR + SHIF(jrho)
    q_as(jrho)    = 1./MU(jrho)
    ptot_as(jrho) = NE(jrho)*TE(jrho) + ni_main_as(jrho)*TI(jrho) + NIZ1(jrho)*TI(jrho) + pfast(jrho) + 0.5*(pblon(jrho) + pbper(jrho))
    bpolz = BTOR*AMETR(jrho)*MU(jrho)/RTOR
    bmod = sqrt(BTOR**2 + bpolz**2)
    vper_as(jrho) = ER(jrho)/(RTOR*bpolz) ! vexb in m/s --> Omega_E
    vpar_as(jrho) = VTOR(jrho) * BTOR/bmod + c_vpol* VPOL(jrho) * bpolz/bmod
    vexb_as(jrho) = -ER(jrho)/bmod ! vexb in m/s (vperp = vexb since the diamagnetic velocity is the curvature drift ac
enddo

call qinterp(RHO(1:NA1), ni_main_as(1:NA1), NA1, rho_m, ni_m(1, :), nrho_m)
call qinterp(RHO(1:NA1),    rmaj_as(1:NA1), NA1, rho_m,     rmaj_m, nrho_m)
call qinterp(RHO(1:NA1),       q_as(1:NA1), NA1, rho_m,        q_m, nrho_m)
call qinterp(RHO(1:NA1),    ptot_as(1:NA1), NA1, rho_m,     ptot_m, nrho_m)
call qinterp(RHO(1:NA1),    vpar_as(1:NA1), NA1, rho_m,     vpar_m, nrho_m)
call qinterp(RHO(1:NA1),    vper_as(1:NA1), NA1, rho_m,     vper_m, nrho_m)
call qinterp(RHO(1:NA1),    vexb_as(1:NA1), NA1, rho_m,     vexb_m, nrho_m)

! Electrons and main ions
tglf_zs_in(1) = -1.
tglf_zs_in(2) = ZMJ

tglf_mass_in(1) = 5.4447e-4/AMJ
tglf_mass_in(2) = AMJ/AMJ  ! AMJ is reference mass
tglf_mass_in(3) = AIM1/AMJ
tglf_mass_in(4) = AIM2/AMJ
tglf_mass_in(5) = AIM3/AMJ
do jr=1, nrho_m
    ni_m(2, jr) = max(1.e-9, ni_m(2, jr))
    ni_m(3, jr) = max(1.e-9, ni_m(3, jr))
    ni_m(4, jr) = max(1.e-9, ni_m(4, jr))
enddo

m0 = AMJ*mp          ! Ref. mass = D ion mass [g]
a0_m  = AMETR(NA1)
a0_cm = 1.d2*a0_m    ! length scale used by GYRO, m -> cm

e_pflux_as = 0.
i_pflux_as = 0.
chie_as  = 0.
chii_as  = 0.
mtori_as = 0.
exchi_as = 0.

! Number of species

tglf_ns_in = nsm

! These will be reset locally in the radial loop
tglf_zs_in(3) = MAXVAL(ZIM1(1:NA1))
tglf_zs_in(4) = MAXVAL(ZIM2(1:NA1))
tglf_zs_in(5) = MAXVAL(ZIM3(1:NA1))

if (tglf_zs_in(5) >= 1.) then
    tglf_ns_in = 5
else
    tglf_ns_in = 4
endif
if (tglf_zs_in(4) < 1.) tglf_ns_in = 3
if (tglf_zs_in(5) >= 1. .and. tglf_ns_in == 3) then
    tglf_ns_in = 4
endif
if (tglf_zs_in(3) < 1.) tglf_ns_in = 2
if (tglf_zs_in(4) >= 1. .and. tglf_ns_in == 2) then
    tglf_ns_in = 3
endif

kygrid_model_tg = 4 !1 Email Angioni Aug 1st 2023

sat_rule = 2
geom_flag = 1

if (geom_flag == 3) then
    tglf_n_elite_in  = nthe_elite - 1
    nrho_equ = SIZE(equil_now%coord_sys%position%r, dim=1)
    nthe_equ = SIZE(equil_now%coord_sys%position%r, dim=2)
    allocate(pfn_equ(nrho_equ))
    allocate(theta_equ(nthe_equ))
    allocate(RR_tg(nrho_m, nthe_equ), ZZ_tg(nrho_m, nthe_equ), Bp_tg(nrho_m, nthe_equ))

! Interpolation on TGLF rho-grid

    pfn_equ = (equil_now%profiles_1d%psi - equil_now%profiles_1d%psi(1))/(equil_now%profiles_1d%psi(nrho_equ) - equil_now%profiles_1d%psi(1))

! Interpolation on TGLF rho-grid
    do jthe=1, nthe_equ
        call qinterp(pfn_equ, equil_now%coord_sys%position%r(:, jthe), nrho_equ, pfn_m, RR_tg(:, jthe), nrho_m)
        call qinterp(pfn_equ, equil_now%coord_sys%position%z(:, jthe), nrho_equ, pfn_m, Zz_tg(:, jthe), nrho_m)
        call qinterp(pfn_equ, equil_now%coord_sys%bpcell    (:, jthe), nrho_equ, pfn_m, Bp_tg(:, jthe), nrho_m)
    enddo

    deallocate(pfn_equ)

    theta_equ = equil_now%coord_sys%position%theta2d
    dtheta_elite = GP2/dble(nthe_elite-1)
    theta_elite = (/ ((jthe - 1.)*dtheta_elite, jthe=1, nthe_elite) /)
endif

write(6, '(A, 5i4)') 'Call TGLF...', NA1, nrho_m, sat_rule, geom_flag, tglf_ns_in

SELECT CASE(sat_rule)
CASE(0)
    nmodes_tg = 2
    xnu_model_tg = 2
    wdia_trap_tg = 0.
    alpha_zf_in  = 0.
CASE(1)
    nmodes_tg = tglf_ns_in + 2
    xnu_model_tg = 2
    wdia_trap_tg = 0.
    alpha_zf_in  = 1.
CASE(2)
    nmodes_tg = tglf_ns_in + 2
    xnu_model_tg = 3
    wdia_trap_tg = 1.
    alpha_zf_in  = 1.
END SELECT

!-------------------------
! General TGLF settings

tglf_find_width_in     = .True.
tglf_iflux_in          = .True.
tglf_use_bper_in       = .True.
tglf_use_bpar_in       = .False.
tglf_use_mhd_rule_in   = .False.
tglf_use_bisection_in  = .True.
tglf_use_inboard_detrapped_in = .False.
tglf_new_eikonal_in    = .True.
tglf_adiabatic_elec_in = .False.
tglf_ibranch_in    = -1
tglf_nmodes_in     = nmodes_tg
tglf_nbasis_max_in = 6 ! email Angioni Aug 1st 2023, 4 old default
tglf_nbasis_min_in = 2
tglf_nxgrid_in     = 24 ! (24 email Gary, for ELITE), 16 old default
tglf_nky_in        = 19
tglf_units_in   = 'CGYRO'
tglf_path_in = '../tglf/'
! Want fluxes from TGLF
tglf_use_transport_model_in = .true.
tglf_use_ave_ion_grid_in    = .true. ! Email Angioni Aug 1st 2023

tglf_sign_Bt_in = 1
tglf_sign_It_in = 1

tglf_ky_in        = 0.3
tglf_width_in     = 1.65
tglf_width_min_in = 0.3
tglf_nwidth_in    = 21

tglf_geometry_flag_in = geom_flag
tglf_dump_flag_in     = .False.   ! Dumps input file
tglf_test_flag_in     = 0
tglf_nn_max_error_in  = 0
tglf_write_wavefunction_flag_in = 0 ! Writes eigenfunction

tglf_theta_trapped_in  = 0.7
tglf_wdia_trapped_in   = wdia_trap_tg
tglf_park_in           = 1.
tglf_ghat_in           = 1.
tglf_gchat_in          = 1.
tglf_wd_zero_in        = 0.1
tglf_linsker_factor_in = 0.
tglf_gradB_factor_in   = 0.
tglf_filter_in         = 2.
tglf_damp_psi_in       = 0.
tglf_damp_sig_in       = 0.
tglf_kx0_loc_in        = 0.

tglf_alpha_e_in          = 1.
tglf_alpha_p_in          = 1.
tglf_alpha_mach_in       = 0.
tglf_alpha_quench_in     = 0.
tglf_alpha_zf_in         = alpha_zf_in
tglf_xnu_factor_in       = 1.
tglf_debye_factor_in     = 1.
tglf_etg_factor_in       = 1.25
tglf_sat_rule_in         = sat_rule
tglf_kygrid_model_in     = kygrid_model_tg
tglf_xnu_model_in        = xnu_model_tg
tglf_vpar_model_in       = 0
tglf_vpar_shear_model_in = 1

tglf_b_model_sa_in  = 1
tglf_ft_model_sa_in = 1

! local field averages
! Initialise to zero for non-calculated species

tglf_as_in   = 0.
tglf_taus_in = 0.
tglf_rlns_in = 0.
tglf_rlts_in = 0.
tglf_vpar_in = 0.
tglf_vpar_shear_in = 0.

allocate(gamma(tglf_nky_in))
allocate(omega(tglf_nky_in))
allocate(kyspectrum(tglf_nky_in))
allocate(efluxspectrum(tglf_nky_in))

radial_loop: do jr=1, nrho_m

!thermal impurities

    tglf_zs_in(3) = max(1., zimp_m(1, jr))
    tglf_zs_in(4) = zimp_m(2, jr)
    tglf_zs_in(5) = zimp_m(3, jr)

    if (tglf_zs_in(5) >= 1. .and. tglf_ns_in == 3) then
        tglf_zs_in(4) = tglf_zs_in(5)
        tglf_mass_in(4) = tglf_mass_in(5)
        ni_m(3, :) = ni_m(4, :)
    endif
    if (tglf_zs_in(4) >= 1. .and. tglf_ns_in == 2) then
        tglf_zs_in(3) = tglf_zs_in(4)
        tglf_mass_in(3) = tglf_mass_in(4)
        ni_m(2, :) = ni_m(3, :)
    endif

!    tglf_ns_in = 3

    tglf_mass_in(tglf_ns_in + 1:) = 0.
    tglf_zs_in  (tglf_ns_in + 1:) = 0.

! Differentials

    jr_r = jr + 1
    jr_l = jr - 1
    if (jr == 1) then
        jr_l = jr
    else if (jr == nrho_m) then
        jr_r = jr
    endif
    dstep = 1./dble(jr_r - jr_l)  ! 0.5 in between, 1 at the edges

    drmin  = dstep*(ametr_m(jr_r) - ametr_m(jr_l))
    drmaj  = dstep*( rmaj_m(jr_r) -  rmaj_m(jr_l))
    drho   = dstep*(  rho_m(jr_r) -   rho_m(jr_l))
    delong = dstep*( elon_m(jr_r) -  elon_m(jr_l))
    dtrian = dstep*( tria_m(jr_r) -  tria_m(jr_l))
    dptot  = dstep*( ptot_m(jr_r) -  ptot_m(jr_l)) * 1E3*1E13
    dte = dstep*(te_m(jr_r) - te_m(jr_l))
    dne = dstep*(ne_m(jr_r) - ne_m(jr_l))
    dq    = dstep*(q_m(jr_r) - q_m(jr_l))
    dvper = dstep*(vper_m(jr_r) - vper_m(jr_l))
    dti = dstep*(ti_m(jr_r) - ti_m(jr_l))
    do jion=1, tglf_ns_in-1
        dni(jion) = dstep*(ni_m(jion, jr_r) - ni_m(jion, jr_l))
    enddo
    dv_r = dstep* &
        (vpar_m(jr_r)/(rmaj_m(jr_r) + ametr_m(jr_r)) - &
         vpar_m(jr_l)/(rmaj_m(jr_l) + ametr_m(jr_l)))
    dr = drmin/a0_m    ! gradients w.r.t. minor radius even for s-alpha geometry
    drhodr = drho/drmin

! local field averages

    tglf_as_in(1)   = 1. ! ne is ref
    tglf_taus_in(1) = 1. ! Te is ref

! Log derivatives
    tglf_rlns_in(1) = -dne/(dr*ne_m(jr))
    tglf_rlts_in(1) = -dte/(dr*te_m(jr))

    do jion=1, tglf_ns_in-1
        tglf_as_in(jion+1)   = ni_m(jion, jr)/ne_m(jr)
        tglf_taus_in(jion+1) = ti_m(jr)/te_m(jr)
        tglf_rlns_in(jion+1) = -dni(jion)/(dr*ni_m(jion, jr))
        tglf_rlts_in(jion+1) = -dti/(dr*ti_m(jr))
    enddo
! Restore quasi-neutrality via main ions

    tglf_as_in(2)   = -1./tglf_zs_in(2)*(SUM(tglf_as_in*tglf_zs_in) - tglf_as_in(2)*tglf_zs_in(2))
    tglf_rlns_in(2) = -1./(tglf_as_in(2)*tglf_zs_in(2))*(SUM(tglf_rlns_in*tglf_as_in*tglf_zs_in) - tglf_rlns_in(2)*tglf_as_in(2)*tglf_zs_in(2))

! GYRO conventions

    N0  = 1E13*ne_m(jr)   ! density scale used by GYRO [1/cm**3]
    T0  = 1E3 *te_m(jr)   ! temperature scale used by GYRO
    Bunit_T = BTOR*drhodr*rho_m(jr)/ametr_m(jr)  ! Miller geometry magnetic field unit [gauss]
    Bunit_gauss = 1.d4*Bunit_T 

! derived units for the plasma

    cs0 = SQRT(k0*T0/m0)          ! thermal velocity unit cm/sec
    cs00 = SQRT(e00*T0/(AMJ*mpp)) ! thermal velocity unit m/sec
    omega0 = e0*Bunit_gauss/(m0*c0)     ! gyrofrequency unit 1/sec
    rhos0 = cs0/omega0            ! gyroradius unit cm
    lnlamda = 24.0 -0.5*LOG(tglf_as_in(1)*N0) + LOG(tglf_taus_in(1)*T0)
    taue = 3.44E5 * (tglf_taus_in(1)*T0)**1.5 / (tglf_as_in(1)*N0*lnlamda)  !  sec

    rmin_m = ametr_m(jr)/a0_m
    cexb = ametr_m(jr)/q_m(jr)

    tglf_vpar_shear_in(2) = -1E2*rmaj_m(jr)*dv_r/(dr*cs0)  !From m/s to cm/s for vpar
    tglf_vpar_shear_in(1) = tglf_vpar_shear_in(2)

    tglf_vpar_in(2) = 1E2*vpar_m(jr)/cs0
    tglf_vpar_in(1) = tglf_vpar_in(2)

    if (tglf_ns_in >= 3) then
        tglf_vpar_shear_in(3: tglf_ns_in) = tglf_vpar_shear_in(2)
        tglf_vpar_in      (3: tglf_ns_in) = tglf_vpar_in(2)
    endif

! local magnetic geometry

    rhostar2 = (rhos0/a0_cm)**2
    drho_cs = drhodr**2*a0_m*rhostar2*cs00
    drho_nt = drhodr*a0_m*rhostar2*N0*e00*T0/(1.e13*mpp)
    nt_cs = 0.001602*N0/1.e13*T0/1.e3*cs00/a0_m*rhostar2

! Share variables with tglf_run via module tglf_interface

    tglf_vexb_shear_in = -1E2*cexb*dvper/(dr*cs0) ! Waltz-Miller definition     !From m/s to cm/s for vexb  -dVexb/dr

! Initialise

    tglf_vexb_in  = 1E2*vexb_m(jr)/cs0
    tglf_betae_in = 8.0*pi*k0*N0*T0/Bunit_gauss**2
    tglf_xnue_in  = 0.75*SQRT(pi)*a0_cm/(taue*cs0)
    tglf_zeff_in  = zef_m(jr)
    tglf_debye_in = SQRT(k0*T0/(4.0*pi*N0*e0**2))/rhos0

    tglf_rmin_loc_in    = rmin_m
    tglf_rmaj_loc_in    = rmaj_m(jr)/a0_m
    tglf_zmaj_loc_in    = 0.
    tglf_drmindx_loc_in = 1.
    tglf_drmajdx_loc_in = drmaj/(dr*a0_m)
    tglf_dzmajdx_loc_in = 0.
    tglf_kappa_loc_in   = elon_m(jr)
    tglf_s_kappa_loc_in = ametr_m(jr)*delong/(drmin*elon_m(jr))
    tglf_delta_loc_in   = tria_m(jr)
    tglf_s_delta_loc_in = ametr_m(jr)*dtrian/drmin
    tglf_zeta_loc_in    = 0.
    tglf_s_zeta_loc_in  = 0.
    tglf_q_loc_in       = q_m(jr)
    tglf_q_prime_loc_in = (q_m(jr)/rmin_m)*dq/dr
    tglf_p_prime_loc_in = (k0/Bunit_gauss**2)*(q_m(jr)/rmin_m)*dptot/dr

    tglf_q_ELITE_in       = tglf_q_loc_in
    tglf_q_prime_ELITE_in = tglf_q_prime_loc_in
    tglf_p_prime_ELITE_in = tglf_p_prime_loc_in

    tglf_rmin_sa_in     = rmin_m
    tglf_rmaj_sa_in     = rmaj_m(jr)/a0_m
    tglf_q_sa_in        = q_m(jr)
    tglf_shat_sa_in     = (ametr_m(jr)/q_m(jr))*dq/drmin
    tglf_alpha_sa_in    = -(8.0*pi*k0/Bunit_gauss**2)*q_m(jr)**2 * rmaj_m(jr)*dptot/drmin
    tglf_xwell_sa_in    = 0.
    tglf_theta0_sa_in   = 0.

! Settings
    if (tglf_dump_flag_in) then
        write(file_dump_local, '(A11, I0)') 'input.tglf_', jr
    endif

    if (geom_flag == 3) then ! R, Z contours for ELITE
! Interpolation on ELITE theta-grid
        call qinterp(theta_equ, RR_tg(jr, :), nthe_equ, theta_elite, RR_elite, nthe_elite)
        call qinterp(theta_equ, Zz_tg(jr, :), nthe_equ, theta_elite, ZZ_elite, nthe_elite)
        call qinterp(theta_equ, Bp_tg(jr, :), nthe_equ, theta_elite, Bp_elite, nthe_elite)
        RR_elite(tglf_n_elite_in+1) = RR_elite(1)
        ZZ_elite(tglf_n_elite_in+1) = ZZ_elite(1)
        Bp_elite(tglf_n_elite_in+1) = Bp_elite(1)
        tglf_R_elite_in(1:tglf_n_elite_in+1)  = RR_elite/a0_m
        tglf_Z_elite_in(1:tglf_n_elite_in+1)  = ZZ_elite/a0_m
        tglf_Bp_elite_in(1:tglf_n_elite_in+1) = Bp_elite/Bunit_T

        if (debug_elite) then
            write(f_elite, '(A11, I0)') 'elite4tglf_', jr
            open(31, FILE=f_elite)
            do jthe=1, nthe_elite
                write(31, '(3F)') tglf_R_elite_in(jthe), tglf_Z_elite_in(jthe), Bp_elite(jthe)
            enddo
            close(31)
        endif
    endif

! ------ Call TGLF ------
    write(*, *) 'Calling TGLF jrho=', jr
    call tglf_run
! -----------------------

! Transport coefficients

    ion_eflux = SUM(tglf_ion_eflux_out(1: tglf_ns_in-1))
    ion_eflux = ion_eflux/(tglf_taus_in(2) * 1e13*ni_m(1, jr)/N0)
    ion_mflux = SUM(tglf_ion_mflux_out(1: tglf_ns_in-1))
    chii (jr) = ion_eflux          /(1e-4 + abs(tglf_rlts_in(2))) *drho_cs
    chie (jr) = tglf_elec_eflux_out/(1e-4 + abs(tglf_rlts_in(1))) *drho_cs
    mtori(jr) = ion_mflux*drho_nt
    epflux(jr) = tglf_elec_pflux_out/drhodr *drho_cs         ! particle flux
    do jion=1, tglf_ns_in-1
        i_pflux(jion, jr) = tglf_ion_pflux_out(jion)/drhodr *drho_cs  !ion particle flux
    enddo
    exchi(jr) = tglf_elec_expwd_out * nt_cs                ! Equipartition
    do kyloop=1, tglf_nky_in
        gamma(kyloop) = get_eigenvalue_spectrum_out(1, kyloop, 1)
        omega(kyloop) = get_eigenvalue_spectrum_out(2, kyloop, 1)
        kyspectrum(kyloop) = get_ky_spectrum_out(kyloop)
        efluxspectrum(kyloop) = get_flux_spectrum_out(2, 1, 1, kyloop, 1)
    enddo
    jgamma_max = maxloc(efluxspectrum(1:tglf_nky_in), 1)
    gamma_max(jr) = gamma(jgamma_max)*(cs0/a0_cm)
    omega_max(jr) = omega(jgamma_max)*(cs0/a0_cm)
    kymax(jr) = kyspectrum(jgamma_max)
enddo radial_loop

! Interpolate back to ASTRA radial grid

call qinterp(rho_m, chii     , nrho_m, RHO(1:NA1),    chii_as(1:NA1), NA1)
call qinterp(rho_m, chie     , nrho_m, RHO(1:NA1),    chie_as(1:NA1), NA1)
call qinterp(rho_m, mtori    , nrho_m, RHO(1:NA1), i_mflux_as(1:NA1), NA1)
call qinterp(rho_m, epflux   , nrho_m, RHO(1:NA1), e_pflux_as(1:NA1), NA1)
call qinterp(rho_m, exchi    , nrho_m, RHO(1:NA1),   exchi_as(1:NA1), NA1)
call qinterp(rho_m, gamma_max, nrho_m, RHO(1:NA1),   gamma_as(1:NA1), NA1)
call qinterp(rho_m, omega_max, nrho_m, RHO(1:NA1),   omega_as(1:NA1), NA1)
do jion=1, tglf_ns_in-1
    call qinterp(rho_m, i_pflux(jion, 1:nrho_m), nrho_m, RHO(1:NA1), i_pflux_as(jion, 1:NA1), NA1)
enddo

do jrho=1, NA1
    gradrhosq_inv = VRS(jrho)/G11(jrho)
    chi_i(jrho) = chii_as(jrho)*gradrhosq_inv ! \chi_i, m^2/s
    chi_e(jrho) = chie_as(jrho)*gradrhosq_inv ! \chi_e, m^2/s
    e_pflux(jrho) = e_pflux_as(jrho)*gradrhosq_inv/a0_m ! D flux
    vimp1(jrho) = i_pflux_as(2, jrho)*gradrhosq_inv/a0_m/(NIZ1(jrho)/NE(jrho))  ! 1st imp convection
    vimp2(jrho) = i_pflux_as(3, jrho)*gradrhosq_inv/a0_m/(NIZ2(jrho)/NE(jrho))  ! 2nd imp convection
enddo

call SYSTEM_CLOCK(t_wall2, rate)
print*, "TGLF_serial wall time", dble(t_wall2 - t_wall1)/dble(rate)

return
END subroutine tglf_serial
