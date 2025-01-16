subroutine tglf_serial(CHI, CHE, VIN, DPH, DPL, DPR, XTB, GM1, OM1)

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
    tglf_q_fourier_in, tglf_q_prime_fourier_in, tglf_p_prime_fourier_in, &
    tglf_nfourier_in, tglf_fourier_in, &
    tglf_q_elite_in, tglf_q_prime_elite_in, tglf_p_prime_elite_in, &
    tglf_n_ELITE_in, tglf_R_ELITE_in, tglf_Z_ELITE_in, tglf_Bp_ELITE_in
use tglf_pkg, only: get_eigenvalue_spectrum_out, get_ky_spectrum_out, &
     get_flux_spectrum_out

use parameter_inc, only: NRD
use const_inc, only: BTOR, RTOR, GP, GP2, ABC, &
    AMJ, AIM1, AIM2, AIM3, ZMJ, PSIAX, PSIBO, &
    NA1, NA1N, NA1E, NA1I
use status_inc, only: NE, TE, NI, TI, &
    ZEF, ZIM1, ZIM2, ZIM3, PBLON, PBPER, &
    PFAST, NIZ3, AMAIN, ER, MU, FP_NORM, &
    RHO, AMETR, SHIF, ELON, &
    NDEUT, NIZ1, NTRIT, NIZ2, NHE3, &
    TRIA, VTOR, NIBM, G11, VPOL, VRS, SHEAR
use parameters_a2equil, only: equil_now

implicit none

logical, parameter :: debug=.false.
integer, parameter :: nrho_tg=40, nthe_elite=400, mpol=6
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

double precision, dimension(*), intent(out) :: CHI, CHE, VIN, DPH, DPL, DPR, XTB, GM1, OM1

!----------------------------------------------------------------------

integer :: jr, jrho, jr_r, jr_l, jgamma_max, jspec, kyloop, mom_order
integer :: sat_rule           ! Saturation rule
integer :: geom_flag          ! 1: Miller; 2: Fourier; 3: ELITE
integer :: nmodes_tg          ! number of unstable modes to use in computing fluxes (max=4)
integer :: kygrid_model_tg    ! select version of ky-grid to use 1
integer :: xnu_model_tg       ! select version of trapped-passing 2
integer :: jthe, jthe_rev, nrho_equ, nthe_equ ! for ELITE geometry

double precision :: bmod, bpolz, alpha_zf_in, ion_eflux, ion_mflux, xstep, rho_min, rho_max
double precision :: dtheta_elite, drmin, drmaj, drho, dte, dne, dq, dptot, &
        delong, dtrian, dvper, drhodr, dstep, dr, dv_r
double precision :: Bunit_gauss, Bunit_T, cs0, cs00, rhos0, omega0, rhostar2, lnlamda, taue, cexb
double precision :: a0_cm, a0_m, T0, N0, m0, rmin_tg, drho_cs, drho_nt, nt_cs
double precision :: wdia_trap_tg          ! parameter for trapped fraction model

double precision, dimension(NRD) :: gradrhosq_exp, rmaj_exp, q_exp, &
    vexb_exp, vpar_exp, vper_exp, mtori_m, &
    chie_m, chii_m, elec_pflux_m, exchi_m, ptot_exp, gamma_m, omega_m
double precision, dimension(nrho_tg) :: mtori, chie, chii, exchi, elec_pflux, rho_tg, &
    gamma_max, omega_max, kymax, te_tg, ne_tg, vpar_tg, vper_tg, vexb_tg, &
    ametr_tg, pf_tg, pfn_tg, elon_tg, tria_tg, rmaj_tg, ptot_tg, q_tg, zef_tg
double precision, dimension(nthe_elite) :: theta_elite, cos_mthe, sin_mthe, da_elite, drda, dzda
double precision, allocatable, dimension(:) :: gamma, omega, kyspectrum, efluxspectrum, ifluxspectrum, pfluxspectrum
double precision, dimension(nsm-1) :: dti, dni
double precision, dimension(nsm-1, nrho_tg) :: ni_tg, ti_tg, z_tg, ion_pflux
double precision, dimension(nsm-1, NRD) :: ni_exp, ion_pflux_m
! Fourier moments
double precision, dimension(mpol) :: xrc, xrs, xzc, xzs
! ELITE
double precision, allocatable, dimension(:) :: pfn_equ, theta_equ, Bp_elite, RR_elite, ZZ_elite
double precision, allocatable, dimension(:, :) :: Bp_tg, RR_tg, ZZ_tg, da_tg
character(len=10) :: time_loc
character(len=120) :: f_fourier, f_elite

call DATE_AND_TIME(TIME=time_loc)
write(*, *) time_loc, ' BEGIN tglf_serial.f90'

! Interpolate from ASTRA grid to TGLF grid
rho_min = RHO(1)
!rho_max = RHO(NA1)
rho_max = max(RHO(NA1I), RHO(NA1E), RHO(NA1N))
xstep = (rho_max - rho_min)/(nrho_tg - 1.)
rho_tg = (/ (rho_min + (jr - 1.)*xstep, jr=1, nrho_tg) /)

call qinterp(RHO(1:NA1),     TI(1:NA1), NA1, rho_tg, ti_tg(1, :), nrho_tg)
call qinterp(RHO(1:NA1),     TE(1:NA1), NA1, rho_tg,       te_tg, nrho_tg)
call qinterp(RHO(1:NA1),   ZIM1(1:NA1), NA1, rho_tg,  z_tg(2, :), nrho_tg)
call qinterp(RHO(1:NA1),   ZIM2(1:NA1), NA1, rho_tg,  z_tg(3, :), nrho_tg)
call qinterp(RHO(1:NA1),   ZIM3(1:NA1), NA1, rho_tg,  z_tg(4, :), nrho_tg)
call qinterp(RHO(1:NA1),     NE(1:NA1), NA1, rho_tg,       ne_tg, nrho_tg)
call qinterp(RHO(1:NA1),    ZEF(1:NA1), NA1, rho_tg,      zef_tg, nrho_tg)
call qinterp(RHO(1:NA1),  AMETR(1:NA1), NA1, rho_tg,    ametr_tg, nrho_tg)
call qinterp(RHO(1:NA1),   ELON(1:NA1), NA1, rho_tg,     elon_tg, nrho_tg)
call qinterp(RHO(1:NA1),   TRIA(1:NA1), NA1, rho_tg,     tria_tg, nrho_tg)
call qinterp(RHO(1:NA1),FP_NORM(1:NA1), NA1, rho_tg,      pfn_tg, nrho_tg)
ti_tg(2, :) = ti_tg(1, :)
ti_tg(3, :) = ti_tg(1, :)
ti_tg(4, :) = ti_tg(1, :)

do jrho=1, NA1
    if (NDEUT(jrho) >= 0.01*NE(jrho)) then
        ni_exp(1, jrho) = NDEUT(jrho)
    else ! likely: NDEUT not defined in equ file, hence zero
        ni_exp(1, jrho) = NI(jrho)
    endif
    ni_exp(2, jrho) = NIZ1(jrho)
    ni_exp(3, jrho) = NIZ2(jrho)
    ni_exp(4, jrho) = NIZ3(jrho)
    rmaj_exp(jrho) = RTOR + SHIF(jrho)
    q_exp(jrho)    = 1./MU(jrho)
    ptot_exp(jrho) = NE(jrho)*TE(jrho) + ni_exp(1, jrho)*TI(jrho) + ni_exp(2, jrho)*TI(jrho) + pfast(jrho) + 0.5*(pblon(jrho) + pbper(jrho))
    bpolz = BTOR*AMETR(jrho)*MU(jrho)/RTOR
    bmod = sqrt(BTOR**2 + bpolz**2)
    gradrhosq_exp(jrho) = G11(jrho)/VRS(jrho)
    vper_exp(jrho) = ER(jrho)/(RTOR*bpolz) ! vexb in m/s --> Omega_E
!    vpar_m(jrho) = ER(jrho)/(RTOR*bpolz)*(RTOR+SHIF(jrho)+AMETR(jrho))  !--> R*Omega_E , no neoclassical terms
    vpar_exp(jrho) = VTOR(jrho) * BTOR/bmod + c_vpol* VPOL(jrho) * bpolz/bmod
    vexb_exp(jrho)  = -ER(jrho)/bmod ! vexb in m/s (vperp = vexb since the diamagnetic velocity is the curvature drift ac
enddo

call qinterp(RHO(1:NA1), ni_exp(1, 1:NA1), NA1, rho_tg, ni_tg(1, :), nrho_tg)
call qinterp(RHO(1:NA1), ni_exp(2, 1:NA1), NA1, rho_tg, ni_tg(2, :), nrho_tg)
call qinterp(RHO(1:NA1), ni_exp(3, 1:NA1), NA1, rho_tg, ni_tg(3, :), nrho_tg)
call qinterp(RHO(1:NA1), ni_exp(4, 1:NA1), NA1, rho_tg, ni_tg(4, :), nrho_tg)
call qinterp(RHO(1:NA1), rmaj_exp(1:NA1), NA1, rho_tg,  rmaj_tg, nrho_tg)
call qinterp(RHO(1:NA1),    q_exp(1:NA1), NA1, rho_tg,     q_tg, nrho_tg)
call qinterp(RHO(1:NA1), ptot_exp(1:NA1), NA1, rho_tg,  ptot_tg, nrho_tg)
call qinterp(RHO(1:NA1), vpar_exp(1:NA1), NA1, rho_tg,  vpar_tg, nrho_tg)
call qinterp(RHO(1:NA1), vper_exp(1:NA1), NA1, rho_tg,  vper_tg, nrho_tg)
call qinterp(RHO(1:NA1), vexb_exp(1:NA1), NA1, rho_tg,  vexb_tg, nrho_tg)

! Electrons and main ions
tglf_zs_in(1) = -1.
tglf_zs_in(2) = ZMJ

tglf_mass_in(1) = 5.4447e-4/AMJ
tglf_mass_in(2) = AMJ/AMJ  ! AMJ is reference mass
tglf_mass_in(3) = AIM1/AMJ
tglf_mass_in(4) = AIM2/AMJ
tglf_mass_in(5) = AIM3/AMJ
do jr=1, nrho_tg
    ni_tg(2, jr) = max(1.e-9, ni_tg(2, jr))
    ni_tg(3, jr) = max(1.e-9, ni_tg(3, jr))
    ni_tg(4, jr) = max(1.e-9, ni_tg(4, jr))
enddo

m0 = AMJ*mp          ! Ref. mass = D ion mass [g]
a0_m  = AMETR(NA1)
a0_cm = 1.d2*a0_m    ! length scale used by GYRO, m -> cm

elec_pflux_m = 0.
ion_pflux_m  = 0.
chie_m  = 0.
chii_m  = 0.
mtori_m = 0.
exchi_m = 0.

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
geom_flag = 3

if (geom_flag == 2 .or. geom_flag == 3) then
    tglf_nfourier_in = mpol
    tglf_n_elite_in  = nthe_elite - 1
    nrho_equ = SIZE(equil_now%coord_sys%position%r, dim=1)
    nthe_equ = SIZE(equil_now%coord_sys%position%r, dim=2)
    allocate(theta_equ(nthe_equ))
    allocate(pfn_equ(nrho_equ))
    allocate(RR_tg(nrho_tg, nthe_equ), ZZ_tg(nrho_tg, nthe_equ), Bp_tg(nrho_tg, nthe_equ), da_tg(nrho_tg, nthe_equ))
    allocate(RR_elite(nthe_elite), ZZ_elite(nthe_elite), Bp_elite(nthe_elite))

    theta_equ = equil_now%coord_sys%position%teta2d
    dtheta_elite = GP2/dble(nthe_elite-1)
    theta_elite = (/ ((jthe - 1.)*dtheta_elite, jthe=1, nthe_elite) /)

    pfn_equ = (equil_now%profiles_1d%psi - equil_now%profiles_1d%psi(1))/(equil_now%profiles_1d%psi(nrho_equ) - equil_now%profiles_1d%psi(1))

! Interpolation on TGLF rho-grid
    do jthe=1, nthe_equ
        call qinterp(pfn_equ, equil_now%coord_sys%position%r   (:, jthe), nrho_equ, pfn_tg, RR_tg(:, jthe), nrho_tg)
        call qinterp(pfn_equ, equil_now%coord_sys%position%z   (:, jthe), nrho_equ, pfn_tg, ZZ_tg(:, jthe), nrho_tg)
        call qinterp(pfn_equ, equil_now%coord_sys%bpcell       (:, jthe), nrho_equ, pfn_tg, Bp_tg(:, jthe), nrho_tg)
        call qinterp(pfn_equ, equil_now%coord_sys%position%rmin(:, jthe), nrho_equ, pfn_tg, da_tg(:, jthe), nrho_tg)
    enddo
    da_tg(1, :) = da_tg(2, :) 
    da_tg(nrho_tg, :) = da_tg(nrho_tg-1, :) 
endif

write(6, '(A, 4i4)') 'Call TGLF...', NA1, nrho_tg, sat_rule, tglf_ns_in

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
allocate(ifluxspectrum(tglf_nky_in))
allocate(pfluxspectrum(tglf_nky_in))

radial_loop: do jr=1, nrho_tg

!thermal impurities

    tglf_zs_in(3) = max(1., z_tg(2, jr))
    tglf_zs_in(4) = z_tg(3, jr)
    tglf_zs_in(5) = z_tg(4, jr)

    if (tglf_zs_in(5) >= 1. .and. tglf_ns_in == 3) then
        tglf_zs_in(4) = tglf_zs_in(5)
        tglf_mass_in(4) = tglf_mass_in(5)
        ni_tg(3, :) = ni_tg(4, :)
        ti_tg(3, :) = ti_tg(4, :)
    endif
    if (tglf_zs_in(4) >= 1. .and. tglf_ns_in == 2) then
        tglf_zs_in(3) = tglf_zs_in(4)
        tglf_mass_in(3) = tglf_mass_in(4)
        ni_tg(2, :) = ni_tg(3, :)
        ti_tg(2, :) = ti_tg(3, :)
    endif

!    tglf_ns_in = 3

    tglf_mass_in(tglf_ns_in + 1:) = 0.
    tglf_zs_in  (tglf_ns_in + 1:) = 0.

! Differentials

    jr_r = jr + 1
    jr_l = jr - 1
    if (jr == 1) then
        jr_l = jr
    else if (jr == nrho_tg) then
        jr_r = jr
    endif
    dstep = 1./dble(jr_r - jr_l)  ! 0.5 in between, 1 at the edges

    drmin  = dstep*(ametr_tg(jr_r) - ametr_tg(jr_l))
    drmaj  = dstep*( rmaj_tg(jr_r) -  rmaj_tg(jr_l))
    drho   = dstep*(  rho_tg(jr_r) -   rho_tg(jr_l))
    delong = dstep*( elon_tg(jr_r) -  elon_tg(jr_l))
    dtrian = dstep*( tria_tg(jr_r) -  tria_tg(jr_l))
    dptot  = dstep*( ptot_tg(jr_r) -  ptot_tg(jr_l)) * 1E3*1E13
    dte = dstep*(te_tg(jr_r) - te_tg(jr_l))
    dne = dstep*(ne_tg(jr_r) - ne_tg(jr_l))
    dq    = dstep*(q_tg(jr_r) - q_tg(jr_l))
    dvper = dstep*(vper_tg(jr_r) - vper_tg(jr_l))
    do jspec=1, tglf_ns_in-1
        dti(jspec) = dstep*(ti_tg(jspec, jr_r) - ti_tg(jspec, jr_l))
        dni(jspec) = dstep*(ni_tg(jspec, jr_r) - ni_tg(jspec, jr_l))
    enddo
    dv_r = dstep* &
        (vpar_tg(jr_r)/(rmaj_tg(jr_r) + ametr_tg(jr_r)) - &
         vpar_tg(jr_l)/(rmaj_tg(jr_l) + ametr_tg(jr_l)))
    dr = drmin/a0_m    ! gradients w.r.t. minor radius even for s-alpha geometry
    drhodr = drho/drmin

! local field averages

    tglf_as_in(1)   = 1. ! ne is ref
    tglf_taus_in(1) = 1. ! Te is ref

! Log derivatives
    tglf_rlns_in(1) = -dne/(dr*ne_tg(jr))
    tglf_rlts_in(1) = -dte/(dr*te_tg(jr))

    do jspec=2, tglf_ns_in
        tglf_as_in(jspec)   = ni_tg(jspec-1, jr)/ne_tg(jr)
        tglf_taus_in(jspec) = ti_tg(jspec-1, jr)/te_tg(jr)
        tglf_rlns_in(jspec) = -dni(jspec-1)/(dr*ni_tg(jspec-1, jr))
        tglf_rlts_in(jspec) = -dti(jspec-1)/(dr*ti_tg(jspec-1, jr))
    enddo
! Restore quasi-neutrality via main ions

    tglf_as_in(2)   = -1./tglf_zs_in(2)*(SUM(tglf_as_in*tglf_zs_in) - tglf_as_in(2)*tglf_zs_in(2))
    tglf_rlns_in(2) = -1./(tglf_as_in(2)*tglf_zs_in(2))*(SUM(tglf_rlns_in*tglf_as_in*tglf_zs_in) - tglf_rlns_in(2)*tglf_as_in(2)*tglf_zs_in(2))

! GYRO conventions

    N0  = 1E13*ne_tg(jr)   ! density scale used by GYRO [1/cm**3]
    T0  = 1E3 *te_tg(jr)   ! temperature scale used by GYRO
    Bunit_T = BTOR*drhodr*rho_tg(jr)/ametr_tg(jr)  ! Miller geometry magnetic field unit [gauss]
    Bunit_gauss = 1.d4*Bunit_T 

! derived units for the plasma

    cs0 = SQRT(k0*T0/m0)          ! thermal velocity unit cm/sec
    cs00 = SQRT(e00*T0/(AMJ*mpp)) ! thermal velocity unit m/sec
    omega0 = e0*Bunit_gauss/(m0*c0)     ! gyrofrequency unit 1/sec
    rhos0 = cs0/omega0            ! gyroradius unit cm
    lnlamda = 24.0 -0.5*LOG(tglf_as_in(1)*N0) + LOG(tglf_taus_in(1)*T0)
    taue = 3.44E5 * (tglf_taus_in(1)*T0)**1.5 / (tglf_as_in(1)*N0*lnlamda)  !  sec

    rmin_tg = ametr_tg(jr)/a0_m
    cexb = ametr_tg(jr)/q_tg(jr)

    tglf_vpar_shear_in(2) = -1E2*rmaj_tg(jr)*dv_r/(dr*cs0)  !From m/s to cm/s for vpar
    tglf_vpar_shear_in(1) = tglf_vpar_shear_in(2)

    tglf_vpar_in(2) = 1E2*vpar_tg(jr)/cs0
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

    tglf_vexb_in  = 1E2*vexb_tg(jr)/cs0
    tglf_betae_in = 8.0*pi*k0*N0*T0/Bunit_gauss**2
    tglf_xnue_in  = 0.75*SQRT(pi)*a0_cm/(taue*cs0)
    tglf_zeff_in  = zef_tg(jr)
    tglf_debye_in = SQRT(k0*T0/(4.0*pi*N0*e0**2))/rhos0

    tglf_rmin_loc_in    = rmin_tg
    tglf_rmaj_loc_in    = rmaj_tg(jr)/a0_m
    tglf_zmaj_loc_in    = 0.
    tglf_drmindx_loc_in = 1.
    tglf_drmajdx_loc_in = drmaj/(dr*a0_m)
    tglf_dzmajdx_loc_in = 0.
    tglf_kappa_loc_in   = elon_tg(jr)
    tglf_s_kappa_loc_in = ametr_tg(jr)*delong/(drmin*elon_tg(jr))
    tglf_delta_loc_in   = tria_tg(jr)
    tglf_s_delta_loc_in = ametr_tg(jr)*dtrian/drmin
    tglf_zeta_loc_in    = 0.
    tglf_s_zeta_loc_in  = 0.
    tglf_q_loc_in       = q_tg(jr)
    tglf_q_prime_loc_in = (q_tg(jr)/rmin_tg)*dq/dr
    tglf_p_prime_loc_in = (k0/Bunit_gauss**2)*(q_tg(jr)/rmin_tg)*dptot/dr

    tglf_q_fourier_in       = tglf_q_loc_in
    tglf_q_prime_fourier_in = tglf_q_prime_loc_in
    tglf_p_prime_fourier_in = tglf_p_prime_loc_in
    tglf_q_ELITE_in       = tglf_q_loc_in
    tglf_q_prime_ELITE_in = tglf_q_prime_loc_in
    tglf_p_prime_ELITE_in = tglf_p_prime_loc_in

    tglf_rmin_sa_in     = rmin_tg
    tglf_rmaj_sa_in     = rmaj_tg(jr)/a0_m
    tglf_q_sa_in        = q_tg(jr)
    tglf_shat_sa_in     = (ametr_tg(jr)/q_tg(jr))*dq/drmin
    tglf_alpha_sa_in    = -(8.0*pi*k0/Bunit_gauss**2)*q_tg(jr)**2 * rmaj_tg(jr)*dptot/drmin
    tglf_xwell_sa_in    = 0.
    tglf_theta0_sa_in   = 0.

! Settings
    if (tglf_dump_flag_in) then
        write(file_dump_local, '(A11, I0)') 'input.tglf_', jr
    endif

    if (geom_flag == 2) then ! Fourier moments
        call scrunch2d(mpol, nthe_elite, RR_tg(jr, :), ZZ_tg(jr, :), xrc, xrs, xzc, xzs)
        tglf_fourier_in(1, 0:mpol-1) = xrc/a0_m
        tglf_fourier_in(2, 0:mpol-1) = xrs/a0_m
        tglf_fourier_in(3, 0:mpol-1) = xzc/a0_m
        tglf_fourier_in(4, 0:mpol-1) = xzs/a0_m
! Radial derivatives       
        call qinterp(theta_equ, da_tg(jr, :), nthe_equ, theta_elite, da_elite, nthe_elite)
        drda = (RR_tg(jr_r, :) - RR_tg(jr_l, :))/da_elite
        dzda = (ZZ_tg(jr_r, :) - ZZ_tg(jr_l, :))/da_elite
        do mom_order=0, mpol-1
            cos_mthe(:) = cos(dble(mom_order)*theta_elite(:))
            sin_mthe(:) = sin(dble(mom_order)*theta_elite(:))
            tglf_fourier_in(5, mom_order) = sum(drda*cos_mthe)  ! <Zsin>
            tglf_fourier_in(6, mom_order) = sum(drda*sin_mthe)  ! <Rsin>
            tglf_fourier_in(7, mom_order) = sum(dzda*cos_mthe)  ! <Zsin>
            tglf_fourier_in(8, mom_order) = sum(dzda*sin_mthe)  ! <Rsin>
        enddo
        tglf_fourier_in(5:8, :) = 2.*tglf_fourier_in(5:8, :)/dble(nthe_elite - 1)
    else if (geom_flag == 3) then ! R, Z contours for ELITE
! Interpolation on ELITE theta-grid
        call qinterp(theta_equ, RR_tg(jr, :), nthe_equ, theta_elite, RR_elite, nthe_elite)
        call qinterp(theta_equ, ZZ_tg(jr, :), nthe_equ, theta_elite, ZZ_elite, nthe_elite)
        call qinterp(theta_equ, Bp_tg(jr, :), nthe_equ, theta_elite, Bp_elite, nthe_elite)
        RR_elite(tglf_n_elite_in+1) = RR_elite(1)
        ZZ_elite(tglf_n_elite_in+1) = ZZ_elite(1)
        Bp_elite(tglf_n_elite_in+1) = Bp_elite(1)
        tglf_R_elite_in  = RR_elite/a0_m
        tglf_Z_elite_in  = ZZ_elite/a0_m
        tglf_Bp_elite_in = Bp_elite/Bunit_T

        if (debug) then
            write(f_fourier, '(A11, I0)') 'four4tglf_', jr
            open(21, FILE=f_fourier)
            do jthe=1, nthe_elite
                write(21, '(8F)') tglf_fourier_in(1:8, jthe)
            enddo
            close(21)

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
    ion_eflux = ion_eflux/(tglf_taus_in(2) * 1e13*ni_tg(1, jr)/N0)
    ion_mflux = SUM(tglf_ion_mflux_out(1: tglf_ns_in-1))
    chii (jr) = ion_eflux          /(1e-4 + abs(tglf_rlts_in(2))) *drho_cs
    chie (jr) = tglf_elec_eflux_out/(1e-4 + abs(tglf_rlts_in(1))) *drho_cs
    mtori(jr) = ion_mflux*drho_nt
    elec_pflux(jr) = tglf_elec_pflux_out/drhodr *drho_cs         ! particle flux
    do jspec=1, tglf_ns_in-1
        ion_pflux(jspec, jr) = tglf_ion_pflux_out(jspec)/drhodr *drho_cs  !ion particle flux
    enddo
    exchi(jr) = tglf_elec_expwd_out * nt_cs                ! Equipartition
    do kyloop=1, tglf_nky_in
        gamma(kyloop) = get_eigenvalue_spectrum_out(1, kyloop, 1)
        omega(kyloop) = get_eigenvalue_spectrum_out(2, kyloop, 1)
        kyspectrum(kyloop) = get_ky_spectrum_out(kyloop)
        efluxspectrum(kyloop) = get_flux_spectrum_out(2, 1, 1, kyloop, 1)
        ifluxspectrum(kyloop) = get_flux_spectrum_out(2, 2, 1, kyloop, 1)
        pfluxspectrum(kyloop) = get_flux_spectrum_out(1, 1, 1, kyloop, 1)
    enddo

    jgamma_max = maxloc(efluxspectrum(1:tglf_nky_in), 1)
    gamma_max(jr) = gamma(jgamma_max)
    omega_max(jr) = omega(jgamma_max)
    kymax(jr) = kyspectrum(jgamma_max)

enddo radial_loop

! Interpolate back to ASTRA radial grid

call qinterp(rho_tg, chii      , nrho_tg, RHO(1:NA1), chii_m(1:NA1)      , NA1)
call qinterp(rho_tg, chie      , nrho_tg, RHO(1:NA1), chie_m(1:NA1)      , NA1)
call qinterp(rho_tg, mtori     , nrho_tg, RHO(1:NA1), mtori_m(1:NA1)     , NA1)
call qinterp(rho_tg, elec_pflux, nrho_tg, RHO(1:NA1), elec_pflux_m(1:NA1), NA1)
call qinterp(rho_tg, exchi     , nrho_tg, RHO(1:NA1), exchi_m(1:NA1)     , NA1)
call qinterp(rho_tg, gamma_max , nrho_tg, RHO(1:NA1), gamma_m(1:NA1)     , NA1)
call qinterp(rho_tg, omega_max , nrho_tg, RHO(1:NA1), omega_m(1:NA1)     , NA1)
do jspec=1, tglf_ns_in-1
    call qinterp(rho_tg, ion_pflux(jspec, 1:nrho_tg), nrho_tg, RHO(1:NA1), ion_pflux_m(jspec, 1:NA1), NA1)
enddo

chii_m (1:2) = chii_m (3)
chie_m (1:2) = chie_m (3)
mtori_m(1:2) = mtori_m(3)
elec_pflux_m(1:2) = elec_pflux_m(3)
exchi_m(1:2) = exchi_m(3)
omega_m(1:2) = omega_m(3)
gamma_m(1:2) = gamma_m(3)

do jrho=1, NA1
    CHI(jrho) = chii_m(jrho)/gradrhosq_exp(jrho) ! \chi_i, m^2/s
    CHE(jrho) = chie_m(jrho)/gradrhosq_exp(jrho) ! \chi_e, m^2/s
    VIN(jrho) = elec_pflux_m(jrho)/a0_m/gradrhosq_exp(jrho) ! D flux
    DPR(jrho) = mtori_m(jrho)
!First impurity only, index 2 of ion species
    if (tglf_ns_in >= 3) then
        DPL(jrho) = ion_pflux_m(2, jrho)/a0_m/gradrhosq_exp(jrho)/(ni_exp(2, jrho)/NE(jrho))  ! 1st imp convection
    endif
    if (tglf_ns_in >= 4) then
        DPH(jrho) = ion_pflux_m(3, jrho)/a0_m/gradrhosq_exp(jrho)/(ni_exp(3, jrho)/NE(jrho))  ! 2nd imp convection
    endif
    XTB(jrho) = exchi_m(jrho) ! turbulent e-i equipartition in MW/m^3
    T0  = 1E3 *TE(jrho)       ! temperature scale used by GYRO
    cs0 = SQRT(k0*T0/m0)      ! thermal velocity unit cm/sec
    GM1(jrho) = gamma_m(jrho)*(cs0/a0_cm)
    OM1(jrho) = omega_m(jrho)*(cs0/a0_cm)
enddo

call DATE_AND_TIME(TIME=time_loc)
write(*, *) time_loc, ' END tglf_serial.f90'

return
END subroutine tglf_serial
