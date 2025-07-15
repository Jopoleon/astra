program tglf_chunk

use mpi

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
    tglf_elec_pflux_out, tglf_ion_pflux_out, tglf_elec_expwd_out
use tglf_pkg, only: get_eigenvalue_spectrum_out, get_ky_spectrum_out, &
     get_flux_spectrum_out
  
implicit none

integer, parameter :: nky=19, nspec_max=5, n_dims=8
double precision, parameter :: &
   k0   = 1.6022E-12, &       ! erg/ev
   e0   = 4.8032E-10, &       ! elementary charge (statcoulombs)
   e00  = 1.6020e-19, &       ! elementary charge (C)
   c0   = 2.9979E+10, &       ! speed of light (cm/sec)
   mp   = 1.6726E-24, &       ! proton mass (g)
   mpp  = 1.6726E-27, &       ! proton mass (kg)
   pi   = 3.141592653589793

integer :: ierr, parent, rank, status(MPI_STATUS_SIZE), i1, i2, j, k
integer :: chunk, nprocs, n_inputs, n_outputs, n_scalars, dims(n_dims)
integer :: sat_rule, jr, jgamma_max, jspec, kyloop
double precision :: Bunit_gauss, Bunit_T, cs0, cs00, rhos0, omega0, rhostar2, lnlamda, taue, cexb
double precision :: a0_cm, a0_m, T0, N0, m0, rmin, drho_cs, drho_nt, nt_cs
double precision :: AMJ, BTOR
double precision :: ion_eflux, ion_mflux
double precision, allocatable, dimension(:) :: inputs, output, scalars
double precision, allocatable, dimension(:) :: mtori, chie, chii, exchi, elec_pflux, rho, &
    gamma_max, omega_max, kymax, te, ne, vpar, vper, vexb, &
    ametr, elon, tria, rmaj, ptot, q, zef, pfn, &
    drmin, drmaj, drho, delong, dtrian, dr, dne, dte, dq, dptot, dvpar, dvper, dv_r, drhodr
double precision, allocatable, dimension(:, :) :: dti, dni, ni, ti, zimp, ion_pflux
double precision, dimension(nky) :: gamma, omega, kyspectrum, efluxspectrum, &
    ifluxspectrum, pfluxspectrum

call MPI_Init(ierr)
call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
call MPI_Comm_get_parent(parent, ierr)

if (parent == MPI_COMM_NULL) then
    if (rank == 0) then
        print *, "No parent communicator!"
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    endif
endif

! Receive dimensions from parent
dims = 0
call MPI_Recv(dims, n_dims, MPI_INTEGER, 0, 100+rank, parent, status, ierr)
n_scalars = dims(2)
n_inputs  = dims(3)
n_outputs = dims(4)
i1 = dims(7)
i2 = dims(8)
chunk = i2 + 1 - i1
allocate(scalars(n_scalars))
allocate(inputs(chunk*n_inputs))
scalars = 0.d0
inputs = 0.d0

! Receive TGLF input scalars and profiles from parent
call MPI_Recv(scalars,       n_scalars, MPI_DOUBLE_PRECISION, 0, 101+rank, parent, status, ierr)
call MPI_Recv(inputs, chunk * n_inputs, MPI_DOUBLE_PRECISION, 0, 102+rank, parent, status, ierr)
call MPI_Barrier(parent, ierr)

allocate(output(chunk*n_outputs))
allocate( mtori(chunk), chie(chunk), chii(chunk), exchi(chunk), elec_pflux(chunk), &
    rho(chunk), gamma_max(chunk), omega_max(chunk), kymax(chunk), &
    te(chunk), ne(chunk), vpar(chunk), vper(chunk), vexb(chunk), &
    ametr(chunk), elon(chunk), tria(chunk), rmaj(chunk), &
    ptot(chunk), q(chunk), zef(chunk), pfn(chunk) )
allocate( drmin(chunk), drmaj(chunk), drho(chunk), delong(chunk), dtrian(chunk), &
    dr(chunk), dne(chunk), dte(chunk), dq(chunk), dptot(chunk), dvpar(chunk), &
    dvper(chunk), dv_r(chunk), drhodr(chunk) )
allocate( dti(nspec_max-1, chunk), dni(nspec_max-1, chunk), ni(nspec_max-1, chunk), &
    ti(nspec_max-1, chunk), zimp(nspec_max-2, chunk), ion_pflux(nspec_max-1, chunk) )

tglf_ns_in = dims(5)

! Initialise to zero for non-calculated species

tglf_zs_in   = 0.
tglf_as_in   = 0.
tglf_mass_in = 0.
tglf_taus_in = 0.
tglf_rlns_in = 0.
tglf_rlts_in = 0.
tglf_vpar_in = 0.
tglf_vpar_shear_in = 0.

! Scalars
BTOR = scalars(2)
a0_m = scalars(3)
AMJ  = scalars(5)
tglf_mass_in(1: 5) = scalars(4:  8)/AMJ
tglf_zs_in(1: 5)   = scalars(9: 13)
m0 = AMJ*mp          ! Ref. mass = D ion mass [g]
a0_cm = 1.d2*a0_m    ! length scale used by GYRO, m -> cm

! Profiles
k = 0
rho    = inputs(k+1:k+chunk);   k = k + chunk
ametr  = inputs(k+1:k+chunk);   k = k + chunk
rmaj   = inputs(k+1:k+chunk);   k = k + chunk
elon   = inputs(k+1:k+chunk);   k = k + chunk
tria   = inputs(k+1:k+chunk);   k = k + chunk
q      = inputs(k+1:k+chunk);   k = k + chunk
pfn    = inputs(k+1:k+chunk);   k = k + chunk
ptot   = inputs(k+1:k+chunk);   k = k + chunk
ne     = inputs(k+1:k+chunk);   k = k + chunk
te     = inputs(k+1:k+chunk);   k = k + chunk
zef    = inputs(k+1:k+chunk);   k = k + chunk
vpar   = inputs(k+1:k+chunk);   k = k + chunk
vper   = inputs(k+1:k+chunk);   k = k + chunk
vexb   = inputs(k+1:k+chunk);   k = k + chunk
drmin  = inputs(k+1:k+chunk);   k = k + chunk
drmaj  = inputs(k+1:k+chunk);   k = k + chunk
drho   = inputs(k+1:k+chunk);   k = k + chunk
delong = inputs(k+1:k+chunk);   k = k + chunk
dtrian = inputs(k+1:k+chunk);   k = k + chunk
dptot  = inputs(k+1:k+chunk);   k = k + chunk
dte    = inputs(k+1:k+chunk);   k = k + chunk
dne    = inputs(k+1:k+chunk);   k = k + chunk
dq     = inputs(k+1:k+chunk);   k = k + chunk
dvper  = inputs(k+1:k+chunk);   k = k + chunk
dv_r   = inputs(k+1:k+chunk);   k = k + chunk
dr     = inputs(k+1:k+chunk);   k = k + chunk
drhodr = inputs(k+1:k+chunk);   k = k + chunk
do j=1, 4
    ti(j, :)  = inputs(k+1:k+chunk);    k = k + chunk
    ni(j, :)  = inputs(k+1:k+chunk);    k = k + chunk
    dti(j, :) = inputs(k+1:k+chunk);    k = k + chunk
    dni(j, :) = inputs(k+1:k+chunk);    k = k + chunk
enddo
do j=1, 3
    zimp(j, :) = inputs(k+1:k+chunk); k = k + chunk
enddo

! TGLF settings
sat_rule = 2
SELECT CASE(sat_rule)
CASE(0)
    tglf_nmodes_in = 2
    tglf_xnu_model_in = 2
    tglf_wdia_trapped_in = 0.
    tglf_alpha_zf_in  = 0.
CASE(1)
    tglf_nmodes_in = tglf_ns_in + 2
    tglf_xnu_model_in = 2
    tglf_wdia_trapped_in = 0.
    tglf_alpha_zf_in  = 1.
CASE(2)
    tglf_nmodes_in = tglf_ns_in + 2
    tglf_xnu_model_in = 3
    tglf_wdia_trapped_in = 1.
    tglf_alpha_zf_in  = 1.
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
tglf_nbasis_max_in = 6 ! email Angioni Aug 1st 2023, 4 old default
tglf_nbasis_min_in = 2
tglf_nxgrid_in     = 24 ! (24 email Gary, for ELITE), 16 old default
tglf_nky_in        = nky
tglf_units_in   = 'CGYRO'
tglf_path_in = 'tglf/'
! Want fluxes from TGLF
tglf_use_transport_model_in = .true.
tglf_use_ave_ion_grid_in    = .true. ! Email Angioni Aug 1st 2023

tglf_sign_Bt_in = 1
tglf_sign_It_in = 1

tglf_ky_in        = 0.3
tglf_width_in     = 1.65
tglf_width_min_in = 0.3
tglf_nwidth_in    = 21

tglf_geometry_flag_in = 1
tglf_dump_flag_in     = .False.   ! Dumps input file
tglf_test_flag_in     = 0
tglf_nn_max_error_in  = 0
tglf_write_wavefunction_flag_in = 0 ! Writes eigenfunction

tglf_theta_trapped_in  = 0.7
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
tglf_xnu_factor_in       = 1.
tglf_debye_factor_in     = 1.
tglf_etg_factor_in       = 1.25
tglf_sat_rule_in         = sat_rule
tglf_kygrid_model_in     = 4
tglf_vpar_model_in       = 0
tglf_vpar_shear_model_in = 1

tglf_b_model_sa_in  = 1
tglf_ft_model_sa_in = 1

radial_loop: do jr=1, chunk

!thermal impurities

    tglf_zs_in(3) = max(1., zimp(1, jr))
    tglf_zs_in(4) = zimp(2, jr)
    tglf_zs_in(5) = zimp(3, jr)

    if (tglf_zs_in(5) >= 1. .and. tglf_ns_in == 3) then
        tglf_zs_in(4) = tglf_zs_in(5)
        tglf_mass_in(4) = tglf_mass_in(5)
        ni(3, :) = ni(4, :)
        ti(3, :) = ti(4, :)
    endif
    if (tglf_zs_in(4) >= 1. .and. tglf_ns_in == 2) then
        tglf_zs_in(3) = tglf_zs_in(4)
        tglf_mass_in(3) = tglf_mass_in(4)
        ni(2, :) = ni(3, :)
        ti(2, :) = ti(3, :)
    endif

!    tglf_ns_in = 3

    tglf_mass_in(tglf_ns_in + 1:) = 0.
    tglf_zs_in  (tglf_ns_in + 1:) = 0.

    tglf_as_in(1)   = 1. ! ne is ref
    tglf_taus_in(1) = 1. ! Te is ref

! Log derivatives
    tglf_rlns_in(1) = -dne(jr)/(dr(jr)*ne(jr))
    tglf_rlts_in(1) = -dte(jr)/(dr(jr)*te(jr))

    do jspec=2, tglf_ns_in
        tglf_as_in(jspec)   = ni(jspec-1, jr)/ne(jr)
        tglf_taus_in(jspec) = ti(jspec-1, jr)/te(jr)
        tglf_rlns_in(jspec) = -dni(jspec-1, jr)/(dr(jr)*ni(jspec-1, jr))
        tglf_rlts_in(jspec) = -dti(jspec-1, jr)/(dr(jr)*ti(jspec-1, jr))
    enddo
! Restore quasi-neutrality via main ions

    tglf_as_in(2)   = -1./tglf_zs_in(2)*(SUM(tglf_as_in*tglf_zs_in) - tglf_as_in(2)*tglf_zs_in(2))
    tglf_rlns_in(2) = -1./(tglf_as_in(2)*tglf_zs_in(2))*(SUM(tglf_rlns_in*tglf_as_in*tglf_zs_in) - tglf_rlns_in(2)*tglf_as_in(2)*tglf_zs_in(2))

! GYRO conventions

    N0  = 1E13*ne(jr)   ! density scale used by GYRO [1/cm**3]
    T0  = 1E3 *te(jr)   ! temperature scale used by GYRO
    Bunit_T = BTOR*drhodr(jr)*rho(jr)/ametr(jr)  ! Miller geometry magnetic field unit [gauss]
    Bunit_gauss = 1.d4*Bunit_T 

! derived units for the plasma

    cs0 = SQRT(k0*T0/m0)          ! thermal velocity unit cm/sec
    cs00 = SQRT(e00*T0/(AMJ*mpp)) ! thermal velocity unit m/sec
    omega0 = e0*Bunit_gauss/(m0*c0)     ! gyrofrequency unit 1/sec
    rhos0 = cs0/omega0            ! gyroradius unit cm
    lnlamda = 24.0 -0.5*LOG(tglf_as_in(1)*N0) + LOG(tglf_taus_in(1)*T0)
    taue = 3.44E5 * (tglf_taus_in(1)*T0)**1.5 / (tglf_as_in(1)*N0*lnlamda)  !  sec

    rmin = ametr(jr)/a0_m
    cexb = ametr(jr)/q(jr)

    tglf_vpar_shear_in(2) = -1E2*rmaj(jr)*dv_r(jr)/(dr(jr)*cs0)  !From m/s to cm/s for vpar
    tglf_vpar_shear_in(1) = tglf_vpar_shear_in(2)

    tglf_vpar_in(2) = 1E2*vpar(jr)/cs0
    tglf_vpar_in(1) = tglf_vpar_in(2)

    if (tglf_ns_in >= 3) then
        tglf_vpar_shear_in(3: tglf_ns_in) = tglf_vpar_shear_in(2)
        tglf_vpar_in      (3: tglf_ns_in) = tglf_vpar_in(2)
    endif

! local magnetic geometry

    rhostar2 = (rhos0/a0_cm)**2
    drho_cs = drhodr(jr)**2*a0_m*rhostar2*cs00
    drho_nt = drhodr(jr)*a0_m*rhostar2*N0*e00*T0/(1.e13*mpp)
    nt_cs = 0.001602*N0/1.e13*T0/1.e3*cs00/a0_m*rhostar2

! Share variables with tglf_run via module tglf_interface

    tglf_vexb_shear_in = -1E2*cexb*dvper(jr)/(dr(jr)*cs0) ! Waltz-Miller definition     !From m/s to cm/s for vexb  -dVexb/dr

! Initialise

    tglf_vexb_in  = 1E2*vexb(jr)/cs0
    tglf_betae_in = 8.0*pi*k0*N0*T0/Bunit_gauss**2
    tglf_xnue_in  = 0.75*SQRT(pi)*a0_cm/(taue*cs0)
    tglf_zeff_in  = zef(jr)
    tglf_debye_in = SQRT(k0*T0/(4.0*pi*N0*e0**2))/rhos0

    tglf_rmin_loc_in    = rmin
    tglf_rmaj_loc_in    = rmaj(jr)/a0_m
    tglf_zmaj_loc_in    = 0.
    tglf_drmindx_loc_in = 1.
    tglf_drmajdx_loc_in = drmaj(jr)/(dr(jr)*a0_m)
    tglf_dzmajdx_loc_in = 0.
    tglf_kappa_loc_in   = elon(jr)
    tglf_s_kappa_loc_in = ametr(jr)*delong(jr)/(drmin(jr)*elon(jr))
    tglf_delta_loc_in   = tria(jr)
    tglf_s_delta_loc_in = ametr(jr)*dtrian(jr)/drmin(jr)
    tglf_zeta_loc_in    = 0.
    tglf_s_zeta_loc_in  = 0.
    tglf_q_loc_in       = q(jr)
    tglf_q_prime_loc_in = (q(jr)/rmin)*dq(jr)/dr(jr)
    tglf_p_prime_loc_in = (k0/Bunit_gauss**2)*(q(jr)/rmin)*dptot(jr)/dr(jr)

    tglf_rmin_sa_in     = rmin
    tglf_rmaj_sa_in     = rmaj(jr)/a0_m
    tglf_q_sa_in        = q(jr)
    tglf_shat_sa_in     = (ametr(jr)/q(jr))*dq(jr)/drmin(jr)
    tglf_alpha_sa_in    = -(8.0*pi*k0/Bunit_gauss**2)*q(jr)**2 * rmaj(jr)*dptot(jr)/drmin(jr)
    tglf_xwell_sa_in    = 0.
    tglf_theta0_sa_in   = 0.

! Settings
    if (tglf_dump_flag_in) then
        write(file_dump_local, '(A11, I0)') 'input.tglf_', jr + i1 
    endif
! -----------------------
    call tglf_run
! -----------------------

! Transport coefficients

    ion_eflux = SUM(tglf_ion_eflux_out(1: tglf_ns_in-1))
    ion_eflux = ion_eflux/(tglf_taus_in(2) * 1e13*ni(1, jr)/N0)
    ion_mflux = SUM(tglf_ion_mflux_out(1: tglf_ns_in-1))
    chii (jr) = ion_eflux          /(1e-4 + abs(tglf_rlts_in(2))) *drho_cs
    chie (jr) = tglf_elec_eflux_out/(1e-4 + abs(tglf_rlts_in(1))) *drho_cs
    mtori(jr) = ion_mflux*drho_nt
    elec_pflux(jr) = tglf_elec_pflux_out/drhodr(jr) *drho_cs         ! particle flux
    do jspec=1, tglf_ns_in-1
        ion_pflux(jspec, jr) = tglf_ion_pflux_out(jspec)/drhodr(jr) *drho_cs  !ion particle flux
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

! Simulated TGLF computation:
output = 0.d0
output(        1:   chunk) = chii
output(  chunk+1: 2*chunk) = chie
output(2*chunk+1: 3*chunk) = mtori
output(3*chunk+1: 4*chunk) = elec_pflux
output(4*chunk+1: 5*chunk) = exchi
output(5*chunk+1: 6*chunk) = gamma_max
output(6*chunk+1: 7*chunk) = omega_max
do jspec=1, tglf_ns_in-1
    output((6+jspec)*chunk+1: (7+jspec)*chunk) = ion_pflux(jspec, :)
enddo

call MPI_Send(output, chunk * n_outputs, MPI_DOUBLE_PRECISION, 0, 103+rank, parent, ierr)
call MPI_Barrier(parent, ierr)
call MPI_Finalize(ierr)

end program tglf_chunk
