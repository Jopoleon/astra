!----------------------------------------------------------------------|
program main

use parameter_inc, only: n_ql

implicit none

integer :: iargc, mampid, mamkey, eignr
character(len=132) :: STRING, eigpath, mampath

if (iargc() /= 4) then
    write(6, '(A)') "Error"
    call a_stop
endif
call getarg(0, STRING)
eigpath = STRING(1:len_trim(STRING))//char(0)
call getarg(1, STRING)
mampath = STRING(1:len_trim(STRING))//char(0)
call getarg(2, STRING)
read(STRING, *) mampid
call getarg(3, STRING)
read(STRING, *) mamkey
call getarg(4, STRING)
read(STRING, *) eignr
call sbp2shm(eigpath, mampath, n_ql, mampid, mamkey, eignr)

end program main

!----------------------------------------------------------------------|
subroutine tglf_interf(jr1_in, jr2_in, nrho, NA1N, NA1E, NA1I, &
    BTOR, RTOR, AMJ, ZMJ, AIM1, AIM2, AIM3, &
    NE, TE, NI, NDEUT, NTRIT, NIZ1, NIZ2, TI, ZEF, ZIM1, AMAIN, &
    MU, RHO, AMETR, SHIF, ELON, TRIA, ER, NIBM, G11, &
    VPOL, VRS, VTOR, SHEAR, PBLON, PBPER, PFAST, NIZ3, ZIM2, ZIM3, &
    ZIMPT, NIMPT, AIMPT, &
! output
    mem_out)

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

implicit none

logical, parameter :: debug_elite=.false.
integer, parameter :: nradial=5, unit_elite=11, nthe_elite=400
double precision, parameter :: &
   k0   = 1.6022d-12, &       ! erg/ev
   e0   = 4.8032d-10, &       ! elementary charge (statcoulombs)
   e00  = 1.6020d-19, &       ! elementary charge (C)
   c0   = 2.9979d+10, &       ! speed of light (cm/sec)
   mp   = 1.6726d-24, &       ! proton mass (g)
   mpp  = 1.6726d-27, &       ! proton mass (kg)
   pi   = 3.141592653589793 
double precision, parameter :: c_vpol=1.d0

integer, intent(in) :: nrho, jr1_in, jr2_in, NA1N, NA1E, NA1I

double precision, intent(in) :: BTOR, RTOR, AMJ, AIM1, AIM2, AIM3, ZMJ
double precision, intent(in), dimension(*) :: NE, TE, NI, TI, &
    ZEF, ZIM1, ZIM2, ZIM3, PBLON, PBPER, PFAST, NIZ3, AMAIN, &
    ER, MU, RHO, AMETR, SHIF, ELON, NDEUT, NIZ1, NTRIT, &
    NIZ2, TRIA, NIBM, G11, VPOL, VRS, VTOR, SHEAR, &
    ZIMPT, NIMPT, AIMPT

double precision, intent(out), dimension(nrho, 15) :: mem_out

!----------------------------------------------------------------------

integer :: i_ion, n, ios
integer :: jr_min, jr_max, jrho, j0, j01, j02, jgamma_max, n_radial, jthe, nthe_equ, nrho_as
integer :: j, jradial, jjgrid(nradial), jspec, kyloop
integer :: sat_rule           ! Saturation rule
integer :: geom_flag          ! 1: Miller; 2: Fourier; 3: ELITE
integer :: nmodes_tg          ! number of unstable modes to use in computing fluxes (max=4)
integer :: kygrid_model_tg    ! select version of ky-grid to use 1
integer :: xnu_model_tg       ! select version of trapped-passing 2

double precision :: bmod, bpolz, alpha_zf_in, ion_eflux, ion_mflux
double precision :: drmin, drmaj, drho, dte, dne, dq, dptot, &
    delong, dtrian, dvper, drhodr, dstep, dr, dv_r
double precision :: Bunit, cs0, cs00, rhos0, omega0, rhostar2, lnlamda, taue, cexb
double precision :: a0_m, a0_cm, T0, N0, m0, rmin_tg, drho_cs, drho_nt, nt_cs
double precision :: wdia_trap_tg          ! parameter for trapped fraction model

double precision, dimension(nrho) :: gradrhosq_exp, rmaj_exp, q_exp, &
    rho_m, vexb2, vpar_m, vper_m, mtori_m, &
    chie_m, chii_m, elec_pflux_m, exchi_m, ptot, gamma_m, omega_m
double precision, dimension(nradial) :: mtori, chie, chii, exchi, elec_pflux, rho_tg, &
    gamma_max, omega_max, kymax
double precision, dimension(nsm-1, nrho) :: ion_pflux_m
double precision, dimension(nsm-1, nradial) :: ion_pflux
double precision, allocatable, dimension(:) :: gamma, omega, kyspectrum, efluxspectrum, ifluxspectrum, pfluxspectrum
double precision, dimension(nsm-1) :: dti, dni
double precision, dimension(nsm-1, nrho) :: ni_m, ti_m
double precision :: dtheta_elite, xstep
double precision, dimension(nthe_elite) :: theta_elite, RR_elite, ZZ_elite, Bp_elite
double precision, allocatable, dimension(:) :: theta_equ
double precision, allocatable, dimension(:, :) :: RR_as, ZZ_as, Bp_as
character(len=120) :: f_elite

!-----------------
! Radial subdomain
!-----------------

jr_min = max(1, jr1_in)
jr_max = min(nrho, jr2_in)
if (jr_min > jr_max) then
    write(*, *) 'Error! jr_min > jr_max', jr_max
    return
endif

xstep = dble(jr_max - jr_min)/(nradial - 1.)
if (xstep <= 1.) then
    do jradial=1, nradial
        jjgrid(jradial) = jr_min + jradial - 1
        if (jjgrid(jradial) == jr_max) EXIT
    enddo
    n_radial = jradial
else
    n_radial = nradial
    do jradial = 1, n_radial-1
        jjgrid(jradial) = jr_min + nint((jradial-1)*xstep)
    enddo
    jjgrid(n_radial) = jr_max
endif

! Electrons and main ions
tglf_zs_in(1) = -1.
tglf_zs_in(2) = ZMJ

tglf_mass_in(1) = 5.4447e-4/AMJ
tglf_mass_in(2) = AMJ/AMJ  ! AMJ is reference mass
tglf_mass_in(3) = AIM1/AMJ
tglf_mass_in(4) = AIM2/AMJ
tglf_mass_in(5) = AIM3/AMJ

do jrho=1, nrho
    rho_m(jrho) = RHO(jrho)
    ti_m(1:4, jrho) = TI(jrho)
    if (NDEUT(jrho) >= 0.01*NE(jrho)) then
        ni_m(1, jrho) = NDEUT(jrho)
    else ! likely: NDEUT not defined in equ file, hence zero
        ni_m(1, jrho) = NI(jrho)
    endif
    ni_m(2, jrho) = max(1.e-9, NIZ1(jrho))
    ni_m(3, jrho) = max(1.e-9, NIZ2(jrho))
    ni_m(4, jrho) = max(1.e-9, NIZ3(jrho))
    rmaj_exp(jrho) = RTOR + SHIF(jrho)
    q_exp(jrho)    = 1./MU(jrho)
    ptot(jrho) = NE(jrho)*TE(jrho) + ni_m(1, jrho)*ti_m(1, jrho) + ni_m(2, jrho)*ti_m(2, jrho) + pfast(jrho) + 0.5*(pblon(jrho) + pbper(jrho))
    bpolz = BTOR*AMETR(jrho)*MU(jrho)/RTOR
    bmod = sqrt(BTOR**2 + bpolz**2)
    gradrhosq_exp(jrho) = G11(jrho)/VRS(jrho)

    vper_m(jrho) = ER(jrho)/(RTOR*bpolz) ! vexb in m/s --> Omega_E    
!    vpar_m(jrho) = ER(jrho)/(RTOR*bpolz)*(RTOR+SHIF(jrho)+AMETR(jrho))  !--> R*Omega_E , no neoclassical terms
    vpar_m(jrho) = VTOR(jrho) * BTOR/bmod + c_vpol* VPOL(jrho) * bpolz/bmod
    vexb2(jrho)  = -ER(jrho)/bmod ! vexb in m/s (vperp = vexb since the diamagnetic velocity is the curvature drift actually

enddo

m0 = AMJ*mp             ! Ref. mass = D ion mass [g]
a0_m = AMETR(nrho)    ! length scale used by GYRO from AMETR (meters) to cm
a0_cm = 1E2*a0_m    ! length scale used by GYRO from AMETR (meters) to cm

elec_pflux_m = 0.
ion_pflux_m  = 0.
chie_m  = 0.
chii_m  = 0.
mtori_m = 0.
exchi_m = 0.

! Number of species

tglf_ns_in = nsm
! These will be reset locally in the radial loop
tglf_zs_in(3) = MAXVAL(ZIM1(1:nrho))
tglf_zs_in(4) = MAXVAL(ZIM2(1:nrho))
tglf_zs_in(5) = MAXVAL(ZIM3(1:nrho))

if (tglf_zs_in(5) >= 1.) tglf_ns_in = 5
if (tglf_zs_in(5) <  1.) tglf_ns_in = 4
if (tglf_zs_in(4) <  1.) tglf_ns_in = 3
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

write(6, '(A, 9i4)') 'Call TGLF...', jjgrid(1: nradial), nrho, sat_rule, geom_flag, tglf_ns_in

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
! Fourier and ELITE coefficients

if (geom_flag == 3) then ! Read ELITE file
    tglf_n_elite_in  = nthe_elite - 1
    dtheta_elite = 2.*pi/dble(nthe_elite-1)
    theta_elite = (/ ((jthe - 1.)*dtheta_elite, jthe=1, nthe_elite) /)
    f_elite = '../tglf/tglf4elite.dat'
    open(unit_elite, file=TRIM(f_elite), iostat=ios, status='OLD')
    if (ios == 0) then
        read(unit_elite, '(i)') nthe_equ
        if (.not. allocated(theta_equ)) allocate(theta_equ(nthe_equ))
        do jthe=1, nthe_equ
            read (unit_elite, '(e14.6)') theta_equ(jthe)
        enddo
        read(unit_elite, '(2i)') nrho_as, nthe_equ
        if (.not. allocated(RR_as)) then
            allocate(RR_as(nrho_as, nthe_equ))
            allocate(ZZ_as(nrho_as, nthe_equ))
            allocate(Bp_as(nrho_as, nthe_equ))
        endif

        do jrho=1, nrho_as
            do jthe=1, nthe_equ
                read(unit_elite, '(3e14.6)') RR_as(jrho, jthe), ZZ_as(jrho, jthe), Bp_as(jrho, jthe)
            enddo
        enddo
    else
        write(*, *) 'Seeting geom_flag=1 only for this time step', ios
        geom_flag = 1
    endif
    close(unit_elite)
endif

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
tglf_nbasis_max_in = 6 ! 4 email Angioni Aug 1st 2023
tglf_nbasis_min_in = 2
tglf_nxgrid_in     = 16
tglf_nky_in        = 19
tglf_units_in   = 'CGYRO'
tglf_path_in = '../tglf/'
! Want fluxes from TGLF
tglf_use_transport_model_in = .true.
tglf_use_ave_ion_grid_in    = .true. ! Email Angioni Aug 1st 2023

tglf_sign_Bt_in = 1
tglf_sign_It_in = 1

tglf_ky_in = 0.3
tglf_width_in      = 1.65
tglf_width_min_in  = 0.3
tglf_nwidth_in     = 21

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

radial_loop: do jradial=1, n_radial

    j0 = jjgrid(jradial)

!thermal impurities

    tglf_zs_in(3) = max(1., ZIM1(j0))
    tglf_zs_in(4) = ZIM2(j0)
    tglf_zs_in(5) = ZIM3(j0)

    if (tglf_zs_in(5) >= 1. .and. tglf_ns_in == 3) then 
        tglf_zs_in(4) = tglf_zs_in(5)
        tglf_mass_in(4) = tglf_mass_in(5)
        ni_m(3, :) = ni_m(4, :)
        ti_m(3, :) = ti_m(4, :)
    endif
    if (tglf_zs_in(4) >= 1. .and. tglf_ns_in == 2) then 
        tglf_zs_in(3) = tglf_zs_in(4)
        tglf_mass_in(3) = tglf_mass_in(4)
        ni_m(2, :) = ni_m(3, :)
        ti_m(2, :) = ti_m(3, :)
    endif

!    tglf_ns_in = 3

    tglf_mass_in(tglf_ns_in + 1:) = 0.
    tglf_zs_in  (tglf_ns_in + 1:) = 0.

! Selected grid for TGLF computation

    rho_tg(jradial) = rho_m(j0)

! Differentials

    j01 = j0 + 1
    j02 = j0 - 1
    if (j0 == 1) then
        j02 = j0
    else if (j0 == nrho) then
        j01 = j0
    endif
    dstep = 1./dble(j01 - j02)

    drmin  = dstep*(AMETR(j01) - AMETR(j02))
    drmaj  = dstep*(rmaj_exp(j01) - rmaj_exp(j02))
    drho   = dstep*(rho(j01) - rho(j02))
    delong = dstep*(ELON(j01) - ELON(j02))
    dtrian = dstep*(TRIA(j01) - TRIA(j02))
    dptot  = dstep*(ptot(j01) - ptot(j02)) * 1E3*1E13
    if (j0 == NA1E .and. NA1E /= nrho) then
        dte = 0.5*dstep*(TE(j0) - TE(j02)) ! Left derivative
    else
        dte = dstep*(TE(j01) - TE(j02))
    endif
    if (j0 == NA1N .and. NA1N /= nrho) then
        dne = 0.5*dstep*(NE(j0) - NE(j02)) ! Left derivative
    else
        dne = dstep*(NE(j01) - NE(j02))
    endif
    dq    = dstep*(q_exp(j01) - q_exp(j02))
    dvper = dstep*(vper_m(j01) - vper_m(j02))
    do jspec=1, tglf_ns_in-1
        if (j0 == NA1I .and. NA1I /= nrho) then
            dti(jspec) = 0.5*dstep*(ti_m(jspec, j0) - ti_m(jspec, j02))
            dni(jspec) = 0.5*dstep*(ni_m(jspec, j0) - ni_m(jspec, j02))
        else
            dti(jspec) = dstep*(ti_m(jspec, j01) - ti_m(jspec, j02))
            dni(jspec) = dstep*(ni_m(jspec, j01) - ni_m(jspec, j02))
        endif
    enddo
    dv_r = dstep* &
        (vpar_m(j01)/(rmaj_exp(j01) + AMETR(j01)) - &
         vpar_m(j02)/(rmaj_exp(j02) + AMETR(j02)))
    dr = drmin/a0_m    ! gradients w.r.t. minor radius even for s-alpha geometry
    drhodr = drho/drmin

! local field averages

    tglf_as_in(1)   = 1. ! ne is ref
    tglf_taus_in(1) = 1. ! Te is ref

! Log derivatives
    tglf_rlns_in(1) = -dne/(dr*NE(j0))
    tglf_rlts_in(1) = -dte/(dr*TE(j0))

    do jspec=2, tglf_ns_in
        tglf_as_in(jspec)   = ni_m(jspec-1, j0)/NE(j0)
        tglf_taus_in(jspec) = ti_m(jspec-1, j0)/TE(j0)
        tglf_rlns_in(jspec) = -dni(jspec-1)/(dr*ni_m(jspec-1, j0))
        tglf_rlts_in(jspec) = -dti(jspec-1)/(dr*ti_m(jspec-1, j0))
    enddo
! Restore quasi-neutrality via main ions

    tglf_as_in(2)   = -1./tglf_zs_in(2)*(SUM(tglf_as_in*tglf_zs_in) - tglf_as_in(2)*tglf_zs_in(2))
    tglf_rlns_in(2) = -1./(tglf_as_in(2)*tglf_zs_in(2))*(SUM(tglf_rlns_in*tglf_as_in*tglf_zs_in) - tglf_rlns_in(2)*tglf_as_in(2)*tglf_zs_in(2))

! GYRO conventions

    N0  = 1E13*NE(j0)   ! density scale used by GYRO
    T0  = 1E3 *TE(j0)   ! temperature scale used by GYRO
    Bunit = 1E4*BTOR*drhodr*rho_m(j0)/AMETR(j0)  ! Miller geometry magnetic field unit

! derived units for the plasma

    cs0 = SQRT(k0*T0/m0)          ! thermal velocity unit cm/sec
    cs00 = SQRT(e00*T0/(AMJ*mpp)) ! thermal velocity unit m/sec
    omega0 = e0*Bunit/(m0*c0)     ! gyrofrequency unit 1/sec
    rhos0 = cs0/omega0            ! gyroradius unit cm
    lnlamda = 24.0 -0.5*LOG(tglf_as_in(1)*N0) + LOG(tglf_taus_in(1)*T0)       
    taue = 3.44E5 * (tglf_taus_in(1)*T0)**1.5 / (tglf_as_in(1)*N0*lnlamda)  !  sec

    rmin_tg = AMETR(j0)/a0_m
    cexb = AMETR(j0)/q_exp(j0)

    tglf_vpar_shear_in(2) = -1E2*rmaj_exp(j0)*dv_r/(dr*cs0)  !From m/s to cm/s for vpar
    tglf_vpar_shear_in(1) = tglf_vpar_shear_in(2)

    tglf_vpar_in(2) = 1E2*vpar_m(j0)/cs0
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

    tglf_vexb_in  = 1E2*vexb2(j0)/cs0
    tglf_betae_in = 8.0*pi*k0*N0*T0/Bunit**2
    tglf_xnue_in  = 0.75*SQRT(pi)*a0_cm/(taue*cs0)
    tglf_zeff_in  = ZEF(j0) 
    tglf_debye_in = SQRT(k0*T0/(4.0*pi*N0*e0**2))/rhos0

    tglf_rmin_loc_in    = rmin_tg
    tglf_rmaj_loc_in    = rmaj_exp(j0)/a0_m
    tglf_zmaj_loc_in    = 0.
    tglf_drmindx_loc_in = 1.
    tglf_drmajdx_loc_in = drmaj/(dr*a0_m)
    tglf_dzmajdx_loc_in = 0.
    tglf_kappa_loc_in   = ELON(j0)
    tglf_s_kappa_loc_in = AMETR(j0)*delong/(drmin*ELON(j0))
    tglf_delta_loc_in   = TRIA(j0)
    tglf_s_delta_loc_in = AMETR(j0)*dtrian/drmin
    tglf_zeta_loc_in    = 0.
    tglf_s_zeta_loc_in  = 0.
    tglf_q_loc_in       = q_exp(j0)
    tglf_q_prime_loc_in = (q_exp(j0)/rmin_tg)*dq/dr
    tglf_p_prime_loc_in = (k0/Bunit**2)*(q_exp(j0)/rmin_tg)*dptot/dr

    tglf_q_ELITE_in       = tglf_q_loc_in
    tglf_q_prime_ELITE_in = tglf_q_prime_loc_in
    tglf_p_prime_ELITE_in = tglf_p_prime_loc_in

    tglf_rmin_sa_in     = rmin_tg
    tglf_rmaj_sa_in     = rmaj_exp(j0)/a0_m
    tglf_q_sa_in        = q_exp(j0)
    tglf_shat_sa_in     = (AMETR(j0)/q_exp(j0))*dq/drmin
    tglf_alpha_sa_in    = -(8.0*pi*k0/Bunit**2)*q_exp(j0)**2 * rmaj_exp(j0)*dptot/drmin
    tglf_xwell_sa_in    = 0.
    tglf_theta0_sa_in   = 0.

! Settings
    if (tglf_dump_flag_in) then
        write(file_dump_local, '(A11, I0)') 'input.tglf_', j0
    endif

    if (geom_flag == 3) then ! R, Z contours for ELITE
! Interpolation on ELITE theta-grid
        call qinterp(theta_equ, RR_as(j0, :), nthe_equ, theta_elite, RR_elite, nthe_elite)
        call qinterp(theta_equ, ZZ_as(j0, :), nthe_equ, theta_elite, ZZ_elite, nthe_elite)
        call qinterp(theta_equ, Bp_as(j0, :), nthe_equ, theta_elite, Bp_elite, nthe_elite)
        RR_elite(tglf_n_elite_in+1) = RR_elite(1)
        ZZ_elite(tglf_n_elite_in+1) = ZZ_elite(1)
        Bp_elite(tglf_n_elite_in+1) = Bp_elite(1)
        tglf_R_elite_in(1:tglf_n_elite_in+1)  = RR_elite/a0_m
        tglf_Z_elite_in(1:tglf_n_elite_in+1)  = ZZ_elite/a0_m
        tglf_Bp_elite_in(1:tglf_n_elite_in+1) = 1.d4*Bp_elite/Bunit

        if (debug_elite) then
            write(f_elite, '(A, I0)') '../tglf/elite4tglf_', j0
            open(31, FILE=f_elite)
            do jthe=1, nthe_elite
                write(31, '(3F)') tglf_R_elite_in(jthe), tglf_Z_elite_in(jthe), Bp_elite(jthe)
            enddo
            close(31)
        endif
    endif

    call tglf_run

! Transport coefficients

    ion_eflux = SUM(tglf_ion_eflux_out(1: tglf_ns_in-1))
    ion_eflux = ion_eflux/(tglf_taus_in(2) * 1e13*NI(j0)/N0)
    ion_mflux = SUM(tglf_ion_mflux_out(1: tglf_ns_in-1))
    chii (jradial) = ion_eflux          /(1e-4 + abs(tglf_rlts_in(2))) *drho_cs
    chie (jradial) = tglf_elec_eflux_out/(1e-4 + abs(tglf_rlts_in(1))) *drho_cs
    mtori(jradial) = ion_mflux*drho_nt
    elec_pflux(jradial) = tglf_elec_pflux_out/drhodr *drho_cs         ! particle flux
    do jspec=1, tglf_ns_in-1
        ion_pflux(jspec, jradial) = tglf_ion_pflux_out(jspec)/drhodr *drho_cs  !ion particle flux
    enddo
    exchi(jradial) = tglf_elec_expwd_out * nt_cs                ! Equipartition
    do kyloop=1, tglf_nky_in
        gamma(kyloop) = get_eigenvalue_spectrum_out(1, kyloop, 1)
        omega(kyloop) = get_eigenvalue_spectrum_out(2, kyloop, 1)
        kyspectrum(kyloop) = get_ky_spectrum_out(kyloop)
        efluxspectrum(kyloop) = get_flux_spectrum_out(2, 1, 1, kyloop, 1)
        ifluxspectrum(kyloop) = get_flux_spectrum_out(2, 2, 1, kyloop, 1)
        pfluxspectrum(kyloop) = get_flux_spectrum_out(1, 1, 1, kyloop, 1)
    enddo

    jgamma_max = maxloc(efluxspectrum(1:tglf_nky_in), 1)
    gamma_max(jradial) = gamma(jgamma_max)
    omega_max(jradial) = omega(jgamma_max)
    kymax(jradial) = kyspectrum(jgamma_max)

enddo radial_loop

call qinterp(rho_tg, chii      , n_radial, rho_m(jr_min:jr_max), chii_m(jr_min:jr_max)      , jr_max-jr_min+1)
call qinterp(rho_tg, chie      , n_radial, rho_m(jr_min:jr_max), chie_m(jr_min:jr_max)      , jr_max-jr_min+1)
call qinterp(rho_tg, mtori     , n_radial, rho_m(jr_min:jr_max), mtori_m(jr_min:jr_max)     , jr_max-jr_min+1)
call qinterp(rho_tg, elec_pflux, n_radial, rho_m(jr_min:jr_max), elec_pflux_m(jr_min:jr_max), jr_max-jr_min+1)
call qinterp(rho_tg, exchi     , n_radial, rho_m(jr_min:jr_max), exchi_m(jr_min:jr_max)     , jr_max-jr_min+1)
call qinterp(rho_tg, gamma_max , n_radial, rho_m(jr_min:jr_max), gamma_m(jr_min:jr_max)     , jr_max-jr_min+1)
call qinterp(rho_tg, omega_max , n_radial, rho_m(jr_min:jr_max), omega_m(jr_min:jr_max)     , jr_max-jr_min+1)
do jspec=1, tglf_ns_in-1
    call qinterp(rho_tg, ion_pflux(jspec, 1:n_radial), n_radial, rho_m(jr_min:jr_max), ion_pflux_m(jspec, jr_min:jr_max), jr_max-jr_min+1)
enddo

chii_m (1:2) = chii_m (3)
chie_m (1:2) = chie_m (3)
mtori_m(1:2) = mtori_m(3)
elec_pflux_m(1:2) = elec_pflux_m(3)
exchi_m(1:2) = exchi_m(3)
omega_m(1:2) = omega_m(3)
gamma_m(1:2) = gamma_m(3)

mem_out = 0.d0

do j=jr_min, jr_max
    mem_out(j, 1) = chii_m(j)/gradrhosq_exp(j) ! \chi_i, m^2/s
    mem_out(j, 2) = chie_m(j)/gradrhosq_exp(j) ! \chi_e, m^2/s
    mem_out(j, 4) = elec_pflux_m(j)/AMETR(nrho)/gradrhosq_exp(j) ! D flux
    mem_out(j, 7) = mtori_m(j)
!First impurity only, index 2 of ion species
    if (tglf_ns_in >= 3) then
        mem_out(j, 6) = ion_pflux_m(2, j)/AMETR(nrho)/gradrhosq_exp(j)/(ni_m(2, j)/NE(j))  ! 1st imp convection
    endif
    if (tglf_ns_in >= 4) then
        mem_out(j, 5) = ion_pflux_m(3, j)/AMETR(nrho)/gradrhosq_exp(j)/(ni_m(3, j)/NE(j))  ! 2nd imp convection
    endif
    mem_out(j,  8) = exchi_m(j)  ! turbulent e-i equipartition in MW/m^3
    mem_out(j, 11) = gamma_m(j)*(cs0/a0_cm)
    mem_out(j, 13) = omega_m(j)*(cs0/a0_cm)
enddo

mem_out(1, 6) = 0.d0 !ensure NIZ1 convection equal to zero on axis
mem_out(1, 5) = 0.d0 !ensure NIZ2 convection equal to zero on axis (already 0 otherwise)

deallocate(gamma)
deallocate(omega)
deallocate(kyspectrum)
deallocate(efluxspectrum)
deallocate(ifluxspectrum)
deallocate(pfluxspectrum)

return
END subroutine tglf_interf
