!----------------------------------------------------------------------|
program main

implicit none

integer :: iargc, mampid, mamkey, eignr
character*132 :: STRING, eigpath, mampath

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
call sbp2shm(eigpath, mampath, mampid, mamkey, eignr)

end program main

!----------------------------------------------------------------------|
subroutine tglf_interf(jr1_in, jr2_in, nrho, NA1N, NA1E, NA1I, &
    BTOR, RTOR, AMJ, ZMJ, AIM1, AIM2, AIM3, &
    NE, TE, NI, NDEUT, NTRIT, NIZ1, NIZ2, TI, ZEF, ZIM1, AMAIN, &
    MU, RHO, AMETR, SHIF, ELON, TRIA, ER, NIBM, G11, &
    VPOL, VRS, VTOR, SHEAR, PBLON, PBPER, PFAST, NIZ3, ZIM2, ZIM3, &
! output
    CHI, CHE, DIF, VIN, DPH, DPL, DPR, XTB, EGM, GAM, GM1, GM2, OM1, OM2, FR1)

!----------------------------------------------------------------------|
! based on stand-alone driver for the GLF23 model
!       "testglf.f" 18-fev-03 version 1.61
!       written by Jon Kinsey, General Atomics
!----------------------------------------------------------------------|
! WORK(1:NA1,1:13) array is used for output
!                              (when i_delay=0 and egamma_d is not used)
!----------------------------------------------------------------------|

USE tglf_interface
use tglf_pkg, only: get_eigenvalue_spectrum_out, get_ky_spectrum_out, &
    get_flux_spectrum_out

implicit none

integer, parameter :: jpd=700, nradial=9

double precision, parameter :: &
   k0   = 1.6022d-12, &       ! erg/ev
   e0   = 4.8032d-10, &       ! elementary charge (statcoulombs)
   e00  = 1.6020d-19, &       ! elementary charge (C)
   c0   = 2.9979d+10, &       ! speed of light (cm/sec)
   mp   = 1.6726d-24, &       ! proton mass (g)
   mpp  = 1.6726d-27, &       ! proton mass (kg)
   pi   = 3.141592653589793 
double precision, parameter :: c_vpol = 1.d0

integer, intent(in) :: nrho, jr1_in, jr2_in, NA1N, NA1E, NA1I

double precision, intent(in) :: BTOR, RTOR, &
    AMJ, AIM1, AIM2, AIM3, ZMJ
double precision, intent(in), dimension(*) :: NE, TE, NI, TI, &
    ZEF, ZIM1, ZIM2, ZIM3, PBLON, PBPER, PFAST, NIZ3, AMAIN, &
    ER, MU, RHO, AMETR, SHIF, ELON, NDEUT, NIZ1, NTRIT, &
    NIZ2, TRIA, NIBM, G11, VPOL, VRS, VTOR, SHEAR

double precision, intent(out), dimension(*) :: CHI, CHE, DIF, VIN, &
    DPH, DPL, DPR, XTB, EGM, GAM, GM1, GM2, OM1, OM2, FR1

!----------------------------------------------------------------------

integer :: jna, jinterval, jr_min, jr_max, jrho, j0, j01, j02, jgamma_max
integer :: j, jradial, j3r, jjgrid(nradial), j1, j2, jspec, kyloop
integer :: sat_rule           ! Saturation rule
integer :: nmodes_tg          ! number of unstable modes to use in computing fluxes (max=4)
integer :: kygrid_model_tg    ! select version of ky-grid to use 1
integer :: xnu_model_tg       ! select version of trapped-passing 2

real :: bmod, bpolz, alpha_zf_in, ion_eflux
real :: drmin, drmaj, drho, dte, dne, dq, dptot, &
        delong, dtrian, dvper, drhodr, dstep, dr, dv_r
real :: Bunit, cs0, cs00, rhos0, omega0, rhostar2, lnlamda, taue, cexb
real :: a0, T0, N0, m0, rmin_tg, drho_cs, drho_nt, nt_cs
real :: wdia_trap_tg          ! parameter for trapped fraction model

real, dimension(jpd) :: rho_m, vexb2, vpar_m, vper_m, &
    gradrhosq_exp, rmaj_exp, q_exp, &
    chie_m, chii_m, elec_pflux_m, exchi_m, ptot, gamma_m, omega_m
real, dimension(nsm-1,jpd) :: ion_pflux_m

real, dimension(nradial) :: chie, chii, exchi, elec_pflux, rho_tg, gamma_max, omega_max, kymax
real, dimension(nsm-1, nradial) :: ion_pflux
real, allocatable, dimension(:) :: gamma, omega, kyspectrum, efluxspectrum
real, dimension(nsm-1) :: dti, dni
real, dimension(nsm-1, jpd) :: ni_m, ti_m

!--------------

jna = max(NA1E, NA1I, NA1N)
if (jna == 0) jna = nrho
if (jna > jpd) then
    return
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

m0 = AMJ*mp                ! Ref. mass = D ion mass [g]
a0 = 1E2*AMETR(nrho)    ! length scale used by GYRO from AMETR (meters) to cm

elec_pflux_m = 0.0
ion_pflux_m  = 0.0
chie_m  = 0.0
chii_m  = 0.0
exchi_m = 0.0

jr_min = max(1, jr1_in)
jr_max = min(nrho, jr2_in)
if (jr_min >= jr_max) return

jinterval = 1 + (jr_max - jr_min)
j3r = max(1, INT(jinterval/nradial))
do jradial = 1, nradial-1
    jjgrid(jradial) = jr_min + (jradial - 1)*j3r
enddo
jjgrid(nradial) = min(nrho, jr_max)

! Number of species

tglf_ns_in = nsm
! These will be reset locally in the radial loop
tglf_zs_in(3) = MAXVAL(ZIM1(1:nrho))
tglf_zs_in(4) = MAXVAL(ZIM2(1:nrho))
tglf_zs_in(5) = MAXVAL(ZIM3(1:nrho))

if (tglf_zs_in(5) >= 1.) tglf_ns_in = 5
if (tglf_zs_in(5) < 1.) tglf_ns_in = 4
if (tglf_zs_in(4) < 1.) tglf_ns_in = 3
if (tglf_zs_in(5) >= 1. .and. tglf_ns_in == 3) then 
    tglf_ns_in = 4
endif
if (tglf_zs_in(3) < 1.) tglf_ns_in = 2
if (tglf_zs_in(4) >= 1. .and. tglf_ns_in == 2) then 
    tglf_ns_in = 3
endif

kygrid_model_tg = 1

sat_rule = 2
write(6, '(A, 6i4)') 'Call TGLF...', jr1_in, jr2_in, nrho, jna, sat_rule, tglf_ns_in

if (sat_rule == 0) then
    nmodes_tg = 2
    xnu_model_tg = 2
    wdia_trap_tg = 0.
    alpha_zf_in  = 0.
endif
if (sat_rule == 1) then
    nmodes_tg = tglf_ns_in + 2
    xnu_model_tg = 2
    wdia_trap_tg = 0.
    alpha_zf_in  = 1.
endif
if (sat_rule == 2) then
    nmodes_tg = tglf_ns_in + 2
    xnu_model_tg = 3
    wdia_trap_tg = 1.
    alpha_zf_in  = 1.
endif

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
tglf_nbasis_max_in = 4
tglf_nbasis_min_in = 2
tglf_nxgrid_in     = 16
tglf_nky_in        = 19
tglf_units_in   = 'CGYRO'
tglf_path_in = '../tglf/'
! Want fluxes from TGLF
tglf_use_transport_model_in = .true.

tglf_sign_Bt_in = 1
tglf_sign_It_in = 1

tglf_ky_in = 0.3
tglf_width_in      = 1.65
tglf_width_min_in  = 0.3
tglf_nwidth_in     = 21

tglf_geometry_flag_in = 1
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

radial_loop: do jradial=1, nradial
   
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

    j01 = j0+1
    j02 = j0-1
    if (j0 == 1) then
        j02 = j0
    else if (j0 == nrho) then
        j01 = j0
    endif
    dstep = 1./float(j01 - j02)

    drmin  = dstep*(AMETR(j01) - AMETR(j02))
    drmaj  = dstep*(rmaj_exp(j01) - rmaj_exp(j02))
    drho   = dstep*(rho(j01) - rho(j02))
    delong = dstep*(ELON(j01) - ELON(j02))
    dtrian = dstep*(TRIA(j01) - TRIA(j02))
    dptot  = dstep*(ptot(j01) - ptot(j02)) * 1E3*1E13
    dte    = dstep*(TE(j01) - TE(j02))
    dne    = dstep*(NE(j01) - NE(j02))
    dq     = dstep*(q_exp(j01) - q_exp(j02))
    dvper  = dstep*(vper_m(j01) - vper_m(j02))
    do jspec=1, tglf_ns_in-1
        dti(jspec) = dstep*(ti_m(jspec, j01) - ti_m(jspec, j02))
        dni(jspec) = dstep*(ni_m(jspec, j01) - ni_m(jspec, j02))
    enddo
    dv_r = dstep* &
        (vpar_m(j01)/(rmaj_exp(j01) + AMETR(j01)) - &
         vpar_m(j02)/(rmaj_exp(j02) + AMETR(j02)))
    dr = 1E2*drmin/a0    ! gradients w.r.t. minor radius even for s-alpha geometry
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

    rmin_tg = 1E2*AMETR(j0)/a0
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

    rhostar2 = (rhos0/a0)**2
    drho_cs = drhodr**2*a0/1e2*rhostar2*cs00
    drho_nt = drhodr**2*a0/1e2*rhostar2*N0*e00*T0/(1.e13*mpp)
    nt_cs = 0.0016*N0/1.e13*T0/1.e3*cs00/(a0/1.e2)*rhostar2

! Share variables with tglf_run via module tglf_interface

    tglf_vexb_shear_in = -1E2*cexb*dvper/(dr*cs0) ! Waltz-Miller definition     !From m/s to cm/s for vexb  -dVexb/dr

! Initialise

    tglf_vexb_in  = 1E2*vexb2(j0)/cs0
    tglf_betae_in = 8.0*pi*k0*N0*T0/Bunit**2
    tglf_xnue_in  = 0.75*SQRT(pi)*a0/(taue*cs0)
    tglf_zeff_in  = ZEF(j0) 
    tglf_debye_in = SQRT(k0*T0/(4.0*pi*N0*e0**2))/rhos0

    tglf_rmin_loc_in    = 1E2*AMETR(j0)/a0
    tglf_rmaj_loc_in    = 1E2*rmaj_exp(j0)/a0
    tglf_zmaj_loc_in    = 0.
    tglf_drmindx_loc_in = 1.
    tglf_drmajdx_loc_in = 1E2*drmaj/(dr*a0)
    tglf_dzmajdx_loc_in = 0.
    tglf_kappa_loc_in   = ELON(j0)
    tglf_s_kappa_loc_in = AMETR(j0)*delong/drmin
    tglf_delta_loc_in   = TRIA(j0)
    tglf_s_delta_loc_in = AMETR(j0)*dtrian/drmin
    tglf_zeta_loc_in    = 0.
    tglf_s_zeta_loc_in  = 0.
    tglf_q_loc_in       = q_exp(j0)
    tglf_q_prime_loc_in = (q_exp(j0)/rmin_tg)*dq/dr
    tglf_p_prime_loc_in = (k0/Bunit**2)*(q_exp(j0)/rmin_tg)*dptot/dr

    tglf_rmin_sa_in     = 1E2*AMETR(j0)/a0
    tglf_rmaj_sa_in     = 1E2*rmaj_exp(j0)/a0
    tglf_q_sa_in        = q_exp(j0)
    tglf_shat_sa_in     = (AMETR(j0)/q_exp(j0))*dq/drmin
    tglf_alpha_sa_in    = -(8.0*pi*k0/Bunit**2)*q_exp(j0)**2 * rmaj_exp(j0)*dptot/drmin
    tglf_xwell_sa_in    = 0.
    tglf_theta0_sa_in   = 0.

! Settings
    if (tglf_dump_flag_in) then
        write(file_dump_local, '(A11, I0)') 'input.tglf_', j0
    endif

    call tglf_run

! Transport coefficients

    ion_eflux = SUM(tglf_ion_eflux_out(1: tglf_ns_in-1))
    ion_eflux = ion_eflux/(tglf_taus_in(2) * 1e13*NI(j0)/N0)
    chii(jradial) = ion_eflux          /(1e-4 + abs(tglf_rlts_in(2))) *drho_cs
    chie(jradial) = tglf_elec_eflux_out/(1e-4 + abs(tglf_rlts_in(1))) *drho_cs
    elec_pflux(jradial) = tglf_elec_pflux_out/drhodr *drho_cs         ! particle flux
    do jspec=1, tglf_ns_in-1
        ion_pflux(jspec, jradial) = tglf_ion_pflux_out(jspec)/drhodr *drho_cs  !ion particle flux
    enddo
    exchi(jradial) = tglf_ion_expwd_out(1) * nt_cs                ! Equipartition
    do kyloop=1, tglf_nky_in 
        gamma(kyloop) = get_eigenvalue_spectrum_out(1, kyloop, 1)
        omega(kyloop) = get_eigenvalue_spectrum_out(2, kyloop, 1)
        kyspectrum(kyloop) = get_ky_spectrum_out(kyloop)
        efluxspectrum(kyloop) = get_flux_spectrum_out(2, 1, 1, kyloop, 1)
    enddo

    jgamma_max = maxloc(efluxspectrum(1:tglf_nky_in), 1)
    gamma_max(jradial) = gamma(jgamma_max)
    omega_max(jradial) = omega(jgamma_max)
    kymax(jradial) = kyspectrum(jgamma_max)
!    write(6, '(A,i3,2e10.2)') 'tglf freq', jradial, gamma_max(jradial), omega_max(jradial)

enddo radial_loop

j1 = jjgrid(1)
j2 = jjgrid(nradial)

call qinterp(rho_tg, chii  , nradial, rho_m(j1:j2), chii_m(j1:j2)  , j2-j1+1)
call qinterp(rho_tg, chie  , nradial, rho_m(j1:j2), chie_m(j1:j2)  , j2-j1+1)
call qinterp(rho_tg, elec_pflux, nradial, rho_m(j1:j2), elec_pflux_m(j1:j2), j2-j1+1)
call qinterp(rho_tg, exchi , nradial, rho_m(j1:j2), exchi_m(j1:j2) , j2-j1+1)
call qinterp(rho_tg, gamma_max, nradial, rho_m(j1:j2), gamma_m(j1:j2) , j2-j1+1)
call qinterp(rho_tg, omega_max, nradial, rho_m(j1:j2), omega_m(j1:j2) , j2-j1+1)
do jspec=1, tglf_ns_in-1
    call qinterp(rho_tg, ion_pflux(jspec, 1:nradial), nradial, rho_m(j1:j2), ion_pflux_m(jspec, j1:j2), j2-j1+1)
enddo

chii_m  (1:2) = chii_m(3)
chie_m  (1:2) = chie_m(3)
elec_pflux_m(1:2) = elec_pflux_m(3)
exchi_m (1:2) = exchi_m(3)
omega_m (1:2) = omega_m(3)
gamma_m (1:2) = gamma_m(3)

DIF(1:nrho) = 0.d0       ! D, electron diffusivity, m^2/s
DPH(1:nrho) = 0.d0       ! D, impurity diffusivity, m^2/s
DPL(1:nrho) = 0.d0       ! impurity convection
DPR(1:nrho) = 0.d0       ! tor. stress
EGM(1:nrho) = 0.d0
GAM(1:nrho) = 0.d0
GM1(1:nrho) = 0.d0
GM2(1:nrho) = 0.d0
OM1(1:nrho) = 0.d0
OM2(1:nrho) = 0.d0
FR1(1:nrho) = 0.d0

do j=jr_min, jr_max
    CHI(j) = chii_m(j)/gradrhosq_exp(j) ! \chi_i, m^2/s -> work(21,:) 
    CHE(j) = chie_m(j)/gradrhosq_exp(j) ! \chi_e, m^2/s
    VIN(j) = elec_pflux_m(j)/AMETR(nrho)/gradrhosq_exp(j) ! D flux
!First impurity only, index 2 of ion species
    if (tglf_ns_in >= 3) then
       DPL(j) = ion_pflux_m(2, j)/AMETR(nrho)/gradrhosq_exp(j)/(ni_m(2, j)/NE(j))  ! 1st imp convection
    endif
    if (tglf_ns_in >= 4) then
       DPH(j) = ion_pflux_m(3, j)/AMETR(nrho)/gradrhosq_exp(j)/(ni_m(3, j)/NE(j))  ! 2nd imp convection
    endif
    XTB(j) = exchi_m(j)  ! turbulent e-i equipartition in MW/m^3

    CHI(j) = max( -10., min(10., CHI(j)) )
    CHE(j) = max( -10., min(10., CHE(j)) )
    VIN(j) = max( -20., min(20., VIN(j)) )
    GM1(j) = gamma_m(j)*(cs0/a0)
    OM1(j) = omega_m(j)*(cs0/a0)
enddo

DPL(1) = 0.d0 !ensure NIZ1 convection equal to zero on axis
DPH(1) = 0.d0 !ensure NIZ2 convection equal to zero on axis (already 0 otherwise)

if (jr_max < nrho-1) return

do j=jna-2, jna-1
    CHI(j) = CHI(jna-3)
    CHE(j) = CHE(jna-3)
    VIN(j) = VIN(jna-3)
    XTB(j) = XTB(jna-3)
enddo

if (jna < nrho) then
    do j=jna-1, nrho
        CHI(j) = 0.d0
        CHE(j) = 0.d0
        VIN(j) = 0.d0
        XTB(j) = 0.d0
    enddo
endif  

return
END subroutine tglf_interf
