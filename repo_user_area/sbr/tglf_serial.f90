!----------------------------------------------------------------------|
subroutine tglf_serial(CHI, CHE, VIN, XTB)

!----------------------------------------------------------------------|
! based on stand-alone driver for the GLF23 model
!       "testglf.f" 18-fev-03 version 1.61
!       written by Jon Kinsey, General Atomics
!----------------------------------------------------------------------|
! WORK(1:NA1,1:13) array is used for output
!                              (when i_delay=0 and egamma_d is not used)
!----------------------------------------------------------------------|

use parameter_inc, only: NRD
USE tglf_interface

use const_inc, only: BTOR, RTOR, AMJ, AIM1, AIM2, AIM3, ZMJ, &
    NA1, NA1N, NA1E, NA1I
use status_inc, only: NE, TE, NI, TI, &
    ZEF, ZIM1, ZIM2, ZIM3, PBLON, PBPER, &
    PFAST, NIZ3, AMAIN, ER, MU, &
    RHO, AMETR, SHIF, ELON, &
    NDEUT, NIZ1, NTRIT, NIZ2, NHE3, &
    TRIA, VTOR, NIBM, G11, VRS, SHEAR

implicit none

integer, parameter :: nradial = 100
real, parameter :: &
   k0   = 1.6022E-12, &       ! erg/ev
   e0   = 4.8032E-10, &       ! elementary charge (statcoulombs)
   e00  = 1.6020e-19, &       ! elementary charge (C)
   c0   = 2.9979E+10, &       ! speed of light (cm/sec)
   mp   = 1.6726E-24, &       ! proton mass (g)
   mpp  = 1.6726E-27, &       ! proton mass (kg)
   pi   = 3.141592653589793 

! Parameters used:
! In A_vars & A_arrs:

double precision, dimension(*), intent(out) :: CHI, CHE, VIN, XTB

!----------------------------------------------------------------------
integer :: jna, jrho, j0, j01, j02, n_red
integer :: j, jradial, j1, j2, jspec, n_step
integer :: jjgrid(100)
integer :: sat_rule           ! Saturation rule
integer :: nmodes_tg          ! number of unstable modes to use in computing fluxes (max=4)
integer :: kygrid_model_tg    ! select version of ky-grid to use 1
integer :: xnu_model_tg       ! select version of trapped-passing 2

real :: bmod, bpolz, alpha_zf_in
real :: drmin, drmaj, drho, dte, dne, dq, dptot, &
        delong, dtrian, dvper, drhodr, dstep, dr, dv_r
real :: Bunit, cs0, cs00, rhos0, omega0, rhostar2, lnlamda, taue, cexb
real :: a0, T0, N0, m0, rmin_tg, drho_cs
real :: wdia_trap_tg          ! parameter for trapped fraction model

real, dimension(NRD) :: rho_m, vexb2, vpar_m, vper_m, &
    gradrhosq_exp, rmaj_exp, q_exp, &
    chie_m, chii_m, pfluxi_m, exchi_m, ptot

real, dimension(nradial) :: chie, chii, exchi, pfluxi, rho_tg
real, dimension(nsm-1) :: dti, dni
real, dimension(nsm-1, NRD) :: ni_m, ti_m
character(10) :: time_loc

!--------------

call DATE_AND_TIME(TIME=time_loc)

jna = max(NA1E, NA1I, NA1N)
if (jna .eq. 0) jna = NA1

! Electrons and main ions
tglf_zs_in(1) = -1.
tglf_zs_in(2) = ZMJ

tglf_mass_in(1) = 5.4447e-4/AMJ
tglf_mass_in(2) = AMJ/AMJ  ! AMJ is reference mass
tglf_mass_in(3) = AIM1/AMJ
tglf_mass_in(4) = AIM2/AMJ
tglf_mass_in(5) = AIM3/AMJ

do jrho=1, NA1
    rho_m(jrho) = RHO(jrho)
    ti_m(1:4, jrho) = TI(jrho)
    ni_m(1, jrho) = NI(jrho)
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
    vpar_m(jrho) = ER(jrho)/(RTOR*bpolz)*(RTOR+SHIF(jrho)+AMETR(jrho))  !--> R*Omega_E , no neoclassical terms
    vexb2(jrho)  = -ER(jrho)/bmod ! vexb in m/s (vperp = vexb since the diamagnetic velocity is the curvature drift actually

enddo

m0 = AMJ*mp                ! Ref. mass = D ion mass [g]
a0 = 1E2*AMETR(NA1)    ! length scale used by GYRO from AMETR (meters) to cm

pfluxi_m  = 0.0
chie_m  = 0.0
chii_m  = 0.0
exchi_m = 0.0
n_red = 25
n_step = int(jna/n_red) + 1
do jradial = 1, n_red
    jjgrid(jradial) = 1 + (jradial - 1)*n_step
    if (jjgrid(jradial) >= jna) EXIT
enddo
n_red = jradial
jjgrid(n_red) = jna
write(*, *) 'GIT tglf', jjgrid(1:n_red)

! Number of species

tglf_ns_in = nsm
! These will be reset locally in the radial loop
tglf_zs_in(3) = MAXVAL(ZIM1(1:NA1))
tglf_zs_in(4) = MAXVAL(ZIM2(1:NA1))
tglf_zs_in(5) = MAXVAL(ZIM3(1:NA1))

if (tglf_zs_in(5) .ge. 1.) tglf_ns_in = 5
if (tglf_zs_in(5) .lt. 1.) tglf_ns_in = 4
if (tglf_zs_in(4) .lt. 1.) tglf_ns_in = 3
if (tglf_zs_in(5) .ge. 1. .and. tglf_ns_in .eq. 3) then 
    tglf_ns_in = 4
endif
if (tglf_zs_in(3) .lt. 1.) tglf_ns_in = 2
if (tglf_zs_in(4) .ge. 1. .and. tglf_ns_in .eq. 2) then 
    tglf_ns_in = 3
endif

kygrid_model_tg = 1

sat_rule = 2

write(6, *) time_loc, ' BEGIN tglf_serial'
write(6, *) jna, sat_rule, tglf_ns_in, NA1N, NA1E, NA1I

if (sat_rule .eq. 0) then
    nmodes_tg = 2
    xnu_model_tg = 2
    wdia_trap_tg = 0.
    alpha_zf_in  = 0.
endif
if (sat_rule .eq. 1) then
    nmodes_tg = tglf_ns_in + 2
    xnu_model_tg = 2
    wdia_trap_tg = 0.
    alpha_zf_in  = 1.
endif
if (sat_rule .eq. 2) then
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
tglf_nxgrid_in     = 12
tglf_nky_in        = 19
tglf_units_in   = 'CGYRO'
tglf_path_in = '/toks/work/git/'
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

radial_loop: do jradial=1, n_red
   
    j0 = jjgrid(jradial)

!thermal impurities

    tglf_zs_in(3) = max(1., ZIM1(j0))
    tglf_zs_in(4) = ZIM2(j0)
    tglf_zs_in(5) = ZIM3(j0)

    if (tglf_zs_in(5) .ge. 1. .and. tglf_ns_in .eq. 3) then 
        tglf_zs_in(4) = tglf_zs_in(5)
        tglf_mass_in(4) = tglf_mass_in(5)
        ni_m(3, :) = ni_m(4, :)
        ti_m(3, :) = ti_m(4, :)
    endif
    if (tglf_zs_in(4) .ge. 1. .and. tglf_ns_in .eq. 2) then 
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
    else if (j0 == NA1) then
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

    if (tglf_ns_in .ge. 3) then
        tglf_vpar_shear_in(3: tglf_ns_in) = tglf_vpar_shear_in(2)
        tglf_vpar_in      (3: tglf_ns_in) = tglf_vpar_in(2)
    endif

! local magnetic geometry

    rhostar2 = (rhos0/a0)**2
    drho_cs = drhodr**2*a0/1e2*rhostar2*cs00

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

    call tglf_run

! Transport coefficients

!    chii(jradial)   = tglf_ion_eflux_low_out(1)/(1e-4 + abs(tglf_rlts_in(2))) *drho_cs
!    chie(jradial)   = tglf_elec_eflux_low_out  /(1e-4 + abs(tglf_rlts_in(1))) *drho_cs
    chii(jradial)   = tglf_ion_eflux_out(1)/(1e-4 + abs(tglf_rlts_in(2))) *drho_cs
    chie(jradial)   = tglf_elec_eflux_out  /(1e-4 + abs(tglf_rlts_in(1))) *drho_cs
    pfluxi(jradial) = tglf_ion_pflux_out(1)/drhodr *drho_cs         ! particle flux
    exchi(jradial)  = tglf_ion_expwd_out(1) *drho_cs                ! Equipartition

enddo radial_loop

call qinterp(rho_tg, chii  , n_red, rho_m(1:jna), chii_m(1:jna)  , jna)
call qinterp(rho_tg, chie  , n_red, rho_m(1:jna), chie_m(1:jna)  , jna)
call qinterp(rho_tg, pfluxi, n_red, rho_m(1:jna), pfluxi_m(1:jna), jna)
call qinterp(rho_tg, exchi , n_red, rho_m(1:jna), exchi_m(1:jna) , jna)

chii_m  (1:2) = chii_m(3)
chie_m  (1:2) = chie_m(3)
pfluxi_m(1:2) = pfluxi_m(3)
exchi_m (1:2) = exchi_m(3)

do j=1, jna
    CHI(j) = chii_m(j)/gradrhosq_exp(j) ! \chi_i, m^2/s : starts from work(21,:) 
    CHE(j) = chie_m(j)/gradrhosq_exp(j) ! \chi_e, m^2/s
    VIN(j) = min(20., pfluxi_m(j)/AMETR(NA1)/gradrhosq_exp(j)) ! D flux
    VIN(j) = max(-20., VIN(j)) ! D flux
!    VIN(j) = pfluxi_m(j)/AMETR(NA1)/gradrhosq_exp(j)
    XTB(j) = exchi_m(j)  ! turbulent e-i equipartition in MW/m^3
enddo   ! End of main loop

do j=jna-2, jna-1
    CHI(j) = CHI(jna-3)
    CHE(j) = CHE(jna-3)
    VIN(j) = VIN(jna-3)
    XTB(j) = XTB(jna-3)
enddo

if (jna .lt. NA1) then
    do j=jna-1, NA1
        CHI(j) = 0.d0
        CHE(j) = 0.d0
        VIN(j) = 0.d0
        XTB(j) = 0.d0
    enddo
endif  

call DATE_AND_TIME(TIME=time_loc)
write(*, *) time_loc, ' END tglf_serial.f90'

return
END subroutine tglf_serial
