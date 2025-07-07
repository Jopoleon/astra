program qlk_chunk

use mpi
USE mod_qualikiz, only: qualikiz
USE kind, only: qlk_output_meth_0, qlk_output_meth_0_sep_0, &
    qlk_primi_meth_0, qlk_sizes, qlk_in_regular, qlk_in_newt
USE mod_make_io, only: zeroout_qlk_in_newt, zeroout_qlk_in_regular, &
    allocate_qlk_in_regular, allocate_qlk_in_newt
USE nanfilter

implicit none
type(qlk_sizes)      :: sizes
type(qlk_in_regular) :: in_regular
type(qlk_in_newt)    :: in_newt

type(qlk_output_meth_0)       :: output_meth_0
type(qlk_output_meth_0_sep_0) :: output_meth_0_sep_0_SI , output_meth_0_sep_0_GB
type(qlk_primi_meth_0)        :: primi_meth_0


integer, parameter :: ntheta=64, numecoefs=13, numicoefs=7, dimx=1, dimn=16, numsols=3, phys_meth=0, nradial=5, nspec_max=7

double precision, parameter :: &
    k0  = 1.6022E-12, &       ! erg/ev
    e0  = 4.8032E-10, &       ! elementary charge (statcoulombs)
    e00 = 1.6020e-19, &       ! elementary charge (C)
    c0  = 2.9979E+10, &       ! speed of light (cm/sec)
    mp  = 1.6726E-24, &       ! proton mass (g)
    mpp = 1.6726E-27, &       ! proton mass (kg)
    pi  = 3.141592653589793 

!-----------------------------------------

integer :: simple_mpi_only_in, maxpts_in, maxruns_in, runcounter_in
integer :: nions, coll_flag_in, rot_flag_in, verbose_in, el_type_in, &
    integration_routine_in, separateflux_in
integer, dimension(dimx, nspec_max-1) :: ion_type_in

double precision :: relacc1_in, relacc2_in, absacc1_in, absacc2_in, R0_in, &
    ETGmultin, collmultin, timeout_in, rhomin, rhomax, rhoscale, xstep
double precision, dimension(dimn) :: kthetarhos_in
double precision, dimension(dimx) :: x_in, rho_in, Ro_in, Rmin_in, Bo_in, &
    qx_in, smag_in, alphax_in, Tex_in, Nex_in, Ate_in, Ane_in, anise_in, &
    danisedr_in, Machtor_in, Autor_in, Machpar_in, Aupar_in, gammaE_in, &
    epf_GB_out, eef_GB_out
double precision, dimension(dimx, nspec_max-1) :: Tix_in, ninorm_in, &
    Ati_in, Ani_in, anis_in, danisdr_in, Ai_in, Zi_in, ipf_GB_out, ief_GB_out
double precision, dimension(dimx, dimn, numsols) :: gam_GB_out, ome_GB_out
double precision, dimension(dimx, nspec_max - 1, numicoefs) :: cftrans_out

! Old solution for non-reset runs
double precision, DIMENSION(:, :, :), ALLOCATABLE :: oldrsol, oldisol, oldrfdsol, oldifdsol

!-----------------------------------------
LOGICAL :: exist1, exist2, exist3, exist4, exist5 !used for checking for existence of files

!MPI variables:
INTEGER :: mpi_ierr, nproc, myrank
INTEGER :: myunit=700, i_mpic
integer :: jr_min, jr_max, jrho, j0, j01, j02, n_radial
integer :: i, j, k, jradial, jjgrid(nradial), jion

double precision :: bmod, bpolz
double precision :: drmin, drmaj, drho, dte, dne, dq, dptot, &
    delong, dtrian, dvper, drhodr, dstep, dr, dv_r
double precision :: Bunit, cs00, rhos00, omega0, rhostar2
double precision :: T0, m0, drho_cs
! Shifted cicle geometry inputs
double precision :: gamma_e_tg, mach_fac, ql_fac

double precision, dimension(nrho) :: vexb2, vpar_m, vper_m, &
    gradrhosq_exp, rmaj_exp, q_exp, &
    chie_m, chii_m, pfluxi_m, exchi_m, ptot

double precision, dimension(nradial) :: chie, chii, exchi, pfluxi, rho_tg
double precision, dimension(nspec_max-1) :: dti, dni
double precision, dimension(nspec_max-1, nrho) :: ni_m, ti_m
double precision :: vpar_in, vpar_shear_in, cexb

COMPLEX(kind=DBL), DIMENSION(:, :, :), ALLOCATABLE :: oldsol_in, oldfdsol_in
complex(kind=DBL), dimension(dimx, dimn, numsols) :: sol_out, fdsol_out

CHARACTER(len=20) :: fmtn
character(len=80) :: fname1, prim_dir


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
call MPI_Recv(dims, 6, MPI_INTEGER, 0, 0, parent, status, ierr)
nrho_tg   = dims(1)
n_scalars = dims(2)
n_inputs  = dims(3)
n_outputs = dims(4)
tglf_ns_in  = dims(5)
tglf_nky_in = dims(6)

chunk = nrho_tg / nprocs  ! Safe here: we now know nrho_tg
allocate(scalars(n_scalars))
allocate(inputs(chunk, n_inputs))
allocate(output(chunk, n_outputs))
allocate( mtori(chunk), chie(chunk), chii(chunk), exchi(chunk), elec_pflux(chunk), &
    rho_tg(chunk), gamma_max(chunk), omega_max(chunk), kymax(chunk), &
    te_tg(chunk), ne_tg(chunk), vpar_tg(chunk), vper_tg(chunk), vexb_tg(chunk), &
    ametr_tg(chunk), elon_tg(chunk), tria_tg(chunk), rmaj_tg(chunk), &
    ptot_tg(chunk), q_tg(chunk), zef_tg(chunk), pfn_tg(chunk) )
allocate( drmin(chunk), drmaj(chunk), drho(chunk), delong(chunk), dtrian(chunk), &
    dr(chunk), dne(chunk), dte(chunk), dq(chunk), dptot(chunk), dvpar(chunk), &
    dvper(chunk), dv_r(chunk), drhodr(chunk) )
allocate( dti(4, chunk), dni(4, chunk), ni_tg(4, chunk), ti_tg(4, chunk), &
    zimp_tg(3, chunk), ion_pflux(4, chunk) )
allocate( gamma(tglf_nky_in), omega(tglf_nky_in), kyspectrum(tglf_nky_in), efluxspectrum(tglf_nky_in), &
    ifluxspectrum(tglf_nky_in), pfluxspectrum(tglf_nky_in) )
! Receive TGLF input scalars and profiles from parent
call MPI_Recv(scalars     , n_scalars, MPI_DOUBLE_PRECISION, 0, 0, parent, status, ierr)
call MPI_Recv(inputs, chunk * n_inputs, MPI_DOUBLE_PRECISION, 0, 0, parent, status, ierr)
rho_tg      = inputs(:,  1)
ametr_tg    = inputs(:,  2)
rmaj_tg     = inputs(:,  3)
elon_tg     = inputs(:,  4)
tria_tg     = inputs(:,  5)
q_tg        = inputs(:,  6)
pfn_tg      = inputs(:,  7)
ptot_tg     = inputs(:,  8)
ne_tg       = inputs(:,  9)
te_tg       = inputs(:, 10)
zef_tg      = inputs(:, 11)
vpar_tg     = inputs(:, 12)
vper_tg     = inputs(:, 13)
vexb_tg     = inputs(:, 14)
ti_tg(1, :) = inputs(:, 15)
ti_tg(2, :) = inputs(:, 16)
ti_tg(3, :) = inputs(:, 17)
ti_tg(4, :) = inputs(:, 18)
ni_tg(1, :) = inputs(:, 19)
ni_tg(2, :) = inputs(:, 20)
ni_tg(3, :) = inputs(:, 21)
ni_tg(4, :) = inputs(:, 22)
zimp_tg(1, :) = inputs(:, 23)
zimp_tg(2, :) = inputs(:, 24)
zimp_tg(3, :) = inputs(:, 25)
drmin       = inputs(:, 26)
drmaj       = inputs(:, 27)
drho        = inputs(:, 28)
delong      = inputs(:, 29)
dtrian      = inputs(:, 30)
dptot       = inputs(:, 31)
dte         = inputs(:, 32)
dne         = inputs(:, 33)
dq          = inputs(:, 34)
dvper       = inputs(:, 35)
dv_r        = inputs(:, 36)
dr          = inputs(:, 37)
drhodr      = inputs(:, 38)
dti(1, :)   = inputs(:, 39)
dti(2, :)   = inputs(:, 40)
dti(3, :)   = inputs(:, 41)
dti(4, :)   = inputs(:, 42)
dni(1, :)   = inputs(:, 43)
dni(2, :)   = inputs(:, 44)
dni(3, :)   = inputs(:, 45)
dni(4, :)   = inputs(:, 46)

AMJ  = scalars(1)
BTOR = scalars(2)
a0_m = scalars(3)
m0 = AMJ*mp          ! Ref. mass = D ion mass [g]
a0_cm = 1.d2*a0_m    ! length scale used by GYRO, m -> cm

! Electrons and main ions
Zi_in(1, 1) = ZMJ

Ai_in(1, 1) = AMJ
Ai_in(1, 2) = AMJ*scalars(6)
Ai_in(1, 3) = AMJ*scalars(7)
Ai_in(1, 4) = AMJ*scalars(8)

radial_loop: do jr=1, chunk

!thermal impurities

    tglf_zs_in(3) = max(1., zimp_tg(1, jr))
    tglf_zs_in(4) = zimp_tg(2, jr)
    tglf_zs_in(5) = zimp_tg(3, jr)

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

    tglf_as_in(1)   = 1. ! ne is ref
    tglf_taus_in(1) = 1. ! Te is ref

! Log derivatives
    tglf_rlns_in(1) = -dne(jr)/(dr(jr)*ne_tg(jr))
    tglf_rlts_in(1) = -dte(jr)/(dr(jr)*te_tg(jr))

    do jspec=2, tglf_ns_in
        tglf_as_in(jspec)   = ni_tg(jspec-1, jr)/ne_tg(jr)
        tglf_taus_in(jspec) = ti_tg(jspec-1, jr)/te_tg(jr)
        tglf_rlns_in(jspec) = -dni(jspec-1, jr)/(dr(jr)*ni_tg(jspec-1, jr))
        tglf_rlts_in(jspec) = -dti(jspec-1, jr)/(dr(jr)*ti_tg(jspec-1, jr))
    enddo
! Restore quasi-neutrality via main ions

    tglf_as_in(2)   = -1./tglf_zs_in(2)*(SUM(tglf_as_in*tglf_zs_in) - tglf_as_in(2)*tglf_zs_in(2))
    tglf_rlns_in(2) = -1./(tglf_as_in(2)*tglf_zs_in(2))*(SUM(tglf_rlns_in*tglf_as_in*tglf_zs_in) - tglf_rlns_in(2)*tglf_as_in(2)*tglf_zs_in(2))

! GYRO conventions

    N0  = 1E13*ne_tg(jr)   ! density scale used by GYRO [1/cm**3]
    T0  = 1E3 *te_tg(jr)   ! temperature scale used by GYRO
    Bunit_T = BTOR*drhodr(jr)*rho_tg(jr)/ametr_tg(jr)  ! Miller geometry magnetic field unit [gauss]
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

    tglf_vpar_shear_in(2) = -1E2*rmaj_tg(jr)*dv_r(jr)/(dr(jr)*cs0)  !From m/s to cm/s for vpar
    tglf_vpar_shear_in(1) = tglf_vpar_shear_in(2)

    tglf_vpar_in(2) = 1E2*vpar_tg(jr)/cs0
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

    tglf_vexb_in  = 1E2*vexb_tg(jr)/cs0
    tglf_betae_in = 8.0*pi*k0*N0*T0/Bunit_gauss**2
    tglf_xnue_in  = 0.75*SQRT(pi)*a0_cm/(taue*cs0)
    tglf_zeff_in  = zef_tg(jr)
    tglf_debye_in = SQRT(k0*T0/(4.0*pi*N0*e0**2))/rhos0

    tglf_rmin_loc_in    = rmin_tg
    tglf_rmaj_loc_in    = rmaj_tg(jr)/a0_m
    tglf_zmaj_loc_in    = 0.
    tglf_drmindx_loc_in = 1.
    tglf_drmajdx_loc_in = drmaj(jr)/(dr(jr)*a0_m)
    tglf_dzmajdx_loc_in = 0.
    tglf_kappa_loc_in   = elon_tg(jr)
    tglf_s_kappa_loc_in = ametr_tg(jr)*delong(jr)/(drmin(jr)*elon_tg(jr))
    tglf_delta_loc_in   = tria_tg(jr)
    tglf_s_delta_loc_in = ametr_tg(jr)*dtrian(jr)/drmin(jr)
    tglf_zeta_loc_in    = 0.
    tglf_s_zeta_loc_in  = 0.
    tglf_q_loc_in       = q_tg(jr)
    tglf_q_prime_loc_in = (q_tg(jr)/rmin_tg)*dq(jr)/dr(jr)
    tglf_p_prime_loc_in = (k0/Bunit_gauss**2)*(q_tg(jr)/rmin_tg)*dptot(jr)/dr(jr)

    tglf_rmin_sa_in     = rmin_tg
    tglf_rmaj_sa_in     = rmaj_tg(jr)/a0_m
    tglf_q_sa_in        = q_tg(jr)
    tglf_shat_sa_in     = (ametr_tg(jr)/q_tg(jr))*dq(jr)/drmin(jr)
    tglf_alpha_sa_in    = -(8.0*pi*k0/Bunit_gauss**2)*q_tg(jr)**2 * rmaj_tg(jr)*dptot(jr)/drmin(jr)
    tglf_xwell_sa_in    = 0.
    tglf_theta0_sa_in   = 0.

! Settings
    if (tglf_dump_flag_in) then
        write(file_dump_local, '(A11, I0)') 'input.tglf_', jr
    endif
! -----------------------
    call tglf_run
! -----------------------

! Transport coefficients

    ion_eflux = SUM(tglf_ion_eflux_out(1: tglf_ns_in-1))
    ion_eflux = ion_eflux/(tglf_taus_in(2) * 1e13*ni_tg(1, jr)/N0)
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
output(:, 1) = chii
output(:, 2) = chie
output(:, 3) = mtori
output(:, 4) = elec_pflux
output(:, 5) = exchi
output(:, 6) = gamma_max
output(:, 7) = omega_max
do jspec=1, tglf_ns_in-1
    output(:, 7+jspec) = ion_pflux(jspec, :)
enddo


call MPI_Send(output, chunk * n_outputs, MPI_DOUBLE_PRECISION, 0, 1, parent, ierr)
call MPI_Finalize(ierr)

end program qlk_chunk
