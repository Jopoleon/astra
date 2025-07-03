subroutine tglf_parent

use mpi

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

implicit none

integer, parameter :: n_scalars=10, n_inputs=55, n_outputs=15, nrho_tg=60, nsm=7, nky_in=19
double precision, parameter :: c_vpol=1.d0

integer :: ierr, intercomm, errcodes(100), status(MPI_STATUS_SIZE)
integer :: jr, jrho, jr_r, jr_l, jgamma_max, jspec, kyloop, mom_order
integer :: sat_rule           ! Saturation rule
integer :: geom_flag          ! 1: Miller; 2: Fourier; 3: ELITE
integer :: nmodes_tg          ! number of unstable modes to use in computing fluxes (max=4)
integer :: kygrid_model_tg    ! select version of ky-grid to use 1
integer :: xnu_model_tg       ! select version of trapped-passing 2
integer :: ns_in              ! Number of species, including electrons
integer :: i, i1, i2, chunk, nprocs, nworkers, dims(6)

double precision, dimension(n_scalars) :: scal_in_tg
double precision, dimension(nrho_tg, n_inputs ) :: prof_in_tg
double precision, dimension(nrho_tg, n_outputs) :: prof_out_tg
double precision :: bmod, bpolz, alpha_zf_in, ion_eflux, ion_mflux, xstep, rho_min, rho_max, dstep
double precision, dimension(nrho_tg) :: drmin, drmaj, drho, dte, dne, dq, dptot, &
    delong, dtrian, dvper, drhodr, dr, dv_r
double precision :: Bunit_gauss, Bunit_T, cs0, cs00, rhos0, omega0, rhostar2, lnlamda, taue, cexb
double precision :: a0_cm, a0_m, T0, N0, m0, rmin_tg, drho_cs, drho_nt, nt_cs
double precision :: wdia_trap_tg          ! parameter for trapped fraction model
double precision, dimension(NRD) :: gradrhosq_exp, rmaj_exp, q_exp, &
    vexb_exp, vpar_exp, vper_exp, mtori_m, &
    chie_m, chii_m, elec_pflux_m, exchi_m, ptot_exp, gamma_m, omega_m
double precision, dimension(nrho_tg) :: mtori, chie, chii, exchi, elec_pflux, rho_tg, &
    gamma_max, omega_max, kymax, te_tg, ne_tg, vpar_tg, vper_tg, vexb_tg, &
    ametr_tg, elon_tg, tria_tg, rmaj_tg, ptot_tg, q_tg, zef_tg, pfn_tg
double precision, dimension(nsm) :: mass_in, zs_in
double precision, dimension(nsm-1, nrho_tg) :: dti, dni, ni_tg, ti_tg, z_tg, ion_pflux
double precision, dimension(nsm-1, NRD) :: ni_exp, ion_pflux_m
double precision, dimension(nky_in) :: gamma, omega, kyspectrum, efluxspectrum, &
    ifluxspectrum, pfluxspectrum
character(len=256) :: worker_exe

worker_exe = "/shares/departments/AUG/users/git/a8/xpr/tglf.x"

! Interpolate from ASTRA grid to TGLF grid
rho_min = RHO(1)
!rho_max = RHO(NA1)
rho_max = max(RHO(NA1I), RHO(NA1E), RHO(NA1N))
xstep = (rho_max - rho_min)/(nrho_tg - 1.)
rho_tg = (/ (rho_min + (jr - 1.)*xstep, jr=1, nrho_tg) /)

z_tg(1, :) = 1.
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
call qinterp(RHO(1:NA1),  FP_NORM(1:NA1), NA1, rho_tg,   pfn_tg, nrho_tg)

! Reference length
a0_m = AMETR(NA1)

! Species cmassses and charges
mass_in(1) = 5.4447e-4/AMJ
mass_in(2) = AMJ/AMJ  ! AMJ is reference mass
mass_in(3) = AIM1/AMJ
mass_in(4) = AIM2/AMJ
mass_in(5) = AIM3/AMJ
do jr=1, nrho_tg
    ni_tg(2, jr) = max(1.e-9, ni_tg(2, jr))
    ni_tg(3, jr) = max(1.e-9, ni_tg(3, jr))
    ni_tg(4, jr) = max(1.e-9, ni_tg(4, jr))
enddo

elec_pflux_m = 0.
ion_pflux_m  = 0.
chie_m  = 0.
chii_m  = 0.
mtori_m = 0.
exchi_m = 0.

! Number of species
ns_in = nsm

! These will be reset locally in the radial loop
zs_in(1) = -1.
zs_in(2) = ZMJ
zs_in(3) = MAXVAL(ZIM1(1:NA1))
zs_in(4) = MAXVAL(ZIM2(1:NA1))
zs_in(5) = MAXVAL(ZIM3(1:NA1))

if (zs_in(5) >= 1.) then
    ns_in = 5
else
    ns_in = 4
endif
if (zs_in(4) < 1.) ns_in = 3
if (zs_in(5) >= 1. .and. ns_in == 3) then
    ns_in = 4
endif
if (zs_in(3) < 1.) ns_in = 2
if (zs_in(4) >= 1. .and. ns_in == 2) then
    ns_in = 3
endif

!--------------
! Differentials

do jr=1, nrho_tg
    jr_r = jr + 1
    jr_l = jr - 1
    if (jr == 1) then
        jr_l = jr
    else if (jr == nrho_tg) then
        jr_r = jr
    endif
    dstep = 1./dble(jr_r - jr_l)  ! 0.5 in between, 1 at the edges
    drmin(jr)  = dstep*(ametr_tg(jr_r) - ametr_tg(jr_l))
    drmaj(jr)  = dstep*( rmaj_tg(jr_r) -  rmaj_tg(jr_l))
    drho(jr)   = dstep*(  rho_tg(jr_r) -   rho_tg(jr_l))
    delong(jr) = dstep*( elon_tg(jr_r) -  elon_tg(jr_l))
    dtrian(jr) = dstep*( tria_tg(jr_r) -  tria_tg(jr_l))
    dptot(jr)  = dstep*( ptot_tg(jr_r) -  ptot_tg(jr_l)) * 1E3*1E13
    dte(jr)    = dstep*(te_tg(jr_r) - te_tg(jr_l))
    dne(jr)    = dstep*(ne_tg(jr_r) - ne_tg(jr_l))
    dq(jr)     = dstep*(q_tg(jr_r) - q_tg(jr_l))
    dvper(jr)  = dstep*(vper_tg(jr_r) - vper_tg(jr_l))
    do jspec=1, ns_in-1
        dti(jspec, jr) = dstep*(ti_tg(jspec, jr_r) - ti_tg(jspec, jr_l))
        dni(jspec, jr) = dstep*(ni_tg(jspec, jr_r) - ni_tg(jspec, jr_l))
    enddo
    dv_r(jr) = dstep* &
        (vpar_tg(jr_r)/(rmaj_tg(jr_r) + ametr_tg(jr_r)) - &
         vpar_tg(jr_l)/(rmaj_tg(jr_l) + ametr_tg(jr_l)))
    dr(jr) = drmin(jr)/a0_m    ! gradients w.r.t. minor radius even for s-alpha geometry
    drhodr(jr) = drho(jr)/drmin(jr)
enddo

!--------------------
! Populate prof_in_tg

scal_in_tg(1) = AMJ
scal_in_tg(2) = BTOR
scal_in_tg(3) = a0_m

prof_in_tg(:,  1) = rho_tg
prof_in_tg(:,  2) = ametr_tg
prof_in_tg(:,  3) = rmaj_tg
prof_in_tg(:,  4) = elon_tg
prof_in_tg(:,  5) = tria_tg
prof_in_tg(:,  6) = q_tg
prof_in_tg(:,  7) = pfn_tg
prof_in_tg(:,  8) = ptot_tg
prof_in_tg(:,  9) = ne_tg
prof_in_tg(:, 10) = te_tg
prof_in_tg(:, 11) = zef_tg
prof_in_tg(:, 12) = vpar_tg
prof_in_tg(:, 13) = vper_tg
prof_in_tg(:, 14) = vexb_tg
prof_in_tg(:, 15) = ti_tg(1, :)
prof_in_tg(:, 16) = ti_tg(2, :)
prof_in_tg(:, 17) = ti_tg(3, :)
prof_in_tg(:, 18) = ti_tg(4, :)
prof_in_tg(:, 19) = ni_tg(1, :)
prof_in_tg(:, 20) = ni_tg(2, :)
prof_in_tg(:, 21) = ni_tg(3, :)
prof_in_tg(:, 22) = ni_tg(4, :)
prof_in_tg(:, 23) = z_tg(1, :)
prof_in_tg(:, 24) = z_tg(2, :)
prof_in_tg(:, 25) = z_tg(3, :)
prof_in_tg(:, 26) = z_tg(4, :)
prof_in_tg(:, 27) = drmin
prof_in_tg(:, 28) = drmaj
prof_in_tg(:, 29) = drho
prof_in_tg(:, 30) = delong
prof_in_tg(:, 31) = dtrian
prof_in_tg(:, 32) = dptot
prof_in_tg(:, 33) = dte
prof_in_tg(:, 34) = dne
prof_in_tg(:, 35) = dq
prof_in_tg(:, 36) = dvper
prof_in_tg(:, 37) = dv_r
prof_in_tg(:, 38) = dr
prof_in_tg(:, 39) = drhodr
prof_in_tg(:, 40) = dti(1, :)
prof_in_tg(:, 41) = dti(2, :)
prof_in_tg(:, 42) = dti(3, :)
prof_in_tg(:, 43) = dti(4, :)
prof_in_tg(:, 44) = dni(1, :)
prof_in_tg(:, 45) = dni(2, :)
prof_in_tg(:, 46) = dni(3, :)
prof_in_tg(:, 47) = dni(4, :)


!--------------
! Send MPI jobs
!--------------

call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

nworkers = 15  ! A submultiple of nrho_tg!
chunk = nrho_tg / nworkers

call MPI_Comm_spawn(worker_exe, MPI_ARGV_NULL, nworkers, MPI_INFO_NULL, 0, MPI_COMM_SELF, intercomm, errcodes, ierr)
print *, "MPI workers = ", nworkers, nrho_tg

dims(1) = nrho_tg
dims(2) = n_scalars
dims(3) = n_inputs
dims(4) = n_outputs
dims(5) = ns_in
dims(6) = nky_in

! Send dimensions and data to workers
do i=0, nworkers-1
    call MPI_Send(dims, 6, MPI_INTEGER, i, 0, intercomm, ierr)
enddo
do i=0, nworkers-1
    call MPI_Send(scal_in_tg, n_scalars, MPI_DOUBLE_PRECISION, i, 0, intercomm, ierr)
enddo
do i=0, nworkers-1
    i1 = i * chunk + 1
    i2 = (i + 1) * chunk
    call MPI_Send(prof_in_tg(i1:i2, :), chunk * n_inputs, MPI_DOUBLE_PRECISION, i, 0, intercomm, ierr)
enddo

! Receive results from each worker
do i=0, nworkers-1
    i1 = i * chunk + 1
    i2 = (i + 1) * chunk
    call MPI_Recv(prof_out_tg(i1:i2, :), chunk * n_outputs, MPI_DOUBLE_PRECISION, i, 1, intercomm, status, ierr)
enddo

write(*, *) 'Out', prof_out_tg(:, 1)

return
end subroutine tglf_parent
