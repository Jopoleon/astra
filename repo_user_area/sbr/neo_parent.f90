subroutine neo_parent(CHI, CHE, VIN, DPH, DPL, DPR, XTB, GM1, OM1)

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

integer, parameter :: n_scalars=20, n_inputs=55, n_outputs=15, nrho_neo=80, nsm=7, nky_in=19
double precision, parameter :: c_vpol=1.d0
double precision, parameter :: &
   k0 = 1.6022E-12, &    ! erg/ev
   mp = 1.6726E-24       ! proton mass (g)

double precision, dimension(*), intent(out) :: CHI, CHE, VIN, DPH, DPL, DPR, XTB, GM1, OM1

integer :: ierr, intercomm, errcodes(100), status(MPI_STATUS_SIZE)
integer :: jr, jrho, jr_r, jr_l, jgamma_max, jspec
integer :: ns_in              ! Number of species, including electrons
integer :: i, i1, i2, chunk, nprocs, nworkers, dims(8)

double precision, dimension(n_scalars) :: scal_in_neo
double precision, dimension(nrho_neo, n_inputs ) :: prof_in_neo
double precision, dimension(nrho_neo, n_outputs) :: prof_out_neo
double precision :: bmod, bpolz, xstep, rho_min, rho_max, dstep, T0, m0, a0_m, a0_cm, cs0
double precision, dimension(nrho_neo) :: drmin, drmaj, drho, dte, dne, dq, dptot, &
    delong, dtrian, dvpar, drhodr, dr
double precision, dimension(NRD) :: gradrhosq_exp, rmaj_exp, q_exp, &
    vexb_exp, vpar_exp, vper_exp, mtori_m, &
    chie_m, chii_m, elec_pflux_m, exchi_m, ptot_exp, gamma_m, omega_m
double precision, dimension(nrho_neo) :: mtori, chie, chii, exchi, elec_pflux, rho_neo, &
    gamma_max, omega_max, kymax, te_neo, ne_neo, vpar_neo, vper_neo, vexb_neo, &
    ametr_neo, elon_neo, tria_neo, rmaj_neo, ptot_neo, q_neo, zef_neo, pfn_neo
double precision, dimension(nsm) :: mass_in, zs_in
double precision, dimension(nsm-1, nrho_neo) :: dti, dni, ni_neo, ti_neo, ion_pflux
double precision, dimension(nsm-2, nrho_neo) :: zimp_neo 
double precision, dimension(nsm-1, NRD) :: ni_exp, ion_pflux_m
character(len=256) :: worker_exe

worker_exe = "xpr/neo.x"

! Interpolate from ASTRA grid to TGLF grid
rho_min = RHO(1)
rho_max = RHO(NA1)
!rho_max = max(RHO(NA1I), RHO(NA1E), RHO(NA1N))
xstep = (rho_max - rho_min)/(nrho_neo - 1.)
rho_neo = (/ (rho_min + (jr - 1.)*xstep, jr=1, nrho_neo) /)

call qinterp(RHO(1:NA1),     TI(1:NA1), NA1, rho_neo,   ti_neo(1, :), nrho_neo)
call qinterp(RHO(1:NA1),     TE(1:NA1), NA1, rho_neo,         te_neo, nrho_neo)
call qinterp(RHO(1:NA1),   ZIM1(1:NA1), NA1, rho_neo, zimp_neo(1, :), nrho_neo)
call qinterp(RHO(1:NA1),   ZIM2(1:NA1), NA1, rho_neo, zimp_neo(2, :), nrho_neo)
call qinterp(RHO(1:NA1),   ZIM3(1:NA1), NA1, rho_neo, zimp_neo(3, :), nrho_neo)
call qinterp(RHO(1:NA1),     NE(1:NA1), NA1, rho_neo,         ne_neo, nrho_neo)
call qinterp(RHO(1:NA1),    ZEF(1:NA1), NA1, rho_neo,        zef_neo, nrho_neo)
call qinterp(RHO(1:NA1),  AMETR(1:NA1), NA1, rho_neo,      ametr_neo, nrho_neo)
call qinterp(RHO(1:NA1),   ELON(1:NA1), NA1, rho_neo,       elon_neo, nrho_neo)
call qinterp(RHO(1:NA1),   TRIA(1:NA1), NA1, rho_neo,       tria_neo, nrho_neo)

ti_neo(2, :) = ti_neo(1, :)
ti_neo(3, :) = ti_neo(1, :)
ti_neo(4, :) = ti_neo(1, :)

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
    vpar_exp(jrho) = VTOR(jrho) * BTOR/bmod + c_vpol* VPOL(jrho) * bpolz/bmod
    vexb_exp(jrho) = -ER(jrho)/bmod ! vexb in m/s (vperp = vexb since the diamagnetic velocity is the curvature drift ac
enddo

call qinterp(RHO(1:NA1), ni_exp(1, 1:NA1), NA1, rho_neo, ni_neo(1, :), nrho_neo)
call qinterp(RHO(1:NA1), ni_exp(2, 1:NA1), NA1, rho_neo, ni_neo(2, :), nrho_neo)
call qinterp(RHO(1:NA1), ni_exp(3, 1:NA1), NA1, rho_neo, ni_neo(3, :), nrho_neo)
call qinterp(RHO(1:NA1), ni_exp(4, 1:NA1), NA1, rho_neo, ni_neo(4, :), nrho_neo)
call qinterp(RHO(1:NA1),  rmaj_exp(1:NA1), NA1, rho_neo,  rmaj_neo, nrho_neo)
call qinterp(RHO(1:NA1),     q_exp(1:NA1), NA1, rho_neo,     q_neo, nrho_neo)
call qinterp(RHO(1:NA1),  ptot_exp(1:NA1), NA1, rho_neo,  ptot_neo, nrho_neo)
call qinterp(RHO(1:NA1),  vpar_exp(1:NA1), NA1, rho_neo,  vpar_neo, nrho_neo)
call qinterp(RHO(1:NA1),  vper_exp(1:NA1), NA1, rho_neo,  vper_neo, nrho_neo)
call qinterp(RHO(1:NA1),  vexb_exp(1:NA1), NA1, rho_neo,  vexb_neo, nrho_neo)
call qinterp(RHO(1:NA1),   FP_NORM(1:NA1), NA1, rho_neo,   pfn_neo, nrho_neo)

! Reference length
a0_m = AMETR(NA1)

! Species cmassses and charges
mass_in(1) = 5.4447e-4
mass_in(2) = AMJ
mass_in(3) = AIM1
mass_in(4) = AIM2
mass_in(5) = AIM3
do jr=1, nrho_neo
    ni_neo(2, jr) = max(1.e-9, ni_neo(2, jr))
    ni_neo(3, jr) = max(1.e-9, ni_neo(3, jr))
    ni_neo(4, jr) = max(1.e-9, ni_neo(4, jr))
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

do jr=1, nrho_neo
    jr_r = jr + 1
    jr_l = jr - 1
    if (jr == 1) then
        jr_l = jr
    else if (jr == nrho_neo) then
        jr_r = jr
    endif
    dstep = 1./dble(jr_r - jr_l)  ! 0.5 in between, 1 at the edges
    drmin(jr)  = dstep*(ametr_neo(jr_r) - ametr_neo(jr_l))
    drmaj(jr)  = dstep*( rmaj_neo(jr_r) -  rmaj_neo(jr_l))
    drho(jr)   = dstep*(  rho_neo(jr_r) -   rho_neo(jr_l))
    delong(jr) = dstep*( elon_neo(jr_r) -  elon_neo(jr_l))
    dtrian(jr) = dstep*( tria_neo(jr_r) -  tria_neo(jr_l))
    dptot(jr)  = dstep*( ptot_neo(jr_r) -  ptot_neo(jr_l)) * 1E3*1E13
    dte(jr)    = dstep*(te_neo(jr_r) - te_neo(jr_l))
    dne(jr)    = dstep*(ne_neo(jr_r) - ne_neo(jr_l))
    dq(jr)     = dstep*(q_neo(jr_r) - q_neo(jr_l))
    dvpar(jr)  = dstep*(vpar_neo(jr_r) - vpar_neo(jr_l))
    do jspec=1, ns_in-1
        dti(jspec, jr) = dstep*(ti_neo(jspec, jr_r) - ti_neo(jspec, jr_l))
        dni(jspec, jr) = dstep*(ni_neo(jspec, jr_r) - ni_neo(jspec, jr_l))
    enddo
    dr(jr) = drmin(jr)/a0_m    ! gradients w.r.t. minor radius even for s-alpha geometry
    drhodr(jr) = drho(jr)/drmin(jr)
enddo

!--------------------
! Populate prof_in_neo

dims(1) = nrho_neo
dims(2) = n_scalars
dims(3) = n_inputs
dims(4) = n_outputs
dims(5) = ns_in
dims(6) = nky_in

scal_in_neo(1) = RTOR
scal_in_neo(2) = BTOR
scal_in_neo(3) = a0_m
scal_in_neo(4:  8) = mass_in(1:5)
scal_in_neo(9: 13) = zs_in(1:5)

prof_in_neo(:,  1) = rho_neo
prof_in_neo(:,  2) = ametr_neo
prof_in_neo(:,  3) = rmaj_neo
prof_in_neo(:,  4) = elon_neo
prof_in_neo(:,  5) = tria_neo
prof_in_neo(:,  6) = q_neo
prof_in_neo(:,  7) = pfn_neo
prof_in_neo(:,  8) = ptot_neo
prof_in_neo(:,  9) = ne_neo
prof_in_neo(:, 10) = te_neo
prof_in_neo(:, 11) = zef_neo
prof_in_neo(:, 12) = vpar_neo
prof_in_neo(:, 13) = vper_neo
prof_in_neo(:, 14) = vexb_neo
prof_in_neo(:, 15) = ti_neo(1, :)
prof_in_neo(:, 16) = ti_neo(2, :)
prof_in_neo(:, 17) = ti_neo(3, :)
prof_in_neo(:, 18) = ti_neo(4, :)
prof_in_neo(:, 19) = ni_neo(1, :)
prof_in_neo(:, 20) = ni_neo(2, :)
prof_in_neo(:, 21) = ni_neo(3, :)
prof_in_neo(:, 22) = ni_neo(4, :)
prof_in_neo(:, 23) = zimp_neo(1, :)
prof_in_neo(:, 24) = zimp_neo(2, :)
prof_in_neo(:, 25) = zimp_neo(3, :)
prof_in_neo(:, 26) = drmin
prof_in_neo(:, 27) = drmaj
prof_in_neo(:, 28) = drho
prof_in_neo(:, 29) = delong
prof_in_neo(:, 30) = dtrian
prof_in_neo(:, 31) = dptot
prof_in_neo(:, 32) = dte
prof_in_neo(:, 33) = dne
prof_in_neo(:, 34) = dq
prof_in_neo(:, 35) = dvpar

prof_in_neo(:, 38) = dr
prof_in_neo(:, 39) = drhodr
prof_in_neo(:, 40) = dti(1, :)
prof_in_neo(:, 41) = dti(2, :)
prof_in_neo(:, 42) = dti(3, :)
prof_in_neo(:, 43) = dti(4, :)
prof_in_neo(:, 44) = dni(1, :)
prof_in_neo(:, 45) = dni(2, :)
prof_in_neo(:, 46) = dni(3, :)
prof_in_neo(:, 47) = dni(4, :)

!--------------
! Send MPI jobs
!--------------

call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

nworkers = 40  ! A submultiple of nrho_neo!
chunk = nrho_neo / nworkers

call MPI_Comm_spawn(worker_exe, MPI_ARGV_NULL, nworkers, MPI_INFO_NULL, 0, MPI_COMM_SELF, intercomm, errcodes, ierr)
print *, "MPI workers = ", nworkers, nrho_neo

! Send dimensions and data to workers
do i=0, nworkers-1
    i1 = i * chunk + 1
    dims(7) = i1
    call MPI_Send(dims, SIZE(dims), MPI_INTEGER, i, 0, intercomm, ierr)
enddo
do i=0, nworkers-1
    call MPI_Send(scal_in_neo, n_scalars, MPI_DOUBLE_PRECISION, i, 0, intercomm, ierr)
enddo
do i=0, nworkers-1
    i1 = i * chunk + 1
    i2 = (i + 1) * chunk
    call MPI_Send(prof_in_neo(i1:i2, :), chunk * n_inputs, MPI_DOUBLE_PRECISION, i, 0, intercomm, ierr)
enddo

! Receive results from each worker
do i=0, nworkers-1
    i1 = i * chunk + 1
    i2 = (i + 1) * chunk
    call MPI_Recv(prof_out_neo(i1:i2, :), chunk * n_outputs, MPI_DOUBLE_PRECISION, i, 1, intercomm, status, ierr)
enddo

! Interpolate back to ASTRA radial grid

call qinterp(rho_neo, prof_out_neo(:, 1), nrho_neo, RHO(1:NA1), chii_m(1:NA1)      , NA1)
call qinterp(rho_neo, prof_out_neo(:, 2), nrho_neo, RHO(1:NA1), chie_m(1:NA1)      , NA1)
call qinterp(rho_neo, prof_out_neo(:, 3), nrho_neo, RHO(1:NA1), mtori_m(1:NA1)     , NA1)
call qinterp(rho_neo, prof_out_neo(:, 4), nrho_neo, RHO(1:NA1), elec_pflux_m(1:NA1), NA1)
call qinterp(rho_neo, prof_out_neo(:, 5), nrho_neo, RHO(1:NA1), exchi_m(1:NA1)     , NA1)
call qinterp(rho_neo, prof_out_neo(:, 6), gamma_max , nrho_neo, RHO(1:NA1), gamma_m(1:NA1)     , NA1)
call qinterp(rho_neo, prof_out_neo(:, 7), nrho_neo, RHO(1:NA1), omega_m(1:NA1)     , NA1)
do jspec=1, ns_in-1
    call qinterp(rho_neo, prof_out_neo(7+jspec, 1:nrho_neo), nrho_neo, RHO(1:NA1), ion_pflux_m(jspec, 1:NA1), NA1)
enddo

chii_m (1:2) = chii_m (3)
chie_m (1:2) = chie_m (3)
mtori_m(1:2) = mtori_m(3)
elec_pflux_m(1:2) = elec_pflux_m(3)
exchi_m(1:2) = exchi_m(3)
omega_m(1:2) = omega_m(3)
gamma_m(1:2) = gamma_m(3)

m0 = AMJ*mp          ! Ref. mass = D ion mass [g]
a0_cm = 1.d2*a0_m    ! length scale used by GYRO, m -> cm

do jrho=1, NA1
    CHI(jrho) = chii_m(jrho)/gradrhosq_exp(jrho) ! \chi_i, m^2/s
    CHE(jrho) = chie_m(jrho)/gradrhosq_exp(jrho) ! \chi_e, m^2/s
    VIN(jrho) = elec_pflux_m(jrho)/a0_m/gradrhosq_exp(jrho) ! D flux
    DPR(jrho) = mtori_m(jrho)
!First impurity only, index 2 of ion species
    if (ns_in >= 3) then
        DPL(jrho) = ion_pflux_m(2, jrho)/a0_m/gradrhosq_exp(jrho)/(ni_exp(2, jrho)/NE(jrho))  ! 1st imp convection
    endif
    if (ns_in >= 4) then
        DPH(jrho) = ion_pflux_m(3, jrho)/a0_m/gradrhosq_exp(jrho)/(ni_exp(3, jrho)/NE(jrho))  ! 2nd imp convection
    endif
    XTB(jrho) = exchi_m(jrho) ! turbulent e-i equipartition in MW/m^3
    T0  = 1E3 *TE(jrho)       ! temperature scale used by GYRO
    cs0 = SQRT(k0*T0/m0)      ! thermal velocity unit cm/sec
    GM1(jrho) = gamma_m(jrho)*(cs0/a0_cm)
    OM1(jrho) = omega_m(jrho)*(cs0/a0_cm)
enddo

return
end subroutine neo_parent
