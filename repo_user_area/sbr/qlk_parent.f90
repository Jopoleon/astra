subroutine qlk_parent(CHI, CHE, VIN, DPH, DPL, DPR, XTB, GM1, OM1)

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

integer, parameter :: n_scalars=20, n_inputs=55, n_outputs=15, nrho_qlk=80, nsm=7, nky_in=19
double precision, parameter :: c_vpol=1.d0
double precision, parameter :: &
   k0 = 1.6022E-12, &    ! erg/ev
   mp = 1.6726E-24       ! proton mass (g)

double precision, dimension(*), intent(out) :: CHI, CHE, VIN, DPH, DPL, DPR, XTB, GM1, OM1

integer :: ierr, intercomm, errcodes(100), status(MPI_STATUS_SIZE)
integer :: jr, jrho, jr_r, jr_l, jgamma_max, jspec
integer :: ns_in              ! Number of species, including electrons
integer :: i, i1, i2, chunk, nprocs, nworkers, dims(8)

double precision, dimension(n_scalars) :: scal_in_qlk
double precision, dimension(nrho_qlk, n_inputs ) :: prof_in_qlk
double precision, dimension(nrho_qlk, n_outputs) :: prof_out_qlk
double precision :: bmod, bpolz, xstep, rho_min, rho_max, dstep, T0, m0, a0_m, a0_cm, cs0
double precision, dimension(nrho_qlk) :: drmin, drmaj, drho, dte, dne, dq, dptot, &
    delong, dtrian, dvper, drhodr, dr, dv_r
double precision, dimension(NRD) :: gradrhosq_exp, rmaj_exp, q_exp, &
    vexb_exp, vpar_exp, vper_exp, mtori_m, &
    chie_m, chii_m, elec_pflux_m, exchi_m, ptot_exp, gamma_m, omega_m
double precision, dimension(nrho_qlk) :: mtori, chie, chii, exchi, elec_pflux, rho_qlk, &
    gamma_max, omega_max, kymax, te_qlk, ne_qlk, vpar_qlk, vper_qlk, vexb_qlk, &
    ametr_qlk, elon_qlk, tria_qlk, rmaj_qlk, ptot_qlk, q_qlk, zef_qlk, pfn_qlk
double precision, dimension(nsm) :: mass_in, zs_in
double precision, dimension(nsm-1, nrho_qlk) :: dti, dni, ni_qlk, ti_qlk, ion_pflux
double precision, dimension(nsm-2, nrho_qlk) :: zimp_qlk 
double precision, dimension(nsm-1, NRD) :: ni_exp, ion_pflux_m
character(len=256) :: worker_exe

worker_exe = "xpr/qlk.x"

! Interpolate from ASTRA grid to QuaLiKiZ grid
rho_min = RHO(1)
rho_max = RHO(NA1)
!rho_max = max(RHO(NA1I), RHO(NA1E), RHO(NA1N))
xstep = (rho_max - rho_min)/(nrho_qlk - 1.)
rho_qlk = (/ (rho_min + (jr - 1.)*xstep, jr=1, nrho_qlk) /)

call qinterp(RHO(1:NA1),     TI(1:NA1), NA1, rho_qlk,   ti_qlk(1, :), nrho_qlk)
call qinterp(RHO(1:NA1),     TE(1:NA1), NA1, rho_qlk,         te_qlk, nrho_qlk)
call qinterp(RHO(1:NA1),   ZIM1(1:NA1), NA1, rho_qlk, zimp_qlk(1, :), nrho_qlk)
call qinterp(RHO(1:NA1),   ZIM2(1:NA1), NA1, rho_qlk, zimp_qlk(2, :), nrho_qlk)
call qinterp(RHO(1:NA1),   ZIM3(1:NA1), NA1, rho_qlk, zimp_qlk(3, :), nrho_qlk)
call qinterp(RHO(1:NA1),     NE(1:NA1), NA1, rho_qlk,         ne_qlk, nrho_qlk)
call qinterp(RHO(1:NA1),    ZEF(1:NA1), NA1, rho_qlk,        zef_qlk, nrho_qlk)
call qinterp(RHO(1:NA1),  AMETR(1:NA1), NA1, rho_qlk,      ametr_qlk, nrho_qlk)
call qinterp(RHO(1:NA1),   ELON(1:NA1), NA1, rho_qlk,       elon_qlk, nrho_qlk)
call qinterp(RHO(1:NA1),   TRIA(1:NA1), NA1, rho_qlk,       tria_qlk, nrho_qlk)

ti_qlk(2, :) = ti_qlk(1, :)
ti_qlk(3, :) = ti_qlk(1, :)
ti_qlk(4, :) = ti_qlk(1, :)

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

call qinterp(RHO(1:NA1), ni_exp(1, 1:NA1), NA1, rho_qlk, ni_qlk(1, :), nrho_qlk)
call qinterp(RHO(1:NA1), ni_exp(2, 1:NA1), NA1, rho_qlk, ni_qlk(2, :), nrho_qlk)
call qinterp(RHO(1:NA1), ni_exp(3, 1:NA1), NA1, rho_qlk, ni_qlk(3, :), nrho_qlk)
call qinterp(RHO(1:NA1), ni_exp(4, 1:NA1), NA1, rho_qlk, ni_qlk(4, :), nrho_qlk)
call qinterp(RHO(1:NA1),  rmaj_exp(1:NA1), NA1, rho_qlk,  rmaj_qlk, nrho_qlk)
call qinterp(RHO(1:NA1),     q_exp(1:NA1), NA1, rho_qlk,     q_qlk, nrho_qlk)
call qinterp(RHO(1:NA1),  ptot_exp(1:NA1), NA1, rho_qlk,  ptot_qlk, nrho_qlk)
call qinterp(RHO(1:NA1),  vpar_exp(1:NA1), NA1, rho_qlk,  vpar_qlk, nrho_qlk)
call qinterp(RHO(1:NA1),  vper_exp(1:NA1), NA1, rho_qlk,  vper_qlk, nrho_qlk)
call qinterp(RHO(1:NA1),  vexb_exp(1:NA1), NA1, rho_qlk,  vexb_qlk, nrho_qlk)
call qinterp(RHO(1:NA1),   FP_NORM(1:NA1), NA1, rho_qlk,   pfn_qlk, nrho_qlk)

! Reference length
a0_m = AMETR(NA1)

! Species cmassses and charges
mass_in(1) = 5.4447e-4
mass_in(2) = AMJ  ! AMJ is reference mass
mass_in(3) = AIM1
mass_in(4) = AIM2
mass_in(5) = AIM3
do jr=1, nrho_qlk
    ni_qlk(2, jr) = max(1.e-9, ni_qlk(2, jr))
    ni_qlk(3, jr) = max(1.e-9, ni_qlk(3, jr))
    ni_qlk(4, jr) = max(1.e-9, ni_qlk(4, jr))
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

do jr=1, nrho_qlk
    jr_r = jr + 1
    jr_l = jr - 1
    if (jr == 1) then
        jr_l = jr
    else if (jr == nrho_qlk) then
        jr_r = jr
    endif
    dstep = 1./dble(jr_r - jr_l)  ! 0.5 in between, 1 at the edges
    drmin(jr)  = dstep*(ametr_qlk(jr_r) - ametr_qlk(jr_l))
    drmaj(jr)  = dstep*( rmaj_qlk(jr_r) -  rmaj_qlk(jr_l))
    drho(jr)   = dstep*(  rho_qlk(jr_r) -   rho_qlk(jr_l))
    delong(jr) = dstep*( elon_qlk(jr_r) -  elon_qlk(jr_l))
    dtrian(jr) = dstep*( tria_qlk(jr_r) -  tria_qlk(jr_l))
    dptot(jr)  = dstep*( ptot_qlk(jr_r) -  ptot_qlk(jr_l))
    dte(jr)    = dstep*(te_qlk(jr_r) - te_qlk(jr_l))
    dne(jr)    = dstep*(ne_qlk(jr_r) - ne_qlk(jr_l))
    dq(jr)     = dstep*(q_qlk(jr_r) - q_qlk(jr_l))
    dvper(jr)  = dstep*(vper_qlk(jr_r) - vper_qlk(jr_l))
    do jspec=1, ns_in-1
        dti(jspec, jr) = dstep*(ti_qlk(jspec, jr_r) - ti_qlk(jspec, jr_l))
        dni(jspec, jr) = dstep*(ni_qlk(jspec, jr_r) - ni_qlk(jspec, jr_l))
    enddo
    dv_r(jr) = dstep* &
        (vpar_qlk(jr_r)/(rmaj_qlk(jr_r) + ametr_qlk(jr_r)) - &
         vpar_qlk(jr_l)/(rmaj_qlk(jr_l) + ametr_qlk(jr_l)))
    dr(jr) = drmin(jr)/a0_m    ! gradients w.r.t. minor radius even for s-alpha geometry
    drhodr(jr) = drho(jr)/drmin(jr)
enddo

!--------------------
! Populate prof_in_qlk

dims(1) = nrho_qlk
dims(2) = n_scalars
dims(3) = n_inputs
dims(4) = n_outputs
dims(5) = ns_in
dims(6) = nky_in

scal_in_qlk(1) = BTOR
scal_in_qlk(2) = a0_m
scal_in_qlk(3:  7) = mass_in(1:5)
scal_in_qlk(8: 12) = zs_in(1:5)

prof_in_qlk(:,  1) = rho_qlk
prof_in_qlk(:,  2) = ametr_qlk
prof_in_qlk(:,  3) = rmaj_qlk
prof_in_qlk(:,  4) = elon_qlk
prof_in_qlk(:,  5) = tria_qlk
prof_in_qlk(:,  6) = q_qlk
prof_in_qlk(:,  7) = pfn_qlk
prof_in_qlk(:,  8) = ptot_qlk
prof_in_qlk(:,  9) = ne_qlk
prof_in_qlk(:, 10) = te_qlk
prof_in_qlk(:, 11) = zef_qlk
prof_in_qlk(:, 12) = vpar_qlk
prof_in_qlk(:, 13) = vper_qlk
prof_in_qlk(:, 14) = vexb_qlk
prof_in_qlk(:, 15) = ti_qlk(1, :)
prof_in_qlk(:, 16) = ti_qlk(2, :)
prof_in_qlk(:, 17) = ti_qlk(3, :)
prof_in_qlk(:, 18) = ti_qlk(4, :)
prof_in_qlk(:, 19) = ni_qlk(1, :)
prof_in_qlk(:, 20) = ni_qlk(2, :)
prof_in_qlk(:, 21) = ni_qlk(3, :)
prof_in_qlk(:, 22) = ni_qlk(4, :)
prof_in_qlk(:, 23) = zimp_qlk(1, :)
prof_in_qlk(:, 24) = zimp_qlk(2, :)
prof_in_qlk(:, 25) = zimp_qlk(3, :)
prof_in_qlk(:, 26) = drmin
prof_in_qlk(:, 27) = drmaj
prof_in_qlk(:, 28) = drho
prof_in_qlk(:, 29) = delong
prof_in_qlk(:, 30) = dtrian
prof_in_qlk(:, 31) = dptot
prof_in_qlk(:, 32) = dte
prof_in_qlk(:, 33) = dne
prof_in_qlk(:, 34) = dq
prof_in_qlk(:, 35) = dvper
prof_in_qlk(:, 36) = dv_r
prof_in_qlk(:, 37) = dr
prof_in_qlk(:, 38) = drhodr
prof_in_qlk(:, 39) = dti(1, :)
prof_in_qlk(:, 40) = dti(2, :)
prof_in_qlk(:, 41) = dti(3, :)
prof_in_qlk(:, 42) = dti(4, :)
prof_in_qlk(:, 43) = dni(1, :)
prof_in_qlk(:, 44) = dni(2, :)
prof_in_qlk(:, 45) = dni(3, :)
prof_in_qlk(:, 46) = dni(4, :)

!--------------
! Send MPI jobs
!--------------

call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

nworkers = 40  ! A submultiple of nrho_qlk!
chunk = nrho_qlk / nworkers

call MPI_Comm_spawn(worker_exe, MPI_ARGV_NULL, nworkers, MPI_INFO_NULL, 0, MPI_COMM_SELF, intercomm, errcodes, ierr)
print *, "MPI workers = ", nworkers, nrho_qlk

! Send dimensions and data to workers
do i=0, nworkers-1
    i1 = i * chunk + 1
    dims(7) = i1
    call MPI_Send(dims, SIZE(dims), MPI_INTEGER, i, 0, intercomm, ierr)
enddo
do i=0, nworkers-1
    call MPI_Send(scal_in_qlk, n_scalars, MPI_DOUBLE_PRECISION, i, 0, intercomm, ierr)
enddo
do i=0, nworkers-1
    i1 = i * chunk + 1
    i2 = (i + 1) * chunk
    call MPI_Send(prof_in_qlk(i1:i2, :), chunk * n_inputs, MPI_DOUBLE_PRECISION, i, 0, intercomm, ierr)
enddo

! Receive results from each worker
do i=0, nworkers-1
    i1 = i * chunk + 1
    i2 = (i + 1) * chunk
    call MPI_Recv(prof_out_qlk(i1:i2, :), chunk * n_outputs, MPI_DOUBLE_PRECISION, i, 1, intercomm, status, ierr)
enddo

! Interpolate back to ASTRA radial grid

call qinterp(rho_qlk, prof_out_qlk(:, 1), nrho_qlk, RHO(1:NA1), chii_m(1:NA1)      , NA1)
call qinterp(rho_qlk, prof_out_qlk(:, 2), nrho_qlk, RHO(1:NA1), chie_m(1:NA1)      , NA1)
call qinterp(rho_qlk, prof_out_qlk(:, 3), nrho_qlk, RHO(1:NA1), mtori_m(1:NA1)     , NA1)
call qinterp(rho_qlk, prof_out_qlk(:, 4), nrho_qlk, RHO(1:NA1), elec_pflux_m(1:NA1), NA1)
call qinterp(rho_qlk, prof_out_qlk(:, 5), nrho_qlk, RHO(1:NA1), exchi_m(1:NA1)     , NA1)
call qinterp(rho_qlk, prof_out_qlk(:, 6), gamma_max , nrho_qlk, RHO(1:NA1), gamma_m(1:NA1)     , NA1)
call qinterp(rho_qlk, prof_out_qlk(:, 7), nrho_qlk, RHO(1:NA1), omega_m(1:NA1)     , NA1)
do jspec=1, ns_in-1
    call qinterp(rho_qlk, prof_out_qlk(7+jspec, 1:nrho_qlk), nrho_qlk, RHO(1:NA1), ion_pflux_m(jspec, 1:NA1), NA1)
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
end subroutine qlk_parent
