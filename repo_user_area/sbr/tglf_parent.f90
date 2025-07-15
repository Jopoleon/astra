subroutine tglf_parent(CHI, CHE, VIN, DPH, DPL, DPR, exchi_out, GM1, OM1)

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

integer, parameter :: n_scalars=13, n_inputs=55, n_outputs=15, n_dims=8, nrho_m=80, nspec_max=5
double precision, parameter :: c_vpol=1.d0
double precision, parameter :: &
   k0 = 1.6022E-12, &    ! erg/ev
   mp = 1.6726E-24       ! proton mass (g)

double precision, dimension(*), intent(out) :: CHI, CHE, VIN, DPH, DPL, DPR, exchi_out, GM1, OM1

integer :: ierr, info, intercomm, errcodes(100), status(MPI_STATUS_SIZE)
integer :: jr, jrho, jr_r, jr_l, jgamma_max, jspec
integer :: ns_in              ! Number of species, including electrons
integer :: i, j, k, i1, i2, chunk, nworkers, dims(n_dims)

double precision, dimension(n_scalars) :: scal_in_m
double precision, dimension(n_outputs, nrho_m) :: prof_out_m
double precision :: bmod, bpolz, xstep, rho_min, rho_max, dstep, T0, m0, a0_m, a0_cm, cs0
double precision, dimension(nrho_m) :: drmin, drmaj, drho, dte, dne, dq, dptot, &
    delong, dtrian, dvpar, dvper, drhodr, dr, dv_r
double precision, dimension(NRD) :: gradrhosq_exp, rmaj_exp, q_exp, &
    vexb_exp, vpar_exp, vper_exp, mtori_m, &
    chie_m, chii_m, elec_pflux_m, exchi_m, ptot_exp, gamma_m, omega_m
double precision, dimension(nrho_m) :: mtori, chie, chii, exchi, elec_pflux, rho_m, &
    gamma_max, omega_max, kymax, te_m, ne_m, vpar_m, vper_m, vexb_m, &
    ametr_m, elon_m, tria_m, rmaj_m, ptot_m, q_m, zef_m, pfn_m
double precision, dimension(nspec_max) :: mass_in, zs_in
double precision, dimension(nspec_max-1, nrho_m) :: dti, dni, ni_m, ti_m, zi_m, ion_pflux
double precision, dimension(nspec_max-1, NRD) :: ni_exp, ion_pflux_m
double precision, allocatable, dimension(:, :) :: send_buffer
character(len=256) :: worker_exe

worker_exe = "xpr/tglf.x"

! Interpolate from ASTRA grid to TGLF grid
rho_min = RHO(1)
rho_max = RHO(NA1)
!rho_max = max(RHO(NA1I), RHO(NA1E), RHO(NA1N))
xstep = (rho_max - rho_min)/(nrho_m - 1.)
rho_m = (/ (rho_min + (jr - 1.)*xstep, jr=1, nrho_m) /)

call qinterp(RHO(1:NA1),     TI(1:NA1), NA1, rho_m, ti_m(1, :), nrho_m)
call qinterp(RHO(1:NA1),     TE(1:NA1), NA1, rho_m,       te_m, nrho_m)
call qinterp(RHO(1:NA1),   ZIM1(1:NA1), NA1, rho_m, zi_m(2, :), nrho_m)
call qinterp(RHO(1:NA1),   ZIM2(1:NA1), NA1, rho_m, zi_m(3, :), nrho_m)
call qinterp(RHO(1:NA1),   ZIM3(1:NA1), NA1, rho_m, zi_m(4, :), nrho_m)
call qinterp(RHO(1:NA1),     NE(1:NA1), NA1, rho_m,       ne_m, nrho_m)
call qinterp(RHO(1:NA1),    ZEF(1:NA1), NA1, rho_m,      zef_m, nrho_m)
call qinterp(RHO(1:NA1),  AMETR(1:NA1), NA1, rho_m,    ametr_m, nrho_m)
call qinterp(RHO(1:NA1),   ELON(1:NA1), NA1, rho_m,     elon_m, nrho_m)
call qinterp(RHO(1:NA1),   TRIA(1:NA1), NA1, rho_m,     tria_m, nrho_m)

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

call qinterp(RHO(1:NA1), ni_exp(1, 1:NA1), NA1, rho_m, ni_m(1, :), nrho_m)
call qinterp(RHO(1:NA1), ni_exp(2, 1:NA1), NA1, rho_m, ni_m(2, :), nrho_m)
call qinterp(RHO(1:NA1), ni_exp(3, 1:NA1), NA1, rho_m, ni_m(3, :), nrho_m)
call qinterp(RHO(1:NA1), ni_exp(4, 1:NA1), NA1, rho_m, ni_m(4, :), nrho_m)
call qinterp(RHO(1:NA1),  rmaj_exp(1:NA1), NA1, rho_m,  rmaj_m, nrho_m)
call qinterp(RHO(1:NA1),     q_exp(1:NA1), NA1, rho_m,     q_m, nrho_m)
call qinterp(RHO(1:NA1),  ptot_exp(1:NA1), NA1, rho_m,  ptot_m, nrho_m)
call qinterp(RHO(1:NA1),  vpar_exp(1:NA1), NA1, rho_m,  vpar_m, nrho_m)
call qinterp(RHO(1:NA1),  vper_exp(1:NA1), NA1, rho_m,  vper_m, nrho_m)
call qinterp(RHO(1:NA1),  vexb_exp(1:NA1), NA1, rho_m,  vexb_m, nrho_m)
call qinterp(RHO(1:NA1),   FP_NORM(1:NA1), NA1, rho_m,   pfn_m, nrho_m)

! Reference length
a0_m = AMETR(NA1)

! Species cmassses and charges
mass_in(1) = 5.4447e-4
mass_in(2) = AMJ
mass_in(3) = AIM1
mass_in(4) = AIM2
mass_in(5) = AIM3
do jr=1, nrho_m
    ni_m(2, jr) = max(1.e-9, ni_m(2, jr))
    ni_m(3, jr) = max(1.e-9, ni_m(3, jr))
    ni_m(4, jr) = max(1.e-9, ni_m(4, jr))
enddo

elec_pflux_m = 0.
ion_pflux_m  = 0.
chie_m  = 0.
chii_m  = 0.
mtori_m = 0.
exchi_m = 0.

! Number of species
ns_in = nspec_max

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

ti_m(2, :) = ti_m(1, :)
ti_m(3, :) = ti_m(1, :)
ti_m(4, :) = ti_m(1, :)
zi_m(1, :) = ZMJ

!--------------
! Differentials

do jr=1, nrho_m
    jr_r = jr + 1
    jr_l = jr - 1
    if (jr == 1) then
        jr_l = jr
    else if (jr == nrho_m) then
        jr_r = jr
    endif
    dstep = 1./dble(jr_r - jr_l)  ! 0.5 in between, 1 at the edges
    drmin(jr)  = dstep*(ametr_m(jr_r) - ametr_m(jr_l))
    drmaj(jr)  = dstep*( rmaj_m(jr_r) -  rmaj_m(jr_l))
    drho(jr)   = dstep*(  rho_m(jr_r) -   rho_m(jr_l))
    delong(jr) = dstep*( elon_m(jr_r) -  elon_m(jr_l))
    dtrian(jr) = dstep*( tria_m(jr_r) -  tria_m(jr_l))
    dptot(jr)  = dstep*( ptot_m(jr_r) -  ptot_m(jr_l)) * 1E3*1E13
    dte(jr)    = dstep*(te_m(jr_r) - te_m(jr_l))
    dne(jr)    = dstep*(ne_m(jr_r) - ne_m(jr_l))
    dq(jr)     = dstep*(q_m(jr_r) - q_m(jr_l))
    dvper(jr)  = dstep*(vper_m(jr_r) - vper_m(jr_l))
    do jspec=1, ns_in-1
        dti(jspec, jr) = dstep*(ti_m(jspec, jr_r) - ti_m(jspec, jr_l))
        dni(jspec, jr) = dstep*(ni_m(jspec, jr_r) - ni_m(jspec, jr_l))
    enddo
    dv_r(jr) = dstep* &
        (vpar_m(jr_r)/(rmaj_m(jr_r) + ametr_m(jr_r)) - &
         vpar_m(jr_l)/(rmaj_m(jr_l) + ametr_m(jr_l)))
    dr(jr) = drmin(jr)/a0_m    ! gradients w.r.t. minor radius even for s-alpha geometry
    drhodr(jr) = drho(jr)/drmin(jr)
enddo

!--------------------
! Populate prof_in_m

dims(1) = nrho_m
dims(2) = n_scalars
dims(3) = n_inputs
dims(4) = n_outputs
dims(5) = ns_in

scal_in_m(1) = RTOR
scal_in_m(2) = BTOR
scal_in_m(3) = a0_m
scal_in_m(4:  8) = mass_in(1:5)
scal_in_m(9: 13) = zs_in(1:5)

!--------------
! Send MPI jobs
!--------------

nworkers = 40  ! A submultiple of nrho_m!
chunk = nrho_m / nworkers
print *, "MPI workers = ", nworkers, nrho_m, chunk*n_inputs
allocate(send_buffer(n_inputs, chunk))

call MPI_Info_create(info, ierr)
call MPI_Info_set(info, "host", "localhost", ierr)
call MPI_Info_set(info, "oversubscribe", "false", ierr)
call MPI_Comm_spawn(worker_exe, MPI_ARGV_NULL, nworkers, MPI_INFO_NULL, 0, MPI_COMM_SELF, intercomm, errcodes, ierr)

! Send dimensions and data to workers
do i=0, nworkers-1
    send_buffer = 0.d0
    i1 = i * chunk + 1
    i2 = (i + 1) * chunk
    dims(7) = i1
    dims(8) = i2
    call MPI_Send(dims, n_dims, MPI_INTEGER, i, 100 + i, intercomm, ierr)
    call MPI_Send(scal_in_m, n_scalars, MPI_DOUBLE_PRECISION, i, 101 + i, intercomm, ierr)

    send_buffer( 1, :) = rho_m(i1:i2)
    send_buffer( 2, :) = ametr_m(i1:i2)
    send_buffer( 3, :) = rmaj_m(i1:i2)
    send_buffer( 4, :) = elon_m(i1:i2)
    send_buffer( 5, :) = tria_m(i1:i2)
    send_buffer( 6, :) = q_m(i1:i2) 
    send_buffer( 7, :) = pfn_m(i1:i2)
    send_buffer( 8, :) = ptot_m(i1:i2)
    send_buffer( 9, :) = ne_m(i1:i2) 
    send_buffer(10, :) = te_m(i1:i2) 
    send_buffer(11, :) = zef_m(i1:i2)
    send_buffer(12, :) = vpar_m(i1:i2)
    send_buffer(13, :) = vper_m(i1:i2)
    send_buffer(14, :) = vexb_m(i1:i2)
    send_buffer(15, :) = drmin(i1:i2)
    send_buffer(16, :) = drmaj(i1:i2)
    send_buffer(17, :) = drho(i1:i2) 
    send_buffer(18, :) = delong(i1:i2)
    send_buffer(19, :) = dtrian(i1:i2)
    send_buffer(20, :) = dptot(i1:i2)
    send_buffer(21, :) = dte(i1:i2)
    send_buffer(22, :) = dne(i1:i2)
    send_buffer(23, :) = dq(i1:i2)
    send_buffer(24, :) = dvper(i1:i2)
    send_buffer(25, :) = dv_r(i1:i2)
    send_buffer(26, :) = dr(i1:i2)
    send_buffer(27, :) = drhodr(i1:i2)
    send_buffer(28, :) = ti_m(1, i1:i2)
    send_buffer(29, :) = ti_m(2, i1:i2)
    send_buffer(30, :) = ti_m(3, i1:i2)
    send_buffer(31, :) = ti_m(4, i1:i2)
    send_buffer(32, :) = ni_m(1, i1:i2)
    send_buffer(33, :) = ni_m(2, i1:i2)
    send_buffer(34, :) = ni_m(3, i1:i2)
    send_buffer(35, :) = ni_m(4, i1:i2)
    send_buffer(36, :) = zi_m(1, i1:i2)
    send_buffer(37, :) = zi_m(2, i1:i2)
    send_buffer(38, :) = zi_m(3, i1:i2)
    send_buffer(39, :) = zi_m(4, i1:i2)
    send_buffer(40, :) = dti(1, i1:i2)
    send_buffer(41, :) = dti(2, i1:i2)
    send_buffer(42, :) = dti(3, i1:i2)
    send_buffer(43, :) = dti(4, i1:i2)
    send_buffer(44, :) = dni(1, i1:i2)
    send_buffer(45, :) = dni(2, i1:i2)
    send_buffer(46, :) = dni(3, i1:i2)
    send_buffer(47, :) = dni(4, i1:i2)
    call MPI_Send(send_buffer, chunk * n_inputs, MPI_DOUBLE_PRECISION, i, 102+i, intercomm, ierr)
enddo
call MPI_Barrier(intercomm, ierr)  ! Optional: ensure child finished before next step

! Receive results from each worker
do i=0, nworkers-1
    i1 = i * chunk + 1
    i2 = (i + 1) * chunk
    call MPI_Recv(prof_out_m(:, i1:i2), chunk * n_outputs, MPI_DOUBLE_PRECISION, i, 103+i, intercomm, status, ierr)
enddo
call MPI_Barrier(intercomm, ierr)  ! Optional: ensure child finished before next step

! Interpolate back to ASTRA radial grid

call qinterp(rho_m, prof_out_m(1, :), nrho_m, RHO(1:NA1), chii_m(1:NA1)      , NA1)
call qinterp(rho_m, prof_out_m(2, :), nrho_m, RHO(1:NA1), chie_m(1:NA1)      , NA1)
call qinterp(rho_m, prof_out_m(3, :), nrho_m, RHO(1:NA1), mtori_m(1:NA1)     , NA1)
call qinterp(rho_m, prof_out_m(4, :), nrho_m, RHO(1:NA1), elec_pflux_m(1:NA1), NA1)
call qinterp(rho_m, prof_out_m(5, :), nrho_m, RHO(1:NA1), exchi_m(1:NA1)     , NA1)
call qinterp(rho_m, prof_out_m(6, :), nrho_m, RHO(1:NA1), gamma_m(1:NA1)     , NA1)
call qinterp(rho_m, prof_out_m(7, :), nrho_m, RHO(1:NA1), omega_m(1:NA1)     , NA1)
do jspec=1, nspec_max-1
    call qinterp(rho_m, prof_out_m(7+jspec, :), nrho_m, RHO(1:NA1), ion_pflux_m(jspec, 1:NA1), NA1)
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
    exchi_out(jrho) = exchi_m(jrho) ! turbulent e-i equipartition in MW/m^3
    T0  = 1E3 *TE(jrho)       ! temperature scale used by GYRO
    cs0 = SQRT(k0*T0/m0)      ! thermal velocity unit cm/sec
    GM1(jrho) = gamma_m(jrho)*(cs0/a0_cm)
    OM1(jrho) = omega_m(jrho)*(cs0/a0_cm)
enddo

return
end subroutine tglf_parent
