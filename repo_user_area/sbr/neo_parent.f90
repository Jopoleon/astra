subroutine neo_parent(chii_m, chie_m)

use mpi

use parameter_inc, only: NRD
use const_inc, only: BTOR, RTOR, AMJ, AIM1, AIM2, AIM3, ZMJ, NA1
use status_inc, only: NE, TE, NI, TI, ER, MU, FP_NORM, &
    ZIM1, ZIM2, ZIM3, NDEUT, NIZ1, NIZ2, NIZ3, &
    RHO, AMETR, SHIF, ELON, TRIA, VTOR, VPOL

implicit none

integer, parameter :: n_dims=8, n_scalars=20, n_inputs=32, n_outputs=8, nrho_m=80, nsm=7
double precision, parameter :: c_vpol=1.d0

double precision, dimension(NRD), intent(out) :: chii_m, chie_m

integer :: ierr, intercomm, errcodes(100), status(MPI_STATUS_SIZE)
integer :: jr, jrho, jr_r, jr_l, jion
integer :: ns_in              ! Number of species, including electrons
integer :: i, i1, i2, chunk, nprocs, nworkers, dims(8)

double precision, dimension(n_scalars) :: scal_in_m
double precision, dimension(nrho_m, n_inputs ) :: prof_in_m
double precision, dimension(nrho_m, n_outputs) :: prof_out_m
double precision :: bmod, bpolz, xstep, rho_min, rho_max, dstep, T0, m0, a0_m, a0_cm, cs0, drho
double precision, dimension(nrho_m) :: drmin, drmaj, dti, dte, dne, dq, &
    delong, dtrian, dvpar, drhodr, dr
double precision, dimension(NRD) :: rmaj_as, q_as, ni_main_as, vpar_as, &
    vippd_m, vittd_m, vippi1_m, vitti1_m, j_boot, elec_pflux_m
double precision, dimension(nrho_m) :: rho_m,  ti_m, te_m, ne_m, vpar_m, &
    ametr_m, elon_m, tria_m, rmaj_m, q_m
double precision, dimension(nsm) :: mass_in, zs_in
double precision, dimension(nsm-1, nrho_m) :: dni, ni_m
double precision, dimension(nsm-2, nrho_m) :: zimp_m 
character(len=256) :: worker_exe

worker_exe = "xpr/neo.x"

! Interpolate from ASTRA grid to TGLF grid
rho_min = RHO(1)
rho_max = RHO(NA1)
xstep = (rho_max - rho_min)/(nrho_m - 1.)
rho_m = (/ (rho_min + (jr - 1.)*xstep, jr=1, nrho_m) /)

call qinterp(RHO(1:NA1),  ZIM1(1:NA1), NA1, rho_m, zimp_m(1, :), nrho_m)
call qinterp(RHO(1:NA1),  ZIM2(1:NA1), NA1, rho_m, zimp_m(2, :), nrho_m)
call qinterp(RHO(1:NA1),  ZIM3(1:NA1), NA1, rho_m, zimp_m(3, :), nrho_m)
call qinterp(RHO(1:NA1),  NIZ1(1:NA1), NA1, rho_m,   ni_m(2, :), nrho_m)
call qinterp(RHO(1:NA1),  NIZ2(1:NA1), NA1, rho_m,   ni_m(3, :), nrho_m)
call qinterp(RHO(1:NA1),  NIZ3(1:NA1), NA1, rho_m,   ni_m(4, :), nrho_m)
call qinterp(RHO(1:NA1),    TI(1:NA1), NA1, rho_m,         ti_m, nrho_m)
call qinterp(RHO(1:NA1),    TE(1:NA1), NA1, rho_m,         te_m, nrho_m)
call qinterp(RHO(1:NA1),    NE(1:NA1), NA1, rho_m,         ne_m, nrho_m)
call qinterp(RHO(1:NA1), AMETR(1:NA1), NA1, rho_m,      ametr_m, nrho_m)
call qinterp(RHO(1:NA1),  ELON(1:NA1), NA1, rho_m,       elon_m, nrho_m)
call qinterp(RHO(1:NA1),  TRIA(1:NA1), NA1, rho_m,       tria_m, nrho_m)

do jrho=1, NA1
    if (NDEUT(jrho) >= 0.01*NE(jrho)) then
        ni_main_as(jrho) = NDEUT(jrho)
    else ! likely: NDEUT not defined in equ file, hence zero
        ni_main_as(jrho) = NI(jrho)
    endif
    rmaj_as(jrho) = RTOR + SHIF(jrho)
    q_as(jrho)    = 1./MU(jrho)
    bpolz = BTOR*AMETR(jrho)*MU(jrho)/RTOR
    bmod = sqrt(BTOR**2 + bpolz**2)
    vpar_as(jrho) = VTOR(jrho) * BTOR/bmod + c_vpol* VPOL(jrho) * bpolz/bmod
enddo

call qinterp(RHO(1:NA1), ni_main_as(1:NA1), NA1, rho_m, ni_m(1, :), nrho_m)
call qinterp(RHO(1:NA1), rmaj_as(1:NA1), NA1, rho_m, rmaj_m, nrho_m)
call qinterp(RHO(1:NA1),    q_as(1:NA1), NA1, rho_m,    q_m, nrho_m)
call qinterp(RHO(1:NA1), vpar_as(1:NA1), NA1, rho_m, vpar_m, nrho_m)

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
chie_m  = 0.
chii_m  = 0.

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
    drho       = dstep*(rho(jr_r) - rho(jr_l))
    delong(jr) = dstep*( elon_m(jr_r) -  elon_m(jr_l))
    dtrian(jr) = dstep*( tria_m(jr_r) -  tria_m(jr_l))
    dti(jr)    = dstep*(ti_m(jr_r) - ti_m(jr_l))
    dte(jr)    = dstep*(te_m(jr_r) - te_m(jr_l))
    dne(jr)    = dstep*(ne_m(jr_r) - ne_m(jr_l))
    dq(jr)     = dstep*(q_m(jr_r) - q_m(jr_l))
    dvpar(jr)  = dstep*(vpar_m(jr_r) - vpar_m(jr_l))
    do jion=1, ns_in-1
        dni(jion, jr) = dstep*(ni_m(jion, jr_r) - ni_m(jion, jr_l))
    enddo
    dr(jr) = drmin(jr)/a0_m    ! gradients w.r.t. minor radius even for s-alpha geometry
    drhodr(jr) = drho/drmin(jr)
enddo

!--------------------
! Populate prof_in_m

dims = 0
dims(1) = n_scalars
dims(2) = n_inputs
dims(3) = n_outputs
dims(4) = ns_in

scal_in_m(1) = RTOR
scal_in_m(2) = BTOR
scal_in_m(3) = a0_m
scal_in_m(4:  8) = mass_in(1: 5)
scal_in_m(9: 13) = zs_in(1: 5)

prof_in_m(:,  1) = rho_m
prof_in_m(:,  2) = ametr_m
prof_in_m(:,  3) = rmaj_m
prof_in_m(:,  4) = elon_m
prof_in_m(:,  5) = tria_m
prof_in_m(:,  6) = q_m
prof_in_m(:,  7) = ti_m
prof_in_m(:,  8) = te_m
prof_in_m(:,  9) = ne_m
prof_in_m(:, 10) = vpar_m
prof_in_m(:, 11) = ni_m(1, :)
prof_in_m(:, 12) = ni_m(2, :)
prof_in_m(:, 13) = ni_m(3, :)
prof_in_m(:, 14) = ni_m(4, :)
prof_in_m(:, 15) = zimp_m(1, :)
prof_in_m(:, 16) = zimp_m(2, :)
prof_in_m(:, 17) = zimp_m(3, :)
prof_in_m(:, 18) = drmin
prof_in_m(:, 19) = drmaj
prof_in_m(:, 20) = delong
prof_in_m(:, 21) = dtrian
prof_in_m(:, 22) = dti
prof_in_m(:, 23) = dte
prof_in_m(:, 24) = dne
prof_in_m(:, 25) = dq
prof_in_m(:, 26) = dvpar
prof_in_m(:, 27) = dr
prof_in_m(:, 28) = drhodr
prof_in_m(:, 29) = dni(1, :)
prof_in_m(:, 30) = dni(2, :)
prof_in_m(:, 31) = dni(3, :)
prof_in_m(:, 32) = dni(4, :)

!--------------
! Send MPI jobs
!--------------

call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

nworkers = 40  ! A submultiple of nrho_m!
chunk = nrho_m / nworkers

call MPI_Comm_spawn(worker_exe, MPI_ARGV_NULL, nworkers, MPI_INFO_NULL, 0, MPI_COMM_SELF, intercomm, errcodes, ierr)
print *, "MPI workers = ", nworkers, nrho_m

! Send dimensions and data to workers
do i=0, nworkers-1
    i1 = i * chunk + 1
    i2 = (i + 1) * chunk
    dims(5) = i1
    dims(6) = i2
    call MPI_Send(dims, SIZE(dims), MPI_INTEGER, i, 0, intercomm, ierr)
enddo
do i=0, nworkers-1
    call MPI_Send(scal_in_m, n_scalars, MPI_DOUBLE_PRECISION, i, 0, intercomm, ierr)
enddo
do i=0, nworkers-1
    i1 = i * chunk + 1
    i2 = (i + 1) * chunk
    call MPI_Send(prof_in_m(i1:i2, :), chunk * n_inputs, MPI_DOUBLE_PRECISION, i, 0, intercomm, ierr)
enddo

! Receive results from each worker
do i=0, nworkers-1
    i1 = i * chunk + 1
    i2 = (i + 1) * chunk
    call MPI_Recv(prof_out_m(i1:i2, :), chunk * n_outputs, MPI_DOUBLE_PRECISION, i, 1, intercomm, status, ierr)
enddo

! Interpolate back to ASTRA radial grid

call qinterp(rho_m, prof_out_m(:, 1), nrho_m, rho_m(1:NA1),       chii_m(1:NA1), NA1)
call qinterp(rho_m, prof_out_m(:, 2), nrho_m, rho_m(1:NA1),       chie_m(1:NA1), NA1)
call qinterp(rho_m, prof_out_m(:, 4), nrho_m, rho_m(1:NA1), elec_pflux_m(1:NA1), NA1)
call qinterp(rho_m, prof_out_m(:, 3), nrho_m, rho_m(1:NA1),      vippd_m(1:NA1), NA1)
call qinterp(rho_m, prof_out_m(:, 4), nrho_m, rho_m(1:NA1),      vittd_m(1:NA1), NA1)
call qinterp(rho_m, prof_out_m(:, 5), nrho_m, rho_m(1:NA1),     vippi1_m(1:NA1), NA1)
call qinterp(rho_m, prof_out_m(:, 6), nrho_m, rho_m(1:NA1),     vitti1_m(1:NA1), NA1)
call qinterp(rho_m, prof_out_m(:, 7), nrho_m, rho_m(1:NA1),       j_boot(1:NA1), NA1)

return
end subroutine neo_parent
