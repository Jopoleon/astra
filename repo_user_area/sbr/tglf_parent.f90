subroutine tglf_parent(chi_i, chi_e, e_pflux, vimp1, vimp2, i_mflux_as, exchi_as, gamma_as, omega_as)

use mpi

use parameter_inc, only: NRD
use const_inc, only: BTOR, RTOR, AMJ, AIM1, AIM2, AIM3, ZMJ, NA1
use status_inc, only: NE, TE, NI, TI, ZEF, PBLON, PBPER, PFAST, &
    ZIM1, ZIM2, ZIM3, NIZ1, NIZ2, NIZ3, ER, MU, FP_NORM, &
    RHO, AMETR, SHIF, ELON, NDEUT, TRIA, VTOR, G11, VPOL, VRS

implicit none

integer, parameter :: n_scalars=10, n_inputs=40, n_outputs=15, n_dims=8, nrho_m=80, nspec_max=5
double precision, parameter :: c_vpol=1.d0

double precision, dimension(NRD), intent(out) :: chi_i, chi_e, e_pflux, vimp2, vimp1, &
     i_mflux_as, exchi_as, gamma_as, omega_as

integer :: ierr, info, intercomm, errcodes(100), status(MPI_STATUS_SIZE)
integer :: jr, jrho, jr_r, jr_l, jgamma_max, jion
integer :: ns_in              ! Number of species, including electrons
integer :: i, i1, i2, chunk, nworkers, dims(n_dims)

double precision, dimension(n_scalars) :: scal_in_m
double precision, dimension(n_outputs, nrho_m) :: prof_out_m
double precision :: bmod, bpolz, xstep, rho_min, rho_max, dstep, a0_m, gradrhosq_inv
double precision, dimension(nrho_m) :: drmin, drmaj, drho, dti, dte, dne, dq, &
    dptot, delong, dtrian, dvpar, dvper, drhodr, dr, dv_r
double precision, dimension(NRD) :: rmaj_as, q_as, ni_main_as, &
    vexb_as, vpar_as, vper_as, chie_as, chii_as, e_pflux_as, ptot_as
double precision, dimension(nrho_m) :: rho_m, gamma_max, omega_max, kymax, &
    ti_m, te_m, ne_m, vpar_m, vper_m, vexb_m, &
    ametr_m, elon_m, tria_m, rmaj_m, ptot_m, q_m, zef_m, pfn_m
double precision, dimension(nspec_max) :: mass_in
double precision, dimension(nspec_max-2) :: zimp_max
double precision, dimension(nspec_max-1, nrho_m) :: dni, ni_m, i_pflux
double precision, dimension(nspec_max-2, nrho_m) :: zimp_m 
double precision, dimension(nspec_max-1, NRD) :: i_pflux_as
double precision, allocatable, dimension(:, :) :: send_buffer
character(len=256) :: worker_exe
character(len=1), dimension(1) :: args
args(1) = ''

worker_exe = "xpr/tglf.x"

! Interpolate from ASTRA grid to TGLF grid
rho_min = RHO(1)
rho_max = RHO(NA1)
xstep = (rho_max - rho_min)/(nrho_m - 1.)
rho_m = (/ (rho_min + (jr - 1.)*xstep, jr=1, nrho_m) /)

call qinterp(RHO(1:NA1), ZIM1(1:NA1), NA1, rho_m, zimp_m(1, :), nrho_m)
call qinterp(RHO(1:NA1), ZIM2(1:NA1), NA1, rho_m, zimp_m(2, :), nrho_m)
call qinterp(RHO(1:NA1), ZIM3(1:NA1), NA1, rho_m, zimp_m(3, :), nrho_m)
call qinterp(RHO(1:NA1), NIZ1(1:NA1), NA1, rho_m,   ni_m(2, :), nrho_m)
call qinterp(RHO(1:NA1), NIZ2(1:NA1), NA1, rho_m,   ni_m(3, :), nrho_m)
call qinterp(RHO(1:NA1), NIZ3(1:NA1), NA1, rho_m,   ni_m(4, :), nrho_m)
call qinterp(RHO(1:NA1),      TI(1:NA1), NA1, rho_m,    ti_m, nrho_m)
call qinterp(RHO(1:NA1),      TE(1:NA1), NA1, rho_m,    te_m, nrho_m)
call qinterp(RHO(1:NA1),      NE(1:NA1), NA1, rho_m,    ne_m, nrho_m)
call qinterp(RHO(1:NA1),     ZEF(1:NA1), NA1, rho_m,   zef_m, nrho_m)
call qinterp(RHO(1:NA1),   AMETR(1:NA1), NA1, rho_m, ametr_m, nrho_m)
call qinterp(RHO(1:NA1),    ELON(1:NA1), NA1, rho_m,  elon_m, nrho_m)
call qinterp(RHO(1:NA1),    TRIA(1:NA1), NA1, rho_m,  tria_m, nrho_m)
call qinterp(RHO(1:NA1), FP_NORM(1:NA1), NA1, rho_m,   pfn_m, nrho_m)

do jrho=1, NA1
    if (NDEUT(jrho) >= 0.01*NE(jrho)) then
        ni_main_as(jrho) = NDEUT(jrho)
    else ! likely: NDEUT not defined in equ file, hence zero
        ni_main_as(jrho) = NI(jrho)
    endif
    rmaj_as(jrho) = RTOR + SHIF(jrho)
    q_as(jrho)    = 1./MU(jrho)
    ptot_as(jrho) = NE(jrho)*TE(jrho) + ni_main_as(jrho)*TI(jrho) + NIZ1(jrho)*TI(jrho) + pfast(jrho) + 0.5*(pblon(jrho) + pbper(jrho))
    bpolz = BTOR*AMETR(jrho)*MU(jrho)/RTOR
    bmod = sqrt(BTOR**2 + bpolz**2)
    vper_as(jrho) = ER(jrho)/(RTOR*bpolz) ! vexb in m/s --> Omega_E
    vpar_as(jrho) = VTOR(jrho) * BTOR/bmod + c_vpol* VPOL(jrho) * bpolz/bmod
    vexb_as(jrho) = -ER(jrho)/bmod ! vexb in m/s (vperp = vexb since the diamagnetic velocity is the curvature drift ac
enddo

call qinterp(RHO(1:NA1), ni_main_as(1:NA1), NA1, rho_m, ni_m(1, :), nrho_m)
call qinterp(RHO(1:NA1),    rmaj_as(1:NA1), NA1, rho_m,     rmaj_m, nrho_m)
call qinterp(RHO(1:NA1),       q_as(1:NA1), NA1, rho_m,        q_m, nrho_m)
call qinterp(RHO(1:NA1),    ptot_as(1:NA1), NA1, rho_m,     ptot_m, nrho_m)
call qinterp(RHO(1:NA1),    vpar_as(1:NA1), NA1, rho_m,     vpar_m, nrho_m)
call qinterp(RHO(1:NA1),    vper_as(1:NA1), NA1, rho_m,     vper_m, nrho_m)
call qinterp(RHO(1:NA1),    vexb_as(1:NA1), NA1, rho_m,     vexb_m, nrho_m)

! Reference length
a0_m = AMETR(NA1)

! Species cmassses and charges
mass_in = 0.d0
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

! Number of species
ns_in = nspec_max

! These will be reset locally in the radial loop
zimp_max(1) = MAXVAL(ZIM1(1:NA1))
zimp_max(2) = MAXVAL(ZIM2(1:NA1))
zimp_max(3) = MAXVAL(ZIM3(1:NA1))
if (zimp_max(3) >= 1.) then
    ns_in = 5
else
    ns_in = 4
endif
if (zimp_max(2) < 1.) ns_in = 3
if (zimp_max(3) >= 1. .and. ns_in == 3) then
    ns_in = 4
endif
if (zimp_max(1) < 1.) ns_in = 2
if (zimp_max(2) >= 1. .and. ns_in == 2) then
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
    drho(jr)   = dstep*(  rho_m(jr_r) -   rho_m(jr_l))
    delong(jr) = dstep*( elon_m(jr_r) -  elon_m(jr_l))
    dtrian(jr) = dstep*( tria_m(jr_r) -  tria_m(jr_l))
    dptot(jr)  = dstep*( ptot_m(jr_r) -  ptot_m(jr_l)) * 1E3*1E13
    dti(jr)    = dstep*(ti_m(jr_r) - ti_m(jr_l))
    dte(jr)    = dstep*(te_m(jr_r) - te_m(jr_l))
    dne(jr)    = dstep*(ne_m(jr_r) - ne_m(jr_l))
    dq(jr)     = dstep*(q_m(jr_r) - q_m(jr_l))
    dvper(jr)  = dstep*(vper_m(jr_r) - vper_m(jr_l))
    do jion=1, nspec_max-1
        dni(jion, jr) = dstep*(ni_m(jion, jr_r) - ni_m(jion, jr_l))
    enddo
    dv_r(jr) = dstep* &
        (vpar_m(jr_r)/(rmaj_m(jr_r) + ametr_m(jr_r)) - &
         vpar_m(jr_l)/(rmaj_m(jr_l) + ametr_m(jr_l)))
    dr(jr) = drmin(jr)/a0_m    ! gradients w.r.t. minor radius even for s-alpha geometry
    drhodr(jr) = drho(jr)/drmin(jr)
enddo

dims(1) = nrho_m
dims(2) = n_scalars
dims(3) = n_inputs
dims(4) = n_outputs
dims(5) = ns_in
dims(6) = nspec_max

scal_in_m(1) = RTOR
scal_in_m(2) = BTOR
scal_in_m(3) = a0_m
scal_in_m(4) = ZMJ
scal_in_m(5: 4+nspec_max) = mass_in(1: nspec_max)

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
call MPI_Comm_spawn(TRIM(worker_exe), args, nworkers, MPI_INFO_NULL, 0, MPI_COMM_SELF, intercomm, errcodes, ierr)

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
    send_buffer( 9, :) = ti_m(i1:i2)
    send_buffer(10, :) = te_m(i1:i2) 
    send_buffer(11, :) = ne_m(i1:i2) 
    send_buffer(12, :) = zef_m(i1:i2)
    send_buffer(13, :) = vpar_m(i1:i2)
    send_buffer(14, :) = vper_m(i1:i2)
    send_buffer(15, :) = vexb_m(i1:i2)
    send_buffer(16, :) = drmin(i1:i2)
    send_buffer(17, :) = drmaj(i1:i2)
    send_buffer(18, :) = drho(i1:i2) 
    send_buffer(19, :) = delong(i1:i2)
    send_buffer(20, :) = dtrian(i1:i2)
    send_buffer(21, :) = dptot(i1:i2)
    send_buffer(22, :) = dti(i1:i2)
    send_buffer(23, :) = dte(i1:i2)
    send_buffer(24, :) = dne(i1:i2)
    send_buffer(25, :) = dq(i1:i2)
    send_buffer(26, :) = dvper(i1:i2)
    send_buffer(27, :) = dv_r(i1:i2)
    send_buffer(28, :) = dr(i1:i2)
    send_buffer(29, :) = drhodr(i1:i2)
    send_buffer(30, :) = ni_m(1, i1:i2)
    send_buffer(31, :) = ni_m(2, i1:i2)
    send_buffer(32, :) = ni_m(3, i1:i2)
    send_buffer(33, :) = ni_m(4, i1:i2)
    send_buffer(34, :) = zimp_m(1, i1:i2)
    send_buffer(35, :) = zimp_m(2, i1:i2)
    send_buffer(36, :) = zimp_m(3, i1:i2)
    send_buffer(37, :) = dni(1, i1:i2)
    send_buffer(38, :) = dni(2, i1:i2)
    send_buffer(39, :) = dni(3, i1:i2)
    send_buffer(40, :) = dni(4, i1:i2)
    call MPI_Send(send_buffer, chunk * n_inputs, MPI_DOUBLE_PRECISION, i, 102+i, intercomm, ierr)
enddo

! Receive results from each worker
do i=0, nworkers-1
    i1 = i * chunk + 1
    i2 = (i + 1) * chunk
    call MPI_Recv(prof_out_m(:, i1:i2), chunk * n_outputs, MPI_DOUBLE_PRECISION, i, 103+i, intercomm, status, ierr)
enddo
call MPI_Barrier(intercomm, ierr)  ! Optional: ensure child finished before next step

! Interpolate back to ASTRA radial grid
e_pflux_as = 0.
i_pflux_as = 0.
i_mflux_as = 0.
chie_as  = 0.
chii_as  = 0.
exchi_as = 0.
gamma_as = 0.
omega_as = 0.
call qinterp(rho_m, prof_out_m(1, :), nrho_m, RHO(1:NA1),    chii_as(1:NA1), NA1)
call qinterp(rho_m, prof_out_m(2, :), nrho_m, RHO(1:NA1),    chie_as(1:NA1), NA1)
call qinterp(rho_m, prof_out_m(3, :), nrho_m, RHO(1:NA1), i_mflux_as(1:NA1), NA1)
call qinterp(rho_m, prof_out_m(4, :), nrho_m, RHO(1:NA1), e_pflux_as(1:NA1), NA1)
call qinterp(rho_m, prof_out_m(5, :), nrho_m, RHO(1:NA1),   exchi_as(1:NA1), NA1)
call qinterp(rho_m, prof_out_m(6, :), nrho_m, RHO(1:NA1),   gamma_as(1:NA1), NA1)
call qinterp(rho_m, prof_out_m(7, :), nrho_m, RHO(1:NA1),   omega_as(1:NA1), NA1)
do jion=1, nspec_max-1
    call qinterp(rho_m, prof_out_m(7+jion, :), nrho_m, RHO(1:NA1), i_pflux_as(jion, 1:NA1), NA1)
enddo

do jrho=1, NA1
    gradrhosq_inv = VRS(jrho)/G11(jrho)
    chi_i(jrho) = chii_as(jrho)*gradrhosq_inv ! \chi_i, m^2/s
    chi_e(jrho) = chie_as(jrho)*gradrhosq_inv ! \chi_e, m^2/s
    e_pflux(jrho) = e_pflux_as(jrho)*gradrhosq_inv/a0_m ! D flux
    vimp1(jrho) = i_pflux_as(2, jrho)*gradrhosq_inv/a0_m/(NIZ1(jrho)/NE(jrho))  ! 1st imp convection
    vimp2(jrho) = i_pflux_as(3, jrho)*gradrhosq_inv/a0_m/(NIZ2(jrho)/NE(jrho))  ! 2nd imp convection
enddo

return
end subroutine tglf_parent
