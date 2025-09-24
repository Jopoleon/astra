subroutine tglf_ipc

use parameter_inc, only: NRD
use io_mod, only: equ_file, exp_file
use ipc_mod, only: mem_tglf, n_sbp_arr_out
use const_inc, only: NA1, GP2, BTOR, RTOR, ABC, AMJ, AIM1, AIM2, AIM3, ZMJ
use status_inc, only: NE, TE, NI, TI, ZEF, PBLON, PBPER, PFAST, &
    ZIM1, ZIM2, ZIM3, NIZ1, NIZ2, NIZ3, ER, MU, FP_NORM, &
    RHO, AMETR, SHIF, ELON, NDEUT, NTRIT, TRIA, VTOR, G11, VPOL, VRS
use parameters_a2equil, only: equil_now

implicit none

logical, parameter :: debug_elite=.false.
integer, parameter :: n_dims=7, n_scalars=8, n_inputs=41, nrho_m=80, nworkers=40, nspec_max=5, nthe_elite=400, mpol=6
double precision, parameter :: c_vpol=1.d0

logical :: first_call=.True.
integer :: i, jr, jrho, jr_r, jr_l, jgamma_max, jion, nchunk
integer :: ns_in, geom_flag=1              ! Number of species, including electrons
integer :: jthe, jthe_rev, nrho_equ, nthe_equ ! for ELITE geometry
integer :: t_wall1, t_wall2, rate
integer, dimension(n_dims) :: dims_in

double precision, dimension(n_sbp_arr_out, nrho_m) :: prof_out
double precision :: bmod, bpolz, xstep, rho_min, rho_max, dstep, a0_m, gradrhosq_inv, dtheta_elite
double precision, dimension(n_scalars) :: scal_in
double precision, dimension(nrho_m) :: drmin, drmaj, drho, dti, dte, dne, dq, &
    dptot, delong, dtrian, dvpar, dvper, drhodr, dr, dv_r
double precision, dimension(NRD) :: rmaj_as, q_as, ni_main_as, &
    vexb_as, vpar_as, vper_as, chie_as, chii_as, e_pflux_as, i_mflux_as, ptot_as
double precision, dimension(nrho_m) :: rho_m, gamma_max, omega_max, kymax, &
    ti_m, te_m, ne_m, vpar_m, vper_m, vexb_m, &
    ametr_m, elon_m, tria_m, rmaj_m, ptot_m, q_m, zef_m, pfn_m
double precision, dimension(nspec_max) :: mass_in
double precision, dimension(nspec_max-1) :: zi_max
double precision, dimension(nspec_max-1, nrho_m) :: dni, ni_m, i_pflux
double precision, dimension(nspec_max-1, nrho_m) :: zi_m 
double precision, dimension(nspec_max-1, NRD) :: i_pflux_as
double precision, dimension(n_inputs, nrho_m) :: prof_in
character(len=64), dimension(nworkers) :: SBP_NAMES
! ELITE
double precision, allocatable, dimension(:) :: theta_equ, pfn_equ
double precision, allocatable, dimension(:, :) :: RR_tg, ZZ_tg, Bp_tg
double precision, dimension(nthe_elite) :: theta_elite, RR_elite, ZZ_elite, Bp_elite
character(len=120) :: f_elite

call SYSTEM_CLOCK(t_wall1, rate)

! Interpolate from ASTRA grid to TGLF grid
rho_min = RHO(1)
rho_max = RHO(NA1)
xstep = (rho_max - rho_min)/(nrho_m - 1.)
rho_m = (/ (rho_min + (jr - 1.)*xstep, jr=1, nrho_m) /)

zi_m(1, :) = ZMJ
call qinterp(RHO(1:NA1), ZIM1(1:NA1), NA1, rho_m, zi_m(2, :), nrho_m)
call qinterp(RHO(1:NA1), ZIM2(1:NA1), NA1, rho_m, zi_m(3, :), nrho_m)
call qinterp(RHO(1:NA1), ZIM3(1:NA1), NA1, rho_m, zi_m(4, :), nrho_m)
call qinterp(RHO(1:NA1), NIZ1(1:NA1), NA1, rho_m, ni_m(2, :), nrho_m)
call qinterp(RHO(1:NA1), NIZ2(1:NA1), NA1, rho_m, ni_m(3, :), nrho_m)
call qinterp(RHO(1:NA1), NIZ3(1:NA1), NA1, rho_m, ni_m(4, :), nrho_m)
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
zi_max(1) = ZMJ
zi_max(2) = MAXVAL(ZIM1(1:NA1))
zi_max(3) = MAXVAL(ZIM2(1:NA1))
zi_max(4) = MAXVAL(ZIM3(1:NA1))
if (zi_max(3) >= 1.) then
    ns_in = 5
else
    ns_in = 4
endif
if (zi_max(3) < 1.) ns_in = 3
if (zi_max(4) >= 1. .and. ns_in == 3) then
    ns_in = 4
endif
if (zi_max(2) < 1.) ns_in = 2
if (zi_max(3) >= 1. .and. ns_in == 2) then
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

!------
! ELITE
!------

if (geom_flag == 3) then
    nrho_equ = SIZE(equil_now%coord_sys%position%r, dim=1)
    nthe_equ = SIZE(equil_now%coord_sys%position%r, dim=2)
    allocate(pfn_equ(nrho_equ))
    allocate(theta_equ(nthe_equ))
    allocate(RR_tg(nrho_m, nthe_equ), ZZ_tg(nrho_m, nthe_equ), Bp_tg(nrho_m, nthe_equ))

! Interpolation on TGLF rho-grid

    pfn_equ = (equil_now%profiles_1d%psi - equil_now%profiles_1d%psi(1))/(equil_now%profiles_1d%psi(nrho_equ) - equil_now%profiles_1d%psi(1))

! Interpolation on TGLF rho-grid
    do jthe=1, nthe_equ
        call qinterp(pfn_equ, equil_now%coord_sys%position%r(:, jthe), nrho_equ, pfn_m, RR_tg(:, jthe), nrho_m)
        call qinterp(pfn_equ, equil_now%coord_sys%position%z(:, jthe), nrho_equ, pfn_m, Zz_tg(:, jthe), nrho_m)
        call qinterp(pfn_equ, equil_now%coord_sys%bpcell    (:, jthe), nrho_equ, pfn_m, Bp_tg(:, jthe), nrho_m)
    enddo

    deallocate(pfn_equ)

    theta_equ = equil_now%coord_sys%position%teta2d
    dtheta_elite = GP2/dble(nthe_elite - 1)
    theta_elite = (/ ((jthe - 1.)*dtheta_elite, jthe=1, nthe_elite) /)
endif

!--------------------
! IPC parallelisation
!--------------------

nchunk = nrho_m / nworkers

dims_in(1) = nchunk
dims_in(2) = n_inputs
dims_in(3) = n_sbp_arr_out
dims_in(4) = nrho_m
dims_in(5) = nspec_max
dims_in(6) = ns_in
dims_in(7) = geom_flag

scal_in(1) = BTOR
scal_in(2) = RTOR
scal_in(3) = ABC
scal_in(4) = AMJ
scal_in(5) = AIM1
scal_in(6) = AIM2
scal_in(7) = AIM3
scal_in(8) = ZMJ

prof_in = 0.d0
prof_in( 1, :) = rho_m
prof_in( 2, :) = ametr_m
prof_in( 3, :) = rmaj_m
prof_in( 4, :) = elon_m
prof_in( 5, :) = tria_m
prof_in( 6, :) = q_m 
prof_in( 7, :) = pfn_m
prof_in( 8, :) = ptot_m
prof_in( 9, :) = ti_m
prof_in(10, :) = te_m 
prof_in(11, :) = ne_m 
prof_in(12, :) = zef_m
prof_in(13, :) = vpar_m
prof_in(14, :) = vper_m
prof_in(15, :) = vexb_m
prof_in(16, :) = drmin
prof_in(17, :) = drmaj
prof_in(18, :) = drho 
prof_in(19, :) = delong
prof_in(20, :) = dtrian
prof_in(21, :) = dptot
prof_in(22, :) = dti
prof_in(23, :) = dte
prof_in(24, :) = dne
prof_in(25, :) = dq
prof_in(26, :) = dvper
prof_in(27, :) = dv_r
prof_in(28, :) = dr
prof_in(29, :) = drhodr
prof_in(30, :) = ni_m(1, :)
prof_in(31, :) = ni_m(2, :)
prof_in(32, :) = ni_m(3, :)
prof_in(33, :) = ni_m(4, :)
prof_in(34, :) = zi_m(1, :)
prof_in(35, :) = zi_m(2, :)
prof_in(36, :) = zi_m(3, :)
prof_in(37, :) = zi_m(4, :)
prof_in(38, :) = dni(1, :)
prof_in(39, :) = dni(2, :)
prof_in(40, :) = dni(3, :)
prof_in(41, :) = dni(4, :)

SBP_NAMES = "xpr/tglfi"//char(0)
if (first_call) then
    call initialise_ipc(nrho_m, n_dims, n_scalars, n_inputs, n_sbp_arr_out, nworkers, equ_file, exp_file)
    call fill_dim2shm(dims_in)
    call send_ipc_jobs(nworkers, nchunk, 64, SBP_NAMES)
    first_call = .False.
endif

! **** Fill shared memory segments
call fill_var2shm(scal_in)
call fill_arr2shm(prof_in)

! **** Free each semaphore
do i=1, nworkers
    call unlock_sbp(i)
enddo

! **** Synchronisation point
call wait4all

! **** Collect data from ShMem
do i=1, nworkers
    call sbp2astra(i, prof_out(1, 1))
enddo

! Interpolate back to ASTRA radial grid
e_pflux_as = 0.
i_pflux_as = 0.
i_mflux_as = 0.
chie_as = 0.
chii_as = 0.

call qinterp(rho_m, prof_out(1, :), nrho_m, RHO(1:NA1),      chii_as(1:NA1), NA1) ! chi_i
call qinterp(rho_m, prof_out(2, :), nrho_m, RHO(1:NA1),      chie_as(1:NA1), NA1) !chi_e
call qinterp(rho_m, prof_out(3, :), nrho_m, RHO(1:NA1),   i_mflux_as(1:NA1), NA1)
call qinterp(rho_m, prof_out(4, :), nrho_m, RHO(1:NA1),   e_pflux_as(1:NA1), NA1) ! Electron flux
call qinterp(rho_m, prof_out(5, :), nrho_m, RHO(1:NA1), mem_tglf(1:NA1,  8), NA1) ! Turb. equip.
call qinterp(rho_m, prof_out(6, :), nrho_m, RHO(1:NA1), mem_tglf(1:NA1, 11), NA1) ! gamma
call qinterp(rho_m, prof_out(7, :), nrho_m, RHO(1:NA1), mem_tglf(1:NA1, 12), NA1) ! omega
do jion=1, nspec_max-1
    call qinterp(rho_m, prof_out(7+jion, :), nrho_m, RHO(1:NA1), i_pflux_as(jion, 1:NA1), NA1)
enddo

do jrho=1, NA1
    gradrhosq_inv = VRS(jrho)/G11(jrho)
    mem_tglf(jrho, 1) = chii_as(jrho)*gradrhosq_inv ! \chi_i, m^2/s
    mem_tglf(jrho, 2) = chie_as(jrho)*gradrhosq_inv ! \chi_e, m^2/s
    mem_tglf(jrho, 4) = e_pflux_as(jrho)*gradrhosq_inv/a0_m ! e flux
    mem_tglf(jrho, 13) = i_pflux_as(2, jrho)*gradrhosq_inv/a0_m/(NIZ1(jrho)/NE(jrho))  ! 1st imp convection
    mem_tglf(jrho, 14) = i_pflux_as(3, jrho)*gradrhosq_inv/a0_m/(NIZ2(jrho)/NE(jrho))  ! 2nd imp convection
enddo

call SYSTEM_CLOCK(t_wall2, rate)
print*, "XPR wall time", dble(t_wall2 - t_wall1)/dble(rate)

return
end subroutine tglf_ipc
