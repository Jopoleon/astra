module a2tglf

implicit none

integer, parameter :: nspec_max=5

type tglf_output
    double precision, allocatable, dimension(:) :: chi_i, chi_e, e_pflux, equipart, &
        mom_flux, gamma, omega
    double precision, allocatable, dimension(:, :) :: ion_pflux
endtype tglf_output
type(tglf_output) :: tglf_out 

contains

!---------------------------------------------------------------------
    subroutine tglf_alloc()

    use scalars, only: NA1

    if (.not. allocated(tglf_out%chi_i)) then
        allocate(tglf_out%chi_i(NA1), tglf_out%chi_e(NA1), tglf_out%e_pflux(NA1), &
             tglf_out%equipart(NA1), tglf_out%mom_flux(NA1), &
             tglf_out%gamma(NA1), tglf_out%omega(NA1))
        allocate(tglf_out%ion_pflux(nspec_max-1, NA1))
    endif

    end subroutine tglf_alloc
  
!---------------------------------------------------------------------
    subroutine tglf_ipc(rho_norm_max)

    use omp_lib
    use pi_const, only: GP2  
    use read_input, only: equ_file, exp_file, awd, astra_exe
    use scalars, only: NA1, BTOR, RTOR, ROC, AMJ, AIM1, AIM2, AIM3, ZMJ
    use status, only: NE, TE, NI, TI, ZEF, PBLON, PBPER, PFAST, &
        ZIM1, ZIM2, ZIM3, NIZ1, NIZ2, NIZ3, ER, MU, FP_NORM, &
        RHO, AMETR, SHIF, ELON, NDEUT, NTRIT, TRIA, VTOR, G11, VPOL, VRS
    use parameters_a2equil, only: equil_now
    use numerical_tools, only: qinterp

    logical, parameter :: debug_elite=.false.
    integer, parameter :: n_arr_out=15, nrho_m=64, nworkers=64, ipcId=0, &
        n_dims=7, n_scalars=8, n_inputs=47, nthe_elite=400, mpol=6
    double precision, parameter :: c_vpol=1.d0

    double precision, intent(in), optional :: rho_norm_max

    logical :: first_call=.True.
    integer :: i, jr, jrho, jr_r, jr_l, jgamma_max, jion, nchunk
    integer :: ns_in, geom_flag=1              ! Number of species, including electrons
    integer :: jthe, jthe_rev, nrho_equ, nthe_equ ! for ELITE geometry
    integer :: t_wall1, t_wall2, rate, max_nworkers, stat, &
        semID, shmID_dims, shmID_vars, shmID_arrs
    integer, dimension(n_dims) :: dims_in

    double precision, dimension(n_arr_out, nrho_m) :: prof_out
    double precision :: bmod, bpolz, xstep, rho_min, rho_max, dstep, a0_m, gradrhosq_inv, dtheta_elite
    double precision, dimension(n_scalars) :: scal_in
    double precision, dimension(nrho_m) :: drmin, drmaj, drho, dte, dne, dq, &
        dptot, delong, dtrian, dvpar, dvper, drhodr, dr, dv_r
    double precision, dimension(NA1) :: rmaj_as, q_as, ni_main_as, &
        vexb_as, vpar_as, vper_as, chie_as, chii_as, e_pflux_as, ptot_as
    double precision, dimension(nrho_m) :: rho_m, gamma_max, omega_max, kymax, &
        te_m, ne_m, vpar_m, vper_m, vexb_m, &
        ametr_m, elon_m, tria_m, rmaj_m, ptot_m, q_m, zef_m, pfn_m
    double precision, dimension(nspec_max-1) :: zi_max
    double precision, dimension(nrho_m, nspec_max-1) :: dti, dni, ni_m, ti_m, zi_m 
    double precision, dimension(NA1, nspec_max-1) :: i_pflux_as
    double precision, dimension(n_inputs, nrho_m) :: prof_in
    character(len=32) :: str_nworkers
    character(len=64) :: SBP_NAME
! ELITE
    double precision, allocatable, dimension(:) :: theta_equ, pfn_equ
    double precision, allocatable, dimension(:, :) :: RR_tg, ZZ_tg, Bp_tg
    double precision, dimension(nthe_elite) :: theta_elite, RR_elite, ZZ_elite, Bp_elite
    character(len=128) :: f_elite, ipc_file, astra_task

    save semID, shmID_dims, shmID_vars, shmID_arrs, first_call

    call SYSTEM_CLOCK(t_wall1, rate)

    write(ipc_file, '(5A, i0, 2A)') TRIM(awd), '/tmp/', TRIM(exp_file), &
        TRIM(equ_file), '-', ipcId, '.ipc', char(0)

    if (first_call) then
        call get_environment_variable("MAX_NWORKERS", str_nworkers, status=stat)
        if (stat /= 0) then
! Not set -> use a fallback (e.g., half of logical CPUs)
            max_nworkers = omp_get_num_procs() / 2
            print *, "MAX_NWORKERS not set, using fallback:", max_nworkers
        else ! Convert string to integer safely
            read(str_nworkers, *) max_nworkers
            print *, "MAX_NWORKERS from environment:", max_nworkers
        endif
        if (nworkers > max_nworkers) then
            write(*, '(A, i3, A, i3)') '>>> ERROR nworkers=', nworkers, ' exceeds #physical CPUs=', max_nworkers
            print*, 'Quitting ASTRA'
            stop
        endif
        if (mod(nrho_m, nworkers) /= 0) then
            write(*, '(A, i3, A, i3)') '>>> ERROR nrho_m=', nrho_m, ' is no multiple of nworkers=', nworkers
            print*, 'Quitting ASTRA'
            stop
        endif
        if (nrho_m > NA1) then
            write(*, '(A, i3, A, i3)') '>>> Warning nrho_m=', nrho_m, ' larger than NA1=', NA1
            print*, 'Possible profile overfit on TGLF grid'
        endif
        SBP_NAME = "xpr/tglfi" // char(0)
        astra_task = TRIM(astra_exe) // char(0)
        call initialise_ipc(nrho_m, n_dims, n_scalars, n_inputs, n_arr_out, &
            nworkers, SBP_NAME, ipc_file, astra_task, semID, shmID_dims, &
            shmID_vars, shmID_arrs, ipcId)
    endif
 
    nchunk = nrho_m / nworkers

! Interpolate from ASTRA grid to TGLF grid
    rho_min = RHO(1)
    if (present(rho_norm_max)) then
        rho_max = rho_norm_max*ROC
    else
        rho_max = ROC
    endif
    xstep = (rho_max - rho_min)/(nrho_m - 1.)
    rho_m = (/ (rho_min + (jr - 1.)*xstep, jr=1, nrho_m) /)

    zi_m(1, :) = ZMJ
    call qinterp(RHO(1:NA1), ZIM1(1:NA1), NA1, rho_m, zi_m(:, 2), nrho_m)
    call qinterp(RHO(1:NA1), ZIM2(1:NA1), NA1, rho_m, zi_m(:, 3), nrho_m)
    call qinterp(RHO(1:NA1), ZIM3(1:NA1), NA1, rho_m, zi_m(:, 4), nrho_m)
    call qinterp(RHO(1:NA1), NIZ1(1:NA1), NA1, rho_m, ni_m(:, 2), nrho_m)
    call qinterp(RHO(1:NA1), NIZ2(1:NA1), NA1, rho_m, ni_m(:, 3), nrho_m)
    call qinterp(RHO(1:NA1), NIZ3(1:NA1), NA1, rho_m, ni_m(:, 4), nrho_m)
    call qinterp(RHO(1:NA1),   TI(1:NA1), NA1, rho_m, ti_m(:, 1), nrho_m)
    call qinterp(RHO(1:NA1),      TE(1:NA1), NA1, rho_m,    te_m, nrho_m)
    call qinterp(RHO(1:NA1),      NE(1:NA1), NA1, rho_m,    ne_m, nrho_m)
    call qinterp(RHO(1:NA1),     ZEF(1:NA1), NA1, rho_m,   zef_m, nrho_m)
    call qinterp(RHO(1:NA1),   AMETR(1:NA1), NA1, rho_m, ametr_m, nrho_m)
    call qinterp(RHO(1:NA1),    ELON(1:NA1), NA1, rho_m,  elon_m, nrho_m)
    call qinterp(RHO(1:NA1),    TRIA(1:NA1), NA1, rho_m,  tria_m, nrho_m)
    call qinterp(RHO(1:NA1), FP_NORM(1:NA1), NA1, rho_m,   pfn_m, nrho_m)

    ti_m(:, 2) = ti_m(:, 1)
    ti_m(:, 3) = ti_m(:, 1)
    ti_m(:, 4) = ti_m(:, 1)

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

    call qinterp(RHO(1:NA1), ni_main_as(1:NA1), NA1, rho_m, ni_m(:, 1), nrho_m)
    call qinterp(RHO(1:NA1),    rmaj_as(1:NA1), NA1, rho_m,     rmaj_m, nrho_m)
    call qinterp(RHO(1:NA1),       q_as(1:NA1), NA1, rho_m,        q_m, nrho_m)
    call qinterp(RHO(1:NA1),    ptot_as(1:NA1), NA1, rho_m,     ptot_m, nrho_m)
    call qinterp(RHO(1:NA1),    vpar_as(1:NA1), NA1, rho_m,     vpar_m, nrho_m)
    call qinterp(RHO(1:NA1),    vper_as(1:NA1), NA1, rho_m,     vper_m, nrho_m)
    call qinterp(RHO(1:NA1),    vexb_as(1:NA1), NA1, rho_m,     vexb_m, nrho_m)

! Reference length
    a0_m = AMETR(NA1)

    do jr=1, nrho_m
        ni_m(jr, 2) = max(1.e-9, ni_m(jr, 2))
        ni_m(jr, 3) = max(1.e-9, ni_m(jr, 3))
        ni_m(jr, 4) = max(1.e-9, ni_m(jr, 4))
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
            jr_l = 1
        else if (jr == nrho_m) then
            jr_r = nrho_m
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
        do jion=1, nspec_max-1
            dni(jr, jion) = dstep*(ni_m(jr_r, jion) - ni_m(jr_l, jion))
            dti(jr, jion) = dstep*(ti_m(jr_r, jion) - ti_m(jr_l, jion))
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

        theta_equ = equil_now%coord_sys%position%theta2d
        dtheta_elite = GP2/dble(nthe_elite - 1)
        theta_elite = (/ ((jthe - 1.)*dtheta_elite, jthe=1, nthe_elite) /)
    endif

!--------------------
! IPC parallelisation
!--------------------

    dims_in(1) = nchunk
    dims_in(2) = n_inputs
    dims_in(3) = n_arr_out
    dims_in(4) = nrho_m
    dims_in(5) = nspec_max
    dims_in(6) = ns_in
    dims_in(7) = geom_flag

    scal_in(1) = BTOR
    scal_in(2) = RTOR
    scal_in(3) = AMETR(NA1)
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
    prof_in( 9, :) = te_m 
    prof_in(10, :) = ne_m 
    prof_in(11, :) = zef_m
    prof_in(12, :) = vpar_m
    prof_in(13, :) = vper_m
    prof_in(14, :) = vexb_m
    prof_in(15, :) = drmin
    prof_in(16, :) = drmaj
    prof_in(17, :) = drho 
    prof_in(18, :) = delong
    prof_in(19, :) = dtrian
    prof_in(20, :) = dptot
    prof_in(21, :) = dte
    prof_in(22, :) = dne
    prof_in(23, :) = dq
    prof_in(24, :) = dvper
    prof_in(25, :) = dv_r
    prof_in(26, :) = dr
    prof_in(27, :) = drhodr
    prof_in(28, :) = ni_m(:, 1)
    prof_in(29, :) = ni_m(:, 2)
    prof_in(30, :) = ni_m(:, 3)
    prof_in(31, :) = ni_m(:, 4)
    prof_in(32, :) = ti_m(:, 1)
    prof_in(33, :) = ti_m(:, 2)
    prof_in(34, :) = ti_m(:, 3)
    prof_in(35, :) = ti_m(:, 4)
    prof_in(36, :) = zi_m(:, 1)
    prof_in(37, :) = zi_m(:, 2)
    prof_in(38, :) = zi_m(:, 3)
    prof_in(39, :) = zi_m(:, 4)
    prof_in(40, :) = dni(:, 1)
    prof_in(41, :) = dni(:, 2)
    prof_in(42, :) = dni(:, 3)
    prof_in(43, :) = dni(:, 4)
    prof_in(44, :) = dti(:, 1)
    prof_in(45, :) = dti(:, 2)
    prof_in(46, :) = dti(:, 3)
    prof_in(47, :) = dti(:, 4)

    if (first_call) then
        call fill_int_shm(n_dims, dims_in, shmID_dims)
        first_call = .False.
    endif

! **** Fill shared memory segments
    call fill_dbl_shm(n_scalars, scal_in, shmID_vars)
    call fill_dbl_shm(nrho_m*n_inputs, prof_in, shmID_arrs)

! **** Free each semaphore
    do i=1, nworkers
        call unlock_sbp(i, semID)
    enddo

! **** Synchronisation point
    call wait4all(semID)

! **** Collect data from ShMem
    do i=1, nworkers
        call sbp2astra(i, nchunk, n_arr_out, ipc_file, prof_out(1, 1))
    enddo

! Interpolate back to ASTRA radial grid
    e_pflux_as = 0.
    i_pflux_as = 0.
    chie_as = 0.
    chii_as = 0.

    call qinterp(rho_m, prof_out(1, :), nrho_m, RHO(1:NA1),    chii_as(1:NA1), NA1, extrap_right=0.d0)
    call qinterp(rho_m, prof_out(2, :), nrho_m, RHO(1:NA1),    chie_as(1:NA1), NA1, extrap_right=0.d0)
    call qinterp(rho_m, prof_out(4, :), nrho_m, RHO(1:NA1), e_pflux_as(1:NA1), NA1, extrap_right=0.d0)
    call qinterp(rho_m, prof_out(3, :), nrho_m, RHO(1:NA1), tglf_out%mom_flux, NA1, extrap_right=0.d0)
    call qinterp(rho_m, prof_out(5, :), nrho_m, RHO(1:NA1), tglf_out%equipart, NA1, extrap_right=0.d0)
    call qinterp(rho_m, prof_out(6, :), nrho_m, RHO(1:NA1),    tglf_out%gamma, NA1, extrap_right=0.d0)
    call qinterp(rho_m, prof_out(7, :), nrho_m, RHO(1:NA1),    tglf_out%omega, NA1, extrap_right=0.d0)
    do jion=1, nspec_max-1
        call qinterp(rho_m, prof_out(7+jion, :), nrho_m, RHO(1:NA1), i_pflux_as(1:NA1, jion), NA1, extrap_right=0.d0)
    enddo

    do jrho=1, NA1
        gradrhosq_inv = VRS(jrho)/G11(jrho)
        tglf_out%chi_i(jrho)   = chii_as(jrho)*gradrhosq_inv ! m^2/s
        tglf_out%chi_e(jrho)   = chie_as(jrho)*gradrhosq_inv ! m^2/s
        tglf_out%e_pflux(jrho) = e_pflux_as(jrho)*gradrhosq_inv/a0_m
        tglf_out%ion_pflux(1, jrho) = i_pflux_as(jrho, 1)*gradrhosq_inv/a0_m/(ni_main_as(jrho)/NE(jrho))  ! main ion particle flux
        tglf_out%ion_pflux(2, jrho) = i_pflux_as(jrho, 2)*gradrhosq_inv/a0_m/(NIZ1(jrho)/NE(jrho))  ! 1st imp particle flux
        tglf_out%ion_pflux(3, jrho) = i_pflux_as(jrho, 3)*gradrhosq_inv/a0_m/(NIZ2(jrho)/NE(jrho))  ! 2nd imp particle flux
        tglf_out%ion_pflux(4, jrho) = i_pflux_as(jrho, 4)*gradrhosq_inv/a0_m/(NIZ3(jrho)/NE(jrho))  ! 2nd imp particle flux
    enddo

    call SYSTEM_CLOCK(t_wall2, rate)
    print*, "XPR wall time", dble(t_wall2 - t_wall1)/dble(rate)

    end subroutine tglf_ipc

end module a2tglf
