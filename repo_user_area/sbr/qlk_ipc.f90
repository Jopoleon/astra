module a2qlk

implicit none

integer, parameter :: nspec_max=7

type qlk_output
    double precision, allocatable, dimension(:) :: chi_i, chi_e, e_pflux, equipart
endtype qlk_output
type(qlk_output) :: qlk_out 

contains

!---------------------------------------------------------------------
    subroutine qlk_alloc

    use const_inc, only: NA1

    if (.not. allocated(qlk_out%chi_i)) then
        allocate(qlk_out%chi_i(NA1), qlk_out%chi_e(NA1), qlk_out%e_pflux(NA1), qlk_out%equipart(NA1))
    endif

    return
    end subroutine qlk_alloc

!---------------------------------------------------------------------
    subroutine qlk_ipc(rho_norm_max)

    use omp_lib
    use parameter_inc, only: NRD
    use io_mod, only: equ_file, exp_file
    use const_inc, only: NA1, BTOR, RTOR, ROC, AMJ, AIM1, AIM2, AIM3, ZMJ
    use status_inc, only: NE, TE, NI, TI, ZIM1, ZIM2, ZIM3, PBLON, PBPER, &
        PFAST, NIZ3, AMAIN, ER, MU, FP_NORM, RHO, AMETR, SHIF, &
        NDEUT, NIZ1, NIZ2, VTOR, NIBM, G11, VPOL, VRS, SHEAR
    use numerical_tools, only: qinterp
    use debugger, only: markloc

    integer, parameter :: n_arr_out=15, nrho_m=64, nworkers=64, n_dims=5, &
         n_scalars=8, n_inputs=39

    double precision, intent(in), optional :: rho_norm_max

    logical :: first_call=.True.
    integer :: nchunk
    integer :: i, j, jr, jrho, jr_r, jr_l, jgamma_max, jspec
    integer :: ns_in              ! Number of species, including electrons
    integer :: t_wall1, t_wall2, rate, max_nworkers, stat
    integer, dimension(n_dims) :: dims_in

    double precision, dimension(n_inputs, nrho_m) :: prof_in
    double precision, dimension(n_arr_out, nrho_m) :: prof_out
    double precision :: bpolz, xstep, rho_min, rho_max, dstep, a0_m
    double precision, dimension(n_scalars) :: scal_in
    double precision, dimension(nrho_m) :: drmin, drmaj, drho, dte, dne, dq, dptot, &
        dvper, drhodr, dr, dv_r
    double precision, dimension(NRD) :: gradrhosq_as, rmaj_as, q_as, &
        vpar_as, vper_as, &
        chie_m, chii_m, elec_pflux_m, exchi_m, ptot_as, gamma_m, omega_m
    double precision, dimension(nrho_m) :: chie, chii, exchi, elec_pflux, rho_m, &
        gamma_max, omega_max, kymax, te_m, ne_m, vpar_m, vper_m, &
        ametr_m, rmaj_m, ptot_m, q_m
    double precision, dimension(nspec_max) :: zs_in
    double precision, dimension(nspec_max-1, nrho_m) :: dti, dni, ni_m, ti_m, ion_pflux
    double precision, dimension(nspec_max-1, nrho_m) :: zi_m 
    double precision, dimension(nspec_max-1, NRD) :: ni_as, ion_pflux_m
    character(len=32) :: str_nworkers
    character(len=64), dimension(:), allocatable :: SBP_NAMES

    call SYSTEM_CLOCK(t_wall1, rate)

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
            write(*, '(A, i3, A, i3)') '>>> WARNING n_workers=', nworkers, ' exceeds #physical CPUs=', max_nworkers
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
        allocate(SBP_NAMES(nworkers))
    endif

    nchunk = nrho_m / nworkers

    ! Interpolate from ASTRA grid to QuaLiKiZ grid
    rho_min = RHO(1)
    if (present(rho_norm_max)) then
        rho_max = rho_norm_max*ROC
    else
        rho_max = ROC
    endif
    xstep = (rho_max - rho_min)/(nrho_m - 1.)
    rho_m = (/ (rho_min + (jr - 1.)*xstep, jr=1, nrho_m) /)

    zi_m(1, :) = ZMJ
    call qinterp(RHO(1:NA1),     TI(1:NA1), NA1, rho_m, ti_m(1, :), nrho_m)
    call qinterp(RHO(1:NA1),     TE(1:NA1), NA1, rho_m,       te_m, nrho_m)
    call qinterp(RHO(1:NA1),   ZIM1(1:NA1), NA1, rho_m, zi_m(2, :), nrho_m)
    call qinterp(RHO(1:NA1),   ZIM2(1:NA1), NA1, rho_m, zi_m(3, :), nrho_m)
    call qinterp(RHO(1:NA1),   ZIM3(1:NA1), NA1, rho_m, zi_m(4, :), nrho_m)
    call qinterp(RHO(1:NA1),     NE(1:NA1), NA1, rho_m,       ne_m, nrho_m)
    call qinterp(RHO(1:NA1),  AMETR(1:NA1), NA1, rho_m,    ametr_m, nrho_m)

    ti_m(2, :) = ti_m(1, :)
    ti_m(3, :) = ti_m(1, :)
    ti_m(4, :) = ti_m(1, :)

    do jrho=1, NA1
        if (NDEUT(jrho) >= 0.01*NE(jrho)) then
            ni_as(1, jrho) = NDEUT(jrho)
        else ! likely: NDEUT not defined in equ file, hence zero
            ni_as(1, jrho) = NI(jrho)
        endif
        ni_as(2, jrho) = NIZ1(jrho)
        ni_as(3, jrho) = NIZ2(jrho)
        ni_as(4, jrho) = NIZ3(jrho)
        rmaj_as(jrho) = RTOR + SHIF(jrho)
        q_as(jrho)    = 1./MU(jrho)
        ptot_as(jrho) = NE(jrho)*TE(jrho) + ni_as(1, jrho)*TI(jrho) + ni_as(2, jrho)*TI(jrho) + pfast(jrho) + 0.5*(pblon(jrho) + pbper(jrho))
        bpolz = BTOR*AMETR(jrho)*MU(jrho)/RTOR
        gradrhosq_as(jrho) = G11(jrho)/VRS(jrho)
        vper_as(jrho) = ER(jrho)/(RTOR*bpolz) ! vexb in m/s --> Omega_E
        vpar_as(jrho) = ER(jrho)/(RTOR*bpolz)*(RTOR+SHIF(jrho) + AMETR(jrho))
    enddo

    call qinterp(RHO(1:NA1), ni_as(1, 1:NA1), NA1, rho_m, ni_m(1, :), nrho_m)
    call qinterp(RHO(1:NA1), ni_as(2, 1:NA1), NA1, rho_m, ni_m(2, :), nrho_m)
    call qinterp(RHO(1:NA1), ni_as(3, 1:NA1), NA1, rho_m, ni_m(3, :), nrho_m)
    call qinterp(RHO(1:NA1), ni_as(4, 1:NA1), NA1, rho_m, ni_m(4, :), nrho_m)
    call qinterp(RHO(1:NA1),  rmaj_as(1:NA1), NA1, rho_m,  rmaj_m, nrho_m)
    call qinterp(RHO(1:NA1),     q_as(1:NA1), NA1, rho_m,     q_m, nrho_m)
    call qinterp(RHO(1:NA1),  ptot_as(1:NA1), NA1, rho_m,  ptot_m, nrho_m)
    call qinterp(RHO(1:NA1),  vpar_as(1:NA1), NA1, rho_m,  vpar_m, nrho_m)
    call qinterp(RHO(1:NA1),  vper_as(1:NA1), NA1, rho_m,  vper_m, nrho_m)

    ! Reference length
    a0_m = AMETR(NA1)

    do jr=1, nrho_m
        ni_m(2, jr) = max(1.e-9, ni_m(2, jr))
        ni_m(3, jr) = max(1.e-9, ni_m(3, jr))
        ni_m(4, jr) = max(1.e-9, ni_m(4, jr))
    enddo

    elec_pflux_m = 0.
    ion_pflux_m  = 0.
    chie_m  = 0.
    chii_m  = 0.
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

    !--------------
    ! Differentials

    dti = 0.d0
    dni = 0.d0
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
        dptot(jr)  = dstep*( ptot_m(jr_r) -  ptot_m(jr_l))
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
    ! IPC parallelisation
    !--------------------

    nchunk = nrho_m / nworkers

    dims_in(1) = nchunk
    dims_in(2) = n_inputs
    dims_in(3) = n_arr_out
    dims_in(4) = nrho_m
    dims_in(5) = ns_in

    scal_in(1) = BTOR
    scal_in(2) = RTOR
    scal_in(3) = AMETR(NA1)
    scal_in(4) = AMJ
    scal_in(5) = AIM1
    scal_in(6) = AIM2
    scal_in(7) = AIM3
    scal_in(8) = ZMJ

    prof_in( 1, :) = rho_m
    prof_in( 2, :) = ametr_m
    prof_in( 3, :) = rmaj_m
    prof_in( 4, :) = q_m
    prof_in( 5, :) = ne_m
    prof_in( 6, :) = te_m
    prof_in( 7, :) = vpar_m
    prof_in( 8, :) = vper_m
    prof_in( 9, :) = ti_m(1, :)
    prof_in(10, :) = ti_m(2, :)
    prof_in(11, :) = ti_m(3, :)
    prof_in(12, :) = ti_m(4, :)
    prof_in(13, :) = ni_m(1, :)
    prof_in(14, :) = ni_m(2, :)
    prof_in(15, :) = ni_m(3, :)
    prof_in(16, :) = ni_m(4, :)
    prof_in(17, :) = zi_m(1, :)
    prof_in(18, :) = zi_m(2, :)
    prof_in(19, :) = zi_m(3, :)
    prof_in(20, :) = zi_m(4, :)
    prof_in(21, :) = drmin
    prof_in(22, :) = drmaj
    prof_in(23, :) = drho
    prof_in(24, :) = dptot
    prof_in(25, :) = dte
    prof_in(26, :) = dne
    prof_in(27, :) = dq
    prof_in(28, :) = dvper
    prof_in(29, :) = dv_r
    prof_in(30, :) = dr
    prof_in(31, :) = drhodr
    prof_in(32, :) = dti(1, :)
    prof_in(33, :) = dti(2, :)
    prof_in(34, :) = dti(3, :)
    prof_in(35, :) = dti(4, :)
    prof_in(36, :) = dni(1, :)
    prof_in(37, :) = dni(2, :)
    prof_in(38, :) = dni(3, :)
    prof_in(39, :) = dni(4, :)

    if (first_call) then
        SBP_NAMES = "xpr/qlki"//char(0)
        call initialise_ipc(nrho_m, n_dims, n_scalars, n_inputs, n_arr_out, nworkers, equ_file, exp_file)
        call fill_dim2shm(dims_in)
        call send_ipc_jobs(nworkers, nchunk, 64, SBP_NAMES)
        first_call = .False.
    endif

    ! **** Fill shared memory segments
    call fill_var2shm(scal_in)
    call fill_arr2shm(prof_in)

    ! **** Free each semaphore
    do j=1, nworkers
        call unlock_sbp(j)
    enddo

    ! **** Synchronisation point
    call wait4all

    ! **** Collect data from ShMem
    do i=1, nworkers
        call sbp2astra(i, prof_out(1, 1))
    enddo

    ! Interpolate back to ASTRA radial grid

    call qinterp(rho_m, prof_out(1, :), nrho_m, RHO(1:NA1), chii_m(1:NA1)      , NA1, extrap_right=0.)
    call qinterp(rho_m, prof_out(2, :), nrho_m, RHO(1:NA1), chie_m(1:NA1)      , NA1, extrap_right=0.)
    call qinterp(rho_m, prof_out(4, :), nrho_m, RHO(1:NA1), elec_pflux_m(1:NA1), NA1, extrap_right=0.)
    call qinterp(rho_m, prof_out(5, :), nrho_m, RHO(1:NA1), exchi_m(1:NA1)     , NA1, extrap_right=0.)

    chii_m (1:2) = chii_m (3)
    chie_m (1:2) = chie_m (3)
    elec_pflux_m(1:2) = elec_pflux_m(3)
    exchi_m(1:2) = exchi_m(3)

    do jrho=1, NA1
        qlk_out%chi_i(jrho)    = chii_m(jrho)/gradrhosq_as(jrho) ! \chi_i, m^2/s
        qlk_out%chi_e(jrho)    = chie_m(jrho)/gradrhosq_as(jrho) ! \chi_e, m^2/s
        qlk_out%e_pflux(jrho)  = elec_pflux_m(jrho)/a0_m/gradrhosq_as(jrho) ! D flux
        qlk_out%equipart(jrho) = exchi_m(jrho) ! turbulent e-i equipartition in MW/m^3
    enddo

    call SYSTEM_CLOCK(t_wall2, rate)
    print*, "XPR wall time", dble(t_wall2 - t_wall1)/dble(rate)

    return
    end subroutine qlk_ipc

end module a2qlk
