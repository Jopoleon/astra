module a2neo

implicit none

integer, parameter :: nspec_max=7
double precision, allocatable, dimension(:, :) :: mem_neo

type neo_output
    double precision, allocatable, dimension(:) :: chi_i, chi_e, e_pflux
endtype neo_output
type(neo_output) :: neo_out 

contains

!---------------------------------------------------------------------
    subroutine neo_alloc

    use const_inc, only: NA1

    if (.not. allocated(neo_out%chi_i)) then
        allocate(neo_out%chi_i(NA1), neo_out%chi_e(NA1), neo_out%e_pflux(NA1))
    endif

    return
    end subroutine neo_alloc

!---------------------------------------------------------------------
    subroutine neo_ipc(rho_norm_max)

    use omp_lib
    use parameter_inc, only: NRD
    use io_mod, only: equ_file, exp_file, awd
    use const_inc, only: NA1, BTOR, RTOR, ROC, AMJ, AIM1, AIM2, AIM3, ZMJ
    use status_inc, only: NE, TE, NI, TI, ER, MU, FP_NORM, &
        ZIM1, ZIM2, ZIM3, NDEUT, NIZ1, NIZ2, NIZ3, &
        RHO, AMETR, SHIF, ELON, TRIA, VTOR, VPOL, G11, VRS
    use numerical_tools, only: qinterp
    use debugger, only: markloc

    integer, parameter :: n_arr_out=15, nrho_m=64, nworkers=64, &
         n_dims=5, n_scalars=8, n_inputs=33
    double precision, parameter :: c_vpol=1.d0

    double precision, intent(in), optional :: rho_norm_max

    logical :: first_call=.True.
    integer :: jr, jrho, jr_r, jr_l, jion
    integer :: ns_in              ! Number of species, including electrons
    integer :: i, j, nchunk, jout, offset
    integer :: t_wall1, t_wall2, rate, max_nworkers, stat
    integer, dimension(n_dims) :: dims_in

    double precision, dimension(n_inputs, nrho_m) :: prof_in
    double precision, dimension(n_arr_out, nrho_m) :: prof_out
    double precision :: bmod, bpolz, xstep, rho_min, rho_max, dstep, T0, m0, a0_m, gradrhosq_inv, drho
    double precision, dimension(n_scalars) :: scal_in
    double precision, dimension(NRD) :: rmaj_as, q_as, ni_main_as, vpar_as, &
        j_boot, e_pflux_as, chii_as, chie_as
    double precision, dimension(nrho_m) :: rho_m, ti_m, te_m, ne_m, vpar_m, &
        ametr_m, elon_m, tria_m, rmaj_m, q_m, &
        drmin, drmaj, dti, dte, dne, dq, delong, dtrian, dvpar, drhodr, dr
    double precision, dimension(nspec_max) :: zs_in
    double precision, dimension(nspec_max-1, nrho_m) :: dni, ni_m
    double precision, dimension(nspec_max-1, nrho_m) :: zi_m 
    character(len=32) :: str_nworkers
    character(len=64) :: SBP_NAME
    character(len=128) :: ipc_file, astra_task

    call SYSTEM_CLOCK(t_wall1, rate)

    ipc_file   = TRIM(awd) // '/tmp/' // TRIM(exp_file) // TRIM(equ_file) // '.ipc'
    astra_task = TRIM(awd) // '/bin/' // TRIM(equ_file) // '.exe'

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
    endif

! Interpolate from ASTRA grid to NEO grid
    rho_min = RHO(1)
    if (present(rho_norm_max)) then
        rho_max = rho_norm_max*ROC
    else
        rho_max = ROC
    endif
    xstep = (rho_max - rho_min)/(nrho_m - 1.)
    rho_m = (/ (rho_min + (jr - 1.)*xstep, jr=1, nrho_m) /)

    zi_m(1, :) = ZMJ
    call qinterp(RHO(1:NA1),  ZIM1(1:NA1), NA1, rho_m, zi_m(2, :), nrho_m)
    call qinterp(RHO(1:NA1),  ZIM2(1:NA1), NA1, rho_m, zi_m(3, :), nrho_m)
    call qinterp(RHO(1:NA1),  ZIM3(1:NA1), NA1, rho_m, zi_m(4, :), nrho_m)
    call qinterp(RHO(1:NA1),  NIZ1(1:NA1), NA1, rho_m, ni_m(2, :), nrho_m)
    call qinterp(RHO(1:NA1),  NIZ2(1:NA1), NA1, rho_m, ni_m(3, :), nrho_m)
    call qinterp(RHO(1:NA1),  NIZ3(1:NA1), NA1, rho_m, ni_m(4, :), nrho_m)
    call qinterp(RHO(1:NA1),    TI(1:NA1), NA1, rho_m,       ti_m, nrho_m)
    call qinterp(RHO(1:NA1),    TE(1:NA1), NA1, rho_m,       te_m, nrho_m)
    call qinterp(RHO(1:NA1),    NE(1:NA1), NA1, rho_m,       ne_m, nrho_m)
    call qinterp(RHO(1:NA1), AMETR(1:NA1), NA1, rho_m,    ametr_m, nrho_m)
    call qinterp(RHO(1:NA1),  ELON(1:NA1), NA1, rho_m,     elon_m, nrho_m)
    call qinterp(RHO(1:NA1),  TRIA(1:NA1), NA1, rho_m,     tria_m, nrho_m)

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
    call qinterp(RHO(1:NA1),    rmaj_as(1:NA1), NA1, rho_m,     rmaj_m, nrho_m)
    call qinterp(RHO(1:NA1),       q_as(1:NA1), NA1, rho_m,        q_m, nrho_m)
    call qinterp(RHO(1:NA1),    vpar_as(1:NA1), NA1, rho_m,     vpar_m, nrho_m)

! Reference length
    a0_m = AMETR(NA1)

    do jr=1, nrho_m
        ni_m(2, jr) = max(1.e-9, ni_m(2, jr))
        ni_m(3, jr) = max(1.e-9, ni_m(3, jr))
        ni_m(4, jr) = max(1.e-9, ni_m(4, jr))
    enddo

    e_pflux_as = 0.
    chie_as = 0.
    chii_as = 0.

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
        drho       = dstep*(  rho_m(jr_r) -   rho_m(jr_l))
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

    prof_in = 0.d0
    prof_in( 1, :) = rho_m
    prof_in( 2, :) = ametr_m
    prof_in( 3, :) = rmaj_m
    prof_in( 4, :) = elon_m
    prof_in( 5, :) = tria_m
    prof_in( 6, :) = q_m
    prof_in( 7, :) = ti_m
    prof_in( 8, :) = te_m
    prof_in( 9, :) = ne_m
    prof_in(10, :) = vpar_m
    prof_in(11, :) = ni_m(1, :)
    prof_in(12, :) = ni_m(2, :)
    prof_in(13, :) = ni_m(3, :)
    prof_in(14, :) = ni_m(4, :)
    prof_in(15, :) = zi_m(1, :)
    prof_in(16, :) = zi_m(2, :)
    prof_in(17, :) = zi_m(3, :)
    prof_in(18, :) = zi_m(4, :)
    prof_in(19, :) = drmin
    prof_in(20, :) = drmaj
    prof_in(21, :) = delong
    prof_in(22, :) = dtrian
    prof_in(23, :) = dti
    prof_in(24, :) = dte
    prof_in(25, :) = dne
    prof_in(26, :) = dq
    prof_in(27, :) = dvpar
    prof_in(28, :) = dr
    prof_in(29, :) = drhodr
    prof_in(30, :) = dni(1, :)
    prof_in(31, :) = dni(2, :)
    prof_in(32, :) = dni(3, :)
    prof_in(33, :) = dni(4, :)

    if (first_call) then
        SBP_NAME = "xpr/neo"//char(0)
        offset = 10
        call initialise_ipc(nrho_m, n_dims, n_scalars, n_inputs, n_arr_out, &
            nchunk, nworkers, offset, SBP_NAME, ipc_file, astra_task)
        call fill_dim2shm(n_dims, dims_in)
        first_call = .False.
    endif

! **** Fill shared memory segments
    call fill_var2shm(n_scalars, scal_in)
    call fill_arr2shm(nrho_m, n_inputs, prof_in)

! **** Free each semaphore
    do i=1, nworkers
        call unlock_sbp(i)
    enddo

! **** Synchronisation point
    call wait4all

! **** Collect data from ShMem
    do i=1, nworkers
        call sbp2astra(i, nchunk, n_arr_out, prof_out(1, 1))
    enddo

! Interpolate back to ASTRA radial grid

    call qinterp(rho_m, prof_out(1, :), nrho_m, RHO(1:NA1),    chii_as, NA1, extrap_right=0.)
    call qinterp(rho_m, prof_out(2, :), nrho_m, RHO(1:NA1),    chie_as, NA1, extrap_right=0.)
    call qinterp(rho_m, prof_out(4, :), nrho_m, RHO(1:NA1), e_pflux_as, NA1, extrap_right=0.)
    do jrho=1, NA1
        gradrhosq_inv = VRS(jrho)/G11(jrho)
        neo_out%chi_i(jrho)   = chii_as(jrho)*gradrhosq_inv ! m^2/s
        neo_out%chi_e(jrho)   = chie_as(jrho)*gradrhosq_inv ! m^2/s
        neo_out%e_pflux(jrho) = e_pflux_as(jrho)*gradrhosq_inv/a0_m
    enddo
    call SYSTEM_CLOCK(t_wall2, rate)
    print*, "XPR wall time", dble(t_wall2 - t_wall1)/dble(rate)

    return
    end subroutine neo_ipc

end module a2neo
