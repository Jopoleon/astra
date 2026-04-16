module set_x_data

implicit none

contains

!---------------------------------------------------------------------
    subroutine get_coil(tim_in, raw_coil, coil_curr)

    use read_input, only: rawCoils 

! Get the control quantities from the exp data at the present time slice

    double precision, intent(in) :: tim_in
    type(rawCoils), intent(in) :: raw_coil
    double precision, intent(out), dimension(raw_coil%ncoils) :: coil_curr

    integer :: j, j1, j2, jt, nt, n_coil
    double precision :: ydt, yd1, yd2

    nt = raw_coil%nt
    n_coil = raw_coil%ncoils
    if (nt == 0) then  
        coil_curr = 0.0
        return
    endif

! Input order:
! t_1  t_2  t_3
! coil currents(t_1)  
! coil currents(t_2) 
    if (tim_in <= raw_coil%time(1)) then
        do j=1, n_coil
            coil_curr(j) = raw_coil%current(j)
        enddo
        return
    endif
    if (tim_in >= raw_coil%time(nt)) then
        j1 = n_coil*(nt - 1)
        do j=1, n_coil
            coil_curr(j) = raw_coil%current(j1+j)
        enddo
        return
    endif

    do j=1, nt
        if (raw_coil%time(j) <= tim_in) jt = j
    enddo
    ydt = raw_coil%time(jt+1) - raw_coil%time(jt)
    yd1 = (tim_in - raw_coil%time(jt))/ydt  ! > 0
    yd2 = (tim_in - raw_coil%time(jt+1))/ydt  ! < 0
    do j=1, n_coil
        j1 = n_coil*(jt-1) + j
        j2 = n_coil*jt + j
        coil_curr(j) = yd1*raw_coil%current(j2) - yd2*raw_coil%current(j1)
    enddo

    end subroutine get_coil

!--------------------------------------------------------------------
    subroutine set_x_arrays(ICALL)
!--------------------------------------------------------------------
! The time evolution of the input data is taken from
!  - raw_profiles  |  for arrays
!
! Then it is stored for the current time in the arrays
!--------------------------------------------------------------------

    use pi_const, only: GP
    use scalars, only: TIME, BTOR, AB, ABC, ROC, VOLUME, NA1, NAB, PSIAX
    use status, only: NRD, AMETR, RHO, FP, VOLUM, profiles_x, rho_pol
    use numerical_tools, only: qinterp, sortab, smooth
    use debugger, only: markloc, astra_stop
    use read_input, only: raw_profiles, jbeg_arrx, IFDFAX, &
        XAXES, DATAX, NPTM, TOUTX
    use standard_functions, only: RZ2A

    integer, intent(in) :: ICALL

    integer :: jprof, n_grid, gridtype, jtarr, jt_start, jto, jt_end, kn, j3, &
        jtn, jt, jt0, jx, jy, NP1, N11
    double precision :: RORZ, YDT, YDTA, YDTB, Y, Y1, dxl, dxr
    double precision, dimension(:), allocatable :: x_grid, dat_exp
    double precision, dimension(NRD) :: XA, DA
    character(len=132) :: err_msg, err_msg_grid
!--------------------------------------------------------------------
! Input
! ICALL
!  = 1 - time interpolation off
!  = 2 - time interpolation on
! raw_profiles%data - data array
! jbeg_arrx(kn)  - pointer to a position in the array TIMEX
! Output
! IFDFAX(kn)   - current pointer to data set in raw_profiles%data
! NPTM(kn)     - number of data points within a<=AB
! XAXES(jprof, kn) - "radial" grid for displayed data
! DATAX(jprof, kn) - array for displayed data
! profiles_x(jprof, kn)   - smoothed curve
!--------------------------------------------------------------------

    call markloc('set_x_arrays')

    var_loop: do jtarr=1, raw_profiles%n_groups
        if (raw_profiles%arr_index(jtarr) == 0) EXIT
        jprof = 0
        if (raw_profiles%arr_index(jtarr+1) == 0) then
            jprof  = jtarr
            jt_end = jtarr
        else
            if (raw_profiles%arr_index(jtarr+1) /= raw_profiles%arr_index(jtarr)) then ! Label change
                jprof = jtarr ! jprof -> group end
                jt_end = jbeg_arrx(raw_profiles%arr_index(jtarr+1)) - 1 ! jt_end -> last time
            endif
        endif
        if (jprof == 0) CYCLE var_loop

        KN = raw_profiles%arr_index(jprof)
        err_msg = 'Quantity  ' // raw_profiles%label(jprof) // ' Input times '
        jt_start = jbeg_arrx(KN)
        jto = jt_start
! Check whether the input time-array is monotonic
        do j3=jt_start+1, jt_end
            if (raw_profiles%time(j3)  < raw_profiles%time(j3-1)) call astra_stop(err_msg // 'not ascending')
            if (raw_profiles%time(j3) == raw_profiles%time(j3-1)) call astra_stop(err_msg // 'repeated')
            if (raw_profiles%time(j3) <= TIME) jto = j3
        enddo

        IFDFAX(KN) = jto

        if (TIME <= raw_profiles%time(jt_start)) then
            jtn = jt_start
        else
            jtn = min(jto+1, jt_end)
        endif

        if (2.*TIME > raw_profiles%time(jtn) + raw_profiles%time(jto)) then
            if (ICALL == 1) then
                jt = jtn
            else
                jt = jto
            endif
        else
            if (ICALL == 1) then
                jt = jto
            else
                jt = jtn
            endif
        endif

        jt0 = 0

        time_loop: do

!--------------------------------------------------------------------
! The following is done below:
! (1) The grid in "a", XA(NP1), and the data DA(NP1) on this grid
!      are defined by
!  (i)  mapping the original grid to the "a" grid x_grid(N11)
!  (ii) transfer (SMOOTH) from {x_grid(N11), dat_exp(N11)} to {XA, DA}
! (2) XAXES(n_grid, KN) is defined which is as "a" grid for exp-dot plots
!      DATAX(n_grid, KN) data on this grid
! (3) profiles_x(NRD, KN) smoothed input arrays interpolated in time

            n_grid = raw_profiles%nrho(jt)
            if (allocated(x_grid)) then
                deallocate(x_grid)
                deallocate(dat_exp)
            endif
            allocate(x_grid(n_grid+1), dat_exp(n_grid+1))

            gridtype = raw_profiles%grid_type(jt)
            write(err_msg_grid, '(A, i, A)')  'Option GRIDTYPE=', gridtype, ' not implemented, exiting'
            jx = raw_profiles%jbeg_grid(jt)
            jy = raw_profiles%jbeg_data(jt)

            N11 = n_grid
            dxl = 1. - 0.5/n_grid
            dxr = 1. + 0.5/n_grid

            do j3=1, n_grid
                dat_exp(j3) = raw_profiles%data(jy + j3 - 1)
                DATAX(j3, KN) = dat_exp(j3)
            enddo
            if (gridtype < 10) then
                do j3=1, n_grid
                    x_grid(j3) = (j3 - 1.)/(n_grid - 1.)
                enddo
            else if (gridtype < 20) then
                if (gridtype == 18 .or. gridtype == 19) then
                    RORZ = raw_profiles%data(jx)
                    jx = jx + 1
                endif
                x_grid(: n_grid) = raw_profiles%data(jx: jx + n_grid - 1)
                XAXES(: n_grid, KN) = x_grid(: n_grid)
            else if (gridtype > 20) then
                CYCLE var_loop
            endif

            SELECT CASE(gridtype)

            CASE(0)
                NP1 = NAB
                XA(: NP1) = AMETR(: NP1)/AB
                XAXES(: n_grid, KN) = AB*x_grid(: n_grid)

            CASE(1)
                NP1 = NA1
                XA(: NP1) = AMETR(: NP1)/ABC
                XAXES(: n_grid, KN) = ABC*x_grid(: n_grid)

            CASE(2)
                NP1 = NA1
                XA(: NP1) = RHO(: NP1)/ROC
                call qinterp(XA(1: NP1), AMETR(1: NP1), NP1, x_grid(1: n_grid), XAXES(1: n_grid, KN), n_grid)
! AMETR(XA(1:NP1)) is given; AMETR(x_grid(j)) is returned;

            CASE(3)
                NP1 = NA1
                XA(: NP1) = rho_pol(: NP1)
                call qinterp(XA, AMETR, NP1, x_grid(1: n_grid), XAXES(1: n_grid, KN), n_grid)

            CASE(4)
                call astra_stop(err_msg_grid)

            CASE(5)
                call astra_stop(err_msg_grid)

            CASE(6)
                call astra_stop(err_msg_grid)

! Normalised grids

            CASE(10)
                NP1 = NAB
                XA(: NP1) = AMETR(: NP1)/AB
                do N11=n_grid, 1, -1
                    if (x_grid(N11) <= dxr*AB) EXIT
                enddo
                if (x_grid(N11) < dxl*AB) N11 = n_grid + 1
                x_grid(: N11-1) = x_grid(: N11-1)/AB
                call qinterp(XA(1: NP1), AMETR(1: NP1), NP1, x_grid(1: n_grid), XAXES(1: n_grid, KN), n_grid)
                x_grid(N11) = 1.
                dat_exp(N11) = DATAX(min(n_grid, N11), KN)

            CASE(11)
                NP1 = NA1
                XA(: NP1) = AMETR(: NP1)/ABC
                do N11=n_grid, 1, -1
                    if (x_grid(N11) <= dxr*ABC) EXIT
                enddo
                if (x_grid(N11)  < dxl*ABC) N11 = n_grid+1
                x_grid(: N11-1) = x_grid(: N11-1)/ABC
                call qinterp(XA(1: NP1), AMETR(1: NP1), NP1, x_grid(1: n_grid), XAXES(1: n_grid, KN), n_grid)
                x_grid(N11) = 1.
                dat_exp(N11) = DATAX(min(n_grid, N11), KN)

            CASE(12)
                NP1 = NA1
                XA(: NP1) = RHO(: NP1)/ROC
                XA = XA / XA(NP1)    ! renormalize in case rho(NP1) is not equal to ROC (can happen if equilibrium is not consistent and/or not calculated at every time step)
                do N11=n_grid, 1, -1
                    if (x_grid(N11) <= dxr) EXIT
                enddo
                if (x_grid(N11) < dxl) N11 = n_grid+1
                call qinterp(XA(1: NP1), AMETR(1: NP1), NP1, x_grid(1: n_grid), XAXES(1: n_grid, KN), n_grid)
                x_grid(N11) = 1.
                dat_exp(N11) = DATAX(min(n_grid, N11), KN)

            CASE(13)
                NP1 = NA1
                XA(: NP1) = rho_pol(: NP1)
                do N11=n_grid, 1, -1
                    if (x_grid(N11) <= dxr) EXIT
                enddo
                if (x_grid(N11) < dxl) N11 = n_grid+1
                call qinterp(XA(1: NP1), AMETR(1: NP1), NP1, x_grid(1: n_grid), XAXES(1: n_grid, KN), n_grid)
                x_grid(N11) = 1.
                dat_exp(N11) = DATAX(min(n_grid, N11), KN)

            CASE(14)
                NP1 = NA1
                XA(: NP1) = sqrt((VOLUM(: NP1) - VOLUM(1))/VOLUME)
                do N11=n_grid, 1, -1
                    if (x_grid(N11) <= dxr) EXIT
                enddo
                if (x_grid(N11) < dxl) N11 = n_grid+1
                call qinterp(XA(1: NP1), AMETR(1: NP1), NP1, x_grid(1: n_grid), XAXES(1: n_grid, KN), n_grid)
                x_grid(N11) = 1.
                dat_exp(N11) = DATAX(min(n_grid, N11), KN)

            CASE(15)
                NP1 = NA1
                XA(: NP1) = (FP(: NP1) - PSIAX)/(FP(j3) - FP(NP1))
                do N11=n_grid, 1, -1
                    if (x_grid(N11) <= dxr) EXIT
                enddo
                if (x_grid(N11) < dxl) N11 = n_grid+1
                call qinterp(XA(1: NP1), AMETR(1: NP1), NP1, x_grid(1: n_grid), XAXES(1: n_grid, KN), n_grid)
                x_grid(N11) = 1.
                dat_exp(N11) = DATAX(min(n_grid, N11), KN)

            CASE(16)
                NP1 = NA1
                XA(: NP1) = RHO(: NP1)/ROC
                Y = 1./(GP*BTOR)
                x_grid(: n_grid) = sqrt(x_grid(: n_grid)*Y)
                do N11=n_grid, 1, -1
                    if (x_grid(N11) <= ROC*dxr) EXIT
                enddo
                if (x_grid(N11) < ROC*dxl) N11 = n_grid+1
                x_grid(: N11-1) = x_grid(: N11-1)/ROC
                call qinterp(XA(1: NP1), AMETR(1: NP1), NP1, x_grid(1: n_grid), XAXES(1: n_grid, KN), n_grid)
                x_grid(N11) = 1.
                dat_exp(N11) = DATAX(min(n_grid, N11), KN)

            CASE(17)
                call astra_stop(err_msg_grid)

            CASE(18)
                NP1 = NA1
                do j3=1, n_grid
                    XAXES(j3, KN) = RZ2A(RORZ, x_grid(j3), NP1)
                    x_grid(j3) = XAXES(j3, KN)
                enddo
                call SORTAB(x_grid, dat_exp, n_grid)
                do N11=n_grid, 1, -1
                    if (x_grid(N11) <= dxr*AB) EXIT
                enddo
                if (x_grid(N11)  < dxl*AB) N11 = n_grid+1
                dat_exp(N11) = dat_exp(min(n_grid, N11))
                x_grid(: N11-1) = x_grid(: N11-1)/AB
                x_grid(N11) = 1.
                NP1 = NAB
                do j3 = 1, NP1
                    XA(j3) = AMETR(j3)/AB
                enddo

            CASE(19)
                NP1 = NA1
                do j3=1, n_grid
                    XAXES(j3, KN) = RZ2A(x_grid(j3), RORZ, NP1)
                    x_grid(j3) = XAXES(j3, KN)
                enddo
                call SORTAB(x_grid, dat_exp, n_grid)
                do N11=n_grid, 1, -1
                    if (x_grid(N11) <= dxr*AB) EXIT
                enddo
                if (x_grid(N11) < dxl*AB) N11 = n_grid+1
                dat_exp(N11) = dat_exp(min(n_grid, N11))
                x_grid(: N11-1) = x_grid(: N11-1)/AB
                x_grid(N11) = 1.
                NP1 = NAB
                do j3 = 1, NP1
                    XA(j3) = AMETR(j3)/AB
                enddo

            CASE(20)
! Input data are given on {r, z} plane
                NP1 = NAB
                do j3=1, n_grid
                    Y  = raw_profiles%data(jx + j3 - 1)
                    Y1 = raw_profiles%data(jx + n_grid + j3 - 1)
                    XAXES(j3, KN) = RZ2A(Y, Y1, NP1)
                    DATAX(j3, KN) = raw_profiles%data(jy + j3 - 1)
                    x_grid(j3) = XAXES(j3, KN)
                    dat_exp(j3) = raw_profiles%data(jy + j3 - 1)
                enddo

            END SELECT

            NPTM(KN) = min(n_grid, N11)
            TOUTX(KN) = raw_profiles%time(jt)

! This is added to avoid too long extrapolation to the magnetic axis
            if ( N11 > 1 ) then
                if ( x_grid(2) - x_grid(1) < x_grid(1) ) x_grid(1) = 0.
            endif

! All input data are mapped to the grid XA(1:NP1) in the variable "a"

            call SMOOTH(raw_profiles%filter(jt), x_grid(1:N11), dat_exp(1:N11), N11, XA(1:NP1), DA(1:NP1), NP1)

!     data interpolation
! jto - pointer to the previous time
! jtn - pointer to the subsequent time
! jt  - pointer to the current time
! if jto=/=jtn two runs are accomplished: with jt=jtn and jt=jto
! the order of runs depends on the time selected for exp output
! namely, the last run determines current XAXES and DATAX

            if (jto == jtn) then ! no time dependence
                profiles_x(: NRD, KN) = DA(: NRD)
                CYCLE var_loop
            endif

            if (jt0 /= 0) EXIT time_loop

            jt0 = 1
            if (jt /= jto) then
                jt = jto
            else
                jt = jtn
            endif
            profiles_x(: NRD, KN) = DA(: NRD)
        enddo time_loop

        ydt  = (raw_profiles%time(jtn) - raw_profiles%time(jto))
        ydta = (raw_profiles%time(jtn) - TIME)/ydt
        ydtb = (TIME - raw_profiles%time(jto))/ydt
        if (jt == jto) then
            do j3=1, NRD
                profiles_x(j3, KN) = profiles_x(j3, KN)*ydtb + DA(j3)*ydta
            enddo
        else
            do j3=1, NRD
                profiles_x(j3, KN) = profiles_x(j3, KN)*ydta + DA(j3)*ydtb
            enddo
        endif

    enddo var_loop

    if (allocated(x_grid)) then
        deallocate(x_grid)
        deallocate(dat_exp)
    endif

    end subroutine set_x_arrays

!-----------------------------------------------------------------------
    subroutine set_x_scalars()
!-----------------------------------------------------------------------
! Time evolution of the scalar input data
! For the current time, a value is stored in the array
! varxValues(NTVAR) - (description in the file astra_variables.json)
!
! Output:
!    varxValues, varValues
!
! IFDFVX(N)  - type of variable
! IFDFVX:
!    = 0 - determined by the data file (independent on time),
!    = 1 - determined by the data file (dependent on time),
!    = 2 - determined by the MODEL,
!    = 3 - Keyboard
!    = 4 - any change of the variable is forbidden 
!          eg. (AB, RTOR, ELONM, TRICH or set interactively)
!-----------------------------------------------------------------------

    use scalars, only: varxValues, varValues, TIME
    use read_input, only: raw_scalars, IFDFVX
    use debugger, only: markloc

    integer :: jtvar, N1, N2
    double precision :: ydt, ydtr, ydtl

    call markloc('set_x_scalars')

    N1 = 0
    N2 = 0

    do jtvar=1, raw_scalars%nt_all
        if (raw_scalars%var_index(jtvar) == 0) EXIT
        N2 = N1
        N1 = raw_scalars%var_index(jtvar)
        if (IFDFVX(N1) >= 0) then
            if (IFDFVX(N1) == 0 .or. N1 /= N2) varxValues(N1) = raw_scalars%data(jtvar)
            if (N1 == N2) then
                if (TIME >= raw_scalars%time(jtvar-1)) then
                    if (TIME < raw_scalars%time(jtvar)) then
                        YDT = raw_scalars%time(jtvar) - raw_scalars%time(jtvar-1)
                        YDTR = (raw_scalars%time(jtvar) - TIME)/YDT
                        YDTL = (TIME - raw_scalars%time(jtvar-1))/YDT
                        varxValues(N1) = raw_scalars%data(jtvar)*YDTL + raw_scalars%data(jtvar-1)*YDTR
                    else
                        varxValues(N1) = raw_scalars%data(jtvar)
                    endif
                endif
            endif
            if (IFDFVX(N1) <= 1) varValues(N1) = varxValues(N1)
        endif
    enddo

    end subroutine set_x_scalars

!----------------------------------------------------------------------
    subroutine astra_assignments()

    use pi_const, only: GP, GP2
    use scalars, only: NA1, NA, NB1, NAB, AB, ABC, AWALL, TIME, TSTART, TPAUSE, &
         TAUMIN, TAUPRP, VOLUME, IPL, IPLN, HRO, HROX, ROC, ROCO, BTN, FTO, FTN, &
         PSIAX, PSIBO, RTOR, BTOR, SHIFT, ROWALL, ELONM, ELONG, TRICH, TRIAN, &
         TINIT, TSCALE, TIMEQL, DTEQL
    use status, only: NRD, XRHO, SXHO, RHO, SRHO, AMETR, &
        G11, G22, VR, VRO, VRS, VOLUM, &
        FP, FPO, FP_NORM, rho_pol, NE, NEO, TE, TEO, UPAR, UPARO, MRHO, &
        AMAIN, UPS0, UPS0O
    use json_vars, only: varNames, n_profx, profxNames, n_var
    use numerical_tools, only: EXTRAP, INTEGR
    use debugger, only: astra_stop
    use read_input, only: raw_boundary, exp_file, IFDFVX
    use metrics, only: setgeo, new_grid, roc3a
    use char_manip, only: to_upper, str_in_list
 
    integer :: j, jt, jthe, jvar, KAWALL
    double precision :: YTP=-1.d9
    double precision, dimension(:), allocatable :: bnd_r, bnd_z
    character(len=132) :: err_msg

    if (NA1 > NRD) then
        write(err_msg, '(2A, i)') '>>> FATAL ERROR: The radial grid size out of range.\n', &
            '                 Parameter "NA1" cannot exceed', NRD
        call astra_stop(err_msg)
    endif

    TIME = TSTART
    if (YTP > -1.d8) TPAUSE = YTP

    call set_x_scalars()
    do j=1, n_var
        jvar = str_in_list(varNames(j), (/ 'AB    ', 'AWALL ', 'RTOR  ', 'ELONM ', 'TRICH ' /))
        if (jvar > 0) IFDFVX(j) = 4
        if (varNames(j) == 'AWALL ') KAWALL = j
    enddo
 
    if (AWALL < AB) then
         if (IFDFVX(KAWALL) >= 0) write(*, *) '>>> Warning: AWALL < AB.  Setting AWALL = AB'
         AWALL = AB
    endif
    if (IFDFVX(KAWALL) < 0 .and. (AWALL < AB .or. AWALL > 1.2*AB)) AWALL = AB

    if (ABC > AB) then
        err_msg = '>>> Error: ABC cannot exceed AB. Check your file exp/' // TRIM(exp_file)
        call astra_stop(TRIM(err_msg))
    endif

    if (AWALL > 1.2*AB .or. AWALL > RTOR) then
        write(*, *) '>>> Warning: AWALL is set unreasonably large'
        write(*, *) '    Check settings in data and log files'
    endif

! If boundary is given, calculates initial geometry from that
    if (raw_boundary%nt > 0) then
! Find time index of most proximum boundary
        j=1
        do jt=1, raw_boundary%nt
            if (raw_boundary%time(jt) <= TSTART) j = jt
        enddo
        jt = j
        allocate(bnd_r(raw_boundary%n_theta), bnd_z(raw_boundary%n_theta))
        do jthe=1, raw_boundary%n_theta
            bnd_r(jthe) = raw_boundary%R((jthe-1)*raw_boundary%nt + jt)
            bnd_z(jthe) = raw_boundary%Z((jthe-1)*raw_boundary%nt + jt)
        enddo
! Calculate ABC
        ABC = (maxval(bnd_r) - minval(bnd_r))/2.
! Calculate elong
        ELONG = (maxval(bnd_z) - minval(bnd_z))/(2.*ABC)
        ELONG = max(ELONG, 1.d0)
        deallocate(bnd_r)
        deallocate(bnd_z)
    endif

! Assign variables here for initialization:

    VOLUME = GP2*GP*RTOR*AB**2 * ELONG

    ROC  = ROC3A(RTOR, SHIFT, ABC, ELONG, TRIAN)
    ROCO = ROC

! Initialization of magnetic quantities
    BTN = BTOR
    FTO = GP*BTOR*ROC**2
    FTN = GP*BTN*ROC**2
    IPLN = IPL

    if (AWALL > RTOR+SHIFT) then
        ROWALL = AWALL*SQRT(max(ELONM, ELONG))
    elseif (ELONM > ELONG) then
        ROWALL = ROC3A(RTOR, SHIFT, AWALL, ELONM, TRICH)
    else
        ROWALL = ROC3A(RTOR, 0.d0, AWALL, ELONG, TRIAN)
    endif

    HROX  = 1.0/(NA1 - 0.5)

    do j=1, NRD
        XRHO(j) = (j - 0.5)*HROX
        SXHO(j) = j*HROX
! real space grids, these are function of ROC
        RHO (j) = XRHO(j)*ROC
        SRHO(j) = SXHO(j)*ROC
    enddo
    HRO  = HROX*ROC

! Compute NB1
    NB1 = NA1
    NA  = NA1 - 1

    call SETGEO()
    call NEW_GRID()

    do J=1, NB1
        G22(J) = RHO(J)
        VR(J)  = (GP2*(RTOR + SHIFT))**2*RHO(J)/RTOR
        VRS(J) = (GP2*(RTOR + SHIFT))**2*J*HRO/RTOR
        G11(J) = VRS(J)
    enddo

    call INTEGR(RHO, 1, VR, VOLUM, NA1)

    PSIBO = FP(NA1)
    PSIAX = EXTRAP(XRHO(1: NA1), FP(1: NA1), 0.d0, NA1, 2, .true.)

    FP_NORM = (FP - PSIAX)/(PSIBO - PSIAX)
    rho_pol = SQRT(FP_NORM)

    VOLUME = VOLUM(NA1)

    NAB = NA1
    if (NA1 < NB1 .and. AB > ABC) then
        do j=NA1+1, NB1
            if (AMETR(j) < AB) NAB = j
        enddo
        if (NAB < NB1) NAB = NAB + 1
    endif

    AMETR(NAB) = AB
    AMETR(NA1) = ABC

    NEO   = NE
    TEO   = TE
    FPO   = FP
    UPARO = UPAR
    VRO   = VR
    MRHO  = AMAIN*NE
    UPS0  = MRHO*RTOR
    UPS0O = UPS0

    TIMEQL = TIME - DTEQL - 1.d-7
    TAUPRP = TAUMIN
    if (TIME > TINIT + 1.025*abs(TSCALE)) TINIT = TSTART

    end subroutine astra_assignments

end module set_x_data
