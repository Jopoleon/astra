subroutine SETARX(ICALL)
!--------------------------------------------------------------------
! The time evolution of the input data is taken from
!  - raw_profiles  |  for arrays
!
! Then it is stored for the current time in the arrays
! profiles_x(NRD, NARRX) - (description in astra_variables.json)
!   
!--------------------------------------------------------------------

use parameter_inc, only: NRD, NRDX, NTARR
use const_inc, only: TIME, BTOR, GP, AB, ABC, ROC, VOLUME, NA1, NAB, PSIAX
use status_inc, only: AMETR, RHO, FP, VOLUM, profiles_x, rho_pol
use numerical_tools, only: qinterp, sortab, smooth
use io_mod, only: jbeg_arrx, IFDFAX, XAXES, DATAX, NPTM, TOUTX
use debugger, only: markloc, astra_stop
use exp_data, only: raw_profiles

implicit  none

integer, intent(in) :: ICALL

integer :: jprof, n_grid, gridtype, jtarr, jt_start, jto, jt_end, kn, j3, &
    jtn, jt, jt0, jx, jy, NP1, N11
double precision :: RORZ, RZ2A, YDT, YDTA, YDTB, Y, Y1, dxl, dxr
double precision, dimension(NRDX) :: x_grid, dat_exp
double precision, dimension(NRD) :: XA, DA
character(len=132) :: err_msg, err_msg_grid
!--------------------------------------------------------------------
!  NARRX    maximal number of arrays readable from a data file 
!  NTARR    maximal number of time slices for all arrays (total)
!--------------------------------------------------------------------
! Input
! ICALL
!  = 1 - time interpolation off
!  = 2 - time interpolation on
! raw_profiles%data(NRDX*NTARR) - data array
! jbeg_arrx(kn)  - pointer to a position in the array TIMEX
! Output
! IFDFAX(kn)   - current pointer to data set in raw_profiles%data
! NPTM(kn)     - number of data points within a<=AB
! XAXES(jprof, kn) - "radial" grid for displayed data
! DATAX(jprof, kn) - array for displayed data
! profiles_x(jprof, kn)   - smoothed curve
!--------------------------------------------------------------------

call markloc('SETARX')

var_loop: do jtarr=1, NTARR
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

        n_grid   = raw_profiles%nrho(jt)
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
            if (x_grid(N11)  < dxl*AB) N11 = n_grid + 1
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
            if (x_grid(N11)  < dxl*AB) N11 = n_grid+1
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

return
end subroutine SETARX
