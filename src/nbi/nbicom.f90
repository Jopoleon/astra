module nbicom

use nbstatus, only: n_rho, n_theta, n_energy, nspec_max

implicit none

integer, parameter :: ndim1=82, ndim2=201
double precision, dimension(n_theta), parameter :: PLEJ2 = (/ &
    -5.000000d-01, -4.994000d-01, -4.976000d-01, -4.946000d-01, &
    -4.904000d-01, -4.850000d-01, -4.784000d-01, -4.706000d-01, &
    -4.616000d-01, -4.514000d-01, -4.400000d-01, -4.274000d-01, &
    -4.136000d-01, -3.986000d-01, -3.824000d-01, -3.650000d-01, &
    -3.464000d-01, -3.266000d-01, -3.056000d-01, -2.834000d-01, &
    -2.600000d-01, -2.354000d-01, -2.096000d-01, -1.826000d-01, &
    -1.544000d-01, -1.250000d-01, -9.440002d-02, -6.260002d-02, &
    -2.960002d-02,  4.599977d-03,  3.999998d-02,  7.659998d-02, &
     1.144000d-01,  1.534000d-01,  1.936000d-01,  2.350000d-01, &
     2.776000d-01,  3.214000d-01,  3.664000d-01,  4.126000d-01, &
     4.599999d-01,  5.085999d-01,  5.584000d-01,  6.094000d-01, &
     6.615999d-01,  7.150000d-01,  7.695999d-01,  8.253999d-01, &
     8.823999d-01,  9.405999d-01,  9.999999d-01/)

logical :: perp_flag
integer :: IHB, n1, ntet1, II(ndim2)
double precision :: aspect_ratio, amin_cm, HB, BZ, RBmin_cm, RBmax_cm, &
    EB, QB, PB, RJ, RJT, YDYH
double precision, dimension(n_energy) :: SVII, VNB, YAQB
double precision, dimension(n_rho) :: YDZ, YTRAP, X, XJ, DX, ELON1, TRIA1
double precision, dimension(ndim1) :: az, YDRY
double precision, dimension(4, 4) :: SVEX
double precision, dimension(n_energy, n_rho) :: ARD, F, YANBA, &
    YAQBA, YACBA, DTCX, YATBA, YANBA1, YSLEJ0, YSLEJ2
double precision, dimension(n_energy, n_rho, n_theta) :: YASBA, YASBA1
double precision, dimension(nspec_max) :: ZB, RMB
double precision :: YSIMPI(n_energy, nspec_max)
double precision, dimension(n_rho) :: stnbdp, sdnbtp, sdnbdp1, sdnbdp2

contains

!---------------------------------------------------------------------
    double precision function power_dist(x_in, C1, C2) 
!----- Vertical/horizontal beam power distribution in the footprint -----
! The input parameters X1 = CVER1/CHOR1, C2 = CVER2/CHOR2 are given in
!     the NBI configuration file for each NBI source (last line)
! x_in is the normalized height/width in respect of the footprint center
! H_half_width = CBMS4*(RBMAX-RBMIN)/2
! power_dist [a.u.] is renormalized in the calling routine
!---------------------------------------------------------------------

    double precision, intent(in) :: x_in, C1, C2

    power_dist = EXP(-C1*abs(x_in)**C2)

    end function power_dist

!---------------------------------------------------------------------
    subroutine nb0(jE_min, n_pencil, jsrc, adqb, yhm, ch1, ch2, cr1, cr2, cs3)

    integer, intent(in) :: jE_min, n_pencil, jsrc
    double precision, intent(in) :: yhm, ch1, ch2, cr1, cr2, cs3, adqb(n_energy)

    integer :: j, jh, jpencil, jhh, jE, n, jh1, jhb, jhb05, jedge
    double precision :: ys, z, zy, y1, zy12, yr1, yr2, yhm_n, yd, ye, &
        ydr, yr0, yz1, yds, ysb, yhbd2, yhocu, ynocu, width, coeff, &
        yhtop, yhtopa, yhedge, yhedga, yd05h, ydzout, ydhz, ydzjh

    n = n1 - 1
    yr1 = RBmin_cm/amin_cm
    yr2 = RBmax_cm/amin_cm

    if (yr2 <= yr1) then
        write(*, *) 'for nbi source n ', jsrc, &
                   ' rmax = ', RBmax_cm, ' < rmin = ', RBmin_cm
        write(*, *) 'continue with (rmax - rmin)/a = 1.d-5'
        yr2 = yr1 + 1.d-5
    endif

    yhbd2 = abs(hb)/2.d0

    if (yhbd2 == 0.d0) then
        write(*, *) 'nbi source n ', jsrc, ' vertical size = 0'
        write(*, *) 'continue with vertical size = 1.e-5'
        yhbd2 = 1.d-5
    endif

    ydr = (yr2 - yr1) / n_pencil
    yd  = 2.d0 / abs(yr2 - yr1)
    yr0 = 0.5d0 * (yr1 + yr2)
    do jh=1, n1
        ydz(jh) = 0.d0
    enddo

! perpendicular nbi: calculation of the nbi height

    yhm_n = yhm/amin_cm

    if (yr0 < (aspect_ratio + dx(1))) then
        do j=1, n1
            yhocu = yhm_n + cs3 * sqrt((aspect_ratio + dx(j) + yr0) * (aspect_ratio + dx(j) - yr0))
            if (abs(yhocu) >= x(j) * elon1(j)) EXIT
        enddo
    endif

    ys = 0.d0

    do jpencil=1, n_pencil
        az(jpencil) = (jpencil - 0.5d0) * ydr + yr1
        width = (abs(az(jpencil) - yr0) + 1.d-7) * yd
        ydry(jpencil) = power_dist(width, cr1, cr2)
        ys = ys + ydry(jpencil)
    enddo

    ys = n_pencil / ys
    do jpencil=1, n_pencil
        ydry(jpencil) = ydry(jpencil) * ys
    enddo

    yhocu  = abs(amin_cm*yhocu)
    yhtop  = yhocu + yhbd2
    yhtopa = yhtop/amin_cm
    ynocu  = yhocu/yhbd2
    yhedga = x(n1)*elon1(n1)
    yhedge = amin_cm*yhedga
    jhb05  = 0

    if (yhedge < (yhocu - yhbd2)) then ! nbi out of plasma
        ihb = 0
        return
    endif

    yd05h = amin_cm/yhbd2
    ys = 0.d0

    if (yhtop > yhedge) then ! nbi partly out of plasma
        jedge = (yhtop / yhedge - 1.d0) * (n1 - 1)
        if (jedge < 1) jedge = 1
        ydzout = (yhtop - yhedge) / (jedge*amin_cm)
        zy = yhtopa

        do jhh=1, jedge
            yz1 = zy
            zy  = yz1 - ydzout
            zy12 = 0.5d0 * (zy + yz1)
            z = abs(zy12 * yd05h - ynocu)
            ydhz = 0.d0
            if (z <= 1.d0) then
                ydhz = power_dist(z, ch1, ch2)
                z = abs(zy12 * yd05h + ynocu)
                if (z <= 1.d0) then
                    ydhz = ydhz + power_dist(z, ch1, ch2)
                endif
            endif
            ydzjh = (yz1 - zy) * ydhz
            ys = ys + ydzjh
        enddo
        jhb   = n1
        jhb05 = n1 - 1
    endif

! nbi completely inside plasma

    yz1 = yhedga
    zy  = yz1
    jhb = n1
    do while (jhb > 1)
        yz1 = x(jhb) * elon1(jhb)
        if (yz1 <= yhtopa) EXIT
        zy  = yz1
        jhb = jhb - 1
    enddo

    jhb05 = jhb

    if (jhb > 1) then
        do jhh=2, jhb05+1
            jh1 = jhb05 - jhh + 1
            yz1 = zy
            jh  = jh1 + 1
            zy  = x(jh) * elon1(jh)
            zy12 = 0.5d0 * (zy + yz1)
            z = abs(zy12 * yd05h - ynocu)
            ydhz = 0.d0
            if (z <= 1.d0) then
                ydhz = power_dist(z, ch1, ch2)
                z = abs(zy12 * yd05h + ynocu)
                if (z <= 1.d0) then
                    ydhz = ydhz + power_dist(z, ch1, ch2)
                endif
            endif
            ydz(jh) = (yz1 - zy) * ydhz
            ys = ys + ydz(jh)
        enddo
    endif

    if (ys <= 0.d0) then
        z  = 0.d0
        jh = jhb05
        ydhz = ydhz + power_dist(z, ch1, ch2)
        ydz(jh) = ydhz * 2.d0 * yhbd2
        ys = ys + ydz(jh)
    endif

    yds = yhbd2 / (amin_cm*ys)
    do jh=1, jhb05
        ydz(jh) = ydz(jh) * yds
    enddo

    ihb = 2*jhb05
    ysb = (yr2 - yr1)*amin_cm*2.d0*yhbd2
    coeff = 2.d0*451.9d0*qb/ysb
    y1 = eb/pb

    do jE=jE_min, n_energy
        ye = y1/(n_energy - jE + 1)
        yaqb(jE) = coeff*adqb(jE)/(ye*pb*sqrt(ye))
    enddo

    end subroutine nb0

!---------------------------------------------------------------------
    subroutine compute_dists(slej0, slej2, YVE, YCU, jE_min, nrho_in, nthe_in)

    integer, intent(in) :: jE_min, nrho_in, nthe_in
    double precision, intent(in), dimension(n_energy) :: YVE, YCU
    double precision, intent(out), dimension(n_energy, n_rho) :: slej0, slej2

    integer :: jE, jrho, jthe
    double precision :: Y, YDEDJ

    YDEDJ = 1.d6/(1.6d0*EB)

    do jE=jE_min, n_energy
        do jrho=1, nrho_in
            YAQBA (jE, jrho) = YAQBA (jE, jrho) * YVE(jE)
            YANBA (jE, jrho) = YANBA (jE, jrho) * YVE(jE)
            YANBA1(jE, jrho) = YANBA1(jE, jrho) * YVE(jE)
            YACBA (jE, jrho) = YACBA (jE, jrho) * YCU(jE)
            YATBA (jE, jrho) = YATBA (jE, jrho) * YCU(jE)
        enddo
    enddo

    slej0 = 0.0d0
    slej2 = 0.0d0

    do jE=jE_min, n_energy
        Y = (n_energy - jE + 1) * YDEDJ * YVE(jE)
        do jrho=1, nrho_in
            do jthe=1, nthe_in
                if (YASBA(jE, jrho, jthe) > 0.d0) then
                    YASBA(jE, jrho, jthe) = YASBA(jE, jrho, jthe) * Y
                    slej2(jE, jrho) = slej2(jE, jrho) + YASBA(jE, jrho, jthe) * 2.5d0 * PLEJ2(jthe)
                    slej0(jE, jrho) = slej0(jE, jrho) + 0.5d0 * YASBA(jE, jrho, jthe)
                endif

                if (YASBA1(jE, jrho, jthe) > 0.d0) then
                    YASBA1(jE, jrho, jthe) = YASBA1(jE, jrho, jthe) * Y
                    slej2(jE, jrho) = slej2(jE, jrho) + YASBA1(jE, jrho, jthe) * 2.5d0 * PLEJ2(jthe)
                    slej0(jE, jrho) = slej0(jE, jrho) + 0.5d0 * YASBA1(jE, jrho, jthe)
                endif
            enddo
        enddo
    enddo

    end subroutine compute_dists

!---------------------------------------------------------------------
    subroutine beam_step_pos(r_jn, az_jpencil, rcr, y1, y2, yc2, j_the1, j_the2)

    double precision, intent(in) :: r_jn, az_jpencil
    integer, intent(out) :: j_the1, j_the2  
    double precision, intent(out) :: rcr, y1
    double precision, intent(inout) :: y2, yc2

    integer :: jt1, jt2, ntet
    double precision :: y12, yc1

    y12 = (r_jn - az_jpencil) * (r_jn + az_jpencil)
    y1  = y2
    yc1 = yc2
    y2  = sqrt(y12)
    yc2 = az_jpencil / r_jn
    rcr = y12 / r_jn

    ntet = ntet1 - 1
    jt1 = ntet * yc1
    jt2 = ntet * yc2
    jt1 = min(jt1 + 1, ntet)
    jt2 = min(jt2 + 1, ntet)
    j_the1 = min(jt1, jt2)
    j_the2 = max(jt1, jt2)

    end subroutine beam_step_pos

!---------------------------------------------------------------------
    subroutine beam_step_neg(r_jn, az_jpencil, rcr, y1, y2, yc2, j_the1, j_the2)

    double precision, intent(in) :: r_jn, az_jpencil
    integer, intent(out) :: j_the1, j_the2  
    double precision, intent(out) :: rcr, y1
    double precision, intent(inout) :: y2, yc2

    integer :: jt1, jt2, ntet
    double precision :: y12, yc1

    y12 = (r_jn - az_jpencil) * (r_jn + az_jpencil)
    y1  = y2
    yc1 = yc2
    if (y12 > 0.d0) then
        y2  = -sqrt(y12)
        yc2 = az_jpencil / r_jn
        rcr = y12 / r_jn
    else
        y2  = 0.d0
        yc2 = 1.d0
        rcr = 0.d0
    endif

    ntet = ntet1 - 1
    jt1 = ntet * yc1
    jt2 = ntet * yc2
    jt1 = min(jt1 + 1, ntet)
    jt2 = min(jt2 + 1, ntet)
    j_the1 = min(jt1, jt2)
    j_the2 = max(jt1, jt2)

    end subroutine beam_step_neg

!---------------------------------------------------------------------
    subroutine nbctr(jjn, jE_min, re, ri, az_jpencil, yve, contr, yct2, yaqbp)
! Replacing "counter" part of height_loop

    integer, intent(in) :: jjn, jE_min
    double precision, intent(in) :: contr, az_jpencil
    double precision, intent(in), dimension(n_energy) :: yve
    double precision, intent(in), dimension(n_rho) :: re, ri
    double precision, intent(inout), dimension(n_energy) :: yct2, yaqbp

    integer :: jj, jn, jn1, jnn, j_jn, jhp1, jt1, jt2, jh, n, ntet
    double precision :: y1, y2, yc2, rcr

    ntet = ntet1 - 1
    n = n1 - 1

! loop1 from R = re(A) to R = RC (re(0)) or to R(Y=0)
    do jj=1, jjn
        jn = n1 - jj
        if (y2 == 0.d0) EXIT  ! R = R(Y=0)
        call beam_step_neg(re(jn), az_jpencil, rcr, y1, y2, yc2, jt1, jt2)
        call nbctr3(jn, jt1, jt2, n1, jE_min, rcr, y1, y2, yc2, az_jpencil, contr, yct2, yaqbp)
    enddo

    if (jj == jjn + 1) jn = jh - 1

    if (y2 /= 0.d0) then  ! beam touch R = R(0) (R=RC)
! loop from R = R(0) (R=RC) to R = R(Y=0) or to R = ri(A)
        do jn=jh, N
            jn1 = jn + 1
            if (y2 == 0.d0) EXIT
            call beam_step_neg(ri(jn1), az_jpencil, rcr, y1, y2, yc2, jt1, jt2)
            call nbctr3(jn, jt1, jt2, n1, jE_min, rcr, y1, y2, yc2, az_jpencil, contr, yct2, yaqbp)
        enddo

        if (jn == n1 .and. perp_flag) then ! beam to internal wall
            return
        endif
        jnn  = jn1 - 1
        jhp1 = jh + 1
        do j_jn=jhp1, jnn
            jn = jnn - j_jn + jhp1 - 1
            call beam_step_pos(ri(jn), az_jpencil, rcr, y1, y2, yc2, jt1, jt2)
            call nbctr3(jn, jt1, jt2, n1, jE_min, rcr, y1, y2, yc2, az_jpencil, contr, yct2, yaqbp)
        enddo

        jn = jh - 1
    endif

    jnn = jn + 1
! loop from R = R(0) (R=RC) or from R = R(Y=0) to R = re(A)
    do jn=jnn, N
        jn1 = jn + 1
        call beam_step_pos(re(jn1), az_jpencil, rcr, y1, y2, yc2, jt1, jt2)
        call nbctr3(jn, jt1, jt2, n1, jE_min, rcr, y1, y2, yc2, az_jpencil, contr, yct2, yaqbp)
    enddo

    end subroutine nbctr

!---------------------------------------------------------------------
    subroutine nb1tr(jE_min, n_pencil, yaqbp, contr, ar)
!---- ripple normalized radius of ripple boundary
!---- banana with rtrap > yriplr is lost

    integer, intent(in) :: jE_min, n_pencil
    double precision, intent(in) :: contr
    double precision, intent(in) :: ar(*)
    double precision, intent(inout) :: yaqbp(n_energy)

    integer :: ntet, ji, jn, jt, jpencil, jjn, j_jn, jj, jhb05, jn1, in, in1, &
        jjh, jh, jh1, jjr, jnn, jhp1, jE, jt1, jt2, n
    double precision :: rcr, y1, y2, y, ye, yh, ds, dv, &
        p_shine, yloss, yvede, ye3, yr, ydh, xjh, &
        zjh, yzjh, y12, yr1, yr2, yc2
    double precision, dimension(n_rho) :: re, ri, rc
    double precision, dimension(n_energy) :: drl, ycu, yct2, yve

    p_shine = 0.d0
    ntet = ntet1 - 1
    n    = n1 - 1
    ds = abs(RBmax_cm - RBmin_cm)*amin_cm / (n*n_pencil)
    dv = ds*amin_cm
    yloss = 0.d0

! -------- clean the sources --------
    do jn=1, n1
        rc(jn) = aspect_ratio + dx(jn)
    enddo
    yanba  = 0.d0
    yanba1 = 0.d0
    yacba  = 0.d0
    yaqba  = 0.d0
    yatba  = 0.d0
    YASBA  = 0.d0
    YASBA1 = 0.d0

    if (ihb == 0) return

! nbi out of plasma

    jhb05 = ihb / 2
    ye3   = eb / pb

! for calculation of surface index ii = jn(x(rcrit)) in trapping analysis
! yr = (rj - rii/ab) normalised distance from ext. boundary rj

    rj   = aspect_ratio + 1.d0
    rjt  = aspect_ratio - 1.d0
    ydyh = 1.d2
    yh   = 1.d-2
    jn1 = n1

    do ji=1, ndim2-1
        yr = yh * (ji - 1)
        do while (jn1 > 1)
            jn = jn1 - 1
            if (yr >= (x(n1) + dx(n1) - dx(jn1) - x(jn1)) .and. &
                yr <  (x(n1) + dx(n1) - dx(jn)  - x(jn))) then
                II(ji) = jn
                EXIT
            else
                if (jn1 == 1) jn1 = -1
                do
                    in1 = -jn1
                    in  = in1 + 1
                    if (yr >= (x(n1) + dx(n1) + x(in1) - dx(in1)) .and. &
                        yr <  (x(n1) + dx(n1) + x(in)  - dx(in))) then
                        II(ji) = in1
                        EXIT
                    endif
                    jn1 = jn1 - 1
                enddo
            endif
            jn1 = jn
        enddo
    enddo

    II(ndim2) = n
    do jE=jE_min, n_energy
        drl(jE) = 1.d0 / (0.0144d0*sqrt(eb*pb/(n_energy - jE + 1)) / (bz*amin_cm))
        do jn=1, n1
            ard(jE, jn) = ar(jn)*drl(jE)
        enddo
    enddo

! === end trapping arrays

    do jE=jE_min, n_energy
        YE = YE3/(n_energy - jE + 1)
        YVEDE = DS*sqrt(YE)/451.9d0
        YVE(jE) = YVEDE*PB*YE*1.d-3
        YCU(jE) = -CONTR*VNB(jE)*YVEDE*5.d-6
    enddo

! loop beam`s height H = X(jn)+dX/2
    height_loop: do jjh=1, jhb05
        jh  = jhb05 - jjh + 1
        jh1 = jh + 1
        y = elon1(jh1) * x(jh1) - elon1(jh) * x(jh)
        if (jh == jhb05) then
            ydh = n * ydz(jh)
            xjh = x(jh) + 0.5d0 * ydh / (n * elon1(jh))
            re(jh) = rc(jh) - (rc(jh1) - rc(jh)) * (xjh - x(jh)) / (x(jh1) - x(jh)) - tria1(jh)
        else
            xjh = x(jh) + 0.5d0 * y / elon1(jh)
            ydh = n * ydz(jh)
            re(jh) = 0.5d0 * (rc(jh1) + rc(jh) - tria1(jh) - tria1(jh1))
        endif

        zjh = xjh * elon1(jh)
        ri(jh) = re(jh)

        do jjr=jh, n
            jn = n1 - jjr + jh
            yzjh = zjh / elon1(jn)
            y = sqrt((x(jn) - yzjh) * (x(jn) + yzjh))
            re(jn) = rc(jn) + y - tria1(jn) * (yzjh / x(jn))**2
            ri(jn) = rc(jn) - y - tria1(jn) * (yzjh / x(jn))**2
        enddo

        do jn=n, jh+1, -1
            jj = jn - 1
            if (re(jn) < re(jj)) re(jj) = re(jn) - (re(jj) - re(jn)) / n
            if (ri(jn) > ri(jj)) ri(jj) = ri(jn) - (ri(jj) - ri(jn)) / n
        enddo

        do jn=jh+1, n1
            if (jn > 1 .and. ri(jn-1) < ri(jn)) then
                write(*, *) jn, jh, jhb05, n1
            endif
        enddo

! loop beam radius R = Z(jr)

        do jpencil=1, n_pencil
            if (yaqb(n_energy) <= 0.d0) EXIT ! no power in this pencil
            do jE=jE_min, n_energy
                yaqbp(jE) = yaqb(jE) * ydh * ydry(jpencil)
            enddo

! integral along beam Y

            y12 = (re(n1) - az(jpencil)) * (re(n1) + az(jpencil))
            if (y12 <= 0.d0) EXIT ! beam beyond tokamak
            y2  = -sqrt(y12)
            yc2 = az(jpencil) / re(n1)
            rcr = y12 / re(n1)

            jjn = n1 - jh
            yr2 = contr * sqrt(rj * (rj - rcr))

            if (rcr < rjt) then
                yr1 = contr * sqrt(rjt * (rjt - rcr))
            else
                yr1 = 0.d0
            endif

            do jE=jE_min, n_energy
                yct2(jE) = 0.d0
            enddo

            if (contr < 0.d0) then
                do jj=1, jjn  ! loop from R = re(A) to R = RC (re(0)) or to R(Y=0)
                    jn  = n1 - jj
                    if (y2 == 0.d0) EXIT ! R = R(Y=0)
                    call beam_step_neg(re(jn), az(jpencil), rcr, y1, y2, yc2, jt1, jt2)
                    call nbco3(jn, jt1, jt2, n1, jE_min, rcr, y1, y2, yc2, az(jpencil), contr, yct2, yaqbp)
                enddo

                if (jj == jjn + 1) jn = jh - 1

                if (y2 /= 0.d0) then

                    do jn=jh, n
                        jn1 = jn + 1
                        if (y2 == 0.d0) EXIT
                        call beam_step_neg(ri(jn1), az(jpencil), rcr, y1, y2, yc2, jt1, jt2)
                        call nbco3(jn, jt1, jt2, n1, jE_min, rcr, y1, y2, yc2, az(jpencil), contr, yct2, yaqbp)
                    enddo

                    if (jn == n + 1 .and. perp_flag) then
                        p_shine = p_shine + SUM(yaqbp(jE_min: n_energy) * yve(jE_min: n_energy))
                        CYCLE
                    endif
                    jnn  = jn1 - 1
                    jhp1 = jh + 1

                    do j_jn=jhp1, jnn
                        jn  = jnn - j_jn + jhp1 - 1
                        y1  = y2
                        call beam_step_pos(ri(jn), az(jpencil), rcr, y1, y2, yc2, jt1, jt2)
                        call nbco3(jn, jt1, jt2, n1, jE_min, rcr, y1, y2, yc2, az(jpencil), contr, yct2, yaqbp)
                    enddo

                    jn = jh - 1
                endif

                jnn = jn + 1

                do jn=jnn, n
                    jn1 = jn + 1
                    call beam_step_pos(re(jn1), az(jpencil), rcr, y1, y2, yc2, jt1, jt2)
                    call nbco3(jn, jt1, jt2, n1, jE_min, rcr, y1, y2, yc2, az(jpencil), contr, yct2, yaqbp)
                enddo
            else
                call nbctr(jjn, jE_min, re, ri, az(jpencil), yve, contr, yct2, yaqbp)
            endif

            write(*, *) 'passed'
            p_shine = p_shine + SUM(yaqbp(jE_min: n_energy) * yve(jE_min: n_energy))
        enddo
    enddo height_loop

    call compute_dists(YSLEJ0, YSLEJ2, YVE, YCU, jE_min, n1, ntet1)

    end subroutine nb1tr

!---------------------------------------------------------------------
    subroutine nbtrag(jE_min, n_pencil, YAQBP, CONTR, AR)
!--------- Neutral beam ionization -----------
! coinjection  with ion's trapping and orbital losses
!  AQBA[MW],ACBA[MA m/s],ANBA[10*13 prtcls],
!  ASBA[10#19 prtcl/s]*Dcos(JT)
!  P_Shine [MW] shine through power
!-------------------------------------- Polevoy

    integer, intent(in) :: jE_min, n_pencil
    double precision, intent(in) :: CONTR, AR(*)
    double precision, intent(out) :: yaqbp(n_energy)

    integer :: je, ntet, n, jn, jt, jhb05, jn1, ji, in, in1, jh, jh1, jjh, &
        jjr, jpencil, jjn, j_jn, jj, jt1, jt2, jnn, jhp1
    double precision :: P_SHINE, DS, DV, YLOSS, YE3, YH, YR, YRN, &
        YE, YVEDE, Y, YDH, XJH, ZJH, YZJH, y12, y2, yc2, YRN1, YR2, YR1, &
        rcr, y1, yc1, YDYS, YDEX
    double precision, dimension(n_energy) :: DRL, YCU, YCT2, YVE
    double precision, dimension(n_rho) :: re, ri, rc

    P_SHINE = 0.d0
    ntet = ntet1 - 1
    N    = n1 - 1
    DS = abs(RBmax_cm - RBmin_cm)*amin_cm/(n*n_pencil)
    DV = DS*amin_cm
    YLOSS = 0.d0

! clean the sources
    do jn=1, n1
        rc(jn) = aspect_ratio + DX(jn)
    enddo
    YANBA  = 0.d0
    YANBA1 = 0.d0
    YACBA  = 0.d0
    YATBA  = 0.d0
    YAQBA  = 0.d0
    YASBA  = 0.d0
    YASBA1 = 0.d0

    if (IHB == 0) return ! NBI out of plasma

    JHB05 = IHB/2
    YE3 = EB/PB

! for calculation of surface index II=jn(X(Rcrit)) in trapping analysis
! YR=(RJ-Rii/AB) normalised distance from ext. boundary RJ
! Trapping analysis with gyroradius and surface averaging
! Rlm = Gyrorad(in the midplain)*YDYH/sqrt(rcr)
    RJ   = aspect_ratio + 1.d0
    RJT  = aspect_ratio - 1.d0
    YDYH = 1.d2
    YH   = 1.d-2
    jn1  = n1
    do jE=jE_min, n_energy
        DRL(jE) = 1.d0/(0.0144d0*sqrt(EB*PB/(n_energy - jE + 1))/(BZ*amin_cm))
        do jn=1, n1
            ARD(jE, jn) = AR(jn)*DRL(jE)
        enddo
    enddo

    do ji=1, ndim2-1
        yr = yh * (ji - 1)
        do while (jn1 > 1)
            jn  = jn1 - 1
            yrn1 = aspect_ratio + dx(jn1) + x(jn1)
            yrn  = aspect_ratio + dx(jn)  + x(jn)
            if (yr >= (x(n1)+dx(n1)-dx(jn1)-x(jn1)) .and. &
                yr <  (x(n1)+dx(n1)-dx(jn) -x(jn))) then
                II(ji) = jn
                EXIT
            endif
            jn1 = jn
        enddo
! --- Second branch (only if first didn’t assign) ---
        if (jn1 <= 1) then
            if (jn1 == 1) jn1 = -1
            do
                in1 = -jn1
                in  = in1 + 1
                yrn1 = aspect_ratio + dx(in1) - x(in1)
                yrn  = aspect_ratio + dx(in)  - x(in)
                if (yr >= (x(n1) + dx(n1) + x(in1) - dx(in1)) .and. &
                    yr <  (x(n1) + dx(n1) + x(in)  - dx(in))) then
                    II(ji) = in1
                    EXIT
                endif
                jn1 = jn1 - 1
            enddo
        endif
    enddo

    II(ndim2) = N

    do jE=jE_min, n_energy
        YE = YE3/(n_energy - jE + 1)
        YVEDE = DS*sqrt(YE)/451.9
        YVE(jE) = YVEDE*PB*YE*1.d-3 
        YCU(jE) = -CONTR*VNB(jE)*YVEDE*5.d-6
    enddo

! loop beam`s height H = X(jn)+dX/2
    height_loop: do jjh=1, jhb05
        jh  = jhb05 - jjh + 1
        jh1 = jh + 1
        y = elon1(jh1) * x(jh1) - elon1(jh) * x(jh)
        if (jh == jhb05) then
            ydh = n * ydz(jh)
            xjh = x(jh) + 0.5d0 * ydh / (n * elon1(jh))
            re(jh) = rc(jh) - (rc(jh1) - rc(jh)) * (xjh - x(jh)) / (x(jh1) - x(jh)) - tria1(jh)
        else
            xjh = x(jh) + 0.5d0 * y / elon1(jh)
            ydh = n * ydz(jh)
            re(jh) = 0.5d0 * (rc(jh1) + rc(jh) - tria1(jh) - tria1(jh1))
        endif

        zjh = xjh * elon1(jh)
        ri(jh) = re(jh)

        do jjr=jh, n
            jn = n1 - jjr + jh
            yzjh = zjh / elon1(jn)
            y = sqrt((x(jn) - yzjh) * (x(jn) + yzjh))
            re(jn) = rc(jn) + y - tria1(jn) * (yzjh / x(jn))**2
            ri(jn) = rc(jn) - y - tria1(jn) * (yzjh / x(jn))**2
        enddo

        do jn=n, jh+1, -1
            jj = jn - 1
            if (re(jn) < re(jj)) re(jj) = re(jn) - (re(jj) - re(jn)) / n
            if (ri(jn) > ri(jj)) ri(jj) = ri(jn) - (ri(jj) - ri(jn)) / n
        enddo

        do jn=jh+1, n1
            if (jn > 1 .and. ri(jn-1) < ri(jn)) then
                write(*, *) jn, jh, jhb05, n1
            endif
        enddo

! loop beam radius R = Z(jpencil)
        do jpencil=1, n_pencil
            if (yaqb(n_energy) <= 0.d0) EXIT ! no power in this pencil
            do je=jE_min, n_energy
                yaqbp(jE) = yaqb(jE) * ydh * ydry(jpencil)
            enddo

! integral along beam Y
            y12 = (re(n1) - az(jpencil)) * (re(n1) + az(jpencil))
            if (y12 <= 0.d0) EXIT ! beam beyond tokamak
            y2  = -sqrt(y12)
            yc2 = az(jpencil) / re(n1)
            jjn = n1 - jh
            rcr = y12 / re(n1)
            yr2 = contr * sqrt(rj * (rj - rcr))

            if (rcr < rjt) then
                yr1 = contr * sqrt(rjt * (rjt - rcr))
            else
                yr1 = 0.d0
            endif

            do je=jE_min, n_energy
                yct2(jE) = 0.d0
            enddo

            if (contr < 0.d0) then
                do jj=1, jjn   ! loop from R = re(A) to R = RC (re(0)) or to R(Y=0)
                    jn  = n1 - jj
                    if (y2 == 0.d0) EXIT ! R = R(Y=0)
                    call beam_step_neg(re(jn), az(jpencil), rcr, y1, y2, yc2, jt1, jt2)
                    call nbco2gc(jn, n1, jE_min, rcr, y1, y2, yc2, az(jpencil), contr, yct2, yaqbp, ydys)
                enddo

                if (jj == jjn + 1) jn = jh - 1

                if (y2 /= 0.d0) then ! beam touch R = R(0) (R=RC)

                    do jn=jh, N ! loop from R = R(0) (R=RC) to R = R(Y=0) or to R = ri(A)
                        jn1 = jn + 1
                        if (y2 == 0.d0) EXIT
                        call beam_step_neg(ri(jn1), az(jpencil), rcr, y1, y2, yc2, jt1, jt2)
                        call nbco2gc(jn, n1, jE_min, rcr, y1, y2, yc2, az(jpencil), contr, yct2, yaqbp, ydys)
                    enddo

                    if (jn == n1) then ! beam to internal wall
                        p_shine = p_shine + SUM(yaqbp(jE_min: n_energy) * yve(jE_min: n_energy))
                        CYCLE
                    endif

! R increases
                    jnn  = jn1 - 1
                    jhp1 = jh + 1

                    do j_jn=jhp1, jnn   ! loop from R = R(Y=0) to R =R(0) (R=RC)
                        jn  = jnn - j_jn + jhp1 - 1
                        call beam_step_pos(ri(jn), az(jpencil), rcr, y1, y2, yc2, jt1, jt2)
                        call nbco2gc(jn, n1, jE_min, rcr, y1, y2, yc2, az(jpencil), contr, yct2, yaqbp, ydys)
                    enddo

                    jn = jh - 1
                endif

                jnn = jn + 1

                do jn=jnn, n  ! loop from R = R(0) (R=RC) or from R = R(Y=0) to R = re(A)
                    jn1 = jn + 1
                    call beam_step_pos(re(jn1), az(jpencil), rcr, y1, y2, yc2, jt1, jt2)
                    call nbco2gc(jn, n1, jE_min, rcr, y1, y2, yc2, az(jpencil), contr, yct2, yaqbp, ydys)
                enddo

            else
                call nbctr(jjn, jE_min, re, ri, az(jpencil), yve, contr, yct2, yaqbp)
            endif

            p_shine = p_shine + SUM(yaqbp(jE_min: n_energy) * yve(jE_min: n_energy)) ! Shinethrough power distribution R[m], Z[m], Q/S[W/cm2]
        enddo  ! loop on beam`s radius R = Z(jpencil)
    enddo height_loop

    call compute_dists(YSLEJ0, YSLEJ2, YVE, YCU, jE_min, n1, ntet1)

    end subroutine nbtrag

!---------------------------------------------------------------------
    subroutine nbsrsr(j_nbi, YCONTR, NA1, RTOR, SHIFT, AB, BTOR, &
        HRO, YHM, CBMH1, CBMH2, CBMS1, CBMS2, CBMS3, CBMS4, &
        CBMR1, CBMR2, dn_rho, fp_flag, EBEAM, power_frac, &
        ABEAM, QBEAM, RBMAX, RBMIN, JSRREC, YEXTARR, calc_fus, jE_min)

!---------------------------------------------------------------------
! fast ion's sourses (for multi sources) + ripple
! Suzuki stopping crossections
! New SCUBM =Nt/m3					     22-APR-13
!--------------------------------------------------------------Polevoy
!	entry:	AMETR,SHIF,NA1,RTOR,AB,BTOR,NI,HBEAM,RBMIN,RBMAX,
!     		CBMS2,dn_rho,EBEAM,QBEAM,ABEAM,power_frac,
!	     	TE,TI,NE,NN,ZEF,AMAIN,PBEAM,SCUBM,CONTR,ELON,
!	CBMS2	number of 'pencils'
!	dn_rho	switch for number of internal mesh points
!			 41 (dn_rho=1),21 (dn_rho=2)
!	power_frac3,2,1 power fraction of energy comps.
!		    	 3(EB,EB/2,EB/3),2(EB,EB/2),1(EB)
!	CONTR	Qcontr/Qbeam
!	j_nbi 	Number of the current hot ion source
!	JSRREC 	Length of the hot ion source record
!	exit:	PBEAM,SCUBM,SNEBM,SNNBM		for MAIN
!	SCUBM	Toroidal pulse [kg*m/s2/m3]	05-AUG-96
!---------------------------------------------------------------------

    use nbstatus, only: ISPE, ISPEND, &
        AMETR, VR, RHO, FP, AMAIN, ELON, TRIA, SHIF, NE, TE, TI, &
        SCUBM, PBEAM, NNBM1, NNBM2, NNBM3, SNNBM, SNEBM, &
        SNIBM1, SNIBM2, SNIBM3
    use cross_sections, only: svdtbp, svddnp1, svddnp2, seiv

    integer, parameter :: JDBL=2
    double precision, parameter :: GP2=6.283185d0

    integer, intent(in) :: j_nbi, NA1, dn_rho, fp_flag, calc_fus
    double precision, intent(in) :: YCONTR, RTOR, SHIFT, AB, BTOR, &
        HRO, YHM, CBMH1, CBMH2, CBMS1, CBMS2, CBMS3, CBMS4, &
        CBMR1, CBMR2, EBEAM, power_frac(n_energy), &
        ABEAM, QBEAM, RBMAX, RBMIN, yEXTARR(n_rho, nspec_max)
    integer, intent(out) :: JSRREC, jE_min

    double precision :: AR(n_rho), Rmaj_cm, YE2, YCX
    double precision, dimension(n_energy) :: AQBP, ADQB
    double precision :: YSCU1, YFVDA, YNHDT, &
        YDBM, Yjn, YPOW, YD, YDV, YDDV, YjE, YE1, YEPS
    integer :: N, J, ntet, jn, jn1, jnA, jnAX, jE, JT, jnA1, jnAC, JS, JSP, n_pencil

    n1 = (NA1 - 1)/dn_rho + 1
    n = n1 - 1
    if (RBMAX <= RBMIN) then
        write(*, *) 'ILLEGAL: NBI source N', j_nbi, ' RBMAX <= RBMIN !!!'
        return
    endif
    if (EBEAM*ABEAM == 0.d0) then
        write(*, *) 'ILLEGAL: NBI source N', j_nbi, ' EBEAM*ABEAM = 0 !!!'
        return
    endif

    n_pencil = CBMS2
    YDBM = SUM(power_frac)

    do j=1, n_energy
        ADQB(j) = 0.d0
    enddo
    if (power_frac(1) <= 0.d0) then
        write(*, *) 'ILLEGAL: NBI source N', j_nbi, ' power_frac(1) <= 0 !!!'
        return
    endif
    ADQB(n_energy) = power_frac(1)/YDBM
    jE_min = n_energy
    if (power_frac(2) > 0.d0) then
        ADQB(2) = power_frac(2)/YDBM
        jE_min = n_energy - 1
    endif
    if (power_frac(n_energy) > 0.d0) then
        ADQB(1) = power_frac(n_energy)/YDBM
        jE_min = n_energy - 2
    endif
    if (fp_flag /= 1) then
        ntet1 = 50/2 + 1
    else
        ntet1 = 50 + 1
    endif

    ntet = ntet1 - 1

    if (fp_flag /= 1) then
! N- number of surfaces,  2*ntet (cntr+co angle)
        JSRREC = JDBL*4*(1 + 3*N*2*ntet)
    endif

    X(1)    = 0.d0
    XJ(1)   = 0.d0
    XJ(NA1) = 1.d0
    AR(1)   = 0.d0
    EB = EBEAM*1.d3
    QB = QBEAM*1.d3
    PB = ABEAM
    Rmaj_cm = (RTOR + SHIFT)*1.d2
    amin_cm = AMETR(NA1)*1.d2
    aspect_ratio = Rmaj_cm/amin_cm
    perp_flag = (((RBMAX + RBMIN)/2.d0) > (RTOR - AB)) ! Perpendicular NBI

    BZ = BTOR/(1. + SHIFT/RTOR)
    HB = CBMS4*(RBMAX - RBMIN)*1.d2*sqrt(1.d0 + CBMS3*CBMS3)
    if (HB <= 0.d0) then
        write(*, *) 'NBI input for the source ', j_nbi, ' is not correct'
        write(*, *) 'please use: RBMAX > RBMIN and CBMS4 > 0'
    endif
    RBmin_cm = RBMIN*1.d2
    RBmax_cm = RBMAX*1.d2
    YD = 1.d4/(BZ*GP2*amin_cm**2)
    do jn=2, NA1
        jn1 = jn - 1
        AR(jn) = (FP(jn) - FP(1))*YD
        XJ(jn) = 5.d1*(AMETR(jn) + AMETR(jn1))/amin_cm
    enddo
    Yjn = 0.d0
    do jn=1, n1
        jnAX = 1 + dn_rho*(jn - 1)
        jnA  = jnAX - dn_rho*Yjn
        if (dn_rho > 1) Yjn  = .5d0
        ELON1(jn) = max(ELON(jnAX), 1.d0)
        TRIA1(jn) = TRIA(jnAX)*XJ(jnAX)
        AR(jn)    = AR(jnA)
        X(jn)     = XJ(jnAX)
        DX(jn)    = 1.d2*(SHIF(jnAX) - SHIFT)/amin_cm
        YEPS      = AMETR(jnAX)/(RTOR + SHIF(jnAX))
        YTRAP(jn) = sqrt(2.*YEPS/(1. + YEPS))
        do jE=1, n_energy
            YAQBA (jE, jn) = 0.d0
            YACBA (jE, jn) = 0.d0
            YATBA (jE, jn) = 0.d0
            YANBA (jE, jn) = 0.d0
            YANBA1(jE, jn) = 0.d0
            YSLEJ0(jE, jn) = 0.d0
            YSLEJ2(jE, jn) = 0.d0
            do JT=1, ntet1
                YASBA (jE, jn, JT) = 0.d0
                YASBA1(jE, jn, JT) = 0.d0
            enddo
        enddo
    enddo

    X(n1)  = 1.d0
    DX(n1) = 0.d0

    if (QB <= 0.d0) return

    call NBSISN(NA1, EBEAM, ABEAM, dn_rho, YEXTARR)
    call NB0(jE_min, n_pencil, j_nbi, ADQB, YHM, CBMH1, CBMH2, CBMR1, CBMR2, CBMS3)

    if (CBMS1 < 1.d0) then
        call NB1TR(jE_min, n_pencil, AQBP, YCONTR, AR)  ! Updating also YSLEJ0, YSLEJ2
    else
        call NBTRAG(jE_min, n_pencil, AQBP, YCONTR, AR) ! Updating also YSLEJ0, YSLEJ2
    endif
    YSCU1 = ABEAM*0.0209d0*0.5d-2*amin_cm
    YPOW = 0.d0
    jnA1 = 1
    do jn=2, n1
        jn1 = jn - 1
        if (jn > 2) jnA1 = 2 + dn_rho*(jn - 2)
        jnA = jnA1 + 1
        jnAC = 1 + dn_rho*(jn - 1)
        if (jn == n1)then
            jnAC = NA1 - 1
            YDV = VR(jnAC)*(RHO(NA1) - jnAC*HRO)
        else
            YDV = 0.d0
        endif
        if (jna1 > jnAC) jnA1 = jnAC
        jnAX = jnAC

        do J=jnA1, jnAC
            YDV = VR(j)*HRO + YDV
        enddo
        YDDV = 1./YDV
        do jE=jE_min, n_energy
            YSLEJ0(jE, jn1) = YSLEJ0(jE, jn)*YDDV
            YSLEJ2(jE, jn1) = YSLEJ2(jE, jn)*YDDV
            YAQBA (jE, jn1) = YAQBA (jE, jn)*YDDV
            YACBA (jE, jn1) = YACBA (jE, jn)*YDDV
            YATBA (jE, jn1) = YATBA (jE, jn)*YDDV
            YjE = EBEAM/(n_energy - jE + 1)
            if (fp_flag > 0) then ! convert prtcl/s => prtcl/s/m^3
                do JT=1, ntet1
                    YASBA (jE, jn1, JT) = YASBA (jE, jn, JT)*YDDV
                    YASBA1(jE, jn1, JT) = YASBA1(jE, jn, JT)*YDDV
                enddo
            endif
            do j=jnA1, jnAC
                PBEAM(j) = PBEAM(j) + YAQBA(jE, jn1)
                if (calc_fus > 0) then ! beam-plasma fusion
                    if (ABEAM == 3.d0) stnbdp(j) = stnbdp(j) + YAQBA(jE, jn1) * &
                        svdtbp(YjE, ABEAM, NE(j), TE(j), TI(j), AMAIN(j), calc_fus)*(625.d0/YjE)
                    if (ABEAM == 2.d0) then
                        sdnbtp(j)  = sdnbtp(j)  + YAQBA(jE, jn1)*(625.d0/YjE) * svdtbp (YjE, ABEAM, NE(j), TE(j), TI(j), AMAIN(j), calc_fus)
                        sdnbdp2(j) = sdnbdp2(j) + YAQBA(jE, jn1)*(625.d0/YjE) * svddnp2(YjE, ABEAM, NE(j), TE(j), TI(j), AMAIN(j), calc_fus)
                        sdnbdp1(j) = sdnbdp1(j) + YAQBA(jE, jn1)*(625.d0/YjE) * svddnp1(YjE, ABEAM, NE(j), TE(j), TI(j), AMAIN(j), calc_fus)
                    endif
                endif
                SCUBM(j) = SCUBM(j) + YATBA(jE, jn1)*YSCU1
            enddo

            YANBA (jE, jn1) = YANBA (jE, jn)*YDDV*625.d0/YjE
            YANBA1(jE, jn1) = YANBA1(jE, jn)*YDDV*625.d0/YjE
        enddo

        do J=jnA1, jnAC
            YFVDA = F(1, jn1)*VNB(1)/amin_cm
            if (YFVDA /= 0.d0) NNBM1(j) = NNBM1(j) + YANBA(3, jn1)/YFVDA
            YFVDA = F(2, jn1)*VNB(2)/amin_cm
            if (YFVDA /= 0.d0) NNBM2(j) = NNBM3(j) + YANBA(2, jn1)/YFVDA
            YFVDA = F(3, jn1)*VNB(3)/amin_cm
            if (YFVDA /= 0.d0) NNBM3(j) = NNBM3(j) + YANBA(1, jn1)/YFVDA
            YNHDT = 0.d0   !  total proton content
            do JS=2, ISPEND
                JSP = ISPE(JS)
                if (ZB(JS) == 1.) YNHDT = YNHDT + yEXTARR(jnA, JSP)
            enddo
! Calc. total particle sourse: SNEBM
            YE1 = NE(jnA)*SEIV(TE(jnA))
            do jE=jE_min, n_energy
                YE2 = 0.d0
                do JS=2, ISPEND
                    JSP = ISPE(JS)
                    YE2 = YE2 + yEXTARR(jnA, JSP)*YSIMPI(jE, JS)*VNB(jE)
                enddo
! Calc. thermal neutral  sourCe: SNNBM
                YCX = YNHDT*SVEX(4, jE)
                SNNBM(j) = SNNBM(j) + YCX/(YCX + YE1 + YE2)*YANBA1(jE, jn1)
                SNEBM(j) = SNEBM(j) + (1.d0 - YCX/(YCX + YE1 + YE2))*YANBA1(jE, jn1)
                if (ABEAM < 1.5d0) SNIBM1(j) = SNIBM1(j) + (1.d0 - YCX/(YCX + YE1 + YE2))*YANBA(jE, jn1) !H ion source
                if (ABEAM < 2.5d0 .and. ABEAM > 1.5d0) SNIBM2(j) = SNIBM2(j) + (1.d0 - YCX/(YCX + YE1 + YE2))*YANBA(jE, jn1) !D ion source
                if (ABEAM > 2.5d0) SNIBM3(j) = SNIBM3(j) + (1.d0 - YCX/(YCX + YE1 + YE2))*YANBA(jE, jn1) !T ion source
            enddo
        enddo
    enddo

    end subroutine nbsrsr

!---------------------------------------------------------------------
    subroutine nbion0(jE_min, NA1, ABEAM, EBEAM, RTOR, dn_rho, yEXTARR)
!---------------------------------------------------------------------
! Steady State (1+2D:(x+MU, V)) Fokker-Plank Solver
!	PEBM, PIBM(X)-power to electrons, ions [MW/m3]
!	CUFI, CUBM(X)- fast ion's and NBI driven currents [MA/m2]
!	NIBM(x) - fast ion's density [10^19/m3]
!-------------------------------------------------------------- Polevoi

    use nbstatus, only: ISPE, ISPEND, AMAIN, NE, TE, TI, ZEF, &
        AMETR, SHIF, PBPER, PBLON, NIBM, PEBM, PIBM, CUFI, CUBM
    use cross_sections, only: fnbf, fnb2, fnbp, fnbi1

    integer, intent(in) :: jE_min, NA1, dn_rho
    double precision, intent(in) :: ABEAM, EBEAM, RTOR, yEXTARR(n_rho, nspec_max)

    integer :: jn, JS, JSP, jn1, jnA, jnA1, jnAC, J, J2, jE
    double precision :: yc1, YS, YSTE, YEPS, YDN, Y, X1, X2, X3, &
        YA, YB, YC, YC0, y1, y12, YA1, YD1, YI0, YI2, YP11, YCRNT, YPIDPB, &
        STSD3
    double precision, dimension(n_rho) :: YZ2D3, YTSE, YFCUR, YLNI, YLNE, YLNZ, YEBDEC

    do jn=1, NA1
        YSTE = sqrt(TE(jn))
        if (EBEAM > 1.d2*ABEAM) then
            YLNI(jn) = 23.7d0 + LOG(AMAIN(jn)/(AMAIN(jn) + ABEAM)*sqrt(1.d-3*ABEAM*EBEAM*TE(jn)/NE(jn)))
        else
            YLNI(jn) = 25.4d0 + LOG(1.d-3*EBEAM*AMAIN(jn)/(AMAIN(jn) + ABEAM)*sqrt(TE(jn)/NE(jn)))
        endif
        YLNE(jn) = 15.85d0 + log(TE(jn)/sqrt(NE(jn)))
        YLNZ(jn) = YLNI(jn)
        YS = 0.d0
        y1 = 0.d0

        do JS=2, ISPEND
            JSP = ISPE(JS)
            if (EBEAM > 1.d2*ABEAM) then
                YLNZ(jn) = 23.7d0 + log(RMB(JS)/(RMB(JS) + ABEAM)*sqrt(1.d-3*ABEAM*EBEAM*TE(jn)/NE(jn)))
            else
                YLNZ(jn) = 25.4d0 + log(1.d-3*EBEAM*RMB(JS)/(RMB(JS) + ABEAM)*sqrt(TE(jn)/NE(jn)))
            endif
            Y = yEXTARR(jn, JSP)*ZB(JS)**2*YLNZ(jn)
            y1 = y1 + Y
            YS = YS + Y/RMB(JS)
        enddo

        YEBDEC(jn) = EBEAM/(TE(jn)*ABEAM*14.8d0*(YS/(NE(jn)*YLNE(jn)))**0.6667d0)
        YZ2D3(jn)  = y1/(YS*3.0d0)
        YTSE(jn)   = 2.d0*ABEAM*YSTE*TE(jn)/(NE(jn)*YLNE(jn))
        YEPS = AMETR(jn)/(RTOR + SHIF(jn))
        YFCUR(jn) = (1.d0 - FNBF(ZEF(jn), YEPS)/ZEF(jn))
        YDN = 0.d0
    enddo

    do jn=2, n1
        jn1  = jn - 1
        jnA1 = 1 + dn_rho*(jn1 - 1)
        jnA  = jnA1 + 1
        jnAC = jnA1 - 1 + dn_rho
        if (jn == n1) jnAC = NA1 - 1
        J2 = 1 + dn_rho*(jn1-1) - dn_rho*YDN
        if (dn_rho > 1) YDN = 0.5d0
        do jE=jE_min, n_energy
            Y = 1.d0/(n_energy - jE + 1)
            do J=jnA1, jnAC
                X2 = YEBDEC(J2)*Y
                X1 = sqrt(X2)
                X3 = X2*X1
                YPIDPB = 2.0d0*FNB2(X2)/X2
                STSD3 = 208.33d0*YTSE(J2)*YAQBA(jE, jn1)
! Tail correction for current
                YB  = 1.d0 + 1.d0/X3
                YA  = 0.5d0*TE(J2)/(EBEAM*Y)*(1.d0 + TI(J2)/(TE(J2)*X3))
                YC0 = YTSE(J2)*DTCX(jE, J2) - 3.d0
                YC  = YZ2D3(J2)*3.d0/X3
                yc1 = YC0 + YC
                y1  = 4.d0*YA*yc1/YB**2
                y12 = 1.d0 + sqrt(1.d0 + y1)
                YA1 = 0.5d0*YB*y12/YA
                YD1 = 2.d0/(y12 + 0.5d0*y1)
                YI0 = (0.5d0 - FNB2(X2)/X2)
                YI2 = FNBP(X1, YZ2D3(J2))
! Beam energy Wbeam[10#19keV/m3] = Pbper + Pblon/2
                YP11 = 2.d0/3.d0*(YI0*YSLEJ0(jE, jn1) + 0.4d0*YI2*YSLEJ2(jE, jn1))
! Beam perpendicular pressure <Mb Vper2>/2 [10*19 keV/m3]
                PBPER(j) = PBPER(j) + EBEAM*Y*YTSE(J2)*(2.d0*YI0*YSLEJ0(jE, jn1) - YP11)
! Beam parallel pressure    <Mb Vpar2> [10*19 keV/m3]
                PBLON(j) = PBLON(j) + EBEAM*Y*YTSE(J2)*2.d0*YP11
! Beam density
                NIBM(j) = NIBM(j) + STSD3*LOG(1.d0 + X3)/(EBEAM*Y)
! Fast ion current without trapping correction for Eb
                YCRNT = YACBA(jE, jn1)*YTSE(J2)*YD1*FNBI1(X2, YZ2D3(J2))
                CUFI(j) = CUFI(j) + YCRNT
                PIBM(j) = PIBM(j) + YAQBA(jE, jn1)*YPIDPB
                PEBM(j) = PEBM(j) + YAQBA(jE, jn1)*(1.0d0 - YPIDPB)
                CUBM(j) = CUBM(j) + YCRNT*YFCUR(J2)
            enddo
        enddo
    enddo

    end subroutine nbion0

!---------------------------------------------------------------------
    subroutine nbsisn(NA1, YEBEAM, YABEAM, dn_rho, YEXTARR)
!---------------------------------------------------------------------
!			25-MAY-11 StotQ->StotQ1 Janev->Suzuki
!	Fij = A/L(ri, Ej) - local normalized inverse Neutral Beam
!	mean free path for Ej energy beam component,  Ej beam power
!	absorption along the beam line :
!	Ij(x)	=Ij0*exp(-S Fj dx), 	dx=dL/A
!------------------------------------------------------------- Polevoy

    use nbstatus, only: ISPE, ISPEND, NE, TE, AMETR
    use cross_sections, only: STOTQ1, seiv, spii, spex, simpi

    integer, intent(in) :: NA1, dn_rho
    double precision, intent(in) :: YEBEAM, YABEAM, YEXTARR(n_rho, nspec_max)

    integer :: J, jn, JSP, jE, JS
    double precision :: YY, y1, y2, y3, y12, y13, y23, YjE

    y3 = yEBEAM/yABEAM
    y1 = y3/3.d0
    y2 = y3/2.d0

    do jE=1, n_energy
        YjE = y3/(n_energy - jE + 1)
        VNB(jE) = 1.383d6*sqrt(1.d3*YjE)
        SVII(jE) = SPII(YjE)
        SVEX(4, jE) = SPEX(YjE)
        do JS=2, ISPEND
            if (ZB(JS) > 1.d0) then
                YSIMPI(jE, JS) = SIMPI(YjE, ZB(JS))
            else
                YSIMPI(jE, JS) = SVEX(4, jE) + SVII(jE)
            endif
        enddo
        SVEX(4, jE) = SPEX(YjE)*VNB(jE)
    enddo

    SVEX(4, 4) = SPEX(1.d0)
    y12 = y1 + y2 - 2.d0 * sqrt(y1*y2)
    y13 = y1 + y3 - 2.d0 * sqrt(y1*y3)
    y23 = y2 + y3 - 2.d0 * sqrt(y2*y3)
    SVEX(1, 2) = SPEX(y12)
    SVEX(1, 3) = SPEX(y13)
    SVEX(2, 3) = SPEX(y23)

    do jE=1, n_energy
        YY = y3/(n_energy - jE + 1)
        do jn =1, n1
            J = (jn - 1) * dn_rho + 1
! Janev, Boley (for Eb>0.1 MeV)
            if (YY >= 1.d2) then
                F(jE, jn) = 0.d0
                do JS=2, ISPEND
                    JSP = ISPE(JS)
                    F(jE, jn) = F(jE, jn) + yEXTARR(J, JSP) * &
                         STOTQ1(yABEAM, YY, NE(j), TE(j), ZB(JS), RMB(JS)) !by Suzuki
                enddo
            else
                F(jE, jn) = NE(j) * SEIV(TE(j)) / VNB(jE)
                do JS=2, ISPEND
                    JSP = ISPE(JS)
                    F(jE, jn) = F(jE, jn) + yEXTARR(J, JSP)*YSIMPI(jE, JS)
                enddo
            endif
            F(jE, jn) = F(jE, jn)*(AMETR(NA1-1) + AMETR(NA1))*5.d1
! Charge-exchange transparancy
        enddo
    enddo

    end subroutine nbsisn

!---------------------------------------------------------------------
    subroutine nbco3(jn, jt1, jt2, n1, jE_min, rcr, y1, y2, yc2, az_jpencil, contr, yct2, yaqbp)

    use nbstatus, only: yriplr

    integer, intent(in) :: jn, jt1, jt2, n1, jE_min
    double precision, intent(in) :: rcr, y1, y2, yc2, az_jpencil, contr
    double precision, dimension(n_energy), intent(inout) :: yct2, yaqbp

    integer :: jii, je, ji, jn1, jtrap, j_the
    double precision :: dthe, yr1, yr2, dy, yd, ydcos, ytcos, ydys, yqbp, yf0, yfa, yddd, dy_dthe

    jn1 = jn + 1
    dthe = 1./(jt2 - jt1 + 1)
    YR2 = -sqrt(RJ*(RJ - rcr))

    if (rcr < RJT) then
        YR1 = -sqrt(RJT*(RJT - rcr))
        JI = n1
    else
        YR1 = 0.d0
        JII = (RJ - rcr)*YDYH
        JI  = II(JII)
    endif

    dy = y2 - y1

    do jE=jE_min, n_energy
        YDYS  = F(jE, jn) * dy
        YQBP  = YAQBP(jE) * exp(-YDYS)
        YF0   = ARD(jE, jn) - az_jpencil
        YFA   = YF0 - ARD(jE, n1)
        YDDD  = YAQBP(jE) - YQBP
        YD    = 0.d0
        YDCOS = 0.d0
        YTCOS = yc2 + YCT2(jE)
        YCT2(jE) = yc2

        if (ARD(jE, JI) < YF0 .and. ARD(jE, n1) > YF0) then
            JTRAP = 0
            if ( (yfa > 0.d0 .or. YFA < YR2) .and. rcr >= yriplr(ji)) then ! kept
                YD = YDDD
            endif
        else     ! passing
            JTRAP = 1
            if (YFA > YR1 .or. YFA < YR2) then ! kept
                YD = YDDD
                YDCOS = YD * (yc2 + yct2(jE))
            endif
        endif

        YANBA(jE, jn1) = YANBA(jE, jn1) + YDDD
        YATBA(jE, jn1) = YATBA(jE, jn1) - 2.d0 * CONTR * YDDD * az_jpencil
        YAQBA(jE, jn1) = YAQBA(jE, jn1) + YD

        if (YTRAP(jn1) <= YTCOS) YACBA(jE, jn1) = YACBA(jE, jn1) + YDCOS

        YAQBP(jE) = YQBP
        dy_dthe = YD * dthe

        if (JTRAP /= 0) then   ! passing
            do j_the=jt1, jt2
                YASBA(jE, jn1, j_the) = YASBA(jE, jn1, j_the) + dy_dthe
            enddo
        else   ! trapped
            do j_the=jt1, jt2
                YASBA (jE, jn1, j_the) = YASBA (jE, jn1, j_the) + 0.5d0 * dy_dthe
                YASBA1(jE, jn1, j_the) = YASBA1(jE, jn1, j_the) + 0.5d0 * dy_dthe
            enddo
        endif

    enddo

    end subroutine nbco3

!---------------------------------------------------------------------
    subroutine nbctr3(jn, jt1, jt2, n1, jE_min, rcr, y1, y2, yc2, az_jpencil, contr, yct2, yaqbp)

    use nbstatus, only: yriplr

    integer, intent(in) :: jn, jt1, jt2, n1, jE_min
    double precision, intent(in) :: rcr, y1, y2, yc2, az_jpencil, contr
    double precision, dimension(n_energy), intent(inout) :: yct2, yaqbp

    integer :: jii, je, ji, jn1, jtrap, j_the
    double precision :: dthe, yr1, yr2, dy, yd, ydcos, ytcos, ydys, yqbp, yf0, yfa, yddd, dy_dthe

    jn1 = jn + 1
    dthe = 1./(jt2 - jt1 + 1)
    YR2 = sqrt(RJ*(RJ - rcr))

    if (rcr < RJT) then
        YR1 = sqrt(RJT*(RJT - rcr))
        JI = n1
    else
        YR1 = 0.d0
        JII = (RJ - rcr)*YDYH
        JI  = II(JII)
    endif

    dy = y2 - y1

    do jE=jE_min, n_energy
        YDYS  = F(jE, jn) * dy
        YQBP  = YAQBP(jE) * exp(-YDYS)
        YF0   = ARD(jE, jn) + az_jpencil
        YFA   = YF0 - ARD(jE, n1)
        YDDD  = YAQBP(jE) - YQBP
        YD    = 0.d0
        YDCOS = 0.d0
        YTCOS = yc2 + yct2(jE)
        yct2(jE) = yc2

        if (ARD(jE, JI) < YF0) then  ! trapped
            JTRAP = 0
            if ( (YFA > YR2 .or. YFA < -YR2) .and. rcr < yriplr(ji) ) then ! kept
                YD = YDDD
            endif
        else     ! passing
            if (YFA < YR1 .or. YFA > YR2) then ! kept
                YD = YDDD
                YDCOS = YD * (yc2 + yct2(jE))
            endif
        endif
 
        YANBA(jE, jn1) = YANBA(jE, jn1) + YDDD
        YATBA(jE, jn1) = YATBA(jE, jn1) + 2.d0 * CONTR * YDDD * az_jpencil
        YAQBA(jE, jn1) = YAQBA(jE, jn1) + YD

        if (YTRAP(jn1) <= YTCOS) YACBA(jE, jn1) = YACBA(jE, jn1) + YDCOS

        YAQBP(jE) = YQBP
        dy_dthe = YD * dthe

        if (JTRAP /= 0) then   ! passing (contr)
            do j_the=jt1, jt2
                YASBA1(jE, jn1, j_the) = YASBA1(jE, jn1, j_the) + dy_dthe
            enddo
        else   ! trapped
            do j_the=jt1, jt2
                YASBA (jE, jn1, j_the) = YASBA (jE, jn1, j_the) + 0.5d0 * dy_dthe
                YASBA1(jE, jn1, j_the) = YASBA1(jE, jn1, j_the) + 0.5d0 * dy_dthe
            enddo
        endif

    enddo

    end subroutine nbctr3

!---------------------------------------------------------------------
    subroutine nbco2gc(jn, n1, jE_min, rcr, y1, y2, yc2, az_jpencil, contr, yct2, yaqbp, ydys)

    use nbstatus, only: yriplr

    integer, intent(in) :: jn, n1, jE_min
    double precision, intent(in) :: rcr, y1, y2, yc2, az_jpencil, contr
    double precision, intent(inout) :: ydys
    double precision, intent(inout), dimension(n_energy) :: yct2, yaqbp

    integer :: jii, je, ji, jn1, jtrap, j_the, ntet, jloss, j_jn, jnR, jnL
    double precision :: yr2, yr1, dy, yd, yqbp, &
         yf0, yfa, yddd, YJSERF, ydex
    double precision, dimension(n_rho) :: ycos

    jn1 = jn + 1
    ntet = ntet1 - 1
    YR2 = -sqrt(RJ*(RJ - rcr))
    if (rcr < RJT) then
        YR1 = -sqrt(RJT*(RJT - rcr))
        JI  = n1
    else
        YR1 = 0.d0
        JII = (RJ - rcr)*YDYH
        JI  = II(JII)
    endif

    dy = y2 - y1
    if (dy <= 0.d0) then
        WRITE(*, *) 'NBTRAG: warning: overlapping of surfaces'
        dy = -dy
    endif

    YDEX = exp(-YDYS) ! Using old YDYS, corresponds to original Polevoi's passing assumption
    do jE=jE_min, n_energy
        YDYS  = F(jE, jn)*dy
        YQBP  = YAQBP(jE)*YDEX
        YF0   = ARD(jE, jn) - az_jpencil
        YFA   = YF0 - ARD(jE, n1)
        YDDD  = YAQBP(jE) - YQBP
        YD    = 0.d0
        jnL   = jn1
        jnR   = jn1

        if (ARD(jE, JI) <= YF0 .and. ARD(jE, n1) >= YF0) then
            JTRAP = 0
            if (YFA <= 0. .and. YFA >= YR2 .or. rcr >= YRIPLR(JI)) then
                JLOSS = 1
            else
                CALL NBORBCO(jE, jn, rcr, RJ, RJT, yc2, YF0, jnR, jnL, JTRAP, JLOSS, ycos)
                if (JLOSS == 0) then
                    YJSERF = 1./(jnR - jnL + 1.D0)
                    YD     = YDDD * YJSERF
                endif
            endif
        else
            JTRAP = 1
            if (YFA <= YR1 .and. YFA >= YR2) then
                JLOSS = 0
            else
                CALL NBORBCO(jE, jn, rcr, RJ, RJT, yc2, YF0, jnR, jnL, JTRAP, JLOSS, ycos)
                if (JLOSS == 0) then
                    YJSERF = 1./(jnR - jnL + 1)
                    YD     = YDDD * YJSERF
                endif
            endif
        endif

        YCT2(jE) = yc2
        YDDD  = -2.d0*CONTR*YDDD*az_jpencil*YJSERF

        do j_jn=jnL, jnR
            YATBA(jE, j_jn) = YATBA(jE, j_jn) + YDDD
            YANBA(jE, j_jn) = YANBA(jE, j_jn) + YD
            YANBA1(jE, jn1) = YANBA1(jE, jn1) + YD
            YAQBA(jE, j_jn) = YAQBA(jE, j_jn) + YD
            if (YTRAP(j_jn) <= YCOS(j_jn)) then
                YACBA(jE, j_jn) = YACBA(jE, j_jn) + YD * 2.d0 * YCOS(j_jn)
            endif
            YAQBP(jE) = YQBP
            if (YCOS(j_jn) > 0.d0) then
                j_the = ntet * YCOS(j_jn)
                if (j_the > ntet) j_the = ntet
                if (j_the < 1)    j_the = 1
                YASBA(jE, j_jn, j_the) = YASBA(jE, j_jn, j_the) + YD
            else
                j_the = -ntet * YCOS(j_jn)
                if (j_the > ntet) j_the = ntet
                if (j_the < 1)    j_the = 1
                YASBA1(jE, j_jn, j_the) = YASBA1(jE, j_jn, j_the) + YD
            endif
        enddo

    enddo

    end subroutine nbco2gc

!---------------------------------------------------------------------
    subroutine nborbco(jE, jn, rcr, rj, rjt, yc2, yf0, jn_r, jn_l, itrap, iloss, ycos)

!----- First orbit analysis taking account ----
! the input source redistribution
! due to the first orbit deviation (guide centre apprx.)
! co_injection at birth
! input:
!    jE, jn, yc2, YF0	starting point of
!    velocity index, surface index, yc2=<v.B>/B
!    rcr, RJ, RJT     = critical, external, internal boundary radii
!    YF0            = F value at the birth position
!    RMAJ_LFS/RMAJ_HFS  = major radius of (j) surface with mid plane LFS/HFS
! output:
!    jn_l, jn_r   = minimum/maximum surface index of the orbit with finit RLM
!    ILOSS      = 0/1 if particle is kept/lost
!    ITRAP      = 0/1 for banana/passing orbits
!    YCOS(jjn)  = <v.B>/vB (x(jjn))

    integer, intent(in) :: jE, jn
    double precision, intent(in) :: rcr, rj, rjt, yc2, yf0
    integer, intent(out) :: itrap, iloss, jn_r, jn_l
    double precision, intent(out), dimension(n_rho) :: ycos

    logical :: found_root
    integer :: j
    double precision :: Y, YY, Rmaj_lfs, Rmaj_hfs, yrcr, yrcr_sq

    jn_r = jn
    jn_l = jn

    itrap = 1 !(def.passing)
    iloss = 0 !(def. kept)
    YCOS(jn) = yc2
    if (jn > n1) then
        write(*, *) 'error in nbrbco: jn > n1: ', jn, '>', n1
        return
    endif
 
    yrcr = 0.5d0 * rcr
    yrcr_sq = yrcr**2

    found_root = .false.
    do j=jn, n1
        Rmaj_lfs = X(j) + aspect_ratio + DX(j)
        y = ard(jE, j) - yf0  ! motion to LFS (R->RJ)
        yy = yrcr + sqrt(yrcr_sq + y**2) ! R=Rc/2+sqrt((Rc/2)^2+y^2)
        if (yy > Rmaj_lfs) then
           found_root = .true.
           EXIT
        endif
        if (yy < rjt .or. yy > rj) then ! lost
            iloss = 1
            return
        endif
        ycos(j) = max(-1.d0, y/yy)
        ycos(j) = min( 1.d0, ycos(j))
        jn_r = j
    enddo
    if (.not. found_root) then
        iloss = 1
        return
    endif

    found_root = .false.
    do j=jn, 1, -1
        y  = ard(jE, j) - yf0
        yy = yrcr + sqrt(yrcr_sq + y**2)
        Rmaj_lfs =  X(j) + aspect_ratio + DX(j)
        Rmaj_hfs = -X(j) + aspect_ratio + DX(j)
        if (yy > Rmaj_lfs .or. yy < Rmaj_hfs) then
            found_root = .true.
            EXIT
        endif
        if (y < 0.d0) itrap = 0
        if (yy < rjt .or. yy > rj) then
            iloss = 1
            return
        endif
        ycos(j) = max(-1.d0, y / yy)
        ycos(j) = min( 1.d0, ycos(j))
        jn_l = j
    enddo
    if (.not. found_root) then
        iloss = 0
    endif

    end subroutine nborbco

end module nbicom
