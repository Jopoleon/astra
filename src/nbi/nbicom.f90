module nbicom

use nbstatus, only: n_rho, n_theta, n_energy

implicit none

integer, parameter :: ndim1=82, ndim2=201
double precision, dimension(n_theta), parameter :: PLEJ2 = (/ &
    -5.000000E-01, -4.994000E-01, -4.976000E-01, -4.946000E-01, &
    -4.904000E-01, -4.850000E-01, -4.784000E-01, -4.706000E-01, &
    -4.616000E-01, -4.514000E-01, -4.400000E-01, -4.274000E-01, &
    -4.136000E-01, -3.986000E-01, -3.824000E-01, -3.650000E-01, &
    -3.464000E-01, -3.266000E-01, -3.056000E-01, -2.834000E-01, &
    -2.600000E-01, -2.354000E-01, -2.096000E-01, -1.826000E-01, &
    -1.544000E-01, -1.250000E-01, -9.440002E-02, -6.260002E-02, &
    -2.960002E-02,  4.599977E-03,  3.999998E-02,  7.659998E-02, &
     1.144000E-01,  1.534000E-01,  1.936000E-01,  2.350000E-01, &
     2.776000E-01,  3.214000E-01,  3.664000E-01,  4.126000E-01, &
     4.599999E-01,  5.085999E-01,  5.584000E-01,  6.094000E-01, &
     6.615999E-01,  7.150000E-01,  7.695999E-01,  8.253999E-01, &
     8.823999E-01,  9.405999E-01,  9.999999E-01/)

integer :: IHB, IRB, IEB, JEB, N1, JTANG, ntet1, JNR, JNL, II(ndim2)
double precision :: Rmaj, amin, BZ, HB, RBMIN1, RBMAX1, EB, QB, PB, Hocu, &
    RCR, RJ, RJT, YDYH
double precision, dimension(n_energy) :: SVII, VNB, YAQB
double precision, dimension(n_rho) :: RC, RE, RI, DRE, DRI, YCOS, YBTDB, &
    REJ, RIJ, YDZ, YTRAP, X, XJ, DX, ELON1, TRIA1
double precision, dimension(ndim1) :: AZ, YDRY
double precision, dimension(4, 4) :: SVEX
double precision, dimension(n_energy, n_rho) :: ARD, YFI, F, YFE, YANBA, &
    YAQBA, YACBA, DTCX, YATBA, YANBA1, YSLEJ0, YSLEJ2
double precision, dimension(n_energy, n_rho, n_theta) :: YASBA, YASBA1
double precision, dimension(9) :: ZB, RMB
double precision :: YSIMPI(3, 9)
double precision, dimension(n_rho) :: stnbdp, sdnbtp, sdnbdp1, sdnbdp2

contains

!---------------------------------------------------------------------
    double precision function power_dist(x_in, C1, C2) 
!----- Vertical/horizontal beam power distribution in the footprint -----
! The input parameters X1 = CVER1/CHOR1, C2 = CVER2/CHOR2 are given in
!     the NBI configuration file for each NBI source (last line)
! x_in is the normalized height/width in respect of the footprint center
! x_in = (H - HBEAM)/H_half_width,   |Z|<=1
! x_in = (R - RBtang)/R_half_width,  |Y|<=1
! H_half_width = CBMS4*(RBMAX-RBMIN)/2
! power_dist [a.u.] is renormalized in the calling routine
!---------------------------------------------------------------------

    double precision, intent(in) :: x_in, C1, C2

    power_dist = EXP(-C1*abs(x_in)**C2)

    end function power_dist

!---------------------------------------------------------------------
    subroutine nb0(adqb, jsrc, yhm, ch1, ch2, cr1, cr2, cs3)

    use nbstatus, only: n_rho

    integer, intent(in) :: jsrc
    double precision, intent(in) :: yhm, ch1, ch2, cr1, cr2, cs3, adqb(*)

    integer :: j, jh, jr, jhh, je, n, jh1, jhb, jhb05, jedge
    double precision :: ys, z, zy, y1, zy12, yr1, yr2, y, yc, yd, ye, &
        ydr, yr0, yz1, yds, ysb, yhbd2, yhocu, ynocu, &
        yhtop, yhtopa, yhedge, yhedga, yd05h, ydzout, ydhz, ydzjh

    n = n1 - 1
    yr1 = rbmin1/amin
    yr2 = rbmax1/amin

    if (yr2 <= yr1) then
        write(*, *) 'for nbi source n ', jsrc, &
                   ' rmax = ', rbmax1, ' < rmin = ', rbmin1
        write(*, *) 'continue with (rmax - rmin)/a = 1.d-5'
        yr2 = yr1 + 1.d-5
    endif

    yc = Rmaj/amin
    y  = yc + 1.d0
    yhbd2 = abs(hb)/2.d0

    if (yhbd2 == 0.d0) then
        write(*, *) 'nbi source n ', jsrc, ' vertical size = 0'
        write(*, *) 'continue with vertical size = 1.e-5'
        yhbd2 = 1.d-5
    endif

    ydr = (yr2 - yr1) / irb
    yd  = 2.d0 / abs(yr2 - yr1)
    yr0 = 0.5d0 * (yr1 + yr2)
    hocu = yhm
    do jh=1, n1
        ydz(jh) = 0.d0
    enddo

! perpendicular nbi: calculation of the nbi height

    y = hocu/amin
    yhocu = y

    if (yr0 < (yc + dx(1))) then
        do j=1, n1
            yhocu = y + cs3 * sqrt((yc + dx(j) + yr0) * (yc + dx(j) - yr0))
            if (abs(yhocu) >= x(j) * elon1(j)) exit
        enddo
    endif

    hocu = amin*yhocu
    ys = 0.d0

    do jr=1, irb
        az(jr) = (jr - 0.5d0) * ydr + yr1
        y = (abs(az(jr) - yr0) + 1.d-7) * yd
        ydry(jr) = power_dist(y, cr1, cr2)
        ys = ys + ydry(jr)
    enddo

    ys = irb / ys
    do jr=1, irb
        ydry(jr) = ydry(jr) * ys
    enddo

    yhocu  = abs(hocu)
    yhtop  = yhocu + yhbd2
    yhtopa = yhtop/amin
    ynocu  = yhocu/yhbd2
    yhedga = x(n1)*elon1(n1)
    yhedge = amin*yhedga
    jhb05  = 0

    if (yhedge < (yhocu - yhbd2)) then ! nbi out of plasma
        ihb = 0
        return
    endif

    yd05h = amin/yhbd2
    ys = 0.d0

    if (yhtop > yhedge) then ! nbi partly out of plasma
        jedge = (yhtop / yhedge - 1.d0) * (n1 - 1)
        if (jedge < 1) jedge = 1
        ydzout = (yhtop - yhedge) / (jedge*amin)
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
        if (yz1 <= yhtopa) exit
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

    yds = yhbd2 / (amin*ys)
    do jh=1, jhb05
        ydz(jh) = ydz(jh) * yds
    enddo

    ihb = 2*jhb05
    ysb = (yr2 - yr1)*amin*2.d0*yhbd2
    y  = 2.d0*451.9d0*qb/ysb
    y1 = eb/pb

    do je=jeb, ieb
        ye = y1/(ieb - je + 1)
        yaqb(je) = y*adqb(je)/(ye*pb*sqrt(ye))
    enddo

    end subroutine nb0

!---------------------------------------------------------------------
    subroutine nb1tr(yaqbp, contr, ar)
!---- ripple normalized radius of ripple boundary
!---- banana with rtrap > yriplr is lost

    use nbstatus, only: n_energy

    double precision, intent(in) :: contr
    double precision, intent(in) :: ar(*)
    double precision, intent(inout) :: yaqbp(*)

    integer :: ntet, ji, jn, jt, jr, jjn, jj, jhb05, jn1, in, in1, &
        jjh, jh, jh1, jjr, jbb, jnn, jhp1, je, jeth, jt1, jt2, n
    double precision :: ydedj, y0, y1, y2, y, ye, yh, ds, dv, &
        yqshth, yloss, yvede, ye3, rj0, yr, ydh, xjh, &
        zjh, yzjh, y12, yr1, yr2, yc1, yc2
    double precision, dimension(n_energy) :: drl, ycu, yct2, yve

    yqshth = 0.d0
    ntet = ntet1 - 1
    n    = n1 - 1
    y0 = Rmaj/amin
    ds = abs(rbmax1 - rbmin1)*amin/(n*irb)
    dv = ds*amin
    yloss = 0.d0

! -------- clean the sources --------
    do jn=1, n1
        do je=jeb, ieb
            yanba(je, jn)  = 0.d0
            yanba1(je, jn) = 0.d0
            yacba(je, jn)  = 0.d0
            yaqba(je, jn)  = 0.d0
            yatba(je, jn)  = 0.d0
            do jt=1, ntet1
                yasba(je, jn, jt)  = 0.d0
                yasba1(je, jn, jt) = 0.d0
            enddo
        enddo
        rc(jn) = y0 + dx(jn)
    enddo

    if (ihb == 0) return

! nbi out of plasma

    jhb05 = ihb / 2
    ye3   = eb / pb

! for calculation of surface index ii = jn(x(rcrit)) in trapping analysis
! yr = (rj - rii/ab) normalised distance from ext. boundary rj

    rj   = y0 + 1.d0
    rjt  = y0 - 1.d0
    rj0  = y0 + dx(1)
    ydyh = 100.0d0
    yh   = 0.01d0
    jn1 = n1

    do ji=1, 200
        yr = yh * (ji - 1)

 81     continue
        if (jn1 > 1) then
            jn = jn1 - 1
            if (yr >= (x(n1) + dx(n1) - dx(jn1) - x(jn1)) .and. yr < (x(n1) + dx(n1) - dx(jn) - x(jn))) then
                ii(ji) = jn
            else
                jn1 = jn
                goto 81
            endif
        else
            if (jn1 == 1) jn1 = -1

 82         continue
            in1 = -jn1
            in  = in1 + 1
            if (yr >= (x(n1) + dx(n1) + x(in1) - dx(in1)) .and. yr < (x(n1) + dx(n1) + x(in) - dx(in))) then
                ii(ji) = in1
            else
                jn1 = jn1 - 1
                goto 82
            endif

        endif
    enddo

    ii(ndim2) = n
    do je=jeb, ieb
        drl(je) = 1.d0 / (0.0144d0*sqrt(eb*pb/(ieb - je + 1)) / (bz*amin))
        do jn = 1, n1
            ard(je, jn) = ar(jn)*drl(je)
        enddo
    enddo

! === end trapping arrays

    do je=jeb, ieb
        YE = YE3/(IEB - JE + 1)
        YVEDE = DS*sqrt(YE)/451.9d0
        YVE(JE) = YVEDE*PB*YE*1.d-3
        YCU(JE) = -CONTR*VNB(JE)*YVEDE*5.d-6
    enddo

! loop beam`s height H = X(JN)+dX/2
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
            re(jh) = 0.5d0 * ( rc(jh1) + rc(jh) - tria1(jh) - tria1(jh1) )
        endif

        zjh = xjh * elon1(jh)
        ri(jh) = re(jh)

        do jjr = jh, n
            jn = n1 - jjr + jh
            yzjh = zjh / elon1(jn)
            y = sqrt((x(jn) - yzjh) * (x(jn) + yzjh))
            re(jn) = rc(jn) + y - tria1(jn) * (yzjh / x(jn))**2
            ri(jn) = rc(jn) - y - tria1(jn) * (yzjh / x(jn))**2
        enddo

        do jn = n, jh + 1, -1
            jj = jn - 1
            if (re(jn) < re(jj)) re(jj) = re(jn) - (re(jj) - re(jn))/n
            if (ri(jn) > ri(jj)) ri(jj) = ri(jn) - (ri(jj) - ri(jn))/n
        enddo

        do jn = jh + 1, n1
            if (jn > 1 .and. ri(jn-1) < ri(jn)) then
                write(*,*) jn, jh, jhb05, n1
            endif
        enddo

        do jn = jh, n1
            dre(jn) = 1.d0 / re(jn)
            dri(jn) = 1.d0 / ri(jn)
        enddo

! ZZ 22 - loop beam radius R = Z(JR)

        do jr=1, irb
            if (yaqb(3) <= 0.d0) EXIT
            jbb = jeb
            do je=jbb, ieb
                yaqbp(je) = yaqb(je) * ydh * ydry(jr)
            enddo

! integral along beam Y

            y12 = (re(n1) - az(jr)) * (re(n1) + az(jr))

            if (y12 <= 0.d0) EXIT

            y2  = -sqrt(y12)
            yc2 = az(jr) * dre(n1)
            jjn = n1 - jh
            rcr = y12 * dre(n1)
            yr2 = contr * sqrt(rj * (rj - rcr))

            if (rcr < rjt) then
                yr1 = contr * sqrt(rjt * (rjt - rcr))
            else
                yr1 = 0.d0
            endif

            do je = jeb, ieb
                yct2(je) = 0.d0
            enddo

            if (contr < 0.d0) then

                do jj=1, jjn
                    jn  = n1 - jj
                    jn1 = jn + 1
                    if (y2 == 0.d0) goto 226
                    y12 = (re(jn) - az(jr)) * (re(jn) + az(jr))
                    y1  = y2
                    y2  = 0.d0
                    yc1 = yc2
                    if (y12 > 0.d0) then
                        y2  = -sqrt(y12)
                        yc2 = az(jr) * dre(jn)
                        rcr = y12 * dre(jn)
                    else
                        rcr = 0.d0
                        yc2 = 1.d0
                    endif
                    jt1 = ntet * yc1 + 1
                    jt2 = ntet * yc2 + 1
                    jt1 = min(jt1, ntet)
                    jt2 = min(jt2, ntet)
                    call nbco3(jn, jn1, JT1, JT2, N1, jbb, ieb, y1, y2, yc2, AZ(jr), contr, yct2, yaqbp)
                enddo

                jn = jh - 1

                if (y2 == 0.d0) goto 226

                do jn=jh, n
                    jn1 = jn + 1
                    if (y2 == 0.d0) goto 224
                    y12 = (ri(jn1) - az(jr)) * (ri(jn1) + az(jr))
                    y1  = y2
                    y2  = 0.d0
                    yc1 = yc2
                    if (y12 > 0.d0) then
                        y2  = -sqrt(y12)
                        yc2 = az(jr) * dri(jn1)
                        rcr = y12 * dri(jn1)
                    else
                        yc2 = 1.d0
                        rcr = 0.d0
                    endif
                    jt1 = ntet * yc1 + 1
                    jt2 = ntet * yc2 + 1
                    jt1 = min(jt1, ntet)
                    jt2 = min(jt2, ntet)
                    call nbco3(jn, jn1, JT1, JT2, N1, jbb, ieb, y1, y2, yc2, AZ(jr), contr, yct2, yaqbp)
                    if (yc1 > yc2) then
                        write(*, *) 'YC1, YC2, 231', yc1, yc2, jn, jh
                    endif
                enddo

                if (jn1 == n1 .and. jtang == 0) goto 220

224             continue

                jnn  = jn1 - 1
                jhp1 = jh + 1

                do jjn=jhp1, jnn
                    jn1 = jnn - jjn + jhp1
                    jn  = jn1 - 1
                    y1  = y2
                    y12 = (ri(jn) - az(jr)) * (ri(jn) + az(jr))
                    y2  = sqrt(y12)
                    yc1 = yc2
                    yc2 = az(jr) * dri(jn)
                    rcr = y12 * dri(jn)
                    jt2 = ntet * yc1 + 1
                    jt1 = ntet * yc2 + 1
                    jt1 = min(jt1, ntet)
                    jt2 = min(jt2, ntet)
                    call nbco3(jn, jn1, JT1, JT2, N1, jbb, ieb, y1, y2, yc2, AZ(jr), contr, yct2, yaqbp)
                    if (yc1 < yc2) then
                        write(*,*) 'YC1, YC2, 251', yc1, yc2, jn
                    endif
                enddo

                jn = jh - 1

 226            continue

                jnn = jn + 1

                do jn=jnn, n
                    jn1 = jn + 1
                    y1 = y2
                    y12 = (re(jn1) - az(jr)) * (re(jn1) + az(jr))
                    y2 = sqrt(y12)
                    yc1 = yc2
                    yc2 = az(jr) * dre(jn1)
                    rcr = y12 * dre(jn1)
                    jt2 = ntet * yc1 + 1
                    jt1 = ntet * yc2 + 1
                    jt1 = min(jt1, ntet)
                    jt2 = min(jt2, ntet)
                    call nbco3(jn, jn1, JT1, JT2, N1, jbb, ieb, y1, y2, yc2, AZ(jr), contr, yct2, yaqbp)
                enddo
            else
! loop from R = RE(A) to R = RC (RE(0)) or to R(Y=0)
                do JJ=1, JJN 
                    JN = N1 - JJ
                    JN1 = JN + 1
                    if (Y2 == 0.) goto 1226  ! R = R(Y=0)
                    Y12 =(RE(JN) - AZ(JR))*(RE(JN) + AZ(JR))
                    Y1 = Y2
                    Y2 = 0.
                    YC1 = YC2
                    if (Y12 > 0.) then
                        Y2 = -SQRT(Y12)
                        YC2 = AZ(JR)*DRE(JN)
                        RCR = Y12*DRE(JN)
                    else
                        RCR = 0.
                        YC2 = 1.
                    endif
                    JT1 = min(int(ntet*YC1) + 1, ntet)
                    JT2 = min(int(ntet*YC2) + 1, ntet)
                    call nbctr3(jn, jn1, JT1, JT2, n1, jbb, ieb, y1, y2, yc2, AZ(jr), contr, yct2, yaqbp)
                enddo

                JN = JH - 1
                if (Y2 == 0.) goto 1226  ! beam touch R = R(0) (R=RC)
! loop from R = R(0) (R=RC) to R = R(Y=0) or to R = RI(A)
                do JN=JH, N
                    JN1 = JN + 1
                    if (Y2 == 0.) goto 1224
                    Y12 = (RI(JN1) - AZ(JR))*(RI(JN1) + AZ(JR))
                    Y1 = Y2
                    Y2 = 0.
                    YC1 = YC2
                    if (Y12 > 0.) then
                        Y2 = -SQRT(Y12)
                        YC2 = AZ(JR)*DRI(JN1)
                        RCR = Y12*DRI(JN1)
                    else
                        YC2 = 1.
                        RCR = 0.
                    endif
                    JT1 = min(int(ntet*YC1) + 1, ntet)
                    JT2 = min(int(ntet*YC2) + 1, ntet)
                    call nbctr3(jn, jn1, JT1, JT2, n1, jbb, ieb, y1, y2, yc2, AZ(jr), contr, yct2, yaqbp)
                enddo

                if (JN1 == N1 .and. JTANG == 0) goto 220 ! beam to internal wall

1224            continue

! R increases
! loop from R = R(Y=0) to R =R(0) (R=RC)
                JNN = JN1 - 1
                JHP1 = JH + 1
                do JJN=JHP1, JNN
                    JN1 = JNN - JJN + JHP1
                    JN = JN1 - 1
                    Y1 = Y2
                    Y12 = (RI(JN) - AZ(JR))*(RI(JN) + AZ(JR))
                    Y2 = SQRT(Y12)
                    YC1 = YC2
                    YC2 = AZ(JR)*DRI(JN)
                    RCR = Y12*DRI(JN)
                    JT2 = min(int(ntet*YC1) + 1, ntet)
                    JT1 = min(int(ntet*YC2) + 1, ntet)
                    call nbctr3(jn, jn1, JT1, JT2, n1, jbb, ieb, y1, y2, yc2, AZ(jr), contr, yct2, yaqbp)
                enddo

                JN = JH - 1

1226            continue

                JNN = JN + 1
! loop from R = R(0) (R=RC) or from R = R(Y=0) to R = RE(A)
                do JN=JNN, N
                    JN1 = JN + 1
                    Y1 = Y2
                    Y12 = (RE(JN1) - AZ(JR))*(RE(JN1) + AZ(JR))
                    Y2 = SQRT(Y12)
                    YC1 = YC2
                    YC2 = AZ(JR)*DRE(JN1)
                    RCR = Y12*DRE(JN1)
                    JT2 = min(int(ntet*YC1) + 1, ntet)
                    JT1 = min(int(ntet*YC2) + 1, ntet)
                    call nbctr3(jn, jn1, JT1, JT2, n1, jbb, ieb, y1, y2, yc2, AZ(jr), contr, yct2, yaqbp)
                enddo

            endif

            write(*, *) 'passed'

220         continue

            do jeth=jeb, ieb
                yqshth = yqshth + yaqbp(jeth) * yve(jeth)
            enddo

        enddo

    enddo height_loop

    do JE=JEB, IEB
        do JN=1, N1
            YAQBA (JE, JN) = YAQBA (JE, JN)*YVE(JE)
            YANBA (JE, JN) = YANBA (JE, JN)*YVE(JE)
            YANBA1(JE, JN) = YANBA1(JE, JN)*YVE(JE)
            YACBA (JE, JN) = YACBA (JE, JN)*YCU(JE)
            YATBA (JE, JN) = YATBA (JE, JN)*YCU(JE)
        enddo
    enddo
    YDEDJ = 1.d6/(1.6d0*EB)

    do JE=JEB, IEB
        Y = (IEB - JE + 1)*YDEDJ*YVE(JE)
        do JN=1, N1
            YSLEJ0(JE, JN) = 0.0d0
            YSLEJ2(JE, JN) = 0.0d0
            do JT=1, ntet1
                if (YASBA(JE, JN, JT) > 0.) then
                    YASBA(JE, JN, JT) = YASBA(JE, JN, JT)*Y
                    YSLEJ2(JE, JN) = YSLEJ2(JE, JN) + YASBA(JE, JN, JT)*2.5d0*PLEJ2(JT)
                    YSLEJ0(JE, JN) = YSLEJ0(JE, JN) + 0.5d0*YASBA(JE,JN,JT)
                endif
                if (YASBA1(JE, JN, JT) > 0.) then
                    YASBA1(JE, JN, JT) = YASBA1(JE, JN, JT)*Y
                    YSLEJ2(JE, JN) = YSLEJ2(JE, JN) + YASBA1(JE, JN, JT)*2.5d0*PLEJ2(JT)
                    YSLEJ0(JE, JN) = YSLEJ0(JE, JN) + 0.5d0*YASBA1(JE, JN, JT)
                endif
            enddo
        enddo
    enddo

    end subroutine nb1tr

!---------------------------------------------------------------------
    subroutine nbtrag(YAQBP, CONTR, AR)
!--------- Neutral beam ionization -----------
! coinjection  with ion's trapping and orbital losses
!  AQBA[MW],ACBA[MA m/s],ANBA[10*13 prtcls],
!  ASBA[10#19 prtcl/s]*Dcos(JT)
!  YQSHth [MW] shine through power
!-------------------------------------- Polevoy

    double precision, intent(in) :: CONTR, AR(*)
    double precision, intent(out) :: YAQBP(*)

    integer :: je, ntet, n, jn, jt, jhb05, jn1, ji, in, in1, jh, jh1, jjh, &
        jjr, jr, jbb, jjn, jj, jt1, jt2, jnn, jhp1, jeth
    double precision :: YQSHTH, Y0, DS, DV, YLOSS, YE3, RJ0, YH, YR, YRN, &
        YE, YVEDE, Y, YDH, XJH, ZJH, YZJH, Y12, Y2, YC2, YRN1, YR2, YR1, &
        YDEDJ, YRBJN, Y1, YC1, YDYS, YDEX
    double precision, dimension(3) :: DRL, YCU, YCT2, YVE

    YQSHTH = 0.
    ntet = ntet1 - 1
    N    = N1 - 1
    Y0 = Rmaj/amin
    DS = abs(RBMAX1 - RBMIN1)*amin/(N*IRB)
    DV = DS*amin
    YLOSS = 0.

! clean the sources
    do JN=1, N1
        do JE=JEB, IEB
            YANBA (JE, JN) = 0.
            YANBA1(JE, JN) = 0.
            YACBA (JE, JN) = 0.
            YATBA (JE, JN) = 0.
            YAQBA (JE, JN) = 0.
            do JT=1, ntet1
                YASBA (JE, JN, JT) = 0.d0
                YASBA1(JE, JN, JT) = 0.d0
            enddo
        enddo
        RC(JN) = Y0 + DX(JN)
    enddo

    if (IHB == 0) return ! NBI out of plasma

    JHB05 = IHB/2
    YE3 = EB/PB

! for calculation of surface index II=JN(X(Rcrit)) in trapping analysis
! YR=(RJ-Rii/AB) normalised distance from ext. boundary RJ
! Trapping analysis with gyroradius and surface averaging
! Rlm = Gyrorad(in the midplain)*YDYH/sqrt(RCR)
    RJ   = Y0 + 1.d0
    RJT  = Y0 - 1.d0
    RJ0  = Y0 + DX(1)
    YDYH = 100.0
    YH   = 0.01
    JN1  = N1
    do JE=JEB, IEB
        DRL(JE) = 1./(0.0144d0*sqrt(EB*PB/(IEB - JE + 1))/(BZ*amin))
        do JN=1, N1
            ARD(JE, JN) = AR(JN)*DRL(JE)
        enddo
    enddo

    do JI =1, 200
        YR = YH*(JI - 1)
81      continue
        if (JN1 > 1) then
            JN = JN1 - 1
            YRN1 = Y0 + DX(JN1) + X(JN1)
            YRN  = Y0 + DX(JN)  + X(JN)
            if (YR >= (X(N1)+DX(N1)-DX(JN1)-X(JN1)) .and.YR < (X(N1)+DX(N1)-DX(JN)-X(JN))) then
                II(JI) =JN
            else
                JN1 =JN
                goto 81
            endif

! For orbit averaging and gyrolosses
        else
            if (JN1 == 1) then
                JN1 = -1
            endif
82          continue
            IN1 = -JN1
            IN = IN1 + 1

            YRN1 = Y0 + DX(IN1) - X(IN1)
            YRN  = Y0 + DX(IN)  - X(IN)
            if (YR >= (X(N1)+DX(N1)+X(IN1)-DX(IN1)) .and. YR < (X(N1)+DX(N1)+X(IN)-DX(IN))) then
                II(JI) = IN1
            else
                JN1 = JN1 - 1
                goto 82
            endif
! For orbit averaging and gyrolosses
        endif
    enddo

    II(ndim2) = N

    do JE=JEB, IEB
        YE = YE3/(IEB - JE + 1)
        YVEDE = DS*SQRT(YE)/451.9
        YVE(JE) = YVEDE*PB*YE*1.E-3 
        YCU(JE) = -CONTR*VNB(JE)*YVEDE*5.E-6
    enddo

! loop beam`s height H = X(JN)+dX/2
    height_loop: do JJH=1, JHB05
        JH = JHB05 - JJH + 1
        JH1 = JH + 1
        Y = ELON1(JH1)*X(JH1) - ELON1(JH)*X(JH)
        if (JH == JHB05) then
            YDH = N*YDZ(JH)
            XJH = X(JH) + 0.5*YDH/(N*ELON1(JH))
            RE(JH) = RC(JH) - (RC(JH1) - RC(JH))*(XJH - X(JH))/(X(JH1) - X(JH)) - TRIA1(JH)
        else
            XJH = X(JH) + 0.5*Y/ELON1(JH)
            YDH = N*YDZ(JH)
            RE(JH) = 0.5*(RC(JH1) + RC(JH) - TRIA1(JH) - TRIA1(JH1)) 
        endif
        ZJH = XJH*ELON1(JH)
        RI(JH) = RE(JH)
        do JJR=JH, N
            JN = N1 - JJR + JH
            YZJH = ZJH/ELON1(JN)
            Y = SQRT((X(JN) - YZJH)*(X(JN) + YZJH))
            RE(JN) = RC(JN) + Y - TRIA1(JN)*(YZJH/X(JN))**2
            RI(JN) = RC(JN) - Y - TRIA1(JN)*(YZJH/X(JN))**2
        enddo
        do JN=N, JH+1, -1
            JJ = JN - 1
            if (RE(JN) < RE(JJ)) RE(JJ) = RE(JN) - (RE(JJ) - RE(JN))/N
            if (RI(JN) > RI(JJ)) RI(JJ) = RI(JN) - (RI(JJ) - RI(JN))/N
        enddo

        do JN=JH+1, N1    
            if (jn > 1 .and. RI(JN-1) < RI(JN)) then
                write(*, *) JN
                write(*, *) JH
                write(*, *) JHB05
                write(*, *) N1
            endif
        enddo

        do JN =JH, N1
            DRE(JN) = 1./RE(JN)
            DRI(JN) = 1./RI(JN)
        enddo
 
!  loop beam`s radius R = Z(JR)
        do JR=1, IRB
            if (YAQB(3) <= 0.) EXIT ! no power in this pencil
            JBB = JEB
            do JE=JBB, IEB
                YAQBP(JE) = YAQB(JE)*YDH*YDRY(JR)
            enddo

!  R  decreases
            Y12 = (RE(N1) - AZ(JR))*(RE(N1) + AZ(JR))
            if (Y12 <= 0.) EXIT ! beam beyond tokamak
            Y2 = -SQRT(Y12)
            YC2 = AZ(JR)*DRE(N1)
            JJN = N1 - JH
            RCR = (RE(N1) - AZ(JR))*(RE(N1) + AZ(JR))*DRE(N1)
            YR2 = CONTR*sqrt(RJ*(RJ - RCR))
            if (RCR < RJT) then
                YR1 = CONTR*sqrt(RJT*(RJT - RCR))
            else
                YR1 = 0.
            endif
            do JE=JEB, IEB
                YCT2(JE) = 0.
            enddo

            if (CONTR < 0.) then
               do JJ=1, JJN   ! loop from R = RE(A) to R = RC (RE(0)) or to R(Y=0)
                   JN  = N1 - JJ
                   JN1 = JN + 1
                   if (Y2 == 0.) goto 226 ! R = R(Y=0)
                   Y12 = (RE(JN) - AZ(JR))*(RE(JN) + AZ(JR))
                   Y1  = Y2
                   Y2  = 0.
                   YC1 = YC2
                    if (Y12 > 0.) then
                        Y2   = -SQRT(Y12)
                        YC2  = AZ(JR)/RE(JN)
                        YRBJN = RE(JN)
                        RCR = (RE(JN) - AZ(JR)*YBTDB(JN)) * (RE(JN) + AZ(JR)*YBTDB(JN))/RE(JN)
                    else
                        RCR = 0.
                        YC2 = 1.
                    endif
                    JT1 = ntet*YC1 + 1
                    JT2 = ntet*YC2 + 1
                    YDEX = EXP(-YDYS)
                    YFE(JE, JN) = YDEX
                    call nbco2gc(jn, jn1, N1, jbb, ieb, y1, y2, yc2, ydex, AZ(jr), contr, yct2, yaqbp, ydys)
                enddo

                JN = JH - 1
                if (Y2 == 0.) goto 226 ! beam touch R = R(0) (R=RC)

                do JN=JH, N ! loop from R = R(0) (R=RC) to R = R(Y=0) or to R = RI(A)
                    JN1 = JN + 1
                    if (Y2 == 0.) goto 224
                    Y12 = (RI(JN1) - AZ(JR))*(RI(JN1) + AZ(JR))
                    Y1  = Y2
                    Y2  = 0.
                    YC1 = YC2
                    if (Y12 > 0.) then
                        Y2   = -SQRT(Y12)
                        YC2  = AZ(JR)*DRI(JN1)
                        YRBJN = RI(JN1)
                        RCR = (RI(JN1) - AZ(JR)) * (RI(JN1) + AZ(JR))/RI(JN1)
                    else
                        YC2 = 1.
                        RCR = 0.
                    endif
                    JT1 = ntet*YC1 + 1
                    JT2 = ntet*YC2 + 1
                    YDEX = EXP(-YDYS)
                    YFI(JE,JN) = YDEX
                    call nbco2gc(jn, jn1, N1, jbb, ieb, y1, y2, yc2, ydex, AZ(jr), contr, yct2, yaqbp, ydys)
                    if (YC1 > YC2) write(*, *) 'YC1,YC2,231', YC1, YC2, JN
                enddo

                if (JN1 == N1) goto 220 ! beam to internal wall : go to 220

224             continue

! R increases
                JNN  = JN1 - 1
                JHP1 = JH + 1

                do JJN=JHP1, JNN   ! loop from R = R(Y=0) to R =R(0) (R=RC)
                    JN1 = JNN - JJN + JHP1
                    JN  = JN1 - 1
                    Y1  = Y2
                    Y12 = (RI(JN) - AZ(JR))*(RI(JN) + AZ(JR))
                    Y2  = SQRT(Y12)
                    YC1 = YC2
                    YC2 = AZ(JR)*DRI(JN)
                    YRBJN = RI(JN)
                    RCR = (RI(JN) - AZ(JR))*(RI(JN) + AZ(JR))/RI(JN)
                    JT2 = ntet*YC1 + 1
                    JT1 = ntet*YC2 + 1
                    YDEX = YFI(JE, JN)
                    call nbco2gc(jn, jn1, N1, jbb, ieb, y1, y2, yc2, ydex, AZ(jr), contr, yct2, yaqbp, ydys)
                    if (YC1 < YC2) write(*,*) 'YC1,YC2,251', YC1, YC2, JN
                enddo

                JN = JH - 1

226             continue
                JNN = JN + 1

                do JN=JNN, N  ! loop from R = R(0) (R=RC) or from R = R(Y=0) to R = RE(A)
                    JN1 = JN + 1
                    Y1 = Y2
                    Y12 = (RE(JN1) - AZ(JR))*(RE(JN1) + AZ(JR))
                    Y2 = SQRT(Y12)
                    YC1 = YC2
                    YC2 = AZ(JR)*DRE(JN1)
                    YRBJN = RE(JN)
                    RCR = (RE(JN1) - AZ(JR))*(RE(JN1) + AZ(JR))/RE(JN1)
                    JT2 = ntet*YC1 + 1
                    JT1 = ntet*YC2 + 1
                    YDEX = YFE(JE,JN)
                    call nbco2gc(jn, jn1, N1, jbb, ieb, y1, y2, yc2, ydex, AZ(jr), contr, yct2, yaqbp, ydys)
                    if (YC1 < YC2) write(*,*) 'YC1,YC2,271', YC1, YC2, JN
                enddo

            else
! loop from R = RE(A) to R = RC (RE(0)) or to R(Y=0)
                do JJ=1, JJN 
                    JN = N1 - JJ
                    JN1 = JN + 1
                    if (Y2 == 0.) goto 1226  ! R = R(Y=0)
                    Y12 =(RE(JN) - AZ(JR))*(RE(JN) + AZ(JR))
                    Y1 = Y2
                    Y2 = 0.
                    YC1 = YC2
                    if (Y12 > 0.) then
                        Y2 = -SQRT(Y12)
                        YC2 = AZ(JR)*DRE(JN)
                        RCR = Y12*DRE(JN)
                    else
                        RCR = 0.
                        YC2 = 1.
                    endif
                    JT1 = min(int(ntet*YC1) + 1, ntet)
                    JT2 = min(int(ntet*YC2) + 1, ntet)
                    call nbctr3(jn, jn1, JT1, JT2, n1, jbb, ieb, y1, y2, yc2, AZ(jr), contr, yct2, yaqbp)
                enddo

                JN = JH - 1
                if (Y2 == 0.) goto 1226  ! beam touch R = R(0) (R=RC)
! loop from R = R(0) (R=RC) to R = R(Y=0) or to R = RI(A)
                do JN=JH, N
                    JN1 = JN + 1
                    if (Y2 == 0.) goto 1224
                    Y12 = (RI(JN1) - AZ(JR))*(RI(JN1) + AZ(JR))
                    Y1 = Y2
                    Y2 = 0.
                    YC1 = YC2
                    if (Y12 > 0.) then
                        Y2 = -SQRT(Y12)
                        YC2 = AZ(JR)*DRI(JN1)
                        RCR = Y12*DRI(JN1)
                    else
                        YC2 = 1.
                        RCR = 0.
                    endif
                    JT1 = min(int(ntet*YC1) + 1, ntet)
                    JT2 = min(int(ntet*YC2) + 1, ntet)
                    call nbctr3(jn, jn1, JT1, JT2, n1, jbb, ieb, y1, y2, yc2, AZ(jr), contr, yct2, yaqbp)
                enddo

                if (JN1 == N1 .and. JTANG == 0) goto 220 ! beam to internal wall

1224            continue

! R increases
! loop from R = R(Y=0) to R =R(0) (R=RC)
                JNN = JN1 - 1
                JHP1 = JH + 1
                do JJN=JHP1, JNN
                    JN1 = JNN - JJN + JHP1
                    JN = JN1 - 1
                    Y1 = Y2
                    Y12 = (RI(JN) - AZ(JR))*(RI(JN) + AZ(JR))
                    Y2 = SQRT(Y12)
                    YC1 = YC2
                    YC2 = AZ(JR)*DRI(JN)
                    RCR = Y12*DRI(JN)
                    JT2 = min(int(ntet*YC1) + 1, ntet)
                    JT1 = min(int(ntet*YC2) + 1, ntet)
                    call nbctr3(jn, jn1, JT1, JT2, n1, jbb, ieb, y1, y2, yc2, AZ(jr), contr, yct2, yaqbp)
                enddo

                JN = JH - 1

1226            continue

                JNN = JN + 1
! loop from R = R(0) (R=RC) or from R = R(Y=0) to R = RE(A)
                do JN=JNN, N
                    JN1 = JN + 1
                    Y1 = Y2
                    Y12 = (RE(JN1) - AZ(JR))*(RE(JN1) + AZ(JR))
                    Y2 = SQRT(Y12)
                    YC1 = YC2
                    YC2 = AZ(JR)*DRE(JN1)
                    RCR = Y12*DRE(JN1)
                    JT2 = min(int(ntet*YC1) + 1, ntet)
                    JT1 = min(int(ntet*YC2) + 1, ntet)
                    call nbctr3(jn, jn1, JT1, JT2, n1, jbb, ieb, y1, y2, yc2, AZ(jr), contr, yct2, yaqbp)
                enddo

            endif
 
220         continue ! There's a goto in the inc files!
! Shinethrough power distribution R[m], Z[m], Q/S[W/cm2] 
            do JETH=JEB, IEB
                YQSHTH = YQSHTH + YAQBP(JETH)*YVE(JETH)
            enddo
        enddo  ! loop on beam`s radius R = Z(JR)
    enddo height_loop

! end of loop on beam`s height
    do JE=JEB, IEB
        do JN=1, N1
            YAQBA (JE, JN) = YAQBA (JE, JN)*YVE(JE)
            YANBA (JE, JN) = YANBA (JE, JN)*YVE(JE)
            YANBA1(JE, JN) = YANBA1(JE, JN)*YVE(JE)
            YACBA (JE, JN) = YACBA (JE, JN)*YCU(JE)
            YATBA (JE, JN) = YATBA (JE, JN)*YCU(JE)
        enddo
    enddo

    YDEDJ = 1.E6/(1.6*EB)

    do JE=JEB, IEB
        Y = (IEB - JE + 1)*YDEDJ*YVE(JE)
        do JN=1, N1
            YSLEJ0(JE, JN) = 0.0
            YSLEJ2(JE, JN) = 0.0
            do JT=1, ntet1
                if (YASBA(JE, JN, JT) > 0.) then
                    YASBA(JE, JN, JT) = YASBA(JE, JN, JT)*Y
                    YSLEJ2(JE, JN) = YSLEJ2(JE, JN) + YASBA(JE, JN, JT)*2.5*PLEJ2(JT)
                    YSLEJ0(JE, JN) = YSLEJ0(JE, JN) + 0.5*YASBA(JE, JN, JT)
                endif
                if (YASBA1(JE, JN, JT) > 0.) then
                    YASBA1(JE, JN, JT) = YASBA1(JE, JN, JT)*Y
                    YSLEJ2(JE, JN) = YSLEJ2(JE, JN) + YASBA1(JE, JN, JT)*2.5*PLEJ2(JT)
                    YSLEJ0(JE, JN) = YSLEJ0(JE, JN) + 0.5*YASBA1(JE, JN, JT)
                endif
            enddo
        enddo
    enddo

    end subroutine nbtrag

!---------------------------------------------------------------------
    subroutine nbsrsr(j_nbi, YCONTR, NA1, RTOR, SHIFT, AB, BTOR, &
        HRO, YHM, CBMH1, CBMH2, CBMS1, CBMS2, CBMS3, CBMS4, &
        CBMR1, CBMR2, dn_rho, fp_flag, EBEAM, power_frac, &
        ABEAM, QBEAM, RBMAX, RBMIN, JSRREC, YEXTARR, calc_fus)

!---------------------------------------------------------------------
! fast ion's sourses (for multi sources) + ripple
! REJ,RIJ						     12-JAN-11
! YZERO -> YASBA1					     12-JUL-12
! Suzuki stopping crossections
! New SCUBM =Nt/m3					     22-APR-13
!--------------------------------------------------------------Polevoy
!	entry:	AMETR,SHIF,NA1,RTOR,AB,BTOR,NI,HBEAM,RBMIN,RBMAX,
!     		CBMS2,dn_rho,EBEAM,QBEAM,ABEAM,power_frac,
!	     	TE,TI,NE,NN,ZEF,AMAIN,PBEAM,SCUBM,CONTR,ELON,
!	CBMS2	number of 'pencils'
!	dn_rho	number of internal mesh points
!			 41 (dn_rho=1),21 (dn_rho=2)
!	power_frac3,2,1 power fraction of energy comps.
!		    	 3(EB,EB/2,EB/3),2(EB,EB/2),1(EB)
!	CONTR	Qcontr/Qbeam
!	j_nbi 	Number of the current hot ion source
!	JSRREC 	Length of the hot ion source record
!	exit:	PBEAM,SCUBM,SNEBM,SNNBM		for MAIN
!	SCUBM	Toroidal pulse [kg*m/s2/m3]	05-AUG-96
!---------------------------------------------------------------------

    use nbstatus, only: n_rho, n_energy, n_theta, ISPE, ISPEND, &
        AMETR, VR, RHO, FP, AMAIN, ELON, TRIA, SHIF, NE, TE, TI, &
        SCUBM, PBEAM, NNBM1, NNBM2, NNBM3, SNNBM, SNEBM, &
        SNIBM1, SNIBM2, SNIBM3
    use cross_sections, only: svdtbp, svddnp1, svddnp2, seiv

    integer, parameter :: JDBL=2
    double precision, parameter :: GP2=6.283185d0

    integer, intent(in) :: j_nbi, NA1, dn_rho, fp_flag, calc_fus
    double precision, intent(in) :: YCONTR, RTOR, SHIFT, AB, BTOR, &
        HRO, YHM, CBMH1, CBMH2, CBMS1, CBMS2, CBMS3, CBMS4, &
        CBMR1, CBMR2, EBEAM, power_frac(3), &
        ABEAM, QBEAM, RBMAX, RBMIN, yEXTARR(n_rho, 9)
    integer, intent(out) :: JSRREC

    double precision :: AR(n_rho), YE2, YCX
    double precision, dimension(n_energy) :: AQBP, ADQB
    double precision :: YZERO(n_energy, n_rho, n_theta), YSCU1, YFVDA, YNHDT, &
        YDBM, YJN, YSCU, YPOW, YD, YDV, YDDV, YJE, YE1, Y, YEPS
    integer :: N, J, ntet, JN, JN1, JNA, JNAX, JE, JT, JNA1, JNAC, JS, JSP

    N1 = (NA1 - 1)/dn_rho + 1
    N = N1 - 1
    if (RBMAX <= RBMIN) then
        write(*, *) 'ILLEGAL: NBI source N', j_nbi, ' RBMAX <= RBMIN !!!'
        return
    endif
    if (EBEAM*ABEAM == 0.) then
        write(*, *) 'ILLEGAL: NBI source N', j_nbi, ' EBEAM*ABEAM = 0 !!!'
        return
    endif

    IEB = 3
    IRB = CBMS2
    YDBM = SUM(power_frac)

    do j=1,  3
        ADQB(j) = 0.d0
    enddo
    if (power_frac(1) <= 0.d0) then
        write(*, *) 'ILLEGAL: NBI source N', j_nbi, ' power_frac(1) <= 0 !!!'
        return
    endif
    ADQB(3) = power_frac(1)/YDBM
    JEB = IEB
    if (power_frac(2) > 0.d0) then
        ADQB(2) = power_frac(2)/YDBM
        JEB = IEB - 1
    endif
    if (power_frac(3) > 0.d0) then
        ADQB(1) = power_frac(3)/YDBM
        JEB = IEB - 2
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
    EB = EBEAM*1000.d0
    QB = QBEAM*1000.d0
    PB = ABEAM
    Rmaj = (RTOR + SHIFT)*100.d0
    amin = AMETR(NA1)*100.d0
    if (((RBMAX + RBMIN)/2.d0) > (RTOR - AB)) then ! Tangentional NBI
        JTANG = 1
    else ! Perpendicular NBI
        JTANG = 0
    endif
    BZ = BTOR/(1. + SHIFT/RTOR)
    HB = CBMS4*(RBMAX - RBMIN)*100.*sqrt(1. + CBMS3*CBMS3)
    if (HB <= 0.)then
        write(*, *) 'NBI input for the source ', j_nbi, ' is not correct'
        write(*, *) 'please use: RBMAX > RBMIN and CBMS4 > 0'
    endif
    RBMIN1 = RBMIN*100.
    RBMAX1 = RBMAX*100.
    YD = 10000.d0/(BZ*GP2*amin**2)
    YBTDB(1) = 1.d0
    do JN=2, NA1
        JN1 = JN - 1
        AR(JN) = (FP(JN) - FP(1))*YD
        YBTDB(JN) = 1.d0
        XJ(JN) = 50.d0*(AMETR(JN) + AMETR(JN1))/amin
    enddo
    YJN = 0.d0
    Y = Rmaj/amin
    do JN=1, N1
        JNAX = 1 + dn_rho*(JN-1)
        JNA  = JNAX - dn_rho*YJN
        if (dn_rho > 1) YJN  = .5d0
        ELON1(JN) = max(ELON(JNAX), 1.d0)
        TRIA1(JN) = TRIA(JNAX)*XJ(JNAX)
        X(JN)     = XJ(JNAX)
        AR(JN)    = AR(JNA)
        DX(JN)    = 100.d0*(SHIF(JNAX) - SHIFT)/amin
        REJ(JN)   =  X(JN) + Y + DX(JN)
        RIJ(JN)   = -X(JN) + Y + DX(JN)
        YBTDB(JN) = YBTDB(JNAX)
        YEPS      = AMETR(JNAX)/(RTOR + SHIF(JNAX))
        YTRAP(JN) = sqrt(2.*YEPS/(1. + YEPS))
        do JE=1, n_energy
            YAQBA (JE, JN) = 0.d0
            YACBA (JE, JN) = 0.d0
            YATBA (JE, JN) = 0.d0
            YANBA (JE, JN) = 0.d0
            YANBA1(JE, JN) = 0.d0
            YSLEJ0(JE, JN) = 0.0d0
            YSLEJ2(JE, JN) = 0.0d0
            do JT=1, ntet1
                YASBA (JE, JN, JT) = 0.d0
                YASBA1(JE, JN, JT) = 0.d0
                YZERO (JE, JN, JT) = 0.d0
            enddo
        enddo
    enddo

    REJ(N1) = Y + 1.d0
    RIJ(N1) = Y - 1.d0

    JN = N1
    X(JN) = 1.0d0
    DX(JN) = 0.d0
    REJ(JN) =  X(JN) + Y + DX(JN)
    RIJ(JN) = -X(JN) + Y + DX(JN)

    if (QB <= 0.d0) goto 999

    call NBSISN(NA1, EBEAM, ABEAM, dn_rho, YEXTARR)
    call NB0(ADQB, j_nbi, YHM, CBMH1, CBMH2, CBMR1, CBMR2, CBMS3)

    if (CBMS1 < 1.d0) then
        call NB1TR(AQBP, YCONTR, AR)
    else
        call NBTRAG(AQBP, YCONTR, AR)
    endif
    YSCU = 2.d3/(9.79d0*sqrt(2000.d0*EBEAM/ABEAM))
    YSCU1 = ABEAM*0.0209d0*0.5d-2*amin
    YPOW = 0.d0
    JNA1 = 1
    do JN=2, N1
        JN1 = JN - 1
        if (JN > 2) JNA1 = 2 + dn_rho*(JN - 2)
        JNA = JNA1 + 1
        JNAC = 1 + dn_rho*(JN - 1)
        if (JN == N1)then
            JNAC = NA1 - 1
            YDV = VR(JNAC)*(RHO(NA1) - JNAC*HRO)
        else
            YDV = 0.d0
        endif
        if (jna1 > JNAC) JNA1 = JNAC
        JNAX = JNAC

        do J=JNA1, JNAC
            YDV = VR(j)*HRO + YDV
        enddo
        YDDV = 1./YDV
        do JE=JEB, IEB
            YSLEJ0(JE, JN1) = YSLEJ0(JE, JN)*YDDV
            YSLEJ2(JE, JN1) = YSLEJ2(JE, JN)*YDDV
            YAQBA (JE, JN1) = YAQBA (JE, JN)*YDDV
            YACBA (JE, JN1) = YACBA (JE, JN)*YDDV
            YATBA (JE, JN1) = YATBA (JE, JN)*YDDV
            YJE = EBEAM/(IEB - JE + 1)
            if (fp_flag > 0) then ! convert prtcl/s => prtcl/s/m^3
                do JT=1, ntet1
                    YASBA (JE, JN1, JT) = YASBA (JE, JN, JT)*YDDV
                    YASBA1(JE, JN1, JT) = YASBA1(JE, JN, JT)*YDDV
                enddo
            endif
            do j=JNA1, JNAC
                PBEAM(j) = PBEAM(j) + YAQBA(JE, JN1)
                if (calc_fus > 0) then ! beam-plasma fusion
                    if (ABEAM == 3.d0) stnbdp(j) = stnbdp(j) + YAQBA(JE, JN1) * &
                        svdtbp(YJE, ABEAM, NE(j), TE(j), TI(j), AMAIN(j), calc_fus)*(625.d0/YJE)
                    if (ABEAM == 2.d0) then
                        sdnbtp(j)  = sdnbtp(j)  + YAQBA(JE, JN1)*(625.d0/YJE) * svdtbp (YJE, ABEAM, NE(j), TE(j), TI(j), AMAIN(j), calc_fus)
                        sdnbdp2(j) = sdnbdp2(j) + YAQBA(JE, JN1)*(625.d0/YJE) * svddnp2(YJE, ABEAM, NE(j), TE(j), TI(j), AMAIN(j), calc_fus)
                        sdnbdp1(j) = sdnbdp1(j) + YAQBA(JE, JN1)*(625.d0/YJE) * svddnp1(YJE, ABEAM, NE(j), TE(j), TI(j), AMAIN(j), calc_fus)
                    endif
                endif
                SCUBM(j) = SCUBM(j) + YATBA(JE, JN1)*YSCU1
            enddo

            YANBA (JE, JN1) = YANBA (JE, JN)*YDDV*625.d0/YJE
            YANBA1(JE, JN1) = YANBA1(JE, JN)*YDDV*625.d0/YJE
        enddo

        do J=JNA1, JNAC
            YFVDA = F(1, JN1)*VNB(1)/amin
            if (YFVDA /= 0.d0) NNBM1(j) = NNBM1(j) + YANBA(3, JN1)/YFVDA
            YFVDA = F(2, JN1)*VNB(2)/amin
            if (YFVDA /= 0.d0) NNBM2(j) = NNBM3(j) + YANBA(2, JN1)/YFVDA
            YFVDA = F(3, JN1)*VNB(3)/amin
            if (YFVDA /= 0.d0) NNBM3(j) = NNBM3(j) + YANBA(1, JN1)/YFVDA
            YNHDT = 0.d0   !  total proton content
            do JS=2, ISPEND
                JSP = ISPE(JS)
                if (ZB(JS) == 1.) YNHDT = YNHDT + yEXTARR(JNA, JSP)
            enddo
! Calc. total particle sourse: SNEBM
            YE1 = NE(JNA)*SEIV(TE(JNA))
            do JE=JEB, IEB
                YE2 = 0.d0
                do JS=2, ISPEND
                    JSP = ISPE(JS)
                    YE2 = YE2 + yEXTARR(JNA, JSP)*YSIMPI(JE, JS)*VNB(JE)
                enddo
! Calc. thermal neutral  sourCe: SNNBM
                YCX = YNHDT*SVEX(4, JE)
                SNNBM(j) = SNNBM(j) + YCX/(YCX + YE1 + YE2)*YANBA1(JE, JN1)
                SNEBM(j) = SNEBM(j) + (1.d0 - YCX/(YCX + YE1 + YE2))*YANBA1(JE, JN1)
                if (ABEAM < 1.5d0) SNIBM1(j) = SNIBM1(j) + (1.d0 - YCX/(YCX + YE1 + YE2))*YANBA(JE, JN1) !H ion source
                if (ABEAM < 2.5d0 .and. ABEAM > 1.5d0) SNIBM2(j) = SNIBM2(j) + (1.d0 - YCX/(YCX + YE1 + YE2))*YANBA(JE, JN1) !D ion source
                if (ABEAM > 2.5d0) SNIBM3(j) = SNIBM3(j) + (1.d0 - YCX/(YCX + YE1 + YE2))*YANBA(JE, JN1) !T ion source
            enddo
        enddo
    enddo

 999 continue

    end subroutine nbsrsr

!---------------------------------------------------------------------
    subroutine nbion0(NA1, ABEAM, EBEAM, RTOR, dn_rho, yEXTARR)
!---------------------------------------------------------------------
! Steady State (1+2D:(x+MU, V)) Fokker-Plank Solver
!	PEBM, PIBM(X)-power to electrons, ions [MW/m3]
!	CUFI, CUBM(X)- fast ion's and NBI driven currents [MA/m2]
!	NIBM(x) - fast ion's density [10^19/m3]
!-------------------------------------------------------------- Polevoi

    use nbstatus, only: n_rho, ISPE, ISPEND, AMAIN, NE, TE, TI, ZEF, &
        AMETR, SHIF, PBPER, PBLON, NIBM, PEBM, PIBM, CUFI, CUBM
    use cross_sections, only: fnbf, fnb2, fnbp, fnbi1

    integer, intent(in) :: NA1, dn_rho
    double precision, intent(in) :: ABEAM, EBEAM, RTOR, yEXTARR(n_rho, 9)

    integer :: JN, JS, JSP, JN1, JNA, JNA1, JNAC, J, J2, JE
    double precision :: YC1, YC2, YS, YSTE, YEPS, YDN, Y, X1, X2, X3, &
        YA, YB, YC, YC0, Y1, Y12, YA1, YD1, YI0, YI2, YP11, YCRNT, YPIDPB, &
        STSD3
    double precision, dimension(n_rho) :: YZ2D3, YTSE, YFCUR, YLNI, YLNE, YLNZ, YEBDEC

    do JN=1, NA1
        YSTE = sqrt(TE(JN))
        if (EBEAM > 100.*ABEAM) then
            YLNI(JN) = 23.7d0 + LOG(AMAIN(JN)/(AMAIN(JN) + ABEAM)*sqrt(1.d-3*ABEAM*EBEAM*TE(JN)/NE(JN)))
        else
            YLNI(JN) = 25.4d0 + LOG(1.d-3*EBEAM*AMAIN(JN)/(AMAIN(JN) + ABEAM)*sqrt(TE(JN)/NE(JN)))
        endif
        YLNE(JN) = 15.85d0 + log(TE(JN)/sqrt(NE(JN)))
        YLNZ(JN) = YLNI(JN)
        YS = 0.
        Y1 = 0.

        do JS=2, ISPEND
            JSP = ISPE(JS)
            if (EBEAM > 100.*ABEAM) then
                YLNZ(JN) = 23.7d0 + log(RMB(JS)/(RMB(JS) + ABEAM)*sqrt(1.d-3*ABEAM*EBEAM*TE(JN)/NE(JN)))
            else
                YLNZ(JN) = 25.4d0 + log(1.d-3*EBEAM*RMB(JS)/(RMB(JS) + ABEAM)*sqrt(TE(JN)/NE(JN)))
            endif
            Y = yEXTARR(JN, JSP)*ZB(JS)**2*YLNZ(JN)
            Y1 = Y1 + Y
            YS = YS + Y/RMB(JS)
        enddo

        YEBDEC(JN) = EBEAM/(TE(JN)*ABEAM*14.8d0*(YS/(NE(JN)*YLNE(JN)))**0.6667)
        YZ2D3(JN)  = Y1/(YS*3.0d0)
        YTSE(JN)   = 2.d0*ABEAM*YSTE*TE(JN)/(NE(JN)*YLNE(JN))
        YEPS = AMETR(JN)/(RTOR + SHIF(JN))
        YFCUR(JN) = (1.d0 - FNBF(ZEF(JN), YEPS)/ZEF(JN))
        YDN = 0.d0
    enddo

    do JN=2, N1
        JN1  = JN - 1
        JNA1 = 1 + dn_rho*(JN1 - 1)
        JNA  = JNA1 + 1
        JNAC = JNA1 - 1 + dn_rho
        if (JN == N1) JNAC = NA1 - 1
        J2 = 1 + dn_rho*(JN1-1) - dn_rho*YDN
        if (dn_rho > 1) YDN = 0.5d0
        do JE=JEB, IEB
            Y = 1.d0/(IEB - JE + 1)
            do J=JNA1, JNAC
                X2 = YEBDEC(J2)*Y
                X1 = sqrt(X2)
                X3 = X2*X1
                YPIDPB = 2.0d0*FNB2(X2)/X2
                STSD3 = 208.33d0*YTSE(J2)*YAQBA(JE, JN1)
! Tail correction for current
                YB  = 1.d0 + 1.d0/X3
                YA  = 0.5*TE(J2)/(EBEAM*Y)*(1.d0 + TI(J2)/(TE(J2)*X3))
                YC0 = YTSE(J2)*DTCX(JE, J2) - 3.d0
                YC  = YZ2D3(J2)*3.d0/X3
                YC1 = YC0 + YC
                YC2 = YC1 + 2.d0*YC
                Y1  = 4.d0*YA*YC1/YB**2
                Y12 = 1.d0 + sqrt(1.d0 + Y1)
                YA1 = 0.5d0*YB*Y12/YA
                YD1 = 2.d0/(Y12 + 0.5d0*Y1)
                YI0 = (0.5d0 - FNB2(X2)/X2)
                YI2 = FNBP(X1, YZ2D3(J2))
! Beam energy Wbeam[10#19keV/m3] = Pbper + Pblon/2
                YP11 = 2.d0/3.d0*(YI0*YSLEJ0(JE, JN1) + 0.4d0*YI2*YSLEJ2(JE, JN1))
! Beam perpendicular pressure <Mb Vper2>/2 [10*19 keV/m3]
                PBPER(j) = PBPER(j) + EBEAM*Y*YTSE(J2)*(2.d0*YI0*YSLEJ0(JE, JN1) - YP11)
! Beam parallel pressure    <Mb Vpar2> [10*19 keV/m3]
                PBLON(j) = PBLON(j) + EBEAM*Y*YTSE(J2)*2.d0*YP11
! Beam density
                NIBM(j) = NIBM(j) + STSD3*LOG(1.d0 + X3)/(EBEAM*Y)
! Fast ion current without trapping correction for Eb
                YCRNT = YACBA(JE, JN1)*YTSE(J2)*YD1*FNBI1(X2, YZ2D3(J2))
                CUFI(j) = CUFI(j) + YCRNT
                PIBM(j) = PIBM(j) + YAQBA(JE, JN1)*YPIDPB
                PEBM(j) = PEBM(j) + YAQBA(JE, JN1)*(1.0d0 - YPIDPB)
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

    use nbstatus, only: n_rho, ISPE, ISPEND, NE, TE, AMETR
    use cross_sections, only: STOTQ1, seiv, spii, spex, simpi

    integer, intent(in) :: NA1, dn_rho
    double precision, intent(in) :: YEBEAM, YABEAM, YEXTARR(n_rho, 9)

    integer :: J, JN, JSP, JE, JS
    double precision :: YY, Y1, Y2, Y3, Y12, Y13, Y23, YJE

    Y3 = yEBEAM/yABEAM
    Y1 = Y3/3.d0
    Y2 = Y3/2.d0

    do JE=1, IEB
        YJE = Y3/(IEB - JE + 1)
        VNB(JE) = 1.383d6*sqrt(1.d3*YJE)
        SVII(JE) = SPII(YJE)
        SVEX(4, JE) = SPEX(YJE)
        do JS=2, ISPEND
            if (ZB(JS) > 1.) then
                YSIMPI(JE, JS) = SIMPI(YJE, ZB(JS))
            else
                YSIMPI(JE, JS) = SVEX(4, JE) + SVII(JE)
            endif
        enddo
        SVEX(4, JE) = SPEX(YJE)*VNB(JE)
    enddo

    SVEX(4, 4) = SPEX(1.d0)
    Y12 = Y1 + Y2 - 2.*sqrt(Y1*Y2)
    Y13 = Y1 + Y3 - 2.*sqrt(Y1*Y3)
    Y23 = Y2 + Y3 - 2.*sqrt(Y2*Y3)
    SVEX(1, 2) = SPEX(Y12)
    SVEX(1, 3) = SPEX(Y13)
    SVEX(2, 3) = SPEX(Y23)

    do JE=1, IEB
        YY = Y3/(IEB - JE + 1)
        do JN =1, N1
            J = (JN - 1)*dn_rho + 1
! Janev, Boley (for Eb>0.1 MeV)
            if (YY >= 100.d0) then
                F(JE, JN) = 0.d0
                do JS=2, ISPEND
                    JSP =ISPE(JS)
                    F(JE, JN) = F(JE, JN) + yEXTARR(J, JSP) * &
                         STOTQ1(yABEAM, YY, NE(j), TE(j), ZB(JS), RMB(JS)) !by Suzuki
                enddo
            else
                F(JE, JN) = NE(j)*SEIV(TE(j))/VNB(JE)
                do JS=2, ISPEND
                    JSP = ISPE(JS)
                    F(JE, JN) = F(JE, JN) + yEXTARR(J, JSP)*YSIMPI(JE, JS)
                enddo
            endif
            F(JE, JN) = F(JE, JN)*(AMETR(NA1-1) + AMETR(NA1))*50.d0
! Charge-exchange transparancy
        enddo
    enddo

    end subroutine nbsisn

!---------------------------------------------------------------------
    subroutine nbspec(AIM1, AIM2, AIM3, AMJ, ZMJ, NA1, &
        yEXTARR, NE, NI, NHYDR, NDEUT, NTRIT, NHE3, NALF, NIZ1, NIZ2, NIZ3,  &
        ZIM1, ZIM2, ZIM3, RMB, ZB, ISPEND, ISPE, IFLAG)
!---------------------------------------------------------------------
! Input (arrays only)
!   NHYDR,  NDEUT,  NTRIT,  NHE3,  NALF,  NIZ1,  NIZ2,  NIZ3,  NI
!   EXTARR(, )
! Output:
!   ISPEND    = 1 + total number of ion species
! Arrays: 1(e),  p,  d,  t,  He3,  He4,  Imp1,  Imp2,  Imp3
!   RMB(9)    m/m_p
!   ZB(9)
!   ISPE(9)   number of ion species in the EXTARR
!   EXTARR ! One of arrays for (p, d, t, He3) is spoiled !
!---------------------------------------------------------------------

    use nbstatus, only: n_rho

    double precision, intent(in) :: AMJ, ZMJ, AIM1, AIM2, AIM3
    double precision, intent(in), dimension(*) :: NE, NHYDR, NDEUT, NTRIT, NHE3, NALF, NI, &
        ZIM1, ZIM2, ZIM3, NIZ1, NIZ2, NIZ3
    double precision, intent(out), dimension(*) :: ZB, RMB
    double precision, intent(out) :: yEXTARR(n_rho, 9)

    integer :: ISPEND, ISPE(*), JIHYDR, JIDEUT, JITRIT, JIHE3, JIALF, JN, &
        JIZ1, JIZ2, JIZ3, NA1, JENE, JINI, IFLAG

! Identification of plasma species for NBI
! Polevoy A.R. 18.03.91
    JIHYDR = 0
    JIDEUT = 0
    JITRIT = 0
    JIHE3  = 0
    JIALF  = 0
    JIZ1   = 0
    JIZ2   = 0
    JIZ3   = 0
    JENE   = 0
    JINI   = 0
    ISPEND = 0
    IFLAG  = 0

    do JN=1, NA1 ! check for non zero species

        if (NE(JN) > 0.d0) then
            if (JENE < 1) then
                JENE = 1
                ISPEND = ISPEND + 1
                RMB(ISPEND) = 1./1836.d0
                ZB(ISPEND) = 1.d0
                ISPE(ISPEND) = 1
            endif
            yEXTARR(JN, 1) = NE(JN)
        endif

        if (NHYDR(JN) > 0.d0) then
            if (JIHYDR < 1) then
                JIHYDR = 1
                ISPEND = ISPEND + 1
                RMB(ISPEND) = 1.d0
                ZB(ISPEND)  = 1.d0
                ISPE(ISPEND) = 2
            endif
            yEXTARR(JN, 2) = NHYDR(JN)
        endif
 
        if (NDEUT(JN) > 0.d0) then
            if (JIDEUT < 1) then
                JIDEUT = 1
                ISPEND = ISPEND + 1
                RMB(ISPEND) = 2.d0
                ZB(ISPEND)  = 1.d0
                ISPE(ISPEND) = 3
            endif
            yEXTARR(JN, 3) = NDEUT(JN)
        endif

        if (NTRIT(JN) > 0.d0) then
            if (JITRIT < 1) then
                JITRIT = 1
                ISPEND = ISPEND + 1
                RMB(ISPEND) = 3.d0
                ZB(ISPEND)  = 1.d0
                ISPE(ISPEND) = 4
            endif
            yEXTARR(JN, 4) = NTRIT(JN)
        endif

        if (NHE3(JN) > 0.d0) then
            if (JIHE3 < 1) then
                JIHE3 = 1
                ISPEND = ISPEND + 1
                RMB(ISPEND) = 3.d0
                ZB(ISPEND)  = 2.d0
                ISPE(ISPEND) = 5
            endif
            yEXTARR(JN, 5) = NHE3(JN)
        endif

        if (NALF(JN) > 0.d0) then
            if (JIALF < 1) then
                JIALF = 1
                ISPEND = ISPEND + 1
                RMB(ISPEND) = 4.d0
                ZB(ISPEND)  = 2.d0
                ISPE(ISPEND) = 6
            endif
            yEXTARR(JN, 6)=NALF(JN)
        endif

        if (NIZ1(JN) > 0.d0) then
            if (JIZ1 < 1) then
                JIZ1 = 1
                ISPEND = ISPEND + 1
                RMB(ISPEND) = AIM1
                ZB(ISPEND)  = ZIM1(1)
                ISPE(ISPEND) =7
            endif
            yEXTARR(JN, 7) = NIZ1(JN)
        endif

        if (NIZ2(JN) > 0.d0) then
            if (JIZ2 < 1) then
                JIZ2 = 1
                ISPEND = ISPEND + 1
                RMB(ISPEND) = AIM2
                ZB(ISPEND)  = ZIM2(1)
                ISPE(ISPEND) = 8
            endif
            yEXTARR(JN, 8) = NIZ2(JN)
        endif

        if (NIZ3(JN) > 0.d0) then
            if (JIZ3 < 1) then
                JIZ3 = 1
                ISPEND = ISPEND + 1
                RMB(ISPEND) = AIM3
                ZB(ISPEND)  = ZIM3(1)
                ISPE(ISPEND) = 9
            endif
            yEXTARR(JN, 9) = NIZ3(JN)
        endif

        if (NI(JN) > 0.d0) then
            if (JINI < 1) JINI = 1
        endif
    enddo

!  No plasma
    if (ispe(1) < 1) then
        write(*, *) 'quit NBI: no ionised plasma,  Ne is not specified'
        IFLAG = 1
        return
    endif

! Main ion specia identification
    if (ISPE(2) < 1) then ! no ion species
        ISPEND = 2
        ISPE(ISPEND) =ISPEND
        if ((AMJ*ZMJ*JINI) > 0.d0) then
            do JN=1, NA1
                yEXTARR(JN, ISPEND) = NI(JN)
                RMB(ISPEND) = AMJ
                ZB(ISPEND)  = ZMJ
            enddo
        else
            do JN=1, NA1
                yEXTARR(JN, ISPEND) = NI(JN)
                RMB(ISPEND) = 1.d0
                ZB(ISPEND)  = 1.d0
            enddo
            write(*, *) 'NBI: no ion specie is specified. Continue with hydrogen'
            IFLAG = 0
        endif
    endif

    end subroutine nbspec

!---------------------------------------------------------------------
    subroutine nbco3(jn, jn1, JT1, JT2, n1, jbb, ieb, y1, y2, yc2, az_jr, contr, yct2, yaqbp)

    use nbstatus, only: yriplr, n_energy

    integer, intent(in) :: jn, jn1, JT1, JT2, N1, jbb, ieb
    double precision, intent(in) :: y1, y2, yc2, az_jr, contr
    double precision, dimension(n_energy), intent(inout) :: yct2, yaqbp

    integer :: jii, je, ji, jtrap, j_the
    double precision :: dt, yr1, yr2, dy, yd, ydcos, ytcos, ydys, yqbp, yf0, yfa, yddd, dy_dt

    dt = 1./(JT2 - JT1 + 1)
    YR2 = -sqrt(RJ*(RJ - RCR))

    if (RCR < RJT) then
        YR1 = -sqrt(RJT*(RJT - RCR))
        JI = N1
    else
        YR1 = 0.
        JII = (RJ - RCR)*YDYH
        JI = II(JII)
    endif

    dy = Y2 - Y1

    do JE=JBB, IEB
        YDYS = F(JE, JN) * dy
        YQBP = YAQBP(JE) * exp(-YDYS)
        YF0 = ARD(JE, JN) - az_jr
        YFA = YF0 - ARD(JE, N1)
        YDDD = (YAQBP(JE) - YQBP)

        if (ARD(JE, JI) < YF0 .and. ARD(JE, N1) > YF0) then
            JTRAP = 0
            if (YFA <= 0. .and. YFA >= YR2 .or. rcr >= yriplr(ji)) then ! lost
                YD = 0.
                YDCOS = 0.
                YTCOS = YC2 + yct2(je)
                yct2(je) = YC2
            else   ! kept
                YD = YDDD
                YDCOS = 0.d0
                YTCOS = YC2 + yct2(JE)
                YCT2(JE) = YC2
            endif
        else     ! passing
            JTRAP = 1
            if (YFA <= YR1 .and. YFA >= YR2) then ! lost
                YD = 0.
                YDCOS = 0.
                YTCOS = YC2 + YCT2(JE)
                YCT2(JE) = YC2
            else     ! kept
                YD = YDDD
                YDCOS = YD * (YC2 + yct2(JE))
                YTCOS = YC2 + YCT2(JE)
                YCT2(JE) = YC2
            endif
        endif

        YANBA(JE, JN1) = YANBA(JE, JN1) + YDDD
        YATBA(JE, JN1) = YATBA(JE, JN1) - 2.d0 * CONTR * YDDD * az_jr
        YAQBA(JE, JN1) = YAQBA(JE, JN1) + YD

        if (YTRAP(JN1) <= YTCOS) YACBA(JE, JN1) = YACBA(JE, JN1) + YDCOS

        YAQBP(JE) = YQBP
        dy_dt = YD * dt

        if (JTRAP /= 0) then   ! passing
            do j_the=JT1, JT2
                YASBA(JE, JN1, j_the) = YASBA(JE, JN1, j_the) + dy_dt
            enddo
        else   ! trapped
            do j_the=JT1, JT2
                YASBA (JE, JN1, j_the) = YASBA (JE, JN1, j_the) + 0.5d0 * dy_dt
                YASBA1(JE, JN1, j_the) = YASBA1(JE, JN1, j_the) + 0.5d0 * dy_dt
            enddo
        endif

    enddo

    end subroutine nbco3

!---------------------------------------------------------------------
    subroutine nbctr3(jn, jn1, JT1, JT2, n1, jbb, ieb, y1, y2, yc2, az_jr, contr, yct2, yaqbp)

    use nbstatus, only: yriplr, n_energy

    integer, intent(in) :: jn, jn1, JT1, JT2, N1, jbb, ieb
    double precision, intent(in) :: y1, y2, yc2, az_jr, contr
    double precision, dimension(n_energy), intent(inout) :: yct2, yaqbp

    integer :: jii, je, ji, jtrap, j_the
    double precision :: dt, yr1, yr2, dy, yd, ydcos, ytcos, ydys, yqbp, yf0, yfa, yddd, dy_dt

    dt = 1./(JT2 - JT1 + 1)
    YR2 = sqrt(RJ*(RJ - RCR))

    if (RCR < RJT) then
        YR1 = sqrt(RJT*(RJT - RCR))
        JI = N1
    else
        YR1 = 0.
        JII = (RJ - RCR)*YDYH
        JI = II(JII)
    endif

    dy = Y2 - Y1

    do JE=JBB, IEB
        YDYS = F(JE, JN) * dy
        YQBP = YAQBP(JE) * exp(-YDYS)
        YF0 = ARD(JE, JN) + az_jr
        YFA = YF0 - ARD(JE, N1)
        YDDD = (YAQBP(JE) - YQBP)

        if (ARD(JE, JI) < YF0) then  ! trapped
            JTRAP = 0
            if (YFA <= YR2 .and. YFA >= -YR2 .or. rcr >= yriplr(ji)) then ! lost
                YD = 0.
                YDCOS = 0.
                YTCOS = YC2 + yct2(je)
                yct2(je) = YC2
            else   ! kept
                YD = YDDD
                YDCOS = 0.d0
                YTCOS = YC2 + yct2(JE)
                YCT2(JE) = YC2
            endif
        else     ! passing
            if (YFA >= YR1 .and. YFA <= YR2) then ! lost
                YD = 0.
                YDCOS = 0.
                YTCOS = YC2 + YCT2(JE)
                YCT2(JE) = YC2
            else     ! kept
                YD = YDDD
                YDCOS = YD * (YC2 + yct2(JE))
                YTCOS = YC2 + YCT2(JE)
                YCT2(JE) = YC2
            endif
        endif

        YANBA(JE, JN1) = YANBA(JE, JN1) + YDDD
        YATBA(JE, JN1) = YATBA(JE, JN1) + 2.d0 * CONTR * YDDD * az_jr
        YAQBA(JE, JN1) = YAQBA(JE, JN1) + YD

        if (YTRAP(JN1) <= YTCOS) YACBA(JE, JN1) = YACBA(JE, JN1) + YDCOS

        YAQBP(JE) = YQBP
        dy_dt = YD * dt

        if (JTRAP /= 0) then   ! passing (contr)
            do j_the=JT1, JT2
                YASBA1(JE, JN1, j_the) = YASBA1(JE, JN1, j_the) + dy_dt
            enddo
        else   ! trapped
            do j_the=JT1, JT2
                YASBA (JE, JN1, j_the) = YASBA (JE, JN1, j_the) + 0.5d0 * dy_dt
                YASBA1(JE, JN1, j_the) = YASBA1(JE, JN1, j_the) + 0.5d0 * dy_dt
            enddo
        endif

    enddo

    end subroutine nbctr3

!---------------------------------------------------------------------
    subroutine nbco2gc(jn, jn1, N1, jbb, ieb, y1, y2, yc2, ydex, &
        az_jr, contr, yct2, yaqbp, ydys)

    use nbstatus, only: yriplr, n_energy

    integer, intent(in) :: jn, jn1, N1, jbb, ieb
    double precision, intent(in) :: y1, y2, yc2, ydex, az_jr, contr
    double precision, intent(out) :: ydys
    double precision, dimension(n_energy), intent(inout) :: yct2, yaqbp

    integer :: jii, je, ji, jtrap, j_the, ntet, jloss, j_jn
    double precision :: yr2, yr1, dy, yd, ydcos, ytcos, yqbp, yf0, yfa, yddd, &
        YJSERF

    ntet = ntet1 - 1
    YR2 = -SQRT(RJ*(RJ - RCR))
    if (RCR < RJT) then
        YR1 = -SQRT(RJT*(RJT - RCR))
        JI  = N1
    else
        YR1 = 0.
        JII = (RJ - RCR)*YDYH
        JI  = II(JII)
    endif

    dy = Y2 - Y1
    if (dy <= 0.D0) then
        WRITE(*, *) 'NBTRAG: warning: overlapping of surfaces'
        dy = -dy
    endif

    do JE=JBB, IEB
        YDYS  = F(JE, JN)*dy
        YQBP = YAQBP(JE)*YDEX
        YF0  = (ARD(JE, JN) - az_jr)
        YFA  = YF0 - ARD(JE, N1)
        YDDD = (YAQBP(JE) - YQBP)

        if (ARD(JE, JI) <= YF0 .and. ARD(JE, N1) >= YF0) then
            JTRAP = 0
            if (YFA <= 0. .and. YFA >= YR2 .or. RCR >= YRIPLR(JI)) then
                YD    = 0.
                YDCOS = 0.
                YTCOS = YC2 + YCT2(JE)
                YCT2(JE) = YC2
                JNL   = JN1
                JNR   = JN1
                JLOSS = 1
            else
                CALL NBORBCO(JE, JN, RCR, RJ, RJT, YC2, YF0, JNR, JNL, JTRAP, JLOSS)
                if (JLOSS == 0) then
                    YJSERF = 1./(JNR - JNL + 1.D0)
                    YD     = YDDD*YJSERF
                    YDCOS  = YD*2.D0
                    YTCOS  = YC2 + YCT2(JE)
                    YCT2(JE) = YC2
                else
                    YD    = 0.
                    YDCOS = 0.
                    YTCOS = YC2 + YCT2(JE)
                    YCT2(JE) = YC2
                    JNL = JN1
                    JNR = JN1
                endif
            endif
        else
            JTRAP = 1
            if (YFA <= YR1 .and. YFA >= YR2) then
                YD    = 0.
                YDCOS = 0.
                YTCOS = YC2 + YCT2(JE)
                YCT2(JE) = YC2
                JNL   = JN1
                JNR   = JN1
                JLOSS = 0
            else
                CALL NBORBCO(JE, JN, RCR, RJ, RJT, YC2, YF0, JNR, JNL, JTRAP, JLOSS)
                if (JLOSS == 0) then
                    YJSERF = 1./(JNR - JNL + 1)
                    YD     = YDDD*YJSERF
                    YDCOS  = YD*2.D0
                    YTCOS  = YC2 + YCT2(JE)
                    YCT2(JE) = YC2
                else
                    YD    = 0.
                    YDCOS = 0.
                    YTCOS = YC2 + YCT2(JE)
                    YCT2(JE) = YC2
                    JNL = JN1
                    JNR = JN1
                endif
            endif
        endif

        YDDD = -2.D0*CONTR*YDDD*az_jr*YJSERF

        do j_jn=JNL, JNR
            YATBA(JE, j_jn) = YATBA(JE, j_jn) + YDDD
            YANBA(JE, j_jn) = YANBA(JE, j_jn) + YD
            YANBA1(JE, JN1) = YANBA1(JE, JN1) + YD
            YAQBA(JE, j_jn) = YAQBA(JE, j_jn) + YD
            if (YTRAP(j_jn) <= YCOS(j_jn)) then
                YACBA(JE, j_jn) = YACBA(JE, j_jn) + YDCOS*YCOS(j_jn)
            endif
            YAQBP(JE) = YQBP
            if (YCOS(j_jn) > 0.D0) then
                j_the = ntet*YCOS(j_jn)*YBTDB(JN)
                if (j_the > ntet) j_the = ntet
                if (j_the < 1)    j_the = 1
                YASBA(JE, j_jn, j_the) = YASBA(JE, j_jn, j_the) + YD
            else
                j_the = (-ntet*YCOS(j_jn)*YBTDB(JN))
                if (j_the > ntet) j_the = ntet
                if (j_the < 1)    j_the = 1
                YASBA1(JE, j_jn, j_the) = YASBA1(JE, j_jn, j_the) + YD
            endif
        enddo

    enddo

    end subroutine nbco2gc

!---------------------------------------------------------------------
    subroutine nborbco(JE, JN, RCR, RJ, RJT, YC2, YF0, JNR, JNL, ITRAP, ILOSS)

!----- First orbit analysis taking account ----
! the input source redistribution
! due to the first orbit deviation (guide centre apprx.)
! co_injection at birth
! input:
!    JE, JN, YC2, YF0	starting point of
!    velocity index, surface index, YC2=<v.B>/B
!    RCR, RJ, RJT     = critical, external, internal boundary radii
!    YF0            = F value at the birth position
!    REJ(j)/RIJ(j)  = major radius of (j) surface with mid plane LFS/HFS
!    RC(j   )       = major radius of mag. axis
! output:
!    JNL, JNR   = minimum/maximum surface index of the orbit with finit RLM
!    ILOSS      = 0/1 if particle is kept/lost
!    ITRAP      = 0/1 for banana/passing orbits
!    YCOS(JJN)  = <v.B>/vB (x(JJN))

    integer, intent(in) :: je, jn
    double precision, intent(in) :: RCR, RJ, RJT, YC2, YF0
    integer, intent(out) :: ITRAP, ILOSS, jnr, jnl

    integer :: J, J1
    double precision :: Y, YY, YY1

    JNR  = JN
    JNL  = JN

    ITRAP = 1 !(def.passing)
    ILOSS = 0 !(def. kept)
    YCOS(JN) = YC2
    J = JN
    if (JN > n1) then
        write(*, *) 'error in nbrbco: JN > N1: ', JN, '>', N1
        return
    endif

    J1 = 1 ! sign of index increment

1    continue

    if (J > n1) goto 3   ! lost

    Y = ARD(JE, J) - YF0  ! motion to LFS (R->RJ)
    yy = RCR/2.d0  !*YBTDB(J)
    yy = yy + sqrt(yy**2 + y**2) ! R=Rc/2+sqrt((Rc/2)^2+y^2)
    if (j == jn) yy1 = yy
    if (yy <= REJ(j)) then !Ri(Z=0)<= R<=Re(Z=0) real root
        if (yy < RJT .or. yy > RJ) goto 3 !lost
        ycos(j) = y/yy
        ycos(j) = max(-1.d0, ycos(j))
        ycos(j) = min( 1.d0, ycos(j))
        JNR  = J
        J = J + J1 ! motion to LFS (R->RJ)
        goto 1
    endif ! exit: end of cycle to LFS: fake root

    J = JN
    J1 = -1

2   continue

    if (J < 1) then ! reaches plasma ceter
        ILOSS = 0 ! kept in plasma
        return
    endif
    Y = ARD(JE, J) - YF0  ! motion to HFS (R->RJT)
    yy = RCR/2.d0  !*YBTDB(J)
    yy = yy + sqrt(yy**2 + y**2) ! R=Rc/2+sqrt((Rc/2)^2+y^2)
    if (yy <= REJ(j) .and. yy >= RIJ(j)) then ! real root
        if (y < 0.) itrap = 0 !trapped to well
        if (yy < RJT .or. yy > RJ) goto 3  !lost
        ycos(j) = y/yy
        ycos(j) = max(-1.d0, ycos(j))
        ycos(j) = min( 1.d0, ycos(j))
        JNL  = J
        J =J + J1 ! motion to HFS (R->RJT)
        yy1 = yy
        goto 2
    endif ! exit: end of cycle to HFS: fake root
    return

3   continue
    ILOSS = 1 ! lost from plasma

    end subroutine nborbco

end module nbicom
