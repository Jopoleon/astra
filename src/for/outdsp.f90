!----------------------------------------------------------------------|
subroutine OUTDSP(MARK, JIFNEW, IYO, ITIMES, TTOUT, TOUT)
!----------------------------------------------------------------------|
! Drawing options:
!   X - axis
!      0 <= rho <= ROC=RHO(NA1)
!      0 <= a <= ABC
!      0 <= a <= AB
!
! MARK = 1     Put marks
! MARK = 0     Solid lines
! MARK =-1     Dashed lines
! JIFNEW  =  0 Re-draw (erase) the previous curves
! JIFNEW =/= 0 New curves only
! JIFNEW < 0   Don't mark resonances q=m/n
! JIFNEW > 10  Call from Review. (JIFNEW-10) is used to control erasing
!----------------------------------------------------------------------|

use parameter_inc, only: NRD, NRDX, NRW
use status_inc, only: AMETR, MU, SHIF, ELON, TRIA
use const_inc, only: XOUT, NAB, NA1, ABC, TINIT, TSCALE, RTOR, &
    MEQUIL, LEQ, UPDWN, TIME
use outcmn_inc, only: frame_wid, frame_hei, canv_hei, canv_wid, nx_canvas, ny_canvas, &
    curves_per_frame, active_tab, MOD10, NWIND1, NWIND3, NWINDX, &
    IFDFAX, IY0, IYM, KPRI, DXLET, DYLET, NPTM, ICVMX, &
    NROUT, ROUT, OSHIFR, NAMER, SCALER, &
    NTOUT, TOUT, OSHIFT, NAMET, SCALET, &
    NXOUT, NAMEX, NARRX, EXARNM, DATAX, TOUTX, LTOUT, &
    XAXES, GRAL, GRAP, pixel_ymid, meter2pixel, LineWidth, &
    Black, WarningColor, EraseColor, &
    equ_file, NBFILE, NBFLAG
use timeoutput_inc, only: NTIMES
use expdat, only: raw_profile_map, DATARR, BNDR, BNDZ
use ac_neg1, only: NUM, NKL1, NKL2, JMIN, JMAX, MODK
use dbl2char, only: fmt_xf, fmt4
use char_manip, only: len_trim_tab, str_in_list
use debugger, only: markloc, debug, astra_stop

implicit none

integer, parameter :: jzero=0, fshift=10
integer, intent(in) :: MARK, JIFNEW, ITIMES
integer, intent(inout) :: IYO(ITIMES,*)
double precision, intent(in) :: TTOUT(ITIMES)
double precision, intent(inout) :: TOUT(ITIMES, NRW)

integer :: PTM(2), PTMO(2, NRDX, NRW), &
    IWN(16), EQOLD(260, 11), &
    IST, IQ1, IQ2, text_posy, jt_old, JS, MODEX, &
    IYM0, LTOUT1, LTOUT2, JFNEW, STYL, x_shift, y_shift, jx_canv, jy_canv, JY, jxout, &
    JW, j_curve, j_canv, &
    IYMN, IYMX, JDSP, test_posx, NPTMO(NRW), jlx(8), &
    NP1, i, j, half_wid, &
    j1, jj, jsco, jn, jpnt, jsc, jposy, jarr, jtyp, n_canvas, &
    jplot_in_canv, jcol, jcol2, jprof, jtrace
double precision :: SC(NRW), YX, r_out, YA, YL, YR, YQ1, YQ2, &
     XROUT, YZ, ABSC, ymin, ymax, px_rmag
double precision, dimension(NRD) :: xplot, yplot, xtrace, ytrace, xtrace_old
double precision, dimension(NRD, ICVMX) :: xold, yold, ytrace_old
character(len=80) :: STRI
character(len=5 ) :: XF4
character(len=6 ) :: CHAR6

save EQOLD, PTMO, NPTMO, YQ1, YQ2, IQ1, IQ2, IWN, xold, yold, xtrace_old, ytrace_old

!----------------------------------------------------------------------|
call markloc('OUTDSP')

half_wid = frame_wid/2
YA = 0.
JFNEW = JIFNEW
if (JFNEW >= 10) JFNEW = JFNEW - 10
if (JFNEW /= 0) then
    do J=1, 16
        IWN(J) = 0
    enddo
endif

! Different options for abscissa. 
! Implemented for modes 0, 1, 2, 3
MODEX = XOUT + 0.49
SELECT CASE(MODEX)
    CASE(0)        ! Draw up to AB against "a"
        NP1 = NAB
    CASE(1:3)      ! Draw up to ABC against "a", "rho", "Psi"
        NP1 = NA1
    CASE DEFAULT ! Unknown option
        return
END SELECT

if (MOD10 == 3) then
    NP1 = NA1
endif

if (MARK ==  0) IST = 1
if (MARK ==  1) IST = 8
if (MARK == -1) IST = 15
IYMN = frame_hei - IYM
IYM0 = IYM - canv_hei
IYMX = frame_hei - IY0
JY = 10*frame_hei

ymin = dble(frame_hei - IYM)
ymax = dble(frame_hei - IY0)

!----------------------------------------------------------------------|
n_canvas = nx_canvas*ny_canvas

SELECT CASE(MOD10)

!------------
CASE(1: 3)  ! Profiles

    call markloc('Drawing mode 1-3', debug_lev=debug*2)
    do JPROF=1, NROUT
        do J=1, NP1
            ROUT(j, jprof) = ROUT(j, jprof) + OSHIFR(jprof)
        enddo
    enddo
    call SCAL(NROUT, SC, SCALER, ROUT(4, 1), NP1 - 3, NRD)

! Plot curves

    j_curve = 0
    plot_prof: do jprof=1, NROUT
        JW = NWIND1(jprof) - curves_per_frame(MOD10)*active_tab(MOD10) ! 1-16 for mode '1'
        if (NAMER(jprof) == '    ') JW = 0
        if (JW <= 0 .or. JW > curves_per_frame(MOD10)) CYCLE plot_prof
        j_canv = MOD(JW - 1, n_canvas) + 1        ! 1-8 for mode '1'
        jx_canv = MOD(j_canv - 1, nx_canvas)      ! 0-3 for mode '1'
        jy_canv = (j_canv - 1)/nx_canvas          ! 0-1
        jplot_in_canv = (JW - 1)/n_canvas
        jcol = jplot_in_canv + 2
        y_shift = canv_hei * (ny_canvas - 1 - jy_canv)
        x_shift = canv_wid*jx_canv

        YL = max(0.d0, ABSC(GRAL(jprof)))
        YR = ABSC(GRAP(jprof))
        if (YL >= YR) YL = 0.d0

! Translate curve into pixel
        jxout = 0
        do j=1, NP1
            YX = ABSC(AMETR(j)) 
            if (YX >= YL .and. YX <= YR) then
                jxout = jxout + 1
                if (jxout == 1 .and. j > 1) then ! left edge interpolation
                    YA = ROUT(J, jprof) + (ROUT(J-1, jprof) - ROUT(J, jprof))*(YL - YX)/(YA - YX)
                    r_out = min(max(YA/SC(jprof), -7.d0), 7.d0)
                    JDSP  = 10*(canv_hei*r_out + IYMN + y_shift)
                    xplot(jxout) = dble(x_shift)
                    yplot(jxout) = frame_hei - min(max(dble(canv_hei)*r_out + ymin + dble(y_shift), ymin), ymax)
                    jxout = jxout + 1
                endif
                r_out = min(max(ROUT(J, jprof)/SC(jprof), -7.d0), 7.d0)
                JDSP  = 10*(canv_hei*r_out + IYMN + y_shift)
                xplot(jxout) = dble(x_shift) + dble(frame_wid)/dble(nx_canvas)*(YX - YL)/(YR - YL)
                yplot(jxout) = frame_hei - min(max(dble(canv_hei)*r_out + ymin + dble(y_shift), ymin), ymax)
            endif
            if (YA <= YR .and. YX > YR) then ! right edge interpolation
                jxout = jxout + 1
                YA = ROUT(J, jprof) + (ROUT(J-1, jprof) - ROUT(J, jprof))*(YR - YX)/(YA - YX)
                r_out = min(max(YA/SC(jprof), -7.d0), 7.d0)
                JDSP  = 10*(canv_hei*r_out + IYMN + y_shift)
                xplot(jxout) = dble(x_shift) + dble(frame_wid)/dble(nx_canvas)
                yplot(jxout) = dble(frame_hei) - min(max(dble(canv_hei)*r_out + ymin + dble(y_shift), ymin), ymax)
            endif
            YA = YX
        enddo
        if (KPRI >= 1 .and. KPRI <= 2) then
            write(STRI, '(1A6, 1A4, 1A1)') 'Plot "', NAMER(jprof), '"'
            j = len_trim_tab(STRI)
            call pscom(STRI, j)
        endif
        STYL = (jcol - 1)*MARK
        j_curve = j_curve + 1
        if (j_curve <= ICVMX) then
            call update_curve(jxout, IWN(JW), jcol, STYL, xold(1:, j_curve), yold(1:, j_curve), xplot, yplot)
            xold(1: jxout, j_curve) = xplot(1: jxout)
            yold(1: jxout, j_curve) = yplot(1: jxout)
        endif
        IWN(JW) = jxout
        do J=1, NP1
            ROUT(j, jprof) = ROUT(j, jprof) - OSHIFR(jprof)
        enddo
! Variable labels
        call colovm(jcol)
        test_posx = x_shift + jplot_in_canv*canv_wid*n_canvas/curves_per_frame(mod10)
        text_posy = DYLET + FSHIFT + (IYM - IY0 + 2 + DYLET)*jy_canv + 2
        call CMARK(text_posy, test_posx, SC(jprof), OSHIFR(jprof), NAMER(jprof), STYL)

    enddo plot_prof

! Plot dots
    do j=1, n_canvas
        jlx(j) = 0
    enddo
    jsco = 0
    jn = 0

    plot_profx: do jxout=1, NXOUT
        if (NWINDX(jxout) == 0) CYCLE plot_profx ! NWINDX set in ininam.tmp
        CHAR6 = NAMEX(jxout)
        if (CHAR6(1: 1) == ' ') CYCLE plot_profx
        do jprof=1, NARRX
            if (EXARNM(jprof) == CHAR6) jn = jprof 
        enddo
        if (jn == 0) then
            write(*, *)
            write(*, *) '>>> Error in the model "', TRIM(equ_file), '"'
            call astra_stop('>>> Unknown data set "' // CHAR6 // '" is requested')
        endif

        jpnt = NPTM(jn)
        if (jpnt <= 0) CYCLE plot_profx

        jsc = NWINDX(jxout)
        if (abs(SC(jsc)) < 1.1E-7) call SCAL(1, SC(jsc), SCALER(jsc), DATAX(1, jn), jpnt, NRDX)
        JW = NWIND1(jsc) - curves_per_frame(MOD10)*active_tab(MOD10)
        if (JW <= 0 .or. JW > curves_per_frame(MOD10)) CYCLE plot_profx
        j_canv = MOD(JW - 1, n_canvas) + 1        ! 1-8 for mode '1'
        jx_canv = MOD(j_canv - 1, nx_canvas)      ! 0-3 for mode '1'
        jy_canv = (j_canv - 1)/nx_canvas          ! 0-1
        jcol = (JW - 1)/n_canvas + 2
        y_shift = canv_hei * (ny_canvas - 1 - jy_canv)
        x_shift = canv_wid*jx_canv

        call colovm(jcol)
        if (jsc /= jsco) then
            jcol2 = 5  ! 13
        else
            jcol2 = jcol2 + 1
        endif

        if (JFNEW == 0) then
            call colovm(EraseColor)
            do j=1, NPTMO(jxout)
                call NMARK(PTMO(1, j, jxout), jcol2)
            enddo
        endif

        if (KPRI >= 1 .and. KPRI <= 2) then
            write(STRI, '(1A6, 1A6, 1A1)')'Dots "', EXARNM(jn), '"'
            j = len_trim_tab(STRI)
            call pscom(STRI, j)
        endif
        YL = max(0.d0, ABSC(GRAL(jsc)))
        YR = ABSC(GRAP(jsc)) 
        if (YL >= YR) YL = 0.d0
        j1 = 0
        do j=1, jpnt
            YX = XAXES(j, jn)
            if (MODEX >= 1 .and. YX > ABC) CYCLE
            YA = ABSC(YX)
            if (YA > 1. .or. YA < YL .or. YA > YR) CYCLE
            j1 = j1 + 1
            PTM(1) = x_shift + frame_wid/nx_canvas*(YA-YL)/(YR-YL)
            r_out= max((DATAX(j, jn) + OSHIFR(jsc))/SC(jsc), -7.d0)
            r_out= min(r_out, 7.d0)
            JDSP = canv_hei*r_out + IYMN + y_shift
            PTM(2) = frame_hei - min(max(JDSP, IYMN), IYMX)
            call colovm(jcol, 2)
            call NMARK(PTM, jcol2)
            PTMO(1, j1, jxout) = PTM(1)
            PTMO(2, j1, jxout) = PTM(2)
        enddo
        NPTMO(jxout) = j1

        if (KPRI >= 1 .and. KPRI <= 2) then
            write(STRI, '(1A)')'Mark & time for dots'
            j = len_trim_tab(STRI)
            call pscom(STRI, j)
        endif
        jlx(j_canv) = jlx(j_canv) + 1
        test_posx = x_shift + canv_wid - 6*DXLET
        text_posy = (1 + jlx(j_canv))*DYLET + FSHIFT + (IYM - IY0 - canv_hei)*(jy_canv) + 3
        XF4 = fmt_xf(TOUTX(jn), 4)
        call textvm(test_posx, text_posy, XF4, 5)
        PTM(1) = test_posx + 5.5*DXLET
        PTM(2) = text_posy - DYLET/2 + 1
        jcol = (JW - 1)/n_canvas + 2
        call NMARK(PTM, jcol2)

        jsco = jsc
        if (jcol2 == 7) jcol2 = 1
    enddo plot_profx

    if (JIFNEW >= 10) return
   
! Erase/put resonance radii
    call colovm(EraseColor)
    do J=1, frame_wid, frame_wid/nx_canvas
        if (IQ1 /= 0) then
            test_posx = J + frame_wid/nx_canvas*YQ1 - 1
            call drawvm(0, test_posx, IYM0, test_posx, IYM0 - 4)
        endif
        if(IQ2 == 0) CYCLE
        test_posx = J + frame_wid/nx_canvas*YQ2 - 1
        call drawvm(0, test_posx, IYM0, test_posx, IYM0 - 4)
        test_posx = test_posx + 2
        call drawvm(0, test_posx, IYM0, test_posx, IYM0 - 4)
    enddo

    call colovm(WarningColor)
    IQ1 = 0
    IQ2 = 0
    do J=1, NAB-1
        if (MU(J) > 1. .and. MU(J+1) <= 1.) IQ1 = J
        if (MU(J) > .5 .and. MU(J+1) <= .5) IQ2 = J
    enddo
    if (KPRI >= 1 .and. KPRI <= 2) then
        write(STRI, '(1A)') 'Mark resonance radii'
        j = len_trim_tab(STRI)
        call pscom(STRI, j)
    endif
    if (IQ1 /= 0) YQ1 = ABSC(AMETR(IQ1)) 
    if (IQ2 /= 0) YQ2 = ABSC(AMETR(IQ2)) 

    do J=1, frame_wid, frame_wid/nx_canvas
        if (IQ1 /= 0) then
            test_posx = J + frame_wid/nx_canvas*YQ1 - 1
            call drawvm(0, test_posx, IYM0, test_posx, IYM0 - 4)
        endif
        if (IQ2 == 0) CYCLE
        test_posx = J + frame_wid/nx_canvas*YQ2 - 1
        call drawvm(0, test_posx, IYM0, test_posx, IYM0 - 4)
        test_posx = test_posx + 2
        call drawvm(0, test_posx, IYM0, test_posx, IYM0 - 4)
    enddo

!----------------
CASE(4:5)
    call markloc('Drawing mode 4/5', debug_lev=2*debug)

!----------------
CASE(6)  ! Time traces

    call markloc('Drawing mode 6', debug_lev=2*debug)
    LTOUT1 = 1
    LTOUT2 = LTOUT

! right_label_position=JDX*JDMX=23*5*5=575 (see typdsp.f)
    do J=1, LTOUT-1
        r_out = (TTOUT(J) - TINIT)*575/abs(TSCALE)
        IYO(J, ICVMX+1) = 6*DXLET + r_out
        xtrace(J) = 6*DXLET + r_out
        if (r_out < 0)  LTOUT1 = J + 1
        if (IYO(J, ICVMX+1) <= frame_wid - 1) LTOUT2 = J - 1
    enddo
    LTOUT2 = LTOUT2 - LTOUT1 + 1
    if (LTOUT2 < 3) then ! No plot for small time
        return
    endif
    do jj=1, min(NTOUT, NRW)
        do j=LTOUT1, LTOUT1 + LTOUT2
            TOUT(j, jj) = TOUT(j, jj) + OSHIFT(jj)
        enddo
    enddo
    call SCAL(NTOUT, SC, SCALET, TOUT(LTOUT1, 1), LTOUT2, ITIMES)

    j_curve = 0
    plot_traces: do jtrace=1, min(NRW, NTOUT)
        JW = NWIND3(jtrace) - curves_per_frame(MOD10)*active_tab(MOD10)
        if (NAMET(jtrace) == '    ') CYCLE plot_traces
        if (JW <= 0 .or. JW > curves_per_frame(MOD10)) CYCLE plot_traces
        jplot_in_canv = (JW - 1)/n_canvas       ! <-> color
        j_canv = MOD(JW - 1, n_canvas) + 1      ! 1-8 for mode '1'
        do J=1, LTOUT
            r_out = max(TOUT(J, jtrace)/SC(jtrace), -7.d0)
            r_out = min(r_out, 7.d0)
            JDSP  = 10*(canv_hei*r_out + IYMN + (n_canvas - j_canv)*canv_hei)
            JDSP  = max(JDSP, 10*IYMN)
            IYO(J, ICVMX+2) = JY - min(JDSP, 10*IYMX)
            ytrace(J) = dble(frame_hei) - min(max(canv_hei*r_out + ymin + (n_canvas - j_canv)*canv_hei, ymin), ymax)
        enddo

        if (JFNEW == 0) then
            jt_old = LTOUT2
        else
            jt_old = 0
        endif
        if (KPRI >= 1 .and. KPRI <= 2) then
            write(STRI, '(1A6, 1A4, 1A1)')'Plot "', NAMET(jj), '"'
            j = len_trim_tab(STRI)
            call pscom(STRI, j)
        endif
        jcol = jplot_in_canv + 1 
        STYL = jcol*MARK
        jcol  = jcol  + 1
        j_curve = j_curve + 1
        if (j_curve <= ICVMX) then
            call update_curve(LTOUT2 + 1, jt_old, jcol, STYL, xtrace_old(1:jt_old), ytrace_old(1: jt_old, j_curve), &
                xtrace(1), ytrace(1) )
            xtrace_old(1: jt_old+1) = xtrace(1: jt_old+1)
            ytrace_old(1: jt_old+1, j_curve) = ytrace(1: jt_old+1)
        endif
        test_posx = 0
        text_posy = FSHIFT + DYLET*(2*jplot_in_canv + 1) + (j_canv - 1)*canv_hei
        call CMARKT(text_posy, test_posx, SC(jtrace), OSHIFT(jtrace), NAMET(jtrace), STYL)

        do j=LTOUT1, LTOUT1 + LTOUT2
            TOUT(j, jtrace) = TOUT(j, jtrace) - OSHIFT(jtrace)
        enddo
    enddo plot_traces

!----------------
CASE(7)
    call markloc('Drawing mode 7', debug_lev=2*debug)

CASE(8)

    call markloc('Drawing mode 8', debug_lev=2*debug)
    px_rmag = (RTOR + SHIF(1))*meter2pixel
    call colovm(Black)
    call drawcurve(0, 2, LineWidth, (/0., px_rmag/), (/pixel_ymid, pixel_ymid/))
! Plot the complete wall structure (Pixmap # 1)
    call plot_wall

    if (NBFLAG /= 0 .and. NBFILE(1:1) /= '*' .and. jifnew /= 0) then ! Create/update NBI Pixmap # 2
        call drawfoot(jifnew)
        call redraw(2)
        NBFLAG = 0
    endif

! if (data file includes NAMEXP BND) then (NBND > 0);
! or (NBND == 8) after calling equil with no boundary points provided;
!     NBND == 0 otherwise

    SELECT CASE(LEQ(5))
    CASE(:1)
        call DRAW3M(jifnew, DYLET, EQOLD, &
            NA1, NAB, RTOR, AMETR, SHIF, UPDWN, ELON, TRIA)
    CASE(3)
        call bnd_draw(JIFNEW, IYO, TIME)
        call DRAWSPFLUX
    CASE(4: 5)
        if (MEQUIL == 0) then
            SHIF = 0.0
            ELON = 1.0
            TRIA = 0.0 
            call DRAW3M(jifnew, DYLET, EQOLD, &
                NA1, NAB, RTOR, AMETR, SHIF, UPDWN, ELON, TRIA)
        else
            call bnd_draw(JIFNEW, IYO, TIME)
            call DRAWSPFLUX
        endif
    END SELECT

! Plot dots of R, z if an exp profile is input as function of R, z (GRIDTYPE 18-20)
    jcol = -1
    loop8: do jxout=1, NXOUT
        if (NWINDX(jxout) == 0) CYCLE loop8
        CHAR6 = NAMEX(jxout)
        if (CHAR6(1:1) == ' ') CYCLE loop8
        if (jxout > 1) then
            do jj=jxout-1, 1, -1
                if (NAMEX(jj) == CHAR6) CYCLE loop8 ! double name?
            enddo
        endif

        jn = str_in_list(CHAR6, EXARNM)
    
        jarr = IFDFAX(jn)
        if (jarr <= 0)  CYCLE loop8
        jtyp = raw_profile_map%grid_type(jarr)
        if (jtyp < 18)  CYCLE loop8
        jpnt = raw_profile_map%nrho(jarr)
        if (jpnt <= 0)  CYCLE loop8
        jcol = jcol+1
        js = raw_profile_map%jbeg_grid(jarr)

        if (JFNEW == 0) then
            call colovm(EraseColor)
            do j=1, NPTMO(jxout)
                call NMARK(PTMO(1, j, jxout), 7)
            enddo
        endif
        call colovm(jcol+2)
        do j=1, jpnt
            if (jtyp == 18)  then
                YR = DATARR(js)
                YZ = DATARR(js+j)
            elseif (jtyp == 19)  then
                YR = DATARR(js+j)
                YZ = DATARR(js)
            elseif (jtyp == 20) then
                YR = DATARR(js-1+j)
                YZ = DATARR(jpnt+js-1+j)
            else
                write(*, *) 'Unknown input-grid type'
            endif
            PTM(1) = YR*meter2pixel
            PTM(2) = pixel_ymid - YZ*meter2pixel
            call NMARK(PTM, 7)
            PTMO(1, j, jxout) = PTM(1)
            PTMO(2, j, jxout) = PTM(2)
        enddo
! git hardcoded 510
        call textvm(510, 100 + 2*DYLET*jcol, NAMEX(jxout), 6)
        NPTMO(jxout) = jpnt
    enddo loop8

!----------------
CASE(9)
    call markloc('User drawing mode', debug_lev=2*debug)

END SELECT

return
end subroutine outdsp

!======================================================================|
subroutine bnd_draw(ifnew, IYO, time_in)
!----------------------------------------------------------------------|
! IFNEW  =  0 Re-draw (erase) the previous curves
! IFNEW =/= 0 New curves only
! IFNEW < 0 Don't mark resonances q=m/n
! IFNEW  > 10 Call from Review. (JIFNEW-10) is used to control erasing

use outcmn_inc, only: Red, NBNT, pixel_ymid, meter2pixel
use expdat, only: BNDTIM, BNDR, BNDZ
use const_inc, only: NBND

implicit none

integer, intent(in) :: ifnew
integer, intent(inout) :: IYO(2, *)
double precision, intent(in) :: time_in

integer :: j, j1, j2, jj, PTM(2)
double precision :: YS, YX, YXL, YXR, YZ

j2 = 1 + NBND/32

if (IFNEW == 0) then
    j = 0   ! Erase the old points
! j = 31   ! Draw the old points with shadow color
    call colovm(j)
    do j=1, NBND, j2
        PTM(1) = IYO(1, j)
        PTM(2) = IYO(2, j)
        call NMARK(PTM, 4)
    enddo
endif

! Boundary points
call colovm(Red)

if (NBNT <= 1) then
    jj = 1
    do j=1, NBND, j2
        j1 = max(1, NBNT + (j - 1)*jj)
        PTM(1) = BNDR(j1)*meter2pixel
        PTM(2) = pixel_ymid - BNDZ(j1)*meter2pixel
        call NMARK(PTM, 4)   !Use (PTM, 4) for *
        IYO(1, j) = PTM(1)
        IYO(2, j) = PTM(2)
    enddo
else if (time_in <= BNDTIM(1) .or. time_in >= BNDTIM(NBNT)) then ! extrapolate flat
    jj = NBNT
    do j=1, NBND, j2
        j1 = NBNT + (j - 1)*jj
        PTM(1) = BNDR(j1)*meter2pixel
        PTM(2) = pixel_ymid - BNDZ(j1+jj)*meter2pixel
        call NMARK(PTM, 4)   !Use (PTM, 4) for *
        IYO(1, j) = PTM(1)
        IYO(2, j) = PTM(2)
    enddo
else               ! interpolate linearly
    do j=1, NBNT   ! Find current time
        if (time_in > BNDTIM(j)) jj = j
    enddo
    YS  = BNDTIM(jj+1) - BNDTIM(jj)
    YXL = (time_in - BNDTIM(jj  ))/YS
    YXR = (time_in - BNDTIM(jj+1))/YS
    do j=1, NBND, j2 ! Time differentiation
        j1 = jj + (j - 1)*NBNT
        YX = YXL*BNDR(j1+1) - YXR*BNDR(j1)
        YZ = YXL*BNDZ(j1+1) - YXR*BNDZ(j1)
        PTM(1) = YX*meter2pixel
        PTM(2) = pixel_ymid - YZ*meter2pixel
        call NMARK(PTM, 4)
        IYO(1, j) = PTM(1)
        IYO(2, j) = PTM(2)
    enddo
endif

return
end subroutine bnd_draw

!======================================================================|
subroutine DRAW3M(jifnew, DYLET, EQOLD, &
    NA1, NAB, RTOR, AMETR, SHIF, UPDWN, ELON, TRIA)
!----------------------------------------------------------------------|

use outcmn_inc, only: Magenta, Pink, EraseColor, Red, frame_hei, &
    LineWidth, pixel_ymid, meter2pixel

use const_inc, only: GP

implicit none

integer, intent(in) :: jifnew, DYLET, NA1, NAB
integer, intent(inout) :: EQOLD(260, 11)
double precision, intent(in) :: RTOR, AMETR(*), ELON(*), TRIA(*), UPDWN, SHIF(*)

integer, parameter :: n_theta=64
integer :: j, jj, jn, JX
integer, dimension(260) :: plot_arr
double precision :: YR, YZ, YFI

!----------------------------------------------------------------------|
! Redraw magnetic surfaces:
do J=1, 10
! Erase:
    if (jifnew == 0) then
        call colovm(EraseColor)
        do JJ=1, 130
            plot_arr(2*JJ-1) = EQOLD(2*JJ-1, J)
            plot_arr(2*JJ  ) = EQOLD(2*JJ, J)
        enddo
        jj = 0
        jn = 130
        call d2polyline(jj, jn, LineWidth, plot_arr)
    endif

! New configuration:
    JX = max(1., 0.1*NA1*J - 1)
    JX = min(NA1, JX)
    if (J == 10) then
        call colovm(Red)
        JX = NA1
        if (NA1 /= NAB) then
            call drawvm(0, frame_hei, 55, 400, 55)
            call textvm(410, 55+DYLET/2, 'Transport boundary', 18)
        endif
    else
        call colovm(Magenta)  ! Outermost flux surface
    endif
    do JJ=1, 130
        YFI = GP*(JJ - 1)/64.
        YZ = UPDWN + AMETR(JX)*ELON(JX)*sin(YFI)
        YR = RTOR + SHIF(JX) + AMETR(JX)*(cos(YFI) + 0.5*TRIA(JX)*(cos(2.*YFI) - 1.))
        plot_arr(2*JJ-1) = 10.*YR*meter2pixel
        plot_arr(2*JJ) = 10.*(pixel_ymid - YZ*meter2pixel)
! Save picture:
        EQOLD(2*JJ-1, J) = plot_arr(2*JJ-1)
        EQOLD(2*JJ  , J) = plot_arr(2*JJ)
    enddo
    jj = 0
    jn = 130
    call d2polyline(jj, jn, LineWidth, plot_arr)
enddo

return
end subroutine DRAW3M

!======================================================================|
subroutine plot_wall
!----------------------------------------------------------------------|

use const_inc, only: AB, ELONM, RTOR, TRICH, GP2
use outcmn_inc, only: wall_gc_file, Blue, White, LineWidth, pixel_ymid, &
    meter2pixel
use debugger, only: debug

implicit none

integer, parameter :: n_theta=64, ngc_max=750
integer :: j, j1, jgc, ios, nSHOT, Ndim, NGC, plot_arr(100), ndim_gc
integer, dimension(40) :: ixbeg, lenix, valix
double precision :: pol_ang, Rwall, Zwall
double precision, dimension(n_theta) :: xwall, ywall
double precision, dimension(ngc_max) :: xGC, yGC
double precision, dimension(ngc_max, 2) :: xyGC
character(len=64) :: STRI

call colovm(Blue)
open(7, FILE=TRIM(wall_gc_file), iostat=ios)
if (ios /= 0) then
! plot the AWALL boundary in blue instead
    do j=1, n_theta
        pol_ang = GP2*(j - 1)/float(n_theta)
        Zwall = AB*ELONM*SIN(pol_ang)
        Rwall = RTOR + AB*(COS(pol_ang) + 0.5*TRICH*(COS(2.*pol_ang) - 1.))
        xwall(j) = Rwall*meter2pixel
        ywall(j) = pixel_ymid - Zwall*meter2pixel
    enddo
    call plot_curve(jgc, 0, LineWidth, xwall, ywall)
    write(*, *) '>>> plot_wall: problems opening file ' // TRIM(wall_gc_file)
else
    read(7, *) STRI
    read(7, *) nSHOT
    read(7, *) Ndim
    if (Ndim <= 750) then
        read(7, *) NGC
    endif
    if (Ndim <= 750 .and. NGC <= 40) then
        read(7, *) ((xyGC(j, j1), j1=1, 2), j=1, Ndim)
        read(7, *) STRI
        read(7, *) (ixbeg(j), j=1, NGC)
        read(7, *) STRI
        read(7, *) (lenix(j), j=1, NGC)
        read(7, *) STRI
        read(7, *) (valix(j), j=1, NGC)
        close(7)
        if (debug > 0) then
             write(*, *) "Plotting device wall contour from " // TRIM(wall_gc_file), meter2pixel
        endif
        do j1=1, NGC
            if (valix(j1) /= White) then
                ndim_gc = lenix(j1)
                do j=1, ndim_gc
                    xgc(j) = meter2pixel*xyGC(ixbeg(j1)+j-1, 1)
                    ygc(j) = pixel_ymid - meter2pixel*xyGC(ixbeg(j1)+j-1, 2)
                enddo
                call plot_curve(ndim_gc, 0, xgc(1:ndim_gc), ygc(1:ndim_gc))
            endif
        enddo
    else
        close(7)
        write(*, *) "Configuration file is too long"
    endif
endif

return
end subroutine plot_wall

!======================================================================|
double precision function ABSC(YIN)
!----------------------------------------------------------------------|
! Input: MODEX, YIN, FP
! Output: Value a=YIN mapped to the current abscissa
!----------------------------------------------------------------------|

use outcmn_inc, only: MOD10
use status_inc, only: AMETR, FP_NORM
use const_inc, only: XOUT, AB, ABC, ROC, NA1, PSIAX, PSIBO
use numerical_tools, only: QUADIN

implicit none

double precision, intent(in) :: YIN

integer :: MODEX, j
double precision :: YAB, RFA, Y

MODEX = XOUT + 0.49

SELECT CASE(MODEX)
CASE(0) ! Draw up to AB against "a"
    YAB = AB
CASE(1) ! Draw up to ABC against "a"
    YAB = ABC
CASE(2) ! Draw up to ROC against "rho"
    YAB = ROC
CASE(3) ! Draw up to ROC against "Psi"
    YAB = 1.
CASE DEFAULT ! Unknown option
    return
END SELECT

if (MOD10 == 3 .or. MODEX == 3) then
    ABSC = QUADIN(NA1, AMETR, FP_NORM, YIN)
    return
endif

if (MODEX == 0 .or. MODEX == 1) then
    ABSC = min(1.d0, max(0.d0, YIN/YAB))
elseif (MODEX == 2) then
    ABSC = RFA(YIN)/YAB
else
    write(*, *) "MODEX is neither 0, nor 1, nor 2, it cannot be"
endif

end function ABSC

!======================================================================|
subroutine DRAWFOOT(jifnew)
!----------------------------------------------------------------------|
! YRBMN maximum radius of the footprint 
! YRBMX minimum radius of the footprint 
! YHBM  the upshift of the beam footprint
! YASP  the aspect ratio of the beam footprint
! YQ    the beam power
!----------------------------------------------------------------------|

use parameter_inc, only: NRD
use const_inc, only: CNB1
use outcmn_inc, only: NBFILE, XWH, XWW, Magenta, DYLET, meter2pixel, pixel_ymid
use debugger, only: markloc, astra_stop

implicit none

integer, intent(in) :: jifnew

integer :: j, jj, JL, JN, ERCODE
double precision :: YRBMN, YRBMX, YHBM, YASP, YH, YQ
double precision, dimension(10) :: plot_arr
double precision, dimension(15, 2*NRD+7) :: work_nbi
character(len=16) :: STRI
character(len=132) :: err_msg
!----------------------------------------------------------------------|

call markloc('DRAWFOOT')

call createpixmap(2) ! All calls except the first are ignored

open(2, file=TRIM(NBFILE), status='OLD')

if (anint(CNB1) < 1) return

do JN=1, anint(CNB1)
! Read a record for one beam source
    call STREAD(2, 20, work_nbi(1, JN), ERCODE)
    if (ERCODE /= 0) EXIT
enddo

SELECT CASE(ERCODE)
CASE(1)
    err_msg = '>>> NB drawing >>> Error in file "' // TRIM(NBFILE) // &
         '": unrecognized variable name'
    call astra_stop(err_msg)
CASE(2)
    call astra_stop('>>> NB file "' // TRIM(NBFILE) // '" read error')
CASE(3)
    write(*, *) '>>> NB file "' // TRIM(NBFILE) // '" wrong format:'
    call astra_stop('            More records expected than available.')
CASE(4) 
    call astra_stop('>>> NB drawing STREAD: array out of limits')
CASE(5)
    write(*, *) '>>> NB file "', TRIM(NBFILE), '" wrong format:'
    write(err_msg, '(A, I2, A, I2, A)')'     ', JN-1, ' beam records available, ', &
         int(CNB1), ' records required.'
    call astra_stop(err_msg)
END SELECT

close(2)
if (jifnew == 0) call redraw(2) ! Erase previous
call cleare(2, 0, 0, XWW - 1, XWH - 1)
call colovm(Magenta)

JL = 0

do JN=1, anint(CNB1)

    YQ = work_nbi(1, JN)

    if (YQ >= 1.d-2) then
        YHBM = work_nbi(11, JN)
        YASP = work_nbi(15, JN)
        YRBMX = work_nbi(12, JN)
        YRBMN = work_nbi(13, JN)
! Beam footprint drawing:
        JL = JL + 1
        write(STRI(1: 2), '(I2)') JN
        STRI(3: 4) = '  '
        write(STRI(5: 10), '(1F6.3)') YQ

        call textnb(16, 35 - 3 + JL*DYLET, STRI(1: 10), 10)
        call textnb(10, 35 - 3, 'Beam  Power', 11)

        YH = 0.5*YASP*(YRBMX - YRBMN)
        plot_arr(1) = max(YRBMN, 0.d0)*meter2pixel
        plot_arr(2) = pixel_ymid + (-YHBM - YH)*meter2pixel
        plot_arr(3) = YRBMX*meter2pixel
        plot_arr(4) = plot_arr(2)
        plot_arr(5) = plot_arr(3)
        plot_arr(6) = pixel_ymid + (-YHBM + YH)*meter2pixel
        plot_arr(7) = plot_arr(1)
        plot_arr(8) = plot_arr(6)
        plot_arr(9) = plot_arr(1)
        plot_arr(10) = plot_arr(2)

        j = 2
        jj = 5
        call drawcurve(j, plot_arr, jj)

    endif

enddo

return
end subroutine DRAWFOOT

!======================================================================|
subroutine DRAWSPFLUX
!----------------------------------------------------------------------|
! Redraw magnetic surfaces:

use outcmn_inc, only: Magenta, pixel_ymid, meter2pixel
use const_inc, only: NEQUIL, MEQUIL
use parameters_a2equil, only: equil_now

implicit none

integer, parameter :: n_surf=556, nrho_plot=12
integer :: jrho, jr, nskip, jrho_loc, n_theta, n_theta1, n_rho_surf
double precision, dimension(n_surf) :: xplot, yplot
double precision, dimension(nrho_plot+1, n_surf) :: xplot_old, yplot_old

save xplot_old, yplot_old

n_rho_surf = NINT(NEQUIL)
n_theta = NINT(MEQUIL)
n_theta1 = n_theta + 1
nskip = 1 + n_rho_surf/nrho_plot

jr = 1
do jrho=1, n_rho_surf + nskip - 1, nskip
    jrho_loc = min(jrho, n_rho_surf)
    xplot(1: n_theta) = meter2pixel*equil_now%coord_sys%position%r(jrho_loc, 1: n_theta)
    xplot(n_theta1)   = meter2pixel*equil_now%coord_sys%position%r(jrho_loc, 1)  ! Close polygon
    yplot(1: n_theta) = pixel_ymid - meter2pixel*equil_now%coord_sys%position%z(jrho_loc, 1:n_theta)
    yplot(n_theta1)   = pixel_ymid - meter2pixel*equil_now%coord_sys%position%z(jrho_loc, 1)
    call update_curve(n_theta1, n_theta1, Magenta, 0, xplot_old(jr, 1:n_theta1), &
         yplot_old(jr, 1:n_theta1), xplot(1:n_theta1), yplot(1:n_theta1))
    xplot_old(jr, 1:n_theta1) = xplot(1:n_theta1)
    yplot_old(jr, 1:n_theta1) = yplot(1:n_theta1)
    jr = jr + 1
enddo

return
end subroutine DRAWSPFLUX
