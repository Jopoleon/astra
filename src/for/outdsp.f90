subroutine initMainWindow

use outcmn_inc, only: astra_gui_ref, astra_gui, plot_area_ref, plot_area, &
    RUNID, NST, MOD10, MODEY, NTOUT, resizeGraph
use const_inc, only: XOUT, NA
use status_inc, only: MU

implicit none

integer :: j, jj, plot_mode
double precision :: CHORDN
double precision, external :: LINEAV
integer, external :: plotMode

call get_runid()
! Resize
jj = max(0, (15 + NTOUT - 64)/16)

astra_gui%LineWidth = int(0.85*resizeGraph) + astra_gui_ref%LineWidth
astra_gui%dxlet    = resizeGraph*astra_gui_ref%dxlet
astra_gui%dylet    = resizeGraph*astra_gui_ref%dylet
astra_gui%yMessage = resizeGraph*astra_gui_ref%yMessage + 135
astra_gui%Width    = resizeGraph*astra_gui_ref%width
astra_gui%Height   = resizeGraph*(astra_gui_ref%Height + 2*jj*resizeGraph*(astra_gui_ref%dylet + 2))
astra_gui%Xpos  = astra_gui_ref%Xpos
astra_gui%Ypos  = astra_gui_ref%Ypos
astra_gui%title = astra_gui_ref%title
astra_gui%resizeGraph = resizeGraph

plot_area%width  = resizeGraph*plot_area_ref%width
plot_area%height = resizeGraph*plot_area_ref%height

!call resizewindow(astra_gui%width, astra_gui%height)
call initvm(astra_gui%xpos, astra_gui%ypos, astra_gui%Width, astra_gui%Height, astra_gui%LineWidth, astra_gui%title, LEN(astra_gui%title)) ! Initialise graphic window
plot_mode = 1
NST = 0
MOD10 = 1
plot_mode = plotMode(MOD10, MODEY)
call set_plot_area(plot_mode)
call set_plot(plot_mode)

j = XOUT + 0.49

call taskmenu(j) ! Task menu
call textbf(0, astra_gui%Height - int(104*resizeGraph), RUNID, 80) ! Task ID

CHORDN = LINEAV()
call up_label(CHORDN, 1./MU(NA))
  
return
end subroutine initMainWindow

!---------------------------------------------------------------------
subroutine OUTDSP(MARK, JIFNEW, IYO, TT_out, t_out)

!---------------------------------------------------------------------
! Drawing options:
!   X - axis
!      0 <= rho <= ROC=RHO(NA1)
!      0 <= a <= ABC
!      0 <= a <= AB
!
! MARK = 1     Put marks
! MARK = 0     Solid lines
! MARK =-1     Dashed lines
! JIFNEW =  0 Re-draw (erase) the previous curves
! JIFNEW = 1 New curves only
!---------------------------------------------------------------------

use parameter_inc, only: NRD, NRDX, NARRX
use status_inc, only: AMETR, MU, SHIF, ELON, TRIA
use const_inc, only: XOUT, NAB, NA1, NA1E, ABC, TINIT, TSCALE, RTOR, &
    MEQUIL, LEQ, TIME
use io_mod, only: IFDFAX, NPTM, XAXES, DATAX, equ_file, TOUTX
use outcmn_inc, only: astra_gui, plot_area, &
    curves_per_frame, active_tab, MOD10, NWIND1, NWIND3, NWINDX, &
    KPRI, nplots_max, NTIMES, NRW, &
    NROUT, ROUT, OSHIFR, NAMER, SCALER, &
    NTOUT, OSHIFT, NAMET, SCALET, &
    NXOUT, NAMEX, LTOUT, &
    GRAL, GRAP, pixel_ymid, meter2pixel, &
    Black, WarningColor, EraseColor, Red, Blue, Green, White
use expdat, only: raw_profile_map, DATARR
use dbl2char, only: fmt_xf
use char_manip, only: len_trim_tab, str_in_list
use debugger, only: markloc, debug, astra_stop
use json_vars, only: profxNames

implicit none

integer, parameter :: jzero=0
integer, intent(in) :: MARK, JIFNEW
integer, intent(inout) :: IYO(NTIMES,*)
double precision, intent(in) :: TT_out(NTIMES)
double precision, intent(inout) :: t_out(NTIMES, NRW)

integer :: PTM(2), PTMO(2, NRDX, NRW), &
    IWN(16), fshift, &
    IST, text_posx, text_posy, JS, MODEX, &
    IYM0, LTOUT1, LTOUT2, STYL, x_shift, y_shift, jx_canv, jy_canv, JY, jxout, &
    jplot_in_tab, j_curve, j_canv, &
    IYMN, IYMX, JDSP, NPTMO(NRW), jlx(8), &
    NP1, j, half_wid, &
    j1, jj, jsco, jn, jpnt, jsc, jarr, jtyp, n_canvas, &
    jplot_in_canv, jcol, jsym, jprof, jtrace
double precision :: SC(NRW), YX, r_out, YA, YL, YR, &
     YZ, abscissa, ymin, ymax, px_rmag, yq1, xq1, xte, te_bc
double precision ,dimension(2) :: xbar, xbar_old, ybar, x8bar, y8bar
double precision, dimension(16) :: xq1_old, xte_old
double precision, dimension(NRD) :: xplot, yplot
double precision, dimension(NTIMES) :: xtrace, ytrace, xtrace_old
double precision, dimension(NRD, nplots_max) :: xold, yold
double precision, dimension(NTIMES, nplots_max) :: ytrace_old
double precision, external :: AFVAL
character(len=80) :: STRI
character(len=5 ) :: XF4
character(len=6 ) :: CHAR6

save PTMO, NPTMO, IWN, xq1_old, xte_old, xold, yold, xtrace_old, ytrace_old

!---------------------------------------------------------------------

call markloc('OUTDSP')

fshift = 12
half_wid = plot_area%width/2
YA = 0.
if (JIFNEW /= 0) then
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
IYMN = plot_area%height - plot_area%ymin
IYM0 = plot_area%ymin - plot_area%canvas_height
IYMX = plot_area%height - plot_area%ymax
JY = 10*plot_area%height

ymin = dble(plot_area%height - plot_area%ymin)
ymax = dble(plot_area%height - plot_area%ymax)

n_canvas = plot_area%nx_canvas*plot_area%ny_canvas

SELECT CASE(MOD10)

!-----------------
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
        jplot_in_tab = NWIND1(jprof) - curves_per_frame(MOD10)*active_tab(MOD10) ! 1-16 for mode '1'
        if (NAMER(jprof) == '    ') jplot_in_tab = 0
        if (jplot_in_tab <= 0 .or. jplot_in_tab > curves_per_frame(MOD10)) CYCLE plot_prof
        j_canv = MOD(jplot_in_tab - 1, n_canvas) + 1        ! 1-8 for mode '1'
        jx_canv = MOD(j_canv - 1, plot_area%nx_canvas)      ! 0-3 for mode '1'
        jy_canv = (j_canv - 1)/plot_area%nx_canvas          ! 0-1
        jplot_in_canv = (jplot_in_tab - 1)/n_canvas
        jcol = jplot_in_canv + 2
        y_shift = plot_area%canvas_height * (plot_area%ny_canvas - 1 - jy_canv)
        x_shift = plot_area%canvas_width*jx_canv

        YL = max(0.d0, abscissa(GRAL(jprof)))
        YR = abscissa(GRAP(jprof))
        if (YL >= YR) YL = 0.d0

! Translate curve into pixel
        jxout = 0
        do j=1, NP1
            YX = abscissa(AMETR(j))
            if (YX >= YL .and. YX <= YR) then
                jxout = jxout + 1
                if (jxout == 1 .and. j > 1) then ! left edge interpolation
                    YA = ROUT(J, jprof) + (ROUT(J-1, jprof) - ROUT(J, jprof))*(YL - YX)/(YA - YX)
                    r_out = min(max(YA/SC(jprof), -7.d0), 7.d0)
                    JDSP  = 10*(plot_area%canvas_height*r_out + IYMN + y_shift)
                    xplot(jxout) = dble(x_shift)
                    yplot(jxout) = plot_area%height - min(max(dble(plot_area%canvas_height)*r_out + ymin + dble(y_shift), ymin), ymax)
                    jxout = jxout + 1
                endif
                r_out = min(max(ROUT(J, jprof)/SC(jprof), -7.d0), 7.d0)
                JDSP  = 10*(plot_area%canvas_height*r_out + IYMN + y_shift)
                xplot(jxout) = dble(x_shift) + dble(plot_area%width)/dble(plot_area%nx_canvas)*(YX - YL)/(YR - YL)
                yplot(jxout) = plot_area%height - min(max(dble(plot_area%canvas_height)*r_out + ymin + dble(y_shift), ymin), ymax)
            endif
            if (YA <= YR .and. YX > YR) then ! right edge interpolation
                jxout = jxout + 1
                YA = ROUT(J, jprof) + (ROUT(J-1, jprof) - ROUT(J, jprof))*(YR - YX)/(YA - YX)
                r_out = min(max(YA/SC(jprof), -7.d0), 7.d0)
                JDSP  = 10*(plot_area%canvas_height*r_out + IYMN + y_shift)
                xplot(jxout) = dble(x_shift) + dble(plot_area%width)/dble(plot_area%nx_canvas)
                yplot(jxout) = dble(plot_area%height) - min(max(dble(plot_area%canvas_height)*r_out + ymin + dble(y_shift), ymin), ymax)
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
        if (j_curve <= nplots_max) then
            call update_curve(jxout, IWN(jplot_in_tab), jcol, STYL, xold(1:, j_curve), yold(1:, j_curve), xplot, yplot)
            xold(1: jxout, j_curve) = xplot(1: jxout)
            yold(1: jxout, j_curve) = yplot(1: jxout)
        endif
        IWN(jplot_in_tab) = jxout
        do J=1, NP1
            ROUT(j, jprof) = ROUT(j, jprof) - OSHIFR(jprof)
        enddo
! Variable labels
        call setColor(jcol)
        text_posx = x_shift
        text_posy = astra_gui%dylet + FSHIFT + (plot_area%ymin - plot_area%ymax + 3 + astra_gui%dylet)*jy_canv
        call CMARK(text_posx, text_posy, SC(jprof), OSHIFR(jprof), NAMER(jprof), STYL, jplot_in_canv)

    enddo plot_prof

! Plot dots
    do j=1, n_canvas
        jlx(j) = 0
    enddo
    jsco = 0
    jn = 0

    plot_profx: do jxout=1, NXOUT
        if (NWINDX(jxout) == 0) CYCLE plot_profx ! NWINDX set in ininam.f90
        CHAR6 = NAMEX(jxout)
        if (CHAR6(1: 1) == ' ') CYCLE plot_profx
        do jprof=1, NARRX
            if (profxNames(jprof) == CHAR6) jn = jprof 
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
        jplot_in_tab = NWIND1(jsc) - curves_per_frame(MOD10)*active_tab(MOD10)
        if (jplot_in_tab <= 0 .or. jplot_in_tab > curves_per_frame(MOD10)) CYCLE plot_profx
        j_canv = MOD(jplot_in_tab - 1, n_canvas) + 1        ! 1-8 for mode '1'
        jx_canv = MOD(j_canv - 1, plot_area%nx_canvas)      ! 0-3 for mode '1'
        jy_canv = (j_canv - 1)/plot_area%nx_canvas          ! 0-1
        jcol = (jplot_in_tab - 1)/n_canvas + 2
        y_shift = plot_area%canvas_height * (plot_area%ny_canvas - 1 - jy_canv)
        x_shift = plot_area%canvas_width*jx_canv

        call setColor(jcol)
        if (jsc /= jsco) then
            jsym = 5  ! 13
        else
            jsym = jsym + 1
        endif

        if (JIFNEW == 0) then
            call setColor(EraseColor)
            do j=1, NPTMO(jxout)
                call NMARK(PTMO(1, j, jxout), jsym)
            enddo
        endif

        if (KPRI >= 1 .and. KPRI <= 2) then
            write(STRI, '(1A6, 1A6, 1A1)')'Dots "', profxNames(jn), '"'
            j = len_trim_tab(STRI)
            call pscom(STRI, j)
        endif
        YL = max(0.d0, abscissa(GRAL(jsc)))
        YR = abscissa(GRAP(jsc)) 
        if (YL >= YR) YL = 0.d0
        j1 = 0
        do j=1, jpnt
            YX = XAXES(j, jn)
            if (MODEX >= 1 .and. YX > ABC) CYCLE
            YA = abscissa(YX)
            if (YA > 1. .or. YA < YL .or. YA > YR) CYCLE
            j1 = j1 + 1
            PTM(1) = x_shift + plot_area%width/plot_area%nx_canvas*(YA-YL)/(YR-YL)
            r_out= max((DATAX(j, jn) + OSHIFR(jsc))/SC(jsc), -7.d0)
            r_out= min(r_out, 7.d0)
            JDSP = plot_area%canvas_height*r_out + IYMN + y_shift
            PTM(2) = plot_area%height - min(max(JDSP, IYMN), IYMX)
            call setColor(jcol, 2)
            call NMARK(PTM, jsym)
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
        text_posx = x_shift + plot_area%canvas_width - astra_gui%dxlet - 45
        text_posy = (1 + jlx(j_canv))*astra_gui%dylet + FSHIFT + (plot_area%ymin - plot_area%ymax - plot_area%canvas_height)*(jy_canv) + 3
        XF4 = fmt_xf(TOUTX(jn), 4)
        call textvm(text_posx, text_posy, XF4, 5) ! Text (time) -> plot legend
        PTM(1) = text_posx + astra_gui%dxlet + 37 ! 12 is fixed, as the font size does not scale
        PTM(2) = text_posy - 0.3*astra_gui%dylet
        jcol = (jplot_in_tab - 1)/n_canvas + 2
        call NMARK(PTM, jsym) ! Marker symbol -> plot legend

        jsco = jsc
        if (jsym == 7) jsym = 1
    enddo plot_profx

! Erase/put q=1 radius, BC for Te
    yq1   = abscissa(AFVAL(MU, 1.0))
    te_bc = abscissa(AMETR(max(NA1E, 1)))
    ymax = dble(IYM0) - 0.8*plot_area%canvas_height
    do j_canv=1, plot_area%nx_canvas
        xq1 = plot_area%canvas_width*(j_canv -1 + YQ1)
        xte = plot_area%canvas_width*(j_canv -1 + te_bc)
        do jy=1, plot_area%ny_canvas
            ybar = (/ dble(IYM0) - (jy-2)*plot_area%canvas_height, ymax - (jy-2)*plot_area%canvas_height/)
            if (yq1 > 1.d-3 .and. yq1 < 0.999) then
                xbar_old = xq1_old(j_canv)
                xbar = xq1
                call update_curve(2, 2,   Red, 0, xbar_old, ybar, xbar, ybar)
            endif
            if (te_bc > 1.d-2 .and. te_bc < 0.999) then
                xbar_old = xte_old(j_canv)
                xbar = xte
                call update_curve(2, 2, Green, 0, xbar_old, ybar, xbar, ybar)
            endif
        enddo
        xq1_old(j_canv) = xq1
        xte_old(j_canv) = xte
    enddo

!-------------------
CASE(4:5)
    call markloc('Drawing mode 4/5', debug_lev=2*debug)

!-------------------
CASE(6)  ! Time traces

    call markloc('Drawing mode 6', debug_lev=2*debug)
    LTOUT1 = 1
    LTOUT2 = LTOUT
    if (LTOUT < 2) return
    ! right_label_position=JDX*JDMX=23*5*5=575 (see typdsp.f)
    do J=1, LTOUT-1
        r_out = (TT_out(J) - TINIT)*575/abs(TSCALE)
        IYO(J, nplots_max+1) = 6*astra_gui%dxlet + r_out
        xtrace(J) = 6*astra_gui%dxlet + r_out
        if (r_out < 0)  LTOUT1 = J + 1
        if (IYO(J, nplots_max+1) <= plot_area%width - 1) LTOUT2 = J - 1
    enddo
    LTOUT2 = LTOUT2 - LTOUT1 + 1
    if (LTOUT2 < 3) then ! No plot for small time
        return
    endif
    do jj=1, min(NTOUT, NRW)
        do j=LTOUT1, LTOUT1 + LTOUT2
            t_out(j, jj) = t_out(j, jj) + OSHIFT(jj)
        enddo
    enddo
    call SCAL(NTOUT, SC, SCALET, t_out(LTOUT1, 1), LTOUT2, NTIMES)

    j_curve = 0
    plot_traces: do jtrace=1, min(NRW, NTOUT)
        jplot_in_tab = NWIND3(jtrace) - curves_per_frame(MOD10)*active_tab(MOD10)
        if (NAMET(jtrace) == '    ') CYCLE plot_traces
        if (jplot_in_tab <= 0 .or. jplot_in_tab > curves_per_frame(MOD10)) CYCLE plot_traces
        jplot_in_canv = (jplot_in_tab - 1)/n_canvas       ! <-> color
        j_canv = MOD(jplot_in_tab - 1, n_canvas) + 1      ! 1-8 for mode '1'
        do J=1, LTOUT
            r_out = max(t_out(J, jtrace)/SC(jtrace), -7.d0)
            r_out = min(r_out, 7.d0)
            JDSP  = 10*(plot_area%canvas_height*r_out + IYMN + (n_canvas - j_canv)*plot_area%canvas_height)
            JDSP  = max(JDSP, 10*IYMN)
            IYO(J, nplots_max+2) = JY - min(JDSP, 10*IYMX)
            ytrace(J) = dble(plot_area%height) - min(max(plot_area%canvas_height*r_out + ymin + (n_canvas - j_canv)*plot_area%canvas_height, ymin), ymax)
        enddo

        if (KPRI >= 1 .and. KPRI <= 2) then
            write(STRI, '(1A6, 1A4, 1A1)') 'Plot "', NAMET(jj), '"'
            j = len_trim_tab(STRI)
            call pscom(STRI, j)
        endif
        jcol = jplot_in_canv + 1 
        STYL = jcol*MARK
        jcol  = jcol  + 1
        j_curve = j_curve + 1
        if (j_curve <= nplots_max) then
            call setColor(White)
            call plot_curve(LTOUT2, STYL, xtrace_old(1:LTOUT2), ytrace_old(1:LTOUT2, j_curve))
            call setColor(jcol)
            call plot_curve(LTOUT2+1, STYL, xtrace(1:LTOUT2+1), ytrace(1:LTOUT2+1) )
            xtrace_old(1: LTOUT2+1) = xtrace(1: LTOUT2+1)
            ytrace_old(1: LTOUT2+1, j_curve) = ytrace(1: LTOUT2+1)
        endif
        text_posx = 0
        text_posy = (j_canv-1)*plot_area%canvas_height + astra_gui%dylet*(jplot_in_canv*2 + 2) + 2*jplot_in_canv
        call CMARKT(text_posx, text_posy, SC(jtrace), OSHIFT(jtrace), NAMET(jtrace), STYL)

        do j=LTOUT1, LTOUT1 + LTOUT2
            t_out(j, jtrace) = t_out(j, jtrace) - OSHIFT(jtrace)
        enddo
    enddo plot_traces

!-------------------
CASE(7)
    call markloc('Drawing mode 7', debug_lev=2*debug)

!-------------------
CASE(8)

    call markloc('Drawing mode 8', debug_lev=2*debug)
    px_rmag = (RTOR + SHIF(1))*meter2pixel
    call setColor(Black)
    x8bar = (/0., px_rmag/)
    y8bar = (/pixel_ymid, pixel_ymid/)
    call drawcurve(0, 2, x8bar, y8bar)
! Plot the complete wall structure (Pixmap # 1)
    call plot_wall

! if (data file includes NAMEXP BND) then (NBND > 0);
! or (NBND == 8) after calling equil with no boundary points provided;
!     NBND == 0 otherwise

    SELECT CASE(LEQ(5))
    CASE(3)
        call plot_lcfs(JIFNEW, IYO, TIME)
    CASE(4: 5)
        if (MEQUIL == 0) then
            SHIF = 0.0
            ELON = 1.0
            TRIA = 0.0 
        else
            call plot_lcfs(JIFNEW, IYO, TIME)
        endif
    END SELECT
    call plot_flux_surfaces

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

        jn = str_in_list(CHAR6, profxNames)
    
        jarr = IFDFAX(jn)
        if (jarr <= 0)  CYCLE loop8
        jtyp = raw_profile_map%grid_type(jarr)
        if (jtyp < 18)  CYCLE loop8
        jpnt = raw_profile_map%nrho(jarr)
        if (jpnt <= 0)  CYCLE loop8
        jcol = jcol+1
        js = raw_profile_map%jbeg_grid(jarr)

        if (JIFNEW == 0) then
            call setColor(EraseColor)
            do j=1, NPTMO(jxout)
                call NMARK(PTMO(1, j, jxout), 7)
            enddo
        endif
        call setColor(jcol+2)
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
        call textvm(510, 100 + 2*astra_gui%dylet*jcol, NAMEX(jxout), 6)
        NPTMO(jxout) = jpnt
    enddo loop8

!-------------------
CASE(9)
    call markloc('User drawing mode', debug_lev=2*debug)

END SELECT

return
end subroutine outdsp

!---------------------------------------------------------------------
subroutine plot_lcfs(ifnew, IYO, time_in)

! Scatter plot of the LCFS
! IFNEW  =  0 Re-draw (erase) the previous curves
! IFNEW =/= 0 New curves only
! IFNEW < 0 Don't mark resonances q=m/n
! IFNEW > 10 Call from Review. (JIFNEW-10) is used to control erasing

use io_mod, only: NBNT
use outcmn_inc, only: Red, EraseColor, pixel_ymid, meter2pixel
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
    call setColor(EraseColor)
    do j=1, NBND, j2
        PTM(1) = IYO(1, j)
        PTM(2) = IYO(2, j)
        call NMARK(PTM, 4)
    enddo
endif

! Boundary points
call setColor(Red)

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
end subroutine plot_lcfs

!---------------------------------------------------------------------
subroutine plot_wall
!---------------------------------------------------------------------
! Plot vessel components reading them from json machine file

use const_inc, only: AB, ELONM, RTOR, TRICH, GP2
use outcmn_inc, only: Blue, White, pixel_ymid, meter2pixel
use debugger, only: debug
use machine_config, only: config, json_cfg, cfg_exists

implicit none

integer, parameter :: n_theta=64, ngc_max=750
integer :: j, j1, jgc, jbeg, ios, nSHOT, NGC, ndim_gc
integer, allocatable, dimension(:) :: contour_len, contour_color
double precision :: pol_ang, Rwall, Zwall
double precision, dimension(n_theta) :: xwall, ywall
double precision, dimension(ngc_max) :: xGC, yGC
double precision, allocatable, dimension(:) :: rGC, zGC
character(len=64) :: STRI

call setColor(Blue)

if (cfg_exists) then
    call config%get('contour_len', contour_len)
    call config%get('contour_color', contour_color)
    call config%get('Rvessel', rGC)
    call config%get('Zvessel', zGC)
    NGC = SIZE(contour_len)
    jbeg = 0
    do j1=1, NGC
        if (contour_color(j1) /= White) then
            ndim_gc = contour_len(j1)
            do j=1, ndim_gc
                xgc(j) = meter2pixel*rGC(jbeg+j)
                ygc(j) = pixel_ymid - meter2pixel*zGC(jbeg+j)
            enddo
            call plot_curve(ndim_gc, 0, xgc(1:ndim_gc), ygc(1:ndim_gc))
        endif
        jbeg = jbeg + contour_len(j1) 
    enddo
else
    do j=1, n_theta
        pol_ang = GP2*(j - 1)/float(n_theta)
        Zwall = AB*ELONM*SIN(pol_ang)
        Rwall = RTOR + AB*(COS(pol_ang) + 0.5*TRICH*(COS(2.*pol_ang) - 1.))
        xwall(j) = Rwall*meter2pixel
        ywall(j) = pixel_ymid - Zwall*meter2pixel
    enddo
    call plot_curve(jgc, 0, xwall, ywall)
    write(*, *) '>>> plot_wall: problems opening file ' // TRIM(json_cfg)
endif

return
end subroutine plot_wall

!---------------------------------------------------------------------
double precision function abscissa(YIN)
!---------------------------------------------------------------------
! Input: MODEX, YIN, FP
! Output: Value a=YIN mapped to the current abscissa

use outcmn_inc, only: MOD10
use status_inc, only: AMETR, FP_NORM
use const_inc, only: XOUT, AB, ABC, ROC, NA1
use numerical_tools, only: QUADIN

implicit none

double precision, intent(in) :: YIN

integer :: MODEX
double precision :: YAB
double precision, external :: RFA

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
    abscissa = QUADIN(NA1, AMETR, FP_NORM, YIN)
    return
endif

if (MODEX == 0 .or. MODEX == 1) then
    abscissa = min(1.d0, max(0.d0, YIN/YAB))
elseif (MODEX == 2) then
    abscissa = RFA(YIN)/YAB
else
    write(*, *) "MODEX is neither 0, nor 1, nor 2, it cannot be"
endif

return
end function abscissa

!---------------------------------------------------------------------
subroutine plot_flux_surfaces
!---------------------------------------------------------------------
! Redraw magnetic surfaces:

use outcmn_inc, only: Magenta, pixel_ymid, meter2pixel
use const_inc, only: NEQUIL, MEQUIL
use parameters_a2equil, only: equil_now

implicit none

integer, parameter :: n_surf=556, nrho_plot=12
integer :: jrho, jr, nskip, n_theta, n_theta1, n_rho_surf
double precision, dimension(n_surf) :: xplot, yplot
double precision, dimension(nrho_plot+1, n_surf) :: xplot_old, yplot_old

save xplot_old, yplot_old

n_rho_surf = NINT(NEQUIL)
n_theta    = NINT(MEQUIL)
n_theta1 = n_theta + 1
nskip = 1 + n_rho_surf/nrho_plot

jr = 1
do jrho=1, n_rho_surf + nskip - 1, nskip
    if (jrho > n_rho_surf) EXIT
    xplot(1: n_theta) = meter2pixel*equil_now%coord_sys%position%r(jrho, 1: n_theta)
    xplot(n_theta1)   = meter2pixel*equil_now%coord_sys%position%r(jrho, 1)  ! Close polygon
    yplot(1: n_theta) = pixel_ymid - meter2pixel*equil_now%coord_sys%position%z(jrho, 1:n_theta)
    yplot(n_theta1)   = pixel_ymid - meter2pixel*equil_now%coord_sys%position%z(jrho, 1)
    call update_curve(n_theta1, n_theta1, Magenta, 0, xplot_old(jr, 1:n_theta1), &
         yplot_old(jr, 1:n_theta1), xplot(1:n_theta1), yplot(1:n_theta1))
    xplot_old(jr, 1:n_theta1) = xplot(1:n_theta1)
    yplot_old(jr, 1:n_theta1) = yplot(1:n_theta1)
    jr = jr + 1
enddo

return
end subroutine plot_flux_surfaces
