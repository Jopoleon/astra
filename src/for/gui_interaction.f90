module gui_interaction

implicit none
contains

!---------------------------------------------------------------------
    integer function if_key(IFKL)
!---------------------------------------------------------------------
! IFKL = 256 call from initial iteration loop,
!         move to KEY analysis skipping D[PRT]OUT checks,
!   no drawing unless in "DSP" mode
! IFKL = 0 standard call
!
! IFKL = -1 set PAUSE mode
!
! IFKL = ichar('G') equivalent to pressing "G" in interactive mode
!
! The following logical units are used
! 1 - for files equ/log/model (once on entry and on request "I")
!  tmp/model.txt (on request "L")
! 12,13 - for equ/model.log file (once on entry)
! 3 - for post-viewer file (first on entry, then periodically)

    use parameter_inc, only: NRD
    use status_inc, only: MU, AMETR, SHIF, ELON, TRIA, EQFF, EQPF, FP, RHO
    use const_inc, only: KEY, ITREQ, DROUT, DTOUT, DPOUT, exp_header, &
       NA, NB1, NA1, NAB, LEQ, TIME, TAU, TINIT, TSCALE, &
       TSTART, TPAUSE, TEQ, DTEQ, HRO, AB, ABC, ROC, XOUT, RTOR, &
       BTOR, IPL, constValues, varValues, internValues
    use graph_utils, only: astra_gui, astra_gui_ref, plot_area, &
        Black, Blue, Magenta, WarningColor, &
        active_tab, curves_per_frame, MOD10, LTOUT, IPOUT, MODEY, &
        NWINDX, NWIND1, NWIND3, NWIND4, NWIND7, &
        NROUT, NTOUT, NXOUT, NST, NDTNAM, NRW, &
        NAMER, NAMET, NAMEX, SCALER, SCALET, ROUT, OSHIFR, OSHIFT, &
        DTNAME, runid, VERSION, AVERS, ARLEAS, AEDIT, &
        GRAP, GRAL, TIM7, NAM7, KPRI, nplots_max, &
        NTIMES, TTOUT, TOUT, ASTWIN, ASXWIN, ASKINT, MENUTABLE, &
        set_plot_area, set_plot, plotMode
    use io_mod, only: n_sbr, equ_file, exp_file, TASK, jbeg_arrx, IFDFVX
    use char_manip, only: null_ch, beep_ch
    use debugger, only: markloc, debug, astra_stop
    use json_vars, only: internNames, constNames, varNames, n_const, n_var, n_intern
    use cpu_usage, only: cpu_report
    use auxiliary, only: lineav

    integer, parameter :: n_portrait=0, n_landscape=1
    integer, intent(in) :: IFKL

    logical :: skip_poll
    integer :: POLLEVENT, WAITEVENT, KIBM, KASCII
    integer :: MARK, J, JJ, NNN, LTOUTO, JTOUT, IDSP, &
        IFLAG, INT4, IRET, plot_mode, &
        MODEX, IX, IY, NU1, j2, J1, ios, &
        YEAR, MONTH, DAY, HOUR, MINUTE, time_arr(8)
! plot_arr dimension: 4*NRD(Mode 5, 8) 320(7) 2*NTIMES(Mode 6) 2*NRD(Modes 1-4)
    integer :: ITO(NTIMES, nplots_max+2)
    double precision :: CHORDN, ABD, ALFA, TIMEB, TROUT, TPOUT=0.d0
    double precision, allocatable :: varValues_old(:) 
    double precision, dimension(NTIMES) :: PRMARK
    character(len=6) :: NAMEP(NTIMES)
    character(len=10), dimension(NRW) :: UNAMES
    character(len=40) :: CNSFIL
    character(len=80) :: HELP(28), STR, STRB
    character(len=132) :: STRI, ps_root, PSNAME

    save ITO, IFLAG, TROUT, MARK, LTOUTO, IDSP
    save NAMEP
    data PRMARK/NTIMES*0./  TROUT/-99999./ &
         IFLAG/0/  &
         JTOUT/0/ LTOUTO/0/ MARK /0/       IDSP/0/

! ASCII codes: ^C 3  <Esc>27 <Space>32  % 37  * 42  . 46  / 47  ? 63
!   0 48  1 49  2 50  3 51  4 52  5 53  6 54  7 55  8 56  9 57
!   A 65  B 66  C 67  D 68  E 69  F 70  G 71  H 72  I 73  J 74
!   K 75 L 76  M 77  N 78  O 79  P 80  Q 81  R 82  S 83  T 84
!   U 85  V 86  W 87  X 88  Y 89  Z 90
    data (HELP(j), j=1, 24)/  &
          'Esc or / - STOP', &
          '<space>  - Set Pause mode, one-step advance', &
          ' <CR>    - Return to Run mode', &
          '  0      - No graphic output', &
          '  1,2,3  - Radial profile graphs', &
          '  4,5    - Time evolution of radial profiles', &
          '  6      - Time evolution of local/global quantities', &
          '  7      - Trace of the discharge in a phase space', &
          '  8      - Magnetic flux surfaces', &
          '  .      - Curve style', &
          '  A      - Adjust colors (for a color monitor only)', &
          '  B & N  - Backward & forward screen scan', &
          '  C & V  - Constants & main Variables control', &
          '  D      - Time step, output times & subroutine call control', &
          '  F & P  - Store displayed curves in a file', &
          '  G & Q  - Save the screen as a PostScript file', &
          '  H & ?  - Show operative keys', &
          '  I      - Save the model constants for the next run', &
          '  L      - Model listing', &
          '  M      - Edit plot layout', &
          '  W      - Re-arrange windows', &
          "  X      - What's X-axis/grid?", &
          '  Y      - Shift curve up/down', &
          ' '/
!----------------------------------------------------------------------|

    call markloc('IF_KEY', debug_lev=2*debug)

    IF_KEY = 0

    allocate(varValues_old(n_const))

    CHORDN = lineav()

    if (IFKL == -1) then
! This sets "pause" mode each time when IF_KEY(-1) is called
        TASK(1:3) = 'DSP'
    endif

    if (IFKL < 0 .or. IFKL > 256) then
        write(*, *) '>>> IF_KEY: wrong input parameter. Call ignored.'
        return
    endif

    skip_poll = .False.
    if (IFKL > 0 .and. IFKL < 256) then
        KEY = IFKL
        skip_poll = .True.
        goto 1
    elseif (IFKL == 256) then
        write(STRI, '(a, i3)') "Iteration #", ITREQ
        call setColor(Magenta) ! Iterations
        call textvm(astra_gui%width-18*astra_gui_ref%dxlet, 2, "equil iterations", 16)
        call setColor(Blue) ! Iteration #
        call textvm(astra_gui%width-17*astra_gui_ref%dxlet, astra_gui_ref%dylet+1, STRI(1:14), 14)
        TROUT = TIME
        call graph_output(MARK, ITO)
    endif

    if (IFKL /= 256 .and. TASK(4:4) /= 'B') call time_label(TIME, 1000.*TAU)

    if (LTOUT > 1) then
        call markloc(str_in='IF_KEY (saving time traces)')
        if (LTOUT >= NTIMES) then
            do J=1, NTIMES-1
                do JJ=1, NTOUT
                    TOUT(J, JJ) = TOUT(J+1, JJ)
                enddo
                TTOUT(J) = TTOUT(J+1)
            enddo
            LTOUT = NTIMES - 1
        endif
    endif

    call TIMOUT

    TTOUT(LTOUT) = TIME
    LTOUT = LTOUT + 1
    JTOUT = JTOUT + 1

! Radial output
    if (MOD10 <= 3 .or. MOD10 >= 8) then
        if (TIME + .5*TAU >= TROUT + DROUT) then
            TROUT = TIME
            call graph_output(MARK, ITO)
        endif
    endif

! Time output
    if (MOD10 == 6 .or. MOD10 == 7) then
        call graph_output(MARK, ITO)
    endif

!-------------
! Key analysis
!-------------

     1 continue

    do while(.True.)
        if (.not. skip_poll) then
            KEY = 0
            if (TASK(4:4) /= 'B') call redraw

! Check Pause time condition
            if (TIME >= TPAUSE .and. IDSP == 0) then
                IDSP = 1
                KEY = 32
            else ! Polling events
                KEY = 0
                if (TASK(1:3) == 'RUN') then
                    KIBM = pollevent(KEY)
! KIBM = 1 - <Ctrl> was pressed
! KIBM = 2 - <Alt>  was pressed
! KEY=318 - the root window was closed
! KEY=322 and then KEY=319 - the root window is opened

                    if (KEY == 0) return
                    if (KIBM == 1 .and. (KEY == 99 .or. KEY == 67)) then ! <Ctrl>+C
                        if (TASK(4:4) /= 'B') call Close_Screen
                        call cpu_report('>>> ASTRA <Ctrl>+C exit >>>')
                        call astra_stop
                    endif

                    if (KIBM == 65006) then
                        TASK = 'RUN '
                        KIBM  = 0
                        return
                    endif
                    if (KIBM == 65005) then
                        TASK = 'RUNB'
                        return
                    else
                        TASK = 'RUN '
                    endif
                    if (KIBM == 1 .or. KIBM == 2) goto 49
                    if (KEY >= 97 .and. KEY <= 122) KEY = KEY - 32
                endif

!-------------------
! Waiting for events

                if (TASK(1:3) == 'DSP') then
                    if (IFLAG == 1) then
! IFLAG is equal 1 if
! (1) <Space> is pressed 
! (2) when a subroutine-call KEY is pressed.
! One time step and re-drawing is done
                        KEY = 0
                        IFLAG = 0
                        call graph_output(MARK, ITO)
                    endif

                    call PUTXY(IX, IY)
                    KASCII = 0
                    KIBM = waitevent(KASCII, ix, iy)
                    KEY  = KASCII
                    if (KEY == 0) CYCLE
                    if (KEY  < 127 .and. (KIBM == 1 .or. KIBM == 2)) goto 49
                    if (KIBM < 65000) CYCLE
                    KIBM = KIBM - 65000
                    if (KIBM >= 361 .and. KIBM <= 364) then
                        call mvcursor(KIBM, ix, iy)
                        CYCLE
                    endif
                    if (KEY >= 97) KEY = KEY - 32
                endif
            endif ! Poll key event
        endif

        SELECT CASE(KEY)

        CASE(13) ! 'ESC'
            TASK = 'RUN '
            call rcurso
            call ERASXY()

        CASE(32) ! 'space'
            KEY = 0
            if (TASK(1:3) == 'DSP') then
                IFLAG = 1    ! for DSP mode only
                IF_KEY = 0
                skip_poll = .False.
                return
            endif
            if (TASK(1:3) == 'RUN') then
                TASK = 'DSP '
                ix = 0
                iy = 0
                call pcurso
                skip_poll = .False.
            endif

        CASE(37) ! '%'
            call cpu_report(null_ch)

        CASE(46) ! '.'
            MARK = MARK + 1
            if (MARK == 2) MARK = -1
            call refresh_plot(IFKL, MARK, PSNAME)

        CASE(47) ! '/'
            if (TASK(4:4) /= 'B') call Close_Screen
            call cpu_report('>>> ASTRA / or "Quit" button exit >>>')
            call astra_stop

        CASE(48: 57) ! '0', '1', '2', '3', '4', '5', '6', '7', '8', '9'
            if (MOD10 >= 2 .and. MOD10 <= 5) then
                if (MOD10 == KEY - 48) then
                    MODEY = -MODEY
                else
                    MODEY = 1
                endif
            endif
            if (MOD10 == 1 .and. KEY == 49) active_tab(MOD10) = 0
            if (MOD10 == 6) then
                if (KEY == 54) then
                    MODEY = MODEY + 1
                    if (MODEY == 2) MODEY = -1
                else
                    MODEY = 1
                endif
            endif
            if (MOD10 /= KEY - 48) then ! Just changed plotting mode
                MOD10 = KEY - 48
                call erasrw
                plot_mode = 1
                NST = 0
                plot_mode = plotMode(MOD10, MODEY)
                call set_plot_area(plot_mode)
                call set_plot(plot_mode)
            endif
            call refresh_plot(IFKL, MARK, PSNAME)

        CASE(63, 72) ! 'H', '?'
            write(*, *)
            write(*, *) "The following keys are operable in this mode:"
            do J=1, 28
                write(*, '(1X, A)') TRIM(HELP(J))
            enddo

        CASE(66) ! 'B'
            if (MOD10 >= 1 .and. MOD10 <= 8) then
                active_tab(MOD10) = active_tab(MOD10) - 1
                if (active_tab(MOD10) < 0) then
                    JJ = curves_per_frame(MOD10)
                    if (MOD10 == 6 .or. MOD10 == 7) then
                        J = NTOUT
                    else
                        J = NROUT
                    endif
                    if (JJ == 0) JJ = 1
                    active_tab(MOD10) = (J - 1)/JJ
                endif
            endif
            call refresh_plot(IFKL, MARK, PSNAME)

        CASE(67) ! 'C'
            call MENUTABLE(n_const, constValues, constNames, 2)

        CASE(68) ! 'D'
            NDTNAM = n_intern + 4*n_sbr
            TIMEB = TIME
            MODEX = XOUT + 0.49
            call MENUTABLE(NDTNAM, internValues, DTNAME, 3) ! Only place requiring internValues(j>44)
            j = XOUT + 0.49
            if (j < 0 .or. j > 3) then
                write(*, *) ">>> Unknown X-axis. Redefinition ignored"
                j = MODEX
                XOUT = MODEX
            endif
            if (j /= MODEX) call xaxis(j)
            if (TIME >= TIMEB) then
                call refresh_plot(IFKL, MARK, PSNAME)
            else
                TROUT = TIME
                TTOUT(LTOUT-1) = TIME
                TPOUT = TIME
                do J=1, n_sbr
                    TEQ(J) = TIME
                enddo
            endif

        CASE(70) ! 'F'
            call TIMOUT
            call writeData(CHORDN)

         CASE(71, 81) ! 71:'G'=portrait, 81:'Q'=landscape
            ps_root = 'dat/' // TRIM(exp_file) // '-' // TRIM(equ_file) // '-'
            PSNAME = set_filename(ps_root, '.ps')
            if (KEY == 71) INT4 = n_portrait
            if (KEY == 81) INT4 = n_landscape
            call PSOPEN(TRIM(PSNAME) // null_ch, INT4, IRET)

            if (IRET == 0) then
                if (KEY == 71) KPRI = 1
                if (KEY == 81) KPRI = 2
                call refresh_plot(IFKL, MARK, PSNAME)
                if (IFKL == KEY) return
            else
                if (IRET == 1) then
                    STRI = '>>>  Can not open file: ' // TRIM(PSNAME)
                    call setColor(WarningColor)
                    call textvm(0, astra_gui%yMessage, STRI, 24+LEN_TRIM(PSNAME))
                endif
                KEY = 0
            endif

        CASE(73) ! 'I'
            CNSFIL = 'equ/log/' // TRIM(equ_file)
            open(1, file=TRIM(CNSFIL), iostat=ios)
            if (ios /= 0) then
                write(*, *) '>>> IF_KEY: file "', TRIM(CNSFIL), '" open error'
            else
! New format of the equ/MODEL.log file for versions => 5.3
               write(1, '(3(1A, 1I1))') ' Start file for version ', AVERS, '.', ARLEAS, '.', AEDIT
               write(1, *) 'Variables:'
               do J=1, n_var
                   if (varNames(J) == 'ZRD1  ') EXIT
                   write(1, '(1A6, 1A2, 1P, 8E11.3)') varNames(J), ' =', varValues(J)
               enddo

               write(1, '(A)')' Constants:'
               do J=1, n_const
                   write(1, '(1A6, 1A2, 1P, 8E11.3)') constNames(J), ' =', constValues(J)
               enddo
               write(1, '(A, I2)') ' Control parameters:', 22
               do J=1, 22   ! Don't save TPAUSE and TEND
                   write(1, '(1A6, 1A2, 1P, 8E11.3)') internNames(J), ' =', internValues(J)
               enddo
               close (1)
               write(*, *) "Default start file is modified"
           endif

! Test field
        CASE(74) ! 'J'
            call system("ipcs -s") ! Report active semaphore sets
            call system("ipcs -m") ! Report active shared memory segments
            call refresh_plot(IFKL, MARK, PSNAME)

        CASE(76) ! 'L'
            CNSFIL = 'src/tmp/model.txt'
            open(1, file=TRIM(CNSFIL), iostat=ios)
            if (ios /= 0) then
                write(*, *) '>>> IF_KEY: "', TRIM(CNSFIL), '" file error'
                if (KEY == 27)  then
                    write(*, '(/2A)') 'Use key "/" for exit', beep_ch ! Beep
                elseif (KEY /= 0 .and. KIBM == 0) then
                    write(*, *) 'Unrecognized key: "', char(KEY), '"', KEY, beep_ch
                endif
                KEY = 0
                CYCLE
            endif
            j2 = 0

            do
                if (j2 < 0) EXIT
                do J1=1, int(1.333*plot_area%height/astra_gui%dylet) - 1
                    if (j2 >= 0) then
                        read(1, '(1A80)', iostat=ios) STR
                        if (ios < 0) j2 = -1
                    endif
                    NNN = (J1 - 1)*astra_gui%dylet + 1
                    if (j2 < 0) then
                        STRB = repeat(' ', 35)
                        write(*, '(1X, A)') TRIM(STRB)
                    else
                        write(*, '(1X, A)') TRIM(STR)
                    endif
                enddo
            enddo
            close (1)

        CASE(77) ! 'M'
            if (MOD10 == 1 .or. MOD10 == 2 .or. MOD10 == 3) call ASXWIN(AB, NROUT, &
                NWIND1, NAMER, SCALER, OSHIFR, GRAL, GRAP, MOD10, MODEY)
            if (MOD10 == 6) call ASTWIN(NTOUT, NWIND3, NAMET, SCALET, &
                OSHIFT, MOD10, MODEY)
            if (MOD10 == 7) then
                call MENUTABLE(4, TIM7, NAM7, 5)
            endif
            if (MOD10 == 4 .or. MOD10 == 5) then
                INT4 = -MAX(4, IPOUT-1)
            endif
            if (MOD10 <= 7) then
                call refresh_plot(IFKL, MARK, PSNAME)
            endif

        CASE(78) ! 'N'
            if (MOD10 >= 0 .and. MOD10 < 7) then
                active_tab(MOD10) = active_tab(MOD10) + 1
                JJ = curves_per_frame(MOD10)
                if (MOD10 == 6 .or. MOD10 == 7) then
                    J = NTOUT
                else
                    J = NROUT
                endif
                if (J <= JJ*active_tab(MOD10)) active_tab(MOD10) = 0
            endif
            call refresh_plot(IFKL, MARK, PSNAME)

        CASE(86) ! 'V'
            do J=1, n_var
                varValues_old(J) = varValues(J)
            enddo
            INT4 = n_var - 96  ! INT4 = n_var - No. of ZRDs
            call MENUTABLE(INT4, varValues, varNames, 1)
            do J=1, n_var
                if (IFDFVX(J) > 3) varValues(J) = varValues_old(J)
                if (ABS(varValues(J)-varValues_old(J)) > 1.d-6*ABS(varValues(J))) IFDFVX(J) = 3
            enddo

        CASE(87) ! 'W'
            if (MOD10 == 1 .or. MOD10 == 2 .or. MOD10 == 3) call ASKINT(NROUT, NWIND1, NAMER)
            if (MOD10 == 4 .or. MOD10 == 5) call ASKINT(NROUT, NWIND4, NAMER)
            if (MOD10 == 6) call ASKINT(NTOUT, NWIND3, NAMET)
            if (MOD10 == 7) call ASKINT(NTOUT, NWIND7, NAMET)
            call refresh_plot(IFKL, MARK, PSNAME)

        CASE(88) ! 'X'
            MODEX = XOUT + 0.49
            if (MOD10 == 0) then
                write(*, *) 'X-axis:   none'
            elseif (MOD10 == 3) then
                write(*, '(A, I)') 'X-axis:   poloidal flux,  1 < j < NA1 =', NA1
            elseif (MOD10 == 4) then
                write(*, '(A, F5.2, A, I)') 'X-axis:   0 < a < AB =', AB, 'm     1 < j < NAB =', NAB
            elseif (MOD10 == 5) then
                write(*, '(A, I)') 'X-axis:   major radius in the mid-plane'
            elseif (MOD10 == 6) then
                write(*, *) 'X-axis:   time [s]'
            elseif (MOD10 == 7) then
                write(*, *) 'X-axis:   phase space'
            elseif (MOD10 == 8) then
                write(*, *) 'X-axis:   major radius [m]'
            elseif (MOD10 == 9) then
                write(*, *) "User's plot"
            elseif (MODEX == 0) then
                write(*, '(A, F5.2, A, I)') 'X-axis:   0 < a < AB =', AB, 'm,     1 < j < NAB =', NAB
            elseif (MODEX == 1) then
                write(*, '(A, F5.2, A, I)') 'X-axis:   0 < a < ABC =', ABC, 'm,    1 < j < NA1 =', NA1
            elseif (MODEX == 2) then
                write(*, '(A, F5.2, A, I)') 'X-axis:   0 < rho < ROC =', ROC, 'm,    1 < j < NA1 =', NA1
            elseif (MODEX == 3) then
                write(*, '(A, I)') 'X-axis:   0 < Psi < FP(NA1),  1 < j < NA1 =', NA1
            else
                write(*, *) 'X-axis:   Unknown option'
            endif

            call refresh_plot(IFKL, MARK, PSNAME)

        END SELECT

        if (KEY >= 46 .and. KEY <= 90) then ! ASCII keys
            if (IFKL == KEY) return ! Important in call from c (tglf-like models)
        endif

        CYCLE

    49 continue

        if (KIBM == 2) then !-------- <Alt> pressed'
            if (KEY == 77 .or. KEY == 109) then
                 if (KIBM == 0) then
                     write(*, *) 'Unrecognized key: "', char(KEY), '"', KEY, beep_ch
                 endif
                 KEY = 0
                 CYCLE
            endif
            if (KEY == 47) then ! <Alt>+/
                if (TASK(4:4) /= 'B') call Close_Screen
                call cpu_report('>>> ASTRA <Alt>+/ exit >>>')
                call astra_stop
            endif
            if (KIBM == 2 .and. (KEY >= 32 .and. KEY <= 126) ) then
                write(*, *) '  "<Alt>+<', char(KEY), '>"  pressed'
            endif
        endif

        if (KIBM == 1 .and. (KEY >= 32 .and. KEY <= 126) ) then
            write(*, *) '  "<Ctrl>+<', char(KEY), '>" pressed'
        endif
        jj = 0
        if (KEY > 90) KEY = KEY - 32
        KEY = KEY - 64

        do j=1, n_sbr
            if (ABS(KEY-DTEQ(4, j)) < 0.1) jj = 1
        enddo
        if (TASK(1:3) == 'DSP' .and. jj == 1) then
            IFLAG = 1     ! for DSP mode only
            IF_KEY = 0
            return
        endif
        if (TASK(1:3) == 'DSP') CYCLE
        if (jj == 1) then
            IF_KEY = 0
            return
        endif

        if (KEY == 27)  then
            write(*, '(/2A)') 'Use key "/" for exit', beep_ch ! Beep
        elseif (KEY /= 0 .and. KIBM == 0) then
            write(*, *) 'Unrecognized key: "', char(KEY), '"', KEY, beep_ch
        endif
        KEY = 0

    enddo

! Exit ASTRA

    if (TASK(4:4) /= 'B') call Close_Screen

    call cpu_report('>>> ASTRA exit: reached END time >>>')
    call astra_stop

    return
    end function if_key

!---------------------------------------------------------------------
    subroutine graph_output(MARK, ITO)

    use const_inc, only: NA
    use status_inc, only: MU
    use io_mod, only: TASK
    use graph_utils, only: MOD10, nplots_max, NTIMES, TTOUT, TOUT, show_plots
    use debugger, only: markloc, debug
    use auxiliary, only: lineav

    integer, intent(in) :: MARK
    integer, intent(inout) :: ITO(NTIMES, nplots_max+2)

    integer :: jt
    double precision :: CHORDN

    call markloc('graph_output', debug_lev=2*debug)

    if (TASK(1:3) == 'BGD' .or. TASK(4:4) == 'B') return

    call TIMOUT
    call RADOUT

    if (MOD10 == 4 .or. MOD10 == 5) then
        print*, 'Plot modes 4-5 not available'
    else
        call show_plots(MARK, 0, ITO, TTOUT, TOUT)
    endif
    CHORDN = lineav()
    call up_label(CHORDN, 1./MU(NA))
    jt = 0
    if (MOD10 <= 7) call down_label(jt, TOUT)
    call redraw

    return
    end subroutine graph_output

!---------------------------------------------------------------------
    subroutine refresh_plot(IFKL, MARK, PSNAME)
! Corresponds to block from statement 201

    use io_mod, only: TASK
    use graph_utils, only: astra_gui, KPRI, MOD10, MODEY, RUNID, &
        WarningColor, nplots_max, NTIMES, TOUT, TTOUT, show_plots, &
        set_plot_area, set_plot, plotMode
    use const_inc, only: XOUT, TIME, TAU, NA
    use status_inc, only: MU
    use debugger, only: markloc, debug
    use auxiliary, only: lineav

    integer, intent(in) :: IFKL, MARK
    character(len=*), intent(in) :: PSNAME

    integer :: plot_mode, j
    integer :: ITO(NTIMES, nplots_max+2)
    double precision :: CHORDN
    character(len=132) :: STRI

    call markloc('refresh_plot', debug_lev=2*debug)

    call erasrw

    plot_mode = plotMode(MOD10, MODEY)
    call set_plot_area(plot_mode)
    call set_plot(plot_mode)

    j = XOUT + 0.49
    call TaskMenu(j)
    call textbf(0, astra_gui%Height - int(104.*astra_gui%resizeGraph), RUNID, 80) ! Task ID
    if (IFKL == 256 .and. TASK(1: 3) /= 'DSP') then
        call PSCLOSE
        return
    endif
    call RADOUT
    call TIMOUT

    if (MOD10 == 4 .or. MOD10 == 5) then
        print*, 'Plot mode 4-5 not available anymore'
    else
        call show_plots(MARK, 1, ITO, TTOUT, TOUT)
    endif
    CHORDN = lineav()
    call up_label(CHORDN, 1./MU(NA))
    if (IFKL /= 256) call time_label(TIME, 1000.*TAU) ! 256 <-> initial iterations
    j = 0
    if (MOD10 <= 5 .or. MOD10 == 7) call down_label(j, TOUT)
    if (MOD10 == 6 .and. KPRI == 0) call down_label(j, TOUT)
    if (KPRI == 1 .or. KPRI == 2) then
        write(*, '(/A/)') '>>>  The figure is stored in the file: ' // TRIM(PSNAME)
        call setColor(WarningColor)
        if (KPRI == 1) call const2ps
        STRI = 'The figure is stored in the file: ' // TRIM(PSNAME)
        call textvm(10, astra_gui%yMessage, TRIM(STRI), LEN_TRIM(STRI))
        call PSCLOSE
        KPRI = 0
    endif

    return
    end subroutine refresh_plot
  
!---------------------------------------------------------------------
    subroutine ERASXY()

    use graph_utils, only: astra_gui, White

    integer :: JN0, JLR
    character(len=118) :: STRI

    JLR = astra_gui%Height - 125
    JN0 = 0
    STRI = repeat(' ', 80)

    call setColor(White)
    call textvm(JN0, JLR, STRI(1: 15), 15)
    call textvm(astra_gui%Width - 83*astra_gui%dxlet + 2, JLR, STRI(1: 80), 80)

    end subroutine ERASXY

!---------------------------------------------------------------------
    subroutine writeData(CHORDN)

    use const_inc, only: XOUT, NAB, NA, NA1, exp_header, RTOR, ABC, BTOR, IPL, TIME, TAU, constValues
    use status_inc, only: MU, AMETR, RHO, FP
    use io_mod, only: AWD, equ_file, exp_file
    use graph_utils, only: LTOUT, NTOUT, NROUT, MOD10, NAMER, NAMET, &
        WarningColor, ROUT, RUNID, NRW, NTIMES, TTOUT, TOUT
    use dbl2char, only: fmt_smart
    use json_vars, only: n_const

    double precision, intent(in) :: CHORDN

    integer :: NCH=0, NP1, MODEX, JBE, JEND, J, JEN, JJ, J1, ios
    double precision :: YQ
    character(len=6) :: CH6
    character(len=118) :: STRI
    character(len=132) :: FNAME, dat_dir, file_out

    MODEX = XOUT + 0.49
    ! MODEX = 0 [0, AB]  against "a"
    ! MODEX = 1 [0, ABC] against "a"
    ! MODEX = 2 [0, ROC] against "rho"
    ! MODEX = 3 [FP(1), FP(NA1)] against "psi"
    ! otherwise Unknown option => MODEX=0
    NP1 = NAB
    if (MODEX >= 1 .and. MODEX <= 3 .or. MOD10 == 3) NP1 = NA1

    dat_dir = TRIM(AWD) // '/dat/'
    call system('mkdir -p ' // TRIM(dat_dir))
    FNAME = TRIM(dat_dir) // TRIM(exp_file) // '-' // TRIM(equ_file) // '-'
    file_out = set_filename(FNAME, '.dat')

    write(*, *) '>>>  Data are written into the file: ' // TRIM(file_out)

    open(7, file=TRIM(file_out), iostat=ios)

    if (ios /= 0) then
        write(*, *) '>>> TYPDSP: Output file error'
        stop
    endif

    YQ = 1./MU(NA)
    STRI = ' '
    STRI(1: 77) = exp_header(1: 16) // upperLabel(CHORDN, YQ) // timeLabel(TIME, 1000.*TAU)

    write(7, 104) STRI

    if (NCH == 0) then

        if (MOD10 <= 5) then   ! Writing radial data
            JBE = 1
            JEND = 16
            do
                JEN = MIN0(NTOUT, JEND)
                write(7, 102) (NAMET(J), J=JBE, JEN)
                STRI = ' '
                STRI(1:5) = fmt_smart(TIME, 4)
                do J=JBE, JEN
                    JJ = 7*(J - JBE) + 8
                    STRI(JJ: JJ+5) = fmt_smart(TOUT(LTOUT, J), 5)
                enddo
                write(7, 104) STRI
                if (JEN == NTOUT) EXIT
                JBE  = JEN + 1
                JEND = JEN + 16
            enddo

            JBE  = 1
            JEND = 16

            do
                JEN = MIN0(NROUT, JEND)
                if (MODEX == 0 .or. MODEX == 1) then
                    write(7, '("     a  ", 16(3X, 1A4))') (NAMER(J), J=JBE, JEN)
                elseif (MODEX == 2) then
                    write(7, '("     rho", 16(3X, 1A4))') (NAMER(J), J=JBE, JEN)
                elseif (MODEX == 3 .or. MOD10 == 3) then
                    write(7, '("     psi", 16(3X, 1A4))') (NAMER(J), J=JBE, JEN)
                else
                    write(7, '("     ???", 16(3X, 1A4))') (NAMER(J), J=JBE, JEN)
                endif
                do J=1, NP1
                    STRI = ' '
                    do JJ=JBE, JEN
                        J1 = 7*(JJ - JBE + 1) + 1
                        STRI(J1: J1+5) = fmt_smart(ROUT(J, JJ), 5)
                    enddo
    ! Different options for a radial variable 
                    if (MODEX == 0) then
                        STRI(1: 5) = fmt_smart(AMETR(j), 4)
                    elseif (MODEX == 1) then
                        STRI(1: 5) = fmt_smart(AMETR(j), 4)
                    elseif (MODEX == 2) then
                        STRI(1: 5) = fmt_smart(RHO(j), 4)
                    elseif (MODEX == 3 .or. MOD10 == 3) then
                        STRI(1: 5) = fmt_smart(FP(j), 4)
                    else
                        STRI(1: 5) = fmt_smart(AMETR(j), 4)
                    endif
                    write(7, 104) STRI
                enddo
                if (JEN == NROUT) EXIT
                JBE  = JEN + 1
                JEND = JEN + 16
            enddo

        else if (MOD10 == 6) then ! Writing time data

            JBE  = 1
            JEND = 16
            do
                JEN = MIN0(NTOUT, JEND)
                write(7, 102) (NAMET(J), J=JBE, JEN)
                do J1=1, LTOUT - 1
                    STRI = ' '
                    STRI(1: 5) = fmt_smart(TTOUT(J1), 4)
                    do J=JBE, JEN
                        JJ = 7*(J - JBE) + 8
                        STRI(JJ: JJ+5) = fmt_smart(TOUT(J1, J), 5)
                    enddo
                    write(7, 104) STRI
                enddo
                if (JEN == NTOUT) EXIT
                JBE  = JEN + 1
                JEND = JEN + 16
            enddo
        endif

    endif !NCH=0

! Writing constants

    write(7, '(10X, A)') TRIM(RUNID)
    write(7, '(A)') 'Constants'
    J1 = 0
    do JEN=1, 11
        STRI = ' '
        do J=1, 16
            J1 = J1 + 1
            if (J1 > n_const) EXIT
            CH6 = fmt_smart(constValues(J1), 5)
            JJ = 7*(J - 1) + 1
            STRI(JJ: JJ+5) = CH6
        enddo
        write(7, '(A)') TRIM(STRI)
    enddo

    if (NCH == 1) then

        if (MOD10 <= 5) then   ! Writing radial data
            JBE  = 1
            JEND = 16
            do
                JEN = MIN0(NTOUT, JEND)
                write(7, '(3X, "Time", 16(3X, 1A4))') (NAMET(J), J=JBE, JEN)
                STRI = ' '
                STRI(1: 5) = fmt_smart(TIME, 4)
                do J=JBE, JEN
                    JJ = 7*(J - JBE) + 8
                    STRI(JJ: JJ+5) = fmt_smart(TOUT(LTOUT, J), 5)
                enddo
                write(7, 104) STRI
                if (JEN == NTOUT) EXIT
                JBE  = JEN + 1
                JEND = JEN + 16
            enddo

            JBE  = 1
            JEND = NRW
            JEN  = MIN0(NROUT, JEND)
    ! Different options for radial variable 
            if (MODEX == 2) then
                write(7, '(8X, "rho ", 64(8X, 1A4))') (NAMER(J), J=JBE, JEN)
                do j=1, NP1
                    write(7, 408) RHO(j), (ROUT(J, JJ), JJ=JBE, JEN)
                enddo
            elseif (MODEX == 3 .or. MOD10 == 3) then
                write(7, '(8X, "psi ", 64(8X, 1A4))') (NAMER(J), J=JBE, JEN)
                do j=1, NP1
                    write(7, 408) FP(j), (ROUT(J, JJ), JJ=JBE, JEN)
                enddo
            else
                write(7, '(8X, "a   ", 64(8X, 1A4))') (NAMER(J), J=JBE, JEN)
                do j=1, NP1
                    write(7, 408) AMETR(j), (ROUT(J, JJ), JJ=JBE, JEN)
                enddo
            endif

        else if (MOD10 == 6) then  ! Writing time data

            STRI=' '
            write(7, 104) STRI
            write(7, 104) STRI
            write(7, 104) STRI
            write(7, 104) STRI
            JBE  = 1
            JEND = MIN(NTOUT, NRW)
            do while(JBE < JEND)
                JEN = JBE + 7
                write(7, '(8X, "Time", 64(8X, 1A4))') (NAMET(J), J=JBE, JEN)
                do J1=1, LTOUT-1
                    STRI = ' '
                    STRI(1: 5) = fmt_smart(TTOUT(J1), 4)
                    write(7, 408) TTOUT(J1), (TOUT(J1, J), J=JBE, min(JEN, JEND))
                enddo
                JBE = JBE + 8
            enddo
        endif
    endif

    close(7)

    102 format('   Time', 16(3X, A))
    104 format(X, A)
    408 format(PE12.3, 64(PE12.3))

    return
    end subroutine writeData

!---------------------------------------------------------------------
    integer function TimeIndex(time_in, time_arr, ntim)
! TimeIndex = 1st index of the array time_arr(1:ntim) element >= time_in
! Returns 1 if time_arr(1) >= time_in, ntim if time_arr(ntim) < time_in
! time_arr has to be monotonically increasing  

    integer, intent(in) :: ntim
    double precision, intent(in) :: time_in, time_arr(*)

    integer :: pos

    pos = findloc(time_arr(1: ntim) >= time_in, .true., dim=1)
    TimeIndex = merge(pos, ntim, pos /= 0)

    return
    end function TimeIndex

!---------------------------------------------------------------------
    subroutine PUTXY(IX, IY)
! Prints x, y coordinates on GUI in "Step" mode

    use graph_utils, only: astra_gui, plot_area, MOD10, &
        scale_bnd, NTIMES, TOUT, TTOUT, &
        IDT, IDX, MODEY, LTOUT, NRW, NTOUT, active_tab, &
        NWIND3, NAMET, White, Red, Blue
    use status_inc, only: AMETR, SHIF, ELON, TRIA, FP, RHO
    use const_inc, only: TIME, TINIT, TSCALE, NA, NA1, NAB, XOUT, AB, ABC, ROC, HRO
    use dbl2char, only: fmt_smart
    use numerical_tools, only: QUADIN
    use standard_functions, only: RZ2A

    implicit none

    integer, parameter :: JN0=0
    integer, intent(in) :: IX, IY

    integer :: JX, JY, JLR, j, j1, JC, JL, MODEX, JW, JN2
    double precision :: DX, DY, YX, YX1, YY, YY1, YA, YA1, YD, YE, YT, &
        YRHO, YFP, YFPC
    character(len=80) :: STRI

    JLR = astra_gui%Height - int(125*astra_gui%resizeGraph)
    if (MOD10 <= 0) return

    STRI = repeat(' ', 80)
    STRI(7:25) = '(x, y)=(     ,     )'
    JX = IX - 10
    JY = IY - 10
    if (plot_area%xmin > JX .or. JX > plot_area%xmax .or. plot_area%ymax > JY .or. JY > plot_area%ymin) then
        STRI = repeat(' ', 80)
        call setColor(White)
        call textvm(JN0, JLR, "               ", 15)
        call textvm(astra_gui%Width - 83*astra_gui%dxlet + 2, JLR, STRI(1: 80), 80)
        return
    endif

    if (MOD10 == 7) plot_area%nx_canvas = 2
    DX = 1./plot_area%nx_canvas
    DY = 1./plot_area%ny_canvas

    YX1 =      (JX - plot_area%xmin + 0.)/(plot_area%xmax - plot_area%xmin)
    YY1 = 1. - (JY - plot_area%ymax + 0.)/(plot_area%ymin - plot_area%ymax)

    if (MOD10 == 6) then
       ! (window_width)/(step=IDX=23)/(n_labels)=592/23/25=1.0295652
        YX = TINIT + 1.029565*YX1*abs(TSCALE)
        do j=1, plot_area%ny_canvas
            YY1 = YY1 - DY
            if (YY1 < 0) EXIT
        enddo
        YY = (YY1 + DY)/DY
        JN2 = plot_area%height + 15 + 3*astra_gui%dylet
        JC = 0
        JL = 0
        YY1 = max(TIME, TTOUT(LTOUT-1), TTOUT(LTOUT))
        YY1 = min(YX, YY1)
        YY1 = max(TTOUT(1), YY1)
        j = TimeIndex(YX, TTOUT, LTOUT)
        YX = TTOUT(j)
        do J1=1, NTOUT
            JW = NWIND3(J1) - 8*active_tab(MOD10)
            if (NAMET(J1) == '    ') JW = 0
            if (JW > 0 .and. JW <= 8) call down_label(j, TOUT) ! for all plotting modes
        enddo
        STRI(1: 11) = 'Time=' // fmt_smart(YY1, 5) // 's'
        call setColor(Red)
        call textvm(astra_gui%dxlet, JN2 - int(2.5*astra_gui%dylet), STRI, 11)
        return
    else if (MOD10 == 8) then
        YX = YX1*scale_bnd*plot_area%canvas_width/IDT/IDX
        YY = (YY1 - 0.5)*scale_bnd*plot_area%canvas_height/IDT/IDX
        YA1= RZ2A(YX, YY, NAB)
        YD = QUADIN(NAB, AMETR, SHIF, YA1)
        YE = QUADIN(NAB, AMETR, ELON, YA1)
        YT = QUADIN(NAB, AMETR, TRIA, YA1)
        STRI( 7: 12) = "(r, z)="
        STRI(32: 37) = "(a, S)="
        STRI(38: 50) = '(     ,     )'
        STRI(39: 43) = fmt_smart(YA1, 5)
        STRI(45: 49) = fmt_smart(YD, 5)
        STRI(57: 62) = "(E, T)="
        STRI(63: 75) = '(     ,     )'
        STRI(64: 68) = fmt_smart(YE, 5)
        STRI(70: 74) = fmt_smart(YT, 5)
    else
        do j=1, plot_area%nx_canvas
            YX1 = YX1 - DX
            if (YX1 < 0) EXIT
        enddo
        YX = (YX1 + DX)/DX
        if (MOD10 < 5) then
            YA = 0.
            YFP = 0.
            YRHO = 0.
            YFPC = 1.125*FP(1) - 0.125*FP(2)
            MODEX = XOUT + 0.49
            if (MOD10 == 3) MODEX = 3
            if (MOD10 == 4) MODEX = 0

            SELECT CASE(modex)
            CASE(0)
                YA = YX*AB
                YRHO = QUADIN(NA1, AMETR, RHO, YA)
                YFP  = QUADIN(NA1, AMETR, FP , YA)
            CASE(1)
                YA = YX*ABC
                YRHO = QUADIN(NA1, AMETR, RHO, YA)
                YFP  = QUADIN(NA1, AMETR, FP , YA)
            CASE(2)
                YRHO = YX*ROC
                YA   = QUADIN(NA1, RHO, AMETR, YRHO)
                YFP  = QUADIN(NA1, RHO, FP   , YRHO)
            CASE(3)
                YFP  = YFPC + (FP(NA1) - YFPC)*YX
                YA   = QUADIN(NA1, FP, AMETR, YFP)
                YRHO = QUADIN(NA1, FP, RHO  , YFP)
            END SELECT

            j = YRHO/HRO + 1
            if (YRHO > 0.5*(RHO(NA) + ROC)) j = NA1

            STRI(32: 33) = "a="
            STRI(34: 38) = fmt_smart(YA, 5)
            YRHO = YRHO/ROC
            STRI(39: 49) = 'm,   rho_t='
            STRI(50: 54) = fmt_smart(YRHO, 5)
            if (YFP > YFPC) then
                YFP = sqrt((YFP - YFPC)/(FP(NA1) - YFPC))
            else
                YFP = 0.
            endif
            STRI(55: 64) = ",   rho_p="
            STRI(65: 69) = fmt_smart(YFP, 5)
            STRI(70: 77) = ",  Node:"
            write(STRI(78: 80), '(1I3)')j
        endif

        do j=1, plot_area%ny_canvas
            YY1 = YY1 - DY
            if (YY1 < 0) EXIT
        enddo
        YY = (YY1 + DY)/DY
        if (MOD10 > 1 .and. MOD10 < 6 .and. MODEY == -1) then
            if (JY - plot_area%ymax > (plot_area%ymin - plot_area%ymax)/plot_area%ny_canvas) YY = YY - 1.
        endif

    endif

    STRI(14: 18) = fmt_smart(YX, 5)
    STRI(20: 24) = fmt_smart(YY, 5)

    call setColor(Blue)
    call textvm(JN0, JLR, "               ", 15)
    call textvm(astra_gui%Width - 83*astra_gui%dxlet, JLR, STRI(1: 80), 80)

    return
    end subroutine putxy

    
!---------------------------------------------------------------------
    subroutine down_label(jt_in, TOUT)
!---------------------------------------------------------------------
! Time dependences for radial output
! Curve to digit conversion
!  Input:  jt defines the current time
!
! if mod10 != 6 or call from run then jt = LTOUT
!---------------------------------------------------------------------

    use graph_utils, only: astra_gui, plot_area, LTOUT, MOD10, NTOUT, &
        NRW, NTIMES, NWIND3, active_tab, NAMET, Black, Blue, curves_per_frame
    use dbl2char, only: fmt_smart

    implicit none

    integer, intent(in) :: jt_in
    double precision, intent(in) :: TOUT(NTIMES, *)

    integer :: jt, JN2, JN0, JEND, JB, JL, JC, JW, JJ, J, fshift
    character(len=5) :: XF4
    character(len=7) :: XF7
    character(len=80) :: STRI, STRIN

    if (jt_in == 0) then
        jt = LTOUT
        call setColor(Black)
    else
        jt = jt_in
        call setColor(Blue)
    endif

    fshift = int(10*astra_gui%resizeGraph)

    if (MOD10 == 6) then

        JN0 = 6*astra_gui%dxlet
        JN2 = plot_area%height + FSHIFT + 3*astra_gui%dylet + 5
        JC = 0
        JL = 0
        write(STRIN, '(79X, 1A1)') ' '
        write(STRI , '(79X, 1A1)') ' '
        do j=1, NTOUT
            JW = NWIND3(j) - curves_per_frame(MOD10)*active_tab(MOD10)
            if (NAMET(j) == '    ') JW = 0
            if (JW <= 0 .or. JW > curves_per_frame(MOD10)) CYCLE
            JC = JC + 1   ! Actual curve number in the window
            jj = 8*JC - 6
            if (jj > 74) CYCLE
            if (jj >= 66) JN0 = 5*astra_gui%dxlet
            if (jj == 74) JN0 = -astra_gui%dxlet
    ! curve #, win #, chan #, mode 6, screen #
            XF7 = fmt_smart(TOUT(jt, j), 6)
            STRI (jj: jj+6) = XF7
            STRIN(jj: jj+6) = '  ' // NAMET(J) // ' '
            JL = max(JL, jj + 6)
        enddo
        call textvm(JN0, JN2, STRI, JL)
        JN2 = JN2 - astra_gui%dylet + 1
        call textvm(JN0, JN2, STRIN, JL)
    else
        JB = 1
        JN0 = 0
        JN2 = plot_area%height - 5*astra_gui%dylet + FSHIFT + 2
        do
            JEND = MIN(JB + 15, NTOUT)
            do J=JB, JEND
                XF4 = fmt_smart(TOUT(jt, J), 4)
                if (NAMET(J) == ' ') XF4 = '    '
                JJ = 5*(J - JB + 1) - 4
                STRI(JJ: JJ+4) = XF4
            enddo
            JN2 = JN2 + 2*astra_gui%dylet + 2
            JL = 5*(JEND - JB + 1)
            call textvm(JN0, JN2, STRI, JL)
            write(STRI, '(16(1X, 1A4))') (NAMET(J), J=JB, JEND)
            JN2 = JN2 - astra_gui%dylet + 1
            call textvm(JN0, JN2, STRI, JL)
            JN2 = JN2 + astra_gui%dylet - 1
            if (JEND == NTOUT .or. JEND == NRW) EXIT
            JB = JB + 16
        enddo
    endif

    return
    end subroutine down_label

!-----------------------------
    function upperLabel(ne_av, q95) result(upper_label)

    use const_inc, only: RTOR, BTOR, IPL, ABC
    use dbl2char, only: fmt_smart

    implicit none

    double precision, intent(in) :: ne_av, q95
    character(len=42) :: upper_label

    upper_label = ' R=' // fmt_smart(RTOR, 4) // ' a=' // fmt_smart(ABC, 4) // ' b=' // &
        fmt_smart(BTOR, 4) // ' I=' // fmt_smart(IPL, 4) // ' q=' // fmt_smart(q95, 4) // &
        ' n=' // fmt_smart(ne_av, 4)

    return
    end function upperLabel

!-----------------------------
    function timeLabel(time_in, dt_in) result(time_lbl)

    use dbl2char, only: fmt_smart

    implicit none

    double precision, intent(in) :: time_in, dt_in
    character(len=19) :: time_lbl

    time_lbl = 'Time=' // fmt_smart(time_in, 5) // ' dt=' // fmt_smart(dt_in, 5)

    return
    end function timeLabel

!---------------------------------------------------------------------
    subroutine up_label(YN, YQ)
! Upper string of the Astra graphic window

    use const_inc, only: exp_header
    use graph_utils, only: astra_gui, active_tab, MOD10, Black, Blue

    implicit none

    double precision, intent(in) :: YN, YQ

    character(len=2) :: CHR

    call setColor(Black)
    call rectvm(0, 0, 0, astra_gui%Width - 1, astra_gui%Height - 1) ! Draw outer frame
    call textvm(0, 2, exp_header(1: 16) // upperLabel(YN, YQ), 58)  ! Type upper label

    write(CHR, '(1I2)') active_tab(MOD10) + 1
    call setColor(Blue)
    call textvm(astra_gui%width - 2*astra_gui%dxlet, astra_gui%dylet + 1, CHR, 2) ! Type screen No.

    return
    end subroutine up_label

!---------------------------------------------------------------------
    subroutine time_label(time_in, dt_in)

    use graph_utils, only: astra_gui, astra_gui_ref, Black

    implicit none

    integer, parameter :: fshift=2, str_len=19
    double precision, intent(in) :: time_in, dt_in

    call setColor(Black)
    call textvm(astra_gui%width - (str_len+3)*astra_gui_ref%dxlet, fshift, timeLabel(time_in, dt_in), str_len)

    return
    end subroutine time_label

!---------------------------------------------------------------------
    subroutine const2ps
! Appending the list of constants to a PS file

    use const_inc, only: constValues, varValues
    use char_manip, only: null_ch
    use io_mod, only: resize
    use dbl2char, only: fmt_smart
    use json_vars, only: n_const, n_var, varNames

    implicit none

    character(len=6), dimension(22), parameter :: CONN = (/ &
        'CF1-> ', 'CF5-> ', 'CF9-> ', 'CF13->', &
        'CV1-> ', 'CV5-> ', 'CV9-> ', 'CV13->', &
        'CHE   ', 'CHI   ', 'CNB   ', 'CNBI  ', &
        'CCD   ', 'CRF   ', 'CNEUT ', 'CPEL  ', &
        'CBND  ', 'CFUS  ', 'CIMP  ', 'CMHD  ', &
        'CRAD  ', 'CSOL  ' /)

    integer :: J, J1, J2, JJ, JDUM, JNY, JNX
    character(len=6 ) :: CH6
    character(len=80) :: STRI

! Writing constants

    JNY = 540*resize
    JNX = 10*resize

    ps_loop: do J2=1, 11

        JDUM = 2*(J2 - 1) + 1
        STRI = CONN(JDUM)

        do J=1, 4
            J1 = J1 + 1
            if (J1 > n_const) EXIT ps_loop
            CH6 = fmt_smart(constValues(J1), 5)
            JJ = 7*(J - 1) + 8
        enddo

        STRI(JJ: JJ+5) = CH6
        STRI(JJ+12: JJ+18) = CONN(JDUM+1)

        do J=5, 8
            J1 = J1 + 1
            if (J1 > n_const) EXIT ps_loop
            CH6 = fmt_smart(constValues(J1), 5)
            JJ = 7*(J - 1) + 20
        enddo

        STRI(JJ: JJ+5) = CH6
        JNY = JNY + 17

        call textvm(JNX, JNY, STRI, 75)

    enddo ps_loop

! Writing variables
    j1 = 1
    JNY = JNY + 20
    JDUM = n_var - 48

    do j=1, JDUM
        CH6 = fmt_smart(varValues(j), 5)
        STRI(j1: j1+19) = varNames(j) // '=' // CH6 // '     '
        j1 = j1+20
        if (j1 > 70 .or. j == JDUM) then
            STRI(j1-5:) = null_ch
            JNY = JNY + 17
            call textvm(JNX, JNY, STRI, j1-5)
            j1 = 1
        endif
    enddo

    return
    end subroutine const2ps

!---------------------------------------------------------------------
    function set_filename(fname_in, ext_in) result(fname_out)
        implicit none

        character(len=*), intent(in) :: fname_in, ext_in
        character(len=128) :: fname_out

        integer :: jext
        logical :: fileExists
        character(len=:), allocatable :: filename
        character(len=16) :: ext   ! long enough for big integers

        fileExists = .true.
        jext = 0

        do while (fileExists)
            jext = jext + 1
            write(ext, '(i0, A)') jext, ext_in
            filename = trim(fname_in) // trim(ext)
            inquire(file=TRIM(filename), exist=fileExists)
        enddo

        fname_out = filename
    end function set_filename

end module gui_interaction

!---------------------------------------------------------------------
integer function ifkey(key)

use gui_interaction, only: if_key

implicit none

integer, intent(in) :: key

ifkey = if_key(key)

return
end function ifkey
