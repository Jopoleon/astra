! 10.03.95 G.P.
! 19.11.20 GIT major f90 cleaning
!---------------------------------------------------------------------
integer function IFKEY_(IFKL)
! IFKEY for calls from c

implicit none

integer, intent(in) :: IFKL
integer :: IFKEY

IFKEY_ = IFKEY(IFKL)

end function IFKEY_

!---------------------------------------------------------------------
integer function IFKEY(IFKL)
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

use parameter_inc, only: NRD, NARRX
use status_inc, only: MU, AMETR, SHIF, ELON, TRIA, EQFF, EQPF, FP, RHO
use const_inc, only: KEY, ITREQ, DROUT, DTOUT, DPOUT, exp_header, &
   NA, NB1, NA1, NAB, NUF, LEQ, TIME, TAU, TINIT, TSCALE, &
   TSTART, TPAUSE, TEQ, DTEQ, HRO, AB, ABC, ROC, XOUT, RTOR, &
   BTOR, IPL, constValues, varValues, internValues, &
   NUFR, NBNDR
use outcmn_inc, only: astra_gui, astra_gui_ref, plot_area, resizeGraph, &
    Black, Blue, Magenta, WarningColor, &
    active_tab, curves_per_frame, MOD10, LTOUT, IPOUT, MODEY, &
    NWINDX, NWIND1, NWIND3, NWIND4, NWIND7, &
    NROUT, NTOUT, NXOUT, NST, NDTNAM, NRW, &
    NAMER, NAMET, NAMEX, SCALER, SCALET, ROUT, OSHIFR, OSHIFT, &
    DTNAME, runid, VERSION, AVERS, ARLEAS, AEDIT, &
    GRAP, GRAL, TIM7, NAM7, KPRI, nplots_max, &
    NTIMES, TTOUT, TOUT
use io_mod, only: n_sbr, n_bnd, NGR, equ_file, exp_file, TASK, jbeg_arrx, IFDFVX
use dbl2char, only: fmt6
use char_manip, only: str_in_list, null_ch, beep_ch
use debugger, only: markloc, debug, astra_stop
use json_vars, only: internNames, constNames, varNames, n_const, n_var, n_intern
use cpu_usage, only: cpu_report
use numerical_tools, only: smooth

implicit none

integer, parameter :: n_portrait=0, n_landscape=1
integer, intent(in) :: IFKL
character(len=10), parameter :: DEFUNA='      .tmp'

logical :: MODADD, skip_poll, ps_exists
integer :: POLLEVENT, WAITEVENT, KIBM, KASCII, jpos, jps
integer :: MARK, J, JJ, NNN, LTOUTO, JTOUT, IDSP, &
    IFLAG, INT4, IRET, plot_mode, &
    MODEX, IX, IY, NU1, j2, J1, ios, &
    YEAR, MONTH, DAY, HOUR, MINUTE, time_arr(8)
! plot_arr dimension: 4*NRD(Mode 5, 8) 320(7) 2*NTIMES(Mode 6) 2*NRD(Modes 1-4)
integer :: ITO(NTIMES, nplots_max+2)
double precision :: LINEAV, CHORDN, ABD, ALFA, TIMEB, TROUT, TPOUT=0.d0
double precision, allocatable :: varValues_old(:) 
double precision, dimension(1) :: rescale_array
double precision, dimension(NTIMES) :: PRMARK, TIMOD4
double precision, dimension(NRD) :: YWA, YWB, YWC
character(len=6) :: NAMEP(NTIMES)
character(len=10), dimension(NRW) :: UNAMES
character(len=40) :: CNSFIL
character(len=80) :: HELP(28), STR, STRB
character(len=132) :: STRI, ps_root, PSNAME
character(len=7), dimension(1), parameter :: rescale_label = (/ 'Rescale' /)
integer, external :: plotMode

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
data (HELP(j), j=1, 10)/  &
      'Esc or / - STOP', &
      '<space>  - Set Pause mode, one-step advance', &
      ' <CR>    - Return to Run mode', &
      '  0      - No graphic output', &
      '  1,2,3  - Radial profile graphs', &
      '  4,5    - Time evolution of radial profiles', &
      '  6      - Time evolution of local/global quantities', &
      '  7      - Trace of the discharge in a phase space', &
      '  8      - Magnetic flux surfaces', &
      '  .      - Curve style'/
data (HELP(j), j=11, 20)/ &
      '  A      - Adjust colors (for a color monitor only)', &
      '  B & N  - Backward & forward screen scan', &
      '  C & V  - Constants & main Variables control', &
      '  D      - Time step, output times & subroutine call control', &
      '  F & P  - Store displayed curves in a file', &
      '  G & Q  - Save the screen as a PostScript file', &
      '  H & ?  - Show operative keys', &
      '  I      - Save the model constants for the next run', &
      '  L      - Model listing', &
      '  M      - Mark time slice[s] in the modes 4, 5, 7'/
data (HELP(j), j=21, 28)/ &
      '  R      - Refresh screen', &
      '  S      - New scales', &
      '  T      - Type numerical values of the current curves', &
      '  U      - Write U-file; (1D in modes 1-3,6; 2D in modes 4,5)', &
      '  W      - Re-arrange windows', &
      "  X      - What's X-axis/grid?", &
      '  Y      - Shift curve up/down', &
      ' '/
!----------------------------------------------------------------------|

call markloc('IFKEY', debug_lev=2*debug)

allocate(varValues_old(n_const))

CHORDN = lineav()

IFKEY = 0

if (IFKL == -1) then
! This sets "pause" mode each time when IFKEY(-1) is called
    TASK(1:3) = 'DSP'
endif

if (IFKL < 0 .or. IFKL > 256) then
    write(*, *) '>>> IFKEY: wrong input parameter. Call ignored.'
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
    call graph_output(MARK, PRMARK, NAMEP, ITO)
endif

if (IFKL /= 256 .and. TASK(4:4) /= 'B') call TIMEDT(TIME, 1000.*TAU)

if (LTOUT > 1) then
    call markloc(str_in='IFKEY (saving time traces)')
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
        call graph_output(MARK, PRMARK, NAMEP, ITO)
    endif
endif

! Time output
if (MOD10 == 6 .or. MOD10 == 7) then
    call graph_output(MARK, PRMARK, NAMEP, ITO)
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
                    call graph_output(MARK, PRMARK, NAMEP, ITO)
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
        call ERASXY() ! git (IX, IY)
        return

    CASE(32) ! 'space'
        KEY = 0
        if (TASK(1:3) == 'DSP') then
            IFLAG = 1    ! for DSP mode only
            IFKEY = 0
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
        call refresh_plot(IFKL, MARK, PRMARK, PSNAME)

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
        call refresh_plot(IFKL, MARK, PRMARK, PSNAME)

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
        call refresh_plot(IFKL, MARK, PRMARK, PSNAME)

    CASE(67) ! 'C'
        call MENUTABLE(n_const, constValues, constNames, 2)

    CASE(68) ! 'D'
        NDTNAM = n_intern + 4*n_sbr
        TIMEB = TIME
        MODEX = XOUT + 0.49
        call MENUTABLE(NDTNAM, internValues, DTNAME, 3) ! Only place requiring internValues(j>44)
        NUF   = int(NUFR)
        n_bnd = int(NBNDR)
        j = XOUT + 0.49
        if (j < 0 .or. j > 3) then
            write(*, *) ">>> Unknown X-axis. Redefinition ignored"
            j = MODEX
            XOUT = MODEX
        endif
        if (j /= MODEX) call xaxis(j)
        if (TIME >= TIMEB) then
            call refresh_plot(IFKL, MARK, PRMARK, PSNAME)
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
        ps_exists = .True.
        jps = 0
        do while(ps_exists)
            jps = jps + 1
            write(PSNAME, '(A, i0, A)') TRIM(ps_root), jps, '.ps'
            inquire(file=TRIM(PSNAME), exist=ps_exists)
        enddo
        if (KEY == 71) INT4 = n_portrait
        if (KEY == 81) INT4 = n_landscape
        call PSOPEN(TRIM(PSNAME) // null_ch, INT4, IRET)

        if (IRET == 0) then
            if (KEY == 71) KPRI = 1
            if (KEY == 81) KPRI = 2
            call refresh_plot(IFKL, MARK, PRMARK, PSNAME)
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
            write(*, *) '>>> IFKEY: file "', TRIM(CNSFIL), '" open error'
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
        call refresh_plot(IFKL, MARK, PRMARK, PSNAME)

    CASE(76) ! 'L'
        CNSFIL = 'src/tmp/model.txt'
        open(1, file=TRIM(CNSFIL), iostat=ios)
        if (ios /= 0) then
            write(*, *) '>>> IFKEY: "', TRIM(CNSFIL), '" file error'
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
        if (MOD10 == 1 .or. MOD10 == 2 .or. MOD10 == 3) call ASXWIN(NROUT, NWIND1, NAMER, SCALER, &
            OSHIFR, GRAL, GRAP, MOD10, MODEY)
        if (MOD10 == 6) call ASTWIN(NTOUT, NWIND3, NAMET, SCALET, &
            OSHIFT, MOD10, MODEY)
        if (MOD10 == 7) then
            call MENUTABLE(4, TIM7, NAM7, 5)
        endif
        if (MOD10 == 4 .or. MOD10 == 5) then
            INT4 = -MAX(4, IPOUT-1)
            call MENUTABLE(INT4, PRMARK, NAMEP, 6)
        endif
        if (MOD10 <= 7) then
            call refresh_plot(IFKL, MARK, PRMARK, PSNAME)
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
        call refresh_plot(IFKL, MARK, PRMARK, PSNAME)

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
        call refresh_plot(IFKL, MARK, PRMARK, PSNAME)

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

        call refresh_plot(IFKL, MARK, PRMARK, PSNAME)

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
        IFKEY = 0
        return
    endif
    if (TASK(1:3) == 'DSP') CYCLE
    if (jj == 1) then
        IFKEY = 0
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
end function ifkey

!---------------------------------------------------------------------
subroutine graph_output(MARK, PRMARK, NAMEP, ITO)

use const_inc, only: NA
use status_inc, only: MU
use io_mod, only: TASK
use outcmn_inc, only: MOD10, nplots_max, NTIMES, TTOUT, TOUT
use debugger, only: markloc, debug

implicit none

integer, intent(in) :: MARK
integer, intent(inout) :: ITO(NTIMES, nplots_max+2)
double precision, intent(in), dimension(NTIMES) :: PRMARK
character(len=6) , intent(in) :: NAMEP(NTIMES)

integer :: jt
double precision :: CHORDN, lineav

call markloc('graph_output', debug_lev=2*debug)

if (TASK(1:3) == 'BGD' .or. TASK(4:4) == 'B') return

call TIMOUT
call RADOUT

if (MOD10 == 4 .or. MOD10 == 5) then
    print*, 'Plot modes 4-5 not available'
else
    call OUTDSP(MARK, 0, ITO, TTOUT, TOUT)
endif
CHORDN = lineav()
call up_label(CHORDN, 1./MU(NA))
jt = 0
if (MOD10 <= 7) call down_label(jt, TOUT)
call redraw

end subroutine graph_output

!---------------------------------------------------------------------
subroutine refresh_plot(IFKL, MARK, PRMARK, PSNAME)
! Corresponds to block from statement 201

use io_mod, only: TASK
use outcmn_inc, only: astra_gui, KPRI, MOD10, MODEY, RUNID, &
    WarningColor, nplots_max, resizeGraph, NTIMES, TOUT, TTOUT
use const_inc, only: XOUT, TIME, TAU, NA
use status_inc, only: MU
use debugger, only: markloc, debug

implicit none

integer, parameter :: nn80=80

integer, intent(in) :: IFKL, MARK
double precision, intent(in), dimension(NTIMES) :: PRMARK
character(len=*) :: PSNAME

integer :: plot_mode, NST, j
integer :: ITO(NTIMES, nplots_max+2)
double precision :: CHORDN, lineav
character(len=6) :: NAMEP(NTIMES)
character(len=132) :: STRI
integer, external :: plotMode

call markloc('refresh_plot', debug_lev=2*debug)

call erasrw

plot_mode = 1
NST = 0

plot_mode = plotMode(MOD10, MODEY)
call set_plot_area(plot_mode)
call set_plot(plot_mode)

j = XOUT + 0.49
call TaskMenu(j)
call textbf(0, astra_gui%Height - int(104.*resizeGraph), RUNID, 80) ! Task ID
if (IFKL == 256 .and. TASK(1: 3) /= 'DSP') then
    call PSCLOSE
    return
endif
call RADOUT
call TIMOUT

if (MOD10 == 4 .or. MOD10 == 5) then
    print*, 'Plot mode 4-5 not available anymore'
else
    call OUTDSP(MARK, 1, ITO, TTOUT, TOUT)
endif
CHORDN = lineav()
call up_label(CHORDN, 1./MU(NA))
if (IFKL /= 256) call TIMEDT(TIME, 1000.*TAU) ! 256 <-> initial iterations
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
double precision function LINEAV

! LINEAV [10#19/m#3]: Horizontal chord average density (r) [m]
! Integral {0, r} ( NE ) dl / a

use const_inc, only: NA, ABC, NA1
use status_inc, only: AMETR, NE

implicit none

integer j

LINEAV = 2.*AMETR(1)*NE(1)
do j=2, NA
    LINEAV = LINEAV + (AMETR(j) - AMETR(j-1))*(NE(j) + NE(j-1))
enddo
LINEAV = 0.5*(LINEAV + (ABC - AMETR(NA))*(NE(NA1) + NE(NA)))/ABC

return
end function lineav

!---------------------------------------------------------------------
subroutine menutable(arr_size, array_in, var_names, id)

use char_manip, only: null_ch

implicit none

integer, intent(in) :: arr_size, id
double precision, intent(in), dimension(arr_size) :: array_in
character(len=6), intent(in), dimension(arr_size) :: var_names

integer :: nameLength, editable=1
character(len=70), dimension(10), parameter :: titles = (/ &
    'Variable control', 'Constant control', 'Times & Grids', 'Sequence control', &
    'Time interval', 'Mark times:  < 0 - skip,  0 - dim,  > 0 - color #', &
    'Equilibrium control', '1D_Ufile', '2D_Ufile', 'NBI const for beam No' /)

if (id == 4) then
    nameLength = 4
else
    namelength = 6
endif

call menubox(TRIM(titles(id)) // null_ch, arr_size, array_in, var_names, nameLength, id, editable)

return
end subroutine menutable
