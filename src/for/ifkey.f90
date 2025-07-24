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

use parameter_inc, only: NRD, NRW, NARRX, NCONST, NSDELOUT
use status_inc, only: MU, AMETR, SHIF, ELON, TRIA, EQFF, EQPF, FP, RHO
use const_inc, only: KEY, ITREQ, DROUT, DTOUT, DPOUT, exp_header, &
   NA, NB1, NA1, NAB, NUF, LEQ, NBND, TIME, TAU, TINIT, TSCALE, &
   TSTART, TPAUSE, TEQ, DTEQ, HRO, AB, ABC, ROC, XOUT, RTOR, &
   BTOR, IPL, CONSTF, DEVAR, DELOUT, XFLAG
use outcmn_inc, only: astra_gui, astra_gui_ref, plot_area, resizeGraph, &
    Black, Blue, Magenta, WarningColor, &
    active_tab, curves_per_frame, MOD10, LTOUT, IPOUT, MODEY, &
    NWINDX, NWIND1, NWIND3, NWIND4, NWIND7, &
    NROUT, NTOUT, NXOUT, NST, NDTNAM, &
    NAMER, NAMET, NAMEX, SCALER, SCALET, ROUT, OSHIFR, OSHIFT, &
    rev_file, DTNAME, runid, VERSION, AVERS, ARLEAS, AEDIT, &
    GRAP, GRAL, TIM7, NAM7, KPRI, ICVMX
use io_mod, only: null_ch, NSBR, NGR, equ_file, exp_file, TASK, jbeg_arrx, IFDFVX
use expdat, only: raw_profile_map, DATARR
use timeoutput_inc, only: NTIMES, TTOUT, TPOUT, TOUT
use dbl2char, only: fmt6
use char_manip, only: str_in_list
use debugger, only: markloc, debug, astra_stop
use json_vars, only: internNames, constNames, varNames, n_const, n_var
use cpu_usage, only: cpu_report

implicit none

integer, parameter :: n_portrait=0, n_landscape=1
integer, intent(in) :: IFKL
character(len=10), parameter :: DEFUNA='      .tmp'

logical :: MODADD, skip_poll, ps_exists
integer*2, dimension(NRD) :: YWD
integer :: POLLEVENT, WAITEVENT, KIBM, KASCII, jpos, jps
integer :: MARK, J, JJ, NNN, LTOUTO, JTOUT, IDSP, &
    IFLAG, INT4, IRET, plot_mode, &
    MODEX, IX, IY, NU1, j2, J1, ios, &
    YEAR, MONTH, DAY, HOUR, MINUTE, time_arr(8)
! plot_arr dimension: 4*NRD(Mode 5, 8) 320(7) 2*NTIMES(Mode 6) 2*NRD(Modes 1-4)
integer :: ITO(NTIMES, ICVMX+2)
double precision :: DEVARO(NCONST), LINEAV, CHORDN, ABD, ALFA, TIMEB, TROUT
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

!-----------------------------------------------
! Append data to post-view file, 2D Radial/Time output 

! The next line suppresses writing a view file during the iteration loop
if (TPOUT + DPOUT < TSTART .or. (IFKL /= 256 .and. TIME + 0.5*TAU >= TPOUT + DPOUT)) then

    call markloc(str_in='RADOUT|1 call from IFKEY')
    call RADOUT

    if (LTOUTO /= 0)  then
        open(3, FILE=rev_file, STATUS='OLD', iostat=ios, &
!        ACCESS='APPEND', FORM='UNFORMATTED')   ! SUN, Alpha
             POSITION='APPEND', FORM='UNFORMATTED') ! Intel
        if (ios /= 0) then
            write(*, '( // A // )') '>>> IFKEY: Review file append error'
        endif
    else
! The file rev_file='.res/profil.dat' (default name) is used in 3 places: 
!  here (w), ifkey.smode5 - modes 4, 5 (r), typdsp - writing 2D U-file (r) 
        open(unit=12, file='equ/' // TRIM(equ_file), iostat=ios)
        if (ios /= 0) write(*, *) '>>> IFKEY: Model file "equ/', TRIM(equ_file), '" open error'
        open(3, file=TRIM(rev_file), iostat=ios, form='unformatted')
        if (ios /= 0) write(*, *) '>>> IFKEY: Review file open error ' // TRIM(rev_file)
        CNSFIL = 'equ/log/' // TRIM(equ_file)
        inquire(FILE=TRIM(CNSFIL), EXIST=MODADD)
        if (MODADD) then
            open(1, file=TRIM(CNSFIL), iostat=ios)
            if (ios /= 0) write(*, *) '>>> IFKEY: File "', TRIM(CNSFIL), '" open error'
            call ADDMOD(3, 12, 1)
            close(1)
        else
            call ADDMOD(3, 12, 0)
        endif
        close(12)

        call date_and_time(VALUES=time_arr)
        YEAR   = time_arr(1)
        MONTH  = time_arr(2)
        DAY    = time_arr(3)
        HOUR   = time_arr(5)
        MINUTE = time_arr(6)

        write(3) exp_file, equ_file, VERSION, exp_header, &
            YEAR, MONTH, DAY, HOUR, MINUTE, n_const, n_var, &
            NROUT, (NAMER(J), J=1, NROUT), (SCALER(J), J=1, NROUT), &
            NTOUT, (NAMET(J), J=1, NTOUT), (SCALET(J), J=1, NTOUT), &
            HRO, NB1, NSBR, NGR, NXOUT, (LEQ(j), j=1, 7)
! Note Change the cycle in NEQNS, (LEQ(j), j=1, NEQNS)
! Presently LEQ is not used by review.f and need not be stored
!     . , HRO, NB1, NSBR, NGR, NXOUT
        if (NXOUT > 0 .and. NGR > 0) then
! Total length: 3*NGR*int+(3*NGR+GDEY(NGR)+NGRIDX(NGR)-1)*real+3*NARRX*int
            write(3) &
                (raw_profile_map%arr_index(j), j=1, NGR), (raw_profile_map%nrho(j), j=1, NGR), (raw_profile_map%grid_type(j), j=1, NGR), &
                (raw_profile_map%time(j), j=1, NGR), (raw_profile_map%jbeg_grid(j), j=1, NGR), (raw_profile_map%jbeg_data(j), j=1, NGR), &
                (DATARR(j), j=1, raw_profile_map%jbeg_data(NGR) + raw_profile_map%nrho(NGR) - 1), &
                (NAMEX(j), j=1, NARRX), (NWINDX(j), j=1, NARRX), &
                (jbeg_arrx(j), j=1, NARRX)
        endif
    endif

    LTOUTO = LTOUT - JTOUT
    write(3) JTOUT
    if (JTOUT /= 0) then
        write(3) (TTOUT(J), (TOUT(J, JJ), JJ=1, NTOUT), J=LTOUTO, LTOUT-1)
        JTOUT = 0
    endif
    write(3) TIME

    write(3) (CONSTF(J), J=1, n_const), (DEVAR(J), J=1, n_var), ABC, ROC, CHORDN, 1./MU(NA)
    write(3) NA1, NAB, (0, j=1, 10), (0.d0, j=1, 10)

    if (LEQ(5) /= 5) call RHSEQ !call this only if equil is not active

    call STUFF(3, AMETR, NAB, YWD)
    call STUFF(3, SHIF , NAB, YWD)
    call STUFF(3, ELON , NAB, YWD)
    call STUFF(3, TRIA , NAB, YWD)
    call STUFF(3, EQFF , NAB, YWD)
    call STUFF(3, EQPF , NAB, YWD)
    call STUFF(3, FP   , NAB, YWD)
    do J=1, NROUT
        call STUFF(3, ROUT(1, J), NAB, YWD)
    enddo
    close(3)

    TPOUT = TIME
    TIMOD4(IPOUT) = TPOUT
    NAMEP(IPOUT) = fmt6(TPOUT)
    if (IPOUT < NTIMES) IPOUT = IPOUT + 1

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
                    call cpu_report('>>> ASTRA <Ctrl>+C exit >>>' // char(0))
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
        call cpu_report(char(0))

    CASE(46) ! '.'
        MARK = MARK + 1
        if (MARK == 2) MARK = -1
        call refresh_plot(IFKL, MARK, PRMARK, PSNAME)

    CASE(47) ! '/'
        if (TASK(4:4) /= 'B') call Close_Screen
        call cpu_report('>>> ASTRA / or "Quit" button exit >>>' // char(0))
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
        call MENUTABLE(n_const, CONSTF, constNames, 2)

    CASE(68) ! 'D'
        NDTNAM = NSDELOUT + 4*NSBR
        TIMEB = TIME
        MODEX = XOUT + 0.49
        call MENUTABLE(NDTNAM, DELOUT, DTNAME, 3) ! Only place requiring DELOUT(j>44)
        if (int(DELOUT(13)) /= NA1) then
            write(*, *)">>> NA1 re-definition ignored"
        endif
        DELOUT(13) = NA1
        NUF   = DELOUT(14)
        NBND  = DELOUT(19)
        XFLAG = DELOUT(20)
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
            do J=1, NSBR
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
        call PSOPEN(TRIM(PSNAME) // char(0), INT4, IRET)

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
               write(1, '(1A6, 1A2, 1P, 8E11.3)') varNames(J), ' =', DEVAR(J)
           enddo

           write(1, '(A)')' Constants:'
           do J=1, n_const
               write(1, '(1A6, 1A2, 1P, 8E11.3)') constNames(J), ' =', CONSTF(J)
           enddo
           write(1, '(A, I2)') ' Control parameters:', 22
           do J=1, 22   ! Don't save TPAUSE and TEND
               write(1, '(1A6, 1A2, 1P, 8E11.3)') internNames(J), ' =', DELOUT(J)
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
        CNSFIL = 'tmp/model.txt'
        open(1, file=TRIM(CNSFIL), iostat=ios)
        if (ios /= 0) then
            write(*, *) '>>> IFKEY: "', TRIM(CNSFIL), '" file error'
            if (KEY == 27)  then
                write(*, '(/2A)') 'Use key "/" for exit', char(7) ! Beep
            elseif (KEY /= 0 .and. KIBM == 0) then
                write(*, *) 'Unrecognized key: "', char(KEY), '"', KEY, char(7)
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

    CASE(82) ! 'R'
        call refresh_plot(IFKL, MARK, PRMARK, PSNAME)

    CASE(83) ! 'S'
        rescale_array(1) = resizeGraph
        call MENUTABLE(1, rescale_array, rescale_label, 4)
        resizeGraph = rescale_array(1)
!    call initMainWindow

    CASE(84) ! 'T'
        call TIMOUT
        call TYPDSP

    CASE(85) ! 'U'
        if (NUF > NRD) then
            write(*, *) '>>> No. of radial points in a U-file is too large:'
            write(*, *) '           Decrease NUF or increase NB1'
            CYCLE
        else
            SELECT CASE(MOD10)
! Different options for radial array 
            CASE(4)
                call ufileBox(NROUT, UNAMES, NAMER, DEFUNA)
                jpos = str_in_list(DEFUNA, UNAMES(1: NROUT), 10)
                if (jpos > 0) then
                    call RADOUT
                    MODEX = XOUT + 0.49
                    if (MODEX == 1) then
! Write up to ABC against "a"
                        NU1 = NA1
                        ABD = ABC
                        do j = 1, NA1
                            YWC(j) = AMETR(j)/ABC
                        enddo
                    elseif (MODEX == 2) then
! Write up to ROC against "rho"
                        NU1 = NA1
                       ABD = 1.
                        do j = 1, NA1
                            YWC(j) = RHO(j)/ROC
                        enddo
                    else
! Write up to AB against "a" (default)
                        if (MODEX /= 0) write(*, *) ">>>  Warning: Unknown radial mode.  Writing anyway [0, AB]"
                        NU1 = NAB
                        ABD = AB
                        do j = 1, NAB
                            YWC(j) = AMETR(j)/AB
                        enddo
                    endif
! Radial coordinate in a U-file in [m] (presently)
                    do jj=1, NROUT
                        if (UNAMES(jj) /= DEFUNA) then
                            ALFA = 1.d-4
! Transfer to an equidistant radial grid
                            do j=1, NUF
                                YWA(j) = (j - 1.)/(NUF - 1.)
                            enddo
                            call SMOOTH(ALFA, NU1, ROUT(1, jj), YWC, NUF, YWB, YWA)
                            do j=1, NUF
                                YWA(j) = YWA(j)*ABD
                            enddo
                            call UF1DWA('AUGD', RUNID, UNAMES(jj), TIME, NAMER(jj), NUF, 1, 4, &
                                YWA, YWB, RTOR, AB, BTOR, IPL, CHORDN, MODEX)
                        endif
                    enddo
!  Write up to AB against "a" (default)
!                do jj=1, NROUT
!                    if (UNAMES(jj) /= DEFUNA) then
!                        call UF2DWA('AUGD', UNAMES(jj), NAMER(jj), jj, NUF, 0, 4, PRMARK, &
!                            TIMOD4, RTOR, AB, BTOR, IPL, CHORDN, YWA, YWB, YWC)
!                    endif
!                enddo
                endif

            CASE(6)
                call TIMOUT
                call ufileBox(NTOUT, UNAMES, NAMET, DEFUNA)
                do jj=1, NTOUT
                    if (UNAMES(jj) /= DEFUNA) then
                        call UF1DWA('AUGD', RUNID, UNAMES(jj), TIME, NAMET(jj), LTOUT-1, 0, &
                            4, TTOUT(1), TOUT(1, jj), RTOR, AB, BTOR, IPL, CHORDN, MODEX)
                    endif
                enddo
            CASE DEFAULT
                write(*, *) '>>>  WARNING: U-file writing is not implemented in this mode'
            END SELECT
        endif

    CASE(86) ! 'V'
        do J=1, n_var
            DEVARO(J) = DEVAR(J)
        enddo
        INT4 = n_var - 96  ! INT4 = n_var - No. of ZRDs
        call MENUTABLE(INT4, DEVAR, varNames, 1)
        do J=1, n_var
            if (IFDFVX(J) > 3) DEVAR(J) = DEVARO(J)
            if (ABS(DEVAR(J)-DEVARO(J)) > 1.d-6*ABS(DEVAR(J))) IFDFVX(J) = 3
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
                 write(*, *) 'Unrecognized key: "', char(KEY), '"', KEY, char(7)
             endif
             KEY = 0
             CYCLE
        endif
        if (KEY == 47) then ! <Alt>+/
            if (TASK(4:4) /= 'B') call Close_Screen
            call cpu_report('>>> ASTRA <Alt>+/ exit >>>' // char(0))
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

    do j=1, NSBR
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
        write(*, '(/2A)') 'Use key "/" for exit', char(7) ! Beep
    elseif (KEY /= 0 .and. KIBM == 0) then
        write(*, *) 'Unrecognized key: "', char(KEY), '"', KEY, char(7)
    endif
    KEY = 0

enddo

! Exit ASTRA

if (TASK(4:4) /= 'B') call Close_Screen

call cpu_report('>>> ASTRA exit: reached END time >>>' // char(0))
call astra_stop

return
end function ifkey

!---------------------------------------------------------------------
subroutine graph_output(MARK, PRMARK, NAMEP, ITO)

use parameter_inc, only: NRD
use const_inc, only: NA
use status_inc, only: MU
use io_mod, only: TASK, null_ch
use outcmn_inc, only: MOD10, ICVMX
use timeoutput_inc, only: NTIMES, TTOUT, TOUT
use debugger, only: markloc, debug

implicit none

integer, intent(in) :: MARK
integer, intent(inout) :: ITO(NTIMES, ICVMX+2)
double precision, intent(in), dimension(NTIMES) :: PRMARK
character(len=6) , intent(in) :: NAMEP(NTIMES)

integer :: jt
double precision :: CHORDN, lineav

call markloc('graph_output', debug_lev=2*debug)

if (TASK(1:3) == 'BGD' .or. TASK(4:4) == 'B') return

call TIMOUT
call RADOUT

if (MOD10 == 4 .or. MOD10 == 5) then
    call SMODE5(MARK, PRMARK, NAMEP, NTIMES)
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

use parameter_inc, only: NRD
use io_mod, only: TASK, null_ch  
use outcmn_inc, only: astra_gui, KPRI, MOD10, MODEY, RUNID, &
    WarningColor, ICVMX, resizeGraph
use timeoutput_inc, only: NTIMES, TOUT, TTOUT
use const_inc, only: XOUT, TIME, TAU, NA
use status_inc, only: MU
use debugger, only: markloc, debug

implicit none

integer, parameter :: nn80=80

integer, intent(in) :: IFKL, MARK
double precision, intent(in), dimension(NTIMES) :: PRMARK
character(len=*) :: PSNAME

integer :: plot_mode, NST, j
integer :: ITO(NTIMES, ICVMX+2)
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
    call SMODE5(MARK, PRMARK, NAMEP, NTIMES)
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
subroutine SMODE5(MARK, PRMARK, NAMEP)

use parameter_inc, only: NRD, NRW
use io_mod, only:  NGR, null_ch
use outcmn_inc, only: astra_gui, plot_area, NROUT, ICVMX, SCALER, &
    ROUT, rev_file, NXOUT, NAMER, NWIND4, active_tab, OSHIFR, &
    MOD10, GRAL, GRAP, MODEY, KPRI, Black, Red
use const_inc, only: AB, NAB
use dbl2char, only: fmt4
use char_manip, only: len_trim_tab
use debugger, only: markloc, debug
use timeoutput_inc, only: NTIMES

implicit none

integer, parameter :: fshift=10
integer, parameter, dimension(6) :: symbols=(/ 1, 111, 43, 42, 120, 36 /)

integer, intent(in) :: MARK
double precision, intent(in) :: PRMARK(*)
character(len=6), intent(in) :: NAMEP(*)

integer*2 :: INTY(NRD)
integer :: JTIM, ios, &
       j, jj, int2, jab, i, is, NP, NP1, jxout, jx, jnl, &
       jc, IYMN, STYL, jpos, jk(5), SKIPM, half_wid
integer, dimension(3) :: plot_arr
double precision :: &
       SC(NRW), TEMPR, SCL, DOWN, YWA(NRD), YWB(NRD), YS, YL, YR, YX, &
       YROUT, YA, YQ1, YQ2, YXR, YXL, yloc, ymin, ymax
double precision, dimension(2*NRD) :: xplot, yplot
character(len=4) :: CHAR4
character(len=5) :: axis_label
character(len=9) :: ST
character(len=80) :: STRI

call markloc('SMODE5', debug_lev=2*debug)

IYMN = plot_area%height - plot_area%ymin
ymin = dble(plot_area%height - plot_area%ymin)
ymax = dble(plot_area%height - plot_area%ymax)
half_wid = plot_area%width/2

!-------------------------------------------------------
! Mode 4 & 5:
! Plots against (a/AB) (scale does not change with time)
! Scale is determined by the current radial distribution

call SCAL(NROUT, SC, SCALER, ROUT(4, 1), NAB - 3, NRD)
JTIM = 0
open(3, file=TRIM(rev_file), iostat=ios, form='unformatted')
if(ios /= 0) then
    write(*, *) ' >>> ERROR >>> Cannot open profile file'
    stop
endif
j = 0
j = SKIPM(3, j)   ! Returned value is not used
read(3, ERR=38, END=39) CHAR4
if (NXOUT > 0 .and. NGR > 0) read(3, ERR=38, END=39) INT2

read_loop: do

    read(3, ERR=38, END=39) INT2
    if (INT2 /= 0) read(3, ERR=38) TEMPR
    read(3, END=39) TEMPR
    read(3, END=39) TEMPR
    if (JTIM > NTIMES) then
        write(*, *) ">>> Too many time slices. Data skipped"
        write(*, *) "    Use Astra post-viewer or reduce DPOUT"
        close(3)
    return
    endif

    JTIM = JTIM + 1
    read(3, END=39) JAB, JAB, (INT2, I=1, 10), (TEMPR, I=1, 10)

    read(3) SCL, DOWN, (INTY(j), j=1, JAB)
    do j=1, JAB
        YWA(j) = DOWN + (32768 + INTY(j))*SCL/65535.
    enddo

    read(3) SCL, DOWN, (INTY(j), j=1, JAB)
    do j=1, JAB
        YWB(j) = DOWN + (32768 + INTY(j))*SCL/65535.
    enddo

    do jj=1, 5
        read(3) SCL, DOWN, (INTY(j), j=1, JAB)
    enddo

    do jj=1, NROUT
        read(3, ERR=38) SCL, DOWN, (INTY(j), j=1, JAB)
        if(NAMER(jj) == '    ' .or. PRMARK(JTIM) < 0) CYCLE
        NP = NWIND4(jj) - 2*active_tab(MOD10)
        JX = (NP - 1)*plot_area%canvas_width
        if (NP /= 1 .and. NP /= 2) CYCLE
        DOWN = DOWN + OSHIFR(jj)
        jnl = astra_gui%dylet + FSHIFT
        do IS=1, 5
            jk(IS) = 0
        enddo

        ST(1:4) = fmt4(SC(jj))
        ST(5: 9) = NAMER(jj)
        call setColor(Black)
        call textvm((NP - 1)*plot_area%canvas_width, jnl, ST, 9)
        YS = OSHIFR(jj)
        if (YS < 0) then
            YS = -YS
            CHAR4 = fmt4(YS)
            ST = '-' // CHAR4
            call textvm((NP - 1)*plot_area%canvas_width + 4*astra_gui%dxlet, jnl + astra_gui%dylet + 2, ST, 5)
        else
            CHAR4 = fmt4(YS)
            ST = '+' // CHAR4
            call textvm((NP - 1)*plot_area%canvas_width + 4*astra_gui%dxlet, jnl + astra_gui%dylet + 2, ST, 5)
        endif

        if (MOD10 == 4) then
            NP1 = JAB
            jxout = NP1
            YL = GRAL(jj)/YWA(NP1)
            YR = GRAP(jj)/YWA(NP1)
            if (YL >= YR) YL = 0.d0
            if (YL < YR .and. (YL > 0.001 .or. YR < 0.999)) then
                plot_arr(1) = NP + half_wid*YL
                plot_arr(2) = plot_area%height - IYMN
                plot_arr(3) = NP + half_wid*YR
                call setColor(Red)
                call drawvm(0, plot_arr(1), plot_arr(2)    , plot_arr(3), plot_arr(2))
                call drawvm(0, plot_arr(1), plot_arr(2) + 1, plot_arr(3), plot_arr(2) + 1)
                call drawvm(0, plot_arr(1), plot_arr(2) + 2, plot_arr(3), plot_arr(2) + 2)
            endif
            jxout = 0
	    YA = 0.
            do j=1, NP1
                YX = YWA(j)/YWA(NP1)
                if (YX > YL .and. YX < YR) then
                    jxout = jxout + 1
                    if (jxout == 1 .and. j > 1) then ! left edge interpolation
                        xplot(jxout) = dble(JX) + 1.
                        YQ1 = DOWN + (32768 + INTY(j-1))*SCL/65535.
                        YQ2 = DOWN + (32768 + INTY(j))  *SCL/65535.
                        YA = YQ2 + (YQ1 - YQ2)*(YL - YX)/(YA - YX)
                        YROUT = min(max(YA/SC(jj), -7.d0), 7.d0)
                        yloc = plot_area%canvas_height*YROUT + IYMN
                        if (MODEY == -1) then
                            yloc = yloc + dble(plot_area%canvas_height)
                        endif
                        yplot = plot_area%height - min(max(yloc, ymin), ymax)
                        jxout = jxout + 1
                    endif
                    xplot(jxout) = dble(JX) + dble(half_wid)*(YX - YL)/(YR - YL) + 1.
                    YROUT = (DOWN + (32768 + INTY(j))*SCL/65535.)/SC(jj)
                    YROUT = min(max(YROUT, -7.d0), 7.d0)
                    yloc = plot_area%canvas_height*YROUT + IYMN
                    if (MODEY == -1) yloc = yloc + plot_area%canvas_height
                    yplot(jxout) = plot_area%height - min(max(yloc, ymin), ymax)
                endif

                if (YA <= YR .and. YX > YR) then ! right edge interpolation
                     jxout = jxout + 1
                     xplot(jxout) = dble(JX + half_wid) + 1. ! git + 1
                     YQ1 = DOWN + (32768 + INTY(j-1))*SCL/65535.
                     YQ2 = DOWN + (32768 + INTY(j))  *SCL/65535.
                     YA = YQ2 + (YQ1 - YQ2)*(YR - YX)/(YA - YX)
                     YROUT = min(max(YA/SC(jj), -7.d0), 7.d0)
                     yloc = plot_area%canvas_height*YROUT + IYMN
                     if (MODEY == -1) yloc = yloc + plot_area%canvas_height
                     yplot(jxout) = plot_area%height - min(max(yloc, ymin), ymax)
                endif
                YA = YX
            enddo
            NP1 = jxout
        else
! Mode 5
! Take AMETR(x, t) and SHIF(x, t) from "profile.dat"
            NP1 = 2*JAB
            jxout = NP1
            do j=1, JAB
                YROUT = (DOWN + (32768 + INTY(j))*SCL/65535.)/SC(jj)
                YROUT = min(max(YROUT, -7.d0), 7.d0)
                yloc = plot_area%canvas_height*YROUT + IYMN
                if (MODEY == -1) then
                   yloc = yloc + plot_area%canvas_height
                endif
                yplot(j) = plot_area%height - min(max(yloc, ymin), ymax)
                yplot(JAB+j) = yplot(j)
            enddo
            do j=1, JAB
                YXR = (YWB(j) + YWA(j))/AB
                YXL = (YWB(j) - YWA(j))/AB
                xplot(JAB+j)   = dble(JX) + 0.5*half_wid*(1. + min( 1.d0, YXR))
                xplot(JAB+1-j) = dble(JX) + 0.5*half_wid*(1. + max(-1.d0, YXL))
                yplot(JAB+1-j) = yplot(JAB+j)
            enddo
        endif
        jc = 31  ! non-marked profiles (shadow color)
        do IS=1, 5
            if (PRMARK(JTIM) == IS) then
                if (jk(IS) /= 0) EXIT
                jk(IS) = 1
                jc = IS + 1
            endif
        enddo

        STYL = (jc - 1)*MARK
        if (PRMARK(JTIM) == 0) STYL=0
        if (KPRI == 1 .or. KPRI == 2) then
            write(STRI, '(1A6, 1A6, 1A1)')'Plot "', NAMEP(JTIM), '"'
            j = len_trim_tab(STRI)
            call pscom(STRI, j)
        endif
        call setColor(jc)
        call plot_curve(NP1, STYL, xplot, yplot)
        if (PRMARK(JTIM) /= 0) then
            JPOS = 3*astra_gui%dxlet + PRMARK(JTIM)*plot_area%canvas_width/6.5 + plot_area%canvas_width*(NP - 1)
            axis_label(2: 5) = NAMEP(JTIM)
            if (STYL >= 7 .and. STYL <= 12) then
                axis_label(1: 1) = char(symbols(STYL - 6))
            else
                axis_label(1: 1) = ' '
            endif
            call textvm(JPOS, jnl, axis_label, 5)
        endif
    enddo

enddo read_loop
   
38 write(*, *)' >>> ERROR >>> Cannot read profile file'
stop

39 close(3)

return
end subroutine SMODE5

!---------------------------------------------------------------------
subroutine ADDMOD(NCHW, NCHM, NCHL)
! Write model & model.log records in the header of a post-view file
!   Both are preceded by one line 32*"^" 
!             and two such lines are written after model.log
! Unit NCHW (profile data file) must be open
! Unit NCHM (model) must be open
! Unit NCHL (model.log) must be open if nonzero
! Note:   1) NCHL =/= NCHM; 2) Empty lines are skipped.

use io_mod, only: null_ch
use char_manip, only: len_trim_tab
use debugger, only: markloc

implicit none

integer, intent(in) :: NCHW, NCHM, NCHL

integer :: NCH, NCHI, LSTRI, j, LSTR, ierr, ios
character(len=80) :: INCNAM
character(len=133) :: STRI, STR

call markloc('ADDMOD')

write(NCHW) "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
NCH = NCHM
NCHI = 0

read_loop: do

    write(STR, '(133A1)')(' ', j=1, 132), null_ch
    read(NCH, FMT='(1A132)', iostat=ios) STR

    if (ios > 0) then
        write(*, *) '>>> Error reading model file'
        EXIT read_loop
    else if (ios < 0) then
        write(NCHW) char(32), "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
        if (NCHL == 0 .or. NCHL == NCH) then
            write(NCHW) char(32), "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
            EXIT read_loop
        endif
        NCH = NCHL
        CYCLE read_loop
    endif

    LSTR = len_trim_tab(STR)
! '#' - pre-processor command
    if (STR(1: 8) /= "#include") then
        write(NCHW) char(LSTR), STR(1: LSTR)
        CYCLE read_loop
    endif

    LSTRI = index(STR(1:), ';') - 1
    if (LSTRI == -1) LSTRI = LSTR
    INCNAM = STR(9: LSTRI)

    NCHI = NCHM + 1
    open(NCHI, file="equ/" // TRIM(INCNAM), iostat=ios)
    if (ios /= 0) then
        write(*, *) '>>> Cannot open include file ', &
                '"equ/' // TRIM(INCNAM), '". Saving model failed'
        write(*, *)
        NCHI = 0
        if (LSTRI < LSTR) write(NCHW) char(LSTR - LSTRI), STR(LSTRI+1: LSTR)
        CYCLE read_loop
    endif

    do
        read(NCHI, FMT='(1A132)', iostat=ios) STRI
        if (ios > 0) then
            write(*, *) '>>> Cannot read include file ', &
                '"equ/' // TRIM(INCNAM), '". Saving model failed'
        endif
        if (ios /= 0) EXIT
        if (STRI(1: 8) == "#include") then
            write(*, *) '>>> Error in include file: "', INCNAM(1:j), '"'
            write(*, *) '    Recurrent "#include" is not allowed. Saving model failed'
            write(*, *)
            NCHI = 0
            if (LSTRI < LSTR) write(NCHW) char(LSTR - LSTRI), STR(LSTRI+1: LSTR)
            CYCLE read_loop
        endif
        ierr = len_trim_tab(STRI)
        write(NCHW) char(ierr), STRI(1:ierr)
    enddo

    close(NCHI)
    NCHI = 0
    if (LSTRI < LSTR) write(NCHW) char(LSTR - LSTRI), STR(LSTRI+1: LSTR)

enddo read_loop

return
end subroutine ADDMOD

!---------------------------------------------------------------------
subroutine STUFF(NUNIT, ARRAY, NA, BITARR)

use debugger, only: markloc

implicit none

integer, intent(in) :: NUNIT, NA
integer*2, intent(out) :: BITARR(NA)
double precision, intent(in) :: ARRAY(NA)

integer :: JJ
double precision :: YUP, YDN, YSC

call markloc('STUFF')

YUP = MAXVAL(ARRAY(1: NA))
YDN = MINVAL(ARRAY(1: NA))
YSC = YUP - YDN
if (YSC /= 0.) then
    do JJ=1, NA
        YUP = (ARRAY(JJ) - YDN)/YSC
        BITARR(JJ) = 65535*YUP - 32768
    enddo
endif
write(NUNIT) YSC, YDN, BITARR

return
end subroutine STUFF

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

use io_mod, only: null_ch

implicit none

integer, intent(in) :: arr_size, id
double precision, intent(in), dimension(arr_size) :: array_in
character(len=6), intent(in), dimension(arr_size) :: var_names

integer :: nameLength, editable=1
character(len=70), dimension(10), parameter :: titles = (/ &
    'Variable control', 'Constant control', 'Times & Grids', 'Scale control', &
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
