program astra

! Copyright (C) 2026 Institut fuer Plasmaphysik - Boltzmannstrasse 2, 85748 Garching (Germany)
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation;
! version 2.1 of the License

use graph_utils, only: astra_gui, astra_gui_ref, gui_init
use cpu_usage, only: cpu_init, cpu_start, wall_start, cpu_report
use scalars, only: IPART, scalars_init, RTOR, UPDWN, SHIFT, PSIAX, PSIBO, &
    TIME, TSTART, TEND, TPAUSE, DPOUT, TAU, ATREQ, IFBEY, NITOT, IPEQL
use status, only: status_init, defarr, setvar
use debugger, only: astra_stop, markloc
use transport2fbe, only: transport2fbe_init
use json_vars, only: read_metadata
use json_rw, only: read_ajson, write_ajson
use read_input, only: readInput, raw_cCoil, TASK, MACHINE, &
    restart, tend_nml, tpause_nml
use auxiliary, only: IFTREQ
use set_x_data, only: set_x_scalars, set_x_arrays, astra_assignments
use metrics, only: eqguess, metric, CCOIL, VCOIL
use gui_interaction, only: if_key

implicit none

!-------------------------------------------
! Find self-consistent initial configuration
!-------------------------------------------

logical :: gui_on
integer :: jt_req, jkey, jt_out, rate
double precision :: t_stop
character(len=16) :: str_iterations

!-------------------- Initial settings --------------------------------|

call CPU_TIME(cpu_start)
call SYSTEM_CLOCK(wall_start, rate)

call read_metadata      ! reads variables/profiles lists from astra_variables.json
call associate_pointers ! associates vars to pointers consistently with lists in astra_variables.json
call cpu_init           ! starting walltime diagnostic for each sbr
gui_on = (TASK(1: 3) /= 'BGD') ! Graphic window yes/no

call scalars_init ! Fallback default values for scalars
call ininam       ! Sets DTEQ, DTNAME and plot labels (from equ file); call after scalars_init!
call set_graph_names

call readInput    ! "restart" is set inside readInput
allocate(CCOIL(raw_cCoil%ncoils), VCOIL(raw_cCoil%ncoils))
if (gui_on) then
    call gui_init ! Start GUI
endif

if (restart > 0) then ! Initial condition from output json file
    call read_ajson(restart) ! Read the desired json file
    call set_x_arrays(2) ! To plot also exp raw data
    if (IPEQL == 5) then ! FEQIS
        call transport2fbe_init(TAU, TSTART, RTOR, UPDWN, SHIFT, PSIAX, PSIBO, &
            MACHINE, raw_cCoil%ncoils, raw_cCoil%current) ! Transfer ASTRA pars to FEQIS
    endif
    tend = tend_nml ! TEND is read from json, TEND_NML from command line "-e TBEG_NML"
    tpause = tpause_nml
    if (gui_on) jkey = if_key(0) ! Plot right now
else ! Iterations for initial convergence
    call status_init       ! Fallback default values for profiles
    call astra_assignments ! ASTRA default assignments
    IPART = 1              ! Mark initial iteration section
    call set_x_arrays(1)   ! Set TEX, CAR1X at current time from raw_*
    call INIVAR            ! User-defined (equ) initialisation of NE, TE, TI, UPAR, FJ, MU/CU
    call SETVAR            ! Set fallback NI, ZEF and a sanity check if ABC > AB
    call DETVAR
    call eqguess
    call INIVAR
    if (IPEQL == 5) then ! FEQIS
        call transport2fbe_init(TAU, TSTART, RTOR, UPDWN, SHIFT, PSIAX, PSIBO, &
            MACHINE, raw_cCoil%ncoils, raw_cCoil%current) ! Transfer ASTRA pars to FEQIS
    endif

    jt_req = 0
    do while (jt_req == 0) ! Till convergence (jt_req /= 0). Max #iterations is set in IFTREQ (status:defarr)
        if (gui_on) jkey = if_key(256) ! Plot right now
        call set_x_scalars    ! Set exp scalars
        call DETVAR
        call DEFARR
        call set_x_arrays(1)  ! Set X-data w/o time interpolation
        call INIVAR
        call markloc("init")
        NITOT = NITOT + 1
        call INIT_CONVERGE_STEP
        call markloc("init done")
        IFBEY = 0 ! no fbe possible here
        call METRIC
        jt_req = IFTREQ(ATREQ)     ! ++ITREQ; Convergence check; 1 - yes
    enddo
endif

if (gui_on) then
    str_iterations(1:16) = ' ' ! Erase iteration number, iterations label top right
    call textvm(astra_gui%width-18*astra_gui_ref%dxlet, 2, str_iterations(1:16), 16)
    call textvm(astra_gui%width-17*astra_gui_ref%dxlet, astra_gui_ref%dylet + 1, str_iterations(1:14), 14)
endif

!---------------
! Time step loop
!---------------

jt_out = 0
t_stop = TEND + max(DPOUT, TAU)
do while (TIME < t_stop)
    if ((TIME - TSTART + 1.E-8)/DPOUT >= jt_out) then
        call write_ajson
        jt_out = jt_out + 1
    endif
    call STEPUP ! Time-dependent evolution
enddo

call CPU_report('>>> ASTRA normal exit >>>')
call astra_stop

end program astra
