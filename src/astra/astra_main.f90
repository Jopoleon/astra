program astra

! Copyright (C) 2026 Institut fuer Plasmaphysik - Boltzmannstrasse 2, 85748 Garching (Germany)
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation;
! version 2.1 of the License

use json_module, only: json_file
use graph_utils, only: astra_gui, astra_gui_ref, gui_init
use io_mod, only: TASK, io_init, MACHINE, awd, restart, exp_file, equ_file
use cpu_usage, only: cpu_init, cpu_start, wall_start, cpu_report
use scalars, only: IPART, const_init, RTOR, UPDWN, SHIFT, PSIAX, PSIBO, &
    TIME, TINIT, TSTART, TEND, DPOUT, TAU, ATREQ, IFBEY, NITOT, NA1, &
    constValues, varValues, varxValues, internValues, intern2Values
use status, only: status_init, defarr, profiles, profiles_x
use debugger, only: astra_stop, markloc
use transport2fbe, only: transport2fbe_init
use json_vars, only: read_metadata, n_intern
use json_rw, only: json_load, write_json, &
    read_scalar_block, read_array_block, read_equil
use read_input, only: readInput, raw_cCoil
use auxiliary, only: IFTREQ
use set_x_data, only: set_x_scalars, set_x_arrays, astra_assignments
use metrics, only: eqguess, metric, CCOIL, VCOIL
use gui_interaction, only: if_key

implicit none

!-------------------------------------------
! Find self-consistent initial configuration
!-------------------------------------------

logical :: gui_on
integer :: j, jj, IM, ios, XSC0, XSC, jt1, jt2, jt3, jt_req, jkey, ierr, jt_out, rate
double precision :: t_stop
double precision, dimension(:), allocatable :: internVal
double precision, dimension(:, :), allocatable :: profs, profs_x
character(len=132) :: STRI, f_json
type(json_file) :: fjson

!-------------------- Initial settings --------------------------------|

call CPU_TIME(cpu_start)
call SYSTEM_CLOCK(wall_start, rate)

call read_metadata
call associate_pointers

call cpu_init
call const_init
call status_init

call ininam
call io_init
call readInput
call astra_assignments ! ASTRA default assignments

IPART = 1   ! Mark initial iteration section

allocate(CCOIL(raw_cCoil%ncoils), VCOIL(raw_cCoil%ncoils))
CCOIL = 0.
VCOIL = 0.

!--------------------
! ASTRA graphic frame
!--------------------

call set_x_arrays(1)
call INIVAR
call SETVAR
call DETVAR
call eqguess
call INIVAR

call transport2fbe_init(TAU, TSTART, RTOR, UPDWN, SHIFT, PSIAX, PSIBO, MACHINE, &
    raw_cCoil%ncoils, raw_cCoil%current)

gui_on = (TASK(1: 3) /= 'BGD')
if (gui_on) then
    call gui_init
endif

if (restart > 0) then
    print*, 'astra_main1, TIME=', TIME
    write(f_json, '(5A, i0, A)') TRIM(awd), '/ncdf_out/', TRIM(exp_file), &
        TRIM(equ_file), '-', restart, '.json'
    call json_load(f_json, fjson)
    call read_scalar_block(fjson, "constants", constValues)
    call read_scalar_block(fjson, "variables", varValues)
    call read_scalar_block(fjson, "variables_x", varxValues)
    call read_scalar_block(fjson, "internal", internVal)
    call read_scalar_block(fjson, "intern2", intern2Values)
    internValues(1: n_intern) = internVal(1: n_intern)
    write(*, '(A, 2f8.4, i0)') 'astra_main2, TIME=', TIME, internValues(2), n_intern
    call read_array_block(fjson, "profiles", profs)
    call read_array_block(fjson, "profiles_x", profs_x)
    profiles(  1:NA1, :) = profs
    profiles_x(1:NA1, :) = profs_x
    call read_equil(fjson)
    call fjson%destroy()
else ! Iterations for initial convergence
    jt_req = 0
    do while (jt_req == 0) ! Till convergence (jt_req /= 0). Max #iterations is set in IFTREQ (status:defarr)
        if (gui_on) jkey = if_key(256)
        call set_x_scalars    ! Set exp scalars
        call DETVAR
        call DEFARR
        call set_x_arrays(1)  ! Set X-data w/o time interpolation
        call INIVAR
        call markloc("init")
        NITOT = NITOT + 1
        call INIT_CONVERGE_STEP
        call markloc("init done")
        IFBEY = 0. ! no fbe possible here
        call METRIC
        jt_req = IFTREQ(ATREQ)     ! ++ITREQ; Convergence check; 1 - yes
    enddo
 
    if (gui_on) then
        STRI(1:16) = ' ' ! Erase iteration number, iterations label top right
        call textvm(astra_gui%width-18*astra_gui_ref%dxlet, 2, STRI(1:16), 16)
        call textvm(astra_gui%width-17*astra_gui_ref%dxlet, astra_gui_ref%dylet + 1, STRI(1:14), 14)
    endif
endif

!---------------
! Time step loop
!---------------

jt_out = 0
t_stop = TEND + max(DPOUT, TAU)
do while (TIME < t_stop)
    if ((TIME - TSTART + 1.E-8)/DPOUT >= jt_out) then
        call write_json
        jt_out = jt_out + 1
    endif
    call STEPUP ! Time-dependent evolution
enddo

call CPU_report('>>> ASTRA normal exit >>>')
call astra_stop

end program astra
