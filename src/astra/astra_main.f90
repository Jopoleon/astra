program astra

! Copyright (C) 2026 Institut fuer Plasmaphysik - Boltzmannstrasse 2, 85748 Garching (Germany)
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation;
! version 2.1 of the License

use json_module, only: json_file
use graph_utils, only: astra_gui, astra_gui_ref, gui_init
use io_mod, only: TASK, io_init, MACHINE, awd, restart
use cpu_usage, only: cpu_start, wall_start, cpu_report
use const_inc, only: IPART, const_init, RTOR, UPDWN, SHIFT, PSIAX, PSIBO, &
    TIME, TINIT, TSTART, TEND, DPOUT, TAU, ATREQ, IFBEY, NITOT, NA1, &
    constValues, varValues, varxValues, internValues, intern2Values
use status_inc, only: status_init, defarr, profiles, profiles_x
use debugger, only: astra_stop, markloc
use ext_bnd, only: use_ext_bnd
use transport2fbe, only: transport2fbe_init
use json_vars, only: read_metadata, n_intern
use json_rw, only: json_load, write_json, read_scalar_block, read_array_block
use read_input, only: readInput, raw_cCoil
use plasma_state, only: plasma_up
use auxiliary, only: IFTREQ
use set_x_data, only: set_x_scalars, set_x_arrays, astra_assignments
use metrics, only: eqguess, metric

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
integer, external :: IFKEY

!-------------------- Initial settings --------------------------------|

call CPU_TIME(cpu_start)
call SYSTEM_CLOCK(wall_start, rate)

call read_metadata
call associate_pointers

call const_init
call status_init

call ininam
call io_init
plasma_up = 1  ! plasma is up by default, can be set to 0 for breakdown by the user in a user-defined sbr called with "<"
call readInput
call astra_assignments ! ASTRA default assignments

use_ext_bnd = 0
IPART = 1   ! Mark initial iteration section

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

if (restart > 0) then
    print*, 'astra_main1, TIME=', TIME
    write(f_json, '(2A, i0, A)') TRIM(awd), '/ncdf_out/aug34954fluxes-', restart, '.json'
    call json_load(f_json, fjson)
    call read_scalar_block(fjson, "constants", constValues)
    call read_scalar_block(fjson, "variables", varValues)
    call read_scalar_block(fjson, "variables_x", varxValues)
    call read_scalar_block(fjson, "internal", internVal)
    call read_scalar_block(fjson, "intern2", intern2Values)
    internValues(1: n_intern) = internVal(1: n_intern)
    write(*, '(A, 5f8.4)') 'astra_main2, TIME=', TIME
    call read_array_block(fjson, "profiles", profs)
    call read_array_block(fjson, "profiles_x", profs_x)
    profiles(  1:NA1, :) = profs
    profiles_x(1:NA1, :) = profs_x
endif

gui_on = (TASK(1: 3) /= 'BGD')
if (gui_on) then
    call gui_init
endif
jt_req = 0
do while (jt_req == 0) ! Till convergence (jt_req /= 0). Max #iterations is set in IFTREQ (for/defarr.f90)

    if (gui_on) jkey = IFKEY(256)
    call set_x_scalars   ! Set exp scalars
    call DETVAR
    call DEFARR
    call set_x_arrays(1)   ! Set X-data w/o time interpolation
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
    call STEPUP
enddo

call CPU_report('>>> ASTRA normal exit >>>')
call astra_stop

end program astra
