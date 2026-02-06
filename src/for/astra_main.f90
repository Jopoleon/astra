program astra

! Copyright (C) 2026 Institut fuer Plasmaphysik - Boltzmannstrasse 2, 85748 Garching (Germany)
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation;
! version 2.1 of the License

use outcmn_inc, only: astra_gui, astra_gui_ref, outcmn_init
use io_mod, only: TASK, io_init
use cpu_usage, only: cpu_start, wall_start, cpu_report
use const_inc, only: IPART, const_init, &
    TIME, TSTART, TEND, DPOUT, TAU, ATREQ, IFBEY, NITOT
use status_inc, only: status_init, defarr
use debugger, only: astra_stop, markloc
use ext_bnd, only: use_ext_bnd
use transport2fbe, only: transport2fbe_init
use json_vars, only: read_metadata
use json_write, only: write_json, write_jsonx
use exp_data, only: read_input

implicit none

!-------------------------------------------
! Find self-consistent initial configuration
!-------------------------------------------

integer :: j, jj, IM, ios, XSC0, XSC, jt1, jt2, jt3, jt_req, jkey, ierr, jt_out=0, rate
character(len=132) :: STRI
integer, external :: IFKEY, IFTREQ

save jt_out

!-------------------- Initial settings --------------------------------|

call CPU_TIME(cpu_start)
call SYSTEM_CLOCK(wall_start, rate)

call read_metadata
call associate_pointers

call const_init
call status_init

call ininam
call io_init
call read_input
call write_jsonx

use_ext_bnd = 0
IPART = 1   ! Mark initial iteration section

!--------------------
! ASTRA graphic frame
!--------------------

if (TASK(1: 3) /= 'BGD') then
    call outcmn_init
    call initMainWindow
endif

call SETARX(1)
call INIVAR
call SETVAR
call DETVAR
call EQGUESS
call INIVAR

call transport2fbe_init

jt_req = 0
do while (jt_req == 0) ! Till convergence (jt_req /= 0). Max #iterations is set in IFTREQ (for/defarr.f90)

    if (TASK(1:3) /= 'BGD') jkey = IFKEY(256)
    call INTVAR      ! Set exp scalars
    call DETVAR
    call DEFARR
    call SETARX(1)   ! Set X-data w/o time interpolation
    call INIVAR
    call markloc("init")
    NITOT = NITOT + 1

    call INIT_CONVERGE_STEP
    call markloc("init done")

    IFBEY = 0. ! no fbe possible here
    call METRIC
    jt_req = IFTREQ(ATREQ)     ! ++ITREQ; Convergence check; 1 - yes
enddo

if (TASK(1:3) /= 'BGD') then
    STRI(1:16) = ' ' ! Erase iteration number, iterations label top right
    call textvm(astra_gui%width-18*astra_gui_ref%dxlet, 2, STRI(1:16), 16)
    call textvm(astra_gui%width-17*astra_gui_ref%dxlet, astra_gui_ref%dylet + 1, STRI(1:14), 14)
endif

!---------------
! Time step loop
!---------------

do while (TIME - TEND + 1.E-8 < DPOUT + TAU)
    if ((TIME - TSTART + 1.E-8)/DPOUT >= jt_out) then
        call write_json
        jt_out = jt_out + 1
    endif
    call STEPUP
enddo

call CPU_report('>>> ASTRA normal exit >>>')
call astra_stop

end program astra
