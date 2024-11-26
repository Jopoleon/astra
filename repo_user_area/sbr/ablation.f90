subroutine ABLATION(trace, pel_prof)

!----------------------------------------------------------------------|
! Pellet ablation routine
!----------------------------------------------------------------------|

use parameter_inc, only: NRD
use const_inc, only: TIME, NA1
use outcmn_inc, only: AWD, nml_file
use parse_utils, only: ufheader, ufrd
use status_inc, only: XRHO

implicit none

integer, parameter :: unit=1032

double precision, intent(out) :: trace, pel_prof(NRD)

double precision :: tpel_next=0.d0
character(len=30) :: rholbl
character(len=120) :: as_nml, rho_abl_file, time_abl_file
integer :: ios=0, jt, jr, nt_u, nx_u, ndim_u
double precision :: dt, rho_abl, mass
double precision, dimension(:), allocatable :: t_u, x_u, var_u

NAMELIST / pellet / rho_abl_file, time_abl_file, mass

! Read pellet data

as_nml = TRIM(awd) // TRIM(nml_file)
write(*, *) 'Reading namelist ', TRIM(as_nml)

open(57, FILE=TRIM(as_nml), delim='apostrophe')
read(57, nml=pellet, iostat=ios)
close(57)

call ufheader(TRIM(rho_abl_file), ndim_u, nt_u, nx_u, rholbl)
allocate(t_u(nt_u))
allocate(x_u(nx_u))
allocate(var_u(nt_u*nx_u))
call ufrd(TRIM(rho_abl_file), ndim_u, nt_u, nx_u, t_u, x_u, var_u)
call uf1dr(time_abl_file, TIME, dt)

dt = dt*10
do jt=1, nt_u
    tpel_next = t_u(jt)
    trace = 0.d0
    pel_prof(1:NA1) = 1.d0 ! Gets to zero due to trace mutliplier
    if (tpel_next > TIME) EXIT
    if (tpel_next + dt >= TIME) then
        trace = mass/dt
        rho_abl = var_u(jt)
        do jr=1, NA1
            if (XRHO(jr) > rho_abl) then
                pel_prof(jr) = 1.d0 - XRHO(jr)
            else
                pel_prof(jr) = 0.d0
            endif
        enddo
        EXIT
    endif
enddo


return
end subroutine ablation
