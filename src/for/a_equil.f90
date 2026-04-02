subroutine A_EQUIL( &
! Input
    equil_in, equil_solver, &
    nr_equ, n_theta, iter_step, &
    ncoils, ccoils, vcoils, tau_step, time_a, &
    ipsibcf, key_no_refits, &
    icircq, ipctrl, &
    iter_itreq, ifbey, inume_3, &
! Output:
    key_start, PSIEXT, PSPLEX, keyplc, equil_out)

use imas_ids, only: type_equilibrium
use parameters_a2equil, only: type_parameters, fix_adapgrid, GP, GP2, s_fazt
use const_inc, only : rtor,shift, updwn
use feqis_circuit, only: psib_ext_feqis

use io_mod, only: MACHINE, nml_file

implicit none

integer, intent(in) :: equil_solver, nr_equ, n_theta, iter_step, ncoils, &
    ipsibcf, key_no_refits, icircq, ipctrl, iter_itreq, ifbey, inume_3
double precision, intent(in) :: tau_step, time_a
double precision, intent(in), dimension(ncoils) :: ccoils, vcoils
type(type_equilibrium), intent(in) :: equil_in

integer, intent(out) :: key_start, keyplc
double precision, intent(out) :: PSIEXT, PSPLEX
type(type_equilibrium), intent(out) :: equil_out

logical :: file_existence
integer :: nstep, i, key_equil, nrp, &
    toric_fourc, toric_file, strahl_file, strahl_fourc, &
    write_coils_diagn, key_plcs, kprs, k_grids, &
    kprs2, fixadapgrid

double precision :: dampfacpsplex, psplexold, epsros, enelss, k_filessss, ipl
double precision, dimension(ncoils) :: t_currents, ucoils
double precision :: psplexavg, psplexavgexp
character(len=120) :: fname

type(type_parameters) :: parameters_equil

data psplexold /2./
save toric_fourc, toric_file
save strahl_file, strahl_fourc, write_coils_diagn
save kprs, k_grids, epsros, enelss, key_plcs, k_filessss
save psplexavg, psplexold, kprs2, psplexavgexp

namelist / spider / kprs, k_grids, epsros, enelss, key_plcs, &
    toric_fourc, toric_file, strahl_file, strahl_fourc, write_coils_diagn, &
    k_filessss, psplexavg, psplexavgexp, fixadapgrid

!for PBE , use p and cu, key_equil=key_dmf=-10, nstep = 0 only at first iteration
key_equil = 0
nrp = 256

nstep = max(0, ifbey-1)

! grids
parameters_equil%dt     = tau_step
parameters_equil%time   = time_a
parameters_equil%neql   = nr_equ
parameters_equil%ntheta = n_theta + 2
parameters_equil%prename = trim(parameters_equil%prename) // trim(MACHINE) // '/'

!defaults
if (nstep == 0) then
    kprs = 0
    k_grids = 0
    epsros = 1.e-9
    enelss = 1.e-9
    key_plcs = 1
    toric_fourc = 4
    toric_file = 1
    strahl_file = 1
    strahl_fourc = 3
    write_coils_diagn = 0
    k_filessss = 0.
    psplexavg = 0.
    psplexavgexp = 0.
    fix_adapgrid = 0
    INQUIRE( FILE=TRIM(nml_file), EXIST=file_existence)
    if (file_existence) then
        open(53, FILE=TRIM(nml_file))
        read(53, nml=spider)
        close(53)
        kprs2 = kprs
    endif
!    fix_adapgrid = fixadapgrid
endif

if (s_fazt == 0) then
    kprs = kprs2
else
    kprs = -2
endif
k_grids = 0

parameters_equil%kpr     = kprs
parameters_equil%k_grid  = k_grids
parameters_equil%epsro   = epsros
parameters_equil%enels   = enelss
parameters_equil%key_plc = key_plcs
parameters_equil%key_dmf = 0

parameters_equil%key_plc = 1    !force plc = 1 if current diffusion is solved
keyplc = parameters_equil%key_plc
parameters_equil%key_out = 0

if (iter_step == 1) parameters_equil%k_fixfree = 0 !astra initialization, no fbe
if (parameters_equil%k_fixfree == 0) nstep = 0  !no fbe, nstep=0
parameters_equil%no_circuit_eq = 0
if (icircq == 0) parameters_equil%no_circuit_eq = 1    !no circuit equations, only static fbe

if (parameters_equil%k_fixfree == 1) then
    SELECT CASE(ipctrl)
    CASE(1, -4, -5)  !controller, refit currents, coil.dat untouched
        parameters_equil%key_start = 1
    CASE(-3: 0)  !no controller, coil.dat untouched
        parameters_equil%key_start = 0
    END SELECT

    key_start = parameters_equil%key_start
    parameters_equil%key_out = 0  ! keep this and use spidupdate call instead

    if (nstep >= 1) parameters_equil%key_start = 0 !fbe with circuit equations, no refit
endif

!Coil currents
if (ncoils > 0) then
    if (parameters_equil%k_fixfree == 1) then
        SELECT CASE(ipctrl)
        CASE(-5, -3, -2)
!            call coil2spider(ccoils, ncoils, parameters_equil)   !Write coil currents from CCOIL in astra to   coil.dat file only for fbe without controller (otherwise CCOIL is reserved for target coil currents and coil.dat is written elsewhere)
        END SELECT
    endif
endif

!use refits currents in coil.dat, only for nitreq >1
if (key_no_refits == 1) then
    if (key_start == 1 .and. iter_itreq > 0) then
        fname = trim(parameters_equil%prename) // 'tcurrs.wr'
        open(1, file=TRIM(fname))
        do i=1, ncoils
            read(1, *) t_currents(i)
        enddo
        close(1)
        write(*, *) 'rewriting coil.dat with new fitted currents'
!        call coil2spider(t_currents*1.e3, ncoils, parameters_equil)  !Write coil currents in coil.dat when key_start inside iterations NITREQ
    endif
endif

if (parameters_equil%k_fixfree == 1) then
    ucoils(1:ncoils)  = vcoils(1:ncoils)
    parameters_equil%nstep = nstep
endif

if (equil_solver == 101) then
    call feqis_main(ncoils, ucoils, parameters_equil, 1, equil_in, equil_out)
else
    call spider_run(ncoils, ucoils, equil_in, equil_out, parameters_equil)
endif
if (ipsibcf /= 0) parameters_equil%key_psibcf = 1

!output from equil_out structure

!psifb = psifb_in
dampfacpsplex = 0.

if (parameters_equil%k_fixfree == 1) then
    ipl = 1.e-6*equil_in%global_param%i_plasma
    if (parameters_equil%k_grid == 0) then
        if (equil_solver == 101) then
            PSIEXT = psib_ext_feqis()
        else
            call psib_ext(PSIEXT)
        endif
    endif
    if (ipsibcf >= 0) then     ! case with PSI_B and dPSI_B implicit
        PSIEXT = -GP2*PSIEXT
        PSPLEX = equil_out%global_param%psplex
        PSPLEX = ((psplexavg*ipl**psplexavgexp)/tau_step*psplexold + PSPLEX) / &
            (1. + (psplexavg*ipl**psplexavgexp)/tau_step)
        psplexold = PSPLEX
    else
        if (parameters_equil%k_grid == 0) then
            if (equil_solver==101) then
                PSIEXT = psib_ext_feqis()
            else
                call psib_ext(PSIEXT)
            endif
        endif
        PSIEXT = -GP2*PSIEXT
        PSPLEX = (dampfacpsplex*PSPLEX + equil_out%global_param%psplex)/(1. + dampfacpsplex)
        PSPLEX = ((psplexavg*ipl**psplexavgexp)/tau_step*psplexold + PSPLEX) / &
            (1. + (psplexavg*ipl**psplexavgexp)/tau_step)
        psplexold = PSPLEX
    endif
endif

!PSPLEX in FEQIS is the Lext already. In SPIDER NOT.
PSPLEX = equil_out%global_param%psplex
PSIEXT = -GP2*equil_out%global_param%psiext

return
end subroutine A_equil

!---------------------------------------------------------------------
subroutine A_equil_2(ncoils, ifbey, time_a, tau_step, vcoils, eq_solver, iplas_vac)

use imas_ids, only: type_equilibrium
use parameters_a2equil, only: type_parameters
use io_mod, only: MACHINE

implicit none

integer, intent(in) :: ifbey, eq_solver, ncoils
double precision, intent(in) :: tau_step, time_a, iplas_vac
double precision, dimension(ncoils), intent(in) :: vcoils

integer :: nstep, key_equil
double precision, dimension(ncoils) :: ucoils

type(type_parameters) :: parameters_equil
type(type_equilibrium) :: equil_in, equil_out

key_equil = 0
nstep = max(0, ifbey - 1)

equil_in%global_param%i_plasma = iplas_vac*1e6   !itm is in A

parameters_equil%dt      = tau_step
parameters_equil%time    = time_a
parameters_equil%prename = 'exp/equ/'//trim(MACHINE)//'/'
parameters_equil%kpr     = -2
parameters_equil%k_grid  = 1
parameters_equil%epsro   = 1.d-9
parameters_equil%enels   = 1.d-9
parameters_equil%key_plc = 1

parameters_equil%key_out   = 0
parameters_equil%k_fixfree = 1
parameters_equil%key_start = 0    !controller, refit currents, coil.dat untouched

ucoils(1:ncoils) = vcoils(1:ncoils)

parameters_equil%nstep = nstep

if (eq_solver == 101) then
    call feqis_main(ncoils, ucoils, parameters_equil, 0, equil_in, equil_out)
endif

return
end subroutine a_equil_2
