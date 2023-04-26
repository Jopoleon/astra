subroutine A_SPIDER( &
! Input
    equil_in, equil_solver, &
    nr_equ, n_theta, iter_step, &
    ncoils, ccoils, vcoils, tau_step, time_a, &
    ipsibcf, key_no_refits, &
    icircq, ipctrl, &
    iter_itreq, machine_name, ifbey, inume_3, &
! Output:
    key_start, PSIEXT, PSPLEX, keyplc, equil_out)

use imas_ids, only: type_equilibrium
use fenix_params, only: s_adapt, s_fazt
use parameters_a2spider, only: type_parameters, fix_adapgrid, GP, GP2

implicit none

integer, intent(in) :: equil_solver, nr_equ, n_theta, iter_step, ncoils, &
    ipsibcf, key_no_refits, icircq, ipctrl, iter_itreq, ifbey, inume_3
real*8, intent(in) :: tau_step, time_a
real*8, intent(in), dimension(ncoils) :: ccoils, vcoils
character(len=4), intent(in) :: machine_name
type(type_equilibrium), intent(in) :: equil_in

integer, intent(out) :: key_start, keyplc
real*8, intent(out) :: PSIEXT, PSPLEX
type(type_equilibrium), intent(out) :: equil_out

logical :: file_existence 
integer :: nstep, i, j, key_equil, nrp, &
    toric_fourc, toric_file, strahl_file, strahl_fourc, &
    write_coils_diagn, key_plcs, kprs, k_grids, &
    kprs2, fixadapgrid

double precision :: psplexavg, psplexavgexp

real*8 :: dampfacpsplex, psplexold, epsros, enelss, k_filessss, ipl
real*8, dimension(ncoils) :: t_currents, ucoils

character(len=80) :: fname

type(type_parameters) :: parameters_spider

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
s_adapt = 0
s_fazt = 0

nstep = max(0, ifbey-1)

! grids
parameters_spider%dt    = tau_step
parameters_spider%time  = time_a
parameters_spider%neql  = nr_equ
parameters_spider%nteta = n_theta + 2
parameters_spider%prename = trim(parameters_spider%prename) // trim(machine_name) // '/'

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
    fname = trim(parameters_spider%prename) // 'namelist_astra.txt'
    INQUIRE( FILE=trim(fname), EXIST=file_existence) 
    if (file_existence) then
        open(53, FILE=fname)
        read(53, nml=spider)
        close(53)
        kprs2 = kprs
    endif
    fix_adapgrid = fixadapgrid
endif

if (s_fazt == 0) then
    kprs = kprs2
else
    kprs = -2
endif
if (s_adapt == 0) then
    k_grids = 0
else
    k_grids = 1
endif

parameters_spider%kpr     = kprs
parameters_spider%k_grid  = k_grids
parameters_spider%epsro   = epsros
parameters_spider%enels   = enelss
parameters_spider%key_plc = key_plcs      
parameters_spider%key_dmf = 0

if (inume_3 /= -1) parameters_spider%key_plc = 1    !force plc = 1 if current diffusion is solved
keyplc = parameters_spider%key_plc
parameters_spider%key_out = 0

if (icircq == 0) nstep = 0    !no circuit equations, only static fbe
if (iter_step == 1) parameters_spider%k_fixfree = 0 !astra initialization, no fbe
if (parameters_spider%k_fixfree == 0) nstep = 0  !no fbe, nstep=0

if (parameters_spider%k_fixfree == 1) then
    SELECT CASE(ipctrl)
    CASE(1, -4, -5)  !controller, refit currents, coil.dat untouched
        parameters_spider%key_start = 1
    CASE(-3: 0)  !no controller, coil.dat untouched
        parameters_spider%key_start = 0
    END SELECT

    key_start = parameters_spider%key_start
    parameters_spider%key_out = 0  ! keep this and use spidupdate call instead

    if (nstep >= 1) parameters_spider%key_start = 0 !fbe with circuit equations, no refit
endif      

!Coil currents
if (ncoils > 0) then
    if (parameters_spider%k_fixfree == 1) then
        SELECT CASE(ipctrl)
        CASE(-5, -3, -2)
            call coil2spider(ccoils, ncoils, parameters_spider)   !Write coil currents from CCOIL in astra to   coil.dat file only for fbe without controller (otherwise CCOIL is reserved for target coil currents and coil.dat is written elsewhere)
        END SELECT
    endif
endif

!use refits currents in coil.dat, only for nitreq >1
if (key_no_refits == 1) then
    if (key_start == 1 .and. iter_itreq > 0) then
        fname = trim(parameters_spider%prename) // 'tcurrs.wr' 
        open(1, file=TRIM(fname))
        do i=1, ncoils
            read(1, *) t_currents(i)  
        enddo
        close(1)
        write(*, *) 'rewriting coil.dat with new fitted currents'
        call coil2spider(t_currents*1.e3, ncoils, parameters_spider)  !Write coil currents in coil.dat when key_start inside iterations NITREQ
    endif
endif 

if (parameters_spider%k_fixfree == 1) then
    if (machine_name == 'aug_'.or.machine_name == 'aug '.or. machine_name == 'aug5') then
        ucoils(1)  = vcoils(1) - vcoils(7)
        ucoils(2)  = vcoils(7) - vcoils(8)
        ucoils(3)  = vcoils(8)
        ucoils(4)  = vcoils(6)
        ucoils(5)  = vcoils(5)
        ucoils(6)  = vcoils(4)
        ucoils(7)  = vcoils(2)
        ucoils(8)  = vcoils(3)
        ucoils(9)  = vcoils(9)
        ucoils(10) = vcoils(10)
        ucoils(11:12) = 0.
    elseif (machine_name == 'dem_') then !DEMO free boundary, to recheck
        ucoils(1:11)  = vcoils(1:11)
    endif
    parameters_spider%nstep = nstep
endif

if (equil_solver == 101) then
    call feqis_main(equil_in, equil_out)
else
    call spider_run(ncoils, ucoils, equil_in, equil_out, parameters_spider)     
endif

if (ipsibcf /= 0) parameters_spider%key_psibcf = 1

if (parameters_spider%k_fixfree == 1 .and. nstep == 0) then
    open(32, file=TRIM(parameters_spider%prename) // 'data.dat')
    read(32, *) j
    write(*, *) j
    if (j == 0) stop
    close(32)
endif

!output from equil_out structure 

!psifb = psifb_in
dampfacpsplex = 0.

if (parameters_spider%k_fixfree == 1) then
    ipl = 1.e-6*equil_in%global_param%i_plasma
    if (parameters_spider%k_grid == 0) then
        call psib_ext(PSIEXT)
    endif
    if (parameters_spider%k_grid == 1) then
        call f_psib_ext(PSIEXT)
    endif
    if (ipsibcf >= 0) then     ! case with PSI_B and dPSI_B implicit 
        PSIEXT = -GP2*PSIEXT
        PSPLEX = equil_out%global_param%psplex
        PSPLEX = ((psplexavg*ipl**psplexavgexp)/tau_step*psplexold + PSPLEX) / &
            (1. + (psplexavg*ipl**psplexavgexp)/tau_step)
        psplexold = PSPLEX
    else
        if (parameters_spider%k_grid == 0) then
            call psib_ext(PSIEXT)
        endif
        if (parameters_spider%k_grid == 1) then
            call f_psib_ext(PSIEXT)
        endif
        PSIEXT = -GP2*PSIEXT
        PSPLEX = (dampfacpsplex*PSPLEX + equil_out%global_param%psplex)/(1. + dampfacpsplex)
        PSPLEX = ((psplexavg*ipl**psplexavgexp)/tau_step*psplexold + PSPLEX) / &
            (1. + (psplexavg*ipl**psplexavgexp)/tau_step)
        psplexold = PSPLEX
    endif
endif

return
end subroutine A_SPIDER
