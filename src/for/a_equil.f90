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
use fenix_params, only: s_adapt, s_fazt
use parameters_a2equil, only: type_parameters, fix_adapgrid, GP, GP2
use const_inc, only : rtor,shift,updwn
use feqis_tools, only: psib_ext_efff, get_zccurb_efff, find_demo_gaps_efff

use flight_sim_geometrics, only: geom1d
use outcmn_inc, only: MACHINE, nml_file
use astra2fbe

implicit none

integer, intent(in) :: equil_solver, nr_equ, n_theta, iter_step, ncoils, &
    ipsibcf, key_no_refits, icircq, ipctrl, iter_itreq, ifbey, inume_3
real*8, intent(in) :: tau_step, time_a
real*8, intent(in), dimension(ncoils) :: ccoils, vcoils
type(type_equilibrium), intent(in) :: equil_in

integer, intent(out) :: key_start, keyplc
real*8, intent(out) :: PSIEXT, PSPLEX
type(type_equilibrium), intent(out) :: equil_out

logical :: file_existence 
integer :: nstep, i, j, key_equil, nrp, nz, &
    toric_fourc, toric_file, strahl_file, strahl_fourc, &
    write_coils_diagn, key_plcs, kprs, k_grids, &
    kprs2, fixadapgrid, &
    jzmin, jzmax, jdemogaps, i_gaps

real*8 :: dampfacpsplex, psplexold, epsros, enelss, k_filessss, ipl
real*8, dimension(ncoils) :: t_currents, ucoils
double precision :: psplexavg, psplexavgexp, Rmag, Zmag, Rgeo, Zgeo, &
    rcurr, zcurr, rgeoc, zgeoc, ahorc, zsquad, psi_sep, psi_axis, &
    Rin, Raus, zoben, zunten, elong, &
    R_strike_in_aug, R_strike_out_aug, delr_oben, amin
double precision, dimension(100, 4) :: demo_gaps
double precision, dimension(n_theta) :: Rbnd, Zbnd
character(len=80) :: fname

type(type_parameters) :: parameters_equil

data jdemogaps /0/
data psplexold /2./
save toric_fourc, toric_file
save strahl_file, strahl_fourc, write_coils_diagn
save kprs, k_grids, epsros, enelss, key_plcs, k_filessss
save psplexavg, psplexold, kprs2, psplexavgexp, jdemogaps, i_gaps
save demo_gaps

namelist / spider / kprs, k_grids, epsros, enelss, key_plcs, &
    toric_fourc, toric_file, strahl_file, strahl_fourc, write_coils_diagn, &
    k_filessss, psplexavg, psplexavgexp, fixadapgrid

!for PBE , use p and cu, key_equil=key_dmf=-10, nstep = 0 only at first iteration
key_equil = 0
nrp = 256

nstep = max(0, ifbey-1)

! grids
parameters_equil%dt    = tau_step
parameters_equil%time  = time_a
parameters_equil%neql  = nr_equ
parameters_equil%nteta = n_theta + 2
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
    raxis_astra = rtor + shift
    zaxis_astra = updwn
    fname = trim(nml_file)
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

parameters_equil%kpr     = kprs
parameters_equil%k_grid  = k_grids
parameters_equil%epsro   = epsros
parameters_equil%enels   = enelss
parameters_equil%key_plc = key_plcs      
parameters_equil%key_dmf = 0

parameters_equil%key_plc = 1    !force plc = 1 if current diffusion is solved
keyplc = parameters_equil%key_plc
parameters_equil%key_out = 0

!if (icircq == 0) nstep = 0    !no circuit equations, only static fbe !obsolete
if (iter_step == 1) parameters_equil%k_fixfree = 0 !astra initialization, no fbe
if (parameters_equil%k_fixfree == 0) nstep = 0  !no fbe, nstep=0

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
            call psib_ext_efff(PSIEXT)
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
                call psib_ext_efff(PSIEXT)
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

! for any machine, geom1d(299) and geom1d(300) are respecetively li3 and betapol from equil

geom1d(299) = equil_out%global_param%li3
geom1d(300) = equil_out%global_param%betpol
do j=1, n_theta
    Rbnd(j) = equil_out%coord_sys%position%r(nr_equ, j)
    Zbnd(j) = equil_out%coord_sys%position%z(nr_equ, j)
enddo
Rmag  = equil_out%coord_sys%position%r(1, 1)
Zmag  = equil_out%coord_sys%position%z(1, 1)
elong = equil_out%profiles_1d%elongation(nr_equ)
nz    = equil_out%eqgeometry%rectgrid%npointsz

jzmin  = minloc(Zbnd, 1)
jzmax  = maxloc(Zbnd, 1)
Rin    = MINVAL(Rbnd)
Raus   = MAXVAL(Rbnd)
zoben  = Zbnd(jzmax)
zunten = Zbnd(jzmin)
Rgeo = 0.5*(Raus + Rin)
Zgeo = 0.5*(zoben + zunten)
amin = 0.5*(Raus - Rin)
delr_oben = (Rgeo - Rbnd(jzmax))/amin
if (equil_solver == 101) then
    call get_zccurb_efff(Rcurr, Zcurr, Zsquad, rgeoc, zgeoc, ahorc)
else if (equil_solver == 3) then
    call get_zccurb(Rcurr, Zcurr, Zsquad, rgeoc, zgeoc, ahorc)
endif

!below are valid for any machine bcs of prescribed boundary already has these values, but can be overwritten more below by specific machines
geom1d(51) = Zcurr ! For now Zsquad = Zcurr
geom1d(52) = zoben
geom1d(53) = rgeoc !should be rgeo, should go somewhere else
geom1d(54) = zgeoc !should be zgeo, should go somewhere else
geom1d(55) = ahorc !minor radius
geom1d(56) = elong
geom1d(57) = Rin
geom1d(58) = Raus
geom1d(59) = Rmag
geom1d(60) = Zmag
geom1d(70) = delr_oben
geom1d(78) = zunten
geom1d(80) = Rbnd(jzmin) ! Xpoint position R
geom1d(81) = Rcurr
geom1d(82) = Zcurr

if (parameters_equil%k_fixfree == 1) then
    if (MACHINE(1:3) == 'aug') then
        R_strike_in_aug  = 1.27
        R_strike_out_aug = 1.72
! inner strike point position
        i = minloc(abs(equil_out%eqgeometry%rectgrid%r2d - R_strike_in_aug), 1)  
        call find_in_vec_spid(nz, -gp2*equil_out%eqgeometry%rectgrid%psirz2d(i, :), &
            equil_out%global_param%psibound, -1, j)
        geom1d(66) = equil_out%eqgeometry%rectgrid%z2d(j)
        if (j < nz .and. j > 1) then
            geom1d(66) = equil_out%eqgeometry%rectgrid%z2d(j) - &
                (-gp2*equil_out%eqgeometry%rectgrid%psirz2d(i, j) - equil_out%global_param%psibound)/  &
                ( -gp2*equil_out%eqgeometry%rectgrid%psirz2d(i, j+1) + &
                   gp2*equil_out%eqgeometry%rectgrid%psirz2d(i, j-1) ) * &
                (equil_out%eqgeometry%rectgrid%z2d(j+1) - equil_out%eqgeometry%rectgrid%z2d(j-1))
        endif
        psi_sep  = equil_out%global_param%psibound
        psi_axis = equil_out%global_param%psiaxis
        i = minloc(abs(equil_out%eqgeometry%rectgrid%r2d - R_strike_out_aug), 1)  
        call find_in_vec_spid(nz, -gp2*equil_out%eqgeometry%rectgrid%psirz2d(i, :), &
            equil_out%global_param%psibound, -1, J) 
        geom1d(67) = equil_out%eqgeometry%rectgrid%z2d(j)
        if (j < nz .and. j > 1) then
            geom1d(67) = equil_out%eqgeometry%rectgrid%z2d(j) - & 
                (-gp2*equil_out%eqgeometry%rectgrid%psirz2d(i,   j) - equil_out%global_param%psibound)/   &
                (-gp2*equil_out%eqgeometry%rectgrid%psirz2d(i, j+1) + &
                  gp2*equil_out%eqgeometry%rectgrid%psirz2d(i, j-1)) * &
                (equil_out%eqgeometry%rectgrid%z2d(j+1) - equil_out%eqgeometry%rectgrid%z2d(j-1))
        endif
    else if (MACHINE(1:3) == 'dem') then
        if (jdemogaps == 0) then
            write(fname, '(a)') TRIM(parameters_equil%prename) // 'demo_gaps.data'
            open(32, file=fname)
            read(32, *) i_gaps
            do i=1, i_gaps
                read(32, *) demo_gaps(i, 1), demo_gaps(i, 2), demo_gaps(i, 3), demo_gaps(i, 4) !1-R,  2-Z,  3-angle,  4-0 if do both sides,  1 if only positive side
            enddo
            close(32)
            jdemogaps = 1
        endif
        if (equil_solver == 101) then
            call find_demo_gaps_efff(i_gaps, demo_gaps(1:i_gaps, 1:4), geom1d(94-i_gaps+1:94))
            call get_zccurb_efff(Rcurr, Zcurr, Zsquad, rgeoc, zgeoc, ahorc)
        else if (equil_solver == 3) then
            call find_demo_gaps(i_gaps, demo_gaps(1:i_gaps, 1:4), geom1d(94-i_gaps+1:94))
            call get_zccurb(Rcurr, Zcurr, Zsquad, rgeoc, zgeoc, ahorc)
        endif
        geom1d(95)  = zunten
        geom1d(96)  = zoben
        geom1d(97)  = Rcurr
        geom1d(98)  = Zcurr
        geom1d(99)  = Rin
        geom1d(100) = Raus
    endif
endif

return
end subroutine A_equil

!---------------------------------------------------------------------
subroutine A_equil_2(ncoils, ifbey, time_a, tau_step, vcoils, eq_solver)

use imas_ids, only: type_equilibrium
use parameters_a2equil, only: type_parameters
use outcmn_inc, only: MACHINE

implicit none

integer, intent(in) :: ifbey, eq_solver, ncoils
double precision, intent(in) :: tau_step, time_a
real*8, dimension(ncoils), intent(in) :: vcoils

integer :: nstep, key_equil
real*8, dimension(ncoils) :: ucoils

type(type_parameters) :: parameters_equil
type(type_equilibrium) :: equil_in, equil_out

key_equil = 0
nstep = max(0, ifbey - 1)

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

if (MACHINE(1:3) == 'aug') then
    ucoils(1) = vcoils(1) - vcoils(2)
    ucoils(2) = vcoils(2) - vcoils(3)
    ucoils(3) = vcoils(3)
    ucoils(4) = vcoils(4)
    ucoils(5) = vcoils(5)
    ucoils(6) = vcoils(6)
    ucoils(7) = vcoils(7)
    ucoils(8) = vcoils(8)
    ucoils(9) = vcoils(9)
    ucoils(10) = vcoils(10)
    ucoils(11:12) = 0.
elseif (MACHINE(1:3) == 'dem') then
    ucoils(1:ncoils) = vcoils(1:ncoils)
else
    ucoils(1:ncoils) = vcoils(1:ncoils)
endif
parameters_equil%nstep = nstep

if (eq_solver == 101) then
    call feqis_main(ncoils, ucoils, parameters_equil, 0, equil_in, equil_out)     
else if (eq_solver == 3) then
    call spider_run_2(ncoils, ucoils, parameters_equil)     
endif

return
end subroutine a_equil_2

!---------------------------------------------------------------------
subroutine find_in_vec_spid(n, y, y0, is, iv)

implicit none

integer, intent(in) :: n, is
double precision, intent(in) :: y0
double precision, intent(in), dimension(n) :: y
integer, intent(out) :: iv

integer :: i, j
integer, dimension(n) :: ic
double precision :: dum1, dum2
double precision, dimension(n) :: yy

yy = abs(y - y0)
j = 0
do i=2, n-1
    dum1 = yy(i) -   yy(i-1)
    dum2 = yy(i+1) - yy(i)
    if (dum1*dum2 <= 0) then
        j = j + 1
        ic(j) = i
    endif
enddo

if (is ==  1) iv = ic(j)
if (is == -1) iv = ic(1)

return
end subroutine find_in_vec_spid
