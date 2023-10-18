subroutine RABBIT

use mod_rabbit_lib, only: do_dump, rabbit_lib_init, rabbit_lib_set_dump_dir, &
    rabbit_lib_dump_beams, rabbit_lib_set_sp_plasma_ratio, rabbit_lib_step, &
    rabbit_lib_get_dv_darea, rabbit_lib_get_wfi
use rabbit_variables, only: fusion_power, neutron_power

use outcmn_inc, only: AWD, exp_file
use const_inc, only: GP2, AIM1, TIME, TAU, QNBI, ROC, &
   RTOR, BTOR, NA1
use status_inc, only: FP, AMAIN, ZMAIN, ZIM1, NE, TE, TI, &
   XRHO, VOLUM, IPOL, PEBM, PIBM, NIBM, CUBM, SNEBM, SCUBM, &
   PBLON, PBPER, MU, VTOR, ZEF, NI, NHYDR, NDEUT, NTRIT
 
use fs_coupling_variables, only: fs_pow_NB
use debugger, only: flightsim
use parameters_a2equil, only : equil_now

implicit none

integer, parameter :: Nrrect=64, Nzrect=64, nnb_max=30, nspc=3, nrhoout=21, unit_lim=11
double precision, parameter :: ALFA=1.d-5

integer, dimension(nnb_max) :: ierr
integer :: n_Rrect, n_Zrect, n_nbi, dum, n_lim, jumpcor, torqjxb_model
integer :: pdim, ldim
integer :: i, j, jlim, jnb, ios, nrho_surf, nthe_surf

double precision :: aimp, zimp, p_i, p_e, tq_i, fi, i_cd, src, nfi
double precision, allocatable, dimension(:) :: aplasma, zplasma, species_plasma_ratio
double precision, dimension(NA1) :: psi1d_n, rhotor1d
double precision, allocatable, dimension(:, :), save :: powe, powi, &
      press, bdep, bdens, jfi, jnbcd,  wfi_par, wfi_perp, wfi_par_lab, &
      torqe, torqi, torqjxb, torqth, torqthcxloss, torqdepo
     
double precision, dimension(nnb_max) :: a_beam, z_beam, pinj,  &
      einj, prot, powe_tot, powi_tot, pshine, porbloss, pcxloss,  &
      Inbcd
double precision, dimension(3, nnb_max) :: start_pos, unit_vec, width_poly

double precision, allocatable, dimension(:, :) :: PSI_rect
double precision, allocatable, dimension(:) :: Rrect, zrect
double precision :: psi_sep, psi_axis, rmag, zmag
double precision :: R_max, R_min, z_max, z_min, dr, dz, drho_eq
double precision :: part_mix(nspc, nnb_max), dt_in, output_timing 
double precision :: tim_prev=-1.d0, dumba1, dumba2

double precision, dimension(:), allocatable :: pf_eq, rho_eq
double precision, allocatable, dimension(:, :) :: r_surf, z_surf
double precision, dimension(NA1) :: rho_interp_plasma, rho_interp_eq, &
   ti_interp, te_interp, ne_interp, omg_interp, zef_interp,  &
   iota, area, vol, ffp, psi_n
double precision, dimension(nrhoout) :: rho_rab_out, &
    bdens_in, pi_rb, pe_rb, dvol, darea, &
    nfi_rb, jcd_rb, src_rb, tq_rb, pfi_par, pfi_perp, nrate
double precision, dimension(:), allocatable :: r_lim, z_lim

character(len=120) :: as_nml, pinj_file, pinj_file2, limiter_file, table_path

double precision, external :: VINT , IINT
namelist / rabbit_beam_geo / start_pos, unit_vec, width_poly
namelist / partmix / part_mix
namelist / nbi_par / n_nbi, a_beam, z_beam, einj, pinj_file
namelist / physics / jumpcor, table_path, limiter_file, torqjxb_model
save tim_prev, bdens_in, n_nbi, einj, part_mix, pinj_file, Rrect, zrect

write(6, *) 'Calling RABBIT...'

rhotor1d(1:NA1) = XRHO(1:NA1)
rhotor1d(NA1) = 1.d0
rhotor1d(1) = 0.d0

pdim = NA1    !Rabbit input plasma grid size
ldim = NA1    !Rabbit input 1D EQ grid size

call GET_NRHO_NTHETA(nrho_surf, nthe_surf)
allocate(pf_eq(nrho_surf), rho_eq(nrho_surf))
allocate(r_surf(nrho_surf, nthe_surf), z_surf(nrho_surf, nthe_surf))
call SURF_CTR(nrho_surf, nthe_surf, r_surf, z_surf)

psi_axis = FP(1)/GP2
psi_sep  = FP(NA1)/GP2
rmag = r_surf(1, 1)
zmag = z_surf(1, 1)

drho_eq = 1./(nrho_surf - 1.d0)
rho_eq = (/ (drho_eq*(i - 1.d0), i=1, nrho_surf) /)

psi1d_n = (FP(1:NA1) - FP(1))/(FP(NA1) - FP(1))
psi1d_n(1) = 1.d-8
psi1d_n(NA1) = 1.d0

if (.not. allocated(aplasma)) then
    allocate(aplasma(0))
    allocate(zplasma(0))
endif
if (maxval(NHYDR(1:NA1)) .gt. 0.) then
    aplasma = [aplasma, 1.]
    zplasma = [zplasma, 1.]
endif
if (maxval(NDEUT(1:NA1)) .gt. 0.) then
    aplasma = [aplasma, 2.]
    zplasma = [zplasma, 1.]
endif
if (maxval(NTRIT(1:NA1)) .gt. 0.) then
    aplasma = [aplasma, 3.]
    zplasma = [zplasma, 1.]
endif
if (.not. allocated(species_plasma_ratio)) allocate(species_plasma_ratio(size(Aplasma)))

do i=1, size(Aplasma)
   if (nint(Aplasma(i)) == 1) species_plasma_ratio(i) = sum(nhydr)
   if (nint(Aplasma(i)) == 2) species_plasma_ratio(i) = sum(ndeut)
   if (nint(Aplasma(i)) == 3) species_plasma_ratio(i) = sum(ntrit)   
enddo
species_plasma_ratio = species_plasma_ratio / sum(species_plasma_ratio)

if (flightsim == 0) then
    n_Rrect = Nrrect
    n_Zrect = Nzrect
else
    n_Rrect = SIZE(equil_now%eqgeometry%rectgrid%r2d)
    n_Zrect = SIZE(equil_now%eqgeometry%rectgrid%z2d)
endif
if (.not. allocated(psi_rect)) allocate(psi_rect(n_Rrect, n_Zrect))
if (.not. allocated(Rrect)) allocate(Rrect(n_Rrect), Zrect(n_Zrect))

if (tim_prev == -1.d0) then  ! --- RABBIT Initialization ---       
    as_nml = TRIM(AWD) // 'exp/nml/' // TRIM(exp_file)

    write(6, *) 'Parsing namelist ' // TRIM(as_nml)
    open(53, FILE=TRIM(as_nml), delim='apostrophe', iostat=ios)
    if (ios < 0) then
        write(*, *) 'RABBIT namelist ' // trim(as_nml) // ' not found, returning'
        return
    endif 
    read(53, nml=rabbit_beam_geo, iostat=ios)
    read(53, nml=partmix, iostat=ios)
    read(53, nml=nbi_par, iostat=ios)
    read(53, nml=physics, iostat=ios)
    close(53)
    write(6, *) '#NBI', n_nbi
    write(6, *) 'Limiter file', TRIM(limiter_file)
    open(unit_lim, file=TRIM(limiter_file), iostat=ios)
    read(unit_lim, '(2i)') dum, n_lim

    allocate(r_lim(n_lim), z_lim(n_lim))

    do jlim=1, n_lim
        read(unit_lim, *) r_lim(jlim), z_lim(jlim)
    enddo
    close(unit_lim)

    if (flightsim == 0) then
        R_min = MINVAL(R_lim) - 0.05
        R_max = MAXVAL(R_lim) + 0.05
        z_min = MINVAL(z_lim) - 0.05
        z_max = MAXVAL(z_lim) + 0.05
        dr = (R_max - R_min)/(n_Rrect - 1.d0)
        dz = (z_max - z_min)/(n_Zrect - 1.d0)
        Rrect = (/ (R_min + dr*(i - 1.d0), i=1, n_Rrect) /)
        Zrect = (/ (z_min + dz*(i - 1.d0), i=1, n_Zrect) /)
    else
        Rrect = equil_now%eqgeometry%rectgrid%r2d
        Zrect = equil_now%eqgeometry%rectgrid%z2d
    endif

!    aplasma = AMAIN(1)
!    zplasma = ZMAIN(1)
    aimp = AIM1
    zimp = ZIM1(1)
 
    if (zimp /= 4 .AND. zimp /= 5 .AND. zimp /= 6  .AND. zimp /= 7 .AND. zimp /= 10 .AND. zimp /= 18 .AND. zimp /= 28) then
        write(6, *) 'No cross-sections for Zimp other than 4, 5, 6, 7, 10, 18, 28'
        write(6, *) 'Forcing Zimp=6'
        zimp = 6
    endif

    call rabbit_lib_init(aplasma, zplasma, aimp, zimp,            & ! plasma species
        size(Aplasma),                                           &  !number of main ion species
        a_beam(1:n_nbi), z_beam(1:n_nbi), start_pos(:, 1:n_nbi), &
        unit_vec(:, 1:n_nbi), width_poly(:, 1:n_nbi),            & ! beam
        nspc, n_nbi,                                             & ! beam
        nrhoout,                                                 & ! output grid dimension
        Rrect, zrect,                                            & ! eq flux matrix grid
        N_Rrect, N_Zrect,                                          & ! eq grid dimensions
        ldim, pdim,                                              & ! plasma grid dimension
        TRIM(as_nml), LEN_TRIM(as_nml), ierr(1:n_nbi))

    if (do_dump) then
       !activate this to enable dumping of Rabbit inputs (for debbuging)
        call rabbit_lib_set_dump_dir(TRIM(awd), LEN_TRIM(awd))
        call rabbit_lib_dump_beams(TRIM(awd), LEN_TRIM(awd), einj, part_mix)
    endif

    tim_prev = max(0.d0, TIME-TAU)
    bdens_in(:) = 0.

! below arrays are saved, i.e. need to be allocated only once.
    allocate(powe(nrhoout, n_nbi))
    allocate(powi(nrhoout, n_nbi))
    allocate(press(nrhoout, n_nbi))
    allocate(bdep(nrhoout, n_nbi))
    allocate(bdens(nrhoout, n_nbi))
    allocate(jfi(nrhoout, n_nbi))
    allocate(jnbcd(nrhoout, n_nbi))
    allocate(torqe(nrhoout, n_nbi))
    allocate(torqi(nrhoout, n_nbi))
    allocate(torqjxb(nrhoout, n_nbi))
    allocate(torqth(nrhoout, n_nbi))
    allocate(torqthcxloss(nrhoout, n_nbi))
    allocate(torqdepo(nrhoout, n_nbi))
    allocate(wfi_par(nrhoout, n_nbi))
    allocate(wfi_perp(nrhoout, n_nbi))
    allocate(wfi_par_lab(nrhoout, n_nbi))
endif

if (flightsim == 1) then
    pinj(1:8) = fs_pow_NB(1:8)*1d6
else
    pinj_file2 = TRIM(awd) // TRIM(pinj_file)
    call uf2dr(pinj_file2, TIME, pinj(1:n_nbi))
endif

QNBI = sum(pinj(1:n_nbi))*1d-6

output_timing = 0.5d0
dt_in = max(TIME - tim_prev, 1.d-6)

rho_interp_plasma = rhotor1d
ne_interp = 1e19*NE(1:NA1)
te_interp = 1e3*TE(1:NA1)
ti_interp = 1e3*TI(1:NA1)
omg_interp = VTOR(1:NA1)/RTOR
zef_interp = ZEF(1:NA1)

rho_interp_eq = rhotor1d
iota = MU(1:NA1)
vol = VOLUM(1:NA1)
ffp = IPOL(1:NA1)*BTOR*RTOR
psi_n = psi1d_n
vol(1) = 0.d0
psi_n(1) = 0.d0
area = vol/(GP2*RTOR)

call qinterp(rhotor1d(1: NA1), FP(1: NA1), NA1, &
      rho_eq(1: nrho_surf), pf_eq(1: nrho_surf), nrho_surf)

!------------
! RABBIT call
!------------

write(6, *) 'Call rabbit_lib_step'

if (flightsim == 0) then
    call ctr2rz_fun(nrho_surf, nthe_surf, pf_eq(1: nrho_surf)/GP2, &
         r_surf(1: nrho_surf, 1: nthe_surf),  z_surf(1: nrho_surf, 1: nthe_surf), &
         n_Rrect, n_Zrect, Rrect, zrect, PSI_rect)
else
    psi_rect = equil_now%eqgeometry%rectgrid%psirz2d(1:n_Rrect, 1:n_Zrect)
    dumba1 = equil_now%eqgeometry%rectgrid%psi_axis
    dumba2 = equil_now%eqgeometry%rectgrid%psi_boundary
    psi_rect = (psi_rect - dumba1)/(dumba2 - dumba1)
    psi_rect = (psi_sep - psi_axis)*psi_rect + psi_axis
endif

call rabbit_lib_set_sp_plasma_ratio(species_plasma_ratio, size(species_plasma_ratio))

call rabbit_lib_step(                                     & ! input
    rho_interp_plasma, ne_interp, te_interp, ti_interp,    &
    zef_interp, omg_interp, pdim,                          & ! Kin profiles & their dim
    PSI_rect, psi_n, vol, area, rho_interp_eq, iota, ffp,  & ! eq
    psi_sep, psi_axis, rmag, zmag,                         & ! eq scalars
    N_Rrect, N_Zrect, ldim,                                  & ! eq dimensions
    pinj(1: n_nbi), einj(1: n_nbi), part_mix(: , 1: n_nbi),    &
    nspc, n_nbi, bdens_in, dt_in, output_timing,           & ! Output
    powe, powi, press, bdep, bdens, jfi, jnbcd,            &
    torqe, torqi, torqjxb, torqth, torqthcxloss, torqdepo, &
    nrate, rho_rab_out, nrhoout,     &
    powe_tot(1: n_nbi), powi_tot(1: n_nbi), pshine(1: n_nbi), & 
    prot(1: n_nbi), porbloss(1: n_nbi), pcxloss(1: n_nbi),    &
    Inbcd(1: n_nbi), ierr(1: n_nbi))

call rabbit_lib_get_dV_dArea(dvol, darea, nrhoout)
call rabbit_lib_get_Wfi(wfi_par, wfi_perp, wfi_par_lab, n_nbi, nrhoout)
     
! Sum over all NBI sources

pe_rb   = sum(powe (: , 1: n_nbi), 2)/1.d6
pi_rb   = sum(powi (: , 1: n_nbi), 2)/1.d6
tq_rb   = sum(torqi(: , 1: n_nbi), 2)
nfi_rb  = sum(bdens(: , 1: n_nbi), 2)/1.d19
jcd_rb  = sum(jnbcd(: , 1: n_nbi), 2)/1.d6
src_rb  = sum(bdep (: , 1: n_nbi), 2)/1.d19

!      pfi_par = sum(press(: , 1: n_nbi), 2)/1602.d0
pfi_par = 2*sum(wfi_par_lab(: , 1: n_nbi), 2)/1602.d0
pfi_perp = sum(wfi_perp(: , 1: n_nbi), 2)/1602.d0
bdens_in = sum(bdens(: , 1: n_nbi), 2)   !store for next call (bdens_in is saved).

! Integrals

p_i  = sum(pi_rb*dvol)
p_e  = sum(pe_rb*dvol)
tq_i = sum(tq_rb*dvol)
src  = sum(src_rb*dvol)
nfi  = sum(nfi_rb*dvol)
i_cd = sum(jcd_rb*darea)

tim_prev = TIME

call qinterp(rho_rab_out, pi_rb   , nrhoout, XRHO(1: NA1), PIBM( 1: NA1), NA1)
call qinterp(rho_rab_out, pe_rb   , nrhoout, XRHO(1: NA1), PEBM( 1: NA1), NA1)
call qinterp(rho_rab_out, nfi_rb  , nrhoout, XRHO(1: NA1), NIBM( 1: NA1), NA1)
call qinterp(rho_rab_out, tq_rb   , nrhoout, XRHO(1: NA1), SCUBM(1: NA1), NA1)
call qinterp(rho_rab_out, jcd_rb  , nrhoout, XRHO(1: NA1), CUBM( 1: NA1), NA1)
call qinterp(rho_rab_out, src_rb  , nrhoout, XRHO(1: NA1), SNEBM(1: NA1), NA1)
call qinterp(rho_rab_out, pfi_par , nrhoout, XRHO(1: NA1), PBLON(1: NA1), NA1)
call qinterp(rho_rab_out, pfi_perp, nrhoout, XRHO(1: NA1), PBPER(1: NA1), NA1)

if (ALFA > 0.) then 
    call smearr(ALFA, PIBM , PIBM )
    call smearr(ALFA, PEBM , PEBM )
    call smearr(ALFA, NIBM , NIBM )
    call smearr(ALFA, CUBM , CUBM )
    call smearr(ALFA, SCUBM, SCUBM)
    call smearr(ALFA, SNEBM, SNEBM)
endif

write(*, * ) 'Pi  integral', p_i/VINT(PIBM, ROC)
write(*, * ) 'Pe  integral', p_e/VINT(PEBM, ROC)
write(*, * ) 'Tq  integral', tq_i/VINT(SCUBM, ROC)
write(*, * ) 'Icd integral', i_cd/IINT(CUBM, ROC)
if (p_i > 0.d0) then
    PIBM(1: NA1)  = PIBM( 1: NA1)*p_i/VINT(PIBM, ROC)
    PEBM(1: NA1)  = PEBM( 1: NA1)*p_e/VINT(PEBM, ROC)
    NIBM(1: NA1)  = NIBM( 1: NA1)*nfi/VINT(NIBM, ROC)
    if (src > 0.d0) SNEBM(1: NA1) = SNEBM(1: NA1)*src/VINT(SNEBM, ROC)
    SCUBM(1: NA1) = SCUBM(1: NA1)*tq_i/VINT(SCUBM, ROC)
    CUBM(1: NA1)  = CUBM( 1: NA1)*i_cd/IINT(CUBM, ROC)
    PBLON(1: NA1) = PBLON(1: NA1)*sum(pfi_par *dvol)/VINT(PBLON, ROC)
    PBPER(1: NA1) = PBPER(1: NA1)*sum(pfi_perp*dvol)/VINT(PBLON, ROC)
endif

write(6, *) 'Done RABBIT'

return
end subroutine RABBIT
