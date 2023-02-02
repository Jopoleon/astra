subroutine A2STRAHL(tau_start, zneocl, dzneocl, dimpsol, ydimp, yvimp, &
    y_zcharge, ymimp0, ynimp, yzeff, ydneo, yvneo, rrates_in, &
    prate1, prate2, shot_in)

!----------------------------------------------------------------------
!    - G. Tardini, Sep 2022
!    - E. Fable, Feb 2012 -   CCCs
!
!  Coupling between ASTRA and STRAHL (R. Dux)
!
!  Inputs to STRAHL:
!         rhop = sqrt(psi/psib), Ne, Te, Ti,
!         r[cm] = sqrt(V/(2*pi^2*Rgeo)), function of rhop
!         R(z = 0)LFS [cm], function of rhop
!         R(z = 0)HFS [cm], function of rhop
!
!         Vloop, q, fc (passing frac), int(dl/Bp)
!
!         <B>, <B^2>, <1/B^2>, F, <R^2*Bp^2/B^2>, <1/R^2>
!
!         <cos(mt)*B^2>, <sin(mt)*B^2>
!         <cos(mt)*B*log(B)>, <sin(mt)*B*log(B)>     where m is given as input
!
!         Impurity species n
!
!         ydimp = D_anom [m^2/s], yvimp = V_anom [m/s]
!
!         tau_start : beginning of STRAHL
!                zneocl : beginning of NEOART
!         dzneocl: dt of NEOART
!         dimpsol: SOL diffusivity
!
! Output: yzeff, yprad [MW/m^2, +]
!         ynimp [10^19 m^-3]
!         y_zcharge [e] (radial function)
!
!         n_e = (main ions) + n_imp
!         prad_sep are work_strahl(:, j), j = 4...24
!----------------------------------------------------------------------|

use parameter_inc, only: NRD
use const_inc, only: TIME, TSTART, TAUPRP, NA1, PSIAX, GP, GP2, RTOR, NA, HRO, IPART
use status_inc, only: FP, UPL, VOLUM, SHIF, NE, TE, TI, AMAIN, ZMAIN
use outcmn_inc, only: machine, awd, exp_file
use strahl, only: profiles_file_write_strahl, grid_write_strahl, &
    prad_tot, ne_source, nneut_imp, prad_strahl, nimp_strahl, nesrc_strahl

implicit none

integer, parameter :: ngmax=140, n_o_max=512, nch_r1=31, nch_r2=32, nch_w1=41, nch_w4=44
character(len=120), parameter :: strahl_output='results.txt', strahl_param_in='param_files/sparams.dat'

double precision, intent(in) :: tau_start, zneocl, dzneocl, dimpsol, &
    rrates_in, prate1, prate2, shot_in
double precision, dimension(NRD), intent(in) :: ydneo, yvneo
double precision, intent(out) :: ymimp0
double precision, dimension(NRD), intent(out) :: ydimp, yvimp, y_zcharge, ynimp, yzeff

integer :: i, j, k, shotn, indexx, nimp_touse, nfour_c, n_grids, &
    Nr_o, ineocl, ineocla, ineocli, i_stepst, ios
integer, dimension(10) :: irecycl

real*8 :: tneocl, tneocl0, rneocl
double precision :: dum1, tau_strahl, &
    ne_decayl, te_decayl, ti_decayl, z_K, zdr_0, zdr_1, rbrlcfs, rlimrlcfs, &
    tolimiter, solflow, solrout1, solrout2, solrout3, solrout4, todivert, addsheathvoltage
double precision, dimension(NRD) :: yprad, rhopol, rhovol, y_mions, ymimp, r_rho, dneo_o, vneo_o
double precision, dimension(10) :: aweight, eneutr, rsources, rrates, trates, ridecay, &
    wrecycl, divpuff, swincm, swoutcm, promptredep, taudiv, taupump
double precision, dimension(n_o_max) :: rpol_o, zeff_o, prad_o, nimp_o, zimp_o, &
    mimp_o, mion_o, nimpneutr_o
double precision, dimension(n_o_max, 10) :: nimpsep_o
double precision, dimension(n_o_max, 11) :: pradsep_o, ne_source_o
double precision, dimension(ngmax) :: rhopolg, neg, teg, tig

character(len=160) :: strahl_dir, cmd_cmd, as_nml
character(len=20) :: rho_coord, elements_touse(10)
character(len=6) :: diffname1_s

data i_stepst /0/
data ineocli /0/
data tneocl0 /0./
save i_stepst, tneocl0, ineocli

NAMELIST / strahl /  tau_strahl, rho_coord, ne_decayl, te_decayl, ti_decayl, &
    nfour_c, nimp_touse, elements_touse, aweight, eneutr, ridecay, irecycl, &
    wrecycl, diffname1_s, z_K, n_grids, zdr_0, zdr_1, rsources, rrates, trates, &
    ineocla, rneocl, divpuff, swincm, swoutcm, promptredep, taudiv, taupump, &
    rbrlcfs, rlimrlcfs, todivert, tolimiter, solflow, addsheathvoltage, &
    solrout1, solrout2, solrout3, solrout4

!--------------------------------------------------------------------------

as_nml = TRIM(awd) // 'exp/nml/' // TRIM(exp_file)
write(*, *) 'Reading namelist ', TRIM(as_nml)

open(nch_r1, FILE=TRIM(as_nml), delim='apostrophe')
read(nch_r1, nml=strahl, iostat=ios)
close(nch_r1)

ineocl = 0
tneocl = TIME - TSTART - tneocl0
if (TIME > zneocl+TAUPRP) then
    if (tneocl > dzneocl) then
        ineocl = 1
        ineocli = 1
        tneocl0 = TIME - TSTART
    endif
endif

if (TIME-TSTART < tau_start) then
    y_zcharge = 1.0
    ymimp0 = 1.0
    ynimp = 0.0
    yzeff = 1.0
    yprad = 0.0
    return
endif

diffname1_s = ''

shotn = nint(shot_in)
strahl_dir = TRIM(awd) // 'strahl/'

call chdir(TRIM(strahl_dir))      ! cdir
call system('mkdir -p result')
call system('mkdir -p nete')
call system('mkdir -p param_files')

indexx = 1

!rhopol, rhovol
do j=1, NA1
    rhopol(j) = (FP(J) - PSIAX)**0.5/(FP(NA1) - PSIAX)**0.5
    rhovol(j) = (VOLUM(j)/(GP2*GP*(RTOR + SHIF(1))))**0.5
enddo

do j=1, NA
    r_rho(j) = (rhovol(j+1) - rhovol(j))/HRO
enddo
r_rho(NA1) = r_rho(NA)

!conversion from rho to rhovol of STRAHL for diffusion and convection
ydimp = ydimp*r_rho
yvimp = yvimp*r_rho

rhopolg(1) = 0.
do j=2, ngmax
    rhopolg(j) = rhopolg(j-1) + 1./(ngmax - 1.)
enddo

if (NA1 > ngmax) then
! Interpolate to astra grid
    call qinterp_metric(rhopol(1:NA1), ydimp(1:NA1), NA1, rhopolg(1:ngmax), ydimp(1:ngmax), ngmax)
    call qinterp_metric(rhopol(1:NA1), yvimp(1:NA1), NA1, rhopolg(1:ngmax), yvimp(1:ngmax), ngmax)
    call qinterp_metric(rhopol(1:NA1), NE(1:NA1), NA1, rhopolg(1:ngmax), neg(1:ngmax), ngmax)
    call qinterp_metric(rhopol(1:NA1), TE(1:NA1), NA1, rhopolg(1:ngmax), teg(1:ngmax), ngmax)
    call qinterp_metric(rhopol(1:NA1), TI(1:NA1), NA1, rhopolg(1:ngmax), tig(1:ngmax), ngmax)
endif

rrates(1)=rrates_in

tau_strahl = min(TAUPRP, tau_strahl)
i = nint(TAUPRP/tau_strahl)
tau_strahl = TAUPRP/i

! Produce STRAHL input files

!Start with main parameter file
open(nch_w1, file=TRIM(strahl_dir)//TRIM(strahl_param_in))
write(nch_w1, '(A)') &
    '               M A I N  I O N ', &
    '   ', &
    'cv    background ion:  atomic weight    charge  '
!Deuterium: A=2.0    Z=1.0
write(nch_w1, '(A, F12.4, A, F12.4)') '          ', amain(1), '   ', zmain(1)
! remember to put variable instead of 2.some for weight and charge
write(nch_w1, '(A)') &
    '   ', &
    '   			G R I D - F I L E', &
    'cv    shot      index'
write(nch_w1, 109) ' 	', 11111, '   ', indexx
write(nch_w1, '(A)') '   ', &
    '       G R I D   P O I N T S   A N D  I T E R A T I O N', &
    'cv     rho = r**K (->K)       number of grid points      dr_0       dr_1'
write(nch_w1, 104) ' ', z_K, ' ', n_grids, ' ', zdr_0, ' ', zdr_1
write(nch_w1, '(A)') &
    '   ', &
    'cv      max. iterations at fixed time  stop iteration if change below(%)', &
    '  1000          -1.		', &
    '   ', &
    '      S T A R T   C O N D I T I O N S', &
    '   ', &
    'cv    start new=0/from old impurity   distribution=1     shot   at    time  index'
! Notice that start new is used if time below 1e-3. should be changed...
if (i_stepst == 0) then
    write(nch_w1, 109) '   0   ', 11111, '   0.0    ', indexx
endif
if (i_stepst == 1) then
    write(nch_w1, 109) '   2  ', 11111, '  0.0  ', indexx
endif
if (i_stepst > 1 .and. ineocli == 0) then
    write(nch_w1, 109) '   2  ', 11111, '  0.0  ', indexx
endif
if (i_stepst > 1 .and. ineocli == 1) then
    if (ineocl == 0 .and. ineocla == 1) then
        write(nch_w1, 109) '   3  ', 11111, '  0.0  ', indexx
    else
        write(nch_w1, 109) '   2  ', 11111, '  0.0  ', indexx
    endif
endif
write(nch_w1, '(A)') &
    '   ', &
    '           O U T P U T   ', &
    '   ', &
    'cv    save all cycles = 1,   save final and start distribution = 0', &
    ' 		           2', &
    '   ', &
    '                 T I M E S T E P S   ', &
    '   ', &
    'cv    number of changes   (start-time+... +stop-time)', &
    '                  2', &
    '   ', &
    '   ', &
    'cv    time    dt at start     increase of dt after cycle      steps per cycle'
write(nch_w1, 105) ' ', TIME, ' ', tau_strahl, ' 1.0  1'
write(nch_w1, 105) ' ', TIME+TAUPRP, ' ', TAUPRP, '  1.0   1'
write(nch_w1, '(A)') &
    '   ', &
    '   ', &
    '            S T A R T     I M P U R I T Y   E L E M E N T S', &
    '       (for each impurity  one input line needed in this block)	', &
    '   ', &
    'cv     number of impurities'
write(nch_w1, 109) '  ', nimp_touse
write(nch_w1, '(/A)') 'cv     element   atomic weight   energy of neutrals(eV)'
do i=1, nimp_touse
    write(nch_w1, '(3A, 2F12.4)') ' ', elements_touse(i), '  ', aweight(i), eneutr(i)
enddo
write(nch_w1, '(A)') &
    '                            ', &
    '    		      S O  U R C E  ', &
    '   ', &
    'cv  r_source-r_lcfs(cm)   constant rate(1/s)    time dependent rate from file(1/0)'
do i=1, nimp_touse
    write(nch_w1, 108) '  ', rsources(i), '   ', rrates(i), '  ', trates(i)
enddo
write(nch_w1, '(/A)') &
    'cv    divertor puff    source_width_in(cm)     source_width_out(cm)   prompt redep'
do i=1, nimp_touse
    write(nch_w1, 112) '	  ', divpuff(i), '   ', swincm(i), '  ', swoutcm(i), '   ', promptredep(i)
enddo
write(nch_w1, '(A)') &
    '', &
    '                    E D G E ,   R E C Y C L I N G', &
    '', &
    'cv    decay length of  impurity outside last grid point(cm)'
do i=1, nimp_touse
    write(nch_w1, *) '                           ', ridecay(i), '           '
enddo
write(nch_w1, '(A)') &
    '                          ', &
    '', &
    'cv    Rec.:ON=1/OFF=0    wall-rec.  Tau-div->SOL(ms)    Tau-pump(ms) '

! reciclying has to be for each imp species
do i=1, nimp_touse
    write(nch_w1, '(A, I, A, F12.4, A, E14.5, A, F12.4)') &
        '   ', irecycl(i), '               ', wrecycl(i), '        ', taudiv(i), '  ', taupump(i)
enddo
write(nch_w1, '(A)') &
    '', &
    ' ', &
    '          E N D    I M P U R I T Y    E L E M E N T S  ', &
    '          ', &
    '                                               Connection lenghts  [m]   Mach #', &
    'cv    r_bound-r_lcfs (cm)   r_lim-r_lcfs(cm)   to divertor    to limiter    SOL Flow'
write(nch_w1, 145) '         ', rbrlcfs, '   ', rlimrlcfs, '     ', todivert, &
    '   ', tolimiter, '   ', solflow
write(nch_w1, '(/A)') 'cv   additional sheath voltage [V]'
write(nch_w1, *) addsheathvoltage
write(nch_w1, '(A)') &
    '', &
    '', &
    '  D E N S I T Y, T E M P E R A T U R E  AND N E U T R A L  H Y D R O G E N  F O R  CX ', &
    ' ', &
    'cv    take from file with:     shot        index'
write(nch_w1, *) '                               11111           ', indexx
write(nch_w1, '(A)') &
    '', &
    '', &
    '                    N E O C L A S S I C A L     T R A N S P O R T ', &
    '                                    NEOART with', &
    '     0 = off,  >0 = %  of Drift,    2= one stage      no BP      max        min', &
    'cv  <0 =figure out, but  dont use   3= all stages    contrib   rho_pol   rho_pol'
if (i_stepst <= 1) then
    write(nch_w1, *) '  0    3   0   1.0  ', rneocl
endif
if (i_stepst > 1 .and. ineocli == 0) then
    write(nch_w1, *) '  0    3   0   1.0  ', rneocl
endif
if (i_stepst > 1 .and. ineocli == 1) then
    write(nch_w1, *) '  100    3   0   1.0  ', rneocl
endif
write(nch_w1, '(A)') &
    '                ', &
    '                    A N O M A L O U S     T R A N S P O R T ', &
    ' ', &
    'cv   # of changes  for transport', &
    '                   1  ', &
    ' ', &
    'cv    time-vector', &
    '              0.00', &
    ' ', &
    'cv      Diffusion  [m^2/s]', &
    "        'interp'", &
    '', &
    ' ', &
    'cv   # of interpolation points'
write(nch_w1, *) '          ', min(NA1, ngmax)-1+4
write(nch_w1, '(/A)') ''
write(nch_w1, '(A)') 'cv   rho poloidal grid for interpolation'
do i=1, min(NA1, ngmax)-1
    if (NA1 > ngmax) then
        write(nch_w1, 101) '     ', rhopolg(i)
    else
        write(nch_w1, 101) '     ', rhopol(i)
    endif
enddo
write(nch_w1, 101) '     ', solrout1
write(nch_w1, 101) '     ', solrout2
write(nch_w1, 101) '     ', solrout3
write(nch_w1, 101) '     ', solrout4
write(nch_w1, '(A)') ''
write(nch_w1, '(A)') 'cv    D[m**2/s]'
do i=1, min(NA1, ngmax)-1
    write(nch_w1, 101) '     ', ydimp(i)
enddo
write(nch_w1, 101) '     ', ydimp(min(NA1, ngmax) - 1) + (1. - rhopol(min(NA1, ngmax) - 1)) * &
    (dimpsol - ydimp(min(NA1, ngmax) - 1))/(1.05 - rhopol(min(NA1, ngmax) - 1))
write(nch_w1, 101) '     ', dimpsol
write(nch_w1, 101) '     ', dimpsol
write(nch_w1, 101) '     ', dimpsol
write(nch_w1, '(A)') &
    '               ', &
    '', &
    '', &
    'cv    Drift function        Drift Parameter/Velocity', &
    "      'interp'                 'velocity'   ", &
    ' ', &
    ' ', &
    'cv   # of interpolation points'
write(nch_w1, *) '         ', min(NA1, ngmax)+4
write(nch_w1, '(/A)') 'cv   rho poloidal grid for      interpolation'
do i=1, min(NA1, ngmax)
    if(NA1 > ngmax) then
        write(nch_w1, 101) '     ', rhopolg(i)
    else
        write(nch_w1, 101) '     ', rhopol(i)
    endif
enddo
write(nch_w1, 101) '     ', solrout1
write(nch_w1, 101) '     ', solrout2
write(nch_w1, 101) '     ', solrout3
write(nch_w1, 101) '     ', solrout4
write(nch_w1, *) ''
write(nch_w1, '(A)') 'cv    V[m/s]'
do i=1, min(NA1, ngmax)
    write(nch_w1, 101) '     ', yvimp(i)
enddo
write(nch_w1, 101) '     ', 0.0
write(nch_w1, 101) '     ', 0.0
write(nch_w1, 101) '     ', 0.0
write(nch_w1, 101) '     ', 0.0
write(nch_w1, '(A)') &
    '', &
    '', &
    'cv    # of sawteeth          inversion radius ', &
    '     0        25.        ', &
    ' ', &
    ' ', &
    'cv    times of sawteeth ', &
    '            0.0'

close(nch_w1)

101 format(A, F12.4)
104 format(A, F12.4 , A, I, A, F12.4, A, F12.4)
105 format(A, E25.11, A, E25.11, 10A)
108 format(A, F12.4 , A, E16.4, A, E16.4)
109 format(A, I, A , I, A, I, A)
112 format(A, F12.4 , A, F12.4, A, F12.4, A, F12.4)
145 format(A, F12.4 , A, F12.4, A, F12.4, A, F12.4, A, F12.4)
! end call params_file_write_strahl

call profiles_file_write_strahl(strahl_dir, &
    rho_coord, rhopol(1:NA1), &
    NA1, NE(1:NA1), TE(1:NA1), TI(1:NA1), &
    ne_decayl, te_decayl, ti_decayl, TIME-TIME, &
    rhopolg(1:ngmax), teg(1:ngmax) , neg(1:ngmax), tig(1:ngmax), ngmax)

call grid_write_strahl(strahl_dir, nfour_c, RTOR+SHIF(1), &
    rhovol(NA1), UPL(NA1), TIME-TIME, machine)

! Main STRAHL call

cmd_cmd = '/afs/ipp/home/r/rld/strahl/amd64_sles11/strahl a q v'
write(*, '(A)') 'Executing', cmd_cmd 
call system(cmd_cmd)      ! run strahl

cmd_cmd = 'rm -f results.txt'
call system(cmd_cmd)    ! rm old results, if existing

cmd_cmd = '/afs/ipp/home/r/rld/strahl/amd64_sles11/result_to_astra '//trim(elements_touse(1))//' > ' // TRIM(strahl_dir) // 'results.txt'
write(*, '(A)') 'Executing', cmd_cmd 
call system(cmd_cmd)   ! produce new result file

call chdir(TRIM(awd))      ! cdir

! Extract results: rho poloidal, zeff, prad_tot, n_main, n_Z_tot, z_av, m_av, Dneo, Vneo, prad_Z1, prad_Z2, ... prad_Zn, prad_main, nZ1, nZ2, ... nZn, ne_source_Z1, ne_source_Z2, ... ne_source_Zn

open(nch_r2, file=TRIM(strahl_dir)//TRIM(strahl_output))
    read(nch_r2, *) Nr_o
    read(nch_r2, *) cmd_cmd
    read(nch_r2, *) (rpol_o(i), i=1, Nr_o)
    read(nch_r2, *) cmd_cmd
    read(nch_r2, *) (zeff_o(i), i=1, Nr_o)
    read(nch_r2, *) cmd_cmd
    read(nch_r2, *) (prad_o(i), i=1, Nr_o)
    read(nch_r2, *) cmd_cmd
    read(nch_r2, *) (mion_o(i), i=1, Nr_o)
    read(nch_r2, *) cmd_cmd
    read(nch_r2, *) (nimp_o(i), i=1, Nr_o)  !these are only ionized imps, no neutrals
    read(nch_r2, *) cmd_cmd
    read(nch_r2, *) (zimp_o(i), i=1, Nr_o)  !these are only ionized imps, no neutrals
    read(nch_r2, *) cmd_cmd
    read(nch_r2, *) (mimp_o(i), i=1, Nr_o)
    read(nch_r2, *) cmd_cmd
    read(nch_r2, *) (dneo_o(i), i=1, Nr_o)
    read(nch_r2, *) cmd_cmd
    read(nch_r2, *) (vneo_o(i), i=1, Nr_o)
    do j=1, nimp_touse+1
        read(nch_r2, *) cmd_cmd
        read(nch_r2, *) (pradsep_o(i, j), i=1, Nr_o)
    enddo
    do j=1, nimp_touse
        read(nch_r2, *) cmd_cmd
        read(nch_r2, *) (nimpsep_o(i, j), i=1, Nr_o)       !these are only ionized imps, no neutrals
    enddo
    do j=1, nimp_touse
        read(nch_r2, *) cmd_cmd
        read(nch_r2, *) (ne_source_o(i, j), i=1, Nr_o)
    enddo
    read(nch_r2, *) cmd_cmd
    read(nch_r2, *) (nimpneutr_o(i), i=1, Nr_o)
close(nch_r2)

! Interpolate to astra grid
call qinterp_metric(rpol_o(1:Nr_o), max(1., zeff_o(1:Nr_o)), Nr_o, &
    rhopol(1:NA1), yzeff(1:NA1), NA1)

call qinterp_metric(rpol_o(1:Nr_o), mimp_o(1:Nr_o), Nr_o, &
    rhopol(1:NA1), ymimp(1:NA1), NA1)

call qinterp_metric(rpol_o(1:Nr_o), 1.E-19*abs(mion_o(1:Nr_o)), Nr_o, &
    rhopol(1:NA1), y_mions(1:NA1), NA1)

call qinterp_metric(rpol_o(1:Nr_o), nimpneutr_o(1:Nr_o)/1.e19, Nr_o, &
    rhopol(1:NA1), nneut_imp(1:NA1), NA1)

call qinterp_metric(rpol_o(1:Nr_o), 1.E-19*abs(nimp_o(1:Nr_o)), Nr_o, &
    rhopol(1:NA1), ynimp(1:NA1), NA1)

call qinterp_metric(rpol_o(1:Nr_o), prad_o(1:Nr_o)/1.E6, Nr_o, &
    rhopol(1:NA1), yprad(1:NA1), NA1)

open(nch_w4, file='fort.224')
    do j=1, NA1
        prad_tot(j) = yprad(j)
        write(nch_w4, '(6E25.11)') rhopol(j), TE(j), NE(j), yprad(j), ynimp(j), nneut_imp(3)
    enddo
close(nch_w4)

call qinterp_metric(rpol_o(1:Nr_o), zimp_o(1:Nr_o), Nr_o, &
    rhopol(1:NA1), y_zcharge(1:NA1), NA1)

if (ineocl == 1) then
    call qinterp_metric(rpol_o(1:Nr_o), dneo_o(1:Nr_o), Nr_o, rhopol(1:NA1), ydneo(1:NA1), NA1)
    call qinterp_metric(rpol_o(1:Nr_o), vneo_o(1:Nr_o), Nr_o, rhopol(1:NA1), yvneo(1:NA1), NA1)
endif

!conversion from strahl rvol to astra rho
dneo_o = dneo_o/r_rho
vneo_o = vneo_o/r_rho

do i=1, nimp_touse+1
    call qinterp_metric(rpol_o(1:Nr_o), pradsep_o(1:Nr_o, i), Nr_o, &
        rhopol(1:NA1), prad_strahl(1:NA1, i), NA1)
enddo
do i=1, nimp_touse
    call qinterp_metric(rpol_o(1:Nr_o), nimpsep_o(1:Nr_o, i), Nr_o, &
        rhopol(1:NA1), nimp_strahl(1:NA1, i), NA1)
    call qinterp_metric(rpol_o(1:Nr_o), ne_source_o(1:Nr_o, i), Nr_o, &
        rhopol(1:NA1), nesrc_strahl(1:NA1, i), NA1)
enddo
do j=1, NA1
    ne_source(j) = sum(nesrc_strahl(j, 1: nimp_touse)) /1.E+19
enddo

!neutrals from Zimp

ymimp0 = sum(ymimp(1:NA1))/NA1

i_stepst = i_stepst+1

if (i_stepst > 5) i_stepst=5

201 format(1F12.4)
203 format(1I8)

return
end subroutine a2strahl
