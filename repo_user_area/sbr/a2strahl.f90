module strahl_mod

use status, only: NRD

implicit none

double precision, dimension(NRD) :: zeff_strahl, prad_tot_strahl, &
    nmain_strahl, prad_main_strahl
double precision, dimension(NRD, 11) :: prad_strahl, nimp_strahl, &
    zavg_strahl, nesrc_strahl, Dneo_strahl, Vneo_strahl, &
    Dz_in_strahl, Vz_in_strahl
double precision, dimension(11) :: rrates_in_strahl 

contains

!----------------------------------------------------------------------
subroutine A2STRAHL(tau_start, zneocl, dzneocl, dimpsol, shot_in)

!============================================================================================!
!    - D. Fajardo, Feb 2023
!    - G. Tardini, Sep 2022
!    - E. Fable, Feb 2012 -   CCCs
!
!  Coupling between ASTRA and STRAHL (R. Dux)
!
! INPUTS
! ======
! * To subroutine called from equ file:
! -------------------------------------
! - tau_start -> ASTRA time at which STRAHL starts
! - zneocl ----> time at which NEOART starts within STRAHL. Set to large number to skip NEOART
! - dzneocl ---> time step for NEOART
! - dimpsol ---> scrape-off layer diffusivity for the impurities, in [m^2/s]
! - shot_in ---> shot number (can be set to zero)
!
! * Given in equ file to array names from src/for/strahl_mod.f90,
!   with isp = 1...nsp being the index of each impurity species:
! --------------------------------------------------------------
! - Dz_in_strahl(:,isp) ---> anomalous diffusion coefficient profile [m^2/s] of species <sp>
! - Vz_in_strahl(:,isp) ---> anomalous convection velocity profile [m/s] of <sp>
! - rrates_in_strahl(isp) -> edge source of impurities [particles/s] of <sp>
!
! OUTPUTS
! =======
! * In src/for/strahl_mod.f90:
! --------------------------------------------------------
! - zeff_strahl(:)------> Zeff [-], effective charge profile
! - prad_tot_strahl(:)--> Prad_tot [MW/m^3], total radiated power density profile
! - nmain_strahl(:)-----> nmain [10^19/m^3], main ion density profile
! - prad_main_strahl(:)-> Prad_main [MW/m^3], radiated power density from bulk plasma 
! - prad_strahl(:,isp)--> Prad_<sp> [MW/m^3], radiated power density profile of <sp>
! - nimp_strahl(:,isp)--> nimp_<sp> [10^19/m^3], impurity density profile of <sp>
! - zavg_strahl(:,isp)--> Zavg_<sp> [-], average charge profile of <sp>
! - nesrc_strahl(:,isp)-> nesrc_<sp> [10^19/m^3/s], electron density source due to <sp>
! - Dneo_strahl(:,isp)--> Dneo_avg_<sp> [m^2/s], avg. neocl. diffusion from NEOART
! - Vneo_strahl(:,isp)--> Vneo_avg_<sp> [m/s], avg. neocl. convection from NEOART
!  
!============================================================================================!

use pi_const, only: GP, GP2
use scalars, only: TIME, TSTART, TAUPRP, NA1, PSIAX, RTOR, NA, HRO, IPART
use status, only: FP, UPL, VOLUM, SHIF, NE, TE, TI, AMAIN, ZMAIN
use read_input, only: machine, awd, nml_file, astra_ext
use numerical_tools, only: qinterp

implicit none

integer, parameter :: ngmax=140, n_o_max=512, nch_r1=31, nch_r2=32, nch_w1=41, nch_w4=44

character(len=120), parameter :: strahl_output='results.txt', strahl_param_in='param_files/sparams.dat'

double precision, intent(in) :: tau_start, zneocl, dzneocl, dimpsol, shot_in

integer :: i, isp, j, k, shotn, indexx, nimp_touse, nfour_c, n_grids, &
    Nr_o, ineocl, ineocla, ineocli, i_stepst, ios

integer, dimension(10) :: irecycl

real*8 :: tneocl, tneocl0, rneocl

double precision :: dum1, tau_strahl, &
    ne_decayl, te_decayl, ti_decayl, z_K, zdr_0, zdr_1, rbrlcfs, rlimrlcfs, &
    tolimiter, solflow, solrout1, solrout2, solrout3, solrout4, todivert, addsheathvoltage

double precision, dimension(NA1) :: rhopol, rhovol, r_rho

double precision, dimension(10) :: aweight, eneutr, rsources, rrates, trates, ridecay, &
    wrecycl, divpuff, swincm, swoutcm, promptredep, taudiv, taupump

double precision, dimension(n_o_max) :: rpol_o, zeff_o, pradtot_o, nmain_o, pradmain_o 

double precision, dimension(n_o_max, 10) :: pradsp_o, nimpsp_o, zavgsp_o, nesrcsp_o, &
     dneosp_o, vneosp_o

double precision, dimension(ngmax) :: rhopolg, neg, teg, tig

double precision, dimension(NA1,10) :: Dzin, Vzin, Dz_anom, Vz_anom

character(len=160) :: strahl_dir, cmd_cmd, as_nml
character(len=20) :: rho_coord, elements_touse(10)
character(len=6) :: diffname1_s

data i_stepst /0/
data ineocli /0/
data tneocl0 /0./
save i_stepst, tneocl0, ineocli

NAMELIST / strahl_par /  tau_strahl, rho_coord, ne_decayl, te_decayl, ti_decayl, &
    nfour_c, nimp_touse, elements_touse, aweight, eneutr, ridecay, irecycl, &
    wrecycl, diffname1_s, z_K, n_grids, zdr_0, zdr_1, rsources, rrates, trates, &
    ineocla, rneocl, divpuff, swincm, swoutcm, promptredep, taudiv, taupump, &
    rbrlcfs, rlimrlcfs, todivert, tolimiter, solflow, addsheathvoltage, &
    solrout1, solrout2, solrout3, solrout4

!--------------------------------------------------------------------------

as_nml = TRIM(awd) // '/' // TRIM(nml_file)
write(*, *) 'Reading namelist ', TRIM(as_nml)

open(nch_r1, FILE=TRIM(as_nml), delim='apostrophe')
read(nch_r1, nml=strahl_par, iostat=ios)
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
    zeff_strahl = 1.0
    prad_tot_strahl = 0.0
    nmain_strahl(1:NA1) = NE
    prad_main_strahl = 0.0
    prad_strahl = 0.0
    nimp_strahl = 0.0
    zavg_strahl = 1.0
    nesrc_strahl = 0.0
    Dneo_strahl = 0.0
    Vneo_strahl = 0.0
    return
endif

print *, "Calling STRAHL..."
print *, "Number of impurities:", nimp_touse
print *, "Species: ", (elements_touse(isp), isp=1, nimp_touse)

diffname1_s = ''

shotn = nint(shot_in)
strahl_dir = TRIM(awd) // '/strahl/'

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

do isp=1, nimp_touse
   Dzin(:,isp) = Dz_in_strahl(:,isp)*r_rho(:)
   Vzin(:,isp) = Vz_in_strahl(:,isp)*r_rho(:)
enddo

rhopolg(1) = 0.
do j=2, ngmax
    rhopolg(j) = rhopolg(j-1) + 1./(ngmax - 1.)
enddo

if (NA1 > ngmax) then
! Interpolate to astra grid
    do isp=1, nimp_touse
       call qinterp(rhopol(1:NA1), Dzin(1:NA1,isp), NA1, rhopolg(1:ngmax), Dz_anom(1:ngmax,isp), ngmax)
       call qinterp(rhopol(1:NA1), Vzin(1:NA1,isp), NA1, rhopolg(1:ngmax), Vz_anom(1:ngmax,isp), ngmax)
    enddo
    call qinterp(rhopol(1:NA1), NE(1:NA1), NA1, rhopolg(1:ngmax), neg(1:ngmax), ngmax)
    call qinterp(rhopol(1:NA1), TE(1:NA1), NA1, rhopolg(1:ngmax), teg(1:ngmax), ngmax)
    call qinterp(rhopol(1:NA1), TI(1:NA1), NA1, rhopolg(1:ngmax), tig(1:ngmax), ngmax)
else
   Dz_anom = 1.*Dzin
   Vz_anom = 1.*Vzin
endif

do isp=1,nimp_touse
   rrates(isp)=rrates_in_strahl(isp)
enddo

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
    'cv     finite vol=1, finite diff=0', &
    '         1', &
    '   ', & 
!    'cv      max. iterations at fixed time  stop iteration if change below(%)', &
!    '  1000          -1.		', &
!    '   ', &
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
do isp=1, nimp_touse
    write(nch_w1, '(3A, 2F12.4)') ' ', elements_touse(isp), '  ', aweight(isp), eneutr(isp)
enddo
write(nch_w1, '(A)') &
    '                            ', &
    '    		      S O  U R C E  ', &
    '   ', &
    'cv  r_source-r_lcfs(cm)   constant rate(1/s)    time dependent rate from file(1/0)'
do isp=1, nimp_touse
    write(nch_w1, 108) '  ', rsources(isp), '   ', rrates(isp), '  ', trates(isp)
enddo
write(nch_w1, '(/A)') &
    'cv    divertor puff    source_width_in(cm)     source_width_out(cm)   prompt redep'
do isp=1, nimp_touse
    write(nch_w1, 112) '	  ', divpuff(isp), '   ', swincm(isp), '  ', swoutcm(isp), '   ', promptredep(isp)
enddo
write(nch_w1, '(A)') &
    '', &
    '                    E D G E ,   R E C Y C L I N G', &
    '', &
    'cv    decay length of  impurity outside last grid point(cm)'
do isp=1, nimp_touse
    write(nch_w1, *) '                           ', ridecay(isp), '           '
enddo
write(nch_w1, '(A)') &
    '                          ', &
    '', &
    'cv    Rec.:ON=1/OFF=0    wall-rec.  Tau-div->SOL(ms)    Tau-pump(ms) '

! reclying has to be for each imp species
do isp=1, nimp_touse
    write(nch_w1, '(A, I, A, F12.4, A, E14.5, A, E12.4)') &
        '   ', irecycl(isp), '               ', wrecycl(isp), '        ', taudiv(isp), '  ', taupump(isp)
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
    "        'minter'", &
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

do isp=1,nimp_touse
   do i=1, min(NA1, ngmax)-1
      write(nch_w1, 101) '     ', Dz_anom(i,isp)
   enddo
   write(nch_w1, 101) '     ', Dz_anom(min(NA1, ngmax) - 1,isp) + (1. - rhopol(min(NA1, ngmax) - 1)) * &
        (dimpsol - Dz_anom(min(NA1, ngmax) - 1,isp))/(1.05 - rhopol(min(NA1, ngmax) - 1))
   write(nch_w1, 101) '     ', dimpsol
   write(nch_w1, 101) '     ', dimpsol
   write(nch_w1, 101) '     ', dimpsol
enddo

write(nch_w1, '(A)') &
    '               ', &
    '', &
    '', &
    'cv    Drift function        Drift Parameter/Velocity', &
    "      'minter'                 'velocity'   ", &
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
do isp=1,nimp_touse
   do i=1, min(NA1, ngmax)
      write(nch_w1, 101) '     ', Vz_anom(i,isp)
   enddo
   write(nch_w1, 101) '     ', 0.0
   write(nch_w1, 101) '     ', 0.0
   write(nch_w1, 101) '     ', 0.0
   write(nch_w1, 101) '     ', 0.0
enddo

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

101 format(A, F15.8)
104 format(A, F15.8 , A, I, A, F15.8, A, F15.8)
105 format(A, E25.11, A, E25.11, 10A)
108 format(A, F15.8 , A, E16.8, A, F15.8)
109 format(A, I, A , I, A, I, A)
112 format(A, F15.8 , A, F15.8, A, F15.8, A, F15.8)
145 format(A, F15.8 , A, F15.8, A, F15.8, A, F15.8, A, F15.8)
! end call params_file_write_strahl

call profiles_file_write_strahl(strahl_dir, &
    rho_coord, rhopol(1:NA1), &
    NA1, NE(1:NA1), TE(1:NA1), TI(1:NA1), &
    ne_decayl, te_decayl, ti_decayl, TIME-TIME, &
    rhopolg(1:ngmax), teg(1:ngmax) , neg(1:ngmax), tig(1:ngmax), ngmax)

call grid_write_strahl(strahl_dir, nfour_c, RTOR+SHIF(1), &
    rhovol(NA1), UPL(NA1), TIME-TIME, machine)

! Main STRAHL call

cmd_cmd = TRIM(astra_ext) // '/strahl/sep23/bin/strahl a q'

write(*, '(A)') 'Executing', cmd_cmd 
call system(cmd_cmd)      ! run strahl

cmd_cmd = 'rm -f results.txt'
call system(cmd_cmd)    ! rm old results, if existing

cmd_cmd = TRIM(astra_ext) // '/strahl/sep23/bin/result_to_astra '// TRIM(elements_touse(1))//' > ' // TRIM(strahl_dir) // 'results.txt'

write(*, '(A)') 'Executing', cmd_cmd 
call system(cmd_cmd)   ! produce new result file

call chdir(TRIM(awd))      ! cdir

! Extract results from stahl/results.txt:

open(nch_r2, file=TRIM(strahl_dir)//TRIM(strahl_output))
    read(nch_r2, *) Nr_o
    read(nch_r2, *) cmd_cmd
    read(nch_r2, *) (rpol_o(i), i=1, Nr_o) ! rho poloidal
    read(nch_r2, *) cmd_cmd
    read(nch_r2, *) (zeff_o(i), i=1, Nr_o) ! effective charge profile
    read(nch_r2, *) cmd_cmd
    read(nch_r2, *) (pradtot_o(i), i=1, Nr_o) ! total radiated power
    read(nch_r2, *) cmd_cmd
    read(nch_r2, *) (nmain_o(i), i=1, Nr_o) ! main ion density
    do isp=1, nimp_touse
        read(nch_r2, *) cmd_cmd
        read(nch_r2, *) (pradsp_o(i, isp), i=1, Nr_o) ! radiated power of each impurity
    enddo
    read(nch_r2, *) cmd_cmd
    read(nch_r2, *) (pradmain_o(i), i=1, Nr_o)  ! radiated power of main plasma
    do isp=1, nimp_touse
        read(nch_r2, *) cmd_cmd
        read(nch_r2, *) (nimpsp_o(i, isp), i=1, Nr_o) ! density of each impurity
    enddo
    do isp=1, nimp_touse
        read(nch_r2, *) cmd_cmd
        read(nch_r2, *) (zavgsp_o(i, isp), i=1, Nr_o) ! average charge of each impurity
    enddo
    do isp=1, nimp_touse
        read(nch_r2, *) cmd_cmd
        read(nch_r2, *) (dneosp_o(i, isp), i=1, Nr_o) ! NEOART avg. diff. coeff. of each impurity
    enddo
    do isp=1, nimp_touse
        read(nch_r2, *) cmd_cmd
        read(nch_r2, *) (vneosp_o(i, isp), i=1, Nr_o) ! NEOART avg. convection of each impurity
    enddo
    do isp=1, nimp_touse
        read(nch_r2, *) cmd_cmd
        read(nch_r2, *) (nesrcsp_o(i, isp), i=1, Nr_o) ! electron source of each impurity
    enddo
close(nch_r2)


! Interpolate to astra grid

! species-independent quantities
! Zeff
call qinterp(rpol_o(1:Nr_o), max(1., zeff_o(1:Nr_o)), Nr_o, &
     rhopol(1:NA1), zeff_strahl(1:NA1), NA1)

! Prad tot
call qinterp(rpol_o(1:Nr_o), pradtot_o(1:Nr_o)/1.E6, Nr_o, &
    rhopol(1:NA1), prad_tot_strahl(1:NA1), NA1)

! nmain
call qinterp(rpol_o(1:Nr_o), 1.E-19*max(nmain_o(1:Nr_o),0.0), Nr_o, &
    rhopol(1:NA1), nmain_strahl(1:NA1), NA1)

! Prad main
call qinterp(rpol_o(1:Nr_o), max(pradmain_o(1:Nr_o),0.0)/1.E6, Nr_o, &
    rhopol(1:NA1), prad_main_strahl(1:NA1), NA1)

! species-dependent quantities
! Prad species
do isp=1, nimp_touse
    call qinterp(rpol_o(1:Nr_o), max(pradsp_o(1:Nr_o, isp),0.0)/1.E6, Nr_o, &
         rhopol(1:NA1), prad_strahl(1:NA1, isp), NA1)
enddo

! nimp species
do isp=1, nimp_touse
    call qinterp(rpol_o(1:Nr_o), 1.E-19*max(nimpsp_o(1:Nr_o, isp),0.0), Nr_o, &
         rhopol(1:NA1), nimp_strahl(1:NA1, isp), NA1)
enddo

! zavg species
do isp=1, nimp_touse
    call qinterp(rpol_o(1:Nr_o), max(zavgsp_o(1:Nr_o, isp),0.0), Nr_o, &
         rhopol(1:NA1), zavg_strahl(1:NA1, isp), NA1)
enddo
 
! ne source species
do isp=1, nimp_touse
    call qinterp(rpol_o(1:Nr_o), 1.E-19*max(nesrcsp_o(1:Nr_o, isp),0.0), Nr_o, &
         rhopol(1:NA1), nesrc_strahl(1:NA1, isp), NA1)
enddo

!conversion from strahl rvol to astra rho

do isp=1,nimp_touse
   dneosp_o(1:NA1,isp) = dneosp_o(1:NA1,isp)/r_rho(1:NA1)
   vneosp_o(1:NA1,isp) = vneosp_o(1:NA1,isp)/r_rho(1:NA1)
enddo

! Dneo species
do isp=1, nimp_touse
    call qinterp(rpol_o(1:Nr_o), dneosp_o(1:Nr_o, isp), Nr_o, &
         rhopol(1:NA1), Dneo_strahl(1:NA1, isp), NA1)
enddo

! Vneo species
do isp=1, nimp_touse
    call qinterp(rpol_o(1:Nr_o), vneosp_o(1:Nr_o, isp), Nr_o, &
         rhopol(1:NA1), Vneo_strahl(1:NA1, isp), NA1)
enddo

i_stepst = i_stepst+1

if (i_stepst > 5) i_stepst=5

print *, "Finished STRAHL call"

201 format(1F15.8)
203 format(1I8)

end subroutine a2strahl

!----------------------------------------------------------------------
subroutine profiles_file_write_strahl(strahl_dir, rho_coord, rhopol, &
    ngrid, ne, te, ti, ne_decayl, te_decayl, ti_decayl, time, &
    rhopolg, teg, neg, tig, ngmax)

implicit none

integer, parameter :: nch_w2=42
character(len=14), parameter :: strahl_prof_in='nete/pp11111.1'

integer, intent(in) :: ngrid, ngmax
double precision, intent(in) :: time,  ne_decayl, te_decayl, ti_decayl
double precision, intent(in), dimension(ngrid) :: ne, te, ti, rhopol
double precision, intent(in), dimension(ngmax) :: rhopolg, teg, tig, neg
character(len=20), intent(in) :: rho_coord
character(len=80), intent(in) :: strahl_dir

integer :: i

!Start with main parameter file
open(nch_w2, file=TRIM(strahl_dir)//strahl_prof_in)
write(nch_w2, '(A)') &
    '          ******************** ', &
    '          **** from ASTRA **** ', &
    '          ******************** ', &
    '   ', &
    'cv    time-vector  ', &
    '   1'
write(nch_w2, 103) '   ', time
write(nch_w2, '(A)') &
    '   ', &
    '   ', &
    'cv      Function', &
    "      'interp'   ", &
    '   ', &
    '   ', &
    'cv    x-coordinate', &
    '      ' // rho_coord, &
    '   ', &
    '   ', &
    'cv   # of radial points'
write(nch_w2, *) '       ', min(ngmax, ngrid)
write(nch_w2, '(A)') &
    '   ', &
    '   ', &
    'cv   x-grid'
do i=1, min(ngmax, ngrid)
    if (ngmax < ngrid) then
        write(nch_w2, 101) rhopolg(i)
    else
        write(nch_w2, 101) rhopol(i)
    endif
enddo
write(nch_w2, '(A)') &
    '   ', &
    '   ', &
    'cv ELECTRON DENSITY (cm**-3)'
if (ngmax < ngrid) then
    write(nch_w2, 107) neg(1)*1.E+13
    do i=1, min(ngmax, ngrid)
        write(nch_w2, 101) neg(i)/neg(1)
    enddo
else
    write(nch_w2, 107) ne(1)*1.E+13
    do i=1, ngrid
        write(nch_w2, 101) ne(i)/ne(1)
    enddo
endif
write(nch_w2, '(A)') &
    '   ', &
    '   ', &
    'cv ne decay length [cm] in rho_volume'
write(nch_w2, 105) '   ', ne_decayl
write(nch_w2, '(A)') &
    '   ', &
    '   ', &
    'cv    time-vector  ', &
    '   1'
write(nch_w2, 103) '   ', time
write(nch_w2, '(A)') &
    '   ', &
    '   ', &
    'cv      Function', &
    "      'interp'   ", &
    '   ', &
    '   ', &
    'cv    x-coordinate', &
    "      " // rho_coord, &
    '   ', &
    '   ', &
    'cv   # of radial points'
write(nch_w2, *) '       ', min(ngmax, ngrid)
write(nch_w2, '(A)') &
    '   ', &
    '   ', &
    'cv   x-grid'
do i=1, min(ngmax, ngrid)
    if (ngmax < ngrid) then
        write(nch_w2, 101) rhopolg(i)
    else
        write(nch_w2, 101) rhopol(i)
    endif
enddo
write(nch_w2, '(A)') &
    '   ', &
    '   ', &
    'cv ELECTRON TEMPERATURE (eV)'
if (ngmax < ngrid) then
    write(nch_w2, 107) teg(1)*1.E+3
    do i=1, min(ngmax, ngrid)
        write(nch_w2, 101) teg(i)/teg(1)
    enddo
else
    write(nch_w2, 107) te(1)*1.E+3
    do i=1, ngrid
        write(nch_w2, 101) te(i)/te(1)
    enddo
endif
write(nch_w2, '(A)') &
    '   ', &
    '   ', &
    'cv te decay length [cm] in rho_volume'
write(nch_w2, 105) '   ', te_decayl
write(nch_w2, '(A)') &
    '   ', &
    '   ', &
    'cv    time-vector  ', &
    '   1'
write(nch_w2, 103) '   ', time
write(nch_w2, '(A)') &
    '   ', &
    '   ', &
    'cv      Function', &
    "      'interp'   ", &
    '   ', &
    '   ', &
    'cv    x-coordinate', &
    "      " // rho_coord, &
    '   ', &
    '   ', &
    'cv   # of radial points'
write(nch_w2, *) '       ', min(ngmax, ngrid)
write(nch_w2, '(A)') &
    '   ', &
    '   ', &
    'cv   x-grid'
do i=1, min(ngmax, ngrid)
    if (ngmax < ngrid) then
        write(nch_w2, 101) rhopolg(i)
    else
        write(nch_w2, 101) rhopol(i)
    endif
enddo
write(nch_w2, '(A)') &
    '   ', &
    '   ', &
    'cv ION TEMPERATURE (eV)'
if (ngmax < ngrid) then
    write(nch_w2, 107) tig(1)*1.E+3
    do i=1, min(ngmax, ngrid)
        write(nch_w2, 101) tig(i)/tig(1)
    enddo
else
    write(nch_w2, 107) ti(1)*1.E+3
    do i=1, ngrid
        write(nch_w2, 101) ti(i)/ti(1)
    enddo
endif
write(nch_w2, '(A)') &
    '   ', &
    '   ', &
    'cv ti decay length [cm] in rho_volume'
write(nch_w2, 105) '   ', ti_decayl
close(nch_w2)

101 format(F15.8)
103 format(A, F15.8, A, F15.8, A, F15.8, A, F15.8)
105 format(A, F15.8)
107 format(E16.8)

end subroutine profiles_file_write_strahl

!----------------------------------------------------------------------
subroutine grid_write_strahl(strahl_dir, nfour_c, Raxis, &
  Rvoltot, Vloop, time, machine)

use parameters_a2equil, only: equil_now

implicit none

integer, parameter :: nch_r3=33, nch_w3=43, nr_max=100
character(len=17) :: strahl_grid_in='nete/grid_11111.1'

integer, intent(in) :: nfour_c
double precision, intent(in) :: Raxis, Rvoltot, Vloop, time
character(len=80), intent(in) :: strahl_dir
character(len=4) , intent(in) :: machine

integer :: i, j, nequil
double precision :: rpol, rvol, sigma, dum1, bdb0, bdb02, bmaxt, btor, fc

nequil = SIZE(equil_now%profiles_1d%volume)

open(nch_w3, file=TRIM(strahl_dir)//strahl_grid_in)
write(nch_w3, '(A)') &
    '   ', &
    'cv  rho volume(LCFS)[cm]  R_axis[cm]   U_loop[V]    time[s] '
write(nch_w3, 103) '  ', Rvoltot*100., '  ', Raxis*100., '   ', Vloop, '  ', time
write(nch_w3, '(A)') &
    '   ', &
    '   ', &
    'cv  number of grid points  points up to separtrix  fourier coefficients'
write(nch_w3, *) ' ', nequil, '  ', nequil, '  ', nfour_c
write(nch_w3, '(A)') &
    '   ', &
    '   ', &
    'cv  sqrt( (Psi-Psi_ax) / (Psi_sep - Psi_ax) )   '
do i=1, nequil
    rpol = ((equil_now%profiles_1d%psi(i)      - equil_now%profiles_1d%psi(1)) / &
            (equil_now%profiles_1d%psi(nequil) - equil_now%profiles_1d%psi(1)))**0.5
    write(nch_w3, 101) rpol
enddo
write(nch_w3, '(A)') &
    '   ', &
    '   ', &
    'cv     rho volume / rho_volume(LCFS)    '

do i=1, nequil
    rvol = sqrt(equil_now%profiles_1d%volume(i)/equil_now%profiles_1d%volume(nequil))
    write(nch_w3, 101) rvol
enddo
write(nch_w3, '(A)') &
    '   ', &
    '   ', &
    'cv   large radius low field side / R_axis '
do i=1, nequil
    write(nch_w3, 101) equil_now%profiles_1d%r_outboard(i)/equil_now%profiles_1d%r_outboard(1)
enddo
write(nch_w3, '(A)') &
    '   ', &
    '   ', &
    'cv   large radius high field side / R_axis '
do i=1, nequil
    write(nch_w3, 101) equil_now%profiles_1d%r_inboard(i)/equil_now%profiles_1d%r_inboard(1)
enddo
write(nch_w3, '(A)') &
    '   ', &
    '   ', &
    'cv  safety factor '
do i=1, nequil
    write(nch_w3, 101) equil_now%profiles_1d%q(i)
enddo
write(nch_w3, '(A)') &
    '   ', &
    '   ', &
    'cv  fraction of circulating particles  '

btor = equil_now%global_param%toroid_field%b0
do i=1, nequil
    BDB0  = equil_now%profiles_1d%bdb0(i)
    BMAXT = equil_now%profiles_1d%bmaxt(i)
    BDB02 = equil_now%profiles_1d%gm4(i)
    dum1  = min(.999999d0, BDB0/BMAXT)
    sigma = 1. - BDB02/BDB0**2 * (1. - SQRT(1. - dum1)*(1. + 0.5*dum1))
    dum1  = 1. - (BDB02/btor**2)*equil_now%profiles_1d%fofb(i)
    fc = 1. - 0.75*sigma - 0.25*dum1
    if (i == 1) fc = 1.
    write(nch_w3, 101) fc
enddo
write(nch_w3, '(A)') &
    '   ', &
    '   ', &
    'cv  Integral( dl_p / B_p) [m/T]  '
do i=1, nequil
    write(nch_w3, 101) 1./equil_now%profiles_1d%dPSIdV(i)
enddo
write(nch_w3, '(A)') &
    '   ', &
    '   ', &
    'cv  < B_total > [T]  '
do i=1, nequil
    write(nch_w3, 101) equil_now%profiles_1d%bdb0(i)
enddo
write(nch_w3, '(A)') &
    '   ', &
    '   ', &
    'cv  < B_total**2 > [T**2]  '
do i=1, nequil
    write(nch_w3, 101) equil_now%profiles_1d%gm4(i)
enddo
write(nch_w3, '(A)') &
    '   ', &
    '   ', &
    'cv  < 1./B_total**2 > [1/T**2]  '
do i=1, nequil
    write(nch_w3, 101) equil_now%profiles_1d%gm5(i)
enddo
write(nch_w3, '(A)') &
    '   ', &
    '   ', &
    'cv  < R B_T > [m*T]  '
do i=1, nequil
    write(nch_w3, 101) equil_now%profiles_1d%F_dia(i)
enddo
write(nch_w3, '(A)') &
    '   ', &
    '   ', &
    'cv  < R**2 B_p**2/B**2 > [m**2]  '
do i=1, nequil
    write(nch_w3, 101) equil_now%profiles_1d%rbp_b2(i)
enddo
write(nch_w3, '(A)') &
    '   ', &
    '   ', &
    'cv  < 1/R**2 > [1/m**2]  '
do i=1, nequil
    write(nch_w3, 101) equil_now%profiles_1d%gm1(i)
enddo
write(nch_w3, '(A)') &
    '   ', &
    '   ', &
    'cv  <cos (m theta) B_total**2> [T**2]'
do j=1, nfour_c
    do i=1, nequil
        write(nch_w3, 101) equil_now%profiles_1d%acosB2a(i, j)
    enddo
enddo
write(nch_w3, '(A)') &
    '   ', &
    '   ', &
    'cv  <sin (m theta) B_total**2> [T**2]'
do j=1, nfour_c
    do i=1, nequil
        write(nch_w3, 101) equil_now%profiles_1d%asinB2a(i, j)
    enddo
enddo
write(nch_w3, '(A)') &
    '   ', &
    '   ', &
    'cv  <cos (m theta) B_total ln(B_total)> [Tln(T)]'
do j=1, nfour_c
    do i=1, nequil
        write(nch_w3, 101) equil_now%profiles_1d%acosBlnBa(i, j)
    enddo
enddo
write(nch_w3, '(A)') &
    '   ', &
    '   ', &
    'cv  <sin (m theta) B_total ln(B_total)> [Tln(T)]'
do j=1, nfour_c
    do i=1, nequil
        write(nch_w3, 101) equil_now%profiles_1d%asinBlnBa(i, j)
    enddo
enddo
write(nch_w3, '(A)') &
    '   ', &
    '   ', &
    'cv  Bp at LFS equator [T]  '
do i=1, nequil
    write(nch_w3, 101) equil_now%profiles_1d%bplfs(i)
enddo
close(nch_w3)

101 format(F15.8)
103 format(A, F15.8, A, F15.8, A, F15.8, A, F15.8)

end subroutine grid_write_strahl

end module strahl_mod
