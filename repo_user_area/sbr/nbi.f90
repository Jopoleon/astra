subroutine NBI
!----------------------------------------------------------------------|
! Interface to the Neutral Beam Injection package by A.R.Polevoi
!      (Edition 19-APR-2000)
! 
! This interface is compatible with the Astra version 5.3 and later
!      (Pereverzev 27-09-01)
!----------------------------------------------------------------------|
! The subroutine NBI is called from ASTRA model and provides
! the multisource neutral beam power, momentum and driven current
! deposition including ripple losses.
! Optionally, the time dependent Fokker-Planck solver (NBIBCE.f) 
!       can be used
!----------------------------------------------------------------------|
! Input:
! The input is split in three parts. (See also description at the end)
!   (1) The tokamak parameters:
!      (i)  The NBI subroutine input parameters. These parameters: 
!           CNB1, CNB2, CNB3, CNB4, CNBI1, CNBI2, CNBI3, CNBI4 either
!           have to be set explicitly in this subroutine (see
!           an example below) or by the data file "exp/DATA_FILE"
!           or in the calling model. 
!      (ii) Parameters to user's function RIPRAD (sbr/nbuser.f).
!   (2) Beam box parameters:
!      (i)  The input via the file FILENA (the file name is
!           constructed as "exp/DATA_FILE.nbi"), 
!      (ii) User's functions NBFRY, NBFHZ (see file sbr/nbuser.f).
!   (3) Input plasma profiles as TE, TI, NE and so on are transferred 
!       through common blocks.
!----------------------------------------------------------------------|
! Output: QNBI [MW] NBI power
!  PBEAM, PEBM, PIBM, NIBM, CUFI, CUBM, PBLON, PBPER, 
!  SCUBM, SNEBM, SNNBM, NNBM1, 2, 3 for MAIN
!----------------------------------------------------------------------|

use scalars
use status
use read_input, only: NBFILE
use nbstatus, only: nbstatus_io

implicit none

integer :: JINOUT

if (NBFILE(1:1).eq.'*') then
    write(*, *)'>>> NBI Error >>> Configuration file not found'
    stop
endif

! CNB1  = 8
! CNB2  = 1   ? Explicit form of CX losses ?
 JINOUT=0 ! ASTRA->NBI
 call nbstatus_io(JINOUT, NB1, NA1, NE, NHYDR, NDEUT, NTRIT, NHE3, &
     NALF, NI, NIZ1, NIZ2, NIZ3, ZIM1, ZIM2, ZIM3, TE, TI, &
     VR, SHIF, SHIV, ELON, TRIA, AMETR, RHO, FP, MU, AMAIN, &
     NN, TN, ZEF, G33, IPOL, NIBM, PIBM, PEBM, PBLON, PBPER, & 
     PBEAM, SNEBM, SNNBM, CUFI, CUBM, SCUBM, &
     SNIBM1, SNIBM2, SNIBM3, NNBM1, NNBM2, NNBM3)

 call NBINJ( &
     trim(NBFILE), BTOR, RTOR, ABC, AB, ROC, SHIFT, UPDWN, &
    HRO, TAU, NA1, NB1, AIM1, AIM2, AIM3, AMJ, ZMJ, NNCL, NNWM, QNBI, &
! Any allowed control parameters can be used below
    CNB1, & ! No. of NB sources |CNB1| if<0 no FP, srs only  (8)
    CNB2, & ! bulk-ion CX losses:  1/0 on/off                (1)
    CNB3, & ! NBI space grid: N1=(NA1-1)/CNB3+1              (2)
    CNB4, & ! 0/1/2 no solver/steady st. FP/time dep. FP solver (1)
! the next 3 parameters are effective if CNB4=2 only&
    CNBI1, & ! fast-ion CX losses due to cold neutrals&
    CNBI2, & ! FP solver time step: TAUNBI=CNBI2*TAU &
    CNBI3, & ! fast-ion CX losses due to NBI  neutrals&
    CNBI4)   ! control of beam-plasma fusion 1/0/2 full/off/Ti=0
!  0 (no fusion), 1 (finite Te, Ti), 2 (finite Te, Ti=0) 

JINOUT=1 ! NBI->ASTRA

call nbstatus_io(JINOUT, NB1, NA1, NE, NHYDR, NDEUT, NTRIT, NHE3, &
    NALF, NI, NIZ1, NIZ2, NIZ3, ZIM1, ZIM2, ZIM3, TE, TI, &
    VR, SHIF, SHIV, ELON, TRIA, AMETR, RHO, FP, MU, AMAIN, &
    NN, TN, ZEF, G33, IPOL, NIBM, PIBM, PEBM, PBLON, PBPER, & 
    PBEAM, SNEBM, SNNBM, CUFI, CUBM, SCUBM, &
    SNIBM1, SNIBM2, SNIBM3, NNBM1, NNBM2, NNBM3)

end subroutine nbi

!======================================================================|
!  Control parameters: (as described by A.R.Polevoi)
! CNB1 abs(CNB1) is  number of PINIs with different geometry.  (Def=1)
! When CNB1 is negative or changes its value during a run then 
! a pop-up menu appears for interactive setting PINI parameters
! 
!
! CNB2 Control of CX bulk-ion heat losses due to NBI  (Def=1)
! CNB2 > 0 (.and. CNB4.eq.1) explicit treatment of CX heat losses
! that is formula PBICX is included in PIBM. This corresponds 
! to the explicit calculations of the heat sink from the bulk 
! ion component at the moment when NBI is called.
! = 0 No losses from bulk ions due to NBI charge exchange
! This setting supposes that the implicit form is used:
!    CNB2=0; PIT=...-.0024*SNNBM
! CNB3 =1, 2... NBI X-mesh points number N1=(NA1-1)/CNBI3+1
! NOTE: if the fast ion gyroradius Rg(EBEAM, ABEAM)
! exceeds the Astra space step, i.e. Rg > ABC/NA1, 
! then the internal NBI solver space step will inrcease 
! automatically: CNB3 = Rg*NA1/ABC: NA1/CNB3, max = 5
!
!  Fast ion solver control:
! CNBI1, 2, 3, 4 are active for CNB4=2 only
! CNB4 =0/1/2 fast ion source only/ steady /time dependent Fokker-Planck Solver (Def=1)
!
! CNBI1 =1 losses (=0 No losses) from fast ions due to CX (Def=1)
! with the cold neutrals
!
! CNBI2 =TAUNBI/TAU ratio of the NBI time-step to TAU (Def=1)
! WARNING: TAUNBI does not coincide automatically with the 
! interval between sequental calls of NBI (that saves
! the calculation time at the NBI-steady state phase) 
! CNBI3 =1 losses (=0 No losses) from fast ions due to CX
! with the NBI neutrals itself
! The same CX cross-section as for cold neutrals is used.
!
!
! CNBI4 = control of beam-plasma fusion              (Def=1)
!       !  1/0/2 finite Te, Ti/off/Ti=0, finite Te 
!
!----------------------------------------------------------------------|
! Unlike the previous versions a number of parameters have to be
! defined in the input file exp/data_file_name.nbi (FILENA)
!       This file provides a description of each beam line and 
!       includes N groups, where N is the total number of beam lines.
!
! Each group consists of 21 records. 
! The first record is the ordinal group number in a separate line.
! The rest of a group consisits of 20 character*12 fields.
! Every field can be either a real number or a name of variable.
! The allowed variables are: 
!   "ZRDn" (or "ZRDnX"), where n stands for an integer number 1<=n<=48
!   "Astra_Constant" (see the full list in "for/const.inc" 
!                     or just press "C" in the run mode).
! 
! Significance of input parameters in each group:
!
! Ordinal_beam_number
! QBEAM  [MW] Beam power
! CONTR       Counter injection fraction (0 co, 1 counter)
! ABEAM [m_p] Mass of beam ions in the proton mass
! ZBEAM       Beam ion charge in the proton charge units
! EBEAM [keV] Beam energy
! DBM1        EBEAM power fraction
! DBM2        EBEAM/2 power fraction
! DBM3        EBEAM/3 power fraction
! Orb_av      Type of averaging over ion orbits
!    =0 No averaging (deposition at the birth point)
!    =1 Averaging with a finite orbit width
!        =2 Averaging with zero orbit width
! Penc_num    Number of pencils in the horizontal plane
! HBEAM   [m] Beam footprint center height
! RBMAX   [m] Beam footprint maximum radius
! RBMIN   [m] Beam footprint minimum radius
! tg(A)       A is the angle between the beam and the midplane
! Aspect      Footprint aspect ratio: Beam_height=Aspect*(RBMAX-RBMIN)
! Cver1       These parameters describe exponential (or any other)
! Cver2         beam power distribution across the beam cross-section
! Chor1         as described by the user functions NBFRY & NBFHZ
! Chor2         (see file sbr/nbuser.f)
!
! Parameter meaning:
! QBEAM [MW] NBI power      (Def=0)
! CONTR=0-1 fraction of contr-inj. power   (Def=1)
! ABEAM=1, 2, 3 Hydrogen, Deuterium, Tritium NBI only  (Def=1)
! EBEAM [keV] the main component energy   (Def=1)
! DBM1=0-1 ~ full energy component EBEAM power fraction  (Def=1)
! DBM2=0-1 ~ half energy component EBEAM/2 power fraction (Def=0)
! DBM3=0-1 ~third energy component EBEAM/2 power fraction (Def=0) 
! Warning: DBM1 > 0 must be used;    SUM(DBMi).ne.1. is allowed
! Power is renormalized withinNBI: Power(i)~DBMi/SUM(DBMi) 
! HBEAM [m] NBI footprint center height in respect to Z=0 (Def=0)
! RBMAX [m] max. major radius of the NBI footprint (Def=0)
! RBMIN [m] min. major radius of the NBI footprint (Def=0)
!
!omments: the footprint position is interpreted as:
!
! for tangential NBI: (RBMAX+RBMIN)/2 > RTOR - AB
! as NBI cross-section by the meridianal plane perpendicular 
! to the vertical plane of the central NBI pencil beam
! with tangential radius of (RBMAX+RBMIN)/2 
!
! for perpendicular NBI: (RBMAX+RBMIN)/2 < RTOR - AB 
! as NBI cross-section with the vertical plane perpendicular to 
! the central NBI pencil beam (RBMAX+RBMIN)/2 vertical plane, 
! placed at R = RTOR (the corrections, connected with 
! Shafranov's shift are calculated into the NBI block
!
! Penc_num  Each PINI is approxomated by Nh x Nr pencils (Def=1)
!  where Nr = Penc_num (horizontal direction), 
!  Nh (vertical direction) is automatically 
!  calculated as a number of magn. surfaces
!  in the vertical aperture of each PINI in 
!  the footprint crossection.
! 
! tg(A) tg(angle between central pencil and midplane) (Def=1)
!  is used for perpendicular NBI only  
! Aspect =(HBmax-HBmin)/(RBMAX-RBMIN): HBEAM=(HBmax+HBmin)/2(Def=1)
!
! User's functions
!
! RIPRAD(ZUPDWN, J) [m] - Ripple loss cone boundary (major 
!  radius) for each magnetic surface: ZUPDWN [m]- 
!  shift in respect to the midplane of ripple 
!  simmetry, J - surface index
! RIPRAD -depends on tokamak mag. field coils' and 
!  plasma configurations  (Def: No Ripple losses) (Def=999)
!======================================================================|
! Parameter mapping: 
!      External name (Astra_5.3) <-> Internal name (NBI and old Astra)
!   CNB1  <-> CBM1
!   CNB2  <-> CBM2
!   CNB3  <-> CBMI3
!   CNB4  <-> CBMI1
!   CNBI1  <-> CBM4
!   CNBI2  <-> CBMI2
!   CNBI3  <-> CBM3
!   CNBI4  <-> CBMI4
