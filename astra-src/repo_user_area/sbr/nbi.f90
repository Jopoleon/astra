subroutine NBI()
!----------------------------------------------------------------------|
! Interface to the Neutral Beam Injection package by A.R.Polevoi
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
use nbstatus, only: set_input, get_output
use nb_injection, only: nbinj

implicit none

integer :: n_nbi, cx_flag, fp_flag, cx_cold, calc_fus, dn_rho

if (NBFILE(1:1).eq.'*') then
    write(*, *)'>>> NBI Error >>> Configuration file not found'
    stop
endif

! CNB1  = 8
! CNB2  = 1   ? Explicit form of CX losses ?
call set_input(NE, NHYDR, NDEUT, NTRIT, NHE3, &
    NALF, NI, NIZ1, NIZ2, NIZ3, ZIM1, ZIM2, ZIM3, TE, TI, &
    VR, SHIF, SHIV, ELON, TRIA, AMETR, RHO, FP, MU, AMAIN, &
    NN, TN, ZEF, G33, IPOL, NIBM, PIBM, PEBM, PBLON, PBPER, & 
    PBEAM, SNEBM, SNNBM, CUFI, CUBM, SCUBM, &
    SNIBM1, SNIBM2, SNIBM3, NNBM1, NNBM2, NNBM3)

n_nbi    = nint(CNB1)
cx_flag  = nint(CNB2)
fp_flag  = nint(CNB4)
cx_cold  = nint(CNBI1)
dn_rho   = nint(CNB3)
calc_fus = nint(CNBI4)
call NBINJ( trim(NBFILE), BTOR, RTOR, ABC, AB, ROC, SHIFT, UPDWN, &
    HRO, TAU, NA1, NB1, AIM1, AIM2, AIM3, AMJ, ZMJ, NNCL, NNWM, QNBI, &
! Any allowed control parameters can be used below
    n_nbi  , & ! No. of NB sources |CNB1| if<0 no FP, srs only     (8)
    cx_flag, & ! bulk-ion CX losses:  1/0 on/off                   (1)
    dn_rho , & ! NBI space grid: N1=(NA1-1)/CNB3+1                 (2)
    fp_flag, & ! 0/1/2 no solver/steady st. FP/time dep. FP solver (1)
! the next 3 parameters are effective if CNB4=2 only
    cx_cold, & ! fast-ion CX losses due to cold neutrals
    CNBI2  , & ! FP solver time step: TAUNBI=CNBI2*TAU
    CNBI3  , & ! fast-ion CX losses due to NBI  neutrals
    calc_fus)  ! control of beam-plasma fusion 1/0/2 full/off/Ti=0
!  0 (no fusion), 1 (finite Te, Ti), 2 (finite Te, Ti=0) 

call get_output(NIBM, PIBM, PEBM, PBLON, PBPER, & 
    PBEAM, SNEBM, SNNBM, CUFI, CUBM, SCUBM, &
    SNIBM1, SNIBM2, SNIBM3, NNBM1, NNBM2, NNBM3)

end subroutine nbi
