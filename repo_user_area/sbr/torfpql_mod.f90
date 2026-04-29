module torfpql_mod

!***********************************************************************
!                                                                      *
!  Wrapper of TORIC-SSFPQL for ASTRA 8                                 *
!                                                                      *
!  Public routines                                                     *
!  ---------------                                                     *
!    toric : subroutine to call TORIC-SSFPQL                           *
!                                                                      *
!  Author/Contact                                                      *
!  --------------                                                      *
!  roberto.bilato@ipp.mpg.de                                           *
!                                                                      *
!  Todo                                                                *
!  ----                                                                *
!***********************************************************************
!
!-----------------------------------------------------------------------
! ASTRA modules
!-----------------------------------------------------------------------
use read_input, only: awd, nml_file
use status, only: nhydr, ndeut, ntrit, nhe3, nalf, &
    zim1, zim2, zim3, niz1, niz2, niz3, ne, te, ti, &
    fp, piicr, peicr, pifw, pefw, cufw, cuicr

use scalars, only: na1, aim1, aim2, aim3, psiax, psibo, roc
     
!-----------------------------------------------------------------------
! Use TORIC-SSFPQL modules
!-----------------------------------------------------------------------
use tor_mod_public, only:  r8,                                       &
     toricmode, inputpath, profnt_file, equil_file,                  &
     freqcy, nphi, hf_pwgoal, accur_goal,                            &
     max_iter, nspele, nspec, isol,                                  &
     lun5, lun6, lun8, lun9, lun17, lun21, lun22                    
use fpq_mod_public, only:                                            &
     nrdpsi, redpsi, pwnet_prof, tpwnet, fpwnet,                     &
     lun18
use tor_mod_deallocate, only: torfpql_deallocate_all
use tor_mod_nmlst, only: tor_open_nmlst, tor_close_nmlst, tor_run_nmlst

!-----------------------------------------------------------------------
! iso_c_binding for portable chdir
!-----------------------------------------------------------------------
use iso_c_binding, only: c_int, c_char, c_null_char

!-----------------------------------------------------------------------
! Module variables
! lun_scr  : logical unit for screen output
! tordir   : working TORIC-SSFPQL directory
! exedir   : directory where TORIC-SSFPQL is executed
! profnt   : filename for plasma profiles
! equigs   : filename for equilibrium
! rho_max  : maximum rho considered in TORIC
! npt_ts  : number of ASTRA radial points within rho_max
! rpol_as  : poloidal rho array of ASTRA
! cmd_str  : buffer for shell command strings
!
!-----------------------------------------------------------------------
implicit none

integer,          save :: lun_scr = 6
character(len=120)     :: tordir, exedir, cmd_str
character(len=20), parameter :: &
     profnt = 'profnt.dat',     &
     equigs = 'toric.eqdsk.equigs'
real(r8)               :: rho_max = 0.99_r8

integer                        :: npt_ts
real(r8), dimension(:), allocatable :: rpol_as

real(r8), parameter :: ths       = 1.e-5_r8

! Scale factor: ASTRA ne [10^19 m^-3] -> TORIC ne [cm^-3]
real(r8), parameter :: ne_scale = 1.0e13_r8

!-----------------------------------------------------------------------
! Derived type to describe a single plasma species (section 3.1)
!-----------------------------------------------------------------------
type :: species_t
   real(r8), pointer :: dens(:) => null()     ! pointer to ASTRA density array
   real(r8)          :: chg  = 0._r8          ! charge number
   real(r8)          :: mass = 0._r8          ! mass number
   real(r8)          :: conc = 0._r8
   character(5)      :: name = '     '
   logical           :: has_var_chg = .false. ! true for impurities with Z array
   real(r8), pointer :: chg_arr(:) => null()  ! pointer to ASTRA Z array
   real(r8)          :: pwic
end type species_t

 ! Species descriptors
integer, parameter :: nspec_max = 8
type(species_t) :: specs(nspec_max)

private

!***********************************************************************
! C interface for chdir
!***********************************************************************
  interface
     function c_chdir(path) bind(C, name='chdir') result(res)
       import :: c_int, c_char
       character(kind=c_char), intent(in) :: path(*)
       integer(c_int)                     :: res
     end function c_chdir
  end interface
  
!-----------------------------------------------------------------------
! Public procedures
!-----------------------------------------------------------------------
  public :: toric

contains

!=======================================================================
!=======================================================================
  subroutine toric(pwic, frq, ntor, toll, nmax, debug)
!
!***********************************************************************
!                                                                      *
!  Purpose                                                             *
!  -------                                                             *
!  Main driver: sets up directories, writes input files, runs          *
!  TORIC-SSFPQL, and broadcasts results back to ASTRA.                 *
!                                                                      *
!***********************************************************************
!
!-----------------------------------------------------------------------
! Optional input parameters
!-----------------------------------------------------------------------
! pwic   : [real: MW] ICRF coupled power          (default: from namelist)
! frq    : [real: Hz] RF frequency
! ntor   : [int] toroidal wave number
! toll   : [real] convergence tolerance        (default: from namelist)
! nmax   : [int]  max TORIC-SSFPQL iterations  (default: from namelist)
! debug  : [int]  0 = remove run directory on exit; 1 = keep it (default = 0)
!-----------------------------------------------------------------------
    use tor_mod_looptorql, only: tor_toricql

    external eqdsk

    real(r8), intent(in), optional :: pwic, frq, toll
    integer,  intent(in), optional :: ntor, nmax, debug

    integer :: ios
    
    real(r8) :: pwic_lc=-1._r8, frq_lc=-1._r8, toll_lc=-1._r8
    integer  :: ntor_lc=-99999, nmax_lc=-1

    if (present(pwic)) pwic_lc = pwic
    if (present(frq))  frq_lc  = frq
    if (present(toll)) toll_lc = toll
    if (present(ntor)) ntor_lc = ntor
    if (present(nmax)) nmax_lc = nmax
    
    write(lun_scr,1100)
1100 format(/,2X,'----- TORIC-SSFPQL starts ---------------------------------')


!-----------------------------------------------------------------------
! Generating the running directory
!-----------------------------------------------------------------------
    call mk_rundir

!-----------------------------------------------------------------------
! Move into the running directory
!-----------------------------------------------------------------------
    call safe_chdir(exedir, ios)

!-----------------------------------------------------------------------
! Namelist file
!-----------------------------------------------------------------------
    call write_nml(pwic_lc, frq_lc, ntor_lc, toll_lc, nmax_lc)

!-----------------------------------------------------------------------
! Plasma profiles
!-----------------------------------------------------------------------
    call profnt4toric

!-----------------------------------------------------------------------
! Equilibrium
!-----------------------------------------------------------------------
    call equigs4toric

!-----------------------------------------------------------------------
! Execute TORIC-SSFPQL loop
!-----------------------------------------------------------------------

!----- Open namelist file ----------------------------------------------    
    call tor_open_nmlst    

!----- Read the namelist for the type of run ---------------------------
    call tor_run_nmlst

!----- TORIC-SSFPQL loop manager ---------------------------------------
    call tor_toricql

!----- Close namelist file ----------------------------------------------
    call tor_close_nmlst

!-----------------------------------------------------------------------
! Broadcast TORIC-SSFPQL results into ASTRA arrays
!-----------------------------------------------------------------------
    call broadcast2astra

!-----------------------------------------------------------------------
! Clean up memory
!-----------------------------------------------------------------------
    call torfpql_deallocate_all

!-----------------------------------------------------------------------
! Return to home directory
!-----------------------------------------------------------------------
    call safe_chdir(awd, ios)

!-----------------------------------------------------------------------
! Remove run directory unless debug==1
!-----------------------------------------------------------------------
    if (present(debug)) then
       if (debug == 1) return   ! keep the directory
    end if
    call clean_up
    
      write(lun_scr,1100)
1101 format(2X,'----- TORIC-SSFPQL ends ---------------------------------',/)
  
    return
  end subroutine toric

!=======================================================================
!=======================================================================
  subroutine mk_rundir
!
!***********************************************************************
!                                                                      *
!  Purpose                                                             *
!  -------                                                             *
!  Create the base toric directory and a timestamped run subdirectory. *
!                                                                      *
!***********************************************************************

    character(len=8)  :: date_str   ! ccyymmdd
    character(len=10) :: time_str   ! hhmmss.sss

!----- Base directory --------------------------------------------------
    tordir = trim(awd) // '/toric'

    cmd_str = 'mkdir -p ' // trim(tordir)
    call run_cmd(cmd_str, 'mkdir tordir')

!----- Timestamped execution subdirectory ------------------------------
    call date_and_time(date=date_str, time=time_str)
    exedir = trim(tordir) // '/torfpql_' // date_str // '_' // time_str(1:6)

    cmd_str = 'mkdir -p ' // trim(exedir)
    call run_cmd(cmd_str, 'mkdir exedir')

    return
  end subroutine mk_rundir

!=======================================================================
!=======================================================================
  subroutine write_nml(pwicrf, frq, ntor, toll, nmax)
!
!***********************************************************************
!                                                                      *
!  Purpose                                                             *
!  -------                                                             *
!  Read the TORIC-SSFPQL namelist template, override with the supplied *
!  optional arguments, and write the updated namelist to ./torica.inp. *
!                                                                      *
!***********************************************************************
    use tor_mod_nmlst, only:                          &
         tor_run_nmlst, tor_loop_nmlst, tor_hf_nmlst, &
         tor_antenna_nmlst, fpq_fpql_nmlst

    real(r8), intent(in) :: pwicrf, frq, toll
    integer,  intent(in) :: ntor, nmax

    integer            :: ios
    character(len=120) :: as_nml

    namelist /toric/ freqcy, hf_pwgoal, accur_goal, max_iter, nphi
    
!----- Locate and parse namelist template ------------------------------
    as_nml = trim(awd) // '/exp/icrh/torica.inp'
    write(lun_scr, '(3A)') 'Parsing namelist ', trim(as_nml), ' for TORIC-SSFPQL'

    open(lun5, file=trim(as_nml), delim='apostrophe', iostat=ios)
    if (ios /= 0) then
       write(lun_scr, '(2A)') 'ERROR: namelist not found: ', trim(as_nml)
       return
    end if

    call tor_run_nmlst('r')
    call tor_loop_nmlst('r')
    call tor_hf_nmlst('r')
    call tor_antenna_nmlst('r')
    call fpq_fpql_nmlst('r')
    
    close(lun5)

!----- Read from the namelist of the experiment ------------------------
    as_nml = trim(awd) // '/' // trim(nml_file)
    open(lun5, file=trim(as_nml), delim='apostrophe', iostat=ios)
    rewind lun5
    read(lun5,nml=toric)
    close(lun5)

!----- Override with supplied values -----------------------------------
    if (pwicrf > 0._r8) hf_pwgoal  = pwicrf
    if (frq > 0._r8)    freqcy     = frq
    if (toll > 0._r8)   accur_goal = toll
    if (nmax > 0)       max_iter   = nmax
    if (ntor > -1000)   nphi       = ntor    
    toricmode   = 'toric_fpql'
    isol        = 1
    profnt_file = profnt
    equil_file  = equigs
    inputpath   = ''

!----- Write updated namelist ------------------------------------------
    ios = 0
    open(lun5, file=trim(exedir)//'/torica.inp', status='new', iostat=ios)
    
    if (ios /= 0) then
       write(lun_scr, *) 'ERROR: ', ios, lun5
       return
    end if
    call tor_run_nmlst('w')
    call tor_loop_nmlst('w')
    call tor_hf_nmlst('w')
    call tor_antenna_nmlst('w')    
    call fpq_fpql_nmlst('w')

    close(lun5)

  end subroutine write_nml

!=======================================================================
!=======================================================================
  subroutine broadcast2astra
!
!***********************************************************************
!                                                                      *
!  Purpose                                                             *
!  -------                                                             *
!  Broadcast TORIC-SSFPQL results to ASTRA arrays.                     *
!                                                                      *
!***********************************************************************
    use numerical_tools, only: qinterp
    use standard_functions, only: VINT

    integer :: i, isp
    real(r8) :: pic_as, pec_as, ptot_as, pic_ts, pec_ts, ptot_ts

!----- Initialized to zero ASTRA arrays --------------------------------
    pifw (1:na1) = 0._r8    ! ion heating from fast wave
    pefw (1:na1) = 0._r8    ! electron heating from fast wave
    cufw (1:na1) = 0._r8    ! fast-wave driven current
    piicr(1:na1) = 0._r8    ! ICRH ion power
    peicr(1:na1) = 0._r8    ! ICRH electron power
    cuicr(1:na1) = 0._r8    ! ICRH driven current

!----- Broadcast/Interpolate to ASTRA radial grid -----------------------
    call qinterp(redpsi *  rpol_as(npt_ts), sum(pwnet_prof(1:nspec, :), dim=1), &
         nrdpsi, rpol_as(1:npt_ts), piicr(1:npt_ts), npt_ts)
    call qinterp(redpsi *  rpol_as(npt_ts), pwnet_prof(nspele, :),              &
         nrdpsi, rpol_as(1:npt_ts), peicr(1:npt_ts), npt_ts)

!----- Update the species structure with their absorbed power -----------    
    i = 0
    do isp = 1, nspec_max
       if (specs(isp)%conc > ths) then
          i = i+1
          specs(isp)%pwic = tpwnet(i)
       endif
    end do
    
!----- Output to screen ---------------------------------------------------
    pic_as = VINT(PIICR, ROC)
    pec_as = VINT(PEICR, ROC)
    ptot_as = pic_as + pec_as

    pic_ts = sum(tpwnet(1:nspec))
    pec_ts = tpwnet(nspele)
    ptot_ts = pic_ts + pec_ts
    
    write(lun_scr,1100) 
    write(lun_scr,1105) ptot_as, ptot_ts, pec_as, &
         pec_as/ptot_as*100._r8, pic_as, pic_as/ptot_as*100._r8
    do isp = 1, nspec_max
       if (specs(isp)%conc > ths) then
          write(lun_scr,1110) trim(specs(isp)%name), specs(isp)%pwic, &
               specs(isp)%pwic/ptot_ts*100._r8
       endif
    enddo
    write(lun_scr,1100) 

1100 format(/,                                                                &
             2X,'--------------------------------------------------------',/)
1105 format( 2X,'Summary TORIC-SSFPQL (ICRF only)',/,                         &
             6X,'Total absorbed power ', 1P,E8.2,' MW, (',1P,E8.2,' MW)',/,   &
             4X,'Power Repartiton',/,                                         &
             6X,'Electrons  ',1P,E8.2,' MW (',0P,F5.2,' %)',/,                 &
             6X,'Ions       ',1P,E8.2,' MW (',0P,F5.2,' %)',/,                 &
             4X,'Ion power repartition')             
1110 format(  5X,'Ion species',A6,2X,1P,E8.2,' MW (',0P,F5.2,' %)' ) 

    return
  end subroutine broadcast2astra

!=======================================================================
!=======================================================================
  subroutine equigs4toric
!
!***********************************************************************
!                                                                      *
!  Purpose                                                             *
!  -------                                                             *
!  Create the TORIC equilibrium file from the ASTRA equilibrium.       *
!                                                                      *
!***********************************************************************

    integer, parameter :: cocos=11
    character(len=80) :: rho_str

!----- Generate raw eqdsk ----------------------------------------------
    call eqdsk(cocos, trim(exedir) // '/toric.eqdsk')

!----- Convert to TORIC inverse-equilibrium format ---------------------
    write(rho_str, '(F8.3)') rho_max
    cmd_str = 'python3 ' //  &
         '../../python/eqdsk_to_TORIC_gs.py toric.eqdsk'  //        &
         ' --rho_max ' // trim(adjustl(rho_str))
    call run_cmd(cmd_str, 'eqdsk_to_TORIC.py')

    return
  end subroutine equigs4toric

!=======================================================================
!=======================================================================
  subroutine profnt4toric
!
!***********************************************************************
!                                                                      *
!  Purpose                                                             *
!  -------                                                             *
!  Write the plasma-profile file for TORIC from ASTRA arrays.          *
!                                                                      *
!***********************************************************************
!-----------------------------------------------------------------------
! Local variables
! specs        : array of species descriptors
! nspc         : number of active ion species found
! mspc         : index of the main (most-abundant) species
! conc(isp)    : flux-surface-averaged relative concentration n_isp/ne
! chg(isp)     : charge number of species isp
! mass(isp)    : mass number of species isp
! max_conc     : running maximum of conc, used to find mspc
! mask         : logical work array
! rne          : 1/ne (zero where ne=0) for concentration computation
!-----------------------------------------------------------------------

    ! Named constants for output
    integer, parameter :: kdiff_idens = 0
    integer, parameter :: kdiff_itemp = 0
    character(10), parameter :: descript = 'torifpql'

    real(r8) :: chg
    integer  :: nspc, mspc, isp, i, lun_out
    real(r8) :: tsum, max_conc
    real(r8) :: rne(na1)
    logical  :: mask(na1) 
    real(r8) :: sum_others

!--- Build species descriptors -------------------------------------------
    ! Fixed-charge species
    specs(1) = species_t(nhydr(1:na1), 1._r8, 1._r8, 0._r8, 'HYDR ', .false., null(), 0._r8)
    specs(2) = species_t(ndeut(1:na1), 1._r8, 2._r8, 0._r8, 'DEUT ', .false., null(), 0._r8)
    specs(3) = species_t(ntrit(1:na1), 1._r8, 3._r8, 0._r8, 'TRIT ', .false., null(), 0._r8)
    specs(4) = species_t(nalf(1:na1),  2._r8, 4._r8, 0._r8, 'He4  ', .false., null(), 0._r8)
    specs(5) = species_t(nhe3(1:na1),  2._r8, 3._r8, 0._r8, 'He3  ', .false., null(), 0._r8)
    ! Variable-charge impurities
    specs(6) = species_t(niz1(1:na1),  0._r8, aim1,  0._r8, 'IMP1 ', .true.,  zim1, 0._r8)
    specs(7) = species_t(niz2(1:na1),  0._r8, aim2,  0._r8, 'IMP2 ', .true.,  zim2, 0._r8)
    specs(8) = species_t(niz3(1:na1),  0._r8, aim3,  0._r8, 'IMP3 ', .true.,  zim3, 0._r8)

!--- initialise output arrays --------------------------------------------
    nspc     = 0
    mspc     = 1
    max_conc = 0._r8

!--- 1/ne, guarded against division by zero ------------------------------
    rne = merge(1._r8 / ne(1:na1), 0._r8, ne(1:na1) > ths)

!-------------------------------------------------------------------------
! Species loop                   
!-------------------------------------------------------------------------
    do isp = 1, nspec_max
       
       mask = (specs(isp)%dens > ths)
       if (.not. any(mask)) cycle

       tsum      = sum(specs(isp)%dens* rne, mask=mask)
       specs(isp)%conc = tsum / real(count(mask), r8)
       nspc      = nspc + 1

       if (specs(isp)%has_var_chg) then
          ! Impurity: average charge from the Z array
          mask = (specs(isp)%chg_arr > ths)
          specs(isp)%chg  = sum(specs(isp)%chg_arr, mask=mask) &
               / real(count(mask), r8)
       end if

       if (specs(isp)%conc > max_conc) then
          max_conc = specs(isp)%conc
          mspc     = isp
       end if

    end do

!----- Charge neutrality -----------------------------------------------
    sum_others = sum(specs(:)%chg * specs(:)%conc) - specs(mspc)%chg * specs(mspc)%conc
    specs(mspc)%conc = (1._r8 - sum_others) / specs(mspc)%chg

!-------------------------------------------------------------------------
! Radial grid
!-------------------------------------------------------------------------
    if (abs(psibo - psiax) < epsilon(psiax)) then
       write(lun_scr, *) 'ERROR: psibo == psiax, cannot compute rho_pol'
       stop 'TORIC: profnt4toric'
    end if

    npt_ts = 0
    if(.not. allocated(rpol_as)) allocate(rpol_as(na1))
    do i = 1, na1
       rpol_as(i) = sqrt((fp(i) - psiax) / (psibo - psiax))
       if (rpol_as(i) > rho_max) exit
       npt_ts = i
    end do

    if (npt_ts == 0) then
       write(lun_scr, *) 'ERROR: no radial points within rho_max'
       stop 'TORIC: profnt4toric'
    end if

!-----------------------------------------------------------------------
! Write profile file in TORIC format
!-----------------------------------------------------------------------
    open(newunit=lun_out, file=profnt)

    write(lun_out, 101) descript, npt_ts, nspc, 1, kdiff_idens, kdiff_itemp

    ! Mass and charge header lines
    write(lun_out, 106) int(specs(mspc)%mass), int(specs(mspc)%chg)
    do isp = 1, nspec_max
       if (specs(isp)%conc > ths .and. isp /= mspc) &
            write(lun_out, 106) int(specs(isp)%mass), int(specs(isp)%chg)
    end do

    ! Normalised poloidal rho grid
    write(lun_out, 102) 'rho pol'
    write(lun_out, 111) (rpol_as(i) / rpol_as(npt_ts), i = 1, npt_ts)

    ! Electron density [cm^-3] 
    write(lun_out, 102) 'ne [cm-3]'
    write(lun_out, 111) (ne(i) * ne_scale, i = 1, npt_ts)

    ! Electron temperature [keV]
    write(lun_out, 102) 'te [keV]'
    write(lun_out, 111) (te(i), i = 1, npt_ts)

    ! Main species concentration (scalar)
    write(lun_out, 102) specs(mspc)%name
    write(lun_out, 107) specs(mspc)%conc

    ! Ion temperature [keV]
    write(lun_out, 102) 'ti [keV]'
    write(lun_out, 111) (ti(i), i = 1, npt_ts)

    ! Minority / impurity concentrations
    do isp = 1, nspec_max
       if (specs(isp)%conc > ths .and. isp /= mspc) then
          write(lun_out, 102) specs(isp)%name
          write(lun_out, 107) specs(isp)%conc
       end if
    end do

101 format(A10, 5I4)
102 format(A10)
106 format(2I4)
107 format(1P, E16.9, 0P)
111 format(1P, 5E16.9, 0P)

    close(lun_out)

    return
  end subroutine profnt4toric

!=======================================================================
!=======================================================================

!***********************************************************************
! Portable chdir wrapper
!***********************************************************************
  subroutine safe_chdir(path, ios)
    character(len=*), intent(in)  :: path
    integer,          intent(out) :: ios
    integer(c_int) :: cret
    cret = c_chdir(trim(path) // c_null_char)
    ios  = int(cret)
    if (ios /= 0) then
       write(lun_scr, '(2A)') 'ERROR: cannot chdir to ', trim(path)
       stop 'TORIC: safe_chdir'
    end if
    return
  end subroutine safe_chdir

!***********************************************************************
! Portable execute_command_line wrapper with consistent error checking
! always checks both cmdstat and exitstat
!***********************************************************************
  subroutine run_cmd(cmd, label)
    character(len=*), intent(in) :: cmd, label
    integer :: exit_ios, cmd_ios
!    call execute_command_line('sleep 1')
    call execute_command_line(trim(cmd), exitstat=exit_ios, cmdstat=cmd_ios)
    if (cmd_ios /= 0 .or. exit_ios /= 0) then
       write(lun_scr, '(3A,2(A,I0))') &
            'ERROR: command failed [', trim(label), ']: ', &
            'cmdstat=', cmd_ios, '  exitstat=', exit_ios
       write(lun_scr, '(2A)') '  cmd: ', trim(cmd)
       stop 'TORIC: run_cmd'
    end if
    return
  end subroutine run_cmd
  
!***********************************************************************
! Remove the local TORIC files
!***********************************************************************
  subroutine clean_up   

    logical :: is_open

!----- Cose all possible files opened by TORIC-SSFPQL ------------------
    inquire(unit=lun5, opened=is_open)
    if (is_open) close(lun5)
    
    inquire(unit=lun6, opened=is_open)
    if (is_open) close(lun6)

    inquire(unit=lun8, opened=is_open)
    if (is_open) close(lun8)

    inquire(unit=lun9, opened=is_open)
    if (is_open) close(lun9)
    
    inquire(unit=lun17, opened=is_open)
    if (is_open) close(lun17)

    inquire(unit=lun18, opened=is_open)
    if (is_open) close(lun18)

    inquire(unit=lun21, opened=is_open)
    if (is_open) close(lun21)

    inquire(unit=lun22, opened=is_open)
    if (is_open) close(lun22)

    inquire(unit=lun6, opened=is_open)
    if (is_open) close(lun6)
    
!----- Remove the directory    
    cmd_str = '/bin/rm -rvf ' // trim(exedir)
    call run_cmd(cmd_str, 'cleanup exedir')
    return
  end subroutine clean_up
!
!=======================================================================
!=======================================================================
end module torfpql_mod
