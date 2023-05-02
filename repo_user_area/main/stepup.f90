subroutine STEPUP
!-------------------------------------------------------------------
! Perform one time step
! Note that now time step is updated at the end of a full time cycle
!-------------------------------------------------------------------

use parameter_inc, only: NRD
use const_inc, only: IPSMK, IPART, ITFBE, IFBEY, IPLFBE, IFBEG, &
    IPCTRL, NCNB, ICIRCQ, ITFBP, ITREQ, UPDWN, FTN, FTO, BTN, BTOR, HRO, ROC, NA1, &
    TAU, TAUMIN, TAUMAX, TAUPRP, TIME, TSTART, ATREQ, LEQ, & 
    PSIFBO, PSIFB, PSIEXO, PSIEXT, PSPLXO, PSPLEX, ADCMPF, RBDOT, BBDOT
use status_inc, only: TE, TI, NE, NI, NIO, FP
use outcmn_inc, only: CCOIL, CCOILO, DUMCT, DUMCTP, CTRLM, VCOIL, MACHINE
use plasma_state, only: plasma_up
use debugger, only: markloc, flightsim
use fs_coupling_variables, only: fs_dt_smlk,fs_dt_tctrl

implicit none

integer :: IFKEY, IFSUB, &
    ibcpsi_fb, jreadd, icurradj, bc_type_for_fp, jkey, &
    IFTREQ, IFSTEP

double precision :: tau_new, updwno, zipctrl, &
    iplfbeo, Apsibcfac, Bpsibcfac, dfpdrbm12, time_ext, &
    tau_temp_smlk, tau_old, dt_smlk

double precision, dimension(NRD) :: dummycoils

data ibcpsi_fb /0/
data jreadd /0/
data tau_temp_smlk /0./

save ibcpsi_fb, jreadd  ! counter to use psi as bc stuff
save dt_smlk
save time_ext
save tau_temp_smlk

!Initialize a few variables for toroidal field
BTN = BTOR
FTN = FTO

if (flightsim < 1) plasma_up = 1

!tau treatment to avoid machine precision errors
if (flightsim >= 1) then
    tau = 1.d-6*nint(tau*1.d6)
    tau_temp_smlk = 1.d-6*nint(tau_temp_smlk*1.d6)
endif

dfpdrbm12 = 0.

tau_old = tau
tau_new = tau

! wait until constants file is read and read control file
if (flightsim == 1) then
    if (jreadd == 0) call read_input_constant_file(time_ext, dt_smlk)
    if (TIME-TSTART == 0) then
        tau     = taumin
        tau_old = tau
        tau_new = tau
        tau_temp_smlk = tau
    endif
    if (tau_temp_smlk == 0.) tau_temp_smlk = tau
else if (flightsim >= 1) then
    if (jreadd == 0) then
        time_ext = TIME
        dt_smlk = fs_dt_smlk
    endif
    if (tau_temp_smlk == 0.) tau_temp_smlk = tau
endif

IPART = 2             ! Mark time evolution section

! reset initial condition
				
NIO = NI  !moved from OLDNEW here to maintain it correctly. detvar goes before oldnew to maintain time derivative computations. 
						
!switch from pbe to fbe gs solver
if (ITFBE < 0.) IFBEY = 0.
if (ITFBE > 0.) then
    if (TIME < ITFBE) then
        IFBEY = 0.
    else
        IFBEY = IFBEY + 1.
    endif
endif

!include 'tmp/detvar.tmp' ! GIT 
call detvar
! Subroutines with the "<" symbol are put here
! here it computes the new NI also

if (plasma_up == 1 .or. ifbey == 0) then
    call OLDNEW           ! Time advance: F(t-tau):=F(t) neo=ne, etc,except ni
    call INTVAR           ! Set exp scalars, moved here for btor consistency
endif

! Update time at the end of everything

iplfbeo = iplfbe

!reset some quantities
CCOILO = CCOIL
PSIFBO = PSIFB          ! reset boundary flux
RBDOT = 0.              ! reset boundary adiabatic factor
BBDOT = 0.              ! reset boundary adiabatic factor
PSIEXO = PSIEXT	      ! reset also external flux from fbe and ce, this is for test!
PSPLXO = PSPLEX  	      ! reset also green function flux from fbe and ce, this is for test!
						
!Get target quantities for control
if (plasma_up == 1 .or. ifbey == 0) then
    if (nint(IPCTRL) /= 0 .and. nint(IPCTRL) > -2) then
        call GETDUMCT(DUMCT)   ! in this case VCOIL and CCOIL are not taken from experimental file, rather from controller
        call GETCTRLMS(CTRLM)   
    else
! Get coils currents from experimental file if no control
	
        call GETDUMCT(DUMCT)   
        call GETCTRLMS(CTRLM)   
        DUMCTP = DUMCT

! do this only if flightsim = 0, so that with -1 it doesnt do this.
! if (flightsim==0) then
        call GETCOILS(VCOIL(1:NCNB), dummycoils(1:NCNB))

        if (ICIRCQ == 0.) then
            call GETCOILS(VCOIL(1:NCNB), CCOIL(1:NCNB))
        endif
    endif
endif

! counter for psi bc = -1 
if (ITFBP == 0.0) ibcpsi_fb = 0

if (IFBEY >= 1.) then
    if (ITFBP < 0.0) then
        ibcpsi_fb = ibcpsi_fb + 1
        ibcpsi_fb = min(ibcpsi_fb, 3)
    else
        ibcpsi_fb = 0
    endif
endif

if (flightsim >= 1) then
    tau = tau_temp_smlk
endif

updwno = updwn          ! for fsim

time_step_accuracy: do

    ITREQ = 0               ! Start Tr-Eq loop

!reset some quantities
    IPLFBE = iplfbeo        ! reset plasma current for fbe
    PSIFB  = PSIFBO         ! reset boundary flux
    RBDOT  = 0.             ! reset boundary adiabatic factor
    BBDOT  = 0.             ! reset boundary adiabatic factor
    PSIEXT = PSIEXO         ! reset also external flux from fbe and ce, this is for test!
    PSPLEX = PSPLXO         ! reset also green function flux from fbe and ce, this is for test!

! start Tr-EQ loop
    jkey = 0
    tr_eq_loop: do while (jkey == 0)

! NITREQ regulates this. 
! 1 - no iterations, 2 - yes. is 1 by default

        if (plasma_up == 0 .or. ifbey == 0) then
            call SETARX(2)       ! Update exp-data with a new metric
        endif

        call METRIC          ! Equilibrium call, compute IPL from dfpdrb, compute PSIEXT, shape, psplex, and metric coefficients, update ROC, FTN

	if (plasma_up == 1) then
            RBDOT = (FTO  - FTN)/(FTO  + FTN)/TAU     !New rbdot for adiabatic compression
            BBDOT = (BTOR - BTN)/(BTOR + BTN)/TAU     !New bbdot for adiabatic compression

            if (nint(ADCMPF) == 2) then
                RBDOT = 0.   !no adiabatic compression whatsoever
                BBDOT = 0.   !no adiabatic compression whatsoever
            endif

! here it should go the correction after 1st free boundary call since geometry changes abruptly
            if (IFBEY == 1.) then
                RBDOT = 0.   !also set compression to zero to avoid jumps
                BBDOT = 0.   !also set compression to zero to avoid jumps
            endif
        endif

        tau_old = tau !remember old tau for later
        if (isnan(hro)) then
            write(*, *) 'hro is nan'
            call err_catch_a
        endif

        if (plasma_up == 1) then
            call eqns_inc(ibcpsi_fb, bc_type_for_fp, dfpdrbm12)
        endif

! catching errors: infinite or nan profiles
        if (sum(abs(te(1:na1)))/na1 > 1e8) then
            write(*, *) 'te isinf'
            call err_catch_a
        endif
        if (sum(abs(ti(1:na1)))/na1 > 1e8) then
            write(*, *) 'ti isinf'
            call err_catch_a
        endif
        if (sum(abs(ne(1:na1)))/na1 > 1e8) then
            write(*, *) 'ne isinf'
            call err_catch_a
        endif
        if (sum(abs(fp(1:na1)))/na1 > 1e8)  then
            write(*, *) 'fp isinf'
            call err_catch_a
        endif
        if (isnan(sum(te(1:na1))))  then
            write(*, *) 'te isnan'
            call err_catch_a
        endif
        if (isnan(sum(ti(1:na1))))  then
            write(*, *) 'ti isnan'
            call err_catch_a
        endif
        if (isnan(sum(ne(1:na1))))  then
            write(*, *) 'ne isnan'
            call err_catch_a
        endif
        if (isnan(sum(fp(1:na1))))  then
            write(*, *) 'fp isnan'
            call err_catch_a
        endif

        jkey = IFKEY(0)                 ! Enables ITREQ iteration control 
        jkey = IFTREQ(ATREQ)            ! ++ITREQ; Tr-Eq loop converged?  

! some options to avoid NITREQ when IFBEY = 1, IPCTR = X.1  --> does not do NITREQ
        if (IFBEY == 1.) then
            zipctrl = IPCTRL - nint(IPCTRL)
            if (zipctrl > 1.e-10) jkey = 1
        endif

!quantitites for psi b.c.
 	if (plasma_up == 1) then
            Apsibcfac = dfpdrbm12
            Bpsibcfac = PSIEXT - PSPLEX*ROC*Apsibcfac
            PSIFB = Bpsibcfac

            icurradj = 0
            if (ITFBP < 0.0) then
                if (IFBEY >= 1.) then
                    if (ibcpsi_fb == 1) then
                        FP = FP - FP(NA1) + PSIFB
                        bc_type_for_fp = 3
                        icurradj = 1
                    endif
                endif
            endif
        endif

    enddo tr_eq_loop

    if (plasma_up == 1 .or. IFBEY == 0) then
        call DEFARR                  ! F(t)>0? Define F(t) outside ABC
    endif
    tau_old = tau

    if (IFSTEP(jkey, updwno) == 0 .and. IFBEY /= 1) then
! Time step accuracy accepted? No(0)
! note that here TAU is modified and TIME updated with time_new = TIME+TAU !
        if (flightsim >= 1) then
            TAU = max(taumin, TAU_old - fs_dt_tctrl) ! correct TAU not to exceed time_ext
            TAU = max(taumin, TAU)
        endif
    else
        if (flightsim >= 1) then
            TAU = min(taumax, TAU_old + fs_dt_tctrl)
        endif
        tau_new = tau
        tau = tau_old
        EXIT
    endif                        

enddo time_step_accuracy

! When circuit equations are used do this
if (IFBEY >= 1.) then         ! is doing free boundary
    if (ICIRCQ > 0.) then    ! circuit equations are solved with whatever code
        if (LEQ(5) == 4) then ! SPIDER
            if (nint(IFBEG) == 0) then
                call SPIDUPDATE(machine, CCOIL(1:NCNB), time, ncnb)    ! Update circuit stuff which has to be outside the iterations of course
            else if (nint(IFBEG) == 1) then
                call f_SPIDUPDATE(machine, CCOIL(1:NCNB), time, ncnb)  ! Update circuit stuff which has to be outside the iterations of course
            endif
        else if (LEQ(5) == 5) then ! FEQIS
!            call FEQISUPDATE(machine, CCOIL(1:NCNB), time, ncnb)    ! Update circuit stuff which has to be
        endif
    endif
endif

!adiabatic compression, not more useful...
if (nint(ADCMPF) == 0 .and. icurradj == 0) then
    call ADCMP(bc_type_for_fp, dfpdrbm12)  ! Do AdComp once per time step, since if 
endif

TIME = TIME + TAU
TAUPRP = tau_new
call POSTEP

! Special Sbrs call (METRIC?)  Subroutines with the ">" symbol in their 
! call are put here, so outside the convergence or time step check loop
! note that in postep if one wants to modify tau, like in tsctrl, better to do it in tauprp
                                                                      
tau_new = tauprp

if (flightsim >= 1) then
    tau_temp_smlk = TAU_new
    taumin = min(tau_temp_smlk, taumin)
    if (TIME-TSTART >= time_ext+dt_smlk-1.e-8) then
        call write_output_diag_file
        jreadd = 0
    else
        jreadd = 1
        if (TIME-TSTART > time_ext+dt_smlk+1.d-6) then
            write(*, *) 'problems with sync, stopping', tau, tau_temp_smlk, &
                time, time_ext + dt_smlk + 1.d-6
            call err_catch_a
        endif
    endif
endif

tau = min(taumax, tau_new)
if (flightsim >= 1) then
    tau     = min(taumax, tau_temp_smlk)
    tau_new = tau
endif
tau = max(taumin, tau)

return
end subroutine STEPUP
