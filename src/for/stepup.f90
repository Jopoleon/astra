subroutine STEPUP
!-------------------------------------------------------------------
! Perform one time step
! Note that now time step is updated at the end of a full time cycle
!-------------------------------------------------------------------

use const_inc, only: IPART, ITFBE, IFBEY, IPLFBE, IFBEG, &
    IPCTRL, ICIRCQ, ITFBP, ITREQ, FTN, FTO, BTN, BTOR, HRO, ROC, NA1, &
    TAU, TAUMIN, TAUMAX, TAUPRP, TIME, TSTART, ATREQ, LEQ, & 
    PSIFBO, PSIFB, PSIEXO, PSIEXT, PSPLXO, PSPLEX, RBDOT, BBDOT
use status_inc, only: TE, TI, NE, NI, NIO, FP, defarr
use io_mod, only: CCOIL, VCOIL, MACHINE, TASK
use read_input, only: raw_cCoil, raw_vCoil
use plasma_state, only: plasma_up

implicit none

integer :: IFKEY, IFSUB, ibcpsi_fb, bc_type_for_fp, jkey, &
    IFTREQ, IFSTEP, n_coils
double precision :: zipctrl, iplfbeo, Apsibcfac, Bpsibcfac, dfpdrbm12, &
    tau_old, tau_new
double precision, dimension(raw_cCoil%ncoils) :: yccoil
double precision, dimension(raw_vCoil%ncoils) :: yvcoil

data ibcpsi_fb /0/

save ibcpsi_fb  ! counter to use psi as bc stuff
save tau_old, tau_new

!Initialize a few variables for toroidal field
BTN = BTOR
FTN = FTO

n_coils = raw_cCoil%ncoils

if (TIME <= TSTART) then
    tau_old = tau
    tau_new = tau
endif
 
! tau treatment to avoid machine precision errors
tau = tau_old ! still using old one up to equations
dfpdrbm12 = 0.
IPART = 2     ! Mark time evolution section

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

call detvar
! Subroutines with the "<" symbol are put here
! here it computes the new NI also. These are run with tau_old

if (plasma_up == 1 .or. ifbey == 0) then
    call OLDNEW           ! Time advance: F(t-tau):=F(t) neo=ne, etc,except ni
    call INTVAR           ! Set exp scalars, moved here for btor consistency
endif

! Update time at the end of everything

iplfbeo = iplfbe

!reset some quantities
PSIFBO = PSIFB       ! reset boundary flux
RBDOT  = 0.          ! reset boundary adiabatic factor
BBDOT  = 0.          ! reset boundary adiabatic factor
PSIEXO = PSIEXT      ! reset also external flux from fbe and ce, this is for test!
PSPLXO = PSPLEX      ! reset also green function flux from fbe and ce, this is for test!

if (TIME > ITFBE) then ! if free boundary, solve circuit equations, ccoil comes from there
    yccoil = CCOIL(1: raw_cCoil%ncoils)
    yvcoil = VCOIL(1: raw_vCoil%ncoils)
else ! if time <= ITFBE, ccoil and vcoil comes from experimental traces in exp file
    call get_coil(TIME, raw_cCoil, yccoil)
    call get_coil(TIME, raw_vCoil, yvcoil)
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

time_step_accuracy: do

    tau = tau_new !use new tau now

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

! here it should go the correction after 1st free boundary call since geometry changes abruptly
            if (IFBEY == 1.) then
                RBDOT = 0.   !also set compression to zero to avoid jumps
                BBDOT = 0.   !also set compression to zero to avoid jumps
            endif
        endif

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

        if (TASK(1:3) /= 'BGD') then
            jkey = IFKEY(0)                 ! Enables ITREQ iteration control 
        endif
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

            if (ITFBP < 0.0) then
                if (IFBEY >= 1.) then
                    if (ibcpsi_fb == 1) then
                        FP = FP - FP(NA1) + PSIFB
                        bc_type_for_fp = 3
                    endif
                endif
            endif
        endif

    enddo tr_eq_loop

    if (plasma_up == 1 .or. IFBEY == 0) then
        call DEFARR                  ! F(t)>0? Define F(t) outside ABC
    endif
    
    tau_old = tau !store old tau before changing it

    if (IFSTEP() == 0 .and. IFBEY /= 1) then
! Time step accuracy accepted? No(0)
! note that here TAU is modified and TIME updated with time_new = TIME+TAU !
        tau_new = tau
    else
        EXIT
    endif                        

enddo time_step_accuracy

tau = max(taumin, tau) !just reforce it to be between limits
tau = min(taumax, tau) !just reforce it to be between limits

tau_new = tau !store new tau
tau = tau_old !reuse old for postep routines

! When circuit equations are used do this
if (IFBEY >= 1.) then         ! is doing free boundary
    if (ICIRCQ > 0.) then    ! circuit equations are solved with whatever code
        if (LEQ(5) == 4) then ! SPIDER
            call SPIDUPDATE(machine, CCOIL(1:n_coils), time, n_coils)    ! Update circuit stuff which has to be outside the iterations of course
        else if (LEQ(5) == 5) then ! FEQIS
            call FEQISUPDATE(CCOIL(1:n_coils), n_coils)    ! Update circuit stuff which has to be
        endif
    endif
endif

TIME = TIME + TAU
TAUPRP = tau
call POSTEP

if (TAU /= TAUPRP) then
    tau_new = tau !store new tau in case it has been changed in postep
    tau = tau_old !in case tau has been changed
    tauprp = tau
endif

! Special Sbrs call.  Subroutines with the ">" symbol in their 
! call are put here, so outside the convergence or time step check loop
! note that in postep if one wants to modify tau, like in tsctrl, better to do it in tauprp
! call TSCTRL at the end of all other subroutines

return
end subroutine STEPUP
