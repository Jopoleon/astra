module metrics

use pi_const, only: GP, GP2, GP2_sq

implicit none

logical :: use_ext_bnd=.false.
logical :: plasma_up=.true.  ! plasma is up by default, can be set to False for breakdown by the user in a user-defined sbr called with "<"
integer :: naxis
real*8, allocatable :: raxiscc(:), zaxiscc(:)
double precision, dimension(:), allocatable :: CCOIL, VCOIL


contains

!---------------------------------------------------------------------
    subroutine METRIC()

    use cpu_usage, only: wallTime_equ, cpuTime_equ
    use status, only: VRO, VR, SHIF, AMETR, ELON, TRIA, XRHO, FP, IPOL
    use scalars, only: IPART, FTO, FTN, ROC, &
        BTOR, ROCO, RTOR, SHIFT, &
        ABC, ELONG, TRIAN, UPDWN, NA1, NB1, MEQUIL, NEQUIL, &
        IPEQL, vmec_option, TIME, TSTART, TIMEQL, DTEQL, BTN
    use debugger, only: markloc
    use parameters_a2equil, only: equil_now
    use numerical_tools, only: qinterp
    use scalars, only: tau

    integer :: i, jexit, NDTEQUILMY, equil_solver, jthe, nrho_surf, nthe_surf
    integer :: t_wall1, t_wall2, rate
    double precision :: theta
    integer :: vmec_vacuum, vmec_tau, vmec_dteq, yes_boozer
    integer :: i_vmec_options_choose
    real :: t_cpu1, t_cpu2
    double precision, allocatable, dimension(:) :: prof_as, prof_eq
    character(len=120) :: err_msg

    data i_vmec_options_choose/0/
    save vmec_vacuum, vmec_tau, vmec_dteq, yes_boozer, i_vmec_options_choose

    call markloc('METRIC')

    if (IPART == 1 .or. IPART == 3) then ! do only at initiation
        FTN = FTO
        BTN = BTOR
        ROC = sqrt(FTO/GP/BTOR)
        ROCO = ROC
        VRO(1: NB1) = VR(1: NB1)
! if IPEQL = 70, 71, 72, ..., 79
! vmec_option between 0 and 9, set IPEQL to 7
    endif

    if (i_vmec_options_choose == 0) then
        SELECT CASE(vmec_option)
        CASE(0)
            vmec_vacuum = 1
            vmec_tau    = 1
            vmec_dteq   = 1
            yes_boozer  = 1
        CASE(1)
            vmec_vacuum = 2
            vmec_tau    = 1
            vmec_dteq   = 1
            yes_boozer  = 1
        CASE(2)
            vmec_vacuum = 2
            vmec_tau    = 2
            vmec_dteq   = 1
            yes_boozer  = 1
        CASE(3)
            vmec_vacuum = 2
            vmec_tau    = 2
            vmec_dteq   = 0
            yes_boozer  = 1
        CASE(4)
            vmec_vacuum = 1
            vmec_tau    = 1
            vmec_dteq   = 1
            yes_boozer  = 0
        CASE(5)
            vmec_vacuum = 2
            vmec_tau    = 1
            vmec_dteq   = 1
            yes_boozer  = 0
        CASE(6)
            vmec_vacuum = 2
            vmec_tau    = 2
            vmec_dteq   = 1
            yes_boozer  = 0
        CASE(7)
            vmec_vacuum = 2
            vmec_tau    = 2
            vmec_dteq   = 0
            yes_boozer  = 0
        CASE(8)
            vmec_vacuum = 2
            vmec_tau    = 2
            vmec_dteq   = 0
            yes_boozer  = 1
        CASE(9)
            vmec_vacuum = 2
            vmec_tau    = 2
            vmec_dteq   = 0
            yes_boozer  = 0
        END SELECT
        i_vmec_options_choose = 1
    endif

    if (IPART == 3) then
        IPART = 1
        return
    endif

    call CPU_TIME(t_cpu1)
    call SYSTEM_CLOCK(t_wall1, rate)

    SELECT CASE(IPEQL)

    CASE(-2)  ! Cylindircal case, No equilibrium solver. No toroidicity
        call EQCYL()
        call RHSEQ()

    CASE(-1)  ! Take metric from exp/data_file
        ROC = ROC3A(RTOR, SHIFT, ABC, ELONG, TRIAN)
        FTO = GP*BTOR*ROC**2
        call set_external_metric()! Main grid: (jj-0.5)*h
        call RHSEQ()! this computes ffprime and pprime

    CASE(-6)  ! Take metric internally subroutine external file (so not X)
        ROC = ROC3A(RTOR, SHIFT, ABC, ELONG, TRIAN)
        FTO = GP*BTOR*ROC**2
        call set_external_metric_2()! Main grid: (jj-0.5)*h
        call RHSEQ()! this computes ffprime and pprime

    CASE(6)  ! Take stellarator metric internally subroutine external file (so not X)
       ROC = ROC3A(RTOR, SHIFT, ABC, ELONG, TRIAN)
       FTO = GP*BTOR*ROC**2
       call set_external_metric_2() ! Main grid: (jj-0.5)*h
       call RHSEQ()  ! this computes ffprime and pprime

    CASE(0) ! No equilibrium solver (NEQUIL=0) .or. data initiation @ 1st entry
        call EQGUESS()
    CASE(1)  ! EMEQ
        if (TIME == TSTART) NDTEQUILMY = 0
        if (TIME >  TSTART) NDTEQUILMY = 1
        if (TIME-TIMEQL >= NDTEQUILMY*DTEQL) then
            call RHSEQ()! Define p', FF', j_tor=CUTOR
            call A2EMEQ(jexit)
            if (jexit /= 0) then
                err_msg = 'Equilibrium problem at the initial iterations'
                if (IPART == 1) then
                    write(*, '(A)') err_msg
                    ERROR STOP
                endif
            endif
            TIMEQL = TIME
        endif

    CASE(3)  ! equil iterations
        if (TIME == TSTART) NDTEQUILMY = 0
        if (TIME >  TSTART) NDTEQUILMY = 1
        if (TIME-TIMEQL >= NDTEQUILMY*DTEQL) then
            call RHSEQ()
            TIMEQL = TIME
        endif

    CASE(4: 5)  ! SPIDER, FEQIS
        if (IPEQL == 4) then
            equil_solver = 3
        else
            equil_solver = 101
        endif
        if (TIME == TSTART) NDTEQUILMY = 0
        if (TIME >  TSTART) NDTEQUILMY = 1
        if (TIME-TIMEQL >= NDTEQUILMY*DTEQL) then
            call RHSEQ2()! Define p', FF', j_tor=CUTOR, but using the gssolver definitions
            call A2GSSOLVER(equil_solver)
            TIMEQL = TIME
        endif

    CASE(7)   !stellarator equilibrium solver
        nrho_surf = abs(NEQUIL)
        nthe_surf = abs(MEQUIL)
        if (associated(equil_now%coord_sys%position%r)) then
            deallocate(equil_now%coord_sys%position%r)
            deallocate(equil_now%coord_sys%position%z)
            deallocate(equil_now%coord_sys%position%rmin)
            deallocate(equil_now%coord_sys%position%psirz)
            deallocate(equil_now%coord_sys%position%theta2d)
            deallocate(equil_now%eqgeometry%rectgrid%psirz2d)
            deallocate(equil_now%eqgeometry%rectgrid%fdia2d)
            deallocate(equil_now%eqgeometry%rectgrid%r2d)
            deallocate(equil_now%eqgeometry%rectgrid%z2d)
        endif
        allocate(equil_now%coord_sys%position%r(nrho_surf, nthe_surf))
        allocate(equil_now%coord_sys%position%z(nrho_surf, nthe_surf))
        allocate(equil_now%coord_sys%position%rmin(nrho_surf, nthe_surf))
        allocate(equil_now%coord_sys%position%psirz(nrho_surf, nthe_surf))
        allocate(equil_now%coord_sys%position%theta2d(nthe_surf))
        allocate(equil_now%eqgeometry%rectgrid%psirz2d(1, 1))
        allocate(equil_now%eqgeometry%rectgrid%fdia2d(1, 1))
        allocate(equil_now%eqgeometry%rectgrid%r2d(1))
        allocate(equil_now%eqgeometry%rectgrid%z2d(1))
        allocate(equil_now%coord_sys%gradvcell(nrho_surf, nthe_surf))
        allocate(equil_now%coord_sys%bpcell(nrho_surf, nthe_surf))
        allocate(equil_now%coord_sys%bcell(nrho_surf, nthe_surf))
        allocate(equil_now%coord_sys%rcell(nrho_surf, nthe_surf))
        allocate(equil_now%coord_sys%darea(nrho_surf, nthe_surf))
        allocate(equil_now%coord_sys%jphi(nrho_surf, nthe_surf))

        if (.not. associated(equil_now%profiles_1d%rho_tor_norm)) then
            allocate(equil_now%profiles_1d%areat  (nrho_surf))
            allocate(equil_now%profiles_1d%bdb0   (nrho_surf))
            allocate(equil_now%profiles_1d%bmaxt  (nrho_surf))
            allocate(equil_now%profiles_1d%bmint  (nrho_surf))
            allocate(equil_now%profiles_1d%dpsidv (nrho_surf))
            allocate(equil_now%profiles_1d%F_dia  (nrho_surf))
            allocate(equil_now%profiles_1d%ffprime(nrho_surf))
            allocate(equil_now%profiles_1d%fofb   (nrho_surf))
            allocate(equil_now%profiles_1d%g1     (nrho_surf))
            allocate(equil_now%profiles_1d%g2     (nrho_surf))
            allocate(equil_now%profiles_1d%ggradro(nrho_surf))
            allocate(equil_now%profiles_1d%gm1    (nrho_surf))
            allocate(equil_now%profiles_1d%gm4    (nrho_surf))
            allocate(equil_now%profiles_1d%gm41   (nrho_surf))
            allocate(equil_now%profiles_1d%gm5    (nrho_surf))
            allocate(equil_now%profiles_1d%perim  (nrho_surf))
            allocate(equil_now%profiles_1d%phi    (nrho_surf))
            allocate(equil_now%profiles_1d%pprime (nrho_surf))
            allocate(equil_now%profiles_1d%pressure(nrho_surf))
            allocate(equil_now%profiles_1d%psi    (nrho_surf))
            allocate(equil_now%profiles_1d%q      (nrho_surf))
            allocate(equil_now%profiles_1d%rho_tor(nrho_surf))
            allocate(equil_now%profiles_1d%rho_tor_norm(nrho_surf))
            allocate(equil_now%profiles_1d%shif   (nrho_surf))
            allocate(equil_now%profiles_1d%surface(nrho_surf))
            allocate(equil_now%profiles_1d%volume (nrho_surf))
            allocate(equil_now%profiles_1d%elongation(nrho_surf))
            allocate(equil_now%profiles_1d%r_inboard (nrho_surf))
            allocate(equil_now%profiles_1d%r_outboard(nrho_surf))
        endif

        if (TIME == TSTART) NDTEQUILMY = 0
        if (TIME >  TSTART) NDTEQUILMY = 1

        if (TIME == TSTART .and. vmec_vacuum > 0) then
            call EQCYL_STELLA
            call RHSEQ
            call VMEC2ASTRA(vmec_vacuum, vmec_tau, vmec_dteq, yes_boozer)
            vmec_vacuum = 0
            return
        endif
        if (TIME-TSTART < 2*TAU .and. vmec_tau == 1) then ! call again vmec and vmec2astra.py
            call VMEC2ASTRA(vmec_vacuum, vmec_tau, vmec_dteq, yes_boozer)
            TIMEQL = TIME
            return
        endif
        if (TIME-TSTART < 2*TAU .and. vmec_tau == 2) then ! call only vmec2astra.py but no vmec
            call VMEC2ASTRA(vmec_vacuum, 0, 0, yes_boozer)
            TIMEQL = TIME
            vmec_tau = 0
            return
        endif
        if (TIME > TSTART .and. TIME-TIMEQL >= NDTEQUILMY*DTEQL) then !call vmec+vmec2astra.py or only vmec2astra.py
            call VMEC2ASTRA(vmec_vacuum, vmec_tau, vmec_dteq, yes_boozer)
            TIMEQL = TIME
            return
        endif
    END SELECT

    call CPU_TIME(t_cpu2)
    call SYSTEM_CLOCK(t_wall2, rate)
    cpuTime_equ = cpuTime_equ + dble(t_cpu2) - dble(t_cpu1)
    wallTime_equ = wallTime_equ + t_wall2 - t_wall1

    if (IPEQL < 3) then
        nrho_surf = abs(NEQUIL)
        nthe_surf = abs(MEQUIL)
        if (nrho_surf == 0) nrho_surf = NA1 + 1
        if (nthe_surf == 0) nthe_surf = 41
        if (.not. associated(equil_now%coord_sys%position%r)) then
            allocate(equil_now%coord_sys%position%r(nrho_surf, nthe_surf))
            allocate(equil_now%coord_sys%position%z(nrho_surf, nthe_surf))
            allocate(equil_now%coord_sys%position%rmin(nrho_surf, nthe_surf))
            allocate(equil_now%coord_sys%position%psirz(nrho_surf, nthe_surf))
            allocate(equil_now%coord_sys%position%theta2d(nthe_surf))
        endif
        if (.not. associated(equil_now%profiles_1d%rho_tor_norm)) then
            allocate(equil_now%profiles_1d%areat  (nrho_surf))
            allocate(equil_now%profiles_1d%bdb0   (nrho_surf))
            allocate(equil_now%profiles_1d%bmaxt  (nrho_surf))
            allocate(equil_now%profiles_1d%bmint  (nrho_surf))
            allocate(equil_now%profiles_1d%dpsidv (nrho_surf))
            allocate(equil_now%profiles_1d%F_dia  (nrho_surf))
            allocate(equil_now%profiles_1d%ffprime(nrho_surf))
            allocate(equil_now%profiles_1d%fofb   (nrho_surf))
            allocate(equil_now%profiles_1d%g1     (nrho_surf))
            allocate(equil_now%profiles_1d%g2     (nrho_surf))
            allocate(equil_now%profiles_1d%ggradro(nrho_surf))
            allocate(equil_now%profiles_1d%gm1    (nrho_surf))
            allocate(equil_now%profiles_1d%gm4    (nrho_surf))
            allocate(equil_now%profiles_1d%gm41   (nrho_surf))
            allocate(equil_now%profiles_1d%gm5    (nrho_surf))
            allocate(equil_now%profiles_1d%perim  (nrho_surf))
            allocate(equil_now%profiles_1d%phi    (nrho_surf))
            allocate(equil_now%profiles_1d%pprime (nrho_surf))
            allocate(equil_now%profiles_1d%pressure(nrho_surf))
            allocate(equil_now%profiles_1d%psi    (nrho_surf))
            allocate(equil_now%profiles_1d%q      (nrho_surf))
            allocate(equil_now%profiles_1d%rho_tor(nrho_surf))
            allocate(equil_now%profiles_1d%rho_tor_norm(nrho_surf))
            allocate(equil_now%profiles_1d%shif   (nrho_surf))
            allocate(equil_now%profiles_1d%surface(nrho_surf))
            allocate(equil_now%profiles_1d%volume (nrho_surf))
            allocate(equil_now%profiles_1d%elongation(nrho_surf))
            allocate(equil_now%profiles_1d%r_inboard (nrho_surf))
            allocate(equil_now%profiles_1d%r_outboard(nrho_surf))
        endif
        if (.not. allocated(prof_as)) then
            allocate(prof_as(NA1))
            allocate(prof_eq(nrho_surf))
        endif
! Fill array values
        equil_now%profiles_1d%rho_tor_norm =  (/ ((i - 1.d0)/(nrho_surf - 1.d0), i=1, nrho_surf) /)
        equil_now%profiles_1d%rho_tor = equil_now%profiles_1d%rho_tor_norm*ROC
        prof_as = FP(1: NA1)
        call qinterp(XRHO(1:NA1), prof_as, NA1, equil_now%profiles_1d%rho_tor_norm, prof_eq, nrho_surf)
        equil_now%profiles_1d%psi = prof_eq
        prof_as = IPOL(1: NA1)*RTOR*BTOR
        call qinterp(XRHO(1:NA1), prof_as, NA1, equil_now%profiles_1d%rho_tor_norm, prof_eq, nrho_surf)
        equil_now%profiles_1d%F_dia = prof_eq
        do jthe=1, nthe_surf
            theta = GP2*DBLE(jthe - 1)/DBLE(nthe_surf - 1)
            prof_as = RTOR + SHIF(1:NA1) + AMETR(1:NA1) *  ( COS(theta) + 0.5*TRIA(1:NA1) * (COS(2.*theta) - 1.))
            call qinterp(XRHO(1:NA1), prof_as, NA1, equil_now%profiles_1d%rho_tor_norm, prof_eq, nrho_surf)
            equil_now%coord_sys%position%r(:, jthe) = prof_eq
            prof_as = UPDWN + AMETR(1:NA1)*ELON(1:NA1)*SIN(theta)
            call qinterp(XRHO(1:NA1), prof_as, NA1, equil_now%profiles_1d%rho_tor_norm, prof_eq, nrho_surf)
            equil_now%coord_sys%position%z(:, jthe) = prof_eq
            equil_now%coord_sys%position%theta2d(jthe) = theta
        enddo
    endif

    end subroutine METRIC

!---------------------------------------------------------------------
    subroutine EQCYL()

!---------------------------------------------------------------------
! Quasi-cylindrical assignment: Called if IPEQL==-2
!
! In: RTOR, SHIFT, ABC, ELONG, TRIAN, NA1, NB1
! Out: NA, HRO, ROC, RHO(j), DRODA, VOLUM, IPOL, G33, GRADRO, G11, G22, SLAT
!  BDB02, B0DB2, BDB0, FOFB, BMAXT, BMINT, DRODA, GRADRO
!---------------------------------------------------------------------

    use status, only: RHO, XRHO, VR, VRS, AMETR, SHIF, SHIV, &
        ELON, TRIA, SLAT, G11, G22, G33, G41, G42, G43, G44, G45, &
        BDB0, BDB02, B0DB2, IPOL, MU, &
        FOFB, BMAXT, BMINT, DRODA, GRADRO, VOLUM
    use scalars, only: VOLUME, RTOR, BTOR, &
        ABC, HRO, ROC, FTO, ROWALL, NA, NA1, NB1
    use numerical_tools, only: integr
    use debugger, only: markloc, debug

    integer :: j

    call markloc('EQCYL', debug_lev=3*debug)

    VOLUME = GP2*RTOR*GP*ABC**2
    ROC = ABC
    FTO = GP*BTOR*ROC**2
    HRO = ABC/(NA1 - 0.5)

    do j=1, NB1
        RHO(J) = XRHO(J)*ROC
        if (RHO(j) <= ROWALL) NA = j
        G22(J)   = J*HRO
        VR(J)    = GP2_sq * RTOR*RHO(J)
        VRS(j)   = GP2_sq * RTOR*G22(J)
        AMETR(J) = RHO(J)
        SHIF(J)  = 0.
        SHIV(J)  = 0.
        ELON(J)  = 1.
        TRIA(J)  = 0.
        IPOL(J)  = 1.
        G33(J)   = 1.
        G11(J)   = VRS(j)
        SLAT(J)  = VRS(j)
        BDB02(j) = 1. + (RHO(j)*MU(j)/RTOR)**2
        B0DB2(j) = 1./BDB02(j)
        BDB0(j)  = sqrt(BDB02(j))
        FOFB(j)  = 1.
        BMAXT(j) = BTOR*BDB0(j)
        BMINT(j) = BMAXT(j)
        DRODA(j) = 1.
        GRADRO(j)= 1.
        G41(J)   = 1.
        G42(J)   = GRADRO(J)
        G43(J)   = GRADRO(J)
        G44(J)   = G11(J)/VRS(J)
        G45(J)   = G11(J)/VRS(J)
    enddo
    if (NA < NB1) then
        NB1 = NA
    endif
    NA = NA1 - 1
    call INTEGR(RHO, 1, VR, VOLUM, NB1)
    VOLUME = VOLUM(NA1)

    end subroutine EQCYL

!---------------------------------------------------------------------
    subroutine extrap_fields_flat()

    use scalars, only: NA1, NAB
    use status, only: SHEAR, BDB02, B0DB2, BMAXT, BMINT, BDB0, FOFB

    integer :: j

    do J=NA1+1, NAB
        SHEAR(j) = SHEAR(NA1)
        BDB02(j) = BDB02(NA1)
        B0DB2(j) = B0DB2(NA1)
        BMAXT(j) = BMAXT(NA1)
        BMINT(j) = BMINT(NA1)
        BDB0(j)  = BDB0(NA1)
        FOFB(j)  = FOFB(NA1)
    enddo

    end subroutine extrap_fields_flat

!---------------------------------------------------------------------
    subroutine EQGUESS()

!---------------------------------------------------------------------
! Guessed equibrium: Called if IPEQL==0 or data_initiation @ 1st_entry
!
! In: RTOR, SHIFT, ABC, ELONG, TRIAN, NA1, NB1, HRO
! Out: NA, RHO(j), DRODA, VOLUM, IPOL, G33, GRADRO, G11, G22, SLAT,
!   G41 G42 G43 G44 G45
!
! Calling:
!         SETGEO -> SHIF, SHIV, ELONG, TRIA, DRODA, AMETR
! new_grid -> Takes ROC, HRO, NB1, NA1, NA=NA1-1, AB, ABC, AMETR(NA1)
!     Returns:
!     NA1, NA=NA1-1, NAB, RHO(NA1)=ROC, AMETR(j>NA1)
!---------------------------------------------------------------------

    use scalars, only: ROC, RTOR, SHIFT, ABC, ELONG, TRIAN, &
        FTO, BTOR, NA, NA1, NB1, HRO, VOLUME
    use status, only: RHO, VR, VRS, AMETR, SHIF, &
        ELON, TRIA, SLAT, G11, G22, G33, G41, G42, G43, G44, G45, &
        BDB0, BDB02, B0DB2, IPOL, MU, FP, SHEAR, &
        FOFB, BMAXT, BMINT, DRODA, GRADRO, VOLUM
    use numerical_tools, only: integr
    use debugger, only: markloc, debug

    integer :: j, j1
    double precision :: YA, YAS, YES, YDS, YDV, YR1

    call markloc('EQGUESS', debug_lev=3*debug)

    ROC = ROC3A(RTOR, SHIFT, ABC, ELONG, TRIAN)
    FTO = GP*BTOR*ROC**2
    NA  = NA1 - 1

! Compute NB1
    call new_grid()
! Define AMETR, SHIF, ELON, TRIA, SHIV

    call SETGEO()
    YDV = 0.

    do J=1, NB1
        if(j < nb1) then
            j1 = j + 1
        else
            j1 = j
        endif
        YA  = 0.5*(AMETR(J) + AMETR(J1))
        YES = 0.5*(ELON(J)  + ELON(J1))
        YDS = 0.5*(TRIA(J)  + TRIA(J1))
        YAS = 0.5*(SHIF(J)  + SHIF(J1))
        VOLUM(J) = GP*GP2*YA**2 * YES*(RTOR + YAS - 0.25*YA*YDS)
        VR(J) = (VOLUM(J) - YDV)/HRO
        YDV = VOLUM(J)
        YR1 = min(AMETR(j)/2., 1.d-2)
        YAS = ROC3A(RTOR, SHIF(j), AMETR(j),       ELON(j), TRIA(j))
        YDS = ROC3A(RTOR, SHIF(j), AMETR(j) - YR1, ELON(j), TRIA(j))
        DRODA(J) = (YAS - YDS)/YR1
    enddo
    VOLUM(NA1) = GP*GP2*ABC**2 * ELONG*(RTOR + SHIFT - 0.25*ABC*TRIAN)

! Input:  ROC, HRO, NB1, NA1, NA=NA1-1, AB, ABC, AMETR(NA1)
    call new_grid()! Output: NAB, RHO(NA1)=ROC, AMETR(j>NA1)

!----- Definition --------------------------- Approximation ----------
! gradRHO = DRODA
! <|grad(a)|> = sqrt(G1)/DRODA
! G1=<(gradRHO)**2>;     G1=DRODA**2
! G11=VR*G1=VR*<(gradRHO)**2>    G11=VR*G1
! G2=<(gradRHO/R)**2>*VR/4/pi**2;    G22=R*G2/J;
! G22=VR*R*<(gradRHO/R)**2>/(GP2)**2/IPOL; C22=G11/R/(GP2)**2/IPOL
! SLAT=VR*DRODA*<|gradA|>;    SLAT=VR*sqrt(G1)

    do J=1, NB1
        IPOL(J) = 1.
        G33(J) = (RTOR/(RTOR + SHIF(J)))**2
        GRADRO(J) = DRODA(J)
        if (j < NB1) VRS(j) = 0.5*(VR(J + 1) + VR(j))
        SLAT(J) = VRS(J)*DRODA(J)
        G11(J)  = VRS(J)*DRODA(J)**2
        G22(J)  = G11(J)/GP2_sq/(RTOR + SHIF(J))
        G41(J)  = 1.0
        G42(J)  = GRADRO(J)
        G43(J)  = GRADRO(J)
        G44(J)  = G11(J)/VRS(J)
        G45(J)  = G11(J)/VRS(J)
    enddo

!---------------------------------------------------------------------
! BDB02 - <B**2/B0**2>
! B0DB2 - <B0**2/B**2>    <(R/R0)^2>
! BMAXT - BMAXT
! BMINT - BMINT
! BDB0  - <B/BTOR>
! FOFB  - <(BTOR/B)**2*(1.-SQRT(1-B/Bmax)*(1+.5B/Bmax))>
! SHEAR -  d[ln(q)]/d[ln(rho)] (to replace fml/shear)
! Alternative definition for BDB02 (G.Pereverzev 10.02.99)
!---------------------------------------------------------------------

    do J=1, NA1
        if (j == 1) then
            SHEAR(J) = (FP(2) - FP(1))/(2.*MU(1) + 0.333*(MU(1) - MU(2)))
        elseif (j <= NA) then
            SHEAR(J) = (FP(j+1) - 2.*FP(j) + FP(j-1))/(MU(j+1) + MU(j))
        endif
        SHEAR(J) = 1. - SHEAR(J)/(GP*BTOR*HRO**2)
        YR1      = RHO(j)*G22(j)*(MU(j)/RTOR)**2
        BDB02(j) = (1. + YR1)*G33(J)*IPOL(J)**2
        B0DB2(j) = ((RTOR + SHIF(J))/RTOR)**2 + 0.75*(AMETR(j)/RTOR)**2
        BMAXT(j) = BTOR*RTOR/(RTOR + SHIF(j) - AMETR(j))
        BMINT(j) = BTOR*RTOR/(RTOR + SHIF(j) + AMETR(j))
        BDB0(j)  = (RTOR/(RTOR + SHIF(j)))
!  - <(BTOR/B)**2*(1.-SQRT(1-B/Bmax)*(1+.5B/Bmax))>
        YR1      = (RTOR + SHIF(j) - AMETR(j))/(RTOR + SHIF(j))
        FOFB(j)  = 1. - sqrt(YR1)*(1. + 0.5*YR1)
    enddo
    SHEAR(NA1) = SHEAR(NA)

    call extrap_fields_flat()
    call INTEGR(RHO, 1, VR, VOLUM, NA1)

    VOLUME = VOLUM(NA1)

    end subroutine EQGUESS

!----------------------------------------------------------
    logical function IFDEFX(XARNAM)
! Name exists in profxNames, and the array is defined

    use read_input, only: IFDFAX
    use json_vars, only: profxNames
    use parse_utils, only: str_in_list

    character(len=6), intent(in) :: XARNAM

    integer :: j

    j = str_in_list(XARNAM, profxNames)
    if (j > 0) then
        if (IFDFAX(j) > 0) IFDEFX = .true. ! True (X-array is defined)
    endif

    end function IFDEFX

!---------------------------------------------------------------------
    subroutine set_external_metric()

! Set external metric  (Pereverzev 10.02.2005)

    use scalars, only: RTOR, BTOR, ABC, ROC, HRO, HROX, &
        SHIFT, ELONG, TRIAN, VOLUME, NA, NA1, NAB, updwn, ipart
    use status, only: SHIF, ELON, TRIA, SHX, ELX, TRX, &
        G11, G22, G33, G11X, G22X, G33X, GRADRO, DRODA, DRODAX, &
        IPOL, IPOLX, VR, VRS, VRX, RHO, XRHO, AMETR, SLAT, SLATX, &
        BDB0, BDB02, B0DB2, BMINT, BMAXT, FOFB, VOLUM, SHEAR, FP, MU, &
        shiv, squarn, shivx, squax
    use debugger, only: markloc, debug
    use numerical_tools, only: integr
    use read_input, only: flightsim

    integer :: j
    double precision :: YNF, YR1

    call markloc('set_external_metric', debug_lev=3*debug)

    if (.not. flightsim) then
        YNF = RTOR*GP2_sq
        do J=1, NA1
            if (IFDEFX('SHX   ')) then
                SHIF(J) = SHX(j)
            else
                SHIF(J) = SHIFT
            endif
            if (IFDEFX('SHIVX ')) then
                SHIV(J) = SHIVX(j)
            else
                SHIV(J) = 0.
            endif
            if (IFDEFX('SQUAX ')) then
                SQUARN(J) = SQUAX(j)
            else
                SQUARN(J) = 0.
            endif
            if (IFDEFX('ELX   ')) then
                ELON(J) = ELX(j)
            else
                ELON(J) = 1.
            endif
            if (IFDEFX('TRX   ')) then
                TRIA(J) = TRX(j)
            else
                TRIA(J) = 0.
            endif
            if (IFDEFX('G33X  ')) then
                G33(J) = G33X(j)
            else
                G33(J) = (RTOR/(RTOR + SHIFT))**2
            endif
            if (IFDEFX('IPOLX ')) then
                IPOL(J) = IPOLX(j)
            else
                IPOL(J) = 1.
            endif
            if (IFDEFX('VRX   ')) then
                VR(J) = VRX(j)
            else
                VR(J) = YNF*RHO(j)/(IPOL(j)*G33(j))
            endif
        enddo
    else if (ipart == 1) then   ! for initialisation
        call eqguess()
        return
    endif

! Compute new ROC

    ROC = VR(NA1)/GP2_sq * G33(NA1)/RTOR
    RHO(1: NA1) = XRHO(1: NA1)*ROC

    if (debug > 0) then
        write(*, *) 'metric', RHO(1: 10)
        write(*, *) VR(1: 10)
        write(*, *) G33(1: 10)
        write(*, *) IPOL(1: 10)
    endif
    HRO  =  RHO(2) - RHO(1)
    HROX = (RHO(2) - RHO(1))/ROC

    if (flightsim) then
        do j=1, NA1
            VRS(j)   = 0.5*(VR(J+1) + VR(j))
            SLAT(J)  = 0.5*(SLAT(J+1) + SLAT(j))
            G11(J)   = 0.5*(G11(j) + G11(j+1))
            G22(J)   = 0.5*(G22(j) + G22(j+1))
            DRODA(J) = 0.5*(DRODA(j) + DRODA(j+1))
         enddo
    else
! Flux grid: j*h
        do J=1, NA
            VRS(j) = 0.5*(VR(J+1) + VR(j))
            if (IFDEFX('SLATX ')) then
                SLAT(J) = 0.5*(SLATX(J+1) + SLATX(j))
            else
                SLAT(J) = VRS(j)
            endif
            if (IFDEFX('G11X  ')) then
                G11(J) = 0.5*(G11X(j) + G11X(j+1))
            else
                G11(J) = VRS(j)
            endif
            if (IFDEFX('G22X  ')) then
                G22(J) = 0.5*(G22X(j) + G22X(j+1))
            else
                G22(J) = RTOR*VRS(j)/(GP2*(RTOR + SHIFT))**2
            endif
            if (IFDEFX('DRODAX')) then
                DRODA(J) = 0.5*(DRODAX(j) + DRODAX(j+1))
            else
                DRODA(J) = 1.
            endif
        enddo
    endif

! Linear extrapolation
    SLAT(NA1)  = 1.5*SLAT(NA)  - 0.5*SLAT(NA-1)
    VRS(NA1)   = 1.5*VRS(NA)   - 0.5*VRS(NA-1)
    G11(NA1)   = 1.5*G11(NA)   - 0.5*G11(NA-1)
    G22(NA1)   = 1.5*G22(NA)   - 0.5*G22(NA-1)
    DRODA(NA1) = 1.5*DRODA(NA) - 0.5*DRODA(NA-1)

! Compute new minor radius => better just take from data?
    call INTEGR(RHO(1: NA1), 1, 1./DRODA(1: NA1), AMETR(1: NA1), NA1)
    ABC   = AMETR(NA1)
    ELONG = ELON(NA1)
    TRIAN = TRIA(NA1)
    SHIFT = SHIF(NA1)
    UPDWN = SHIV(NA1)

    do J=1, NA1
        if (j == 1) then
            SHEAR(J) = (FP(2) - FP(1))/(2.*MU(1) + 0.333*(MU(1) - MU(2)))
        elseif (j <= NA) then
            SHEAR(J) = (FP(j+1) - 2.*FP(j) + FP(j-1))/(MU(j+1) + MU(j))
        endif
        SHEAR(J) = 1. - SHEAR(J)/(GP*BTOR*HRO**2)
        GRADRO(j) = DRODA(J)
        YR1 = RHO(j)*G22(j)*(MU(j)/RTOR)**2
        BDB02(j) = (1. + YR1)*G33(J)*IPOL(J)**2
        B0DB2(j) = ((RTOR + SHIF(J))/RTOR)**2 + 0.75*(AMETR(j)/RTOR)**2
        BMAXT(j) = BTOR*RTOR/(RTOR + SHIF(j) - AMETR(j))
        BMINT(j) = BTOR*RTOR/(RTOR + SHIF(j) + AMETR(j))
        BDB0(j)  = (RTOR/(RTOR + SHIF(j)))
        YR1      = (RTOR + SHIF(j) - AMETR(j))/(RTOR + SHIF(j))
        FOFB(j)  = 1. - sqrt(YR1)*(1. + 0.5*YR1)
    enddo
    SHEAR(NA1) = SHEAR(NA)
    if (NA1 < NAB) then
        do J=NA1+1, NAB
            SHIF(J) = SHIFT
            SHIV(J) = 0.
            SQUARN(J) = 0.
            ELON(J) = 1.
            TRIA(J) = 0.
            G33(J) = (RTOR/(RTOR + SHIFT))**2
            IPOL(J) = 1.
            VR(J) = YNF*RHO(j)/(IPOL(j)*G33(j))
            AMETR(J) = RHO(J)
            VRS(j) = 0.5*(VR(J+1) + VR(j))
            G11(J) = VRS(j)
            G22(J) = RTOR*VRS(j)/(RTOR + SHIFT)**2
            DRODA(J) = 1.
            SLAT(J)  = VRS(J)*DRODA(J)
        enddo
        call extrap_fields_flat()
        call INTEGR(RHO, 1, VR, VOLUM, NA1) ! Compute volume(rho) on shifted grid integrating VR
        call new_grid()! The RHO-grid and NA, NA1 are updated
        VOLUME = VOLUM(NA1)
    endif

    end subroutine set_external_metric

!---------------------------------------------------------------------
    subroutine extmetric_input()

    use scalars, only: NA1
    use status, only: SHIF, ELON, TRIA, G33, IPOL, VR, SLAT, G11, G22, &
        DRODA, SHIV, SQUARN

    open(32, file='input_metric.dat')
    read(32, '(5555E25.11)') SHIF(1:NA1), elon(1:NA1), tria(1:NA1), &
        g33(1:NA1), ipol(1:NA1), vr(1:NA1), slat(1:NA1), g11(1:NA1), &
        g22(1:NA1), droda(1:NA1), shiv(1:NA1), squarn(1:NA1)
    close(32)

    end subroutine extmetric_input

!---------------------------------------------------------------------
    subroutine set_external_metric_2()

! Set external metric

    use scalars, only: RTOR, BTOR, ABC, ROC, HRO, HROX, &
        SHIFT, ELONG, TRIAN, VOLUME, NA, NA1, NAB, updwn
    use status, only: SHIF, ELON, TRIA, &
        G11, G22, G33, GRADRO, DRODA, &
        IPOL, VR, VRS, RHO, XRHO, AMETR, SLAT, &
        BDB0, BDB02, B0DB2, BMINT, BMAXT, FOFB, VOLUM, SHEAR, FP, MU, SHIV, SQUARN
    use debugger, only: markloc, debug
    use numerical_tools, only: integr
    use read_input, only: flightsim

    integer :: j
    double precision :: YNF, YR1

    call markloc('set_external_metric_2', debug_lev=3*debug)

    if (.not. flightsim) then
        call extmetric_input()
    endif

    YNF = RTOR*GP2_sq

! Compute new roc

    ROC = VR(NA1)/GP2_sq * G33(NA1)/RTOR
    RHO(1: NA1) = XRHO(1: NA1)*ROC

    if (debug > 0) then
        write(*, *) 'metric', RHO(1: 10), ROC
        write(*, *) VR(1: 10)
        write(*, *) G33(1: 10)
        write(*, *) IPOL(1: 10)
    endif
    HRO  =  RHO(2) - RHO(1)
    HROX = (RHO(2) - RHO(1))/ROC

! Flux grid: j*h
    do J=1, NA
        VRS(j) = 0.5*(VR(J+1) + VR(j))
    enddo

! Linear extrapolation
    SLAT(NA1)  = 1.5*SLAT(NA)  - 0.5*SLAT(NA-1)
    VRS(NA1)   = 1.5*VRS(NA)   - 0.5*VRS(NA-1)
    G11(NA1)   = 1.5*G11(NA)   - 0.5*G11(NA-1)
    G22(NA1)   = 1.5*G22(NA)   - 0.5*G22(NA-1)
    DRODA(NA1) = 1.5*DRODA(NA) - 0.5*DRODA(NA-1)

! Compute new minor radius
    call INTEGR(RHO(1: NA1), 1, 1./DRODA(1: NA1), AMETR(1: NA1), NA1)
    ABC   = AMETR(NA1)
    ELONG = ELON(NA1)
    TRIAN = TRIA(NA1)
    SHIFT = SHIF(NA1)
    UPDWN = SHIV(NA1)

    do J=1, NA1
        if (j == 1) then
            SHEAR(J) = (FP(2) - FP(1))/(2.*MU(1) + 0.333*(MU(1) - MU(2)))
        elseif (j <= NA) then
            SHEAR(J) = (FP(j+1) - 2.*FP(j) + FP(j-1))/(MU(j+1) + MU(j))
        endif
        SHEAR(J) = 1. - SHEAR(J)/(GP*BTOR*HRO**2)
        GRADRO(j) = DRODA(J)
        YR1 = RHO(j)*G22(j)*(MU(j)/RTOR)**2
        BDB02(j) = (1. + YR1)*G33(J)*IPOL(J)**2
        B0DB2(j) = ((RTOR + SHIF(J))/RTOR)**2 + 0.75*(AMETR(j)/RTOR)**2
        BMAXT(j) = BTOR*RTOR/(RTOR + SHIF(j) - AMETR(j))
        BMINT(j) = BTOR*RTOR/(RTOR + SHIF(j) + AMETR(j))
        BDB0(j)  = (RTOR/(RTOR + SHIF(j)))
        YR1      = (RTOR + SHIF(j) - AMETR(j))/(RTOR + SHIF(j))
        FOFB(j)  = 1. - sqrt(YR1)*(1. + 0.5*YR1)
    enddo
    SHEAR(NA1) = SHEAR(NA)
    if (NA1 < NAB) then
        do J=NA1+1, NAB
            SHIF(J) = SHIFT
            SHIV(J) = 0.
            SQUARN(J) = 0.
            ELON(J) = 1.
            TRIA(J) = 0.
            G33(J) = (RTOR/(RTOR + SHIFT))**2
            IPOL(J) = 1.
            VR(J) = YNF*RHO(j)/(IPOL(j)*G33(j))
            AMETR(J) = RHO(J)
            VRS(j) = 0.5*(VR(J+1) + VR(j))
            G11(J) = VRS(j)
            G22(J) = RTOR*VRS(j)/(RTOR + SHIFT)**2
            DRODA(J) = 1.
            SLAT(J)  = VRS(J)*DRODA(J)
        enddo
        call extrap_fields_flat()
        call INTEGR(RHO, 1, VR, VOLUM, NA1) ! Compute volume(rho) integrating VR
        call new_grid() ! The RHO-grid and NA, NA1 are updated
        VOLUME = VOLUM(NA1)
    endif

    end subroutine set_external_metric_2

!---------------------------------------------------------------------
    double precision function ROC3A(Rmaj, shaf_shift, a_min, elongation, triangularity)

!---------------------------------------------------------------------
! ROC3A [m]: Analytical formula for the "rho_tor" in vacuum
!      (Pereverzev 24.03.00)
! Input: Rmaj - Major radius [m]
!  shaf_shift - Shafranov shift [m]
!  a_min - Minor radius [m]
!  elongation - Elongation [d/l]
!  triangularity - Triangularity [d/l]
! Output:
!  ROC3A - Dimensional toroidal "rho" [m]
!---------------------------------------------------------------------

    use debugger, only: markloc, debug

    double precision, intent(in) :: Rmaj, shaf_shift, a_min, elongation, triangularity

    double precision :: YD1, YT2, YGE, Y1, Y2, Y3, Y4, Y5, Y6

    call markloc('ROC3A', debug_lev=2*debug)

    YT2 = 2.*triangularity
    YGE = a_min/(Rmaj + shaf_shift)
    if (abs(YGE) > 1.d0) then
        write(*, *) " >>> Error >>> ROC3A >>> Illegal input: R+Delta < a"
        write(*, '(1P, 16X, 2(A, E10.3))') "R+Delta =", Rmaj + shaf_shift, ",    a =", a_min
        write(*, *) Rmaj, shaf_shift, a_min, elongation, triangularity
        STOP
     endif

    if (YGE < 0) then
        write(*, *) " >>> Error >>> ROC3A: Illegal input"
        write(*, *) Rmaj, shaf_shift, a_min, elongation, triangularity
        STOP
    else if (YGE == 0) then
        ROC3A = 0.
    else
        YD1 = 1.+YT2*(YT2-2./YGE)
        if (YD1 >= 0) then
            YD1 = sqrt(YD1)
            Y1 = (2. - YT2*YGE)*YD1/(1. + YD1)
            Y2 = YGE - YT2
            Y3 = sqrt(1. - YT2*YGE + YGE*YD1)
            Y4 = sqrt(1. - YT2*YGE - YGE*YD1)
            Y5 = sqrt(1. + YGE)
            Y6 = sqrt(1. - YGE)
            Y1 = (Y1 + Y2)/(Y3 + (1. - YT2)*Y5) - (Y1 - Y2)/(Y4 + (1. + YT2)*Y6)
            Y1 = Y1/sqrt(2.*(1. - YT2*YGE + Y5*Y6))
        else
            Y5 = sqrt(1. + YGE)
            Y6 = sqrt(1. - YGE)
            Y1 = (1. - YT2)*Y5 + (1. + YT2)*Y6
            Y1 = (1. - Y1/sqrt(2.*(1. - YT2*YGE + Y5*Y6)))/YT2
        endif
        ROC3A = 2.*sqrt(Rmaj*a_min*elongation*Y1)
    endif

    end function ROC3A

!---------------------------------------------------------------------
    subroutine A2EMEQ(jexit)

!---------------------------------------------------------------------
! Module call sequence in A2EMEQ
!
!    EQ -> mapping -> EMEQ -> new_grid -> mapping
!
!                                         | -> EQGB3
!                                         |                 | -> EQLVU3
!    inside EMEC:          EMEQ -> | -> EQAB3 -> | -> EQK3
!                                         |                 | -> EQC1
!                                         | -> EQPPAB
!---------------------------------------------------------------------
! This solver is called if IPEQL == 3
! EMEQ grid is defined as min(NA1, NEQUIL, NP)
! For (NEQUIL = 0 ) metric is prescribed by a simple formula
! For (NEQUIL = -1) metric is taken from a data file
! For (NEQUIL = 1 ) metric is frozen (can be used interactively)
!---------------------------------------------------------------------

    use emeq_mod, only: NP, emeq
    use scalars, only: HRO, ABC, ROC, RTOR, BTOR, IPL, &
         TIME, ELONG, TRIAN, SHIFT, UPDWN, VOLUME, &
         NA1, NA, NAB, NEQUIL
    use status, only: NRD, TE, TI, CU, SHEAR, SHIV, &
         RHO, AMETR, EQPF, EQFF, IPOL, MU, FP, SXHO, &
         SLAT, VOLUM, SHIF, ELON, TRIA, DRODA, GRADRO, VR, VRS, XRHO, &
         G11, G22, G33, G41, G42, G43, G44, G45, &
         BDB0, BDB02, B0DB2, BMAXT, BMINT, FOFB
    use numerical_tools, only: integr, qinterp, smooth
    use debugger, only: markloc
    use standard_functions, only: IINT

    double precision, parameter :: ACEQLB=1.d-6

    integer, intent(out) :: jexit

    integer :: NR_EQU, j, jp, jt, jcall
    double precision :: ALFA, Y1, Y2, YDA, YRO, YCB, TRIABC, BTOOO
    double precision, dimension(NRD) :: BA, BB, GR, GBD, GL, GSD, &
        A, B, C, D, BC, BD, XTR, BMOD_EQU, FOFB_EQU, GRDA_EQU, &
        X_EQU, B2B0_EQU, B0B2_EQU, BMAX_EQU, BMIN_EQU, VR_EQU, VRS_EQU, &
        G11_EQU, G22_EQU, G33_EQU, IPOL_EQU, DRODA_EQU, GRADRO_EQU

    save jcall
    data jcall/0/

!--------------------------------------------------
! Prepare input data for the 3M equilibrium solver:

    call markloc('A2EMEQ')

    jexit = 0

! Define NR_EQU as min(NA1, NEQIL, NP)
    NR_EQU = NEQUIL
    if (NR_EQU > NP) then
        write(*, '(A, I4)') " >>> Warning >>> Maximum size of the equilibrium grid is", NP
        write(*, '(17X, A, I4)') "The grid will be reduced to NP =", NP
        NR_EQU = NP
    endif

! Rise triangularity (initiation stage only)
    jt = 5
    TRIABC = TRIAN*ABC

! Rise pressure (initiation stage only)
    jp = 5
    if (jcall <= jp) then
        Y1 = jcall    ! 0 <= jcall <= jp
        Y1 = Y1/jp
        do J=1, NA1
            A(J)  = TE(J)
            B(J)  = TI(J)
            C(J)  = CU(J)
            TE(J) = Y1*A(J)
            TI(J) = Y1*B(J)
            if (jcall == 0) then
                CU(J) = (1 - (RHO(J)/ROC)**2)
            else
                CU(J) = Y1*C(J) + (1. - Y1)*(1. - (RHO(J)/ROC)**2)
            endif
        enddo
        YCB = IINT(CU, ROC)
        do J=1, NA1
            CU(J) = CU(J)*IPL/YCB
        enddo
        call RHSEQ()
        do J=1, NA1
            TE(J) = A(J)
            TI(J) = B(J)
            CU(J) = C(J)
        enddo
        jcall = jcall + 1
    endif

    YCB = RTOR/(RTOR + SHIFT)
    do J=1, NA1
        XTR(J) = AMETR(J)/AMETR(NA1)
        B(J) = EQPF(J)/YCB
        A(J) = B(J) + YCB*EQFF(J)
    enddo

    do J=1, NR_EQU
        X_EQU(J) = (j - 1.)/(NR_EQU - 1.)
    enddo
! From transport grid in "a" to equidistant grid in "a"
    ALFA = 0.001
    call SMOOTH(ALFA, XTR, A, NA1, X_EQU, BA, NR_EQU)
    call SMOOTH(ALFA, XTR, B, NA1, X_EQU, BB, NR_EQU)

    call EMEQ( &
! Input:
        BA, BB, RTOR + SHIFT, ABC, ELONG, TRIABC, NR_EQU, ACEQLB, &  ! relative accuracy
        BTOR*RTOR/(RTOR + SHIFT), IPL, &  ! Total plasma current
! output
        GR, GBD, GL, GSD, A, BD, B, BA, BB, BC, C, D, &
        B2B0_EQU, B0B2_EQU, BMAX_EQU, BMIN_EQU, BMOD_EQU, FOFB_EQU, GRDA_EQU, &
! input
        TIME)
    BTOOO = BTOR*RTOR/(RTOR + SHIFT)

!----------------------------------------------------------------------|
! Output: GR  - rho
!  GBD - shift
!  GL  - elongation
!  GSD - triangularity = \delta_Zakh = a*\delta_Astra
!  A    - <g^{11}>  = <[nabla(a)]^2>
!  BD  - <sqrt[g^{11}]> = <|nabla(a)|> flux surface area
!  B    - <g^{11}g^{33}> = <[nabla(a)/r]^2>
!  BA  - <g^{33}>  = <1/r^2>=G33/RTOR**2
!  BB  - IPOL*RTOR*BTOR = I
!  BC  - d\rho/da
!  C    - (dV/da)/(4\pi^2)
!  D    - V(a)
!                    B2B0_EQU - <B**2/B0**2>
!                    B0B2_EQU - <B0**2/B**2>
!                    BMAX_EQU - BMAXT
!                    BMIN_EQU - BMINT
!                    BMOD_EQU - <B/BTOR>
!                    FOFB_EQU - <(BTOR/B)**2*(1.-SQRT(1-B/Bmax)*(1+.5B/Bmax))>
!                    GRDA_EQU - <grad a>
!----------------------------------------------------------------------|

    call markloc('3-moment solver')

! NR_EQU <= 1 can be returned by EQAB3 via EMEQ
    if (NR_EQU > 10) then
        if (isnan(GR(NR_EQU))) then
            jexit = 2
            return
        endif
    endif

!---------------------------------------
! Define a new RHO-grid:
    YRO = sqrt(RTOR/(RTOR + SHIFT))
    ROC = YRO*GR(NR_EQU)  ! Define a new RHO_edge
! FTN = GP*BTN*ROC*ROC

    call new_grid() ! The RHO-grid and NA, NA1, HRO are updated
!---------------------------------------
! Define a new auxiliary (shifted) grid:
    Y2 = 0.5d0/ROC

    do J=1, NR_EQU
        DRODA_EQU(J)  = YRO*BC(J)
        X_EQU(J)      = GR(J)/GR(NR_EQU)
        G11_EQU(J)    = A(J)*DRODA_EQU(J)**2
        G22_EQU(J)    = B(J)*DRODA_EQU(J)**2
        G33_EQU(J)    = BA(J)*RTOR*RTOR
        VRS_EQU(J)    = GP2_sq*C(J)/DRODA_EQU(J)
        IPOL_EQU(J)   = BB(J)/RTOR/BTOR
        GRADRO_EQU(J) = BD(J)*DRODA_EQU(J)
        VR_EQU(J)     = VRS_EQU(j)
    enddo
    call qinterp(X_EQU(1:NR_EQU), VRS_EQU(1:NR_EQU), NR_EQU, SXHO(1: NA1), VRS(1: NA1), NA1)
    call SMOOTH(ALFA, X_EQU,    G11_EQU, NR_EQU, SXHO,    G11, NA1)
    call SMOOTH(ALFA, X_EQU,    G22_EQU, NR_EQU, SXHO,    G22, NA1)
    call SMOOTH(ALFA, X_EQU,    G33_EQU, NR_EQU, SXHO,    G33, NA1)
    call SMOOTH(ALFA, X_EQU,   IPOL_EQU, NR_EQU, SXHO,   IPOL, NA1)
    call SMOOTH(ALFA, X_EQU,  DRODA_EQU, NR_EQU, SXHO,  DRODA, NA1)
    call SMOOTH(ALFA, X_EQU, GRADRO_EQU, NR_EQU, SXHO, GRADRO, NA1)

! Multiply above quantities by linear in rho factors (i.e. f(0)=0)
    do J=1, NA1
        G11(J) = G11(J)*VRS(j)
        G22(J) = G22(J)/G33(j)*(RTOR/IPOL(j))**2
        if (j < NA1) then
            G22(J) = G22(J)*SXHO(j)*ROC
        else
            G22(J) = G22(J)*(NA*HRO + 0.5*HRO)
        endif
        SLAT(J) = GRADRO(J)*VRS(J)
    enddo

! Define VR, G33 and IPOL on the main transport grid:
    call qinterp(X_EQU(1:NR_EQU), VR_EQU(1:NR_EQU), NR_EQU, XRHO(1: NA1), VR(1: NA1), NA1)
    call SMOOTH(ALFA, X_EQU,  G33_EQU, NR_EQU, XRHO(1: NA1),  G33(1: NA1), NA1)
    call SMOOTH(ALFA, X_EQU, IPOL_EQU, NR_EQU, XRHO(1: NA1), IPOL(1: NA1), NA1)
    call SMOOTH(ALFA, X_EQU,      GBD, NR_EQU, XRHO(1: NA1), SHIF(1: NA1), NA1)
    call SMOOTH(ALFA, X_EQU,       GL, NR_EQU, XRHO(1: NA1), ELON(1: NA1), NA1)

! GL -> a
! GSD -> \delta^{ASTRA} (dimensionless)
    YDA = ABC/(NR_EQU - 1.)
    do j=2, NR_EQU
        A(j) = YDA*(J - 1.)
        B(J) = GSD(J)/A(J)
    enddo
    A(1) = 0.
    B(1) = 0.
    call qinterp(X_EQU(1:NR_EQU), B(1:NR_EQU), NR_EQU, XRHO(1: NA1),  TRIA(1: NA1), NA1)
    call qinterp(X_EQU(1:NR_EQU), A(1:NR_EQU), NR_EQU, XRHO(1: NA1), AMETR(1: NA1), NA1)

    do J=1, NA
        SHIF(J) = SHIFT + SHIF(J)
        if (j == 1) then
            SHEAR(J) = (FP(2) - FP(1))/(2.*MU(1) + 0.333*(MU(1) - MU(2)))
        else
            SHEAR(J) = (FP(j+1) - 2.*FP(j) + FP(j-1))/(MU(j+1) + MU(j))
        endif
        SHEAR(J) = 1. - SHEAR(J)/(GP*BTOR*HRO**2)
    enddo
    SHEAR(NA1) = SHEAR(NA)

    SHIV(1: NAB) = UPDWN

    call qinterp(X_EQU(1:NR_EQU), B2B0_EQU(1:NR_EQU), NR_EQU, XRHO(1: NA1), BDB02(1: NA1), NA1)
    call qinterp(X_EQU(1:NR_EQU), B0B2_EQU(1:NR_EQU), NR_EQU, XRHO(1: NA1), B0DB2(1: NA1), NA1)
    call qinterp(X_EQU(1:NR_EQU), BMAX_EQU(1:NR_EQU), NR_EQU, XRHO(1: NA1), BMAXT(1: NA1), NA1)
    call qinterp(X_EQU(1:NR_EQU), BMIN_EQU(1:NR_EQU), NR_EQU, XRHO(1: NA1), BMINT(1: NA1), NA1)
    call qinterp(X_EQU(1:NR_EQU), BMOD_EQU(1:NR_EQU), NR_EQU, XRHO(1: NA1), BDB0( 1: NA1), NA1)
    call qinterp(X_EQU(1:NR_EQU), FOFB_EQU(1:NR_EQU), NR_EQU, XRHO(1: NA1), FOFB( 1: NA1), NA1)

    BDB02 = BDB02*BTOOO**2 / BTOR**2
    BDB0  = BDB0*BTOOO/BTOR
    B0DB2 = B0DB2/BTOOO**2 * BTOR**2

    call extrap_fields_flat()
    call INTEGR(RHO, 1, VR, VOLUM, NA1)
    VOLUME = VOLUM(NA1)

!Efable add G41-G45 : these are on the shifted grid.
    do J=1, NA1
        G41(J) = 1.0
        G42(J) = GRADRO(J)
        G43(J) = GRADRO(J)
        G44(J) = G11(J)/VRS(J)
        G45(J) = G11(J)/VRS(J)
    enddo

    end subroutine A2EMEQ

!---------------------------------------------------------------------
    subroutine A2GSSOLVER(equil_solver)

    use read_input, only: machine
    use scalars, only: NEQUIL, MEQUIL, IPART, IPCTRL, TAU, NA, NA1, &
        RTOR, BTOR, IPL, HRO, ROC, ABC, &
        VOLUME, SHIFT, ELONG, UPDWN, TRIAN, &
        ITFBP, IPLFBE, IFBEY, ITREQ, ICIRCQ, ITFBE, &
        NB2EQL, TIME, LEQ, PSIFB, PSPLEX, PSIEXT, IPEQL, IPROT
    use status, only: G11, G22, G22E, G33, G33E, G41, G42, G43, G44, G45, &
        FP, IPOL, MU, SHEAR, &
        AMETR, VR, VRS, SLAT, GRADRO, DRODA, &
        NE, TE, NI, TI, MRHO, PBLON, PBPER, PFAST, EQPF, EQFF, &
        BMAXT, BMINT, BDB02, BDB0, B0DB2, FOFB, &
        VOLUM, SHIF, ELON, TRIA, XRHO, AREAT, PERIM, SHIV, SQUARN, VTOR
    use debugger, only: markloc
    use imas_ids, only: type_equilibrium
    use read_input, only: raw_boundary, raw_cCoil
    use gs_solver, only: gssolver

    integer, parameter :: itfbe_ctrl=0

    integer, intent(in) :: equil_solver

    integer :: i, j, jneql, jntheta, jnbnd, j_rotation, n_coils, &
        neql, k_fixfree, no_circuit_eq
    double precision :: yrocnew, iplnew, ychipfp, yipl, yupdwn
    double precision, dimension(NA1) :: yg11, yg22, yg33, yvr, yvrs, yslat, yg41, &
        ygradro, yipol, ydroda, ypres, ybmaxt, ybmint, yfp, &
        ybdb02, ybdb0, yb0db2, yvolum, yametr, yshif, yelon, &
        ytria, yfofb, yeqpf, yeqff, yshiv, ysquare, omega_rot
    double precision, dimension(raw_cCoil%ncoils) :: yccoil, yvcoil
    double precision, dimension(1000) :: rbnd, zbnd

    type(type_equilibrium) :: equil_in, equil_out

    call markloc('A2GSSOLVER')

    jneql   = abs(NEQUIL)
    jntheta = abs(MEQUIL)
    n_coils = raw_cCoil%ncoils

    if (raw_boundary%n_theta == 0) then
        jnbnd = jntheta ! "NAMEXP BND" not found
        raw_boundary%n_theta = jnbnd
        allocate(raw_boundary%R(1, jnbnd), raw_boundary%Z(1, jnbnd))
    else
        jnbnd = raw_boundary%n_theta
    endif

! provide grid for t=TIME+TAU
    call BNDRY(rbnd(1: jnbnd), zbnd(1: jnbnd))

    if (raw_boundary%n_theta == 0) then
        raw_boundary%R(1, 1: jnbnd) = rbnd(1: jnbnd)
        raw_boundary%Z(1, 1: jnbnd) = zbnd(1: jnbnd)
    endif

    do j=1, n_coils
        yccoil(j) = CCOIL(j)
        yvcoil(j) = VCOIL(j)
    enddo

    iplnew = G22(NA)/RTOR/0.4/GP * (FP(NA1) - FP(NA))/HRO * IPOL(NA1)

    if (ITFBP /= 0) IPLFBE = iplnew      ! current for free boundary equilibrium
    if (IPART == 1) then
        iplnew = IPL            ! if in initialization mode, use plasma current
    endif

    if (.not. plasma_up) then
        iplnew = IPL
        IPLFBE = IPL
    endif

    if (LEQ(4) <= 0) iplnew = IPL !if CU:AS, current is assigned from model file

    yipl  = iplnew

    j_rotation = 0
    omega_rot = 0.
    if (abs(IPROT) == 3 .or. abs(IPROT) == 4) then
        j_rotation = 1
        omega_rot(1:NA1) = VTOR(1:NA1)/(RTOR+SHIF(1:NA1) + AMETR(1:NA1)) ! Flux function Omega from Vtor_LFS / R_LFS
    endif

    do j=1, NA1
        yfp(j) = FP(j)
        yametr(j) = AMETR(j)

!Efable Reput values since now they are used for ff' computation
        yg11(j)    = G11(j)*VRS(j)                                    !g11 = <(grad(V)^2)>
        yg22(j)    = G22(j)*GP2_sq*IPOL(j)/RTOR*VRS(j)    !g22 = <(grad(V)/R)^2>
        yg33(j)    = G33(j)/(RTOR**2)                             !g33 = <1/R^2>
        yvr(j)     = VR(j)
        yvrs(j)    = VRS(j)
        yslat(j)   = SLAT(j)
        ygradro(j) = GRADRO(j)*VRS(j)                !gradro = <(grad(V))>
        yipol(j)   = IPOL(j)*RTOR*BTOR                 !ipol = R*Bphi
        ypres(j)   = 1.60218E-3*((NE(j)*TE(j) + NI(j)*TI(j)) + &
            NB2EQL*(0.5*PBLON(j) + 0.5*PBPER(j)) + PFAST(j))  !thermal + fast ions  from NBI + fast alpha in keV/m^3 *1e19 to MJ/m^3
        yeqpf(j)   = EQPF(j)
        yeqff(j)   = EQFF(j)
        yvolum(j)  = VOLUM(j)
        yshif(j)   = SHIF(j)
        yelon(j)   = ELON(j)
        ytria(j)   = TRIA(j)
    enddo
    ychipfp = FP(NA1)

! call to equil
    i = 1    !fbe is off
    if (IFBEY >= 1) i = 2    !fbe is on
    if (IPART == 1 ) i = 1    !fbe is off

    if (ifbey > 0 .and. .not. plasma_up) then
        equil_in%global_param%i_plasma = IPLFBE*1e6   !itm is in A
        k_fixfree = 1
        neql = 100
        no_circuit_eq = 0
        if (equil_solver == 101) then
            call feqis_main(n_coils, vcoil(1:n_coils), neql, k_fixfree, no_circuit_eq, 0, machine, equil_in, equil_out)
        endif
    else
        call GSSOLVER( &
! Input:
            equil_solver, jneql, jntheta, jnbnd, NA1, rbnd, zbnd, XRHO(1: NA1), RTOR, BTOR, &
            ROC, yfp, ypres, VOLUME, n_coils, yvcoil, i, ITREQ, &
            TAU, ITFBP, ICIRCQ, IPCTRL, IFBEY, &
            TIME, ychipfp, PSIFB, &
            omega_rot, j_rotation, TI(1: NA1), NI(1: NA1), MRHO(1: NA1), &
! Output:
            yrocnew, yipl, yg11, yg41, yg22, yg33, G22E(1: jneql), G33E(1: jneql), &
            yeqpf, yeqff, yvr, yvrs, yslat, ygradro, yipol, ybmaxt, ybmint, &
            ybdb02, ybdb0, yb0db2, ydroda, yvolum, yametr, yupdwn, yshif, yelon, ytria, &
            yfofb, AREAT(1: NA1), PERIM(1: NA1), yshiv, ysquare)

        ROC  = YROCNEW  ! Define a new RHO_edge

        do j=1, NA1
            VR(j)    = yvr(j)
            VRS(j)   = yvrs(j)
            SLAT(j)  = yslat(j)
            BMAXT(j) = ybmaxt(j)
            BMINT(j) = ybmint(j)
            BDB02(j) = ybdb02(j)
            BDB0(j)  = ybdb0(j)
            B0DB2(j) = yb0db2(j)
            DRODA(j) = ydroda(j)
            IPOL(j)  = yipol(j)
            G11(j)   = yg11(j)
            G22(j)   = yg22(j)
            G33(j)   = yg33(j)
            GRADRO(j)= ygradro(j)
            VOLUM(j) = yvolum(j)
            AMETR(J) = yametr(J)
            SHIF(J)  = yshif(J)
            ELON(J)  = yelon(J)
            TRIA(J)  = ytria(J)
            G41(J)   = yg41(J)
            G42(J)   = GRADRO(J)
            G43(J)   = GRADRO(J)
            G44(J)   = G11(J)/VRS(J)
            G45(J)   = G11(J)/VRS(J)
            FOFB(J)  = yfofb(J)
            FP(J)    = yfp(J)      ! due to adiabatic compression done in the code
            EQPF(J)  = yeqpf(J)    ! due to adiabatic compression done in the code
            EQFF(J)  = yeqff(J)    ! due to adiabatic compression done in the code
        enddo

        if (raw_boundary%nt > 0 .or. TIME >= ITFBE .or. use_ext_bnd) then
            UPDWN = yupdwn
            ABC   = yametr(NA1)
            ELONG = ELON(NA1)
            TRIAN = TRIA(NA1)
            SHIFT = SHIF(NA1)
        endif
        if (itfbe_ctrl > 0) then
            UPDWN = yupdwn
            ABC   = yametr(NA1)
            ELONG = ELON(NA1)
            TRIAN = TRIA(NA1)
            SHIFT = SHIF(NA1)
        endif

! Deallocate equil_out%metric_coefs%g1 & co
        call new_grid() ! The RHO-grid and NA, NA1, HRO are updated, also AMETR(NA1) = ABC is done there

        VOLUM(NA1) = yvolum(NA1)
        G22 = G22/VRS*RTOR/GP2_sq/IPOL
        G11 = G11/VRS
        GRADRO = GRADRO/VRS
        DRODA  = DRODA/VRS
        PSPLEX = equil_out%global_param%psplex
        PSIEXT = -GP2*equil_out%global_param%psiext

        if (IPEQL == 5) then  ! FEQIS
            PSPLEX = PSPLEX/(0.4*GP*RTOR*ROC)*0.5*(G22(NA) + G22(NA1)) ! If LEXT only
        endif

        do J=1, NA
            if (j == 1) then
                SHEAR(J) = (FP(2) - FP(1))/(2.*MU(1) + 0.333*(MU(1) - MU(2)))
            else
                SHEAR(J) = (FP(j+1) - 2.*FP(j) + FP(j-1))/(MU(j+1) + MU(j))
            endif
            SHEAR(J) = 1. - SHEAR(J)/(GP*BTOR*HRO**2)
        enddo
        SHEAR(NA1) = SHEAR(NA)
        do J=1, NA1
            SHIV(J) = yshiv(J)
            SQUARN(J) = ysquare(J)
        enddo

        call extrap_fields_flat()
        VOLUME = VOLUM(NA1)
    endif

    end subroutine A2GSSOLVER

!---------------------------------------------------------------------
    subroutine BNDRY(RPB, ZPB)

!---------------------------------------------------------------------
! If 3M solver is used the subroutine is not called.
! Otherwise, if a general equilibrium solver, equil code, is called
! then
! 1) In case of the plasma boundary defined by 3 moments,
!     this subroutine writes 8 points on the boundary into arrays BNDR, BNDZ
!     and into arrays RPB(1:n_bnd), ZPB(1:n_bnd)
! 2) If the plasma boundary is defined by a data file then
!     this subroutine uses the arrays BNDR, BNDZ as an input and
!     produces output in [time dependent] arrays RPB, ZPB
!---------------------------------------------------------------------
!  call from ESC:
!  call BNDRY(RPB, ZPB)
!  call from equil:
!  call BNDRY(RZPB, RZPB(n_bnd+1))
!---------------------------------------------------------------------

    use read_input, only: raw_boundary
    use scalars, only: TIME, RTOR, SHIFT, ABC, TRIAN, UPDWN, ELONG

    double precision, intent(out) :: RPB(*), ZPB(*)

    integer :: j, j1, jt, nt_bnd, n_bnd
    double precision :: ydt, yd1, yd2, yfi
    double precision, dimension(:, :), allocatable :: ext_bnd_in ! 50 , 2 boundary values R, Z

    nt_bnd = raw_boundary%nt
    n_bnd  = raw_boundary%n_theta

    if (nt_bnd <= 1) then
        if (nt_bnd == 0) then  ! No input group "NAMEXP BND" found
            if (n_bnd == 0) then
                n_bnd = 8       ! call from ESC
                yd1 = 0.75      ! sin^2(pi/3)
                yd2 = 0.5       ! cos(pi/3)
                ydt = sqrt(yd1) ! sin(pi/3)
                RPB(1) = RTOR + SHIFT - ABC*TRIAN
                ZPB(1) = UPDWN + ABC*ELONG
                RPB(2) = RTOR + SHIFT - ABC*TRIAN
                ZPB(2) = UPDWN - ABC*ELONG
                RPB(3) = RTOR + SHIFT - ABC
                ZPB(3) = UPDWN
                RPB(4) = RTOR + SHIFT + ABC
                ZPB(4) = UPDWN
                RPB(5) = RTOR + SHIFT - ABC*(TRIAN*yd1 + yd2)
                ZPB(5) = UPDWN + ABC*ELONG*ydt
                RPB(6) = RTOR + SHIFT - ABC*(TRIAN*yd1 - yd2)
                ZPB(6) = UPDWN + ABC*ELONG*ydt
                RPB(7) = RTOR + SHIFT - ABC*(TRIAN*yd1 - yd2)
                ZPB(7) = UPDWN - ABC*ELONG*ydt
                RPB(8) = RTOR + SHIFT - ABC*(TRIAN*yd1 + yd2)
                ZPB(8) = UPDWN - ABC*ELONG*ydt
            else
                if (use_ext_bnd) then
                    if (.not. allocated(ext_bnd_in)) allocate(ext_bnd_in(n_bnd, 2))
                    do j=1, n_bnd
                        ZPB(j) = ext_bnd_in(j, 2)
                        RPB(j) = ext_bnd_in(j, 1)
                    enddo
                else
                    do j=1, n_bnd
                        YFI = GP2*(j - 1.)/n_bnd
                        YD1 = sin(YFI)
                        ZPB(j) = UPDWN + ABC*ELONG*YD1
                        RPB(j) = RTOR + SHIFT + ABC*(cos(YFI) - TRIAN*YD1**2)
                    enddo
                endif
            endif
        else ! nt_nbd = 1
            do j=1, n_bnd
                RPB(j) = raw_boundary%R(1, j)
                ZPB(j) = raw_boundary%Z(1, j)
            enddo
        endif
        return
    endif

    if (TIME <= raw_boundary%time(1)) then ! Take bnd at time=t1
        do j=1, n_bnd
            RPB(j) = raw_boundary%R(1, j)
            ZPB(j) = raw_boundary%Z(1, j)
        enddo
        return
    endif
    if (TIME >= raw_boundary%time(nt_bnd)) then ! Take bnd at time=t_nt_bnd
        do j=1, n_bnd
            RPB(j) = raw_boundary%R(nt_bnd, j)
            ZPB(j) = raw_boundary%Z(nt_bnd, j)
        enddo
        return
    endif

! Time interpolation of boundary R(t, theta), Z(t, theta)
    do j=1, nt_bnd
        if (TIME > raw_boundary%time(j)) jt = j
    enddo
    if (jt == nt_bnd) write(*, *) "OGOGO"
    ydt = raw_boundary%time(jt+1) - raw_boundary%time(jt)
    yd1 = (TIME - raw_boundary%time(jt))/ydt
    yd2 = (TIME - raw_boundary%time(jt+1))/ydt
    do j=1, n_bnd
        RPB(j) = yd1*raw_boundary%R(jt+1, j) - yd2*raw_boundary%R(jt, j)
        ZPB(j) = yd1*raw_boundary%Z(jt+1, j) - yd2*raw_boundary%Z(jt, j)
    enddo
    if (n_bnd > 12) return

!----------------------------------------------------------------------|
! The order is essential
!  Top(1), bottom(2), inward(3), outward(4), ...
    YDT = -1.
    do j=1, n_bnd
        if (ZPB(j) > YDT) then
            YDT = ZPB(j)
            j1 = j
        endif
    enddo
    YD1 = RPB(1)
    YD2 = ZPB(1)
    RPB(1)  = RPB(j1)
    ZPB(1)  = ZPB(j1)
    RPB(j1) = YD1
    ZPB(j1) = YD2
! Bottom(2)
    YDT = 1.
    do j=2, n_bnd
        if (ZPB(j) < YDT) then
            YDT = ZPB(j)
            j1 = j
        endif
    enddo
    YD1 = RPB(2)
    YD2 = ZPB(2)
    RPB(2) = RPB(j1)
    ZPB(2) = ZPB(j1)
    RPB(j1) = YD1
    ZPB(j1) = YD2
! Inward(3)
    YDT = 1000.
    do j=3, n_bnd
        if (RPB(j) < YDT) then
            YDT = RPB(j)
            j1 = j
        endif
    enddo
    YD1 = RPB(3)
    YD2 = ZPB(3)
    RPB(3)  = RPB(j1)
    ZPB(3)  = ZPB(j1)
    RPB(j1) = YD1
    ZPB(j1) = YD2
! Outward(4)
    YDT = -1.
    do j=4, n_bnd
        if (RPB(j) > YDT) then
            YDT = RPB(j)
            j1 = j
        endif
    enddo
    YD1 = RPB(4)
    YD2 = ZPB(4)
    RPB(4)  = RPB(j1)
    ZPB(4)  = ZPB(j1)
    RPB(j1) = YD1
    ZPB(j1) = YD2

    end subroutine BNDRY

!---------------------------------------------------------------------
    subroutine RHSEQ()

!---------------------------------------------------------------------
! Input: RTOR, BTOR, NA, NA1, HRO, NB2EQL,
!  NE, NI, TE, TI, MU, CU, AMETR, RHO, PBLON, PBPER, G22, G33, IPOL
! Output:
!         EQPF
!         EQFF
!         CUTOR
! Both quantities EQPF (~p') and EQFF (~II') are given in [MA/m^2]
! EQPF = -1.6E-3*(2*\pi*R_0)\prti{n_13*T_keV}{\psi[Vs=T*m^2]}
!        = -1.E-6*/(2*\pi*R_0)\prti{p[J/m^3=Pascal]}{\psi[Vs]}
! EQFF = -1.E-6*2*\pi/(R_0*\mu_0)*I*\prti{I}{\psi}
!        = -5./R_0*I*\prti{I}{\psi}
! Local toroidal current density j[MA/m^2] is EQPF*r/R_0+EQFF*R_0/r, i.e.
!         j(r, z) = r*(\vec j\cdot\nabla\zeta) = EQPF*r/R_0+EQFF*R_0/r ,
! ASTRA average toroidal current density is
!    R_0*<\vec j\cdot\nabla\zeta> = EQPF+EQFF*<R_0^2/r^2>
!---------------------------------------------------------------------

    use scalars, only: RTOR, BTOR, HRO, NA, NA1, NB2EQL
    use status, only: EQPF, EQFF, NE, TE, NI, TI, PBLON, PBPER, PFAST, &
        RHO, AMETR, CU, CUTOR, G22, G33, MU, IPOL
    use debugger, only: markloc, debug

    integer :: j
    double precision :: YCB, YG, YTH2

    call markloc('RHSEQ', debug_lev=3*debug)

! Preparing input for the 3M equilibrium solver:

    YCB = 1.6E-3*RTOR/(BTOR*HRO**2)
    do J=2, NA
        EQFF(J) = ( (NE(J+1)*TE(J+1) - NE(J)*TE(J)) + &
                    (NI(J+1)*TI(J+1) - NI(J)*TI(J)) )/J
        EQFF(J) = EQFF(J) + 0.5*NB2EQL * &
            (PBLON(J+1) - PBLON(J) + PBPER(J+1) - PBPER(J))/J
        EQFF(J) = EQFF(J) + (PFAST(J+1) - PFAST(J))/J
        EQFF(J) = -YCB*EQFF(J)/(MU(J))
    enddo
    EQFF(1) = EQFF(2)
    EQFF(NA1) = EQFF(NA) + (EQFF(NA) - EQFF(NA-1)) * &
        (AMETR(NA) - AMETR(NA-1))/(AMETR(NA1) - AMETR(NA))
    EQFF(NA1) = EQFF(NA) + (EQFF(NA) - EQFF(NA-1))
    do J=1, NA1
        EQPF(j) = EQFF(j)
        YTH2 = RHO(j)*G22(J)*(MU(J)/RTOR)**2
        YG = (1. + YTH2)*G33(J)
        EQFF(J)  = (CU(J)/IPOL(J) - EQPF(j))/YG
        CUTOR(J) = (CU(J)/IPOL(J) + YTH2*EQPF(j))/(1. + YTH2)
    enddo
!    else ! calculates only cutor
!        do J=1, NA1
!            YTH2 = RHO(j)*G22(J)*(MU(J)/RTOR)**2
!            CUTOR(J) = (CU(J)/IPOL(J) + YTH2*EQPF(j))/(1. + YTH2)
!        enddo
!    endif

    end subroutine RHSEQ

!---------------------------------------------------------------------
    subroutine RHSEQ2()

!---------------------------------------------------------------------
! Input: RTOR, BTOR, NA, NA1, HRO, NB2EQL,
!  NE, NI, TE, TI, MU, CU, AMETR, RHO, PBLON, PBPER, G22, G33, IPOL
! Output:
!         EQPF
!         EQFF
!         CUTOR
! Both quantities EQPF (~p') and EQFF (~II') are given in [MA/m^2]
! EQPF = -1.6E-3*(2*\pi*R_0)\prti{n_13*T_keV}{\psi[Vs=T*m^2]}
!        = -1.E-6*/(2*\pi*R_0)\prti{p[J/m^3=Pascal]}{\psi[Vs]}
! EQFF = -1.E-6*2*\pi/(R_0*\mu_0)*I*\prti{I}{\psi}
!        = -5./R_0*I*\prti{I}{\psi}
! Local toroidal current density j[MA/m^2] is EQPF*r/R_0+EQFF*R_0/r, i.e.
!         j(r, z) = r*(\vec j\cdot\nabla\zeta) = EQPF*r/R_0+EQFF*R_0/r ,
! ASTRA average toroidal current density is
!    R_0*<\vec j\cdot\nabla\zeta> = EQPF+EQFF*<R_0^2/r^2>
!---------------------------------------------------------------------

    use scalars, only: RTOR, BTOR, NA, NA1, NB2EQL
    use status, only: EQPF, EQFF, NE, TE, NI, TI, PBLON, PBPER, PFAST, &
        RHO, AMETR, CU, CUTOR, G22, MU, IPOL, FP
    use debugger, only: markloc, debug

    integer :: j
    double precision :: YTH2, press
    double precision :: z1

    call markloc('RHSEQ2', debug_lev=3*debug)

! Preparing input for the 3M equilibrium solver:
    do J=2, NA
        press = ( (NE(J+1)*TE(J+1) - NE(J)*TE(J)) + (NI(J+1)*TI(J+1) - NI(J)*TI(J)) )
        press = press + 0.5*NB2EQL * (PBLON(J+1) - PBLON(J) + PBPER(J+1) - PBPER(J))
        press = press + (PFAST(J+1) - PFAST(J))
        EQPF(J) = 1602.*press/(FP(J+1)-FP(J))
    enddo
    EQPF(1) = EQPF(2)
    EQPF(NA1) = EQPF(NA) + (EQPF(NA) - EQPF(NA-1)) * (AMETR(NA) - AMETR(NA-1))/(AMETR(NA1) - AMETR(NA))
    EQPF(NA1) = EQPF(NA) + (EQPF(NA) - EQPF(NA-1))
    do J=2, NA
        press = 0.5*(IPOL(J+1)**2 - IPOL(J)**2)/(FP(J+1)-FP(J))
        EQFF(J) = press * (RTOR*BTOR)**2
    enddo
    EQFF(1) = EQFF(2)
    EQFF(NA1) = EQFF(NA) + (EQFF(NA) - EQFF(NA-1)) * (AMETR(NA) - AMETR(NA-1))/(AMETR(NA1) - AMETR(NA))
    EQFF(NA1) = EQFF(NA) + (EQFF(NA) - EQFF(NA-1))
    do J=1, NA1
        z1 = 1.e-6/(GP2*RTOR)*EQPF(j)
        YTH2 = RHO(j)*G22(J)*(MU(J)/RTOR)**2
        CUTOR(J) = (CU(J)/IPOL(J) + YTH2*z1)/(1. + YTH2)
    enddo
!    else ! calculates only cutor
!        do J=1, NA1
!            YTH2 = RHO(j)*G22(J)*(MU(J)/RTOR)**2
!            CUTOR(J) = (CU(J)/IPOL(J) + YTH2*EQPF(j))/(1. + YTH2)
!        enddo
!    endif

    end subroutine RHSEQ2

!---------------------------------------------------------------------
    subroutine CUOFMU()
!---------------------------------------------------------------------
! Compute CU(rho) and FP(rho) from MU(rho)
!---------------------------------------------------------------------
! Input: RHO - radial grid (m)
!          NA1 - number of grid points
!          GP2 - 2\pi
!          RTOR - R_0 [m]
!          BTOR - B_0 [T]
!          G22(1:NA1-1) - <g22/g>*............
!          G33(1:NA1-1) -
!          IPOL(1:NA1-1)-
!          CD(1:NA1-1)  - external (+bootstrap) current
!          MU(1:NA1)     - (1/rho)dF/d(rho) rotational transform
! Output: CU(1:NA1) - (1/rho)d{K*dF/d(rho)}/d(rho) current density
!            FP(1:NA1) - poloidal flux [Vs]
!---------------------------------------------------------------------

    use scalars, only: RTOR, BTOR, NA1, NA, HRO, IPEQL
    use status, only: SRHO, XRHO, CU, MU, FP, G22, G33, IPOL, &
    SG11, SG12, SG21, SG22, VR, RHO
    use numerical_tools, only: extrap, integr, deriv, qinterp

    integer :: j
    double precision :: YAJ, YCJ
    double precision, dimension(NA1) :: YAR, dumI, dumF, YAR1

    do j=1, NA1
        YAR(j) = GP2*BTOR*MU(j)*SRHO(j)
    enddo

    call INTEGR(SRHO(1:NA1), 2, YAR(1:NA1), FP(1:NA1), NA1)

    YAJ = 0.
    do J=1, NA
        YCJ = YAJ
        YAJ = (FP(j+1) - FP(j))/HRO**2
        YAJ = G22(j)*YAJ
        CU(j) = (YAJ - YCJ)/(HRO*(j - 0.5))
    enddo

    CU(NA1) = EXTRAP(XRHO(1:NA), CU(1:NA), XRHO(NA1), NA, 2, .false.)
    YCJ = 1.25/(GP**2 * RTOR)
    do J=1, NA1
        CU(j) = YCJ*CU(j)*G33(J)*IPOL(J)**3
    enddo

! Computation of CU for stellarators (FP is equivalent for tokamaks)
    if (IPEQL == 9 .or. IPEQL == 6 .or. IPEQL == 7) then
        call DERIV(RHO(1:NA1), SRHO(NA1), 1, FP(1:NA1), YAR(1:NA1), 1, NA1, 1)
! Calculate plasma current and poloidal current via SGij (defined on shifted grid)
        dumI(1:NA1) = (SG11(1:NA1)*YAR(1:NA1) + GP2*BTOR*SRHO(1:NA1)*SG12(1:NA1))/(0.4*GP)
        dumF(1:NA1) = (SG21(1:NA1)*YAR(1:NA1) + GP2*BTOR*SG22(1:NA1))/(0.4*GP)
! Calculate deriv(I/F) and F**2 on main grid
        call DERIV(SRHO(1:NA1), RHO(NA1), 2, dumI(1:NA1)/dumF(1:NA1), YAR(1:NA1), 1, NA1, 0)
        call qinterp(RHO(1:NA1), dumF(1:NA1)**2, NA1, srho(1:NA1), YAR1(1:NA1), NA1)
! Calculate CU
        CU(1:NA1) = 0.4*GP*YAR1(1:NA1)*YAR(1:NA1)/(VR(1:NA1)*BTOR)
        CU(NA1) = EXTRAP(XRHO(1:NA) , CU(1:NA ), XRHO(NA1), NA, 2, .false.)
        CU(1)   = EXTRAP(XRHO(2:NA1), CU(2:NA1), XRHO(1  ), NA, 2, .true. )
    endif

    end subroutine CUOFMU

!---------------------------------------------------------------------
    subroutine CUOFP()
!---------------------------------------------------------------------
! Compute CU(rho) and MU(rho) from FP(rho)
!---------------------------------------------------------------------
! Input: HRO - radial step (m)
!  NA1 - number of grid points
!  G22(1:NA) - <g22/g>*............
!  G33(1:NA) -
!  IPOL(1:NA) -
!  CD(1:NA) - external (+bootstrap) current
!  FP(1:NA1) - poloidal flux
! Output: CU(1:NA1) - (1/rho)d{K*dF/d(rho)}/d(rho) current density
!  MU(1:NA1) - (1/rho)dF/d(rho)      rotational transform
!---------------------------------------------------------------------

    use status, only: XRHO, FP, MU, CU, IPOL, G22, G33, SXHO, &
        SG11, SG12, SG21, SG22, VR, SXHO, MV, FV, RHO, SRHO
    use scalars, only: RTOR, HRO, BTOR, NA, NA1, IPEQL
    use numerical_tools, only: extrap, deriv, qinterp

    integer :: j
    double precision :: YAJ, YCJ
    double precision, dimension(NA1) :: YAR, YAR1, dumI, dumF

    YAJ = 0.
    CU  = 0.
    MU  = 0.

    do J=1, NA
        YCJ = YAJ
        YAJ = (FP(j+1) - FP(j))/HRO**2
        MU(j) = YAJ/j
        YAJ = G22(j)*YAJ
        CU(j) = (YAJ - YCJ)/HRO
        CU(j) = CU(j)/(j - 0.5)
    enddo
    CU(NA1) = EXTRAP(XRHO(1:NA), CU(1:NA), XRHO(NA1), NA, 2, .false.)
    MU(NA1) = EXTRAP(SXHO(1:NA), MU(1:NA), SXHO(NA1), NA, 2, .false.)
    YCJ = 1.25/(GP**2 * RTOR)
    YAJ = 0.5/(GP*BTOR)
    do J=1, NA1
        CU(j) = YCJ*CU(j)*G33(J)*IPOL(J)**3
        MU(J) = YAJ*MU(j)
    enddo

! Computation of CU for stellarators (MU is equivalent to tokamak)
    if (IPEQL == 9 .or. IPEQL == 6 .or. IPEQL == 7) then
        call DERIV(RHO(1:NA1), SRHO(NA1), 1, FP(1:NA1), YAR(1:NA1), 1, NA1, 1)		
! Calculate plasma current and poloidal current via SGij (defined on shifted grid)
        dumI(1:NA1) = (SG11(1:NA1)*YAR(1:NA1) + GP2*BTOR*SRHO(1:NA1)*SG12(1:NA1))/(0.4*GP)
! Calculate deriv(I/F) and F**2 on main grid
        call DERIV(SRHO(1:NA1), RHO(NA1), 2, dumI(1:NA1)/IPOL(1:NA1), YAR(1:NA1), 1, NA1, 0)
        call qinterp(RHO(1:NA1), IPOL(1:NA1)**2, NA1, SRHO(1:NA1), YAR1(1:NA1), NA1)
! Calculate CU
        CU(1:NA1) = GP2*RTOR*YAR1(1:NA1)*YAR(1:NA1)/VR(1:NA1)
        CU(NA1) = EXTRAP(XRHO(1:NA ), CU(1:NA ), XRHO(NA1), NA, 2, .false.)
        CU(1)   = EXTRAP(XRHO(2:NA1), CU(2:NA1), XRHO(  1), NA, 2, .true.)
    endif

    end subroutine CUOFP

!---------------------------------------------------------------------
    subroutine FPMUOFCU()
!---------------------------------------------------------------------
! Compute FP(rho) and MU(rho) from CU(rho)

    use scalars, only:  RTOR, BTOR, NA1, NA, HRO, IPEQL, IPL
    use status, only: RHO, SRHO, XRHO, CU, MU, FP, G22, G33, IPOL, &
        SG11, SG12, SG21, SG22, VR, CUBS, CD, FV, MV, SXHO
    use numerical_tools, only: extrap, integr, deriv, derivcc, integrcc

    integer :: j
    double precision :: YAJ, YCJ
    double precision, dimension(NA1) :: YAR, YAR1, dumI, dumF, ywa, phi, dumFp
    double precision ym, ymcd, yioh, yicd, yc, ym1, yf, yu, yj_cu, y_ipl

!computation of FP and MU for stellarators from CU (jpar)
    if (IPEQL == 9 .or. IPEQL == 6 .or. IPEQL == 7) then
        phi = GP*BTOR*SRHO(1:NA1)**2
        yar(1:NA1) = CU(1:NA1)*VR(1:NA1)/IPOL(1:NA1)**2/(GP2*RTOR)  ! F I' - F' I
        call integr(rho(1:NA1), 1, yar(1:NA1), dumF(1:NA1), NA1)
        dumI(1:NA1) = dumF(1:NA1)*IPOL(1:NA1)
        yc = IPL/dumI(NA1) ! the user needs to not input a cu that integrates to zero current. IPL can be zero
        CU = CU * yc
        dumI = dumI * yc

! Calculate full mu
        MU(1:NA1) = (0.4*GP/(GP2*BTOR*SRHO(1:NA1))*dumI(1:NA1) - SG12(1:NA1))/SG11(1:NA1)

        call integr(srho(1:NA1), 2, GP2*BTOR*srho(1:NA1)*mv(1:NA1), FV(1:NA1), NA1)
        call integr(srho(1:NA1), 2, GP2*BTOR*srho(1:NA1)*mu(1:NA1), FP(1:NA1), NA1)
        return
    endif

    YF = GP2*HRO**2 * BTOR
    YC = 0.4*GP*RTOR/BTOR
    YM = 0.
    YMCD = 0.
    do j=1, NA1
        if (j == NA1) CYCLE
        YWA(J) = CUBS(J)+CD(J)
        YM   = YM   +  CU(J)*RHO(J)/(G33(J)*IPOL(J)**3)
        YMCD = YMCD + YWA(J)*RHO(J)/(G33(J)*IPOL(J)**3)
    enddo
    YIOH = GP2*YM*HRO*IPOL(NA1)
    YICD = GP2*YMCD*HRO*IPOL(NA1)
    CU(1: NA1) = CU(1: NA1)*(IPL - YICD)/YIOH
    YM = 0.
    do j=1, NA
        YM = YM + YC*CU(J)*RHO(J)/(G33(J)*IPOL(J)**3)
        YM1 = YM/G22(J)
        MU(J) = YM1/J
        FP(J+1) = FP(J) + YF*YM1
    enddo
    mu(NA1) = EXTRAP(SXHO(1:NA), MU(1:NA), SXHO(NA1), NA, 2, .false.)
    YU = GP2*RTOR
    YJ_CU = (FP(NA1) - FP(NA))/HRO * IPOL(NA1) * G22(NA)/(0.4*GP*RTOR)
    YJ_CU = IPL/YJ_CU
    do j=1, NA1
        CU(J) = YJ_CU*CU(J)
        MU(J) = YJ_CU*MU(J)
        FP(J) = YJ_CU*(FP(J) - FP(1)) + FP(1)
    enddo

    end subroutine FPMUOFCU

!---------------------------------------------------------------------
    subroutine new_grid()
!---------------------------------------------------------------------
! Input:  XRHO, SXHO, HROX, ROC, NA1, AB, ABC, AMETR(NA1)
! Output: NB1, RHO, SRHO, HRO, AMETR(j>NA1)
!---------------------------------------------------------------------

    use status, only: NRD, RHO, XRHO, SRHO, SXHO, AMETR
    use scalars, only: HRO, HROX, AB, ABC, ROC, ROWALL, &
        FTO, BTOR, NA, NA1, NB1, NAB

    integer :: j
    double precision :: YDA

    HRO = HROX*ROC
    do j=1, NRD
        RHO(j)  = XRHO(j)*ROC
        SRHO(j) = SXHO(j)*ROC
    enddo
    FTO = GP*BTOR*ROC**2

    do j=1, NRD
        if (RHO(j) >= ROWALL) EXIT
    ENDDO
    NB1 = min(NRD, j)

    AMETR(NA1) = ABC

    if ( 2.*abs(AB - ABC) < AMETR(NA1) - AMETR(NA) ) then
        AB = ABC
        return
    endif

    NAB = min(NRD, nint(NA1/ABC*AB))

    if (NA1 + 1 > NB1 .or. NA1 == NAB) return

    YDA = (AB - ABC)/(NAB - NA1)
!the values outside NA1 need to be controlled if NAB is not equal to NA1/ABC*AB
    if (NA1 < NB1) then
        do j=NA1+1, NB1
            AMETR(j) = ABC + YDA*(j - NA1)
            if (AMETR(j) < AB) NAB = j
        enddo
        if (NAB < NB1) NAB = NAB + 1
    endif
! the values at NAB maybe wrong, do not use them!
    NAB = min(NRD, NAB)
    AMETR(NAB) = AB

    end subroutine new_grid

!---------------------------------------------------------------------
    subroutine SETGEO()
!---------------------------------------------------------------------
! input:  NB1, HRO, ROC, AB, RHO(j)
! Output: SHIF(1:NB1), ELON(1:NB1), TRIA(1:NB1),
!  AMETR(1:NB1), DRODA(1:NB1)

    use scalars, only: NB1, AB, ROC, RTOR, SHIFT, UPDWN, ELONG, TRIAN
    use status, only: RHO, SHIF, SHIV, ELON, TRIA, AMETR, DRODA

    integer :: j
    double precision :: YDA, YA, YR1, YR2, rho_n_sq

    YDA = 0.1*AB/NB1
    YR1 = 0.
    YA  = 0.
    YR2 = 0.

    do j=1, NB1
        rho_n_sq = min(1.d0, (RHO(J)/ROC)**2)
        SHIF(J) = SHIFT
        SHIV(J) = UPDWN
        ELON(J) = 0.5*(1. + ELONG + (ELONG - 1.)*rho_n_sq)
        TRIA(J) = TRIAN*rho_n_sq
        do while(YR2 <= RHO(j))
            YR1 = YR2
            YA = YA + YDA
            YR2 = ROC3A(RTOR, SHIF(j), YA, ELON(j), TRIA(j))
        enddo
        DRODA(j) = YDA/(YR2 - YR1)
        AMETR(j) = YA - YDA + (RHO(j) - YR1)*DRODA(j)
    enddo

    end subroutine SETGEO

!---------------------------------------------------------------------
    subroutine fourier_expansion(m_equil, mnmax, brangle, &
        xm, xnn, rmnc_lcfs, rmns_lcfs, zmnc_lcfs, zmns_lcfs, bndr, bndz)
! Construct flux surfaces from Fourier moments

    integer, intent(in) :: m_equil, mnmax
    double precision, intent(in) :: brangle
    integer, intent(in), dimension(*) :: xm, xnn
    real*8, intent(in), dimension(*) :: rmnc_lcfs, rmns_lcfs, zmnc_lcfs, zmns_lcfs
    real*8, intent(out), dimension(m_equil) :: bndr, bndz

    integer :: k1, k
    double precision :: theta, angle

    bndr = 0.
    bndz = 0.
    do k1=1, m_equil
        theta = GP2*(k1 - 1.)/(0. + m_equil)
        do k=1, mnmax
            angle = xm(k)*theta - xnn(k)*brangle
            bndr(k1) = bndr(k1) + rmnc_lcfs(k)*cos(angle) + rmns_lcfs(k)*sin(angle)
            bndz(k1) = bndz(k1) + zmns_lcfs(k)*sin(angle) + zmnc_lcfs(k)*cos(angle)
        enddo
    enddo
    end subroutine fourier_expansion

!---------------------------------------------------------------------
    subroutine VMEC2ASTRA(vmec_vacuum, vmec_tau, vmec_dteq, yes_boozer)
!---------------------------------------------------------------------
! Reads VMEC output and converts it to ASTRA quantities

    use scalars, only: NEQUIL, MEQUIL, NA1, HROX, HRO, ROC, RTOR, ABC, BTOR, &
        VOLUME, GVAC, FTO, UPDWN
    use status
    use stella_module, only: dphidsb_stella, dphidvpb_ip, dphidvpb_f0, &
        stella_which_surf
    use parameters_a2equil, only: equil_now
    use numerical_tools, only: qinterp

    integer, intent(in):: vmec_vacuum, vmec_tau, vmec_dteq, yes_boozer
    double precision :: Rmaj(NA1), theta, angle, f_boundary, nfperiods, brangle
    double precision :: dummo1, dummo2, phi_full_surfaces
    character(len=128) :: cmd, str_NA1, mpi_command
    integer ivmec, k, k1
    integer :: mnmax
    integer, allocatable :: xm(:), xnn(:)
    integer :: ns_temp, k2
    logical :: lasym
    real*8, allocatable :: rmnc_lcfs(:), zmns_lcfs(:)
    real*8, allocatable :: rmnc_all(:, :), zmns_all(:, :), bndr(:), bndz(:)
    real*8, allocatable :: rmnc_interp(:, :), zmns_interp(:, :), xrho_eq(:), phi_temp(:)
    real*8, allocatable :: rmns_lcfs(:), zmnc_lcfs(:)
    real*8, allocatable :: rmns_all(:, :), zmnc_all(:, :)
    real*8, allocatable :: rmns_interp(:, :), zmnc_interp(:, :)
    integer :: nt_bnd, n_phase

    save f_boundary

    if (.not.associated(equil_now%coord_sys%position%r)) then
        allocate(equil_now%coord_sys%position%r(NEQUIL, MEQUIL))
        allocate(equil_now%coord_sys%position%z(NEQUIL, MEQUIL))
    endif
    if (.not.allocated(bndr)) then
        allocate(bndr(MEQUIL))
        allocate(bndz(MEQUIL))
    endif

    call a2vmec(vmec_vacuum, vmec_tau, vmec_dteq, yes_boozer, f_boundary, phi_full_surfaces)

    write(str_NA1, '(I0)') NA1

    call get_environment_variable("PYTHON_BIN", mpi_command)
    cmd = trim(mpi_command)//' python/VMEC2ASTRA.py ' // trim(str_NA1)
    write(*, *) cmd
    call execute_command_line(trim(cmd), wait=.true.)

    open(unit=10, file='dat/VMEC2ASTRA.bin', form='unformatted', access='stream')
    read(10) HROX, HRO, ROC, RTOR, ABC, BTOR, volume, GVAC, f_boundary, nfperiods, &
        dphidsb_stella, dummo1, dummo2
    read(10) RHO(1:NA1), SRHO(1:NA1), SG11(1:NA1), SG12(1:NA1), &
        SG21(1:NA1), SG22(1:NA1), MV(1:NA1), VR(1:NA1), &
        VRS(1:NA1), GRADRO(1:NA1), G11(1:NA1), Rmaj(1:NA1), &
        VOLUM(1:NA1), AMETR(1:NA1), SLAT(1:NA1), FTPT(1:NA1)
    read(10) mnmax

    AREAT = 0.*SLAT  ! to be implemented

    FTO = GP*BTOR*ROC**2 !needed

! Compute F = IPOL*RTOR*BTOR --> IPOL = F/RTOR/BTOR
    IPOL = (SG21*MU  + SG22/SRHO)*SRHO/RTOR   ! From F. Solfronk et al., PPCF 2026

    allocate(xm(mnmax), xnn(mnmax))
    allocate(rmnc_lcfs(mnmax), zmns_lcfs(mnmax))
    allocate(rmns_lcfs(mnmax), zmnc_lcfs(mnmax))
    rmns_lcfs = 0.0
    zmnc_lcfs = 0.0

    read(10) xm
    read(10) xnn
    read(10) rmnc_lcfs
    read(10) zmns_lcfs
    read(10) lasym
    if (lasym) then
        read(10) rmns_lcfs
        read(10) zmnc_lcfs
    endif
    read(10) naxis
    if (allocated(raxiscc)) deallocate(raxiscc)
    if (allocated(zaxiscc)) deallocate(zaxiscc)
    allocate(raxiscc(naxis), zaxiscc(naxis))
    read(10) raxiscc
    read(10) zaxiscc
! --- Read number of flux surfaces ---
    read(10) ns_temp

! --- Allocate full Fourier coefficient arrays ---
    allocate(rmnc_all(mnmax, ns_temp))
    allocate(zmns_all(mnmax, ns_temp))
    allocate(rmns_all(mnmax, ns_temp))
    allocate(zmnc_all(mnmax, ns_temp))
    allocate(rmnc_interp(ns_temp, mnmax))
    allocate(zmns_interp(ns_temp, mnmax))
    allocate(rmns_interp(ns_temp, mnmax))
    allocate(zmnc_interp(ns_temp, mnmax))
    allocate(xrho_eq(ns_temp))
    allocate(phi_temp(ns_temp))

! --- Read full arrays ---
    read(10) rmnc_all
    read(10) zmns_all
    rmns_all = 0.0
    zmnc_all = 0.0
    if (lasym) then
        read(10) rmns_all
        read(10) zmnc_all
    endif
    read(10) phi_temp

    close(10)

    if (phi_full_surfaces < 0.) then ! plot LCFS for various toroidal angles if phi_full_surfaces < 0
        stella_which_surf = 0
        nt_bnd = 1
! phi = 0, pi/4, pi/2, 3*pi/4, pi
        do n_phase=1, 4
            brangle = GP2*dble(n_phase)/5./nfperiods
            call fourier_expansion(MEQUIL, mnmax, brangle, &
                xm, xnn, rmnc_lcfs, rmns_lcfs, zmnc_lcfs, zmns_lcfs, &
                equil_now%coord_sys%position%r(NEQUIL-n_phase, 1:MEQUIL), &
                equil_now%coord_sys%position%z(NEQUIL-n_phase, 1:MEQUIL))
        enddo
    else ! all flux surfaces at 1 toroidal angle given by phi_full_surfaces >=0.
        stella_which_surf=1
        do k1=1, ns_temp
            xrho_eq(k1) = sqrt(phi_temp(k1)/phi_temp(ns_temp))
        enddo

        do k1=1, mnmax
            call qinterp(xrho_eq(1:ns_temp), rmnc_all(k1, 1:ns_temp), ns_temp, xrho_eq(1:ns_temp), rmnc_interp(1:ns_temp, k1), ns_temp)
            call qinterp(xrho_eq(1:ns_temp), zmns_all(k1, 1:ns_temp), ns_temp, xrho_eq(1:ns_temp), zmns_interp(1:ns_temp, k1), ns_temp)
            call qinterp(xrho_eq(1:ns_temp), rmns_all(k1, 1:ns_temp), ns_temp, xrho_eq(1:ns_temp), rmns_interp(1:ns_temp, k1), ns_temp)
            call qinterp(xrho_eq(1:ns_temp), zmnc_all(k1, 1:ns_temp), ns_temp, xrho_eq(1:ns_temp), zmnc_interp(1:ns_temp, k1), ns_temp)
        enddo

        brangle = phi_full_surfaces
        nt_bnd = 1

! Construct flux surfaces, phi = 0
        do k2=1, ns_temp
            call fourier_expansion(MEQUIL, mnmax, brangle, &
                xm, xnn, rmnc_interp(k2, :), rmns_interp(k2, :), zmnc_interp(k2, :), zmns_interp(k2, :), &
                equil_now%coord_sys%position%r(k2, 1:MEQUIL), &
                equil_now%coord_sys%position%z(k2, 1:MEQUIL))
        enddo
    endif

! Calculate additional quantities
    SHIF(1:NA1) = RMAJ(1:NA1) - RTOR
    UPDWN = 0.0

    deallocate(xm, xnn, rmnc_lcfs, zmns_lcfs)
    deallocate(rmns_lcfs, zmnc_lcfs)
    deallocate(rmnc_all)
    deallocate(zmns_all)
    deallocate(rmnc_interp)
    deallocate(zmns_interp)
    deallocate(rmns_all)
    deallocate(zmnc_all)
    deallocate(rmns_interp)
    deallocate(zmnc_interp)
    deallocate(xrho_eq)
    deallocate(phi_temp)
    deallocate(bndr, bndz)

    end subroutine VMEC2ASTRA

!---------------------------------------------------------------------
    subroutine A2VMEC(vmec_vacuum, vmec_tau, vmec_dteq, yes_boozer, f_boundary, phi_full_surfaces)

    use status, only: NRD
    use scalars, only: NA1, RTOR, BTOR, TIME, ROC, SGNIP, SGNBT, ABC, IPL, PHIEDG, &
        SGNBT, IPART, TSTART, NEQUIL
    use status, only: TE, NE, FP, XRHO, ZEF, MU, ELON, SHif , IPOL, &
        AMETR, VOLUM, PEECR, CUECR, AREAT, rho_pol, FP_NORM, &
        PBLON, PBPER, PFAST, TI, NI, CU, SG11, SG12, MV
    use read_input, only: AWD, astra_ext, nml_file
    use stella_module, only: phi_edge_total, &
        dphidsb_stella, dphidvpb_ip, dphidvpb_f0
    use numerical_tools, only: qinterp, integr, derivcc
    use parameters_a2equil, only : equil_now

    integer, intent(in):: vmec_vacuum, vmec_tau, vmec_dteq, yes_boozer
    double precision, intent(out) :: f_boundary, phi_full_surfaces

    logical :: file_exists
    integer :: ios, i, ivmec, vac_phase_stel, count_rate, start_count, end_count, &
         mboz, nboz, init_vmecco
    double precision :: phi_edgehog, t1, t2, dt_boozer, t_boozero
    double precision, dimension(NRD) :: pressure, svmec, curtorprof, dum1

    character(len=32) :: s_curtor, s_phi, n_nodes
    character(len=256) :: path_to_vmec, stellopt_dir, &
        dat_in_file, mpi_command, iomsg
    character(len=1000) :: command_line, raxis_str, zaxis_str, f_true_surf='vmec_io/true_surfaces.txt'
    character(len=120) :: s_nequil, as_nml, PMASS_FILE="dat/vmecp.dat", &
        PIOTA_FILE="dat/vmeci.dat", PCURR_FILE="dat/vmecc.dat"

    data init_vmecco/0/
    data t_boozero/0./
    save init_vmecco, t_boozero

    NAMELIST / vmec / phi_full_surfaces, dt_boozer, vac_phase_stel, mboz, nboz

    CALL getenv('STELLOPT_PATH', stellopt_dir)
    path_to_vmec = TRIM(stellopt_dir) // '/VMEC2000/Release/'

    as_nml = TRIM(awd) // '/' // TRIM(nml_file)
    write(*, *) 'Reading namelist ', TRIM(as_nml)
    open(57, FILE=TRIM(as_nml), delim='apostrophe')
    read(57, nml=vmec, iostat=ios)
    close(57)

    dat_in_file = TRIM(awd) // '/dat/vmecinput.dat'

    call execute_command_line("nproc")
    call get_environment_variable("MPI_COMMAND", mpi_command)
    call get_environment_variable("N_NODES", n_nodes)

    if (vmec_dteq == 1 .and. vmec_vacuum < 2) then  ! run VMEC
        if (vmec_vacuum == 1) phi_edgehog = PHIEDG  ! use user guess for vacuum
        if (vmec_vacuum == 0) call compute_phi_edgehog(phi_edgehog, f_boundary)  !compute it with plasma
! this for now while the compute is fixed
        phi_edge_total = phi_edgehog
        svmec = XRHO**2  ! s vmec is phi/phi_b
        pressure = NE*TE + NI*TI + 0.5*(PBLON + PBPER) + PFAST
        pressure = 1602.*pressure  ! Pascal

        curtorprof = 0.
! calculate actual curtorprof for vmec
        curtorprof = (MU - MV)*SG11/(0.4*GP)*GP2*BTOR*XRHO*ROC

        if (vac_phase_stel == 1) curtorprof = curtorprof/curtorprof(NA1-1)*IPL ! ATTENTION
        curtorprof(NA1) = IPL

!the input namelist to run vmec using curtor and phiedge are:
! NCURR = 1   ! given toroidal current IPL
! IMATCH_PHIEDGE = 1    ! phiedge is enforced to be matched exactly in vmec
! PMASS_TYPE, PCURR_TYPE and files PMASS_FILE, PCURR_FILE. Iota not used as input
! PRES_SCALE = 1.
! SPRES_PED = 1.
! BLOAT = 1.
! lfreeb = .false.
! RAXIS and ZAXIS guessues

! python sub curtor and phiedge in the rigjht place
        if (vmec_vacuum == 0 .or. vmec_vacuum == 1) then
            call execute_command_line('rm -f ' // TRIM(dat_in_file)) ! Clean, to raise errors
            write(s_curtor, '(F)') vac_phase_stel*IPL*1.e6
            write(s_phi   , '(F)') SGNBT*phi_edgehog
            write(s_nequil, '(I0)') NEQUIL
            command_line = "python python/vmecmodin.py " // trim(s_curtor) // &
                " " // trim(s_phi) // " " // trim(s_nequil)
            if (vmec_vacuum == 0) then
                raxis_str = '"'
                zaxis_str = '"'
                do i=1, naxis
                    write(raxis_str(len_trim(raxis_str)+1:), '(F10.6)') raxiscc(i)
                    write(zaxis_str(len_trim(zaxis_str)+1:), '(F10.6)') zaxiscc(i)
                    if (i < naxis) then
                        raxis_str(len_trim(raxis_str)+1: len_trim(raxis_str)+3) = ", "
                        zaxis_str(len_trim(zaxis_str)+1: len_trim(zaxis_str)+3) = ", "
                    endif
                enddo
                raxis_str = trim(raxis_str) // '"'
                zaxis_str = trim(zaxis_str) // '"'
                command_line = TRIM(command_line) // ' -r ' // TRIM(raxis_str) // &
                                                     ' -z ' // TRIM(zaxis_str)
            endif
!            pause 
            call execute_command_line(command_line)
!            pause
        endif

! Write pressure, iota and current files

        open(132, file=TRIM(PMASS_FILE))
        write(132, *) NA1
        do i=1, NA1
            write(132, *) svmec(i), vac_phase_stel*pressure(i)
        enddo
        close(132)

        open(133, file=TRIM(PIOTA_FILE))
        write(133, *) NA1
        do i=1, NA1
            write(133, *) svmec(i), MU(i)
        enddo
        close(133)
 
        open(134, file=TRIM(PCURR_FILE))
        write(134, *) NA1
        do i=1, NA1
            write(134, *) svmec(i), vac_phase_stel*curtorprof(i)*1.e6  ! A/m^2
        enddo
        write(134, *) vac_phase_stel*IPL*1.e6, SGNBT*phi_edgehog
        close(134)

        call system_clock(count_rate=count_rate)
        call system_clock(start_count)

        command_line = trim(mpi_command) // ' ' // trim(n_nodes) // ' ' // &
            trim(path_to_vmec) // '/xvmec2000' // ' ' // trim(dat_in_file)
        if (init_vmecco == 0) then
            init_vmecco = 1
        else
            command_line = trim(command_line) // ' reset=dat/wout_VMECoutput.nc'
        endif
        write(*, *) command_line
        call execute_command_line(trim(command_line))

        call system_clock(end_count)
        write(*, *) 'time spent on vmec : ', real(end_count-start_count, 8)/real(count_rate, 8)

        command_line = 'mv wout_dat.nc dat/wout_VMECoutput.nc'
        call execute_command_line(command_line)
        write(*, *) 'end vmec'
    endif  ! vmec dteq command

    if (yes_boozer == 1) then !run boozer after vmec only if ivmec ==2, vmec run with 1 or 2
        open(32, file='dat/inboozer.in')
        write(32, '(33333I8)') mboz, nboz
        write(32, *) ' VMECoutput '
        write(32, '(33333I8)') [(i, i=1, NEQUIL)]
        close(32)

        command_line = 'cd dat && ' // trim(mpi_command) // ' ' // trim(n_nodes) // ' ' // &
            trim(stellopt_dir) // '/BOOZ_XFORM/Release/xbooz_xform inboozer.in ../' // TRIM(f_true_surf) // ' cd ..'

        if (TIME >= TSTART+t_boozero .or. TIME <= TSTART) then
            call execute_command_line(command_line)
            t_boozero = time - tstart + dt_boozer
        endif
    endif

! ---------------------------------------------------------------------
! extract_boozer_data.py is called EVERY time .
! It auto-detects the Boozer file format from &VMEC_TO_ASTRA_INPUTS
! BOOZER_FILE:
!   * text  Boozer (archive)  -> reads the full-surface archive that made the
!                                reused DKES table (no Boozer run needed);
!   * NetCDF boozmn (generated)-> reads the fresh 7-surface transform just
!                                produced above when yes_boozer==1.
! Either way it writes dat/b00_profile_boozer.txt and dat/minorradiusW7AS.txt
! for the DKES interface.  NA1 is passed so the B00 profile length matches.
! ---------------------------------------------------------------------
    write(s_curtor, '(I0)') NA1
    command_line = 'python python/extract_boozer_data.py ' // trim(s_curtor)
    call execute_command_line(command_line)

    end subroutine a2vmec

!---------------------------------------------------------------------
    subroutine compute_phi_edgehog(phi_edgehog, f_boundary)

    use status, only: NRD
    use scalars, only: btor, roc, NA1, phiedg, time, volume, IPL, rtor
    use status, only: FP, MU, SG11, SG12, SG21, SG22, &
     PFAST, PBLON, PBPER, NE, TE, TI, NI, VOLUM, XRHO, rho, srho, vr, vrs, sxho, MV
    use numerical_tools, only: derivcc, integrcc
    use stella_module, only: dphidsb_stella, dphidvpb_ip, dphidvpb_f0, f_vacuum_rb

    double precision, intent(out) :: phi_edgehog
    double precision phi(nrd), pressure(nrd), dphidv(nrd)
    double precision pprimv(nrd), muhat(nrd), a(nrd), b(nrd), c(nrd), dum(nrd), d(nrd)
    double precision c1(nrd), c2(nrd), c3(nrd), c4(nrd), dphidv1, F_0, dphidv0(nrd)
    double precision f_boundary, errtol, tolerr, v_temp, i1(nrd), f1(nrd)
    double precision H_b, detS, phipvb, f_vacuum, rocco, s22_temp(nrd), mu_temp(nrd)
    integer niter, i_initto
    data i_initto/0/
    save i_initto

! here insert calculation of ROC
! should solve an equation for roc similar as for the tokamak case

! H = Phi_v

! the equation is:   I_v Psi_v + F_v Phi_v = 4*pi2 *mu_0 P_v
! which becomes H*I_v muhat/Phib  + H * F_v = Pprimv
! I = S11 Psi_v + S12 Phi_v = (S11*muhat / Phib + S12) H = i1 H
! F = S21 Psi_v + S22 Phi_v = (S21*muhat/ Phib + S22) H = f1 H
! Need to solve for Phi given Psi and mu and Pprime
!
! H(i1 H)v * muhat/Phib + H (f1 H)v = Pv
! muhat/phib*H2*i1v + muhat/phib*i1/2*(H2)v + f1v H2 + f1/2 (H2)v = Pv

! c1 H2 + c2 (H2)v = Pv

! phiedg is the vacuum phiedge from user
    phi_edgehog = GP*BTOR*ROC**2
    muhat = mu*phi_edgehog
    tolerr = 1.e-12
    errtol = 100.
    v_temp = volume
    s22_temp = SG22/SRHO

    pressure = NE*TE + NI*TI + 0.5*(PBLON+PBPER)+PFAST
    pressure = 1602.*pressure  !Pascal
    call derivcc(NA1, rho(1:NA1), pressure(1:NA1), pprimv(1:NA1), 2)
    pprimv = -1.*4.*GP**2*0.4*GP*1.e-6*pprimv ! mu0 pprime

    if (i_initto == 0) then
        call integrcc(NA1, srho(1:NA1), 1./(SG21(1:NA1)*MV(1:NA1)+s22_temp(1:NA1)), dum(1:NA1))
        f_vacuum_rb = phiedg/(GP2*dum(NA1))
        f_vacuum = f_vacuum_rb
        i_initto = 1
    else
        f_vacuum = f_vacuum_rb
    endif

    dum(NA1) = GP2*f_vacuum/(SG21(NA1)*MV(NA1) + s22_temp(NA1))
    detS = s22_temp(NA1)*SG11(NA1) - SG12(NA1)*SG21(NA1)
    f_boundary = f_vacuum   ! f_plasma is zero by topological argument ???
    call derivcc(NA1, rho(1:NA1), fp(1:NA1), c1(1:NA1), 2)
    call derivcc(NA1, rho(1:NA1), gp*btor*rho(1:NA1)**2, c2(1:NA1), 2)
    niter = 0

    do while (errtol > tolerr)
        niter = niter + 1
        rocco = sqrt(phi_edgehog/GP/BTOR)
        mu_temp = muhat/phi_edgehog

        phi =  phi_edgehog*xrho**2
        i1 = (SG11*mu_temp + SG12)
        f1 = (SG21*mu_temp + s22_temp)
        H_b = GP2*f_boundary/f1(NA1)
        call derivcc(NA1, rocco*sxho(1:NA1), i1(1:NA1), c3(1:NA1), 2)
        call derivcc(NA1, rocco*sxho(1:NA1), f1(1:NA1), c4(1:NA1), 2)
        c1 = (c3*mu_temp + c4)/(0.5*(i1*mu_temp + f1))
        c2 = pprimv/(0.5*(i1*mu_temp + f1))

        call integrcc(NA1, rocco*xrho(1:NA1), c1(1:NA1), c3(1:NA1))
        call integrcc(NA1, rocco*xrho(1:NA1), c2(1:NA1)*exp(c3(1:NA1) - c3(NA1)), c4(1:NA1))
        c4 = H_b**2 + c4 - c4(NA1)

        dphidv = sqrt(c4*exp(-(c3 - c3(NA1))))/vrs

        call integrcc(NA1, sxho(1:NA1)**2*phi_edgehog, 1./dphidv(1:NA1), dum(1:NA1))

        v_temp = 0.9*v_temp + 0.1*dum(NA1-1)
        errtol = abs(v_temp - volume)/volume*100.
        phi_edgehog = phi_edgehog*(1. - 2.*(v_temp - volume)/(v_temp + volume))

        if (niter > 1000000 ) errtol = 0.
    enddo !phi loop

    write(*, *) 'New phiedge : ', phi_edgehog, GP*BTOR*ROC**2, phiedg

    end subroutine compute_phi_edgehog

!---------------------------------------------------------------------
    subroutine EQCYL_STELLA
!---------------------------------------------------------------------
! Quasi-cylindrical assignment: Called if IPEQL==-2
!
! In: RTOR, SHIFT, ABC, ELONG, TRIAN, NA1, NB1
! Out: NA, HRO, ROC, RHO(j), DRODA, VOLUM, IPOL, G33, GRADRO, G11, G22, SLAT
!  BDB02, B0DB2, BDB0, FOFB, BMAXT, BMINT, DRODA, GRADRO
!---------------------------------------------------------------------

    use status, only: RHO, XRHO, VR, VRS, AMETR, SHIF, SHIV, &
        ELON, TRIA, SLAT, G11, G22, G33, G41, G42, G43, G44, G45, &
        BDB0, BDB02, B0DB2, IPOL, MU, &
        FOFB, BMAXT, BMINT, DRODA, GRADRO, VOLUM
    use scalars, only: VOLUME, RTOR, BTOR,  &
        ABC, HRO, ROC, FTO, ROWALL, NA, NA1, NB1, PHIEDG
    use numerical_tools, only: integr
    use debugger, only: markloc, debug

    integer :: j

    call markloc('EQCYL_STELLA', debug_lev=3*debug)

    FTO = PHIEDG
    VOLUME = GP2*RTOR*GP*ABC**2
    ROC = sqrt(FTO/(GP*BTOR))
    HRO = ROC/(NA1 - 0.5)

    do j=1, NB1
        RHO(J) = XRHO(J)*ROC
        if (RHO(j) <= ROWALL) NA = j
        G22(J)   = J*HRO
        VR(J)    = GP2**2 * RTOR*RHO(J)
        VRS(j)   = GP2**2 * RTOR*G22(J)
        AMETR(J) = RHO(J)/ROC*ABC
        SHIF(J)  = 0.
        SHIV(J)  = 0.
        ELON(J)  = 1.
        TRIA(J)  = 0.
        IPOL(J)  = 1.
        G33(J)   = 1.
        G11(J)   = VRS(j)
        SLAT(J)  = VRS(j)
        BDB02(j) = 1. + (RHO(j)*MU(j)/RTOR)**2
        B0DB2(j) = 1./BDB02(j)
        BDB0(j)  = sqrt(BDB02(j))
        FOFB(j)  = 1.
        BMAXT(j) = BTOR*BDB0(j)
        BMINT(j) = BMAXT(j)
        DRODA(j) = 1.
        GRADRO(j)= 1.
        G41(J)   = 1.
        G42(J)   = GRADRO(J)
        G43(J)   = GRADRO(J)
        G44(J)   = G11(J)/VRS(J)
        G45(J)   = G11(J)/VRS(J)
    enddo
    if (NA < NB1) then
        NB1 = NA
    endif
    NA = NA1 - 1
    call INTEGR(RHO, 1, VR, VOLUM, NB1)
    VOLUME = VOLUM(NA1)

    end subroutine EQCYL_STELLA

end module metrics
