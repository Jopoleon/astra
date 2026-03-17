module gs_solver

implicit none
contains

!---------------------------------------------------------------------
    subroutine gssolver( &
! GSSOLVEREF solves for the fix boundary GSE
! Input:
        equil_solver, &
        nr_equ, &  ! radial grid
        ntheta, &   ! poloidal grid
        nbnd, &    ! # boundary points
        jna1, &    ! radial grid for input (dimension of all arrays below)
        rbnd, zbnd, &
        xrho, rtor, btor, &
        roc, fp, pres_in, volume, &
        ncoils, yvcoil, iter_step, iter_part, iter_itreq, &
        inume_3, tau_step, ipsibcf, icircq, ipctrl, ifbey, time_a, &
        psifb_in, psifb, &
        omega_rot, i_rotation, ion_temp, ion_dens, plasma_mass, &
! Output:
        rocnew, ipl, g11, g41, g22, g33, g22e, g33e, eqpf, eqff, &
        vr, vrs, slat, gradro, ipol, bmaxt, bmint, bdb02, bdb0, b0db2, &
        droda, volum, ametr, updwn, shif, elon, tria, fofb, areat, perim, shiv, square)

    use imas_ids, only: type_equilibrium
    use numerical_tools, only: reinterp_back, reinterp_back_quad, qinterp, &
        derivcc, integrcc
    use pi_const, only: GP, GP2, GP2_sq, muvac
    use parameters_a2equil, only: time_fix_eqpff, fix_eqpf_eqff, &
        cheb_degree, spidat_yes, iter_one_only_fbe, advanced_methods, &
        i3method, diagnostic_gsef, do_adcmp, urelax, urelax2, &
        murelax2, ydiff, ydiff2, max_iter, miter_ext, interp_routine, &
        interp_method_rect, epsf_tol, epss_tol, epsv_tol, epsg_tol, &
        key_no_startz, key_no_refits, equil_now
    use io_mod, only: nml_file
    use machine_config, only: json_cfg, config, cfg_exists

    integer, parameter :: nbtabp=1000

    integer, intent(in) :: equil_solver, ntheta, nr_equ, jna1, nbnd, ncoils, &
        iter_step, iter_part, ipsibcf, icircq, ipctrl, &
        iter_itreq, inume_3, ifbey, i_rotation

    double precision, intent(in) :: tau_step, psifb_in, rtor, btor, roc, time_a
    double precision, intent(in), dimension(ncoils) :: yvcoil
    double precision, intent(in), dimension(nbtabp) :: rbnd, zbnd
    double precision, intent(in), dimension(jna1) :: xrho, pres_in, fp, & 
        omega_rot, ion_temp, ion_dens, plasma_mass

    double precision, intent(out) :: rocnew, updwn, psifb
    double precision, intent(out), dimension(jna1) :: ametr, vr, vrs, &
       slat, gradro, shif, tria, elon, ipol, bmaxt, bmint, bdb02, &
       bdb0, b0db2, droda, fofb, areat, perim, volum, &
       g11, g41, g22, g33, shiv, square
    double precision, intent(inout), dimension(jna1) :: eqpf, eqff
    double precision, intent(out), dimension(nr_equ) :: g22e, g33e
    double precision, intent(inout) :: ipl, volume

    logical :: file_existence, found
    integer :: i, j, n_theta, i_call_gsss, k, k1, key_start, keyplc, &
        jiter, p, jveps, jr, jz, nr, nz
    double precision :: dum1r, R0, Z0, Fvacuum, dxrho_sp, dx, &
        phib, PSIb, deltaPSI, PSI0, phibm, phibl, IPLX, Vtemp, Veps, &
        zfuncb, errG, roc_sp, g2ediff, errght, ybound, Rmag, vtemp_counter, &
        Rmin, Rmax, Zmin, Zmax
    double precision, dimension(nbnd) :: Rb, Zb
    double precision, dimension(jna1) :: dpsi_ad, dp_ad, pres, sxho, vxho, xrho_sq, sxho_sq, vxho_sq
    double precision, dimension(nr_equ) :: volum_in, PSI, psi_minus, PRESS, xrho_sp, xrho_sp_sq, &
        GG2, GG3, g11_sp, g41_sp, gradro_sp, xrho_roc_sp, &
        bdb02_sp, ipol_sp, &
        bdb0_sp, b0db2_sp, droda_sp, vr_sp, ametr_sp, &
        volum_sp, tria_sp, &
        eqpf_sp, eqff_sp, dPSI_adcmp, dP_adcmp, &
        PSIn_grid, G2f, pprimp, pprimx, ffprimp, yprimp, &
        psipx, Chat, betahat, lhs, PHI, qqsg, qg3s, As, Bm, AAs, &
        expAA, expAAm, y, H, Hinv, dPSIdV, zfunc, &
        G2m, G3m, G2p, G3p, G2mt, G2pt, dum1, dum2, dum3, &
        Hout, Houtt, hin1, hin2, hout1, hout2, &
        G2tild1, Htild1, G2tild2, Htild2, G2corr2, Hcorr2, & 
        o_rot, i_temp, i_dens, i_mass
    type(type_equilibrium) :: equil_in

!----------------------------------------------------------------------

    data i_call_gsss /0/
    save i_call_gsss

    data vtemp_counter/0.5/
    save vtemp_counter

    namelist / equil_settings / time_fix_eqpff, fix_eqpf_eqff, &
        cheb_degree, spidat_yes, iter_one_only_fbe, advanced_methods, &
        i3method, diagnostic_gsef, do_adcmp, urelax, urelax2, &
        murelax2, ydiff, ydiff2, max_iter, miter_ext, interp_routine, &
        interp_method_rect, epsf_tol, epss_tol, epsv_tol, epsg_tol, &
        key_no_startz, key_no_refits

!-----------------------------

    if (ifbey == 0 .and. i_call_gsss == 0) then
!defaults
        time_fix_eqpff = 0.
        fix_eqpf_eqff = 0
        cheb_degree = 5
        spidat_yes = 0
        iter_one_only_fbe = 0
        advanced_methods = 3
        i3method = 3
        diagnostic_gsef = 0
        do_adcmp = 0
        urelax = 0.5
        urelax2 = 0.1
        murelax2 = 0.2
        ydiff = 1.E+2
        ydiff2 = 1.
        max_iter = 10000
        miter_ext = 1
        interp_routine = 2
        interp_method_rect = 4
        epsf_tol = 1.E-9
        epss_tol = 1.E-5
        epsv_tol = 1.E-5
        epsg_tol = 1.E-8
        key_no_startz = 0
        key_no_refits = 0
        INQUIRE(FILE=TRIM(nml_file), EXIST=file_existence)
        if (file_existence) then
            open(53, FILE=TRIM(nml_file))
            read(53, nml=equil_settings)
            close(53)
        endif
    endif

    n_theta = ntheta
    if (n_theta == 1) n_theta = nbnd
    if (n_theta == 0) n_theta = 1

    dum1r = (nbnd + 1.e-9)/(n_theta + 1.e-9)
    do i=1, n_theta
        j = nint((i-1)*dum1r) + 1
        Rb(i) = rbnd(j)
        Zb(i) = zbnd(j)
    enddo

    pres = 1.E+06*pres_in
    R0   = rtor
    Z0   = 0

    psifb = fp(jna1)
    Fvacuum = rtor*btor
    phib = GP*btor * roc**2

    dx = xrho(2) - xrho(1)
    do j=1, jna1
        sxho(j) = xrho(j) + 0.5*dx
    enddo
    vxho = sxho
    vxho(jna1) = 1.

    dxrho_sp = 1./(nr_equ - 1.)
    do j=1, nr_equ
        xrho_sp(j) = (j - 1.)*dxrho_sp
        volum_in(j) = 0.0
    enddo

!go from astra grid to equilibrium radial grid
    xrho_sq = xrho**2
    sxho_sq = sxho**2
    vxho_sq = vxho**2
    xrho_sp_sq = xrho_sp**2
    call reinterp_back(xrho_sq, fp   , jna1, xrho_sp_sq, PSI     , nr_equ, interp_routine)
    call reinterp_back(xrho_sq, pres , jna1, xrho_sp_sq, PRESS   , nr_equ, interp_routine)
    call reinterp_back(sxho_sq, g22  , jna1, xrho_sp_sq, GG2     , nr_equ, interp_routine)
    call reinterp_back(xrho_sq, g33  , jna1, xrho_sp_sq, GG3     , nr_equ, interp_routine)
    call reinterp_back(vxho_sq, volum, jna1, xrho_sp_sq, volum_in, nr_equ, interp_routine)
    if (fix_eqpf_eqff == -1) then !use eqpf and eqff from metric or user defined
        call reinterp_back(xrho_sq, eqpf, jna1, xrho_sp_sq, eqpf_sp, nr_equ, interp_routine)
        call reinterp_back(xrho_sq, eqff, jna1, xrho_sp_sq, eqff_sp, nr_equ, interp_routine)
    endif

    if (i_rotation == 1) then
        call reinterp_back(xrho_sq, omega_rot  , jna1, xrho_sp_sq, o_rot , nr_equ, interp_routine)
        call reinterp_back(xrho_sq, ion_temp   , jna1, xrho_sp_sq, i_temp, nr_equ, interp_routine)
        call reinterp_back(xrho_sq, ion_dens   , jna1, xrho_sp_sq, i_dens, nr_equ, interp_routine)
        call reinterp_back(xrho_sq, plasma_mass, jna1, xrho_sp_sq, i_mass, nr_equ, interp_routine)
    endif

    GG2(1) = 0.0
    volum_in(1) = 0.0

    if (i_call_gsss /= 0) then
        GG2 = g22e
        GG3 = g33e
    endif

!---------------------------------------
! call EQUIL_CALL_SPID
!---------------------------------------

    IPLX = IPL

    key_start = 0

    G3m = GG3
    G2m = GG2
    G3p = GG3
    G2p = GG2
    G2f = GG2

    phibm = phib
    phibl = phibm

    PSI0 = PSI(1)
    PSIb = PSI(nr_equ)
    deltaPSI = PSIb - PSI0

!here uses as inputs PSI, PRESS, xrho_sp,
!here goes the coupling scheme that computes y (y is F)

! this scheme uses the analytical solution and finite differences for derivatives and integrals. Not very accurate

!compute pprime P', pprimx P_x, psipx dPsi/dx
    call derivcc(nr_equ, PSI    , PRESS, pprimp, 2)
    call derivcc(nr_equ, xrho_sp, PRESS, pprimx, 1)
    call derivcc(nr_equ, xrho_sp, PSI  , psipx , 1)

    zfunc = 0.
    zfuncb = 0.
    dPSI_adcmp = 0.
    dP_adcmp = 0.
    errG = 0.
    g2ediff = 1.

    G2mt = G2m
    G2pt = G2m

    hin1  = 0.
    hin2  = 0.
    hout1 = 0.
    hout2 = 0.
    Houtt = 0.

    allocate(equil_in%profiles_1d%psi(nr_equ))
    allocate(equil_in%profiles_1d%pprime(nr_equ))
    allocate(equil_in%profiles_1d%ffprime(nr_equ))
    allocate(equil_in%profiles_1d%pressure(nr_equ))
    allocate(equil_in%profiles_1d%F_dia(nr_equ))
    allocate(equil_in%eqgeometry%boundary%r(n_theta))
    allocate(equil_in%eqgeometry%boundary%z(n_theta))

    equil_in%eqgeometry%boundary%npoints = n_theta    !one periodic point

    if (i_call_gsss == 0) then
        if (cfg_exists) then
            call config%get('nR', nr, found)
            call config%get('nZ', nz, found)
            call config%get('Rmin', Rmin, found)
            call config%get('Rmax', Rmax, found)
            call config%get('Zmin', Zmin, found)
            call config%get('Zmax', Zmax, found)
        else
            nr = 65
            nz = 65
            Rmin = MINVAL(Rbnd) - 0.05
            Rmax = MAXVAL(Rbnd) + 0.05
            Zmin = MINVAL(Zbnd) - 0.05
            Zmax = MAXVAL(Zbnd) + 0.05
        endif
        allocate(equil_in%eqgeometry%rectgrid%r2d(nr))
        allocate(equil_in%eqgeometry%rectgrid%z2d(nz))
        do jr=1, nr
            equil_in%eqgeometry%rectgrid%r2d(jr) = Rmin + (jr - 1.)*(Rmax - Rmin)/(nr - 1.)
        enddo
        do jz=1, nz
            equil_in%eqgeometry%rectgrid%z2d(jz) = Zmin + (jz - 1.)*(Zmax - Zmin)/(nz - 1.)
        enddo
    endif

!Iteration cycle
    iter_loop: do jiter=1, miter_ext
        p = 1
        Vtemp = volume
        Veps = 1
        phibm = phib
        phibl = phibm
        G2f = G2m

!Iteration cycle

        fsa_gse_loop: do jveps=1, max_iter

            p = p + 1

!start of algorithm
            PHI = phibm*(xrho_sp_sq)
            call derivcc(nr_equ, PSI, PHI, qqsg, 2)
            qg3s = qqsg/G3m

            AAs = G2f/qg3s**2 + (GP2_sq**2)*G3m
            Bm = -GP2_sq*muvac*pprimx/AAs
            dum3 = G2f/qg3s

            call derivcc(nr_equ, xrho_sp, dum3, dum2, 1)
            As = -(dum2/qg3s)/AAs
            ybound = Fvacuum/GP2
            zfuncb = 0.5 * ybound**2

!Solve 1st ODE for zfunc:   z'  - 2*As*z = Bm
! Bm = C/A, As = B/A

            call integrcc(nr_equ, xrho_sp, As, dum2)
            expAA = exp(2.*dum2)  ! e^(A)
            dum2 = expAA/expAA(nr_equ)  ! e^(A)
            expAA = dum2
            dum2 = Bm/expAA
            call integrcc(nr_equ, xrho_sp, dum2, expAAm)  ! int(B*exp(-A))
            expAAm = expAAm - expAAm(nr_equ)
            zfunc = expAA*(zfuncb + expAAm)

            dum2 = G2f/qg3s
            call derivcc(nr_equ, PSI, dum2, dum3, 2)
            betahat = (dum3/qg3s)/AAs
            Chat = -GP2_sq*muvac*pprimp/AAs

            y = sqrt(2.*zfunc)
            H = y/qg3s
            dPSIdV = H
            Hinv = 1./H

! Compute new volume
            call integrcc(nr_equ, PSI, Hinv, dum1)
!end of algorithm

            Vtemp = vtemp_counter*vtemp+(1. - vtemp_counter)*dum1(nr_equ)

            Veps = abs(Vtemp - volume)/volume

! new method
            if (advanced_methods > 0) then
                dum2 = G2f
                if (Veps < epsf_tol) then
                    call qinterp(volum_in, G2m, nr_equ, dum1, G2f, nr_equ)
                    G2f(1) = 0.
                    lhs(1) = sum(abs(G2f(2:nr_equ) - dum2(2:nr_equ))/dum2(2:nr_equ), 1)/nr_equ
                else
                    lhs(1) = 100.
                endif
            else ! old method
                lhs(1) = 0.
                G2f = G2m
            endif

            phibm = phibm*(1. - 2.*(Vtemp - volume)/(Vtemp + volume))

            if (Veps + lhs(1) <= epsf_tol) EXIT

        enddo fsa_gse_loop

        vtemp_counter = 0.

! G2f = G2m   ! Be careful this is now done because numerically I have to find a good way to integrate G2f...

!Update
!Rescale on current
        IPL = 1.E-06*1./(GP2_sq*muvac)*dPSIdV(nr_equ)*G2f(nr_equ)
        H = H/IPL*IPLX
        dPSIdV = H

        y = H*qg3s
        yprimp = Chat - betahat * y**2

        phib = phibm
        roc_sp =  sqrt(phib/(GP*btor))
        PHI = phib * xrho_sp_sq

        call derivcc(nr_equ, PSI, PHI, qqsg, 2)

        dum3 = G2f

! Put correction for error control
        hin2 = hin1
        hin1 = H

        if (advanced_methods > 1 .and. advanced_methods < 4) then
            if (jiter > 1) then
                LHS(4) = exp(min(10., ydiff*errG))
!Alternative method 2
                do j=1, nr_equ
                    LHS(3) = (xrho_sp(j)**4)*exp(-3.*(1. - xrho_sp(j)**4))
                    LHS(2) = LHS(3)*LHS(4)*errG/(1. + LHS(4)*errG)
                    LHS(1) = 1. - LHS(2)
                    dum3(j) = dum3(j)*LHS(1) + LHS(2)*G2f(nr_equ)*H(nr_equ)/H(j)
                enddo
            endif
        endif

        ffprimp = GP2_sq*yprimp
! Skipping g2 preconditioner for now, and for some reasons...

        IPL = 1.E-06*1./(GP2_sq*muvac)*dPSIdV(nr_equ)*G2f(nr_equ)

        PSIn_grid = sqrt((PSI - PSI(1))/(PSI(nr_equ) - PSI(1)))

        if (fix_eqpf_eqff == -1) then !use eqpf and eqff from metric or user defined
            ffprimp(1:nr_equ) = eqff_sp(1:nr_equ)
            pprimp(1:nr_equ)  = eqpf_sp(1:nr_equ)
        else
            eqpf_sp(1:nr_equ) = pprimp(1:nr_equ)
            eqff_sp(1:nr_equ) = ffprimp(1:nr_equ)
        endif
    		
! Compute F according to newfound dPSIdV
        call integrcc(nr_equ, PSI, ffprimp, dum3)
        do j=1, nr_equ
            ipol_sp(j) = (2.**0.5)* ( 0.5*Fvacuum**2 - dum3(nr_equ) + dum3(j) )**0.5
            ipol_sp(j) = ipol_sp(j)/Fvacuum
        enddo

        PSIn_grid = sqrt((PSI - PSI(1))/(PSI(nr_equ) - PSI(1)))

        equil_in%global_param%toroid_field%r0 = R0
        equil_in%global_param%toroid_field%b0 = btor
        equil_in%global_param%i_plasma = IPLX*1e6   !itm is in A
        do i=1, nr_equ
            equil_in%profiles_1d%psi(i)      = PSI(i)
            equil_in%profiles_1d%pprime(i)   = eqpf_sp(i)
            equil_in%profiles_1d%ffprime(i)  = eqff_sp(i)
            equil_in%profiles_1d%pressure(i) = PRESS(i)
            equil_in%profiles_1d%F_dia(i)    = ipol_sp(i)*btor*r0
        enddo

! for current diffusion equation in equil
        do i=1, n_theta
            equil_in%eqgeometry%boundary%r(i) = Rb(i)
            equil_in%eqgeometry%boundary%z(i) = Zb(i)
        enddo

        call a_equil( &
! Inputs
            equil_in, equil_solver, &
            nr_equ, n_theta, iter_step, &
            ncoils, yvcoil, tau_step, time_a, &
            ipsibcf, key_no_refits, &
            icircq, ipctrl, iter_itreq, ifbey, inume_3, &
! Outputs
            key_start, keyplc, equil_now)
        
        !reassign input profiles to output    
        equil_now%profiles_1d%pprime(1:nr_equ)   = eqpf_sp(1:nr_equ)
        equil_now%profiles_1d%ffprime(1:nr_equ)  = eqff_sp(1:nr_equ)
        equil_now%profiles_1d%psi(1:nr_equ)  = PSI(1:nr_equ)

        psifb = psifb_in
        Rmag  = equil_now%coord_sys%position%r(1, 1)
        updwn = equil_now%profiles_1d%shiv(nr_equ)
        Hout(     1:nr_equ) = equil_now%profiles_1d%dPSIdV(1:nr_equ)
        g11_sp(   2:nr_equ) = equil_now%profiles_1d%g1(2:nr_equ)
        G2p(      2:nr_equ) = equil_now%profiles_1d%g2(2:nr_equ)
        G3p(      2:nr_equ) = equil_now%profiles_1d%gm1(2:nr_equ)
        g41_sp(   2:nr_equ) = equil_now%profiles_1d%gm41(2:nr_equ)
        gradro_sp(2:nr_equ) = equil_now%profiles_1d%ggradro(2:nr_equ)
        tria_sp(  1:nr_equ) = 0.5*(equil_now%profiles_1d%tria_upper(1:nr_equ) + &
                                   equil_now%profiles_1d%tria_lower(1:nr_equ))
        volum_sp( 1:nr_equ) = equil_now%profiles_1d%volume(1:nr_equ)
        ametr_sp( 1:nr_equ) = 0.5*(equil_now%profiles_1d%r_outboard(1:nr_equ) - &
                                   equil_now%profiles_1d%r_inboard( 1:nr_equ))

        g11_sp(1)    = 0.
        g41_sp(1)    = Rmag**2
        G2p(1)       = 0.
        G3p(1)       = 1./Rmag**2
        gradro_sp(1) = 0.
        psi_minus = -equil_now%profiles_1d%psi(1:nr_equ)
        
        call derivcc(nr_equ, volum_sp, psi_minus, Hout, 2)
        call derivcc(nr_equ, ametr_sp, volum_sp, droda_sp, 1)

        droda_sp(1)  = 0. 
        ametr_sp(1)  = 0.

        g2ediff = 2.*(G2p(nr_equ) - G2m(nr_equ))/(G2p(nr_equ) + G2m(nr_equ))
        hout2   = 0.
        hout1   = 0.
        G2tild1 = 0.
        Htild1  = 0.
        G2tild2 = 0.
        Htild2  = 0.
        G2corr2 = 0.
        Hcorr2  = 0.
        if (advanced_methods > 2 .and. errG < 1./ydiff2) then
            SELECT CASE(i3method)
            CASE(1)
                errght = 1.
                G2tild1 = G2m
                G2tild2 = G2p
                Htild1 = Houtt
                Htild2 = Hout
                if (jiter > 2) then
                    iter_loop_small: do k=1, miter_ext
! apply correction
                        do k1=1, nr_equ
                            hin1(k1) = min(murelax2, log(Htild2(k1)/Htild1(k1)))
                            hin2(k1) = max(-murelax2, hin1(k1))
                        enddo
                        hin1 = hin2
                        hin2 = G2tild2*hin1
                        G2corr2 = G2tild2 + urelax2*hin2   ! new G2
                        do k1=1, nr_equ
                            hin1(k1) = min(murelax2, log((1.E-3 + G2tild2(k1)) &
                                 /(1.E-3 + G2tild1(k1))))
                            hin2(k1) = max(-murelax2, hin1(k1))
                        enddo
                        hin1 = hin2
                        hin2 = Htild2*hin1
                        Hcorr2 = Htild2 + urelax2*hin2  ! new H
                        hin2(nr_equ)   = 1.E-06*1./(GP2_sq*muvac) *G2corr2(nr_equ)*Hcorr2(nr_equ)
                        hin2(nr_equ-1) = 1.E-06*1./(GP2_sq*muvac) *G2tild2(nr_equ)*Htild2(nr_equ)
                        G2corr2 = G2corr2*hin2(nr_equ-1)/hin2(nr_equ)
                        errght = 2.*abs(G2tild2(nr_equ-5) - G2tild1(nr_equ-5)) &
                                      /(G2tild1(nr_equ-5) + G2tild1(nr_equ-5))
                        G2tild1 = G2tild2
                        G2tild2 = 0.5*(G2tild2 + G2corr2)
                        G2tild2 = 0.5*(G2tild2 + G2corr2)
                        Htild1 = Htild2
                        Htild2 = 0.5*(Htild2 + Hcorr2)
                        if (errght < epsg_tol) EXIT
                    enddo iter_loop_small
                    G2p = G2tild2
                endif
                Houtt = Hout

            CASE(2)
                errght = 0.
                G2tild1 = G2m
                G2tild2 = G2p
                Htild1 = Houtt
                Htild2 = Hout
                if (jiter > 2)  then
                    do k1=1, nr_equ
                        hin1(k1) = min(murelax2, log(Htild2(k1)/Htild1(k1)))
                        hin2(k1) = max(-murelax2, hin1(k1))
                    enddo
                    hin1 = hin2
                    hin2 = G2tild2*hin1
                    G2corr2 = G2tild2 + urelax2*hin2   ! new G2
                    Hcorr2 = Htild2 ! new H
                    hin2(nr_equ)   = 1.E-06*1./(GP2_sq*muvac) *G2corr2(nr_equ)*Hcorr2(nr_equ)
                    hin2(nr_equ-1) = 1.E-06*1./(GP2_sq*muvac) *G2tild2(nr_equ)*Htild2(nr_equ)
                    G2corr2 = G2corr2*hin2(nr_equ-1)/hin2(nr_equ)
                    G2tild1 = G2tild2
                    G2tild2 = 0.5*(G2tild2 + G2corr2)
                    G2tild2 = 0.5*(G2tild2 + G2corr2)
                    Htild1 = Htild2
                    Htild2 = 0.5*(Htild2 + Hcorr2)
                    G2p = G2tild2
                endif
                Houtt = Hout

            CASE(3)
                errght = 0.
                G2tild1 = G2m
                G2tild2 = G2p
                Htild1 = Houtt
                Htild2 = Hout
                if (jiter > 2) then
                    call derivcc(nr_equ, xrho_sp, Htild1, hout1, 1)
                    call derivcc(nr_equ, xrho_sp, Htild2, hout2, 1)
                    hin1 = G2tild2*(hout2/Htild2 - hout1/Htild1)
                    call integrcc(nr_equ, xrho_sp, hin1, hin2)  ! correction to H
                    hin1 = hin2
                    hin2 = hin1
                    G2corr2 = G2tild2 + urelax2*hin2   ! new G2
                    Hcorr2 = Htild2 ! new H
                    G2tild1 = G2tild2
                    Htild1 = Htild2
                    G2p = G2corr2
                endif
                Houtt = Hout
            END SELECT

        endif

        if (iter_one_only_fbe == 1 .and. ifbey > 1) then
            G2p = ((1. - urelax)*G2p + urelax*G2m)
            G3p = ((1. - urelax)*G3p + urelax*G3m)
            errG = 0.0
            if (errG <= epss_tol) then
                EXIT iter_loop
            endif
            G2m = G2p; !on sxrho_sp
            G3m = G3p; !on xrho_sp
        endif

        if (iter_one_only_fbe == -1 .and. ifbey == 1) then
            G2p = ((1. - urelax)*G2p + urelax*G2m)
            G3p = ((1. - urelax)*G3p + urelax*G3m)
            errG = 0.0
            if (errG <= epss_tol) then
                EXIT iter_loop
            endif
            G2m = G2p; !on sxrho_sp
            G3m = G3p; !on xrho_sp
        endif

        if (key_start == 1 .and. key_no_startz == 0) then
            G2p = ((1. - urelax)*G2p + urelax*G2m)
            G3p = ((1. - urelax)*G3p + urelax*G3m)
            errG = 0.0
            if (errG <= epss_tol) then
                EXIT iter_loop
            endif
            G2m = G2p; !on sxrho_sp
            G3m = G3p; !on xrho_sp
        endif

!update G2

        G2p = ((1. - urelax)*G2p + urelax*G2m)
        G3p = ((1. - urelax)*G3p + urelax*G3m)

        errG = 0.0
        do j=1, nr_equ
            errG = max(errG, (abs(G2p(j) - G2m(j))/(abs(G2p(j) + G2m(j)) + 1.E-3)*2.0))
        enddo

        if (errG <= epss_tol) then
            EXIT iter_loop
        endif

        G2m = G2p; !on sxrho_sp
        G3m = G3p; !on xrho_sp
        volume = ((1. - urelax)*volum_sp(nr_equ) + urelax*volum_in(nr_equ))
        volum_sp = volum_sp/volum_sp(nr_equ)*volume

    enddo iter_loop

    deallocate(equil_in%profiles_1d%psi)
    deallocate(equil_in%profiles_1d%pprime)
    deallocate(equil_in%profiles_1d%ffprime)
    deallocate(equil_in%profiles_1d%pressure)
    deallocate(equil_in%profiles_1d%F_dia)
    deallocate(equil_in%eqgeometry%boundary%r)
    deallocate(equil_in%eqgeometry%boundary%z)

    GG2 = G2p
    GG3 = G3p

    PHI = phib * xrho_sp_sq

    call derivcc(nr_equ, PSI, PHI, qqsg, 2)
    xrho_roc_sp = roc_sp*xrho_sp
    call derivcc(nr_equ, xrho_roc_sp, volum_sp, vr_sp, 1)

    IPL = IPLX

! Conversion to ASTRA conventions
    do j=1, nr_equ
        GG3(j)  = GG3(j) *(R0**2)
        g41_sp(j) = g41_sp(j)/(R0**2)
        bdb02_sp(j) = equil_now%profiles_1d%gm4(j)/btor**2
        b0db2_sp(j) = equil_now%profiles_1d%gm5(j)*btor**2
        bdb0_sp(j)  = equil_now%profiles_1d%bdb0(j)/btor
    enddo

    GG2(1) = 0.0
    ametr_sp(1) = 0.0
    droda_sp(1) = 0.0

!------------------------------------------------
! END call EQUIL_CALL_SPID
!------------------------------------------------

    call reinterp_back_quad(xrho_sp, GG2      , nr_equ, sxho, g22    , jna1)
    call reinterp_back_quad(xrho_sp, GG3      , nr_equ, xrho, g33    , jna1)
    call reinterp_back_quad(xrho_sp, vr_sp    , nr_equ, xrho, vr     , jna1)
    call reinterp_back_quad(xrho_sp, vr_sp    , nr_equ, sxho, vrs    , jna1)
    call reinterp_back_quad(xrho_sp, volum_sp , nr_equ, sxho, volum  , jna1)
    call reinterp_back_quad(xrho_sp, g11_sp   , nr_equ, sxho, g11    , jna1)
    call reinterp_back_quad(xrho_sp, g41_sp   , nr_equ, sxho, g41    , jna1)
    call reinterp_back_quad(xrho_sp, droda_sp , nr_equ, sxho, droda  , jna1)
    call reinterp_back_quad(xrho_sp, gradro_sp, nr_equ, sxho, gradro , jna1)
    call reinterp_back_quad(xrho_sp, ipol_sp  , nr_equ, xrho, ipol   , jna1)
    call reinterp_back_quad(xrho_sp, ametr_sp , nr_equ, xrho, ametr  , jna1)
    call reinterp_back_quad(xrho_sp, tria_sp  , nr_equ, xrho, tria   , jna1)
    call reinterp_back_quad(xrho_sp, eqpf_sp  , nr_equ, xrho, eqpf   , jna1)
    call reinterp_back_quad(xrho_sp, eqff_sp  , nr_equ, xrho, eqff   , jna1)
    call reinterp_back_quad(xrho_sp, bdb02_sp , nr_equ, xrho, bdb02  , jna1)
    call reinterp_back_quad(xrho_sp, bdb0_sp  , nr_equ, xrho, bdb0   , jna1)
    call reinterp_back_quad(xrho_sp, b0db2_sp , nr_equ, xrho, b0db2  , jna1)
    call reinterp_back_quad(xrho_sp, dPSI_adcmp,nr_equ, xrho, dpsi_ad, jna1)
    call reinterp_back_quad(xrho_sp, dP_adcmp , nr_equ, xrho, dp_ad  , jna1)
    call reinterp_back_quad(xrho_sp, equil_now%profiles_1d%elongation(1:nr_equ), nr_equ, xrho, elon, jna1)
    call reinterp_back_quad(xrho_sp, equil_now%profiles_1d%surface   (1:nr_equ), nr_equ, xrho, slat, jna1)
    call reinterp_back_quad(xrho_sp, equil_now%profiles_1d%shif (1:nr_equ), nr_equ, xrho, shif , jna1)
    call reinterp_back_quad(xrho_sp, equil_now%profiles_1d%shiv (1:nr_equ), nr_equ, xrho, shiv , jna1)
    call reinterp_back_quad(xrho_sp, equil_now%profiles_1d%squareness(1:nr_equ), nr_equ, xrho, square , jna1)
    call reinterp_back_quad(xrho_sp, equil_now%profiles_1d%fofb (1:nr_equ), nr_equ, xrho, fofb , jna1)
    call reinterp_back_quad(xrho_sp, equil_now%profiles_1d%areat(1:nr_equ), nr_equ, xrho, areat, jna1)
    call reinterp_back_quad(xrho_sp, equil_now%profiles_1d%perim(1:nr_equ), nr_equ, xrho, perim, jna1)
    call reinterp_back_quad(xrho_sp, equil_now%profiles_1d%bmaxt(1:nr_equ), nr_equ, xrho, bmaxt, jna1)
    call reinterp_back_quad(xrho_sp, equil_now%profiles_1d%bmint(1:nr_equ), nr_equ, xrho, bmint, jna1)

    volum(jna1) = volum_sp(nr_equ)
    rocnew = sqrt(phib/(GP*btor))

    g22e = GG2
    g33e = GG3/R0**2

    end subroutine gssolver

!---------------------------------------------------------------------
    subroutine A_EQUIL( &
! Input
        equil_in, equil_solver, &
        nr_equ, n_theta, iter_step, &
        ncoils, yvcoils, tau_step, time_a, &
        ipsibcf, key_no_refits, &
        icircq, ipctrl, &
        iter_itreq, ifbey, inume_3, &
! Output:
        key_start, keyplc, equil_out)

    use imas_ids, only: type_equilibrium
    use spider_params, only: type_parameters
    use pi_const, only: GP, GP2
    use parameters_a2equil, only: fix_adapgrid, s_fazt
    use scalars, only : rtor,shift, updwn
    use io_mod, only: nml_file, machine

    integer, intent(in) :: equil_solver, nr_equ, n_theta, iter_step, ncoils, &
        ipsibcf, key_no_refits, icircq, ipctrl, iter_itreq, ifbey, inume_3
    double precision, intent(in) :: tau_step, time_a
    double precision, intent(in), dimension(ncoils) :: yvcoils
    type(type_equilibrium), intent(in) :: equil_in

    integer, intent(out) :: key_start, keyplc
    type(type_equilibrium), intent(out) :: equil_out

    logical :: file_existence
    integer :: nstep, i, key_equil, nrp, &
        toric_fourc, toric_file, strahl_file, strahl_fourc, &
        write_coils_diagn, key_plcs, kprs, k_grids, &
        kprs2, fixadapgrid

    double precision :: epsros, enelss, k_filessss, ipl, psplexavg, psplexavgexp
    double precision, dimension(ncoils) :: t_currents, ucoils
    character(len=120) :: fname

    type(type_parameters) :: parameters_equil

    save toric_fourc, toric_file
    save strahl_file, strahl_fourc, write_coils_diagn
    save kprs, k_grids, epsros, enelss, key_plcs, k_filessss
    save psplexavg, kprs2, psplexavgexp

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
    parameters_equil%ntheta = n_theta + 2

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
        INQUIRE( FILE=TRIM(nml_file), EXIST=file_existence)
        if (file_existence) then
            open(53, FILE=TRIM(nml_file))
            read(53, nml=spider)
            close(53)
            kprs2 = kprs
        endif
!    fix_adapgrid = fixadapgrid
    endif

    if (s_fazt == 0) then
        kprs = kprs2
    else
        kprs = -2
    endif
    k_grids = 0

    parameters_equil%kpr     = kprs
    parameters_equil%k_grid  = k_grids
    parameters_equil%epsro   = epsros
    parameters_equil%enels   = enelss
    parameters_equil%key_plc = key_plcs
    parameters_equil%key_dmf = 0

    parameters_equil%key_plc = 1    !force plc = 1 if current diffusion is solved
    keyplc = parameters_equil%key_plc
    parameters_equil%key_out = 0

    if (iter_step == 1) parameters_equil%k_fixfree = 0 !astra initialization, no fbe
    if (parameters_equil%k_fixfree == 0) nstep = 0  !no fbe, nstep=0
    parameters_equil%no_circuit_eq = 0
    if (icircq == 0) parameters_equil%no_circuit_eq = 1    !no circuit equations, only static fbe

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
        ucoils(1:ncoils)  = yvcoils(1:ncoils)
        parameters_equil%nstep = nstep
    endif

    if (equil_solver == 101) then
        call feqis_main(ncoils, ucoils, parameters_equil%neql, parameters_equil%k_fixfree, parameters_equil%no_circuit_eq, 1, machine, equil_in, equil_out)
    else
        call spider_run(ncoils, ucoils, equil_in, equil_out, parameters_equil)
    endif
    if (ipsibcf /= 0) parameters_equil%key_psibcf = 1

!output from equil_out structure

    end subroutine A_equil

end module GS_SOLVER
