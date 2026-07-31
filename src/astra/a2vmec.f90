module a2vmec

use pi_const, only: GP, GP2
use read_input, only: AWD

implicit none

integer :: naxis
real*8, allocatable :: raxiscc(:), zaxiscc(:)
character(len=256) :: f_vmec_wout, f_boozer, f_vmec_in, f_vmec_metric, &
    f_boozer_in, vmec_work_dir

contains

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
    subroutine vmec_io_files

    integer :: ios
    character(len=256) :: f_stella_nml
    character(len=256) :: VMEC_WOUT_FILE, BOOZER_FILE, VMEC_IN_FILE, &
        VMEC2A_METRIC, BOOZER_INFILE, VMEC_WD

    f_stella_nml = TRIM(awd) // '/vmec_io/stell_files.nml'
    NAMELIST / vmec_to_astra_inputs / VMEC_WOUT_FILE, BOOZER_FILE, VMEC_IN_FILE, &
        VMEC2A_METRIC, BOOZER_INFILE, VMEC_WD

    write(*, *) 'Reading namelist ', TRIM(f_stella_nml)
    open(58, FILE=TRIM(f_stella_nml), delim='apostrophe')
    read(58, nml=vmec_to_astra_inputs, iostat=ios)
    close(58)

    vmec_work_dir = TRIM(awd) // '/' // TRIM(VMEC_WD) // '/'
    f_vmec_wout   = TRIM(vmec_work_dir) // TRIM(VMEC_WOUT_FILE)
    f_boozer      = TRIM(vmec_work_dir) // TRIM(BOOZER_FILE)
    f_vmec_in     = TRIM(vmec_work_dir) // TRIM(VMEC_IN_FILE)
    f_vmec_metric = TRIM(vmec_work_dir) // TRIM(VMEC2A_METRIC)
    f_boozer_in   = TRIM(vmec_work_dir) // TRIM(BOOZER_INFILE)

    end subroutine vmec_io_files

!---------------------------------------------------------------------
    subroutine vmec_interface(vmec_vacuum, vmec_dteq, yes_boozer)
!---------------------------------------------------------------------
! Reads VMEC output and converts it to ASTRA quantities

    use scalars, only: NA1

    integer, intent(in):: vmec_vacuum, vmec_dteq, yes_boozer

    double precision :: f_boundary, phi_full_surfaces
    character(len=128) :: cmd, py_exe

    save f_boundary

! Set I/O paths for VMEDC I/O files

    call vmec_io_files()

! Run VMEC stand-alone
    call astra2vmec(vmec_vacuum, vmec_dteq, yes_boozer, f_boundary, phi_full_surfaces)

! Collect VMEC output and store it into ASTRA arrays
    call get_environment_variable("PYTHON_BIN", py_exe)
    write(cmd, '(A, A, I0)') trim(py_exe), ' python/vmec2bin.py ', NA1
    write(*, *) cmd
    call execute_command_line(trim(cmd), wait=.true.)

    call vmecBin2astra(phi_full_surfaces)

    end subroutine vmec_interface

!---------------------------------------------------------------------
    subroutine vmecBin2astra(phi_full_surfaces)

    use scalars, only: NEQUIL, MEQUIL, NA1, HROX, HRO, ROC, RTOR, ABC, BTOR, &
        VOLUME, GVAC, FTO, UPDWN
    use status, only: RHO, SRHO, SG11, SG12, SG21, SG22, MV, VR, VRS, GRADRO, &
        G11, VOLUM, AMETR, SLAT, FTPT, AREAT, IPOL, MU, SHIF
    use parameters_a2equil, only: equil_now
    use stella_module, only: dphidsb_stella, stella_which_surf
    use numerical_tools, only: qinterp

    double precision, intent(in) :: phi_full_surfaces
    logical :: lasym
    integer :: k, k1, k2, mnmax, ns_temp, nt_bnd, n_phase
    integer, allocatable, dimension(:) :: xm, xnn
    real*8, allocatable, dimension(:) :: rmnc_lcfs, rmns_lcfs, zmnc_lcfs, zmns_lcfs, &
        xrho_eq, phi_temp
    real*8, allocatable, dimension(:, :) :: rmnc_all, rmns_all, zmnc_all, zmns_all, &
        rmnc_interp, rmns_interp, zmnc_interp, zmns_interp
    double precision :: f_boundary, nfperiods, brangle, dum1, dum2
    double precision :: Rmaj(NA1)

    save f_boundary

    open(unit=10, file=TRIM(f_vmec_metric), form='unformatted', access='stream')
    read(10) HROX, HRO, ROC, RTOR, ABC, BTOR, volume, GVAC, f_boundary, nfperiods, &
        dphidsb_stella, dum1, dum2
    read(10) RHO(1:NA1), SRHO(1:NA1), SG11(1:NA1), SG12(1:NA1), &
        SG21(1:NA1), SG22(1:NA1), MV(1:NA1), VR(1:NA1), &
        VRS(1:NA1), GRADRO(1:NA1), G11(1:NA1), Rmaj(1:NA1), &
        VOLUM(1:NA1), AMETR(1:NA1), SLAT(1:NA1), FTPT(1:NA1)
    read(10) mnmax

    AREAT = 0.*SLAT  ! to be implemented
    FTO = GP*BTOR*ROC**2 ! needed

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

! Clean working files to make sure errors are raised
    call execute_command_line('rm ' // TRIM(vmec_work_dir) // '/*.dat')
    call execute_command_line('rm ' // TRIM(vmec_work_dir) // '/*.txt')

    if (.not.associated(equil_now%coord_sys%position%r)) then
        allocate(equil_now%coord_sys%position%r(NEQUIL, MEQUIL))
        allocate(equil_now%coord_sys%position%z(NEQUIL, MEQUIL))
    endif

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
        stella_which_surf = 1
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

    end subroutine vmecBin2astra

!---------------------------------------------------------------------
    subroutine ASTRA2VMEC(vmec_vacuum, vmec_dteq, yes_boozer, f_boundary, phi_full_surfaces)

    use status, only: NRD
    use scalars, only: NA1, RTOR, BTOR, TIME, ROC, SGNIP, SGNBT, ABC, IPL, PHIEDG, &
        SGNBT, IPART, TSTART, NEQUIL
    use status, only: TE, NE, FP, XRHO, ZEF, MU, ELON, SHif , IPOL, &
        AMETR, VOLUM, PEECR, CUECR, AREAT, rho_pol, FP_NORM, &
        PBLON, PBPER, PFAST, TI, NI, CU, SG11, SG12, MV
    use read_input, only: nml_file
    use stella_module, only: phi_edge_total

    integer, intent(in):: vmec_vacuum, vmec_dteq, yes_boozer
    double precision, intent(out) :: f_boundary, phi_full_surfaces

    integer :: ios, i, vac_phase_stel, count_rate, start_count, end_count, &
         mboz, nboz, init_vmecco
    double precision :: phi_edgehog, dt_boozer, t_boozero
    double precision, dimension(NRD) :: pressure, svmec, curtorprof

    character(len=32) :: n_nodes
    character(len=256) :: stellopt_dir, mpi_command, f_true_surf
    character(len=1000) :: command_line, raxis_str, zaxis_str
    character(len=120) :: as_nml, PMASS_FILE="vmecp.dat", &
        PIOTA_FILE="vmeci.dat", PCURR_FILE="vmecc.dat"

    data init_vmecco/0/
    data t_boozero/0./
    save init_vmecco, t_boozero

    NAMELIST / vmec / phi_full_surfaces, dt_boozer, vac_phase_stel, mboz, nboz

    CALL getenv('STELLOPT_PATH', stellopt_dir)
    f_true_surf = TRIM(awd) // '/vmec_io/true_surfaces.txt'
    as_nml = TRIM(awd) // '/' // TRIM(nml_file)

    write(*, *) 'Reading namelist ', TRIM(as_nml)
    open(57, FILE=TRIM(as_nml), delim='apostrophe')
    read(57, nml=vmec, iostat=ios)
    close(57)

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

!--------------------------------------------------------------
! Replace CURTOR, PHIEDGE and NEQUIL -> write VMEC_WD/vmecinput.dat

        if (vmec_vacuum == 0 .or. vmec_vacuum == 1) then
            call execute_command_line('rm -f ' // TRIM(f_vmec_in)) ! Clean, to raise errors
            write(command_line, '(A, F, 1X, F, 1X, I0)') "python python/vmecmodin.py ", &
                vac_phase_stel*IPL*1.e6, SGNBT*phi_edgehog, NEQUIL
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
            call execute_command_line(TRIM(command_line))
        endif

!-----------------------------------------------------------------
! Write pressure, iota and current files (vmecp, vmeci, vmecc.dat)

        open(132, file=TRIM(vmec_work_dir) // TRIM(PMASS_FILE))
        write(132, *) NA1
        do i=1, NA1
            write(132, *) svmec(i), vac_phase_stel*pressure(i)
        enddo
        close(132)

        open(133, file=TRIM(vmec_work_dir) // TRIM(PIOTA_FILE))
        write(133, *) NA1
        do i=1, NA1
            write(133, *) svmec(i), MU(i)
        enddo
        close(133)
 
        open(134, file=TRIM(vmec_work_dir) // TRIM(PCURR_FILE))
        write(134, *) NA1
        do i=1, NA1
            write(134, *) svmec(i), vac_phase_stel*curtorprof(i)*1.e6  ! A/m^2
        enddo
        write(134, *) vac_phase_stel*IPL*1.e6, SGNBT*phi_edgehog
        close(134)

!---------------
! VMEC execution
!---------------
        call system_clock(count_rate=count_rate)
        call system_clock(start_count)

        command_line = trim(mpi_command) // ' -n ' // trim(n_nodes) // ' ' // &
              TRIM(stellopt_dir) // '/VMEC2000/Release/xvmec2000 ' // trim(f_vmec_in)
        if (init_vmecco == 0) then
            init_vmecco = 1
        else
            command_line = trim(command_line) // ' reset=' // TRIM(f_vmec_wout)
        endif
        write(*, *) TRIM(command_line)
        call execute_command_line(trim(command_line))

! Output files: parvmecinfo.txt, jxbout_dat.nc, mercier.dat, wout_dat.nc, 
!     threed1.dat, timings.txt

        call system_clock(end_count)
        write(*, *) 'time spent on vmec : ', real(end_count-start_count, 8)/real(count_rate, 8)

        command_line = 'mv wout_dat.nc ' // TRIM(f_vmec_wout)
        write(*, *) TRIM(command_line)
        call execute_command_line(TRIM(command_line))
    endif  ! vmec dteq command

!------------------------------
! Dump boozer geometry for DKES

    if (yes_boozer == 1) then
        open(32, file=TRIM(f_boozer_in))
        write(32, '(33333I8)') mboz, nboz
        write(32, *) ' VMECoutput '
        write(32, '(33333I8)') [(i, i=1, NEQUIL)]
        close(32)

        if (TIME >= TSTART+t_boozero .or. TIME <= TSTART) then
            command_line = 'cd ' // TRIM(vmec_work_dir) // ' && ' // &
                TRIM(mpi_command)  // ' -n 1 ' // &
                TRIM(stellopt_dir) // '/BOOZ_XFORM/Release/xbooz_xform ' // &
                TRIM(f_boozer_in)  // ' ' // TRIM(f_true_surf)
            write(*, *) TRIM(command_line)
            call execute_command_line(TRIM(command_line)) ! Input: VMEC_WD/inboozer.in, vmec_io/true_surfaces.txt, ./wout_VMECoutput.nc; Output: ./boozmn_VMECoutput.nc
            t_boozero = time - tstart + dt_boozer
        endif
    endif

! ---------------------------------------------------------------------
! extract_boozer_data.py is called EVERY time .
! It auto-detects the Boozer file-format from vmec_io/stell_files.nml, &VMEC_TO_ASTRA_INPUTS
! BOOZER_FILE:
!   * text  Boozer (archive)  -> reads the full-surface archive that made the
!                                reused DKES table (no Boozer run needed);
!   * NetCDF boozmn (generated)-> reads the fresh 7-surface transform just
!                                produced above when yes_boozer==1.
! Either way it writes VMEC_WD/b00_profile_boozer.txt and VMEC_WD/minorradiusW7AS.txt
! for the DKES interface.  NA1 is passed so the B00 profile length matches.
! ---------------------------------------------------------------------
    write(command_line, '(A, I0)') 'python python/extract_boozer_data.py ', NA1
    write(*, *) TRIM(command_line)
    call execute_command_line(command_line) ! Input: VMEC_WD/boozmn_VMECoutput.nc; Output VMEC_WD/b00_profile_boozer.txt, VMEC_WD/minorradiusW7AS.txt

    end subroutine astra2vmec

!---------------------------------------------------------------------
    subroutine compute_phi_edgehog(phi_edgehog, f_boundary)

    use status, only: NRD
    use scalars, only: btor, roc, NA1, phiedg, volume
    use status, only: FP, MU, SG11, SG12, SG21, SG22, &
     PFAST, PBLON, PBPER, NE, TE, TI, NI, VOLUM, XRHO, rho, srho, vr, vrs, sxho, MV
    use numerical_tools, only: derivcc, integrcc
    use stella_module, only: f_vacuum_rb

    double precision, intent(out) :: phi_edgehog

    integer :: niter, i_initto
    double precision :: f_boundary, errtol, tolerr, v_temp, H_b, &
        detS, f_vacuum, rocco
    double precision, dimension(NRD) :: phi, pressure, dphidv, pprimv, muhat, &
        a, b, c, dum, d, c1, c2, c3, c4, i1, f1, s22_temp, mu_temp
 
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
    muhat  = mu*phi_edgehog
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

        phi = phi_edgehog * xrho**2
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

end module a2vmec
