subroutine ABLATION_NGS(src, ndep)

!----------------------------------------------------------------------|
! Stellarator pellet ablation model (Neutral Gas Shielding).
!
! Port of the THRIFT/STELLOPT NGS pellet model and of F. Wilms'
! NGS_pellet_model.py (D. Kulla 2026) to ASTRA's 1D flux-surface grid.
!
! A frozen hydrogenic pellet is launched at the plasma edge and traced
! across the flux surfaces.  On every surface the local electron density
! / temperature (and, for the McClenaghan law, the local field B) drive
! the ablation rate that shrinks the pellet; the ablated atoms are
! deposited as a volumetric electron source [10^19/(m^3 s)] (electrons
! follow by quasineutrality, Z=1 hydrogenic) that can be added into SN.
!
! ---------------------------------------------------------- ablation law
! PEL_MODEL selects the scaling:
!
!   'MCL23' (default) -- McClenaghan et al 2023 NF 63 036015, Eq. (1),
!       a MASS ablation rate with an explicit field dependence:
!         dM/dt [g/s] = 39 * (2/B)^0.35 * (A/A_D)^(2/3)
!                          * (r_cm/0.2)^(4/3) * (n_e/1e20)^(1/3)
!                          * (T_keV/2)^1.64
!         dN/dt       = (dM/dt)/A * N_A            [atoms/s]
!       A = molecular weight of the (D,T) mix, A_D = 2.0141.
!
!   'NGS' / 'NGPS' -- Parks neutral-gas-shielding scaling (no B):
!         dN/dt = C * n_e[cm^-3]^0.333 * T_e[eV]^1.64
!                   * r_p[cm]^1.333 * (2/amu)^(1/3)
!       C = 5.2e16 (NGS) / 6.1e16 (NGPS).
!
! Both share the same n_e^1/3, T^1.64, r^4/3 exponents; MCL23 adds the
! (2/B)^0.35 field factor and a physically normalised prefactor.  The
! absolute prefactor is first-order; calibrate through PEL_SCALE.
!
! ------------------------------------------------------------ geometry
! Two geometry modes, selected by whether PEL_CHORD_FILE is present:
!
!   (1) 1D flux-surface-averaged (no chord file)
!       Traced in the flux label rho; residence time in a rho-step is
!       dt = ds/(vp*<|grad rho|>) via GRADRO.  The field for MCL23 is the
!       reference BTOR (no local variation).  (phi,theta) do not enter.
!
!   (2) 3D straight-line chord (PEL_CHORD_FILE set)
!       The pellet flies along the real straight line fired at the VMEC
!       equilibrium.  src/for/VMEC2ASTRA.py traces that line through the
!       3D metric and tabulates, per real path length l, the flux label
!       rho(l) and the LOCAL field modB(l) into PEL_CHORD_FILE.  Here the
!       pellet is stepped along the tabulated path (dt = dl/vp) using the
!       local rho and B, so both penetration and the B-dependent ablation
!       follow the true 3D geometry.  The table is the FULL path (edge ->
!       deepest point -> exit); a pellet that survives to the far edge
!       deposits only what it ablated (partial ablation).
!
! Usage in an equ model (postep hook, like ABLATION):
!     ABLATION_NGS(CAR60, CF15)>::::;
!     SN = SNEBM + CAR60;          ! CAR60 = pellet source, CF15 = atoms
!
! &pellet_ngs namelist (exp nml file):
!     npel        number of pellets
!     pel_model   'MCL23' (default) | 'NGS' | 'NGPS'
!     pel_scale   calibration multiplier
!     pel_amu     pellet atomic mass [amu]      (Parks laws)
!     pel_fd      deuterium fraction f_D         (MCL23 material model)
!     pel_pdens   solid atom density [m^-3]; <=0 -> compute from f_D
!     pel_radius  pellet radius [m]              (per pellet)
!     pel_vp      pellet speed  [m/s]            (per pellet)
!     pel_tinject injection time [s]             (per pellet)
!     pel_window  deposition window [s]
!     pel_chord_file  rho(l)/modB(l) table -> 3D straight-line mode
!
! Timing model (pel_tinject, pel_window): a single pellet's whole ablated
! inventory is deposited as a constant-in-time source rate spread over the
! window [pel_tinject, pel_tinject+pel_window]:  src = natoms/dvol/pel_window.
! Integrated over the window this returns exactly natoms, so pel_window is a
! NUMERICAL smear (it sets the source *rate*, not the total), chosen for a
! gentle density source rather than the physical ~L/v (sub-ms) ablation time.
! This matches the window-smear representation used in PELPOL.  A repetitive
! pellet train is not modelled as a frequency here; inject several discrete
! pellets via per-pellet pel_tinject instead.
!----------------------------------------------------------------------|

use status, only: NRD
use scalars, only: NA1, HRO, ROC, TIME, TAU, BTOR
use status, only: XRHO, RHO, VR, NE, TE, GRADRO
use read_input, only: AWD, nml_file

implicit none

double precision, intent(out) :: src(NRD), ndep

integer, parameter :: maxpel = 16
double precision, parameter :: pi = 3.141592653589793d0

! physical constants for the McClenaghan material / particle model
double precision, parameter :: A_D   = 2.014101778d0   ! D atomic weight [g/mol]
double precision, parameter :: A_T   = 3.0160492779d0  ! T atomic weight [g/mol]
double precision, parameter :: RHO_D = 200.d0          ! solid D density [kg/m^3]
double precision, parameter :: RHO_T = 318.d0          ! solid T density [kg/m^3]
double precision, parameter :: N_AVO = 6.02214076d23   ! Avogadro [1/mol]

! --- namelist pellet parameters --------------------------------------
integer          :: npel
character(len=8) :: pel_model
double precision :: pel_scale, pel_amu, pel_pdens, pel_window, pel_fd
double precision :: pel_radius(maxpel), pel_vp(maxpel), pel_tinject(maxpel)
double precision :: pel_phi(maxpel), pel_theta(maxpel)
logical          :: pel_cyl
character(len=160) :: pel_chord_file

namelist /pellet_ngs/ npel, pel_model, pel_scale, pel_amu, pel_pdens, &
                      pel_window, pel_fd, pel_cyl, pel_radius, pel_vp, &
                      pel_tinject, pel_phi, pel_theta, pel_chord_file

! --- 3D chord table: l, rho(l), modB(l) from PEL_CHORD_FILE -----------
integer, parameter :: maxchord = 8192
integer :: nchord
double precision :: chord_l(maxchord), chord_rho(maxchord), chord_b(maxchord)
double precision :: chord_len
logical :: use3d

! --- locals ----------------------------------------------------------
integer, parameter :: ksub = 50
character(len=160) :: as_nml
integer :: ios, ip, j, istep, nstep, ibin
logical :: use_mcl
double precision :: rp, vp, cpref, inv0, dt_step, ds, step_len
double precision :: natoms(NRD), dvol, dtwin, gr_l, rpos, pathpos, Bloc
double precision :: molwt, pdens, massdens, fd, ft, amol, fmass

! NGS / NGPS Parks reference prefactors (ne in cm^-3, Te in eV, rp in cm)
double precision, parameter :: C_NGS  = 5.2d16
double precision, parameter :: C_NGPS = 6.1d16

src  = 0.d0
ndep = 0.d0

! ---- defaults, then read the namelist -------------------------------
npel        = 0
pel_model   = 'MCL23'
pel_scale   = 1.d0
pel_amu     = 2.d0
pel_fd      = 0.5d0
pel_pdens   = -1.d0            ! <=0 -> compute from the material model
pel_window  = 1.0d-3
pel_radius  = 0.d0
pel_vp      = 0.d0
pel_tinject = -1.d0
pel_phi     = 0.d0
pel_theta   = 0.d0
pel_cyl     = .false.
pel_chord_file = ''

as_nml = TRIM(AWD) // '/' // TRIM(nml_file)
open(57, FILE=TRIM(as_nml), delim='apostrophe', status='old', iostat=ios)
if (ios == 0) then
    read(57, nml=pellet_ngs, iostat=ios)
    close(57)
endif
if (npel <= 0) return
npel = min(npel, maxpel)

! ---- ablation-law selection + material / particle model -------------
use_mcl = .not. (TRIM(pel_model) == 'NGS'  .or. TRIM(pel_model) == 'ngs' .or. &
                 TRIM(pel_model) == 'NGPS' .or. TRIM(pel_model) == 'ngps')

if (TRIM(pel_model) == 'NGPS' .or. TRIM(pel_model) == 'ngps') then
    cpref = C_NGPS
else
    cpref = C_NGS
endif

! (D,T) mixture: molecular weight, solid mass & number density
fd    = min(max(pel_fd, 0.d0), 1.d0)
ft    = 1.d0 - fd
molwt = fd*A_D + ft*A_T
massdens = (1.d0 - fd + fd*(A_D/A_T)) / &
           ((1.d0 - fd)/RHO_T + fd*(A_D/A_T)/RHO_D)          ! kg/m^3
if (pel_pdens > 0.d0) then
    pdens = pel_pdens                                        ! explicit override
else
    pdens = massdens / 1.d-3 / molwt * N_AVO                 ! m^-3
endif

! atomic weight used by each law: (D,T) molecular weight for MCL23,
! the namelist pel_amu for the Parks (2/amu) factor
if (use_mcl) then
    amol = molwt
else
    amol = max(pel_amu, 1.d-3)
endif
fmass = (2.d0/max(amol,1.d-3))**(1.d0/3.d0)                  ! Parks mass factor
dtwin = max(pel_window, TAU)

! ---- loop over pellets, deposit those active in this window ---------
do ip = 1, npel

    if (pel_tinject(ip) < 0.d0) cycle
    if (TIME < pel_tinject(ip) .or. TIME >= pel_tinject(ip) + dtwin) cycle

    rp = pel_radius(ip)
    vp = pel_vp(ip)
    if (rp <= 0.d0 .or. vp <= 0.d0) cycle

    ! A4: pel_radius is the sphere-equivalent radius the NGS model assumes.
    ! If pel_cyl, interpret it as the radius of an equal-aspect (length =
    ! diameter) cylindrical pellet and convert to the volume-equivalent
    ! sphere: (4/3)pi r_s^3 = pi/4 d^2 l, d = l = 2 r_cyl -> r_s = 1.5^(1/3) r_cyl.
    if (pel_cyl) rp = rp * 1.5d0**(1.d0/3.d0)

    ! B5: load this pellet's straight-line chord (per-pellet launch).  Tries
    ! <stem>_<ip>.<ext> first, then the plain PEL_CHORD_FILE (shared/global).
    use3d   = .false.
    nchord  = 0
    chord_len = 0.d0
    if (LEN_TRIM(pel_chord_file) > 0) call load_chord(ip)

    inv0   = 1.d0 / ((4.d0/3.d0) * pi * pdens)   ! d(rp^3) = -dN*inv0
    natoms = 0.d0

    if (use3d) then
        ! ---- 3D: march the tabulated straight-line path ONCE --------
        do istep = 2, nchord
            if (rp <= 0.d0) exit
            ds      = chord_l(istep) - chord_l(istep-1)          ! dl [m]
            if (ds <= 0.d0) cycle
            dt_step = ds / vp
            rpos    = 0.5d0*(chord_rho(istep)+chord_rho(istep-1)) * ROC
            Bloc    = 0.5d0*(chord_b(istep)  +chord_b(istep-1))
            call ablate_at(rpos, Bloc, dt_step)
        enddo
    else
        ! ---- 1D: flux-averaged trace, edge -> axis -> edge ----------
        step_len = HRO / dble(ksub)
        nstep    = 2 * (NA1 - 1) * ksub
        do istep = 1, nstep
            if (rp <= 0.d0) exit
            pathpos = dble(istep) * step_len
            rpos    = abs(ROC - pathpos)
            if (rpos > ROC) rpos = ROC
            gr_l    = gradro_at(rpos)
            dt_step = step_len / (vp * max(gr_l, 1.d-3))
            call ablate_at(rpos, BTOR, dt_step)
        enddo
    endif

    ! --- turn deposited atoms into a volumetric source rate ----------
    do j = 1, NA1
        if (natoms(j) <= 0.d0) cycle
        dvol = VR(j) * HRO
        if (dvol <= 0.d0) cycle
        src(j) = src(j) + natoms(j) / dvol / dtwin * 1.d-19
        ndep   = ndep + natoms(j)
    enddo

enddo

return

contains

    !------------------------------------------------------------------|
    ! Ablate the pellet over one path element of duration dt_step at
    ! flux position rpos (ASTRA rho) in local field Bloc; shrink rp and
    ! bin the ablated atoms onto the transport grid.  Uses host rp,
    ! natoms, inv0, use_mcl, cpref, pel_scale, amol, fmass.
    !------------------------------------------------------------------|
    subroutine ablate_at(rpos, Bloc, dt_step)
    double precision, intent(in) :: rpos, Bloc, dt_step
    integer :: jl
    double precision :: fj, frac, ne_l, te_l, ne_cm3, te_ev, ne20, tkev
    double precision :: gdot, dMdt, dN, rp3, Bf

    ! linear interpolation of NE/TE on the uniform RHO grid
    fj   = rpos / HRO
    jl   = min(max(int(fj) + 1, 1), NA1 - 1)
    frac = min(max(rpos / HRO - dble(jl - 1), 0.d0), 1.d0)
    ne_l = NE(jl) + (NE(jl+1) - NE(jl)) * frac
    te_l = TE(jl) + (TE(jl+1) - TE(jl)) * frac
    if (ne_l <= 0.d0 .or. te_l <= 0.d0) return

    if (use_mcl) then
        ! McClenaghan 2023: mass rate [g/s] -> particle rate [1/s]
        Bf   = max(abs(Bloc), 0.1d0)
        ne20 = ne_l * 0.1d0                     ! 1e19 m^-3 -> 1e20 m^-3
        tkev = te_l                             ! keV
        dMdt = 39.d0 * (2.d0/Bf)**0.35d0 &
                     * (amol/A_D)**(2.d0/3.d0) &
                     * (rp*1.d2/0.2d0)**(4.d0/3.d0) &
                     * ne20**(1.d0/3.d0) &
                     * (tkev/2.d0)**1.64d0
        dN   = pel_scale * (dMdt/amol) * N_AVO * dt_step   ! atoms
    else
        ! Parks NGS/NGPS: particle rate [1/s]
        ne_cm3 = ne_l * 1.d13                   ! 1e19 m^-3 -> cm^-3
        te_ev  = te_l * 1.d3                     ! keV -> eV
        gdot   = pel_scale * cpref * ne_cm3**0.333d0 * te_ev**1.64d0 &
                            * (rp*1.d2)**1.333d0 * fmass
        dN     = gdot * dt_step
    endif
    if (dN <= 0.d0) return

    ! shrink the pellet, limited by remaining inventory
    rp3 = rp*rp*rp - dN * inv0
    if (rp3 <= 0.d0) then
        dN = (rp*rp*rp) / inv0
        rp = 0.d0
    else
        rp = rp3**(1.d0/3.d0)
    endif

    ! Deposit the ablated atoms with linear (cloud-in-cell) weighting between
    ! the two bracketing cells jl, jl+1, reusing the NE/TE interpolation
    ! indices.  This conserves the atom count (weights sum to 1) and avoids the
    ! nearest-grid-point aliasing that made the deposition profile spiky when
    ! the ~500 chord sub-steps are binned onto the ~NA1 transport cells.
    natoms(jl)   = natoms(jl)   + dN * (1.d0 - frac)
    natoms(jl+1) = natoms(jl+1) + dN * frac
    end subroutine ablate_at

    !------------------------------------------------------------------|
    ! <|grad rho|> interpolated on the RHO grid at ASTRA rho = rpos.
    !------------------------------------------------------------------|
    double precision function gradro_at(rpos)
    double precision, intent(in) :: rpos
    integer :: jl
    double precision :: frac
    jl   = min(max(int(rpos/HRO) + 1, 1), NA1 - 1)
    frac = min(max(rpos/HRO - dble(jl - 1), 0.d0), 1.d0)
    gradro_at = GRADRO(jl) + (GRADRO(jl+1) - GRADRO(jl)) * frac
    end function gradro_at

    !------------------------------------------------------------------|
    ! Load the 3D straight-line chord written by VMEC2ASTRA.py.
    ! Format (comment lines start with '#'):
    !     # header ...
    !     <nchord>
    !     l[m]  rho[-]  s[-]  modB[T]  R[m]  Z[m]     (l ascending)
    ! Reads l, rho, modB (columns 1,2,4); s,R,Z ignored.  On success
    ! sets use3d=.true.
    !------------------------------------------------------------------|
    subroutine load_chord(ip)
    integer, intent(in) :: ip
    character(len=160) :: chord_path, perpel, line
    character(len=16)  :: ipstr
    integer :: k, ic, dotpos
    double precision :: xl, xr, xs, xb

    ! per-pellet file:  <stem>_<ip>.<ext>  (insert before the extension)
    write(ipstr, '(I0)') ip
    dotpos = index(TRIM(pel_chord_file), '.', back=.true.)
    if (dotpos > 1) then
        perpel = pel_chord_file(1:dotpos-1)//'_'//TRIM(ipstr)//pel_chord_file(dotpos:)
    else
        perpel = TRIM(pel_chord_file)//'_'//TRIM(ipstr)
    endif

    ! resolve to abs path and try the per-pellet file, then the plain stem
    call resolve_chord(perpel, chord_path)
    open(58, FILE=TRIM(chord_path), status='old', iostat=ic)
    if (ic /= 0) then
        call resolve_chord(pel_chord_file, chord_path)
        open(58, FILE=TRIM(chord_path), status='old', iostat=ic)
    endif
    if (ic /= 0) then
        write(*,*) 'ABLATION_NGS: no chord file for pellet ', ip, &
                   ', using 1D mode'
        return
    endif

    nchord = 0
    do
        read(58, '(A)', iostat=ic) line
        if (ic /= 0) exit
        line = ADJUSTL(line)
        if (LEN_TRIM(line) == 0) cycle
        if (line(1:1) == '#' .or. line(1:1) == '!') cycle
        read(line, *, iostat=ic) nchord
        exit
    enddo
    if (nchord <= 1 .or. nchord > maxchord) then
        write(*,*) 'ABLATION_NGS: bad chord count, using 1D mode:', nchord
        close(58); nchord = 0; return
    endif

    k = 0
    do
        read(58, '(A)', iostat=ic) line
        if (ic /= 0) exit
        line = ADJUSTL(line)
        if (LEN_TRIM(line) == 0) cycle
        if (line(1:1) == '#' .or. line(1:1) == '!') cycle
        read(line, *, iostat=ic) xl, xr, xs, xb
        if (ic /= 0) cycle
        k = k + 1
        chord_l(k)   = xl
        chord_rho(k) = xr
        chord_b(k)   = xb
        if (k >= nchord) exit
    enddo
    close(58)

    if (k < 2) then
        write(*,*) 'ABLATION_NGS: too few chord points, using 1D mode'
        nchord = 0; return
    endif
    nchord    = k
    chord_len = chord_l(nchord)
    use3d     = .true.
    write(*,'(A,I0,A,I5,A,F7.4,A,F6.3,2A)') &
          ' ABLATION_NGS: pellet ', ip, ' 3D straight-line chord, ', &
          nchord, ' pts, length ', chord_len, ' m, deepest rho ', &
          minval(chord_rho(1:nchord)), '  from ', TRIM(chord_path)
    end subroutine load_chord

    !------------------------------------------------------------------|
    ! Resolve a chord filename to an absolute path (prepend AWD unless
    ! already absolute).
    !------------------------------------------------------------------|
    subroutine resolve_chord(fname, path)
    character(len=*), intent(in)  :: fname
    character(len=*), intent(out) :: path
    if (fname(1:1) == '/') then
        path = TRIM(fname)
    else
        path = TRIM(AWD) // '/' // TRIM(fname)
    endif
    end subroutine resolve_chord

end subroutine ABLATION_NGS
