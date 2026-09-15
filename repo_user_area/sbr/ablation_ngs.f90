subroutine ABLATION_NGS(src, ndep)

!----------------------------------------------------------------------|
! Stellarator pellet ablation model (Neutral Gas Shielding).
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
!       equilibrium.  vmec/python/pellet.py traces that line through the
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
! ------------------------------------------ normalised output (pel_norm)
! pel_norm = .true. makes src the deposition SHAPE [m^-3] instead of a
! source, normalised so that VINT(CAR60B) = 1.  Multiplied by a fuelling
! rate in 10^19 particles/s it is a source in SN's units:
!     SN = CF16*CAR60;             ! CF16 = fuelling rate [10^19 part/s]
! and VINT(SN) = CF16.  Nothing is deposited: ndep comes back zero, and
! pel_tinject / pel_window play no part, the rate carries the timing.  The
! shape is the first pellet with a radius and a speed, traced on every call,
! 0.091 ms in 1D mode at NA1 = 101 (0.19 us per idle call with pel_norm off).
! Nothing ablated leaves CAR60 zero, with one warning.
!
! ------------------------------------------ plasmoid drift (pel_cdrift)
! The ionised cloud polarises in the 1/R field and drifts toward the low-field
! side before it spreads over the flux surface, so the atoms land away from
! where they were ablated, by
!     d = |pel_cdrift| * r_p[mm]^(4/3) * T_e[keV]^(1/2)   [m]
! (Lazerson, An Analytical Model of Pellet Fueling and Plasmoid Drift, 2026).
! pel_cdrift = 0, the default, deposits at the ablation point.  In 3D chord
! mode the direction is geometric, the projection of +R on grad rho there, so
! an inboard point drifts inward and an outboard one outward and the sign of
! pel_cdrift does not matter.  1D mode has no inboard/outboard, so there the
! sign carries the launch side: negative inward for HFS, positive outward for
! LFS.  A landing point outside [0, ROC] is clipped, never dropped.
!
! &pellet_ngs namelist (exp nml file):
!     npel        number of pellets
!     pel_model   'MCL23' (default) | 'NGS' | 'NGPS'
!     pel_scale   calibration multiplier
!     pel_amu     pellet atomic mass [amu]      (Parks laws)
!     pel_fh      protium (H) fraction f_H       (MCL23 material model)
!     pel_fd      deuterium fraction f_D         (MCL23 material model)
!                 tritium fraction f_T = 1 - f_H - f_D.  f_H = 0 is a pure
!                 (D,T) pellet; for a hydrogen pellet set pel_fh = 1.0
!                 and pel_fd = 0.0.
!     pel_pdens   solid atom density [m^-3]; <=0 -> compute from f_D
!     pel_radius  pellet radius [m]              (per pellet)
!     pel_vp      pellet speed  [m/s]            (per pellet)
!     pel_tinject injection time [s]             (per pellet)
!     pel_window  deposition window [s]
!     pel_cdrift  plasmoid drift scale [m]; 0 = deposit where it ablates
!     pel_norm    .true. -> src is the shape, not a source; see above
!     pel_ux      launch direction, Cartesian, per pellet.  Scaled to unit
!     pel_uy      length on read, so any length will do; all three zero
!     pel_uz      means unset.
!     pel_x0      injector axis, Cartesian [m], per pellet: pel_?0 is the
!     pel_y0      start point and pel_?1 the tip.  The line through them is
!     pel_z0      the path, so the pair replaces pel_theta/pel_phi and
!     pel_x1      pel_ux/uy/uz.  Coincident points mean unset.
!     pel_y1
!     pel_z1
!     pel_chord_file  rho(l)/modB(l) table -> 3D straight-line mode
!
! Launch, two ways: pel_theta and pel_phi place it on the LCFS with pel_ux/uy/uz
! for the direction, or the two injector points fix both.  pel_vp is the speed.
! vmec/python/pellet.py reads this group and writes the chord pel_chord_file
! names; it finds where the injector line crosses the LCFS.
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
!
! Call cadence: a postep hook's output PERSISTS until the routine is called
! again, so an interval-called hook (e.g. >:0.5:2.*TAU::;) holds src for the
! whole call stride.  The source is therefore normalised to the coming stride,
! estimated by the elapsed one (= TAU for a per-step hook), and a per-pellet
! deposited fraction is carried: each call deposits the part of the window
! that falls before TIME+stride and is not yet deposited, and nothing before
! pel_tinject.  A per-step hook reproduces the natoms/pel_window rate above.
! Under a call interval longer than the window every pellet's atoms are still
! conserved exactly, deposited at the first call after injection and
! time-aggregated over one stride, so at most one stride late.
! Call per step (>::::;) when the train structure matters; idle calls, with
! no pellet depositing, are near-free.  pel_norm traces on every call.
!----------------------------------------------------------------------|

use status, only: NRD
use scalars, only: NA1, HRO, ROC, TIME, TAU, BTOR, RTOR, SHIFT, UPDWN
use status, only: XRHO, RHO, VR, NE, TE, GRADRO
use read_input, only: AWD, nml_file

implicit none

double precision, intent(out) :: src(NRD), ndep

integer, parameter :: maxpel = 16
double precision, parameter :: pi = 3.141592653589793d0

! physical constants for the McClenaghan material / particle model
double precision, parameter :: A_H   = 1.00782503d0    ! H atomic weight [g/mol]
double precision, parameter :: A_D   = 2.014101778d0   ! D atomic weight [g/mol]
double precision, parameter :: A_T   = 3.0160492779d0  ! T atomic weight [g/mol]
double precision, parameter :: RHO_H = 88.d0           ! solid H density [kg/m^3]
double precision, parameter :: RHO_D = 200.d0          ! solid D density [kg/m^3]
double precision, parameter :: RHO_T = 318.d0          ! solid T density [kg/m^3]
double precision, parameter :: N_AVO = 6.02214076d23   ! Avogadro [1/mol]

! --- namelist pellet parameters --------------------------------------
integer          :: npel
character(len=8) :: pel_model
double precision :: pel_scale, pel_amu, pel_pdens, pel_window, pel_fd, pel_fh
double precision :: pel_cdrift
double precision :: pel_radius(maxpel), pel_vp(maxpel), pel_tinject(maxpel)
double precision :: pel_phi(maxpel), pel_theta(maxpel)
double precision :: pel_ux(maxpel), pel_uy(maxpel), pel_uz(maxpel)
double precision :: pel_x0(maxpel), pel_y0(maxpel), pel_z0(maxpel)
double precision :: pel_x1(maxpel), pel_y1(maxpel), pel_z1(maxpel)
logical          :: pel_cyl, pel_norm
character(len=160) :: pel_chord_file

namelist /pellet_ngs/ npel, pel_model, pel_scale, pel_amu, pel_pdens, &
                      pel_window, pel_fd, pel_fh, pel_cdrift, pel_cyl, pel_norm, &
                      pel_radius, pel_vp, pel_tinject, pel_phi, pel_theta, &
                      pel_ux, pel_uy, pel_uz, pel_x0, pel_y0, pel_z0, &
                      pel_x1, pel_y1, pel_z1, pel_chord_file

! the namelist is static during a run: read it once and keep it
save :: npel, pel_model, pel_scale, pel_amu, pel_pdens, pel_window, &
        pel_fd, pel_fh, pel_cdrift, pel_cyl, pel_norm, pel_radius, pel_vp, &
        pel_tinject, pel_phi, pel_theta, pel_ux, pel_uy, pel_uz, &
        pel_x0, pel_y0, pel_z0, pel_x1, pel_y1, pel_z1, pel_chord_file
logical, save :: nml_loaded = .false.

! --- 3D chord table: l, rho(l), modB(l) from PEL_CHORD_FILE -----------
integer, parameter :: maxchord = 8192
integer :: nchord
double precision :: chord_l(maxchord), chord_rho(maxchord), chord_b(maxchord)
double precision :: chord_R(maxchord), chord_Z(maxchord)
double precision :: chord_len
logical :: use3d

! --- chord cache, persists across calls ------------------------------
! The chord table is static while the equilibrium is frozen (vmec_option 7),
! so re-parsing the file for every active pellet on every timestep is wasted
! I/O.  Keyed on PEL_CHORD_FILE, so a changed filename forces a reload.
! NOTE: when VMEC runs live (vmec_option 0) vmec/python/pellet.py rewrites
! the chord on each equilibrium update, and this cache then serves the
! launch-time geometry for the whole run.
integer, save :: nchord_c(maxpel) = 0
double precision, save :: chord_l_c(maxchord, maxpel), &
                          chord_rho_c(maxchord, maxpel), &
                          chord_b_c(maxchord, maxpel), &
                          chord_R_c(maxchord, maxpel), &
                          chord_Z_c(maxchord, maxpel)
double precision, save :: chord_len_c(maxpel) = 0.d0
logical, save :: use3d_c(maxpel) = .false.
logical, save :: chord_ready(maxpel) = .false.
character(len=160), save :: chord_file_c(maxpel) = ''

! --- call-cadence bookkeeping (see "Call cadence" in the header) -------
double precision, save :: tlast_call = -1.d30    ! TIME of the previous call
double precision, save :: fdep(maxpel) = 0.d0    ! deposited window fraction

! --- locals ----------------------------------------------------------
integer, parameter :: ksub = 50
character(len=160) :: as_nml
integer :: ios, ip, j, ibin, ipnorm
logical :: use_mcl
double precision :: rp, vp, cpref, inv0
double precision :: natoms(NRD), dvol, dtwin, unorm
double precision :: molwt, pdens, massdens, fh, fd, ft, amol, fmass
double precision :: stride, f_target, df, shp_tot, drift_proj
logical, save :: warned_noabl = .false.

! NGS / NGPS Parks reference prefactors (ne in cm^-3, Te in eV, rp in cm)
double precision, parameter :: C_NGS  = 5.2d16
double precision, parameter :: C_NGPS = 6.1d16

src  = 0.d0
ndep = 0.d0

! ---- defaults, then read the namelist (first call only) -------------
if (.not. nml_loaded) then
    npel        = 0
    pel_model   = 'MCL23'
    pel_scale   = 1.d0
    pel_amu     = 2.d0
    pel_fh      = 0.d0
    pel_fd      = 0.5d0
    pel_pdens   = -1.d0        ! <=0 -> compute from the material model
    pel_window  = 1.0d-3
    pel_cdrift  = 0.d0         ! plasmoid drift scale [m], 0 = no drift
    pel_radius  = 0.d0
    pel_vp      = 0.d0
    pel_tinject = -1.d0
    pel_phi     = 0.d0
    pel_theta   = 0.d0
    pel_ux      = 0.d0         ! launch direction, all three zero = unset
    pel_uy      = 0.d0
    pel_uz      = 0.d0
    pel_x0      = 0.d0         ! injector start and tip, coincident = unset
    pel_y0      = 0.d0
    pel_z0      = 0.d0
    pel_x1      = 0.d0
    pel_y1      = 0.d0
    pel_z1      = 0.d0
    pel_cyl     = .false.
    pel_norm    = .false.
    pel_chord_file = ''

    as_nml = TRIM(AWD) // '/' // TRIM(nml_file)
    open(57, FILE=TRIM(as_nml), delim='apostrophe', status='old', iostat=ios)
    if (ios == 0) then
        read(57, nml=pellet_ngs, iostat=ios)
        close(57)
    endif
    npel = min(npel, maxpel)

    ! Unit length, so the chord tracer can take it as read.  All three
    ! zero means unset.
    do ip = 1, npel
        unorm = sqrt(pel_ux(ip)**2 + pel_uy(ip)**2 + pel_uz(ip)**2)
        if (unorm <= 0.d0) cycle
        pel_ux(ip) = pel_ux(ip)/unorm
        pel_uy(ip) = pel_uy(ip)/unorm
        pel_uz(ip) = pel_uz(ip)/unorm
    enddo

    nml_loaded = .true.
endif
if (npel <= 0) return

! ---- ablation-law selection + material / particle model -------------
use_mcl = .not. (TRIM(pel_model) == 'NGS'  .or. TRIM(pel_model) == 'ngs' .or. &
                 TRIM(pel_model) == 'NGPS' .or. TRIM(pel_model) == 'ngps')

if (TRIM(pel_model) == 'NGPS' .or. TRIM(pel_model) == 'ngps') then
    cpref = C_NGPS
else
    cpref = C_NGS
endif

! (H,D,T) mixture: molecular weight, solid mass & number density.
! Atom fractions f_H, f_D, f_T = 1-f_H-f_D; the solid mass density is the
! molar-volume mix  rho = (sum x_i A_i) / (sum x_i A_i/rho_i).
fh    = min(max(pel_fh, 0.d0), 1.d0)
fd    = min(max(pel_fd, 0.d0), 1.d0 - fh)
ft    = 1.d0 - fh - fd
molwt = fh*A_H + fd*A_D + ft*A_T
massdens = molwt / (fh*A_H/RHO_H + fd*A_D/RHO_D + ft*A_T/RHO_T)  ! kg/m^3
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
dtwin = max(pel_window, 1.d-9)

! Coming call stride: src written here persists until the NEXT call, i.e.
! for one stride.  Estimate it by the elapsed stride (uniform for an
! interval-called hook; = the timestep for a per-step hook), TAU on the
! first call.
if (tlast_call > -1.d29 .and. TIME > tlast_call) then
    stride = TIME - tlast_call
else
    stride = max(TAU, 1.d-12)
endif
tlast_call = TIME

! ---- normalised output: the deposition SHAPE and nothing else -------
! Nothing is deposited here, so ndep keeps the zero set above: an atom
! count beside a normalised profile is fuel the model could count twice.
! The shape is the first configured pellet at the present ne/Te.
if (pel_norm) then

    ipnorm = 0
    do ip = 1, npel
        if (pel_radius(ip) > 0.d0 .and. pel_vp(ip) > 0.d0) then
            ipnorm = ip
            exit
        endif
    enddo

    shp_tot = 0.d0
    if (ipnorm > 0) then
        call trace_pellet(ipnorm)
        do j = 1, NA1
            if (natoms(j) <= 0.d0) cycle
            if (VR(j)*HRO <= 0.d0) cycle
            shp_tot = shp_tot + natoms(j)
        enddo
    endif

    ! atoms per cell / (total atoms * cell volume), so that
    ! HRO*sum_j src(j)*VR(j) = 1, the sum VINT computes.  Only cells
    ! carrying volume take part, as in the physical branch.
    if (shp_tot > 0.d0) then
        do j = 1, NA1
            if (natoms(j) <= 0.d0) cycle
            dvol = VR(j) * HRO
            if (dvol <= 0.d0) cycle
            src(j) = natoms(j) / shp_tot / dvol
        enddo
    else if (.not. warned_noabl) then
        warned_noabl = .true.
        write(*,*) 'ABLATION_NGS: pel_norm is on and nothing ablates, the', &
                   ' normalised profile is zero -- check pel_radius, pel_vp', &
                   ' and the ne/Te at this time'
    endif

    return
endif

! ---- loop over pellets, deposit the not-yet-deposited window fraction
!      that falls before TIME+stride (see "Call cadence" in the header) --
do ip = 1, npel

    if (pel_tinject(ip) < 0.d0) cycle
    if (fdep(ip) >= 1.d0) cycle
    if (TIME < pel_tinject(ip)) cycle      ! causal: nothing before launch
    f_target = (TIME + stride - pel_tinject(ip)) / dtwin
    f_target = min(max(f_target, 0.d0), 1.d0)
    df = f_target - fdep(ip)
    if (df <= 0.d0) cycle

    if (pel_radius(ip) <= 0.d0 .or. pel_vp(ip) <= 0.d0) cycle

    call trace_pellet(ip)

    ! --- turn this call's share of the inventory into a source rate --
    ! df of the traced atoms, spread over the coming stride: integrated
    ! until the next call this deposits exactly df*natoms.  A per-step
    ! hook (df = TAU/pel_window, stride = TAU) recovers natoms/pel_window.
    do j = 1, NA1
        if (natoms(j) <= 0.d0) cycle
        dvol = VR(j) * HRO
        if (dvol <= 0.d0) cycle
        src(j) = src(j) + natoms(j) * df / dvol / stride * 1.d-19
        ndep   = ndep + natoms(j) * df
    enddo
    fdep(ip) = f_target

    if (df >= 1.d0 .and. stride > 2.d0*dtwin) write(*,'(A,I0,2A,F8.5,A,F8.5,A)') &
        ' ABLATION_NGS: pellet ', ip, ' time-aggregated into one interval call', &
        ' (call stride ', stride, ' s >> window ', dtwin, ' s)'

enddo

return

contains

    !------------------------------------------------------------------|
    ! Trace pellet ip across the plasma at the present ne/Te and bin what
    ! it ablates into natoms.  Sets the host rp, vp, inv0 and the chord
    ! state.  Both outputs come from this one trace.
    !------------------------------------------------------------------|
    subroutine trace_pellet(ip)
    integer, intent(in) :: ip
    integer :: istep, nstep
    double precision :: ds, dt_step, step_len, rpos, pathpos, gr_l, Bloc

    rp = pel_radius(ip)
    vp = pel_vp(ip)

    ! A4: pel_radius is the sphere-equivalent radius the NGS model assumes.
    ! If pel_cyl, interpret it as the radius of an equal-aspect (length =
    ! diameter) cylindrical pellet and convert to the volume-equivalent
    ! sphere: (4/3)pi r_s^3 = pi/4 d^2 l, d = l = 2 r_cyl -> r_s = 1.5^(1/3) r_cyl.
    if (pel_cyl) rp = rp * 1.5d0**(1.d0/3.d0)

    ! B5: load this pellet's straight-line chord (per-pellet launch).  Tries
    ! <stem>_<ip>.<ext> first, then the plain PEL_CHORD_FILE (shared/global),
    ! reusing the cached table if this pellet's filename hasn't changed.
    use3d   = .false.
    nchord  = 0
    chord_len = 0.d0
    if (LEN_TRIM(pel_chord_file) > 0) then
        if (chord_ready(ip) .and. TRIM(chord_file_c(ip)) == TRIM(pel_chord_file)) then
            use3d     = use3d_c(ip)
            nchord    = nchord_c(ip)
            chord_len = chord_len_c(ip)
            if (use3d) then
                chord_l(1:nchord)   = chord_l_c(1:nchord, ip)
                chord_rho(1:nchord) = chord_rho_c(1:nchord, ip)
                chord_b(1:nchord)   = chord_b_c(1:nchord, ip)
                chord_R(1:nchord)   = chord_R_c(1:nchord, ip)
                chord_Z(1:nchord)   = chord_Z_c(1:nchord, ip)
            endif
        else
            call load_chord(ip)
            chord_ready(ip)  = .true.
            chord_file_c(ip) = pel_chord_file
            use3d_c(ip)      = use3d
            nchord_c(ip)     = nchord
            chord_len_c(ip)  = chord_len
            if (use3d) then
                chord_l_c(1:nchord, ip)   = chord_l(1:nchord)
                chord_rho_c(1:nchord, ip) = chord_rho(1:nchord)
                chord_b_c(1:nchord, ip)   = chord_b(1:nchord)
                chord_R_c(1:nchord, ip)   = chord_R(1:nchord)
                chord_Z_c(1:nchord, ip)   = chord_Z(1:nchord)
            endif
        endif
    endif

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
            drift_proj = lfs_proj(0.5d0*(chord_R(istep)+chord_R(istep-1)), &
                                  0.5d0*(chord_Z(istep)+chord_Z(istep-1)))
            call ablate_at(rpos, Bloc, dt_step)
        enddo
    else
        ! ---- 1D: flux-averaged trace, edge -> axis -> edge ----------
        ! No inboard/outboard here, so the user carries the launch side in
        ! the sign of pel_cdrift.
        drift_proj = sign(1.d0, pel_cdrift)
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
    end subroutine trace_pellet

    !------------------------------------------------------------------|
    ! Projection of the outward major-radius direction, the one the
    ! plasmoid drifts along, on grad rho at (Rl, Zl): +1 on the outboard
    ! midplane, -1 on the inboard one, 0 top and bottom.  Flux surfaces
    ! taken as circles about (RTOR+SHIFT, UPDWN).
    !------------------------------------------------------------------|
    double precision function lfs_proj(Rl, Zl)
    double precision, intent(in) :: Rl, Zl
    double precision :: dR, dZ, rpol
    dR   = Rl - (RTOR + SHIFT)
    dZ   = Zl - UPDWN
    rpol = sqrt(dR*dR + dZ*dZ)
    if (rpol <= 1.d-6) then
        lfs_proj = 0.d0
    else
        lfs_proj = dR / rpol
    endif
    end function lfs_proj

    !------------------------------------------------------------------|
    ! Ablate the pellet over one path element of duration dt_step at
    ! flux position rpos (ASTRA rho) in local field Bloc; shrink rp and
    ! bin the ablated atoms onto the transport grid.  Uses host rp,
    ! natoms, inv0, use_mcl, cpref, pel_scale, amol, fmass.
    !------------------------------------------------------------------|
    subroutine ablate_at(rpos, Bloc, dt_step)
    double precision, intent(in) :: rpos, Bloc, dt_step
    integer :: jl, jd
    double precision :: fj, frac, ne_l, te_l, ne_cm3, te_ev, ne20, tkev
    double precision :: gdot, dMdt, dN, rp3, Bf, rdep, wd

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

    ! Plasmoid drift: the ionised cloud moves toward the low-field side before
    ! it spreads over the flux surface, so the atoms land off the ablation
    ! point (see "plasmoid drift" in the header).  Clipped to the plasma, so
    ! the count is conserved either way.
    rdep = rpos
    if (pel_cdrift /= 0.d0) rdep = min(max(rpos + drift_proj * abs(pel_cdrift) &
        * (rp*1.d3)**(4.d0/3.d0) * sqrt(te_l), 0.d0), ROC)

    ! shrink the pellet, limited by remaining inventory
    rp3 = rp*rp*rp - dN * inv0
    if (rp3 <= 0.d0) then
        dN = (rp*rp*rp) / inv0
        rp = 0.d0
    else
        rp = rp3**(1.d0/3.d0)
    endif

    ! Deposit the ablated atoms with linear (cloud-in-cell) weighting between
    ! the two cells bracketing rdep.  This conserves the atom count (weights
    ! sum to 1) and avoids the nearest-grid-point aliasing that made the
    ! deposition profile spiky when the ~500 chord sub-steps are binned onto
    ! the ~NA1 transport cells.
    jd = min(max(int(rdep/HRO) + 1, 1), NA1 - 1)
    wd = min(max(rdep/HRO - dble(jd - 1), 0.d0), 1.d0)
    natoms(jd)   = natoms(jd)   + dN * (1.d0 - wd)
    natoms(jd+1) = natoms(jd+1) + dN * wd
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
    ! Load the 3D straight-line chord written by vmec/python/pellet.py.
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
    double precision :: xl, xr, xs, xb, xmaj, xvrt

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
        read(line, *, iostat=ic) xl, xr, xs, xb, xmaj, xvrt
        if (ic /= 0) cycle
        k = k + 1
        chord_l(k)   = xl
        chord_rho(k) = xr
        chord_b(k)   = xb
        chord_R(k)   = xmaj
        chord_Z(k)   = xvrt
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
