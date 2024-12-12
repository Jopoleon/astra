subroutine qlknn_serial(CHI, CHE, VIN)

use parameter_inc, only: NRD

use const_inc, only: BTOR, RTOR, AMJ, AIM1, AIM2, AIM3, ZMJ, &
    NA1, NA1N, NA1E, NA1I
use status_inc, only: NE, TE, NI, TI, &
    ZEF, ZIM1, ZIM2, ZIM3, PBLON, PBPER, &
    PFAST, NIZ3, AMAIN, ER, MU, &
    RHO, AMETR, SHIF, ELON, &
    NDEUT, NIZ1, NTRIT, NIZ2, NHE3, &
    TRIA, VTOR, NIBM, G11, VRS, SHEAR

USE qlknn_evaluate_nets, only: lli, qlknn_options, qlknn_normpars, evaluate_jetexp_net, &
    default_qlknn_hyper_options, default_qlknn_hornnet_options, &
    default_qlknn_fullflux_options, default_qlknn_jetexp_options, &
    evaluate_qlknn_10d, evaluate_hornnet_constants, &
    hornnet_flux_from_constants, evaluate_fullflux_net
USE qlknn_disk_io, only: nets, blocks, &
   load_qlknn_hyper_nets_from_disk, load_hornnet_nets_from_disk, &
   load_fullflux_nets_from_disk   , load_jetexp_nets_from_disk
USE, INTRINSIC :: IEEE_ARITHMETIC, only: IEEE_IS_NAN, IEEE_IS_FINITE

implicit none

double precision, intent(out), dimension(NRD) :: CHI, CHE, VIN
!--------------------------------

integer, parameter :: dimx=1, nradial=100, nspec_max=7

real, parameter :: &
   e0   = 4.8032E-10,      &    ! elementary charge (statcoulombs)
   e00  = 1.602176462e-19, &    ! elementary charge (C)
   c0   = 2.9979E+10,      &    ! speed of light (cm/sec)
   mpp = 1.67262158E-27,   &
   z_electronmass = 5.44631e-4  ! Mass_e/mass_p

!-----------------------------------------

integer :: nions, rot_flag_in, verbose_in

double precision :: R0_in, collmultin, rhomin, rhomax, rhoscale
double precision, dimension(dimx) :: x_in, rho_in, Ro_in, Rmin_in, Bo_in, &
    qx_in, smag_in, alphax_in, Tex_in, Nex_in, Ate_in, Ane_in, &
    Machtor_in, Autor_in, Machpar_in, Aupar_in, gammaE_in, &
    epf_GB_out, eef_GB_out, x_in_eff, Zeffx, epsilon
double precision, dimension(dimx, nspec_max-1) :: Tix_in, ninorm_in, &
    Ati_in, Ani_in, Ai_in, Zi_in, ipf_GB_out, ief_GB_out
integer :: NN_ETG, NN_part_trans_switch, NN_TEM, NN_MODEL

integer :: n_red, n_step, jna, jrho, jr_l, jr_r
integer :: i, j, k, ihi, ilow, jion

double precision :: bmod, bpolz
double precision :: drmin, drmaj, drho, dte, dne, dq, dptot, &
        delong, dtrian, dvper, drhodr, dstep, dr, dv_r
double precision :: Bunit, cs00, rhos00, omega0, rhostar2
double precision :: T0, m0, drho_cs
! Shifted cicle geometry inputs
double precision :: gamma_e_tg, vpar_tg, mach_fac, ql_fac

double precision, dimension(NRD) :: vexb2, vpar_m, vper_m, &
    gradrhosq_exp, rmaj_exp, q_exp, &
    chie_m, chii_m, pfluxi_m, ptot

double precision, dimension(nradial) :: chie, chii, pfluxi, rho_tg
double precision, dimension(nspec_max-1) :: dti, dni
double precision, dimension(nspec_max-1, NRD) :: ni_m, ti_m
double precision :: vpar_in, vpar_shear_in, cexb

CHARACTER(len=20) :: fmtnets
character(len=120) :: AEXT

TYPE (qlknn_options), SAVE :: qlknn_opts ! for QuaLiKiz Neural Network options
TYPE (qlknn_normpars), SAVE :: qlknn_norms ! for normalisation conversion of gammaE inside QLKNN
REAL, ALLOCATABLE, DIMENSION(:, :) :: qlknn_in
REAL, ALLOCATABLE, DIMENSION(:, :) :: qlknn_out
REAL, DIMENSION(dimx, 15) :: qlknn_hornnet_constants
REAL, DIMENSION(dimx, 13) :: qlknn_eb
LOGICAL, DIMENSION(dimx, 13) :: qlknn_validity
LOGICAL, DIMENSION(13) :: qlknn_validity_mask
LOGICAL, DIMENSION(dimx) :: qlknn_rho_validity
INTEGER, DIMENSION(dimx, 13), SAVE :: n_qlknn_invalid ! counter for NN validity check failures. Writes every maxruns times
INTEGER(lli), SAVE :: qlknn_members = 10
LOGICAL, SAVE :: first_run = .TRUE. ! for initialising NN readout
INTEGER, SAVE :: freqio = 0 ! counter for NN write-to-disk. Writes every maxruns times
LOGICAL :: qlknn_out_of_bounds ! flag for running original QLK model when committee NN variance indicates poor NN predictive capability
INTEGER :: ii
REAL, DIMENSION(dimx) :: Lambe, Nue, Nustar, filterprof
REAL :: epsilonNN ! Epsilon constant of the trained NN

character(len=132) :: qlknn_sets_dir

!---------------------------------

call GETENV('ASTRA_EXT', AEXT)
qlknn_sets_dir = TRIM(AEXT) // '/qlk_nn/nov24/'

nions = nspec_max - 1

jna = max(NA1E, NA1I, NA1N)

! Electrons and main ions
Zi_in(1, 1) = ZMJ

Ai_in(1, 1) = AMJ  ! AMJ is reference mass
Ai_in(1, 2) = AIM1
Ai_in(1, 3) = AIM2
Ai_in(1, 4) = AIM3

do jrho=1, NA1
    ti_m(1:4, jrho) = TI(jrho)
    ni_m(1, jrho) = NDEUT(jrho)
    ni_m(2, jrho) = max(1.e-9, NIZ1(jrho))
    ni_m(3, jrho) = max(1.e-9, NIZ2(jrho))
    ni_m(4, jrho) = max(1.e-9, NIZ3(jrho))
    rmaj_exp(jrho) = RTOR + SHIF(jrho)
    q_exp(jrho)    = 1./MU(jrho)
    ptot(jrho) = NE(jrho)*TE(jrho) + ni_m(1, jrho)*ti_m(1, jrho) + ni_m(2, jrho)*ti_m(2, jrho) + pfast(jrho) + 0.5*(pblon(jrho) + pbper(jrho))
    bpolz = BTOR*AMETR(jrho)*MU(jrho)/RTOR
    bmod = sqrt(BTOR**2 + bpolz**2)
    gradrhosq_exp(jrho) = G11(jrho)/VRS(jrho)

    vper_m(jrho) = ER(jrho)/(RTOR*bpolz) ! vexb in m/s --> Omega_E
    vpar_m(jrho) = ER(jrho)/(RTOR*bpolz)*(RTOR+SHIF(jrho)+AMETR(jrho))  !--> R*Omega_E , no neoclassical terms
    vexb2(jrho)  = -ER(jrho)/bmod ! vexb in m/s (vperp = vexb since the diamagnetic velocity is the curvature drift actually

enddo

! GYRO conventions

m0 = AMJ*mpp*1000         ! Ref. mass = D ion mass [g]

pfluxi_m  = 0.0
chie_m  = 0.0
chii_m  = 0.0

! These will be reset locally in the radial loop
Zi_in(1, 2) = MAXVAL(ZIM1(1:NA1))
Zi_in(1, 3) = MAXVAL(ZIM2(1:NA1))
Zi_in(1, 4) = MAXVAL(ZIM3(1:NA1))

if (Zi_in(1, 4) .ge. 1.) nions = 4
if (Zi_in(1, 4) .lt. 1.) nions = 3
if (Zi_in(1, 3) .lt. 1.) nions = 2
if (Zi_in(1, 4) .ge. 1. .and. nions .eq. 2) then
    nions = 3
endif
if (Zi_in(1, 2) .lt. 1.) nions = 1
if (Zi_in(1, 3) .ge. 1. .and. nions .eq. 1) then
    nions = 2
endif

vpar_in = 0.
vpar_shear_in = 0.
rhomin = 0.
rhomax = 1.
verbose_in  = 0
collmultin = 1.
rot_flag_in = 0  ! git / Originally: 2; 0 suggested by Pierre
NN_ETG = 1
NN_part_trans_switch = 1
freqio   = 20
NN_TEM   = 0
NN_MODEL = 0

R0_in = RTOR + SHIF(NA1)
Bo_in(1) = BTOR
Rmin_in(1) = AMETR(NA1)
rhoscale = rho(NA1)
write(6, *) 'Calling qualikiz NN', nions, NN_MODEL

radial_loop: do jrho=1, NA1

    qx_in(1) = q_exp(jrho)
    rho_in(1) = RHO(jrho)
    x_in(1) = AMETR(jrho)/Rmin_in(1)
    rho_tg(jrho) = rho_in(1)
    Ro_in(1) = RTOR + SHIF(jrho)
    T0  = 1E3 *te(jrho)   ! eV

!thermal impurities

    Zi_in(1, 2) = max(1., ZIM1(jrho))
    Zi_in(1, 3) = ZIM2(jrho)
    Zi_in(1, 4) = ZIM3(jrho)

    if (Zi_in(1, 4) .ge. 1. .and. nions .eq. 2) then
        Zi_in(1, 3) = Zi_in(1, 4)
        Ai_in(1, 3) = Ai_in(1, 4)
        ni_m(3, :) = ni_m(4, :)
        ti_m(3, :) = ti_m(4, :)
    endif
    if (Zi_in(1, 3) .ge. 1. .and. nions .eq. 1) then
        Zi_in(1, 2) = Zi_in(1, 3)
        Ai_in(1, 2) = Ai_in(1, 3)
        ni_m(2, :) = ni_m(3, :)
        ti_m(2, :) = ti_m(3, :)
    endif

    Ai_in(1, nions + 1:) = 0.
    Zi_in(1, nions + 1:) = 0.
    Ati_in(1, nions + 1:) = 0.
    Ani_in(1, nions + 1:) = 0.

! Differentials

    jr_r = jrho+1
    jr_l = jrho-1
    if (jrho == 1) then
        jr_l = jrho
    else if (jrho == NA1) then
        jr_r = jrho
    endif
    dstep = 1./float(jr_r - jr_l)

    drmin  = dstep*(AMETR(jr_r) - AMETR(jr_l))
    drmaj  = dstep*(rmaj_exp(jr_r) - rmaj_exp(jr_l))
    drho   = dstep*(rho(jr_r) - rho(jr_l))
    delong = dstep*(ELON(jr_r) - ELON(jr_l))
    dtrian = dstep*(TRIA(jr_r) - TRIA(jr_l))
    dptot  = dstep*(ptot(jr_r) - ptot(jr_l))
    dte    = dstep*(TE(jr_r) - TE(jr_l))
    dne    = dstep*(NE(jr_r) - NE(jr_l))
    dq     = dstep*(q_exp(jr_r) - q_exp(jr_l))
    dvper  = dstep*(vper_m(jr_r) - vper_m(jr_l))
    do jion=1, nions
        dti(jion) = dstep*(ti_m(jion, jr_r) - ti_m(jion, jr_l))
        dni(jion) = dstep*(ni_m(jion, jr_r) - ni_m(jion, jr_l))
    enddo
    dv_r = dstep* &
        (vpar_m(jr_r)/(rmaj_exp(jr_r) + AMETR(jr_r)) - &
         vpar_m(jr_l)/(rmaj_exp(jr_l) + AMETR(jr_l)))
    dr = drmin/Rmin_in(1)    ! gradients w.r.t. minor radius even for s-alpha geometry
    drhodr = drho/drmin
    smag_in(1) = (x_in(1)/qx_in(1))*dq/dr        ! r/q dq/dr

! local field averages
    Nex_in(1) = NE(jrho)
    Tex_in(1) = TE(jrho)
    ZEFFX(1) = ZEF(jrho)
    epsilon(1) = Rmin_in(1) / Ro_in(1)
    Ate_in(1) = -dte*R0_in/(drmin*Tex_in(1))
    Ane_in(1) = -dne*R0_in/(drmin*Nex_in(1))
    do jion=1, nions
        ninorm_in(1, jion) = ni_m(jion, jrho)/Nex_in(1)
        Tix_in(1, jion) = ti_m(jion, jrho)
        Ati_in(1, jion) = -dti(jion)*R0_in/(drmin*ti_m(jion, jrho))
        Ani_in(1, jion) = -dni(jion)*R0_in/(drmin*ni_m(jion, jrho))
    enddo

! Restore quasi-neutrality via main ions
    ninorm_in(1, 1) = (1 - SUM(ninorm_in(1, 2:nions)*Zi_in(1, 2:nions)))/Zi_in(1, 1)
    Ani_in(1, 1) = (Ane_in(1) - SUM(ninorm_in(1, 2: nions)*Ani_in(1, 2: nions)* &
               Zi_in(1, 2: nions)))/(ninorm_in(1, 1)*Zi_in(1, 1))

! derived units for the plasma

    Bunit = 1E4*BTOR*drhodr*rho_in(1)/AMETR(jrho)  ! Miller geometry magnetic field unit
    cs00 = SQRT(e00*T0/(AMJ*mpp))                ! thermal velocity unit m/sec
    omega0 = e0*Bunit/(m0*c0)                    ! gyrofrequency unit 1/sec
    rhos00 = cs00/omega0                         ! gyroradius unit m

    vpar_shear_in = -rmaj_exp(jrho)*dv_r/(dr*cs00) ! From m/s to cm/s for vpar
    vpar_in = vpar_m(jrho)/cs00

! Local magnetic geometry

    rhostar2 = (rhos00/Rmin_in(1))**2
    drho_cs = drhodr**2 * Rmin_in(1) * rhostar2 * cs00

    alphax_in(1) = -(0.0040267/BTOR**2) * qx_in(1)**2 * R0_in * dptot/(dr*Rmin_in(1))
    if (NN_MODEL == 0) then
        smag_in = smag_in - 0.5*alphax_in
    endif
    cexb = AMETR(jrho)/qx_in(1)         ! r/(q)
    gamma_e_tg = cexb*dvper/(dr*cs00) ! Waltz-Miller definition
    mach_fac = sqrt(Tex_in(1)/AMJ)

    Machtor_in(1) = vpar_in*mach_fac
    Machpar_in(1) = vpar_in*mach_fac
    Autor_in(1)  = vpar_shear_in*R0_in/Rmin_in(1)*mach_fac
    Aupar_in(1)  = vpar_shear_in*R0_in/Rmin_in(1)*mach_fac
    gammaE_in(1) = gamma_e_tg   *R0_in/Rmin_in(1)*mach_fac

    if (FIRST_RUN) then
        SELECT CASE(NN_model)
        CASE(0)
            write(*, *) 'Using QLKNN-hyper as QuaLiKiz surrogate'
        CASE(1)
            write(*, *) 'Using QLKNN-HornNet as QuaLiKiz surrogate'
        CASE(-1)
            write(*, *) 'Using QLKNN-fullflux as QuaLiKiz surrogate. Warning! Developer only!'
        CASE(2)
            write(*, *) 'Using QLKNN-jetexp as QuaLiKiz surrogate'
            n_qlknn_invalid = 0
        CASE DEFAULT
            write(*, *) 'QLKNN NN_model out of bounds [-1, 2]. Stopping'
            stop
        END SELECT
    endif

! Allocate QLKNN in-out
    if(.not. allocated(qlknn_in)) then
        SELECT CASE(NN_model)
        CASE(0)
            ALLOCATE(qlknn_in(11, dimx))
            ALLOCATE(qlknn_out(dimx, 10))
        CASE(1)
            ALLOCATE(qlknn_in(11, dimx))
            ALLOCATE(qlknn_out(dimx, 4))
        CASE(-1)
            ALLOCATE(qlknn_in(9, dimx))
            ALLOCATE(qlknn_out(dimx, 3))
        CASE(2)
            ALLOCATE(qlknn_in(15, dimx))
            ALLOCATE(qlknn_out(dimx, 13))
        END SELECT
    endif

    Lambe(:) = 15.2 - 0.5*LOG(0.1*Nex_in(:)) + LOG(Tex_in(:)) !Coulomb constant and collisionality. Wesson 2nd edition p661-66

    Nue(:) = 1./(1.09d-3) *Zeffx(:)*Nex_in(:)*Lambe(:)/(Tex_in(:))**1.5*collmultin
    Nustar(:) = Nue(:)*qx_in(:)*Ro_in(:)/epsilon(:)**1.5/(SQRT(Tex_in(:)*1.d3*e00/(z_electronmass*mpp)))

    if (first_run) then
        SELECT CASE(NN_model)
        CASE(0)
! data location set as local softlink by rjettov
            call load_qlknn_hyper_nets_from_disk(trim(qlknn_sets_dir) // &
                'qlknn-hyper-namelists', INT(verbose_in, lli))
            if (verbose_in >= 1) then
                write(*, *) 'Hyper neural networks successfully loaded!'
            endif

            call default_qlknn_hyper_options(qlknn_opts)

            qlknn_opts%use_etg=NN_ETG
            qlknn_opts%use_tem=NN_TEM
            qlknn_opts%constrain_outputs=.FALSE.

            SELECT CASE(NN_part_trans_switch)
            CASE(1)
                qlknn_opts%use_ion_diffusivity_networks = .FALSE.
                qlknn_opts%use_effective_diffusivity = .TRUE.
            CASE(2)
                qlknn_opts%use_ion_diffusivity_networks = .FALSE.
                qlknn_opts%use_effective_diffusivity = .FALSE.
            CASE(3)
                qlknn_opts%use_ion_diffusivity_networks = .TRUE.
                qlknn_opts%use_effective_diffusivity = .FALSE.
            CASE DEFAULT
                write(*, *) 'QLKNN NN_part_trans_switch out of bounds (1-3)! Stopping'
                stop
            END SELECT

        CASE(1)
! data location set as local softlink by rjettov
            call load_hornnet_nets_from_disk(trim(qlknn_sets_dir) // &
                'qlknn-hornnet-namelists', INT(verbose_in, lli))
            if (verbose_in >= 1) then
                write(*, *) 'Hornnet neural networks successfully loaded!'
            endif

            call default_qlknn_hornnet_options(qlknn_opts)

            qlknn_opts%use_etg=NN_ETG
            qlknn_opts%use_tem=NN_TEM
            qlknn_opts%constrain_outputs=.FALSE.

            if (NN_part_trans_switch /= 1) then
                write(*, *) 'QLKNN NN_part_trans_switch should be 1 for MegaHornNet! Stopping'
                stop
            endif

        CASE(-1)
            call load_fullflux_nets_from_disk(trim(qlknn_sets_dir) // &
                'qlknn-fullflux-namelists', INT(verbose_in, lli))
            if (verbose_in >= 1) then
                write(*, *) 'Fullflux neural networks successfully loaded!'
            endif

            call default_qlknn_fullflux_options(qlknn_opts)

            if (.not. NN_ETG) then
                write(*, *) 'Fullflux nets are always running with ETG on'
                stop
            endif
            if (.not. NN_TEM) then
                write(*, *) 'Fullflux nets are always running with TEM on'
                stop
            endif

            qlknn_opts%constrain_outputs=.FALSE.

            if (NN_part_trans_switch /= 1) then
                write(*, *) 'Fullflux nets are always running with part_trans_switch == 1'
                stop
            endif

        CASE(2)
            call load_jetexp_nets_from_disk(trim(qlknn_sets_dir) // &
                'qlknn-jetexp-namelists', qlknn_members, INT(verbose_in, lli))
            if (verbose_in >= 1) then
                write(*, *) 'Jetexp neural networks successfully loaded!'
            endif

            call default_qlknn_jetexp_options(qlknn_opts)

            qlknn_opts%constrain_outputs=.FALSE.

        END SELECT

        first_run = .FALSE.
    endif ! if first_run

    SELECT CASE(NN_model)
    CASE(0, 1, -1)
        if (NN_model /= 1) then
            epsilonNN = nets%a / nets%R_0 ! epsilon used to train 10D NN
        else
            epsilonNN = 1. / 3. ! epsilon used to train hornnet NN
        endif

        qlknn_in(1, :) = Zeffx(:)
        qlknn_in(2, :) = Ati_in(:, 1)
        qlknn_in(3, :) = Ate_in(:)
        qlknn_in(4, :) = Ane_in(:)
        qlknn_in(5, :) = qx_in(:)
        qlknn_in(6, :) = smag_in(:)
        qlknn_in(7, :) = x_in(:)
        qlknn_in(8, :) = Tix_in(:, 1) / Tex_in(:)
        qlknn_in(9, :) = LOG10(Nustar(:))

! Set gammaE normalisation quantities
        IF( .not. allocated(qlknn_norms%A1) ) ALLOCATE(qlknn_norms%A1(dimx))
        qlknn_norms%A1 = Ai_in(:, 1)
        qlknn_norms%R0 = R0_in
        qlknn_norms%a  = Rmin_in(1)  ! Is it aref?

    CASE(2) ! QLKNN-jetexp
        qlknn_in(nets%Ane_ind, :)     = Ane_in(:)
        qlknn_in(nets%Ate_ind, :)     = Ate_in(:)
        qlknn_in(nets%Autor_ind, :)   = Autor_in(:)
        qlknn_in(nets%Machtor_ind, :) = Machtor_in(:)
        qlknn_in(nets%x_ind, :)       = x_in(:)
        qlknn_in(nets%Zeff_ind, :)    = Zeffx(:)
        qlknn_in(nets%gammaE_ind, :)  = gammaE_in(:) ! gammaE QLK [cref/R]
        qlknn_in(nets%q_ind, :)       = qx_in(:)
        qlknn_in(nets%smag_ind, :)    = smag_in(:)
        qlknn_in(nets%alpha_ind, :)   = alphax_in(:)
        qlknn_in(nets%Ati0_ind, :)    = Ati_in(:, 1)

        if (nspec_max-1 > 1) then
            qlknn_in(nets%Ani1_ind, :) = Ani_in(:, 2)
            qlknn_in(nets%normni1_ind, :) = ninorm_in(:, 2)
        else
            qlknn_in(nets%Ani1_ind, :) = 0
            qlknn_in(nets%normni1_ind, :) = 0
        endif

        qlknn_in(nets%Ti_Te0_ind, :) = Tix_in(:, 1) / Tex_in
        qlknn_in(nets%logNustar_ind, :) = LOG10(Nustar)
    END SELECT

    SELECT CASE(NN_model)
    CASE(0)
! Modify rotation quantities to be consistent with rot_flag choice
        SELECT CASE(rot_flag_in)
        CASE(0) ! no rotation
            qlknn_opts%apply_victor_rule = .FALSE.
        CASE(1) ! victor rule. Full radius only
            qlknn_opts%apply_victor_rule = .TRUE.
            qlknn_in(10, :) = gammaE_in ! gammaE QLK [cref/R]
            qlknn_in(11, :) = Tex_in
        CASE(2) ! victor rule. Half radius
            qlknn_opts%apply_victor_rule = .TRUE.
            ihi=dimx
            do i=1, dimx ! define filterprof
                if(rho_in(i)/rhoscale < 0.4) then
                    filterprof(i)=0.
                    ilow=i
                elseif (rho_in(i)/rhoscale > 0.6) then
                    filterprof(i)=1.
                    if(i<ihi) ihi=i
                endif
            enddo
            filterprof(ilow:ihi) = (/( REAL((i-ilow))/REAL((ihi-ilow))  , i=ilow, ihi)/)
            qlknn_in(10, :) = gammaE_in*filterprof ! gammaE QLK [cref/R]
            qlknn_in(11, :) = Tex_in
        CASE DEFAULT
            stop 'QLK rot_flag out of bounds (0-2)! Stopping'
        END SELECT

        call evaluate_QLKNN_10D(qlknn_in, nets, qlknn_out, INT(verbose_in-1, lli), qlknn_opts, qlknn_norms)

        call evaluate_QLKNN_10D(qlknn_in, nets, qlknn_out, INT(verbose_in-1, lli), qlknn_opts, qlknn_norms)

        eef_GB_out(:) = qlknn_out(:, 1)
        epf_GB_out(:) = qlknn_out(:, 4)

        do j=1, nspec_max-1
            ief_GB_out(:, j) = qlknn_out(:, 3) !provide all ions with same GB heat flux. Important for multiple-isotope simulations
        enddo

    CASE(1)
        if (rot_flag_in /= 0) then
            write(*, *) 'MegaHornNet cannot run with rotation'
            stop
        endif

        call evaluate_hornnet_constants(qlknn_in, blocks, qlknn_hornnet_constants, INT(verbose_in-1, lli), qlknn_opts, qlknn_norms)
        call hornnet_flux_from_constants(qlknn_in, blocks, qlknn_hornnet_constants, qlknn_out, INT(verbose_in-1, lli), qlknn_opts, qlknn_norms)
        if (verbose_in >= 2) then
            write(*, *) 'MegaHornNet evaluated'
        endif

        eef_GB_out(:) = qlknn_out(:, 1)
        epf_GB_out(:) = qlknn_out(:, 4)

        do j=1, nspec_max-1
            ief_GB_out(:, j) = qlknn_out(:, 3) ! same GB heat flux to all ion species. Important for multiple-isotope simulations
        enddo

    CASE(-1) ! QLKNN-fullflux
! Modify rotation quantities to be consistent with rot_flag choice
        if (rot_flag_in /= 0) then
            write(*, *) 'Fullflux nets cannot run with rotation'
            stop
        endif

        call evaluate_fullflux_net(qlknn_in, nets, qlknn_out, INT(verbose_in-1, lli), qlknn_opts)

        eef_GB_out(:) = qlknn_out(:, 1)
        epf_GB_out(:) = qlknn_out(:, 4)

        do j=1, nspec_max-1
            ief_GB_out(:, j) = qlknn_out(:, 2) !provide all ions with same GB heat flux. Important for multiple-isotope simulations
        enddo

    CASE(2)
! Modify rotation quantities to be consistent with rot_flag choice
        SELECT CASE(rot_flag_in)
        CASE(0) ! no rotation
            qlknn_in(nets%Machtor_ind, :) = Machtor_in*0.
            qlknn_in(nets%Autor_ind, :) = Autor_in*0.
            qlknn_in(nets%gammaE_ind, :) = gammaE_in*0. ! gammaE QLK [cref/R]
        CASE(1) ! apply on full radius, effectively do nothing
        CASE(2) ! apply only on half radius
            ihi=dimx
            do i=1, dimx ! define filterprof
                if(rho_in(i)/rhoscale < 0.4) then
                    filterprof(i)=0.
                    ilow=i
                elseif (rho_in(i)/rhoscale > 0.6) then
                    filterprof(i)=1.
                    if(i<ihi) ihi=i
                endif
            enddo
            filterprof(ilow:ihi) = (/( REAL((i-ilow))/REAL((ihi-ilow))  , i=ilow, ihi)/)
! newest 14D networks trained with inherent rot_flag=2, no need to modify inputs

        CASE DEFAULT
            stop 'QLK rot_flag out of bounds (0-2)! Stopping'
        END SELECT

        call evaluate_jetexp_net(qlknn_in, nets, qlknn_members, qlknn_out, qlknn_eb, INT(verbose_in-1, lli), qlknn_opts, qlknn_validity=qlknn_validity)

        qlknn_validity_mask = .FALSE.  ! Only check validity for values actually being used by JETTO
        qlknn_validity_mask(1) = .TRUE.
        qlknn_validity_mask(2) = .TRUE.
        qlknn_validity_mask(3) = .TRUE.
        SELECT CASE(NN_part_trans_switch)
        CASE(1) ! NN part trans output is electron particle flux
            qlknn_validity_mask(4) = .TRUE.
        CASE(2) ! Gammae and De output
            qlknn_validity_mask(4) = .TRUE.
            qlknn_validity_mask(7) = .TRUE.
        CASE(3) ! Gammai and Di output
            qlknn_validity_mask(4) = .TRUE.
            qlknn_validity_mask(10) = .TRUE.
        CASE(4) ! Ion particle flux (main ion only)
            qlknn_validity_mask(5) = .TRUE.
        CASE(5) ! De and Ve output
            qlknn_validity_mask(7) = .TRUE.
            qlknn_validity_mask(8) = .TRUE.
            qlknn_validity_mask(9) = .TRUE.
        CASE(6) ! Di and Vi output (main ion only)
            qlknn_validity_mask(10) = .TRUE.
            qlknn_validity_mask(11) = .TRUE.
            qlknn_validity_mask(12) = .TRUE.
            qlknn_validity_mask(13) = .TRUE.
        END SELECT

        j = 0
        k = 0
        qlknn_rho_validity = .FALSE.
        do i=1, dimx
! Only flip switch based on quantities used, covered by mask
            WHERE (.not. qlknn_validity_mask) qlknn_validity(i, :) = .TRUE.
            qlknn_rho_validity(i) = ALL(qlknn_validity(i, :))
! Only flip switch based on validity inside QLK boundaries

            if(rho(i)/rhoscale >= rhomin .AND. rho(i)/rhoscale <= rhomax) then
                j = j + 1
                if( .not. qlknn_rho_validity(i) ) then
                    k = k + 1
                endif
            endif
        enddo
        qlknn_out_of_bounds = (k * 2 >= j)

        if(.not. qlknn_out_of_bounds) then
            eef_GB_out(:) = qlknn_out(:, 1)
            epf_GB_out(:) = qlknn_out(:, 4)

            do j=1, nspec_max-1
                ief_GB_out(:, j) = 0.
               ! Provide all ions with Z=1 with same GB heat flux. Important for multiple-isotope simulations
               ! Impurity heat flux negligible due to ni in GB descaling
                WHERE (Zi_in(:, j) == 1.0)
                    ief_GB_out(:, j) = qlknn_out(:, 3)
                ENDWHERE
             enddo

        else
            write(*, *) 'QuaLiKiz neural network is extrapolating! Switching to original QuaLiKiz model...'
        endif

! Counter to keep track of invalid time slices in run (maybe keep track per rho later?)
        do i=1, dimx
            do j=1, size(qlknn_validity, 2)
                if( .not. qlknn_validity(i, j) ) then
                    n_qlknn_invalid(i, j) = n_qlknn_invalid(i, j) + 1
                endif
            enddo
        enddo
        if( (freqio == 0) .AND. (verbose_in >= 1) ) then
            open(unit=702, file="nsteps_NN_invalid.qlk", action="write", status="replace")
            write(fmtnets, '(A, I0, A)') '(', size(n_qlknn_invalid, 1), 'I6)'
            do i=1, dimx
                write(702, fmtnets) (n_qlknn_invalid(i, j), j=1, size(n_qlknn_invalid, 2))
            enddo
            close(702)
        endif

    END SELECT

    ql_fac = drhodr**2 * Rmin_in(1) * rhostar2 * cs00

! Chii

    chii(jrho)   = ql_fac * ief_gb_out(1, 1)/(1e-4 + Rmin_in(1)/R0_in * abs(Ati_in(1, 1)))
    chie(jrho)   = ql_fac * eef_gb_out(1)/(1e-4 + Rmin_in(1)/R0_in * abs(Ate_in(1)))
    pfluxi(jrho) = ql_fac * epf_gb_out(1)

enddo radial_loop

call qinterp(rho_tg, chii  , nradial, RHO(1:jna), chii_m(1:jna)  , jna)
call qinterp(rho_tg, chie  , nradial, RHO(1:jna), chie_m(1:jna)  , jna)
call qinterp(rho_tg, pfluxi, nradial, RHO(1:jna), pfluxi_m(1:jna), jna)

chii_m  (1:2) = chii_m(3)
chie_m  (1:2) = chie_m(3)
pfluxi_m(1:2) = pfluxi_m(3)

do j=1, jna
    CHI(j) = chii_m(j)/gradrhosq_exp(j) ! \chi_i, m^2/s
    CHE(j) = chie_m(j)/gradrhosq_exp(j) ! \chi_e, m^2/s
    VIN(j) = pfluxi_m(j)/AMETR(NA1)/gradrhosq_exp(j) ! D flux
    VIN(j) = min( 20., max(-20., VIN(j)) ) ! D flux
enddo

do j=jna-2, jna-1
    CHI(j) = CHI(jna-3)
    CHE(j) = CHE(jna-3)
    VIN(j) = VIN(jna-3)
enddo

if (jna .lt. NA1) then
    do j=jna-1, NA1
        CHI(j) = 0.d0
        CHE(j) = 0.d0
        VIN(j) = 0.d0
    enddo
endif

return
END subroutine qlknn_serial
