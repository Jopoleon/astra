subroutine qlk_interf(jproc, dims_in, scal_in, prof_in, prof_out)

USE mod_qualikiz, only: qualikiz
USE kind, only: qlk_output_meth_0, qlk_output_meth_0_sep_0, &
    qlk_primi_meth_0, qlk_sizes, qlk_in_regular, qlk_in_newt
USE mod_make_io, only: zeroout_qlk_in_newt, zeroout_qlk_in_regular, &
    allocate_qlk_in_regular, allocate_qlk_in_newt
USE mpi
USE nanfilter

implicit none

type(qlk_sizes)      :: sizes
type(qlk_in_regular) :: in_regular
type(qlk_in_newt)    :: in_newt

type(qlk_output_meth_0)       :: output_meth_0
type(qlk_output_meth_0_sep_0) :: output_meth_0_sep_0_SI , output_meth_0_sep_0_GB
type(qlk_primi_meth_0)        :: primi_meth_0

integer, parameter :: ntheta=64, numecoefs=13, numicoefs=7, dimx=1, dimn=16, numsols=3, phys_meth=0, nspec_max=7

double precision, parameter :: &
    e0  = 4.8032E-10, &       ! elementary charge (statcoulombs)
    e00 = 1.6020e-19, &       ! elementary charge (C)
    c0  = 2.9979E+10, &       ! speed of light (cm/sec)
    mp  = 1.6726E-24, &       ! proton mass (g)
    mpp = 1.6726E-27          ! proton mass (kg)

integer, intent(in) :: jproc, dims_in(*)
double precision, intent(in) :: scal_in(*)
double precision, intent(in), dimension(dims_in(2), dims_in(4)) :: prof_in
double precision, intent(out), dimension(dims_in(3), dims_in(1)) :: prof_out

!-----------------------------------------

integer :: jr1, jr2, n_inputs, n_outputs, nrho, ns_in, chunk
integer :: simple_mpi_only_in, maxpts_in, maxruns_in, runcounter_in
integer :: nions, coll_flag_in, rot_flag_in, verbose_in, el_type_in, &
     integration_routine_in, separateflux_in
integer :: mpi_ierr, nproc, myrank, i_mpic
integer, dimension(dimx, nspec_max-1) :: ion_type_in

double precision :: a0_m, T0, m0, rmin_tg, drho_cs, drho_nt, nt_cs
double precision :: relacc1_in, relacc2_in, absacc1_in, absacc2_in, R0_in, &
    ETGmultin, collmultin, timeout_in, rhomin, rhomax, rhoscale
double precision, dimension(dimn) :: kthetarhos_in
double precision, dimension(dimx) :: x_in, rho_in, Ro_in, Rmin_in, Bo_in, &
    qx_in, smag_in, alphax_in, Tex_in, Nex_in, Ate_in, Ane_in, anise_in, &
    danisedr_in, Machtor_in, Autor_in, Machpar_in, Aupar_in, gammaE_in, &
    epf_GB_out, eef_GB_out
double precision, dimension(dimx, nspec_max-1) :: Tix_in, ninorm_in, &
    Ati_in, Ani_in, anis_in, danisdr_in, Ai_in, Zi_in, ipf_GB_out, ief_GB_out
double precision, dimension(dimx, dimn, numsols) :: gam_GB_out, ome_GB_out
double precision, dimension(dimx, nspec_max - 1, numicoefs) :: cftrans_out

! Old solution for non-reset runs
double precision, DIMENSION(:, :, :), ALLOCATABLE :: oldrsol, oldisol, oldrfdsol, oldifdsol

!-----------------------------------------
LOGICAL :: exist1, exist2, exist3, exist4, exist5 !used for checking for existence of files

integer, parameter :: unit_runc=600, unit_rsol=610, unit_isol=620, unit_rfd=630, unit_ifd=640, myunit=700, unit_input=500
integer :: i, j, k, jr, jion

double precision, allocatable, dimension(:) :: drmin, drmaj, drho, dte, dne, dq, dptot, &
    dvper, drhodr, dstep, dr, dv_r
double precision :: BTOR, RTOR, ABC, AMJ, AIM1, AIM2, AIM3, ZMJ
double precision :: Bunit, cs00, rhos00, omega0, rhostar2
! Shifted cicle geometry inputs
double precision :: gamma_e, mach_fac, ql_fac

double precision, allocatable, dimension(:) :: chie, chii, exchi, pfluxi, rho, &
    te, ne, vpar, vper, ametr, rmaj, q_saf
double precision, allocatable, dimension(:, :) :: dti, dni, ni, ti, zi
double precision :: vpar_in, vpar_shear_in, cexb

COMPLEX(kind=DBL), DIMENSION(:, :, :), ALLOCATABLE :: oldsol_in, oldfdsol_in
complex(kind=DBL), dimension(dimx, dimn, numsols) :: sol_out, fdsol_out

CHARACTER(len=20) :: fmtn
character(len=80) :: fname1, prim_dir

!---------------------------------

save i_mpic
data i_mpic /0/

verbose_in = 0

if (i_mpic == 0) then
    CALL mpi_init(mpi_ierr)
    if (mpi_ierr /= 0) then
        write(6, *) 'GIt qlk_interface mpi_init error=', mpi_ierr
    endif
    CALL mpi_comm_size(mpi_comm_world, nproc,  mpi_ierr)
    CALL mpi_comm_rank(mpi_comm_world, myrank, mpi_ierr)
    i_mpic = 1
endif

verbose_in = 0

chunk     = dims_in(1)
n_inputs  = dims_in(2)
n_outputs = dims_in(3)
nrho      = dims_in(4)
ns_in     = dims_in(5)

jr1 = (jproc - 1)*chunk + 1
jr2 = jproc*chunk

BTOR = scal_in(1)
RTOR = scal_in(2)
ABC  = scal_in(3)
AMJ  = scal_in(4)
AIM1 = scal_in(5)
AIM2 = scal_in(6)
AIM3 = scal_in(7)
ZMJ  = scal_in(8)

nions = ns_in - 1

allocate( chie(chunk), chii(chunk), exchi(chunk), pfluxi(chunk), rho(chunk), &
    te(chunk), ne(chunk), vpar(chunk), vper(chunk), &
    ametr(chunk), rmaj(chunk), q_saf(chunk) )
allocate( drmin(chunk), drmaj(chunk), drho(chunk), &
    dr(chunk), dne(chunk), dte(chunk), dq(chunk), dptot(chunk), &
    dvper(chunk), dv_r(chunk), drhodr(chunk) )
allocate( dti(4, chunk), dni(4, chunk), ni(4, chunk), ti(4, chunk), zi(4, chunk) )

! Receive TGLF input scalars and profiles from parent
rho      = prof_in( 1, jr1:jr2)
ametr    = prof_in( 2, jr1:jr2)
rmaj     = prof_in( 3, jr1:jr2)
q_saf    = prof_in( 4, jr1:jr2)
ne       = prof_in( 5, jr1:jr2)
te       = prof_in( 6, jr1:jr2)
vpar     = prof_in( 7, jr1:jr2)
vper     = prof_in( 8, jr1:jr2)
ti(1, :) = prof_in( 9, jr1:jr2)
ti(2, :) = prof_in(10, jr1:jr2)
ti(3, :) = prof_in(11, jr1:jr2)
ti(4, :) = prof_in(12, jr1:jr2)
ni(1, :) = prof_in(13, jr1:jr2)
ni(2, :) = prof_in(14, jr1:jr2)
ni(3, :) = prof_in(15, jr1:jr2)
ni(4, :) = prof_in(16, jr1:jr2)
zi(1, :) = prof_in(17, jr1:jr2)
zi(2, :) = prof_in(18, jr1:jr2)
zi(3, :) = prof_in(19, jr1:jr2)
zi(4, :) = prof_in(20, jr1:jr2)
drmin      = prof_in(21, jr1:jr2)
drmaj      = prof_in(22, jr1:jr2)
drho       = prof_in(23, jr1:jr2)
dptot      = prof_in(24, jr1:jr2)
dte        = prof_in(25, jr1:jr2)
dne        = prof_in(26, jr1:jr2)
dq         = prof_in(27, jr1:jr2)
dvper      = prof_in(28, jr1:jr2)
dv_r       = prof_in(29, jr1:jr2)
dr         = prof_in(30, jr1:jr2)
drhodr     = prof_in(31, jr1:jr2)
dti(1, :)  = prof_in(32, jr1:jr2)
dti(2, :)  = prof_in(33, jr1:jr2)
dti(3, :)  = prof_in(34, jr1:jr2)
dti(4, :)  = prof_in(35, jr1:jr2)
dni(1, :)  = prof_in(36, jr1:jr2)
dni(2, :)  = prof_in(37, jr1:jr2)
dni(3, :)  = prof_in(38, jr1:jr2)
dni(4, :)  = prof_in(39, jr1:jr2)

a0_m = ABC
R0_in    = prof_in(3, nrho) ! nrho has to be NA1?
rhoscale = prof_in(1, nrho) ! ??
m0 = AMJ*mp          ! Ref. mass = D ion mass [g]

! Electrons and main ions
Ai_in(1, 1) = AMJ
Ai_in(1, 2) = AIM1
Ai_in(1, 3) = AIM2
Ai_in(1, 4) = AIM3
Zi_in(1, 1) = ZMJ

vpar_in = 0.
vpar_shear_in = 0.

kthetarhos_in(1:8) = (/0.1, 0.175, 0.25, 0.325, 0.4, 0.5, 0.7, 1.0/)
if (dimn == 16) then
    kthetarhos_in(9:16) = (/1.8,   3.0,  9.0,  15.0, 21.0, 27.0, 36.0, 45.0/)
endif

coll_flag_in = 1
el_type_in   = 1
separateflux_in    = 0
simple_mpi_only_in = 1
maxpts_in  = 50000000 ! 500000 default
maxruns_in = 10

Bo_in(1) = BTOR
Rmin_in(1) = a0_m
WRITE(fmtn, '(A, I0, A)') '(', dimn, 'G15.7)'

radial_loop: do jr=1, chunk
   
    write(fname1, '(A, i0)') 'qlkzin_' , jr1 - 1 + jr
    qx_in(1)  = q_saf(jr)
    rho_in(1) = rho(jr)
    x_in(1)   = ametr(jr)/Rmin_in(1)
    Ro_in(1)  = rmaj(jr)
    T0 = 1E3*te(jr)   ! temperature scale used by GYRO
    smag_in(1) = (x_in(1)/qx_in(1))*dq(jr)/dr(jr)        ! r/q dq/dr 

!thermal impurities

    Zi_in(1, 2) = max(1., zi(2, jr))
    Zi_in(1, 3) = zi(3, jr)
    Zi_in(1, 4) = zi(4, jr)

    if (Zi_in(1, 4) >= 1. .and. nions == 2) then
        Zi_in(1, 3) = Zi_in(1, 4)
        Ai_in(1, 3) = Ai_in(1, 4)
        ni(3, :) = ni(4, :)
        ti(3, :) = ti(4, :)
    endif
    if (Zi_in(1, 3) >= 1. .and. nions == 1) then
        Zi_in(1, 2) = Zi_in(1, 3)
        Ai_in(1, 2) = Ai_in(1, 3)
        ni(2, :) = ni(3, :)
        ti(2, :) = ti(3, :)
    endif

    Ai_in(1, nions + 1:) = 0.
    Zi_in(1, nions + 1:) = 0.
    Ati_in(1, nions + 1:) = 0.
    Ani_in(1, nions + 1:) = 0.

! local field averages
    Nex_in(1) = ne(jr)
    Tex_in(1) = te(jr)
    Ate_in(1) = -dte(jr)*R0_in/(drmin(jr)*Tex_in(1))
    Ane_in(1) = -dne(jr)*R0_in/(drmin(jr)*Nex_in(1))
    anise_in(1) = 1.0
    danisedr_in(1) = 0.0
    do jion=1, nions  ! deut, impurities
        ninorm_in(1, jion) = ni(jion, jr)/Nex_in(1)
        Tix_in(1, jion) = ti(jion, jr)
        Ati_in(1, jion) = -dti(jion, jr)*R0_in/(drmin(jr)*ti(jion, jr))
        Ani_in(1, jion) = -dni(jion, jr)*R0_in/(drmin(jr)*ni(jion, jr))
        ion_type_in(1, jion) = 1
        anis_in    (1, jion) = 1.
        danisdr_in (1, jion) = 0.
    enddo

! Restore quasi-neutrality via main ions
    ninorm_in(1, 1) = (1 - SUM(ninorm_in(1, 2:nions)*Zi_in(1, 2:nions)))/Zi_in(1, 1)
    Ani_in(1, 1) = (Ane_in(1) - SUM(ninorm_in(1, 2: nions)*Ani_in(1, 2: nions)* &
               Zi_in(1, 2: nions)))/(ninorm_in(1, 1)*Zi_in(1, 1))

! derived units for the plasma

    Bunit = 1E4*BTOR*drhodr(jr)*rho_in(1)/ametr(jr)  ! Miller geometry magnetic field unit
    cs00 = SQRT(e00*T0/(AMJ*mpp))  ! thermal velocity unit m/sec
    omega0 = e0*Bunit/(m0*c0)      ! gyrofrequency unit 1/sec
    rhos00 = cs00/omega0           ! gyroradius unit m

    vpar_shear_in = -rmaj(jr)*dv_r(jr)/(dr(jr)*cs00)  !From m/s to cm/s for vpar
    vpar_in = vpar(jr)/cs00

! local magnetic geometry

    rhostar2 = (rhos00/Rmin_in(1))**2
    drho_cs = drhodr(jr)**2 * Rmin_in(1) * rhostar2 * cs00

! 1.e16 comes from cgs to SI for ptot
    alphax_in(1) = -(0.0040267/BTOR**2) * qx_in(1)**2 * R0_in * dptot(jr)/(dr(jr)*Rmin_in(1))
    cexb = ametr(jr)/qx_in(1)            ! r/(q)
    gamma_e = cexb*dvper(jr)/(dr(jr)*cs00) ! Waltz-Miller definition
    mach_fac = sqrt(Tex_in(1)/AMJ)

    Machtor_in(1) = vpar_in*mach_fac
    Machpar_in(1) = vpar_in*mach_fac
    Autor_in(1)  = vpar_shear_in*R0_in/Rmin_in(1)*mach_fac
    Aupar_in(1)  = vpar_shear_in*R0_in/Rmin_in(1)*mach_fac
    gammaE_in(1) = gamma_e*R0_in/Rmin_in(1)*mach_fac

    absacc2_in = 0.
    absacc1_in = 0.
    relacc1_in = 1.e-3
    relacc2_in = 0.02
    timeout_in = 20.
    ETGmultin  = 1.
    collmultin = 1.  ! git / Originally: 1. ; 0.1 suggested by Jon
    rot_flag_in = 0  ! git / Originally: 2; 0 suggested by Pierre
    integration_routine_in = 1 ! 0 for NAG routines, 1 for Cubature

    if (verbose_in > 0) then
        print*, 'aNe', Ane_in
        print*, 'aNi', Ani_in(1, 1:nions)
        print*, 'ate', Ate_in
        print*, 'ati', Ati_in(1, 1:nions)
        print*, 'Ai', Ai_in(1, 1:nions)
        print*, 'Ni/Ne', ninorm_in(1, 1:nions)
        print*, 'rhos', R0_in, Ro_in(1), x_in(1), Rmin_in(1)
        print*, 'check: ', m0, Ati_in(1, 1), gammaE_in, Autor_in, Machtor_in, alphax_in(1)
    endif

    INQUIRE(file="qualikiz/"//trim(fname1)//"/runcounter.dat", EXIST=exist1)
    INQUIRE(file=trim(prim_dir)//'/rsol.dat', EXIST=exist2)
    INQUIRE(file=trim(prim_dir)//'/isol.dat', EXIST=exist3)
    INQUIRE(file=trim(prim_dir)//'/rfdsol.dat', EXIST=exist4)
    INQUIRE(file=trim(prim_dir)//'/ifdsol.dat', EXIST=exist5)

    IF ( exist1 .AND. exist2 .AND. exist3 .AND. exist4 .AND. exist5 ) THEN
        OPEN(unit=unit_runc, file="qualikiz/"//trim(fname1)//"/runcounter.dat", status="old", action="read")
        READ(unit_runc, *) runcounter_in
        CLOSE(unit_runc)
    ELSE
        prim_dir = 'qualikiz/'//trim(fname1)//'/output/primitive'
        call system('mkdir -p ' // trim(prim_dir))
        runcounter_in = 0 !First run
    ENDIF

    IF (runcounter_in >= maxruns_in) THEN !Reset if we're at our maximum number of runs
        runcounter_in = 0
    ENDIF
!    runcounter_in = 0 ! GIT force calculation from scratch

    IF (runcounter_in == 0) THEN !load old rsol and isol if we're not doing a reset run
        print *,'Qualikiz from scratch'
    ELSE
        ALLOCATE( oldrsol  (dimx, dimn, numsols) )
        ALLOCATE( oldisol  (dimx, dimn, numsols) )
        ALLOCATE( oldrfdsol(dimx, dimn, numsols) )
        ALLOCATE( oldifdsol(dimx, dimn, numsols) )
        if (.not. ALLOCATED(oldsol_in)) then
            ALLOCATE( oldsol_in   (dimx, dimn, numsols) )
            ALLOCATE( oldfdsol_in (dimx, dimn, numsols) )
        endif
        OPEN(unit=unit_rsol, file=trim(prim_dir)//'/rsol.dat', action="read", status="old")
        READ(unit_rsol, fmtn) (((oldrsol(i, j, k), j=1, dimn), i=1, dimx), k=1, numsols)
        CLOSE(unit_rsol)

        OPEN(unit=unit_isol, file=trim(prim_dir)//'/isol.dat', action="read", status="old")
        READ(unit_isol, fmtn) (((oldisol(i, j, k), j=1, dimn), i=1, dimx), k=1, numsols)
        CLOSE(unit_isol)

        OPEN(unit=unit_rfd, file=trim(prim_dir)//'/rfdsol.dat', action="read", status="old")
        READ(unit_rfd, fmtn) (((oldrfdsol(i, j, k), j=1, dimn), i=1, dimx), k=1, numsols)
        CLOSE(unit_rfd)

        OPEN(unit=unit_ifd, file=trim(prim_dir)//'/ifdsol.dat', action="read", status="old")
        READ(unit_ifd, fmtn) (((oldifdsol(i, j, k), j=1, dimn), i=1, dimx), k=1, numsols)
        CLOSE(unit_ifd)

        oldsol_in   = CMPLX(oldrsol  , oldisol)
        oldfdsol_in = CMPLX(oldrfdsol, oldifdsol)

        DEALLOCATE( oldrsol )
        DEALLOCATE( oldisol )
        DEALLOCATE( oldrfdsol )
        DEALLOCATE( oldifdsol )
        print *, 'Call qualikiz from old solution'
    ENDIF

    sizes%dimn = dimn
    sizes%dimx = dimx
    sizes%nions = nions
    sizes%numsols = numsols

! Prepare and clean input structures
  ! Read input and initialize all arrays

    if (.not. ALLOCATED(in_regular%Ane)) then
        CALL allocate_qlk_in_regular(dimx, dimn, nions, in_regular)
        CALL allocate_qlk_in_newt(dimx, dimn, numsols, in_newt)
        CALL zeroout_qlk_in_regular(in_regular)
        CALL zeroout_qlk_in_newt(in_newt)
    endif

    in_regular%Ane(:) = Ane_in(1:dimx)
    in_regular%Ate(:) = Ate_in(1:dimx)
    in_regular%Aupar(:) = Aupar_in(1:dimx)
    in_regular%Autor(:) = Autor_in(1:dimx)
    in_regular%Machpar(:) = Machpar_in(1:dimx)
    in_regular%Machtor(:) = Machtor_in(1:dimx)
    in_regular%x(:) = x_in(1:dimx)
    in_regular%Bo(:) = Bo_in(1:dimx)
    in_regular%gammaE(:) = gammaE_in(1:dimx)
    in_regular%ne(:) = Nex_in(1:dimx)
    in_regular%q(:) = qx_in(1:dimx)
    in_regular%Ro(:) = Ro_in(1:dimx)
    in_regular%Rmin(:) = Rmin_in(1:dimx)
    in_regular%smag(:) = smag_in(1:dimx)
    in_regular%Te(:) = Tex_in(1:dimx)
    in_regular%alpha(:) = alphax_in(1:dimx)
    in_regular%rho(:) = rho_in(1:dimx)/rhoscale

    in_regular%anise(:) = anise_in(1:dimx)
    in_regular%danisedr(:) = danisedr_in(1:dimx)

    in_regular%kthetarhos(:) = kthetarhos_in(1:dimn)

    in_regular%Ai(:,:) = Ai_in(1:dimx, 1:nions)
    in_regular%Ani(:,:) = Ani_in(1:dimx, 1:nions)
    in_regular%Ati(:,:) = Ati_in(1:dimx, 1:nions)
    in_regular%Zi(:,:) = Zi_in(1:dimx, 1:nions)
    in_regular%normni(:,:) = ninorm_in(1:dimx, 1:nions)
    in_regular%Ti(:,:) = Tix_in(1:dimx, 1:nions)

    in_regular%anis(:,:) = anis_in(1:dimx, 1:nions)
    in_regular%danisdr(:,:) = danisdr_in(1:dimx, 1:nions)
    in_regular%ion_type(:,:) = ion_type_in(1:dimx, 1:nions)

    in_regular%el_type = el_type_in
    in_regular%coll_flag = coll_flag_in
    in_regular%maxpts = maxpts_in
    in_regular%maxruns = maxruns_in
    in_regular%separateflux = separateflux_in
    in_regular%phys_meth = phys_meth
    in_regular%verbose = verbose_in
    in_regular%integration_routine = integration_routine_in
    in_regular%simple_mpi_only = simple_mpi_only_in
    in_regular%write_primi = 0
    in_regular%rot_flag = rot_flag_in

!    in_regular%R0 = R0_in
    in_regular%relacc1 = relacc1_in
    in_regular%relacc2 = relacc2_in
    in_regular%absacc2 = absacc2_in
    in_regular%absacc1 = absacc1_in
    in_regular%collmult = collmultin
    in_regular%ETGmult = ETGmultin
    in_regular%timeout = timeout_in
    rhomin = 0.
    rhomax = 1.
    in_regular%rhomin = rhomin !/rhoscale
    in_regular%rhomax = rhomax !/rhoscale

    OPEN(unit=unit_input, file="qualikiz/"//trim(fname1)//"/input.dat", status="replace", action="write")
    WRITE(unit_input, '(6(e14.6))') dq(jr), rmaj(jr), dv_r(jr), dr(jr), cs00
    WRITE(unit_input, '(6(e14.6))') Ane_in(1), Ate_in(1), Aupar_in(1), Autor_in(1), Machpar_in(1), Machtor_in(1)
    WRITE(unit_input, '(6(e14.6))') x_in(1), Bo_in(1), gammaE_in(1), Nex_in(1), qx_in(1), Ro_in(1)
    WRITE(unit_input, '(6(e14.6))') Rmin_in(1), smag_in(1), Tex_in(1), alphax_in(1), rho_in(1)/rhoscale, anise_in(1)
    WRITE(unit_input, '(6(e14.6))') Ai_in(1, 1:nions)
    WRITE(unit_input, '(6(e14.6))') Ani_in(1, 1:nions)
    WRITE(unit_input, '(6(e14.6))') Ati_in(1, 1:nions)
    WRITE(unit_input, '(6(e14.6))') Zi_in(1, 1:nions)
    WRITE(unit_input, '(6(e14.6))') ninorm_in(1, 1:nions)
    WRITE(unit_input, '(6(e14.6))') Tix_in(1, 1:nions)
    CLOSE(unit_input)

    print*, 'Calling qualikiz', nions, jr1 - 1 + jr, runcounter_in, dptot(jr)

    if (runcounter_in == 0) then
        call qualikiz(sizes, in_regular, &
            output_meth_0, &
            output_meth_0_sep_0_SI, output_meth_0_sep_0_GB, &
            primi_meth_0=primi_meth_0, &
            runcounterin=runcounter_in)
    else
        in_newt%oldsol = oldsol_in
        in_newt%oldfdsol = oldfdsol_in
        call qualikiz(sizes, in_regular, &
            output_meth_0, &
            output_meth_0_sep_0_SI, output_meth_0_sep_0_GB, &
            primi_meth_0=primi_meth_0, &
            in_newt=in_newt, runcounterin=runcounter_in)
    endif

    cftrans_out = output_meth_0%cftrans
    gam_GB_out  = output_meth_0_sep_0_GB%gam
    ome_GB_out  = output_meth_0_sep_0_GB%ome
    epf_GB_out  = output_meth_0_sep_0_GB%pfe
    eef_GB_out  = output_meth_0_sep_0_GB%efe
    ipf_GB_out  = output_meth_0_sep_0_GB%pfi
    ief_GB_out  = output_meth_0_sep_0_GB%efi
    sol_out   = primi_meth_0%sol
    fdsol_out = primi_meth_0%fdsol

    OPEN(unit=unit_runc, file="qualikiz/"//trim(fname1)//"/runcounter.dat", status="replace", action="write") !Replace old runcounter with new runcounter
    WRITE(unit_runc, *) runcounter_in + 1
    CLOSE(unit_runc)

    OPEN(unit=unit_rsol, file=trim(prim_dir)//'/rsol.dat', action="write", status="replace")
    WRITE(unit_rsol, fmtn) (((REAL(sol_out(i, j, k)), j=1, dimn), i=1, dimx), k=1, numsols)
    CLOSE(unit_rsol)

    OPEN(unit=unit_isol, file=trim(prim_dir)//'/isol.dat', action="write", status="replace")
    WRITE(unit_isol, fmtn) (((AIMAG(sol_out(i, j, k)), j=1, dimn), i=1, dimx), k=1, numsols)
    CLOSE(unit_isol)

    OPEN(unit=unit_rfd, file=trim(prim_dir)//'/rfdsol.dat', action="write", status="replace")
    WRITE(unit_rfd, fmtn) (((REAL(fdsol_out(i, j, k)), j=1, dimn), i=1, dimx), k=1, numsols)
    CLOSE(unit_rfd)

    OPEN(unit=unit_ifd, file=trim(prim_dir)//'/ifdsol.dat', action="write", status="replace")
    WRITE(unit_ifd, fmtn) (((AIMAG(fdsol_out(i, j, k)), j=1, dimn), i=1, dimx), k=1, numsols)
    CLOSE(unit_ifd)

    ql_fac = drhodr(jr)**2 * Rmin_in(1) * rhostar2 * cs00

    IF (ALLOCATED(oldsol_in)) THEN !Call with optional old solution input
        DEALLOCATE( oldsol_in  )
        DEALLOCATE( oldfdsol_in )
    endif

! Chii
 
    chii(jr)   = ql_fac * ief_gb_out(1, 1)/(1e-4 + Rmin_in(1)/R0_in * abs(Ati_in(1, 1)))
    chie(jr)   = ql_fac * eef_gb_out(1)/(1e-4 + Rmin_in(1)/R0_in * abs(Ate_in(1)))
    pfluxi(jr) = ql_fac * epf_gb_out(1)
    exchi(jr)  = cftrans_out(1, 2, 1)

enddo radial_loop

! Simulated TGLF computation:
prof_out(1, :) = chii
prof_out(2, :) = chie
prof_out(4, :) = pfluxi
prof_out(5, :) = exchi

end subroutine qlk_interf
