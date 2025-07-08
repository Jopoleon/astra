program qlk_chunk

use mpi
USE mod_qualikiz, only: qualikiz
USE kind, only: qlk_output_meth_0, qlk_output_meth_0_sep_0, &
    qlk_primi_meth_0, qlk_sizes, qlk_in_regular, qlk_in_newt
USE mod_make_io, only: zeroout_qlk_in_newt, zeroout_qlk_in_regular, &
    allocate_qlk_in_regular, allocate_qlk_in_newt
USE nanfilter

implicit none

type(qlk_sizes)      :: sizes
type(qlk_in_regular) :: in_regular
type(qlk_in_newt)    :: in_newt

type(qlk_output_meth_0)       :: output_meth_0
type(qlk_output_meth_0_sep_0) :: output_meth_0_sep_0_SI , output_meth_0_sep_0_GB
type(qlk_primi_meth_0)        :: primi_meth_0

integer, parameter :: ntheta=64, numecoefs=13, numicoefs=7, dimx=1, dimn=16, numsols=3, phys_meth=0, nradial=5, nspec_max=7

double precision, parameter :: &
    k0  = 1.6022E-12, &       ! erg/ev
    e0  = 4.8032E-10, &       ! elementary charge (statcoulombs)
    e00 = 1.6020e-19, &       ! elementary charge (C)
    c0  = 2.9979E+10, &       ! speed of light (cm/sec)
    mp  = 1.6726E-24, &       ! proton mass (g)
    mpp = 1.6726E-27, &       ! proton mass (kg)
    pi  = 3.141592653589793

!-----------------------------------------

integer :: simple_mpi_only_in, maxpts_in, maxruns_in, runcounter_in
integer :: nions, coll_flag_in, rot_flag_in, verbose_in, el_type_in, &
     integration_routine_in, separateflux_in
integer :: ierr, parent, rank, nproc, status(MPI_STATUS_SIZE)
integer :: chunk, nprocs, nrho_qlk, n_inputs, n_outputs, n_scalars, dims(8)
integer, dimension(dimx, nspec_max-1) :: ion_type_in

double precision :: a0_m, T0, m0, rmin_tg, drho_cs, drho_nt, nt_cs
double precision :: AMJ, RTOR, BTOR
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

INTEGER :: unit_runc=600, unit_rsol=610, unit_isol=620, unit_rfd=630, unit_ifd=640, myunit=700, i_mpic
integer :: i, j, k, jr, jion, jrho_beg

double precision, allocatable, dimension(:) :: drmin, drmaj, drho, dte, dne, dq, dptot, &
    dvper, drhodr, dstep, dr, dv_r
double precision :: Bunit, cs00, rhos00, omega0, rhostar2
! Shifted cicle geometry inputs
double precision :: gamma_e_qlk, mach_fac, ql_fac

double precision, allocatable, dimension(:) :: chie, chii, exchi, pfluxi, rho_qlk, &
    te_qlk, ne_qlk, vpar_qlk, vper_qlk, vexb_qlk, &
    ametr_qlk, rmaj_qlk, q_qlk
double precision, allocatable, dimension(:, :) :: dti, dni, ni_qlk, ti_qlk, zimp_qlk
double precision :: vpar_in, vpar_shear_in, cexb
double precision, allocatable :: inputs(:, :), output(:, :), scalars(:)

COMPLEX(kind=DBL), DIMENSION(:, :, :), ALLOCATABLE :: oldsol_in, oldfdsol_in
complex(kind=DBL), dimension(dimx, dimn, numsols) :: sol_out, fdsol_out

CHARACTER(len=20) :: fmtn
character(len=80) :: fname1, prim_dir

verbose_in = 1

call MPI_Init(ierr)
call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
call MPI_Comm_get_parent(parent, ierr)

if (parent == MPI_COMM_NULL) then
    if (rank == 0) then
        print *, "No parent communicator!"
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    endif
endif

! Receive dimensions from parent
call MPI_Recv(dims, 8, MPI_INTEGER, 0, 0, parent, status, ierr)
nrho_qlk  = dims(1)
n_scalars = dims(2)
n_inputs  = dims(3)
n_outputs = dims(4)
nions  = dims(5) - 1
jrho_beg = dims(7)

chunk = nrho_qlk / nprocs  ! Safe here: we now know nrho_qlk
allocate(scalars(n_scalars))
allocate(inputs(chunk, n_inputs))
allocate(output(chunk, n_outputs))
allocate( chie(chunk), chii(chunk), exchi(chunk), pfluxi(chunk), rho_qlk(chunk), &
    te_qlk(chunk), ne_qlk(chunk), vpar_qlk(chunk), vper_qlk(chunk), vexb_qlk(chunk), &
    ametr_qlk(chunk), rmaj_qlk(chunk), q_qlk(chunk) )
allocate( drmin(chunk), drmaj(chunk), drho(chunk), &
    dr(chunk), dne(chunk), dte(chunk), dq(chunk), dptot(chunk), &
    dvper(chunk), dv_r(chunk), drhodr(chunk) )
allocate( dti(4, chunk), dni(4, chunk), ni_qlk(4, chunk), ti_qlk(4, chunk), &
    zimp_qlk(3, chunk) )

! Receive TGLF input scalars and profiles from parent
call MPI_Recv(scalars     , n_scalars, MPI_DOUBLE_PRECISION, 0, 0, parent, status, ierr)
call MPI_Recv(inputs, chunk * n_inputs, MPI_DOUBLE_PRECISION, 0, 0, parent, status, ierr)
rho_qlk      = inputs(:,  1)
ametr_qlk    = inputs(:,  2)
rmaj_qlk     = inputs(:,  3)

q_qlk        = inputs(:,  6)

ne_qlk       = inputs(:,  9)
te_qlk       = inputs(:, 10)
vpar_qlk     = inputs(:, 12)
vper_qlk     = inputs(:, 13)
vexb_qlk     = inputs(:, 14)
ti_qlk(1, :) = inputs(:, 15)
ti_qlk(2, :) = inputs(:, 16)
ti_qlk(3, :) = inputs(:, 17)
ti_qlk(4, :) = inputs(:, 18)
ni_qlk(1, :) = inputs(:, 19)
ni_qlk(2, :) = inputs(:, 20)
ni_qlk(3, :) = inputs(:, 21)
ni_qlk(4, :) = inputs(:, 22)
zimp_qlk(1, :) = inputs(:, 23)
zimp_qlk(2, :) = inputs(:, 24)
zimp_qlk(3, :) = inputs(:, 25)
drmin       = inputs(:, 26)
drmaj       = inputs(:, 27)
drho        = inputs(:, 28)

dptot       = inputs(:, 31)
dte         = inputs(:, 32)
dne         = inputs(:, 33)
dq          = inputs(:, 34)

dvper       = inputs(:, 36)
dv_r        = inputs(:, 37)
dr          = inputs(:, 38)
drhodr      = inputs(:, 39)
dti(1, :)   = inputs(:, 40)
dti(2, :)   = inputs(:, 41)
dti(3, :)   = inputs(:, 42)
dti(4, :)   = inputs(:, 43)
dni(1, :)   = inputs(:, 44)
dni(2, :)   = inputs(:, 45)
dni(3, :)   = inputs(:, 46)
dni(4, :)   = inputs(:, 47)

RTOR = scalars(1)
BTOR = scalars(2)
a0_m = scalars(3)
AMJ  = scalars(5)
R0_in = scalars(14)
rhoscale = scalars(15)
m0 = AMJ*mp          ! Ref. mass = D ion mass [g]

! Electrons and main ions
Ai_in(1, 1: 4) = scalars( 5:  8)
Zi_in(1, 1: 4) = scalars(10: 13)

if (Zi_in(1, 4) >= 1.) then
    nions = 4
else
    nions = 3
endif
if (Zi_in(1, 3)  < 1.) nions = 2
if (Zi_in(1, 4) >= 1. .and. nions == 2) then 
    nions = 3
endif
if (Zi_in(1, 2) < 1.) nions = 1
if (Zi_in(1, 3) >= 1. .and. nions == 1) then 
    nions = 2
endif

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
   
    write(fname1, '(A, i0)') 'qlkzin_' , jrho_beg - 1 + jr
    qx_in(1)  = q_qlk(jr)
    rho_in(1) = rho_qlk(jr)
    x_in(1)   = ametr_qlk(jr)/Rmin_in(1)
    Ro_in(1)  = rmaj_qlk(jr)
    T0 = 1E3*te_qlk(jr)   ! temperature scale used by GYRO
    smag_in(1) = (x_in(1)/qx_in(1))*dq(jr)/dr(jr)        ! r/q dq/dr 

!thermal impurities

    Zi_in(1, 2) = max(1., zimp_qlk(1, jr))
    Zi_in(1, 3) = zimp_qlk(2, jr)
    Zi_in(1, 4) = zimp_qlk(3, jr)

    if (Zi_in(1, 4) >= 1. .and. nions == 2) then
        Zi_in(1, 3) = Zi_in(1, 4)
        Ai_in(1, 3) = Ai_in(1, 4)
        ni_qlk(3, :) = ni_qlk(4, :)
        ti_qlk(3, :) = ti_qlk(4, :)
    endif
    if (Zi_in(1, 3) >= 1. .and. nions == 1) then
        Zi_in(1, 2) = Zi_in(1, 3)
        Ai_in(1, 2) = Ai_in(1, 3)
        ni_qlk(2, :) = ni_qlk(3, :)
        ti_qlk(2, :) = ti_qlk(3, :)
    endif

    Ai_in(1, nions + 1:) = 0.
    Zi_in(1, nions + 1:) = 0.
    Ati_in(1, nions + 1:) = 0.
    Ani_in(1, nions + 1:) = 0.

! local field averages
    Nex_in(1) = ne_qlk(jr)
    Tex_in(1) = te_qlk(jr)
    Ate_in(1) = -dte(jr)*R0_in/(drmin(jr)*Tex_in(1))
    Ane_in(1) = -dne(jr)*R0_in/(drmin(jr)*Nex_in(1))
    anise_in(1) = 1.0
    danisedr_in(1) = 0.0
    do jion=1, nions  ! deut, impurities
        ninorm_in(1, jion) = ni_qlk(jion, jr)/Nex_in(1)
        Tix_in(1, jion) = ti_qlk(jion, jr)
        Ati_in(1, jion) = -dti(jion, jr)*R0_in/(drmin(jr)*ti_qlk(jion, jr))
        Ani_in(1, jion) = -dni(jion, jr)*R0_in/(drmin(jr)*ni_qlk(jion, jr))
        ion_type_in(1, jion) = 1
        anis_in    (1, jion) = 1.
        danisdr_in (1, jion) = 0.
    enddo

! Restore quasi-neutrality via main ions
    ninorm_in(1, 1) = (1 - SUM(ninorm_in(1, 2:nions)*Zi_in(1, 2:nions)))/Zi_in(1, 1)
    Ani_in(1, 1) = (Ane_in(1) - SUM(ninorm_in(1, 2: nions)*Ani_in(1, 2: nions)* &
               Zi_in(1, 2: nions)))/(ninorm_in(1, 1)*Zi_in(1, 1))

! derived units for the plasma

    Bunit = 1E4*BTOR*drhodr(jr)*rho_in(1)/ametr_qlk(jr)  ! Miller geometry magnetic field unit
    cs00 = SQRT(e00*T0/(AMJ*mpp))  ! thermal velocity unit m/sec
    omega0 = e0*Bunit/(m0*c0)      ! gyrofrequency unit 1/sec
    rhos00 = cs00/omega0           ! gyroradius unit m

    vpar_shear_in = -rmaj_qlk(jr)*dv_r(jr)/(dr(jr)*cs00)  !From m/s to cm/s for vpar
    vpar_in = vpar_qlk(jr)/cs00

! local magnetic geometry

    rhostar2 = (rhos00/Rmin_in(1))**2
    drho_cs = drhodr(jr)**2 * Rmin_in(1) * rhostar2 * cs00

! 1.e16 comes from cgs to SI for ptot
    alphax_in(1) = -(0.0040267/BTOR**2) * qx_in(1)**2 * R0_in * dptot(jr)/(dr(jr)*Rmin_in(1))
    cexb = ametr_qlk(jr)/qx_in(1)            ! r/(q)
    gamma_e_qlk = cexb*dvper(jr)/(dr(jr)*cs00) ! Waltz-Miller definition
    mach_fac = sqrt(Tex_in(1)/AMJ)

    Machtor_in(1) = vpar_in*mach_fac
    Machpar_in(1) = vpar_in*mach_fac
    Autor_in(1)  = vpar_shear_in*R0_in/Rmin_in(1)*mach_fac
    Aupar_in(1)  = vpar_shear_in*R0_in/Rmin_in(1)*mach_fac
    gammaE_in(1) = gamma_e_qlk  *R0_in/Rmin_in(1)*mach_fac

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

    print*, 'Calling qualikiz', nions, jrho_beg - 1 + jr, runcounter_in, rank, dptot(jr)

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
output(:, 1) = chii
output(:, 2) = chie
output(:, 4) = pfluxi
output(:, 5) = exchi

call MPI_Send(output, chunk * n_outputs, MPI_DOUBLE_PRECISION, 0, 1, parent, ierr)
call MPI_Finalize(ierr)

end program qlk_chunk
