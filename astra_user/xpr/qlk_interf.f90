!----------------------------------------------------------------------|
program main

implicit none

integer :: iargc, mampid, mamkey, eignr
character*132 :: STRING, eigpath, mampath

if (iargc() /= 4) then
    write(6, '(A)') "Error"
    call a_stop
endif
call getarg(0, STRING)
eigpath = STRING(1:len_trim(STRING))//char(0)
call getarg(1, STRING)
mampath = STRING(1:len_trim(STRING))//char(0)
call getarg(2, STRING)
read(STRING, *) mampid
call getarg(3, STRING)
read(STRING, *) mamkey
call getarg(4, STRING)
read(STRING, *) eignr
call sbp2shm(eigpath, mampath, mampid, mamkey, eignr)

end program main


!----------------------------------------------------------------------|
subroutine qlk_interf(jr1_in, jr2_in, nrho, NA1N, NA1E, NA1I, &
    BTOR, RTOR, AMJ, ZMJ, AIM1, AIM2, AIM3, &
    NE, TE, NI, NDEUT, NTRIT, NIZ1, NIZ2, TI, ZEF, ZIM1, AMAIN, &
    MU, RHO, AMETR, SHIF, ELON, TRIA, ER, NIBM, G11, &
    VPOL, VRS, VTOR, SHEAR, PBLON, PBPER, PFAST, NIZ3, ZIM2, ZIM3, &
! output
    CHI, CHE, DIF, VIN, DPH, DPL, DPR, XTB, EGM, GAM, GM1, GM2, OM1, OM2, FR1)
 
USE mod_qualikiz, only: qualikiz
USE kind, only: qlk_output_meth_0, qlk_output_meth_0_sep_0, qlk_output_meth_0_sep_1, &
         qlk_output_meth_1, qlk_output_meth_1_sep_0, qlk_output_meth_1_sep_1, &
         qlk_output_meth_2, qlk_output_meth_2_sep_0, qlk_output_meth_2_sep_1, &
         qlk_primi_meth_0, qlk_primi_meth_1, qlk_primi_meth_2, qlk_sizes, &
         qlk_in_regular, qlk_in_newt
USE mod_make_io, only: zeroout_qlk_in_newt, zeroout_qlk_in_regular, &
         allocate_qlk_in_regular, allocate_qlk_in_newt
USE mpi
USE nanfilter

implicit none

integer, intent(in) :: nrho, jr1_in, jr2_in, NA1N, NA1E, NA1I

double precision, intent(in) :: BTOR, RTOR, &
    AMJ, AIM1, AIM2, AIM3, ZMJ
double precision, intent(in), dimension(*) :: NE, TE, NI, TI, &
    ZEF, ZIM1, ZIM2, ZIM3, PBLON, PBPER, PFAST, NIZ3, AMAIN, &
    ER, MU, RHO, AMETR, SHIF, ELON, NDEUT, NIZ1, NTRIT, &
    NIZ2, TRIA, NIBM, G11, VPOL, VRS, VTOR, SHEAR

double precision, intent(out), dimension(*) :: CHI, CHE, DIF, VIN, &
    DPH, DPL, DPR, XTB, EGM, GAM, GM1, GM2, OM1, OM2, FR1

!--------------------------------
type(qlk_sizes) :: sizes
type(qlk_in_regular) :: in_regular
type(qlk_in_newt) :: in_newt

type(qlk_output_meth_0) :: output_meth_0
type(qlk_output_meth_0_sep_0) :: output_meth_0_sep_0_SI , output_meth_0_sep_0_GB
type(qlk_output_meth_0_sep_1) :: output_meth_0_sep_1_SI, output_meth_0_sep_1_GB
type(qlk_output_meth_1) :: output_meth_1
type(qlk_output_meth_1_sep_0) :: output_meth_1_sep_0_SI, output_meth_1_sep_0_GB
type(qlk_output_meth_1_sep_1) :: output_meth_1_sep_1_SI, output_meth_1_sep_1_GB
type(qlk_output_meth_2)       :: output_meth_2
type(qlk_output_meth_2_sep_0) :: output_meth_2_sep_0_SI, output_meth_2_sep_0_GB
type(qlk_output_meth_2_sep_1) :: output_meth_2_sep_1_SI, output_meth_2_sep_1_GB
type(qlk_primi_meth_0)        :: primi_meth_0
type(qlk_primi_meth_1)        :: primi_meth_1
type(qlk_primi_meth_2)        :: primi_meth_2


integer, parameter :: ntheta=64, numecoefs=13, numicoefs=7, dimx=1, dimn=16, numsols=3, phys_meth=0, jpd=700, nradial=5, nspec_max=7

real, parameter :: &
   k0   = 1.6022E-12, &       ! erg/ev
   e0   = 4.8032E-10, &       ! elementary charge (statcoulombs)
   e00  = 1.6020e-19, &       ! elementary charge (C)
   c0   = 2.9979E+10, &       ! speed of light (cm/sec)
   mp   = 1.6726E-24, &       ! proton mass (g)
   mpp  = 1.6726E-27, &       ! proton mass (kg)
   pi   = 3.141592653589793 

!-----------------------------------------

integer :: simple_mpi_only_in, maxpts_in, maxruns_in, runcounter_in, nspecies
integer :: nions, coll_flag_in, rot_flag_in, verbose_in, el_type_in, &
    integration_routine_in, separateflux_in
integer, dimension(dimx, nspec_max-1) :: ion_type_in

real(kind=DBL) :: relacc1_in, relacc2_in, absacc1_in, absacc2_in, R0_in, &
   ETGmultin, collmultin, timeout_in, rhomin, rhomax, rhoscale
real(kind=DBL), dimension(dimn) :: kthetarhos_in
real(kind=DBL), dimension(dimx) :: x_in, rho_in, Ro_in, Rmin_in, Bo_in, &
    qx_in, smag_in, alphax_in, Tex_in, Nex_in, Ate_in, Ane_in, anise_in, &
    danisedr_in, Machtor_in, Autor_in, Machpar_in, Aupar_in, gammaE_in, &
    epf_GB_out, eef_GB_out
real(kind=DBL), dimension(dimx, nspec_max-1) :: Tix_in, ninorm_in, &
    Ati_in, Ani_in, anis_in, danisdr_in, Ai_in, Zi_in, ipf_GB_out, ief_GB_out
real(kind=DBL), dimension(dimx, dimn, numsols) :: gam_GB_out, ome_GB_out
real(kind=DBL), dimension(dimx, nspec_max - 1, numicoefs) :: cftrans_out

! Old solution for non-reset runs
REAL(KIND=DBL), DIMENSION(:, :, :), ALLOCATABLE :: oldrsol, oldisol, oldrfdsol, oldifdsol

!-----------------------------------------
LOGICAL :: exist1, exist2, exist3, exist4, exist5 !used for checking for existence of files

!MPI variables:
INTEGER :: mpi_ierr, nproc, myrank
INTEGER :: myunit=700, i_mpic
integer :: jna, jinterval, jr_min, jr_max, jrho, j0, j01, j02
integer :: i, j, k, jradial, j3r, jjgrid(nradial), j1, j2, jion

real(kind=DBL) :: bmod, bpolz
real(kind=DBL) :: drmin, drmaj, drho, dte, dne, dq, dptot, &
        delong, dtrian, dvper, drhodr, dstep, dr, dv_r
real(kind=DBL) :: Bunit, cs00, rhos00, omega0, rhostar2
real(kind=DBL) :: T0, m0, drho_cs
! Shifted cicle geometry inputs
real(kind=DBL) :: gamma_e_tg, vpar_tg, mach_fac, ql_fac

real(kind=DBL), dimension(nrho) :: vexb2, vpar_m, vper_m, &
    gradrhosq_exp, rmaj_exp, q_exp, &
    chie_m, chii_m, pfluxi_m, exchi_m, ptot

real(kind=DBL), dimension(nradial) :: chie, chii, exchi, pfluxi, rho_tg
real(kind=DBL), dimension(nspec_max-1) :: dti, dni
real(kind=DBL), dimension(nspec_max-1, nrho) :: ni_m, ti_m
real(kind=DBL) :: vpar_in, vpar_shear_in, cexb

COMPLEX(kind=DBL), DIMENSION(:, :, :), ALLOCATABLE :: oldsol_in, oldfdsol_in
complex(kind=DBL), dimension(dimx, dimn, numsols) :: sol_out, fdsol_out

CHARACTER(len=20) :: fmtn
character(len=80) :: fname1, prim_dir
logical :: oldsol_loaded = .false.

!---------------------------------

save i_mpic
data i_mpic /0/

verbose_in = 0

if (i_mpic == 0) then
    if (verbose_in > 0) then
        write(6, *) 'GIt qlk_interface first call mpi_init'
    endif
    CALL mpi_init(mpi_ierr)
    if (mpi_ierr /= 0) then
        write(6, *) 'GIt qlk_interface mpi_init error=', mpi_ierr
    endif
    CALL mpi_comm_size(mpi_comm_world, nproc,  mpi_ierr)
    if (verbose_in > 0) then
        write(6, *) 'GIt qlk_interface debug 2'
    endif
    CALL mpi_comm_rank(mpi_comm_world, myrank, mpi_ierr)
    i_mpic = 1
endif

nions = nspec_max - 1

jna = max(NA1E, NA1I, NA1N)
if (jna == 0) jna = nrho
if (jna > jpd) then
    return
endif

! Electrons and main ions
Zi_in(1, 1) = ZMJ

Ai_in(1, 1) = AMJ
Ai_in(1, 2) = AIM1
Ai_in(1, 3) = AIM2
Ai_in(1, 4) = AIM3

do jrho=1, nrho
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

m0 = AMJ*mp                ! Ref. mass = D ion mass [g]

pfluxi_m  = 0.0
chie_m  = 0.0
chii_m  = 0.0
exchi_m = 0.0

jr_min = max(1, jr1_in)
jr_max = min(nrho, jr2_in)

jinterval = 1 + (jr_max - jr_min)
j3r = max(1, INT(jinterval/nradial))
do jradial = 1, nradial-1
    jjgrid(jradial) = jr_min + (jradial - 1)*j3r
enddo
jjgrid(nradial) = min(nrho, jr_max)


! These will be reset locally in the radial loop
Zi_in(1, 2) = MAXVAL(ZIM1(1:nrho))
Zi_in(1, 3) = MAXVAL(ZIM2(1:nrho))
Zi_in(1, 4) = MAXVAL(ZIM3(1:nrho))

if (Zi_in(1, 4) >= 1.) nions = 4
if (Zi_in(1, 4) < 1.) nions = 3
if (Zi_in(1, 3) < 1.) nions = 2
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
maxpts_in   = 50000000 ! 500000 default
maxruns_in  = 10

R0_in = RTOR + SHIF(nrho)
Bo_in(1) = BTOR
Rmin_in(1) = AMETR(nrho)
rhoscale = rho(nrho)

WRITE(fmtn, '(A, I0, A)') '(', dimn, 'G15.7)'

radial_loop: do jradial=1, nradial
   
    j0 = jjgrid(jradial)
    write(fname1, '(A, i0)') 'qlkzin_' , j0
    qx_in(1) = q_exp(j0)
    rho_in(1) = RHO(j0)
    x_in(1) = AMETR(j0)/Rmin_in(1)
    rho_tg(jradial) = rho_in(1)

    Ro_in(1) = RTOR + SHIF(j0)
    T0  = 1E3 *te(j0)   ! temperature scale used by GYRO

!thermal impurities

    Zi_in(1, 2) = max(1., ZIM1(j0))
    Zi_in(1, 3) = ZIM2(j0)
    Zi_in(1, 4) = ZIM3(j0)

    if (Zi_in(1, 4) >= 1. .and. nions == 2) then 
        Zi_in(1, 3) = Zi_in(1, 4)
        Ai_in(1, 3) = Ai_in(1, 4)
        ni_m(3, :) = ni_m(4, :)
        ti_m(3, :) = ti_m(4, :)
    endif
    if (Zi_in(1, 3) >= 1. .and. nions == 1) then 
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

    j01 = j0+1
    j02 = j0-1
    if (j0 == 1) then
        j02 = j0
    else if (j0 == nrho) then
        j01 = j0
    endif
    dstep = 1./float(j01 - j02)

    drmin  = dstep*(AMETR(j01) - AMETR(j02))
    drmaj  = dstep*(rmaj_exp(j01) - rmaj_exp(j02))
    drho   = dstep*(rho(j01) - rho(j02))
    delong = dstep*(ELON(j01) - ELON(j02))
    dtrian = dstep*(TRIA(j01) - TRIA(j02))
    dptot  = dstep*(ptot(j01) - ptot(j02))
    dte    = dstep*(TE(j01) - TE(j02))
    dne    = dstep*(NE(j01) - NE(j02))
    dq     = dstep*(q_exp(j01) - q_exp(j02))
    dvper  = dstep*(vper_m(j01) - vper_m(j02))
    do jion=1, nions
        dti(jion) = dstep*(ti_m(jion, j01) - ti_m(jion, j02))
        dni(jion) = dstep*(ni_m(jion, j01) - ni_m(jion, j02))
    enddo
    dv_r = dstep* &
        (vpar_m(j01)/(rmaj_exp(j01) + AMETR(j01)) - &
         vpar_m(j02)/(rmaj_exp(j02) + AMETR(j02)))
    dr = drmin/Rmin_in(1)    ! gradients w.r.t. minor radius even for s-alpha geometry
    drhodr = drho/drmin
    smag_in(1) = (x_in(1)/qx_in(1))*dq/dr        ! r/q dq/dr 

! local field averages
    Nex_in(1) = NE(j0)
    Tex_in(1) = TE(j0)
    Ate_in(1) = -dte*R0_in/(drmin*Tex_in(1))
    Ane_in(1) = -dne*R0_in/(drmin*Nex_in(1))
    anise_in(1) = 1.0
    danisedr_in(1) = 0.0
    do jion=1, nions  ! deut, impurities
        ninorm_in(1, jion) = ni_m(jion, j0)/Nex_in(1)
        Tix_in(1, jion) = ti_m(jion, j0)
        Ati_in(1, jion) = -dti(jion)*R0_in/(drmin*ti_m(jion, j0))
        Ani_in(1, jion) = -dni(jion)*R0_in/(drmin*ni_m(jion, j0))
        ion_type_in(1, jion) = 1
        anis_in    (1, jion) = 1.
        danisdr_in (1, jion) = 0.
    enddo

! Restore quasi-neutrality via main ions
    ninorm_in(1, 1) = (1 - SUM(ninorm_in(1, 2:nions)*Zi_in(1, 2:nions)))/Zi_in(1, 1)
    Ani_in(1, 1) = (Ane_in(1) - SUM(ninorm_in(1, 2: nions)*Ani_in(1, 2: nions)* &
               Zi_in(1, 2: nions)))/(ninorm_in(1, 1)*Zi_in(1, 1))

! derived units for the plasma

    Bunit = 1E4*BTOR*drhodr*rho_in(1)/AMETR(j0)  ! Miller geometry magnetic field unit
    cs00 = SQRT(e00*T0/(AMJ*mpp))  ! thermal velocity unit m/sec
    omega0 = e0*Bunit/(m0*c0)      ! gyrofrequency unit 1/sec
    rhos00 = cs00/omega0           ! gyroradius unit m

    vpar_shear_in = -rmaj_exp(j0)*dv_r/(dr*cs00)  !From m/s to cm/s for vpar
    vpar_in = vpar_m(j0)/cs00

! local magnetic geometry

    rhostar2 = (rhos00/Rmin_in(1))**2
    drho_cs = drhodr**2 * Rmin_in(1) * rhostar2 * cs00

! 1.e16 comes from cgs to SI for ptot
    alphax_in(1) = -(0.0040267/BTOR**2) * qx_in(1)**2 * R0_in * dptot/(dr*Rmin_in(1))
    cexb = AMETR(j0)/qx_in(1)            ! r/(q)
    gamma_e_tg = cexb*dvper/(dr*cs00) ! Waltz-Miller definition
    mach_fac = sqrt(Tex_in(1)/AMJ)

    Machtor_in(1) = vpar_in*mach_fac
    Machpar_in(1) = vpar_in*mach_fac
    Autor_in(1)  = vpar_shear_in*R0_in/Rmin_in(1)*mach_fac
    Aupar_in(1)  = vpar_shear_in*R0_in/Rmin_in(1)*mach_fac
    gammaE_in(1) = gamma_e_tg   *R0_in/Rmin_in(1)*mach_fac

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
        write(6, *) 'aNe', Ane_in
        write(6, *) 'aNi', Ani_in(1, 1:nions)
        write(6, *) 'ate', Ate_in
        write(6, *) 'ati', Ati_in(1, 1:nions)
        write(6, *) 'Ai', Ai_in(1, 1:nions)
        write(6, *) 'Ni/Ne', ninorm_in(1, 1:nions)
        write(6, *) 'check: ', Ati_in(1, 1), gammaE_in, Autor_in, Machtor_in, alphax_in(1), x_in(1)
    endif
 
    write(1331, *) 'inputs,total'
    write(1331, *) dimx, rho_in/rhoscale, dimn, nions, numsols, phys_meth, coll_flag_in, &
         rot_flag_in, verbose_in, 0,  kthetarhos_in, & !general param
         'x', x_in, Ro_in, Rmin_in, R0_in, Bo_in, qx_in, smag_in, alphax_in, & !geometry
         el_type_in, Tex_in, Nex_in, Ate_in, Ane_in, anise_in, danisedr_in, & !electrons
         ion_type_in(1, 1:nions), Ai_in(1, 1:nions), Zi_in(1, 1:nions), Tix_in(1, 1:nions), ninorm_in(1, 1:nions), Ati_in(1, 1:nions), Ani_in(1, 1:nions), anis_in(1, 1:nions), danisdr_in(1, 1:nions), & !ions
         Machtor_in, Autor_in, Machpar_in, Aupar_in, gammaE_in

    INQUIRE(file="../qualikiz/"//trim(fname1)//"/runcounter.dat", EXIST=exist1)
    INQUIRE(file=trim(prim_dir)//'/rsol.dat', EXIST=exist2)
    INQUIRE(file=trim(prim_dir)//'/isol.dat', EXIST=exist3)
    INQUIRE(file=trim(prim_dir)//'/rfdsol.dat', EXIST=exist4)
    INQUIRE(file=trim(prim_dir)//'/ifdsol.dat', EXIST=exist5)

    IF ( exist1 .AND. exist2 .AND. exist3 .AND. exist4 .AND. exist5 ) THEN
        OPEN(unit=700, file="../qualikiz/"//trim(fname1)//"/runcounter.dat", status="old", action="read")
        READ(700,*) runcounter_in
        CLOSE(700)
    ELSE
        prim_dir = '../qualikiz/'//trim(fname1)//'/output/primitive'
        call system('mkdir -p ' // trim(prim_dir))
        runcounter_in = 0 !First run
    ENDIF

    IF (runcounter_in >= maxruns_in) THEN !Reset if we're at our maximum number of runs
        runcounter_in = 0
    ENDIF
    runcounter_in = 0 ! GIT force calculation from scratch

    IF (runcounter_in == 0) THEN !load old rsol and isol if we're not doing a reset run
        write(6, *) 'Qualikiz from scratch'
    ELSE
        ALLOCATE( oldrsol  (dimx, dimn, numsols) )
        ALLOCATE( oldisol  (dimx, dimn, numsols) )
        ALLOCATE( oldrfdsol(dimx, dimn, numsols) )
        ALLOCATE( oldifdsol(dimx, dimn, numsols) )
        if (.not. ALLOCATED(oldsol_in)) then
            ALLOCATE( oldsol_in   (dimx, dimn, numsols) )
            ALLOCATE( oldfdsol_in (dimx, dimn, numsols) )
        endif
        OPEN(unit=myunit, file=trim(prim_dir)//'/rsol.dat', action="read", status="old")
        READ(myunit,fmtn) (((oldrsol(i, j, k), j=1, dimn), i=1, dimx), k=1, numsols)
        CLOSE(myunit)

        OPEN(unit=myunit, file=trim(prim_dir)//'/isol.dat', action="read", status="old")
        READ(myunit,fmtn) (((oldisol(i, j, k), j=1, dimn), i=1, dimx), k=1, numsols)
        CLOSE(myunit)

        OPEN(unit=myunit, file=trim(prim_dir)//'/rfdsol.dat', action="read", status="old")
        READ(myunit,fmtn) (((oldrfdsol(i, j, k), j=1, dimn), i=1, dimx), k=1, numsols)
        CLOSE(myunit)

        OPEN(unit=myunit, file=trim(prim_dir)//'/ifdsol.dat', action="read", status="old")
        READ(myunit,fmtn) (((oldifdsol(i, j, k), j=1, dimn), i=1, dimx), k=1, numsols)
        CLOSE(myunit)

        oldsol_in   = CMPLX(oldrsol  , oldisol)
        oldfdsol_in = CMPLX(oldrfdsol, oldifdsol)

        DEALLOCATE( oldrsol )
        DEALLOCATE( oldisol )
        DEALLOCATE( oldrfdsol )
        DEALLOCATE( oldifdsol )
        write(6, *) 'Call qualikiz from old solution'
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

    in_regular%R0 = R0_in
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

    write(6, *) 'Calling qualikiz', nions, j0, runcounter_in, myrank

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

    OPEN(unit=700, file="../qualikiz/"//trim(fname1)//"/runcounter.dat", status="replace", action="write") !Replace old runcounter with new runcounter
    WRITE(700,*) runcounter_in + 1
    CLOSE(700)

    OPEN(unit=myunit, file=trim(prim_dir)//'/rsol.dat', action="write", status="replace")
    WRITE(myunit, fmtn) (((REAL(sol_out(i, j, k)), j=1, dimn), i=1, dimx), k=1, numsols)
    CLOSE(myunit)

    OPEN(unit=myunit, file=trim(prim_dir)//'/isol.dat', action="write", status="replace")
    WRITE(myunit, fmtn) (((AIMAG(sol_out(i, j, k)), j=1, dimn), i=1, dimx), k=1, numsols)
    CLOSE(myunit)

    OPEN(unit=myunit, file=trim(prim_dir)//'/rfdsol.dat', action="write", status="replace")
    WRITE(myunit, fmtn) (((REAL(fdsol_out(i, j, k)), j=1, dimn), i=1, dimx), k=1, numsols)
    CLOSE(myunit)

    OPEN(unit=myunit, file=trim(prim_dir)//'/ifdsol.dat', action="write", status="replace")
    WRITE(myunit, fmtn) (((AIMAG(fdsol_out(i, j, k)), j=1, dimn), i=1, dimx), k=1, numsols)
    CLOSE(myunit)

    ql_fac = drhodr**2 * Rmin_in(1) * rhostar2 * cs00

    IF (ALLOCATED(oldsol_in)) THEN !Call with optional old solution input
        DEALLOCATE( oldsol_in  )
        DEALLOCATE( oldfdsol_in )
    endif

! Chii
 
    chii(jradial)   = ql_fac * ief_gb_out(1, 1)/(1e-4 + Rmin_in(1)/R0_in * abs(Ati_in(1, 1)))
    chie(jradial)   = ql_fac * eef_gb_out(1)/(1e-4 + Rmin_in(1)/R0_in * abs(Ate_in(1)))
    pfluxi(jradial) = ql_fac * epf_gb_out(1)
    exchi(jradial)  = cftrans_out(1, 2, 1)

enddo radial_loop

j1 = jjgrid(1)
j2 = jjgrid(nradial)

call qinterp(rho_tg, chii  , nradial, RHO(j1:j2), chii_m(j1:j2)  , j2-j1+1)
call qinterp(rho_tg, chie  , nradial, RHO(j1:j2), chie_m(j1:j2)  , j2-j1+1)
call qinterp(rho_tg, pfluxi, nradial, RHO(j1:j2), pfluxi_m(j1:j2), j2-j1+1)
call qinterp(rho_tg, exchi , nradial, RHO(j1:j2), exchi_m(j1:j2) , j2-j1+1)

chii_m  (1:2) = chii_m(3)
chie_m  (1:2) = chie_m(3)
pfluxi_m(1:2) = pfluxi_m(3)
exchi_m (1:2) = exchi_m(3)

DIF(1:nrho) = 0.d0                       ! D, electron diffusivity, m^2/s
DPH(1:nrho) = 0.d0                       ! D, impurity diffusivity, m^2/s
DPL(1:nrho) = 0.d0                       ! impurity convection
DPR(1:nrho) = 0.d0                       ! tor. stress
EGM(1:nrho) = 0.d0
GAM(1:nrho) = 0.d0
GM1(1:nrho) = 0.d0
GM2(1:nrho) = 0.d0
OM1(1:nrho) = 0.d0
OM2(1:nrho) = 0.d0
FR1(1:nrho) = 0.d0
do j=jr_min, jr_max
    CHI(j) = chii_m(j)/gradrhosq_exp(j) ! \chi_i, m^2/s : starts from work(21,:) 
    CHE(j) = chie_m(j)/gradrhosq_exp(j) ! \chi_e, m^2/s
    VIN(j) = min(20., pfluxi_m(j)/AMETR(nrho)/gradrhosq_exp(j)) ! D flux
    VIN(j) = max(-20., VIN(j)) ! D flux
!    VIN(j) = pfluxi_m(j)/AMETR(nrho)/gradrhosq_exp(j)
    XTB(j) = exchi_m(j)  ! turbulent e-i equipartition in MW/m^3
enddo   ! End of main loop

if (jr_max < nrho-1) return

do j=jna-2, jna-1
    CHI(j) = CHI(jna-3)
    CHE(j) = CHE(jna-3)
    VIN(j) = VIN(jna-3)
    XTB(j) = XTB(jna-3)
enddo

if (jna < nrho) then
    do j=jna-1, nrho
        CHI(j) = 0.d0
        CHE(j) = 0.d0
        VIN(j) = 0.d0
        XTB(j) = 0.d0
    enddo
endif  

return
END subroutine qlk_interf
