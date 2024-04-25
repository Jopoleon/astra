!----------------------------------------------------------------------|
program main

implicit none

integer :: iargc, mampid, mamkey, eignr
character*132 :: STRING, eigpath, mampath

if (iargc() .ne. 4) then
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
subroutine neo_interf(jr1_in, jr2_in, nrho, NA1N, NA1E, NA1I, &
    BTOR, RTOR, AMJ, ZMJ, AIM1, AIM2, AIM3, &
    NE, TE, NI, NDEUT, NTRIT, NIZ1, NIZ2, TI, ZEF, ZIM1, AMAIN, &
    MU, RHO, AMETR, SHIF, ELON, TRIA, ER, NIBM, G11, &
    VPOL, VRS, VTOR, SHEAR, PBLON, PBPER, PFAST, NIZ3, ZIM2, ZIM3, &
    ZIMPT, NIMPT, AIMPT, &
! output
    CHI, CHE, DIF, VIN, DPH, DPL, DPR, XTB, EGM, GAM, GM1, GM2, OM1, OM2, FR1)

!----------------------------------------------------------------------|
! WORK(1:NA1,1:13) array is used for output
!                              (when i_delay=0 and egamma_d is not used)
!----------------------------------------------------------------------|

use neo_interface

implicit none

logical, parameter :: verbose=.False.
integer, parameter :: jpd=700, nradial=5, nsm=11

real, parameter :: &
   k0   = 1.6022E-12, &       ! erg/ev
   e00  = 1.6020e-19, &       ! elementary charge (C)
   c0   = 2.9979E+10, &       ! speed of light (cm/sec)
   mpp  = 1.6726E-27, &       ! proton mass (kg)
   pi   = 3.141592653589793 

integer, intent(in) :: nrho, jr1_in, jr2_in, NA1N, NA1E, NA1I

double precision, intent(in) :: BTOR, RTOR, &
    AMJ, AIM1, AIM2, AIM3, ZMJ
double precision, intent(in), dimension(*) :: NE, TE, NI, TI, &
    ZEF, ZIM1, ZIM2, ZIM3, PBLON, PBPER, PFAST, NIZ3, AMAIN, &
    ER, MU, RHO, AMETR, SHIF, ELON, NDEUT, NIZ1, NTRIT, &
    NIZ2, TRIA, NIBM, G11, VPOL, VRS, VTOR, SHEAR, &
    ZIMPT, NIMPT, AIMPT

double precision, intent(out), dimension(*) :: CHI, CHE, DIF, VIN, &
    DPH, DPL, DPR, XTB, EGM, GAM, GM1, GM2, OM1, OM2, FR1

!----------------------------------------------------------------------
integer :: jr_min, jr_max, jrho, j0, j01, j02, n_radial
integer :: j, jradial, jjgrid(nradial), jspec
integer :: i_ion, n_ions
real :: bmod, bpolz, alpha_zf_in, ion_eflux
real :: drmin, drmaj, drho, dte, dne, dq, &
        delong, dtrian, dvpar, dvper, drhodr, dstep, dr, xstep
real :: Bunit, cs0, rhos0, omega0, rhostar2, lnlamda, taue, cexb, xnuei
real :: T0, anorm, mnorm, tnorm, nnorm, vnorm, &
   pflux_e_neo, eflux_e_neo, jboots, tgyro_neo_gv_flag, &
   Gamma_neo_GB, Q_neo_GB, Pi_neo_GB, Jpar_GB

real, dimension(nrho) :: rho_m, vexb2, vpar_m, vper_m, &
    gradrhosq_exp, epar0_in, rmaj_exp, q_exp, &
    chie_m, chii_m, vippd_m, vittd_m, vippi1_m, vitti1_m, j_boot, elec_pflux_m

real, dimension(nradial) :: chie, chii, elec_pflux, rho_tg, &
   vippd, vittd, vippi1, vitti1, jbs
real, dimension(nsm-1) :: dti, dni, pflux_i_neo, eflux_i_neo, vpflux_neo, vtflux_neo
real, dimension(nsm-1, nrho) :: ni_m, ti_m
real, dimension(nsm, 2) :: energy_flux, particle_flux
 
character(len=80) :: path_in

!-----------------
! Radial subdomain
!-----------------

jr_min = max(1, jr1_in)
jr_max = min(nrho, jr2_in)
if (jr_min > jr_max) then
    write(*, *) 'Error! jr_min > jr_max', jr_max
    return
endif

xstep = float(jr_max - jr_min)/(nradial - 1.)
if (xstep <= 1.) then
    do jradial=1, nradial
        jjgrid(jradial) = jr_min + jradial - 1
        if (jjgrid(jradial) == jr_max) EXIT
    enddo
    n_radial = jradial
else
    n_radial = nradial
    do jradial = 1, n_radial-1
        jjgrid(jradial) = jr_min + nint((jradial-1)*xstep)
    enddo
    jjgrid(n_radial) = jr_max
endif

!-----------------
! Flags and inputs

tgyro_neo_gv_flag = 0.

!Main ions
neo_z_in(1) = -1.
neo_z_in(2) = ZMJ

neo_mass_in(1) = 5.4447e-4/AMJ
neo_mass_in(2) = AMJ/AMJ  ! AMJ is reference mass
neo_mass_in(3) = AIM1/AMJ
neo_mass_in(4) = AIM2/AMJ
neo_mass_in(5) = AIM3/AMJ

do jrho=1, nrho
    rho_m(jrho) = RHO(jrho)
    ti_m(1:4, jrho) = TI(jrho)
    if (NDEUT(jrho) >= 0.01*NE(jrho)) then
        ni_m(1, jrho) = NDEUT(jrho)
    else ! likely: NDEUT not defined in equ file, hence zero
        ni_m(1, jrho) = NI(jrho)
    endif
    ni_m(2, jrho) = max(1.e-9, NIZ1(jrho))
    ni_m(3, jrho) = max(1.e-9, NIZ2(jrho))
    ni_m(4, jrho) = max(1.e-9, NIZ3(jrho))
    rmaj_exp(jrho) = RTOR + SHIF(jrho)
    q_exp(jrho)    = 1./MU(jrho)
    bpolz = BTOR*AMETR(jrho)*MU(jrho)/RTOR
    bmod = sqrt(BTOR**2 + bpolz**2)
    gradrhosq_exp(jrho) = G11(jrho)/VRS(jrho)
    epar0_in(jrho) = 0 !NIZ3(jrho)*AMETR(nrho)/(TE(jrho)*1.e3) !NIZ3 is supposed to be Epar*e*a/T

    vper_m(jrho) = ER(jrho)/(RTOR*bpolz) ! vexb in m/s --> Omega_E    
    vpar_m(jrho) = ER(jrho)/(RTOR*bpolz)*(RTOR+SHIF(jrho)+AMETR(jrho))  !--> R*Omega_E , no neoclassical terms
    vexb2(jrho)  = -ER(jrho)/bmod ! vexb in m/s (vperp = vexb since the diamagnetic velocity is the curvature drift actually

enddo

elec_pflux_m  = 0.0
chie_m   = 0.0
chii_m   = 0.0
vippd_m  = 0.0
vittd_m  = 0.0
vippi1_m = 0.0
vitti1_m = 0.0

! Number of species

n_ions = nsm - 1
! These will be reset locally in the radial loop
neo_z_in(3) = MAXVAL(ZIM1(1:nrho))
neo_z_in(4) = MAXVAL(ZIM2(1:nrho))
neo_z_in(5) = MAXVAL(ZIM3(1:nrho))

if (neo_z_in(5) .ge. 1.) n_ions = 4
if (neo_z_in(5) .lt. 1.) n_ions = 3
if (neo_z_in(4) .lt. 1.) n_ions = 2
if (neo_z_in(5) .ge. 1. .and. n_ions .eq. 2) then 
    n_ions = 3
endif
if (neo_z_in(3) .lt. 1.) n_ions = 1
if (neo_z_in(4) .ge. 1. .and. n_ions .eq. 1) then 
    n_ions = 2
endif

! local field averages
! Initialise to zero for non-calculated species

neo_dens_in   = 0.
neo_temp_in   = 0.
neo_dlnndr_in = 0.
neo_dlntdr_in = 0.

radial_loop: do jradial=1, n_radial

    path_in='./'
    call neo_init_serial(path_in)

    j0 = jjgrid(jradial)

!thermal impurities

    neo_z_in(3) = max(1., ZIM1(j0))
    neo_z_in(4) = ZIM2(j0)
    neo_z_in(5) = ZIM3(j0)

    if (neo_z_in(5) .ge. 1. .and. n_ions .eq. 2) then 
        neo_z_in(4) = neo_z_in(5)
        neo_mass_in(4) = neo_mass_in(4)
        ni_m(3, :) = ni_m(4, :)
        ti_m(3, :) = ti_m(4, :)
    endif
    if (neo_z_in(4) .ge. 1. .and. n_ions .eq. 1) then 
        neo_z_in(3) = neo_z_in(4)
        neo_mass_in(3) = neo_mass_in(4)
        ni_m(2, :) = ni_m(3, :)
        ti_m(2, :) = ti_m(3, :)
    endif

!    n_ions = 2 ! git brute force

    neo_mass_in(n_ions+2:) = 0.
    neo_z_in   (n_ions+2:) = 0.

! Selected grid for TGLF computation

    rho_tg(jradial) = rho_m(j0)

! Ref variables for normalisation

    anorm = AMETR(nrho)
    tnorm = TE(j0) 
    nnorm = NE(j0)
    mnorm = AMJ*mpp
    vnorm = sqrt(e00*1.e3*tnorm/mnorm)
    T0 = tnorm*1.e3
    cs0 = vnorm       ! thermal velocity unit m/sec

    if (verbose) then
       write(*, '(A, 5e11.4)') 'Norm', anorm, tnorm, nnorm, mnorm, vnorm
    endif

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
    if (j0 == NA1E .and. NA1E /= nrho) then
        dte = 0.5*dstep*(TE(j0) - TE(j02)) ! Left derivative
    else
        dte = dstep*(TE(j01) - TE(j02))
    endif
    if (j0 == NA1N .and. NA1N /= nrho) then
        dne = 0.5*dstep*(NE(j0) - NE(j02)) ! Left derivative
    else
        dne = dstep*(NE(j01) - NE(j02))
    endif
    dq    = dstep*(q_exp(j01) - q_exp(j02))
    dvpar = dstep*(vpar_m(j01) - vpar_m(j02))
    dvper = dstep*(vper_m(j01) - vper_m(j02))
    do jspec=1, n_ions
       if (j0 == NA1I .and. NA1I /= nrho) then
            dti(jspec) = 0.5*dstep*(ti_m(jspec, j0) - ti_m(jspec, j02))
            dni(jspec) = 0.5*dstep*(ni_m(jspec, j0) - ni_m(jspec, j02))
        else
            dti(jspec) = dstep*(ti_m(jspec, j01) - ti_m(jspec, j02))
            dni(jspec) = dstep*(ni_m(jspec, j01) - ni_m(jspec, j02))
        endif
    enddo
    dr = drmin/anorm    ! gradients w.r.t. minor radius even for s-alpha geometry
    drhodr = drho/drmin

! Log derivatives
    neo_dens_in(1) = 1.
    neo_temp_in(1) = 1.
    neo_dlnndr_in(1) = -dne/(dr*NE(j0))
    neo_dlntdr_in(1) = -dte/(dr*TE(j0))

    do i_ion=1, n_ions
        neo_dens_in(i_ion+1) = ni_m(i_ion, j0)/NE(j0)
        neo_temp_in(i_ion+1) = ti_m(i_ion, j0)/TE(j0)
        neo_dlnndr_in(i_ion+1) = -dni(i_ion)/(dr*ni_m(i_ion, j0))
        neo_dlntdr_in(i_ion+1) = -dti(i_ion)/(dr*ti_m(i_ion, j0))
    enddo
! Restore quasi-neutrality via main ions?

! neo model parameters
    neo_sim_model_in = 2  ! type of NEO calculation: 1 analytic, 2 kinetic
    neo_equilibrium_model_in = 2
    neo_silent_flag_in = 0 ! DUmp file for stand-alone
    neo_test_flag_in = 0

! Resolution 
    neo_n_energy_in = 5  ! number of energy points
    neo_n_xi_in     = 17 ! number of xi points
    neo_n_theta_in  = 29 ! number of theta points
    neo_ipccw_in = -1
    neo_btccw_in = -1
    neo_n_species_in = n_ions + 1

    Bunit = BTOR*drhodr*rho(j0)/AMETR(j0)
! Miller geometry magnetic field unit
    omega0 = e00*Bunit/mnorm    ! gyrofrequency unit 1/sec
    rhos0 = cs0/omega0      ! gyroradius unit cm

    neo_epar0_in = epar0_in(j0)*BTOR/Bunit
    neo_rho_star_in  = rhos0/anorm

!---------
! Species
!---------
  
! Electrons
    lnlamda = 24.0 -0.5*LOG(nnorm*1.e13)+LOG(T0)       
    taue = (3.44E5)*((T0)**1.5)/(nnorm*1e13*lnlamda)  !  sec
    xnuei = 0.75*SQRT(pi)/taue             ! 1/sec 
    neo_nu_1_in = xnuei*anorm/cs0          ! normalized electron-ion collision 

! Geometry
    neo_rmin_over_a_in = AMETR(j0)/anorm
    neo_rmaj_over_a_in = rmaj_exp(j0)/anorm
    neo_q_in           = q_exp(j0)
    neo_shear_in       = AMETR(j0)*dq/(drmin*q_exp(j0))
    neo_shift_in       = drmaj/drmin
    neo_kappa_in       = ELON(j0)
    neo_s_kappa_in     = AMETR(j0)*delong/(drmin*ELON(j0))

    neo_delta_in       = TRIA(j0)
    neo_s_delta_in     = AMETR(j0)*dtrian/drmin

    if (verbose) then
       write(*, '(A, 8f8.4, e11.4)') 'GEO', neo_rmin_over_a_in, neo_rmaj_over_a_in, neo_q_in, &
          neo_shear_in, neo_shift_in, neo_kappa_in, neo_s_kappa_in, &
          neo_delta_in,  neo_s_delta_in
    endif

! Rotation is always *on* in NEO.
! COORDINATES: The signs of all rotation-related quantities below are 
! inherited (unchanged) from input.profiles.  In general NEO expects 
! these to be correctly signed/oriented.
! no rotation at the moment

    neo_rotation_model_in = 2
    neo_omega_rot_in       =  anorm*vpar_m(j0)/(rmaj_exp(j0) * cs0)
    neo_omega_rot_deriv_in = -anorm*dvpar/(drmaj*cs0)
    neo_rmin_over_a_2_in = neo_rmin_over_a_in ! used only for global runs

    Gamma_neo_GB = anorm*vnorm
    Q_neo_GB     = anorm*vnorm
    Pi_neo_GB    = anorm*vnorm
    Jpar_GB      = e00*vnorm*nnorm*1.e19*Bunit/BTOR/1.e6

    pflux_i_neo(:) = 0.0
    pflux_e_neo    = 0.0
    eflux_i_neo(:) = 0.0
    eflux_e_neo    = 0.0

! derived units for the plasma

    SELECT CASE (neo_sim_model_in) 

    CASE(1) ! analytic
        write(*,*) 'run neo analytic'
        call neo_run
        pflux_i_neo(1) = neo_pflux_thHH_out *Gamma_neo_GB
        eflux_i_neo(1) = neo_eflux_thCHi_out*Q_neo_GB
        pflux_e_neo    = neo_pflux_thHH_out *Gamma_neo_GB 
        eflux_e_neo    = neo_eflux_thHHe_out*Q_neo_GB
        jboots = neo_jpar_thS_out*Jpar_GB
        write(*,*) 'end neo analytic', pflux_i_neo(1), eflux_i_neo(1), pflux_e_neo, eflux_e_neo, jboots

    CASE(2) ! kinetic calculation
        write(*,*) 'run neo DKE', n_ions
        call neo_run

        pflux_e_neo = (neo_pflux_dke_out(1) + tgyro_neo_gv_flag*neo_pflux_gv_out(1)) * Gamma_neo_GB
        eflux_e_neo = (neo_efluxncv_dke_out(1) + tgyro_neo_gv_flag*neo_efluxncv_gv_out(1)) * Q_neo_GB

        do i_ion=1, n_ions
            pflux_i_neo(i_ion) = (neo_pflux_dke_out(i_ion+1) + &
                tgyro_neo_gv_flag*neo_pflux_gv_out(i_ion+1))*Gamma_neo_GB/neo_dens_in(i_ion)
            eflux_i_neo(i_ion) = (neo_efluxncv_dke_out(i_ion+1) + &
                tgyro_neo_gv_flag*neo_efluxncv_gv_out(i_ion+1))*Q_neo_GB
            vpflux_neo(i_ion) = neo_vpol_dke_out(i_ion+1)*vnorm  !poloidal flow on outboard mid-plane of main ions
            vtflux_neo(i_ion) = neo_vtor_dke_out(i_ion+1)*vnorm  !toroidal flow on outboard mid-plane of main ions
        enddo

        jboots = neo_jpar_dke_out*Jpar_GB
        write(*, '(A, 2e11.4)') 'end neo DKE', drhodr**2 *pflux_e_neo, drhodr**2 *eflux_e_neo

    END SELECT

    particle_flux(1, 1) = drhodr**2 *pflux_e_neo
    energy_flux  (1, 1) = drhodr**2 *eflux_e_neo
! saves gv on 2nd index
    particle_flux(1, 2) = drhodr**2 *neo_pflux_gv_out(2)*Gamma_neo_GB
    energy_flux  (1, 2) = drhodr**2 *neo_efluxncv_gv_out(2)*Q_neo_GB
    do i_ion=1,n_ions
        particle_flux(i_ion+1, 1) = drhodr**2 *pflux_i_neo(i_ion) !impurity particle flux
        energy_flux  (i_ion+1, 1) = drhodr**2 *eflux_i_neo(i_ion)   !impurity energy flux
        particle_flux(i_ion+1, 2) = drhodr**2 *neo_pflux_gv_out(i_ion+1)*Gamma_neo_GB
        energy_flux  (i_ion+1, 2) = drhodr**2 *neo_efluxncv_gv_out(i_ion+1)*Q_neo_GB
    enddo

    jbs(jradial) = jboots

! END call_ganeo

! Transport coefficients

    ion_eflux = SUM(energy_flux(2: n_ions+1, 1))
    chii(jradial) = ion_eflux        /(1e-4 + abs(neo_dlntdr_in(2)))
    chie(jradial) = energy_flux(1, 1)/(1e-4 + abs(neo_dlntdr_in(1)))
    elec_pflux(jradial) = particle_flux(1, 1)/drhodr         ! particle flux
    vippd(jradial)  = vpflux_neo(1)
    vittd(jradial)  = vtflux_neo(1)
    vippi1(jradial) = vpflux_neo(2)
    vitti1(jradial) = vtflux_neo(2)

enddo radial_loop

call qinterp(rho_tg, chii  , nradial, rho_m(jr_min:jr_max), chii_m  (jr_min:jr_max), jr_max-jr_min+1)
call qinterp(rho_tg, chie  , nradial, rho_m(jr_min:jr_max), chie_m  (jr_min:jr_max), jr_max-jr_min+1)
call qinterp(rho_tg, vippd , nradial, rho_m(jr_min:jr_max), vippd_m (jr_min:jr_max), jr_max-jr_min+1)
call qinterp(rho_tg, vittd , nradial, rho_m(jr_min:jr_max), vittd_m (jr_min:jr_max), jr_max-jr_min+1)
call qinterp(rho_tg, vippi1, nradial, rho_m(jr_min:jr_max), vippi1_m(jr_min:jr_max), jr_max-jr_min+1)
call qinterp(rho_tg, vitti1, nradial, rho_m(jr_min:jr_max), vitti1_m(jr_min:jr_max), jr_max-jr_min+1)
call qinterp(rho_tg, jbs   , nradial, rho_m(jr_min:jr_max), j_boot  (jr_min:jr_max), jr_max-jr_min+1)
call qinterp(rho_tg, elec_pflux, nradial, rho_m(jr_min:jr_max), elec_pflux_m(jr_min:jr_max), jr_max-jr_min+1)

chii_m  (1:2) = chii_m  (3)
chie_m  (1:2) = chie_m  (3)
vippd_m (1:2) = vippd_m (3)
vittd_m (1:2) = vittd_m (3)
vippi1_m(1:2) = vippi1_m(3)
vitti1_m(1:2) = vitti1_m(3)
j_boot  (1:2) = j_boot  (3)
elec_pflux_m(1:2) = elec_pflux_m(3)

DIF(1:nrho) = 0.d0       ! D, electron diffusivity, m^2/s
DPH(1:nrho) = 0.d0       ! D, impurity diffusivity, m^2/s
DPL(1:nrho) = 0.d0       ! impurity convection
GM2(1:nrho) = 0.d0
OM1(1:nrho) = 0.d0
OM2(1:nrho) = 0.d0
FR1(1:nrho) = 0.d0

do j=jr_min, jr_max
    CHI(j) = chii_m(j)/gradrhosq_exp(j) ! \chi_i, m^2/s : starts from work(21,:) 
    CHE(j) = chie_m(j)/gradrhosq_exp(j) ! \chi_e, m^2/s
    VIN(j) = elec_pflux_m(j)/anorm/gradrhosq_exp(j)
    DPR(j) = j_boot(j)    ! bootstrap current
    XTB(j) = vippd_m(j)   ! main ions poloidal flow
    EGM(j) = vippi1_m(j)  ! 1st imp poloidal flow
    GAM(j) = vittd_m(j)   ! main ions toroidal flow
    GM1(j) = vitti1_m(j)  ! 1st imp toroidal flow
enddo

return
END subroutine neo_interf
