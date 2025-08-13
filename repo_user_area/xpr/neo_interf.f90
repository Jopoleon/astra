subroutine neo_interf(jr1_in, jr2_in, n_inputs, n_outputs, nrho, &
    nspec_max, ns_in, BTOR, RTOR, ABC, AMJ, ZMJ, AIM1, AIM2, AIM3, &
    inputs, &
! output
    outputs)

use neo_interface, only: neo_mass_in, neo_z_in, neo_dens_in, neo_temp_in, &
    neo_dlnndr_in, neo_dlntdr_in, neo_sim_model_in, neo_equilibrium_model_in, &
    neo_silent_flag_in, neo_test_flag_in, neo_n_energy_in, neo_n_xi_in, &
    neo_n_theta_in, neo_ipccw_in, neo_btccw_in, neo_n_species_in, &
    neo_epar0_in, neo_rho_star_in, neo_nu_1_in, neo_q_in, neo_shear_in, &
    neo_rmin_over_a_in, neo_rmaj_over_a_in, neo_rmin_over_a_2_in, &
    neo_shift_in, neo_kappa_in, neo_s_kappa_in, neo_delta_in, neo_s_delta_in, &
    neo_rotation_model_in, neo_omega_rot_in, neo_omega_rot_deriv_in, &
    neo_pflux_thHH_out, neo_eflux_thCHi_out, neo_eflux_thHHe_out, &
    neo_jpar_thS_out, neo_pflux_dke_out, neo_pflux_dke_out, neo_efluxncv_dke_out, &
    neo_vpol_dke_out, neo_vtor_dke_out, neo_jpar_dke_out, &
    neo_pflux_gv_out, neo_efluxncv_gv_out

implicit none

logical, parameter :: verbose=.False.
integer, parameter :: jpd=700, nradial=5

double precision, parameter :: &
   e00  = 1.6020e-19, &       ! elementary charge (C)
   mpp  = 1.6726E-27, &       ! proton mass (kg)
   pi   = 3.141592653589793

integer, intent(in) :: n_inputs, n_outputs, nrho, jr1_in, jr2_in, nspec_max, ns_in

double precision, intent(in) :: BTOR, RTOR, ABC, AMJ, AIM1, AIM2, AIM3, ZMJ
double precision, intent(in), dimension(n_inputs, nrho) :: inputs
double precision, intent(out), dimension(jr2_in + 1 - jr1_in, n_outputs) :: outputs

!----------------------------------------------------------------------
integer :: i_ion, n_ions, j, jr, jspec, chunk
double precision :: ion_eflux, drhodr_sq
double precision :: Bunit, cs0, rhos0, omega0, lnlamda, taue, xnuei
double precision :: T0, anorm, mnorm, tnorm, nnorm, vnorm, &
   pflux_e, eflux_e, jboots, tgyro_gv_flag, &
   Gamma_GB, Q_GB, Pi_GB, Jpar_GB

double precision, dimension(jr2_in+1-jr1_in) :: chie, chii, elec_pflux, rho, &
    ametr, rmaj, elon, tria, q, ne, te, ti, vpar, &
    vippd, vittd, vippi1, vitti1, jbs, epar0_in, &
    drmin, drmaj, delong, dtrian, dr, dne, dte, dti, dq, dvpar, drhodr
double precision, dimension(ns_in, jr2_in+1-jr1_in) :: ni, zimp, dni
double precision, dimension(ns_in) :: pflux_i, eflux_i, vpflux, vtflux
double precision, dimension(ns_in, 2):: energy_flux, particle_flux
character(len=80) :: path_in

chunk = jr2_in + 1 - jr1_in

!-----------------
! Get input arrays
!-----------------
rho   = inputs( 1, jr1_in:jr2_in)
ametr = inputs( 2, jr1_in:jr2_in)
rmaj  = inputs( 3, jr1_in:jr2_in)
elon  = inputs( 4, jr1_in:jr2_in)
tria  = inputs( 5, jr1_in:jr2_in)
q     = inputs( 6, jr1_in:jr2_in)
ti    = inputs( 7, jr1_in:jr2_in)
te    = inputs( 8, jr1_in:jr2_in)
ne    = inputs( 9, jr1_in:jr2_in)
vpar  = inputs(10, jr1_in:jr2_in)
ni(1, :)   = inputs(11, jr1_in:jr2_in)
ni(2, :)   = inputs(12, jr1_in:jr2_in)
ni(3, :)   = inputs(13, jr1_in:jr2_in)
ni(4, :)   = inputs(14, jr1_in:jr2_in)
zimp(1, :) = inputs(15, jr1_in:jr2_in)
zimp(2, :) = inputs(16, jr1_in:jr2_in)
zimp(3, :) = inputs(17, jr1_in:jr2_in)
drmin  = inputs(18, jr1_in:jr2_in)
drmaj  = inputs(19, jr1_in:jr2_in)
delong = inputs(20, jr1_in:jr2_in)
dtrian = inputs(21, jr1_in:jr2_in)
dti    = inputs(22, jr1_in:jr2_in)
dte    = inputs(23, jr1_in:jr2_in)
dne    = inputs(24, jr1_in:jr2_in)
dq     = inputs(25, jr1_in:jr2_in)
dvpar  = inputs(26, jr1_in:jr2_in)
dr     = inputs(27, jr1_in:jr2_in)
drhodr = inputs(28, jr1_in:jr2_in)
dni(1, :) = inputs(29, jr1_in:jr2_in)
dni(2, :) = inputs(30, jr1_in:jr2_in)
dni(3, :) = inputs(31, jr1_in:jr2_in)
dni(4, :) = inputs(32, jr1_in:jr2_in)

anorm = ABC
neo_mass_in(1) = 5.4447e-4
neo_mass_in(2) = AMJ
neo_mass_in(3) = AIM1
neo_mass_in(4) = AIM2
neo_mass_in(5) = AIM3
neo_mass_in = neo_mass_in/AMJ
neo_z_in(1) = -1.
neo_z_in(2) = ZMJ

mnorm = AMJ*mpp

tgyro_gv_flag = 0.
epar0_in = 0.

! Number of species
n_ions = ns_in - 1

neo_dens_in   = 0.
neo_temp_in   = 0.
neo_dlnndr_in = 0.
neo_dlntdr_in = 0.

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

if (jr1_in == 1) print*, 'Run NEO', n_ions

radial_loop: do jr=1, chunk

!thermal impurities

    neo_z_in(3) = max(1., zimp(1, jr))
    neo_z_in(4) = zimp(2, jr)
    neo_z_in(5) = zimp(3, jr)

    if (neo_z_in(5) >= 1. .and. ns_in == 3) then
        neo_z_in(4) = neo_z_in(5)
        neo_mass_in(4) = neo_mass_in(5)
        ni(3, :) = ni(4, :)
    endif
    if (neo_z_in(4) >= 1. .and. ns_in == 2) then
        neo_z_in(3) = neo_z_in(4)
        neo_mass_in(3) = neo_mass_in(4)
        ni(2, :) = ni(3, :)
    endif

    neo_mass_in(n_ions+2:) = 0.
    neo_z_in   (n_ions+2:) = 0.

    path_in = 'neo/'
    call neo_init_serial(path_in)

! Ref variables for normalisation

    tnorm = te(jr)
    nnorm = ne(jr)
    vnorm = sqrt(e00*1.e3*tnorm/mnorm)
    T0  = tnorm*1.e3
    cs0 = vnorm       ! thermal velocity unit m/sec

    neo_dens_in(1) = 1.
    neo_temp_in(1) = 1.
    neo_dlnndr_in(1) = -dne(jr)/(dr(jr)*ne(jr))
    neo_dlntdr_in(1) = -dte(jr)/(dr(jr)*te(jr))

    do i_ion=1, n_ions
        neo_dens_in(i_ion+1) = ni(i_ion, jr)/ne(jr)
        neo_temp_in(i_ion+1) = ti(jr)/te(jr)
        neo_dlnndr_in(i_ion+1) = -dni(i_ion, jr)/(dr(jr)*ni(i_ion, jr))
        neo_dlntdr_in(i_ion+1) = -dti(jr)/(dr(jr)*ti(jr))
    enddo

    Bunit = BTOR*drhodr(jr)*rho(jr)/ametr(jr)
! Miller geometry magnetic field unit
    omega0 = e00*Bunit/mnorm    ! gyrofrequency unit 1/sec
    rhos0 = cs0/omega0      ! gyroradius unit cm

    neo_epar0_in = epar0_in(jr)*BTOR/Bunit
    neo_rho_star_in = rhos0/anorm

! Collisionalities
    
    lnlamda = 24.0 - 0.5*LOG(nnorm*1.e13) + LOG(T0)       
    taue = 3.44E5*(T0**1.5)/(nnorm*1e13*lnlamda)  !  sec
    xnuei = 0.75*SQRT(pi)/taue             ! 1/sec 
    neo_nu_1_in = xnuei*anorm/cs0          ! normalized electron-ion collision 

! Geometry
    neo_rmin_over_a_in = ametr(jr)/anorm
    neo_rmaj_over_a_in = rmaj(jr)/anorm
    neo_q_in           = q(jr)
    neo_shear_in       = ametr(jr)*dq(jr)/(drmin(jr)*q(jr))
    neo_shift_in       = drmaj(jr)/drmin(jr)
    neo_kappa_in       = elon(jr)
    neo_s_kappa_in     = ametr(jr)*delong(jr)/(drmin(jr)*elon(jr))
    neo_delta_in       = tria(jr)
    neo_s_delta_in     = ametr(jr)*dtrian(jr)/drmin(jr)

! Rotation is always *on* in NEO.
! COORDINATES: The signs of all rotation-related quantities below are 
! inherited (unchanged) from input.profiles.  In general NEO expects 
! these to be correctly signed/oriented.
! no rotation at the moment

    neo_rotation_model_in = 2
    neo_omega_rot_in       =  anorm*vpar(jr)/(rmaj(jr) * cs0)
    neo_omega_rot_deriv_in = -anorm*dvpar(jr)/(drmaj(jr)*cs0)
    neo_rmin_over_a_2_in = neo_rmin_over_a_in ! used only for global runs

    Gamma_GB = anorm*vnorm
    Q_GB     = anorm*vnorm
    Pi_GB    = anorm*vnorm
    Jpar_GB      = e00*vnorm*nnorm*1.e19*Bunit/BTOR/1.e6

    pflux_i = 0.0
    pflux_e = 0.0
    eflux_i = 0.0
    eflux_e = 0.0
    drhodr_sq = drhodr(jr)**2
! derived units for the plasma

    SELECT CASE (neo_sim_model_in)

    CASE(1) ! analytic
        call neo_run
        pflux_i(1) = neo_pflux_thHH_out *Gamma_GB
        eflux_i(1) = neo_eflux_thCHi_out*Q_GB
        pflux_e    = neo_pflux_thHH_out *Gamma_GB 
        eflux_e    = neo_eflux_thHHe_out*Q_GB
        jboots = neo_jpar_thS_out*Jpar_GB
        print*, 'neo analytic', pflux_i(1), eflux_i(1)

    CASE(2) ! kinetic calculation
        call neo_run

        pflux_e = (neo_pflux_dke_out(1)    + tgyro_gv_flag*neo_pflux_gv_out(1)) *Gamma_GB
        eflux_e = (neo_efluxncv_dke_out(1) + tgyro_gv_flag*neo_efluxncv_gv_out(1)) * Q_GB

        do i_ion=1, n_ions
            pflux_i(i_ion) = (neo_pflux_dke_out(i_ion+1) + &
                tgyro_gv_flag*neo_pflux_gv_out(i_ion+1))*Gamma_GB/neo_dens_in(i_ion)
            eflux_i(i_ion) = (neo_efluxncv_dke_out(i_ion+1) + &
                tgyro_gv_flag*neo_efluxncv_gv_out(i_ion+1))*Q_GB
            vpflux(i_ion) = neo_vpol_dke_out(i_ion+1)*vnorm  !poloidal flow on outboard mid-plane of main ions
            vtflux(i_ion) = neo_vtor_dke_out(i_ion+1)*vnorm  !toroidal flow on outboard mid-plane of main ions
        enddo

        jboots = neo_jpar_dke_out*Jpar_GB
        write(*, '(A, 2e11.4)') 'neo DKE', pflux_e, eflux_e

    END SELECT

    particle_flux(1, 1) = drhodr_sq * pflux_e
    energy_flux  (1, 1) = drhodr_sq * eflux_e
! saves gv on 2nd index
    particle_flux(1, 2) = drhodr_sq * neo_pflux_gv_out(2)*Gamma_GB
    energy_flux  (1, 2) = drhodr_sq * neo_efluxncv_gv_out(2)*Q_GB
    do i_ion=1,n_ions
        particle_flux(i_ion+1, 1) = drhodr_sq * pflux_i(i_ion) ! impurity particle flux
        energy_flux  (i_ion+1, 1) = drhodr_sq * eflux_i(i_ion) ! impurity energy flux
        particle_flux(i_ion+1, 2) = drhodr_sq * neo_pflux_gv_out(i_ion+1)*Gamma_GB
        energy_flux  (i_ion+1, 2) = drhodr_sq * neo_efluxncv_gv_out(i_ion+1)*Q_GB
    enddo

    jbs(jr) = jboots

! END call_ganeo

! Transport coefficients

    ion_eflux = SUM(energy_flux(2: n_ions+1, 1))
    chii(jr) = ion_eflux        /(1e-4 + abs(neo_dlntdr_in(2)))
    chie(jr) = energy_flux(1, 1)/(1e-4 + abs(neo_dlntdr_in(1)))
    elec_pflux(jr) = particle_flux(1, 1)/drhodr(jr)         ! particle flux
    vippd(jr)  = vpflux(1)
    vittd(jr)  = vtflux(1)
    vippi1(jr) = vpflux(2)
    vitti1(jr) = vtflux(2)

enddo radial_loop

! Simulated NEO computation:
outputs(:, 1) = chii
outputs(:, 2) = chie
outputs(:, 3) = jbs
outputs(:, 4) = elec_pflux
outputs(:, 5) = vippd
outputs(:, 6) = vittd
outputs(:, 7) = vippi1
outputs(:, 8) = vitti1

end subroutine neo_interf
