program neo_chunk

use mpi
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

integer, parameter :: n_dims=8
double precision, parameter :: &
   e00  = 1.6020e-19, &       ! elementary charge (C)
   mpp  = 1.6726E-27, &       ! proton mass (kg)
   pi   = 3.141592653589793

integer :: ierr, parent, rank, status(MPI_STATUS_SIZE)
integer :: chunk, n_inputs, n_outputs, n_scalars, dims(n_dims)
integer :: jr, i_ion, nsm, n_ions, i1, i2
double precision :: Bunit, cs0, rhos0, omega0, lnlamda, taue, xnuei
double precision :: anorm, mnorm, tnorm, nnorm, vnorm, T0, nt_cs, &
    pflux_e_neo, eflux_e_neo, jboots, tgyro_neo_gv_flag, &
    Gamma_neo_GB, Q_neo_GB, Pi_neo_GB, Jpar_GB
double precision :: AMJ, BTOR
double precision :: ion_eflux, drhodr_sq
double precision, allocatable :: inputs(:, :), output(:, :), scalars(:)
double precision, allocatable, dimension(:) :: chie, chii, elec_pflux, rho_neo, &
    ametr_neo, rmaj_neo, elon_neo, tria_neo, q_neo, ne_neo, te_neo, ti_neo, vpar_neo, &
    vippd, vittd, vippi1, vitti1, jbs, epar0_in, &
    drmin, drmaj, delong, dtrian, dr, dne, dte, dti, dq, dvpar, drhodr
double precision, allocatable, dimension(:, :) :: ni_neo, zimp_neo, dni
double precision, allocatable, dimension(:) :: pflux_i_neo, eflux_i_neo, vpflux_neo, vtflux_neo
double precision, allocatable, dimension (:, :):: energy_flux, particle_flux
character(len=80) :: path_in

call MPI_Init(ierr)
call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
call MPI_Comm_get_parent(parent, ierr)

if (parent == MPI_COMM_NULL) then
    if (rank == 0) then
        print *, "No parent communicator!"
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    endif
endif

! Receive dimensions from parent
call MPI_Recv(dims, n_dims, MPI_INTEGER, 0, 0, parent, status, ierr)
n_scalars = dims(1)
n_inputs  = dims(2)
n_outputs = dims(3)
nsm       = dims(4)
i1 = dims(5)
i2 = dims(6)
chunk = i2 + 1 - i1
allocate(scalars(n_scalars))
allocate(inputs(chunk, n_inputs))
allocate(output(chunk, n_outputs))
allocate( chie(chunk), chii(chunk), elec_pflux(chunk), &
    rho_neo(chunk), ametr_neo(chunk), rmaj_neo(chunk), elon_neo(chunk), tria_neo(chunk), &
    q_neo(chunk), ne_neo(chunk), te_neo(chunk), ti_neo(chunk),vpar_neo(chunk), &
    vippd(chunk), vittd(chunk), vippi1(chunk), vitti1(chunk), jbs(chunk), &
    drmin(chunk), drmaj(chunk), delong(chunk), dtrian(chunk), dr(chunk), &
    dne(chunk), dte(chunk), dti(chunk), dq(chunk), dvpar(chunk), drhodr(chunk) )
allocate( pflux_i_neo(nsm-1), eflux_i_neo(nsm-1), vpflux_neo(nsm), vtflux_neo(nsm) )
allocate( energy_flux(nsm, 2), particle_flux(nsm, 2) )
allocate( ni_neo(nsm-1, chunk), dni(nsm-1, chunk), zimp_neo(nsm-2, chunk))

! Receive TGLF input scalars and profiles from parent
call MPI_Recv(scalars     , n_scalars, MPI_DOUBLE_PRECISION, 0, 0, parent, status, ierr)
call MPI_Recv(inputs, chunk * n_inputs, MPI_DOUBLE_PRECISION, 0, 0, parent, status, ierr)
rho_neo      = inputs(:,  1)
ametr_neo    = inputs(:,  2)
rmaj_neo     = inputs(:,  3)
elon_neo     = inputs(:,  4)
tria_neo     = inputs(:,  5)
q_neo        = inputs(:,  6)
ti_neo       = inputs(:,  7)
te_neo       = inputs(:,  8)
ne_neo       = inputs(:,  9)
vpar_neo     = inputs(:, 10)
ni_neo(1, :) = inputs(:, 11)
ni_neo(2, :) = inputs(:, 12)
ni_neo(3, :) = inputs(:, 13)
ni_neo(4, :) = inputs(:, 14)
zimp_neo(1, :) = inputs(:, 15)
zimp_neo(2, :) = inputs(:, 16)
zimp_neo(3, :) = inputs(:, 17)
drmin       = inputs(:, 18)
drmaj       = inputs(:, 19)
delong      = inputs(:, 20)
dtrian      = inputs(:, 21)
dti         = inputs(:, 22)
dte         = inputs(:, 23)
dne         = inputs(:, 24)
dq          = inputs(:, 25)
dvpar       = inputs(:, 26)
dr          = inputs(:, 27)
drhodr      = inputs(:, 28)
dni(1, :)   = inputs(:, 29)
dni(2, :)   = inputs(:, 30)
dni(3, :)   = inputs(:, 31)
dni(4, :)   = inputs(:, 32)

BTOR  = scalars(2)
anorm = scalars(3)
AMJ   = scalars(5)
neo_mass_in(1: 5) = scalars(4:  8)/AMJ
neo_z_in(1: 5)    = scalars(9: 13)
mnorm = AMJ*mpp

tgyro_neo_gv_flag = 0.
epar0_in = 0.

! Number of species
n_ions = nsm - 1

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

if (i1 == 1) print*, 'Run NEO', n_ions

radial_loop: do jr=1, chunk

    path_in='./neo/'
    call neo_init_serial(path_in)

!thermal impurities
    neo_z_in(3) = max(1., zimp_neo(1, jr))
    neo_z_in(4) = zimp_neo(2, jr)
    neo_z_in(5) = zimp_neo(3, jr)

    if (neo_z_in(5) >= 1. .and. n_ions == 2) then
        neo_z_in(4) = neo_z_in(5)
        neo_mass_in(4) = neo_mass_in(5)
        ni_neo(3, :) = ni_neo(4, :)
    endif
    if (neo_z_in(4) >= 1. .and. n_ions == 1) then
        neo_z_in(3) = neo_z_in(4)
        neo_mass_in(3) = neo_mass_in(4)
        ni_neo(2, :) = ni_neo(3, :)
    endif

    neo_mass_in(n_ions+2:) = 0.
    neo_z_in   (n_ions+2:) = 0.

! Ref variables for normalisation

    tnorm = te_neo(jr)
    nnorm = ne_neo(jr)
    vnorm = sqrt(e00*1.e3*tnorm/mnorm)
    T0  = tnorm*1.e3
    cs0 = vnorm       ! thermal velocity unit m/sec

! Log derivatives
    neo_dens_in(1) = 1.
    neo_temp_in(1) = 1.
    neo_dlnndr_in(1) = -dne(jr)/(dr(jr)*ne_neo(jr))
    neo_dlntdr_in(1) = -dte(jr)/(dr(jr)*te_neo(jr))

    do i_ion=1, n_ions
        neo_dens_in(i_ion+1) = ni_neo(i_ion, jr)/ne_neo(jr)
        neo_temp_in(i_ion+1) = ti_neo(jr)/te_neo(jr)
        neo_dlnndr_in(i_ion+1) = -dni(i_ion, jr)/(dr(jr)*ni_neo(i_ion, jr))
        neo_dlntdr_in(i_ion+1) = -dti(jr)/(dr(jr)*ti_neo(jr))
    enddo

    Bunit = BTOR*drhodr(jr)*rho_neo(jr)/ametr_neo(jr)
! Miller geometry magnetic field unit
    omega0 = e00*Bunit/mnorm    ! gyrofrequency unit 1/sec
    rhos0 = cs0/omega0      ! gyroradius unit cm

    neo_epar0_in = epar0_in(jr)*BTOR/Bunit
    neo_rho_star_in = rhos0/anorm

!---------
! Species
!---------
  
! Electrons
    lnlamda = 24.0 - 0.5*LOG(nnorm*1.e13) + LOG(T0)       
    taue = 3.44E5*(T0**1.5)/(nnorm*1e13*lnlamda)  !  sec
    xnuei = 0.75*SQRT(pi)/taue             ! 1/sec 
    neo_nu_1_in = xnuei*anorm/cs0          ! normalized electron-ion collision 

! Geometry
    neo_rmin_over_a_in = ametr_neo(jr)/anorm
    neo_rmaj_over_a_in = rmaj_neo(jr)/anorm
    neo_q_in           = q_neo(jr)
    neo_shear_in       = ametr_neo(jr)*dq(jr)/(drmin(jr)*q_neo(jr))
    neo_shift_in       = drmaj(jr)/drmin(jr)
    neo_kappa_in       = elon_neo(jr)
    neo_s_kappa_in     = ametr_neo(jr)*delong(jr)/(drmin(jr)*elon_neo(jr))
    neo_delta_in       = tria_neo(jr)
    neo_s_delta_in     = ametr_neo(jr)*dtrian(jr)/drmin(jr)

! Rotation is always *on* in NEO.
! COORDINATES: The signs of all rotation-related quantities below are 
! inherited (unchanged) from input.profiles.  In general NEO expects 
! these to be correctly signed/oriented.
! no rotation at the moment

    neo_rotation_model_in = 2
    neo_omega_rot_in       =  anorm*vpar_neo(jr)/(rmaj_neo(jr) * cs0)
    neo_omega_rot_deriv_in = -anorm*dvpar(jr)/(drmaj(jr)*cs0)
    neo_rmin_over_a_2_in = neo_rmin_over_a_in ! used only for global runs

    Gamma_neo_GB = anorm*vnorm
    Q_neo_GB     = anorm*vnorm
    Pi_neo_GB    = anorm*vnorm
    Jpar_GB      = e00*vnorm*nnorm*1.e19*Bunit/BTOR/1.e6

    pflux_i_neo = 0.0
    pflux_e_neo = 0.0
    eflux_i_neo = 0.0
    eflux_e_neo = 0.0
    drhodr_sq = drhodr(jr)**2
! derived units for the plasma

    SELECT CASE (neo_sim_model_in)

    CASE(1) ! analytic
        call neo_run
        pflux_i_neo(1) = neo_pflux_thHH_out *Gamma_neo_GB
        eflux_i_neo(1) = neo_eflux_thCHi_out*Q_neo_GB
        pflux_e_neo    = neo_pflux_thHH_out *Gamma_neo_GB 
        eflux_e_neo    = neo_eflux_thHHe_out*Q_neo_GB
        jboots = neo_jpar_thS_out*Jpar_GB
        write(*,*) 'end neo analytic', pflux_i_neo(1), eflux_i_neo(1), pflux_e_neo, eflux_e_neo, jboots

    CASE(2) ! kinetic calculation
        call neo_run

        pflux_e_neo = (neo_pflux_dke_out(1)    + tgyro_neo_gv_flag*neo_pflux_gv_out(1)) *Gamma_neo_GB
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
        write(*, '(A, 2e11.4)') 'neo DKE', pflux_e_neo, eflux_e_neo

    END SELECT

    particle_flux(1, 1) = drhodr_sq * pflux_e_neo
    energy_flux  (1, 1) = drhodr_sq * eflux_e_neo
! saves gv on 2nd index
    particle_flux(1, 2) = drhodr_sq * neo_pflux_gv_out(2)*Gamma_neo_GB
    energy_flux  (1, 2) = drhodr_sq * neo_efluxncv_gv_out(2)*Q_neo_GB
    do i_ion=1,n_ions
        particle_flux(i_ion+1, 1) = drhodr_sq * pflux_i_neo(i_ion) ! impurity particle flux
        energy_flux  (i_ion+1, 1) = drhodr_sq * eflux_i_neo(i_ion) ! impurity energy flux
        particle_flux(i_ion+1, 2) = drhodr_sq * neo_pflux_gv_out(i_ion+1)*Gamma_neo_GB
        energy_flux  (i_ion+1, 2) = drhodr_sq * neo_efluxncv_gv_out(i_ion+1)*Q_neo_GB
    enddo

    jbs(jr) = jboots

! END call_ganeo

! Transport coefficients

    ion_eflux = SUM(energy_flux(2: n_ions+1, 1))
    chii(jr) = ion_eflux        /(1e-4 + abs(neo_dlntdr_in(2)))
    chie(jr) = energy_flux(1, 1)/(1e-4 + abs(neo_dlntdr_in(1)))
    elec_pflux(jr) = particle_flux(1, 1)/drhodr(jr)         ! particle flux
    vippd(jr)  = vpflux_neo(1)
    vittd(jr)  = vtflux_neo(1)
    vippi1(jr) = vpflux_neo(2)
    vitti1(jr) = vtflux_neo(2)

enddo radial_loop

! Simulated TGLF computation:
output(:, 1) = chii
output(:, 2) = chie
output(:, 3) = jbs
output(:, 4) = elec_pflux
output(:, 5) = vippd
output(:, 6) = vittd
output(:, 7) = vippi1
output(:, 8) = vitti1

call MPI_Send(output, chunk * n_outputs, MPI_DOUBLE_PRECISION, 0, 1, parent, ierr)
call MPI_Finalize(ierr)

end program neo_chunk
