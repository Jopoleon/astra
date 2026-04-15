subroutine feqis_main(nucoils, ucoils, neql, k_fixfree, no_circuit_eq, ifplasma, machine_name, & 
    equil_in, equil_out)

use imas_ids, only: type_equilibrium  
use circuit, only: psi_cur_old, psiplasmatoconduc, voltage, &
    psi_mutual_effect_conductors_simple
use fbe_core, only: nr2, nz2, nr1, nz1, Rrect, Zrect, &
    psirz, psiextrz, &
    psi_external_calc, nbnd, i_plasmatype, &
    r_xpoint, z_xpoint, n_of_xpoints, active_x_point
use pbe_core, only: nrho, ntheta, psibndp, psiaxisp
use feqis_scalars, only: psplex, iplasma
use transport2fbe, only: refit_mode, simple_plasma_model_breakdown, &
    plasma_config, x_point_save
use feqis_solvers, only: feqis_init, equil_init_circ, fix_boundary, &
    full_system_advance, convert_boundary_to_pbe, circuit_eq_advance, &
    equil_assignments, solve_gse2d_fbe_full

implicit none

integer, intent(in) :: nucoils, ifplasma, neql, k_fixfree, no_circuit_eq
double precision, intent(in) , dimension(nucoils) :: ucoils
character(len=*) :: machine_name
type(type_equilibrium), intent(in) :: equil_in
type(type_equilibrium), intent(out) :: equil_out

integer :: nrplasma, j_init, j_call, j_vacplas

data j_call/0/
data j_vacplas/0/
data j_init/0/
save j_call, j_init, j_vacplas

call feqis_init(equil_in, neql, k_fixfree, j_init, ifplasma)

nrplasma = nrho
nbnd = ntheta ! for fbe

if (j_call == 0) then
    nr2 = SIZE(equil_in%eqgeometry%rectgrid%r2d)
    if (nr2 > 0) then
        nz2 = SIZE(equil_in%eqgeometry%rectgrid%z2d)
        nr1 = nr2 - 1
        nz1 = nz2 - 1
        if (.not. allocated(Rrect)) then
            allocate(Rrect(nr2))
            allocate(Zrect(nz2))
            allocate(psirz(nr2, nz2))
        endif
        Rrect = equil_in%eqgeometry%rectgrid%r2d
        Zrect = equil_in%eqgeometry%rectgrid%z2d
    endif
    if (k_fixfree == 1 .or. refit_mode == 818) then   ! also if refit mode = 818, initialize free boundary stuff
        call equil_init_circ(machine_name)
        if (refit_mode == 818) then
            j_call = 1 ! this is because if free boundary was never called, it needs to initialize these arrays
            refit_mode = 0 ! this is because if free boundary was never called, it needs to initialize these arrays
        endif
    endif
endif

if (ifplasma == 1) then
    allocate(equil_out%eqgeometry%boundary%r(ntheta))
    allocate(equil_out%eqgeometry%boundary%z(ntheta))
    allocate(equil_out%eqgeometry%rectgrid%r2d(nr2))
    allocate(equil_out%eqgeometry%rectgrid%z2d(nz2))
    allocate(equil_out%eqgeometry%rectgrid%psirz2d(nr2, nz2))
    allocate(equil_out%eqgeometry%rectgrid%fdia2d(nr2, nz2))

    allocate(equil_out%coord_sys%position%r(nrplasma, ntheta))
    allocate(equil_out%coord_sys%position%z(nrplasma, ntheta))    
    allocate(equil_out%coord_sys%position%theta2d(ntheta))    
    allocate(equil_out%coord_sys%position%rmin(nrplasma, ntheta))    
    allocate(equil_out%coord_sys%position%psirz(nrplasma, ntheta))    
    allocate(equil_out%coord_sys%gradvcell(nrplasma, ntheta))
    allocate(equil_out%coord_sys%bpcell(nrplasma, ntheta))
    allocate(equil_out%coord_sys%bcell(nrplasma, ntheta))
    allocate(equil_out%coord_sys%rcell(nrplasma, ntheta))
    allocate(equil_out%coord_sys%darea(nrplasma, ntheta))
    allocate(equil_out%coord_sys%jphi(nrplasma, ntheta))
    
    allocate(equil_out%profiles_1d%psi(nrplasma))
    allocate(equil_out%profiles_1d%pressure(nrplasma))
    allocate(equil_out%profiles_1d%phi(nrplasma))
    allocate(equil_out%profiles_1d%pprime(nrplasma))
    allocate(equil_out%profiles_1d%ffprime(nrplasma))
    allocate(equil_out%profiles_1d%F_dia(nrplasma))
    allocate(equil_out%profiles_1d%q(nrplasma))
    allocate(equil_out%profiles_1d%gm1(nrplasma))
    allocate(equil_out%profiles_1d%gm4(nrplasma))
    allocate(equil_out%profiles_1d%gm5(nrplasma))
    allocate(equil_out%profiles_1d%gm41(nrplasma))
    allocate(equil_out%profiles_1d%rbp_b2(nrplasma))
    allocate(equil_out%profiles_1d%bplfs(nrplasma))
    allocate(equil_out%profiles_1d%rho_tor_norm(nrplasma) )
    allocate(equil_out%profiles_1d%jparallel(nrplasma) )
    allocate(equil_out%profiles_1d%sigmapar%value(nrplasma) )
    allocate(equil_out%profiles_1d%jni%value(nrplasma) )
    allocate(equil_out%profiles_1d%te%value(nrplasma) )
    allocate(equil_out%profiles_1d%acosB2a(nrplasma, 5))
    allocate(equil_out%profiles_1d%asinB2a(nrplasma, 5))
    allocate(equil_out%profiles_1d%acosBlnBa(nrplasma, 5))
    allocate(equil_out%profiles_1d%asinBlnBa(nrplasma, 5))
    allocate(equil_out%profiles_1d%g1(nrplasma))
    allocate(equil_out%profiles_1d%g2(nrplasma))
    allocate(equil_out%profiles_1d%g2int(nrplasma))
    allocate(equil_out%profiles_1d%fofb(nrplasma))
    allocate(equil_out%profiles_1d%areat(nrplasma))
    allocate(equil_out%profiles_1d%perim(nrplasma))
    allocate(equil_out%profiles_1d%ggradro(nrplasma))
    allocate(equil_out%profiles_1d%gdroda(nrplasma))
    allocate(equil_out%profiles_1d%bmaxt(nrplasma))
    allocate(equil_out%profiles_1d%bmint(nrplasma))
    allocate(equil_out%profiles_1d%bdb0(nrplasma))
    allocate(equil_out%profiles_1d%dPSIdV(nrplasma))
    allocate(equil_out%profiles_1d%surface(nrplasma))
    allocate(equil_out%profiles_1d%volume(nrplasma))
    allocate(equil_out%profiles_1d%r_inboard(nrplasma))
    allocate(equil_out%profiles_1d%r_outboard(nrplasma))
    allocate(equil_out%profiles_1d%elongation(nrplasma))
    allocate(equil_out%profiles_1d%tria_upper(nrplasma))
    allocate(equil_out%profiles_1d%tria_lower(nrplasma))
    allocate(equil_out%profiles_1d%shif(nrplasma))
    allocate(equil_out%profiles_1d%shiv(nrplasma))
    allocate(equil_out%profiles_1d%squareness(nrplasma))
endif

if (k_fixfree == 1 .and. refit_mode /= 818) then !any other mode than 818
    voltage = 0.
    voltage(1:nucoils) = ucoils(1:nucoils) ! voltage inputs for active conductors
    if (ifplasma == 0) then       ! only circuit equations solved

        if (simple_plasma_model_breakdown == 0) then
         psi_cur_old = 0.
         psiplasmatoconduc = 0.
        else
         iplasma = equil_in%global_param%i_plasma/1.e6 	
         call psi_mutual_effect_conductors_simple(psiplasmatoconduc)
        endif
        call circuit_eq_advance(j_call)
        
	psi_cur_old = psiplasmatoconduc
        
        call psi_external_calc()
        psirz = psiextrz
        j_vacplas = 0
    else if (ifplasma == 1) then  ! full plasma solved
        if (j_vacplas == 0) j_call = 0
        call full_system_advance(j_call, no_circuit_eq)
        if (j_call == -1) then
            call convert_boundary_to_pbe()
            call fix_boundary(1)
        endif
        j_vacplas = 1
! below 2 diagnostics for configuration and xpoints
        plasma_config = i_plasmatype ! assign configuration of plasma to external variables
        if (n_of_xpoints >= 1) then
            x_point_save(1:min(20, n_of_xpoints), 1) = r_xpoint(1:min(20, n_of_xpoints))
            x_point_save(1:min(20, n_of_xpoints), 2) = z_xpoint(1:min(20, n_of_xpoints))
            x_point_save(20, 1) = r_xpoint(active_x_point)
            x_point_save(20, 2) = z_xpoint(active_x_point)
        endif
    endif
else if (k_fixfree == 0 .and. refit_mode /= 818) then !any other mode than 818
    call fix_boundary(j_init)
    equil_out%global_param%psplex   = psplex
    equil_out%global_param%psibound = psibndp
    equil_out%global_param%psiaxis  = psiaxisp
    call equil_assignments(equil_out)
    j_init = 1
    return
endif

if (refit_mode == 818) then  ! run prescribed boundary but with coil currents fitting in the bakcground
    call psi_external_calc()
    call solve_gse2d_fbe_full(0)
    call fix_boundary(1)
    equil_out%global_param%psplex   = psplex
    equil_out%global_param%psibound = psibndp
    equil_out%global_param%psiaxis  = psiaxisp
    call equil_assignments(equil_out)
    j_call = 1
    j_init = 1
! below 2 diagnostics for configuration and xpoints
    plasma_config = i_plasmatype ! assign configuration of plasma to external variables
    if (n_of_xpoints >= 1) then
        x_point_save(1:min(20,n_of_xpoints), 1) = r_xpoint(1:min(20,n_of_xpoints))
        x_point_save(1:min(20,n_of_xpoints), 2) = z_xpoint(1:min(20,n_of_xpoints))
        x_point_save(20, 1) = r_xpoint(active_x_point)
        x_point_save(20, 2) = z_xpoint(active_x_point)
    endif
    return
endif

j_call = 1
j_init = 1

if (ifplasma == 1) then
    call equil_assignments(equil_out)
endif

end subroutine feqis_main
