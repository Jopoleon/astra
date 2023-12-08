subroutine feqis_main(nucoils, ucoils, parameters_equil, ifplasma, & 
    equil_in, equil_out)

use imas_ids, only: type_equilibrium  
use parameters_a2equil, only: type_parameters
use feqis_circuit, only: ncoils, nrho, nteta, nr2, nz2, &
    psi_cur_old, psiplasmatoconduc, psirz, psiextrz, &
    psplex, psibndp, psiaxisp, &
    ucoils, voltage, &
    psi_external_calc

implicit none

integer, intent(in) :: nucoils, ifplasma
double precision, intent(in) , dimension(nucoils) :: ucoils
type(type_parameters), intent(in) :: parameters_equil
type(type_equilibrium), intent(in) :: equil_in
type(type_equilibrium), intent(out) :: equil_out

integer :: nrplasma, j_init, j_call, j_vacplas

data j_call/0/
data j_vacplas/0/
data j_init/0/
save j_call, j_init, j_vacplas

write(*, *) ifplasma
ncoils = nucoils
call feqis_init(equil_in, parameters_equil, j_init, ifplasma)
voltage(1:ncoils) = ucoils(1:ncoils) ! voltage inputs for active conductors

write(*, *) ifplasma
nrplasma = nrho

! Init coils and grid
if (j_call == 0) then
    if (parameters_equil%k_fixfree == 1) then
        call equil_feqis_init_circ
    endif
endif

if (ifplasma == 1) then
    allocate(equil_out%eqgeometry%boundary%r(nteta))
    allocate(equil_out%eqgeometry%boundary%z(nteta))
    allocate(equil_out%eqgeometry%rectgrid%r2d(nr2))
    allocate(equil_out%eqgeometry%rectgrid%z2d(nz2))
    allocate(equil_out%eqgeometry%rectgrid%psirz2d(nr2, nz2))

    allocate(equil_out%coord_sys%position%r(nrplasma, nteta))
    allocate(equil_out%coord_sys%position%z(nrplasma, nteta))    
    allocate(equil_out%coord_sys%position%teta2d(nteta))    
    allocate(equil_out%coord_sys%position%rmin(nrplasma, nteta))    
    allocate(equil_out%coord_sys%position%psirz(nrplasma, nteta))    
    allocate(equil_out%coord_sys%gradvcell(nrplasma, nteta))
    allocate(equil_out%coord_sys%bpcell(nrplasma, nteta))
    allocate(equil_out%coord_sys%bcell(nrplasma, nteta))
    allocate(equil_out%coord_sys%rcell(nrplasma, nteta))
    allocate(equil_out%coord_sys%darea(nrplasma, nteta))
    allocate(equil_out%coord_sys%jphi(nrplasma, nteta))
    
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
    allocate(equil_out%profiles_1d%rho_tor(nrplasma) )
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

write(*, *) 'fix and nstep', parameters_equil%k_fixfree, parameters_equil%nstep, nrho, nteta

write(*, *) parameters_equil%k_fixfree, j_init

if (parameters_equil%k_fixfree == 1) then
    if (ifplasma == 0) then       ! only circuit equations solved
        write(*, *) 'vacuum'
        psi_cur_old = 0.
        psiplasmatoconduc = 0.
        call circuit_eq_advance_feqis(j_call)
        call psi_external_calc
        psirz = psiextrz
        j_vacplas = 0
    else if (ifplasma == 1) then  ! full plasma solved
        if (j_vacplas == 0) j_call = 0
        call full_system_advance_feqis(j_call)
        if (j_call == -1) then
            call convert_boundary_to_pbe
            call fix_boundary_feqis(1)
        endif
        j_vacplas = 1
    endif
else if (parameters_equil%k_fixfree == 0) then
    write(*, *) 'call fix equil code'
    call fix_boundary_feqis(j_init)
    equil_out%global_param%psplex   = psplex
    equil_out%global_param%psibound = psibndp
    equil_out%global_param%psiaxis  = psiaxisp
    call equil_assignments(equil_out)
    write(*, *) 'end fix equil code'
    j_init = 1
    write(*, *) 'end equil code'
    return
endif

write(*, *) 'end equil code'
j_call = 1
j_init = 1

if (ifplasma == 1) then
    call equil_assignments(equil_out)
endif

return
end subroutine feqis_main
