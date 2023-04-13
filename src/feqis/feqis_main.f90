subroutine feqis_main(nucoils, ucoils, parameters_spider, ifplasma, equil_in, equil_out)

use imas_ids, only: type_equilibrium
use parameters_a2spider, only: type_parameters
use ef_circuit, only: nrho, nteta, nr2, nz2, &
    rho, teta, iplasma, ipol, jrhoteta, btor0, rgeom0, li3, betapol, &
    raxp, zaxp, rbndp, zbndp, rpol, zpol, &
    psiaxisp, psibndp, psigrida, psirhoteta, &
    pressure, pprime, ffprime
use pi_grec_vars, only:  GPI2
use debugger_ef, only: markloc_ef, debug
use exchange_with_astra, only: nonegcurr, &
    raxis_astra, zaxis_astra, psi0_astra, psib_astra

implicit none

integer, intent(in) :: nucoils, ifplasma
double precision, intent(in), dimension(nucoils) :: ucoils
type(type_parameters), intent(in) :: parameters_spider
type(type_equilibrium), intent(in)  :: equil_in
type(type_equilibrium), intent(out) :: equil_out

integer :: j_init, i, j

data j_init/0/
save j_init

debug = 0 
call markloc_ef('feqis_main', debug_lev=debug)

call definitions_feqis(equil_in, parameters_spider, j_init, ifplasma)

if (ifplasma == 1) then
    allocate(equil_out%profiles_1d%psi(nrho))
    allocate(equil_out%profiles_1d%pressure(nrho))
    allocate(equil_out%profiles_1d%phi(nrho))
    allocate(equil_out%profiles_1d%pprime(nrho))
    allocate(equil_out%profiles_1d%ffprime(nrho))
    allocate(equil_out%profiles_1d%F_dia(nrho))
    allocate(equil_out%profiles_1d%q(nrho))
    allocate(equil_out%coord_sys%position%r(nrho, nteta))
    allocate(equil_out%coord_sys%position%z(nrho, nteta))    
    allocate(equil_out%coord_sys%position%teta2d(nteta))    
    allocate(equil_out%coord_sys%position%rmin(nrho, nteta))    
    allocate(equil_out%coord_sys%position%psirz(nrho, nteta))    
    allocate(equil_out%eqgeometry%boundary%r(nteta))
    allocate(equil_out%eqgeometry%boundary%z(nteta))
    allocate(equil_out%profiles_1d%rho_tor(nrho) )
    allocate(equil_out%profiles_1d%jparallel(nrho) )
    allocate(equil_out%profiles_1d%sigmapar%value(nrho) )
    allocate(equil_out%profiles_1d%jni%value(nrho) )
    allocate(equil_out%profiles_1d%te%value(nrho) )
    allocate(equil_out%coord_sys%gradvcell(nrho, nteta))
    allocate(equil_out%coord_sys%bpcell(nrho, nteta))
    allocate(equil_out%coord_sys%bcell(nrho, nteta))
    allocate(equil_out%coord_sys%rcell(nrho, nteta))
    allocate(equil_out%profiles_1d%gm1(nrho))
    allocate(equil_out%profiles_1d%gm4(nrho))
    allocate(equil_out%profiles_1d%gm5(nrho))
    allocate(equil_out%profiles_1d%gm41(nrho))
    allocate(equil_out%profiles_1d%rbp_b2(nrho))
    allocate(equil_out%profiles_1d%bplfs(nrho))
    allocate(equil_out%profiles_1d%acosB2a(nrho, 5))
    allocate(equil_out%profiles_1d%asinB2a(nrho, 5))
    allocate(equil_out%profiles_1d%acosBlnBa(nrho, 5))
    allocate(equil_out%profiles_1d%asinBlnBa(nrho, 5))     
    allocate(equil_out%profiles_1d%g1(nrho))
    allocate(equil_out%profiles_1d%g2(nrho))
    allocate(equil_out%profiles_1d%g2int(nrho))
    allocate(equil_out%profiles_1d%fofb(nrho))
    allocate(equil_out%profiles_1d%areat(nrho))
    allocate(equil_out%profiles_1d%perim(nrho))
    allocate(equil_out%profiles_1d%ggradro(nrho))
    allocate(equil_out%profiles_1d%gdroda(nrho))
    allocate(equil_out%profiles_1d%bmaxt(nrho))
    allocate(equil_out%profiles_1d%bmint(nrho))
    allocate(equil_out%profiles_1d%bdb0(nrho))
    allocate(equil_out%profiles_1d%dPSIdV(nrho))
    allocate(equil_out%profiles_1d%surface(nrho))
    allocate(equil_out%profiles_1d%volume(nrho))
    allocate(equil_out%profiles_1d%r_inboard(nrho))
    allocate(equil_out%profiles_1d%r_outboard(nrho))
    allocate(equil_out%profiles_1d%elongation(nrho))
    allocate(equil_out%profiles_1d%tria_upper(nrho))
    allocate(equil_out%profiles_1d%tria_lower(nrho)) 
    allocate(equil_out%profiles_1d%shif(nrho))
    allocate(equil_out%profiles_1d%shiv(nrho))
    allocate(equil_out%eqgeometry%rectgrid%r2d(nr2))
    allocate(equil_out%eqgeometry%rectgrid%z2d(nz2))
    allocate(equil_out%eqgeometry%rectgrid%psirz2d(nr2, nz2))
endif

write(*, *) 'call fix equil code'

! initial guess
if (j_init == 0) then
    raxp = raxis_astra
    zaxp = zaxis_astra
    psiaxisp = psi0_astra
    psibndp = psib_astra 
!boundary from experiment
endif 

!boundary from previous time step
call PHI_EQ_2d_PBE( &
    nrho, nteta, psigrida(1: nrho), iplasma, pressure(1: nrho), &
    ffprime(1: nrho), pprime(1: nrho), btor0, rgeom0, &
    rbndp(1: nteta), zbndp(1: nteta), Raxp, Zaxp, &
    psibndp, ipol(1: nrho), & 
    equil_out%coord_sys%position%r(1: nrho, 1: nteta), & 
    equil_out%coord_sys%position%z(1: nrho, 1: nteta), & 
    psirhoteta(1: nrho, 1: nteta), & 
! git    equil_out%profiles_1d%dPSIdV(1: nrho), & 
    equil_out%profiles_1d%psi(1: nrho), & 
    equil_out%profiles_1d%g2(1: nrho), &
    equil_out%profiles_1d%gm1(1: nrho), &
    equil_out%profiles_1d%r_outboard(1: nrho), &
    equil_out%profiles_1d%r_inboard(1: nrho), &
    equil_out%profiles_1d%volume(1: nrho), &
    equil_out%profiles_1d%g1(1: nrho), &
    equil_out%profiles_1d%gm41(1: nrho), &
    equil_out%profiles_1d%ggradro(1: nrho), &
    equil_out%profiles_1d%bmaxt(1: nrho), &
    equil_out%profiles_1d%bmint(1: nrho), &
    equil_out%profiles_1d%gm4(1: nrho), &
    equil_out%profiles_1d%bdb0(1: nrho), &
    equil_out%profiles_1d%gm5(1: nrho), &
    equil_out%profiles_1d%fofb(1: nrho), &
    equil_out%profiles_1d%areat(1: nrho), &
    equil_out%profiles_1d%perim(1: nrho), &
    equil_out%profiles_1d%shif(1: nrho), &
    equil_out%profiles_1d%elongation(1: nrho), &
    equil_out%profiles_1d%surface(1: nrho), & ! lateral surface
    equil_out%profiles_1d%tria_upper(1: nrho), &
! git    equil_out%profiles_1d%q(1: nrho), &
    equil_out%coord_sys%position%teta2d(1: nteta), &
    equil_out%coord_sys%position%rmin(1: nrho, 1: nteta), &
    jrhoteta(1: nrho, 1: nteta), &
! git    nonegcurr, &
    li3, betapol)

jrhoteta(1: nrho, nteta+1) = jrhoteta(1: nrho, 1) !periodic j

! additional info from rectangular grid
 
equil_out%global_param%li3 = li3
equil_out%global_param%betpol = betapol 
equil_out%global_param%i_plasma = iplasma*1.e6
equil_out%profiles_1d%tria_lower = equil_out%profiles_1d%tria_upper
equil_out%profiles_1d%ffprime(1: nrho) = 0.
equil_out%profiles_1d%pprime(1: nrho) = 0.
equil_out%profiles_1d%pressure(1: nrho) = 0.
equil_out%profiles_1d%rho_tor(1: nrho) = 0.
equil_out%profiles_1d%F_dia(1: nrho) = 0.

equil_out%coord_sys%position%psirz(1: nrho, 1: nteta) = psirhoteta(1: nrho, 1: nteta)/GPI2
rpol(1: nrho, 1: nteta) = equil_out%coord_sys%position%r(1: nrho, 1: nteta)
zpol(1: nrho, 1: nteta) = equil_out%coord_sys%position%z(1: nrho, 1: nteta)

do j=1, nteta
    do i=1, nrho
        rho(i, j) = sqrt((rpol(i, j) - raxp)**2 + (zpol(i, j) - zaxp)**2)
    enddo
enddo

teta(1: nteta) = equil_out%coord_sys%position%teta2d(1: nteta)
rho(1: nrho, nteta+1) = rho(1: nrho, 1)
teta(nteta+1) = teta(1) + GPI2

equil_out%global_param%psplex = 0.
equil_out%global_param%psibound = psibndp
equil_out%global_param%psiaxis = psiaxisp
write(*, *) 'end fix equil code'
j_init = 1

return
end subroutine feqis_main
