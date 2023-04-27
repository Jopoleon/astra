subroutine feqis_main(equil_in, equil_out)

use imas_ids, only: type_equilibrium
use feqis_geom, only: nrho, nteta, &
    teta, rgeom0, raxp, zaxp

implicit none

double precision, parameter :: GPI=3.141592653589793, &
    GPI2=2.*GPI

type(type_equilibrium), intent(in)  :: equil_in
type(type_equilibrium), intent(out) :: equil_out

integer :: jthe, j_init
double precision :: psi0, psiB
double precision, allocatable, dimension(:) :: psi_norm_in, psi_norm_out
double precision, allocatable, dimension(:, :) :: jrhoteta

data j_init/0/

save j_init

if (.not. allocated(teta)) allocate(teta(nteta))
if (.not. allocated(jrhoteta)) then
    allocate(psi_norm_in(nrho))
    allocate(psi_norm_out(nrho))
    allocate(jrhoteta(nrho, nteta))
endif

if (j_init == 0) then
    nteta = equil_in%eqgeometry%boundary%npoints
    nrho  = SIZE(equil_in%profiles_1d%pressure)
    Rgeom0 = equil_in%global_param%toroid_field%r0
! Initial polar grid: regular from 0 to 2*pi
    do jthe=1, nteta+1
        teta(jthe) = GPI2*(jthe - 1.)/(nteta + 0.)
    enddo
    raxp = Rgeom0
    zaxp = 0.
endif

psi0 = equil_in%profiles_1d%psi(1)
psiB = equil_in%profiles_1d%psi(nrho)
psi_norm_in = (equil_in%profiles_1d%psi - psi0)/(psiB - psi0) ! normalized: 0 axis,  1 sep

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

write(*, '(A, 2f8.4)') 'Call fix equil code', psi0, psiB

!boundary from previous time step
call PHI_EQ_2d_PBE( &
    nrho, nteta, &
    psi_norm_in, &
    equil_in%global_param%i_plasma/1.e6, &
    equil_in%profiles_1d%pressure, &
    equil_in%profiles_1d%ffprime, &
    equil_in%profiles_1d%pprime, &
    equil_in%global_param%toroid_field%b0, &
    rgeom0, &
    equil_in%eqgeometry%boundary%r, &
    equil_in%eqgeometry%boundary%z, &
    Raxp, Zaxp, &
    psi0, psiB, &
    equil_in%profiles_1d%F_dia, & 
    equil_out%coord_sys%position%r, &
    equil_out%coord_sys%position%z, &
    equil_out%coord_sys%position%psirz, &
    psi_norm_out, & 
    equil_out%profiles_1d%g2, &
    equil_out%profiles_1d%gm1, &
    equil_out%profiles_1d%r_outboard, &
    equil_out%profiles_1d%r_inboard, &
    equil_out%profiles_1d%volume, &
    equil_out%profiles_1d%g1, &
    equil_out%profiles_1d%gm41, &
    equil_out%profiles_1d%ggradro, &
    equil_out%profiles_1d%bmaxt, &
    equil_out%profiles_1d%bmint, &
    equil_out%profiles_1d%gm4, &
    equil_out%profiles_1d%bdb0, &
    equil_out%profiles_1d%gm5, &
    equil_out%profiles_1d%fofb, &
    equil_out%profiles_1d%areat, &
    equil_out%profiles_1d%perim, &
    equil_out%profiles_1d%shif, &
    equil_out%profiles_1d%elongation, &
    equil_out%profiles_1d%surface, & ! lateral surface
    equil_out%profiles_1d%tria_upper, &
    equil_out%coord_sys%position%teta2d, &
    equil_out%coord_sys%position%rmin, &
    jrhoteta, &
    equil_out%global_param%li3, &
    equil_out%global_param%betpol, &
    psi0)

write(*, '(A, 2f8.4)') 'Done fix equil code', psi0, psiB

raxp = equil_out%coord_sys%position%r(1, 1)
zaxp = equil_out%coord_sys%position%z(1, 1)

! additional info from rectangular grid
 
equil_out%global_param%i_plasma  = equil_in%global_param%i_plasma
equil_out%profiles_1d%tria_lower = equil_out%profiles_1d%tria_upper
equil_out%profiles_1d%ffprime  = 0.
equil_out%profiles_1d%pprime   = 0.
equil_out%profiles_1d%pressure = 0.
equil_out%profiles_1d%rho_tor  = 0.
equil_out%profiles_1d%F_dia    = 0.
equil_out%coord_sys%position%psirz = equil_out%coord_sys%position%psirz/GPI2

teta(1: nteta) = equil_out%coord_sys%position%teta2d(1: nteta)
teta(nteta+1) = teta(1) + GPI2

equil_out%global_param%psplex   = 0.
equil_out%global_param%psibound = psiB
equil_out%global_param%psiaxis  = psi0
! Unnormalise psi
equil_out%profiles_1d%psi = psi0 + (psiB - psi0)*psi_norm_out
j_init = 1

return
end subroutine feqis_main
