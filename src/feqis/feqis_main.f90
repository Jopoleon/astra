subroutine feqis_main(equil_in, equil_out)

use imas_ids, only: type_equilibrium
use ef_circuit, only: nrho, nteta, nr2, nz2, &
    teta, iplasma, ipol, &
    btor0, rgeom0, li3, betapol, &
    raxp, zaxp, rbndp, zbndp, psiaxisp, psibndp, psigrida, &
    pressure, pprime, ffprime

implicit none

double precision, parameter :: GPI=3.141592653589793, &
    GPI2=2.*GPI, mu0=0.4*GPI

type(type_equilibrium), intent(in)  :: equil_in
type(type_equilibrium), intent(out) :: equil_out

integer :: jthe
integer :: j_init
double precision, allocatable, dimension(:, :) :: jrhoteta

data j_init/0/

save j_init

if (j_init == 0) then
    nteta = equil_in%eqgeometry%boundary%npoints
    nrho  = SIZE(equil_in%profiles_1d%pressure)
    Rgeom0 = equil_in%global_param%toroid_field%r0
! polar grid, goes from 0 to 2*pi
    do jthe=1, nteta+1
        teta(jthe) = GPI2*(jthe - 1.)/(nteta + 0.)
    enddo
endif

btor0   = equil_in%global_param%toroid_field%b0
iplasma = equil_in%global_param%i_plasma/1.e6
pressure(1: nrho) = equil_in%profiles_1d%pressure(1: nrho)
pprime(  1: nrho) = equil_in%profiles_1d%pprime(1: nrho)
ffprime( 1: nrho) = equil_in%profiles_1d%ffprime(1: nrho)
psigrida(1: nrho) = equil_in%profiles_1d%psi(1: nrho) !unnormalized
psigrida(1: nrho) = (psigrida(1: nrho) - psigrida(1)) / &
                    (psigrida(nrho)    - psigrida(1)) ! normalized: 0 axis,  1 sep
ipol(1: nrho)   = equil_in%profiles_1d%F_dia(1: nrho)
rbndp(1: nteta) = equil_in%eqgeometry%boundary%r(1: nteta)
zbndp(1: nteta) = equil_in%eqgeometry%boundary%z(1: nteta)
rbndp(nteta+1)  = rbndp(1)
zbndp(nteta+1)  = zbndp(1)

if (.not. allocated(jrhoteta)) allocate(jrhoteta(nrho, nteta))

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

write(*, *) 'call fix equil code'

! initial guess
if (j_init == 0) then
    raxp = 0.
    zaxp = 0.
    psiaxisp = equil_in%profiles_1d%psi(1)
    psibndp  = equil_in%profiles_1d%psi(nrho)
!boundary from experiment
endif 

!boundary from previous time step
call PHI_EQ_2d_PBE( &
    nrho, nteta, &
    psigrida(1: nrho), &
    iplasma, &
    pressure(1: nrho), &
    ffprime(1: nrho), &
    pprime(1: nrho), &
    btor0, &
    rgeom0, &
    rbndp(1: nteta), &
    zbndp(1: nteta), &
    Raxp, &
    Zaxp, &
    psibndp, &
    ipol(1: nrho), & 
    equil_out%coord_sys%position%r, &
    equil_out%coord_sys%position%z, &
    equil_out%coord_sys%position%psirz, &
    equil_out%profiles_1d%psi, & 
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
    li3, betapol)

raxp = equil_out%coord_sys%position%r(1, 1)
zaxp = equil_out%coord_sys%position%z(1, 1)

! additional info from rectangular grid
 
equil_out%global_param%li3 = li3
equil_out%global_param%betpol = betapol 
equil_out%global_param%i_plasma = iplasma*1.e6
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
equil_out%global_param%psibound = psibndp
equil_out%global_param%psiaxis  = psiaxisp
j_init = 1

return
end subroutine feqis_main
