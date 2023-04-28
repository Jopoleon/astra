subroutine feqis_main(equil_in, equil_out)

use imas_ids, only: type_equilibrium
use feqis_geom, only: theta, raxp, zaxp
use interp_mod, only: polyfitcc
use pi_vars, only: GPI, GPI2

implicit none

type(type_equilibrium), intent(in)  :: equil_in
type(type_equilibrium), intent(out) :: equil_out

integer :: jrho, jthe, j_init, nrho, ntheta, i, i1, i2
double precision :: psi0, psiB, cnorm, X0, Y0, &
    yrr, ya, yrmax, yrmin, yzmax, yzmin, yrzmax, yrzmin, R0
double precision, dimension(3) :: xxxx1, yyyy1, pppp1
double precision, allocatable, dimension(:) :: psi_norm_in, psi_norm_out, &
    elon, tria, shif, r_in, r_out
double precision, allocatable, dimension(:, :) :: jrhotheta, &
    lambda2d, lambda2dp, XX, YY

data j_init/0/

save j_init

ntheta = equil_in%eqgeometry%boundary%npoints
nrho  = SIZE(equil_in%profiles_1d%pressure)

if (.not. allocated(theta)) then
    allocate(theta(ntheta+1))
endif
if (.not. allocated(lambda2d)) then
    allocate(psi_norm_in(nrho), psi_norm_out(nrho), elon(nrho), &
        tria(nrho), shif(nrho), r_in(nrho), r_out(nrho))
    allocate(jrhotheta(nrho, ntheta), XX(nrho, ntheta), YY(nrho, ntheta), &
        lambda2d(nrho, ntheta), lambda2dp(nrho, ntheta))
endif

if (j_init == 0) then
! Initial polar grid: regular from 0 to 2*pi
    do jthe=1, ntheta+1
        theta(jthe) = GPI2*(jthe - 1.)/(ntheta + 0.)
    enddo
    raxp = equil_in%global_param%toroid_field%r0
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
allocate(equil_out%coord_sys%position%r(nrho, ntheta))
allocate(equil_out%coord_sys%position%z(nrho, ntheta))    
allocate(equil_out%coord_sys%position%teta2d(ntheta))    
allocate(equil_out%coord_sys%position%rmin(nrho, ntheta))    
allocate(equil_out%coord_sys%position%psirz(nrho, ntheta))    
allocate(equil_out%eqgeometry%boundary%r(ntheta))
allocate(equil_out%eqgeometry%boundary%z(ntheta))
allocate(equil_out%profiles_1d%rho_tor(nrho) )
allocate(equil_out%profiles_1d%jparallel(nrho) )
allocate(equil_out%profiles_1d%sigmapar%value(nrho) )
allocate(equil_out%profiles_1d%jni%value(nrho) )
allocate(equil_out%profiles_1d%te%value(nrho) )
allocate(equil_out%coord_sys%gradvcell(nrho, ntheta))
allocate(equil_out%coord_sys%bpcell(nrho, ntheta))
allocate(equil_out%coord_sys%bcell(nrho, ntheta))
allocate(equil_out%coord_sys%rcell(nrho, ntheta))
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

write(*, '(A, 2f8.4)') 'Call FIX equil code', psi0, psiB

!boundary from previous time step
call PHI_EQ_2d_PBE( &
    nrho, ntheta, &
    psi_norm_in, &
    equil_in%global_param%i_plasma/1.e6, &
    equil_in%profiles_1d%ffprime, &
    equil_in%profiles_1d%pprime, &
    equil_in%global_param%toroid_field%b0, &
    equil_in%global_param%toroid_field%r0, &
    equil_in%eqgeometry%boundary%r, &
    equil_in%eqgeometry%boundary%z, &
    Raxp, Zaxp, &
    psi0, psiB, &
! Output
    equil_out%coord_sys%position%r, &
    equil_out%coord_sys%position%z, &
    equil_out%coord_sys%position%psirz, &
    psi_norm_out, lambda2d, lambda2dp, & 
    equil_out%coord_sys%position%teta2d, &
    psi0, cnorm, X0, Y0)

call build_2dgrid(nrho, ntheta,  &
    equil_in%eqgeometry%boundary%r, &
    equil_in%eqgeometry%boundary%z, &
    X0, Y0, &
    lambda2d, lambda2dp, psi_norm_out, &
    equil_out%coord_sys%position%psirz, &
    R0, &
    equil_in%profiles_1d%pressure, &
    equil_in%global_param%toroid_field%b0, &
    equil_in%profiles_1d%F_dia, &
    equil_in%global_param%i_plasma/1.e6, &
! Output
    XX, YY, &
    equil_out%coord_sys%position%rmin, &
    equil_out%profiles_1d%g2, &
    equil_out%profiles_1d%gm1, &
    equil_out%profiles_1d%areat, &
    equil_out%profiles_1d%perim, &
    equil_out%profiles_1d%volume, &
    equil_out%profiles_1d%g1, &
    equil_out%profiles_1d%ggradro, &
    equil_out%profiles_1d%bmaxt, &
    equil_out%profiles_1d%bmint, &
    equil_out%profiles_1d%gm4, &
    equil_out%profiles_1d%bdb0, &
    equil_out%profiles_1d%gm5, &
    equil_out%profiles_1d%fofb, &
    equil_out%profiles_1d%surface, &
    equil_out%global_param%li3, &
    equil_out%global_param%betpol)

do jthe=1, ntheta
    jrhotheta(:, jthe) = -GPI2*(equil_in%profiles_1d%ffprime(:) * 1./XX(:, jthe)/0.4/GPI + &
            XX(:, jthe)*1.e-6*equil_in%profiles_1d%pprime(:))/cnorm
enddo

do jrho=2, nrho

    i = minloc(yy(jrho, :), 1)
    i1 = i - 1
    i2 = i + 1
    if (i == 1 ) i1 = ntheta
    if (i == ntheta) i2 = 1
    xxxx1(1) = xx(jrho, i1)
    xxxx1(2) = xx(jrho, i)
    xxxx1(3) = xx(jrho, i2)
    yyyy1(1) = yy(jrho, i1)
    yyyy1(2) = yy(jrho, i)
    yyyy1(3) = yy(jrho, i2)
    call polyfitcc(xxxx1, yyyy1, pppp1)
    yrzmin = -pppp1(2)/(2.*pppp1(1))
    yzmin = pppp1(1)*yrzmin**2 + pppp1(2)*yrzmin + pppp1(3)

    i = maxloc(yy(jrho, :), 1)
    i1 = i - 1
    i2 = i + 1
    if (i ==  1) i1 = ntheta
    if (i == ntheta) i2 = 1
    xxxx1(1) = xx(jrho, i1)
    xxxx1(2) = xx(jrho, i)
    xxxx1(3) = xx(jrho, i2)
    yyyy1(1) = yy(jrho, i1)
    yyyy1(2) = yy(jrho, i)
    yyyy1(3) = yy(jrho, i2)
    call polyfitcc(xxxx1, yyyy1, pppp1)
    yrzmax = -pppp1(2)/(2.*pppp1(1))
    yzmax = pppp1(1)*yrzmax**2 + pppp1(2)*yrzmax + pppp1(3)

    yrmin = MINVAL(xx(jrho, :))
    yrmax = MAXVAL(xx(jrho, :))

    yrr = 0.5*(yrmax + yrmin)
    ya  = 0.5*(yrmax - yrmin)

    SHIF (jrho) = yrr - R0
    ELON (jrho) = (yzmax - yzmin)/(yrmax - yrmin)
    TRIA (jrho) = (yrr - 0.5*(yrzmin + yrzmax))/ya
    r_out(jrho) = yrmax
    r_in (jrho) = yrmin
enddo

r_out(1) = xx(1, 1)
r_in(1)  = xx(1, 1)
ELON(1)  = ELON(2)
TRIA(1)  = 0.d0
SHIF(1)  = XX(1, 1) - R0

equil_out%coord_sys%position%r = XX
equil_out%coord_sys%position%z = YY
equil_out%profiles_1d%r_outboard = r_out
equil_out%profiles_1d%r_inboard  = r_in
equil_out%profiles_1d%shif = SHIF
equil_out%profiles_1d%elongation = ELON
equil_out%profiles_1d%tria_upper = TRIA

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

theta(1: ntheta) = equil_out%coord_sys%position%teta2d(1: ntheta)
theta(ntheta+1) = theta(1) + GPI2

equil_out%global_param%psplex   = 0.
equil_out%global_param%psibound = psiB
equil_out%global_param%psiaxis  = psi0
! Unnormalise psi
equil_out%profiles_1d%psi = psi0 + (psiB - psi0)*psi_norm_out
j_init = 1

return
end subroutine feqis_main
