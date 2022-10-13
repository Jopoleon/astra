subroutine feqis_main(nucoils, ucoils, parameters_spider, &
    ifplasma, equil_in, equil_out)

use schemas, only: type_equilibrium
use parameters, only: type_parameters
use ef_circuit, only: nrho, nteta, ncoils, nr2, nz2, &
    voltage, psi_cur_old, psiplasmatoconduc, psirz, psiextrz, &
    psplex, psibnd, psiaxis, psibndp, psiaxisp, raxp, rax, iaxis, jaxis
use debugger_ef, only: markloc_ef, debug

implicit none

integer, intent(in) :: nucoils, ifplasma
double precision, intent(in), dimension(nucoils) :: ucoils
type(type_parameters), intent(in) :: parameters_spider
type(type_equilibrium), intent(in)  :: equil_in
type(type_equilibrium), intent(out) :: equil_out

integer :: j_init, j_call, j_vacplas

data j_call/0/
data j_vacplas/0/
data j_init/0/
save j_call, j_init, j_vacplas

debug = 0 
call markloc_ef('feqis_main', debug_lev=debug)

call definitions_feqis(equil_in, parameters_spider, j_init, ifplasma)
ncoils = nucoils
voltage(1: ncoils) = ucoils(1: ncoils) ! voltage inputs for active conductors

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

write(*, *) 'fix and nstep', parameters_spider%k_fixfree, & 
    parameters_spider%nstep, ifplasma, nrho, nteta

! init coils and grid
if (j_call == 0) then
    if (parameters_spider%k_fixfree == 1) then
        call feqis_init_circ
    endif
endif

write(*, *) parameters_spider%k_fixfree, j_init
! call fix boundary code, also here boundary comes from experiment

if (parameters_spider%k_fixfree == 1) then
    if (ifplasma == 0) then ! only circuit equations
        write(*, *) 'vacuum'
        psi_cur_old = 0.
        psiplasmatoconduc = 0.
        call circuit_eq_advance_ef(j_call) 
        call psi_external_calc_ef
        psirz = psiextrz
        j_vacplas = 0
    endif

    if (ifplasma == 1) then ! full plasma solved
        write(*, *) 'full plasma'
        if (j_vacplas == 0) j_call = 0
        call full_system_advance_ef(j_call)
        if (j_call == -1) then
            call find_new_boundary_part2
            call convert_boundary_to_pbe
            call fix_boundary_ef(1, equil_out)
            call psiplex_calc_ef
            equil_out%global_param%psplex = psplex
            equil_out%global_param%psibound = psibnd
            equil_out%global_param%psiaxis = psiaxis
        endif
        j_vacplas = 1  
    endif

endif !kfixfree

if (parameters_spider%k_fixfree == 0) then
    write(*, *) 'call fix equil code'
    call fix_boundary_ef(j_init, equil_out)
    equil_out%global_param%psplex = 0.
    equil_out%global_param%psibound = psibndp
    equil_out%global_param%psiaxis = psiaxisp
    write(*, *) 'end fix equil code', rax, raxp, iaxis, jaxis
    j_init = 1
    return
endif

write(*, *) 'end equil code'
j_call = 1
j_init = 1

return
end subroutine feqis_main

!------------------------------------------------------------
subroutine convert_boundary_to_pbe

use pi_grec_vars, only: GPI2
use ef_circuit, only: nbnd, nteta, i_dim5, teta, &
    rax, zax, raxp, zaxp, rbnd, zbnd, rbndp, zbndp, &
    psiaxis, psiaxisp, psibnd, psibndp
use debugger_ef, only: markloc_ef, debug

implicit none

double precision :: teta_fbe(nbnd+1)
integer :: i

call markloc_ef('convert_boundary_to_pbe', debug_lev=debug)

do i=1, nbnd
    call find_angle_ef(rax, zax, rbnd(i), zbnd(i), teta_fbe(i)) 
enddo
do i=1, nbnd-1
    if (teta_fbe(i+1) < teta_fbe(i)) teta_fbe(i) = teta_fbe(i) - GPI2
enddo
teta_fbe(nbnd+1) = teta_fbe(1) + GPI2

if (teta_fbe(1) < 0.) then
    teta_fbe(1: nbnd) = teta_fbe(2: nbnd+1)
    rbnd(1: nbnd) = rbnd(2: nbnd+1)
    zbnd(1: nbnd) = zbnd(2: nbnd+1)
    teta_fbe(nbnd+1) = teta_fbe(1) + GPI2 
    rbnd(nbnd+1) = rbnd(1)
    zbnd(nbnd+1) = zbnd(1)
endif

call linterp_ef_feqis(teta_fbe(1: nbnd+1), rbnd(1: nbnd+1), nbnd+1, &
    teta(1: nteta+1), rbndp(1: nteta+1), nteta+1)
call linterp_ef_feqis(teta_fbe(1: nbnd+1), zbnd(1: nbnd+1), nbnd+1, &
    teta(1: nteta+1), zbndp(1: nteta+1), nteta+1)

raxp = rax
zaxp = zax
psibndp = psibnd
psiaxisp = psiaxis

return
end subroutine convert_boundary_to_pbe

!------------------------------------------------------------
subroutine fix_boundary_ef(j_init, equil_out)

use pi_grec_vars, only:  GPI2
use ef_circuit, only: nrho, nteta, nr2, nz2, &
    rho, teta, &
    iplasma, ipol, jrhoteta, btor0, rgeom0, li3, betapol, &
    raxp, zaxp, rbndp, zbndp, r, z, rpol, zpol, &
    psiaxisp, psibndp, psigrida, psirhoteta, psirz, &
    pressure, pprime, ffprime
use exchange_with_astra, only: nonegcurr, &
    raxis_astra, zaxis_astra, psi0_astra, psib_astra
use schemas, only: type_equilibrium
use debugger_ef, only: markloc_ef, debug

implicit none

integer, intent(in) :: j_init
type(type_equilibrium), intent(inout) :: equil_out

integer :: i, j

call markloc_ef('fix_boundary_ef', debug_lev=debug)

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
    equil_out%profiles_1d%dPSIdV(1: nrho), & 
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
    equil_out%profiles_1d%q(1: nrho), &
    equil_out%coord_sys%position%teta2d(1: nteta), &
    equil_out%coord_sys%position%rmin(1: nrho, 1: nteta), &
    jrhoteta(1: nrho, 1: nteta), &
    nonegcurr, li3, betapol)

jrhoteta(1: nrho, nteta+1) = jrhoteta(1: nrho, 1) !periodic j

! additional info from rectangular grid
 
equil_out%eqgeometry%rectgrid%npointsr = nr2
equil_out%eqgeometry%rectgrid%npointsz = nz2
equil_out%eqgeometry%rectgrid%r2d(1: nr2) = r(1: nr2)
equil_out%eqgeometry%rectgrid%z2d(1: nz2) = z(1: nz2) 
equil_out%eqgeometry%rectgrid%psirz2d(1: nr2, 1: nz2) = -psirz(1: nr2, 1: nz2)/GPI2
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

return
end subroutine fix_boundary_ef

!------------------------------------------------------------
subroutine solve_gs2d(g)

use pi_grec_vars, only:  GPI, mu0
use ef_circuit, only: i_dim2, nr1, nz1, nr, nz, nr2, nz2, &
    dr, dz, r, jrz, rcomp
use debugger_ef, only: markloc_ef, debug

implicit none

double precision, intent(out), dimension(i_dim2, i_dim2) :: g

integer :: i, j, k, k_fourier, j_init
double precision :: x1, x2
double precision, dimension(258) :: A, B, C, z_fourier
double precision, dimension(i_dim2, i_dim2) :: gt, rhs, wrhs

data j_init/0/
save A, B, C, j_init, z_fourier

! g(1, jz), g(nr2, jz), g(jr, 1), gr(jr, nz2) contain the b.c.
! solution is g(2: nr1, 2: nz1) --> nr*nz points
! solution goes back to g
call markloc_ef('solve_gs2d', debug_lev=debug)

do j=1, nz2
    do i=1, nr2
        rhs(i, j) = -mu0*r(i)*jrz(i, j)
    enddo
enddo

rhs(2: nr1, 2)   = rhs(2: nr1, 2) - g(2: nr1, 1)/dz**2
rhs(nr1, 2: nz1) = rhs(nr1, 2: nz1) - g(nr2, 2: nz1)/dr**2*2./(r(nr2) + r(nr1))*r(nr1)
rhs(2: nr1, nz1) = rhs(2: nr1, nz1) - g(2: nr1, nz2)/dz**2
rhs(2, 2: nz1)   = rhs(2, 2: nz1) - g(1, 2: nz1)/dr**2*2./(r(1) + r(2))*r(2)

do i=2, nr1
   call discrete_sine_transform_ef(nz, rhs(i, 2: nz1), wrhs(i, 2: nz1))
enddo

k_fourier = nz
!create inverse matrix for gs2d
if (j_init == 0) then
    A = 0.
    B = 0.
    C = 0.
    do i=2, nr-1  
        x1 = 0.5*(rcomp(i) + rcomp(i-1))  
        x2 = 0.5*(rcomp(i+1) + rcomp(i))  
        B(i) = -rcomp(i)/dr**2*(1./x2 + 1./x1) 
        A(i) =  rcomp(i)/x2/dr**2 
        C(i) =  rcomp(i)/x1/dr**2 
    enddo
    i = 1
    x1 = 0.5*(rcomp(i) + r(1))  
    x2 = 0.5*(rcomp(i+1) + rcomp(i))  
    B(i) = -rcomp(i)/dr**2*(1./x2 + 1./x1) 
    A(i) =  rcomp(i)/x2/dr**2 
    i = nr
    x1 = 0.5*(rcomp(i) + rcomp(i-1))  
    x2 = 0.5*(r(nr+2) + rcomp(i))  
    B(i) = -rcomp(i)/dr**2*(1./x2 + 1./x1) 
    C(i) =  rcomp(i)/x1/dr**2 
    j_init = 1
    do k=2, nz1
        z_fourier(k) = 2./dz**2*(cos((k - 1)*GPI/(nz1)) - 1.)
    enddo
endif

gt = 0.
!solve matrix
do k=2, nz1
    call solve_tridiag_fbe_ef(C(1: nr), B(1: nr) + z_fourier(k), A(1: nr), wrhs(2: nr1, k), gt(2: nr1, k), nr)
enddo

do i=2, nr1
    call discrete_sine_transform_ef(nz, gt(i, 2: nz1), gt(i, 2: nz1))
enddo
g(2: nr1, 2: nz1) = 2./(nz + 1)*gt(2: nr1, 2: nz1)

return
end subroutine solve_gs2d

!------------------------------------------------------------
subroutine solve_tridiag_fbe_ef(A, B, C, R, f, Ngrid)

! Provides solution of the system: 
!
!   Aj fj-1  + Bj fj  + Cj fj+1 = Rj
!   where j = 1..Ngrid
!   bcbound = 1  -> given f_NA1
!   eximp = 2:  implicit

use debugger_ef, only: markloc_ef, debug

implicit none

integer, intent(in) :: Ngrid
double precision, intent(in) , dimension(Ngrid) :: A, B, C, R 
double precision, intent(out), dimension(Ngrid) :: f

integer :: j, k
double precision, dimension(Ngrid) :: alpha, beta

call markloc_ef('solve_tridiag_fbe_ef', debug_lev=debug)

alpha = 0.
beta  = 0.

j = 1
alpha(j) = C(j)/B(j)
beta(j)  = R(j)/B(j)

do j=2, Ngrid-1
    alpha(j) = C(j)/(B(j) - A(j)*alpha(j-1))
enddo
do j=2, Ngrid
    beta(j) = (R(j) - A(j)*beta(j-1))/(B(j) - A(j)*alpha(j-1))
enddo

! Note that boundary value is assumed to be on the main last grid point, so 1-dx/2. This has to be
! corrected later on... CEfable
j = Ngrid

! This one should be appropriate with extrapolation... but now go back to real b.c.

f(j) = beta(j)
do k=1, Ngrid-1
    j = Ngrid - 1 - k + 1
    f(j) = beta(j) - alpha(j)*f(j+1)
enddo

return
end subroutine solve_tridiag_fbe_ef

!------------------------------------------------------------
subroutine find_new_axis_part1

use ef_circuit, only: iaxis, jaxis, rax, zax, trax, tzax, r, z, &
    psiaxis, psirz, derivpsi
use debugger_ef, only: markloc_ef, debug

implicit none

integer :: i1, i2, j1, j2, j, iax, jax
double precision :: errtol, tolerr
double precision, dimension(2) :: ppx
double precision, dimension(6) :: ccc
double precision, dimension(9) :: bub, xub, yub

!the used function is psirz

! 1) find new magnetic axis
!old axis

call markloc_ef('find_new_axis_part1', debug_lev=debug)

iax = iaxis
jax = jaxis
i1 = 10000
j1 = 10000
i2 = 10000
j2 = 10000
errtol = 1.e6
tolerr = 10.
j = 1

do while(errtol > tolerr)
    i1 = i2
    j1 = j2
    i2 = iax
    j2 = jax
    if (psirz(iax+1, jax) > psirz(iax, jax)) then
        iax = iax + 1
        jax = jax
    endif
    if (psirz(iax-1, jax) > psirz(iax, jax)) then
        iax = iax - 1
        jax = jax
    endif
    if (psirz(iax, jax+1) > psirz(iax, jax)) then
        iax = iax
        jax = jax + 1
    endif
    if (psirz(iax, jax-1) > psirz(iax, jax)) then
        iax = iax
        jax = jax - 1
    endif
    if (psirz(iax+1, jax-1) > psirz(iax, jax)) then
        iax = iax + 1
        jax = jax - 1
    endif
    if (psirz(iax-1, jax+1) > psirz(iax, jax)) then
        iax = iax - 1
        jax = jax + 1
    endif
    if (psirz(iax-1, jax-1) > psirz(iax, jax)) then
        iax = iax - 1
        jax = jax - 1
    endif
    if (psirz(iax+1, jax+1) > psirz(iax, jax)) then
        iax = iax + 1
        jax = jax + 1
    endif

    rax = r(iax)
    zax = z(jax)

    j = j + 1
    if ((iax == i1) .and. (jax == j1)) EXIT
    if (j >= 100000) EXIT
enddo

iaxis = iax
jaxis = jax

!find true axis
xub(1) = r(iax-1)
xub(2) = r(iax)
xub(3) = r(iax+1)
xub(4) = r(iax)
xub(5) = r(iax)
xub(6) = r(iax-1)
xub(7) = r(iax-1)
xub(8) = r(iax+1)
xub(9) = r(iax+1)
yub(1) = z(jax)
yub(2) = z(jax)
yub(3) = z(jax)
yub(4) = z(jax-1)
yub(5) = z(jax+1)
yub(6) = z(jax-1)
yub(7) = z(jax+1)
yub(8) = z(jax-1)
yub(9) = z(jax+1)
bub(1) = psirz(iax-1, jax  )
bub(2) = psirz(iax  , jax  )
bub(3) = psirz(iax+1, jax  )
bub(4) = psirz(iax  , jax-1)
bub(5) = psirz(iax  , jax+1)
bub(6) = psirz(iax-1, jax-1)
bub(7) = psirz(iax-1, jax+1)
bub(8) = psirz(iax+1, jax-1)
bub(9) = psirz(iax+1, jax+1)

call least_square_biquad_ef(xub, yub, bub, 9, ccc, ppx(1), ppx(2), psiaxis, derivpsi)
rax = ppx(1) ! r(iaxis)
zax = ppx(2) ! z(jaxis)
trax = rax
tzax = zax

return
end subroutine find_new_axis_part1 

!------------------------------------------------------------
subroutine find_psi_boundary

use pi_grec_vars, only:  GPI
use ef_circuit, only: nr1, nz1, nr2, nz2, nlimiter, n_of_xpoints, &
    i_plasmatype, max_xpoints, psiaxis, psibnd, psirz, &
    zlimpotential, limiterr, limiterz, u_n, &
    rax, zax, dr, dz, r, z, rinner, r_xpoint, z_xpoint, raus, ztop, zbot
use debugger_ef, only: markloc_ef, debug

implicit none

integer :: i, j, k, i4, i5, i_county
double precision :: x1, x2, x5
double precision :: pos_xpoint(2), ddipsi(8), psi_limp(500), &
    psi_xpoint(max_xpoints)
double precision, external :: frlim_ef

data i_county/0/
save i_county

call markloc_ef('find_psi_boundary', debug_lev=debug)

i_plasmatype = 0

if (i_county == 1) then

! go through old x-points and see where they end up
319 continue
    if (n_of_xpoints >= 1) then
        do i=1, n_of_xpoints
            call find_actual_index_ef(r_xpoint(i), z_xpoint(i), j, k)
            call nine_point_regression(r(j), z(k), pos_xpoint, ddipsi, x1)
            if ((pos_xpoint(1) >= r(1)) .and. (pos_xpoint(1) <= r(nr2)) .and. & 
                (pos_xpoint(2) >= z(1)) .and. (pos_xpoint(2) <= z(nz2))) then
                r_xpoint(i) = pos_xpoint(1)
                z_xpoint(i) = pos_xpoint(2)
            else
! xpoint doesnt exist anymore
                r_xpoint(i: n_of_xpoints-1) = r_xpoint(i+1: n_of_xpoints)
                z_xpoint(i: n_of_xpoints-1) = z_xpoint(i+1: n_of_xpoints)
                n_of_xpoints = n_of_xpoints - 1
                goto 319
            endif
        enddo
    endif

! Scan the boundary to find new x-points
    j = 2
    do i=2, nr1
        call nine_point_regression(r(i), z(j), pos_xpoint, ddipsi, x1)
            x5 = (ddipsi(5)**2 - ddipsi(3)*ddipsi(4))
            if ((abs(pos_xpoint(1) - r(i)) <= 2.*dr) .and. &
                (abs(pos_xpoint(2) - z(j)) <= 2.*dz) .and. (x5 >= 0.)) then
                n_of_xpoints = min(max_xpoints, n_of_xpoints+1)
                r_xpoint(n_of_xpoints) = pos_xpoint(1)
                z_xpoint(n_of_xpoints) = pos_xpoint(2)
                do k=1, n_of_xpoints ! check if double counted
                if ((k < n_of_xpoints) .and. &
                    (abs(pos_xpoint(1) - r_xpoint(k)) <= 2.*dr) .and. & 
                    (abs(pos_xpoint(2) - z_xpoint(k)) <= 2.*dz)) then
                    r_xpoint(k) = 0.5*(r_xpoint(n_of_xpoints) + r_xpoint(k)) 
                    z_xpoint(k) = 0.5*(z_xpoint(n_of_xpoints) + z_xpoint(k)) 
                    n_of_xpoints = n_of_xpoints - 1
                endif
            enddo
        endif
    enddo

    j = nz1
    do i=2, nr1
        call nine_point_regression(r(i), z(j), pos_xpoint, ddipsi, x1)
        x5 = (ddipsi(5)**2 - ddipsi(3)*ddipsi(4))
        if ((abs(pos_xpoint(1) - r(i)) <= 2.*dr) .and. &
            (abs(pos_xpoint(2) - z(j)) <= 2.*dz) .and. (x5 >= 0.)) then
            n_of_xpoints = min(max_xpoints, n_of_xpoints+1)
            r_xpoint(n_of_xpoints) = pos_xpoint(1)
            z_xpoint(n_of_xpoints) = pos_xpoint(2)
            do k=1, n_of_xpoints ! check if double counted
                if ((k < n_of_xpoints) .and. &
                    (abs(pos_xpoint(1) - r_xpoint(k)) <= 2.*dr) .and. & 
                    (abs(pos_xpoint(2) - z_xpoint(k)) <= 2.*dz)) then
                    r_xpoint(k) = 0.5*(r_xpoint(n_of_xpoints) + r_xpoint(k)) 
                    z_xpoint(k) = 0.5*(z_xpoint(n_of_xpoints) + z_xpoint(k)) 
                    n_of_xpoints = n_of_xpoints - 1
                endif
            enddo
        endif
    enddo

    i = 2
    do j=2, nz1
        call nine_point_regression(r(i), z(j), pos_xpoint, ddipsi, x1)
        x5 = (ddipsi(5)**2 - ddipsi(3)*ddipsi(4))
        if ((abs(pos_xpoint(1) - r(i)) <= 2.*dr) .and. &
            (abs(pos_xpoint(2) - z(j)) <= 2.*dz) .and. (x5 >= 0.)) then
            n_of_xpoints = min(max_xpoints, n_of_xpoints+1)
            r_xpoint(n_of_xpoints) = pos_xpoint(1)
            z_xpoint(n_of_xpoints) = pos_xpoint(2)
            do k=1, n_of_xpoints ! check if double counted
                if ((k < n_of_xpoints) .and. &
                    (abs(pos_xpoint(1) - r_xpoint(k)) <= 2.*dr) .and. & 
                    (abs(pos_xpoint(2) - z_xpoint(k)) <= 2.*dz)) then
                    r_xpoint(k) = 0.5*(r_xpoint(n_of_xpoints) + r_xpoint(k)) 
                    z_xpoint(k) = 0.5*(z_xpoint(n_of_xpoints) + z_xpoint(k)) 
                    n_of_xpoints = n_of_xpoints - 1
                endif
            enddo
        endif
    enddo

    i = nr1
    do j=2, nz1
        call nine_point_regression(r(i), z(j), pos_xpoint, ddipsi, x1)
        x5 = (ddipsi(5)**2 - ddipsi(3)*ddipsi(4))
        if ((abs(pos_xpoint(1) - r(i)) <= 2.*dr) .and. &
            (abs(pos_xpoint(2) - z(j)) <= 2.*dz) .and. (x5 >= 0.)) then
            n_of_xpoints = min(max_xpoints, n_of_xpoints+1)
            r_xpoint(n_of_xpoints) = pos_xpoint(1)
            z_xpoint(n_of_xpoints) = pos_xpoint(2)
            do k=1, n_of_xpoints ! check if double counted
                if ((k < n_of_xpoints) .and. &
                    (abs(pos_xpoint(1) - r_xpoint(k)) <= 2.*dr) .and. & 
                    (abs(pos_xpoint(2) - z_xpoint(k)) <= 2.*dz)) then
                    r_xpoint(k) = 0.5*(r_xpoint(n_of_xpoints) + r_xpoint(k)) 
                    z_xpoint(k) = 0.5*(z_xpoint(n_of_xpoints) + z_xpoint(k)) 
                    n_of_xpoints = n_of_xpoints - 1
                endif
            enddo
        endif
    enddo

endif

if (i_county == 0) then ! do a full pass to find all X-points
    n_of_xpoints = 0
    do j=2, nz1
        do i=2, nr1
            call nine_point_regression(r(i), z(j), pos_xpoint, ddipsi, x1)
            x5 = (ddipsi(5)**2 - ddipsi(3)*ddipsi(4))
            if ((abs(pos_xpoint(1) - r(i)) <= 2.*dr) .and. &
                (abs(pos_xpoint(2) - z(j)) <= 2.*dz) .and. (x5 >= 0.)) then
                n_of_xpoints = min(max_xpoints, n_of_xpoints+1)
                r_xpoint(n_of_xpoints) = pos_xpoint(1)
                z_xpoint(n_of_xpoints) = pos_xpoint(2)
                do k=1, n_of_xpoints ! check if double counted
                    if ((k < n_of_xpoints) .and. &
                        (abs(pos_xpoint(1) - r_xpoint(k)) <= 2.*dr) .and. & 
                        (abs(pos_xpoint(2) - z_xpoint(k)) <= 2.*dz)) then
                        r_xpoint(k) = 0.5*(r_xpoint(n_of_xpoints) + r_xpoint(k))
                        z_xpoint(k) = 0.5*(z_xpoint(n_of_xpoints) + z_xpoint(k))
                        n_of_xpoints = n_of_xpoints - 1
                    endif
                enddo
            endif
        enddo
    enddo

    i_county = 1
 
 !exclude magnetic axis
    do k=1, n_of_xpoints 
        if (abs(r_xpoint(k) - rax) <= 2.*dr .and. &
            abs(z_xpoint(k) - zax) <= 2.*dz) then
            r_xpoint(k: n_of_xpoints-1) = r_xpoint(k+1: n_of_xpoints)
            z_xpoint(k: n_of_xpoints-1) = z_xpoint(k+1: n_of_xpoints)
            n_of_xpoints = n_of_xpoints - 1
            EXIT
        endif
    enddo
endif

!calculate limiter flux
do i=1, nlimiter
    call find_fields_interp_ef_psionly(limiterR(i), limiterZ(i), &
        psi_limp(i)) !give back psi, br, bz at r0, z0
enddo

! now, remove limiters that are in the shadow of xpoints
if (n_of_xpoints >= 1) then
    i_plasmatype = 1
    do i=1, n_of_xpoints
        call find_actual_index_ef(r_xpoint(i), z_xpoint(i), j, k)
        call nine_point_regression(r(j), z(k), pos_xpoint, ddipsi, psi_xpoint(i))
        if (zlimpotential(j, k) > 0.) then
            do j=1, nlimiter
                x1 = frlim_ef(ddipsi(1: 5), limiterz(j), &
                    pos_xpoint(1), pos_xpoint(2), rax, zax)
                x2 = frlim_ef(ddipsi(1: 5), zax, &
                    pos_xpoint(1), pos_xpoint(2), rax, zax)
                if ((limiterr(j) - x1)*(rax - x2) <= 0.) psi_limp(j) = -1.e6
            enddo
        else
! xpoint is outside of the limiter, doesnt count
            psi_xpoint(i) = -1.e6
        endif
    enddo
!psi boundary is the highest of all the values
    x1 = maxval(psi_xpoint(1: n_of_xpoints), 1)
    x2 = maxval(psi_limp(1: nlimiter), 1)
    i4 = maxloc(psi_xpoint(1: n_of_xpoints), 1)
    i5 = maxloc(psi_limp(1: nlimiter), 1)
    psibnd = max(x1, x2)
    if (x2 > x1) then
        i_plasmatype = 0
    else
        i_plasmatype = 1
    endif
    write(*, *) psibnd, x1, x2, i_plasmatype, psi_xpoint(i4), &
        psi_limp(i5), limiterr(i5), limiterz(i5), r_xpoint(i4), z_xpoint(i4)
else
!no x-points, take largest limiter flux
    psibnd = maxval(psi_limp(1: nlimiter), 1)
    i_plasmatype = 0
    write(*, *) psibnd, i_plasmatype
endif

!normalized flux
u_n(1: nr2, 1: nz2) = (psirz(1: nr2, 1: nz2) - psiaxis)/(psibnd - psiaxis)

!find plasma limits  !improve this using xpoints angles 
ztop =  1.e6
zbot = -1.e6
raus =  1.e6
rinner = 0.
if (n_of_xpoints >= 1) then
    do i=1, n_of_xpoints
        if (psi_xpoint(i) > -1.e5) then
            x1 = ATAN2(z_xpoint(i) - zax, r_xpoint(i) - rax)
            if (r_xpoint(i) > rax .and. (x1 >= -1./4.*GPI .and. x1 <= GPI/4.)) raus = min(raus, r_xpoint(i))
            if (z_xpoint(i) > zax .and. (x1 >= GPI/4. .and. x1 <= 3./4.*GPI)) ztop = min(ztop, z_xpoint(i))
            if (r_xpoint(i) < rax .and. (abs(x1) >= 3./4.*GPI)) rinner = max(rinner, r_xpoint(i))
            if (z_xpoint(i) < zax .and. (x1 >= -3./4.*GPI .and. x1 <= -1./4.*GPI)) zbot = max(zbot, z_xpoint(i))
        endif
    enddo
endif

write(*, *) raus, rinner, ztop, zbot

open(32, file='fort.44442')
    write(32, *) psirz(1: nr2, 1: nz2)
close(32)

return
end subroutine find_psi_boundary

!------------------------------------------------------------
subroutine find_new_boundary_part2

use ef_circuit, only: nr2, nz2, n_of_xpoints, nbnd, i_dim2, i_dim5, &
    iaxis, jaxis, ibnd, jbnd, i_plasmatype, &
    rax, zax, r, z, dr, dz, rbnd, zbnd, raus, rmin, zlimpotential, z_xpoint, &
    psirz, psibnd
use debugger_ef, only: markloc_ef, debug

implicit none

integer :: i, j, iax, jax, i0, i1, i2, i3, i4, i5, i6, j2, j3
integer, dimension(5) :: ipath, jpath
double precision :: errtol, tolerr, br1, bz1, dbrr, dbrz, dbzr, dbzz, &
    rstart, icase, psiloc, rleft, rright, psibtmp, tin, tup, &
    x1, x2, x3, x4, x6, x7, x8, x9, x22, x33, hard_left, hard_right
double precision, dimension(i_dim5) :: xbnd, ybnd
double precision, dimension(i_dim2) :: r_temp, z_temp

call markloc_ef('find_new_boundary_part2', debug_lev=debug)

errtol = 1.e6
tolerr = dr
i0 = 0
icase = 0

iax = iaxis
jax = jaxis

CALL CPU_TIME(tin)

iter_loop: do while(errtol > tolerr) 

    if (i0 == 0) then
        do i=iax, nr2 
            if (psirz(i, jax) >= psibnd .and. r(i) <= raus) i1 = i
        enddo 
        rleft = r(i1) - dr
        rright = r(i1)
        hard_left = rax
        hard_right = r(i1)
    endif
    i0 = i0 + 1
    if (i0 >= 5200) then
        call err_catch_a  
    endif
!start contouring
    x1 = 0.5*(rleft + rright) ! first starting point is zaxis, rmid between axis and limiter
    x2 = z(jaxis) ! first starting point is zaxis, rmid between axis and limiter
    i4 = 0 !counter

    call find_fields_interp_ef_psionly_neg(x1, x2, psibtmp) !give back psi, br, bz at r0, z0
    call find_floor_index_ef(x1, x2, i2, j2)

! Move grid to coincide with boundary psi
    r_temp(1: nr2) = r(1: nr2) - r(i2) + x1
    z_temp(1: nz2) = z(1: nz2)
    i3 = 100
    j3 = 100
    ipath = 1e6
    jpath = 1e6

! Contouring
    do while((i2 /= i3) .or. (j2 /= j3))
        i4 = i4 + 1
        if (i4 == 1) then
            i3 = i2
            j3 = j2
            i5 = 0
        endif

        xbnd(i4) = r_temp(i3)
        ybnd(i4) = z_temp(j3)
        ibnd(i4) = i3
        jbnd(i4) = j3

        call find_fields_interp_ef(r_temp(i3), z_temp(j3), psiloc, &
            br1, bz1, dbrr, dbrz, dbzr, dbzz) !give back psi, br, bz at r0, z0

! Check if X-point

! Vector components of motion
        x1 = sqrt(br1**2 + bz1**2)  ! poloidal field
        x2 = br1/x1
        x3 = bz1/x1
        x22 = abs(br1/x1)
        x33 = abs(bz1/x1)

! Move point
        if (x33 > x22 .and. x3 > 0.) then
            icase = 1
            j3 = j3 + 1
! apply correction if needed
            call find_fields_interp_ef_psionly_neg(r_temp(i3)  , z_temp(j3), x6)
            call find_fields_interp_ef_psionly_neg(r_temp(i3-1), z_temp(j3), x8)
            call find_fields_interp_ef_psionly_neg(r_temp(i3+1), z_temp(j3), x9)
            x7 = x6 - psibtmp
            if (x7 > 0.) then ! at right of boundary
                icase = 11.
                x4 = abs(x7/(x6 - x8))
                if (x4 >= dr/2.) i3 = i3 - nint(x4)
            endif 
            if (x7 < 0.) then !at left of boundary
                icase = 12.
                x4 = abs(x7/(x6 - x9))
                if (x4 >= dr/2.) i3 = i3 + nint(x4)
            endif 
        endif
        if (x33 > x22 .and. x3 < 0.)  then
            icase = 2
            j3 = j3 - 1
            call find_fields_interp_ef_psionly_neg(r_temp(i3)  , z_temp(j3), x6)
            call find_fields_interp_ef_psionly_neg(r_temp(i3-1), z_temp(j3), x8)
            call find_fields_interp_ef_psionly_neg(r_temp(i3+1), z_temp(j3), x9)
! Apply correction if needed
            x7 = x6 - psibtmp
            if (x7 < 0.) then ! at right of boundary
                icase = 21.
                x4 = abs(x7/(x8 - x6))
                if (x4 >= dr/2.) i3 = i3 - nint(x4)
            endif 
            if (x7 > 0.) then !at left of boundary
                icase = 22.
                x4 = abs(x7/(x6 - x9))
                if (x4 >= dr/2.) i3 = i3 + nint(x4)
            endif 
        endif
        if (x22 > x33 .and. x2 > 0.) then
            icase = 3
            i3 = i3 + 1
            call find_fields_interp_ef_psionly_neg(r_temp(i3), z_temp(j3), x6)
            call find_fields_interp_ef_psionly_neg(r_temp(i3), z_temp(j3-1), x8)
            call find_fields_interp_ef_psionly_neg(r_temp(i3), z_temp(j3+1), x9)
! Apply correction if needed
            x7 = x6 - psibtmp
            if (x7 > 0.) then ! under of boundary
                icase = 31.
                x4 = abs(x7/(x9 - x6))
                if (x4 >= dz/2.) j3 = j3 + nint(x4)
            endif 
            if (x7 < 0.) then ! over of boundary
                icase = 32.
                x4 = abs(x7/(x8 - x6))
                if (x4 >= dz/2.) j3 = j3 - nint(x4)
            endif
        endif
        if (x22 > x33 .and. x2 < 0.) then
            icase = 4.
            i3 = i3 - 1
            call find_fields_interp_ef_psionly_neg(r_temp(i3), z_temp(j3)  , x6)
            call find_fields_interp_ef_psionly_neg(r_temp(i3), z_temp(j3-1), x8)
            call find_fields_interp_ef_psionly_neg(r_temp(i3), z_temp(j3+1), x9)

! Apply correction if needed
            x6 = -psirz(i3, j3)
            x7 = x6 - psibtmp
            if (x7 < 0.) then ! under of boundary
                icase = 41.
                x4 = abs(x7/(x9 - x6))
                if (x4 >= dz/2.) j3 = j3 + nint(x4)
            endif 
            if (x7 > 0.) then !over of boundary
                icase = 42.
                x4 = abs(x7/(x8 - x6))
                if (x4 >= dz/2.) j3 = j3 - nint(x4)
            endif 
        endif
        if (x22 == x33) then
            icase = 5.
            if (x2 > 0 .and. x3 > 0.) then
                i3 = i3 + 1
                j3 = j3 + 1
            endif
            if (x2 < 0 .and. x3 > 0.) then
                i3 = i3 - 1
                j3 = j3 + 1
            endif
            if (x2 > 0 .and. x3 < 0.) then
                i3 = i3 + 1
                j3 = j3 - 1
            endif
            if (x2 < 0 .and. x3 < 0.) then
                i3 = i3 - 1
                j3 = j3 - 1
            endif
        endif

        i3 = max(1, i3)
        i3 = min(nr2, i3)
        j3 = max(1, j3)
        j3 = min(nz2, j3)   

! Check if intersects limiter, right boundary
        i6 = nint((r_temp(i3) - rmin)/dr + 1.)
        if (zlimpotential(i6, j3) < 0.) then
            i_plasmatype = 0
            rright = r_temp(i2)
            hard_right = rright
            rleft = max(rleft - dr, hard_left)
            CYCLE iter_loop
        endif

!if it closes on itself, its the left boundary
        if ((i3 == i2) .and. (j3 == j2)) then
            rleft = r_temp(i2)
            hard_left = rleft
            rright = min(hard_right, rright + dr)
        endif

        if (i4 >= 5200) then
            rright = r_temp(i2)
            rleft = max(rleft - dr, hard_left)
            hard_right = rright
            CYCLE iter_loop
        endif

    enddo

    if (abs(rleft - rright) <= 1.e-9) EXIT

enddo iter_loop

rstart = rleft

CALL CPU_TIME(tup)
write(*, *) 'contouring', tup - tin 

nbnd = i4
psibnd = -psibtmp
rbnd(1: i4) = xbnd(1: i4)
zbnd(1: i4) = ybnd(1: i4) 
rbnd(i4+1) = rbnd(1)
zbnd(i4+1) = zbnd(1)

if (n_of_xpoints >= 1) then
!kill boundary points outside of X point
    do i=1, n_of_xpoints
        if (z_xpoint(i) < zax) then
305         continue
            do j=1, i4
                if (zbnd(j) < z_xpoint(i)) then
                    zbnd(j: i4-1) = zbnd(j+1: i4)
                    rbnd(j: i4-1) = rbnd(j+1: i4)
                    ibnd(j: i4-1) = ibnd(j+1: i4)
                    jbnd(j: i4-1) = jbnd(j+1: i4)
                    i4 = i4 - 1
                    goto 305
                endif
            enddo
        endif
    enddo
    nbnd = i4
    rbnd(i4+1) = rbnd(1)
    zbnd(i4+1) = zbnd(1)
    ibnd(i4+1) = ibnd(1)
    jbnd(i4+1) = jbnd(1)
endif

open(32, file = 'fort.607')
    do i=1, nbnd
        write(32, '(8E25.11)') rbnd(i), zbnd(i)
    enddo
close(32)

return
end subroutine find_new_boundary_part2

!------------------------------------------------------------
subroutine new_jrz_ef
! calculate new right hand side given new boundary

use ef_circuit, only: nr2, nz2, nrho, i_dim2, &
    iaxis, jaxis, iplasma, area_eff, u_n, &
    r, z, dr, dz, psiaxis, psibnd, psirz, &
    ppp_2d, ffp_2d, &
    rinner, raus, zbot, ztop, zlimpotential, jrz
use debugger_ef, only: markloc_ef, debug

implicit none

integer :: i, j
double precision :: dum1, dumc(i_dim2, i_dim2)

call markloc_ef('new_jrz_ef', debug_lev=debug)

dumc = 0. 
area_eff = 0.

! sweep from axis to exterior and fill in the current
open(6712, file='fort.6712')

write(*, *) ' 1st quadrant', ztop, rinner
i = iaxis - 1
j = jaxis
558 continue
i = i + 1
if (u_n(i, j) <= 1. .and. zlimpotential(i, j) == 1.) then 
    write(6712, *) r(i), z(j)
    call fill_in_current(r(i), z(j), dr, dz, area_eff(i, j), & 
        nrho, ppp_2d(1: nrho), ffp_2d(1: nrho), dumc(i, j), psibnd, u_n(i, j))
    goto 558
else
    i = iaxis - 1
    j = j + 1
    do
        if (u_n(i+1, j) <= 1. .and. zlimpotential(i, j) == 1.) goto 558
        if (u_n(i+1, j) > 1. .or. zlimpotential(i, j) < 1.) then
            i = i + 1
            if (i /= nr2 .and. j /= nz2 .and. z(j) <= ztop .and. r(i) <= raus) CYCLE
        endif
        EXIT
    enddo
endif

write(*, *) ' 2nd quadrant', ztop, rinner
i = iaxis
j = jaxis
658 continue
i = i - 1
if (u_n(i, j) <= 1. .and. zlimpotential(i, j) == 1.) then 
    write(6712, *) r(i), z(j), u_n(i, j), psirz(i, j)
    call fill_in_current(r(i), z(j), dr, dz, area_eff(i, j), & 
        nrho, ppp_2d(1: nrho), ffp_2d(1: nrho), dumc(i, j), psibnd, u_n(i, j))
    goto 658
else
    i = iaxis
    j = j + 1
    do
        if (u_n(i-1, j) <= 1. .and. zlimpotential(i, j) == 1.) goto 658
        if (u_n(i-1, j) > 1. .or. zlimpotential(i, j) < 1.) then
            i = i - 1
            if (i /= 1 .and. j /= nz2 .and. z(j) <= ztop .and. r(i) >= rinner) CYCLE
        endif
        EXIT
    enddo
endif

write(*, *) ' 3rd quadrant'
i = iaxis
j = jaxis - 1
758 continue
i = i - 1
if (u_n(i, j) <= 1. .and. zlimpotential(i, j) == 1.) then 
    write(6712, *) r(i), z(j)
    call fill_in_current(r(i), z(j), dr, dz, area_eff(i, j), & 
        nrho, ppp_2d(1: nrho), ffp_2d(1: nrho), dumc(i, j), psibnd, u_n(i, j))
    goto 758
else
    i = iaxis
    j = j - 1
    do
        if (u_n(i-1, j) <= 1. .and. zlimpotential(i, j) == 1.) goto 758
        if (u_n(i-1, j) > 1. .or. zlimpotential(i, j) < 1.) then
            i = i - 1
            if (i /= 1 .and. j /= 1 .and. z(j) >= zbot .and. r(i) >= rinner) CYCLE
        endif
    EXIT
    enddo
endif

write(*, *) ' 4th quadrant'
i = iaxis - 1
j = jaxis - 1
858 continue
i = i + 1
if (u_n(i, j) <= 1. .and. zlimpotential(i, j) == 1.) then 
    write(6712, *) r(i), z(j)
    call fill_in_current(r(i), z(j), dr, dz, area_eff(i, j), & 
        nrho, ppp_2d(1: nrho), ffp_2d(1: nrho), dumc(i, j), psibnd, u_n(i, j))
    goto 858
else
    i = iaxis - 1
    j = j - 1
    do
        if (u_n(i+1, j) <= 1. .and. zlimpotential(i, j) == 1.) goto 858
        if (u_n(i+1, j) > 1. .or. zlimpotential(i, j) < 1.) then
            i = i + 1
            if (i /= nr2 .and. j /= 1 .and. z(j) >= zbot .and. r(i) <= raus) CYCLE
        endif
    enddo
endif

close(6712)

!  rescale plasma
jrz(1: nr2, 1: nz2) = dumc(1: nr2, 1: nz2)

dum1 = sum(jrz*area_eff)*dr*dz
jrz = jrz/dum1*iplasma

write(*, *) 'total current ', dum1, iplasma, psiaxis, psibnd

return
end subroutine new_jrz_ef

!------------------------------------------------------------
subroutine fill_in_current(r0, z0, dr, dz, area_eff, nx, &
    ppp_2d, ffp_2d, dumc, psibnd, un)

use debugger_ef, only: markloc_ef, debug

implicit none

integer, intent(in) :: nx
double precision, intent(in) :: r0, z0, dr, dz, psibnd
double precision, intent(in), dimension(nx) :: ppp_2d, ffp_2d
double precision, intent(out) :: area_eff, un, dumc

integer :: k
double precision :: t1, t2, t3, t4, zeta
 
call markloc_ef('fill_in_current', debug_lev=debug)

call find_fields_interp_ef_psionly(r0 + dr/2., z0, t1)
call find_fields_interp_ef_psionly(r0 - dr/2., z0, t2)
call find_fields_interp_ef_psionly(r0, z0 + dz/2., t3)
call find_fields_interp_ef_psionly(r0, z0 - dz/2., t4)

area_eff = abs((max(t1, psibnd) - max(t2, psibnd))*(max(t3, psibnd) - max(t4, psibnd))/((t1 - t2)*(t3 - t4)))

if (un <= 0.) un = 0.

zeta = un*(nx - 1.) + 1.
k = max(1, floor(zeta))
if (k == nx) then
    dumc = ppp_2d(k)*r0 + ffp_2d(k)/r0
else
    dumc = ((ppp_2d(k+1)*(zeta - k) + ppp_2d(k)*(k + 1. - zeta))*r0 + &
       (ffp_2d(k+1)*(zeta - k) + ffp_2d(k)*(k + 1. - zeta))/r0)
endif

return
end subroutine fill_in_current
