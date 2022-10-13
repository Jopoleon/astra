subroutine definitions_feqis(equil_in, params, j_call, ifplasma)

use pi_grec_vars, only: GPI2, mu0
use parameters, only: type_parameters
use schemas, only: type_equilibrium
use ef_circuit
use exchange_with_astra, only: tau_circuit_ef, tau_gseq_ef, &
    activate_coil_ef, current_limit_ef, psi0_astra, psib_astra, &
    raxis_astra, zaxis_astra, nonegcurr
use debugger_ef, only: markloc_ef, debug

implicit none

integer, intent(in) :: j_call, ifplasma
type(type_parameters) , intent(in) :: params
type(type_equilibrium), intent(in) :: equil_in

integer :: i

call markloc_ef('definitions_ef_equil', debug_lev=debug)

if (j_call == 0) then
    data_dir = params%prename(1: params%kname)
    data_dir_k = params%kname
    nteta = equil_in%eqgeometry%boundary%npoints
    nrho = params%neql
    psistabR = 0.
    psistabZ = 0.
!normalized psi from 0 axis to 1 edge,  equispaced
    do i=1, nrho
        psia_2d(i) = (i - 1.)/(nrho - 1.)
    enddo
!constants

    Rgeom0 = equil_in%global_param%toroid_field%r0
    tau_circuit_ef = 0.001 !default value    
    tau_gseq_ef = 0.001     !default value
    activate_coil_ef = 1. ! when 0.,  coil is forced to 0 current
    current_limit_ef(: , 1) =  1e6 ! cant be higher than 1e6 MA
    current_limit_ef(: , 2) = -1e6 ! cant be lower than -1e6 MA
    max_iter = 520 !hardwired
    errgs2dfull = 1.e-7 !hardwired
! teta for polar grid,  goes from 0 to 2*pi-dteta,  but point nt+1 is the periodic one

    do i=1, nteta+1
        teta(i) = GPI2*(i - 1.)/(nteta + 0.)
    enddo
    psi0_astra = equil_in%profiles_1d%psi(1)
    psib_astra = equil_in%profiles_1d%psi(nrho)
    raxp = raxis_astra
    zaxp = zaxis_astra    
    psibnd = -1.e6
endif

if (ifplasma == 1) then
    btor0   = equil_in%global_param%toroid_field%b0
    iplasma = equil_in%global_param%i_plasma/1.e6
    pressure(1: nrho) = equil_in%profiles_1d%pressure(1: nrho)
    pprime(  1: nrho) = equil_in%profiles_1d%pprime(1: nrho)
    ffprime( 1: nrho) = equil_in%profiles_1d%ffprime(1: nrho)
    psigrida(1: nrho) = equil_in%profiles_1d%psi(1: nrho) !unnormalized
    psigrida(1: nrho) = (psigrida(1: nrho) - psigrida(1)) / &
                        (psigrida(nrho)    - psigrida(1)) ! normalized: 0 axis,  1 sep
    call linterp_ef(psigrida(1: nrho), ffprime(1: nrho), nrho, &
        psia_2d(1: nrho), ffp_2d(1: nrho), nrho)
    call linterp_ef(psigrida(1: nrho), pprime(1: nrho), nrho, &
        psia_2d(1: nrho), ppp_2d(1: nrho), nrho)
    call linterp_ef(psigrida(1: nrho), IPOL(1: nrho), nrho, &
        psia_2d(1: nrho), ipol_2d(1: nrho), nrho)
    call linterp_ef(psigrida(1: nrho), pressure(1: nrho), nrho, &
        psia_2d(1: nrho), pres_2d(1: nrho), nrho)
    ffp_2d = -GPI2/mu0*ffp_2d
    ppp_2d = -GPI2*1.e-6*ppp_2d
    if (nonegcurr == 1) then
        do i = 1, nrho
            ffp_2d(i) = max(0., ffp_2d(i))
            ppp_2d(i) = max(0., ppp_2d(i))
        enddo
    endif

    ipol(1: nrho) = equil_in%profiles_1d%F_dia(1: nrho)

    if (params%k_fixfree == 0) then !if 1,  comes from free boundary
        rexp(1: nteta) = equil_in%eqgeometry%boundary%r(1: nteta)
        zexp(1: nteta) = equil_in%eqgeometry%boundary%z(1: nteta)
        rexp(nteta+1)  = rexp(1)
        zexp(nteta+1)  = zexp(1)
        rbndp(1: nteta+1) = rexp(1: nteta+1)
        zbndp(1: nteta+1) = zexp(1: nteta+1)
    endif
endif

return
end subroutine definitions_feqis

!------------------------------------------------------------------
subroutine feqis_init_circ

use pi_grec_vars,  only:  GPI,  GPI2,  mu0
use ef_circuit
use green_matrix, only: greeni
use debugger_ef, only: markloc_ef, debug

implicit none

integer :: j_files, i, j, k, ii, jj, iii, jjj, ielem
integer, dimension(100) :: numeqcump
integer, dimension(300) :: tempcoilturns, tempcoilelem
integer, dimension(5200) :: identcoil, equivtmp, tatmp
integer, dimension(300, 300) :: jjelem
double precision :: dummy1, r1, r2, z1, z2, r3, z3, r4, z4, gtemp, &
    dr1, dz1, x1, x2, x3, x4, x5, x6, x7, x8, x9, greenf
double precision, dimension(300) :: tempcoilr, tempcoilz,  &
    tempcoilangle, tempcoildr, tempcoildz, areactmp
double precision, dimension(5200) :: rcetmp, zcetmp, datmp, drcetmp, dzcetmp

call markloc_ef('feqis_init_circ', debug_lev=debug)

open(32, file=data_dir(1: data_dir_k)//'refit_ef.dat')
    read(32, *) nr_of_fit_parameters   ! # fit parameters from rmag to zxp_fit
    read(32, *) rmag_fit               ! R mag axis
    read(32, *) zmag_fit               ! Z mag axis
    read(32, *) k_fit                  ! elongation at mag axis
    read(32, *) rxp_fit                ! X point R
    read(32, *) zxp_fit                ! X point Z
close(32)

open(32, file=data_dir(1: data_dir_k)//'data.dat')
    read(32, *) j_files
    read(32, *) nr
    read(32, *) nz
    read(32, *) rmin
    read(32, *) rmax
    read(32, *) zmin
    read(32, *) zmax
    read(32, *) alpsep
close(32)

!compatibility with spider
nr = nr - 1  !because for spider its 65, 65 for example,  but here its 64, 64
nz = nz - 1

!define grid
nr2 = nr + 2
nz2 = nz + 2
nr1 = nr + 1
nz1 = nz + 1
do i=1, nr2
     r(i) = rmin + (i - 1.)*(rmax - rmin)/nr1   ! computational domain is r(2: nr+1),  boundaries are r(1) and r(nr+2)
enddo
do i=1, nz2
    z(i) = zmin + (i - 1.)*(zmax - zmin)/nz1
enddo
rcomp(1: nr) = r(2: nr+1)
zcomp(1: nz) = z(2: nz+1)
dr = r(2) - r(1)
dz = z(2) - z(1)

!load coils
r_cond = 0.
z_cond = 0.
numeqcump=0
open(32, file=data_dir(1: data_dir_k)//'coil.dat')
    read(32, *) ncoils
    do i=1, ncoils
        read(32, *) nelemcoil(i)
        read(32, *) rcoil(i), zcoil(i), drcoil(i), dzcoil(i), dummy1, &
            anglecoil(i), curcoil(i), mturns(i), mequivalence(i)
        curconduc(mequivalence(i)) = curcoil(i) !assign current to conductor
        r_cond(mequivalence(i)) = r_cond(mequivalence(i)) + rcoil(i) !assign current to conductor
        z_cond(mequivalence(i)) = z_cond(mequivalence(i)) + zcoil(i) !assign current to conductor
        numeqcump(mequivalence(i)) = numeqcump(mequivalence(i)) + 1
    enddo
    anglecoil = anglecoil/180.*GPI
close(32)
nconduc = maxval(mequivalence(1: ncoils))
r_cond(1: nconduc) = r_cond(1: nconduc)/numeqcump(1: nconduc)
z_cond(1: nconduc) = z_cond(1: nconduc)/numeqcump(1: nconduc)

!load coilres
open(32, file=data_dir(1: data_dir_k)//'coilres.dat')
    read(32, *) nreseqcoil
    do i=1, nreseqcoil
        read(32, *) j
        read(32, *) resconduc(i, 1: nreseqcoil)
    enddo
close(32)
nactive = nconduc
write(*, *) nactive 

!load limiter
open(32, file=data_dir(1: data_dir_k)//'limpnt.dat')
read(32, *) nlimiter
    do i=1, nlimiter
        read(32, *) limiterr(i), limiterz(i)
    enddo
close(32)

!load first wall (not working yet)
open(32, file=data_dir(1: data_dir_k)//'blanfw.dat')
    read(32, *) reswall
    read(32, *) nfirstwall
    if (nfirstwall >= 1) then
        do i=1, nfirstwall
            read(32, *) j, rwall(i)
        enddo
    endif
close(32)

!load blanket (works)
open(32, file=data_dir(1: data_dir_k)//'blanbp.dat')
    read(32, *) resblan
    read(32, *) nblanket
    nelemblanket = 9  !nelemblanket sub element of a blanket element,  for now hardwired width to 10 cm
    if (nblanket >= 1) then
        do i=1, nblanket
            j = 2*i - 1
            read(32, *) jjj, x1, x2, x3, x4, x5, x6
            x7 = x3 - x1
            x8 = x4 - x2
            rblan(j)   = x1 + 1./4.*x7
            rblan(2*i) = x1 + 3./4.*x7
            zblan(j)   = x2 + 1./4.*x8
            zblan(2*i) = x2 + 3./4.*x8
            x9 = sqrt(x7**2. + x8**2.)
            areablan(j)   = 0.1*x9/2
            areablan(2*i) = 0.1*x9/2
        enddo
        write(*, *) resblan, nblanket, x1, x2, x3, x4, x5, x6, x7, x8, &
            rblan(1), zblan(1), areablan(1)
        nblanket = 2*nblanket
        do i=1, nblanket
            nconduc = nconduc + 1
            curconduc(nconduc) = 0.
            resconduc(nconduc, nconduc) = resblan
            r_cond(nconduc) = rblan(i)
            z_cond(nconduc) = zblan(i)
        enddo
    endif
close(32)

!load passive conduc  (works)
open(32, file=data_dir(1: data_dir_k)//'blanbpc.dat')
    read(32, *) nblanketpc
    nelemblanketpc = 9  !nelemblanketpc sub element of a blanket element
    if (nblanketpc >= 1) then
        do i=1, nblanketpc
            read(32, *) rblanpc(i), zblanpc(i), resblanpc(i), &
                areablanpc(i), curblanpc(i)
            nconduc = nconduc + 1
            curconduc(nconduc) = curblanpc(i)
            resconduc(nconduc, nconduc) = resblanpc(i)
            r_cond(nconduc) = rblanpc(i)
            z_cond(nconduc) = zblanpc(i)
        enddo
    endif
close(32)

! currents are in MA!
npassive = nconduc - nactive

if (j_files == 0) then
    
! create circuit structure
    nconduc = 0    
    ielem = 0
    tatmp = 0
! coils first
    identcoil(1: ncoils) = 1
    do i=1, ncoils
        k = 0
        if (identcoil(i) == 1) then
            nconduc = nconduc + 1
            k = k + 1
            tempcoilr(k) = rcoil(i)
            tempcoilz(k) = zcoil(i)
            tempcoilangle(k) = anglecoil(i)
            tempcoildr(k) = drcoil(i)  
            tempcoildz(k) = dzcoil(i)
            tempcoilturns(k) = mturns(i)
            tempcoilelem(k) = nelemcoil(i)
            tatmp(nconduc) = tatmp(nconduc) + mturns(i)
            do j = 1, ncoils
                if ((mequivalence(i) == mequivalence(j)) .and. (i /= j)) then
                    identcoil(j) = 0
                    k = k + 1
                    tempcoilr(k) = rcoil(j)
                    tempcoilz(k) = zcoil(j)
                    tempcoilangle(k) = anglecoil(j)
                    tempcoildr(k) = drcoil(j)
                    tempcoildz(k) = dzcoil(j)
                    tempcoilturns(k) = mturns(j)
                    tempcoilelem(k) = nelemcoil(j)
                    tatmp(nconduc) = tatmp(nconduc) + mturns(j)
                endif
            enddo

! Analysis coil
            do j=1, k
                iii = nint(sqrt(tempcoilelem(j) + 1.d-6))
                dr1 = tempcoildr(j)/iii
                dz1 = tempcoildz(j)/iii
                x1 = tempcoildz(j)/sin(tempcoilangle(j))
                x2 = x1*cos(tempcoilangle(j))
                r1 = tempcoilr(j) - tempcoildr(j)/2. - x2/4.
                r2 = tempcoilr(j) + tempcoildr(j)/2. - x2/4.
                z1 = tempcoilz(j) - tempcoildz(j)/2.
                z2 = tempcoilz(j) - tempcoildz(j)/2.
                r3 = r2 + x2
                r4 = r1 + x2
                z3 = tempcoilz(j) + tempcoildz(j)/2.
                z4 = tempcoilz(j) + tempcoildz(j)/2.
                do jjj=1, iii
                    do jj=1, iii
                        ielem = ielem + 1
                        rcetmp(ielem) = r1 + (jj - 1.)*dr1 + (jjj - 1.)*dz1*cos(tempcoilangle(j)) + dr1/2.     !center
                        zcetmp(ielem) = z1 + (jjj - 1.)*dz1 + dz1/2. !center
                        drcetmp(ielem) = dr1
                        dzcetmp(ielem) = dz1
                        datmp(ielem) = dr1*dz1/sin(tempcoilangle(j)) !area including turns
                        equivtmp(ielem) = nconduc
                    enddo
                enddo !cycles over 1 single coil element
            enddo !cycle over equivalent coils
        endif
    enddo !cycle over equivalent coils

    write(333, *) rcetmp(1: ielem), zcetmp(1: ielem), 0, 0

! Now blanket
    write(*, *) nblanketpc
    if (nblanketpc >= 1) then
        do j=1, nblanketpc
            nconduc = nconduc + 1
            tatmp(nconduc) = 1 !area including turns
            iii = nint(sqrt(nelemblanketpc + 1.d-6))
            dr1 = sqrt(areablanpc(j))/iii
            dz1 = dr1
            x1 = dz1
            x2 = 0
            r1 = rblanpc(j) - dr1/2.
            r2 = rblanpc(j) + dr1/2.
            z1 = zblanpc(j) - dz1/2.
            z2 = zblanpc(j) - dz1/2.
            r3 = r2
            r4 = r1
            z3 = zblanpc(j) + dz1/2.
            z4 = zblanpc(j) + dz1/2.
            do jjj=1, iii
                do jj=1, iii
                    ielem = ielem + 1
                    rcetmp(ielem) = r1 + (jj  - 1./2.)*dr1   ! center
                    zcetmp(ielem) = z1 + (jjj - 1./2.)*dz1   ! center
                    drcetmp(ielem) = dr1
                    dzcetmp(ielem) = dz1
                    datmp(ielem) = dr1*dz1 !area including turns
                    equivtmp(ielem) = nconduc
                enddo
            enddo !cycles over 1 single coil element
        enddo !cycle over equivalent coils
    endif

! Now blanket
    write(*, *) nblanket
    if (nblanket >= 1) then
        do j=1, nblanket
            nconduc = nconduc + 1
            tatmp(nconduc) = 1 !area including turns
            iii = nint(sqrt(nelemblanket + 1.d-6))
            dr1 = sqrt(areablan(j))/iii
            dz1 = dr1
            x1 = dz1
            x2 = 0
            r1 = rblan(j) - dr1/2.
            r2 = rblan(j) + dr1/2.
            z1 = zblan(j) - dz1/2.
            z2 = zblan(j) - dz1/2.
            r3 = r2
            r4 = r1
            z3 = zblan(j) + dz1/2.
            z4 = zblan(j) + dz1/2.
            do jjj=1, iii
                do jj=1, iii
                    ielem = ielem + 1
                    rcetmp(ielem) = r1 + (jj  - 1./2.)*dr1 ! center
                    zcetmp(ielem) = z1 + (jjj - 1./2.)*dz1 ! center
                    drcetmp(ielem) = dr1
                    dzcetmp(ielem) = dz1
                    datmp(ielem) = dr1*dz1 !area including turns
                    equivtmp(ielem) = nconduc
                    write(*, *) 'ielem', ielem, rcetmp(ielem), zcetmp(ielem)
                enddo
            enddo !cycles over 1 single coil element
        enddo !cycle over equivalent coils
    endif

!calculate self-inductances
    write(*, *) 'self', nconduc, npassive, nactive, ielem
    indconduc = 0.
    areactmp = 0.
    identcoil = 1
    jjelem = 0
    do i=1, ielem
        iii = equivtmp(i)
        areactmp(iii) = areactmp(iii) + datmp(i)
        do j=1, ielem
            if ((equivtmp(j) == iii) .and. (i /= j)) then
                call green_function(rcetmp(i), zcetmp(i), rcetmp(j), &
                     zcetmp(j), gtemp)
                indconduc(iii, iii) = indconduc(iii, iii) + mu0/GPI*gtemp
            endif
            if ((equivtmp(j) == iii) .and. (i == j)) then
                jjelem(iii, iii) = jjelem(iii, iii) + 1
                call green_function_identity(rcetmp(i), gtemp, &
                    drcetmp(i), dzcetmp(i))
                indconduc(iii, iii) = indconduc(iii, iii) + mu0/GPI*gtemp/2
            endif
        enddo
    enddo
    do i=1, nconduc
        indconduc(i, i) = GPI2*indconduc(i, i)*tatmp(i)**2./jjelem(i, i)**2.
    enddo

!calculate mutual-inductances
    write(*, *) 'mutual'
    identcoil = 1
    do i=1, ielem
        iii = equivtmp(i)
        do j=1, ielem
            if ((equivtmp(j) /= iii)) then
                jjelem(equivtmp(j), iii) = jjelem(equivtmp(j), iii) + 1
                call green_function(rcetmp(i), zcetmp(i), rcetmp(j), &
                    zcetmp(j), gtemp)
                indconduc(iii, equivtmp(j)) = indconduc(iii, equivtmp(j)) + mu0/GPI*gtemp
                indconduc(equivtmp(j), iii) = indconduc(iii, equivtmp(j))
            endif
        enddo
    enddo

    do j=1, nconduc
        do i=1, nconduc
            if (i /= j) then
                indconduc(i, j) = GPI2*indconduc(i, j)*tatmp(i)*tatmp(j)/2./jjelem(i, j)
            endif
        enddo
    enddo

! calculate grid-inductances
    write(*, *) 'grid'
    greeni = 0.
    do i=1, ielem
        iii = equivtmp(i)
        do jj = 1, nz2
            do ii = 1, nr2
                call green_function(rcetmp(i), zcetmp(i), r(ii), z(jj), gtemp)
                greeni(ii, jj, iii) = greeni(ii, jj, iii) + mu0/GPI*gtemp
            enddo
        enddo
    enddo
    do i=1, nconduc
        greeni(: , : , i) = greeni(: , : , i)*tatmp(i)/jjelem(i, i)
    enddo

!generate zlimpotential
    do j=1, nz2
        do i=1, nr2   
            call generate_zlim_potential(r(i), z(j), nlimiter, &
                limiterR(1: nlimiter), limiterZ(1: nlimiter), zlimpotential(i, j))
        enddo
    enddo

! write everything on file
    open(32, file=data_dir(1: data_dir_k)//'induc_matrix.dat')
        write(32, *) nconduc
        do i=1, nconduc
            write(32, *) indconduc(i, 1: nconduc)
        enddo
    close(32)

    open(32, file = data_dir(1: data_dir_k)//'resistance_matrix.dat')
        write(32, *) nconduc
        do i=1, nconduc
            write(32, *) resconduc(i, 1: nconduc)
        enddo
    close(32)

    open(32, file=data_dir(1: data_dir_k)//'greeni_matrix.dat')
        do i=1, nconduc
            do j=1, nr2
                write(32, *) greeni(j, 1: nz2, i)
            enddo
        enddo
    close(32)

    open(32, file=data_dir(1: data_dir_k)//'zlim_potential.dat')
        do jj=1, nz2
            do ii=1, nr2
                write(32, *) zlimpotential(ii, jj)
            enddo
        enddo
    close(32)

    write(*, *) 'bound'
    open(32, file=data_dir(1: data_dir_k)//'green_boundary.dat')
        write(32, *) nint((2.*nr2 + 2.*nz2)*(2.*nr1 + 2.*nz1))
! lower side
        do i=1, nr2
            do j=1, nr1
                call green_function(r(i), z(1), r(j) + dr/2., z(1), greenf)
                write(32, *) greenf
            enddo
!right side
            do j=1, nz1
                call green_function(r(i), z(1), r(nr2), z(j) + dz/2., greenf)
                write(32, *) greenf
            enddo
! upper side
            do j=1, nr1
                call green_function(r(i), z(1), r(j) + dr/2., z(nz2), greenf)
                write(32, *) greenf
            enddo
!left side
            do j=1, nz1
                call green_function(r(i), z(1), r(1), z(j) + dz/2., greenf)
                write(32, *) greenf
            enddo
        enddo

! right side
        do i=1, nz2
            do j=1, nr1
                call green_function(r(nr2), z(i), r(j) + dr/2., z(1), greenf)
                write(32, *) greenf
            enddo
!right side
            do j=1, nz1
                call green_function(r(nr2), z(i), r(nr2), z(j) + dz/2., greenf)
                write(32, *) greenf
            enddo
! upper side
            do j=1, nr1
                call green_function(r(nr2), z(i), r(j) + dr/2., z(nz2), greenf)
                write(32, *) greenf
            enddo
!left side
            do j=1, nz1
                call green_function(r(nr2), z(i), r(1), z(j) + dz/2., greenf)
                write(32, *) greenf
            enddo
        enddo

! upper side
        do i=1, nr2
            do j=1, nr1
                call green_function(r(i), z(nz2), r(j) + dr/2., z(1), greenf)
                write(32, *) greenf
            enddo
!right side
            do j=1, nz1
                call green_function(r(i), z(nz2), r(nr2), z(j) + dz/2., greenf)
                write(32, *) greenf
            enddo
! upper side
            do j=1, nr1
                call green_function(r(i), z(nz2), r(j) + dr/2., z(nz2), greenf)
                write(32, *) greenf
            enddo
!left side
            do j=1, nz1
                call green_function(r(i), z(nz2), r(1), z(j) + dz/2., greenf)
                write(32, *) greenf
            enddo
        enddo

! left side
        do i=1, nz2
            do j=1, nr1
                call green_function(r(1), z(i), r(j) + dr/2., z(1), greenf)
                write(32, *) greenf
            enddo
!right side
            do j=1, nz1
                call green_function(r(1), z(i), r(nr2), z(j) + dz/2., greenf)
                write(32, *) greenf
            enddo
! upper side
            do j=1, nr1
                call green_function(r(1), z(i), r(j) + dr/2., z(nz2), greenf)
                write(32, *) greenf
            enddo
!left side
            do j=1, nz1
                call green_function(r(1), z(i), r(1), z(j) + dz/2., greenf)
                write(32, *) greenf
            enddo
        enddo
    close(32)

    write(*, *) 'files done, please restart'
    call err_catch_a

else
! load everything from file
    open(32, file=data_dir(1: data_dir_k)//'induc_matrix.dat')
        read(32, *) nconduc
        do i=1, nconduc
            read(32, *) indconduc(i, 1: nconduc)
        enddo
    close(32)
    open(32, file=data_dir(1: data_dir_k)//'resistance_matrix.dat')
        read(32, *) nconduc
        do i=1, nconduc
            read(32, *) resconduc(i, 1: nconduc)
        enddo
    close(32)
    open(32, file=data_dir(1: data_dir_k)//'greeni_matrix.dat')
    do i=1, nconduc
        do j=1, nr2
            read(32, *) greeni(j, 1: nz2, i)
        enddo
    enddo
    close(32)

    open(32, file=data_dir(1: data_dir_k)//'zlim_potential.dat')
        do jj=1, nz2
            do ii=1, nr2
                read(32, *) zlimpotential(ii, jj)
            enddo
        enddo
    close(32)

!use spider values !! temporary,  to improve my calculations
    open(32, file=data_dir(1: data_dir_k)//'ppind_mat.wr')
        read(32, *) nconduc
        read(32, *) ((indconduc(i, j), i=1, nconduc), j=1, nconduc)
    close(32)
    indconduc = indconduc*GPI2
!use spider values
    open(32, file=data_dir(1: data_dir_k)//'res_mat.wr')
        read(32, *) ((resconduc(i, j), i=1, nconduc), j=1, nconduc)
    close(32)
    resconduc = resconduc*GPI2

    open(32, file=data_dir(1: data_dir_k)//'green_boundary.dat')
        read(32, *) ngbnd
        read(32, *) green_bnd_f(1: ngbnd)
    close(32)

endif

return
end subroutine feqis_init_circ
