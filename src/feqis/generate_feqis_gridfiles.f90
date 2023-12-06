subroutine generate_files_feqis(data_dir2, machine)

use feqis_tools, only: green_function, green_function_identity, &
    green_function_non_identity, green_function_includingsamepoint

implicit none

integer, parameter :: ncoils_max=300, nplas_max=300, nblanket_max=200, n_max=5200, nlim_max=500
double precision, parameter :: GPI=3.1415926, GPI2=2.*GPI, mu0=0.4*GPI

character(len=*), intent(in) :: data_dir2, machine

integer :: j_files, i, j, k, ii, jj, iii, jjj, ielem, &
    nr, nr1, nr2, nz, nz1, nz2, nlimiter, &
    ncoils, nreseqcoil, nconduc, nblocks, npassive, nactive, &
    nblanket, nblanketpc, nelemblanketpc, nelemblanket, &
    ilim_minr, ilim_maxr, ilim_minz, ilim_maxz
integer, dimension(ncoils_max) :: tempcoilturns, tempcoilelem, tempnnc, &
    nelemcoil, mturns, mequivalence
integer, dimension(100) :: numeqcump
integer, dimension(n_max) :: identcoil, nctype, equivtmp, equivforce

double precision :: dummy1, r1, r2, z1, z2, r3, z3, r4, z4, gtemp, dr1, dz1, &
    x1, x2, x3, x4, x5, x6, x7, x8, x9, greenf, ssfw, &
    Rmin, Rmax, Zmin, Zmax, alpsep, dR, dZ, &
    lim_minR, lim_maxR, lim_minZ, lim_maxZ, &
    resblan, widthblan
double precision, dimension(ncoils_max) :: tempcoilr, tempcoilz, &
    tempcoilangle, tempcoilangleh, tempcoildr, tempcoildz, areactmp, &
    r_cond, z_cond, curconduc, Rcoil, Zcoil, dRcoil, dZcoil, anglecoil, anglehcoil, &
    Rblan, Zblan, areablan, Rblanpc, Zblanpc, resblanpc, areablanpc, curblanpc
double precision, dimension(nplas_max) :: Rcomp, Zcomp, R, Z
double precision, dimension(nlim_max) :: limiterR, limiterZ
double precision, dimension(nblanket_max) :: dhoriz, dvert !blanket elements lengths
double precision, dimension(n_max) :: rcetmp, zcetmp, datmp, &
    drcetmp, dzcetmp, tatmp
double precision, dimension(ncoils_max, ncoils_max) :: resconduc, indconduc, &
    dgreenirj, dgreenizj
integer, dimension(nplas_max, nplas_max) :: zlimpotential
double precision, dimension(nplas_max, nplas_max, ncoils_max) :: greeni, dgreenirpl, dgreenizpl
double precision, dimension(ncoils_max) :: cos1,cos2,sin1,sin2,Det

character(len=120) :: fname, dumstring1

fname = trim(data_dir2)//'/machine_description_in.'//trim(machine)

open(32, file=trim(fname))

!general grid file
read(32, *) dumstring1
read(32, *) nr2
read(32, *) nz2
read(32, *) rmin
read(32, *) rmax
read(32, *) zmin
read(32, *) zmax
read(32, *) alpsep

!define grid
nr1 = nr2 - 1
nz1 = nz2 - 1
nr  = nr1 - 1
nz  = nz1 - 1
do i=1, nr2
    r(i) = rmin + (i - 1.)*(rmax - rmin)/nr1     ! computational domain is r(2:nr+1), boundaries are r(1) and r(nr+2)
enddo
do i=1, nz2
    z(i) = zmin + (i - 1.)*(zmax - zmin)/nz1
enddo
rcomp(1:nr) = r(2:nr1)
zcomp(1:nz) = z(2:nz1)
dr = r(2) - r(1)
dz = z(2) - z(1)

!load coils
!Definitions:
! rcoil, zcoil : geometric centroid of the coil (R,Z)
! drcoil : length of the "horizontal" side l1
! dzcoil : length of the "vertical" side l2
! Area of romboid is: A = l1*l2*Det    ! Det = cos(angleh)*sin(angle) - cos(angle)*sin(angleh). For angleh=0 and angle = 90, Det = 1.


r_cond = 0.
z_cond = 0.
numeqcump = 0
read(32, *) dumstring1
read(32, *) ncoils
do i=1, ncoils
    read(32, *) rcoil(i), zcoil(i), drcoil(i), dzcoil(i), anglehcoil(i), anglecoil(i), &
         mturns(i), mequivalence(i), nelemcoil(i)
    r_cond(mequivalence(i)) = r_cond(mequivalence(i)) + rcoil(i) !assign current to conductor
    z_cond(mequivalence(i)) = z_cond(mequivalence(i)) + zcoil(i) !assign current to conductor
    numeqcump(mequivalence(i)) = numeqcump(mequivalence(i)) + 1
		cos1(i)=cos(anglehcoil(i))
		cos2(i)=cos(anglecoil(i))
		sin1(i)=sin(anglehcoil(i))
		sin2(i)=sin(anglecoil(i))
		Det(i)=cos1(i)*sin2(i)-cos2(i)*sin1(i)
enddo
anglehcoil = anglehcoil/180.*GPI
anglecoil = anglecoil/180.*GPI
nconduc = maxval(mequivalence(1:ncoils))
r_cond(1:nconduc) = r_cond(1:nconduc)/numeqcump(1:nconduc)
z_cond(1:nconduc) = z_cond(1:nconduc)/numeqcump(1:nconduc)


write(*,*) anglecoil(1:ncoils),anglehcoil(1:ncoils)

nblocks = 0

!load coilres
read(32, *) dumstring1
read(32, *) nreseqcoil
do i=1, nreseqcoil
    read(32, *) resconduc(i, 1:nreseqcoil)
enddo
nactive = nconduc
write(*, *) nactive

!load limiter
read(32, *) dumstring1
read(32, *) nlimiter
do i=1, nlimiter
    read(32, *) limiterr(i), limiterz(i)
enddo
read(32, *) lim_maxR, lim_minR, lim_maxZ, lim_minZ

!forces limiter to adapt to grid points !EFable test, but this should be better than leaving it floating
do i=1, nlimiter
    j = nint((limiterr(i) - rmin)/dr + 1.)
    k = nint((limiterz(i) - zmin)/dz + 1.)
    limiterr(i) = r(j)
    limiterz(i) = z(k)
enddo

ilim_maxR = 1 + nint((lim_maxR - r(1))/dr)
lim_maxR = r(ilim_maxR)

ilim_minR = 1 + nint((lim_minR - r(1))/dr)
lim_minR = r(ilim_minR)

ilim_maxZ = 1 + nint((lim_maxZ - z(1))/dz)
lim_maxZ = z(ilim_maxZ)

ilim_minZ = 1 + nint((lim_minZ - z(1))/dz)
lim_minZ = z(ilim_minZ)

!load blanket (works)
read(32, *) dumstring1
read(32, *) resblan, widthblan
read(32, *) nblanket
nelemblanket = 9  !nelemblanket sub element of a blanket element, for now hardwired width to 10 cm
if (nblanket >= 1) then
    ssfw = 0.
    do i=1, nblanket
        j = 2*i - 1
        read(32, *) jjj, x1, x2, x3, x4, x5, x6
        x7 = x3 - x1
        x8 = x4 - x2
        rblan(j)   = x1 + 1./4.*x7
        rblan(2*i) = x1 + 3./4.*x7
	resblanpc(j)   = x5
	resblanpc(2*i) = x5
	curblanpc(j)   = x6
	curblanpc(2*i) = x6
        zblan(j)   = x2 + 1./4.*x8
        zblan(2*i) = x2 + 3./4.*x8
        x9 = sqrt(x7**2 + x8**2)
        dhoriz(j)   = x7/2.
        dvert(j)    = x8/2.
        dhoriz(2*i) = x7/2.
        dvert(2*i)  = x8/2.
        areablan(j)   = x9*widthblan
        areablan(2*i) = x9*widthblan
        ssfw = ssfw + areablan(j)/rblan(j) + areablan(2*i)/rblan(2*i)
    enddo
    nblanket = 2*nblanket
    do i=1, nblanket
        nconduc = nconduc + 1
        curconduc(nconduc) = curblanpc(i)
        r_cond(nconduc) = rblan(i)
        z_cond(nconduc) = zblan(i)
	if (resblan > 0.) resconduc(nconduc, nconduc) = resblan*r_cond(nconduc)/areablan(i)*ssfw
	if (resblan < 0.) resconduc(nconduc, nconduc) = resblanpc(i)
    enddo
endif
resblanpc = 0.
curblanpc = 0.

!load passive conduc  (works)
read(32, *) dumstring1
read(32, *) nblanketpc
nelemblanketpc = 9  !nelemblanketpc sub element of a blanket element, hardwired to 9 for now
if (nblanketpc >= 1) then
    do i=1, nblanketpc
        read(32, *) rblanpc(i), zblanpc(i), resblanpc(i), areablanpc(i), curblanpc(i)
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
nblocks = ncoils + npassive

!create circuit structure
nconduc = 0
nblocks = 0
ielem = 0
tatmp = 0.
! coils first
identcoil(1:ncoils) = 1
do i=1, ncoils
    nblocks = nblocks + 1
    k = 0
    if (identcoil(i) == 1) then
        nconduc = nconduc + 1
        k = k + 1
        tempcoilr(k) = rcoil(i)
        tempcoilz(k) = zcoil(i)
        tempcoilangleh(k) = anglehcoil(i)
        tempcoilangle(k) = anglecoil(i)
        tempcoildr(k) = drcoil(i)
        tempcoildz(k) = dzcoil(i)
        tempcoilturns(k) = mturns(i)
        tempcoilelem(k) = nelemcoil(i)
        tempnnc(k) = i
        do j=i+1, ncoils
            if ((mequivalence(i) == mequivalence(j))) then
                identcoil(j) = 0
                k = k + 1
                tempcoilr(k) = rcoil(j)
                tempcoilz(k) = zcoil(j)
                tempcoilangleh(k) = anglehcoil(j)
                tempcoilangle(k) = anglecoil(j)
                tempcoildr(k) = drcoil(j)
                tempcoildz(k) = dzcoil(j)
                tempcoilturns(k) = mturns(j)
                tempcoilelem(k) = nelemcoil(j)
                tempnnc(k) = j
            endif
        enddo
! analysis coil
!  double precision rcetmp(1200), zcetmp(1200), datmp(1200), equivtmp(1200)
        do j=1, k
            r1 = tempcoilr(j) - 0.5*(tempcoildr(j)*cos(tempcoilangleh(j)) + tempcoildz(j)*cos(tempcoilangle(j)))
            z1 = tempcoilz(j) - 0.5*(tempcoildr(j)*sin(tempcoilangleh(j)) + tempcoildz(j)*sin(tempcoilangle(j)))
           	x3 = tempcoildr(j)
            x4 = tempcoildz(j)
            x5 = x3*x4*(cos(tempcoilangleh(j))*sin(tempcoilangle(j))-cos(tempcoilangle(j))*sin(tempcoilangleh(j)))
            x6 = x5/tempcoilelem(j)

            iii = nint(x3/sqrt(x6) + 0.5)
            jjj = nint(x4/sqrt(x6) + 0.5)
						tempcoilelem(j)=iii*jjj
            x3 = x3/iii !dl1
            x4 = x4/jjj !dl2
	write(*,*) i,j,r1,z1,iii,jjj,x3,x4,tempcoildr(j),tempcoildz(j)
						
            do jj=1, jjj
                do ii=1, iii
                    ielem = ielem + 1
                    rcetmp(ielem) = r1 + 0.5*(x3*cos(tempcoilangleh(j)) + x4*cos(tempcoilangle(j)))+(ii-1.)*x3*cos(tempcoilangleh(j)) + (jj-1.)*x4*cos(tempcoilangle(j))   !center
                    zcetmp(ielem) = z1 + 0.5*(x3*sin(tempcoilangleh(j)) + x4*sin(tempcoilangle(j)))+(ii-1.)*x3*sin(tempcoilangleh(j)) + (jj-1.)*x4*sin(tempcoilangle(j)) !center
                    drcetmp(ielem) = x3   !dl1
                    dzcetmp(ielem) = x4 !dl2
                    datmp(ielem) = x3*x4*(cos(tempcoilangleh(j))*sin(tempcoilangle(j))-cos(tempcoilangle(j))*sin(tempcoilangleh(j))) !area 
                    equivtmp(ielem) = nconduc
                    equivforce(ielem) = tempnnc(j)
                    tatmp(ielem) = (0.+tempcoilturns(j))/(0.+tempcoilelem(j))+0. ! turns
                    write(*, '(9E25.11)') ielem + 0., rcetmp(ielem), zcetmp(ielem), 0. + equivtmp(ielem), 0. + equivforce(ielem),tempcoilturns(j)+0.,tempcoilelem(j)+0.,tatmp(ielem)
                    nctype(ielem) = 2 !rectangular coil block
                enddo
            enddo !cycles over 1 single coil element
        enddo !cycle over equivalent coils
    endif
enddo !cycle over equivalent coils

!now blanket
if (nblanketpc >= 1) then
    do j=1, nblanketpc
        nconduc = nconduc + 1
        nblocks = nblocks + 1
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
                rcetmp(ielem) = r1 + (jj  - 0.5)*dr1   !center
                zcetmp(ielem) = z1 + (jjj - 0.5)*dz1 !center
                drcetmp(ielem) = dr1   !center
                dzcetmp(ielem) = dz1 !center
                datmp(ielem) = dr1*dz1 !area including turns
                equivtmp(ielem) = nconduc
                equivforce(ielem) = nblocks
                tatmp(ielem) = dr1*dz1/areablanpc(j) !area including turns
                nctype(ielem) = 2 !rectangular passive element
            enddo
        enddo !cycles over 1 single coil element
    enddo !cycle over equivalent coils
endif

!now blanket
if (nblanket >= 1) then
    do j=1, nblanket
        nconduc = nconduc + 1
        nblocks = nblocks + 1
        ielem = ielem + 1
        rcetmp(ielem) = rblan(j)   !center
        zcetmp(ielem) = zblan(j) !center
        drcetmp(ielem) = dhoriz(j)   !center
        dzcetmp(ielem) = dvert(j) !center
        datmp(ielem) = dhoriz(j)*dvert(j) !area including turns
        equivtmp(ielem) = nconduc
        equivforce(ielem) = nblocks
        tatmp(ielem) = 1. !area including turns
        nctype(ielem) = 1 !segment block
    enddo !cycle over equivalent coils
endif

!calculate self-inductances ! Fable, trying to match SPIDER inductances....
write(*, *) 'self', nconduc, npassive, nactive, ielem
indconduc = 0.
areactmp = 0.
identcoil = 1
do i=1, ielem
    iii = equivtmp(i)
    do j=1, ielem
        if (iii == equivtmp(j)) then
            if (i /= j) then
                gtemp = green_function_non_identity(rcetmp(i), zcetmp(i), rcetmp(j), zcetmp(j), &
                    drcetmp(i), dzcetmp(i), drcetmp(j), dzcetmp(j), nctype(i), nctype(j))
                indconduc(iii, iii) = indconduc(iii, iii) + mu0/GPI*gtemp*tatmp(i)*tatmp(j)/1.
            else
                gtemp = green_function_identity(rcetmp(i), zcetmp(i), drcetmp(i), dzcetmp(i), nctype(i))
                indconduc(iii, iii) = indconduc(iii, iii) + mu0/GPI*gtemp*tatmp(i)*tatmp(j)/2.
            endif
        endif
    enddo
enddo
do i=1, nconduc
    indconduc(i, i) = GPI2*indconduc(i, i)
enddo
write(*, *) indconduc(1, 1)/GPI2

!calculate mutual-inductances
write(*, *) 'mutual'
identcoil = 1
do i=1, ielem
    iii = equivtmp(i)
    do j=1, ielem
        if ((equivtmp(j) /= iii)) then
            gtemp = green_function_non_identity(rcetmp(i), zcetmp(i), rcetmp(j), zcetmp(j), &
                drcetmp(i), dzcetmp(i), drcetmp(j), dzcetmp(j), nctype(i), nctype(j))
            indconduc(iii, equivtmp(j)) = indconduc(iii, equivtmp(j)) + mu0/GPI*gtemp*tatmp(i)*tatmp(j)
        endif
    enddo
enddo
do j=1, nconduc
    do i=1, nconduc
        if (i /= j) then
            indconduc(i, j) = GPI2*indconduc(i, j)
        endif
    enddo
enddo

write(*, *) indconduc(1, 2)/GPI2

!calculate grid-inductances
write(*, *) 'grid'
greeni = 0.
do i=1, ielem
    iii = equivtmp(i)
    do jj=1, nz2
        do ii=1, nr2
            gtemp = green_function_non_identity(rcetmp(i), zcetmp(i), r(ii), z(jj), &
                drcetmp(i), dzcetmp(i), dr, dz, nctype(i), 2)
            greeni(ii, jj, iii) = greeni(ii, jj, iii) + mu0/GPI*gtemp*tatmp(i)
        enddo
    enddo
enddo

!calculate inter-blocks forces, to check if these are right btw
dgreenirpl = 0.
dgreenizpl = 0.
dgreenirj = 0.
dgreenizj = 0.

do i=1, ielem
    iii = equivforce(i)
    do j=1, ielem
        if ((equivforce(j) /= iii)) then
            x1 = green_function(rcetmp(i) + dr/2., zcetmp(i), rcetmp(j), zcetmp(j))
            x2 = green_function(rcetmp(i) - dr/2., zcetmp(i), rcetmp(j), zcetmp(j))
            dgreenirj(iii, equivforce(j)) = dgreenirj(iii, equivforce(j)) - &
                mu0*2.*tatmp(i)*tatmp(j)* (x1 - x2)/dr
            x1 = green_function(rcetmp(i), zcetmp(i) + dz/2., rcetmp(j), zcetmp(j))
            x2 = green_function(rcetmp(i), zcetmp(i) - dz/2., rcetmp(j), zcetmp(j))
            dgreenizj(iii, equivforce(j)) = dgreenizj(iii, equivforce(j)) - &
                mu0*2.*tatmp(i)*tatmp(j)* (x1 - x2)/dz
        endif
    enddo
enddo

write(*, *) 'grid forces from plasma to coil'
do i=1, ielem
    iii = equivforce(i)
    do jj=1, nz2
        do ii=1, nr2
            x1 = green_function(rcetmp(i) + dr/2., zcetmp(i), r(ii), z(jj))
            x2 = green_function(rcetmp(i) - dr/2., zcetmp(i), r(ii), z(jj))
            dgreenirpl(ii, jj, iii) = dgreenirpl(ii, jj, iii) - mu0*2.*tatmp(i)* (x1 - x2)/dr
            x1 = green_function(rcetmp(i), zcetmp(i) + dz/2., r(ii), z(jj))
            x2 = green_function(rcetmp(i), zcetmp(i) - dz/2., r(ii), z(jj))
            dgreenizpl(ii, jj, iii) = dgreenizpl(ii, jj, iii) - mu0*2.*tatmp(i)* (x1 - x2)/dz
        enddo
    enddo
enddo


!generate zlimpotential
zlimpotential = 1
do j=1, nz2
    do i=1, nr2
        if (i < ilim_minR) zlimpotential(i, j) = 0
        if (i > ilim_maxR) zlimpotential(i, j) = 0
        if (j < ilim_minZ) zlimpotential(i, j) = 0
        if (j > ilim_maxZ) zlimpotential(i, j) = 0
    enddo
enddo

open(32,file='limatrix.txt')
do i=1, nconduc
    write(32, *) indconduc(i, 1:nconduc)
enddo
close(32)

! write everything on file

fname = trim(data_dir2)//'/machine_description_out.'//trim(machine)
open(32, file=trim(fname))

write(32, *) nr2, nz2
write(32, *) rmin
write(32, *) rmax
write(32, *) zmin
write(32, *) zmax
write(32, *) alpsep

write(32, *) nactive,npassive

write(32, *) ncoils
do i=1, ncoils
    write(32, *) rcoil(i), zcoil(i), drcoil(i), dzcoil(i), anglehcoil(i), anglecoil(i), mequivalence(i)
enddo

write(32, *) nlimiter
do i=1, nlimiter
    write(32, *) limiterr(i),limiterz(i)
enddo
write(32,*) lim_maxR
write(32,*) lim_minR
write(32,*) lim_maxZ
write(32,*) lim_minZ

do i=nactive+1,npassive
    write(32, *) r_cond(i),z_cond(i)
enddo


write(32, *) nconduc
do i=1, nconduc
    write(32, *) indconduc(i, 1:nconduc)
enddo

write(32, *) nactive, nconduc
do i=1, nactive
    write(32, *) (resconduc(i, j), j=1, nactive)
enddo
do i=nactive+1, nconduc
    write(32, *) resconduc(i, i)
enddo

do i=1, nconduc
    do j=1, nr2
        write(32, *) greeni(j, 1:nz2, i)
    enddo
enddo

!force matrix
write(32, *) nblocks
do j=1, nblocks
    do i=1, nblocks
        write(32, *) dgreenirj(i, j), dgreenizj(i, j)
    enddo
enddo
do ii=1, nblocks
    do j=1, nz2
        do i=1, nr2
            write(32, *) dgreenirpl(i, j, ii), dgreenizpl(i, j, ii)
        enddo
    enddo
enddo

do jj=1, nz2
    do ii=1, nr2
        write(32, *) zlimpotential(ii, jj)
    enddo
enddo

write(*, *) 'bound'
write(32, *) nint((2.*nr + 2.*nz)*(2.*nr + 2.*nz))
! lower side
do i=2, nr1
    do j=2, nr1
        greenf = green_function_includingsamepoint(r(i), z(1), r(j), z(1), dr, r(i))
        write(32, *) greenf
    enddo
!right side
    do j=2, nz1
        greenf = green_function_includingsamepoint(r(i), z(1), r(nr2), z(j), dz, r(i))
        write(32, *) greenf
    enddo
! upper side
    do j=2, nr1
        greenf = green_function_includingsamepoint(r(i), z(1), r(j), z(nz2), dr, r(i))
        write(32, *) greenf
    enddo
!left side
    do j=2, nz1
        greenf = green_function_includingsamepoint(r(i), z(1), r(1), z(j), dz, r(i))
        write(32, *) greenf
    enddo
enddo

! right side
do i=2, nz1
    do j=2, nr1
        greenf = green_function_includingsamepoint(r(nr2), z(i), r(j), z(1), dr, r(nr2))
        write(32, *) greenf
    enddo
!right side
    do j=2, nz1
        greenf = green_function_includingsamepoint(r(nr2), z(i), r(nr2), z(j), dz, r(nr2))
        write(32, *) greenf
    enddo
! upper side
    do j=2, nr1
        greenf = green_function_includingsamepoint(r(nr2), z(i), r(j), z(nz2), dr, r(nr2))
        write(32, *) greenf
    enddo
!left side
    do j=2, nz1
        greenf = green_function_includingsamepoint(r(nr2), z(i), r(1), z(j), dz, r(nr2))
        write(32, *) greenf
    enddo
enddo

! upper side
do i=2, nr1
    do j=2, nr1
        greenf = green_function_includingsamepoint(r(i), z(nz2), r(j), z(1), dr, r(i))
        write(32, *) greenf
    enddo
!right side
    do j=2, nz1
        greenf = green_function_includingsamepoint(r(i), z(nz2), r(nr2), z(j), dz, r(i))
        write(32, *) greenf
    enddo
! upper side
    do j=2, nr1
        greenf = green_function_includingsamepoint(r(i), z(nz2), r(j), z(nz2), dr, r(i))
        write(32, *) greenf
    enddo
!left side
    do j=2, nz1
        greenf = green_function_includingsamepoint(r(i), z(nz2), r(1), z(j), dz, r(i))
        write(32, *) greenf
    enddo
enddo

! left side
do i=2, nz1
    do j=2, nr1
        greenf = green_function_includingsamepoint(r(1), z(i), r(j), z(1), dr, r(1))
        write(32, *) greenf
    enddo
!right side
    do j=2, nz1
        greenf = green_function_includingsamepoint(r(1), z(i), r(nr2), z(j), dz, r(1))
        write(32, *) greenf
    enddo
! upper side
    do j=2, nr1
        greenf = green_function_includingsamepoint(r(1), z(i), r(j), z(nz2), dr, r(1))
        write(32, *) greenf
    enddo
!left side
    do j=2, nz1
        greenf = green_function_includingsamepoint(r(1), z(i), r(1), z(j), dz, r(1))
        write(32, *) greenf
    enddo
enddo

close(32)

write(*, *) 'Stored green functions'

return
end subroutine generate_files_feqis
