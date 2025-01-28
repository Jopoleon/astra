subroutine psi_fil(psicon, ncequi)

use comrec, only: zaindk
use comblc, only: nj1, ni1, curf, dri, dzj

implicit none

integer, intent(in) :: ncequi
real*8, intent(out), dimension(*) :: psicon

integer :: iq, i, j
real*8 :: vrftfa, zpscon

do iq=1, ncequi
    zpscon = 0.d0
    do j=2, nj1
        do i=2, ni1
            vrftfa = zaindk(i, j, iq)
            zpscon = zpscon + vrftfa*curf(i, j)*dri(i)*dzj(j)
        enddo
    enddo
    psicon(iq) = zpscon
enddo

return
end subroutine psi_fil

!----------------------------------------------------------------
subroutine extfil(pcequi, ncequi)

use sp_parameters, only: nip, njp, amu0
use iopath, only: path
use comevl, only: nloc, nequi, npfc
use comblc, only: nj, ni, ue, ue_k, rcondzzz

implicit none

integer, intent(in) :: ncequi
real*8, intent(in), dimension(*) :: pcequi

integer :: i, j, iq, ncpfc, nkread
real*8, dimension(nip, njp) :: aindk
character(len=80) :: fname

ncpfc = nloc(npfc)

do j=1, nj
    do i=1, ni
        ue(i, j) = 0.d0
    enddo
enddo

ue_k = 0.

write(fname,'(a,a)') TRIM(path), '/exf.wr'
open(1, file=fname, form='formatted')
read(1,*) nkread, nequi
do iq=1,ncequi
    read(1,*) ((aindk(i, j), j=1, nj), i=1, ni)
    do j=1, nj
        do i=1, ni
            ue(i, j) = ue(i, j) + aindk(i, j)*PCEQui(iq)*amu0
            ue_k(i, j, iq) = ue_k(i, j, iq) + aindk(i,j)*PCEQui(iq)*amu0
        enddo
    enddo
    rcondzzz = iq
enddo
close(1)

return
end subroutine extfil

!----------------------------------------------------------------
subroutine exfmat(rk, zk, nk, rk1, zk1, rk2, zk2, ntipe, NECON, WECON )

use sp_parameters, only: nip, njp, pi
use iopath, only: path
use comevl, only: nequi, npfc, nloc
use comblc, only: ni, nj, r, z

implicit none

integer, intent(in) :: nk
integer, intent(in), dimension(*) :: ntipe, NECON
real*8, intent(in), dimension(*) :: rk, zk, WECON, rk1, zk1, rk2, zk2

integer :: i, j, iq, ik, ncpfc, nves, ncequi, ibeg, ntip, ncam
real*8 :: ddx, zaindk, r1, r2, z1, z2, fint, dlong
real*8, dimension(nip, njp) :: aindk
real*8, external :: greeni
character(len=80) :: fname

ncpfc  = nloc(npfc)
nves   = nk - ncpfc
ncequi = nequi + nves

write(fname,'(a, a)') TRIM(path), '/exf.wr'
open(1, file=fname, form='formatted')
write(1, *) ncequi, nequi
do iq=1, nequi
    do i=1, ni
        do j=1, nj
            aindk(i, j) = 0.d0
        enddo
    enddo
    do i=1, ni
        do j=1, nj
            do ik=1,ncpfc
                if (necon(ik) == iq)  then
                    ddx = sqrt( (r(i) - rk(ik))**2 + (z(j) - zk(ik))**2 )
                    zaindk = greeni(r(i), z(j), rk(ik), zk(ik))/pi
                    aindk(i, j) = aindk(i, j) + zaindk*wecon(ik)
                endif
            enddo
        enddo
    enddo
    write(1, *) ((aindk(i, j), j=1, nj), i=1, ni)
enddo

ibeg = ncpfc + 1

if (ibeg > nk) goto 11

ncam=0.

do ik=ibeg,nk
    ntip=ntipe(ik)
    if (ntip == 1) then
        r1 = rk1(ik)
        z1 = zk1(ik)
        r2 = rk2(ik)
        z2 = zk2(ik)
        dlong = sqrt( (r2 - r1)**2 + (z2 - z1)**2 )
    endif

    do i=1, ni
        do j=1, nj
            ddx = sqrt( (r(i) - rk(ik))**2 + (z(j) - zk(ik))**2 )
            if (ntip == 1) then
                call bint(r(i), z(j), r1, z1, r2, z2, fint, 1)
                aindk(i, j) = fint/dlong
            else
                aindk(i,j)=greeni(r(i),z(j),rk(ik),zk(ik))/pi
            endif
         enddo
      enddo

      write(1,*) ((aindk(i, j), j=1, nj), i=1, ni)
      ncam = ncam + 1
enddo

 11 continue

write(1, *) ncam
close(1)

return
end subroutine exfmat

!----------------------------------------------------------------
subroutine cfr_mat(rk, zk, WECON)

use sp_parameters, only: nip, njp, npfc0, pi
use iopath, only: path
use comevl, only: npfc, nloc
use comblc, only: ni, nj, r, z

implicit none

real*8, intent(in), dimension(*) :: rk, zk, WECON

integer :: i, j, k, m, ik, ib, ie, jb, je, kb, ke
real*8 :: SumFRij, SumFZij, SumFRPijk, SumFZPijk, rkk, zkk, rkm, zkm, dGdr, dGdz
real*8, dimension(npfc0, npfc0) :: fr_indk, fz_indk
real*8, dimension(nip, njp) :: fpr_indk, fpz_indk
character(len=80) :: fname

! coil coil interaction

do i=1, npfc
    do j=1, npfc
        SumFRij = 0.d0
        SumFZij = 0.d0
        if (i /= j) then
            ib = Nloc(i-1) + 1
            ie = Nloc(i)
            jb = Nloc(j-1) + 1
            je = Nloc(j)
            do k=ib, ie
                do m=jb, je
                    rkk = rk(k) !1.d0
                    zkk = zk(k) !1.d-2
                    rkm = rk(m) !1.d0
                    zkm = zk(m) !-1.d-2
                    call grGREN(rkk, zkk, rkm, zkm, dGdr, dGdz)
                    SumFRij = SumFRij + dGdr*wecon(k)*wecon(m)/rk(m)/pi
                    SumFZij = SumFZij - dGdz*wecon(k)*wecon(m)*2.d0
                enddo
            enddo
            fr_indk(i, j) = SumFRij
            fz_indk(i, j) = SumFZij
        endif
    enddo
enddo

write(fname, '(a, a)') TRIM(path), '/forcmat.wr'
open(1, file=fname, form='formatted')
write(1,*) npfc
write(1,*) ((fr_indk(i,j),j=1,npfc),i=1,npfc)
write(1,*) ((fz_indk(i,j),j=1,npfc),i=1,npfc)

! plasma coil interaction
do k=1, npfc
    kb = Nloc(k-1) + 1
    ke = Nloc(k)
    do i=1, ni
        do j=1, nj
            SumFRPijk = 0.d0
            SumFZPijk = 0.d0
            do ik=kb, ke
                call grGREN(r(i), z(j), rk(ik), zk(ik), dGdr, dGdz)
                SumFRPijk = SumFRPijk + dGdr*wecon(ik)/rk(ik)/pi
                SumFZPijk = SumFZPijk - dGdz*wecon(ik)*2.d0
            enddo
            fpr_indk(i, j)= SumFRPijk
            fpz_indk(i, j)= SumFZPijk
        enddo
    enddo
    write(1, *) ((fpr_indk(i, j), j=1, nj), i=1, ni)
    write(1, *) ((fpz_indk(i, j), j=1, nj), i=1, ni)
enddo
close(1)

return
end subroutine cfr_mat

!----------------------------------------------------------------
subroutine ext_fil(pcequi, ncequi)

use sp_parameters, only: pi, amu0
use comrec, only: zaindk
use comblc, only: ni, nj, ue_k, ue, rcondzzz

implicit none

integer, intent(in) :: ncequi
real*8, intent(in), dimension(ncequi) :: pcequi

integer :: iq

iq = 1
ue_k(1:ni, 1:nj, iq) = zaindk(1:ni, 1:nj, iq)*PCEQui(iq)*amu0
ue(1:ni, 1:nj) = ue_k(1:ni, 1:nj, iq)
do iq=2, ncequi
    ue_k(1:ni, 1:nj, iq) = zaindk(1:ni, 1:nj, iq)*PCEQui(iq)*amu0
    ue(1:ni, 1:nj) = ue(1:ni, 1:nj) + ue_k(1:ni, 1:nj, iq)
enddo
rcondzzz = iq

return
end subroutine ext_fil

!----------------------------------------------------------------
subroutine rdexf(ncequi)

use sp_parameters, only: nip, njp
use iopath, only: path
use comevl, only: nequi, nloc, npfc
use comrec, only: zaindk
use comblc, only: ni, nj

implicit none

integer, intent(in) :: ncequi

integer :: i, j, ncpfc, nk, iq
real*8, dimension(nip, njp) :: aindk
character(len=80) :: fname

write(fname, '(a, a)') TRIM(path), '/exf.wr'
open(1, file=fname, form='formatted')
read(1, *) nk, nequi
ncpfc = nloc(npfc)
do iq=1, ncequi
    read(1, *) ((aindk(i, j), j=1, nj), i=1, ni)
    do j=1, nj
        do i=1, ni
            zaindk(i, j, iq)=aindk(i, j)
        enddo
    enddo
enddo
close(1)

return
end subroutine rdexf

!----------------------------------------------------------------
subroutine bndmat

use sp_parameters, only: nbndp
use comblc, only: ni, ni1, nj, nj1, nbnd, r, z, binadg

implicit none

integer :: ib, i, j, ibc
real*8 :: r0, r1, z0, z1, rr, zz, fint
real*8, dimension(nbndp) :: xk, yk

ib = 1

xk(ib) = r(1)
yk(ib) = z(1)

do i=2, ni1
    ib = ib + 1
    xk(ib) = r(i)
    yk(ib) = z(1)
enddo

ib = ib + 1
xk(ib) = r(ni)
yk(ib) = z(1)

do j=2, nj1
    ib = ib + 1
    xk(ib) = r(ni)
    yk(ib) = z(j)
enddo

ib = ib + 1
xk(ib) = r(ni)
yk(ib) = z(nj)

do i=ni1, 2, -1
    ib = ib + 1
    xk(ib) = r(i)
    yk(ib) = z(nj)
enddo

ib = ib + 1
xk(ib) = r(1)
yk(ib) = z(nj)

do j=nj1, 2, -1
    ib = ib + 1
    xk(ib) = r(1)
    yk(ib) = z(j)
enddo

ib = ib + 1
xk(ib) = r(1)
yk(ib) = z(1)

do ib = 1, nbnd
    rr = xk(ib)
    zz = yk(ib)
    do ibc = 1, nbnd
        r0 = xk(ibc)
        z0 = yk(ibc)
        r1 = xk(ibc+1)
        z1 = yk(ibc+1)
        call bint(rr, zz, R0, Z0, r1, z1, Fint, 1)
        binadg(ibc, ib) = fint
    enddo
enddo

return
end subroutine bndmat

!----------------------------------------------------------------
real*8 function blin(i, j, r0, z0)

use comblc, only: r, z, ui

implicit none

integer, intent(in) :: i, j
real*8, intent(in) :: r0, z0

real*8 :: r1, r2, z1, z2, u1, u2, u3, u4, s1, s2, s3, s4

z1 = z(j)
z2 = z(j+1)
r1 = r(i)
r2 = r(i+1)

u1 = ui(i, j)
u2 = ui(i+1, j)
u3 = ui(i+1, j+1)
u4 = ui(i, j+1)

s1 = (r2-r0)*(z2-z0)
s3 = (r0-r1)*(z0-z1)
s2 = (r0-r1)*(z2-z0)
s4 = (r2-r0)*(z0-z1)

blin = (s1*u1+s2*u2+s3*u3+s4*u4)/(s1+s2+s3+s4)

return
end function blin

!----------------------------------------------------------------
real*8 function bline(i, j, r0, z0)

use comblc, only: r, z, ue

implicit none

integer, intent(in) :: i, j
real*8, intent(in) :: r0, z0

real*8 :: r1, r2, z1, z2, u1, u2, u3, u4, s1, s2, s3, s4

z1 = z(j)
z2 = z(j+1)
r1 = r(i)
r2 = r(i+1)

u1 = ue(i, j)
u2 = ue(i+1, j)
u3 = ue(i+1, j+1)
u4 = ue(i, j+1)

s1 = (r2 - r0)*(z2 - z0)
s3 = (r0 - r1)*(z0 - z1)
s2 = (r0 - r1)*(z2 - z0)
s4 = (r2 - r0)*(z0 - z1)

bline = (s1*u1 + s2*u2 + s3*u3 + s4*u4)/(s1 + s2 + s3 + s4)

return
end function bline

!----------------------------------------------------------------
real*8 function blint(i, j, r0, z0)

use comblc, only: r, z, u

implicit none

integer, intent(in) :: i, j
real*8, intent(in) :: r0, z0

real*8 :: r1, r2, z1, z2, u1, u2, u3, u4, s1, s2, s3, s4

z1 = z(j)
z2 = z(j+1)
r1 = r(i)
r2 = r(i+1)

u1 = u(i, j)
u2 = u(i+1, j)
u3 = u(i+1, j+1)
u4 = u(i, j+1)

s1 = (r2 - r0)*(z2 - z0)
s3 = (r0 - r1)*(z0 - z1)
s2 = (r0 - r1)*(z2 - z0)
s4 = (r2 - r0)*(z0 - z1)

blint = (s1*u1 + s2*u2 + s3*u3 + s4*u4)/(s1 + s2 + s3 + s4)

return
end function blint

!----------------------------------------------------------------
subroutine grid_spid

use comblc, only: ni, ni1, nj, nj1, rmin, rmax, zmin, zmax, r, z

implicit none

integer :: i, j
real*8 :: ddr, ddz

ddr = (rmax - rmin)/FLOAT(ni1)
ddz = (zmax - zmin)/FLOAT(nj1)

r(1) = rmin
z(1) = zmin

do i=2, ni
    r(i) = r(i-1) + ddr
enddo

do j=2, nj
    z(j) = z(j-1) + ddz
enddo

return
end subroutine grid_spid

!----------------------------------------------------------------
subroutine geom

use comblc, only: ni1, nj1, nblm, iblm, jblm, r, z, dr, dz, &
    dri, dzj, r12, rblm, zblm

implicit none

integer :: i, j, ilm, ic, jc

do i=1, ni1
    dr(i) = r(i+1) - r(i)
    r12(i) = (r(i+1) + r(i))*0.5d0
    if (i == 1) CYCLE
    dri(i) = (r(i+1) - r(i-1))*0.5d0
enddo

do j=1, nj1
    dz(j) = z(j+1) - z(j)
    if (j == 1) CYCLE
    dzj(j) = (z(j+1) - z(j-1))*0.5d0
enddo

do ilm=1, nblm
    do i=1, ni1
        if (rblm(ilm) < r(i+1) .AND. rblm(ilm) >= r(i) ) then
            ic = i
            EXIT
        endif
    enddo
    iblm(ilm) = ic
    do j = 1, nj1
        if (zblm(ilm) < z(j+1) .AND. zblm(ilm) >= z(j) ) then
            jc = j
            EXIT
        endif
    enddo
    jblm(ilm) = jc
enddo

return
end subroutine geom

!----------------------------------------------------------------
subroutine zero_ax(isol)

use sp_parameters, only: neqp, amu0
use comblc, only: ni, ni1, nj, nj1, imax, jmax, r, z, rm0, zm0, &
    right, curf, g, dri, dzj, ux0, ui

implicit none

integer, intent(in) :: isol

integer :: i, j, ic, jc, il, nlin
real*8 :: r0, z0, tokint, curden

do i=1, ni1
    if (rm0 < r(i+1) .AND. rm0 >= r(i) ) then
        ic = i
        EXIT
    endif
enddo
imax = ic

do j=1, nj1
    if (zm0 < z(j+1) .AND. zm0 >= z(j) ) then
        jc = j
        EXIT
    endif
enddo
jmax = jc

do il=1, neqp
    right(il) = 0.d0
enddo

do i=1, ni
    do j=1, nj
        curf(i, j) = 0.d0
        g(i, j) = 0.d0
    enddo
enddo

tokint = 0.d0
do i=2, ni1
    r0 = r(i)
    do j=2, nj1
        z0 = z(j)
        call cur_map(curden, r0, z0)
        il = nlin(i, j)
        right(il) = -curden*dri(i)*dzj(j)
        curf(i, j) = curden
        tokint = tokint + curden*dri(i)*dzj(j)
    enddo
enddo
tokint = tokint/amu0
ux0 = -10d9

call solve(isol, g)
call bound
call solve(isol, ui)

return
end subroutine zero_ax

!----------------------------------------------------------------
subroutine cur_map(curden, rk, zk)

use sp_parameters, only: twopi
use compol, only: rm, zm, teta, nr, nt, nt1, ro, iplas, cur

implicit none

real*8, intent(in) :: rk, zk
real*8, intent(out) :: curden

integer :: i, j, ic, jc
real*8 :: dr0, dz0, ro0, ro1, ro2, ro3, ro4, ro12, ro34, &
    rob1, rob2, rob12, u1, u2, u3, u4, tetp, tet0, tet1, tet2, blintr

curden = 0.d0
dr0 = rk - rm
dz0 = zk - zm

ro0 = sqrt(dr0**2 + dz0**2)
tetp = ACOS(dr0/ro0)
if (dz0 < 0.d0) then
    tet0 = twopi - tetp
else
    tet0 = tetp
endif

if (tet0 < teta(2)) tet0 = tet0 + twopi
if (tet0 >= teta(nt)) tet0 = tet0 - twopi
jc = 0
do j=2, nt1
    if (tet0 >= teta(j) .AND. tet0 < teta(j+1)) jc = j
enddo

if (jc == 0) write(6, *) 'curmap: jc is not defined'

rob1 = ro(nr, jc)
rob2 = ro(nr, jc+1)
tet1 = teta(jc)
tet2 = teta(jc+1)
rob12 = (rob1*(tet2-tet0)+rob2*(tet0-tet1))/(tet2-tet1)

if (ro0 < rob12) then
    do i=1, iplas-1
        ro1 = ro(i, jc)
        ro2 = ro(i, jc+1)
        ro3 = ro(i+1, jc)
        ro4 = ro(i+1, jc+1)
        ro12 = (ro1*(tet2-tet0)+ro2*(tet0-tet1))/(tet2-tet1)
        ro34 = (ro3*(tet2-tet0)+ro4*(tet0-tet1))/(tet2-tet1)
        if (ro0 > ro12 .AND. ro0 <= ro34) then
            ic = i
            EXIT
        endif
    enddo

    ro1 = ro(ic, jc)
    ro2 = ro(ic+1, jc)
    ro3 = ro(ic+1, jc+1)
    ro4 = ro(ic, jc+1)

    u1 = cur(ic, jc)
    u2 = cur(ic+1, jc)
    u3 = cur(ic+1, jc+1)
    u4 = cur(ic, jc+1)

    curden = blintr(tet0, ro0, tet1, tet2, ro1, ro2, ro3, ro4, u1, u2, u3, u4)
endif

return
end subroutine cur_map

!----------------------------------------------------------------
real*8 function blintr(tet0, ro0, tet1, tet2, ro1, ro2, ro3, ro4, u1, u2, u3, u4)

implicit none

real*8, intent(in) :: tet0, ro0, tet1, tet2, ro1, ro2, ro3, ro4, u1, u2, u3, u4
real*8 :: ro14, ro23, u14, u23

ro14 = (ro1*(tet2 - tet0) + ro4*(tet0 - tet1))/(tet2 - tet1)
ro23 = (ro2*(tet2 - tet0) + ro3*(tet0 - tet1))/(tet2 - tet1)

u14 = (u1*(tet2 - tet0) + u4*(tet0 - tet1))/(tet2 - tet1)
u23 = (u2*(tet2 - tet0) + u3*(tet0 - tet1))/(tet2 - tet1)

blintr = (u14*(ro23 - ro0) + u23*(ro0 - ro14))/(ro23 - ro14)

return
end function blintr
