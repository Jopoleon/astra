subroutine extpol

use comrec, only: x, y, ue
use compol_add, only: psie
use compol, only: nr, nt, r, z

implicit none

integer :: i, j, ic, jc
real*8 :: r0, r1, r2, z0, z1, z2, u1, u2, u3, u4, ddx, ddy
real*8, external :: blin_

ddx = x(2) - x(1)
ddy = y(2) - y(1)

do i=1, nr
    do j=1, nt
        r0 = r(i, j)
        z0 = z(i, j)
        ic = (r0 - x(1))/ddx + 1
        jc = (z0 - y(1))/ddy + 1

        r1 = x(ic)
        r2 = x(ic+1)
        z1 = y(jc)
        z2 = y(jc+1)

        u1 = ue(ic, jc)
        u2 = ue(ic+1, jc)
        u3 = ue(ic+1, jc+1)
        u4 = ue(ic, jc+1)

        psie(i, j) = blin_(r0, z0, r1, r2, z1, z2, u1, u2, u3, u4)
    enddo
enddo

return
end subroutine extpol

!----------------------------------------------------------------
subroutine f_psiful

use compol_add, only: psii, psie
use compol, only: nr, nt, psi

implicit none

integer :: i, j

do i=1, nr
    do j=1, nt
        psi(i, j) = psii(i, j) + psie(i, j)
    enddo
enddo

return
end subroutine f_psiful

!----------------------------------------------------------------
subroutine f_getdrdz(rl_1, zl_1, rr_1, zr_1, drdzl, drdzr)

use compol_add, only: psii, psie
use compol, only: nr, nt, nt1, iplas, r, z, psi

implicit none

integer :: i, i1, j
double precision :: drdzl, drdzr, rl_1, zl_1, rr_1, zr_1, &
    rl_2, zl_2, rr_2, zr_2, dum1, dum2(nt), dum3, dum4, dum5, dum6

dum1 = psi(iplas, nt1)
dum2 = psi(nr, 1:nt)

i = minloc(abs(dum2 - dum1), 1)
i1 = i + 1
if (i == nt) i1 = 1
dum3 = r(nr, i) + (dum1 - psi(nr, i))/ &
    (psi(nr, i) - psi(nr, i1))*(r(nr, i) - r(nr, i1))
dum4 = z(nr, i) + (dum1 - psi(nr, i))/ &
    (psi(nr, i) - psi(nr, i1))*(z(nr, i) - z(nr, i1))

dum2(i) = 1.e6
i = minloc(abs(dum2 - dum1), 1)
i1 = i + 1
if (i == nt) i1 = 1
dum5 = r(nr, i) + (dum1 - psi(nr, i))/ &
    (psi(nr, i) - psi(nr, i1))*(r(nr, i) - r(nr, i1))
dum6 = z(nr, i) + (dum1 - psi(nr, i))/ &
    (psi(nr, i) - psi(nr, i1))*(z(nr, i) - z(nr, i1))

if (dum5 > dum3) then
    rr_1 = dum5
    zr_1 = dum6
    rl_1 = dum3
    zl_1 = dum4
else
    rl_1 = dum5
    zl_1 = dum6
    rr_1 = dum3
    zr_1 = dum4
endif

dum2 = psi(nr-1, 1:nt)

i = minloc(abs(dum2 - dum1), 1)
i1 = i + 1
if (i == nt) i1 = 1
dum3 = r(nr-1, i) + (dum1 - psi(nr-1, i))/ &
    (psi(nr-1, i) - psi(nr-1, i1))*(r(nr-1, i) - r(nr-1, i1))
dum4 = z(nr-1, i) + (dum1 - psi(nr-1, i))/ &
    (psi(nr-1, i) - psi(nr-1, i1))*(z(nr-1, i) - z(nr-1, i1))

dum2(i) = 1.e6
i = minloc(abs(dum2-dum1), 1)
i1 = i + 1
if (i == nt) i1 = 1
dum5 = r(nr-1, i) + (dum1 - psi(nr-1, i))/ &
    (psi(nr-1, i) - psi(nr-1, i1))*(r(nr-1, i) - r(nr-1, i1))
dum6 = z(nr-1, i) + (dum1 - psi(nr-1, i))/ &
    (psi(nr-1, i) - psi(nr-1, i1))*(z(nr-1, i) - z(nr-1, i1))

if (dum5 > dum3) then
    rr_2 = dum5
    zr_2 = dum6
    rl_2 = dum3
    zl_2 = dum4
else
    rl_2 = dum5
    zl_2 = dum6
    rr_2 = dum3
    zr_2 = dum4
endif

drdzl = (zl_2 - zl_1)/(rl_2 - rl_1)
drdzr = (zr_2 - zr_1)/(rr_2 - rr_1)

do i=1, nr
    do j=1, nt
        psi(i, j) = psii(i, j) + psie(i, j)
    enddo
enddo

return
end subroutine f_getdrdz

!----------------------------------------------------------------
subroutine artfil

use sp_parameters, only: ntp
use keys, only: kpr
use compol_add, only: clr, clz
use compol, only: nr, nt, nt1, r, z, psi, rm

implicit none

integer, parameter :: nshp=ntp+1
integer :: i, j, nsh
real*8 :: clrn, clzn, sigma
real*8, dimension(nshp) :: xs, ys, fun
real*8, dimension(5) :: dpm

nsh = 1

xs(nsh) = r(1, 1)
ys(nsh) = z(1, 1)
fun(nsh) = psi(1, 1)
if (kpr == 1) then
    write(*, *) 'artfil:psi(1, 1)', fun(nsh)
endif

do j=2, nt1
    nsh = nsh + 1
    xs(nsh) = r(2, j)
    ys(nsh) = z(2, j)
    fun(nsh) = psi(2, j)
enddo

call deriv5(xs, ys, fun, nsh, 5, dpm)

clrn = -dpm(1)*0.5d0/rm
clzn = -dpm(2)

sigma = 1.d0

clr = sigma*clrn + (1.d0 - sigma)*clr
clz = sigma*clzn + (1.d0 - sigma)*clz

do i=1, nr
    do j=1, nt
        psi(i, j) = psi(i, j) + clz*z(i, j) + clr*r(i, j)*r(i, j)
    enddo
enddo

return
end subroutine artfil

!----------------------------------------------------------------
subroutine f_xpoint(rx, zx, ix, jx, psx, tet0, kodex)

use sp_parameters, only: twopi
use compol, only: nr, nr1, nt, nt1, teta, r, z, psi, rm, zm

implicit none

integer, parameter :: nshp=20

integer, intent(inout) :: ix, jx
integer, intent(out) :: kodex
real*8, intent(inout) :: rx, zx
real*8, intent(out) :: psx, tet0

integer :: i, j, k, l, jc, ixn, jxn, numbit, nsh
real*8 :: det, dro, dromn, dropl, rrx, zzx, xxyy, tetp, &
    rxn, zxn, dr0, dz0, ro0
real*8, dimension(10) :: dp
real*8, dimension(nshp) :: xs, ys, fun

numbit = 0
kodex = 0

 100 continue

numbit = numbit + 1
rxn = rx
zxn = zx
ixn = ix
jxn = jx
dr0 = rx - rm
dz0 = zx - zm
ro0 = sqrt(dr0**2 + dz0**2)

tetp = ACOS(dr0/ro0)
if (dz0 < 0.d0) then
    tet0 = twopi - tetp
else
    tet0 = tetp
endif

if (tet0 < teta(1) ) tet0 = tet0 + twopi
if (tet0 > teta(nt)) tet0 = tet0 - twopi

do j=1, nt1
   if (tet0 >= teta(j) .AND. tet0 < teta(j+1)) jc = j
enddo

dro = 1.d9

do i=3, nr1
    dromn = sqrt((r(i+1, jc)   - rx)**2 + (z(i+1, jc  ) - zx)**2)
    dropl = sqrt((r(i+1, jc+1) - rx)**2 + (z(i+1, jc+1) - zx)**2)
    if (dropl <= dro .AND. dropl <= dromn) then
        dro = dropl
        ix = i + 1
        jx = jc + 1
    elseif (dromn <= dro .AND. dromn <= dropl) then
        dro = dromn
        ix = i + 1
        jx = jc
    endif
enddo

if (jx == nt) jx = 2   
if (jx == 1) jx = nt1   

if (ix == nr) then   !xpoint out of box
   kodex = 1
   return
endif

nsh = 1
xs(nsh) = r(ix, jx)
ys(nsh) = z(ix, jx)
fun(nsh) = psi(ix, jx)

do k=-1, 1
    i = ix+k
    do l=-1, 1
        j = jx + l
        if (k == 0 .AND. l == 0 ) CYCLE
        if (i > nr ) CYCLE
        nsh = nsh + 1
        xs(nsh) = r(i, j)
        ys(nsh) = z(i, j)
        fun(nsh) = psi(i, j)
    enddo
enddo

do k=-2, 2, 2
    i = ix + k
    j = jx
    if (k == 0) CYCLE
    if (i > nr ) CYCLE
    nsh = nsh + 1
    xs(nsh) = r(i, j)
    ys(nsh) = z(i, j)
    fun(nsh) = psi(i, j)
enddo

call deriv5(xs, ys, fun, nsh, 5, dp)

DET = dp(3)*dp(5) - dp(4)**2

Rx = Xs(1) + ( dp(2)*dp(4) - dp(1)*dp(5) )/DET
Zx = Ys(1) + ( dp(1)*dp(4) - dp(2)*dp(3) )/DET
rrx = xs(1)
zzx = ys(1)

psx = fun(1)+ dp(1)*(rx - rrx) + dp(2)*(zx - zzx) + &
        0.5d0*dp(3)*(rx - rrx)*(rx - rrx) + &
              dp(4)*(rx - rrx)*(zx - zzx) + &
        0.5d0*dp(5)*(zx - zzx)*(zx - zzx)

xxyy = dp(3)*dp(5)

if ( (ix /= ixn .OR. jx /= jxn) .AND. (numbit < 10) ) then
    goto 100
else
    xxyy = dp(3)*dp(5)
endif

return
end subroutine f_xpoint
