subroutine flux_r(psitok, ncequi)

use comrec, only: x, y, zaindk
use compol, only: nt1, r, z, rm, zm, cur, sq1, sq2, sq3, sq4, iplas

implicit none

integer, intent(in) :: ncequi
real*8, intent(out), dimension(*) :: psitok

integer :: i, j, k, ic, jc
real*8 :: ddx, ddy, sqcen, sqk, r0, z0, u1, u2, u3, u4, r1, r2, z1, z2
real*8, external :: blin_

ddx = x(2) - x(1)
ddy = y(2) - y(1)

sqcen = 0.d0
do j=2, nt1
    sqcen = sqcen + sq1(1, j) +sq4(1, j)
enddo

r0 = rm
z0 = zm
if (r0 > 100) then
    write(*, *) 'mag axis major radius is > 100 m'
    call err_catch_a
endif
if (abs(z0) > 100)  then
    write(*, *) 'mag axis Z is <> 100 m'
    call err_catch_a
endif
if (isnan(r0))  then
    write(*, *) 'mag axis major radius is NaN'
    call err_catch_a
endif
if (isnan(z0))  then
    write(*, *) 'mag axis Z is NaN'
    call err_catch_a
endif
ic = ((r0 - x(1))/ddx) + 1
jc = ((z0 - y(1))/ddy) + 1

r1 = x(ic)
r2 = x(ic+1)
z1 = y(jc)
z2 = y(jc+1)

do k=1, ncequi
    u1 = zaindk(ic, jc, k)
    u2 = zaindk(ic+1, jc, k)
    u3 = zaindk(ic+1, jc+1, k)
    u4 = zaindk(ic, jc+1, k)
    psitok(k) = blin_(r0, z0, r1, r2, z1, z2, u1, u2, u3, u4)*cur(1, 2)*sqcen
enddo

do i=2, iplas
    do j=2, nt1
        if (i /= iplas) then
            sqk = sq1(i, j) + sq2(i-1, j) + sq3(i-1, j-1) + sq4(i, j-1)
        else
            sqk = sq2(i-1, j) + sq3(i-1, j-1)
        endif
        r0 = r(i, j)
        z0 = z(i, j)
        if (r0 > 100) then
            write(*,*) 'mag axis major radius is > 100 m'
            call err_catch_a
        endif
        if (abs(z0) > 100)  then
            write(*,*) 'mag axis Z is <> 100 m'
            call err_catch_a
        endif
        if (isnan(r0))  then
            write(*,*) 'mag axis major radius is NaN'
            call err_catch_a
        endif
        if (isnan(z0))  then
            write(*,*) 'mag axis Z is NaN'
            call err_catch_a
        endif
        ic = ((r0 - x(1))/ddx) + 1
        jc = ((z0 - y(1))/ddy) + 1
        r1 = x(ic)
        r2 = x(ic+1)
        z1 = y(jc)
        z2 = y(jc+1)
        do k=1, ncequi
            u1 = zaindk(ic, jc, k)
            u2 = zaindk(ic+1, jc, k)
            u3 = zaindk(ic+1, jc+1, k)
            u4 = zaindk(ic, jc+1, k)
            psitok(k) = psitok(k) + blin_(r0, z0, r1, r2, z1, z2, u1, u2, u3, u4)*cur(i, j)*sqk
        enddo
    enddo
enddo

return
end subroutine flux_r

!---------------------------------------------------------------------
subroutine numcel(rrk, zzk, icell, jcell)

use sp_parameters, only: twopi
use compol, only: nr, nt, nt1, teta, r, z, rm, zm

implicit none

real*8, intent(in) :: rrk, zzk
integer, intent(out) :: icell, jcell

integer :: i, j, ic, jc
real*8 :: ro0, dr0, dz0, tetp, tet0, drc, drv, dzc, dzv, &
    drcb, drvb, dzcb, dzvb, vecpro

dr0 = rrk - rm
dz0 = zzk - zm

ro0 = sqrt(dr0**2 + dz0**2)
tetp = ACOS(dr0/ro0)
if (dz0 < 0.d0) then
    tet0 = twopi - tetp
else
    tet0 = tetp
endif

if (tet0 < teta(1) ) tet0 = tet0 + twopi
if (tet0 > teta(nt)) tet0 = tet0 - twopi
if (tet0 < teta(1) ) tet0 = tet0 + twopi
if (tet0 > teta(nt)) tet0 = tet0 - twopi
if (tet0 < teta(1) ) pause 'numcel: tet0 < teta(1)'
if (tet0 > teta(nt)) pause 'numcel: tet0 > teta(nt)' 

do j=1,nt1
    if (tet0 >= teta(j) .AND. tet0 < teta(j+1)) jc = j
enddo

jcell = jc
drvb = rrk - r(nr, jc)
dzvb = zzk - z(nr, jc)
drcb = r(nr, jc+1) -r(nr, jc)
dzcb = z(nr, jc+1) -z(nr, jc)
vecpro = drcb*dzvb - drvb*dzcb
if (vecpro <= 0.d0) then
    icell = nr
else
    do i=nr, 2, -1
        drv = rrk - r(i, jc)
        dzv = zzk - z(i, jc)
        drc = r(i, jc+1) - r(i, jc)
        dzc = z(i, jc+1) - z(i, jc)
        vecpro=drc*dzv-drv*dzc
        if (vecpro > 0.d0) then
            ic=i-1
        else
            EXIT
        endif
    enddo
    icell = ic
endif

return
end subroutine numcel

!---------------------------------------------------------------------
real*8 function blin_tr(tet0, ro0, tet1, tet2, &
    ro1, ro2, ro3, ro4, u1, u2, u3, u4)

implicit none

real*8, intent(in) :: tet0, ro0, tet1, tet2, &
    ro1, ro2, ro3, ro4, u1, u2, u3, u4
real*8 :: ro14, ro23, u14, u23

ro14 = (ro1*(tet2 - tet0) + ro4*(tet0 - tet1))/(tet2 - tet1)
ro23 = (ro2*(tet2 - tet0) + ro3*(tet0 - tet1))/(tet2 - tet1)
u14 = (u1*(tet2 - tet0) + u4*(tet0 - tet1))/(tet2 - tet1)
u23 = (u2*(tet2 - tet0) + u3*(tet0 - tet1))/(tet2 - tet1)
blin_tr = (u14*(ro23 - ro0) + u23*(ro0 - ro14))/(ro23 - ro14)

return
end function blin_tr
