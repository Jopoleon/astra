subroutine renet

use sp_parameters, only: twopi
use comrec
use keys, only: kpr
use compol_add, only: ron_max_g
use compol, only: iplas, nr, nt, nt1, ro, ronor, psia, teta, rm, zm, r, z

implicit none

integer :: i, j, i_bc, key_ang
real*8 :: delro, ropla, robon, rb_min, rb_max, zb_min, zb_max, tet0

do i=iplas, 1, -1
    if (psia(i) >= 0.05d0) then
        i_bc = i
        exit
    endif
enddo

do j=2, nt1
    ronor(nr, j) = ron_max_g(j)
    ro(nr, j) = ro(i_bc, j)*ronor(nr, j)
enddo

do i=2, iplas
    do j=2, nt1
        ro(i, j) = ro(iplas, j)*ronor(i, j)
    enddo
enddo

do i=iplas+1, nr-1
    do j=2, nt1
        ropla = ro(iplas, j)
        robon = ro(   nr, j)
        delro = (robon - ropla)/(nr - iplas)
        ro(i, j) = ropla + delro*(i - iplas)
        ronor(i, j) = ro(i, j)/ropla
    enddo
enddo

do i=2, nr
    do j=2, nt1
        r(i, j) = ro(i, j)*COS(teta(j)) + rm
        z(i, j) = ro(i, j)*SIN(teta(j)) + zm
    enddo
enddo

do i=1, nr
    ro(i,  1) = ro(i, nt1)
    ro(i, nt) = ro(i, 2)
    ronor(i,  1) = ronor(i, nt1)
    ronor(i, nt) = ronor(i, 2)
    r(i,  1) = r(i, nt1)
    r(i, nt) = r(i, 2)
    z(i,  1) = z(i, nt1)
    z(i, nt) = z(i, 2)
enddo

rb_max = rm
rb_min = rm
zb_max = zm
zb_min = zm

do j=2, nt1
    if (r(nr, j) > rb_max) rb_max = r(nr, j)
    if (r(nr, j) < rb_min) rb_min = r(nr, j)
    if (z(nr, j) > zb_max) zb_max = z(nr, j)
    if (z(nr, j) < zb_min) zb_min = z(nr, j)
enddo

if (rb_max > xmax .or. rb_min < xmin .or. &
    zb_max > ymax .or. zb_min < ymin ) then

    if (kpr == 1) then
        write(*, *) 'egg is larger than rectangular box'
        call f_wrd
        write(*, *) 'STOP'
        call err_catch_a
        stop
    endif
endif

! angle theta redefinition
key_ang = 1
if (key_ang == 1) then
    tet0 = teta(2)
    do j=2, nt1
        teta(j) = tet0 + twopi*(j-2)/(nt - 2.d0)
    enddo
    teta( 1) = teta(nt1) - twopi
    teta(nt) = teta(  2) + twopi

    do i=2, nr
        do j=2, nt1
            r(i, j) = ro(i, j)*COS(teta(j)) + rm
            z(i, j) = ro(i, j)*SIN(teta(j)) + zm
        enddo
    enddo

    do i=1, nr
        ro(i,  1) = ro(i, nt1)
        ro(i, nt) = ro(i, 2)
        ronor(i,  1) = ronor(i, nt1)
        ronor(i, nt) = ronor(i, 2)
        r(i,  1) = r(i, nt1)
        r(i, nt) = r(i, 2)
        z(i,  1) = z(i, nt1)
        z(i, nt) = z(i, 2)
    enddo
endif

return
end subroutine renet

!---------------------------------------------------------------------
subroutine f_remesh(erro)

use compol, only: iplas1, nt1, psi

implicit none

real*8, intent(out) :: erro

integer :: i, j, imax, jmax
real*8 :: psimax

psimax = psi(1, 2)
imax = 1
jmax = 2

do i=2, iplas1
    do j=2, nt1
        if (psi(i, j) > psimax) then
            psimax = psi(i, j)
            imax = i
            jmax = j
        endif
    enddo
enddo

if (imax >= 2) then
    call f_prgrid(imax, jmax, erro)
else
    call f_regrid0(erro)
endif

return
end subroutine f_remesh

!---------------------------------------------------------------------
subroutine f_prgrid(imax, jmax, erro)

use sp_parameters, only: twopi, nrp, ntp
use keys, only: kpr
use compol_add, only: nctrl, rx0, rx1, rx2, zx0, zx1, zx2, ixp1, ixp2, &
    jxp1, jxp2, psix0, psix1, alp, numlim, nblm, rblm, zblm
use compol, only: nr, nr1, nt, nt1, iplas, iplas1, itin, r, z, rm, zm, &
    psi, psia, psim, psin, psip, ro, teta

implicit none

integer, parameter :: nshp=ntp+1

integer, intent(in) :: imax, jmax
real*8, intent(out) :: erro

integer :: i, j, k, l, is, isc, isw, ibeg, jj, llim, nsh, nctrli, &
    kodex1, kodex2, kodxp, ierm, jerm, npro
real*8 :: delpsn, det, errox, errpsi, grad, spro, zps, zwgt, &
    rimjm, zimjm, rm0, zm0, rmx, zmx, rxold, zxold, drx, dzx, rnj, znj, &
    psim0, psix2, psipn, psipr, psiblm, psnn, psinm0, pspl, psa1, psa2, &
    psai, ps0, ps1, ps2, tetx1, tetx2, tetp, tetv, tt0, ttpl, ttj, &
    robj, ro2j, ro0, ropl, zro, rodeb, ro1, ro2, ros1, ros2, rosi, &
    row, dro, droj, roerr
real*8, dimension(5) :: dp
real*8, dimension(nshp) :: xs, ys, fun
real*8, dimension(nrp) :: roh, rh, zh
real*8, dimension(ntp) :: ron, rop, rn, zn, tetn, rob, teti, roi
real*8, dimension(nrp, ntp) :: roplt, psiold

nctrli = nctrl
npro = nr - iplas + 1

rimjm = r(imax, jmax)
zimjm = z(imax, jmax)
rm0 = r(1, 1)
zm0 = z(1, 1)
psim0 = psi(1, 1)
nsh = 1
xs(nsh) = r(imax, jmax)
ys(nsh) = z(imax, jmax)
fun(nsh) = psi(imax, jmax)

if (imax == 1) then
    do j=2, nt1
        nsh = nsh + 1
        xs(nsh) = r(2, j)
        ys(nsh) = z(2, j)
        fun(nsh) = psi(2, j)
    enddo
elseif (imax == 2) then
    nsh = nsh + 1
    xs(nsh) = r(1, 2)
    ys(nsh) = z(1, 2)
    fun(nsh) = psi(1, 2)
    do j=2, nt1
        nsh = nsh + 1
        xs(nsh) = r(3, j)
        ys(nsh) = z(3, j)
        fun(nsh) = psi(3, j)
    enddo
else
    do k=-1, 1
        i =  imax + k
        do l=-1, 1
            j =  jmax + l
            if (i /= imax .OR. j /= jmax) then
                nsh = nsh + 1
                xs(nsh) = r(i, j)
                ys(nsh) = z(i, j)
                fun(nsh) = psi(i, j)
            endif
        enddo
    enddo
endif

call deriv5(xs, ys, fun, nsh, 5, dp)

DET = dp(3)*dp(5) - dp(4)**2
Rm = Xs(1) + (dp(2)*dp(4) - dp(1)*dp(5))/DET
Zm = Ys(1) + (dp(1)*dp(4) - dp(2)*dp(3))/DET

erro = SQRT((rm - rm0)**2 + (zm - zm0)**2)

rmx = r(imax, jmax)
zmx = z(imax, jmax)

psim = fun(1) + &
     dp(1)*(rm - rmx) + dp(2)*(zm - zmx) + 0.5d0*dp(3)*(rm - rmx)**2 + &
             dp(4)*(rm - rmx)*(zm - zmx) + 0.5d0*dp(5)*(zm - zmx)**2

if (kpr == 1) then
    write(*, *) 'rm zm psim', rm, zm, psim
    write(*, *) 'prgrid: errma', erro
endif

! new position of x-point
rxold = rx0
zxold = zx0

call f_xpoint(rx1, zx1, ixp1, jxp1, psix1, tetx1, kodex1)
call f_xpoint(rx2, zx2, ixp2, jxp2, psix2, tetx2, kodex2)

if (kpr == 1) then
     write(*, *) 'kodex1, kodex2', kodex1, kodex2
endif

kodxp = kodex1*kodex2
if (kodxp /= 0) then
    if (kpr == 1) then
        write(*, *) ' all x-points out of box'
        write(*, *) ' only limiter case '
    endif
    nctrli = 1
    psip = -1.d12
else
    psipn = psip

    if (kodex1 == 0 .AND. kodex2 == 0) then
        if (psix1 >= psix2) then
            psix0 = psix1
            rx0 = rx1
            zx0 = zx1
        else
            psix0 = psix2
            rx0 = rx2
            zx0 = zx2
        endif
    elseif (kodex1 == 0) then
        psix0 = psix1
        rx0 = rx1
        zx0 = zx1
    elseif (kodex2 == 0) then
        psix0 = psix2
        rx0 = rx2
        zx0 = zx2
    endif
    errox = SQRT((rxold - rx0)**2 + (zxold - zx0)**2)
    if (kpr == 1) then
        write(*, *) 'erro x-point', errox
    endif

    psip = psim - alp*(psim - psix0)
    if (itin > 1) then
        zwgt = 1.0d0
        psipr = zwgt*psip + (1.d0 - zwgt)*psipn
        psip = dmax1(psip, psipr)
    endif
endif

if (nctrl == 1) then
    numlim = 0
    do llim=1, nblm
        call fdefln(psiblm, rblm(llim), zblm(llim))
        if (psiblm > psip ) then
            psip = psiblm
            numlim = llim
        endif
    enddo
endif
if (kpr == 1) then
    write(*, *) 'numlim', numlim
endif

! nomalized poloidal flux definition
errpsi = 0.d0

do i=1, nr
    do j=1, nt
        psnn = psin(i, j)
        psin(i, j) = (psi(i, j) - psip)/(psim - psip)
        delpsn = dabs(psin(i, j) - psnn)
        if (delpsn > errpsi) then
            ierm = i
            jerm = j
            errpsi = delpsn
        endif
    enddo
enddo

do j=1, nt
    tetn(j) = teta(j)
enddo

do j=1, nt
    drx = r(nr, j) - rm
    dzx = z(nr, j) - zm
    robj = SQRT(drx**2 + dzx**2)
    rob(j) = robj
    tetp = dacos(drx/robj)
    if (dzx < 0.d0) then
        teta(j) = -tetp
    else
        teta(j) = tetp
    endif
enddo

do j=2, nt
    if (teta(j) < teta(j-1)) then
        teta(j) = teta(j) + twopi
    endif
    if (teta(j) < teta(j-1)) then
        teta(j) = teta(j) + twopi
    endif
enddo

teta(1) = teta(nt1) - twopi
teta(nt) = teta(2) + twopi

! grid moving along rais

psinm0 = psin(1, 2)

do i=3, iplas1
    if (psia(i) < psinm0) then
        ibeg = i + 1
        EXIT
    endif
enddo

do i=ibeg, iplas
    zps = psia(i)
    do j=1, nt
        do is=1, nr1
            ps0  = psin(is, j)
            pspl = psin(is+1, j)
            ro0  = ro(is, j)
            ropl = ro(is+1, j)
            if (zps <= ps0 .and. zps >= pspl) then
                isc = is
                EXIT
            endif
        enddo
        grad = -(pspl - ps0)/(ropl - ro0)
        zro = ro0 - (psia(i) - ps0)/grad
        ron(j) = zro
    enddo
    do j=1, nt
        rnj = rm0 + ron(j)*COS(tetn(j))
        znj = zm0 + ron(j)*SIN(tetn(j))
        drx = rnj - rm
        dzx = znj - zm
        roi(j) = SQRT(drx**2 + dzx**2)
        tetp = dacos(drx/roi(j))
        if (dzx < 0.d0) then
            teti(j) = -tetp
        else
            teti(j) = tetp
        endif
    enddo
    do j=2, nt
        if (teti(j) < teti(j-1)) then
            teti(j) = teti(j) + twopi
        endif
        if (teti(j) < teti(j-1)) then
            teti(j) = teti(j) + twopi
        endif
    enddo

    teti(1) = teti(nt1) - twopi
    teti(nt) = teti(2) + twopi
    roi(1) = roi(nt1)
    roi(nt) = roi(2)
    do j=2, nt1
        tetv = teta(j)
        if (tetv < teti( 1)) tetv = tetv + twopi
        if (tetv < teti( 1)) tetv = tetv + twopi
        if (tetv > teti(nt)) tetv = tetv - twopi
        if (tetv > teti(nt)) tetv = tetv - twopi
        do jj=1, nt1
            tt0  = teti(jj)
            ttpl = teti(jj+1)
            ro0  = roi(jj)
            ropl = roi(jj+1)
            if (tetv <= ttpl .and. tetv >= tt0) then
                zro = ((ttpl - tetv)*ro0 + (tetv - tt0)*ropl)/(ttpl - tt0)
                EXIT
            endif
        enddo
        rop(j)  = zro
        r(i, j) = rop(j)*COS(teta(j)) + rm
        z(i, j) = rop(j)*SIN(teta(j)) + zm
    enddo
enddo

i = 2
zps = psia(i)
ps0 = psip + zps*(psim - psip)

do j=1, nt
    ttj=teta(j)
    ro2j = (ps0 - psim)/(0.5d0*dp(3)*COS(ttj)**2        + &
                               dp(4)*COS(ttj)*SIN(ttj) + &
                         0.5d0*dp(5)*SIN(ttj)**2)
    rodeb = ro2j
    rop(j) = SQRT(ro2j)
    r(i, j) = rm + rop(j)*COS(teta(j))
    z(i, j) = zm + rop(j)*SIN(teta(j))
enddo

if (ibeg > 3) then
    do j=2, nt1
        ros1 = (r(2, j) - rm)**2 + (z(2, j) - zm)**2
        ros2 = (r(ibeg, j) - rm)**2 + (z(ibeg, j) - zm)**2
        psa1 = psia(2)
        psa2 = psia(ibeg)
        do i=3, ibeg-1
            psai = psia(i)
            rosi = ( (psa2 - psai)*ros1 + (psai - psa1)*ros2 )/(psa2 - psa1)
            rosi = SQRT(rosi)
            r(i, j) = rm + rosi*COS(teta(j))
            z(i, j) = zm + rosi*SIN(teta(j))
        enddo
    enddo
endif

do i=2, iplas
    do j=2, nt1
        ro(i, j) = SQRT((r(i, j) - rm)**2 + (z(i, j) - zm)**2)
    enddo
enddo

do j=2, nt1
    spro = rob(j) - ro(iplas-1, j)
    droj = ro(iplas, j) - ro(iplas-1, j)
    row = ro(iplas, j)
    dro = (rob(j) - ro(iplas, j))/(nr - iplas)
    do i=iplas+1, nr1
        row = row + dro
        rh(i) = rm + row*COS(teta(j))
        zh(i) = zm + row*SIN(teta(j))
        roh(i) = row
        do is=2, nr1
            ps0  = psin(is  , j)
            pspl = psin(is+1, j)
            ro0  = ro(is  , j)
            ropl = ro(is+1, j)
            if (row <= ropl .and. row >= ro0) then
                isw = is
                EXIT
            endif
        enddo
        ro2 = ropl
        ro1 = ro0
        ps2 = pspl
        ps1 = ps0
        psia(i) = ( ps1*(ro2 - row) + ps2*(row - ro1) )/(ro2 - ro1)
    enddo
    do i=iplas+1, nr1
        roerr = dabs(roh(i) - ro(i, j))
        erro = dmax1(erro, roerr)
        ro(i, j) = roh(i)
        r( i, j) = roh(i)*COS(teta(j)) + rm
        z( i, j) = roh(i)*SIN(teta(j)) + zm
        psin(i, j) = psia(i)
    enddo
enddo

do j=2, nt1
    ro(nr, j) = rob(j)
    ro( 1, j) = 0.d0
    r(1, j) = rm
    z(1, j) = zm
enddo

do i=1, nr
    ro(i,  1) = ro(i, nt1)
    ro(i, nt) = ro(i, 2)
    r(i,  1) = r(i, nt1)
    r(i, nt) = r(i, 2)
    z(i,  1) = z(i, nt1)
    z(i, nt) = z(i, 2)
    roplt(i,  1) = roplt(i, nt1)
    roplt(i, nt) = roplt(i, 2)
enddo

do i=1, iplas
    do j=1, nt
        psin(i, j) = psia(i)
    enddo
enddo

do i=1, nr
    do j=1, nt
        psiold(i, j) = psip + psin(i, j)*(psim - psip)
    enddo
enddo

call f_wrd

return
end subroutine f_prgrid

!---------------------------------------------------------------------
subroutine f_regrid0(erro)

use sp_parameters, only: nrp, ntp, twopi
use keys, only: kpr
use compol_add, only: nctrl, rx0, rx1, rx2, zx0, zx1, zx2, ixp1, ixp2, &
    jxp1, jxp2, psix0, psix1, alp, numlim, nblm, rblm, zblm
use compol, only: nr, nr1, nt, nt1, iplas, itin, ngav, &
    r, z, rm, zm, ro, ronor, teta, &
    psi, psia, psim, psin, psip

implicit none

real*8, intent(out) :: erro

integer :: i, j, jj, is, isw, ierm, jerm, nctrli, npro, llim, &
    kodex1, kodex2, kodxp
real*8 :: alfa, alp_pr, delpsn, errpsi, spro, errox, zwgt, zps, &
    tetx1, tetx2, ttj, tt0, ttpl, tetp, tetv, ro2j, robj, &
    ro0, ro1, ro2, ropl, roerr, row, dro, droj, zro, rmold, zmold, &
    rma, zma, rxold, zxold, drx, dzx, rnj, znj, &
    ps0, ps1, ps2, psiblm, psima, psix2, psipn, psipr, psnn, pspl
real*8, dimension(5) :: dp
real*8, dimension(nrp) :: roh, rh, zh
real*8, dimension(ntp) :: ron, rop, rn, zn, tetn, rob, teti, roi
real*8, dimension(nrp, ntp) :: roplt, psiold

save psiold

nctrli = nctrl
npro = nr - iplas + 1

if (ngav == 0) then
    alfa = 1.0d0
else
    alfa = 0.75d0
endif

rmold = rm
zmold = zm

if (kpr == 1) then
    write(*, *) 'rm, zm, psim', rm, zm, psim
endif

call axdef(rma, zma, psima, dp)

if (kpr == 1) then
    write(*, *) 'rma, zma, psima', rma, zma, psima
endif
if (isnan(rma)) then
    write(*, *) 'mag axis major radius is NaN'
    call err_catch_a
endif

rm = rma
zm = zma
psim = psima
erro = SQRT((rmold - rm)**2 + (zmold - zm)**2)

if (kpr == 1) then
    write(*, *) 'erro mag.axis', erro
endif

errpsi = 0.d0
rxold = rx0
zxold = zx0

call f_xpoint(rx1, zx1, ixp1, jxp1, psix1, tetx1, kodex1)
call f_xpoint(rx2, zx2, ixp2, jxp2, psix2, tetx2, kodex2)

if (kpr == 1) then
    write(*, *) 'kodex1, kodex2', kodex1, kodex2
endif

kodxp = kodex1*kodex2
if (kodxp /= 0) then
    if (kpr == 1) then
        write(*, *) ' all x-points out of box'
        write(*, *) ' only limiter case '
    endif
    nctrli = 1
    psip = -1.d12
else
    psipn = psip
    if (kodex1 == 0 .AND. kodex2 == 0) then
        if (psix1 >= psix2) then
            psix0 = psix1
            rx0 = rx1
            zx0 = zx1
        else
            psix0 = psix2
            rx0 = rx2
            zx0 = zx2
        endif
    elseif (kodex1 == 0) then
        psix0 = psix1
        rx0 = rx1
        zx0 = zx1
    elseif (kodex2 == 0) then
        psix0 = psix2
        rx0 = rx2
        zx0 = zx2
    endif

    errox = SQRT((rxold - rx0)**2 + (zxold - zx0)**2)
    if (kpr == 1) then
        write(*, *) 'erro x-point', errox
        write(*, *) 'rx zx', rx0, zx0, psix0
    endif

    psip = psim - alp*(psim - psix0)

    if (itin > 1) then
        zwgt = 1.0d0
        psipr = zwgt*psip + (1.d0 - zwgt)*psipn
        psip = dmax1(psip, psipr)
        alp_pr = (psim - psip)/(psim - psix0)
    endif

endif

if (nctrl == 1) then
    numlim = 0
    do llim=1, nblm
        call fdefln(psiblm, rblm(llim), zblm(llim))
        if (psiblm > psip) then
            psip = psiblm
            numlim = llim
        endif
    enddo
endif
if (kpr == 1) then
    write(*, *) 'numlim', numlim
endif

do i=1, nr
    do j=1, nt
        psnn = psin(i, j)
        psin(i, j) = (psi(i, j) - psip)/(psim - psip)
        delpsn = dabs(psin(i, j) - psnn)
        if (delpsn > errpsi) then
            ierm = i
            jerm = j
            errpsi = delpsn
        endif
    enddo
enddo

if (kpr == 1) then
    write(*, *) 'i, j errpsi', ierm, jerm, errpsi, (psim - psip)
endif

! definition of new angle grid teta(j)

do j=1, nt
    tetn(j) = teta(j)
enddo

do j=1, nt
    drx = r(nr, j) - rm
    dzx = z(nr, j) - zm
    robj = SQRT(drx**2 + dzx**2)
    rob(j) = robj
    tetp = dacos(drx/robj)
    if (dzx < 0.d0) then
        teta(j) = -tetp
    else
        teta(j) = tetp
    endif
enddo

do j=2, nt
    if (teta(j) < teta(j-1)) then
        teta(j) = teta(j) + twopi
    endif
    if (teta(j) < teta(j-1)) then
        teta(j) = teta(j) + twopi
    endif
enddo

teta( 1) = teta(nt1) - twopi
teta(nt) = teta(  2) + twopi

i = 2

zps  = psia(i)
ps0 = psip + zps*(psim - psip)

do j=1, nt
    ttj = teta(j)
    ro2j = (ps0 - psim)/( 0.5d0*dp(3)*COS(ttj)**2 + &
                                dp(4)*COS(ttj)*SIN(ttj) + &
                          0.5d0*dp(5)*SIN(ttj)**2 )
    rop(j) = SQRT(ro2j)
    r(i, j) = rm + rop(j)*COS(teta(j))
    z(i, j) = zm + rop(j)*SIN(teta(j))
enddo

! grid moving along rais

do i=3, iplas
    zps = psia(i)
    do j=1, nt
        do is=1, nr1
            ps0  = psin(is, j)
            pspl = psin(is+1, j)
            ro0  = ro(is, j)
            ropl = ro(is+1, j)
            if (zps <= ps0 .and. zps >= pspl) then
                zro = ((pspl - zps)*ro0 - (ps0 - zps)*ropl)/(pspl - ps0)
                EXIT
            endif
        enddo
        ron(j) = zro
    enddo

    do j=1, nt
        rnj = rmold + ron(j)*COS(tetn(j))
        znj = zmold + ron(j)*SIN(tetn(j))
        drx = rnj - rm
        dzx = znj - zm
        roi(j) = SQRT(drx**2 + dzx**2)
        tetp = dacos(drx/roi(j))
        if (dzx < 0.d0) then
            teti(j) = -tetp
        else
            teti(j) = tetp
        endif
    enddo

    do j=2, nt
        if (teti(j) < teti(j-1)) then
            teti(j) = teti(j) + twopi
        endif
        if (teti(j) < teti(j-1)) then
            teti(j) = teti(j) + twopi
        endif
    enddo

    teti( 1) = teti(nt1) - twopi
    teti(nt) = teti(  2) + twopi
    roi( 1) = roi(nt1)
    roi(nt) = roi(2)

    do j=2, nt1
        tetv = teta(j)
        if (tetv < teti(1) ) tetv = tetv + twopi
        if (tetv < teti(1) ) tetv = tetv + twopi
        if (tetv > teti(nt)) tetv = tetv - twopi
        if (tetv > teti(nt)) tetv = tetv - twopi
        do jj=1, nt1
            tt0  = teti(jj)
            ttpl = teti(jj+1)
            ro0  = roi(jj)
            ropl = roi(jj+1)
            if (tetv <= ttpl .and. tetv >= tt0) then
                zro = ((ttpl - tetv)*ro0 + (tetv - tt0)*ropl)/(ttpl - tt0)
                EXIT
            endif
        enddo
        rop(j) = zro
        roerr = dabs(ron(j) - ro(i, j))
        if (roerr > erro) then
            erro = roerr
            ierm = i
            jerm = j
        endif
        roplt(i, j) = ron(j) - ro(i, j)
        r(i, j) = rop(j)*COS(teta(j)) + rm
        z(i, j) = rop(j)*SIN(teta(j)) + zm
    enddo
enddo

if (kpr == 1) then
    write(6, *) 'erro max.', erro
    write(6, *) 'i, j erro', ierm, jerm, erro
endif

do i=2, iplas
    do j=2, nt1
        ro(i, j) = SQRT((r(i, j) - rm)**2 + (z(i, j) - zm)**2)
        ronor(i, j) = ro(i, j)/rob(j)
        psin(i, j) = psia(i)
    enddo
enddo

do j=2, nt1
    spro = rob(j) - ro(iplas-1, j)
    droj = ro(iplas, j) - ro(iplas-1, j)
    row = ro(iplas, j)
    dro = (rob(j) - ro(iplas, j))/(nr - iplas)
    do i=iplas+1, nr1
        row = row + dro
        rh(i) = rm + row*COS(teta(j))
        zh(i) = zm + row*SIN(teta(j))
        roh(i) = row
        do is=2, nr1
            ps0  = psin(is, j)
            pspl = psin(is+1, j)
            ro0  = ro(is, j)
            ropl = ro(is+1, j)
            if (row <= ropl .and. row >= ro0) then
                isw = is
                EXIT
            endif
        enddo
        ro2 = ropl
        ro1 = ro0
        ps2 = pspl
        ps1 = ps0
        psia(i) = ( ps1*(ro2 - row) + ps2*(row - ro1) )/(ro2 - ro1)
    enddo

    do i=iplas+1, nr1
        roerr = dabs(roh(i) - ro(i, j))
        erro = dmax1(erro, roerr)
        ro(i, j) = roh(i)
        r(i, j) = roh(i)*COS(teta(j)) + rm
        z(i, j) = roh(i)*SIN(teta(j)) + zm
        psin(i, j) = psia(i)
    enddo
enddo

do j=2, nt1
    ro(nr, j) = rob(j)
    ro(1, j) = 0.d0
    r(1, j) = rm
    z(1, j) = zm
enddo

do i=1, nr
    do j=2, nt1
        ronor(i, j) = ro(i, j)/ro(iplas, j)
    enddo
enddo

do i=1, nr
    ro(i, 1) = ro(i, nt1)
    ro(i, nt) = ro(i, 2)
    ronor(i, 1) = ronor(i, nt1)
    ronor(i, nt) = ronor(i, 2)
    r(i, 1) = r(i, nt1)
    r(i, nt) = r(i, 2)
    z(i, 1) = z(i, nt1)
    z(i, nt) = z(i, 2)
    roplt(i, 1) = roplt(i, nt1)
    roplt(i, nt) = roplt(i, 2)
enddo

do i=1, iplas
    do j=1, nt
        psin(i, j) = psia(i)
    enddo
enddo

do i=1, nr
    do j=1, nt
        psiold(i, j) = psip + psin(i, j)*(psim - psip)
    enddo
enddo

return
end subroutine f_regrid0

!---------------------------------------------------------------------
subroutine grid_spdr

use comrec
use keys, only: kpr
use compol_add, only: ron_max_g, rx0, rx1, rx2, zx0, zx1, zx2, ixp1, ixp2, jxp1, jxp2
use compol, only: nr, nt, nt1, iplas, ro, teta, ronor, r, z, rm, zm, &
    psia, psim, psip

implicit none

integer :: i, j, i_bc
real*8 :: robon, delro, ropla, rb_min, rb_max, zb_min, zb_max

rm = r(1, 2)
zm = z(1, 2)

do i=iplas, 1, -1
    if (psia(i) >= 0.05d0) then
        i_bc = i
        exit
    endif
enddo

do j=2, nt1
    ronor(nr, j) = ron_max_g(j)
    ro(nr, j) = ro(i_bc, j)*ronor(nr, j)
enddo

do i=iplas+1, nr-1
    do j=2, nt1
        ropla = ro(iplas, j)
        robon = ro(nr, j)
        delro = (robon - ropla)/(nr - iplas)
        ro(i, j) = ropla + delro*(i - iplas)
        ronor(i, j) = ro(i, j)/ropla
    enddo
enddo

do i=iplas+1, nr
    do j=2, nt1
        r(i, j) = ro(i, j)*COS(teta(j)) + rm
        z(i, j) = ro(i, j)*SIN(teta(j)) + zm
    enddo
enddo

do i=1, nr
    ro(i,  1) = ro(i, nt1)
    ro(i, nt) = ro(i, 2)
    ronor(i,  1) = ronor(i, nt1)
    ronor(i, nt) = ronor(i, 2)
    r(i,  1) = r(i, nt1)
    r(i, nt) = r(i, 2)
    z(i,  1) = z(i, nt1)
    z(i, nt) = z(i, 2)
enddo

rb_max = rm
rb_min = rm
zb_max = zm
zb_min = zm

do j=2, nt1
    if (r(nr, j) > rb_max) rb_max = r(nr, j)
    if (r(nr, j) < rb_min) rb_min = r(nr, j)
    if (z(nr, j) > zb_max) zb_max = z(nr, j)
    if (z(nr, j) < zb_min) zb_min = z(nr, j)
enddo

if (rb_max > xmax .or. rb_min < xmin .or. &
    zb_max > ymax .or. zb_min < ymin ) then
    if (kpr == 1) then
        write(*, *) 'egg is large then rectangular box'
        call f_wrd
        call err_catch_a
    endif
endif

rx1 = xx1  !   x-points from rect. grid
zx1 = yx1

write(*, *) 'jjj ', rx1, zx1, xx2, yx2

rx2 = xx2
zx2 = yx2
rx0 = xx0
zx0 = yx0
ixp1 = 0
jxp1 = 0
ixp2 = 0
jxp2 = 0
psim = um
psip = up

return
end subroutine grid_spdr

!---------------------------------------------------------------------
subroutine ada(erro)

use sp_parameters, only: nrp, twopi
use keys, only: kpr
use compol_add, only: rx0, rx1, rx2, zx0, zx1, zx2, ixp1, ixp2, jxp1, jxp2, &
    psix0, psix1, nctrl, alp, alpnew, jrolim, rolim, numlim, nblm, rblm, zblm
use compol, only: nr, nr1, nt, nt1, iplas, ngav, ro, ronor, teta, r, z, rm, zm, &
    psi, psia, psim, psin, psip

implicit none

real*8, intent(out) :: erro
real*8, dimension(nrp) :: ron, rn, zn

integer :: i, j, is, isc, isw, llim, nctrli, npro, kodex1, kodex2, kodxp
real*8 :: tetx1, tetx2, psix2, psiblm, ps0, pspl, psmn, &
    grad, gradpl, gradmn, ro0, ropl, romn, zps, zro, alfa, &
    spro, dro, droj, row, ro1, ro2, roerr, ps1, ps2

nctrli = nctrl
npro = nr - iplas + 1
erro = 0.d0
psim = psi(1, 1)

call f_xpoint(rx1, zx1, ixp1, jxp1, psix1, tetx1, kodex1)

if (kpr == 1) then
    write(*, *) 'kodex1', kodex1
    write(*, *) 'rx1 zx1', rx1, zx1, psix1
endif

call f_xpoint(rx2, zx2, ixp2, jxp2, psix2, tetx2, kodex2)

if (kpr == 1) then
    write(*, *) 'kodex2', kodex2
    write(*, *) 'rx2 zx2', rx2, zx2, psix2
endif

kodxp = kodex1*kodex2

if (kodxp /= 0) then
    if (kpr == 1) then
        write(*, *) ' all x-points out of box'
        write(*, *) ' only limiter case '
    endif
    nctrli = 1
 psip = -1.d12
else
    if (kodex1 == 0 .AND. kodex2 == 0) then
        if (psix1 > psix2) then
            psix0 = psix1
            rx0 = rx1
            zx0 = zx1
        else
            psix0 = psix2
            rx0 = rx2
            zx0 = zx2
        endif
    elseif (kodex1 == 0) then
        psix0 = psix1
        rx0 = rx1
        zx0 = zx1
    elseif (kodex2 == 0) then
        psix0 = psix2
        rx0 = rx2
        zx0 = zx2
    endif

    if (kpr == 1) then
        write(*, *) 'kodex1, kodex2', kodex1, kodex2
        write(*, *) 'rx zx', rx0, zx0, psix0
    endif

    if (ngav/10*10 == ngav) then
        psip = psim - alp*(psim - psix0)
    else
        psip = psi(iplas, jrolim)
    endif
endif

numlim = 0
if (nctrl == 1) then
    do llim=1, nblm
        call fdefln(psiblm, rblm(llim), zblm(llim))
        if (psiblm > psip) then
            psip = psiblm
            numlim = llim
        endif
    enddo
endif

alpnew = (psim - psip)/(psim - psix0)

do i=1, nr
    do j=1, nt
        psin(i, j) = (psi(i, j) - psip)/(psim - psip)
    enddo
enddo

if (kpr == 1) then
    write(*, *) 'ada:', psim, psip
    write(*, *) 'ada: nctrli, nctrl', nctrli, nctrl
    write(*, *) 'ada: alpnew, numlim', alpnew, numlim, psix0
endif

ron(1) = 0.d0

do j=1, nt
    do i=2, iplas
        zps = psia(i)
        do is=1, nr1
            ps0  = psin(is, j)
            pspl = psin(is+1, j)
            ro0  = ro(is, j)
            ropl = ro(is+1, j)
            if (zps <= ps0 .and. zps >= pspl) then
                isc = is
                EXIT
            endif
        enddo
        if (isc == 1) then
            grad = -(pspl - ps0)/(ropl - ro0)
        else
            psmn = psin(isc-1, j)
            romn = ro(isc-1, j)
            gradpl = -(pspl - ps0)/(ropl - ro0)
            gradmn = -(ps0 - psmn)/(ro0 - romn)
            grad = dmax1(gradpl, gradmn)
        endif
        grad = -(pspl - ps0)/(ropl - ro0)
        zro = ro0 - (zps - ps0)/grad
        ron(i) = zro
        if (ngav == 0) then
            alfa = 0.75d0
        else
            alfa = 0.550d0
        endif
        ron(i) = alfa*zro + (1.d0 - alfa)*ro(i, j)
        if (ron(i) < ron(i-1)) then
            ron(i) = ron(i-1) + 1.d-8
            if (kpr == 1) then
                write(*, *) 'grid crash i, j', i, j
                call err_catch_a
            endif
        endif
    enddo

    if (ngav > 0 .AND. j == jrolim .AND. kpr == 1) then
        write(*, *) 'ro(iplas), rolim', ron(iplas), rolim
    endif
    spro = ro(nr, j) - ron(iplas-1)

    droj = ron(iplas) - ron(iplas-1)
    row = ron(iplas)
    dro = (ro(nr, j) - row)/(nr-iplas)

    do i=iplas+1, nr1
        row = row + dro
        rn(i) = rm + row*COS(teta(j))
        zn(i) = zm + row*SIN(teta(j))
        ron(i) = row
        do is=isc, nr1
            ps0  = psin(is, j)
            pspl = psin(is+1, j)
            ro0  = ro(is, j)
            ropl = ro(is+1, j)
            if (row <= ropl .and. row >= ro0) then
                isw = is
                EXIT
            endif
        enddo
        ro2 = ropl
        ro1 = ro0
        ps2 = pspl
        ps1 = ps0
        psia(i) = ( ps1*(ro2 - row) + ps2*(row - ro1) )/(ro2 - ro1)
    enddo
    do i=2, nr1
        roerr = dabs(ron(i) - ro(i, j))
        erro = dmax1(erro, roerr)
        ro(i, j) = ron(i)
        r(i, j) = ron(i)*COS(teta(j)) + rm
        z(i, j) = ron(i)*SIN(teta(j)) + zm
        psin(i, j) = psia(i)
    enddo
enddo

teta(1) = teta(nt1) - twopi
teta(nt) = teta(2) + twopi
do i=1, nr
    ro(i, 1) = ro(i, nt1)
    ro(i, nt) = ro(i, 2)
enddo
do i=1, nr
    do j=2, nt1
        ronor(i, j) = ro(i, j)/ro(iplas, j)
    enddo
enddo

do i=1, nr
    ronor(i,  1) = ronor(i, nt1)
    ronor(i, nt) = ronor(i, 2)
    r(i,  1) = r(i, nt1)
    r(i, nt) = r(i, 2)
    z(i,  1) = z(i, nt1)
    z(i, nt) = z(i, 2)
    psin(i,  1) = psin(i, nt1)
    psin(i, nt) = psin(i, 2)
enddo

do i=1, nr
    do j=1, nt
        psi(i, j) = psip + psin(i, j)*(psim - psip)
    enddo
enddo

return
end subroutine ada

!---------------------------------------------------------------------
subroutine reform

use sp_parameters, only: ntp, pi, twopi
use compol_add, only: jrolim, rolim
use compol, only: nr, nt, nt1, iplas, ngav, ro, teta, r, z, rm, zm

implicit none

integer :: i, j
real*8 :: drx, dzx, tet0, tetp
real*8, dimension(ntp) :: rob

do j=1, nt
    drx = r(nr, j) - rm
    dzx = z(nr, j) - zm
    rob(j) = SQRT(drx**2 + dzx**2)
    tetp = dacos(drx/rob(j))
    if (dzx < 0.d0) then
        tet0 = twopi - tetp
    else
        tet0 = tetp
    endif
    teta(j) = tet0
enddo

do j=2, nt
    if (teta(j) < teta(j-1)) then
        teta(j) = teta(j) + pi
    endif
    if (teta(j) < teta(j-1)) then
        teta(j) = teta(j) + pi
    endif
enddo
teta(1) = teta(nt1) - twopi
teta(nt) = teta(2) + twopi
do j=1, nt
    do i=1, nr
        ro(i, j) = ro(i, j)*rob(j)/ro(nr, j)
        r(i, j) = rm + ro(i, j)*COS(teta(j))
        z(i, j) = zm + ro(i, j)*SIN(teta(j))
    enddo
enddo

if (ngav/10*10 /= ngav) then
    ro(iplas, jrolim) = rolim
    r(iplas, jrolim) = rm + ro(iplas, jrolim)*COS(teta(jrolim))
    z(iplas, jrolim) = zm + ro(iplas, jrolim)*SIN(teta(jrolim))
endif

return
end subroutine reform

!---------------------------------------------------------------------
subroutine fdefln(psidf, rgv, zgv)

use sp_parameters, only: twopi
use compol, only: nr, nt, nt1, iplas, ro, teta, rm, zm, psi

implicit none

real*8, intent(in) :: rgv, zgv
real*8, intent(out) :: psidf

integer :: i, j, ic, jc, i_rlmp
real*8 :: dr0, dz0, der_psi, tetp, tet0, tet1, tet2, rob1, rob2, rob12, &
    ro0, ro1, ro2, ro3, ro4, ro12, ropl1, ropl2, ropl12, u1, u2, u3, u4
real*8, external :: blin_tr

psidf = 0.d0
dr0 = rgv - rm
dz0 = zgv - zm
ro0 = SQRT(dr0**2 + dz0**2)
tetp = dacos(dr0/ro0)
if (dz0 < 0.d0) then
    tet0 = twopi - tetp
else
    tet0 = tetp
endif

if (tet0 < teta( 1)) tet0 = tet0 + twopi
if (tet0 > teta(nt)) tet0 = tet0 - twopi

do j=1, nt1
    if (tet0 >= teta(j) .AND. tet0 < teta(j+1)) jc = j
enddo

rob1 = ro(nr, jc)
rob2 = ro(nr, jc+1)

ropl1 = ro(iplas, jc)
ropl2 = ro(iplas, jc+1)

tet1 = teta(jc)
tet2 = teta(jc+1)

rob12  = (rob1 *(tet2 - tet0) + rob2 *(tet0 - tet1))/(tet2 - tet1)
ropl12 = (ropl1*(tet2 - tet0) + ropl2*(tet0 - tet1))/(tet2 - tet1)

if (ro0 < rob12) then
    do i=nr, 2, -1
        ro1 = ro(i, jc)
        ro2 = ro(i, jc+1)
        ro12 = (ro1*(tet2 - tet0) + ro2*(tet0 - tet1))/(tet2 - tet1)
        if (ro0 < ro12) ic = i - 1
    enddo

    i_rlmp = 1
    if (ic > iplas) then
        do i=ic+1, iplas-1, -1
            der_psi = psi(i, jc) - psi(i-1, jc)
            if (der_psi > 0.d0) then
                i_rlmp = 0
                psidf = -1.0d12
            endif
        enddo
    endif

    if (i_rlmp == 1) then
        ro1 = ro(ic, jc)
        ro2 = ro(ic+1, jc)
        ro3 = ro(ic+1, jc+1)
        ro4 = ro(ic, jc+1)
        u1 = psi(ic, jc)
        u2 = psi(ic+1, jc)
        u3 = psi(ic+1, jc+1)
        u4 = psi(ic, jc+1)
        psidf = blin_tr(tet0, ro0, tet1, tet2, ro1, ro2, ro3, ro4, u1, u2, u3, u4)
    endif
else
    psidf = -1.0d12
endif

return
end subroutine fdefln
