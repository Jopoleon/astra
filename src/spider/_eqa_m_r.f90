subroutine eqa(keyctr, igdf, platok, psax, i_betp, betplx, &
    rax, zax, rxpnt, zxpnt, psbo, psdel, ncequi, psitok, &
    betpol, nflag, errarr)

use keys, only: kpr
use compol_add, only: rolim, jrolim, ich, g, psii, clr, clz, &
    itrmax, nitmax, rx0, zx0, fpv
use compol, only: nr, nt, iplas, ngav, iter, itin, &
    psi, psiax, psibon, psim, psip, psipla, &
    rm, zm, f, fvac, cnor, ro, tok, tokp
use status_inc, only: error_catch

implicit none

integer, intent(in) :: keyctr, igdf, i_betp, ncequi
real*8, intent(in) :: betplx, betpol
real*8, intent(in), dimension(*) :: psitok

integer, intent(out) :: nflag
real*8, intent(out) :: platok, rax, zax, rxpnt, zxpnt, psax, psbo, psdel
real*8, intent(out), dimension(*) :: errarr

integer :: iplasm, nroi, ntetj, itout, ish, nshift
integer, dimension(4) :: iwrk
real*8 :: epscrz, epspsm, epsfpv, errpsm, errfpv, errpsb, epsro, &
    erro, cab, fvv, fpv0, det, &
    dell, delr, delz, delro, delfpv, drolim, &
    cr0, cr1, cr2, cr3, cr4, cz0, cz1, cz2, cz3, cz4, &
    pm0, pm1, pm2, pm3, pm4, fv0, fv1, fv2, fv3, fv4, &
    rm0, zm0, rolim0, drm, dfpv, dfvdr, dfvdz, dfvdro, dfvdf, &
    dcrdr, dcrdz, dczdz, dczdr, dpmdr, dpmdz, dpmdf, dpmdro, &
    dcrdro, dczdro, dcrdf, dczdf
real*8, dimension(4) :: blm, xlm
real*8, dimension(4, 4) :: alm

cnor = 1.d0
iplasm = iplas
nroi = nr
ntetj = nt
ngav = keyctr
tok = platok

epscrz = 5.d-6
epspsm = 5.d-7
epsfpv = 5.d-4
epsro  = 1.d-7

rolim = ro(iplas, jrolim)

if (ngav /= 2) then
    fpv = 0.d0
endif

itout = 0
iter = 0
itin = 0
ich = 0

call f_bndmat

 1000 continue

do
    iter = iter + 1
    itin = itin + 1
    call f_metric
    call f_matcof
    call matpla 
    call extpol
    call rigext 
    call solext 
    call f_matrix
    call f_rightg

    if (i_betp == 1) then
        if (iter > 4) call skbetp(betplx, betpol)
    endif

    call f_solve(0, g)
    call f_rightp
    call f_solve(1, psii)
    call f_psiful
    call artfil
    if (kpr == 1) then
        write(*, *) 'artfil    ', clr, clz
        write(*, *) 'psiax, psim', psiax, psim
        if (isnan(psi(1, 1))) call error_catch
    endif

    if (ngav == 0 .AND. igdf == 2) then
        call qst_b
        call grdef(igdf)
    endif

    call ada(erro)
    if (kpr == 1) then
        write(*, *) 'ada       '
        write(*, *) 'erro=', erro
    endif

! accuracy test parameters

    cab = ABS(clr) + ABS(clz)

    if (ngav/10*10 == ngav) then
        errpsm = 0.d0
        errfpv = 0.d0
    elseif (ngav == 1) then
        errpsm = ABS((psiax - psim)/psipla)
        errpsb = ABS((psip - psibon)/psipla)
        fvv = sqrt(f(iplas)**2 + fpv)
        errfpv = 0.d0
    elseif (ngav == 2) then
        errpsm = ABS((psiax - psim)/psipla)
        errpsb = ABS((psip - psibon)/psipla)
        fvv = sqrt(f(iplas)**2 + fpv)
        errfpv = ABS((fvac - fvv)/(f(1) - f(iplas)))
    endif

    if (erro < epsro .OR. itin > itrmax) then
        if (ich /= 0) EXIT
        if (kpr == 1) then
            write(*, *) 'itin', itin
            write(*, *) 'crz', cab
            write(*, *) 'errpsm', errpsm
            write(*, *) 'errpsb', errpsb
            write(*, *) 'errfpv', errfpv
        endif
        if ( (cab < epscrz) .AND. (errpsm < epspsm) .AND. &
             (errfpv < epsfpv) ) goto 3000
        if ( itout > Nitmax ) goto 3000
        EXIT
    endif
enddo

 2000 continue
itin = 0

if (ich == 0) then
    cr0 = clr
    cz0 = clz
    pm0 = psim
    fv0 = f(iplas)**2 + fpv
    rm0 = rm
    zm0 = zm
    rolim0 = rolim
    fpv0 = fpv
    ich = 1
    drm = 0.101*(ro(2, 1) - ro(1, 1))
    zm = zm + drm

    call reform
    goto 1000
elseif (ich == 1) then
    cr1 = clr
    cz1 = clz
    pm1 = psim
    fv1 = f(iplas)**2 + fpv
    ich = 2
    rm = rm0 + drm
    zm = zm0
    rolim = rolim0
    fpv = fpv0

    call reform
    goto 1000

elseif (ich == 2 .AND. ngav/10*10 == ngav) then

    cr2 = clr
    cz2 = clz
    dcrdr = (cr2 - cr0)/drm
    dczdr = (cz2 - cz0)/drm
    dcrdz = (cr1 - cr0)/drm
    dczdz = (cz1 - cz0)/drm
    det = dcrdr*dczdz - dczdr*dcrdz
    delr = (cz0*dcrdz - cr0*dczdz)/det
    delz = (cr0*dczdr - cz0*dcrdr)/det

    dell = sqrt(delr**2 + delz**2)
    if (dell > drm*10.d0) then
        nshift = dell/(drm*10.d0)
        delr = drm*10.d0*delr/dell
        delz = drm*10.d0*delz/dell
        if (nshift > 10) nshift = 10
        do ish=1, nshift
            rm = rm0 + delr*ish
            zm = zm0 + delz*ish
            call reform
            call f_metric
            call extpol
            call f_matcof
            call matpla
            call rigext
            call solext
            call f_matrix
            call f_rightg
            call f_solve(0, g)
            call f_rightp
            call f_solve(1, psii)
            call f_psiful
            call artfil
            call ada(erro)
        enddo
    else
        rm = rm0 + delr
        zm = zm0 + delz
        call reform
    endif
    ich = 0
    itout = itout+1
    goto 1000

elseif (ich == 2 .AND. ngav > 0) then

    cr2 = clr
    cz2 = clz
    pm2 = psim
    fv2 = f(iplas)**2 + fpv
    ich = 3
    rm = rm0
    zm = zm0
    drolim = -drm
    rolim = rolim0 + drolim
    fpv = fpv0
    call reform
    goto 1000

elseif (ich == 3 .AnD. ngav == 1) then

    cr3 = clr
    cz3 = clz
    pm3 = psim
    dcrdr = (cr2 - cr0)/drm
    dczdr = (cz2 - cz0)/drm
    dpmdr = (pm2 - pm0)/drm
    dcrdz = (cr1 - cr0)/drm
    dczdz = (cz1 - cz0)/drm
    dpmdz = (pm1 - pm0)/drm
    dcrdro = (cr3 - cr0)/drolim
    dczdro = (cz3 - cz0)/drolim
    dpmdro = (pm3 - pm0)/drolim

    alm(1, 1) = dcrdr
    alm(1, 2) = dcrdz
    alm(1, 3) = dcrdro
    alm(2, 1) = dczdr
    alm(2, 2) = dczdz
    alm(2, 3) = dczdro
    alm(3, 1) = dpmdr
    alm(3, 2) = dpmdz
    alm(3, 3) = dpmdro

    blm(1) = -cr0
    blm(2) = -cz0
    blm(3) =  psiax - pm0

    call ge(3, 4, alm, blm, xlm, iwrk)

    delr = xlm(1)
    delz = xlm(2)
    delro = xlm(3)
    rm = rm0 + delr
    zm = zm0 + delz
    rolim = rolim0 + delro

    call reform

    ich = 0
    itout = itout + 1
    goto 1000

elseif (ich == 3 .AnD. ngav > 1) then

    cr3 = clr
    cz3 = clz
    pm3 = psim
    fv3 = f(iplas)**2 + fpv
    ich = 4
    rm = rm0
    zm = zm0
    rolim = rolim0
    dfpv = (f(1)**2 - fvac**2)*1.d-3
    fpv = fpv0 + dfpv

    call reform
    goto 1000

elseif (ich == 4) then

    cr4 = clr
    cz4 = clz
    pm4 = psim
    fv4 = f(iplas)**2 + fpv

    dcrdr = (cr2 - cr0)/drm
    dczdr = (cz2 - cz0)/drm
    dpmdr = (pm2 - pm0)/drm
    dfvdr = (fv2 - fv0)/drm
    dcrdz = (cr1 - cr0)/drm
    dczdz = (cz1 - cz0)/drm
    dpmdz = (pm1 - pm0)/drm
    dfvdz = (fv1 - fv0)/drm
    dcrdro = (cr3 - cr0)/drolim
    dczdro = (cz3 - cz0)/drolim
    dpmdro = (pm3 - pm0)/drolim
    dfvdro = (fv3 - fv0)/drolim
    dcrdf = (cr4 - cr0)/dfpv
    dczdf = (cz4 - cz0)/dfpv
    dpmdf = (pm4 - pm0)/dfpv
    dfvdf = (fv4 - fv0)/dfpv

    alm(1, 1) = dcrdr
    alm(1, 2) = dcrdz
    alm(1, 3) = dcrdro
    alm(1, 4) = dcrdf

    alm(2, 1) = dczdr
    alm(2, 2) = dczdz
    alm(2, 3) = dczdro
    alm(2, 4) = dczdf

    alm(3, 1) = dpmdr
    alm(3, 2) = dpmdz
    alm(3, 3) = dpmdro
    alm(3, 4) = dpmdf

    alm(4, 1) = dfvdr
    alm(4, 2) = dfvdz
    alm(4, 3) = dfvdro
    alm(4, 4) = dfvdf

    blm(1) = -cr0
    blm(2) = -cz0
    blm(3) =  psiax-pm0
    blm(4) =  fvac**2-fv0

    call ge(4, 4, alm, blm, xlm, iwrk)

    delr = xlm(1)
    delz = xlm(2)
    delro = xlm(3)
    delfpv = xlm(4)

    rm = rm0 + delr
    zm = zm0 + delz
    rolim = rolim0 + delro
    fpv = fpv0 + delfpv

    call reform
 
    ich = 0
    itout = itout + 1
    goto 1000

endif

 3000 continue

nflag = 0

if (itout > nitmax) nflag = 1
errarr(1) = erro    !  out of accuracy tests
errarr(2) = cab
errarr(3) = errpsm
errarr(4) = errfpv

if (itout < 1) then
    go to 2000
endif

if (ngav == 0) then
    call qst_b
    psiax = psim
    psipla = psim-psip
    Fvac  =  f(iplas)
endif

if (kpr == 1) then
    write(*, *) 'iter itout', iter, itout
endif

call flux_r(psitok, ncequi)

platok = tokp
rax    = rm
zax    = zm
rxpnt  = rx0
zxpnt  = zx0
psax   = psim
psbo   = psip
psdel  = psim - psip

call bt_pol(betpol)
call f_wrd

iter = 0
itin = 0

return
end subroutine eqa

!---------------------------------------------------------------------
subroutine eqa_ax(dt, time, keyctr, igdf, nstep, platok, psax, i_betp, betplx, &
    rax, zax, rxpnt, zxpnt, psbo, psdel, pcequi, ncequi, psitok, &
    betpol, nflag, errarr)

use keys, only: kpr, kstep
use tim, only: dtim, ctim
use compol_add, only: itrmax, nitmax, g, psii, rx0, zx0
use compol, only: ngav, nitdel, nitbeg, iter, itin, &
    rm, zm, cnor, tok, tokp, erru, &
    psiax, psim, psip, psipla

implicit none

integer, intent(in) :: keyctr, igdf, nstep, i_betp, ncequi
real*8, intent(in) :: dt, time, betplx, betpol
real*8, intent(in), dimension(*) :: pcequi, psitok

integer, intent(out) :: nflag
real*8, intent(out) :: platok, rax, zax, rxpnt, zxpnt, psax, psbo, psdel
real*8, intent(out), dimension(*) :: errarr

integer :: j_ip_count, max_max_iterj, nstepo
real*8 :: erro, errpsm
double precision :: plappok

data j_ip_count/0/
data plappok/0.1/

save j_ip_count, plappok, nstepO, max_max_iterj

ngav = keyctr
kstep = nstep
dtim = dt
ctim = time

itrmax = 50
Nitmax = 5
nitdel = 7
nitbeg = 5
tok = platok

if (kpr == -2) then
    if (j_ip_count == max_max_iterj) then
        tok = platok
        plappok = 0
        j_ip_count = 0
    else
        j_ip_count = j_ip_count + 1
        if (j_ip_count >= 1) then
            plappok = plappok + platok
        endif
        tok = plappok/j_ip_count
    endif
    write(*, *) 'plap', platok, tok, plappok
else
    plappok = 0   
endif

cnor = 1.d0

call f_ext_fil(pcequi, ncequi)

if (nstep /= nstepO) then 
    itin = 0
    erru = 0.d0
    write(*, *) 'renet', platok !EFable
    call renet
    call f_bndmat
endif

iter =iter + 1
itin =itin + 1

if (kpr == 1) then
    write(*, *) ' '
    write(*, *) 'iter=', iter, itin
endif

call f_metric
call f_matcof
call matpla
call extpol
call rigext
call solext
call f_matrix
call f_rightg

if (i_betp == 1) then
    if (iter > 4) call skbetp(betplx, betpol)
endif

call f_solve(0, g)
call f_rightp
call f_solve(1, psii)
call f_psiful

if (ngav <= 0 .AND. igdf == 2) then
    call qst_b
    call grdef(igdf)
endif

call f_remesh(erro)
erru = erro

if (ngav == 0) then
    errpsm = 0.d0
elseif (ngav == 1) then
    errpsm = ABS((psiax - psim)/psipla)
endif

call flux_r(psitok, ncequi)
if (kpr == -2) then
    tok  = platok
    tokp = platok
endif
platok = tokp
rax    = rm
zax    = zm
rxpnt  = rx0
zxpnt  = zx0
psax   = psim
psbo   = psip
psdel  = psim-psip

if (ngav <= 0) then
    psiax  = psim
    psipla = psim - psip
endif
call bt_pol(betpol)

nflag = 0
errarr(1) = erro    !  out of accuracy tests
nstepO = nstep

if (ngav < 0) call retab_L

return
end subroutine eqa_ax
