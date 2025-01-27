subroutine eqa_in(alf0, alf1, alf2, bet0, bet1, bet2, platok, &
    pcequi, ncequi, b0cen, r0cen)

use parameters_a2equil, only: fix_adapgrid
use sp_parameters, only: nursp
use iopath, only: path
use compol_add, only: alp, nctrl, ron_max, ron_max_g, itrmax, &
    nitmax, jrolim, nblm, rblm, zblm
use compol, only: nr, nr1, nr2, nt, nt1, nt2, iplas, iplas1, &
     nitdel, nitbeg, r0ax, b0ax, tok, cnor, fvac

implicit none

integer, parameter :: nursp4=nursp+4, nursp6=nursp4*6

integer, intent(in) :: ncequi
real*8, intent(in) :: alf0, alf1, alf2, bet0, bet1, bet2, platok, b0cen, r0cen
real*8, intent(in), dimension(*) :: pcequi

logical :: file_existence
integer :: i, i_vac, k_dummy, nlimadap, nroi, ntetj, iplasm
real*8 :: alf0p, alf1p, alf2p, bet0f, bet1f, bet2f
real*8, dimension(500) :: rlimadap, zlimadap
character(len=80) :: fname

save rlimadap, zlimadap, nlimadap

alf0p = alf0
alf1p = alf1
alf2p = alf2
bet0f = bet0
bet1f = bet1
bet2f = bet2
b0ax = b0cen
r0ax = r0cen
tok = platok

!...input initial data

write(fname, '(a, a)') TRIM(path), '/egg.dat'
open(1, file=fname)
read(1, *) i_vac
read(1, *) alp
read(1, *) k_dummy  !(is not used)
read(1, *) nctrl
read(1, *) ron_max
close(1)

ron_max_g = ron_max

if (fix_adapgrid == 1) then		
    write(fname, '(a, a)') TRIM(path), '/egg_g.dat'
    write(*, *) fname
    INQUIRE( FILE=trim(fname), EXIST= file_existence)
    if (file_existence) then		
        write(*, *) 'reading'
        open(53, FILE=fname)
        read(53, *) nlimadap
        do i=1, nlimadap
            read(53, *) rlimadap(i), zlimadap(i)
        enddo
        close(53)
        write(*, *) 'closing'	
    endif
endif

Nr = iplas + i_vac
itrmax = 100
Nitmax = 5
nitdel = 7
nitbeg = 5
cnor   = 1.d0
jrolim = 2  !+(nt-2)/2

if (nctrl == 1) then
    write(fname, '(a, a)') TRIM(path), '/limpnt_d.wr'
    open(1, file=fname)
    read(1, *) nblm
    do i=1, nblm
        read(1, *) rblm(i), zblm(i)
    enddo
    close(1)
endif

call rdrec
call grid_spdr
call f_wrd

fvac = b0cen*r0cen

nr1 = nr - 1
nt1 = nt - 1
nr2 = nr - 2
nt2 = nt - 2
nroi  = nr
ntetj = nt

iplas1 = iplas - 1
iplasm = iplas

call f_rdexf(ncequi)
call f_ext_fil(pcequi, ncequi)

return
end subroutine eqa_in
