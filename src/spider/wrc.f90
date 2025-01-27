SUBROUTINE wrcoil(nk, nkcoil, rk, zk, tk, necon, wecon)

use iopath, only: path

implicit none

integer, intent(in) :: nk, nkcoil, necon(nkcoil)
real*8, intent(in), dimension(nk) :: rk, zk, tk
real*8, intent(in) :: wecon(nkcoil)

integer :: j
character(len=80) :: fname

write(fname, '(a, a)') TRIM(path), '/ecur.wr'
open(1, file=fname, form='formatted')
write(1, *) nk, nkcoil
write(1, *) (rk(j), j=1, nk)
write(1, *) (zk(j), j=1, nk)
write(1, *) (tk(j), j=1, nk)
write(1, *) (wecon(j), j=1, nkcoil)
write(1, *) (necon(j), j=1, nkcoil)
close(1)

return
end subroutine wrcoil 

!---------------------------------------------------------------
SUBROUTINE rdcoil(nk, nkcoil, rk, zk, tk, necon, wecon)

use iopath, only: path

implicit none

integer, intent(out) :: nk, nkcoil, necon(*)
real*8, intent(out) :: rk(*), zk(*), tk(*), wecon(*)

integer :: j
character(len=80) :: fname

write(fname, '(a, a)') TRIM(path), '/ecur.wr'
open(1, file=fname, form='formatted')
read(1, *) nk, nkcoil
read(1, *) (rk(j), j=1, nk)
read(1, *) (zk(j), j=1, nk)
read(1, *) (tk(j), j=1, nk)
read(1, *) (wecon(j), j=1, nkcoil)
read(1, *) (necon(j), j=1, nkcoil)
close(1)

return
end subroutine rdcoil
