subroutine f_bndmat(rk, zk, nk, rlop, zlop, nlop, rprob, zprob, nprob)

use compol_add, only: binadg
use compol, only: nr, nt1, r, z
      
implicit none

integer, intent(in) :: nk, nlop, nprob
real*8, intent(in), dimension(*) :: rk, zk, rlop, zlop, rprob, zprob

integer :: j, jb
real*8 :: r0, r1, z0, z1, rr, zz, fint

do j=2,nt1
    rr = r(nr, j)
    zz = z(nr, j)
    do jb=2, nt1
        r0 = r(nr, jb)
        z0 = z(nr, jb)
        r1 = r(nr, jb+1)
        z1 = z(nr, jb+1)
        call bint(rr, zz, R0, Z0, r1, z1, Fint, 1)
        binadg(jb, j) = fint
    enddo
enddo

return
end subroutine f_bndmat
