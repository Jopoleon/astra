subroutine get_nrho_ntheta(nrho_surf, nthe_surf)

use const_inc, only: NA1, MEQUIL, NEQUIL, LEQ

implicit none

integer, intent(out) :: nrho_surf, nthe_surf

if (LEQ(5) <= 2) then
    nrho_surf = NA1
    nthe_surf = 41
else
! SPIDER, FEQIS
    nrho_surf = abs(nint(NEQUIL)) !neqlp
    nthe_surf = abs(nint(MEQUIL)) !ntetap
endif

return
end subroutine get_nrho_ntheta

!---------------------------------------------------------------------
subroutine SURF_CTR(nrho_surf, nthe_surf, r_surf, z_surf)

use const_inc, only: GP2, RTOR, UPDWN, LEQ
use status_inc, only: SHIF, AMETR, ELON, TRIA
use parameters_a2equil, only: equil_now

implicit none

integer, intent(in) :: nrho_surf, nthe_surf
double precision, dimension(nrho_surf, nthe_surf), intent(out) :: r_surf, z_surf
integer :: jthe, jrho
double precision :: theta

! Surface contours

if (LEQ(5) <= 2) then
! 3 moments equilibrium
    do jthe = 1, nthe_surf
        theta = GP2/DBLE(nthe_surf-1)*(jthe - 1)
        do jrho=1, nrho_surf
            r_surf(jrho, jthe) = RTOR + SHIF(jrho) + AMETR(jrho) * &
                  ( COS(theta) + 0.5*TRIA(jrho) * (COS(2.*theta) - 1.))
            z_surf(jrho, jthe) = UPDWN + AMETR(jrho)*ELON(jrho)*SIN(theta)
        enddo
    enddo
else
! SPIDER, FEQIS
    r_surf = equil_now%coord_sys%position%r
    z_surf = equil_now%coord_sys%position%z
endif

return
end subroutine SURF_CTR
