subroutine SURF_CTR(order, nrho_surf, nthe_surf, r_surf, z_surf)

use const_inc, only: NA1, GP2, RTOR, UPDWN, MEQUIL, NEQUIL, LEQ
use status_inc, only: SHIF, AMETR, ELON, TRIA
use parameters_a2spider, only: equil_now

implicit none

integer, parameter :: n_surf=556
logical, intent(in) :: order
integer, intent(out) :: nrho_surf, nthe_surf
double precision, dimension(n_surf, n_surf), intent(out) :: r_surf, z_surf
integer :: jthe, jrho
double precision :: theta, rtmp, ztmp

! Surface contours

if (MEQUIL == 0. .or. LEQ(5) <= 1) then

! 3 moments equilibrium
   nrho_surf = NA1
   nthe_surf = 41
   do jthe = 1, nthe_surf
      theta = GP2/DBLE(nthe_surf-1)*(jthe - 1)
      do jrho=1, nrho_surf
         rtmp = RTOR + SHIF(jrho) + AMETR(jrho) * &
                  ( COS(theta) + 0.5*TRIA(jrho) * (COS(2.*theta) - 1.))
         ztmp = UPDWN + AMETR(jrho)*ELON(jrho)*SIN(theta)
         if(order) then
            r_surf(jrho, jthe) = rtmp
            z_surf(jrho, jthe) = ztmp
         else
            r_surf(jthe, jrho) = rtmp
            z_surf(jthe, jrho) = ztmp
         endif
      enddo
   enddo

else

! SPIDER
   nrho_surf = abs(nint(NEQUIL)) !neqlp
   nthe_surf = abs(nint(MEQUIL)) !ntetap
   do jrho=1, nrho_surf
      if(order) then
          r_surf(jrho, 1:nthe_surf) = equil_now%coord_sys%position%r(jrho, 1:nthe_surf)
          z_surf(jrho, 1:nthe_surf) = equil_now%coord_sys%position%z(jrho, 1:nthe_surf)
      else
          r_surf(1:nthe_surf, jrho) = equil_now%coord_sys%position%r(jrho, 1:nthe_surf)
          z_surf(1:nthe_surf, jrho) = equil_now%coord_sys%position%z(jrho, 1:nthe_surf)
      endif
   enddo
endif

end subroutine SURF_CTR
