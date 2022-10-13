      subroutine matrix

      use sp_parameters, only: nip
      use comblc, only: ni1, nj1, r, dr, r12, dri, dz

      implicit none

      integer :: i, j
      real*8 :: ri, rpl, rmn, drpl, drmn, dzpl, dzmn
      real*8, dimension(nip) :: axm, ax0, axp
      common /comffs/ axm, ax0, axp

!matrix pakage for fourier solver

      do i=2,ni1
         ri=r(i)
         drpl=dr(i)
         drmn=dr(i-1)
         rpl=r12(i)
         rmn=r12(i-1)
         do j=2,nj1
            dzpl=dz(j)
            dzmn=dz(j-1)
            axm(i)=ri/(drmn*rmn*dri(i))
            axp(i)=ri/(drpl*rpl*dri(i))
            ax0(i)=axm(i)+axp(i)
         enddo
      enddo

      return
      end subroutine matrix

!----------------------------------------------------------------
      integer function nlin(i, j)

      use comblc, only: ni2

      implicit none

      integer, intent(in) :: i, j

      nlin=(j-2)*ni2+(i-1)

      return
      end function nlin
