      integer :: nblm, nctrl, numlim, ich, itrmax, Nitmax, 
     &   ixp1, jxp1, ixp2, jxp2, nk_out,
     &   nlop_out, nprob_out, jrolim, iterbf
      integer, dimension(nblmp) :: iblm, jblm
      integer, dimension(nkp) :: iprcon
      integer, dimension(nloopp) :: iprlop
      integer, dimension(nprobp) :: iprprob

      real*8 :: ublmax, alp, alpnew, pscen, psix0, psix1, psxi2, 
     &    rx0, zx0, rx1, zx1, rx2, zx2, rl, zl, clr, clz,
     &    rolim, ron_max, fpv
      real*8, dimension(500) :: ron_max_g
      real*8, dimension(nblmp) :: rblm, zblm
      real*8, dimension(neqp) :: zpro
      real*8, dimension(ntp) :: dgdn
      real*8, dimension(ntp, ntp) :: binadg
      real*8, dimension(nrp, ntp) :: psie, psii, g, aex
      real*8, dimension(nkp, ntp) :: pinadg
      real*8, dimension(nloopp, ntp) :: adginl
      real*8, dimension(nprobp, ntp) :: adginr, adginz

      common /compol_add/ rblm, zblm, ublmax, 
     & iblm, jblm, nblm, nctrl, numlim,
     & alp, alpnew, pscen,
     & zpro,
     & ich, itrmax, Nitmax,
     % psie, psii, g, 
     & psix0, psix1, psxi2, 
     & rx0, zx0, rx1, zx1, rx2, zx2, 
     & ixp1, jxp1, ixp2, jxp2,
     & dgdn, binadg, pinadg, adginl, adginr, adginz,
     & iprcon , nk_out, iprlop, nlop_out, 
     & iprprob, nprob_out,
     & rl, zl, clr, clz, 
     & rolim, jrolim, iterbf, ron_max, ron_max_g,
     & aex, fpv

!      common /comblm/ rblm, zblm, ublmax, 
!     &                iblm, jblm, nblm, nctrl, numlim
!      common /comcon_add/ alp, alpnew, pscen
!      common /commat_add/ zpro
!      common /comitr_add/ ich, itrmax, Nitmax
!      common /compot_add/ psie, psii, g, 
!     +                   psix0, psix1, psxi2, 
!     +                   rx0, zx0, rx1, zx1, rx2, zx2, 
!     +                   ixp1, jxp1, ixp2, jxp2
!      common /combnd_add/ dgdn, binadg, pinadg, adginl, adginr, adginz,
!     +                    iprcon , nk_out, iprlop, nlop_out, 
!     +                    iprprob, nprob_out
!      common /comstb_add/ rl, zl, clr, clz, 
!     &   rolim, jrolim, iterbf, ron_max, ron_max_g
!      common /compla_add/ aex, fpv
