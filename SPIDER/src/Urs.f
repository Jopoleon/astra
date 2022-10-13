      real*8 function funppp(psi)

      use sp_parameters, only: nursp

      implicit none

      real*8, intent(in) :: psi

      integer :: nurs, i
      real*8 :: dpsip, dpsim
      real*8, dimension(nursp) :: ppp, fff, www, psit, purs, furs, wurs

      common/comppp/ ppp, fff, www
      common/comurs/psit, purs, furs, wurs, nurs

      if(psi.lt.0.d0) then
         funppp =0.d0
         return
      endif
      if(psi.gt.1.d0) then
         funppp =ppp(nurs)
         return
      endif

      i=1+psi*(nurs-1)

      if(i.eq.nurs) i=nurs-1
      if(i.lt.1) i=1

      dpsip=psit(i+1)-psi
      dpsim=psi-psit(i)

      funppp=(ppp(i)*dpsip+ppp(i+1)*dpsim)/(dpsip+dpsim)

      return
      end function funppp

!----------------------------------------------------------------
      real*8 function funcur(ri, psi)

      use sp_parameters, only: nursp

      implicit none

      real*8, intent(in) :: ri, psi

      integer :: nurs, i
      real*8 :: dpsip, dpsim, pp, fp, wp
      real*8, dimension(nursp) :: psit, purs, furs, wurs

      common/comurs/psit, purs, furs, wurs, nurs

      if(psi.lt.0.d0) then
         funcur =0.d0
         return
      endif

      if(psi.gt.1.d0) then
         funcur =purs(nurs)*ri+furs(nurs)/ri+wurs(nurs)*ri**3
         return
      endif

      i=1+psi*(nurs-1)

      if(i.ge.nurs) i=nurs-1
      if(i.lt.1) i=1

      dpsip=psit(i+1)-psi
      dpsim=psi-psit(i)

      pp=(purs(i)*dpsip+purs(i+1)*dpsim)/(dpsip+dpsim)
      fp=(furs(i)*dpsip+furs(i+1)*dpsim)/(dpsip+dpsim)
      wp=(wurs(i)*dpsip+wurs(i+1)*dpsim)/(dpsip+dpsim)

      funcur = ri*pp + fp/ri + wp*ri**3

      return
      end function funcur

!----------------------------------------------------------------
      subroutine taburs(ien, coin, nursb)
 
      use ppf_modul       
      use sp_parameters, only: nursp
      use keys, only: kastr

      implicit none

      integer, parameter :: nursp4=nursp+4,nursp6=nursp4*6

      integer, intent(in) :: ien, nursb
      real*8, intent(in) :: coin

      integer :: nurs, i, j, IFAIL, key_int
      real*8 :: alf0n, dpsip, dpsim, dpsi, zpsi, funfp, funpp, funwp, 
     &   alf0p, alf1p, alf2p, bet0f, bet1f, bet2f, alw0p, alw1p, alw2p
      real*8 :: CWK(4)
      real*8, dimension(nursp) :: psit, purs, furs, wurs, 
     &   ppp, fff, www
      real*8 RRK(nursp4), CCK(nursp4), WRK(nursp6)

      common/comurs/psit, purs, furs, wurs, nurs
      common/comppp/ ppp, fff, www
      common/comsav/ alf0n
      common /comabw/ alf0p, alf1p, alf2p, bet0f, bet1f, 
     &    bet2f, alw0p, alw1p, alw2p

! table for fast calculation P' and FF' as funct. of PSI_nor
! P' - pptab
! FF'- fptab
! PSI_nor - pstab
 
      if(kastr.eq.1) then 
         alw0p=0.d0
         alw1p=0.d0
         alw2p=0.d0
      endif

      if(ien.eq.0) then
         nurs=3999

! if nursb<0 p'and ff'are assumed to be defined as tab.(file 'tabppf.dat'

         psit(1)=0.d0
         dpsi=1.d0/(nurs-1.d0)

         do i=2,nurs-1
            zpsi=psit(i-1)+dpsi
            psit(i)=zpsi
         enddo

         psit(nurs)=1.d0

         if(nursb.gt.0) then
            do i=1,nurs
               zpsi=psit(i)
               purs(i)=funpp(zpsi)
               furs(i)=funfp(zpsi)
               wurs(i)=funwp(zpsi)
            enddo
         else
	    do i=1,nutab
               pstab(i)=pstab(i)/pstab(nutab)
            enddo

            key_int=0
            if(key_int.eq.0) then
               do i=2,nurs-1
                  zpsi=1.d0-psit(i)
	          do j=1,nutab-1
                     if(zpsi.ge.pstab(j) .AnD. zpsi.lt.pstab(j+1)) then
                        dpsip=pstab(j+1)-zpsi
                        dpsim=zpsi-pstab(j)
                        purs(i)=( pptab(j)*dpsip+pptab(j+1)*dpsim )/
     &                          (dpsip+dpsim)
                        furs(i)=( fptab(j)*dpsip+fptab(j+1)*dpsim )/
     &                          (dpsip+dpsim)
                     endif
                  enddo
               enddo
            else
               CALL E01BAF(Nutab,pstab,pptab,RRK,CCK,
     &                     nutab+4,WRK,6*nutab+16,IFAIL)
               if(ifail.ne.0) write(*,*) 'ifail=',ifail
               do i=2,nurs-1
                  zpsi=1.d0-psit(i)
                  CALL E02BCF(Nutab+4,RRK,CCK,zpsi,0,CWk,IFAIL)
                  if(ifail.ne.0) write(*,*) 'ifail=',ifail
                  purs(i)=cwk(1)
               enddo

               CALL E01BAF(Nutab,pstab,fptab,RRK,CCK,
     &                     nutab+4,WRK,6*nutab+16,IFAIL)
               if(ifail.ne.0) write(*,*) 'ifail=',ifail

               do i=2,nurs-1
                  zpsi=1.d0-psit(i)
                  CALL E02BCF(Nutab+4,RRK,CCK,zpsi,0,CWk,IFAIL)
                  if(ifail.ne.0) write(*,*) 'ifail=',ifail
                  furs(i)=cwk(1)
               enddo
            endif

            purs(1)=pptab(nutab)
            purs(nurs)=pptab(1)
            furs(1)=fptab(nutab)
            furs(nurs)=fptab(1)

            do i=1,nurs
               zpsi=psit(i)
               wurs(i)=funwp(zpsi)
            enddo
         endif

      elseif(ien.eq.1) then

         do i=1,nurs
            purs(i)=purs(i)*coin
         enddo

      endif

      ppp(1)=0.d0
      fff(1)=0.d0
      www(1)=0.d0
      dpsi=1.d0/(nurs-1.d0)

      do i=2,nurs
         ppp(i)=ppp(i-1) +(purs(i-1)+purs(i))*dpsi*0.5d0
         fff(i)=fff(i-1) +(furs(i-1)+furs(i))*dpsi*0.5d0
         www(i)=www(i-1) +(wurs(i-1)+wurs(i))*dpsi*0.5d0
      enddo

      return
      end subroutine taburs

!----------------------------------------------------------------
      subroutine tabnor(cnor)

      use sp_parameters, only: nursp

      implicit none

      integer, parameter :: nursp4=nursp+4,nursp6=nursp4*6

      real*8, intent(in) :: cnor

      integer :: nurs, i
      real*8, dimension(nursp) :: psit, purs, furs, wurs, ppp, fff, www

      common/comurs/psit, purs, furs, wurs, nurs
      common/comppp/ ppp, fff, www

      do i=1,nurs
	 purs(i)=cnor*purs(i)
         furs(i)=cnor*furs(i)
         ppp(i)=cnor*ppp(i)
         fff(i)=cnor*fff(i)
      enddo

      return
      end subroutine tabnor

!----------------------------------------------------------------
      real*8 function funpp(psi)

      use sp_parameters, only: nursp

      implicit none

      real*8, intent(in) :: psi

      real*8 :: epss, psm, zpp
      real*8 :: alf0p, alf1p, alf2p, bet0f, bet1f, bet2f, 
     &    alw0p, alw1p, alw2p
      common /comabw/ alf0p, alf1p, alf2p, bet0f, bet1f, bet2f, 
     &    alw0p, alw1p, alw2p

      epss=1.d-10

      if(psi.lt.0.d0) then
         funpp=0.d0
         return
      endif

      psm = 1.d0-psi + epss

      zpp = dabs( 1.d0 -(dabs(psm))**alf1p ) + epss
      funpp  = alf0p*( zpp**alf2p )

      return
      end function funpp

!----------------------------------------------------------------
      real*8 function funfp(psi)

      use sp_parameters, only: nursp

      implicit none

      real*8, intent(in) :: psi

      real*8 :: alf0p, alf1p, alf2p, bet0f, bet1f, bet2f, 
     &    alw0p, alw1p, alw2p
      real*8 :: epss, zzfp, psm
      common /comabw/ alf0p, alf1p, alf2p, bet0f, bet1f, bet2f, 
     &    alw0p, alw1p, alw2p

      epss=1.d-10

      if(psi.lt.0.d0) then
         funfp=0.d0
         return
      endif

      psm = 1.d0- psi + epss

      zzfp = dabs( 1.d0 - psm**bet1f ) + epss
      funfp   = bet0f*( zzfp**bet2f )

      return
      end function funfp

!----------------------------------------------------------------
      real*8 function tabw(psi)

      use sp_parameters, only: nursp

      implicit none

      real*8, intent(in) :: psi

      integer :: nurs, i
      real*8 :: wp, dpsip, dpsim
      real*8, dimension(nursp) :: psit, purs, furs, wurs

      common/comurs/psit, purs, furs, wurs, nurs

      if(psi.lt.0.d0) then
         wp =0.d0
         return
      endif

      if(psi.gt.1.d0) then
         wp =wurs(nurs)
         return
      endif

      i=1+psi*(nurs-1)

      if(i.eq.nurs) i=nurs-1
      if(i.lt.1) i=1

      dpsip=psit(i+1)-psi
      dpsim=psi-psit(i)

      wp=(wurs(i)*dpsip+wurs(i+1)*dpsim)/(dpsip+dpsim)

      tabw = wp

      return
      end function tabw

!----------------------------------------------------------------
      real*8 function tabp(psi)

      use sp_parameters, only: nursp

      implicit none

      real*8, intent(in) :: psi

      integer :: nurs, i
      real*8 :: pp, dpsip, dpsim
      real*8, dimension(nursp) :: psit, purs, furs, wurs

      common/comurs/psit, purs, furs, wurs, nurs

      if(psi.lt.0.d0) then
         pp =0.d0
         return
      endif

      if(psi.gt.1.d0) then
         pp =purs(nurs)
         return
      endif

      i=1+psi*(nurs-1)

      if(i.eq.nurs) i=nurs-1
      if(i.lt.1) i=1

      dpsip=psit(i+1)-psi
      dpsim=psi-psit(i)

      pp=(purs(i)*dpsip+purs(i+1)*dpsim)/(dpsip+dpsim)

      tabp = pp

      return
      end function tabp

!----------------------------------------------------------------
      real*8 function tabf(psi)

      use sp_parameters, only: nursp
      implicit none

      real*8, intent(in) :: psi

      integer :: nurs, i
      real*8 :: fp, dpsip, dpsim
      real*8, dimension(nursp) :: psit, purs, furs, wurs

      common/comurs/psit, purs, furs, wurs, nurs

      if(psi.lt.0.d0) then
         fp =0.d0
         return
      endif

      if(psi.gt.1.d0) then
         fp =furs(nurs)
         return
      endif

      i=1+psi*(nurs-1)

      if(i.eq.nurs) i=nurs-1
      if(i.lt.1) i=1

      dpsip=psit(i+1)-psi
      dpsim=psi-psit(i)

      fp=(furs(i)*dpsip+furs(i+1)*dpsim)/(dpsip+dpsim)

      tabf = fp

      return
      end function tabf

!----------------------------------------------------------------
      real*8 function funwp(psi)

      use sp_parameters, only: pi, nursp

      implicit none

      real*8, intent(in) :: psi

      real*8 :: wp, psm, epss, zwp
      real*8 :: alf0p, alf1p, alf2p, bet0f, bet1f, 
     &    bet2f, alw0p, alw1p, alw2p

      common /comabw/ alf0p, alf1p, alf2p, bet0f, bet1f, 
     &    bet2f, alw0p, alw1p, alw2p

      epss=1.d-10

      if(psi.lt.0.d0) then
         funwp=0.d0
         return
      endif

      psm = 1.d0 - psi + epss

      zwp = dabs( 1.d0 -(dabs(psm))**alw1p ) + epss
      wp  = alw0p*( zwp**alw2p )

      funwp = wp

      return
      end function funwp
