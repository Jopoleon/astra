      subroutine remesh(erro, errpsi, imov)

      use keys, only: kpr
      use compol, only: psi, iplas1, nt1, iswtch, iter, ngav

      implicit none

      integer, intent(in) :: imov
      real*8, intent(out) :: erro, errpsi

      integer :: i, j, imax, jmax
      real*8 :: psimax

      psimax=psi(1,2)
      imax=1
      jmax=2

      do i=2,iplas1
         do j=2,nt1
            if(psi(i,j).gt.psimax) then
               psimax=psi(i,j)
               imax=i
               jmax=j
            endif
         enddo
      enddo

      iswtch=0
	     
      if(iter.ne.1) then
         if(imax.eq.1 .and. erro.lt.5.d-4) iswtch=1
         if(kpr.eq.1) then
            write(*,*) 'mesh:imax jmax iswth',imax,jmax,iswtch
         endif
      else
         iswtch=0       
      endif

      if(ngav.le.0) iswtch=0
      if(imax.ge.2) then
         call prgrid(imax,jmax,erro,imov,errpsi)
      else
         if(iswtch.eq.0) then
            call regrid0(erro,imov,errpsi)
         else
            call regrid(erro,imov,errpsi)
         endif
      endif

      return
      end subroutine remesh

!----------------------------------------------------------------
      subroutine prgrid(imax, jmax, erro, imov, errpsi)

      use sp_parameters, only: ntp, twopi
      use keys, only: kpr
      use compol, only: r, z, rm, zm, nt1, psi, psia, psim, psin, psip, 
     & teta, iplas, iplas1, ro, nr, nt

      implicit none

      integer, parameter :: nshp=ntp+1

      integer, intent(in) :: imax, jmax, imov
      real*8, intent(out) :: erro, errpsi

      integer :: i, j, k, l, nsh, ierm, jerm, ibeg, is, isc, jj
      real*8 :: rmx, zmx, rm0, zm0, drx, dzx, zps,
     &   psim0, psinm0, det, psnn, delpsn, robj, tetp, ps0, pspl, ro0,
     &   ropl, ro2j, rodeb, grad, zro, rnj, znj, tetv, tt0, ttpl, ttj,
     &   ros1, ros2, rosi, psa1, psa2, psai
      real*8, dimension(5) :: dp
      real*8, dimension(nshp) :: xs, ys, fun
      real*8, dimension(ntp) :: ron, rop, tetn, rob, teti, roi

! New position of magnetic axis

      if(kpr.eq.1) then
         write(*,*) 'prgrid:::'
      endif
      rm0=r(1,1)
      zm0=z(1,1)
      psim0=psi(1,1)

      nsh=1
      xs(nsh)=r(imax,jmax)
      ys(nsh)=z(imax,jmax)
      fun(nsh)=psi(imax,jmax)

      if(imax.eq.1) then
         do j=2,nt1
            nsh=nsh+1
            xs(nsh)=r(2,j)
            ys(nsh)=z(2,j)
            fun(nsh)=psi(2,j)
         enddo
      elseif(imax.eq.2) then
         nsh=nsh+1
         xs(nsh)=r(1,2)
         ys(nsh)=z(1,2)
         fun(nsh)=psi(1,2)
         do j=2,nt1
            nsh=nsh+1
            xs(nsh)=r(3,j)
            ys(nsh)=z(3,j)
            fun(nsh)=psi(3,j)
         enddo
      else
         do k=-1,1
            i= imax+k
            do l=-1,1
               j= jmax+l
               if(i.ne.imax .OR. j.ne.jmax) then
                  nsh=nsh+1
                  xs(nsh)=r(i,j)
                  ys(nsh)=z(i,j)
                  fun(nsh)=psi(i,j)
	       endif
            enddo
         enddo  
      endif

      call deriv5(xs,ys,fun,nsh,5,dp)

      DET = dp(3)*dp(5) - dp(4)**2
      Rm = Xs(1) + ( dp(2)*dp(4) - dp(1)*dp(5) )/DET
      Zm = Ys(1) + ( dp(1)*dp(4) - dp(2)*dp(3) )/DET

      erro=SQRT( (rm-rm0)**2+(zm-zm0)**2 )

      rmx=r(imax,jmax)
      zmx=z(imax,jmax)

      psim=fun(1)+ dp(1)*(rm-rmx) + dp(2)*(zm-zmx) +
     &       0.5d0*dp(3)*(rm-rmx)*(rm-rmx) +
     &             dp(4)*(rm-rmx)*(zm-zmx) +
     &       0.5d0*dp(5)*(zm-zmx)*(zm-zmx)

      if(kpr.eq.1) then
         write(6,*) 'rm zm psim',rm,zm,psim
         write(6,*) 'prgrid:errma',erro
      endif

! Nomalized poloidal flux definition

      errpsi=0.d0

      do i=1,iplas
         do j=1,nt
            psnn=psin(i,j)
            psin(i,j)=(psi(i,j)-psip)/(psim-psip)
            delpsn=ABS(psin(i,j)-psnn)
            if(delpsn.gt.errpsi) then
               ierm=i
               jerm=j
               errpsi=delpsn
            endif
         enddo
      enddo
      if(imov.eq.0) return

! Definition of new angle grid teta(j)
      do j=1,nt
         tetn(j)=teta(j)
      enddo

      do j=1,nt
         drx=r(iplas,j)-rm
         dzx=z(iplas,j)-zm
         robj=SQRT(drx**2+dzx**2)
         rob(j)=robj
         tetp=ACOS(drx/robj)
         if(dzx.lt.0.d0) then
            teta(j)=-tetp
         else
            teta(j)=tetp
         endif
      enddo

      do j=2,nt
         if(teta(j).lt.teta(j-1)) then
            teta(j)=teta(j)+twopi
         endif
         if(teta(j).lt.teta(j-1)) then
            teta(j)=teta(j)+twopi
         endif
      enddo

      teta(1)=teta(nt1)-twopi
      teta(nt)=teta(2)+twopi

! grid moving along rais

      psinm0=psin(1,2)

      do i=3,iplas1
         if(psia(i).lt.psinm0) then
            ibeg=i+1
	    EXIT
	 endif
      enddo

      do i=ibeg,iplas-1
         zps =psia(i)
         do j=1,nt
            do is=1,iplas1
               ps0 =psin(is,j)
               pspl=psin(is+1,j)
               ro0 =ro(is,j)
               ropl=ro(is+1,j)
               if(zps.le.ps0 .and. zps.ge.pspl) then
                  isc=is
		  EXIT
               endif
            enddo
            grad=-(pspl-ps0)/(ropl-ro0)
	    zro=ro0-(psia(i)-ps0)/grad
            ron(j)=zro
         enddo
         do j=1,nt
            rnj=rm0+ron(j)*COS(tetn(j))
            znj=zm0+ron(j)*SIN(tetn(j))
            drx=rnj-rm
            dzx=znj-zm
            roi(j)=SQRT(drx**2+dzx**2)
            tetp=ACOS(drx/roi(j))
            if(dzx.lt.0.d0) then
               teti(j)=-tetp
            else
               teti(j)=tetp
            endif
         enddo
         do j=2,nt
            if(teti(j).lt.teti(j-1)) then
               teti(j)=teti(j)+twopi
            endif
            if(teti(j).lt.teti(j-1)) then
               teti(j)=teti(j)+twopi
            endif
         enddo

         teti(1)=teti(nt1)-twopi
         teti(nt)=teti(2)+twopi
         roi(1)=roi(nt1)
         roi(nt)=roi(2)

         do j=2,nt1
            tetv=teta(j)
            if(tetv.lt.teti(1)) tetv=tetv+twopi
            if(tetv.lt.teti(1)) tetv=tetv+twopi
            if(tetv.gt.teti(nt)) tetv=tetv-twopi
            if(tetv.gt.teti(nt)) tetv=tetv-twopi
            do jj=1,nt1
               tt0 =teti(jj)
               ttpl=teti(jj+1)
               ro0 =roi(jj)
               ropl=roi(jj+1)
	       if(tetv.le.ttpl .and. tetv.ge.tt0) then
                  zro=((ttpl-tetv)*ro0+(tetv-tt0)*ropl)/(ttpl-tt0)
                  EXIT
	       endif
            enddo

            rop(j)=zro
            r(i,j)=rop(j)*COS(teta(j))+rm
            z(i,j)=rop(j)*SIN(teta(j))+zm
         enddo

      enddo  ! i-loop

!-----------------------------------------
           
      i=2
      zps =psia(i)
      ps0=psip+zps*(psim-psip)

      do j=1,nt
         ttj=teta(j)
         ro2j=(ps0-psim)/( 0.5d0*dp(3)*COS(ttj)**2 +
     &                           dp(4)*COS(ttj)*SIN(ttj) +
     &                     0.5d0*dp(5)*SIN(ttj)**2 )
         rodeb=ro2j
         rop(j)=SQRT(ro2j)
         r(i,j)=rm+rop(j)*COS(teta(j))
         z(i,j)=zm+rop(j)*SIN(teta(j))
      enddo

      if(ibeg.gt.3) then
         do j=2,nt1
            ros1=(r(2,j)-rm)**2+(z(2,j)-zm)**2
            ros2=(r(ibeg,j)-rm)**2+(z(ibeg,j)-zm)**2
            psa1=psia(2)
            psa2=psia(ibeg)
            do i=3,ibeg-1
               psai=psia(i)
               rosi=( (psa2-psai)*ros1+(psai-psa1)*ros2 )/(psa2-psa1)
               rosi=SQRT(rosi)
               r(i,j)=rm+rosi*COS(teta(j))
               z(i,j)=zm+rosi*SIN(teta(j))
            enddo
         enddo
      endif

      do i=2,iplas1
         do j=2,nt1
            ro(i,j)=SQRT((r(i,j)-rm)**2+(z(i,j)-zm)**2)
         enddo
      enddo

      do j=2,nt1
         ro(iplas,j)=rob(j)
         r(iplas,j)=rob(j)*COS(teta(j))+rm
         z(iplas,j)=rob(j)*SIN(teta(j))+zm
         ro(1,j)=0.d0
         r(1,j)=rm
         z(1,j)=zm
      enddo

      do i=1,iplas
         ro(i,1)=ro(i,nt1)
         ro(i,nt)=ro(i,2)
         r(i,1)=r(i,nt1)
         r(i,nt)=r(i,2)
         z(i,1)=z(i,nt1)
         z(i,nt)=z(i,2)
      enddo

      do i=1,nr
         do j=1,nt
            psin(i,j)=psia(i)
            psi(i,j)=psip+psin(i,j)*(psim-psip)
         enddo
      enddo

      call wrb

      return
      end subroutine prgrid

!----------------------------------------------------------------
      subroutine regrid(erro, imov, errpsi)

      use sp_parameters, only: nrp, ntp, twopi
      use keys, only: kpr
      use compol, only: ngav, rm, zm, r, z, teta, q, errm, 
     &  iplas, iplas1, nt, nt1, nr, ro, ronor,
     &  psi, psia, psin, psip, psim
      use status_inc, only: error_catch

      implicit none

      integer, intent(in) :: imov
      real*8, intent(out) :: erro, errpsi

      integer :: i, j, ierm, jerm, jj
      real*8 :: rmold, zmold, alfa, rma, zma, psima, psnn, delpsn,
     &   drx, dzx, tetp, tetv, ttj, tt0, ttpl, zps, ps0, psmn, pspl, 
     &   grad, gradpl, gradmn, zro, rnj, znj,
     &   robj, ro2j, ro0, romn, ropl, roerr
      real*8, dimension(5) :: dp
      real*8, dimension(ntp) :: ron, rop, tetn, rob, teti, roi
      real*8, dimension(nrp, ntp) :: roplt

      if(kpr.eq.1) then
         write(*,*) 'regrid:::'
      endif
      if(ngav.eq.0) then
         alfa=1.0d0
      else
         alfa=1.0d0
      endif

      rmold=rm
      zmold=zm

      if(kpr.eq.1) then
         write(*,*) 'rm,zm,psim',rm,zm,psim
      endif

      call axdef(rma,zma,psima,dp)

      if(kpr.eq.1) then
         write(*,*) 'rma,zma,psima',rma,zma,psima
      endif
      if (isnan(rma)) then
         write(*,*) 'mag axis major radius is NaN'
         call error_catch
      endif
      rm=rma
      zm=zma
      psim=psima
      erro=SQRT((rmold-rm)**2+(zmold-zm)**2)

      if(kpr.eq.1) then
         write(6,*) 'erro mag.axis',erro
      endif

      errpsi=0.d0

      do i=1,iplas
         do j=1,nt
            psnn=psin(i,j)
            psin(i,j)=(psi(i,j)-psip)/(psim-psip)
            delpsn=ABS(psin(i,j)-psnn)
            if(delpsn.gt.errpsi) then
               ierm=i
               jerm=j
               errpsi=delpsn
            endif
         enddo
      enddo

      if(kpr.eq.1) then
         write(*,*) 'i,j errpsi',ierm,jerm,errpsi,(psim-psip)
      endif
      if(imov.eq.0) return

! definition of new angle grid teta(j)

      do j=1,nt
         tetn(j)=teta(j)
      enddo

      do j=1,nt
         drx=r(iplas,j)-rm
         dzx=z(iplas,j)-zm
         robj=SQRT(drx**2+dzx**2)
         rob(j)=robj
         tetp=ACOS(drx/robj)
         if(dzx.lt.0.d0) then
            teta(j)=-tetp
         else
            teta(j)=tetp
         endif
      enddo

      do j=2,nt
         if(teta(j).lt.teta(j-1)) then
            teta(j)=teta(j)+twopi
         endif
         if(teta(j).lt.teta(j-1)) then
            teta(j)=teta(j)+twopi
         endif
      enddo

      teta(1)=teta(nt1)-twopi
      teta(nt)=teta(2)+twopi

      i=2

      zps =psia(i)
      ps0=psip+zps*(psim-psip)

      do j=1,nt
         ttj=teta(j)
         ro2j=(ps0-psim)/( 0.5d0*dp(3)*COS(ttj)**2 +
     &                           dp(4)*COS(ttj)*SIN(ttj) +
     &                     0.5d0*dp(5)*SIN(ttj)**2 )
         rop(j)=SQRT(ro2j)
         r(i,j)=rm+rop(j)*COS(teta(j))
         z(i,j)=zm+rop(j)*SIN(teta(j))
      enddo

! grid moving along rais

      do i=2,iplas-1
         zps =psia(i)
         do j=1,nt
            ps0 =psin(i,j)
            psmn=psin(i-1,j)
            pspl=psin(i+1,j)
            ro0 =ro(i,j)
            romn=ro(i-1,j)
            ropl=ro(i+1,j)
            gradpl=-(pspl-ps0)/(ropl-ro0)
            gradmn=-(ps0-psmn)/(ro0-romn)
            grad=MAX(gradpl,gradmn)
            if (ngav.eq.0) then
               alfa=1.5d0
            else
               alfa=1.0d0+0.45d0*ABS(q(i)-q(1))/q(1)
            endif
            if(alfa.lt.1.5d0) alfa=1.5d0 
            if(ngav.gt.0) grad=grad*alfa
            zro=ro0-(psia(i)-ps0)/grad
            ron(j)=zro
         enddo

         do j=1,nt
            rnj=rmold+ron(j)*COS(tetn(j))
            znj=zmold+ron(j)*SIN(tetn(j))
            drx=rnj-rm
            dzx=znj-zm
            roi(j)=SQRT(drx**2+dzx**2)
            tetp=ACOS(drx/roi(j))
            if(dzx.lt.0.d0) then
               teti(j)=-tetp
            else
               teti(j)=tetp
            endif
         enddo

         do j=2,nt
            if(teti(j).lt.teti(j-1)) then
               teti(j)=teti(j)+twopi
            endif
            if(teti(j).lt.teti(j-1)) then
               teti(j)=teti(j)+twopi
            endif
	 enddo

         teti(1)=teti(nt1)-twopi
         teti(nt)=teti(2)+twopi
         roi(1)=roi(nt1)
         roi(nt)=roi(2)

         do j=2,nt1
            tetv=teta(j)
            if(tetv.lt.teti(1)) tetv=tetv+twopi
            if(tetv.lt.teti(1)) tetv=tetv+twopi
            if(tetv.gt.teti(nt)) tetv=tetv-twopi
            if(tetv.gt.teti(nt)) tetv=tetv-twopi
            do jj=1,nt1
               tt0 =teti(jj)
               ttpl=teti(jj+1)
               ro0 =roi(jj)
               ropl=roi(jj+1)
	       if(tetv.le.ttpl .and. tetv.ge.tt0) then
                  zro=((ttpl-tetv)*ro0+(tetv-tt0)*ropl)/(ttpl-tt0)
                  EXIT
	       endif
            enddo

            rop(j)=zro
            roerr=ABS(ron(j)-ro(i,j))

            if(roerr.gt.erro) then
               erro=roerr
               ierm=i
               jerm=j
            endif

            roplt(i,j)=ron(j)-ro(i,j)
            r(i,j)=rop(j)*COS(teta(j))+rm
            z(i,j)=rop(j)*SIN(teta(j))+zm
         enddo
      enddo

      if(kpr.eq.1) then
         write(6,*) 'erro max.',erro
         write(6,*) 'i,j erro',ierm,jerm,erro
      endif
      do i=2,iplas1
         do j=2,nt1
            ro(i,j)=SQRT((r(i,j)-rm)**2+(z(i,j)-zm)**2)
            ronor(i,j)=ro(i,j)/rob(j)
         enddo
      enddo

      do j=2,nt1
         ro(iplas,j)=rob(j)
         r(iplas,j)=rob(j)*COS(teta(j))+rm
         z(iplas,j)=rob(j)*SIN(teta(j))+zm
         ro(1,j)=0.d0
         r(1,j)=rm
         z(1,j)=zm
      enddo

      do i=1,iplas
         ro(i,1)=ro(i,nt1)
         ro(i,nt)=ro(i,2)
         ronor(i,1)=ronor(i,nt1)
         ronor(i,nt)=ronor(i,2)
         r(i,1)=r(i,nt1)
         r(i,nt)=r(i,2)
         z(i,1)=z(i,nt1)
         z(i,nt)=z(i,2)
         roplt(i,1)=roplt(i,nt1)
         roplt(i,nt)=roplt(i,2)
      enddo

      errm=erro

      do i=1,nr
         do j=1,nt
            psin(i,j)=psia(i)
            psi(i,j)=psip+psin(i,j)*(psim-psip)
         enddo
      enddo

      return
      end subroutine regrid

!----------------------------------------------------------------
      subroutine regrid0(erro, imov, errpsi)

      use sp_parameters, only: nrp, ntp, twopi
      use keys, only: kpr
      use compol, only: ngav, rm, zm, r, z, teta, errm, 
     & iplas, iplas1, nt, nt1, nr, ro, ronor,
     & psi, psia, psin, psip, psim
      use status_inc, only: error_catch

      implicit none

      integer, intent(in) :: imov
      real*8, intent(out) :: erro, errpsi

      integer :: i, j, is, ierm, jerm, jj
      real*8 :: rmold, zmold, alfa, rma, zma, psima, psnn, delpsn,
     &   drx, dzx, tetp, tetv, ttj, tt0, ttpl, zps, ps0, pspl, 
     &   zro, rnj, znj,
     &   robj, ro2j, ro0, ropl, roerr
      real*8, dimension(5) :: dp
      real*8, dimension(ntp) :: ron, rop, tetn, rob, teti, roi
      real*8, dimension(nrp, ntp) :: roplt

      if(kpr.eq.1) then
         write(*,*) 'regrid0:::'
      endif
      if(ngav.eq.0) then
         alfa=1.0d0
      else
         alfa=1.0d0
      endif

      rmold=rm
      zmold=zm

      if(kpr.eq.1) then
         write(*,*) 'rm,zm,psim',rm,zm,psim
      endif
      call axdef(rma,zma,psima,dp)

      if(kpr.eq.1) then
         write(*,*) 'rma,zma,psima',rma,zma,psima
      endif
      if (isnan(rma)) then
         write(*,*) 'mag axis major radius is NaN'
         call error_catch
      endif
      rm=rma
      zm=zma
      psim=psima
      erro=SQRT((rmold-rm)**2+(zmold-zm)**2)

      if(kpr.eq.1) then
         write(*,*) 'erro mag.axis',erro
      endif

      errpsi=0.d0

      do i=1,iplas
         do j=1,nt
	    psnn=psin(i,j)
            psin(i,j)=(psi(i,j)-psip)/(psim-psip)
            delpsn=ABS(psin(i,j)-psnn)
            if(delpsn.gt.errpsi) then
               ierm=i
               jerm=j
               errpsi=delpsn
            endif
         enddo
      enddo

      if(kpr.eq.1) then
         write(*,*) 'i,j errpsi',ierm,jerm,errpsi,(psim-psip)
      endif
      if(imov.eq.0) return

! definition of new angle grid teta(j)
      do j=1,nt
         tetn(j)=teta(j)
      enddo
      do j=1,nt
         drx=r(iplas,j)-rm
         dzx=z(iplas,j)-zm
         robj=SQRT(drx**2+dzx**2)
         rob(j)=robj
         tetp=ACOS(drx/robj)
         if(dzx.lt.0.d0) then
            teta(j)=-tetp
         else
            teta(j)=tetp
         endif
      enddo

      do j=2,nt
         if(teta(j).lt.teta(j-1)) then
            teta(j)=teta(j)+twopi
         endif
         if(teta(j).lt.teta(j-1)) then
            teta(j)=teta(j)+twopi
         endif
      enddo

      teta(1)=teta(nt1)-twopi
      teta(nt)=teta(2)+twopi

      i=2

      zps =psia(i)
      ps0=psip+zps*(psim-psip)

      do j=1,nt
         ttj=teta(j)
         ro2j=(ps0-psim)/( 0.5d0*dp(3)*COS(ttj)**2  +
     &                           dp(4)*COS(ttj)*SIN(ttj) +
     &                     0.5d0*dp(5)*SIN(ttj)**2 )
         rop(j)=SQRT(ro2j)
         r(i,j)=rm+rop(j)*COS(teta(j))
         z(i,j)=zm+rop(j)*SIN(teta(j))
      enddo

! grid moving along rais

      do i=3,iplas-1
         zps =psia(i)
         do j=1,nt
            do is=1,iplas1
               ps0 =psin(is,j)
               pspl=psin(is+1,j)
               ro0 =ro(is,j)
               ropl=ro(is+1,j)
	       if(zps.le.ps0 .and. zps.ge.pspl) then
                  zro=((pspl-zps)*ro0-(ps0-zps)*ropl)/(pspl-ps0)
                  EXIT
	       endif
            enddo
            ron(j)=zro
         enddo
         do j=1,nt
            rnj=rmold+ron(j)*COS(tetn(j))
            znj=zmold+ron(j)*SIN(tetn(j))
            drx=rnj-rm
            dzx=znj-zm
            roi(j)=SQRT(drx**2+dzx**2)
            tetp=ACOS(drx/roi(j))
            if(dzx.lt.0.d0) then
               teti(j)=-tetp
            else
               teti(j)=tetp
            endif
         enddo
         do j=2,nt
            if(teti(j).lt.teti(j-1)) then
               teti(j)=teti(j)+twopi
            endif
            if(teti(j).lt.teti(j-1)) then
               teti(j)=teti(j)+twopi
            endif
	 enddo

         teti(1)=teti(nt1)-twopi
         teti(nt)=teti(2)+twopi
         roi(1)=roi(nt1)
         roi(nt)=roi(2)

         do j=2,nt1
            tetv=teta(j)
            if(tetv.lt.teti(1)) tetv=tetv+twopi
            if(tetv.lt.teti(1)) tetv=tetv+twopi
            if(tetv.gt.teti(nt)) tetv=tetv-twopi
            if(tetv.gt.teti(nt)) tetv=tetv-twopi
            do jj=1,nt1
               tt0 =teti(jj)
               ttpl=teti(jj+1)
               ro0 =roi(jj)
               ropl=roi(jj+1)
	       if(tetv.le.ttpl .and. tetv.ge.tt0) then
                  zro=((ttpl-tetv)*ro0+(tetv-tt0)*ropl)/(ttpl-tt0)
                  EXIT
	       endif
            enddo
            rop(j)=zro
            roerr=ABS(ron(j)-ro(i,j))
            if(roerr.gt.erro) then
               erro=roerr
               ierm=i
               jerm=j
            endif
            roplt(i,j)=ron(j)-ro(i,j)
            r(i,j)=rop(j)*COS(teta(j))+rm
            z(i,j)=rop(j)*SIN(teta(j))+zm
         enddo
      enddo

      if(kpr.eq.1) then
         write(*,*) 'erro max.',erro
         write(*,*) 'i,j erro',ierm,jerm,erro
      endif
      do i=2,iplas1
         do j=2,nt1
            ro(i,j)=SQRT((r(i,j)-rm)**2+(z(i,j)-zm)**2)
            ronor(i,j)=ro(i,j)/rob(j)
         enddo
      enddo
      do j=2,nt1
         ro(iplas,j)=rob(j)
         r(iplas,j)=rob(j)*COS(teta(j))+rm
         z(iplas,j)=rob(j)*SIN(teta(j))+zm
         ro(1,j)=0.d0
         r(1,j)=rm
         z(1,j)=zm
      enddo
      do i=1,iplas
         ro(i,1)=ro(i,nt1)
         ro(i,nt)=ro(i,2)
         ronor(i,1)=ronor(i,nt1)
         ronor(i,nt)=ronor(i,2)
         r(i,1)=r(i,nt1)
         r(i,nt)=r(i,2)
         z(i,1)=z(i,nt1)
         z(i,nt)=z(i,2)
         roplt(i,1)=roplt(i,nt1)
         roplt(i,nt)=roplt(i,2)
      enddo

      errm=erro

      do i=1,nr
         do j=1,nt
            psin(i,j)=psia(i)
            psi(i,j)=psip+psin(i,j)*(psim-psip)
         enddo
      enddo

      return
      end subroutine regrid0

!----------------------------------------------------------------
      subroutine grid_1

      use bnd_modul, only: nbtab, rbtab, zbtab
      use sp_parameters, only: twopi
      use compol, only: nr, nt, nt1, iplas, iplas1,
     & ro, ronor, rm, zm, r, z, 
     & teta, psia, psin , psi, psip, psim

      implicit none

      integer, parameter :: nbtabp=1000

      integer :: i, j, ib, nbsh, nbn
      real*8 :: rm0, zm0, rc0, zc0, asp0, el_up, el_lw, tr_up, tr_lw,
     &   rbomax, rbomin, zbomax, zbomin, rc0new, zc0new, drc0, dzc0,
     &   drx, dzx, tetp
      real*8, dimension(nbtabp) :: robn, tetbn

      common /combsh/ rm0, zm0, rc0, zc0, asp0, el_up, el_lw, 
     &   tr_up, tr_lw, nbsh

      if(nbsh.eq.1) then

         call arc_x_bnd(nt)

! new position for magn.axis

         rbomax=0.d0
         rbomin=rm
         zbomax=zm
         zbomin=zm

         do ib=1,nbtab
            if(rbtab(ib).gt.rbomax) rbomax=rbtab(ib)
            if(rbtab(ib).lt.rbomin) rbomin=rbtab(ib)
            if(zbtab(ib).gt.zbomax) zbomax=zbtab(ib)
            if(zbtab(ib).lt.zbomin) zbomin=zbtab(ib)
	 enddo

         rc0new=0.5d0*(rbomax+rbomin)
         zc0new=0.5d0*(zbomax+zbomin)

         drc0=rc0new-rc0
         dzc0=zc0new-zc0

         rm=rm+drc0
         zm=zm+dzc0

      endif

      rbomax=0.d0
      rbomin=rm
      zbomax=zm
      zbomin=zm

      do ib=1,nbtab
         drx=rbtab(ib)-rm
         dzx=zbtab(ib)-zm

         robn(ib)=SQRT(drx**2+dzx**2)
         tetp=ACOS(drx/robn(ib))
         if(dzx.lt.0.d0) then
            tetbn(ib)=-tetp
         else
            tetbn(ib)=tetp
         endif
         if(rbtab(ib).gt.rbomax) rbomax=rbtab(ib)
         if(rbtab(ib).lt.rbomin) rbomin=rbtab(ib)
         if(zbtab(ib).gt.zbomax) zbomax=zbtab(ib)
         if(zbtab(ib).lt.zbomin) zbomin=zbtab(ib)
      enddo

      rc0=0.5d0*(rbomax+rbomin)
      zc0=0.5d0*(zbomax+zbomin)

      do ib=2,nbtab
         if(tetbn(ib).lt.tetbn(ib-1)) tetbn(ib)=tetbn(ib)+twopi
      enddo
      do ib=2,nbtab
          if(tetbn(ib).lt.tetbn(ib-1)) tetbn(ib)=tetbn(ib)+twopi
      enddo

      nbn=nbtab
      if(tetbn(1)+twopi-tetbn(nbn) .GT. 1.d-4) then
         nbn=nbn+1
         tetbn(nbn)=tetbn(1)+twopi
         robn(nbn)=robn(1)
      endif

      do j=1,nbn
         teta(j+1)=tetbn(j)
         ro(iplas,j+1)=robn(j)
         ro(1,j)=0.d0
      enddo

      teta(1)=teta(nt1)-twopi
      teta(nt)=teta(2)+twopi
      ro(iplas,nt)=ro(iplas,2)
      ro(iplas,1)=ro(iplas,nt1)

      do j=1,Nt
         tetp=teta(j)
         ro(1,j)=0.d0
         r(iplas,j)=rm+ro(iplas,j)*COS(tetp)
         z(iplas,j)=zm+ro(iplas,j)*SIN(tetp)
      enddo

      do i=2,iplas1
         do j=2,nt
	    ro(i,j)=ro(iplas,j)*ronor(i,j)
         enddo
      enddo

      do i=1,iplas1
         do j=2,nt1
	    r(i,j)=rm+ro(i,j)*COS(teta(j))
	    z(i,j)=zm+ro(i,j)*SIN(teta(j))
            psin(i,j)=psia(i)
         enddo
      enddo

      do j=1,nt
         r(1,j)=rm
         z(1,j)=zm
         ro(1,j)=0.d0
         psin(1,j)=1.d0
      enddo

      do i=1,nr
         r(i,1)=r(i,nt1)
         z(i,1)=z(i,nt1)
         ro(i,1)=ro(i,nt1)
         psin(i,1)=psin(i,nt1)

         ro(i,nt)=ro(i,2)
         r(i,nt)=r(i,2)
         z(i,nt)=z(i,2)
         psin(i,nt)=psin(i,2)
      enddo

      do i=1,nr
         do j=1,nt
            psi(i,j)=psip+psin(i,j)*(psim-psip)
         enddo
      enddo

      return
      end

!----------------------------------------------------------------
      subroutine grid_p0

      use compol, only: erru, nt, psi, iplas

      implicit none

      integer :: j

      erru=1.d0
      do j=1,nt
         psi(iplas, j)=0.d0 
      enddo

      return
      end subroutine grid_p0

!----------------------------------------------------------------
      subroutine grid_p1

      use compol, only: nt, iplas, psi

      implicit none

      integer :: j

      do j=1,nt
         psi(iplas, j)=0.d0
      enddo

      return
      end subroutine grid_p1

!----------------------------------------------------------------
      real*8 function frbon(r0, aa, tr, tet)

      implicit none
      real*8, intent(in) :: r0, aa, tr, tet

      frbon = r0 + (r0/aa)*COS(tet + tr*SIN(tet))

      return
      end function frbon

!----------------------------------------------------------------
      real*8 function fzbon(r0, z0, aa, el, tet)

      implicit none
      real*8, intent(in) :: r0, z0, aa, el, tet
      
      fzbon = z0 + (r0/aa)*el*SIN(tet)

      return
      end function fzbon

!----------------------------------------------------------------
      subroutine grid_0(igdf, nstep)

      use bnd_modul, only: nbtab, rbtab, zbtab
      use sp_parameters, only: twopi
      use compol, only: iplas, iplas1, rm, zm, r, z,
     & nr, nt, nt1, ro, teta,
     & psia, psim, psip, psin, psi, dpsda

      implicit none

      integer, parameter :: nbtabp=1000, nb4=nbtabp+4, nb6=nb4*6

      integer, intent(in) :: igdf, nstep

      integer :: i, j, nbsh, ib, nbn
      real*8 :: rm0, zm0, rc0, zc0, asp0, el_up, el_lw, tr_up, tr_lw,
     &   dtet, tet, tetp, tri, ell, rrr, zzz, roxx, drx, dzx,
     &   rbomin, rbomax, zbomin, zbomax
      real*8, dimension(nbtabp) :: robn, tetbn
      real*8, external :: frbon, fzbon

      common /combsh/ rm0, zm0, rc0, zc0, asp0, el_up, el_lw, 
     &   tr_up, tr_lw, nbsh

      if(nstep.eq.0) then
         if(igdf.eq.0) then
	    do i=1,iplas
               psia(i)=(iplas-i)/(iplas-1.d0)
            enddo
	    do i=1,iplas1
               dpsda(i)=-1.d0
            enddo
         elseif(igdf.eq.1 .OR. igdf.eq.2 .OR. igdf.eq.3) then
	    do i=1,iplas
               psia(i)=1.d0-((i-1)/(iplas-1.d0))**2
            enddo
	    do i=1,iplas1
               dpsda(i)=(1-2*i)/(iplas-1.d0)
            enddo
         endif
      endif

      rm=rm0
      zm=zm0
      psim=1.d0
      psip=0.d0

      if(nbsh.eq.0) then
         dtet=twopi/(nt-2)
         teta(1)=-dtet
         do j=2,nt
            teta(j)=teta(j-1)+dtet
         enddo
         teta(1)=teta(nt1)-twopi
         teta(nt)=teta(2)+twopi
         do j=1,nt
            tet=teta(j)
            tri=0.5d0*(tr_up+tr_lw+(tr_up-tr_lw)*SIN(tet))
            ell=0.5d0*(el_up+el_lw+(el_up-el_lw)*SIN(tet))
            rrr=frbon(rc0,asp0,tri,tet)
            r(iplas,j)=rrr
            zzz=fzbon(rc0,zc0,asp0,ell,tet)
            z(iplas,j)=zzz
            roxx=SQRT((r(iplas,j)-rm)**2+(z(iplas,j)-zm)**2)
            ro(iplas,j)=roxx
            ro(1,j)=0.d0
         enddo
 
	 nbtab=nt1
	 do ib=1,nbtab
            rbtab(ib)= r(iplas,ib)
            zbtab(ib)= z(iplas,ib)
         enddo
      endif

      call arc_x_bnd(nt)

      rbomax=0.d0
      rbomin=rm
      zbomax=zm
      zbomin=zm

      do ib=1,nbtab
         drx=rbtab(ib)-rm
         dzx=zbtab(ib)-zm
         robn(ib)=SQRT(drx**2+dzx**2)
         tetp=ACOS(drx/robn(ib))
         if(dzx.lt.0.d0) then
            tetbn(ib)=-tetp
         else
            tetbn(ib)=tetp
         endif
         if(rbtab(ib).gt.rbomax) rbomax=rbtab(ib)
         if(rbtab(ib).lt.rbomin) rbomin=rbtab(ib)
         if(zbtab(ib).gt.zbomax) zbomax=zbtab(ib)
         if(zbtab(ib).lt.zbomin) zbomin=zbtab(ib)
      enddo

      rc0=0.5d0*(rbomax+rbomin)
      zc0=0.5d0*(zbomax+zbomin)

      do ib=2,nbtab
         if(tetbn(ib).lt.tetbn(ib-1)) tetbn(ib)=tetbn(ib)+twopi
      enddo
      do ib=2,nbtab
         if(tetbn(ib).lt.tetbn(ib-1)) tetbn(ib)=tetbn(ib)+twopi
      enddo

      nbn=nbtab

      if(tetbn(1)+twopi-tetbn(nbn) .GT. 1.d-4) then
         nbn=nbn+1
         tetbn(nbn)=tetbn(1)+twopi
         robn(nbn)=robn(1)
      endif

      do j=1,nbn
         teta(j+1)=tetbn(j)
         ro(iplas,j+1)=robn(j)
         ro(1,j)=0.d0
      enddo

      teta(1)=teta(nt1)-twopi
      teta(nt)=teta(2)+twopi

      do j=2,Nt1
         tetp=teta(j)
         r(iplas,j)=rm+ro(iplas,j)*COS(tetp)
         z(iplas,j)=zm+ro(iplas,j)*SIN(tetp)
      enddo
      do i=2,iplas1
         do j=2,nt
            ro(i,j)=ro(iplas,j)*SQRT(1.d0-psia(i))
         enddo
      enddo

      do i=1,iplas1
         do j=2,nt1
            r(i,j)=rm+ro(i,j)*COS(teta(j))
            z(i,j)=zm+ro(i,j)*SIN(teta(j))
            psin(i,j)=psia(i)
         enddo
      enddo

      do j=1,nt
         r(1,j)=rm
         z(1,j)=zm
         ro(1,j)=0.d0
         psin(1,j)=1.d0
      enddo

      do i=1,nr
         r(i,1)=r(i,nt1)
         z(i,1)=z(i,nt1)
         ro(i,1)=ro(i,nt1)
         psin(i,1)=psin(i,nt1)
         r(i,nt)=r(i,2)
         z(i,nt)=z(i,2)
         psin(i,nt)=psin(i,2)
      enddo

      do i=1,nr
         do j=1,nt
            psi(i,j)=psip+psin(i,j)*(psim-psip)
         enddo
      enddo

      call wrb

      return
      end subroutine grid_0

!----------------------------------------------------------------
      subroutine loopL_b(tetpol, ntet, ro0, u0)

      use sp_parameters, only: nip, njp, nbndp2, pi, twopi
      use keys, only: kpr

      implicit none

      integer, intent(in) :: ntet
      real*8, intent(in) :: u0
      real*8, intent(in) , dimension(*) :: tetpol
      real*8, intent(out), dimension(*) :: ro0

      integer :: ni, nj, ni1, nj1, ni2, nj2, nbnd, nkin, nkout,
     &   i, j, ii, jj, ix1, ix2, jx1, jx2, ic, ic1, jc, jc1, lin, 
     &   nxb, ig, icell, jcell, imax, imax1, jmax, jmax1, jw
      real*8 :: ux0, ux1, ux2, up, um, xm, ym, xx0, yx0, xx1, yx1, 
     &   xx2, yx2, xx10, yx10, xx20, yx20, xm0, ym0, psi_bon,
     &   ropl, romn, tetpl, tetmn, drx, dzx, tetp
      real*8, dimension(nip) :: x, dx, dxi, x12
      real*8, dimension(njp) :: y, dy, dyj
      real*8, dimension(nip, njp) :: u, ue, un, ui, g, ut
      real*8, dimension(nbndp2) :: rxb, zxb, roxb, tetxb

      real*8, external :: xzer

      common /comind/ ni, nj, ni1, nj1, ni2, nj2, nbnd, nkin, nkout
      common /comrz/ x, y, dx, dy, dxi, dyj, x12
      common /compot/ u, ue, un, ui, g, 
     &                ux0, ux1, ux2, up, um, xm, ym, 
     &                xx0, yx0, xx1, yx1, xx2, yx2, imax, jmax, 
     &                ix1, jx1, ix2, jx2, 
     &                xx10, yx10, xx20, yx20, xm0, ym0,
     &                psi_bon

      imax1=imax-1
      jmax1=jmax-1

      do i=imax1,imax
         if(xm.le.x(i+1) .AND. xm.gt.x(i)) icell=i
      enddo

      do j=jmax1,jmax
         if(ym.le.y(j+1) .AND. ym.gt.y(j)) jcell=j
      enddo

      i=icell
      j=jcell+1

      do ii=1,ni
         do jj=1,nj
            ut(ii,jj)=un(ii,jj)-u0
         enddo
      enddo

      ig=1

      do i=imax,ni1
         if(ut(i,j)*ut(i+1,j).le.0.d0) then
            ic=i
            jc=j
            rxb(ig)=xzer(x(i+1),x(i),ut(i+1,j),ut(i,j))
            zxb(ig)=y(j)
            goto 886
         endif
      enddo

      if(kpr.eq.1) then
         write(*,*) 'loopL: first point was not find'
      endif
      stop

 886  continue

      ic1=ic
      jc1=jc
      lin=1

      do

         ig=ig+1
         i=ic
         j=jc

         if(ut(i+1,j)*ut(i,j).le.0.d0 .AND. lin.ne.1) then
            ic=i
            jc=j-1
            lin=3
            rxb(ig)=xzer(x(i+1),x(i),ut(i+1,j),ut(i,j))
            zxb(ig)=y(j)
         elseif(ut(i+1,j+1)*ut(i+1,j).le.0.d0.AND. lin.ne.2) then
            ic=i+1
            jc=j
            lin=4
            rxb(ig)=x(i+1)
            zxb(ig)=xzer(y(j+1),y(j),ut(i+1,j+1),ut(i+1,j))
         elseif(ut(i+1,j+1)*ut(i,j+1).le.0.d0 .AND. lin.ne.3) then
            ic=i
            jc=j+1
            lin=1
            rxb(ig)=xzer(x(i+1),x(i),ut(i+1,j+1),ut(i,j+1))
            zxb(ig)=y(j+1)
         elseif(ut(i,j+1)*ut(i,j).le.0.d0.AND. lin.ne.4) then
            ic=i-1
            jc=j
            lin=2
            rxb(ig)=x(i)
            zxb(ig)=xzer(y(j+1),y(j),ut(i,j+1),ut(i,j))
         endif

         if(jc.eq.jc1 .AND. ic.eq.ic1) then
            ig=ig+1
            rxb(ig)=rxb(2)
            zxb(ig)=zxb(2)
            EXIT
         endif

      enddo

      nxb=ig

      do ig=1,nxb
         drx=rxb(ig)-xm
         dzx=zxb(ig)-ym
         tetxb(ig)=ATAN2(dzx, drx)
         roxb(ig)=SQRT(drx**2+dzx**2)
      enddo

      do j=1,Ntet
         tetp=tetpol(j)
         if(tetp.lt.tetxb(1)) then
            tetp=tetp+twopi
         elseif(tetp.gt.tetxb(nxb)) then
            tetp=tetp-twopi
         endif
         do jw=1,nxb-1
            tetmn=tetxb(jw)
            tetpl=tetxb(jw+1)
            romn=roxb(jw)
            ropl=roxb(jw+1)
            if(tetp.le.tetpl .AnD. tetp.ge.tetmn) EXIT
         enddo
         ro0(j)=(ropl*(tetp-tetmn)+romn*(tetpl-tetp))/(tetpl-tetmn)
      enddo

      return
      end subroutine loopL_b

!----------------------------------------------------------------
      subroutine grid_b(igdf, nstep)

      use sp_parameters, only: nip, njp, ntp, twopi
      use compol, only: iplas, iplas1, rm, zm, r, z,
     & nr, nt, nt1, ro, teta,
     & psia, psim, psip, psin, psi, dpsda

      implicit none

      integer, parameter :: nshp=10, nbtabp=1000, nb4=nbtabp+4, 
     &   nb6=nb4*6

      integer, intent(in) :: igdf, nstep

      integer :: i, j, k, l, ii, jj, ibeg, imax, jmax, nsh,
     &   ix1, ix2, jx1, jx2
      real*8 :: ux0, ux1, ux2, up, um, xm, ym, xx0, yx0, xx1, yx1, 
     &   xx2, yx2, xx10, yx10, xx20, yx20, xm0, ym0, psi_bon
      real*8 :: u0, ur0, stpx, stpy, dtet, ttj, ro2j
      real*8, dimension(5) :: dp
      real*8, dimension(nip) :: x, dx, dxi, x12
      real*8, dimension(njp) :: y, dy, dyj
      real*8, dimension(nip, njp) :: u, ue, un, ui, g
      real*8, dimension(ntp) :: ron, ronm
      real*8, dimension(nshp) :: xs, ys, fun

      common /comrz/ x, y, dx, dy, dxi, dyj, x12
      common /compot/ u, ue, un, ui, g,
     &                ux0, ux1, ux2, up, um, xm, ym, 
     &                xx0, yx0, xx1, yx1, xx2, yx2, imax, jmax, 
     &                ix1, jx1, ix2, jx2, 
     &                xx10, yx10, xx20, yx20, xm0, ym0, 
     &                psi_bon

      if(nstep.eq.0) then
         if(igdf.eq.0)then
	    do i=1,iplas
               psia(i)=(iplas-i)/(iplas-1.d0)
            enddo
	    do i=1,iplas1
	       dpsda(i)=-1.d0
            enddo
         elseif(igdf.eq.1 .OR. igdf.eq.2 .OR. igdf.eq.3) then
	    do i=1,iplas
               psia(i)=1.d0-((i-1)/(iplas-1.d0))**2
            enddo
	    do i=1,iplas1
               dpsda(i)=(1-2*i)/(iplas-1.d0)
            enddo
         endif
      endif

      rm=xm
      zm=ym
      psim=um-up
      psip=0.d0

      dtet=twopi/(nt-2)
      teta(1)=-dtet

      do j=2,nt
         teta(j)=teta(j-1)+dtet
      enddo

      teta(1)=teta(nt1)-twopi
      teta(nt)=teta(2)+twopi
      psia(1)=1.d0

      nsh=1
      xs(nsh)=x(imax)
      ys(nsh)=y(jmax)
      fun(nsh)=u(imax,jmax)

      do k=-1,1
         ii= imax+k
         do l=-1,1
            jj= jmax+l
            if(ii.ne.imax .OR. jj.ne.jmax) then
               nsh=nsh+1
               xs(nsh)=x(ii)
               ys(nsh)=y(jj)
               fun(nsh)=u(ii,jj)
            endif
         enddo
      enddo

      call deriv5(xs,ys,fun,nsh,5,dp)

      stpx=(x(imax+1)-x(imax))
      stpy=(y(jmax+1)-y(jmax))*0.75d0

      do i=2,iplas
         u0=psia(i)
         ur0=up+u0*(um-up)
         do j=1,nt
            ttj=teta(j)
            ro2j=(ur0-um)/( 0.5d0*dp(3)*COS(ttj)**2 +
     &                            dp(4)*COS(ttj)*SIN(ttj) +
     &                      0.5d0*dp(5)*SIN(ttj)**2 )
            ron(j)=SQRT(ro2j)
            ronm(j)=ron(j)    
            r(i,j)=rm+ron(j)*COS(teta(j))
            z(i,j)=zm+ron(j)*SIN(teta(j))
            ro(i,j)=ron(j)
            psin(i,j)=u0
         enddo
	 ibeg=i+1
	 if( ron(2) .GT. MAX(stpx,stpy) ) EXIT
      enddo

      do i=iplas,ibeg,-1
         u0=psia(i)
         call loopL_b(teta,nt,ron,u0)
         do j=1,nt
            if(ron(j) .LT. ronm(j)) then
               ron(j)=(ro(i+1,j)+ronm(j)*(i-ibeg+1.d0))/(i-ibeg+2.d0)
            endif
         enddo
         do j=1,nt
            r(i,j)=rm+ron(j)*COS(teta(j))
            z(i,j)=zm+ron(j)*SIN(teta(j))
            ro(i,j)=ron(j)
            psin(i,j)=u0
         enddo
      enddo

      do j=1,nt
         r(1,j)=rm
         z(1,j)=zm
         ro(1,j)=0.d0
         psin(1,j)=1.d0
      enddo

      do i=1,nr
         r(i,1)=r(i,nt1)
         z(i,1)=z(i,nt1)
         ro(i,1)=ro(i,nt1)
         r(i,nt)=r(i,2)
         z(i,nt)=z(i,2)
         ro(i,nt)=ro(i,2)
      enddo

      teta(1)=teta(nt1)-twopi
      teta(nt)=teta(2)+twopi

      do i=1,nr
         do j=1,nt
            psi(i,j)=psip+psin(i,j)*(psim-psip)
         enddo
      enddo

      return
      end subroutine grid_b

!----------------------------------------------------------------
      subroutine grid_b_ef

      use numerical_tools, only: reinterp_back_quad
      use sp_parameters, only: nip, njp, twopi
      use compol, only: iplas, iplas1, rm, zm, r, z,
     & nr, nt, nt1, ro, ronor, teta,
     & psia, psim, psip, psin, psi, dpsda

      implicit none

      integer, parameter :: nshp=10, nbtabp=1000, nb4=nbtabp+4, 
     &   nb6=nb4*6

      integer :: i, j, k, imax, jmax,
     &   ix1, ix2, jx1, jx2, j1, j2
      real*8 :: ux0, ux1, ux2, up, um, xm, ym, xx0, yx0, xx1, yx1, 
     &   xx2, yx2, xx10, yx10, xx20, yx20, xm0, ym0, psi_bon
      real*8 :: u0, ur0, dtet, xstepp, u002,
     &   um_t, up_t, xm_t, ym_t, dumx, dumy, ronzz, ronzz2
      real*8, dimension(9) :: uuuu
      real*8, dimension(9, 2) :: xxx
      real*8, dimension(500) :: ucaz
      real*8, dimension(nip) :: x, dx, dxi, x12
      real*8, dimension(njp) :: y, dy, dyj
      real*8, dimension(nip, njp) :: u, ue, un, ui, g, utemp
      real*8, dimension(iplas) :: rggr, psitemp

      common /comrz/ x, y, dx, dy, dxi, dyj, x12
      common /compot/ u, ue, un, ui, g,
     &                ux0, ux1, ux2, up, um, xm, ym, 
     &                xx0, yx0, xx1, yx1, xx2, yx2, imax, jmax, 
     &                ix1, jx1, ix2, jx2, 
     &                xx10, yx10, xx20, yx20, xm0, ym0, 
     &                psi_bon

      do i=1,iplas
         psia(i)=1.d0-((i-1)/(iplas-1.d0))**2
      enddo

      do i=1,iplas1
         dpsda(i)=(1-2*i)/(iplas-1.d0)
      enddo
      utemp=u
      up_t=up
      xm_t=xm
      ym_t=ym
      um_t=um
      rm=xm_t
      zm=ym_t
      psim=um_t-up_t
      psip=0.d0

      dtet=twopi/(nt-2)
      teta(1)=-dtet

      do j=2,nt
         teta(j)=teta(j-1)+dtet
      enddo

      teta(1)=teta(nt1)-twopi
      teta(nt)=teta(2)+twopi
      psia(1)=1.d0

! first find boundary
      i=iplas
      xstepp=0.1*SQRT((x(2)-x(1))**2.+(y(2)-y(1))**2.)

      do j=1,nt
         do k=1,nip*njp
            ronzz=xstepp*k
            dumx=rm+ronzz*COS(teta(j))
            dumy=zm+ronzz*SIN(teta(j))
            j1=((dumx-x(1))/(x(2)-x(1)))+1	
            j2=((dumy-y(1))/(y(2)-y(1)))+1	
            uuuu(1)=utemp(j1-1,j2-1)
            uuuu(2)=utemp(j1,j2-1)
            uuuu(3)=utemp(j1-1,j2)
            uuuu(4)=utemp(j1,j2)
            uuuu(5)=utemp(j1-1,j2+1)
            uuuu(6)=utemp(j1,j2+1)
            uuuu(7)=utemp(j1+1,j2+1)
            uuuu(8)=utemp(j1+1,j2)
            uuuu(9)=utemp(j1+1,j2-1)

	    call biquadratic_interp_spidef(x(j1-1), x(j1), x(j1+1),
     &             y(j2-1), y(j2), y(j2+1), uuuu, dumx, dumy, u002)
            ucaz(k)=u002
            if (ucaz(k).eq.up_t) then
               ro(i,j)=ronzz
               EXIT
            endif
            if (ucaz(k).lt.up_t) then
               ronzz2=xstepp*(k-1)
               ro(i,j)=ronzz2+(up_t-ucaz(k-1))/(ucaz(k-1)-ucaz(k))*
     &                (ronzz2-ronzz)
               EXIT
            endif
         enddo

         psitemp(1)=1.
         psitemp(iplas)=0.
         rggr(iplas)=ro(iplas,j)
         rggr(1)=0.

         do i=2,iplas-1
            rggr(i)=ro(iplas,j)*(i-1.)/(iplas-1.)
            dumx=rm+rggr(i)*COS(teta(j))
            dumy=zm+rggr(i)*SIN(teta(j))
            j1=((dumx-x(1))/(x(2)-x(1)))+1	
            j2=((dumy-y(1))/(y(2)-y(1)))+1	

            xxx(1,1)=x(j1-1)
            xxx(1,2)=y(j2-1)
            uuuu(1)=utemp(j1-1,j2-1)
            xxx(2,1)=x(j1)
            xxx(2,2)=y(j2-1)
            uuuu(2)=utemp(j1,j2-1)
            xxx(3,1)=x(j1-1)
            xxx(3,2)=y(j2)
            uuuu(3)=utemp(j1-1,j2)
            xxx(4,1)=x(j1)
            xxx(4,2)=y(j2)
            uuuu(4)=utemp(j1,j2)
            xxx(5,1)=x(j1-1)
            xxx(5,2)=y(j2+1)
            uuuu(5)=utemp(j1-1,j2+1)
            xxx(6,1)=x(j1)
            xxx(6,2)=y(j2+1)
            uuuu(6)=utemp(j1,j2+1)
            xxx(7,1)=x(j1+1)
            xxx(7,2)=y(j2+1)
            uuuu(7)=utemp(j1+1,j2+1)
            xxx(8,1)=x(j1+1)
            xxx(8,2)=y(j2)
            uuuu(8)=utemp(j1+1,j2)
            xxx(9,1)=x(j1+1)
            xxx(9,2)=y(j2-1)
            uuuu(9)=utemp(j1+1,j2-1)

            call biquadratic_interp_spidef(x(j1-1), x(j1), x(j1+1),
     &          y(j2-1), y(j2), y(j2+1), uuuu, dumx, dumy, u002)

            psitemp(i)=(u002-up_t)/(um_t-up_t)
	 enddo

         call reinterp_back_quad((1.-psitemp(1:iplas))**0.5,
     &         rggr(1:iplas), iplas, 
     &         (1.-psia(1:iplas))**0.5,ro(1:iplas,j),iplas)

      enddo

      do i=2,iplas
         u0=psia(i)
         ur0=up_t+u0*(um_t-up_t)
         do j=1,nt
            r(i,j)=rm+ro(i,j)*COS(teta(j))
            z(i,j)=zm+ro(i,j)*SIN(teta(j))
            psin(i,j)=u0
         enddo
      enddo
      do j=1,nt
         r(1,j)=rm
         z(1,j)=zm
         ro(1,j)=0.d0
         psin(1,j)=1.d0
      enddo
      do i=1,nr
         r(i,1)=r(i,nt1)
         z(i,1)=z(i,nt1)
         ro(i,1)=ro(i,nt1)
         r(i,nt)=r(i,2)
         z(i,nt)=z(i,2)
         ro(i,nt)=ro(i,2)
      enddo
      do j=1,nt
         ronor(iplas,j)=1.d0
      enddo
      do j=1,nt 
         do i=1,iplas
            ronor(i,j)=ro(i,j)/ro(iplas,j)
            psin(i,j)=psia(i)
         enddo				  
      enddo				  

      teta(1)=teta(nt1)-twopi
      teta(nt)=teta(2)+twopi

      do i=1,nr
         do j=1,nt
            psi(i,j)=psip+psin(i,j)*(psim-psip)
         enddo
      enddo

      return
      end subroutine grid_b_ef

!-----------------------------------------------------
Cquadratic inerpolation

      subroutine biquadratic_interp_spidef(x1, x2, x3,
     &   z1, z2, z3, uuuu, r2, thet2, y2)

      implicit none

      double precision, intent(in) :: x1, x2, x3, r2, thet2
      double precision, intent(in), dimension(9) :: uuuu
      double precision, intent(out) :: z1, z2, z3, y2

      integer :: j
      double precision :: A, B, C, t1, t2, t3
      double precision, dimension(3) :: yy, thet1, r0, p0
      double precision, dimension(3, 3) :: r1, y1

      r1(1,1)=x1
      r1(2,1)=x2
      r1(3,1)=x3
      r1(:,2)=r1(:,1)
      r1(:,3)=r1(:,1)	
      thet1(1)=z1
      thet1(2)=z2
      thet1(3)=z3
      y1(1,1)=uuuu(1)
      y1(2,1)=uuuu(2)
      y1(3,1)=uuuu(9)
      y1(1,2)=uuuu(3)
      y1(2,2)=uuuu(4)
      y1(3,2)=uuuu(8)
      y1(1,3)=uuuu(5)
      y1(2,3)=uuuu(6)
      y1(3,3)=uuuu(7)

      do j=1,3
         r0=r1(:,j)
         p0=y1(:,j)
         t1=p0(1)
         t2=p0(2)
         t3=p0(3)
         z1=r0(1)
         z2=r0(2)
         z3=r0(3)
         A = (t3-t2-(z3-z2)*(t1-t2)/(z1-z2)) / ((z3-z2)*(z3-z1))
         B = (t1-t2)/(z1-z2) - A*(z1+z2)
         C=t2-A*(z2**2.0)-B*z2
      	 yy(j)=A*(r2**2.0)+B*r2+C
      enddo

      t1=yy(1)
      t2=yy(2)
      t3=yy(3)
      z1=thet1(1)
      z2=thet1(2)
      z3=thet1(3)
      A = (t3-t2-(z3-z2)*(t1-t2)/(z1-z2)) / ((z3-z2)*(z3-z1))
      B = (t1-t2)/(z1-z2) - A*(z1+z2)
      C = t2-A*(z2**2.0)-B*z2
      y2 = A*(thet2**2.0) + B*thet2 + C

      return
      end subroutine biquadratic_interp_spidef

!----------------------------------------------------------------
      subroutine grid_b1

      use sp_parameters, only: nip, njp, ntp, twopi
      use compol, only: iplas, rm, zm, r, z,
     & nr, nt, nt1, ro, ronor, teta,
     & psia, psim, psip, psin, psi

      implicit none

      integer :: i, j, imax, jmax, ix1, jx1, ix2, jx2
      real*8 :: ux0, ux1, ux2, up, um, xm, ym, 
     &          xx0, yx0, xx1, yx1, xx2, yx2, 
     &          xx10, yx10, xx20, yx20, xm0, ym0, psi_bon,
     &   dtet, u0
      real*8, dimension(ntp) :: ron
      real*8, dimension(nip, njp) :: u, ue, un, ui, g

      common /compot/ u, ue, un, ui, g, 
     &                ux0, ux1, ux2, up, um, xm, ym, 
     &                xx0, yx0, xx1, yx1, xx2, yx2, imax, jmax, 
     &                ix1, jx1, ix2, jx2, 
     &                xx10, yx10, xx20, yx20, xm0, ym0,
     &                psi_bon

      rm=xm 
      zm=ym    !  position of magn. axis

      psip=0.d0

      dtet=twopi/(nt-2)
      teta(1)=-dtet

      do j=2,nt
         teta(j)=teta(j-1)+dtet
      enddo

      teta(1)=teta(nt1)-twopi
      teta(nt)=teta(2)+twopi

      psia(1)=1.d0
      u0=0.d0

      call loopL_b(teta,nt,ron,u0)

      do j=1,nt
         ronor(iplas,j)=1.d0
      enddo				  

      do j=1,nt 
         do i=1,iplas
            ro(i,j)=ronor(i,j)*ron(j)
            r(i,j)=rm+ro(i,j)*COS(teta(j))
            z(i,j)=zm+ro(i,j)*SIN(teta(j))
            psin(i,j)=psia(i)
         enddo				  
      enddo				  

      do j=1,nt
         r(1,j)=rm
         z(1,j)=zm
         ro(1,j)=0.d0
         ronor(1,j)=0.d0
         psin(1,j)=1.d0
      enddo

      do i=1,nr
         r(i,1)=r(i,nt1)
         z(i,1)=z(i,nt1)
         ro(i,1)=ro(i,nt1)
         ronor(i,1)=ronor(i,nt1)
         r(i,nt)=r(i,2)
         z(i,nt)=z(i,2)
         ro(i,nt)=ro(i,2)
         ronor(i,nt)=ronor(i,2)
      enddo

      do i=1,nr
         do j=1,nt
            psi(i,j)=psip+psin(i,j)*(psim-psip)
         enddo
      enddo

      return
      end subroutine grid_b1

!----------------------------------------------------------------
      subroutine arc_x_bnd(nteta)

      use sp_parameters, only: pi
      use bnd_modul, only: nbtab, rbtab, zbtab
      use keys, only: kpr

      implicit none

      integer, parameter :: ntz=1000, nsz=6*ntz+16

      integer, intent(in) :: nteta

      integer :: j, ib, m, m11, mx, m1x, m2x, jvmin, jvmin2, nbnd, 
     &   isb, isbo, isbw
      real*8 :: abo, ccur, clock, cmin, culen, epsbo, hbnd, tolen, 
     &   rint, zint, rmin, rmax, zmin, zmax, dr, dz, x2len
      real*8, dimension(ntz) :: uk1, vk1, ukw, vkw, t1, t1n
      real*8, dimension(nsz) :: w1, w2, w3

      do ib=2,nbtab
         dr=ABS(rbtab(ib)-rbtab(1))
         dz=ABS(zbtab(ib)-zbtab(1))
         if(dr.lt.1.d-6 .AND. dz.lt.1.d-6) then
            nbnd=ib-1
            goto 1947
         endif
      enddo
      nbnd=nbtab

 1947 continue

      m=nbnd+1
      do j=1,nbnd
         ukw(j)=rbtab(j)
         vkw(j)=zbtab(j)
      enddo
      ukw(m)=ukw(1)
      vkw(m)=vkw(1)

! x-point as starting
!  internal point as middle
      rmin=ukw(1)
      rmax=ukw(1)
      zmin=vkw(1)
      zmax=vkw(1)
      do j=2,m
         if(ukw(j).lt.rmin) then
            rmin=ukw(j)
         endif
         if(ukw(j).gt.rmax) then
            rmax=ukw(j)
         endif
         if(vkw(j).lt.zmin) then
            zmin=vkw(j)
         endif
         if(vkw(j).gt.zmax) then
            zmax=vkw(j)
         endif
      enddo
      rint=0.5d0*(rmin+rmax)
      zint=0.5d0*(zmin+zmax)
! orientation
      clock=(ukw(1)-rint)*(vkw(2)-zint)-(vkw(1)-zint)*(ukw(2)-rint)
      j=1
      cmin=((ukw(j)-ukw(m-1))*(ukw(j+1)-ukw(j))  +
     &      (vkw(j)-vkw(m-1))*(vkw(j+1)-vkw(j))) /
     &      SQRT(((ukw(j)-ukw(m-1))**2+(vkw(j)-vkw(m-1))**2) *
     &            ((ukw(j+1)-ukw(j))**2+(vkw(j+1)-vkw(j))**2))
      jvmin=1

      do j=2,m-1
         ccur=((ukw(j)-ukw(j-1))*(ukw(j+1)-ukw(j))  +
     &         (vkw(j)-vkw(j-1))*(vkw(j+1)-vkw(j))) /
     &         SQRT(((ukw(j)-ukw(j-1))**2+(vkw(j)-vkw(j-1))**2) *
     &               ((ukw(j+1)-ukw(j))**2+(vkw(j+1)-vkw(j))**2))
! if(clock.lt.0.d0) ccur=-ccur
         if(ccur.lt.cmin) then
            cmin=ccur
            jvmin=j
         endif
      enddo

      do j=jvmin,m-1
         uk1(j-jvmin+1)=ukw(j)
         vk1(j-jvmin+1)=vkw(j)
      enddo
      do j=1,jvmin-1
         uk1(m-jvmin+j)=ukw(j)
         vk1(m-jvmin+j)=vkw(j)
      enddo
      uk1(m)=uk1(1)
      vk1(m)=vk1(1)

      if(kpr.eq.1) then
         write(*,*) ' bound points ',m
         write(*,*) ' bound start point ',uk1(1),vk1(1)
         write(*,*) ' min angle ',
     &      (1.d0-ACOS(cmin)/pi)*180.d0,' degree'
      endif

      if((1.d0-ACOS(cmin)/pi)*180.d0.lt.100.d0) then
         uk1(2)=0.5d0*(uk1(1)+uk1(3))
         uk1(m-1)=0.5d0*(uk1(m)+uk1(m-2))
         vk1(2)=0.5d0*(vk1(1)+vk1(3))
         vk1(m-1)=0.5d0*(vk1(m)+vk1(m-2))
         if(kpr.eq.1) then
            write(*,*) ' linear interpolated near x-point '
	 endif
      endif

      do j=1,m
        ukw(j)=uk1(j)
        vkw(j)=vk1(j)
      enddo

! second x-point search
      j=2
      cmin=((ukw(j)-ukw(j-1))*(ukw(j+1)-ukw(j))  +
     &      (vkw(j)-vkw(j-1))*(vkw(j+1)-vkw(j))) /
     &       SQRT(((ukw(j)-ukw(j-1))**2+(vkw(j)-vkw(j-1))**2) *
     &             ((ukw(j+1)-ukw(j))**2+(vkw(j+1)-vkw(j))**2))
      jvmin2=2
      do j=2,m-1
         ccur=((ukw(j)-ukw(j-1))*(ukw(j+1)-ukw(j))  +
     &         (vkw(j)-vkw(j-1))*(vkw(j+1)-vkw(j))) /
     &        SQRT(((ukw(j)-ukw(j-1))**2+(vkw(j)-vkw(j-1))**2) *
     &              ((ukw(j+1)-ukw(j))**2+(vkw(j+1)-vkw(j))**2))
! if(clock.lt.0.d0) ccur=-ccur
         if(ccur.lt.cmin) then
            cmin=ccur
            jvmin2=j
         endif
      enddo
! check the x-point angle
      if((1.d0-ACOS(cmin)/pi)*180.d0.lt.100.d0) then
         uk1(jvmin2+1)=0.5d0*(uk1(jvmin2)+uk1(jvmin2+2))
         uk1(jvmin2-1)=0.5d0*(uk1(jvmin2)+uk1(jvmin2-2))
         vk1(jvmin2+1)=0.5d0*(vk1(jvmin2)+vk1(jvmin2+2))
         vk1(jvmin2-1)=0.5d0*(vk1(jvmin2)+vk1(jvmin2-2))
         if(kpr.eq.1) then
            write(*,*) ' min angle ',
     &         (1.d0-ACOS(cmin)/pi)*180.d0,' degree'
            write(*,*) ' linear interpolated near the second x-point '
         endif

! arclength from the first to the second x-point
         x2len=0.d0
         tolen=0.d0
         do j=2,m
            culen=SQRT((uk1(j)-uk1(j-1))**2+(vk1(j)-vk1(j-1))**2)
            tolen=tolen+culen
            if(j.le.jvmin2) then
               x2len=x2len+culen
            endif
         enddo

! distribute the point over the two branches
! first branch
        m11=nteta-1
        m1x=(m11-1)*(x2len/tolen)+1
        m2x=m11-m1x+1
        mx=jvmin2
        hbnd=1.d0/(m1x-1)
        t1(1)=0.d0
        do j=2,m1x-1
           t1(j)=t1(j-1)+hbnd
        enddo
        t1(m1x)=1.d0
        call ARCM(MX,UK1,VK1,T1N)
        call SPLNA1(NTZ,MX,T1N,UK1,M1X,T1,NSZ,W1,W2,W3)
        call SPLNA1(NTZ,MX,T1N,VK1,M1X,T1,NSZ,W1,W2,W3)

        ISBO=20
        EPSBO=1.d-10
        do ISB=1,ISBO
           ISBW=ISB
           call ARCM(M1X ,UK1,VK1,T1N)
           call SPLNA1(NTZ,M1X ,T1N,UK1,M1X,T1,NSZ,W1,W2,W3)
           call SPLNA1(NTZ,M1X ,T1N,VK1,M1X,T1,NSZ,W1,W2,W3)
           ABO=0.d0
           do J=1,M1X
              ABO=MAX(ABO,ABS(T1N(J)-T1(J)))
           enddo
           if (ABO.LE.EPSBO) EXIT
        enddo

! second branch
        mx=m-jvmin2+1
        do j=1,mx
           uk1(j+m1x-1)=ukw(j+jvmin2-1)
           vk1(j+m1x-1)=vkw(j+jvmin2-1)
        enddo
        hbnd=1.d0/(m2x-1)
        t1(1)=0.d0
        do j=2,m2x-1
           t1(j)=t1(j-1)+hbnd
        enddo
        t1(m2x)=1.d0
        call ARCM(MX,UK1(m1x),VK1(m1x),T1N)
        call SPLNA1(NTZ,MX,T1N,UK1(m1x),M2X,T1,NSZ,W1,W2,W3)
        call SPLNA1(NTZ,MX,T1N,VK1(m1x),M2X,T1,NSZ,W1,W2,W3)

        ISBO=20
        EPSBO=1.d-10
        do ISB=1,ISBO
           ISBW=ISB
           call ARCM(M2X ,UK1(m1x),VK1(m1x),T1N)
           call SPLNA1(NTZ,M2X ,T1N,UK1(m1x),M2X,T1,NSZ,W1,W2,W3)
           call SPLNA1(NTZ,M2X ,T1N,VK1(m1x),M2X,T1,NSZ,W1,W2,W3)
           ABO=0.d0
           do J=1,M2X
              ABO=MAX(ABO,ABS(T1N(J)-T1(J)))
           enddo
           if (ABO.LE.EPSBO) EXIT
        enddo

        write(*,*) ' '
        write(*,*) ' 2x    ITERATIONS TO MATCH ARCL. DISTR. ',
     &             'AT THE BOUNDARY ', ISBW

      else

! no second x-point
! respline to equal arclength
         m11=nteta-1
         hbnd=1.d0/(m11-1)
         t1(1)=0.d0
         do j=2,m11-1
            t1(j)=t1(j-1)+hbnd
         enddo
         t1(m11)=1.d0
         call ARCM(M ,UK1,VK1,T1N)
         call SPLNA1(NTZ,M ,T1N,UK1,M11,T1,NSZ,W1,W2,W3)
         call SPLNA1(NTZ,M ,T1N,VK1,M11,T1,NSZ,W1,W2,W3)

! ARCLENGTH MESH FOR NEW BO. (ISBO ITERATIONS)
         ISBO=20
         EPSBO=1.d-10
         do ISB=1,ISBO
            ISBW=ISB
            call ARCM(M11 ,UK1,VK1,T1N)
            call SPLNA1(NTZ,M11 ,T1N,UK1,M11,T1,NSZ,W1,W2,W3)
            call SPLNA1(NTZ,M11 ,T1N,VK1,M11,T1,NSZ,W1,W2,W3)
            ABO=0.d0
            do J=1,M11
              ABO=MAX(ABO,ABS(T1N(J)-T1(J)))
            enddo
            if (ABO.LE.EPSBO) EXIT
         enddo

      endif

      nbtab=m11
      deallocate( rbtab, zbtab )
      allocate( rbtab(nbtab), zbtab(nbtab) )

      do j=1,m11
         if(clock.gt.0.) then 
            rbtab(j)=uk1(j)
            zbtab(j)=vk1(j)
         else
            rbtab(j)=uk1(m11-j+1)
            zbtab(j)=vk1(m11-j+1)
         endif
      enddo

      return
      end subroutine arc_x_bnd

!----------------------------------------------------------------
      subroutine ARCM(M, US, VS, AW)
! Arclength mesh for given curve

      implicit none

      integer, intent(in) :: M
      real*8, intent(in), dimension(M) :: US, VS
      real*8, intent(out), dimension(M) :: AW

      integer :: j
      real*8 :: AWJ, RWM
! A.L. computation
      AW(1)=0.d0
      do J=2,M
         AWJ=SQRT( (US(J)-US(J-1))**2+(VS(J)-VS(J-1))**2 )
         AW(J)=AW(J-1)+AWJ
      enddo

      RWM=1.d0/AW(M)
      do J=1,M
         AW(J)=AW(J)*RWM
      enddo

      return
      end subroutine ARCM

!----------------------------------------------------------------
      subroutine SPLNA1(NAZ, N, S, FU, NK, SK, NSZ, W, WW, WWW)
! Spline reconstruction for FU given new mesh SK(NK)

      implicit none

      integer, intent(in) :: NAZ, N, NK, NSZ
      real*8, intent(in), dimension(NAZ) :: S
      real*8, intent(in), dimension(NSZ) :: W, WW, WWW
      real*8, intent(out), dimension(NAZ) :: FU
      real*8, intent(inout), dimension(NAZ) :: SK

      integer :: i, N4, IFAIL
      real*8, dimension(4) :: CW(1:4)

      N4=N+4

! Match ends
      SK(1)=S(1)
      SK(NK)=S(N)

      call E01BAF(N, S, FU, W, WW, NSZ, WWW, NSZ, IFAIL)
      do I=1, NK
         call E02BCF(N4, W, WW, SK(I), 0, CW, IFAIL)
         FU(I)=CW(1)
      enddo

      return
      end subroutine SPLNA1
