      subroutine solint(imov)

      use sp_parameters, only: neqp, nrp, nspp
      use keys, only: kpr
      use compol, only: itin, neqpla, a, nr, nt, nt1, iplas, iplas1,
     &     ia, ja, q, iter, iswtch, ngav, psi, psia, right, errm,
     & nitbeg, nitdel

      implicit none

      integer, intent(in) :: imov

      integer :: i, j, ieq, ipath, flag, esp, numlin
      integer :: p1(neqp), ip1(neqp), isp1(nspp)
      real*8 :: errpss, delpsi
      real*8 :: zw(neqp), rsp1(nspp), wght(nrp)

      equivalence (rsp1(1),isp1(1))

      common /comwrp/ rsp1, p1, ip1

      if(imov.eq.0) then
         if(itin.eq.1) then
            ipath=1
            call odrvd(neqpla,ia,ja,a,p1,ip1,nspp,isp1,1,flag)
         else
            ipath=3
         endif
         call sdrvd(neqpla,p1,ip1,ia,ja,a,right,zw,nspp,
     &              isp1,rsp1,esp,ipath,flag)
      else

         if( itin/nitdel*nitdel+nitbeg .eq. itin 
     & 		.OR. itin.lt.nitbeg 
     & 		.OR. errm.gt.7.7d-2 ) then

            call odrvd(neqpla,ia,ja,a,p1,ip1,nspp,isp1,1,flag)
            ipath=1

            call sdrvd(neqpla,p1,ip1,ia,ja,a,right,zw,nspp,
     &                 isp1,rsp1,esp,ipath,flag)

            if(kpr.eq.1) then
               write(6,*) 'sdrv: flag,esp',flag,esp
            endif
         else
            call solbit(zw)
         endif

      endif

      errpss=0.d0
      do i=1,iplas-1
         if(ngav.eq.0.or.i.eq.1) then
            wght(i)=0.5d0
         else
            wght(i)=1.d0/(2.0d0+sqrt(abs(1.d0-psia(i)))*q(i)/q(1))
            if(iswtch.eq.1)  wght(i)=0.5d0  
         endif

         if(iter.eq.1) wght(i)=1.0d0
      enddo

      do i=1,iplas-1
         do j=2,Nt1
            ieq=numlin(i,j,nr,nt)
            delpsi=dabs(psi(i,j)-zw(ieq))
            errpss=dmax1(errpss,delpsi)
            psi(i,j)=zw(ieq)*wght(i)+(1.d0-wght(i))*psi(i,j)
         enddo
      enddo

      do i=1,iplas1
         psi(i,1)=psi(i,nt1)
         psi(i,nt)=psi(i,2)
      enddo

      return
      end subroutine solint

!----------------------------------------------------------------
      subroutine solbit(zw)

      use sp_parameters, only: neqp, nspp
      use compol, only: nr, nr1, nt, nt1, ia, ja, dapp, app0, right,
     & neqpla, psi

      implicit none

      real*8, intent(out), dimension(neqp) :: zw

      integer :: i, j, i1, i2, ic, ieq, il, im, flag, esp
      integer :: p1(neqp), ip1(neqp), isp1(nspp)
      real*8 :: znes
      real*8, dimension(nspp) :: rsp1
      real*8, dimension(neqp) :: zyy, wrr
      integer, external :: numlin

      equivalence (rsp1(1), isp1(1))

      common /comwrp/ rsp1,p1,ip1

      do i=1,Nr1
         do j=2,Nt1
            ieq=numlin(i,j,nr,nt)
            zw(ieq)=psi(i,j)
         enddo
      enddo

      do il=1,neqpla
         wrr(il)=right(il)
      enddo			     

      do il=1,neqpla
         i1=ia(il)
         i2=ia(il+1)-1
         znes=0.d0
         do im=i1,i2
            ic=ja(im)
            znes=znes+dapp(im)*zw(ic)
         enddo
         zyy(il)=wrr(il)-znes
      enddo

      call sdrvd(neqpla,p1,ip1,ia,ja,app0,zyy,zw,nspp,
     *           isp1,rsp1,esp,3,flag)

      return
      end subroutine solbit
