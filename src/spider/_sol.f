      subroutine f_solve(isol, wdm)

!--ISOL-->if ISOL equal 1 reodering and factorization are asumed
!--       to be done previosly.(PATH=3 will be used in SDRVD)
!--
!--       if ISOL equal 0 reodering and factorization will be done
!--       (PATH=1 will be used in ODRVD and SDRVD)
!--
!--WDM    two dimensional grid array containing rezult of solving

      use sp_parameters

      implicit none

      integer, intent(in) :: isol
      real*8, intent(out) :: wdm(nrp, ntp)

      integer :: i, j, ieq, 
     &   p(neqp), ip(neqp), isp(nspp), ipath, flag, esp
      real*8 zw(neqp), rsp(nspp)
      integer, external :: numlin

      include 'compol.inc'
      include 'compol_add.inc'

      common /comwrc/ rsp,p,ip

      equivalence (rsp(1),isp(1))

      if( itin/nitdel*nitdel+nitbeg .eq. itin 
     * 		.OR. itin.lt.nitbeg ) then
         ipath=3
         if(isol.eq.0) then
            call odrvd(neq,ia,ja,a,p,ip,nspp,isp,1,flag)
            ipath=1
         endif
         call sdrvd(neq,p,ip,ia,ja,a,right,zw,nspp,
     *              isp,rsp,esp,ipath,flag)

      else
         call solvit(isol,zw)
      endif

      do i=1,Nr1
         do j=2,Nt1
            ieq=numlin(i,j,nr,nt)
            wdm(i,j)=zw(ieq)
         enddo
      enddo

      do i=1,Nr
         wdm(i,1)=wdm(i,nt1)
         wdm(i,nt)=wdm(i,2)
      enddo

      return
      end subroutine f_solve

!----------------------------------------------------------------
      subroutine solvit(isol, zw)

      use sp_parameters

      implicit none

      integer, intent(in) :: isol
      real*8, intent(out), dimension(neqp) :: zw

      integer :: i, j, i1, i2, ieq, ic, il, im, icp(neqp), ip(neqp), 
     &   isp(nspp), ipath, flag, esp
      real*8 :: znes, p
      real*8, dimension(nspp) :: rsp
      real*8, dimension(neqp) :: zyy, wpp, wzz, wrr, zuu
      integer, external :: numlin

      include 'compol.inc'
      include 'compol_add.inc'

      common /comwrc/ rsp,p,ip

      equivalence (rsp(1),isp(1))

      do il=1,neq
         wrr(il)=right(il)
      enddo

      if(isol.eq.0) then
         do i=1,Nr1
            do j=2,Nt1
               ieq=numlin(i,j,nr,nt)
               zw(ieq)=g(i,j)
            enddo
         enddo
      else if(isol.eq.1) then
         do i=1,Nr1
            do j=2,Nt1
               ieq=numlin(i,j,nr,nt)          
               zw(ieq)=psii(i,j)
            enddo
         enddo
      endif

      do il=1,neq
         i1=ia(il)
         i2=ia(il+1)-1
         znes=0.d0
         do im=i1,i2
            ic=ja(im)
            znes=znes+daop(im)*zw(ic)
         enddo
         zyy(il)=wrr(il)-znes
      enddo

      call sdrvd(neq,p,ip,ia,ja,aop0,zyy,zw,nspp,
     *           isp,rsp,esp,3,flag)

      return
      end subroutine solvit

!----------------------------------------------------------------
      subroutine solext

      use sp_parameters

      implicit none

      integer :: i, j, ieq, p1(neqp),ip1(neqp),isp1(nspp),
     &   ipath, flag, esp
      real*8 :: zw(neqp), rsp1(nspp)
      integer, external :: numlin

      include 'compol.inc'
      include 'compol_add.inc'

      common /comwrp/ rsp1,p1,ip1

      equivalence (rsp1(1),isp1(1))

      if( itin/nitdel*nitdel+nitbeg .eq. itin 
     * 		.OR. itin.lt.nitbeg ) then
         call odrvd(neqpla,ia,ja,a,p1,ip1,nspp,isp1,1,flag)
         ipath=1
         call sdrvd(neqpla,p1,ip1,ia,ja,a,right,zw,nspp,
     *              isp1,rsp1,esp,ipath,flag)
      else
         call soleit(zw)
      endif

      do i=1,iplas-1
         do j=2,Nt1
            ieq=numlin(i,j,nr,nt)
            psie(i,j)=zw(ieq)
         enddo
      enddo

      do i=1,Nr
         psie(i,1)=psie(i,nt1)
         psie(i,nt)=psie(i,2)
      enddo

      return
      end subroutine solext

!----------------------------------------------------------------
      subroutine soleit(zw)

      use sp_parameters

      implicit none

      real*8, intent(out), dimension(neqp) :: zw

      integer :: i, j, i1, i2, ic, ieq, il, im, p1(neqp), ip1(neqp), 
     &   isp1(nspp), ipath, flag, esp
      real*8 :: znes
      real*8, dimension(nspp) :: rsp1
      real*8, dimension(neqp) :: zyy, wpp, wzz, wrr, zuu
      integer, external :: numlin

      include 'compol.inc'
      include 'compol_add.inc'

      common /comwrp/ rsp1, p1, ip1

      equivalence (rsp1(1),isp1(1))

      do i=1, Nr1
         do j=2, Nt1
            ieq=numlin(i, j, nr, nt)
            zw(ieq)=psie(i, j)
         enddo
      enddo

      do il=1, neqpla
         wrr(il)=right(il)
      enddo

      do il=1, neqpla
         i1=ia(il)
         i2=ia(il+1)-1
         znes=0.d0
         do im=i1,i2
            ic=ja(im)
            znes=znes+dapp(im)*zw(ic)
         enddo
         zyy(il)=wrr(il)-znes
      enddo

      call sdrvd(neqpla, p1, ip1, ia, ja, app0, zyy, zw, nspp,
     *           isp1, rsp1, esp, 3, flag)

      return
      end subroutine soleit
