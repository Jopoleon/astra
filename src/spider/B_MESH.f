      subroutine qst(qcen, cnor, b0ax, r0ax)

      use sp_parameters, only: twopi

      implicit none

      integer, parameter :: nshp=10, n_tp=128, n_rp=64

      real*8, intent(in) :: qcen, cnor, b0ax, r0ax

      integer :: i, j, iplas, nt, iplas1, nt1
      real*8 :: Drr, Drz, Dzz, Dxx, Dyy, tg2a, cos2a, sin2a, bfcen,
     &   fcen, fvac, ffp, dpsi, ps14, 
     &   r0, r1, r2, r3, r4, z1, z2, z3, z4, rm, zm,
     &   sqpol, q2pi, psim, psip, helinp, helout
      real*8 :: funsq, tabf
      real*8, dimension(5) :: dp
      real*8, dimension(n_rp) :: psia, psi, fplr, qplr, dflufi, fluxfi
      real*8, dimension(n_rp, n_tp) :: rplr, zplr

      common /complr/ rplr, zplr, dp, psia, psi, psim, psip, rm, zm, 
     &                fplr, qplr, dflufi, fluxfi, 
     &                iplas, nt, iplas1, nt1
      common /comhel/ helinp, helout

      Drr=dp(3)
      Drz=dp(4)
      Dzz=dp(5)

      tg2a=2.d0*drz/(drr-dzz)
      cos2a=1.d0/sqrt(1.d0+tg2a**2)
      sin2a=cos2a*tg2a

      Dxx=0.5d0*(Drr+Dzz)+0.5d0*cos2a*(Drr-Dzz)+sin2a*Drz
      Dyy=0.5d0*(Drr+Dzz)+0.5d0*cos2a*(Dzz-Drr)-sin2a*Drz

      bfcen=qcen*dsqrt(dxx*dyy)

      fcen=bfcen*rm
      fvac=b0ax*r0ax

      dpsi=(psia(iplas)-psia(iplas-1))
      ps14=-0.25d0*dpsi
      ffp=tabf(ps14)
      fplr(iplas-1)=sqrt(fvac**2-ffp*cnor*(psim-psip)*dpsi)

      do i=iplas-2,1,-1
         dpsi=0.5d0*(psia(i+2)-psia(i))*(psim-psip)
         ffp=tabf(psia(i+1))
         fplr(i)=sqrt(fplr(i+1)**2-2.d0*ffp*dpsi)
      enddo

      fplr(iplas)=fvac

      do i=1,iplas1

         dflufi(i)=0.

         do j=2,nt1
            r1=rplr(i,j)
            r2=rplr(i+1,j)
            r3=rplr(i+1,j+1)
            r4=rplr(i,j+1)

            z1=zplr(i,j)
            z2=zplr(i+1,j)
            z3=zplr(i+1,j+1)
            z4=zplr(i,j+1)

            r0=(r1+r2+r3+r4)*0.25d0
            sqpol=funsq(r1,r2,r3,r4,z1,z2,z3,z4)
            dflufi(i)=dflufi(i)+sqpol*fplr(i)/r0
         enddo

         q2pi=-dflufi(i)/((psia(i+1)-psia(i))*(psim-psip))
         qplr(i)=q2pi/twopi

      enddo

      fluxfi(1)=0.d0
      do i=2,iplas
         fluxfi(i)=fluxfi(i-1) + dflufi(i-1)
      enddo

      helout = - fluxfi(iplas)*psip  
      do i=1,iplas1
         helout = helout + 0.5d0*( psi(i)+psi(i+1) )*dflufi(i) 
      enddo

      return
      end subroutine qst

!----------------------------------------------------------------
      subroutine gridpl

      use sp_parameters, only: twopi
      use comblc, only: um, up, rm, zm, r, z, u, imax, jmax

      implicit none

      integer, parameter :: nshp=10, n_tp=128, n_rp=64

      integer :: i, j, k, l, ii, jj, ibeg, iplas, nt, iplas1, nt1, nsh,
     &   ixp1, ixp2, jxp1, jxp2
      real*8 :: psim, psip, xm, ym, dtet, stpx, stpy, ttj, ro2j, u0, ur0
      real*8, dimension(n_tp) :: ron, teta
      real*8, dimension(nshp) :: xs, ys, fun
      real*8, dimension(5) :: dp
      real*8, dimension(n_rp) :: psia, psi, fplr, qplr, dflufi, fluxfi
      real*8, dimension(n_rp, n_tp) :: rplr, zplr

      common /complr/ rplr, zplr, dp, psia, psi, psim, psip, xm, ym,  
     *             fplr, qplr, dflufi, fluxfi,
     *             iplas, nt, iplas1, nt1

      iplas  = 61
      nt     = 95
      iplas1 = iplas-1
      nt1    = nt-1

      ixp1=0
      jxp1=0

      ixp2=0
      jxp2=0

      psim=um
      psip=up
      xm=rm
      ym=zm

      dtet=twopi/(nt-2)
      teta(1)=-dtet

      do j=2,nt
         teta(j) = teta(j-1) + dtet
      enddo

      teta(1)=teta(nt1)
      teta(nt)=teta(2)

      psia(1)=1.d0

      nsh  = 1
      xs(nsh) = r(imax)
      ys(nsh) = z(jmax)
      fun(nsh) = u(imax,jmax)

      do k=-1,1
         ii = imax + k
         do l=-1,1
            jj = jmax + l
            if(ii.ne.imax .OR. jj.ne.jmax) then
               nsh  = nsh+1
               xs(nsh) = r(ii)
               ys(nsh) = z(jj)
               fun(nsh) = u(ii,jj)
	    endif
         enddo
      enddo

      call deriv5(xs,ys,fun,nsh,5,dp)
      stpx=(r(imax+1)-r(imax))*1.0
      stpy=(z(jmax+1)-z(jmax))*0.7

      do i=2,iplas
         u0=1.d0-((i-1)/(iplas-1.))  !  **2
         psia(i)=u0
         ur0=up+u0*(um-up)
         do j=1,nt
            ttj=teta(j) 
            ro2j=(ur0-um)/( 0.5*dp(3)*dcos(ttj)**2 +
     &                        dp(4)*dcos(ttj)*dsin(ttj) +
     &                    0.5*dp(5)*dsin(ttj)**2 )
            ron(j)=dsqrt(ro2j)     
            rplr(i,j)=rm+ron(j)*dcos(teta(j))
            zplr(i,j)=zm+ron(j)*dsin(teta(j))
         enddo

	 ibeg=i+1
	 if( ron(2) .GT. dmax1(stpx,stpy) ) EXIT

      enddo

      do i=ibeg,iplas
         u0=(iplas-i)/(iplas-1.d0)
         psia(i)=u0

         call loop95(teta,nt,ron,u0)

         do j=1,nt
            rplr(i,j)=rm+ron(j)*dcos(teta(j))
            zplr(i,j)=zm+ron(j)*dsin(teta(j))
         enddo

      enddo

      do j=1,nt
         rplr(1,j) = rm
         zplr(1,j) = zm
      enddo

      do i=1,iplas
         rplr(i,1)=rplr(i,nt1)
         zplr(i,1)=zplr(i,nt1)
         rplr(i,nt)=rplr(i,2)
         zplr(i,nt)=zplr(i,2)
      enddo

      teta(1)=teta(nt1)-twopi
      teta(nt)=teta(2)+twopi

      do i=1,iplas
         do j=1,nt
            psi(i)=psip+psia(i)*(psim-psip)
         enddo
      enddo

      return
      end subroutine gridpl

!----------------------------------------------------------------
      subroutine loop95(tetpol, ntet, ro0, u0)

      use sp_parameters, only: pi, nbndp2, nbndp4, nbndp6, nip, njp 
      use comblc, only: imax, jmax, rm, r, zm, z, ni, nj, un, ni1

      implicit none

      integer, intent(in) :: ntet
      real*8, intent(in) :: tetpol(1)
      real*8, intent(out) :: ro0(1), u0

      integer :: i, j, ic, ic1, jc, jc1, ig, ii, jj, icell, jcell,
     &   imax1, jmax1, lin, nxb, ifail
      real*8 :: drx, dzx , deltet, tetp
      real*8, dimension(nbndp2) :: rxb, zxb, roxb, tetxb
      real*8, dimension(nip, njp) :: ut
      real*8 RRK(nbndp4), CCK(nbndp4), WRK(nbndp6)
      real*8 CWK(4)
      real*8 :: xzer

      if(u0.lt.2.d-5) u0=2.d-5

      imax1=imax-1
      jmax1=jmax-1

      do i=imax1,imax
         if(rm.le.r(i+1) .AND. rm.gt.r(i)) icell=i
      enddo

      do j=jmax1,jmax
         if(zm.le.z(j+1) .AND. zm.gt.z(j)) jcell=j
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
         if(ut(i,j)*ut(i+1,j).le.0.) then
            ic=i
            jc=j
            rxb(ig)=xzer(r(i+1),r(i),ut(i+1,j),ut(i,j))
            zxb(ig)=z(j)
            go to 886
         endif
      enddo
      stop

 886  continue

      ic1=ic
      jc1=jc
      lin=1

      do

         ig=ig+1
         i=ic
         j=jc

         if(ut(i+1,j)*ut(i,j).le.0. .AND. lin.ne.1) then
            ic=i
            jc=j-1
            lin=3
            rxb(ig)=xzer(r(i+1),r(i),ut(i+1,j),ut(i,j))
            zxb(ig)=z(j)
         elseif(ut(i+1,j+1)*ut(i+1,j).le.0..AND. lin.ne.2) then
            ic=i+1
            jc=j
            lin=4
            rxb(ig)=r(i+1)
            zxb(ig)=xzer(z(j+1),z(j),ut(i+1,j+1),ut(i+1,j))
         elseif(ut(i+1,j+1)*ut(i,j+1).le.0. .AND. lin.ne.3) then
            ic=i
            jc=j+1
            lin=1
            rxb(ig)=xzer(r(i+1),r(i),ut(i+1,j+1),ut(i,j+1))
            zxb(ig)=z(j+1)
         elseif(ut(i,j+1)*ut(i,j).le.0..AND. lin.ne.4) then
            ic=i-1
            jc=j
            lin=2
            rxb(ig)=r(i)
            zxb(ig)=xzer(z(j+1),z(j),ut(i,j+1),ut(i,j))
         endif

         if(jc.eq.jc1 .AND. ic.eq.ic1) then
            ig=ig+1
            rxb(ig)=rxb(2 )
            zxb(ig)=zxb(2 )
            EXIT
         endif

      enddo

      nxb=ig

      do ig=1,nxb

         drx=rxb(ig)-rm
         dzx=zxb(ig)-zm
         tetp=datan(dzx/drx)
         if(ig.ne.1) then
            if(drx.lt.0.) tetp=tetp+pi
            deltet=tetp-tetxb(ig-1)
            if(deltet.lt.0.) tetp=tetp+2.*pi
         endif

         tetxb(ig)=tetp
         roxb(ig)=dsqrt(drx**2+dzx**2)

      enddo

      call E01BAF(Nxb,tetxb,roxb,RRK,CCK,nxb+4,WRK,6*nxb+16,IFAIL)

      do j=1,Ntet
         tetp=tetpol(j)
         if(tetp.lt.tetxb(1)) then
            tetp=tetp+2.*pi
         elseif(tetp.gt.tetxb(nxb)) then
            tetp=tetp-2.*pi
         endif
         CALL E02BCF(Nxb+4, RRK, CCK, tetp, 0, CWk, IFAIL)
         ro0(j)=cwk(1)
      enddo

      return
      end subroutine loop95
