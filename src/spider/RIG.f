      subroutine bound

      use sp_parameters, only: nbndp
      use comblc, only: ni, ni1, nj, nj1, nbnd,
     &   binadg, dgdn, g, ui, dr, dz, r, r12

      implicit none

      integer :: i, j, ib, ibc
      real*8 :: zpsi
      real*8, dimension(nbndp) :: psib

! array dgdn initialization (calculation dg/dn)

      ib=1

      dgdn(ib)=0.d0
      g(1,1)=0.d0

      do i=2,ni1
         ib=ib+1
         dgdn(ib) = g(i,2)/dz(1)/r(i)
         g(i,1)=0.d0
      enddo

      ib=ib+1
      dgdn(ib)=0.d0
      g(ni,1)=0.d0

      do j=2,nj1
         ib=ib+1
         dgdn(ib) = g(ni1,j)/( dr(ni1)*r12(ni1))
         g(ni,j)=0.d0
      enddo

      ib=ib+1
      dgdn(ib)=0.d0
      g(ni,nj)=0.d0

      do i=ni1,2,-1
         ib=ib+1
         dgdn(ib) = g(i,nj1)/dz(nj1)/r(i)
         g(i,nj)=0.d0
      enddo

      ib=ib+1
      dgdn(ib)=0.d0
      g(1,nj)=0.d0

      do j=nj1,2,-1
         ib=ib+1
         dgdn(ib) = g(2,j)/(dr(1)*r12(1))
         g(1,j)=0.d0
      enddo

      ib=ib+1
      dgdn(ib)=0.d0

      do ib=1,nbnd
         dgdn(ib)= (dgdn(ib)+dgdn(ib+1))*0.5
      enddo

      do ib=1,nbnd
         zpsi=0.d0
         do ibc=1,nbnd
            zpsi= zpsi+dgdn(ibc)*binadg(ibc,ib)
         enddo
         psib(ib)=zpsi
      enddo

! boundary condition for ui

      ib=0
      do i=1,ni
         ib=ib+1
         ui(i,1)=psib(ib)
      enddo

      do j=2,nj
         ib=ib+1
         ui(ni,j)=psib(ib)
      enddo

      do i=ni1,1,-1
         ib=ib+1
         ui(i,nj)=psib(ib)
      enddo

      do j=nj1,2,-1
         ib=ib+1
         ui(1,j)=psib(ib)
      enddo

      return
      end subroutine bound

!----------------------------------------------------------------
      subroutine right0(ill, jll, icelm, jcelm, ngav1)

      use sp_parameters
      use keys, only: kpr, kxwx, key_plc
      use comblc

      implicit none

      integer, parameter :: nshp=10

      integer, intent(in) :: ill, jll, icelm, jcelm, ngav1

      integer :: i, j, k, l, ic, jc, ix, jx, imp, jmp, il, is, nsh, 
     &   nctrli, kodex1, kodex2, kodxp, iprlm, nulim1, nulim2, neigb,
     &   ilm, ilma, jlm, jlma, icelma, jcelma, llim, isum, nloc,
     &   ipr1, ipr2, ipr3, ipr4, i_rus, iprij
      real*8 :: omg, sigm, helinp, helout, tokp, sij, sij4, wght,
     &   umax, uem, uen, uumm, uxold, uix, ublm, unold, ups, upp,
     &   ulmax1, ulmax2, uartm,
     &   rrmm, zzmm, rmold, zmold, rix, zjx, 
     &   rlim1, rlim2, rmlim1, rmlim2, frlim,
     &   delcr, delcz, delc, delun, zver, zvmin, zvmax, xi, psi,
     &   cur1, cur2, cur3, cur4, funcur, bcur,
     &   x0, x1, x2, x3, x4, x12, x23, x34, x14,
     &   z0, z1, z2, z3, z4, z12, z23, z34, z14,
     &   u0, u1, u2, u3, u4, u12, u23, u34, u14
      real*8, dimension(nshp) :: xs, ys, fun
      real*8, dimension(5) :: dp, dp1, dp2
      real*8, dimension(6) :: clsq
      real*8, dimension(nip, njp) :: curn
      integer, external :: nlin
      real*8, external :: bline, blint

      common /comomg/ omg, sigm
      common /comkjf/ curn, nctrli
      common /comhel/ helinp, helout

      nctrli=nctrl

! right hand side for eq. Lg=Jf in rectangular box.

      do il=1,neqp
         right(il)=0.d0
      enddo
!  definition of Umax- maximum poloidal flux function value on
!   grid in plasma.

      umax=u(imax,jmax)

      if(iter.ge.iterbf) then
         imax=ill
         jmax=jll
         uem=bline(icelm,jcelm,rl,zl)
      else
         do
            ic=imax
            jc=jmax
            do k=-1,1
               imp= imax+k
               do l=-1,1
                  jmp= jmax+l
                  if( u(imp,jmp).gt.umax) then
                     ic=imp
                     jc=jmp
                     umax=u(imp,jmp)
                  endif
               enddo
            enddo
            if(ic.ne.imax .OR. jc.ne.jmax) then
               imax=ic
               jmax=jc
               umax=u(imax,jmax) ! and stay in the loop
            else
               EXIT
            endif
         enddo
      endif

      do

         if(kpr.eq.1)write(*,*) 'imax jmax umax',imax,jmax,umax
         if(kpr.eq.1)write(*,*) 'iclm jclm     ',icelm,jcelm
         if (isnan(umax)) then
            write(*,*) 'nan umax'
	 call err_catch_a
         endif

! definition of Um - poloidal flux function value at magn. axis.

         nsh=1
         xs(nsh)=r(imax)
         ys(nsh)=z(jmax)
         fun(nsh)=u(imax,jmax)

         do k=-1,1
            i= imax+k
            do l=-1,1
               j= jmax+l
               if(i.eq.imax .AND. j.eq.jmax) CYCLE
               nsh=nsh+1
               xs(nsh)=r(i)
               ys(nsh)=z(j)
               fun(nsh)=u(i,j)
            enddo
         enddo

         call lsq_sur6(xs,ys,fun,nsh,clsq,rrmm,zzmm,uumm,dp)
         rmold=rm
         zmold=zm
         rm=rrmm
         zm=zzmm
         um=uumm
         errm=dsqrt( (rm-rmold)**2+(zm-zmold)**2 )

         if(iter.lt.iterbf) then
            clr=0.d0
            clz=0.d0
	 rl=rm
            zl=zm
            call shab(ilma,jlma,icelma,jcelma)
            uem=bline(icelma,jcelma,rl,zl)
            EXIT
         endif

! stabilization by artifical field

         delcr= -(     dp(3)*(rl-rm) + dp(4)*(zl-zm) )*0.5d0/rl
         delcz= -(     dp(5)*(zl-zm) + dp(4)*(rl-rm) )
         clr=clr+delcr
         clz=clz+delcz

         do i=1,ni
            do j=1,nj
               u(i,j)=u(i,j)+delcz*z(j)+delcr*r(i)*r(i)
            enddo
         enddo

         delc = dabs(delcr)+dabs(delcz)
         if(delc.lt.1.d-9) EXIT

      enddo

      if(kpr.eq.1)write(*,*) '*****rm,zm,rl,zl',rm,zm,rl,zl
      if(kpr.eq.1)write(*,*) '*****clr,clz ',clr,clz

!----------------------------------------------------------
! Definition ux- poloidal flux function value on separatrix

      uxold=ux0
      if(nctrl.lt.0 ) then
         rx1=rx10
         zx1=zx10
         rx2=rx20
         zx2=zx20
      endif

      call xpoint(ux1,rx1,zx1,ix1,jx1,dp1,1,kodex1)
      if(kodex1.gt.0) then
         rx1=rx10
         zx1=zx10
      endif
      if(kpr.eq.1)write(*,*) 'xpoint1',kodex1

      call xpoint(ux2,rx2,zx2,ix2,jx2,dp2,2,kodex2)
      if(kodex2.gt.0) then
         rx2=rx20
         zx2=zx20
      endif
      if(kpr.eq.1)write(*,*) 'xpoint2',kodex2

      kodxp=kodex1*kodex2
      if(kodxp.ne.0) then
         if(kpr.eq.1)write(*,*) ' all x-points out of box'
         if(kpr.eq.1)write(*,*) ' only limiter case '
         nctrli=1
      endif
      if(kpr.eq.1)then
         write(*,*)'rx1,zx1,ux1',rx1,zx1,ux1
         write(*,*)'rx2,zx2,ux2',rx2,zx2,ux2
         write(*,*) 'ux1,ux2',ux1,ux2
      endif
      if(ux1.gt.ux2) then
         if(kpr.eq.1) write(*,*) 'xpoint=1'
         ux0= ux1
         rx0= rx1
         zx0= zx1
         ix = ix1
         jx = jx1
         rix = r(ix)
         zjx = z(jx)
         uix = u(ix,jx)
         do is=1,5
            dp(is)=dp1(is)
         enddo
      else
         if(kpr.eq.1) write(*,*) 'xpoint=2'
         ux0= ux2
         rx0= rx2
         zx0= zx2
         ix = ix2
         jx = jx2
         rix = r(ix)
         zjx = z(jx)
         uix = u(ix,jx)
         do is=1,5
            dp(is)=dp2(is)
         enddo
      endif

      if(uxold.ge.ux0) ux0=(1.d0-sigm)*uxold+sigm*ux0
      zvmax=dmax1(zx1,zx2)
      zvmin=dmin1(zx1,zx2)
      zver=dabs(zm-zx0)

! definition up- poloidal flux function value on plasma boundary.

      if(kodxp.eq.0) then
      ups=um - alp*(um-ux0)
      else       !!!!<- no x-points
         ups=-1.d12
      endif

      if(nctrli.eq.0 ) then
         up=um - alp*(um-ux0)
         alpnew=alp
         numlim=0
      elseif(nctrli.eq.-2 .AnD. iter.gt.3 ) then
         rx0= rx1
         zx0= zx1
         ux0= ux1
         upp=um - alp*(um-ux0)
         up=dmax1(psi_bon,upp)
         alpnew=(um-up)/(um-ux0)
         kodex2=1
      elseif(nctrli.lt.0 ) then
         upp=um - alp*(um-ux0)
         if(psi_bon.lt.ux0) then
		  up=upp
            alpnew=alp
         else
            up=psi_bon
            alpnew=(um-up)/(um-ux0)
         endif
         numlim=0
      else   ! if(nctrli.gt.0 ) <-limiter case
         ulmax1=ups
         ulmax2=ups
         nulim1=0
         nulim2=0
         neigb=0
         do llim=1,nblm
            if(kodex1.eq.0) then
               rlim1=frlim(dp1,zblm(llim),rx1,zx1,rm,zm)
               rmlim1=frlim(dp1,zm,rx1,zx1,rm,zm)
            else
               rlim1=rmax
               rmlim1=rmax
            endif

            if(kodex2.eq.0) then
               rlim2=frlim(dp2,zblm(llim),rx2,zx2,rm,zm)
               rmlim2=frlim(dp2,zm,rx2,zx2,rm,zm)
            else
               rlim2=rmax
               rmlim2=rmax
            endif

            if( (rblm(llim)-rlim1)*(rm-rmlim1).le.0.d0) CYCLE
            if( (rblm(llim)-rlim2)*(rm-rmlim2).le.0.d0) CYCLE

            ilm=iblm(llim)
            jlm=jblm(llim)
            ublm=blint(iblm(llim),jblm(llim),rblm(llim),zblm(llim))
            iprlm=ipr(ilm,jlm)+ipr(ilm+1,jlm)
     &              +ipr(ilm,jlm+1)+ipr(ilm+1,jlm+1) 
	          if(iprlm.gt.0) then
               if(ublm.gt.ulmax1)then
                  ulmax1=ublm
                  nulim1=llim
                  neigb=1
               endif
            else
               if(ublm.gt.ulmax2)then
                  ulmax2=ublm
                  nulim2=llim
               endif
            endif

         enddo

         if(neigb.eq.1) then
            ublmax=ulmax1
            numlim=nulim1
         elseif(neigb.eq.0) then
            ublmax=ulmax2
            numlim=nulim2
         endif

         up=dmax1(ups,ublmax)
         if(kodxp.eq.0) alpnew=(um-up)/(um-ux0)

      endif

      if(kpr.eq.1) then
         write(*,*)'rx0,zx0,ux0',rx0,zx0,ux0
         write(*,*) 'up,numlim',up,numlim,neigb
         write(*,*) 'psi_bon',psi_bon
         write(*,*) 'um-up',um-up
         write(*,*) 'alpnew',alpnew
      endif

!  un- normal poloidal flux.

      erru=0.d0

      do i=1,ni
         do j=1,nj
            unold=un(i,j)
            un(i,j)=(u(i,j)-up)/(um-up)
            delun=dabs(un(i,j)-unold)
            erru=dmax1(delun,erru)
         enddo
      enddo

      if(kpr.eq.1) write(*,*) 'erru', erru

!  dimension ipr initialization.
!  if ipr(i,j)=1,then point (i,j) is in plasma.

      do i=1,ni
         do j=1,nj
            ipr(i,j)=0
         enddo
      enddo
      ipr(imax-1,jmax-1)=1

      do j=jmax,nj1
         if(kodex1.eq.0) then
            rlim1=frlim(dp1,z(j),rx1,zx1,rm,zm)
            rmlim1=frlim(dp1,zm,rx1,zx1,rm,zm)
         else
            rlim1=rmax
            rmlim1=rmax
         endif
         if(kodex2.eq.0) then
            rlim2=frlim(dp2,z(j),rx2,zx2,rm,zm)
            rmlim2=frlim(dp2,zm,rx2,zx2,rm,zm)
         else
            rlim2=rmax
            rmlim2=rmax
         endif
         do i=imax,ni1
            if(un(i,j).lt.0.d0) then
               ipr(i,j)=0
            elseif( (r(i)-rlim1)*(rm-rmlim1).le.0.d0) then
               ipr(i,j)=0
            elseif( (r(i)-rlim2)*(rm-rmlim2).le.0.d0) then
               ipr(i,j)=0
            else
               isum=ipr(i-1,j-1)+ipr(i,j-1)+ipr(i+1,j-1)+ipr(i-1,j)
               if(isum.gt.0) ipr(i,j)=1
            endif
         enddo
      enddo

      do j=jmax,nj1
         if(kodex1.eq.0) then
            rlim1=frlim(dp1,z(j),rx1,zx1,rm,zm)
            rmlim1=frlim(dp1,zm,rx1,zx1,rm,zm)
         else
            rlim1=rmin
            rmlim1=rmin
         endif
         if(kodex2.eq.0) then
            rlim2=frlim(dp2,z(j),rx2,zx2,rm,zm)
            rmlim2=frlim(dp2,zm,rx2,zx2,rm,zm)
         else
            rlim2=rmin
            rmlim2=rmin
         endif
         do i=imax,2,-1
            if(un(i,j).lt.0.d0) then
               ipr(i,j)=0
            elseif( (r(i)-rlim1)*(rm-rmlim1).le.0.d0) then
               ipr(i,j)=0
            elseif( (r(i)-rlim2)*(rm-rmlim2).le.0.d0) then
               ipr(i,j)=0
            else
               isum=ipr(i-1,j-1)+ipr(i,j-1)+ipr(i+1,j-1)+ipr(i+1,j)
               if(isum.gt.0.) ipr(i,j)=1
            endif
         enddo
      enddo

      do j=jmax,2,-1
         if(kodex1.eq.0) then
            rlim1=frlim(dp1,z(j),rx1,zx1,rm,zm)
            rmlim1=frlim(dp1,zm,rx1,zx1,rm,zm)
         else
            rlim1=rmax
            rmlim1=rmax
         endif
         if(kodex2.eq.0) then
            rlim2=frlim(dp2,z(j),rx2,zx2,rm,zm)
            rmlim2=frlim(dp2,zm,rx2,zx2,rm,zm)
         else
            rlim2=rmax
            rmlim2=rmax
         endif
         do i=imax,ni1
            if(un(i,j).lt.0.) then
               ipr(i,j)=0
            elseif( (r(i)-rlim1)*(rm-rmlim1).le.0.d0) then
               ipr(i,j)=0
            elseif( (r(i)-rlim2)*(rm-rmlim2).le.0.d0) then
               ipr(i,j)=0
            else
               isum=ipr(i-1,j+1)+ipr(i,j+1)+ipr(i+1,j+1)+ipr(i-1,j)
               if(isum.gt.0.) ipr(i,j)=1
            endif
         enddo
      enddo

      do j=jmax,2,-1
         if(kodex1.eq.0) then
            rlim1=frlim(dp1,z(j),rx1,zx1,rm,zm)
            rmlim1=frlim(dp1,zm,rx1,zx1,rm,zm)
         else
            rlim1=rmin
            rmlim1=rmin
         endif
         if(kodex2.eq.0) then
            rlim2=frlim(dp2,z(j),rx2,zx2,rm,zm)
            rmlim2=frlim(dp2,zm,rx2,zx2,rm,zm)
         else
            rlim2=rmin
            rmlim2=rmin
         endif
         do i=imax,2,-1
            if(un(i,j).lt.0.d0) then
               ipr(i,j)=0
            elseif( (r(i)-rlim1)*(rm-rmlim1).le.0.) then
               ipr(i,j)=0
            elseif( (r(i)-rlim2)*(rm-rmlim2).le.0.) then
               ipr(i,j)=0
            else
               isum=ipr(i-1,j+1)+ipr(i,j+1)+ipr(i+1,j+1)+ipr(i+1,j)
               if(isum.gt.0) ipr(i,j)=1
            endif
         enddo
      enddo

!-----------------------------------
! Toroidal current density in plasma

      do i=1,ni
         do j=1,nj
            curf(i,j)=0.d0
         enddo
      enddo

      tokp=0.d0

      do i=2,ni2
         x1=r(i)
         x2=r(i+1)
         x3=r(i+1)
         x4=r(i)
         x12=0.5d0*(x1+x2)
         x23=0.5d0*(x2+x3)
         x34=0.5d0*(x3+x4)
         x14=0.5d0*(x1+x4)
         x0=0.5d0*(x12+x34)
         do j=2,nj2
            cur1=0.d0
            cur2=0.d0
            cur3=0.d0
            cur4=0.d0
            if (kodex1.eq.0 .and. i.le.ix1 .and. i.ge.ix1-1 .and.
     &                            j.le.jx1 .and. j.ge.jx1-1) then
               nloc=21
               call xdetal(i,j,ix1,jx1,rx1 ,zx1 ,dp1,nloc,
     &                     cur1,cur2,cur3,cur4)
               goto 9010
            endif

            if (kodex2.eq.0 .and. i.le.ix2 .and. i.ge.ix2-1 .and.
     &                            j.le.jx2 .and. j.ge.jx2-1) then
               nloc=21
               call xdetal(i,j,ix2,jx2,rx2 ,zx2 ,dp2,nloc,
     &                     cur1,cur2,cur3,cur4)
               goto 9010
            endif

            ipr1=ipr(i,j)
            ipr2=ipr(i+1,j)
            ipr3=ipr(i+1,j+1)
            ipr4=ipr(i,j+1)
            isum=ipr1+ipr2+ipr3+ipr4
            if(isum.eq.0) CYCLE

            z1=z(j)
            z2=z(j)
            z3=z(j+1)
            z4=z(j+1)
            z12=0.5d0*(z1+z2)
            z23=0.5d0*(z2+z3)
            z34=0.5d0*(z3+z4)
            z14=0.5d0*(z1+z4)
            z0=0.5d0*(z12+z34)

            u1=un(i,j)
            u2=un(i+1,j)
            u3=un(i+1,j+1)
            u4=un(i,j+1)
            if(ipr1.eq.0 .AND. u1.gt.0.d0) u1=0.d0
            if(ipr2.eq.0 .AND. u2.gt.0.d0) u2=0.d0
            if(ipr3.eq.0 .AND. u3.gt.0.d0) u3=0.d0
            if(ipr4.eq.0 .AND. u4.gt.0.d0) u4=0.d0
            u12=0.5d0*(u1+u2)
            u23=0.5d0*(u2+u3)
            u34=0.5d0*(u3+u4)
            u14=0.5d0*(u1+u4)
            u0=0.5d0*(u12+u34)

            sij4=0.25d0*dr(i)*dz(j)
            cur1=0.d0
            cur2=0.d0
            cur3=0.d0
            cur4=0.d0

            if(isum.eq.4) then
               cur1=sij4*0.25d0*(funcur(x1,u1)+funcur(x12,u12) +
     &                           funcur(x0,u0)+funcur(x14,u14) )
               cur2=sij4*0.25d0*(funcur(x2,u2)+funcur(x23,u23) +
     &                           funcur(x0,u0)+funcur(x12,u12) )
               cur3=sij4*0.25d0*(funcur(x3,u3)+funcur(x34,u34) +
     &                           funcur(x0,u0)+funcur(x23,u23) )
               cur4=sij4*0.25d0*(funcur(x4,u4)+funcur(x14,u14) +
     &                           funcur(x0,u0)+funcur(x34,u34) )
               goto 9010
            endif

! quadr..(i,j)...

            if(ipr1.eq.1) then
               cur1=bcur(x1,x12,x0,x14, z1,z12,z0,z14, u1,u12,u0,u14)
            elseif(u12.gt.0.d0) then
               cur1=bcur(x12,x0,x14,x1, z12,z0,z14,z1, u12,u0,u14,u1)
            elseif(u0.gt.0.d0) then
               cur1=bcur(x0,x14,x1,x12, z0,z14,z1,z12, u0,u14,u1,u12)
            elseif(u14.gt.0.d0) then
               cur1=bcur(x14,x1,x12,x0, z14,z1,z12,z0, u14,u1,u12,u0)
            endif

! quadr..(i+1,j)...

            if(ipr2.eq.1) then
               cur2=bcur(x2,x23,x0,x12, z2,z23,z0,z12, u2,u23,u0,u12)
            elseif(u23.gt.0.d0) then
               cur2=bcur(x23,x0,x12,x2, z23,z0,z12,z2, u23,u0,u12,u2)
            elseif(u0.gt.0.d0) then
               cur2=bcur(x0,x12,x2,x23, z0,z12,z2,z23, u0,u12,u2,u23)
            elseif(u12.gt.0) then
               cur2=bcur(x12,x2,x23,x0, z12,z2,z23,z0, u12,u2,u23,u0)
            endif

! quadr..(i+1,j+1)...

            if(ipr3.eq.1) then
               cur3=bcur(x3,x34,x0,x23, z3,z34,z0,z23, u3,u34,u0,u23)
            elseif(u34.gt.0.d0) then
               cur3=bcur(x34,x0,x23,x3, z34,z0,z23,z3, u34,u0,u23,u3)
            elseif(u0.gt.0.d0) then
               cur3=bcur(x0,x23,x3,x34, z0,z23,z3,z34, u0,u23,u3,u34)
            elseif(u23.gt.0.d0) then
               cur3=bcur(x23,x3,x34,x0, z23,z3,z34,z0, u23,u3,u34,u0)
            endif

! quadr..(i,j+1)...

            if(ipr4.eq.1) then
               cur4=bcur(x4,x14,x0,x34, z4,z14,z0,z34, u4,u14,u0,u34)
            elseif(u14.gt.0.d0) then
               cur4=bcur(x14,x0,x34,x4, z14,z0,z34,z4, u14,u0,u34,u4)
            elseif(u0.gt.0.d0) then
               cur4=bcur(x0,x34,x4,x14, z0,z34,z4,z14, u0,u34,u4,u14)
            elseif(u34.gt.0.d0) then
               cur4=bcur(x34,x4,x14,x0, z34,z4,z14,z0, u34,u4,u14,u0)
            endif

 9010       continue

            tokp=tokp+cur1+cur2+cur3+cur4

            curf(i,j)     = curf(i,j)    + cur1
            curf(i+1,j)   = curf(i+1,j)  + cur2
            curf(i+1,j+1) = curf(i+1,j+1)+ cur3
            curf(i,j+1)   = curf(i,j+1)  + cur4

         enddo
      enddo

! -------------------------------------------------------

      tokn=tok

      if(ngav1.eq.2 .OR. ngav1.eq.3) then
         if(iter.ge.iterbf .OR. nnstpp.ge.1) then
            uartm=clr*rl*rl+clz*zl
            tok=tok*(ucen-uem-uartm)/(um-uem-uartm)
         endif
      endif
                 
      wght=0.25d0

      if(ngav1.eq.4 .OR. ngav1.eq.5) then
         if(iter.ge.iterbf .OR. nnstpp.ge.1) then
            tok=tok*helinp/helout
            tok=tok*wght+tokn*(1.d0-wght)
         endif
      endif

      cnor=amu0*tok/tokp
      i_rus=0
      if(i_rus.eq.1)then
         tokp=0.d0
         do i=2,ni1
            do j=2,nj1
               curf(i,j) =  0.d0
               iprij=ipr(i,j)
               xi=r(i)
               psi=un(i,j)
               sij=dri(i)*dzj(j)
               if(iprij.eq.1) then
                  curf(i,j) =  funcur( xi,psi )*sij
                  tokp=tokp+curf(i,j)
               endif
            enddo 
         enddo 
         cnor=amu0*tok/tokp
      endif

      if(key_plc.eq.0) cnor=1.d0

! right hand side definition

      do i=2,ni1
         do j=2,nj1
            curf(i,j) = cnor*curf(i,j)
            il=nlin(i,j)
            if(iter.eq.1) then
               right(il)=-curf(i,j)
            else
               right(il)=-curf(i,j)*omg-(1.d0-omg)*curn(i,j)
            endif
            curn(i,j)=curf(i,j)
            curf(i,j)= curf(i,j)/(dri(i)*dzj(j))
         enddo
      enddo

      return
      end subroutine right0

!----------------------------------------------------------------
      subroutine xpoint(ux, rx, zx, ix, jx, dp, numxp, kodex)

      use comblc, only: ni1, nj1, rmin, rmax, zmin, zmax, r, z, u, dr

      implicit none

      integer, parameter :: nshp=10

      integer, intent(in) :: numxp
      integer, intent(out) :: kodex
      real*8, intent(out), dimension(5) :: dp
      integer, intent(inout) :: ix, jx
      real*8, intent(inout) :: rx, zx

      integer :: i, j, k, l, icx, jcx, ixp, jxp, ix00, jx00, 
     &   numitc, numit, nsh
      real*8 :: ux, sdmin, dlx, det, rr, zz, rrx, zzx, 
     &   dold, dnew, deld01
      real*8, dimension(nshp) :: xs, ys, fun

      ix00=ix
      jx00=jx

      numitc=0
      kodex=0

 444  continue

         numitc=numitc+1
         if(numitc.ge.30) then
            return
         endif

         if(rx.ge.rmax .OR. rx.le.rmin .OR.
     &      zx.ge.zmax .OR. zx.le.zmin) then
            ux=-1.d6
            kodex=1
            return
         endif

! definition of cell, containig x-point

         do i=1,ni1
            icx=i
            if( (rx.lt.r(i+1)) .AND. (rx.ge.r(i)) ) EXIT
         enddo

         do j=1,nj1
            jcx=j
            if( (zx.lt.z(j+1)) .AND. (zx.ge.z(j)) ) EXIT
         enddo

         if(icx.eq.1 .OR. icx.eq.ni1 .OR.
     &      jcx.eq.1 .OR. jcx.eq.nj1) then
            ux=-1.d6
            kodex=1
            return
         endif

         sdmin=rmax
         do k=0,1
            rr=r(icx+k)
            do l=0,1
               zz=z(jcx+l)
               dlx=dsqrt( (rr-rx)**2+(zz-zx)**2 )
               if(dlx.lt.sdmin) then
                  sdmin=dlx
                  ix=icx+k
                  jx=jcx+l
               endif
            enddo
         enddo

         numit=0

 555     continue

         numit=numit+1

         nsh=1
         xs(nsh)=r(ix)
         ys(nsh)=z(jx)
         fun(nsh)=u(ix,jx)

         do k=-1,1
            i= ix+k
            do l=-1,1
               j= jx+l
               if(k.eq.0  .AND. l.eq.0 ) CYCLE
               nsh=nsh+1
               xs(nsh)=r(i)
               ys(nsh)=z(j)
               fun(nsh)=u(i,j)
            enddo
         enddo

         call deriv5(xs,ys,fun,nsh,5,dp)

         DET = dp(3)*dp(5) - dp(4)**2
         Rx = Xs(1) + ( dp(2)*dp(4) - dp(1)*dp(5) )/DET
         Zx = Ys(1) + ( dp(1)*dp(4) - dp(2)*dp(3) )/DET

         rrx=xs(1)
         zzx=ys(1)

         ux = fun(1) + dp(1)*(rx-rrx) + dp(2)*(zx-zzx) +
     &           0.5d0*dp(3)*(rx-rrx)*(rx-rrx)         +
     &                 dp(4)*(rx-rrx)*(zx-zzx)         +
     &           0.5d0*dp(5)*(zx-zzx)*(zx-zzx)

         if( rx.gt.r(ix+1) .OR. rx.lt.r(ix-1) .or.
     &       zx.gt.z(jx+1) .OR. zx.lt.z(jx-1) ) goto 444

      sdmin=rmax
      do k=-1,1
         rr=r(ix+k)
         do l=-1,1
            zz=z(jx+l)
            dlx=dsqrt( (rr-rx)**2+(zz-zx)**2 )
            if(dlx.lt.sdmin) then
               sdmin=dlx
               ixp=ix+k
               jxp=jx+l
            endif
         enddo
      enddo

      if(ixp.ne.ix .OR. jxp.ne.jx) then

         if(numit.ge.20) then
            return
         endif

         ix=ixp
         jx=jxp

         goto 555

      endif

      if(ix.ne.ix00 .OR. jx.ne.jx00) then

         dold=dsqrt( (rx-r(ix00))**2 + (zx-z(jx00))**2 )
         dnew=dsqrt( (rx-r(ix))**2 + (zx-z(jx))**2 )
         deld01=dabs(dold-dnew)

         if(deld01.lt.0.001d0*dr(ix)) then
            ix=ix00
            jx=jx00
            goto 555
         endif

      endif

      return
      end subroutine xpoint

!----------------------------------------------------------------
      real*8 function squtr(x1, x2, x3, y1, y2, y3)

      implicit none

      real*8, intent(in) :: x1, x2, x3, y1, y2, y3

      squtr = 0.5d0*(y3*(x2-x1)+y1*(x3-x2)+y2*(x1-x3))

      return
      end function squtr

!----------------------------------------------------------------
      real*8 function bcur(r1, r2, r3, r4, z1, z2, z3, z4,
     &                     psi1, psi2, psi3, psi4)

      implicit none

      real*8, intent(in) :: r1, r2, r3, r4, z1, z2, z3, z4,
     &                      psi1, psi2, psi3, psi4

      integer :: iflag1, iflag2
      real*8 :: cur1, cur2, cur3, cur4, cur10, cur20, sij, zbcur,
     &   r10, r20, z10, z20, psi10, psi20
      real*8, external :: xzer, squtr, funcur
 
      iflag1=0
      iflag2=0

      if(psi2.gt.0.d0 .AND. psi3.gt.0.d0 .AND. psi4.gt.0.d0) then
         cur1=funcur(r1, psi1)
         cur2=funcur(r2, psi2)
         cur3=funcur(r3, psi3)
         cur4=funcur(r4, psi4)

         sij = squtr(r1, r2, r3, z1, z2, z3) + 
     &         squtr(r1, r3, r4, z1, z3, z4)
         bcur=0.25d0*sij*(cur1+cur2+cur3+cur4)
         return
      endif

! definition of first zero

      if(psi2.le.0.d0) then
         r10=xzer(r1, r2, psi1, psi2)
         z10=xzer(z1, z2, psi1, psi2)
         psi10=0.d0
      elseif(psi3.le.0.d0) then
         r10=xzer(r2, r3, psi2, psi3)
         z10=xzer(z2, z3, psi2, psi3)
         psi10=0.
      elseif(psi4.le.0.d0) then
         iflag1=1
         r10=xzer(r3, r4, psi3, psi4)
         z10=xzer(z3, z4, psi3, psi4)
         psi10=0.d0
      endif

! definition of second zero

      if(psi4.le.0.d0) then
         r20=xzer(r1, r4, psi1, psi4)
         z20=xzer(z1, z4, psi1, psi4)
         psi20=0.d0
      elseif(psi3.le.0.d0) then
         r20=xzer(r4, r3, psi4, psi3)
         z20=xzer(z4, z3, psi4, psi3)
         psi20=0.d0
      elseif(psi2.le.0.d0) then
         iflag2=1
         r20=xzer(r3, r2, psi3, psi2)
         z20=xzer(z3, z2, psi3, psi2)
         psi20=0.d0
      endif

      cur1=funcur(r1, psi1)
      cur2=funcur(r2, psi2)
      cur3=funcur(r3, psi3)
      cur4=funcur(r4, psi4)
      cur10=funcur(r10, psi10)
      cur20=funcur(r20, psi20)

      zbcur = squtr(r1, r20, r4 , z1, z20, z4 )*(cur1+cur20+cur4)  +
     &        squtr(r1, r10, r20, z1, z10, z20)*(cur1+cur10+cur20) +
     &        squtr(r1, r2 , r10, z1, z2 , z10)*(cur1+cur2+cur10)  +
     & iflag1*squtr(r2, r3 , r10, z2, z3 , z10)*(cur2+cur3+cur10)  +
     & iflag2*squtr(r3, r4 , r20, z3, z4 , z20)*(cur3+cur4+cur20)

      bcur=zbcur/3.d0

      return
      end function bcur

!----------------------------------------------------------------
      subroutine xdetal(i, j, ix, jx, rxpn, zxpn, dp, nloc,
     &                  cur1, cur2, cur3, cur4)

      use comblc, only: r, z, rm, zm, u, um, up

      implicit none

      integer, parameter :: nlocp=64

      integer, intent(in) :: i, j, ix, jx, nloc
      real*8, intent(in) :: rxpn, zxpn
      real*8, intent(in), dimension(5) :: dp
      real*8, intent(out) :: cur1, cur2, cur3, cur4

      integer :: is, js, nlc
      real*8 :: drloc, dzloc, ris, rix, zjs, zjx, uisjs, uix, 
     &   sloc, psilim, psiloc, rlim, frlim, rmlim
      real*8, dimension(nlocp) :: rloc, zloc
      real*8, dimension(nlocp, nlocp) :: uloc, curloc
      real*8, external :: funcur

      psilim=0.05d0
      nlc=(nloc+1)/2

      drloc=(r(i+1)-r(i))/(nloc-1.d0)
      dzloc=(z(j+1)-z(j))/(nloc-1.d0)

      rloc(1)=r(i)
      zloc(1)=z(j)

      do is=2,nloc
         rloc(is)=rloc(is-1)+drloc
         zloc(is)=zloc(is-1)+dzloc
      enddo

      uix=u(ix,jx)
      rix=r(ix)
      zjx=z(jx)

      do is=1,nloc
         ris=rloc(is)
         do js=1,nloc
            zjs=zloc(js)
            uisjs = uix + dp(1)*(ris-rix) + dp(2)*(zjs-zjx) +
     &                0.5*dp(3)*(ris-rix)*(ris-rix) +
     &                    dp(4)*(ris-rix)*(zjs-zjx) +
     &                0.5*dp(5)*(zjs-zjx)*(zjs-zjx)

            uloc(is,js)=(uisjs-up)/(um-up)
         enddo
      enddo

      do is=1,nloc
         ris=rloc(is)
         do js=1,nloc
            zjs=zloc(js)
            psiloc=uloc(is,js)
            if(psiloc.le.0.d0) then
               curloc(is,js)=0.d0
               CYCLE
            endif

            rlim=frlim(dp,zjs,rxpn,zxpn,rm,zm)
            rmlim=frlim(dp,zm,rxpn,zxpn,rm,zm)
            if( (ris -rlim)*(rm-rmlim).le.0.) then
               curloc(is,js)=0.d0
               CYCLE
            endif

            curloc(is,js)=funcur(ris,psiloc)

         enddo
      enddo

      sloc=drloc*dzloc

! quadr (i,j)

      do is=1,nlc-1
         do js=1,nlc-1
            cur1 = cur1 + curloc(is  , js)   + curloc(is+1, js)  +
     &                    curloc(is+1, js+1) + curloc(is  , js+1)
         enddo
      enddo
      cur1=cur1*0.25d0*sloc

! quadr (i+1,j)

      do is=nlc,nloc-1
         do js=1,nlc-1
            cur2 = cur2 + curloc(is  , js  ) + curloc(is+1, js)  +
     &                    curloc(is+1, js+1) + curloc(is  , js+1)
         enddo
      enddo
      cur2=cur2*0.25d0*sloc

! quadr (i+1,j+1)

      do is=nlc,nloc-1
         do js=nlc,nloc-1
            cur3 = cur3 + curloc(is  , js ) + curloc(is+1, js)  +
     &                    curloc(is+1,js+1) + curloc(is  , js+1)
         enddo
      enddo
      cur3=cur3*0.25d0*sloc

! quadr (i,j+1)

      do is=1,nlc-1
         do js=nlc,nloc-1
            cur4=cur4+curloc(is,js)+curloc(is+1,js)+
     &                 curloc(is+1,js+1)+curloc(is,js+1)
         enddo
      enddo
      cur4=cur4*0.25d0*sloc

      return
      end subroutine xdetal

!----------------------------------------------------------------
      real*8 function frlim(dp, ylim, rx, zx, rm, zm)

      implicit none

      real*8, intent(in) :: ylim, rx, zx, rm, zm
      real*8, intent(in), dimension(5) :: dp

      real*8 :: Dxx, Dxy, Dyy, disc, cc, cdpls, cdmns, ang1, ang2,
     &   c1, c2, dl2x, dl2y

      Dxx=dp(3)
      Dxy=dp(4)
      Dyy=dp(5)
      disc=(Dxy/Dyy)**2 - Dxx/Dyy

      if(Disc.lt.0.) then
         cc=(zm-zx)/(rm-rx)
         cc=-1.d0/cc
      else

         cdpls = -Dxy/Dyy + dsqrt(disc)
         cdmns = -Dxy/Dyy - dsqrt(disc)
         ang1 = 0.5d0*(datan(cdpls)+datan(cdmns))
         ang2 =-0.5d0*(datan(1.d0/cdpls)+datan(1.d0/cdmns))

! calculation D2u/Dl2(direction ang1)
         c1=dtan(ang1)
         Dl2x=Dyy*c1*c1+2.d0*Dxy*c1+Dxx

! calculation D2u/Dl2(direction ang2)
         c2=dtan(ang2)
         Dl2y=Dyy*c2*c2+2.*Dxy*c2+Dxx

         if(Dl2x.lt.0.) then
            cc=c1
         elseif(Dl2y.lt.0.) then
            cc=c2
         else
            cc=(zm-zx)/(rm-rx)
            cc=-1.d0/cc
         endif

      endif

      frlim=rx+(ylim-zx)/cc

      return
      end function frlim
