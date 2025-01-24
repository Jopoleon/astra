      subroutine rightg

! right-hand side for problem L(g)=J
      use sp_parameters
      use keys
      use compol, only: tok, tokp, tokpp, tokff, tokww, right, ngav,
     & iplas, neqpla, erru, fvac, f, psia, psim, 
     &     nr, nt, nt1, cur, psin, r, rm, dpdpsi, dfdpsi, dwdpsi, cnor,
     & sq1, sq2, sq3, sq4

      implicit none

      integer :: i, j, il
      real*8 :: cur_mu, psn, r0, curp, curf, curw, curcen, sqcen, sqk, 
     &   tok_pl, dfdpsi_sur
      integer, external :: numlin
      real*8, external :: tabf, tabp, tabw

      save cur_mu

      tokp=0.d0
      tokff=0.d0
      tokpp=0.d0
      tokww=0.d0

      do il=1,neqp
         right(il)=0.d0
      enddo

!use the following if p', ff'(normalized poloidal flux) are assumed in ASTRA (kastr=1)
      if(kstep.eq.0 .OR. ngav.eq.0) then       
!!use the following if p', ff'(normalized toroidal flux) are assumed in ASTRA (kastr=1)
         do i=1,iplas
            do j=1,Nt
               cur(i,j)=0.d0
            enddo
         enddo

!! central point

         psn=psin(1,2)
         r0=r(1,2)

         curp=tabp(psn)
         curf=tabf(psn)
         curw=tabw(psn)

         dpdpsi(1)=curp
         dfdpsi(1)=curf
         dwdpsi(1)=curw

         curcen=r0*curp+curf/r0+curw*r0**3
         sqcen=0.d0

         do j=2,nt1
            sqcen=sqcen+sq1(1,j)+sq4(1,j)
         enddo

         tokp=curcen*sqcen
         tokff=curf*sqcen/r0
         tokpp=curp*sqcen*r0
         tokww=curw*sqcen*r0**3

         il=numlin(1,1,nr,nt)

         right(il)=curcen*sqcen

! regular points

         do i=2,Iplas
            do j=2,Nt1
               if(i.ne.iplas) then
                  sqk=sq1(i,j)+sq2(i-1,j)+sq3(i-1,j-1)+sq4(i,j-1)
                  r0=r(i,j)
                  psn=psin(i,j)
               else
                  sqk=sq2(i-1,j)+sq3(i-1,j-1)
                  r0=r(i,j)
                  psn=psin(i,j)
                  r0=(3.d0*r(i,j)+r(i-1,j))*0.25d0
                  psn=(3.d0*psin(i,j)+psin(i-1,j))*0.25d0
               endif

               curp=tabp(psn)
               curf=tabf(psn)
               curw=tabw(psn)
               cur(i,j)=r0*curp+curf/r0+curw*r0**3

               tokp=tokp+cur(i,j)*sqk
               tokff=tokff+curf*sqk/r0
               tokpp=tokpp+curp*sqk*r0
               tokww=tokww+curw*sqk*r0**3

               il=numlin(i,j,nr,nt)
               right(il)=cur(i,j)*sqk
            enddo

            dpdpsi(i)=curp
            dfdpsi(i)=curf
            dwdpsi(i)=curw

         enddo

         dpdpsi(iplas)=tabp(0.d0)
         dfdpsi(iplas)=tabf(0.d0)
         dwdpsi(iplas)=tabw(0.d0)

         cnor=amu0*tok/tokp
         if(key_plc.eq.0) cnor=1.d0
         if(key_0st.eq.1) cnor=1.d0
         tok_pl=0.d0

         do j=1,Nt
            cur(1,j)=curcen*cnor
         enddo

         do i=2,iplas
            do j=2,Nt1
               cur(i,j)=cur(i,j)*cnor
            enddo
         enddo

         do i=1,iplas
            cur(i,1)=cur(i,nt1)
            cur(i,nt)=cur(i,2)
         enddo

         do il=1,neqpla
            right(il)=right(il)*cnor
         enddo

         do i=1,iplas
            dpdpsi(i)=dpdpsi(i)*cnor
            dfdpsi(i)=dfdpsi(i)*cnor
            dwdpsi(i)=dwdpsi(i)*cnor
         enddo

         cur_mu=tokp*cnor
         tokp=tokp/amu0*cnor

      else
         if(kastr.eq.1 .AND. key_prs.eq.1 .AND. erru.lt.5.d-3) then
            call pres_d_psi !
         endif

         call procof(1,cur_mu)

         dfdpsi_sur=-(fvac**2-f(iplas-1)**2)/psia(iplas-1)/psim

! central point

         curcen=rm*dpdpsi(1)+dfdpsi(1)/rm+dwdpsi(1)*rm**3
         sqcen=0.d0

         do j=2,nt1
            sqcen=sqcen+sq1(1,j)+sq4(1,j)
         enddo

         tokp=curcen*sqcen
         il=numlin(1,1,nr,nt)
         right(il)=curcen*sqcen
         tok_pl=0.d0

         do j=1,Nt
            cur(1,j)=curcen
         enddo

! regular points

         do i=2,Iplas
            do j=2,Nt1

               if(i.ne.iplas) then
                  sqk=sq1(i,j)+sq2(i-1,j)+sq3(i-1,j-1)+sq4(i,j-1)
                  r0=r(i,j)
                  curp=dpdpsi(i)
                  curf=dfdpsi(i)
                  curw=dwdpsi(i)
               else
                  sqk=sq2(i-1,j)+sq3(i-1,j-1)
                  r0=(3.d0*r(i,j)+r(i-1,j))*0.25d0
                  curp=dpdpsi(i)
                  curf=dfdpsi(i)
                  curw=dwdpsi(i)
               endif

               cur(i,j)=r0*curp+curf/r0+curw*r0**3
               tokp=tokp+cur(i,j)*sqk
               il=numlin(i,j,nr,nt)
               right(il)=cur(i,j)*sqk
            enddo
         enddo

         do i=1,iplas
            cur(i,1)=cur(i,nt1)
            cur(i,nt)=cur(i,2)
         enddo
         cnor=cur_mu/tokp
         if(key_plc.eq.0) cnor=1.d0
         tokp=tokp/amu0

         if(ngav.eq.0 ) then
            tokp=cur_mu/amu0
            do il=1,neqpla
               right(il)=right(il)*cnor
            enddo
            do i=1,iplas
               do j=1,nt
                  cur(i,j)=cur(i,j)*cnor
               enddo
            enddo
         endif

         if(kpr.eq.1) then
            write(*,*) 'cnor=',cnor
         endif
      endif

      return
      end subroutine rightg

!----------------------------------------------------------------
      subroutine psib_pla(pspl_av)

      use sp_parameters
      use compol, only: iplas, nt, nt1, psi, psip, dlt, cur, sq2, sq3

      implicit none

      real*8, intent(out) :: pspl_av

      integer :: i, j, l, jb
      real*8 :: a1, a2, a3, g1, g2, g3, sqk, dltk, dgdnl, psb
      real*8 :: avr_bnd
      real*8, dimension(ntp) :: pspl_b, dgdn
      real*8, dimension(nrp, ntp) :: a12, a23, a34, a14, a13, a24
      real*8, dimension(ntp, ntp) :: binadg

      common/comaaa/ a12, a23, a34, a14, a13, a24
      common/com_bgr/ binadg, dgdn

      i=iplas

      do j=2,nt1
         a1=a13(i-1,j-1)
         a2=a34(i-1,j-1)+a12(i-1,j)
         a3=a24(i-1,j)
         g1=psi(i-1,j-1)-psip
         g2=psi(i-1,j)-psip
         g3=psi(i-1,j+1)-psip
         sqk=sq2(i-1,j)+sq3(i-1,j-1)
         dltk=(dlt(i,j-1)+dlt(i,j))*0.5d0
         dgdnl=-(a1*g1+a2*g2+a3*g3) + cur(i,j)*sqk !
         dgdn(j)=dgdnl/dltk
      enddo

      dgdn(1)=dgdn(nt1)
      dgdn(nt)=dgdn(2)

      do l=2,nt1
         psb=0.d0
         do jb=2,nt1
            psb=psb+binadg(jb,l)*(dgdn(jb)+dgdn(jb+1))*0.5d0
         enddo
         pspl_b(l)=psb
      enddo

      pspl_b(1)=pspl_b(nt1)
      pspl_b(nt)=pspl_b(2)

      pspl_av = avr_bnd(pspl_b)

      return
      end subroutine psib_pla

!----------------------------------------------------------------
      subroutine psloop_e(rlop, zlop, psilop, psilopk, nlop,
     &   ntp, rcondzzzz)

      use sp_parameters, only: nip, njp

      implicit none

      integer, intent(in) :: nlop, ntp
      real*8, intent(in), dimension(*) :: rlop, zlop
      real*8, intent(out), dimension(*) :: psilop
      real*8, intent(out), dimension(ntp,500) :: psilopk
      double precision, intent(out) :: rcondzzzz

      integer :: i, j, ic, jc, imax, jmax, ix1, jx1, ix2, jx2
      real*8 :: ux0, ux1, ux2, up, um, xm, ym, xx0, yx0, 
     &   xx1, yx1, xx2, yx2, xx10, yx10, xx20, yx20, xm0, ym0, 
     &   psi_bon, rcondzzz, ddx, ddy, u1, u2, u3, u4,
     &   r0, r1, r2, z0, z1, z2
      real*8 :: blin_
      real*8, dimension(nip) :: x, dx, dxi, x12
      real*8, dimension(njp) :: y, dy, dyj
      real*8, dimension(nip, njp) :: u, ue, un, ui, g
      real*8, dimension(nip, njp, 500) :: ue_k
      common /comrz/ x, y, dx, dy, dxi, dyj, x12
      common /compot/ u, ue, un, ui, g, ux0, ux1, ux2, up, um, xm, ym,
     &                xx0, yx0, xx1, yx1, xx2, yx2, imax, jmax, 
     &                ix1, jx1, ix2, jx2, 
     &                xx10, yx10, xx20, yx20, xm0, ym0, 
     &                psi_bon, ue_k, rcondzzz
	
      ddx=x(2)-x(1)
      ddy=y(2)-y(1)

      do i=1, nlop

         r0=rlop(i)
         z0=zlop(i)
         ic=(r0-x(1))/ddx+1
         jc=(z0-y(1))/ddy+1

         r1=x(ic)
         r2=x(ic+1)

         z1=y(jc)
         z2=y(jc+1)

         u1=ue(ic,jc)
         u2=ue(ic+1,jc)
         u3=ue(ic+1,jc+1)
         u4=ue(ic,jc+1)

         psilop(i)=blin_(r0,z0,r1,r2,z1,z2,u1,u2,u3,u4 )

         do j=1,nint(rcondzzz)	
            u1=ue_k(ic,jc,j)
            u2=ue_k(ic+1,jc,j)
            u3=ue_k(ic+1,jc+1,j)
            u4=ue_k(ic,jc+1,j)
            psilopk(i,j)=blin_(r0,z0,r1,r2,z1,z2,u1,u2,u3,u4)
         enddo
	 rcondzzzz=rcondzzz
      enddo

      return
      end subroutine psloop_e

!----------------------------------------------------------------
      subroutine psib_ext(psex_av)

      use sp_parameters
      use compol, only: iplas, nt, nt1, r, z

      implicit none

      real*8, intent(out) :: psex_av

      integer :: i, j
      real*8 :: avr_bnd
      real*8, dimension(ntp) :: psex_b, rbon, zbon
      real*8 :: psexk_b(ntp, 500), psexk_av(500)
      double precision rzzz

      i=iplas

      do j=1,nt
         rbon(j)=r(i,j)
         zbon(j)=z(i,j)
      enddo

      call psloop_e(rbon, zbon, psex_b, psexk_b, nt, ntp, rzzz)

      psex_b(1)=psex_b(nt1)
      psex_b(nt)=psex_b(2)
      psex_av = avr_bnd(psex_b)

      do i=1,nint(rzzz)
         psexk_b(1,i)=psexk_b(nt1,i)
         psexk_b(nt,i)=psexk_b(2,i)
         psexk_av(i) = avr_bnd(psexk_b(1:nt,i))
      enddo

      return
      end subroutine psib_ext
