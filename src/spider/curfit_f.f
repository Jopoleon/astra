      subroutine closest_knot(r0, z0, ik, jk)

      use comblc, only: rmin, rmax, zmin, zmax, dr, dz, r, z

      implicit none

      real*8, intent(in) :: r0, z0
      integer, intent(out) :: ik, jk
      integer :: k, l, ic, jc
      real*8 :: rr, zz, dlx, sdmin

      ic=(r0-rmin)/dr(1)+1
      jc=(z0-zmin)/dz(1)+1

! Looking for nearest knot

      sdmin=rmax

      do k=0,1
         rr=r(ic+k)
         do l=0,1
            zz=z(jc+l)
            dlx=dsqrt( (rr-r0)**2+(zz-z0)**2 )
            if(dlx.lt.sdmin) then
               sdmin=dlx
               ik=ic+k
               jk=jc+l
            endif
         enddo
      enddo

      return
      end subroutine closest_knot

!----------------------------------------------------------------
      subroutine prefit(rk, zk, ncpfc, NECON, WECON, rax, zax, 
     &   alp_b, psi_bnd)
        
      use bnd_modul
      use sp_parameters
      use iopath, only: path
      use parcur
      use comevl, only: nequi
      use keys, only: kastr, ksnf, kxwx

      implicit none

      integer, intent(in) :: ncpfc, necon(*)
      real*8, intent(in) :: rax, zax, alp_b
      real*8, intent(in), dimension(*) :: rk, zk, wecon
      real*8, intent(out) :: psi_bnd

      integer :: i, j, ik, k, l, i_min, i_max, j_min, j_max
      real*8 :: wwl, wdk, wsig, cxcs, rb_min, rb_max, zb_min, zb_max
      character(len=80) :: fname

      alph=alp_b
      r_ax=rax
      z_ax=zax

      if(kastr.eq.0) then 
         write(fname,'(a,a)') TRIM(path), '/bonfit.dat'
         open(1,file=fname,form='formatted')
            read(1,*) wwl,wdk,wsig
            read(1,*) (d_wght(ik),ik=1,NEQUI)
            read(1,*) Lfit          !fitting points for plasma boundary
            if(Lfit.ne.0) then
               do l=1,Lfit
                  read(1,*) rfit(l),zfit(l)
               enddo
            endif
            read(1,*) Lpre   !exact prescribed points for plasma boundary
            if(Lpre.ne.0) then
               do l=1,Lpre
                  read(1,*) rpre(l),zpre(l)
               enddo
            endif
            read(1,*) rx_p,zx_p
            read(1,*) psi_bnd
            if(kxwx.eq.2) then 
               read(1,*) rx2_p,zx2_p
            endif
         close(1)

         if(ksnf.eq.1) then
            write(fname,'(a,a)') TRIM(path), '/matC.wr'
            open(1,file=fname,form='formatted')
               read(1,*) c_wght
               read(1,*) n_cc,m_cc
               read(1,*) (g_wght(i),i=1,n_cc)
               read(1,*) ((Cpir(i,j),i=1,n_cc),j=1,m_cc)
            close(1)
        
            do i=1,n_cc
               do j=1,m_cc
                  Cpir(i,j)=Cpir(i,j)*g_wght(i)
                enddo
            enddo       
        
            do i=1,m_cc
               do j=1,m_cc
                  CxCs=0.d0
                  do k=1,n_cc
                     CxCs=CxCs+Cpir(k,i)*Cpir(k,j)
                  enddo
                  CxC(i,j)=CxCs
               enddo
            enddo
         endif          

      elseif(kastr.eq.1) then 

         wwl= 1.d0
         wdk= 1.d0
         wsig=1.d-2
         write(*,*) nequi

         write(fname,'(a,a)') TRIM(path), '/bonfit.dat'
         open(1,file=fname,form='formatted')
            read(1,*) wwl,wdk,wsig
            read(1,*) (d_wght(ik),ik=1,NEQUI)
         close(1)

         Lfit=nbtab

         rb_max=rbtab(1)
         rb_min=rbtab(1)
         zb_max=zbtab(1)
         zb_min=zbtab(1)

         do l=1,Lfit
            rfit(l)=rbtab(l)
            zfit(l)=zbtab(l)   
            if( rfit(l).ge.rb_max ) then
               i_max=l
               rb_max=rfit(l)
            endif 
            if( rfit(l).le.rb_min ) then 
               i_min=l          
               rb_min=rfit(l)
            endif 
            if( zfit(l).ge.zb_max ) then
               j_max=l          
               zb_max=zfit(l)
            endif
            if( zfit(l).le.zb_min ) then
               j_min=l          
               zb_min=zfit(l)
            endif
         enddo

         Lpre=2

         rpre(1)=rfit(i_min)
         zpre(1)=zfit(i_min)
         rpre(2)=rfit(i_max)
         zpre(2)=zfit(i_max)
           
         rx_p=rax
         zx_p=zax

      endif

      do l=1,Lfit
          w_wght(l)=wwl*2.d0
      enddo
      s_wght=wsig

      call grimat(rk,zk,ncpfc,NECON, WECON)
           
      if(kxwx.eq.2)then
         call grimat2x(rk,zk,ncpfc,NECON, WECON)
      endif

      return
      end subroutine prefit

!----------------------------------------------------------------
      subroutine curfit(rk, zk, nk, NECON, WECON, psi_bnd)

      use sp_parameters
      use comevl, only: nequi, pfceqw
      use keys, only: kxwx, ksnf

      implicit none

      integer, intent(in) :: nk, NECON(*)
      real*8, intent(in), dimension(*) :: rk, zk, wecon
      real*8, intent(out) :: psi_bnd

      integer :: ik

      if(ksnf.eq.1) then
         call psi_cpn
         call cursol_snf(NEQUI,PFCEQW,psi_bnd)
      else        
         if(kxwx.eq.1)then
            call psi_cpn
            call cursol(NEQUI,PFCEQW,psi_bnd)
         elseif(kxwx.eq.2)then
            call psi_cpn
            call psi_cpn_2x
            call cursol_2x(NEQUI,PFCEQW,psi_bnd)
         elseif(kxwx.eq.0)then
            call psi_cpn_wx
            call cursol_wx(NEQUI,PFCEQW,psi_bnd)
         endif
      endif
      write(*,*) 'coil currents' ,nequi

      do ik=1,nequi
         write(*,*) PFCEQW(ik),ik 
      enddo

      return
      end subroutine curfit

!----------------------------------------------------------------
      subroutine curfit_L(rk, zk, nk, NECON, WECON, psi_bnd)

      use sp_parameters
      use comevl, only: nequi, pfceqw

      implicit none

      integer, intent(in) :: nk, NECON(*)
      real*8, intent(in), dimension(*) :: rk, zk, wecon
      real*8, intent(out) :: psi_bnd

      call psi_cpn
      call cursol_(NEQUI, PFCEQW, psi_bnd)

      return
      end subroutine curfit_L

!----------------------------------------------------------------
      subroutine curfit_(rk,zk,nk,NECON,WECON,psi_bnd)

      use sp_parameters
      use comevl, only: nequi, pfceqw

      implicit none

      integer, intent(in) :: nk, NECON(*)
      real*8, intent(in), dimension(*) :: rk, zk, wecon
      real*8, intent(out) :: psi_bnd

      integer :: ik

      call bonpsi
      call cursol_(NEQUI,PFCEQW,psi_bnd)

      write(*,*) 'coil currents',nequi,nequi 

      do ik=1,nequi
         write(*,*) PFCEQW(ik),ik 
      enddo

      return
      end subroutine curfit_

!----------------------------------------------------------------
      subroutine grimat(rk,zk,ncpfc,NECON, WECON )

      use sp_parameters
      use parcur
      use comevl, only: nequi

      implicit none

      integer, intent(in) :: ncpfc, NECON(*)
      real*8, intent(in), dimension(*) :: rk, zk, wecon

      integer :: ip, l, ik, iq
      real*8 :: dgdr, dgdz, d2gdrr, d2gdrz, d2gdzz, zgindk, zgindp
      real*8 :: greeni

      do iq=1,NEQUI

         do l=1,Lfit
            Gindk(iq,l)=0.d0
         enddo
         do l=1,Lpre
            Gindp(iq,l)=0.d0
         enddo

         Gx_r(iq)=0.d0
         Gx_z(iq)=0.d0

         Gx_rz(iq)=0.d0
         Gx_zz(iq)=0.d0
         Gx_rr(iq)=0.d0

         if(Lfit.ne.0) then
            do l=1,Lfit
               do ik=1,ncpfc
                  if( necon(ik) .eq. iq ) then
                     zgindk=greeni(rfit(l),zfit(l),rk(ik),zk(ik))/pi
                     Gindk(iq,l)=Gindk(iq,l)+zgindk*wecon(ik)
                  endif
               enddo
            enddo
         endif

         if(Lpre.ne.0) then
            do l=1,Lpre
               do ik=1,ncpfc
                  if( necon(ik) .eq. iq ) then
                     zgindp=greeni(rpre(l),zpre(l),rk(ik),zk(ik))/pi
                     Gindp(iq,l)=Gindp(iq,l)+zgindp*wecon(ik)
                  endif
               enddo
            enddo
         endif

         do ik=1,ncpfc
            if( necon(ik) .eq. iq ) then
               call grGREN(rk(ik),zk(ik),Rx_p,Zx_p, dGdr,dGdz)
               Gx_r(iq)=Gx_r(iq)+dGdr*wecon(ik)/pi
               Gx_z(iq)=Gx_z(iq)+dGdz*wecon(ik)/pi
               
!!for snowflake
               call d2GREN(rk(ik),zk(ik),Rx_p,Zx_p,d2Gdrz,d2Gdzz,d2Gdrr)
               Gx_rz(iq)=Gx_rz(iq)+d2Gdrz*wecon(ik)/pi
               Gx_zz(iq)=Gx_zz(iq)+d2Gdzz*wecon(ik)/pi
               Gx_rr(iq)=Gx_rr(iq)+d2Gdrr*wecon(ik)/pi
!
            endif
         enddo

      enddo

      return
      end subroutine grimat

!----------------------------------------------------------------
      subroutine grimat2x(rk,zk,ncpfc,NECON, WECON )

      use sp_parameters
      use parcur
      use comevl, only: nequi

      implicit none

      integer, intent(in) :: ncpfc, NECON(*)
      real*8, intent(in), dimension(*) :: rk, zk, wecon

      integer :: iq, ik
      real*8 :: dgdr, dgdz

      do iq=1,NEQUI
         Gx2_r(iq)=0.d0
         Gx2_z(iq)=0.d0
         do ik=1,ncpfc
            if( necon(ik) .eq. iq )  then
               call grGREN(rk(ik),zk(ik),Rx2_p,Zx2_p, dGdr,dGdz)
               Gx2_r(iq)=Gx2_r(iq)+dGdr*wecon(ik)/pi
               Gx2_z(iq)=Gx2_z(iq)+dGdz*wecon(ik)/pi
            endif
         enddo
      enddo

      return
      end subroutine grimat2x

!----------------------------------------------------------------
      subroutine precal(rk, zk, nk, NECON, WECON)

      use sp_parameters, only: pi
      use parcur
      use comevl, only: NEQUI
      use comblc, only: r, z, ui

      implicit none

      integer, parameter :: nshp=10

      real*8, intent(in), dimension(*) :: rk, zk, WECON
      integer, intent(in) :: nk, NECON(*)

      integer :: i, j, k, l, iq, ik, jk, nsh, ncpfc
      real*8, dimension(nshp) :: xs, ys, fun
      real*8 :: dp(5), zginda, zgindx, greeni, rra, zza

      ncpfc=nk

      do iq=1,NEQUI
         Ginda(iq)=0.d0
         Gindx(iq)=0.d0
         do ik=1,ncpfc
            if( necon(ik) .eq. iq )  then
               zginda=greeni(r_ax,z_ax,rk(ik),zk(ik))/pi
               Ginda(iq)=Ginda(iq)+zginda*wecon(ik)
               zgindx=greeni(rx_p,zx_p,rk(ik),zk(ik))/pi
               Gindx(iq)=Gindx(iq)+zgindx*wecon(ik)
            endif
         enddo
      enddo

      call closest_knot(r_ax, z_ax, ik, jk)

      nsh=1
      xs(nsh)=r(ik)
      ys(nsh)=z(jk)
      fun(nsh)=ui(ik,jk)

      do k=-1,1
         i= ik+k
         do l=-1,1

            j= jk+l
            if(k.eq.0  .AND. l.eq.0 ) CYCLE
            nsh=nsh+1
            xs(nsh)=r(i)
            ys(nsh)=z(j)
            fun(nsh)=ui(i,j)
         enddo
      enddo

      rra=xs(1)
      zza=ys(1)

      call deriv5(xs, ys, fun, nsh, 5, dp)

      psip_a=fun(1)+ dp(1)*(r_ax-rra) + dp(2)*(z_ax-zza)
     &          + 0.5d0*dp(3)*(r_ax-rra)*(r_ax-rra)
     &          +       dp(4)*(r_ax-rra)*(z_ax-zza)
     &          + 0.5d0*dp(5)*(z_ax-zza)*(z_ax-zza)

      return
      end subroutine precal

!----------------------------------------------------------------
      subroutine precal_wx(rk, zk, nk, NECON, WECON )

      use sp_parameters, only: pi
      use parcur
      use comevl, only: NEQUI
      use comblc, only: rx0, zx0, r, z, ui

      implicit none

      integer, parameter :: nshp=10

      real*8, intent(in), dimension(*) :: rk, zk, WECON
      integer, intent(in) :: nk, NECON(*)

      integer :: i, j, k, l, iq, ik, jk, nsh, ncpfc
      real*8, dimension(nshp) :: xs, ys, fun
      real*8 :: dp(5), zginda, zgindx, greeni, rra, zza

      ncpfc=nk

      do iq=1,NEQUI
         Ginda(iq)=0.d0
         Gindx(iq)=0.d0
         do ik=1,ncpfc
            if( necon(ik) .eq. iq ) then
               zginda=greeni(r_ax,z_ax,rk(ik),zk(ik))/pi
               Ginda(iq)=Ginda(iq)+zginda*wecon(ik)
               zgindx=greeni(rx0,zx0,rk(ik),zk(ik))/pi
               Gindx(iq)=Gindx(iq)+zgindx*wecon(ik)
            endif
         enddo
      enddo

      call closest_knot(r_ax, z_ax, ik, jk)

      nsh=1
      xs(nsh)=r(ik)
      ys(nsh)=z(jk)
      fun(nsh)=ui(ik,jk)

      do k=-1,1
         i= ik+k
         do l=-1,1
            j= jk+l
            if(k.eq.0  .AND. l.eq.0 ) CYCLE
            nsh=nsh+1
            xs(nsh)=r(i)
            ys(nsh)=z(j)
           fun(nsh)=ui(i,j)
         enddo
      enddo

      rra=xs(1)
      zza=ys(1)

      call deriv5(xs,ys,fun,nsh,5,dp)

      psip_a=fun(1)+ dp(1)*(r_ax-rra) + dp(2)*(z_ax-zza)
     &          + 0.5d0*dp(3)*(r_ax-rra)*(r_ax-rra)
     &          +       dp(4)*(r_ax-rra)*(z_ax-zza)
     &          + 0.5d0*dp(5)*(z_ax-zza)*(z_ax-zza)

      return
      end subroutine precal_wx

!----------------------------------------------------------------
      subroutine cursol(n_equi, pf_ceqw, psi_bnd)

      use sp_parameters, only: ncf_p, amu0
      use parcur

      implicit none

      integer, intent(in) :: n_equi
      real*8, intent(out) :: psi_bnd
      real*8, intent(out), dimension(*) :: pf_ceqw

      integer :: j, k, l, ll, iq
      integer, dimension(ncf_p) :: IP
      real*8 :: asum, gr_xp, gz_xp, ps_ma, ps_xp
      real*8, dimension(ncf_p) :: X, Y, psictr
      real*8, dimension(ncf_p, ncf_p) :: A

      do k=1,n_equi 
         do j=1,n_equi 
            asum=0.d0
            do l=1,Lfit 
               asum=asum+Gindk(k,l)*Gindk(j,l)*w_wght(l)
            enddo
            a(j,k)=asum
         enddo
      enddo

      do j=1,n_equi 
         a(j,j)=a(j,j)+d_wght(j)*s_wght
      enddo

! a(j,n_equi+l)*Lam(l) ,j-number of equation

      do j=1,n_equi 
         do l=1,Lpre 
            a(j,n_equi+l)=Gindp(j,l)
         enddo
      enddo

! a(j,n_equi+Lpre+1)*Lam_r(l) ,j-number of equation

      do j=1,n_equi 
         a(j,n_equi+Lpre+1)=Gx_r(j)
      enddo

! a(j,n_equi+Lpre+2)*Lam_z(l) ,j-number of equation

      do j=1,n_equi 
         a(j,n_equi+Lpre+2)=Gx_z(j)
      enddo

! a(j,n_equi+Lpre+3)*Lamx ,j-number of equation

      do j=1,n_equi 
         a(j,n_equi+Lpre+3)=(1.d0-alph)*Ginda(j)+alph*Gindx(j)
      enddo

! a(j,n_equi+Lpre+4)*Psi_b  ,j-number of equation

      do j=1,n_equi 
         asum=0.d0
         do l=1,Lfit 
            asum=asum+w_wght(l)*Gindk(j,l)
         enddo
         a(j,n_equi+Lpre+4)=-asum
      enddo

! right hand side

      do j=1,n_equi 
         asum=0.d0
         do l=1,Lfit 
            asum=asum+w_wght(l)*Gindk(j,l)*psifit(l)
         enddo
         y(j)=-asum+d_wght(j)*s_wght*curref(j)*amu0
      enddo

!------------------------
! equation (.)*d(Lam_l)=0

! a(j,k)*I(k) ,j-number of equation

      do l=1,Lpre 
         j=n_equi+l
         do k=1,n_equi 
            a(j,k)=Gindp(k,l)
         enddo
      enddo

! a(j,n_equi+l)*Lam(l) ,j-number of equation
           
      do l=1,Lpre 
         j=n_equi+l
         do ll=1,Lpre 
            k=n_equi+ll
            a(j,k)=0.d0
         enddo
      enddo

      do l=1,Lpre 
         j=n_equi+l
         k=n_equi+Lpre+1
         a(j,k)=0.d0
         k=n_equi+Lpre+2
         a(j,k)=0.d0
         k=n_equi+Lpre+3
         a(j,k)=0.d0
      enddo

! a(j,n_equi+Lpre+4)*Psi_b  ,j-number of equation
 
      do l=1,Lpre 
         j=n_equi+l
         k=n_equi+Lpre+4 
         a(j,k)=-1.d0
      enddo
! right hand side
      do l=1,Lpre 
         j=n_equi+l
         y(j)=-psipre(l)
      enddo

!------------------------
! equation (.)*d(Lam_r)=0

! a(j,k)*I(k) ,j-number of equation

      j=n_equi+Lpre+1

      do k=1,n_equi 
         a(j,k)=Gx_r(k)
      enddo

! a(j,n_equi+l)*Lam(l) ,j-number of equation
           
      do ll=1,Lpre 
         k=n_equi+ll
         a(j,k)=0.d0
      enddo

      k=n_equi+Lpre+1
      a(j,k)=0.d0
      k=n_equi+Lpre+2
      a(j,k)=0.d0
      k=n_equi+Lpre+3
      a(j,k)=0.d0

! a(j,n_equi+Lpre+4)*Psi_b  ,j-number of equation
 
      k=n_equi+Lpre+4 
      a(j,k)=0.d0
! right hand side
      y(j)=-dpsx_r

!------------------------
! equation (.)*d(Lam_z)=0

! a(j,k)*I(k) ,j-number of equation

      j=n_equi+Lpre+2

      do k=1,n_equi 
         a(j,k)=Gx_z(k)
      enddo

! a(j,n_equi+l)*Lam(l) ,j-number of equation
           
      do ll=1,Lpre 
         k=n_equi+ll
         a(j,k)=0.d0
      enddo

      k=n_equi+Lpre+1
      a(j,k)=0.d0
      k=n_equi+Lpre+2
      a(j,k)=0.d0
      k=n_equi+Lpre+3
      a(j,k)=0.d0

! a(j,n_equi+Lpre+3)*Psi_b  ,j-number of equation
 
      k=n_equi+Lpre+4 
      a(j,k)=0.d0
! right hand side
      y(j)=-dpsx_z

!-----------------------
! equation (.)*d(Lamx)=0

! a(j,k)*I(k) ,j-number of equation

      j=n_equi+Lpre+3

      do k=1,n_equi 
         a(j,k)=(1.d0-alph)*Ginda(k)+alph*Gindx(k)
      enddo

! a(j,n_equi+l)*Lam(l) ,j-number of equation
           
      do ll=1,Lpre 
         k=n_equi+ll
         a(j,k)=0.d0
      enddo

      k=n_equi+Lpre+1
      a(j,k)=0.d0
      k=n_equi+Lpre+2
      a(j,k)=0.d0
      k=n_equi+Lpre+3
      a(j,k)=0.d0

! a(j,n_equi+Lpre+4)*Psi_b  ,j-number of equation
 
      k=n_equi+Lpre+4 
      a(j,k)=-1.d0

! right hand side
      y(j)=-(1.d0-alph)*psip_a-alph*psip_x

!-------------------------
! equation (.)*d(-psi_b)=0

! a(j,k)*I(k) ,j-number of equation

      j=n_equi+Lpre+4

      do k=1,n_equi 
         asum=0.d0
         do l=1,Lfit 
            asum=asum+w_wght(l)*Gindk(k,l)
         enddo
         a(j,k)=asum
      enddo

! a(j,n_equi+l)*Lam(l) ,j-number of equation

      do l=1,Lpre 
         k=n_equi+l
         a(j,k)=1.d0
      enddo

! a(j,n_equi+Lpre+1)*Lam_r ,j-number of equation
      a(j,n_equi+Lpre+1)=0.d0

! a(j,n_equi+Lpre+2)*Lam_z ,j-number of equation
      a(j,n_equi+Lpre+2)=0.d0

! a(j,n_equi+Lpre+3)*Lamx ,j-number of equation
      a(j,n_equi+Lpre+3)=1.d0

! a(j,n_equi+Lpre+4)*Psi_b  ,j-number of equation

      asum=0.d0
      do l=1,Lfit 
         asum=asum+w_wght(l)
      enddo
      a(j,n_equi+Lpre+4)=-asum

! right hand side
      asum=0.d0
      do l=1,Lfit 
         asum=asum+w_wght(l)*psifit(l)
      enddo
      y(j)=-asum

      call GE(n_equi+Lpre+4,ncf_p,A,Y,X,IP)

      do iq=1,n_equi
         pf_ceqw(iq)=x(iq)
      enddo

! Check
      do l=1,Lfit
         psictr(l)=psifit(l)
         do iq=1,n_equi
            psictr(l)=psictr(l)+Gindk(iq,l)*pf_ceqw(iq)
         enddo
      enddo

      gr_xp=dpsx_r
      gz_xp=dpsx_z
      ps_xp=psip_x
      ps_ma=psip_a
      do iq=1,n_equi
         gr_xp=gr_xp+Gx_r(iq)*pf_ceqw(iq)
         gz_xp=gz_xp+Gx_z(iq)*pf_ceqw(iq)
         ps_xp=ps_xp+Gindx(iq)*pf_ceqw(iq)
         ps_ma=ps_ma+Ginda(iq)*pf_ceqw(iq)
      enddo
      write(*,*) 'cursol:'
      write(*,*) 'gr_xp=',gr_xp
      write(*,*) 'gz_xp=',gz_xp

      do l=1,Lpre
         psictr(Lfit+l)=psipre(l)
         do iq=1,n_equi
            psictr(Lfit+l)=psictr(Lfit+l)+Gindp(iq,l)*pf_ceqw(iq)
         enddo
      enddo

      do iq=1,n_equi
         pf_ceqw(iq)=x(iq)/amu0
      enddo

      psi_bnd=x(n_equi+Lpre+4)

      return
      end subroutine cursol

!----------------------------------------------------------------
      subroutine cursol_snf(n_equi, pf_ceqw, psi_bnd)

      use sp_parameters, only: ncf_p, amu0
      use parcur

      implicit none

      integer, intent(in) :: n_equi
      real*8, intent(out) :: psi_bnd
      real*8, intent(out), dimension(*) :: pf_ceqw

      integer :: j, k, l, ll, iq
      integer, dimension(ncf_p) :: IP
      real*8 :: asum, gsext, d2rr_xp, d2rz_xp, d2zz_xp, gr_xp, gz_xp,
     &   ps_ma, ps_xp
      real*8, dimension(ncf_p) :: X, Y, psictr
      real*8, dimension(ncf_p, ncf_p) :: A

      do k=1,n_equi 
         do j=1,n_equi 
           asum=0.d0
           do l=1,Lfit 
              asum=asum+Gindk(k,l)*Gindk(j,l)*w_wght(l)
           enddo
           a(j,k)=asum
        enddo
      enddo

      do j=1,n_equi 
         a(j,j)=a(j,j)+d_wght(j)*s_wght
      enddo

      do k=1,n_equi 
         do j=1,n_equi 
            a(j,k)=a(j,k)+CxC(j,k)*c_wght
         enddo
      enddo

! a(j,n_equi+l)*Lam(l) ,j-number of equation

      do j=1,n_equi 
         do l=1,Lpre 
            a(j,n_equi+l)=Gindp(j,l)
         enddo
      enddo

! a(j,n_equi+Lpre+1)*Lam_r(l) ,j-number of equation

      do j=1,n_equi 
         a(j,n_equi+Lpre+1)=Gx_r(j)
      enddo

! a(j,n_equi+Lpre+2)*Lam_z(l) ,j-number of equation

      do j=1,n_equi 
         a(j,n_equi+Lpre+2)=Gx_z(j)
      enddo

! a(j,n_equi+Lpre+3)*Lamx ,j-number of equation

      do j=1,n_equi 
         a(j,n_equi+Lpre+3)=(1.d0-alph)*Ginda(j)+alph*Gindx(j)
      enddo

! a(j,n_equi+Lpre+4)*Psi_b  ,j-number of equation

      do j=1,n_equi 
         asum=0.d0
         do l=1,Lfit 
            asum=asum+w_wght(l)*Gindk(j,l)
         enddo
         a(j,n_equi+Lpre+4)=-asum
      enddo

! a(j,n_equi+Lpre+5)*Lam_rz(l) ,j-number of equation

      do j=1,n_equi 
         a(j,n_equi+Lpre+5)=Gx_rz(j)
      enddo

! a(j,n_equi+Lpre+6)*Lam_zz(l) ,j-number of equation

      do j=1,n_equi 
         a(j,n_equi+Lpre+6)=Gx_zz(j)
      enddo

! right hand side

      do j=1,n_equi 
         asum=0.d0
         do l=1,Lfit 
            asum=asum+w_wght(l)*Gindk(j,l)*psifit(l)
         enddo
         y(j)=-asum+d_wght(j)*s_wght*curref(j)*amu0
      enddo

!------------------------
! equation (.)*d(Lam_l)=0

! a(j,k)*I(k) ,j-number of equation

      do l=1,Lpre 
         j=n_equi+l
         do k=1,n_equi 
            a(j,k)=Gindp(k,l)
         enddo
      enddo

! a(j,n_equi+l)*Lam(l) ,j-number of equation
           
      do l=1,Lpre 
         j=n_equi+l
         do ll=1,Lpre 
            k=n_equi+ll
            a(j,k)=0.d0
         enddo
      enddo

      do l=1,Lpre 
         j=n_equi+l
         k=n_equi+Lpre+1
         a(j,k)=0.d0
         k=n_equi+Lpre+2
         a(j,k)=0.d0
         k=n_equi+Lpre+3
         a(j,k)=0.d0
      enddo

! a(j,n_equi+Lpre+4)*Psi_b  ,j-number of equation
 
      do l=1,Lpre 
         j=n_equi+l
         k=n_equi+Lpre+4 
         a(j,k)=-1.d0
      enddo
   
! a(j,n_equi+Lpre+5)*Lam_rz ,j-number of equation
! a(j,n_equi+Lpre+6)*Lam_zz ,j-number of equation

      do l=1,Lpre 
         j=n_equi+l
         k=n_equi+Lpre+5
         a(j,k)=0.d0
         k=n_equi+Lpre+6
         a(j,k)=0.d0
      enddo
  
! right hand side
      do l=1,Lpre 
         j=n_equi+l
         y(j)=-psipre(l)
      enddo

!------------------------
! equation (.)*d(Lam_r)=0

! a(j,k)*I(k) ,j-number of equation

      j=n_equi+Lpre+1

      do k=1,n_equi 
         a(j,k)=Gx_r(k)
      enddo

! a(j,n_equi+l)*Lam(l) ,j-number of equation
           
      do ll=1,Lpre 
         k=n_equi+ll
         a(j,k)=0.d0
      enddo

      k=n_equi+Lpre+1
      a(j,k)=0.d0
      k=n_equi+Lpre+2
      a(j,k)=0.d0
      k=n_equi+Lpre+3
      a(j,k)=0.d0

! a(j,n_equi+Lpre+4)*Psi_b  ,j-number of equation
 
      k=n_equi+Lpre+4 
      a(j,k)=0.d0
          
! a(j,n_equi+Lpre+5)*Lam_rz ,j-number of equation
! a(j,n_equi+Lpre+6)*Lam_zz ,j-number of equation

      k=n_equi+Lpre+5
      a(j,k)=0.d0
      k=n_equi+Lpre+6
      a(j,k)=0.d0
          
! right hand side
      y(j)=-dpsx_r

!-----------------------
! equation (.)*d(Lam_z)=0

! a(j,k)*I(k) ,j-number of equation

      j=n_equi+Lpre+2

      do k=1,n_equi 
         a(j,k)=Gx_z(k)
      enddo

! a(j,n_equi+l)*Lam(l) ,j-number of equation
           
      do ll=1,Lpre 
         k=n_equi+ll
         a(j,k)=0.d0
      enddo

      k=n_equi+Lpre+1
      a(j,k)=0.d0
      k=n_equi+Lpre+2
      a(j,k)=0.d0
      k=n_equi+Lpre+3
      a(j,k)=0.d0

! a(j,n_equi+Lpre+3)*Psi_b  ,j-number of equation
 
      k=n_equi+Lpre+4 
      a(j,k)=0.d0

      k=n_equi+Lpre+5
      a(j,k)=0.d0
      k=n_equi+Lpre+6
      a(j,k)=0.d0
      
! right hand side
      y(j)=-dpsx_z

!-----------------------
! equation (.)*d(Lamx)=0

! a(j,k)*I(k) ,j-number of equation

      j=n_equi+Lpre+3

      do k=1,n_equi 
         a(j,k)=(1.d0-alph)*Ginda(k)+alph*Gindx(k)
      enddo

! a(j,n_equi+l)*Lam(l) ,j-number of equation
           
      do ll=1,Lpre 
         k=n_equi+ll
         a(j,k)=0.d0
      enddo

      k=n_equi+Lpre+1
      a(j,k)=0.d0
      k=n_equi+Lpre+2
      a(j,k)=0.d0
      k=n_equi+Lpre+3
      a(j,k)=0.d0

! a(j,n_equi+Lpre+4)*Psi_b  ,j-number of equation
 
      k=n_equi+Lpre+4 
      a(j,k)=-1.d0

      k=n_equi+Lpre+5
      a(j,k)=0.d0
      k=n_equi+Lpre+6
      a(j,k)=0.d0

! right hand side
      y(j)=-(1.d0-alph)*psip_a-alph*psip_x

!-------------------------
! equation (.)*d(-psi_b)=0

! a(j,k)*I(k) ,j-number of equation

      j=n_equi+Lpre+4

      do k=1,n_equi 
         asum=0.d0
         do l=1,Lfit 
            asum=asum+w_wght(l)*Gindk(k,l)
         enddo
         a(j,k)=asum
      enddo

! a(j,n_equi+l)*Lam(l) ,j-number of equation

      do l=1,Lpre 
         k=n_equi+l
         a(j,k)=1.d0
      enddo

! a(j,n_equi+Lpre+1)*Lam_r ,j-number of equation
      a(j,n_equi+Lpre+1)=0.d0
! a(j,n_equi+Lpre+2)*Lam_z ,j-number of equation
      a(j,n_equi+Lpre+2)=0.d0
! a(j,n_equi+Lpre+3)*Lamx ,j-number of equation
      a(j,n_equi+Lpre+3)=1.d0

! a(j,n_equi+Lpre+4)*Psi_b  ,j-number of equation

      asum=0.d0
      do l=1,Lfit 
         asum=asum+w_wght(l)
      enddo
      a(j,n_equi+Lpre+4)=-asum

! a(j,n_equi+Lpre+5)*Lam_rz ,j-number of equation
      a(j,n_equi+Lpre+5)=0.d0
! a(j,n_equi+Lpre+2)*Lam_zz ,j-number of equation
      a(j,n_equi+Lpre+6)=0.d0

! right hand side
      asum=0.d0
      do l=1,Lfit 
         asum=asum+w_wght(l)*psifit(l)
      enddo
      y(j)=-asum

!-------------------------
! equation (.)*d(Lam_rz)=0

! a(j,k)*I(k) ,j-number of equation

      j=n_equi+Lpre+5

      do k=1,n_equi 
         a(j,k)=Gx_rz(k)
      enddo

! a(j,n_equi+l)*Lam(l) ,j-number of equation
       
      do ll=1,Lpre 
         k=n_equi+ll
         a(j,k)=0.d0
      enddo

      k=n_equi+Lpre+1
      a(j,k)=0.d0
      k=n_equi+Lpre+2
      a(j,k)=0.d0
      k=n_equi+Lpre+3
      a(j,k)=0.d0

! a(j,n_equi+Lpre+4)*Psi_b  ,j-number of equation
 
      k=n_equi+Lpre+4 
      a(j,k)=0.d0
! a(j,n_equi+Lpre+5)*Lam_rz  ,j-number of equation
      k=n_equi+Lpre+5 
      a(j,k)=0.d0
! a(j,n_equi+Lpre+6)*Lam_zz  ,j-number of equation
      k=n_equi+Lpre+6 
      a(j,k)=0.d0
! right hand side
       y(j)=-dpsx_rz

!------------------------
! equation (.)*d(Lam_zz)=0

! a(j,k)*I(k) ,j-number of equation

      j=n_equi+Lpre+6

      do k=1,n_equi 
         a(j,k)=Gx_zz(k)
      enddo

! a(j,n_equi+l)*Lam(l) ,j-number of equation
       
      do ll=1,Lpre 
         k=n_equi+ll
         a(j,k)=0.d0
      enddo

      k=n_equi+Lpre+1
      a(j,k)=0.d0
      k=n_equi+Lpre+2
      a(j,k)=0.d0
      k=n_equi+Lpre+3
      a(j,k)=0.d0

! a(j,n_equi+Lpre+4)*Psi_b  ,j-number of equation
      k=n_equi+Lpre+4 
      a(j,k)=0.d0

! a(j,n_equi+Lpre+5)*Lam_rz  ,j-number of equation
      k=n_equi+Lpre+5 

! a(j,n_equi+Lpre+6)*Lam_zz  ,j-number of equation
      a(j,k)=0.d0
      k=n_equi+Lpre+6 
      a(j,k)=0.d0

! right hand side
      y(j)=-dpsx_zz

      call GE(n_equi+Lpre+6,ncf_p,A,Y,X,IP)

      do iq=1,n_equi
         pf_ceqw(iq)=x(iq)
      enddo

! Check
      do l=1,Lfit
         psictr(l)=psifit(l)
         do iq=1,n_equi
            psictr(l)=psictr(l)+Gindk(iq,l)*pf_ceqw(iq)
         enddo
      enddo
    
      gr_xp=dpsx_r
      gz_xp=dpsx_z
      ps_xp=psip_x
      ps_ma=psip_a
      d2rz_xp=dpsx_rz
      d2zz_xp=dpsx_zz
      d2rr_xp=dpsx_rr

      GSExt=0.d0
          
      do iq=1,n_equi
         gr_xp=gr_xp+Gx_r(iq)*pf_ceqw(iq)
         gz_xp=gz_xp+Gx_z(iq)*pf_ceqw(iq)
         ps_xp=ps_xp+Gindx(iq)*pf_ceqw(iq)
         ps_ma=ps_ma+Ginda(iq)*pf_ceqw(iq)
  
         d2rz_xp=d2rz_xp+Gx_rz(iq)*pf_ceqw(iq)
         d2zz_xp=d2zz_xp+Gx_zz(iq)*pf_ceqw(iq)
         d2rr_xp=d2rr_xp+Gx_rr(iq)*pf_ceqw(iq)
          
         GSExt=GSExt+( Gx_zz(iq)+Gx_rr(iq)-Gx_r(iq)/rx_p )*pf_ceqw(iq)

      enddo

      write(*,*) 'cursol:'
      write(*,*) 'd2rz_xp=',d2rz_xp
      write(*,*) 'd2zz_xp=',d2zz_xp
      write(*,*) 'd2rr_xp=',d2rr_xp
      write(*,*) 'gr_xp=',gr_xp
      write(*,*) 'gz_xp=',gz_xp
      write(*,*) 'GSExt',GSExt

      do l=1,Lpre
         psictr(Lfit+l)=psipre(l)
         do iq=1,n_equi
            psictr(Lfit+l)=psictr(Lfit+l)+Gindp(iq,l)*pf_ceqw(iq)
         enddo
      enddo

      do iq=1,n_equi
         pf_ceqw(iq)=x(iq)/amu0
      enddo

      psi_bnd=x(n_equi+Lpre+4)

      return
      end subroutine cursol_snf

!----------------------------------------------------------------
      subroutine cursol_2x(n_equi, pf_ceqw, psi_bnd)

      use sp_parameters, only: ncf_p, amu0
      use parcur

      implicit none

      integer, intent(in) :: n_equi
      real*8, intent(out) :: psi_bnd
      real*8, intent(out), dimension(*) :: pf_ceqw

      integer :: j, k, l, ll, iq
      integer, dimension(ncf_p) :: IP
      real*8 :: asum, gr_xp, gr_xp2, gz_xp, gz_xp2,
     &   ps_ma, ps_xp
      real*8, dimension(ncf_p) :: X, Y, psictr
      real*8, dimension(ncf_p, ncf_p) :: A

      do k=1,n_equi 
         do j=1,n_equi 
            asum=0.d0
            do l=1,Lfit 
               asum=asum+Gindk(k,l)*Gindk(j,l)*w_wght(l)
            enddo
            a(j,k)=asum
         enddo
      enddo

      do j=1,n_equi 
         a(j,j)=a(j,j)+d_wght(j)*s_wght
      enddo

      do k=1,n_equi 
         do j=1,n_equi 
            a(j,k)=a(j,k)+CxC(j,k)*c_wght
         enddo
      enddo

! a(j,n_equi+l)*Lam(l) ,j-number of equation
      do j=1,n_equi 
         do l=1,Lpre 
            a(j,n_equi+l)=Gindp(j,l)
         enddo
      enddo

! a(j,n_equi+Lpre+1)*Lam_r(l) ,j-number of equation
      do j=1,n_equi 
         a(j,n_equi+Lpre+1)=Gx_r(j)
      enddo

! a(j,n_equi+Lpre+2)*Lam_z(l) ,j-number of equation
      do j=1,n_equi 
         a(j,n_equi+Lpre+2)=Gx_z(j)
      enddo

! a(j,n_equi+Lpre+3)*Lamx ,j-number of equation
      do j=1,n_equi 
         a(j,n_equi+Lpre+3)=(1.d0-alph)*Ginda(j)+alph*Gindx(j)
      enddo

! a(j,n_equi+Lpre+4)*Psi_b  ,j-number of equation
      do j=1,n_equi 
         asum=0.d0
         do l=1,Lfit 
            asum=asum+w_wght(l)*Gindk(j,l)
         enddo
         a(j,n_equi+Lpre+4)=-asum
      enddo

! a(j,n_equi+Lpre+5)*Lam_r2(l) ,j-number of equation
      do j=1,n_equi 
         a(j,n_equi+Lpre+5)=Gx2_r(j)
      enddo

! a(j,n_equi+Lpre+6)*Lam_z2(l) ,j-number of equation
      do j=1,n_equi 
         a(j,n_equi+Lpre+6)=Gx2_z(j)
      enddo

! right hand side
      do j=1,n_equi 
         asum=0.d0
         do l=1,Lfit 
            asum=asum+w_wght(l)*Gindk(j,l)*psifit(l)
         enddo
         y(j)=-asum+d_wght(j)*s_wght*curref(j)*amu0
      enddo

!------------------------
! equation (.)*d(Lam_l)=0

! a(j,k)*I(k) ,j-number of equation

      do l=1,Lpre 
         j=n_equi+l
         do k=1,n_equi 
            a(j,k)=Gindp(k,l)
         enddo
      enddo

! a(j,n_equi+l)*Lam(l) ,j-number of equation
           
      do l=1,Lpre 
         j=n_equi+l
         do ll=1,Lpre 
            k=n_equi+ll
            a(j,k)=0.d0
         enddo
      enddo

! a(j,n_equi+Lpre+1)*Lam_r ,j-number of equation
! a(j,n_equi+Lpre+2)*Lam_z ,j-number of equation
! a(j,n_equi+Lpre+3)*Lamx ,j-number of equation

      do l=1,Lpre 
         j=n_equi+l
         k=n_equi+Lpre+1
         a(j,k)=0.d0
         k=n_equi+Lpre+2
         a(j,k)=0.d0
         k=n_equi+Lpre+3
         a(j,k)=0.d0
      enddo

! a(j,n_equi+Lpre+4)*Psi_b  ,j-number of equation
 
      do l=1,Lpre 
         j=n_equi+l
         k=n_equi+Lpre+4 
         a(j,k)=-1.d0
      enddo
     
! a(j,n_equi+Lpre+5)*Lam_r2 ,j-number of equation
! a(j,n_equi+Lpre+6)*Lam_z2 ,j-number of equation

      do l=1,Lpre 
         j=n_equi+l
         k=n_equi+Lpre+5
         a(j,k)=0.d0
         k=n_equi+Lpre+6
         a(j,k)=0.d0
      enddo
   
! right hand side
      do l=1,Lpre 
         j=n_equi+l
         y(j)=-psipre(l)
      enddo

!------------------------
! equation (.)*d(Lam_r)=0

! a(j,k)*I(k) ,j-number of equation

      j=n_equi+Lpre+1

      do k=1,n_equi 
         a(j,k)=Gx_r(k)
      enddo

! a(j,n_equi+l)*Lam(l) ,j-number of equation
           
      do ll=1,Lpre 
         k=n_equi+ll
         a(j,k)=0.d0
      enddo

! a(j,n_equi+Lpre+1)*Lam_r ,j-number of equation
! a(j,n_equi+Lpre+2)*Lam_z ,j-number of equation
! a(j,n_equi+Lpre+3)*Lamx ,j-number of equation

      k=n_equi+Lpre+1
      a(j,k)=0.d0
      k=n_equi+Lpre+2
      a(j,k)=0.d0
      k=n_equi+Lpre+3
      a(j,k)=0.d0

! a(j,n_equi+Lpre+4)*Psi_b  ,j-number of equation
 
      k=n_equi+Lpre+4 
      a(j,k)=0.d0
        
! a(j,n_equi+Lpre+5)*Lam_r2 ,j-number of equation
! a(j,n_equi+Lpre+6)*Lam_z2 ,j-number of equation

      k=n_equi+Lpre+5
      a(j,k)=0.d0
      k=n_equi+Lpre+6
      a(j,k)=0.d0
        
! right hand side
      y(j)=-dpsx_r

!------------------------
! equation (.)*d(Lam_z)=0

! a(j,k)*I(k) ,j-number of equation

      j=n_equi+Lpre+2

      do k=1,n_equi 
         a(j,k)=Gx_z(k)
      enddo

! a(j,n_equi+l)*Lam(l) ,j-number of equation
         
      do ll=1,Lpre 
         k=n_equi+ll
         a(j,k)=0.d0
      enddo

! a(j,n_equi+Lpre+1)*Lam_r ,j-number of equation
! a(j,n_equi+Lpre+2)*Lam_z ,j-number of equation
! a(j,n_equi+Lpre+3)*Lamx ,j-number of equation

      k=n_equi+Lpre+1
      a(j,k)=0.d0
      k=n_equi+Lpre+2
      a(j,k)=0.d0
      k=n_equi+Lpre+3
      a(j,k)=0.d0

! a(j,n_equi+Lpre+3)*Psi_b  ,j-number of equation

      k=n_equi+Lpre+4 
      a(j,k)=0.d0
  
! a(j,n_equi+Lpre+5)*Lam_r2 ,j-number of equation
! a(j,n_equi+Lpre+6)*Lam_z2 ,j-number of equation

      k=n_equi+Lpre+5
      a(j,k)=0.d0
      k=n_equi+Lpre+6
      a(j,k)=0.d0

! right hand side
      y(j)=-dpsx_z

!------------------------
! equation (.)*d(Lamx)=0

! a(j,k)*I(k) ,j-number of equation

      j=n_equi+Lpre+3

      do k=1,n_equi 
         a(j,k)=(1.d0-alph)*Ginda(k)+alph*Gindx(k)
      enddo

! a(j,n_equi+l)*Lam(l) ,j-number of equation
           
      do ll=1,Lpre 
         k=n_equi+ll
         a(j,k)=0.d0
      enddo

! a(j,n_equi+Lpre+1)*Lam_r ,j-number of equation
! a(j,n_equi+Lpre+2)*Lam_z ,j-number of equation
! a(j,n_equi+Lpre+3)*Lamx ,j-number of equation

      k=n_equi+Lpre+1
      a(j,k)=0.d0
      k=n_equi+Lpre+2
      a(j,k)=0.d0
      k=n_equi+Lpre+3
      a(j,k)=0.d0

! a(j,n_equi+Lpre+4)*Psi_b  ,j-number of equation
 
      k=n_equi+Lpre+4 
      a(j,k)=-1.d0

! a(j,n_equi+Lpre+5)*Lam_r2 ,j-number of equation
! a(j,n_equi+Lpre+6)*Lam_z2 ,j-number of equation

      k=n_equi+Lpre+5
      a(j,k)=0.d0
      k=n_equi+Lpre+6
      a(j,k)=0.d0

! right hand side
      y(j)=-(1.d0-alph)*psip_a-alph*psip_x

!-------------------------
! equation (.)*d(-psi_b)=0

! a(j,k)*I(k) ,j-number of equation

      j=n_equi+Lpre+4

      do k=1,n_equi 
         asum=0.d0
         do l=1,Lfit 
            asum=asum+w_wght(l)*Gindk(k,l)
         enddo
         a(j,k)=asum
      enddo

! a(j,n_equi+l)*Lam(l) ,j-number of equation

      do l=1,Lpre 
         k=n_equi+l
         a(j,k)=1.d0
      enddo

! a(j,n_equi+Lpre+1)*Lam_r ,j-number of equation
      a(j,n_equi+Lpre+1)=0.d0

! a(j,n_equi+Lpre+2)*Lam_z ,j-number of equation
      a(j,n_equi+Lpre+2)=0.d0

! a(j,n_equi+Lpre+3)*Lamx ,j-number of equation
      a(j,n_equi+Lpre+3)=1.d0

! a(j,n_equi+Lpre+4)*Psi_b  ,j-number of equation

      asum=0.d0
      do l=1,Lfit 
         asum=asum+w_wght(l)
      enddo
      a(j,n_equi+Lpre+4)=-asum

! a(j,n_equi+Lpre+5)*Lam_r2 ,j-number of equation
      a(j,n_equi+Lpre+5)=0.d0

! a(j,n_equi+Lpre+2)*Lam_z2 ,j-number of equation
      a(j,n_equi+Lpre+6)=0.d0

! right hand side
      asum=0.d0
      do l=1,Lfit 
         asum=asum+w_wght(l)*psifit(l)
      enddo
      y(j)=-asum

!------------------------
! equation (.)*d(Lam_r2)=0

! a(j,k)*I(k) ,j-number of equation

      j=n_equi+Lpre+5

      do k=1,n_equi 
         a(j,k)=Gx2_r(k)
      enddo

! a(j,n_equi+l)*Lam(l) ,j-number of equation
           
      do ll=1,Lpre 
         k=n_equi+ll
         a(j,k)=0.d0
      enddo

! a(j,n_equi+Lpre+1)*Lam_r ,j-number of equation
! a(j,n_equi+Lpre+2)*Lam_z ,j-number of equation
! a(j,n_equi+Lpre+3)*Lamx ,j-number of equation

      k=n_equi+Lpre+1
      a(j,k)=0.d0
      k=n_equi+Lpre+2
      a(j,k)=0.d0
      k=n_equi+Lpre+3
      a(j,k)=0.d0

! a(j,n_equi+Lpre+4)*Psi_b  ,j-number of equation
      k=n_equi+Lpre+4 
      a(j,k)=0.d0

! a(j,n_equi+Lpre+5)*Lam_r2  ,j-number of equation
      k=n_equi+Lpre+5 
      a(j,k)=0.d0

! a(j,n_equi+Lpre+6)*Lam_z2  ,j-number of equation
      k=n_equi+Lpre+6 
      a(j,k)=0.d0

! right hand side
      y(j)=-dpsx2_r

!------------------------
! equation (.)*d(Lam_z2)=0

! a(j,k)*I(k) ,j-number of equation

      j=n_equi+Lpre+6

      do k=1,n_equi 
         a(j,k)=Gx2_z(k)
      enddo

! a(j,n_equi+l)*Lam(l) ,j-number of equation
       
      do ll=1,Lpre 
         k=n_equi+ll
         a(j,k)=0.d0
      enddo

! a(j,n_equi+Lpre+1)*Lam_r ,j-number of equation
! a(j,n_equi+Lpre+2)*Lam_z ,j-number of equation
! a(j,n_equi+Lpre+3)*Lamx ,j-number of equation

      k=n_equi+Lpre+1
      a(j,k)=0.d0
      k=n_equi+Lpre+2
      a(j,k)=0.d0
      k=n_equi+Lpre+3
      a(j,k)=0.d0

! a(j,n_equi+Lpre+4)*Psi_b  ,j-number of equation
      k=n_equi+Lpre+4 
      a(j,k)=0.d0

! a(j,n_equi+Lpre+5)*Lam_r2  ,j-number of equation
      k=n_equi+Lpre+5 

! a(j,n_equi+Lpre+6)*Lam_z2  ,j-number of equation
      a(j,k)=0.d0
      k=n_equi+Lpre+6 
      a(j,k)=0.d0

! right hand side
      y(j)=-dpsx2_z

      call GE(n_equi+Lpre+6,ncf_p,A,Y,X,IP)

      do iq=1,n_equi
         pf_ceqw(iq)=x(iq)
      enddo

! Check
      do l=1,Lfit
         psictr(l)=psifit(l)
         do iq=1,n_equi
            psictr(l)=psictr(l)+Gindk(iq,l)*pf_ceqw(iq)
         enddo
      enddo

      gr_xp=dpsx_r
      gz_xp=dpsx_z
      ps_xp=psip_x
      ps_ma=psip_a
      gr_xp2=dpsx2_r
      gz_xp2=dpsx2_z

      do iq=1,n_equi 
         gr_xp=gr_xp+Gx_r(iq)*pf_ceqw(iq)
         gz_xp=gz_xp+Gx_z(iq)*pf_ceqw(iq)

         ps_xp=ps_xp+Gindx(iq)*pf_ceqw(iq)
         ps_ma=ps_ma+Ginda(iq)*pf_ceqw(iq)
    
         gr_xp2=gr_xp2+Gx2_r(iq)*pf_ceqw(iq)
         gz_xp2=gz_xp2+Gx2_z(iq)*pf_ceqw(iq)   
      enddo

      write(*,*) 'cursol:'
      write(*,*) 'gr_xp=',gr_xp
      write(*,*) 'gz_xp=',gz_xp
      write(*,*) 'gr_xp2=',gr_xp2
      write(*,*) 'gz_xp2=',gz_xp2

      do l=1,Lpre
         psictr(Lfit+l)=psipre(l)
         do iq=1,n_equi
            psictr(Lfit+l)=psictr(Lfit+l)+Gindp(iq,l)*pf_ceqw(iq)
         enddo
      enddo

      do iq=1,n_equi
         pf_ceqw(iq)=x(iq)/amu0
      enddo

      psi_bnd=x(n_equi+Lpre+4)

      return
      end subroutine cursol_2x

!----------------------------------------------------------------
      subroutine cursol_wx(n_equi, pf_ceqw, psi_bnd)

      use sp_parameters, only: ncf_p, amu0
      use parcur

      implicit none

      integer, intent(in) :: n_equi
      real*8, intent(out) :: psi_bnd
      real*8, intent(out), dimension(*) :: pf_ceqw

      integer :: j, k, l, ll, iq
      integer, dimension(ncf_p) :: IP
      real*8 :: asum, gr_xp, gz_xp, ps_ma, ps_xp
      real*8, dimension(ncf_p) :: X, Y, psictr
      real*8, dimension(ncf_p, ncf_p) :: A

      do k=1,n_equi 
         do j=1,n_equi 
            asum=0.d0
            do l=1,Lfit 
               asum=asum+Gindk(k,l)*Gindk(j,l)*w_wght(l)
            enddo
            a(j,k)=asum
         enddo
      enddo

      do j=1,n_equi 
         a(j,j)=a(j,j)+d_wght(j)*s_wght
      enddo

! a(j,n_equi+l)*Lam(l) ,j-number of equation
      do j=1,n_equi
         do l=1,Lpre 
            a(j,n_equi+l)=Gindp(j,l)
         enddo
      enddo

      do j=1,n_equi 
         a(j,n_equi+Lpre+1)=(1.d0-alph)*Ginda(j)+alph*Gindx(j)
      enddo

! a(j,n_equi+Lpre+4)*Psi_b  ,j-number of equation

      do j=1,n_equi 
         asum=0.d0
         do l=1,Lfit 
            asum=asum+w_wght(l)*Gindk(j,l)
         enddo
         a(j,n_equi+Lpre+2)=-asum
      enddo

! right hand side

      do j=1,n_equi 
         asum=0.d0
         do l=1,Lfit 
            asum=asum+w_wght(l)*Gindk(j,l)*psifit(l)
         enddo
         y(j)=-asum+d_wght(j)*s_wght*curref(j)*amu0
      enddo

!------------------------
! equation (.)*d(Lam_l)=0

! a(j,k)*I(k) ,j-number of equation
      do l=1,Lpre 
         j=n_equi+l
         do k=1,n_equi 
            a(j,k)=Gindp(k,l)
         enddo
      enddo

! a(j,n_equi+l)*Lam(l) ,j-number of equation
           
      do l=1,Lpre 
         j=n_equi+l
         do ll=1,Lpre 
            k=n_equi+ll
            a(j,k)=0.d0
         enddo
      enddo

      do l=1,Lpre 
         j=n_equi+l
         k=n_equi+Lpre+1
         a(j,k)=0.d0
      enddo

! a(j,n_equi+Lpre+2)*Psi_b  ,j-number of equation
      do l=1,Lpre 
         j=n_equi+l
         k=n_equi+Lpre+2 
         a(j,k)=-1.d0
      enddo

! right hand side
      do l=1,Lpre 
         j=n_equi+l
         y(j)=-psipre(l)
      enddo

!------------------------
! equation (.)*d(Lamx)=0

! a(j,k)*I(k) ,j-number of equation
      j=n_equi+Lpre+1
      do k=1,n_equi 
         a(j,k)=(1.d0-alph)*Ginda(k)+alph*Gindx(k)
      enddo

! a(j,n_equi+l)*Lam(l) ,j-number of equation        
      do ll=1,Lpre
         k=n_equi+ll
         a(j,k)=0.d0
      enddo

      k=n_equi+Lpre+1
      a(j,k)=0.d0

! a(j,n_equi+Lpre+2)*Psi_b  ,j-number of equation
      k=n_equi+Lpre+2 
      a(j,k)=-1.d0

! right hand side
      y(j)=-(1.d0-alph)*psip_a-alph*psip_x

!------------------------
! equation (.)*d(-psi_b)=0

! a(j,k)*I(k) ,j-number of equation
      j=n_equi+Lpre+2
      do k=1,n_equi 
         asum=0.d0
         do l=1,Lfit 
            asum=asum+w_wght(l)*Gindk(k,l)
         enddo
         a(j,k)=asum
      enddo

! a(j,n_equi+l)*Lam(l) ,j-number of equation
      do l=1,Lpre
         k=n_equi+l
         a(j,k)=1.d0
      enddo

      a(j,n_equi+Lpre+1)=1.d0

! a(j,n_equi+Lpre+2)*Psi_b  ,j-number of equation
      asum=0.d0
      do l=1,Lfit 
         asum=asum+w_wght(l)
      enddo
      a(j,n_equi+Lpre+2)=-asum

! right hand side
      asum=0.d0
      do l=1,Lfit 
         asum=asum+w_wght(l)*psifit(l)
      enddo
      y(j)=-asum

      call GE(n_equi+Lpre+2,ncf_p,A,Y,X,IP)

      do iq=1,n_equi
         pf_ceqw(iq)=x(iq)
      enddo

      do l=1,Lfit
         psictr(l)=psifit(l)
         do iq=1,n_equi
            psictr(l)=psictr(l)+Gindk(iq,l)*pf_ceqw(iq)
         enddo
      enddo

      gr_xp=dpsx_r
      gz_xp=dpsx_z
      ps_xp=psip_x
      ps_ma=psip_a
      do iq=1,n_equi
         gr_xp=gr_xp+Gx_r(iq)*pf_ceqw(iq)
         gz_xp=gz_xp+Gx_z(iq)*pf_ceqw(iq)
         ps_xp=ps_xp+Gindx(iq)*pf_ceqw(iq)
         ps_ma=ps_ma+Ginda(iq)*pf_ceqw(iq)
      enddo

      do l=1,Lpre
          psictr(Lfit+l)=psipre(l)
          do iq=1,n_equi
             psictr(Lfit+l)=psictr(Lfit+l)+Gindp(iq,l)*pf_ceqw(iq)
          enddo
      enddo

      do iq=1,n_equi
         pf_ceqw(iq)=x(iq)/amu0
      enddo

      psi_bnd=x(n_equi+Lpre+2)

      return
      end subroutine cursol_wx

!------------------------------------------------------
      subroutine psi_cpn_wx

      use parcur
      use comblc, only: rmin, zmin, rx0, zx0, dr, dz, r, z, ui

      implicit none

      integer, parameter :: nshp=10

      integer :: i, j, k, l, ic, jc, ik, jk, nsh
      real*8 :: r0, z0, rrx, zzx
      real*8, dimension(5) :: dp
      real*8, dimension(nshp) :: xs, ys, fun
      real*8, external :: blin

      do l=1,Lfit
         r0=rfit(l)
         z0=zfit(l)
         ic=(r0-rmin)/dr(1)+1
         jc=(z0-zmin)/dz(1)+1
         psifit(l)=blin(ic, jc, r0, z0)
      enddo 

      do l=1,Lpre
         r0=rpre(l)
         z0=zpre(l)
         ic=(r0-rmin)/dr(1)+1
         jc=(z0-zmin)/dz(1)+1
         psipre(l)=blin(ic, jc, r0, z0)
      enddo 

      call closest_knot(rx0, zx0, ik, jk)

      nsh=1
      xs(nsh)=r(ik)
      ys(nsh)=z(jk)
      fun(nsh)=ui(ik,jk)

      do k=-1,1
         i= ik+k
         do l=-1,1
            j= jk+l
            if(k.eq.0  .AND. l.eq.0 ) CYCLE
            nsh=nsh+1
            xs(nsh)=r(i)
            ys(nsh)=z(j)
            fun(nsh)=ui(i,j)
         enddo
      enddo

      rrx=xs(1)
      zzx=ys(1)

      call deriv5(xs, ys, fun, nsh, 5, dp)

      dpsx_r = dp(1) + dp(3)*(rx0-rrx) + dp(4)*(zx0-zzx)
      dpsx_z = dp(2) + dp(5)*(zx0-zzx) + dp(4)*(rx0-rrx)      
      dpsx_rz = dp(4)
      dpsx_zz = dp(5)
      dpsx_rr = dp(3)

      psip_x=fun(1)+ dp(1)*(rx0-rrx) + dp(2)*(zx0-zzx) +
     &         0.5d0*dp(3)*(rx0-rrx)*(rx0-rrx) +
     &               dp(4)*(rx0-rrx)*(zx0-zzx) +
     &         0.5d0*dp(5)*(zx0-zzx)*(zx0-zzx)

      return
      end subroutine psi_cpn_wx

!----------------------------------------------------------------
      subroutine psi_cpn

      use parcur
      use comblc, only: rmin, zmin, dr, dz, r, z, ui, ni, nj

      implicit none

      integer, parameter :: nshp=26

      integer :: i, j, k, l, ic, ik, jc, jk, nsh
      real*8 :: r0, z0, rrx, zzx
      real*8, dimension(9) :: dp
      real*8, dimension(nshp) :: xs, ys, fun
      real*8, external :: blin

      do l=1,Lfit
         r0=rfit(l)
         z0=zfit(l)
         ic=(r0-rmin)/dr(1)+1
         jc=(z0-zmin)/dz(1)+1
         psifit(l)=blin(ic,jc,r0,z0)
      enddo 

      do l=1,Lpre
         r0=rpre(l)
         z0=zpre(l)
         ic=(r0-rmin)/dr(1)+1
         jc=(z0-zmin)/dz(1)+1
         psipre(l)=blin(ic,jc,r0,z0)
      enddo

      call closest_knot(rx_p, zx_p, ik, jk)

      nsh=1
      xs(nsh)=r(ik)
      ys(nsh)=z(jk)
      fun(nsh)=ui(ik,jk)

      do k=-2,2
         i= ik+k
         if(i.le.0 .OR. i.gt.ni) return
         do l=-2,2
            j= jk+l
            if(j.le.0 .OR. j.gt.nj) return
            if(k.eq.0  .AND. l.eq.0 ) CYCLE
            nsh=nsh+1
            xs(nsh)=r(i)
            ys(nsh)=z(j)
            fun(nsh)=ui(i,j)
         enddo
      enddo

      rrx=xs(1)
      zzx=ys(1)

      call deriv9(xs,ys,fun,nsh,9,dp)

      dpsx_r = dp(1) + dp(3)*(rx_p-rrx) + dp(4)*(zx_p-zzx) +
     &         dp(6)*(rx_p-rrx)**2*3+dp(8)*(zx_p-zzx)**2   +
     &         dp(7)*(rx_p-rrx)*(zx_p-zzx)*2
      dpsx_z = dp(2) + dp(5)*(zx_p-zzx) + dp(4)*(rx_p-rrx) +
     &         dp(7)*(rx_p-rrx)**2+dp(9)*(zx_p-zzx)**2*3   +
     &         dp(8)*(rx_p-rrx)*(zx_p-zzx)*2
      dpsx_rz = dp(4) + dp(7)*(rx_p-rrx)*2+dp(8)*(zx_p-zzx)*2
      dpsx_zz = dp(5) + dp(8)*(rx_p-rrx)*2+dp(9)*(zx_p-zzx)*6
      dpsx_rr = dp(3) + dp(7)*(zx_p-zzx)*2+dp(6)*(rx_p-rrx)*6

      psip_x=fun(1)+ dp(1)*(rx_p-rrx) + dp(2)*(zx_p-zzx)    +
     &         0.5d0*dp(3)*(rx_p-rrx)*(rx_p-rrx)            +
     &               dp(4)*(rx_p-rrx)*(zx_p-zzx)            +
     &         0.5d0*dp(5)*(zx_p-zzx)*(zx_p-zzx)            +
     &               dp(6)*(rx_p-rrx)*(rx_p-rrx)*(rx_p-rrx) +
     &               dp(7)*(rx_p-rrx)*(rx_p-rrx)*(zx_p-zzx) +
     &               dp(8)*(rx_p-rrx)*(zx_p-zzx)*(zx_p-zzx) +
     &               dp(9)*(zx_p-zzx)*(zx_p-zzx)*(zx_p-zzx)

      return
      end subroutine psi_cpn

!----------------------------------------------------------------
      subroutine psi_cpn_2x

      use parcur
      use comblc, only: r, z, ui

      implicit none

      integer, parameter :: nshp=26

      integer :: i, j, k, l, ik, jk, nsh
      real*8 :: rrx, zzx
      real*8, dimension(9) :: dp
      real*8, dimension(nshp) :: xs, ys, fun
      real*8, external :: blin

      call closest_knot(rx2_p, zx2_p, ik, jk)

      nsh=1
      xs(nsh)=r(ik)
      ys(nsh)=z(jk)
      fun(nsh)=ui(ik,jk)

      do k=-2,2
         i= ik+k
         do l=-2,2
            j= jk+l
            if(k.eq.0  .AND. l.eq.0 ) CYCLE
            nsh=nsh+1
            xs(nsh)=r(i)
            ys(nsh)=z(j)
            fun(nsh)=ui(i,j)
         enddo
      enddo

      rrx=xs(1)
      zzx=ys(1)

      call deriv9(xs,ys,fun,nsh,9,dp)

      dpsx2_r = dp(1) + dp(3)*(rx2_p-rrx) + dp(4)*(zx2_p-zzx)       +
     &                  dp(6)*(rx2_p-rrx)**2*3+dp(8)*(zx2_p-zzx)**2 +
     &                  dp(7)*(rx2_p-rrx)*(zx2_p-zzx)*2

      dpsx2_z = dp(2) + dp(5)*(zx2_p-zzx) + dp(4)*(rx2_p-rrx)       +
     &                  dp(7)*(rx2_p-rrx)**2+dp(9)*(zx2_p-zzx)**2*3 +
     &                  dp(8)*(rx2_p-rrx)*(zx2_p-zzx)*2

      psip_x2= fun(1)+ dp(1)*(rx2_p-rrx) + dp(2)*(zx2_p-zzx)     +
     &           0.5d0*dp(3)*(rx2_p-rrx)*(rx2_p-rrx)             +
     &                 dp(4)*(rx2_p-rrx)*(zx2_p-zzx)             +
     &           0.5d0*dp(5)*(zx2_p-zzx)*(zx2_p-zzx)             +
     &                 dp(6)*(rx2_p-rrx)*(rx2_p-rrx)*(rx2_p-rrx) +
     &                 dp(7)*(rx2_p-rrx)*(rx2_p-rrx)*(zx2_p-zzx) +
     &                 dp(8)*(rx2_p-rrx)*(zx2_p-zzx)*(zx2_p-zzx) +
     &                 dp(9)*(zx2_p-zzx)*(zx2_p-zzx)*(zx2_p-zzx)

      return
      end subroutine psi_cpn_2x

!----------------------------------------------------------------
      subroutine cursol_(n_equi, pf_ceqw, psi_bnd)

      use sp_parameters, only: ncf_p, amu0
      use parcur

      implicit none

      integer, intent(in) :: n_equi
      real*8, intent(out) :: psi_bnd
      real*8, intent(out), dimension(*) :: pf_ceqw

      integer :: j, k, l, ll, iq
      integer, dimension(ncf_p) :: IP
      real*8 :: asum, gr_xp, gz_xp, ps_ma, ps_xp
      real*8, dimension(ncf_p) :: X, Y, psictr
      real*8, dimension(ncf_p, ncf_p) :: A

!----------------------
! equation (.)*d(I_j)=0

! a(j,k)*I(k) ,j-number of equation

      do k=1,n_equi 
         do j=1,n_equi 
            asum=0.d0
            do l=1,Lfit 
              asum=asum+Gindk(k,l)*Gindk(j,l)*w_wght(l)
            enddo
            a(j,k)=asum
         enddo
      enddo

      do j=1,n_equi 
         a(j,j)=a(j,j)+d_wght(j)*s_wght
      enddo

! a(j,n_equi+l)*Lam(l) ,j-number of equation

      do j=1,n_equi 
         do l=1,Lpre 
            a(j,n_equi+l)=Gindp(j,l)
         enddo
      enddo

! a(j,n_equi+Lpre+1)*Psi_b  ,j-number of equation

      do j=1,n_equi 
         asum=0.d0
         do l=1,Lfit 
            asum=asum+w_wght(l)*Gindk(j,l)
         enddo
         a(j,n_equi+Lpre+1)=-asum
      enddo

      do j=1,n_equi 
         asum=0.d0
         do l=1,Lfit 
            asum=asum+w_wght(l)*Gindk(j,l)*psifit(l)
         enddo
         y(j)=-asum+d_wght(j)*s_wght*curref(j)*amu0
      enddo

!------------------------
! equation (.)*d(Lam_l)=0

! a(j,k)*I(k) ,j-number of equation
      do l=1,Lpre 
         j=n_equi+l
         do k=1,n_equi 
            a(j,k)=Gindp(k,l)
         enddo
      enddo

! a(j,n_equi+l)*Lam(l) ,j-number of equation     
      do l=1,Lpre 
         j=n_equi+l
         do ll=1,Lpre 
            k=n_equi+ll
            a(j,k)=0.d0
         enddo
      enddo

! a(j,n_equi+Lpre+1)*Psi_b  ,j-number of equation
      do l=1,Lpre 
         j=n_equi+l
         k=n_equi+Lpre+1 
         a(j,k)=-1.d0
      enddo

      do l=1,Lpre 
         j=n_equi+l
         y(j)=-psipre(l)
      enddo

!-------------------------
! equation (.)*d(-psi_b)=0
!
! a(j,k)*I(k) ,j-number of equation
      do k=1,n_equi 
         asum=0.d0
         do l=1,Lfit 
            asum=asum+w_wght(l)*Gindk(k,l)
         enddo
         a(n_equi+Lpre+1,k)=asum
      enddo

! a(j,n_equi+l)*Lam(l) ,j-number of equation
      do l=1,Lpre 
         k=n_equi+l
         a(n_equi+Lpre+1,k)=1.d0
      enddo

! a(j,n_equi+Lpre+1)*Psi_b  ,j-number of equation
      asum=0.d0
      do l=1,Lfit 
         asum=asum+w_wght(l)
      enddo
      a(n_equi+Lpre+1,n_equi+Lpre+1)=-asum

      asum=0.d0
      do l=1,Lfit 
         asum=asum+w_wght(l)*psifit(l)
      enddo
      y(n_equi+Lpre+1)=-asum

      call GE(n_equi+Lpre+1,ncf_p,A,Y,X,IP)

      do iq=1,n_equi
         pf_ceqw(iq)=x(iq)
      enddo

      do l=1,Lfit
         psictr(l)=psifit(l)
         do iq=1,n_equi
            psictr(l)=psictr(l)+Gindk(iq,l)*pf_ceqw(iq)
         enddo
      enddo

      do l=1,Lpre
         psictr(Lfit+l)=psipre(l)
         do iq=1,n_equi
            psictr(Lfit+l)=psictr(Lfit+l)+Gindp(iq,l)*pf_ceqw(iq)
         enddo
      enddo

      do iq=1,n_equi
         pf_ceqw(iq)=x(iq)/amu0
      enddo

      psi_bnd=x(n_equi+Lpre+1)

      return
      end subroutine cursol_

!----------------------------------------------------------------
      subroutine bonpsi

      use sp_parameters
      use parcur
      use compol, only: iplas, nt1, psi, dlt, nt, r, z

      implicit none

      integer :: i, j, l, jb
      real*8 :: a1, a2, a3, g1, g2, g3, dltk, dgdnl, psb, fint,
     &   r0, r1, rr, z0, z1, zz
      real*8, dimension(ntp) :: dgdn
      real*8, dimension(ntp, ntp) :: binadg
      real*8, dimension(nrp, ntp) :: a12, a23, a34, a14, a13, a24
      common/comaaa/ a12, a23, a34, a14, a13, a24

      i=iplas

      do j=2,nt1
         a1=a13(i-1,j-1)
         a2=a34(i-1,j-1)+a12(i-1,j)
         a3=a24(i-1,j)

         g1=psi(i-1,j-1)
         g2=psi(i-1,j)
         g3=psi(i-1,j+1)

         dltk=(dlt(i,j-1)+dlt(i,j))*0.5d0
         dgdnl=a1*g1+a2*g2+a3*g3
         dgdn(j)=dgdnl/dltk
      enddo

      dgdn(1)=dgdn(nt1)
      dgdn(nt)=dgdn(2)

      do l=1,Lfit
         rr=rfit(l)
         zz=zfit(l)
         do jb=1,nt1
            r0=r(iplas,jb)
            z0=z(iplas,jb)
            r1=r(iplas,jb+1)
            z1=z(iplas,jb+1)
            call bint(rr,zz,R0,Z0,r1,z1,Fint,1)
            binadg(jb,l)=fint
         enddo
      enddo

      do l=1,Lfit
         psb=0.d0
         do jb=2,nt1
            psb=psb+binadg(jb,l)*(dgdn(jb)+dgdn(jb+1))*0.5d0
         enddo
         psifit(l)=-psb
      enddo

      do l=1,Lpre
         rr=rpre(l)
         zz=zpre(l)
         do jb=1,nt1
            r0=r(iplas,jb)
            z0=z(iplas,jb)
            r1=r(iplas,jb+1)
            z1=z(iplas,jb+1)
            call bint(rr,zz,R0,Z0,r1,z1,Fint,1)
            binadg(jb,l)=fint
         enddo
      enddo

      do l=1,Lpre
         psb=0.d0
         do jb=2,nt1
            psb=psb+binadg(jb,l)*(dgdn(jb)+dgdn(jb+1))*0.5d0
         enddo
         psipre(l)=-psb
      enddo

      return
      end subroutine bonpsi
