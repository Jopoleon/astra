       subroutine spider(nstep, time, dt, key_dmf, k_grid, k_auto,
     &                   k_fixfree, dpsdt, key_start, voltpf)  

!-----------------------------------------------
!   MAIN  PROGRAM  OF  THE EVOLUTION CODE  "PET"
!-----------------------------------------------

      use fenix_params, only: s_adapt, s_fazt
      use keys, only: key_0st, kpr, key_out

      implicit none

      integer, intent(in) :: key_dmf, k_grid, k_fixfree, key_start
      real*8, intent(in) :: time, dt, dpsdt
      integer, intent(out) :: nstep, k_auto
      double precision, intent(out), dimension(*) :: voltpf

      integer :: plasma_up, plasma_trig, j_switch, kluch, yesfitcc, 
     &   k_dmf, nstep_local, nnstep, kkey_dmf

      real*8 :: rax, zax
      real*8, dimension(2500) :: contvals_mat
      real*8, dimension(500)  :: d_pf_mat, d_cam_mat
      double precision :: ibkdw, ifbey

      data j_switch /0/
      save j_switch
      save nstep_local

      common /yes_fit_cc/ yesfitcc
      common /plasma_state/ plasma_up, plasma_trig
      common /yes_bkdw/ ibkdw, ifbey

      s_fazt = 0
      if (ibkdw*ifbey.eq.0) then			 

         call put_tim(dt, time)
         call put_key_fix(k_fixfree)

         if(nstep.eq.0) then
            KLUCH = 0
         else
            KLUCH = 1
         endif

         if(k_fixfree.eq.0) then    
            call B_STEPON( KLUCH, k_auto, nstep, dt, time,
     &                     rax, zax, key_dmf, dpsdt)
            if(key_0st.eq.1 .AnD. nstep.eq.0) then
               nnstep=1
               kkey_dmf=-10        
               KLUCH = 1
               call B_STEPON( KLUCH, k_auto, nnstep, dt, time,
     &                        rax, zax, kkey_dmf, dpsdt)
               nstep=0
            endif

            return
         endif

! k_grid= 0   rect. grid
! k_grid= 1   adap. grid

         if(KLUCH.eq.0) then    !initialization

! basic free bound rectan, equilibrium ( KLUCH=0 )

            if(key_start .le. 0) then
               write(*,*) 'fix grid call'
               call sstepon( KLUCH, k_auto, nstep, dt, time,
     &                       voltpf, d_pf_mat, d_cam_mat, key_dmf)    
               if(key_out.gt.0) then
                  call spidupdate
               endif
            else
               write(*,*) 'fix grid call init'
               call cf_init( k_auto, nstep, dt, time,
     &                       voltpf, d_pf_mat, d_cam_mat )
               if(key_out.gt.0) then
                  call spidupdate
               endif
            endif

            call wrrec

            write(*,*) 'after stepon initial'

            k_auto= 0  ! don't change!

            if(k_grid.eq.1) then
               call f_stepon( KLUCH, k_auto, nstep, dt, time, voltpf,
     &                        d_pf_mat, d_cam_mat, rax, zax, key_dmf)
               if(key_out.gt.0) then
                  call f_spidupdate
               endif
            endif

         elseif(KLUCH.eq.1) then   !time steping

            if(k_grid.eq.0) then

               if (s_fazt.eq.2) kpr=-1
               if (s_fazt.eq.3) kpr=-2
               call sstepon( KLUCH, k_auto, nstep, dt, time, voltpf,
     &                       d_pf_mat, d_cam_mat, key_dmf)
               if(key_out.gt.0) then
                  call spidupdate
               endif
               if(kpr.ge.0) call wrd
          
            elseif(k_grid.eq.1) then
               call f_stepon( KLUCH, k_auto, nstep, dt, time, voltpf,
     &                        d_pf_mat, d_cam_mat, rax, zax, key_dmf)
               if(key_out.gt.0) then
                  call f_spidupdate
               endif
               call cur_avg
               if(kpr.ge.0) call f_wrd
            endif

            if(nstep/20*20.eq.nstep) then 
               call wrd_tim
            endif
         endif        !time steping

      else

!here is with breakdown

         call put_tim(dt,time)
         call put_key_fix(k_fixfree)

         if(nstep.eq.0) then
            KLUCH = 0
         else
            KLUCH = 1
         endif

         if (plasma_up.eq.1.and.plasma_trig.eq.1) then 
            KLUCH=0
      	    j_switch=0
         endif

         if (plasma_up.eq.1.and.plasma_trig.eq.0) KLUCH=1
         if (s_adapt.eq.1.and.j_switch.eq.0) then
            j_switch=1
         endif
      
         if(KLUCH.eq.0) then    !initialization

            if (plasma_up.eq.0) then	
               call fbkdw_stepon( KLUCH, k_auto, nstep, dt, time,
     &                 voltpf, d_pf_mat, d_cam_mat, rax, zax, key_dmf)
               return
            endif

            if (yesfitcc.eq.0) then
               call sstepon_bkdw( KLUCH, k_auto, nstep, dt, time,
     &                            voltpf, d_pf_mat, d_cam_mat, key_dmf)
            else
               if (yesfitcc.eq.1) call cf_init( k_auto, nstep, dt, time,
     &                     voltpf, d_pf_mat,d_cam_mat ) !fit only 12 currents up to psl
      	       if (yesfitcc.eq.2) call cf_init_bkwd( k_auto, nstep, dt,
     &                     time, voltpf, d_pf_mat,d_cam_mat ) !now it works
      	       if (yesfitcc.eq.3) call cf_init_full( k_auto, nstep, dt,
     &                     time, voltpf, d_pf_mat,d_cam_mat ) !not yet implemented 
            endif

            call  wrrec

         elseif(KLUCH.eq.1) then   !time steping

            if (plasma_up.eq.0) then
               call fbkdw_stepon( KLUCH, k_auto, nstep, dt, time,
     &                  voltpf, d_pf_mat, d_cam_mat, rax, zax, key_dmf)
               return
            endif

            if(k_grid.eq.0) then
               if (s_fazt.eq.2) kpr=-1
               if (s_fazt.eq.3) kpr=-2
               call sstepon_bkdw( KLUCH, k_auto, nstep, dt, time,
     &                  voltpf, d_pf_mat, d_cam_mat, rax, zax, key_dmf)
               if(kpr.ge.0) call wrd
            else
               if (j_switch.eq.2) then
                  write(*,*) 'adaptige grid evol', KLUCH, k_auto,
     &               nstep_local, dt, time, s_fazt
                  if (s_fazt.eq.0) kpr=0
                  if (s_fazt.eq.1) kpr=-2
                  call f_stepon_bkdw( 1, 0, 1+nstep_local, dt, time,
     &                  voltpf, d_pf_mat, d_cam_mat, rax, zax, key_dmf)
      	       else
                  write(*,*) 'reset rectangular grid'
                  kpr=0
                  call sstepon_bkdw000( 0, 0, 0, dt, time, voltpf,
     &                  d_pf_mat, d_cam_mat, rax, zax, key_dmf)
                  call wrrec
                  write(*,*) 'switch to adaptive'     
                  call f_stepon_bkdw( 0, 0, 0, dt, time, voltpf,
     &                  d_pf_mat, d_cam_mat, rax, zax, key_dmf)
                  write(*,*) 'adaptive evol'     
                  call f_stepon_bkdw( 1, 0, 1, dt, time, voltpf,
     &                  d_pf_mat, d_cam_mat, rax, zax, key_dmf)
                  j_switch=2
               endif

               if(kpr.ge.0) call f_wrd
               call cur_avg
            endif     

            nstep_local=nstep_local+1
            if(nstep/20*20.eq.nstep) then 
               call wrd_tim
            endif

         endif        !time steping

      endif

      return
      end subroutine spider

!----------------------------------------------------------------
      subroutine get_psib(r_ax, z_ax, psi_b)

      use comblc, only: rm, zm, up

      implicit none

      real*8, intent(out) :: r_ax, z_ax, psi_b

      r_ax  = rm
      z_ax  = zm
      psi_b = up

      return
      end subroutine get_psib

!----------------------------------------------------------------
      subroutine kpr_calc(kpr_xx)

      use keys, only: kpr

      implicit none

      integer, intent(in) :: kpr_xx

      kpr=kpr_xx

      return
      end subroutine kpr_calc

!----------------------------------------------------------------
      subroutine put_name(name)

      use iopath, only: path

      implicit none

      character(len=80), intent(in) :: name
     
      path = name

      return
      end subroutine put_name

!----------------------------------------------------------------
      subroutine put_Ipl(placur)

      implicit none

      real*8, intent(in) :: placur
      real*8 :: cur_pl
      common /com_curpl/ cur_pl         

      cur_pl = placur

      return
      end subroutine put_Ipl

!----------------------------------------------------------------
      subroutine get_Ipl(placur)

      implicit none

      real*8, intent(out) :: placur
      real*8 :: cur_pl

      common /com_curpl/ cur_pl         

      placur=cur_pl

      return
      end subroutine get_Ipl

!----------------------------------------------------------------
      subroutine get_zccur(rc_cur, zc_cur, z2c_cur)

      use comblc, only: ni, nj, ipr, curf, r, z, dri, dzj

      implicit none

      real*8, intent(out) :: rc_cur, zc_cur, z2c_cur

      integer :: i, j
      real*8 :: sum_z2c, sum_zc, sum_rc, sum_Ipla

      sum_z2c=0.d0
      sum_zc=0.d0
      sum_rc=0.d0
      sum_Ipla=0.d0
         
      do i=2,ni-1
         do j=2,nj-1
            if(ipr(i,j) .eq. 1) then
               sum_Ipla=sum_Ipla+dri(i)*dzj(j)*curf(i,j)
               sum_z2c = sum_z2c+dri(i)*dzj(j)*curf(i,j)*z(j)**2
               sum_zc = sum_zc+dri(i)*dzj(j)*curf(i,j)*z(j)
               sum_rc = sum_rc+dri(i)*dzj(j)*curf(i,j)*r(i)
            endif
         enddo
      enddo

      z2c_cur=sum_z2c/sum_Ipla
      zc_cur=sum_zc/sum_Ipla
      rc_cur=sum_rc/sum_Ipla

      return
      end subroutine get_zccur

!----------------------------------------------------------------
      subroutine get_zccurb(rc_cur,zc_cur,z2c_cur,
     & rgeoc,zgeoc,ahorc) 

      use sp_parameters

      implicit none

      real*8, intent(out) :: rc_cur, zc_cur, z2c_cur,
     &   rgeoc, zgeoc, ahorc
      include 'compol.inc'

      integer :: i, j
      real*8 :: perimz, ahorc2, sum_rc, sum_zc, sum_z2c, sum_ipla
      
      sum_z2c=0.d0
      sum_zc=0.d0
      sum_rc=0.d0
      sum_Ipla=0.d0
      		 
      perimz=0.         
      rgeoc=0.
      zgeoc=0.
      ahorc=0.
      ahorc2=0.
      	
      i=iplas

      do j=2,nt1
         perimz=perimz+dlt(i,j)
         rgeoc=rgeoc+r(i,j)*dlt(i,j)
         zgeoc=zgeoc+z(i,j)*dlt(i,j)
      enddo

      rgeoc=rgeoc/perimz
      zgeoc=zgeoc/perimz
      
      do j=2,nt1
         ahorc2=ahorc2+(r(i,j)-rgeoc)**2.*dlt(i,j)
      enddo
      ahorc2=ahorc2/perimz
      ahorc=2.*sqrt(ahorc2)

      zc_cur=sum_zc/sum_Ipla
      z2c_cur=zc_cur
      rc_cur=sqrt(sum_rc/sum_Ipla)

      return
      end subroutine get_zccurb

!----------------------------------------------------------------
      subroutine get_zccur_ada(rc_cur, zc_cur, z2c_cur)

      use sp_parameters

      implicit none

      real*8, intent(out) :: rc_cur, zc_cur, z2c_cur

      integer :: i, j
      real*8 :: sum_z2c, sum_zc, sum_rc, sum_Ipla, sqcen, r0, z0, sqk

      include 'compol.inc'

      sum_z2c=0.d0
      sum_zc=0.d0
      sum_rc=0.d0
      sum_Ipla=0.d0         
      sqcen=0.d0

      do j=2,nt1
         sqcen=sqcen+sq1(1,j)+sq4(1,j)
      enddo
   
      sum_Ipla=sum_Ipla+sqcen*cur(1,2)
      sum_z2c = sum_z2c+sqcen*cur(1,2)*zm**2
      sum_zc = sum_zc+sqcen*cur(1,2)*zm
      sum_rc = sum_rc+sqcen*cur(1,2)*rm
         
      do i=2,iplas
         do j=2,nt1
            if(i.ne.iplas) then
               sqk=sq1(i,j)+sq2(i-1,j)+sq3(i-1,j-1)+sq4(i,j-1)
               r0=r(i,j)
               z0=z(i,j)
            else
               sqk=sq2(i-1,j)+sq3(i-1,j-1)
               r0=r(i,j)
               z0=z(i,j)
            endif
            sum_Ipla=sum_Ipla+sqk*cur(i,j)
            sum_z2c = sum_z2c+sqk*cur(i,j)*z0**2
            sum_zc = sum_zc+sqk*cur(i,j)*z0
            sum_rc = sum_rc+sqk*cur(i,j)*r0
         enddo
      enddo
         
      z2c_cur=sum_z2c/sum_Ipla
      zc_cur=sum_zc/sum_Ipla
      rc_cur=sum_rc/sum_Ipla

      return
      end subroutine get_zccur_ada

!----------------------------------------------------------------
      subroutine find_demo_gaps(ngaps, demo_gaps, geom1d)

      use comblc, only: r, z, u, up

      implicit none

      integer, intent(in) :: ngaps
      double precision, intent(in), dimension(ngaps, 3) :: demo_gaps
      double precision, intent(out), dimension(ngaps) :: geom1d

      integer :: i, j, n_iterz, j1, j2, j3, j4
      double precision :: dumx, dumy, dumz, l_ref, l_gap, dumx2, dumy2,
     &   dumu1, dumu2, dumx3, dumy3, l_gap2, l_gap3, dumx0, dumy0, 
     &   dumxx, dumyy, d_step, gapmin, gapmax, x00, x002, x003, 
     &   y00, y002, y003, u00, u002
      real*8 :: blin_

      d_step=0.01 !advance in 1 cm steps	
      n_iterz=1000
      gapmin=-1.
      gapmax=1.5

      do i=1,ngaps
! first positive gap
         l_ref=d_step

         dumx0=demo_gaps(i,1) ! R0
         dumy0=demo_gaps(i,2) ! Z0
       	 dumz=demo_gaps(i,3)  ! angle	
      	 dumx=dumx0           ! R
      	 dumy=dumy0           ! Z
   	
         do j=1,n_iterz
            j1=((dumx-r(1))/(r(2)-r(1)))+1
            j2=((dumy-z(1))/(z(2)-z(1)))+1
            x00=r(j1)
            y00=z(j2)
            x002=r(j1+1)
            y002=z(j2+1)
            u00 = blin_(dumx, dumy, x00, x002, y00, y002, u(j1,j2), 
     &                u(j1+1,j2), u(j1+1,j2+1), u(j1,j2+1))

            dumx2=dumx+l_ref*cos(dumz)
            dumy2=dumy+l_ref*sin(dumz)

            j3=((dumx2-r(1))/(r(2)-r(1)))+1	
            j4=((dumy2-z(1))/(z(2)-z(1)))+1	

      	    x002=r(j3)
      	    y002=z(j4)
      	    x003=r(j3+1)
      	    y003=z(j4+1)
      	    u002 = blin_(dumx2, dumy2, x002, x003, y002, y003, 
     &          u(j3,j4), u(j3+1,j4), u(j3+1,j4+1), u(j3,j4+1))

            if (u002.eq.up) goto 117
            if (u00.eq.up) goto 118
            if (u002.gt.up.and.u00.lt.up) goto 116
            if (u002.lt.up.and.u00.gt.up) goto 116

!case with no intersection: larger 1. m
            if (l_ref.ge.gapmax.or.l_ref.lt.gapmin) goto 119

            dumx=dumx2
            dumy=dumy2
         enddo

116      continue
         dumxx=(dumx2*abs(u002-up)+dumx*abs(u00-up))/abs(u002-u00)
         dumyy=(dumy2*abs(u002-up)+dumy*abs(u00-up))/abs(u002-u00)
         l_gap=sqrt((dumxx-dumx0)**2.+(dumyy-dumy0)**2.)
         goto 211

117      continue
         l_gap=sqrt((dumx2-dumx0)**2.+(dumy2-dumy0)**2.)
         goto 211

119      continue !case there is no intersection
         l_gap=-50.
         goto 211

118      continue
         l_gap=sqrt((dumx-dumx0)**2.+(dumy-dumy0)**2.)

211      continue

!now negative gap
         l_ref=-d_step 

         dumx0=demo_gaps(i,1) ! R0
         dumy0=demo_gaps(i,2) ! Z0
         dumz=demo_gaps(i,3)  ! angle	
         dumx=dumx0           ! R
         dumy=dumy0           ! Z
      	
         do j=1,n_iterz

      	    j1=((dumx-r(1))/(r(2)-r(1)))+1	
      	    j2=((dumy-z(1))/(z(2)-z(1)))+1	

            x00=r(j1)
            y00=z(j2)
            x002=r(j1+1)
            y002=z(j2+1)
            u00 = blin_(dumx, dumy, x00, x002, y00, y002, u(j1,j2),
     &            u(j1+1,j2), u(j1+1,j2+1), u(j1,j2+1))

            dumx2=dumx+l_ref*cos(dumz)
            dumy2=dumy+l_ref*sin(dumz)

            j3=((dumx2-r(1))/(r(2)-r(1)))+1	
            j4=((dumy2-z(1))/(z(2)-z(1)))+1	

      	    x002=r(j3)
      	    y002=z(j4)
      	    x003=r(j3+1)
      	    y003=z(j4+1)
      	    u002 = blin_(dumx2, dumy2, x002, x003, y002, y003, 
     &             u(j3,j4), u(j3+1,j4), u(j3+1,j4+1), u(j3,j4+1))

            if (u002.eq.up) goto 317
            if (u00.eq.up) goto 318
            if (u002.gt.up.and.u00.lt.up) goto 316
            if (u002.lt.up.and.u00.gt.up) goto 316

!case with no intersection: larger 1. m
            if (l_ref.ge.gapmax.or.l_ref.lt.gapmin) goto 319

            dumx=dumx2
            dumy=dumy2
         enddo

316	 continue
         dumxx=(dumx2*abs(u002-up)+dumx*abs(u00-up))/abs(u002-u00)
         dumyy=(dumy2*abs(u002-up)+dumy*abs(u00-up))/abs(u002-u00)
         l_gap2=-sqrt((dumxx-dumx0)**2.+(dumyy-dumy0)**2.)
         goto 411

317      continue
         l_gap2=-sqrt((dumx2-dumx0)**2.+(dumy2-dumy0)**2.)
         goto 411

319      continue !case there is no intersection
         l_gap2=50.
         goto 411
      
318      continue
         l_gap2=-sqrt((dumx-dumx0)**2.+(dumy-dumy0)**2.)

411      continue

!choose minimum of absolute values
         l_gap3=l_gap
         if (abs(l_gap2).lt.abs(l_gap)) l_gap3=l_gap2

         geom1d(i)=min(gapmax,max(gapmin,l_gap3))

      enddo

      return
      end subroutine find_demo_gaps
