       subroutine spider(nstep,time,dt,key_dmf,k_grid,k_auto,k_fixfree,
     &                   dpsdt,key_start,voltpf)  

!---------------------------------------------------------------
!   MAIN  PROGRAM  OF  THE EVOLUTION CODE  "PET"
!---------------------------------------------------------------

      use parameters, only: s_adapt, s_fazt

      implicit none
	  common /com_0st/ key_0st,key_prs
	  common /c_kpr/ kpr
	  common /com_kout/ key_out
			common /plasma_state/ plasma_up,plasma_trig
	common /yes_bkdw/ ibkdw,ifbey
			integer plasma_up,plasma_trig
	double precision ibkdw,ifbey
	double precision time,dt,voltpf
	integer key_prs
	integer j_switch,k_grid,kpr
	common /yes_fit_cc/ yesfitcc

	real*8 rax,zax
	integer k_fixfree,nstep,kluch,k_auto
	integer yesfitcc,k_dmf

	integer nstep_local,key_out,key_start
	data j_switch /0/
	real*8 contvals_mat,d_pf_mat,d_cam_mat
        dimension  contvals_mat(2500),d_pf_mat(500)
        dimension  voltpf(*)
        dimension  d_cam_mat(500)
	integer key_0st,nnstep,key_dmf,kkey_dmf
	real*8 dpsdt
	save j_switch
	save nstep_local
      s_fazt = 0
	if (ibkdw*ifbey.eq.0) then			 

              call put_tim(dt,time)
              call put_key_fix(k_fixfree)

       if(nstep.eq.0) then
          KLUCH = 0
       else
          KLUCH = 1
       endif

         if(k_fixfree.eq.0) then    
          call B_STEPON( KLUCH, k_auto, nstep, dt, time,
     *                   rax,zax ,key_dmf,dpsdt)
   
        if(key_0st.eq.1 .AnD. nstep.eq.0) then
           nnstep=1        
           kkey_dmf=-10        
           KLUCH = 1
        
          call B_STEPON( KLUCH, k_auto, nnstep, dt, time,
     *                   rax,zax ,kkey_dmf,dpsdt)
        
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

          call  sstepon( KLUCH, k_auto,nstep,dt,time,
     *                  voltpf, d_pf_mat ,d_cam_mat,key_dmf)
     
          if(key_out.gt.0) then
            call spidupdate
          endif

       else
	write(*,*) 'fix grid call init'
          call  cf_init( k_auto, nstep, dt, time,
     *                     voltpf, d_pf_mat,d_cam_mat )
     
          if(key_out.gt.0) then
            call spidupdate
          endif

       endif

          call  wrrec

	write(*,*) 'after stepon initial'

       k_auto= 0  ! don't change!

        if(k_grid.eq.1) then
!	write(*,*) 'adaptive grid call'
          call  f_stepon( KLUCH, k_auto,nstep,dt,time,
     *                  voltpf, d_pf_mat ,d_cam_mat,rax,zax,key_dmf)
     
          if(key_out.gt.0) then
            call f_spidupdate
          endif

        endif

        elseif(KLUCH.eq.1) then   !time steping

        if(k_grid.eq.0) then
!	write(*,*) 'fix grid call evol'
	if (s_fazt.eq.2) kpr=-1
	if (s_fazt.eq.3) kpr=-2
          call  sstepon( KLUCH, k_auto,nstep,dt,time,
     *                  voltpf, d_pf_mat ,d_cam_mat,key_dmf)
     
          if(key_out.gt.0) then
            call spidupdate
          endif

          if(kpr.ge.0) call wrd
          
        elseif(k_grid.eq.1) then

!	write(*,*) 'adaptive grid call evol'
          call  f_stepon( KLUCH, k_auto,nstep,dt,time,
     *                  voltpf, d_pf_mat ,d_cam_mat,rax,zax,key_dmf)
     
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
!	write(*,*) 'fbkdw'
			
	          call  fbkdw_stepon( KLUCH, k_auto,nstep,dt,time,
     *                  voltpf, d_pf_mat ,d_cam_mat,rax,zax,key_dmf)

	return
	endif

	if (yesfitcc.eq.0) then
          call  sstepon_bkdw( KLUCH, k_auto,nstep,dt,time,
     *                  voltpf, d_pf_mat ,d_cam_mat,key_dmf)
	else
		if (yesfitcc.eq.1)	call  cf_init( k_auto, nstep, dt, time,
     *                     voltpf, d_pf_mat,d_cam_mat ) !fit only 12 currents up to psl
		if (yesfitcc.eq.2)	call  cf_init_bkwd( k_auto, nstep, dt, time,
     *                     voltpf, d_pf_mat,d_cam_mat ) !now it works
		if (yesfitcc.eq.3)	call  cf_init_full( k_auto, nstep, dt, time,
     *                     voltpf, d_pf_mat,d_cam_mat ) !not yet implemented 

	endif

          call  wrrec

        elseif(KLUCH.eq.1) then   !time steping

	if (plasma_up.eq.0) then
			
!       k_auto= 0  ! don't change!
	          call  fbkdw_stepon( KLUCH, k_auto,nstep,dt,time,
     *                  voltpf, d_pf_mat ,d_cam_mat,rax,zax,key_dmf)

	return
	endif

!      k_grid=1
	if(k_grid.eq.0) then

! write(*,*) 'call fixed grid evol'  
	if (s_fazt.eq.2) kpr=-1
	if (s_fazt.eq.3) kpr=-2
	        call  sstepon_bkdw( KLUCH, k_auto,nstep,dt,time,
     *                  voltpf, d_pf_mat ,d_cam_mat,rax,zax,key_dmf)
          if(kpr.ge.0) call wrd

				else

	if (j_switch.eq.2) then
	write(*,*) 'adaptige grid evol', KLUCH, k_auto,nstep_local,
     & dt,time,s_fazt
	if (s_fazt.eq.0) kpr=0
	if (s_fazt.eq.1) kpr=-2
          call  f_stepon_bkdw( 1, 0,
     & 1+nstep_local,dt,time,
     *                  voltpf, d_pf_mat ,d_cam_mat,rax,zax,key_dmf)
		else
	write(*,*) 'reset rectangular grid'     

					kpr=0
          call  sstepon_bkdw000( 0, 0,0,dt,time,
     *                  voltpf, d_pf_mat ,d_cam_mat,rax,zax,key_dmf)
					call wrrec
				
					write(*,*) 'switch to adaptive'     
          call  f_stepon_bkdw( 0, 0,0,dt,time,
     *                  voltpf, d_pf_mat ,d_cam_mat,rax,zax,key_dmf)
					write(*,*) 'adaptive evol'     
          call  f_stepon_bkdw( 1, 0,1,dt,time,
     *                  voltpf, d_pf_mat ,d_cam_mat,rax,zax,key_dmf)

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
       end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_psix(r_xp,z_xp,psi_xp)

        include 'double.inc'
        include 'param.inc'
        include 'comblc.inc'

         r_xp=rx0
         z_xp=zx0
         psi_xp=ux0

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_fpsix(r_xp,z_xp,psi_xp)

        include 'double.inc'
        include 'parevo.inc'
        parameter(nkp=njlim)
         include 'dim.inc'
         include 'compol.inc'
         include 'compol_add.inc'

         r_xp=rx0
         z_xp=zx0
         psi_xp=psix0

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_psib(r_ax,z_ax,psi_b)

        include 'double.inc'
        include 'param.inc'
        include 'comblc.inc'

         r_ax=rm
         z_ax=zm
         psi_b=up

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine  kpr_calc(kpr_xx)

      !include 'double.inc'
	common
     *  /c_kpr/kpr
	
	kpr=kpr_xx

	!print *,' kpr  FOR DEBUGING',kpr
	return
	end

	subroutine pau()
	return
	end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         subroutine put_name(name,ksym)
          include 'iopath.inc'
          character*40 name
          integer ksym
            path=name
            kname=ksym
         return
         end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine put_Ipl(placur)

        include 'double.inc'
        common /com_curpl/ cur_pl         

         cur_pl=placur

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_Ipl(placur)

        include 'double.inc'
        common /com_curpl/ cur_pl         

         placur=cur_pl

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_zccur(rc_cur,zc_cur,z2c_cur)

        include 'double.inc'
        include 'param.inc'
        include 'comblc.inc'

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
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_zccurb(rc_cur,zc_cur,z2c_cur,
     & rgeoc,zgeoc,ahorc,bpcell)

        include 'double.inc'
        include 'dim.inc'
        include 'compol.inc'
	real*8 perimz,ahorc2,ahorc,rgeoc,zgeoc
	real*8 bpcell(nt1-1)
	
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
	bpcell(nt1-1)=bpcell(nt1-2)

          do j=2,nt1
          
          
	perimz=perimz+dlt(i,j)
	sum_rc=sum_rc+(r(i,j)**2.)*bpcell(j-1)*dlt(i,j)
	sum_zc=sum_zc+z(i,j)*bpcell(j-1)*dlt(i,j)
	sum_Ipla=sum_Ipla+bpcell(j-1)*dlt(i,j)

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
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine find_demo_gaps(ngaps,demo_gaps,geom1d)

        include 'double.inc'
        include 'param.inc'
        include 'comblc.inc'

	integer ngaps,i,j
	integer n_iterz
	double precision demo_gaps(ngaps,3),geom1d(ngaps)
	double precision dumx,dumy,dumz,l_ref,l_gap
	double precision dumx2,dumy2,dumu1,dumu2
	double precision dumx3,dumy3,l_gap2,l_gap3
	double precision dumx0,dumy0,dumxx,dumyy,d_step
	double precision gapmin,gapmax

	d_step=0.01 !advance in 1 cm steps	
	n_iterz=1000
	gapmin=-1.
	gapmax=1.5

	do i=1,ngaps


!first positive gap!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	l_ref=d_step

			dumx0=demo_gaps(i,1) ! R0
			dumy0=demo_gaps(i,2) ! Z0
			dumz=demo_gaps(i,3) ! angle	 
			
			dumx=dumx0 ! R
			dumy=dumy0 ! Z
		
	do j=1,n_iterz

		j1=((dumx-r(1))/(r(2)-r(1)))+1	
		j2=((dumy-z(1))/(z(2)-z(1)))+1	

		x00=r(j1)
		y00=z(j2)
		x002=r(j1+1)
		y002=z(j2+1)
		u00=blin_(dumx,dumy,x00,x002,y00,y002,u(j1,j2)
     & ,u(j1+1,j2),u(j1+1,j2+1),u(j1,j2+1))

		dumx2=dumx+l_ref*cos(dumz)
		dumy2=dumy+l_ref*sin(dumz)

		j3=((dumx2-r(1))/(r(2)-r(1)))+1	
		j4=((dumy2-z(1))/(z(2)-z(1)))+1	

		x002=r(j3)
		y002=z(j4)
		x003=r(j3+1)
		y003=z(j4+1)
		u002=blin_(dumx2,dumy2,x002,x003,y002,y003,u(j3,j4)
     & ,u(j3+1,j4),u(j3+1,j4+1),u(j3,j4+1))
		
	if (u002.eq.up) goto 117
	if (u00.eq.up) goto 118
	if (u002.gt.up.and.u00.lt.up) goto 116
	if (u002.lt.up.and.u00.gt.up) goto 116

	
!case with no intersection: larger 1. m
	if (l_ref.ge.gapmax.or.l_ref.lt.gapmin)	goto 119

	dumx=dumx2
	dumy=dumy2
	goto 108	

105	continue

108	continue
	
	enddo

116	continue
	dumxx=(dumx2*abs(u002-up)+dumx*abs(u00-up))/abs(u002-u00)
	dumyy=(dumy2*abs(u002-up)+dumy*abs(u00-up))/abs(u002-u00)
	l_gap=sqrt((dumxx-dumx0)**2.+(dumyy-dumy0)**2.)
	goto 211

117	continue
	l_gap=sqrt((dumx2-dumx0)**2.+(dumy2-dumy0)**2.)
	goto 211

119	continue !case there is no intersection
	l_gap=-50.
	goto 211
	
118	continue
	l_gap=sqrt((dumx-dumx0)**2.+(dumy-dumy0)**2.)

211	continue
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!now negative gap
	l_ref=-d_step 

			dumx0=demo_gaps(i,1) ! R0
			dumy0=demo_gaps(i,2) ! Z0
			dumz=demo_gaps(i,3) ! angle	 
			
			dumx=dumx0 ! R
			dumy=dumy0 ! Z
		
	do j=1,n_iterz

		j1=((dumx-r(1))/(r(2)-r(1)))+1	
		j2=((dumy-z(1))/(z(2)-z(1)))+1	

		x00=r(j1)
		y00=z(j2)
		x002=r(j1+1)
		y002=z(j2+1)
		u00=blin_(dumx,dumy,x00,x002,y00,y002,u(j1,j2)
     & ,u(j1+1,j2),u(j1+1,j2+1),u(j1,j2+1))

		dumx2=dumx+l_ref*cos(dumz)
		dumy2=dumy+l_ref*sin(dumz)

		j3=((dumx2-r(1))/(r(2)-r(1)))+1	
		j4=((dumy2-z(1))/(z(2)-z(1)))+1	

		x002=r(j3)
		y002=z(j4)
		x003=r(j3+1)
		y003=z(j4+1)
		u002=blin_(dumx2,dumy2,x002,x003,y002,y003,u(j3,j4)
     & ,u(j3+1,j4),u(j3+1,j4+1),u(j3,j4+1))
		
	if (u002.eq.up) goto 317
	if (u00.eq.up) goto 318
	if (u002.gt.up.and.u00.lt.up) goto 316
	if (u002.lt.up.and.u00.gt.up) goto 316

	
!case with no intersection: larger 1. m
	if (l_ref.ge.gapmax.or.l_ref.lt.gapmin)	goto 319

	dumx=dumx2
	dumy=dumy2
	goto 308	

305	continue
!	l_ref=min(0.5,l_ref*1.1)

308	continue
	
	enddo

316	continue
	dumxx=(dumx2*abs(u002-up)+dumx*abs(u00-up))/abs(u002-u00)
	dumyy=(dumy2*abs(u002-up)+dumy*abs(u00-up))/abs(u002-u00)
	l_gap2=-sqrt((dumxx-dumx0)**2.+(dumyy-dumy0)**2.)
	goto 411

317	continue
	l_gap2=-sqrt((dumx2-dumx0)**2.+(dumy2-dumy0)**2.)
	goto 411

319	continue !case there is no intersection
	l_gap2=50.
	goto 411
	
318	continue
	l_gap2=-sqrt((dumx-dumx0)**2.+(dumy-dumy0)**2.)

411	continue
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!choose minimum of absolute values
	l_gap3=l_gap
	if (abs(l_gap2).lt.abs(l_gap)) l_gap3=l_gap2

	geom1d(i)=min(gapmax,max(gapmin,l_gap3))

	enddo
	
	end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine find_demo_gaps2(ngaps,demo_gaps,geom1d)

        include 'double.inc'
        include 'param.inc'
        include 'comblc.inc'

	integer ngaps,i,j
	integer n_iterz
	double precision demo_gaps(ngaps,3),geom1d(ngaps)
	double precision dumx,dumy,dumz,l_ref,l_gap
	double precision dumx2,dumy2,dumu1,dumu2
	double precision dumx3,dumy3

	do i=1,ngaps

	l_ref=0.1 !advance in 10 cm steps	
	n_iterz=10

			dumx=demo_gaps(i,1) ! R
			dumy=demo_gaps(i,2) ! Z
			dumz=demo_gaps(i,3) ! angle	 
			
		j1=((dumx-r(1))/(r(2)-r(1)))+1	
		j2=((dumy-z(1))/(z(2)-z(1)))+1	

		x00=r(j1)
		y00=z(j2)
		x002=r(j1+1)
		y002=z(j2+1)
		u00=blin_(dumx,dumy,x00,x002,y00,y002,u(j1,j2)
     & ,u(j1+1,j2),u(j1+1,j2+1),u(j1,j2+1))
		
	do j=1,n_iterz

		dumx2=dumx+l_ref*cos(dumz)
		dumy2=dumy+l_ref*sin(dumz)

		j3=((dumx2-r(1))/(r(2)-r(1)))	
		j4=((dumy2-z(1))/(z(2)-z(1)))	

		x002=r(j3)
		y002=z(j4)
		u002=u(j3,j4)
		
		dpsidl=(u002-u00)/l_ref
	if (dpsidl.eq.0.) then
		l_gap=l_ref*1.4
	else
		l_gap=(up-u00)/dpsidl
	endif
		l_ref=0.5*(l_ref+l_gap)

	enddo
	geom1d(i)=l_gap

	write(999,*) geom1d(i)
	enddo

	end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine find_demo_gaps3(ngaps,demo_gaps,geom1d)

        include 'double.inc'
        include 'param.inc'
        include 'comblc.inc'

	integer ngaps,i,j
	integer n_iterz
	double precision demo_gaps(ngaps,3),geom1d(ngaps)
	double precision dumx,dumy,dumz,l_ref,l_gap
	double precision dumx2,dumy2,dumu1,dumu2
	double precision dumx3,dumy3


	do i=1,ngaps

	l_ref=0.1 !advance in 5 cm steps	
	n_iterz=10

			dumx=demo_gaps(i,1) ! R
			dumy=demo_gaps(i,2) ! Z
			dumz=demo_gaps(i,3) ! angle	 
			
		j1=((dumx-r(1))/(r(2)-r(1)))	
		j2=((dumy-z(1))/(z(2)-z(1)))	

		x00=r(j1)
		y00=z(j2)
		u00=u(j1,j2)
		
	do j=1,n_iterz

		dumx2=dumx+l_ref*cos(dumz)
		dumy2=dumy+l_ref*sin(dumz)

		j3=((dumx2-r(1))/(r(2)-r(1)))	
		j4=((dumy2-z(1))/(z(2)-z(1)))	

		x002=r(j3)
		y002=z(j4)
		u002=u(j3,j4)
		
		dpsidl=(u002-u00)/l_ref
	if (dpsidl.eq.0.) then
		l_gap=l_ref*1.4
	else
		l_gap=(up-u00)/dpsidl
	endif
		l_ref=0.5*(l_ref+l_gap)

	enddo
	geom1d(i)=l_gap

	write(999,*) geom1d(i)
	enddo

	end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine find_demo_meas(nmea1,nmea2,
     & bri1,bro1,bzi1,bzo1,fli1,flo1,u_cd)

        include 'double.inc'
        include 'param.inc'
        include 'comblc.inc'

	integer nmea1,nmea2,i,j,j1,j2
	double precision bri1(nmea1,3)
	double precision bro1(nmea2,3)
	double precision bzi1(nmea1,3)
	double precision bzo1(nmea2,3)
	double precision fli1(nmea1,2)
	double precision flo1(nmea2,2)
	double precision u_cd(3*nmea1+3*nmea2)
	double precision dumx,dumy,dumz,l_ref,l_gap
	double precision dumx2,dumy2
	double precision dumx3,dumy3

	do i=1,nmea1

			dumx=bri1(i,1) ! R
			dumy=bri1(i,2) ! Z
			dumz=bri1(i,3) ! angle	 
			
		j1=((dumx-r(1))/(r(2)-r(1)))	
		j2=((dumy-z(1))/(z(2)-z(1)))	

		x00=r(j1)
		y00=z(j2)
		u_cd(i)=-1/r(j1)*(u(j1,j2+1)-u(j1,j2))/
     & (z(2)-z(1))*cos(dumz)
	enddo

	do i=1,nmea2

			dumx=bro1(i,1) ! R
			dumy=bro1(i,2) ! Z
			dumz=bro1(i,3) ! angle	 
			
		j1=((dumx-r(1))/(r(2)-r(1)))	
		j2=((dumy-z(1))/(z(2)-z(1)))	

		x00=r(j1)
		y00=z(j2)
		u_cd(60+i)=-1/r(j1)*(u(j1,j2+1)-u(j1,j2))/
     & (z(2)-z(1))*cos(dumz)
	enddo

	do i=1,nmea1

			dumx=bzi1(i,1) ! R
			dumy=bzi1(i,2) ! Z
			dumz=bzi1(i,3) ! angle	 
			
		j1=((dumx-r(1))/(r(2)-r(1)))	
		j2=((dumy-z(1))/(z(2)-z(1)))	

		x00=r(j1)
		y00=z(j2)
		u_cd(60+68+i)=1/r(j1)*(u(j1+1,j2)-u(j1,j2))/
     & (r(2)-r(1))*sin(dumz)
	enddo

	do i=1,nmea2

			dumx=bzo1(i,1) ! R
			dumy=bzo1(i,2) ! Z
			dumz=bzo1(i,3) ! angle	 
			
		j1=((dumx-r(1))/(r(2)-r(1)))	
		j2=((dumy-z(1))/(z(2)-z(1)))	

		x00=r(j1)
		y00=z(j2)
		u_cd(60+68+60+i)=1/r(j1)*(u(j1+1,j2)-u(j1,j2))/
     & (r(2)-r(1))*sin(dumz)
	enddo

	do i=1,nmea1

			dumx=fli1(i,1) ! R
			dumy=fli1(i,2) ! Z
			
		j1=((dumx-r(1))/(r(2)-r(1)))	
		j2=((dumy-z(1))/(z(2)-z(1)))	

		x00=r(j1)
		y00=z(j2)
		u_cd(60+68+60+68+i)=u(j1,j2)
	enddo

	do i=1,nmea2

			dumx=flo1(i,1) ! R
			dumy=flo1(i,2) ! Z
			
		j1=((dumx-r(1))/(r(2)-r(1)))	
		j2=((dumy-z(1))/(z(2)-z(1)))	

		x00=r(j1)
		y00=z(j2)
		u_cd(60+68+60+68+60+i)=u(j1,j2)
	enddo

	end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_zccur_ada(rc_cur,zc_cur,z2c_cur)

        include 'double.inc'
        include 'dim.inc'
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
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_zbran(r_bra,z_bra)

        include 'double.inc'
        include 'param.inc'
        include 'comblc.inc'

!default min Z

        z_bra = z(1)
        
!definition of cell index i

         do i=1,ni-1
         
          if(r_bra .gt. r(i) .AnD. r_bra .le. r(i+1)) then
           i_cell=i
           exit
          endif
         
         enddo
         
!definition of cell index j
         
         do j=1,nj-1
         
           u1=u(i_cell,j)
           u2=u(i_cell+1,j)
           u3=u(i_cell+1,j+1)
           u4=u(i_cell,j+1)
         
           r1=r(i_cell)
           r2=r(i_cell+1)
         
           u12=( u1*(r2-r_bra)+u2*(r_bra-r1) )/(r2-r1)
           u34=( u4*(r2-r_bra)+u3*(r_bra-r1) )/(r2-r1)
           
            psi_min=dmin1(u12,u34)
            psi_max=dmax1(u12,u34)
         
           if(ux0 .gt. psi_min .AnD. ux0 .le. psi_max) then     
              j_cell=j
              exit
           endif
         
          enddo
         
         
         z1=z(j_cell)
         z3=z(j_cell+1)
         
         z_bra = (ux0*(z3-z1) + u34*z1 - u12*z3)/(u34-u12)

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine f_get_zbran(r_bra,z_bra)

        include 'double.inc'
        include 'param.inc'
        !include 'comblc.inc'
        include 'comrec.inc'
        include 'dimpl1.inc'
        include 'compol_add.inc'
        real*8 ar1(1),az1(1),ar2(1),az2(1),psitok(1)

!default min Z

        z_bra = y(1)
        
!definition of cell index i
!	i_cell=1
         do i=1,ni-1
         
          if(r_bra .gt. x(i) .AnD. r_bra .le. x(i+1)) then
           i_cell=i
           exit
          endif
         
         enddo
         
!definition of cell index j
         
         do j=1,nj-1
         
           u1=ue(i_cell,j)
           u2=ue(i_cell+1,j)
           u3=ue(i_cell+1,j+1)
           u4=ue(i_cell,j+1)
         
           ar1(1)=x(i_cell)
           ar2(1)=x(i_cell+1)
           az1(1)=y(j)
           az2(1)=y(j+1)
           
           call flux_p(psitok,ar1,az1,1)
           u1 = u1 + psitok(1)
           call flux_p(psitok,ar2,az1,1)
           u2 = u2 + psitok(1)
           call flux_p(psitok,ar2,az2,1)
           u3 = u3 + psitok(1)
           call flux_p(psitok,ar1,az2,1)
           u4 = u4 + psitok(1)
         
           r1=ar1(1)
           r2=ar2(1)
           z1=az1(1)
           z2=az2(1)
           
           u12=( u1*(r2-r_bra)+u2*(r_bra-r1) )/(r2-r1)
           u34=( u4*(r2-r_bra)+u3*(r_bra-r1) )/(r2-r1)
           
            psi_min=dmin1(u12,u34)
            psi_max=dmax1(u12,u34)
         
           if(psix0 .gt. psi_min .AnD. psix0 .le. psi_max) then     
              j_cell=j
              exit
           endif
         
          enddo
         
         
         z1=y(j_cell)
         z3=y(j_cell+1)
         
         z_bra = (psix0*(z3-z1) + u34*z1 - u12*z3)/(u34-u12)

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_rmax(r_max,r_min,r_geo,z_max,delta_up)

        include 'double.inc'
        include 'param.inc'
       ! include 'comblc.inc'
        common/comlop/ rxb(nbndp2), zxb(nbndp2), nxb
        real*8 xxxx1(3), yyyy1(3), pppp1(3)

         r_max=0.d0
         
         do i=1,nxb

           if(rxb(i) .gt. r_max) then
            r_max = rxb(i)
           endif

         enddo
         
         r_min=r_max
                  
          do i=1,nxb

           if(rxb(i) .lt. r_min) then
            r_min = rxb(i)
           endif
           
         enddo
        
            r_geo=0.5d0*(r_max+r_min)
            
            
            z_max=zxb(1)
            izmax = 1
            
           do i=1,nxb

           if(zxb(i) .gt. z_max) then
            z_max = zxb(i)
            izmax = i
           endif
           
         enddo
         
         !delta_up = 2*(r_geo - rxb( maxloc(zxb(1:nxb),1)))/(r_max-r_min)
         !second order z(r)
         !izmax = maxloc(zxb(1:nxb),1)
          izmax1 = izmax-1
          if(izmax1.le.0) then
            izmax1 = nxb-1
          endif
          !izmax2 = izmax +2
          izmax2 = izmax +1
          if(izmax2.gt.nxb) then
            izmax2 = 2
          endif

          xxxx1(1)=rxb(izmax1)
	    xxxx1(2)=rxb(izmax)
	    xxxx1(3)=rxb(izmax2)
	    yyyy1(1)=zxb(izmax1)
	    yyyy1(2)=zxb(izmax)
	    yyyy1(3)=zxb(izmax2)
	    call polyfitcc_spid(xxxx1,yyyy1,pppp1)
	    yrzmax=-pppp1(2)/(2.*pppp1(1))
	    !yzmax=pppp1(1)*yrzmax**2.0+pppp1(2)*yrzmax+pppp1(3)
         delta_up = 2*(r_geo - yrzmax)/(r_max-r_min)
	    
      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_rmax_ada(r_max,r_min,r_geo,z_max,delta_up)

        include 'double.inc'
        include 'dim.inc'
        include 'compol.inc'
        real*8 xxxx1(3), yyyy1(3), pppp1(3)

         r_max=0.d0
         
         i=iplas
         
         do j=2,nt1

           if(r(i,j) .gt. r_max) then
            r_max = r(i,j)
           endif

         enddo
         
         r_min=r_max
                  
         do j=2,nt1

           if(r(i,j) .lt. r_min) then
            r_min = r(i,j)
           endif
           
         enddo
        
            r_geo=0.5d0*(r_max+r_min)
            
            z_max=z(i,1)
            jzmax = 1
            
         do j=2,nt1

           if(z(i,j) .gt. z_max) then
            z_max = z(i,j)
            jzmax = j
           endif
           
         enddo

        xxxx1(1)=r(i,jzmax-1)
	    xxxx1(2)=r(i,jzmax)
	    xxxx1(3)=r(i,jzmax+1)
	    yyyy1(1)=z(i,jzmax-1)
	    yyyy1(2)=z(i,jzmax)
	    yyyy1(3)=z(i,jzmax+1)
	    call polyfitcc_spid(xxxx1,yyyy1,pppp1)
	    yrzmax=-pppp1(2)/(2.*pppp1(1))
	    !yzmax=pppp1(1)*yrzmax**2.0+pppp1(2)*yrzmax+pppp1(3)
         delta_up = 2*(r_geo - yrzmax)/(r_max-r_min)
	    
      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_pfccur(pfc,ncoil)

        include 'double.inc'
        include 'prm.inc'
       common/comst1/ PJK(NJLIM),PJKP1(NJLIM),PJKP(NJLIM),PJKD(NJLIM)
       real*8 pfc(ncoil)

         do i=1,ncoil
          pfc(i)=pjk(i)*1.d6
         enddo

      return
      end
