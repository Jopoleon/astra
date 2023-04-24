       subroutine spider(nstep,time,dt,key_dmf,k_grid,k_auto,k_fixfree,
     &                   dpsdt,key_start,voltpf)  

C---------------------------------------------------------------
C   MAIN  PROGRAM  OF  THE EVOLUTION CODE  "PET"
C---------------------------------------------------------------
	use fenix_params
	use plasma_state
			 implicit none
	  common /com_0st/ key_0st,key_prs
	  common /c_kpr/ kpr
	  common /com_kout/ key_out
!	include 'includes_fsim.dat'

	double precision ibkdw,ifbey
	double precision time,dt,voltpf
	integer key_prs
	integer j_switch,k_grid,kpr

	real*8 rax,zax
	integer k_fixfree,nstep,kluch,k_auto
	integer k_dmf

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

	ibkdw=zibkdw
	ifbey=zifbey
!	write(*,*) 'sa',isafazt
!        save i_enter
!        
!        i_enter=i_enter+1
!        
!      if(i_enter .eq. 1) then
!         do i=1,19
!	 voltpf(i)=0.d0
!         enddo
!      endif
        
c  kpr=1 for debugging, kpr=0 no printing

	    !kpr=1
          !call  kpr_calc(kpr)

       !tau_con=2.0d-2
       !k_con=(tau_con+1.d-8)/dt
       !k_con=99999

			 
	if (ibkdw*ifbey.eq.0) then			 

!	write(*,*) 'no breakdown',nstep
!below here is without breakdown

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              call put_tim(dt,time)
              call put_key_fix(k_fixfree)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       if(nstep.eq.0) then
          KLUCH = 0
       else
          KLUCH = 1
       endif
!	write(*,*) 'fdfd',nstep,kluch,k_fixfree

         if(k_fixfree.eq.0) then    
!	write(*,*) 'call fixfree',kpr
          call B_STEPON( KLUCH, k_auto, nstep, dt, time,
     *                   rax,zax ,key_dmf,dpsdt)
         ! call cur_avg
         ! call  wrb
!        write(*,*) '**'
!        write(*,*) 'nstep,time',nstep,time
!        write(*,*) '**'
        
        
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

c---------------------
       !k_auto= 0
       
c----------------
!!!             ! k_grid= 0   rect. grid
!!!             ! k_grid= 1   adap. grid


       if(KLUCH.eq.0) then    !initialization
         
!
!!!!! basic free bound rectan, equilibrium ( KLUCH=0 )
!
        !key_start=1

       if(key_start .le. 0) then
!	write(*,*) 'fix grid call'

          call  sstepon( KLUCH, k_auto,nstep,dt,time,
     *                  voltpf, d_pf_mat ,d_cam_mat,key_dmf)
     
          if(key_out.gt.0) then
!            call spidupdate
          endif

       else
!	write(*,*) 'fix grid call init'
          call  cf_init( k_auto, nstep, dt, time,
     *                     voltpf, d_pf_mat,d_cam_mat )
     
          if(key_out.gt.0) then
!            call spidupdate
          endif

       endif

          !call  wrfb
          call  wrrec

!	write(*,*) 'after stepon initial'
c       stop

       k_auto= 0  ! don't change!

        if(k_grid.eq.1) then
!	write(*,*) 'adaptive grid call'
!          call  f_stepon( KLUCH, k_auto,nstep,dt,time,
!     *                  voltpf, d_pf_mat ,d_cam_mat,rax,zax,key_dmf)
     
          if(key_out.gt.0) then
!            call f_spidupdate
          endif

        endif

	   !call cntrlr_iter(0,k_grid,voltpf)

        elseif(KLUCH.eq.1) then   !time steping



!	print *,' Next step: stepon',nstep

      !   nstep1=i_enter-1
      !if(nstep1/k_con*k_con .eq. nstep1) then
      !   do i=1,19
	! voltpf(i)=0.d0
      !   enddo
	! call cntrlr_iter(1,k_grid,voltpf)
      !endif

        if(k_grid.eq.0) then
	if (isafazt.eq.2) kpr=-1
	if (isafazt.eq.3) kpr=-2
!	write(*,*) 'fix grid call evol',isafazt,kpr

	if (kpr.ge.-1) then
!	write(*,*) 'fix grid call evol slow'
          call  sstepon( KLUCH, k_auto,nstep,dt,time,
     *                  voltpf, d_pf_mat ,d_cam_mat,key_dmf)
	endif
	if (kpr.eq.-2) then
!	write(*,*) 'fix grid call evol fast'
          call  sstepon_bkdw_555( KLUCH, k_auto,nstep,dt,time,
     *                  voltpf, d_pf_mat ,d_cam_mat,key_dmf)
	endif
     
          if(key_out.gt.0) then
!            call spidupdate
          endif

          if(kpr.ge.0) call wrd
          
        elseif(k_grid.eq.1) then

!	write(*,*) 'adaptive grid call evol'
!          call  f_stepon( KLUCH, k_auto,nstep,dt,time,
!     *                  voltpf, d_pf_mat ,d_cam_mat,rax,zax,key_dmf)
     
          if(key_out.gt.0) then
!            call f_spidupdate
          endif
         call cur_avg
         if(kpr.ge.0) call f_wrd
         
        endif

!	print *,' After stepon'


        if(nstep/20*20.eq.nstep) then 
          call wrd_tim
        endif

!        write(*,*) '**'
!        write(*,*) 'nstep,time',nstep,time
!        write(*,*) '**'

        endif        !time steping

















	else








!	write(*,*) 'yes breakdown'


	!here is with breakdown
	


	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              call put_tim(dt,time)
              call put_key_fix(k_fixfree)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
     *   voltpf, d_pf_mat ,d_cam_mat,key_dmf)

	return
	endif

!	write(*,*) 'status sst',nstep,kluch,plasma_up,plasma_trig

!       if(key_start .le. 0) then
!	write(*,*) 'fix grid call bkde kluch 0'
	
	if (yesfitcc.eq.0) then
          call  sstepon_bkdw_3( KLUCH, k_auto,nstep,dt,time,
     *                  voltpf, d_pf_mat ,d_cam_mat,key_dmf)
!check psiext
	else
		if (yesfitcc.eq.1)	call  cf_init( k_auto, nstep, dt, time,
     *        voltpf, d_pf_mat,d_cam_mat ) !fit only 12 currents up to psl
		if (yesfitcc.eq.2)	call  cf_init_bkwd( k_auto, nstep, dt, time,
     *       voltpf, d_pf_mat,d_cam_mat ) !now it works
		if (yesfitcc.eq.3)	call  cf_init_full( k_auto, nstep, dt, time,
     *        voltpf, d_pf_mat,d_cam_mat ) !not yet implemented 

	endif
!          call  sstepon_bkdw( KLUCH, k_auto,nstep,dt,time,
!     *                  voltpf, d_pf_mat ,d_cam_mat,key_dmf)

!	call err_catch_a
     
!       else
!	write(*,*) 'refit currents bkdw'
          call  wrrec

        elseif(KLUCH.eq.1) then   !time steping


	if (plasma_up.eq.0) then
!	write(*,*) 'fbkdw'
			
!       k_auto= 0  ! don't change!
	          call  fbkdw_stepon( KLUCH, k_auto,nstep,dt,time,
     *        voltpf, d_pf_mat ,d_cam_mat,key_dmf)

	return
	endif

!	write(*,*) 'status fstep',nstep,kluch,plasma_up,plasma_trig
!	call err_catch_a
!	k_grid=1
	if(k_grid.eq.0) then

	if (isafazt.eq.2) kpr=-1
	if (isafazt.eq.3) kpr=-2
!					write(*,*) 'call fixed grid evol',isafazt,kpr  
	        call  sstepon_bkdw( KLUCH, k_auto,nstep,dt,time,
     *      voltpf, d_pf_mat ,d_cam_mat,key_dmf)
          if(kpr.ge.0) call wrd
!	write(*,*) 'done stepon',isafazt,kpr
				else

	if (j_switch.eq.2) then
!	write(*,*) 'adaptige grid evol', KLUCH, k_auto,nstep_local,
!     & dt,time,isafazt
	if (isafazt.eq.0) kpr=0
	if (isafazt.eq.1) kpr=-2
!          call  f_stepon_bkdw( 1, 0,
!     & 1+nstep_local,dt,time,
!     *     voltpf, d_pf_mat ,d_cam_mat,key_dmf)
		else
!	write(*,*) 'reset rectangular grid'     

!					kpr=0
!          call  sstepon_bkdw000( 0, 0,0,dt,time,
!     *                  voltpf, d_pf_mat ,d_cam_mat,rax,zax,key_dmf)
!					call wrrec
!	call err_catch_a				
!					write(*,*) 'switch to adaptive'     
!          call  f_stepon_bkdw( 0, 0,0,dt,time,
!     *                  voltpf, d_pf_mat ,d_cam_mat,rax,zax,key_dmf)
!					write(*,*) 'adaptive evol'     
!          call  f_stepon_bkdw( 1, 0,1,dt,time,
!     *                  voltpf, d_pf_mat ,d_cam_mat,rax,zax,key_dmf)

!	call err_catch_a
!	j_switch=2
	endif


         if(kpr.ge.0) call f_wrd
         call cur_avg

	endif     

	nstep_local=nstep_local+1
        if(nstep/20*20.eq.nstep) then 
          call wrd_tim
        endif

!        write(*,*) '**'
!        write(*,*) 'nstep,time',nstep,time
!        write(*,*) '**'

        endif        !time steping

	endif




































       return
       end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_zccurb(rc_cur,zc_cur,z2c_cur,
     & rgeoc,zgeoc,ahorc) !,bpcell)

        include 'double.inc'
        include 'dim.inc'
        include 'compol.inc'
	real*8 perimz,ahorc2,ahorc,rgeoc,zgeoc
	real*8 bpcell(iplas-1,nt1-1),avgelem
	integer i_zcur
	
         sum_z2c=0.d0
         sum_zc=0.d0
         sum_rc=0.d0
         sum_Ipla=0.d0
			  	perimz=0.         
	rgeoc=0.
	zgeoc=0.
	ahorc=0.
	ahorc2=0.
!	i_zcur=nint(0.66*iplas)
	i_zcur=nint(1.*iplas)
		
	do i=1,i_zcur

          do j=2,nt1
	sum_rc=sum_rc+(r(i,j)**2.)*cur(i,j-1)*s(i,j-1)
	sum_zc=sum_zc+z(i,j)*cur(i,j-1)*s(i,j-1)
	sum_Ipla=sum_Ipla+cur(i,j-1)*s(i,j-1)
          enddo
	enddo
	
         zc_cur=sum_zc/sum_Ipla
         z2c_cur=zc_cur
         rc_cur=sqrt(sum_rc/sum_Ipla)

!	write(*,*) 'zcur ',zc_cur,sum_ipla,z(1,2),cur(1,1),s(1,1)

!	open(32,file='fort.344')
!		write(32,*) r(1:iplas,2:nt1),z(1:iplas,2:nt1),
!     & cur(1:iplas,1:nt1-1),s(1:iplas,1:nt1-1),
!     & psi(1:iplas,2:nt1)
!	close(32)


!	write(*,*) 'avgeleem'
	do i=1,iplas
          do j=2,nt1
	avgelem=s(i,j-1) !/bpcell(i,j-1)
!	write(*,*) i,j-1,dlt(i,j-1),bpcell(i,j-1),avgelem
	perimz=perimz+avgelem
	rgeoc=rgeoc+r(i,j)*avgelem
	zgeoc=zgeoc+z(i,j)*avgelem
          enddo
	enddo
	rgeoc=rgeoc/perimz
	zgeoc=zgeoc/perimz
	
	do i=1,iplas
          do j=2,nt1
	avgelem=s(i,j-1) !/bpcell(i,j-1)
	ahorc2=ahorc2+(r(i,j)-rgeoc)**2.*avgelem
          enddo
	enddo
	ahorc2=ahorc2/perimz
	ahorc=2.*sqrt(ahorc2)


!	write(*,*) 'rcur',sum_ipla,sum_zc,zc_cur
!	write(*,*) 'rcur',z(i,2:nt1),'A ',bpcell(1:nt1-1),'A ',dlt(i,2:nt1)
!	call err_catch_a
      return
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine find_demo_gaps(ngaps,demo_gaps,geom1d)

        include 'double.inc'
        include 'param.inc'
        include 'comblc.inc'

	integer ngaps,i,j
	integer n_iterz
	double precision demo_gaps(ngaps,4),geom1d(ngaps)
	double precision dumx,dumy,dumz,l_ref,l_gap
	double precision dumx2,dumy2,dumu1,dumu2
	double precision dumx3,dumy3,l_gap2,l_gap3
	double precision dumx0,dumy0,dumxx,dumyy,d_step
	double precision gapmin,gapmax
	double precision tolez,bolez,dur1,dur2
	integer onlypos
	
	tolez=1.e-5

	d_step=0.1 !advance in 1 cm steps	
	n_iterz=100
	gapmin=-1.5
	gapmax=2.5


	do i=1,ngaps
	onlypos=nint(demo_gaps(i,4))
	
!first positive gap!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	l_ref=d_step
			dumx0=demo_gaps(i,1) ! R0
			dumy0=demo_gaps(i,2) ! Z0
			dumz=demo_gaps(i,3) ! angle	 
			
	dur1=0.
	dur2=l_ref

1781	continue

	dumx=dumx0+dur1*cos(dumz)
	dumy=dumy0+dur1*sin(dumz)
	dumx2=dumx0+dur2*cos(dumz)
	dumy2=dumy0+dur2*sin(dumz)

		j1=floor((dumx-r(1))/(r(2)-r(1)))+1	
		j2=floor((dumy-z(1))/(z(2)-z(1)))+1	
		x00=r(j1)
		y00=z(j2)
		x002=r(j1+1)
		y002=z(j2+1)
		u00=blin_(dumx,dumy,x00,x002,y00,y002,u(j1,j2)
     & ,u(j1+1,j2),u(j1+1,j2+1),u(j1,j2+1))
		j3=floor((dumx2-r(1))/(r(2)-r(1)))+1	
		j4=floor((dumy2-z(1))/(z(2)-z(1)))+1	
		x002=r(j3)
		y002=z(j4)
		x003=r(j3+1)
		y003=z(j4+1)
		u002=blin_(dumx2,dumy2,x002,x003,y002,y003,u(j3,j4)
     & ,u(j3+1,j4),u(j3+1,j4+1),u(j3,j4+1))


	if (abs(l_ref).lt.tolez) goto 131
			
	if (u00.eq.up) goto 118
	if (u002.eq.up) goto 117
	if (u002.gt.up.and.u00.lt.up) goto 115
	if (u002.lt.up.and.u00.gt.up) goto 116

	
!case with no intersection: larger 1. m
	l_gap=0.5*(dur1+dur2)
	if (l_gap.ge.gapmax.or.l_gap.lt.gapmin)	goto 119

	dur1=dur1+l_ref
	dur2=dur2+l_ref

	goto 1781
	
	
116	continue
115	continue

	l_ref=-0.5*l_ref
	dur1=dur2
	dur2=dur1+l_ref
	
	
	goto 1781

118	continue
	l_gap=dur1
	goto 211
	
117	continue
	l_gap=dur2
	goto 211

119	continue !case there is no intersection
	l_gap=-5000.
	goto 211
	
131	continue
	l_gap=0.5*(dur1+dur2)

211	continue
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!now negative gap
	l_ref=-d_step 

			dumx0=demo_gaps(i,1) ! R0
			dumy0=demo_gaps(i,2) ! Z0
			dumz=demo_gaps(i,3) ! angle	 
			
	dur1=0.
	dur2=l_ref

3781	continue

	dumx=dumx0+dur1*cos(dumz)
	dumy=dumy0+dur1*sin(dumz)
	dumx2=dumx0+dur2*cos(dumz)
	dumy2=dumy0+dur2*sin(dumz)

		j1=floor((dumx-r(1))/(r(2)-r(1)))+1	
		j2=floor((dumy-z(1))/(z(2)-z(1)))+1	
		x00=r(j1)
		y00=z(j2)
		x002=r(j1+1)
		y002=z(j2+1)
		u00=blin_(dumx,dumy,x00,x002,y00,y002,u(j1,j2)
     & ,u(j1+1,j2),u(j1+1,j2+1),u(j1,j2+1))
		j3=floor((dumx2-r(1))/(r(2)-r(1)))+1	
		j4=floor((dumy2-z(1))/(z(2)-z(1)))+1	
		x002=r(j3)
		y002=z(j4)
		x003=r(j3+1)
		y003=z(j4+1)
		u002=blin_(dumx2,dumy2,x002,x003,y002,y003,u(j3,j4)
     & ,u(j3+1,j4),u(j3+1,j4+1),u(j3,j4+1))


	if (abs(l_ref).lt.tolez) goto 331
			
	if (u00.eq.up) goto 318
	if (u002.eq.up) goto 317
	if (u002.gt.up.and.u00.lt.up) goto 315
	if (u002.lt.up.and.u00.gt.up) goto 316

	
!case with no intersection: larger 1. m
	l_gap2=0.5*(dur1+dur2)
	if (l_gap2.ge.gapmax.or.l_gap2.lt.gapmin)	goto 319

	dur1=dur1+l_ref
	dur2=dur2+l_ref

	goto 3781
	
	
316	continue
315	continue

	l_ref=-0.5*l_ref
	dur1=dur2
	dur2=dur1+l_ref
	
	
	goto 3781

318	continue
	l_gap2=dur1
	goto 411
	
317	continue
	l_gap2=dur2
	goto 411

319	continue !case there is no intersection
	l_gap2=-5000.
	goto 411
	
331	continue
	l_gap2=0.5*(dur1+dur2)

411	continue
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!choose minimum of absolute values
	if (onlypos.eq.0) then
		l_gap3=l_gap
		if (abs(l_gap2).lt.abs(l_gap)) l_gap3=l_gap2

!	write(*,*) 'gaps ',i,dumx0,dumy0,dumx,dumy, 
!     & dumx2,dumy2,dumxx,dumyy,l_gap2,l_gap,l_gap3
!	write(*,*) 'gaps ',i,dumx0,dumy0,dumx,dumy,
!     & dumx2,dumy2,dumxx,dumyy,l_gap,l_ref
		geom1d(i)=min(gapmax,max(gapmin,l_gap3))
!	write(999,'(3E25.11)') l_gap,l_gap2,geom1d(i)
	else
		geom1d(i)=min(gapmax,max(gapmin,l_gap))
	endif

	enddo
	
!	call err_catch_a
	end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!	write(*,*) x002,y002,u002,up,dpsidl,l_gap
!		dumx3=dumx+l_gap*cos(dumz)
!		dumy3=dumy+l_gap*sin(dumz)
!		j33=((dumx3-r(1))/(r(2)-r(1)))	
!		j43=((dumy3-z(1))/(z(2)-z(1)))	
!		x003=r(j33)
!		y003=z(j43)
!		u003=u(j33,j43)
	enddo
	geom1d(i)=l_gap
!		write(*,*) x00,y00,u00,x002,y002,u002
!		write(*,*) x003,y003,u003,up,l_gap,dpsidl
!		write(*,*) geom1d(i)
	write(999,*) geom1d(i)
	enddo

	end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!	write(*,*) x002,y002,u002,up,dpsidl,l_gap
!		dumx3=dumx+l_gap*cos(dumz)
!		dumy3=dumy+l_gap*sin(dumz)
!		j33=((dumx3-r(1))/(r(2)-r(1)))	
!		j43=((dumy3-z(1))/(z(2)-z(1)))	
!		x003=r(j33)
!		y003=z(j43)
!		u003=u(j33,j43)
	enddo
	geom1d(i)=l_gap
!		write(*,*) x00,y00,u00,x002,y002,u002
!		write(*,*) x003,y003,u003,up,l_gap,dpsidl
!		write(*,*) geom1d(i)
	write(999,*) geom1d(i)
	enddo

	end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine find_demo_meas2021(nmea1,nmea2,
     & bri1,bro1,fli1,flo1)

        include 'double.inc'
        include 'param.inc'
        include 'comblc.inc'

	integer nmea1,nmea2,i,j,j1,j2
	double precision bri1(nmea1,3)
	double precision bro1(nmea1)
	double precision fli1(nmea2,2)
	double precision flo1(nmea2)
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
		bro1(i)=-1/r(j1)*(u(j1,j2+1)-u(j1,j2))/
     & (z(2)-z(1))*cos(dumz)
		bro1(i)=bro1(i)+1/r(j1)*(u(j1+1,j2)-u(j1,j2))/
     & (r(2)-r(1))*sin(dumz)
	enddo

	do i=1,nmea2

			dumx=fli1(i,1) ! R
			dumy=fli1(i,2) ! Z
			
		j1=((dumx-r(1))/(r(2)-r(1)))	
		j2=((dumy-z(1))/(z(2)-z(1)))	

		x00=r(j1)
		y00=z(j2)
		flo1(i)=u(j1,j2)
	enddo
	end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine find_demo_meas2021vacuum(nmea1,nmea2,
     & bri1,bro1,fli1,flo1,pcequi,ncequi,i_read)

        include 'double.inc'
        include 'param.inc'
        include 'comblc.inc'

	integer nmea1,nmea2,i,j,j1,j2
	double precision bri1(nmea1,3)
	double precision bro1(nmea1)
	double precision fli1(nmea2,2)
	double precision flo1(nmea2)
	double precision dumx,dumy,dumz,l_ref,l_gap
	double precision dumx2,dumy2
	double precision dumx3,dumy3
				integer ncequi,i_read
         real*8 pcequi(ncequi)

	write(*,*) 'meas'

	if (i_read.eq.1)	then
        open(1,file='exp/equ/dem_/data.dat')
             read(1,*) i
             read(1,*) ni
             read(1,*) nj
             read(1,*) rmin
             read(1,*) rmax
             read(1,*) zmin
             read(1,*) zmax
        close(1)
		ni1=ni-1
		nj1=nj-1
         ddr=(rmax-rmin)/dfloat(ni1)
         ddz=(zmax-zmin)/dfloat(nj1)
         r(1)=rmin
         z(1)=zmin
         do 10 i=2,ni
 10      r(i)=r(i-1)+ddr
         do 20 j=2,nj
 20      z(j)=z(j-1)+ddz
	call extfil(pcequi,ncequi)
	endif

!	write(*,*) r(1),z(1)
!	write(333,'(4E25.11,2I,1E25.11)') 0
	
	call ext_fil(pcequi,ncequi)

	do i=1,nmea1

			dumx=bri1(i,1) ! R
			dumy=bri1(i,2) ! Z
			dumz=bri1(i,3) ! angle	 
			
		j1=((dumx-r(1))/(r(2)-r(1)))+1	
		j2=((dumy-z(1))/(z(2)-z(1)))+1	

		x00=r(j1)
		y00=z(j2)
		bro1(i)=-1/r(j1)*(ue(j1,j2+1)-ue(j1,j2))/
     & (z(2)-z(1))*cos(dumz)
		bro1(i)=bro1(i)+1/r(j1)*(ue(j1+1,j2)-ue(j1,j2))/
     & (r(2)-r(1))*sin(dumz)
!	write(333,'(4E25.11,2I,1E25.11)') dumx,dumy,x00,y00,j1,j2,ue(j1,j2)
	enddo

	do i=1,nmea2

			dumx=fli1(i,1) ! R
			dumy=fli1(i,2) ! Z
			
		j1=((dumx-r(1))/(r(2)-r(1)))+1	
		j2=((dumy-z(1))/(z(2)-z(1)))+1	

		x00=r(j1)
		y00=z(j2)
!	write(333,'(4E25.11,2I,1E25.11)') dumx,dumy,x00,y00,j1,j2,ue(j1,j2)

		flo1(i)=ue(j1,j2)
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	subroutine sol_current_model(nsol,znocur)

	implicit none
	
	integer nc,nsol,i,j,k
	double precision rc(500),zc(500),rt,zt,znocur(nsol)
	integer jiter
	common /com_ind/ nr,nt,nr1,nt1,nr2,nt2,iplas,iplas1
	common /com_rz/ r,z,ro,teta
	data jiter/0/
	save jiter,rc,zc,nc
	

	integer nr,nt,nr1,nt1,nr2,nt2,iplas,iplas1
	real*8 r(256,1556),z(256,1556)
	real*8 ro(256,1556),teta(256,1556)
	real*8 rgeo,zgeo,anglt(1556),discrim
	real*8 anglf(1556),diffa(1556)
	real*8 anglr(1556),dum1

	if (jiter.eq.0) then
	jiter=1
	open(32,file='exp/equ/aug_/blansoldiv.dat')
	read(32,*) nc
	do i=1,nc
	read(32,*) rc(i),zc(i),dum1,dum1,dum1	
	enddo
	close(32)
	endif

	znocur(1:nsol)=1.0
!	write(888,*) nc,rc(1:nc),'$',zc(1:nc)
!	write(888,*) ' '
	do i=1,nsol
	rt=rc(i)
	zt=zc(i)	

!check if sol current filament is inside plasma
!	rgeo=sum(r(iplas,1:nt))/nt
!	zgeo=sum(z(iplas,1:nt))/nt
	do j=2,nt-1

	if (-(zt-z(iplas,j)).ge.0
     & .and.
     & -(rt-r(iplas,j)).ge.0) then
	anglr(j)=atan((zt-z(iplas,j))
     & /(rt-r(iplas,j)))	
	endif
	if (-(zt-z(iplas,j)).ge.0
     & .and.
     & -(rt-r(iplas,j)).lt.0) then
	anglr(j)=3.141592+atan((zt-z(iplas,j))
     & /(rt-r(iplas,j)))	
	endif
	if (-(zt-z(iplas,j)).lt.0
     & .and.
     & -(rt-r(iplas,j)).lt.0) then
	anglr(j)=3.141592+atan((zt-z(iplas,j))
     & /(rt-r(iplas,j)))	
	endif
	if (-(zt-z(iplas,j)).lt.0
     & .and.
     & -(rt-r(iplas,j)).ge.0) then
	anglr(j)=3.141592*2+atan((zt-z(iplas,j))
     & /(rt-r(iplas,j)))	
	endif
	if (j.gt.2.and.anglr(j).lt.anglr(j-1)) then
		anglr(j)=anglr(j)+2*3.141592
	endif
	enddo

	do j=2,nt-1
	if (j.gt.2)	diffa(j)=(anglr(j)-anglr(j-1))
	if (j.eq.2)	diffa(j)=(anglr(j)-anglr(nt-1)+2*3.141592)
	enddo

	anglt(2:nt-1)=sin(anglr(2:nt-1))*diffa(2:nt-1)
	anglf(2:nt-1)=cos(anglr(2:nt-1))*diffa(2:nt-1)	

	discrim=abs(sum(anglt(2:nt-1))/(nt-2))+
     & 	abs(sum(anglf(2:nt-1))/(nt-2))

	if (i.eq.7) then
	do j=2,nt-1
!	write(885,'(6E25.11)') r(iplas,j),z(iplas,j),
!     & anglr(j),anglt(j),anglf(j),diffa(j)
	enddo
!	write(*,*) 'fdfd',rt,zt
	endif
	



	if (discrim.lt.0.1) znocur(i)=0.

!	write(888,*) rt,zt,discrim,znocur(i)

! end of cycle
	enddo
!	write(888,*) ' '
!	do i=2,nt-1
!		write(889,*) r(iplas,i),z(iplas,i)
!	enddo		
!	call err_catch_a	
	znocur=0.
!
	end




	subroutine limpotential_generate(rcent,zcent,nl,rl,zl,znocur)

	use fenix_params
	implicit none

!	include 'includes_fsim.dat'

	
	integer nc,nsol,i,j,k,nl
	double precision rl(nl),zl(nl),rt,zt,znocur
	double precision rcent,zcent
	integer jiter
	

	integer nr,nt,nr1,nt1,nr2,nt2,iplas,iplas1
	real*8 r(256,1256),z(256,1256)
	real*8 ro(256,1256),teta(256,1256)
	real*8 rgeo,zgeo,anglt(1256),discrim
	real*8 anglf(1256),diffa(1256)
	real*8 anglr(1256),dum1



	znocur=-1.

	rt=rcent
	zt=zcent	

!check if sol current filament is inside plasma
	do j=1,nl

	if (-(zt-zl(j)).ge.0
     & .and.
     & -(rt-rl(j)).ge.0) then
	anglr(j)=atan((zt-zl(j))
     & /(rt-rl(j)))	
	endif
	if (-(zt-zl(j)).ge.0
     & .and.
     & -(rt-rl(j)).lt.0) then
	anglr(j)=3.141592+atan((zt-zl(j))
     & /(rt-rl(j)))	
	endif
	if (-(zt-zl(j)).lt.0
     & .and.
     & -(rt-rl(j)).lt.0) then
	anglr(j)=3.141592+atan((zt-zl(j))
     & /(rt-rl(j)))	
	endif
	if (-(zt-zl(j)).lt.0
     & .and.
     & -(rt-rl(j)).ge.0) then
	anglr(j)=3.141592*2+atan((zt-zl(j))
     & /(rt-rl(j)))	
	endif
	if (j.ge.2) then
		if (anglr(j).lt.anglr(j-1)) then
			anglr(j)=anglr(j)+2*3.141592
		endif
	endif
	
	enddo

	nt=nl
	do j=1,nt
	if (j.gt.1)	then
			diffa(j)=abs(anglr(j)-anglr(j-1))
			if (diffa(j).gt.1.) diffa(j)=abs(diffa(j)-2.*3.141592)
		endif
	if (j.eq.1)	then
	diffa(j)=abs(anglr(j)-anglr(nt))
			if (diffa(j).gt.1.) diffa(j)=abs(diffa(j)-2.*3.141592)
	endif
!	write(555,*) rcent,zcent,j,rl(j),zl(j),anglr(j),diffa(j)
	enddo

	anglt(1:nt)=sin(anglr(1:nt))*diffa(1:nt)
	anglf(1:nt)=cos(anglr(1:nt))*diffa(1:nt)	

	discrim=abs(sum(anglt(1:nt)))/sum(diffa(1:nt))+
     & 	abs(sum(anglf(1:nt)))/sum(diffa(1:nt))


	if (discrim.lt.0.6) znocur=1.

	if (use_zlim_pot.eq.0) znocur=1.


!	write(555,*) 'res',rcent,zcent,discrim,znocur


	end
