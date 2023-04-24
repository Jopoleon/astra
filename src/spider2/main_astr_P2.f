 !uses new interface with cpos
	     
         program main_astr_fixed
 
         use schemas       
         use parameters 
      !   use durs_d_modul
               
      implicit real*8 (a-h,o-z) !!!
       
          integer ni_p, nj_p, np, nb1, nbtabp, nbtabp2
	
          parameter(ni_p=129,nj_p=257)
          parameter(np=1000,nb1=np+1) 
          parameter(nbtabp=1000,nbtabp2=nbtabp*2)
      include 'dimpl1.inc'
       common/savt0/ psi0(nrp),fi0(nrp),f0(nrp),ri0(nrp),q0(nrp),
     *               dpsidt(nrp),dfidt(nrp),rm0,ac0n,skcen0

        real*8 rbtab(nbtabp),zbtab(nbtabp),nbtab
        real*8 Cc(np),CUbs(np),Te(np),Cd(np),pres(np)
        real*8 Cu(np)
        real*8 Cu_out(np)
        real*8 pres_s(np)
        real*8 ymu(np),yfp(np)

        real*8 rzbnd(nbtabp2),eqpf(np),eqff(np),fp(np),ipl,ybetpl,yli3,
     *       rtor,btor,rho(np),roc,eqpf0(np),
     *       g11(np),g22(np),g33(np),ccoils(12),                          
     *       vr(np),vrs(np),slat(np),gradro(np),rocnew, 
     *       mu(np),ipol(np),bmaxt(np),bmint(np),bdb02(np),bdb0(np),
     *	     b0db2(np),droda(np),W_Dj(np),yFOFB(np)
        real*8 roc_input                      

        character*40 prename
        character*40 eqdfn
               
       type(type_equilibrium) equil_in
       type(type_equilibrium) equil_out0
       type(type_equilibrium) equil_out
       type(type_parameters) parameters_spider
       
       !locals
       integer neql,nteta,nbnd,na1,nstep,i,j
       real*8 hro, yreler, platok,time, dt
       integer key_equil
       integer key_plc
       integer key_out
       integer kfixfree
       common/com_flag_ef/kastr2
	integer kastr2       
       !voltages
       real*8 vcoils(20)
       real*8 pfc(20),CCOIL(20)
       integer ncoils,dumi1
       
       integer NCNB,NCTP
       double precision yvcoilx(20) ,yccoilx(20),t_starts, TIMEC, TAU
       double precision DUMCTP(13), VCOIL(20)
       double precision CTRLM(3), DUMCT(13),t_max
       double precision t_cde,iplx,pkin,pkin0
        character*4 machine_name

!	machine_name='aug '
!	machine_name='benc'

	kastr2=1   
!         open(1,file='exp/equ/aug/spidat.dat')  
         open(1,file='spidat.dat')  
	     read(1,*) 
	     read(1,*) neql,nteta,nbnd
	     read(1,*) 
	     read(1,*) (rzbnd(i),i=1,nbnd*2)
	     read(1,*) 
	     read(1,*) na1
	     read(1,*) 
	     read(1,*) (eqpf(i),i=1,na1)
	     read(1,*) 
	     read(1,*) (eqff(i),i=1,na1)
	     read(1,*) 
	     read(1,*) (fp(i),i=1,na1)
	     read(1,*) 
	     read(1,*) (rho(i),i=1,na1)
	     read(1,*) 
	     read(1,*) ipl,rtor,btor,roc,nstep
	     read(1,*) 
 	     read(1,*) (fp(i),i=1,na1)
	     read(1,*) 
 	     read(1,*) (g11(i),i=1,na1)
	     read(1,*) 
 	     read(1,*) (g22(i),i=1,na1)
	     read(1,*) 
 	     read(1,*) (g33(i),i=1,na1)
	     read(1,*) 
 	     read(1,*) (vr(i),i=1,na1)
	     read(1,*) 
 	     read(1,*) (vrs(i),i=1,na1)
	     read(1,*) 
 	     read(1,*) (slat(i),i=1,na1)
	     read(1,*) 
 	     read(1,*) (gradro(i),i=1,na1)
	     read(1,*) 
 	     read(1,*) (mu(i),i=1,na1)
	     read(1,*) 
	     read(1,*) (ipol(j),j=1,na1)

	 read(1,*)                  !    	 
	 read(1,*) (cc(j),j=1,na1)  !      
	 read(1,*) 	              !      
	 read(1,*) (Te(j),j=1,na1)  !      
	 read(1,*) 	              !      
	 read(1,*) (cubs(j),j=1,na1)!      
	 read(1,*) 	              !      
	 read(1,*) (cd(j),j=1,na1)  !        
	 read(1,*) 	              !      
	 read(1,*) (pres(j),j=1,na1)!        
	 read(1,*) 	               !     
	 read(1,*) (cu(j),j=1,na1)   !     
        close(1)

!	ipl=0.09

        nstep=0
        !key_equil=-10
        !key_equil=-2
        key_equil=0

        key_plc=1
        call put_key_plc(key_plc)
        
        key_out=1

        call put_key_out(key_out)

	write(*,*) na1,nbnd,cc(1)

	iplx=ipl

	eqpf0=eqpf

!	nteta=nbnd+2
!			neql=61
!        parameters_spider%neql = neql
!        parameters_spider%nteta = nteta


!        allocate( equil_in%profiles_1d%jni%value(na1) )
!        allocate( equil_in%profiles_1d%te%value(na1) )
!        allocate( equil_in%profiles_1d%sigmapar%value(na1) )
!        allocate( equil_in%profiles_1d%q(na1) )
!        allocate( equil_in%profiles_1d%pressure(na1) )
!        allocate( equil_in%profiles_1d%rho_tor(na1) )
!        equil_in%profiles_1d%jni%value = cubs(1:na1)
!        equil_in%profiles_1d%te%value = Te(1:na1)
!        equil_in%profiles_1d%sigmapar%value = cc(1:na1)
!								equil_in%profiles_1d%pressure= pres(1:na1)
!								equil_in%profiles_1d%rho_tor= rho(1:na1)
!        do i =1,na1
!												equil_in%profiles_1d%q(i) = 1.d0/mu(i)
!        enddo
       parameters_spider%k_fixfree=1
        
          call spidat2itm(neql,nteta,nbnd,rzbnd,
     *                    na1,eqpf,eqff,fp,ipl,
     *                    rtor,btor,rho,roc,mu,
     *                    cc,Te,cubs,cd,pres,cu,
     *                    equil_in,parameters_spider )
!	equil_in%global_param%toroid_field%b0=btor


!	write(*,*) pres(1:10)
	     
!        input from astra to spider
!
!        all values are input parameters 
!        ALL UNITS ARE LIKE IN ASTRA
!
!    nbnd  -  number of boundary points
!    rzbnd  - 1D array of r(1:nbnd) and z(1:nbnd) coordinates of boundary points.
!             rzbnd size is 2*nbnd
!    na1 - astra radial grid size
!    key_equil - control key for spider
!              if key_equil=0 - equilibrium with prescribed p' and ff' profiles and toroidal plasma current; 
!              the following inputs are used (rho,cu,pres are not used): 
!    fp - psi 1D  array,fp size is na1
!    eqpf - p'*rtor 1D array,yeqpf size is na1
!    eqff - ff'/rtor 1D array,yeqff size is na1
!    ipl  -  toroidal plasma current 
!              if key_equil=-10 - equilibrium with prescribed p' and <(j,B)>;
!              the following inputs are used (fp,eqpf,eqpp,ipl are not used): 
!    rtor,btor - btor is vacuum toroidal field at r=rtor
!    rho - normalized toroidal flux 1D array, sqrt(Phi/btor) rho size is na1
!    cu - <(j,B)>/btor 1D array, cu size is na1
!    pres - pressure 1D array,pres size is na1
!    nstep - time step number,
!                              nstep=0  initial equilibrium computed (should always preceed nstep>0 run)
!                              nstep>0   previously computed equilibrium is used as initial guess 
!    parameters_spider - SPIDER parameters (exchange
!              defaults are read and SPIDER  parameters are setup for equilibrium problem to solve
        
!        ncoils = 12
        ncoils = 12
	parameters_spider%prename='./'
	parameters_spider%kname=len_trim(parameters_spider%prename)

	write(*,*) parameters_spider%prename
	write(*,*) parameters_spider%kname
        parameters_spider%key_pres = 0
	
	open(32,file='currents.wr')
	read(32,*) dumi1,dumi1
	read(32,*) (ccoils(i),i=1,ncoils)
	close(32)
	write(*,*) ccoils(1:ncoils)	

!	call coil2spider(ccoils*1.E3,ncoils,parameters_spider)	

	parameters_spider%nstep = nstep	
	parameters_spider%kpr=0
	write(*,*) ipl,roc	,	parameters_spider%kpr				
	parameters_spider%k_grid=0
						  vcoils = 0
       kfixfree=parameters_spider%k_fixfree 

       call spider_run(ncoils,vcoils,
     &                  equil_in,equil_out0,parameters_spider)     




	stop

	pkin0=equil_out0%global_param%wkin
					
	pkin=pkin0					
					

	write(*,*) '1'						

	
!	call spidupdate
C	pause
!        run spider with parameters setup by input2spider
!        input/output equilibrium parameters are to be extended to CPO finally
!
!       equil_in - dummy ASTRA equilibrium (input
!       equil_out - SPIDER equilibrium (toroidal flux phi only by now) (output          
!       parameters_spider - SPIDER parameters (input          

!       rho_boundary=rocnew from the spider equilibrium to be conrolled 
!       and possibly new rho grid specified (e.g. the last mesh step 
!       added/removed/modififed)
!
!	  neql = parameters_spider%neql
!       rocnew = sqrt(equil_out%profiles_1d%phi(neql)/(PI*btor))
!       roc_input = rho(na1)
!       do i=1,na1
!           rho(i) = rho(i)*roc_new/roc_input
!       enddo
	time=0.
	do jj=1,na1
	write(434,'(5E25.11)') time,rho(jj),fp(jj),roc,cu(jj)
	enddo

	write(*,*) '2'						
!        call spider2output(
!     *                     equil_out0,parameters_spider,
!     *                     rtor,btor,rho,roc,na1,
!     *                     g11,g22,g33,vr,vrs,slat,gradro,rocnew,
!     *                     ymu,ipol,bmaxt,bmint,bdb02,b0db2,bdb0,droda,
!     *				       yreler,yli3,platok,cu,fp,
!     *                     pres_s)
!           output from spider to astra


!
!         main astra grid:
!          rho(j)=h*(j-0.5) j=1,2,...,na1-1
!          rho(na1)=roc 
!
!         intermediate astra grid:
!          rho_h(j)=0.5*(rho(j)+rho(j+1)) , j=1,2,...,na1-2
!          rho_h(na)=rho_h(na-1)+h, h=rho(2)-rho(1)
!          rho_h(na1)=roc 
!
!     equil_out -  supposed to be input SPIDER equilibrium but stores the SPIDER magnetic surfaces by now (output
!     parameters_spider - SPIDER parameters (input          
!     rtor,btor -  described above (input
!     rho - normalized toroidal flux 1D array,rho size is na1 (input
!
!     g11,g22,g33,vr,vrs,slat,gradro -  1D arrays of metric coef.size is na1 (output
!     g11,g22,vrs,gradro,slat - are defined on intermediate grid 
!     vr,g33 - are defined on main grid 
!
!     roc - boundary value of normalized toroidal flux (input
!     rocnew - boundary value of normalized toroidal flux (output
!
!     ymu - rotational transform 1D array, size is na1 (output
!     ymu - is defined on intermediate grid
!
!     ipol - =f/rtor/btor - normalized poloidal current 1D array, size is na1 (output 
!     ipol - is defined on main grid
! 
!     bmaxt,bmint,bdb02,b0db2,bdb0,droda-  1D arrays of averaged values, size is na1 (output 
!     bmaxt,bmint,bdb02,b0db2,bdb0 - are defined on intermediate grid 
!
!     yli3 - li3 (output
!     yreler - bet_pol (output
!     platok  -  toroidal plasma current (output
!
!     cu_out  - <(j,B)>/B0 1D array, cu_out size is na1 (output 
!     cu_out - are defined on main grid
! 
!     yfp - psi 1D  array, size is na1 (output
!     yfp - are defined on main grid
!
!     pres_s - pressure.1D array, pres_s size is jna1 (output
!     pres_s - are defined on main grid
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!time !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!! stepping!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	stop	
	write(*,*) '3'						
        parameters_spider%time=0.004
								t_cde= 20.0022
        time = parameters_spider%time
        !dt = parameters_spider%dt
        dt = 1.e-3
								t_max=0.9
        parameters_spider%dt = dt
        
        nstep=0

        NCNB=ncoils
        NCTP=13
        CTRLM = (/3.,3.,101./)
       DUMCT=(/800000.0000, 2.1400, 1.6500, 1.6300,
     &          0.0900, 0.0250, 0.5000, -1.3500,
     &          -1.0500, 0.0000, 1.1200, 0.0000, 0.0000/)
!       CCOIL=(/ -23.9072,  -8.9162 , -0.8602,  7.7685,
!     &   12.8620 , -5.5518 , -10.2462  ,-4.9744 ,
!     &   18.3100,  15.4312 , -14.4453 , -11.8908/)  
        
      if(kfixfree.ne.0) then       
        
!        CTRLM = (/3.,3.,101./)
!       DUMCT=(/800000.0000, 2.1400, 1.6500, 1.6300,
!     &          0.0900, 0.0250, 0.5000, -1.3500,
!     &          -1.0500, 0.0000, 1.1200, 0.0000, 0.0000/)
!        CTRLM = (/3.,6.,101./)
!       DUMCT=(/800000.0000, 2.1400, 1.6500, 1.6300,
!     &          0.0900, 0.0250, 0.85000, -1.3500,
!     &          -1.0500, 0.0000, 1.1200, 0.0000, 0.0000/)

        TIMEC = time
        TAU = dt
        t_starts = 0.d0
        yvcoilx=0.d0
        DUMCTP=DUMCT
	write(*,*) '4'						
        
	yccoilx(1:ncoils)=CCOIL(1:ncoils)
C        call get_pfccur(yccoilx,ncoils)
        !in kA
C        yccoilx = yccoilx*1e-3
        
        open(111, file='exp/equ/benc/tim_tim.wr', status='replace')            
        close(111)
        open(111, file='exp/equ/benc/pfc_tim.wr', status='replace')            
        close(111)
        open(111, file='exp/equ/benc/ctp_tim.wr', status='replace')           
        close(111)
        open(111, file='exp/equ/benc/vol_tim.wr', status='replace')           
        close(111)
        open(111, file='exp/equ/benc/elo_tim.wr', status='replace')           
        close(111)
      endif  
        

	pkin0=equil_out0%global_param%wkin
        
        do i=1,nint((t_max-time)/dt)
!   	  do i=1,3

        parameters_spider%key_dmf = -2
        key_equil=-2
        parameters_spider%key_pres = 0

        key_plc=1
        call put_key_plc(key_plc)
        key_out=1
        call put_key_plc(key_out)

!	if (time.le.t_cde) key_equil=-3
!	if (time.le.t_cde) parameters_spider%key_dmf = -3
	if (time.le.t_cde) key_equil=0
	if (time.le.t_cde) parameters_spider%key_dmf = 0
	if (key_equil.eq.-3 .or. key_equil.eq.0) 
     &  equil_out0%global_param%i_plasma=iplx*1e6
	if (key_equil.eq.-3 .or. key_equil.eq.0) 
     &  ipl=iplx
!	if (time.gt.t_cde) dt=2.e-4

       parameters_spider%dt = dt

					 	write(*,*) 'ipl: ',ipl
        
              nstep=nstep+1
  
              time=time+dt
															write(*,*) 'time = ',time,roc


!								equil_in%profiles_1d%pressure= pres(1:na1)
!								equil_out%profiles_1d%pressure= pres(1:na1)

!	equil_out%profiles_1d%pressure=equil_in%profiles_1d%pressure
      

!           the spider parameters can be changed during the run, in particular in addition to
!              if key_dmf=0 - equilibrium with prescribed p' and ff' profiles and toroidal plasma current; 
!              if key_dmf=-10 - equilibrium with prescribed p and <(j,B)>;
!           other equilibrium problems can be solved (with the profiles taken from the intermediate equilibria) 
!              if key_dmf=2 - equilibrium with prescribed p',q and poloidal flux in plasma 
!              if key_dmf=-2 - equilibrium + magnetic field diffusion with f and psi_ext at boundary
!              if key_dmf=-3 - equilibrium + magnetic field diffusion with f at boundary and toroidal plasma current

        parameters_spider%nstep = nstep
        parameters_spider%kpr = 0
        parameters_spider%time = time
        parameters_spider%key_plc = 1
        parameters_spider%key_out = 1
        parameters_spider%epsro = 1.0d-7
        parameters_spider%enels = 5.e-7
        parameters_spider%k_auto = 0
        parameters_spider%k_grid = 0
        
        ! needed for input2spider setup (dummy equil_in)
        call put_enels(parameters_spider%enels)
        call put_key_plc(parameters_spider%key_plc)
        call put_key_plc(parameters_spider%key_out)
       
      if(kfixfree.ne.0) then       

        call get_rmax(DUMCTP(2),DUMCTP(11),DUMCTP(4),DUMCTP(7),
     &                DUMCTP(10))                                !(r_max,r_min,r_geo,z_max,delta_up)
        call get_zccur(DUMCTP(3),DUMCTP(5),DUMCTP(6))            !(rc_cur,zc_cur,z2c_cur)
        
        !Zsquad = Zcurr
        DUMCTP(6) = sqrt(DUMCTP(6))
        DUMCTP(6) = DUMCTP(5)
        
!        call get_zbran(1.27d0,DUMCTP(9))
!        call get_zbran(1.72d0,DUMCTP(8))
        call get_pfccur(pfc,ncoils)
	write(*,*) pfc(1:ncoils)
C	pause
        DUMCTP(12)=pfc(9)
        DUMCTP(13)=pfc(10)
        
        !plasma current from adaptive grid
        DUMCTP(1)=equil_out0%global_param%i_plasma*1e6

        
        TIMEC = time
               
C        call CTAUG_STD(NCNB,NCTP,CTRLM,
C     &       yvcoilx,yccoilx,t_starts,CCOIL,
C     &       TIMEC,TAU,DUMCT,DUMCTP,VCOIL)
     
        vcoils(1:ncoils) = VCOIL(1:ncoils)
        

	write(*,*) 'vcoils: ',vcoils(1:ncoils),DUMCTP(1)
        !vcoils(4) = 0.
        !vcoils(2:3) = 0.
        !vcoils(1:8) = 0.
	vcoils=0.

    
       endif         


	do jj=1,na1
	write(434,'(5E25.11)') time,rho(jj),fp(jj),roc,cu(jj)
	enddo

	write(*,*) 'stuff from out :',yli3,platok
        
	write(*,*) equil_out0%profiles_1d%psi(1),
     > equil_out0%profiles_1d%psi(neql)
	write(*,*) equil_out0%profiles_1d%pprime(1),
     > equil_out0%profiles_1d%pprime(neql)
	write(*,*) equil_out0%profiles_1d%ffprime(1),
     > equil_out0%profiles_1d%ffprime(neql)
	write(*,*) equil_out0%global_param%toroid_field%r0,
     & equil_out0%global_param%toroid_field%b0,
     & equil_out0%global_param%i_plasma   !itm is in A
								
        call spider_run(ncoils,vcoils,
     &                  equil_out0,equil_out,parameters_spider)
!	stop
!	call f_spidupdate

	write(179,'(6E25.12)') time,equil_out%global_param%li3,
     &  equil_out%global_param%betpol,
     &  equil_out%global_param%wkin,
     &  equil_out%global_param%bpkin,
     &  equil_out%profiles_1d%volume(neql)

	pkin=equil_out%global_param%wkin
	
	write(*,*) 'pkin0/pkin ',pkin0/pkin
	write(3211,'(4E25.11)') time,pkin,equil_out%global_param%li3
     &  ,equil_out%global_param%betpol
	
C	pause          
        call wrd_tim
        open(111, file='exp/equ/benc/tim_tim.wr', access='append')            
        write(111,'(1p12e20.12)') time
        close(111)
        open(111, file='exp/equ/benc/pfc_tim.wr', access='append')            
        write(111,'(1p20e20.12)') pfc(1:ncoils)
        close(111)
        open(111, file='exp/equ/benc/ctp_tim.wr', access='append')            
        write(111,'(1p13e20.12)') DUMCTP(1:13)
        close(111)
        open(111, file='exp/equ/benc/vol_tim.wr', access='append')            
        write(111,'(1p20e20.12)') vcoils(1:ncoils)
        close(111)
        open(111, file='exp/equ/benc/elo_tim.wr', access='append')            
        write(111,'(1p12e20.12)') equil_out%profiles_1d
     &       %elongation(size(equil_out%profiles_1d%elongation))
        close(111)


	write(*,*) 'platok: ',platok
!	do jj=1,nrp
!	write(444,'(3E25.11)') time,fi0(jj),psi0(jj)
!	enddo
!	do jj=1,na1
!	write(445,'(6E25.11)') time,g22(jj),rho(jj),rocnew,vrs(jj),pres_s(jj)
!	enddo
        enddo


        stop
        end        


        
         subroutine spidat2itm(neql,nteta,nbnd,rzbnd,
     *                    na1,eqpf,eqff,fp,ipl,
     *                    rtor,btor,rho,roc,mu,
     *                    cc,Te,cubs,cd,pres,cu,
     *                    equil_in,parameters_spider )

         use schemas       
         use parameters 
               
         implicit real*8(a-h,o-z)
       
        real*8 rzbnd(*),eqpf(*),eqff(*),fp(*),ipl,
     *         rtor,btor,rho(*),roc,mu(*),
     *         cc(*),Te(*),cubs(*),cd(*),pres(*),cu(*) 
               
       type(type_equilibrium) equil_in
       type(type_parameters) parameters_spider
       
       !locals
       integer neql,nteta,nbnd,na1,nstep,i,j
       
      real(8), allocatable :: arrw(:), rhow(:)
      allocate( arrw(na1+1), rhow(na1+1) )
     
         amu0=2.d0*TWOPI/10.d0
 
        parameters_spider%neql = neql
        parameters_spider%nteta = nteta
        equil_in%global_param%toroid_field%r0 = rtor
        equil_in%global_param%toroid_field%b0 = btor
        equil_in%global_param%i_plasma = ipl*1.E6
 
 !! exstrapolation on magn. axis
       
             rh0=0.d0
		              
             rh1=rho(1)           
             rh2=rho(2)           
             rh3=rho(3)           

             fpx1=fp(1)
             fpx2=fp(2)
             fpx3=fp(3)
             
       call EXTRP2(rh0,fpx0, rh1,rh2,rh3, fpx1,fpx2,fpx3)
       
           arrw(1)=0.d0
           rhow(1)=0.d0
          do i=2,na1+1
           arrw(i) =fp(i-1)-fpx0
           rhow(i)=rho(i-1)
          enddo
          
        allocate( equil_in%profiles_1d%psi(na1+1) )
        allocate( equil_in%profiles_1d%rho_tor(na1+1) )
        equil_in%profiles_1d%psi = arrw
        equil_in%profiles_1d%rho_tor = rhow
        		
             ppx1=eqpf(1)
             ppx2=eqpf(2)
             ppx3=eqpf(3)
             
       call EXTRP2(rh0,ppx0, rh1,rh2,rh3, ppx1,ppx2,ppx3)
        		
            arrw(1)=ppx0/Rtor/TWOPI*1.d6
          do i=2,na1+1
           arrw(i) =eqpf(i-1)/Rtor/TWOPI*1.d6
          enddo
       		
         allocate( equil_in%profiles_1d%pprime(na1+1) )
         equil_in%profiles_1d%pprime = -arrw
        		
             ffx1=eqff(1)
             ffx2=eqff(2)
             ffx3=eqff(3)

       call EXTRP2(rh0,ffx0, rh1,rh2,rh3, ffx1,ffx2,ffx3)

           arrw(1)=ffx0*Rtor/TWOPI*amu0
          do i=2,na1+1
           arrw(i)=eqff(i-1)*Rtor/TWOPI*amu0
          enddo		
         allocate( equil_in%profiles_1d%ffprime(na1+1) )
         equil_in%profiles_1d%ffprime = -arrw
 
             ccx1=cc(1)
             ccx2=cc(2)
             ccx3=cc(3)
             
        call EXTRP2(rh0,ccx0, rh1,rh2,rh3, ccx1,ccx2,ccx3)
            
              arrw(1)=ccx0
          do i=2,na1+1
           arrw(i)=cc(i-1)
          enddo		
         allocate( equil_in%profiles_1d%sigmapar%value(na1+1) )
         equil_in%profiles_1d%sigmapar%value = arrw
          
             cbx1=cubs(1)
             cbx2=cubs(2)
             cbx3=cubs(3)
             
       call EXTRP2(rh0,cbx0, rh1,rh2,rh3, cbx1,cbx2,cbx3)
             
              arrw(1)=cbx0
          do i=2,na1+1
           arrw(i)=cubs(i-1)
          enddo		
         allocate( equil_in%profiles_1d%jni%value(na1+1) )
         equil_in%profiles_1d%jni%value = arrw
             
             cex1=te(1)
             cex2=te(2)
             cex3=te(3)
             
        call EXTRP2(rh0,cex0, rh1,rh2,rh3, cex1,cex2,cex3)
            
              arrw(1)=cex0
          do i=2,na1+1
           arrw(i)=te(i-1)
          enddo		
         allocate( equil_in%profiles_1d%te%value(na1+1) )
         equil_in%profiles_1d%te%value = arrw
            
         allocate( equil_in%eqgeometry%boundary%r(nbnd) )
         allocate( equil_in%eqgeometry%boundary%z(nbnd) )
            
	do ib=1,nbnd
	  equil_in%eqgeometry%boundary%r(ib) = rzbnd(ib)
        equil_in%eqgeometry%boundary%z(ib) = rzbnd(ib+nbnd)
      enddo		
        equil_in%eqgeometry%boundary%npoints = nbnd
        
        return
        end
