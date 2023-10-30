      
         program main_astr_fixed
 
         use schemas       
         use parameters 
      !   use durs_d_modul
               
         implicit none      
       
          integer ni_p, nj_p, np, nb1, nbtabp, nbtabp2
          parameter(ni_p=256,nj_p=256)
          parameter(np=1000,nb1=np+1) 
          parameter(nbtabp=1000,nbtabp2=nbtabp*2)

        real*8 rbtab(nbtabp),zbtab(nbtabp),nbtab
        real*8 Cc(np),CUbs(np),Te(np),Cd(np),pres(np)
        real*8 Cu(np)
        real*8 Cu_out(np)
        real*8 pres_s(np)
        real*8 ymu(np),yfp(np)

        real*8 rzbnd(nbtabp2),eqpf(np),eqff(np),fp(np),ipl,ybetpl,yli3,
     *       rtor,btor,rho(np),roc,
     *       g11(np),g22(np),g33(np),                          
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
       
       !voltages
       real*8 vcoils(13)
       real*8 pfc(13)
       integer ncoils
       
       integer NCNB,NCTP
       double precision yvcoilx(12) ,yccoilx(12),t_starts, TIMEC, TAU
       double precision DUMCTP(13), VCOIL(12)
       double precision CTRLM(3), DUMCT(13)
       
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



          call spidat2itm(neql,nteta,nbnd,rzbnd,
     *                    na1,eqpf,eqff,fp,ipl,
     *                    rtor,btor,rho,roc,mu,
     *                    cc,Te,cubs,cd,pres,cu,
     *                    equil_in,parameters_spider )


        nstep=0
        !key_equil=-10
        !key_equil=-2
        key_equil=0

        key_plc=1
        call put_key_plc(key_plc)
        
        key_out=1
        call put_key_out(key_out)


        !parameters_spider%neql = neql
        !parameters_spider%nteta = nteta
        !allocate( equil_in%profiles_1d%jni%value(na1) )
        !allocate( equil_in%profiles_1d%te%value(na1) )
        !allocate( equil_in%profiles_1d%sigmapar%value(na1) )
        !allocate( equil_in%profiles_1d%q(na1) )
        !equil_in%profiles_1d%jni%value = cubs
        !equil_in%profiles_1d%te%value = Te
        !equil_in%profiles_1d%sigmapar%value = cc
        !do i =1,na1
          !equil_in%profiles_1d%q(i) = 1.d0/mu(i)
        !enddo
        
!        call input2spider(nbnd,rzbnd,na1,eqpf,eqff,fp,ipl,
!     *                    rtor,btor,rho,pres,cu,
!     *                    nstep, key_equil,             
!     *                    equil_in,parameters_spider)
     
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
        

        ncoils = 12
        vcoils = 0
        
       parameters_spider%k_grid = 1
       kfixfree = 1
       parameters_spider%k_fixfree = kfixfree
       parameters_spider%key_start = 0
       
       parameters_spider%kpr = 0
       parameters_spider%key_dmf = 0
       parameters_spider%key_0stp = 0
       parameters_spider%key_pres = 0
       parameters_spider%k_auto = 1
       parameters_spider%key_plc = 1
       parameters_spider%nstep = 0
 
       call spider_run(ncoils,vcoils,
     &                  equil_in,equil_out0,parameters_spider)     

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!time !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!! stepping!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        time = parameters_spider%time
        !dt = parameters_spider%dt
        dt = 5e-4
        !dt = 1e-4
        parameters_spider%dt = dt
        
        nstep=0
        
      if(kfixfree.ne.0) then       
        
        NCNB=12
        NCTP=13
        CTRLM = (/0.,3.,101./)
       DUMCT=(/800000.0000, 2.1400, 1.6500, 1.6300,
     &          0.0900, 0.0250, 0.5000, -1.3500,
     &          -1.0500, 0.0000, 1.1200, 0.0000, 0.0000/)
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
        
        call get_pfccur(yccoilx,ncoils)
        !in kA
        yccoilx = yccoilx*1e-3
      endif  
        
        open(111, file='tim_tim.wr', status='replace')            
        close(111)
        open(111, file='pfc_tim.wr', status='replace')            
        close(111)
        open(111, file='ctp_tim.wr', status='replace')           
        close(111)
        open(111, file='vol_tim.wr', status='replace')           
        close(111)
        open(111, file='elo_tim.wr', status='replace')           
        close(111)
        open(111, file='psplex_tim.wr', status='replace')           
        close(111)        
!           the spider parameters can be changed during the run, in particular in addition to
!              if key_dmf=0 - equilibrium with prescribed p' and ff' profiles and toroidal plasma current; 
!              if key_dmf=-10 - equilibrium with prescribed p and <(j,B)>;
!           other equilibrium problems can be solved (with the profiles taken from the intermediate equilibria) 
!              if key_dmf=2 - equilibrium with prescribed p',q and poloidal flux in plasma 
!              if key_dmf=-2 - equilibrium + magnetic field diffusion with f and psi_ext at boundary
!              if key_dmf=-3 - equilibrium + magnetic field diffusion with f at boundary and toroidal plasma current

        parameters_spider%kpr = 0
        parameters_spider%key_dmf = -2
        parameters_spider%epsro = 1.0d-7
        parameters_spider%enels = 5.0d-7
        parameters_spider%k_auto = 0
        !parameters_spider%key_plc = 1
        parameters_spider%key_out = 1
        parameters_spider%key_pres = 0
        parameters_spider%n_dmf =3
        
        do i=1,20
        
              nstep=nstep+1
              time=time+dt
  
!          !first evolution step with fixed profiles (for testing ASTRA behavior)
!              if(nstep .eq. 1) then
!        parameters_spider%key_dmf = -0
!              else
!        parameters_spider%key_dmf = -2
!              endif

        parameters_spider%nstep = nstep
        parameters_spider%time = time
        
!        ! needed for input2spider setup (dummy equil_in)
!        call put_enels(parameters_spider%enels)
!        call put_key_plc(parameters_spider%key_plc)
       
      if(kfixfree.ne.0) then       

             
        if(parameters_spider%k_grid .eq. 0) then
            call get_rmax(DUMCTP(2),DUMCTP(11),DUMCTP(4),DUMCTP(7),
     &                DUMCTP(10))!(r_max,r_min,r_geo,z_max,delta_up)
            call get_zccur(DUMCTP(3),DUMCTP(5),DUMCTP(6))!(rc_cur,zc_cur,z2c_cur)
            !call get_zccur_ada(DUMCTP(3),DUMCTP(5),DUMCTP(6))!(rc_cur,zc_cur,z2c_cur)
            call get_zbran(1.27d0,DUMCTP(9))       
            call get_zbran(1.72d0,DUMCTP(8))
        else if(parameters_spider%k_grid .eq. 1) then
            call get_rmax_ada(DUMCTP(2),DUMCTP(11),DUMCTP(4),DUMCTP(7),
     &                DUMCTP(10))!(r_max,r_min,r_geo,z_max,delta_up)
            call get_zccur_ada(DUMCTP(3),DUMCTP(5),DUMCTP(6))!(rc_cur,zc_cur,z2c_cur)
            call f_get_zbran(1.27d0,DUMCTP(9))
            call f_get_zbran(1.72d0,DUMCTP(8))
       endif
        
        !Zsquad = Zcurr
        DUMCTP(6) = sqrt(DUMCTP(6))
        DUMCTP(6) = DUMCTP(5)
        
        call get_pfccur(pfc,ncoils)
        DUMCTP(12)=pfc(9)
        DUMCTP(13)=pfc(10)
        
        !plasma current from adaptive grid
        if(i .eq. 1) then
            DUMCTP(1)=equil_out0%global_param%i_plasma           
            !for initial equilibrium equil_out0
            open(111, file='tim_tim.wr', access='append')            
            write(111,'(1p12e20.12)') time-dt
            close(111)
            open(111, file='pfc_tim.wr', access='append')            
            write(111,'(1p12e20.12)') pfc(1:12)
            close(111)
            open(111, file='ctp_tim.wr', access='append')            
            write(111,'(1p13e20.12)') DUMCTP(1:13)
            close(111)
            open(111, file='vol_tim.wr', access='append')            
            write(111,'(1p12e20.12)') vcoils(1:12)*0
            close(111)
            open(111, file='elo_tim.wr', access='append')            
            write(111,'(1p12e20.12)') equil_out0%profiles_1d
     &       %elongation(size(equil_out0%profiles_1d%elongation))
            close(111)
            open(111, file='psplex_tim.wr', access='append')            
            write(111,'(1p12e20.12)')
     &       equil_out0%global_param%psplex
!     &       /equil_out0%profiles_1d%q(size(equil_out0%profiles_1d%q))
            close(111)
        else
            DUMCTP(1)=equil_out%global_param%i_plasma
        endif
        
        TIMEC = time
               
        call CTAUG_STD(NCNB,NCTP,CTRLM,
     &       yvcoilx,yccoilx,t_starts,
     &       TIMEC,TAU,DUMCT,DUMCTP,VCOIL)
     
        vcoils(1:12) = VCOIL
        !vcoils(1:12) = 0
        
        !vcoils(4) = 0.
        !vcoils(2:3) = 0.
        !vcoils(1:8) = 0.
        !vcoils=0.
        
       else
        !for initial equilibrium equil_out0
        if(i .eq. 1) then
            open(111, file='tim_tim.wr', access='append')            
            write(111,'(1p12e20.12)') time-dt
            close(111)
            open(111, file='psplex_tim.wr', access='append')            
            write(111,'(1p12e20.12)')
     &       equil_out0%global_param%psplex
!     &       /equil_out0%profiles_1d%q(size(equil_out0%profiles_1d%q))
            close(111)
        endif
       endif         
        
        call spider_run(ncoils,vcoils,
!     &                  equil_in,equil_out,parameters_spider)     
!     &                  equil_out,equil_out,parameters_spider)
     &                  equil_out0,equil_out,parameters_spider)
          
        call wrd_tim
        
      if(kfixfree.ne.0) then       
        open(111, file='tim_tim.wr', access='append')            
        write(111,'(1p12e20.12)') time
        close(111)
        call get_pfccur(pfc,ncoils)
        open(111, file='pfc_tim.wr', access='append')            
        write(111,'(1p12e20.12)') pfc(1:12)
        close(111)
        open(111, file='ctp_tim.wr', access='append')            
        write(111,'(1p13e20.12)') DUMCTP(1:13)
        close(111)
        open(111, file='vol_tim.wr', access='append')            
        write(111,'(1p12e20.12)') vcoils(1:12)
        close(111)
        open(111, file='elo_tim.wr', access='append')            
        write(111,'(1p12e20.12)') equil_out%profiles_1d
     &       %elongation(size(equil_out%profiles_1d%elongation))
        close(111)
            open(111, file='psplex_tim.wr', access='append')            
            write(111,'(1p12e20.12)')
     &       equil_out%global_param%psplex
!     &       /equil_out%profiles_1d%q(size(equil_out%profiles_1d%q))
            close(111)
       else
        open(111, file='tim_tim.wr', access='append')            
        write(111,'(1p12e20.12)') time
        close(111)
            open(111, file='psplex_tim.wr', access='append')            
            write(111,'(1p12e20.12)')
     &       equil_out%global_param%psplex
!     &       /equil_out%profiles_1d%q(size(equil_out%profiles_1d%q))
            close(111)
       endif
            
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
 
 !! extrapolation on magn. axis
       
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
