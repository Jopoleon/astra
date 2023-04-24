      subroutine spider_run(ncoils,vcoils,equil_in,equil_out,params)

      use imas_ids       
      use parameters_a2spider, only: type_parameters, TWOPI
      use vol2d_eq_face
 
      implicit none
 
      type(type_equilibrium) equil_in, equil_out
      type(type_parameters) params
      integer ncoils
      real*8 vcoils(ncoils),startt,fstartt
      integer key_0st,key_prs
	common /com_0st/ key_0st,key_prs
      
      !locals
      integer nstep, key_dmf, k_grid, k_auto, k_fixfree, key_start
      integer kpr, kname,npointzz
      real(8) time, dt, dpsdt
      integer i,neql
      integer itmastra
      real*8 um, psi_boundary,press0,btor
      
!	call cpu_time(startt)
						npointzz=size(equil_in%profiles_1d%psi)
!		write(*,*) 'npointzz ', npointzz
      nstep = params%nstep
      time = params%time
      dt = params%dt
      key_dmf = params%key_dmf
      k_grid = params%k_grid
      k_auto = params%k_auto
      k_fixfree = params%k_fixfree
      dpsdt = params%dpsdt
      key_start = params%key_start
      
      kpr=params%kpr
      call  kpr_calc(kpr)
      kname=params%kname  
      call  put_name(params%prename,kname)

      !set the key_plc for plasma current
      call put_key_plc(params%key_plc)
      !set the key_out for currents and inductive voaltages update in circuit equation
      call put_key_out(params%key_out)
      !set the circuit equation accuracy
      call put_enels(params%enels)
      call put_ndmf(params%n_dmf)
      key_0st = params%key_0stp
      key_prs = params%key_pres
              
      call put_tim(dt,time)
!      call savepsi
      
      itmastra = 0
      if(associated(equil_in%profiles_1d%psi)) then
						 itmastra = 1
      endif
!	write(*,*) itmastra						
      if(itmastra .eq. 1) then
!        write(*,*) 'spider_run: CPO input'
        call itm2spider(equil_in,params)
       endif
!	call cpu_time(fstartt)
!	write(*,*) 'itm2spid ',fstartt-startt       	
!	write(*,*) 'call spider'
!	call cpu_time(startt)
      call spider(nstep,time,dt,key_dmf,k_grid,k_auto,k_fixfree,
     *              dpsdt,key_start,vcoils)  
!	call cpu_time(fstartt)
!	write(*,*) 'spider in interfac ',fstartt-startt       	
!	write(*,*) 'call spider out'
	!call cpu_time(startt)
        !call  pla_volt(dt)
!	write(*,*) 'curavg'
        call cur_avg
        !call  wrb
!	write(*,*) 'wrspik'
        call wr_spik

        !equil_out
        if(k_fixfree .eq. 0) then
          psi_boundary = 0.d0
        else if(k_grid .eq. 0) then
          call get_umup(um,psi_boundary)
        else if(k_grid .eq. 1) then
          call f_get_umup(um,psi_boundary)
        endif
		if(associated(equil_in%profiles_1d%pressure)) then
						press0=equil_in%profiles_1d%pressure(npointzz)
	  endif
								btor=equil_in%global_param%toroid_field%b0

!	write(*,*) 'geteq'
        call get_eq(equil_out,psi_boundary,press0,btor)
        
        ! generic flux surface averaging
!        if(itmastra .eq. 1) then
!	write(*,*) 'puteq'
            call put_eq(equil_out,params)
!            !check averaging
!            allocate(equil_out%profiles_1d%jparallel( params%neql ))
!            call vol2d_eq(equil_out, equil_out%coord_sys%bcell, 
!     &             equil_out%profiles_1d%jparallel)
!            equil_out%profiles_1d%jparallel = 
!     &    equil_out%profiles_1d%jparallel
!     &    /equil_out%global_param%toroid_field%b0
!        endif
	!call cpu_time(fstartt)
	!write(*,*) 'spid2itm ',fstartt-startt       	

      return
      end





      subroutine spider_run_2(ncoils,vcoils,params)

      use imas_ids       
      use parameters_a2spider, only: type_parameters, TWOPI
      use vol2d_eq_face
 
      implicit none
 
      type(type_equilibrium) equil_in, equil_out
      type(type_parameters) params
      integer ncoils
      real*8 vcoils(ncoils)
      integer key_0st,key_prs
	common /com_0st/ key_0st,key_prs
      
      !locals
      integer nstep, key_dmf, k_grid, k_auto, k_fixfree, key_start
      integer kpr, kname,npointzz
      real(8) time, dt, dpsdt
      integer i,neql
      integer itmastra
      real*8 um, psi_boundary,press0,btor
      
      nstep = params%nstep
      time = params%time
      dt = params%dt
      key_dmf = params%key_dmf
      k_grid = params%k_grid
      k_auto = params%k_auto
      k_fixfree = params%k_fixfree
      dpsdt = params%dpsdt
      key_start = params%key_start
      
      kpr=params%kpr
      call  kpr_calc(kpr)
      kname=params%kname  
      call  put_name(params%prename,kname)

      !set the key_plc for plasma current
      call put_key_plc(params%key_plc)
      !set the key_out for currents and inductive voaltages update in circuit equation
      call put_key_out(params%key_out)
      !set the circuit equation accuracy
      call put_enels(params%enels)
      call put_ndmf(params%n_dmf)
      key_0st = params%key_0stp
      key_prs = params%key_pres
              
      call put_tim(dt,time)
!      call savepsi
      
						 itmastra = 1
						
      if(itmastra .eq. 1) then
        call itm2spider_2(params)
       endif
!	write(*,*) 'call spider'

      call spider(nstep,time,dt,key_dmf,k_grid,k_auto,k_fixfree,
     *              dpsdt,key_start,vcoils)  
       
!	write(*,*) 'call spider out'
        !call  pla_volt(dt)
        call cur_avg
        !call  wrb
        call wr_spik

        !equil_out

      return
      end
      
      subroutine itm2spider_2(params)
      use imas_ids       
      use parameters_a2spider, only: type_parameters, TWOPI
      use durs_d_modul       
      use ppf_modul
      use bnd_modul
      
      implicit none
      type(type_equilibrium) equil_in
      type(type_parameters) params
      
      integer na1, nbnd, neql, nteta, nstep, key_ini
      real*8 ipl, rtor, btor
      
      !include 'double.inc'
	integer nrp, ntp, nr1p, nt1p, nr2p, nt2p
	include 'dimpl1.inc'
      integer nbtabp,nursp,nursp4,nursp6
      parameter(nbtabp=1000)
      parameter(nursp=1000,nursp4=nursp+4,nursp6=nursp4*6)
      integer key_0st,key_prs
	common /com_0st/ key_0st,key_prs
	integer i, ib
	real*8 ps0, psb, pscale, fscale
	real*8  rbmax, rbmin, zbmax, zbmin, rc0, zc0
	character*40   eqdfn
      
      !set the key_plc for plasma current
      call put_key_plc(params%key_plc)
      !set the key_out for currents and inductive voaltages update in circuit equation
      call put_key_out(params%key_out)
      !set the circuit equation accuracy
      call put_enels(params%enels)
              
!! input dimension by the size of allocated psi
      neql = params%neql
      nteta = params%nteta
      nstep = params%nstep
      key_ini = params%key_ini
      key_0st = params%key_0stp
      key_prs = params%key_pres
      i_eqdsk= params%i_eqdsk
      eqdfn = params%eqdfn
      
      igdf=2
      nurs=-399
      i_betp=0
      keyctr=0 
      i_betp=0

      ! fixed boundary equilibrium accuracy
      !epsro=1.d-6
      epsro=params%epsro
      
      call put_Ipl(0.d0)
      
      call aspid_flag(1)

  
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine itm2spider(equil_in,params)
      use imas_ids       
      use parameters_a2spider, only: type_parameters, TWOPI
      use durs_d_modul       
      use ppf_modul
      use bnd_modul
      
      implicit none
      type(type_equilibrium) equil_in
      type(type_parameters) params
      
      integer na1, nbnd, neql, nteta, nstep, key_ini
      real*8 ipl, rtor, btor
      
      !include 'double.inc'
	integer nrp, ntp, nr1p, nt1p, nr2p, nt2p
	include 'dimpl1.inc'
      integer nbtabp,nursp,nursp4,nursp6
      parameter(nbtabp=1000)
      parameter(nursp=1000,nursp4=nursp+4,nursp6=nursp4*6)
      integer key_0st,key_prs
	common /com_0st/ key_0st,key_prs
	integer i, ib
	real*8 ps0, psb, pscale, fscale
	real*8  rbmax, rbmin, zbmax, zbmin, rc0, zc0
	character*40   eqdfn
      
      if(.not.associated(equil_in%profiles_1d%psi)) then
        write(*,*) 'itm2astra: INVALID input equilibrium CPO; return'
        pause
        return
      endif
       
      !set the key_plc for plasma current
      call put_key_plc(params%key_plc)
      !set the key_out for currents and inductive voaltages update in circuit equation
      call put_key_out(params%key_out)
      !set the circuit equation accuracy
      call put_enels(params%enels)
              
!! input dimension by the size of allocated psi
      na1 = size(equil_in%profiles_1d%psi)
      nbnd = equil_in%eqgeometry%boundary%npoints
      neql = params%neql
      nteta = params%nteta
      nstep = params%nstep
      key_ini = params%key_ini
      key_0st = params%key_0stp
      key_prs = params%key_pres
      i_eqdsk= params%i_eqdsk
      eqdfn = params%eqdfn
      
!! check dimensions
      if(nbnd.gt.nbtabp) then
        write(*,*) 'spider:nbnd exeeds maximum value',nbtabp		   	
        write(*,*) 'nbnd=',nbnd		   				
        write(*,*) 'program is interrupted'
        stop					   				
      endif				
      if(na1.gt.nursp) then
        write(*,*) 'spider:na1 exeeds maximum value',nursp		   			
        write(*,*) 'na1=',na1		   				
        write(*,*) 'program is interrupted'
        stop					   				
      endif  
      if(neql.gt.nrp .or. nteta.gt.ntp) then
        write(*,*) 'spider: neql or nteta exeeds maximum value'		
        write(*,*) 'neql=',neql,'nrp=',nrp		   				
        write(*,*) 'nteta=',nteta,'ntp=',ntp		   				
        write(*,*) 'program is interrupted'
        stop					   				
       endif
       
!! toroidal field, radius, plasma current
      rtor = equil_in%global_param%toroid_field%r0
      btor = equil_in%global_param%toroid_field%b0
      ipl = equil_in%global_param%i_plasma*1e-6
C      write(*,*), na1, rtor,neql,ipl,btor,nbnd
!! p' and ff' profiles
      nutab = na1     
      if(allocated(pstab))
     %deallocate( pstab, pptab, fptab )
      allocate( pstab(na1), pptab(na1), fptab(na1) )
!        write(fname,'(a,a)') path(1:kname),'tabppf.dat'
!        open(1,file=fname,form='formatted')
!          !open(1,file='tabppf.dat')
!              write(1,*) nutab
!           do i=1,nutab
!              write(1,*) pstab(i),pptab(i),fptab(i)
!           enddo
!          close(1)
      
      ps0=equil_in%profiles_1d%psi(1)
      psb=equil_in%profiles_1d%psi(na1)
      pstab(1)=0.d0
      do i=2,nutab
        pstab(i)=(equil_in%profiles_1d%psi(i)-ps0)/(psb-ps0)
	  if(pstab(i).le.pstab(i-1))   then
          write(*,*) 'psi is nonmonotonic!!! ', pstab(i),i
        endif
      enddo   
      do i=1,nutab
        pptab(i)=equil_in%profiles_1d%pprime(i)
        fptab(i)=equil_in%profiles_1d%ffprime(i)
      enddo
      !normalization
      pscale = 2 * TWOPI * 1.d-7 * TWOPI
      fscale = TWOPI
      do i=1,nutab
        pptab(i)=-pptab(i)*pscale
        fptab(i)=-fptab(i)*fscale
      enddo

!! boundary
      nbtab = nbnd            
      if(allocated(rbtab))
     %deallocate( rbtab, zbtab )
      allocate( rbtab(nbtab), zbtab(nbtab) )
!        write(fname,'(a,a)') path(1:kname),'tab_bnd.dat'
!        open(1,file=fname,form='formatted')
!         !open(1,file='tab_bnd.dat')
!	     write(1,*) nbtab
!	     do ib=1,nbtab
!	        write(*,*) rbtab(ib),zbtab(ib)
!	     enddo
!         close(1)
	do ib=1,nbtab
	  rbtab(ib)=equil_in%eqgeometry%boundary%r(ib)
        zbtab(ib)=equil_in%eqgeometry%boundary%z(ib)
	enddo
      rbmax=rbtab(1)
      rbmin=rbtab(1)
      zbmax=zbtab(1)
      zbmin=zbtab(1)
      do ib=1,nbtab
        if(rbtab(ib).ge.rbmax) rbmax=rbtab(ib)
        if(rbtab(ib).le.rbmin) rbmin=rbtab(ib)
        if(zbtab(ib).ge.zbmax) zbmax=zbtab(ib)
        if(zbtab(ib).le.zbmin) zbmin=zbtab(ib)
	enddo
      rc0=0.5d0*(rbmax+rbmin)           
      zc0=0.5d0*(zbmax+zbmin)           
   
      !if(nstep.ne.0) then
      if(associated(equil_in%profiles_1d%rho_tor)) then
       call S_profro_p(equil_in,params)
       call S_profro_sbc(equil_in,params)
       if(associated(equil_in%profiles_1d%pressure)) then
         call S_proffi_pres(equil_in,params)
       endif
      endif

                                                           
!! parameters durs_d	          
!        open(1,file=fname,form='formatted')
!          !open(1,file='durs_d.dat')
!              write(1,*) n_tht,n_psi,igdf,nurs,keyctr,i_eqdsk,i_betp
!              write(1,*) epsro,betplx,tokf,psax,b0,r0,rax,zax
!              write(1,*) alf0,alf1,alf2,bet0,bet1,bet2
!          close(1)
      n_tht=nteta
      n_psi=neql
      igdf=2
      nurs=-399
      i_betp=0
      keyctr=0 
      i_betp=0

      ! fixed boundary equilibrium accuracy
      !epsro=1.d-6
      epsro=params%epsro
      
      tokf=ipl
      b0 = btor
      r0 = rtor
      rax=rc0
      zax=zc0
      call put_Ipl(tokf)
      
      if(nstep.eq.0 .AnD. key_ini.eq.0) then
          if(i_betp.eq.0) then
             betplx = 0.d0
          endif
          psax = ps0
          if(keyctr.ne.2) then
             psax = 1.d0
           endif
          if(nurs.lt.0) then
            alf0 = 0.d0
            alf1 = 0.d0
            alf2 = 0.d0
            bet0 = 0.d0
            bet1 = 0.d0
            bet2 = 0.d0
          endif
      endif
      
      ! read eqdsk instead
      if(nstep.eq.0 .AnD. key_ini.eq.0) then
        i_eqdsk=1
        call tab_efit(tokf,psax,eqdfn,rax,zax,b0,r0)
        keyctr=0 
        key_0st=0
      endif
						
					 !         if(params%key_dmf.eq.-3) then
	  		 !         call profro_cu(na1,neql,rtor,btor,cu,rho,roc)
    		!         call profro_sbc(na1,neql,rtor,btor,cc,Te,cubs,cd,rho,roc)
			   !       endif
 
         
      if(keyctr .eq. 0) call taburs(0,1.d0,nurs)
      
      call aspid_flag(1)

  
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         subroutine put_key_plc(k)
          common /com_kplc/ key_plc
          integer k,key_plc
            key_plc=k
         return
         end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         subroutine put_key_out(k)
          common /com_kout/ key_out
          integer k,key_out
            key_out=k
         return
         end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         subroutine put_enels(e)
          common /com_enels/ ENELS
          real*8 e,ENELS
            ENELS=e
         return
         end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         subroutine put_ndmf(k_dmf)
          common /com_ndmf/ n_dmf
          integer n_dmf,k_dmf
            n_dmf=k_dmf
         return
         end

        subroutine astra2spider(neql,nteta,nbnd,rzbnd,key_dmf,
     *                    na1,yeqpf,yeqff,fp,ipl,
     *                    rtor,btor,rho,roc,nstep,yreler,mu,
     *                    cc,Te,cubs,cd,key_ini,eqdfn,pres,cu,
     *                    key_0stp,key_pres                     )
                               

         use durs_d_modul       
         use ppf_modul       
         use bnd_modul       

          include 'double.inc'
	    include 'dimpl1.inc'
          parameter(nbtabp=1000)
          parameter(nursp=1000,nursp4=nursp+4,nursp6=nursp4*6)
!          common /com_bas/ rbtab(nbtabp),zbtab(nbtabp),nbtab
!          common/com_pas/ pstab(nursp),pptab(nursp),fptab(nursp),nutab
	    common /creler/ relerr
	    common /com_0st/ key_0st,key_prs

          real*8 rout(nrp,ntp),zout(nrp,ntp)

          real   t_start, t_finish

        real*8 rzbnd(*),yeqpf(*),yeqff(*),fp(*),ipl,ybetpl,yli3,
     *         rtor,btor,rho(*),roc,yreler,mu(*),
     *         cc(*),Te(*),cubs(*),cd(*),pres(*),cu(*) 
                              
          character*40 eqdfn

      real(8), allocatable :: eqpf(:), eqff(:)
     
      allocate( eqpf(na1), eqff(na1) )
      
        pi=3.14159265359d0
        amu0=0.4d0*pi

           do i=1,na1
             eqpf(i)=yeqpf(i)* amu0
             eqff(i)=yeqff(i)* amu0
           enddo


		relerr	=yreler
 				!write(*,*) 'b_eqb, rz', (rzbnd(j),j=1,nbnd)
 				!write(*,*) 'nstep', nstep
             if(nbnd.gt.nbtabp) then
              write(*,*) 'spider:nbnd exeeds maximum value',nbtabp		   	
              write(*,*) 'nbnd=',nbnd		   				
              write(*,*) 'program is interrupted'
              stop					   				
             endif				
				
             if(na1.gt.nursp) then
              write(*,*) 'spider:na1 exeeds maximum value',nursp		   			
              write(*,*) 'na1=',na1		   				
              write(*,*) 'program is interrupted'
              stop					   				
             endif
		   
             if( neql.gt.nrp .or. nteta.gt.ntp) then
              write(*,*) 'spider: neql or nteta exeeds maximum value'		
              write(*,*) 'neql=',neql,'nrp=',nrp		   				
              write(*,*) 'nteta=',nteta,'ntp=',ntp		   				
              write(*,*) 'program is interrupted'
              stop					   				
             endif
				          
             key_0st=key_0stp
             key_prs=key_pres
             nurs=-399
             i_betp=0
             i_bsh=1
             igdf=2
             i_eqdsk=0
             epsro=1.d-6
             maxit=250
            ! nstep=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             n_tht=nteta
             n_psi=neql
             platok=ipl
             tokf=ipl
C	write(*,*) 'platok in a2spid: ',tokf,ipl,neql,na1
           call put_Ipl(tokf)

        if(nstep.eq.0 .AnD. key_ini.eq.0) then
               call tab_efit(tokf,psax,eqdfn,rax,zax,b0,r0)
             keyctr=0 
             key_0st=0
!        write(fname,'(a,a)') path(1:kname),'durs_d.dat'
        if(i_betp.eq.0) then
            betplx = 0.d0
        endif
        if(keyctr.ne.2) then
            psax = 1.d0
        endif
        if(nurs.lt.0) then
            alf0 = 0.d0
            alf1 = 0.d0
            alf2 = 0.d0
            bet0 = 0.d0
            bet1 = 0.d0
            bet2 = 0.d0
        endif
!        open(1,file=fname,form='formatted')
!          !open(1,file='durs_d.dat')
!              write(1,*) n_tht,n_psi,igdf,nurs,keyctr,i_eqdsk,i_betp
!              write(1,*) epsro,betplx,tokf,psax,b0,r0,rax,zax
!              write(1,*) alf0,alf1,alf2,bet0,bet1,bet2
!          close(1)
         return
        endif

         nutab=na1+1

      if(allocated(pstab))
     %deallocate( pstab, pptab, fptab )
      allocate( pstab(na1+1), pptab(na1+1), fptab(na1+1) )
     
         nbtab=nbnd

      if(allocated(rbtab))
     %deallocate( rbtab, zbtab )
      allocate( rbtab(nbtab), zbtab(nbtab) )

	     do ib=1,nbtab

	        rbtab(ib)=rzbnd(ib)
              zbtab(ib)=rzbnd(ib+nbnd)

	     enddo


!        write(fname,'(a,a)') path(1:kname),'tab_bnd.dat'
!        open(1,file=fname,form='formatted')
!         !open(1,file='tab_bnd.dat')
!	     write(1,*) nbtab
!	     do ib=1,nbtab
!	        write(1,*) rbtab(ib),zbtab(ib)
!	     enddo
!         close(1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!__CH

             rbmax=rbtab(1)
             rbmin=rbtab(1)
             zbmax=zbtab(1)
             zbmin=zbtab(1)

	     do ib=1,nbtab

            if(rbtab(ib).ge.rbmax) rbmax=rbtab(ib)
            if(rbtab(ib).le.rbmin) rbmin=rbtab(ib)
            if(zbtab(ib).ge.zbmax) zbmax=zbtab(ib)
            if(zbtab(ib).le.zbmin) zbmin=zbtab(ib)

	     enddo

          rc0=0.5d0*(rbmax+rbmin)           
          zc0=0.5d0*(zbmax+zbmin)           

          rax=rc0
          zax=zc0
!__CH
          b0=btor
          r0=rtor
!__CH


!! p' and ff' profiles

        if(nstep.eq.0) then

          if(key_0st.eq.1) then
         
              key_prs=1
              do i=1,na1
                eqpf(i)=0.d0
                eqff(i)=cu(i)
              enddo
                tokf=1.d0
                
              do i=1,na1
                eqpf(i)=eqpf(i) !*amu0
                eqff(i)=eqff(i) !*amu0
              enddo
              
                !nutab=na1+1
                pstab(1)=0.d0
              do i=2,nutab
                pstab(i)=(rho(i-1)/rho(na1))**2
              enddo
              
            call profro_cu(na1,neql,rtor,btor,cu,rho,roc)
            call proffi_pres(na1,neql,pres,rho,roc)

          elseif(key_0st.eq.0) then
          
              do i=1,na1
                eqpf(i)=eqpf(i) !*amu0
                eqff(i)=eqff(i) !*amu0
              enddo
!! exstrapolation psi on magn. axis
       
             rh0=0.d0
		              
             rh1=rho(1)           
             rh2=rho(2)           
             rh3=rho(3)           

             ps1=fp(1)           
             ps2=fp(2)           
             ps3=fp(3)           


       !call EXTRP2(rh0,ps0, rh1,rh2,rh3, ps1,ps2,ps3)

           ps0=(ps1*rh2**2-ps2*rh1**2)/(rh2**2-rh1**2)


           psb=fp(na1)

                pstab(1)=0.d0


        !nutab=na1+1

              do i=2,nutab
                pstab(i)=(fp(i-1)-ps0)/(psb-ps0)
	       if(pstab(i).le.pstab(i-1))   then
          write(*,*) 'psi is nonmonotonic!!! ', pstab(i),i
             endif
              enddo


          endif  !(key_0st

!! exstrapolation on magn. axis
       
             rh0=0.d0
		              
             rh1=rho(1)           
             rh2=rho(2)           
             rh3=rho(3)           

             ppx1=eqpf(1)
             ppx2=eqpf(2)
             ppx3=eqpf(3)

             ffx1=eqff(1)
             ffx2=eqff(2)
             ffx3=eqff(3)

       call EXTRP2(rh0,ppx0, rh1,rh2,rh3, ppx1,ppx2,ppx3)
       call EXTRP2(rh0,ffx0, rh1,rh2,rh3, ffx1,ffx2,ffx3)

                !pptab(1)=ppx0
                !fptab(1)=ffx0
                pptab(1)=ppx1
                fptab(1)=ffx1

        !nutab=na1+1

              do i=2,nutab
                pptab(i)=eqpf(i-1)
                fptab(i)=eqff(i-1)
              enddo

!normalization        

              do i=1,nutab
                pptab(i)=pptab(i)/rtor
                fptab(i)=fptab(i)*rtor
              enddo
              
!        write(fname,'(a,a)') path(1:kname),'tabppf.dat'
!        open(1,file=fname,form='formatted')
!          !open(1,file='tabppf.dat')

!              write(1,*) nutab


!          close(1)
                                                           
        endif   ! (nstep.eq.0

!__CH

       if(nstep.ne.0) then
        
            !call profro_pres_L(na1,neql,pres,rho,roc)
            !call profro_pres(na1,neql,pres,rho,roc)
            !call proffi_pres(na1,neql,pres,rho,roc)
        
          if(key_dmf.ne.0) then
           
            call profro_cu(na1,neql,rtor,btor,cu,rho,roc)
            
            if(key_prs.eq.0) then
             call profro_p(na1,neql,rtor,eqpf,rho,roc) !eqpf \times mu0
            else 
             call proffi_pres(na1,neql,pres,rho,roc)
            endif
           
            if(key_dmf.eq.-2 .OR. key_dmf.eq.-3) then
             call profro_sbc(na1,neql,rtor,btor,cc,Te,cubs,cd,rho,roc)
            endif
            
          else   !(key_dmf=0)
          
           call profro(na1,neql,rtor,eqpf,eqff,rho,roc)
           !call retab_L_ast(pstab,pptab,fptab,nutab) !
           call retab_L
            
          endif
 
       endif
!__CH
              !do i=1,nutab
              !  pptab(i)=eqpf(i)
              !  fptab(i)=eqff(i)
              !  pstab(i)=fp(i)
              !enddo

          !open(17,file='out.pr')



 1149     format(a40)

        if(nstep.eq.0) then
             keyctr=0 
!        write(fname,'(a,a)') path(1:kname),'durs_d.dat'
        if(i_betp.eq.0) then
            betplx = 0.d0
        endif
        if(keyctr.ne.2) then
            psax = 1.d0
        endif
        if(nurs.lt.0) then
            alf0 = 0.d0
            alf1 = 0.d0
            alf2 = 0.d0
            bet0 = 0.d0
            bet1 = 0.d0
            bet2 = 0.d0
        endif
!        open(1,file=fname,form='formatted')
!          !open(1,file='durs_d.dat')
!              write(1,*) n_tht,n_psi,igdf,nurs,keyctr,i_eqdsk,i_betp
!              write(1,*) epsro,betplx,tokf,psax,b0,r0,rax,zax
!              write(1,*) alf0,alf1,alf2,bet0,bet1,bet2
!          close(1)
        endif


             ! do i=1,na1
             !   eqpf(i)=eqpf(i)/amu0
             !   eqff(i)=eqff(i)/amu0
             ! enddo

         return
         end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine spider2astra(rout,zout,rtor,btor,rho,roc,na1,
     *                    g11,g22,g33,vr,vrs,slat,gradro,rocnew,
     *                    mu,ipol,bmaxt,bmint,bdb02,b0db2,bdb0,droda,
     *				    yreler,yli3,ni_p,nj_p,platok,cu,fp,pres,W_Dj,
     *                  yFOFB )
                               

         use ppf_modul       
         use bnd_modul       

          include 'double.inc'
	    include 'dimpl1.inc'
          parameter(nbtabp=1000)
          parameter(nursp=1000,nursp4=nursp+4,nursp6=nursp4*6)
!          common /com_bas/ rbtab(nbtabp),zbtab(nbtabp),nbtab
!          common/com_pas/ pstab(nursp),pptab(nursp),fptab(nursp),nutab
      	common /creler/ relerr
	integer na1,ni_p,nj_p
          real*8 rout(ni_p,nj_p),zout(ni_p,nj_p)

        real*8 
     *       g11(*),g22(*),g33(*),                          
     *       vr(*),vrs(*),slat(*),gradro(*),rocnew,ybetpl,betpol,platok,
     *       mu(*),ipol(*),bmaxt(*),bmint(*),bdb02(*),bdb0(*),
     *	     b0db2(*),droda(*),rtor,btor,rho(*),roc,cu(*),fp(*),pres(*),
     *	     W_Dj(*),yFOFB(*)
                              
          character*40 eqdfn


         call get_rz(rout,zout,ni_p,nj_p)

cw	write(*,*) 'after eqb'

             ybetpl=betpol

           call comet(na1,platok,
     *                    rtor,btor,rho,roc,nstep,
     *                    g11,g22,g33,vr,vrs,slat,gradro,rocnew,
     *                    mu,ipol,bmaxt,bmint,bdb02,b0db2,bdb0,droda,
     *                    ybetpl,yli3,cu,fp,pres,W_Dj,yFOFB)


!       deallocate( pstab, pptab, fptab )
!       deallocate( rbtab, zbtab )

       return
       end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!last change: 27.05.03
 
           subroutine comet(na1,platok,
     *                      rtor,btor,rho,roc,nstep,
     *                      g11,g22,g33,vr,vrs,slat,gradro,rocnew,
     *                      mu,ipol,bmaxt,bmint,bdb02,b0db2,bdb0,droda,
     *					  ybetpl,yli3,cu,fp,pres,W_Dj,yFOFB)

         include 'double.inc'
         include 'dim.inc'
         parameter(nrpl=nrp+1) 
          parameter(nrpl4=nrpl+4,nrpl6=nrpl4*6)
         include 'compol.inc'
       common /combsh/ rm0,zm0,rc0,zc0,asp0,el_up,el_lw,tr_up,tr_lw,nbsh
        common /com_jb/ BJ_av(nrp),curfi_av(nrp)
        common/com_heat_Dj/ WDj(nrp)
        common /com_trap/ trap(nrp)

          real*8 rhos(nrp),vols(nrpl),g11s(nrpl),g22s(nrpl),g33s(nrpl)
          real*8 gradrs(nrpl),amus(nrpl),fnors(nrpl)
          real*8 dvdro(nrpl),sa(nrp),rhocs(nrpl),rhocn(nrpl),rhosn(nrpl)
          real*8 b_maxt(nrpl),b_mint(nrpl),b_db02(nrpl),b_db0(nrpl),
     *		     b_0db2(nrpl),d_roda(nrpl)
          real*8 psi_1(nrp)
          real*8 prs_1(nrp)
     	
          real*8 RRK(nrpl4),CCK(nrpl4),WRK(nrpl6)
          real*8 CWK(4)
          !real*8 rhowr(1000)
          real(8), allocatable :: rhowr(:), rhowrh(:)
         
        real*8 rtor,btor,rho(*),roc,ybetpl,yli3,
     *       g11(*),g22(*),g33(*),                          
     *       vr(*),vrs(*),slat(*),gradro(*),rocnew,droda(*), 
     *       mu(*),ipol(*),bmaxt(*),bmint(*),bdb02(*),b0db2(*),bdb0(*),
     *       cu(*),fp(*),pres(*),W_Dj(*)
     
          real*8 yFOFB(*)
          real*8 traps(nrpl)
         dimension Btot(nrp,ntp)
          
          !common /fp_sav/ fp_0(1000)
          !common /fp_dot/ dfpdt(1000),nna1
                          
           sqrt(xx)=dsqrt(xx)
       
         !do i=1,na1  
          !fp_0(i)=fp(i)
         !enddo
         !nna1=na1
       
         allocate( rhowr(na1), rhowrh(na1) )
       
         na=na1-1
         platok=tokp       
         do i=1,iplas
          rhos(i)=sqrt(flx_fi(i)/(pi*btor))
         enddo
         
         call get_psibon(psi_bn1)
         
         do i=1,iplas
          psi_1(i)=(psi(i,2)+psi_bn1)*2.d0*pi
         enddo
         
         do i=1,iplas
          psn=psia(i)
          prs_1(i)=funppp(psn)*(psim-psip)*cnor/amu0
         enddo
CP
cp
         rhocs(1)=0.d0
          rhos(1)=0.d0

         do i=2,iplas
          rhocs(i)=0.5d0*(rhos(i-1)+rhos(i))
         enddo

          rocnew=rhos(iplas)
          rhocs(iplas+1)=rhos(iplas)

c=============================andrey=02.03.2003vvvvvvvvvvvvv
          ARC0=rhocs(iplas+1)
          ARC1=rhocs(iplas)
          ARC2=rhocs(iplas-1)
          ARC3=rhocs(iplas-2)
          SPW1=(ARC0-ARC2)*(ARC0-ARC3)
     &        /(ARC1-ARC2)/(ARC1-ARC3)
          SPW2=(ARC0-ARC1)*(ARC0-ARC3)
     &        /(ARC2-ARC1)/(ARC2-ARC3)
          SPW3=(ARC0-ARC1)*(ARC0-ARC2)
     &        /(ARC3-ARC1)/(ARC3-ARC2)
          q(iplas)=SPW1*q(iplas-1)+SPW2*q(iplas-2)+SPW3*q(iplas-3)
!          trap(iplas)=SPW1*trap(iplas-1)+SPW2*trap(iplas-2)+
!     &    SPW3*trap(iplas-3)

c=============================andrey=02.03.2003^^^^^^^^^^^^^^
 


         do i=1,iplas+1
          rhocn(i)=rhocs(i)/rhos(iplas)
         enddo
          rhocn(1)=0.d0

         do i=1,iplas
          rhosn(i)=rhos(i)/rhos(iplas)
	    droda(i)=rhosn(i)
         enddo
          rhosn(1)=0.d0
cp

          rhocn(iplas+1)=1.d0
          rhosn(iplas)=1.d0

         do i=1,na
          rhowr(i)=rho(i)/roc
         enddo
          rhowr(na1)=1.d0-1.d-9
          !assuming main astra grid: rho(j)=h*(j-0.5) j=1,2,...,na1-1
          ! intermediate astra grid: rho_h(j)=h*j
         do i=1,na
          rhowrh(i)=(rhowr(i)+rhowr(i+1))*0.5d0
         enddo
          rhowrh(na)=rhowrh(na-1)+(rhowr(2)-rhowr(1))
          rhowrh(na1)=1.d0-1.d-9

          sa(1)=0.d0
cw	write(*,*) (rhocn(j),j=1,iplas)
cw	write(*,*) (rhosn(j),j=1,iplas)

         do i=1,iplas-1

           dsqi=0.d0
           dvoli=0.d0
           av_r2=0.d0
           av_gr1=0.d0
           av_gr2=0.d0
           av_grr2=0.d0

          u1=rhos(i)
          u2=rhos(i+1)
          u3=rhos(i+1)
          u4=rhos(i)

          do j=2,nt1

           dvoli=dvoli+vol(i,j)
           dsqi=dsqi+sr(i+1,j)

          r1=r(i,j)
          r2=r(i+1,j)
          r3=r(i+1,j+1)
          r4=r(i,j+1)
          r0=(r1+r2+r3+r4)*0.25d0

          z1=z(i,j)
          z2=z(i+1,j)
          z3=z(i+1,j+1)
          z4=z(i,j+1)

           sss=s(i,j)
           
           dadr= 0.5d0*((u1+u2)*(z2-z1)+(u2+u3)*(z3-z2)+
     *                   (u3+u4)*(z4-z3)+(u4+u1)*(z1-z4))/sss

           dadz=-0.5d0*((u1+u2)*(r2-r1)+(u2+u3)*(r3-r2)+
     *                   (u3+u4)*(r4-r3)+(u4+u1)*(r1-r4))/sss

           gr2=dadr**2+dadz**2

           av_gr1=av_gr1+vol(i,j)*sqrt(gr2)
           av_gr2=av_gr2+vol(i,j)*gr2
           av_grr2=av_grr2+vol(i,j)*gr2/r0**2
           av_r2=av_r2+vol(i,j)/r0**2

          enddo
          sa(i+1)=dsqi*2.d0*pi
c          vols(i+1)=dvoli
          gradrs(i+1)=av_gr1/dvoli
          amus(i+1)=(2.d0*pi)/q(i)
          fnors(i+1)=f(i)/(rtor*btor)

          dvdro(i+1)=2.d0*pi*dvoli/(rhos(i+1)-rhos(i))
cc          g11s(i+1)=av_gr2/dvoli
cc          g22s(i+1)=av_grr2/(dvoli*4.d0*pi*pi)*rtor
          g11s(i+1)=av_gr2*dvdro(i+1)/dvoli
        g22s(i+1)=av_grr2*dvdro(i+1)/(dvoli*4.d0*pi*pi)*rtor/fnors(i+1)
          g33s(i+1)=av_r2*rtor**2/dvoli
         enddo

cccc      magnetic field

CW	write(*,*) 'before b'

           bp2av=0.d0
		   bp2=0.d0
cc          do i=2,iplas-1 
          do i=1,iplas-1 

           dvoli=0.d0
           bmx=0.d0
           bmn=1.d9*btor
           bavr=0.d0
           bavr2=0.d0
           bavrd2=0.d0

           do j=2,nt1

           dvoli=dvoli+vol(i,j)

          r1=r(i,j)
          r2=r(i+1,j)
          r3=r(i+1,j+1)
          r4=r(i,j+1)
          r0=(r1+r2+r3+r4)*0.25d0

            bm_pol=psim*(psia(i+1)-psia(i))/st(i,j)
            bp_pol=psim*(psia(i+1)-psia(i))/st(i,j+1)

         if(i.ne.1) then
        bpol2v=bm_pol**2*( vol1(i,j)/sin1(i,j)+vol2(i,j)/sin2(i,j))+ 
     +         bp_pol**2*( vol3(i,j)/sin3(i,j)+vol4(i,j)/sin4(i,j)) 
         else
        bpol2v=bm_pol**2*(                     vol2(i,j)/sin2(i,j))+ 
     +         bp_pol**2*( vol3(i,j)/sin3(i,j)                    ) 
         endif

            bpol2=bpol2v /vol(i,j)
            b_fi2=(f(i)/r0)**2

            b_tot2=bpol2+b_fi2
            b_tot=dsqrt(b_tot2)
!!   aai 06/02/10{{{{{{{{{{{
           Btot(i,j)=b_tot
!!   }}}}}}}}}}}

            bmx=dmax1(b_tot,bmx)
            bmn=dmin1(b_tot,bmn)

           bavr2=bavr2+b_tot2*vol(i,j)
           bavr =bavr +b_tot *vol(i,j)
	   bavrd2=bavrd2+vol(i,j)/(b_tot2+1.d-8)
           bp2av=bp2av+bpol2v
		bp2=bp2+bpol2v
          enddo            
           b_maxt(i+1)=bmx
           b_mint(i+1)=bmn
           b_db02(i+1)=bavr2/(dvoli*btor**2)
           b_0db2(i+1)=(bavrd2*btor**2)/(dvoli)
           b_db0(i+1)=bavr/(dvoli*btor)


        enddo 
           b_maxt(1)=.5d0*(b_maxt(2)+b_mint(2))
           b_mint(1)=b_maxt(1)
           b_db02(1)=b_db02(2)
           b_0db2(1)=b_0db2(2)
           b_db0(1)=b_db0(2)


c            yli3  =4.d0*pi*bp2av/(rtor*tokp*tokp)
            yli3  =  4.d0*pi*bp2/(rtor*(.4d0*pi*tokp)**2)     



!!   aai 06/02/10{{{{{{{{{{{
        do i=1,iplas-1
           Bmax=0.d0
           Svol=0.d0         
          do j=2,nt-1
           B_ij=Btot(i,j)
           Bmax=dmax1(Bmax,B_ij)
           Svol=Svol+vol(i,j)
          enddo

           avr_v=0.d0
          do j=2,nt-1
           Bnor=Btot(i,j)/Bmax
           avr_v=avr_v+
     &          (b0ax/Btot(i,j))**2
     &         *( 1.d0-dsqrt(1.d0-Bnor)*(1.d0+.5d0*Bnor) )*vol(i,j)
          enddo
          
           trap(i)=avr_v/Svol
           traps(i+1)=trap(i)
           
        enddo
!!   }}}}}}}}}}}

C*NEW
 
c	goto 9991 
	ip0=1
	ip1=2
	ip2=3
	ip3=4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          x3=rhocs(ip3)
          x2=rhocs(ip2)
          x1=rhocs(ip1)
          x0=rhocs(ip0)
!!!!!!!!!!!!g11(1)
call	b_extrp(x0, X1,X2,X3, ip1,ip2,ip3,ip0,g11s,jerr)
!!!!!!!!!!!!g22(1)
call	b_extrp(x0,X1,X2,X3, ip1,ip2,ip3,ip0,g22s,jerr)
          
!!!!!!!!!!!!g33(1)
	call	b_extrp(x0,X1,X2,X3, ip1,ip2,ip3,ip0,g33s,jerr)
          
!!!!!!!!!!!!fnors(1)
	call	b_extrp(x0,X1,X2,X3, ip1,ip2,ip3,ip0,fnors,jerr)
 
!!!!!!!!!!!!gradrs(1)
	call	b_extrp(x0,X1,X2,X3, ip1,ip2,ip3,ip0,gradrs,jerr)

!!!!!!!!!!!!amus(1)
	call	b_extrp(x0,X1,X2,X3, ip1,ip2,ip3,ip0,amus,jerr)
!!!!!!!!!!!!tpaps(1)
	call	b_extrp(x0,X1,X2,X3, ip1,ip2,ip3,ip0,traps,jerr)


!          x0=rhocs(1)
!          x1=rhocs(2)
!          x2=rhocs(3)
!          x3=rhocs(4)

!
 9991	continue

	ip0=iplas+1
	ip1=iplas-2
	ip2=iplas-1
	ip3=iplas
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          x1=rhocs(ip3)
          x2=rhocs(ip2)
          x3=rhocs(ip1)
          x0=rhocs(ip0)


cw	write(*,*) 'passed 1'
!!!!!!!!!!!gradrs(iplas+1)
	call	b_extrp(x0, X1,X2,X3, ip3,ip2,ip1,ip0,gradrs,jerr)
 
!!!!!!!!!!!amus(iplas+1)
	call	b_extrp(x0, X1,X2,X3, ip3,ip2,ip1,ip0,amus,jerr)

        
!!!!!!!!!!!dvdro(iplas+1)
	call	b_extrp(x0, X1,X2,X3, ip3,ip2,ip1,ip0,dvdro,jerr)

 
!!!!!!!!!!!!g11(iplas+1)
	call	b_extrp(x0, X1,X2,X3, ip3,ip2,ip1,ip0,g11s,jerr)

!!!!!!!!!!!!g22(iplas+1)
	call	b_extrp(x0, X1,X2,X3, ip3,ip2,ip1,ip0,g22s,jerr)
          
!!!!!!!!!!!!g33(iplas+1)
	call	b_extrp(x0,X1,X2,X3, ip3,ip2,ip1,ip0,g33s,jerr)

!!!!!!!!!!!!b_maxt(iplas+1)
	call	b_extrp(x0, X1,X2,X3, ip3,ip2,ip1,ip0,b_maxt,jerr)

!!!!!!!!!!!!b_mint(iplas+1)
	call	b_extrp(x0, X1,X2,X3, ip3,ip2,ip1,ip0,b_mint,jerr)

!!!!!!!!!!!!b_db02(iplas+1)
	call	b_extrp(x0, X1,X2,X3, ip3,ip2,ip1,ip0,b_db02,jerr)


!!!!!!!!!!!!b_0db2(iplas+1)
	call	b_extrp(x0, X1,X2,X3, ip3,ip2,ip1,ip0,b_0db2,jerr)

 
!!!!!!!!!!!!b_db0(iplas+1)
	call	b_extrp(x0, X1,X2,X3, ip3,ip2,ip1,ip0,b_db0,jerr)

!!!!!!!!!!!!tpaps(iplas+1)
	call	b_extrp(x0, X1,X2,X3, ip3,ip2,ip1,ip0,traps,jerr)

!!!!!!!!!!!!!!!!!!!!!splining to astra grids

	rhocn(1)=0.d0
	rhocn(iplas+1)=1.d0
	dvdro(1)=0.d0
	sa(1)=0.d0
	g22s(1)=0.d0
	g11s(1)=0.d0
        fnors(iplas+1)=fvac/(rtor*btor)

cp        n3spl=iplas
        n3spl=iplas+1

cw	write(*,*) (b_mint(i), i=1,iplas)
!!!!!!!!!!!!!!!!!!!!! bmint !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CALL E01BAF(n3spl,rhocn,b_mint,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na
               !zspl=(rhowr(i)+rhowr(i+1))*0.5d0
               zspl=rhowrh(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               bmint(i)=CWk(1)
           enddo
               bmint(na1)=b_mint(iplas)

cw	write(*,*) 'bmint'

!!!!!!!!!!!!!!!!!!!!! bmaxt !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CALL E01BAF(n3spl,rhocn,b_maxt,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na
               !zspl=(rhowr(i)+rhowr(i+1))*0.5d0
               zspl=rhowrh(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               bmaxt(i)=CWk(1)
           enddo
               bmaxt(na1)=b_maxt(iplas)

cw	write(*,*) 'bmaxt'

!!!!!!!!!!!!!!!!!!!!! bdb02 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CALL E01BAF(n3spl,rhocn,b_db02,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na
               !zspl=(rhowr(i)+rhowr(i+1))*0.5d0
               zspl=rhowrh(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               bdb02(i)=CWk(1)
           enddo
               bdb02(na1)=b_db02(iplas)

cw	write(*,*) 'bdb02'

!!!!!!!!!!!!!!!!!!!!! b0db2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CALL E01BAF(n3spl,rhocn,b_0db2,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na
               !zspl=(rhowr(i)+rhowr(i+1))*0.5d0
               zspl=rhowrh(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               b0db2(i)=CWk(1)
           enddo
               b0db2(na1)=b_0db2(iplas)

cw	write(*,*) 'b0db2'

!!!!!!!!!!!!!!!!!!!!! bdb0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CALL E01BAF(n3spl,rhocn,b_db0,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na
               !zspl=(rhowr(i)+rhowr(i+1))*0.5d0
               zspl=rhowrh(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               bdb0(i)=CWk(1)
          enddo
              bdb0(na1)=b_db0(iplas)

cw	write(*,*) 'bdb0'

!!!!!!!!!!!!!!!!!!!!!! cu !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        n3spl=iplas

        CALL E01BAF(n3spl,rhosn,BJ_av,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na
               zspl=rhowr(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               cu(i)=CWk(1)/btor/amu0
           enddo
               cu(na1)=BJ_av(iplas)/btor/amu0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! pres !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        n3spl=iplas

        CALL E01BAF(n3spl,rhosn,prs_1,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na
               zspl=rhowr(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               pres(i)=CWk(1)
           enddo
               pres(na1)=prs_1(iplas)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! fp !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        n3spl=iplas

        CALL E01BAF(n3spl,rhosn,psi_1,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na
               zspl=rhowr(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               fp(i)=-CWk(1)
          enddo
               fp(na1)=-psi_1(iplas)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! W_Dj !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        n3spl=iplas

        CALL E01BAF(n3spl,rhosn,WDj,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na
               zspl=rhowr(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               W_Dj(i)=-CWk(1)
          enddo
               W_Dj(na1)=WDj(iplas)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          !do i=1,na1
          ! dfpdt(i)=(fp_0(i)-fp(i))/dtim
          !enddo

!!!!!!!!!!!!!!!!!!!!!! slat !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        n3spl=iplas

        CALL E01BAF(n3spl,rhosn,sa,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na
               !zspl=(rhowr(i)+rhowr(i+1))*0.5d0
               zspl=rhowrh(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               slat(i)=CWk(1)
           enddo
               slat(na1)=sa(iplas)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         n3spl=iplas+1

!!!!!!!!!!!!!!!!!!!!!!! ipol !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CALL E01BAF(n3spl,rhocn,fnors,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na1
               zspl=rhowr(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               ipol(i)=CWk(1)
           enddo
               ipol(na1)=fnors(iplas+1)
!!!!!!!!!!!!!!!!!!!!!!! mu !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CALL E01BAF(n3spl,rhocn,amus,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na
		       !zspl=(rhowr(i)+rhowr(i+1))*0.5d0
               zspl=rhowrh(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               mu(i)=CWk(1)
           enddo
               mu(na1)=amus(iplas+1)

!!!!!!!!!!!!!!!!!!!!!!! yFOFB !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CALL E01BAF(n3spl,rhocn,traps,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na
               !zspl=(rhowr(i)+rhowr(i+1))*0.5d0
               zspl=rhowrh(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               yFOFB(i)=CWk(1)
           enddo
               yFOFB(na1)=traps(iplas+1)


!!!!!!!!!!!!!!!! dvdro !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CALL E01BAF(n3spl,rhocn,dvdro,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na1
               zspl=rhowr(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               vr(i)=CWk(1)
           enddo

          do i=1,na
               !zspl=(rhowr(i)+rhowr(i+1))*0.5d0
               zspl=rhowrh(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               vrs(i)=CWk(1)
           enddo
               vrs(na1)=dvdro(iplas)

          
!!!!!!!!!!!!!!!!!! g11 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cw	write(*,*) 'g11s',(g11s(j),j=1,iplas)

        CALL E01BAF(n3spl,rhocn,g11s,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na
               !zspl=(rhowr(i)+rhowr(i+1))*0.5d0
               zspl=rhowrh(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               g11(i)=CWk(1)
           enddo
		g11(na1)=g11s(iplas+1)



!!!!!!!!!!!!!!!!!!!! g22 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CALL E01BAF(n3spl,rhocn,g22s,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na
               !zspl=(rhowr(i)+rhowr(i+1))*0.5d0
               zspl=rhowrh(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
              g22(i)=CWk(1)

          enddo

		g22(na1)=g22s(iplas+1)

!!!!!!!!!!!!!!!!!!!!!! g33 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CALL E01BAF(n3spl,rhocn,g33s,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na1
               zspl=rhowr(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               g33(i)=CWk(1)
           enddo


!!!!!!!!!!!!!!!!!!!!!!! gradro !!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        CALL E01BAF(n3spl,rhocn,gradrs,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na
               !zspl=(rhowr(i)+rhowr(i+1))*0.5d0
               zspl=rhowrh(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               gradro(i)=CWk(1)
           enddo
               gradro(na1)=gradrs(iplas+1)

CW		write(*,*) 'end'


           return
           end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           subroutine comet_(na1,platok,
     *                      rtor,btor,rho,roc,nstep,
     *                      g11,g22,g33,vr,vrs,slat,gradro,rocnew,
     *                      mu,ipol,bmaxt,bmint,bdb02,b0db2,bdb0,droda,
     *					  ybetpl,yli3)

         include 'double.inc'
         include 'dim.inc'
         parameter(nrpl=nrp+1) 
          parameter(nrpl4=nrpl+4,nrpl6=nrpl4*6)
         include 'compol.inc'
       common /combsh/ rm0,zm0,rc0,zc0,asp0,el_up,el_lw,tr_up,tr_lw,nbsh
          real*8 rhos(nrp),vols(nrpl),g11s(nrpl),g22s(nrpl),g33s(nrpl)
          real*8 gradrs(nrpl),amus(nrpl),fnors(nrpl)
          real*8 dvdro(nrpl),sa(nrp),rhocs(nrpl),rhocn(nrpl),rhosn(nrpl)
          real*8 b_maxt(nrpl),b_mint(nrpl),b_db02(nrpl),b_db0(nrpl),
     *		   b_0db2(nrpl),d_roda(nrpl)	
          real*8 RRK(nrpl4),CCK(nrpl4),WRK(nrpl6)
          real*8 CWK(4)
          real*8 rhowr(1000)         
        real*8 rtor,btor,rho(*),roc,ybetpl,yli3,
     *       g11(*),g22(*),g33(*),                          
     *       vr(*),vrs(*),slat(*),gradro(*),rocnew,droda(*), 
     *       mu(*),ipol(*),bmaxt(*),bmint(*),bdb02(*),b0db2(*),bdb0(*)
                          
           sqrt(xx)=dsqrt(xx)

         na=na1-1
       
         platok=tokp       
         do i=1,iplas
          !rhos(i)=sqrt(flucf(i)/(pi*btor))
          rhos(i)=sqrt(flx_fi(i)/(pi*btor))
         enddo
CP
cp
         rhocs(1)=0.d0
          rhos(1)=0.d0

         do i=2,iplas
          rhocs(i)=0.5d0*(rhos(i-1)+rhos(i))
         enddo

          rocnew=rhos(iplas)
          rhocs(iplas+1)=rhos(iplas)

c=============================andrey=02.03.2003vvvvvvvvvvvvv
          ARC0=rhocs(iplas+1)
          ARC1=rhocs(iplas)
          ARC2=rhocs(iplas-1)
          ARC3=rhocs(iplas-2)
          SPW1=(ARC0-ARC2)*(ARC0-ARC3)
     &        /(ARC1-ARC2)/(ARC1-ARC3)
          SPW2=(ARC0-ARC1)*(ARC0-ARC3)
     &        /(ARC2-ARC1)/(ARC2-ARC3)
          SPW3=(ARC0-ARC1)*(ARC0-ARC2)
     &        /(ARC3-ARC1)/(ARC3-ARC2)
          q(iplas)=SPW1*q(iplas-1)+SPW2*q(iplas-2)+SPW3*q(iplas-3)

c=============================andrey=02.03.2003^^^^^^^^^^^^^^
 


         do i=1,iplas+1
          rhocn(i)=rhocs(i)/rhos(iplas)
         enddo
          rhocn(1)=0.d0

         do i=1,iplas
          rhosn(i)=rhos(i)/rhos(iplas)
	    droda(i)=rhosn(i)
         enddo
          rhosn(1)=0.d0
cp

          rhocn(iplas+1)=1.d0
          rhosn(iplas)=1.d0

         do i=1,na
          rhowr(i)=rho(i)/roc
         enddo
          rhowr(na1)=1.d0-1.d-9

          sa(1)=0.d0
cw	write(*,*) (rhocn(j),j=1,iplas)
cw	write(*,*) (rhosn(j),j=1,iplas)

         do i=1,iplas-1

           dsqi=0.d0
           dvoli=0.d0
           av_r2=0.d0
           av_gr1=0.d0
           av_gr2=0.d0
           av_grr2=0.d0

          u1=rhos(i)
          u2=rhos(i+1)
          u3=rhos(i+1)
          u4=rhos(i)

          do j=2,nt1

           dvoli=dvoli+vol(i,j)
           dsqi=dsqi+sr(i+1,j)

          r1=r(i,j)
          r2=r(i+1,j)
          r3=r(i+1,j+1)
          r4=r(i,j+1)
          r0=(r1+r2+r3+r4)*0.25d0

          z1=z(i,j)
          z2=z(i+1,j)
          z3=z(i+1,j+1)
          z4=z(i,j+1)

           sss=s(i,j)
           
           dadr= 0.5d0*((u1+u2)*(z2-z1)+(u2+u3)*(z3-z2)+
     *                   (u3+u4)*(z4-z3)+(u4+u1)*(z1-z4))/sss

           dadz=-0.5d0*((u1+u2)*(r2-r1)+(u2+u3)*(r3-r2)+
     *                   (u3+u4)*(r4-r3)+(u4+u1)*(r1-r4))/sss

           gr2=dadr**2+dadz**2

           av_gr1=av_gr1+vol(i,j)*sqrt(gr2)
           av_gr2=av_gr2+vol(i,j)*gr2
           av_grr2=av_grr2+vol(i,j)*gr2/r0**2
           av_r2=av_r2+vol(i,j)/r0**2

          enddo
          sa(i+1)=dsqi*2.d0*pi
c          vols(i+1)=dvoli
          gradrs(i+1)=av_gr1/dvoli
          amus(i+1)=(2.d0*pi)/q(i)
          fnors(i+1)=f(i)/(rtor*btor)

          dvdro(i+1)=2.d0*pi*dvoli/(rhos(i+1)-rhos(i))
cc          g11s(i+1)=av_gr2/dvoli
cc          g22s(i+1)=av_grr2/(dvoli*4.d0*pi*pi)*rtor
          g11s(i+1)=av_gr2*dvdro(i+1)/dvoli
        g22s(i+1)=av_grr2*dvdro(i+1)/(dvoli*4.d0*pi*pi)*rtor/fnors(i+1)
          g33s(i+1)=av_r2*rtor**2/dvoli
         enddo

cccc      magnetic field

CW	write(*,*) 'before b'

           bp2av=0.d0
		bp2=0.d0
cc          do i=2,iplas-1 
          do i=1,iplas-1 

           dvoli=0.d0
           bmx=0.d0
           bmn=1.d9*btor
           bavr=0.d0
           bavr2=0.d0
           bavrd2=0.d0

           do j=2,nt1

           dvoli=dvoli+vol(i,j)

          r1=r(i,j)
          r2=r(i+1,j)
          r3=r(i+1,j+1)
          r4=r(i,j+1)
          r0=(r1+r2+r3+r4)*0.25d0

            bm_pol=psim*(psia(i+1)-psia(i))/st(i,j)
            bp_pol=psim*(psia(i+1)-psia(i))/st(i,j+1)

        bpol2v=bm_pol**2*( vol1(i,j)/sin1(i,j)+vol2(i,j)/sin2(i,j))+ 
     +         bp_pol**2*( vol3(i,j)/sin3(i,j)+vol4(i,j)/sin4(i,j)) 

            bpol2=bpol2v /vol(i,j)
            b_fi2=(f(i)/r0)**2

            b_tot2=bpol2+b_fi2
            b_tot=dsqrt(b_tot2)

            bmx=dmax1(b_tot,bmx)
            bmn=dmin1(b_tot,bmn)

           bavr2=bavr2+b_tot2*vol(i,j)
           bavr =bavr +b_tot *vol(i,j)
	   bavrd2=bavrd2+vol(i,j)/(b_tot2+1.d-8)
           bp2av=bp2av+bpol2v
		bp2=bp2+bpol2v
          enddo            
           b_maxt(i+1)=bmx
           b_mint(i+1)=bmn
           b_db02(i+1)=bavr2/(dvoli*btor**2)
           b_0db2(i+1)=(bavrd2*btor**2)/(dvoli)
           b_db0(i+1)=bavr/(dvoli*btor)


        enddo 
           b_maxt(1)=.5d0*(b_maxt(2)+b_mint(2))
           b_mint(1)=b_maxt(1)
           b_db02(1)=b_db02(2)
           b_0db2(1)=b_0db2(2)
           b_db0(1)=b_db0(2)


c            yli3  =4.d0*pi*bp2av/(rtor*tokp*tokp)
            yli3  =  4.d0*pi*bp2/(rtor*(.4d0*pi*tokp)**2)     
C*NEW
 
c	goto 9991 
	ip0=1
	ip1=2
	ip2=3
	ip3=4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          x3=rhocs(ip3)
          x2=rhocs(ip2)
          x1=rhocs(ip1)
          x0=rhocs(ip0)
!!!!!!!!!!!!g11(1)
call	b_extrp(x0, X1,X2,X3, ip1,ip2,ip3,ip0,g11s,jerr)
!!!!!!!!!!!!g22(1)
call	b_extrp(x0,X1,X2,X3, ip1,ip2,ip3,ip0,g22s,jerr)
          
!!!!!!!!!!!!g33(1)
	call	b_extrp(x0,X1,X2,X3, ip1,ip2,ip3,ip0,g33s,jerr)
          
!!!!!!!!!!!!fnors(1)
	call	b_extrp(x0,X1,X2,X3, ip1,ip2,ip3,ip0,fnors,jerr)
 
!!!!!!!!!!!!gradrs(1)
	call	b_extrp(x0,X1,X2,X3, ip1,ip2,ip3,ip0,gradrs,jerr)

!!!!!!!!!!!!amus(1)
	call	b_extrp(x0,X1,X2,X3, ip1,ip2,ip3,ip0,amus,jerr)


!          x0=rhocs(1)
!          x1=rhocs(2)
!          x2=rhocs(3)
!          x3=rhocs(4)

!
 9991	continue

	ip0=iplas+1
	ip1=iplas-2
	ip2=iplas-1
	ip3=iplas
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          x1=rhocs(ip3)
          x2=rhocs(ip2)
          x3=rhocs(ip1)
          x0=rhocs(ip0)


cw	write(*,*) 'passed 1'
!!!!!!!!!!!gradrs(iplas+1)
	call	b_extrp(x0, X1,X2,X3, ip3,ip2,ip1,ip0,gradrs,jerr)
 
!!!!!!!!!!!amus(iplas+1)
	call	b_extrp(x0, X1,X2,X3, ip3,ip2,ip1,ip0,amus,jerr)

        
!!!!!!!!!!!dvdro(iplas+1)
	call	b_extrp(x0, X1,X2,X3, ip3,ip2,ip1,ip0,dvdro,jerr)

 
!!!!!!!!!!!!g11(iplas+1)
	call	b_extrp(x0, X1,X2,X3, ip3,ip2,ip1,ip0,g11s,jerr)

!!!!!!!!!!!!g22(iplas+1)
	call	b_extrp(x0, X1,X2,X3, ip3,ip2,ip1,ip0,g22s,jerr)
          
!!!!!!!!!!!!g33(iplas+1)
	call	b_extrp(x0,X1,X2,X3, ip3,ip2,ip1,ip0,g33s,jerr)

!!!!!!!!!!!!b_maxt(iplas+1)
	call	b_extrp(x0, X1,X2,X3, ip3,ip2,ip1,ip0,b_maxt,jerr)

!!!!!!!!!!!!b_mint(iplas+1)
	call	b_extrp(x0, X1,X2,X3, ip3,ip2,ip1,ip0,b_mint,jerr)

!!!!!!!!!!!!b_db02(iplas+1)
	call	b_extrp(x0, X1,X2,X3, ip3,ip2,ip1,ip0,b_db02,jerr)


!!!!!!!!!!!!b_0db2(iplas+1)
	call	b_extrp(x0, X1,X2,X3, ip3,ip2,ip1,ip0,b_0db2,jerr)

 
!!!!!!!!!!!!b_db0(iplas+1)
	call	b_extrp(x0, X1,X2,X3, ip3,ip2,ip1,ip0,b_db0,jerr)


!!!!!!!!!!!!!!!!!!!!!splining to astra grids

	rhocn(1)=0.d0
	rhocn(iplas+1)=1.d0
	dvdro(1)=0.d0
	sa(1)=0.d0
	g22s(1)=0.d0
	g11s(1)=0.d0
        fnors(iplas+1)=fvac/(rtor*btor)

cp        n3spl=iplas
        n3spl=iplas+1

cw	write(*,*) (b_mint(i), i=1,iplas)
!!!!!!!!!!!!!!!!!!!!! bmint !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CALL E01BAF(n3spl,rhocn,b_mint,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na
               zspl=(rhowr(i)+rhowr(i+1))*0.5d0
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               bmint(i)=CWk(1)
           enddo
               bmint(na1)=b_mint(iplas)

cw	write(*,*) 'bmint'

!!!!!!!!!!!!!!!!!!!!! bmaxt !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CALL E01BAF(n3spl,rhocn,b_maxt,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na
               zspl=(rhowr(i)+rhowr(i+1))*0.5d0
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               bmaxt(i)=CWk(1)
           enddo
               bmaxt(na1)=b_maxt(iplas)

cw	write(*,*) 'bmaxt'

!!!!!!!!!!!!!!!!!!!!! bdb02 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CALL E01BAF(n3spl,rhocn,b_db02,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na
               zspl=(rhowr(i)+rhowr(i+1))*0.5d0
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               bdb02(i)=CWk(1)
           enddo
               bdb02(na1)=b_db02(iplas)

cw	write(*,*) 'bdb02'

!!!!!!!!!!!!!!!!!!!!! b0db2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CALL E01BAF(n3spl,rhocn,b_0db2,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na
               zspl=(rhowr(i)+rhowr(i+1))*0.5d0
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               b0db2(i)=CWk(1)
           enddo
               b0db2(na1)=b_0db2(iplas)

cw	write(*,*) 'b0db2'

!!!!!!!!!!!!!!!!!!!!! bdb0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CALL E01BAF(n3spl,rhocn,b_db0,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na
cp               zspl=(rhowr(i)+rhowr(i+1))*0.5d0
               zspl=(rhowr(i)+rhowr(i+1))*0.5d0
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               bdb0(i)=CWk(1)
          enddo
              bdb0(na1)=b_db0(iplas)

cw	write(*,*) 'bdb0'


!!!!!!!!!!!!!!!!!!!!!! slat !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        n3spl=iplas

        CALL E01BAF(n3spl,rhosn,sa,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na
               zspl=(rhowr(i)+rhowr(i+1))*0.5d0
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               slat(i)=CWk(1)
           enddo
               slat(na1)=sa(iplas)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         n3spl=iplas+1

!!!!!!!!!!!!!!!!!!!!!!! ipol !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CALL E01BAF(n3spl,rhocn,fnors,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na1
               zspl=rhowr(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               ipol(i)=CWk(1)
           enddo
               ipol(na1)=fnors(iplas+1)
!!!!!!!!!!!!!!!!!!!!!!! mu !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CALL E01BAF(n3spl,rhocn,amus,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na
		zspl=(rhowr(i)+rhowr(i+1))*0.5d0
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               mu(i)=CWk(1)
           enddo
               mu(na1)=amus(iplas+1)



!!!!!!!!!!!!!!!! dvdro !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CALL E01BAF(n3spl,rhocn,dvdro,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na1
               zspl=rhowr(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               vr(i)=CWk(1)
           enddo

          do i=1,na
               zspl=(rhowr(i)+rhowr(i+1))*0.5d0
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               vrs(i)=CWk(1)
           enddo
               vrs(na1)=dvdro(iplas)

          
!!!!!!!!!!!!!!!!!! g11 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cw	write(*,*) 'g11s',(g11s(j),j=1,iplas)

        CALL E01BAF(n3spl,rhocn,g11s,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na
		zspl=(rhowr(i)+rhowr(i+1))*0.5d0
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               g11(i)=CWk(1)
           enddo
		g11(na1)=g11s(iplas+1)



!!!!!!!!!!!!!!!!!!!! g22 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CALL E01BAF(n3spl,rhocn,g22s,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na
		zspl=(rhowr(i)+rhowr(i+1))*0.5d0
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
              g22(i)=CWk(1)

          enddo

		g22(na1)=g22s(iplas+1)

!!!!!!!!!!!!!!!!!!!!!! g33 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CALL E01BAF(n3spl,rhocn,g33s,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na1
               zspl=rhowr(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               g33(i)=CWk(1)
           enddo


!!!!!!!!!!!!!!!!!!!!!!! gradro !!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        CALL E01BAF(n3spl,rhocn,gradrs,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na
               zspl=(rhowr(i)+rhowr(i+1))*0.5d0
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               gradro(i)=CWk(1)
           enddo
               gradro(na1)=gradrs(iplas+1)

CW		write(*,*) 'end'


           return
           end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       subroutine profro(na1,neql,rtor,eqpf,eqff,rho,roc)
c	29 JAN 2003
           include 'double.inc'
          parameter(nursp=1000,nursp4=nursp+4,nursp6=nursp4*6)
           include 'dim.inc'
           include 'compol.inc'

       real*8 rosn(nrp),ppr(nursp),ffpr(nursp),rhow(nursp)
          real*8 rrk(nursp4),cck(nursp4),wrk(nursp6)
          real*8 cwk(4)

        real*8 eqpf(*),eqff(*),rho(*),rtor,roc

          !return

        iplas=neql
        nspl=na1+1
         do i=1,iplas
          rosn(i)=(i-1.d0)/(iplas-1.d0)
         enddo

!! exstrapolation on magn. axis
       
             rh0=0.d0
		              
             rh1=rho(1)           
             rh2=rho(2)           
             rh3=rho(3)           

             ppx1=eqpf(1)
             ppx2=eqpf(2)
             ppx3=eqpf(3)

             ffx1=eqff(1)
             ffx2=eqff(2)
             ffx3=eqff(3)

       call EXTRP2(rh0,ppx0, rh1,rh2,rh3, ppx1,ppx2,ppx3)
       call EXTRP2(rh0,ffx0, rh1,rh2,rh3, ffx1,ffx2,ffx3)

           ppr(1)=ppx0
           ffpr(1)=ffx0
           rhow(1)=0.d0

          do i=2,na1+1
           ppr(i)=eqpf(i-1)
           ffpr(i)=eqff(i-1)
           rhow(i)=rho(i-1)/roc
          enddo		

           dpdpsi(1)=ppr(1)
           dfdpsi(1)=ffpr(1)

           dpdpsi(iplas)=ppr(na1+1)
           dfdpsi(iplas)=ffpr(na1+1)


           CALL E01BAF(nspl,rhow,ppr,RRK,CCK,
     *                        nspl+4,WRK,6*nspl+16,IFAIL)

           do i=2,iplas-1

              zrho=rosn(i)
              CALL E02BCF(nspl+4,RRK,CCK,zrho,0,CWk,IFAIL)
              dpdpsi(i)=cwk(1)

           enddo


           CALL E01BAF(nspl,rhow,ffpr,RRK,CCK,
     *                        nspl+4,WRK,6*nspl+16,IFAIL)

           do i=2,iplas-1

              zrho=rosn(i)
              CALL E02BCF(nspl+4,RRK,CCK,zrho,0,CWk,IFAIL)
              dfdpsi(i)=cwk(1)

           enddo

              do i=1,iplas
                dpdpsi(i)=dpdpsi(i)/rtor
                dfdpsi(i)=dfdpsi(i)*rtor
              enddo

            return
            end
		 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       subroutine profro_p(na1,neql,rtor,eqpf,rho,roc)
c	march 2007
           include 'double.inc'
          parameter(nursp=1000,nursp4=nursp+4,nursp6=nursp4*6)
           include 'dim.inc'
           include 'compol.inc'

       real*8 rosn(nrp),rhot(nursp),ppr(nursp),ffpr(nursp),rhow(nursp)
          real*8 rrk(nursp4),cck(nursp4),wrk(nursp6)
          real*8 cwk(4)

        real*8 eqpf(*),rho(*),rtor,roc

        iplas=neql
        nspl=na1+1
         do i=1,iplas
          rosn(i)=(i-1.d0)/(iplas-1.d0)
         enddo

!! exstrapolation on magn. axis
       
             rh0=0.d0
		              
             rh1=rho(1)           
             rh2=rho(2)           
             rh3=rho(3)           

             ppx1=eqpf(1)
             ppx2=eqpf(2)
             ppx3=eqpf(3)

             !ffx1=eqff(1)
             !ffx2=eqff(2)
             !ffx3=eqff(3)

       call EXTRP2(rh0,ppx0, rh1,rh2,rh3, ppx1,ppx2,ppx3)
       !call EXTRP2(rh0,ffx0, rh1,rh2,rh3, ffx1,ffx2,ffx3)

           ppr(1)=ppx0
           !ffpr(1)=ffx0
           rhow(1)=0.d0

          do i=2,na1+1
           ppr(i)=eqpf(i-1)
           !ffpr(i)=eqff(i-1)
           rhow(i)=rho(i-1)/roc
          enddo		

           dpdpsi(1)=ppr(1)
           !dfdpsi(1)=ffpr(1)

           dpdpsi(iplas)=ppr(na1+1)
           !dfdpsi(iplas)=ffpr(na1+1)


           CALL E01BAF(nspl,rhow,ppr,RRK,CCK,
     *                        nspl+4,WRK,6*nspl+16,IFAIL)

           do i=2,iplas-1

              zrho=rosn(i)
              CALL E02BCF(nspl+4,RRK,CCK,zrho,0,CWk,IFAIL)
              dpdpsi(i)=cwk(1)

           enddo


!           CALL E01BAF(nspl,rhow,ffpr,RRK,CCK,
!     *                        nspl+4,WRK,6*nspl+16,IFAIL)
!
!           do i=2,iplas-1
!
!              zrho=rosn(i)
!              CALL E02BCF(nspl+4,RRK,CCK,zrho,0,CWk,IFAIL)
!              dfdpsi(i)=cwk(1)
!
!           enddo

              do i=1,iplas
                dpdpsi(i)=dpdpsi(i)/rtor
                !dfdpsi(i)=dfdpsi(i)*rtor
              enddo

            return
            end
		 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       subroutine S_profro_p(equil_in,parameters_spider)
c	apr 2013
         use imas_ids       
         use parameters_a2spider, only: type_parameters, TWOPI 

           include 'double.inc'
          parameter(nursp=1000,nursp4=nursp+4,nursp6=nursp4*6)
           include 'dim.inc'
           include 'compol.inc'

       real*8 rosn(nrp),rhot(nursp),ppr(nursp),ffpr(nursp),rhow(nursp)
          real*8 rrk(nursp4),cck(nursp4),wrk(nursp6)
          real*8 cwk(4)

        !real*8 eqpf(*),rho(*),rtor,roc
       type(type_equilibrium) equil_in
       type(type_parameters) parameters_spider

        iplas=parameters_spider%neql
        nspl=size(equil_in%profiles_1d%psi)

         do i=1,iplas
          rosn(i)=(i-1.d0)/(iplas-1.d0)
         enddo


            roc = equil_in%profiles_1d%rho_tor(nspl)

          do i=1,nspl
           rhow(i)=equil_in%profiles_1d%rho_tor(i)/roc
           ppr(i) = equil_in%profiles_1d%pprime(i)
          enddo	

           dpdpsi(1)=ppr(1)

           dpdpsi(iplas)=ppr(nspl)

           CALL E01BAF(nspl,rhow,ppr,RRK,CCK,
     *                        nspl+4,WRK,6*nspl+16,IFAIL)

           do i=2,iplas-1

              zrho=rosn(i)
              CALL E02BCF(nspl+4,RRK,CCK,zrho,0,CWk,IFAIL)
              dpdpsi(i)=cwk(1)

           enddo


!           CALL E01BAF(nspl,rhow,ffpr,RRK,CCK,
!     *                        nspl+4,WRK,6*nspl+16,IFAIL)
!
!           do i=2,iplas-1
!
!              zrho=rosn(i)
!              CALL E02BCF(nspl+4,RRK,CCK,zrho,0,CWk,IFAIL)
!              dfdpsi(i)=cwk(1)
!
!           enddo
          pscale = 2 * TWOPI * 1.d-7 * TWOPI

              do i=1,iplas
                dpdpsi(i)=-dpdpsi(i)*pscale
                !dfdpsi(i)=dfdpsi(i)*rtor
              enddo

            return
            end
		 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       subroutine profro_p_L(na1,neql,rtor,eqpf,rho,roc)
c	march 2007
           include 'double.inc'
          parameter(nursp=1000,nursp4=nursp+4,nursp6=nursp4*6)
           include 'dim.inc'
           include 'compol.inc'

       real*8 rosn(nrp),rhot(nursp),ppr(nursp),ffpr(nursp),rhow(nursp)
          real*8 rrk(nursp4),cck(nursp4),wrk(nursp6)
          real*8 cwk(4)

        real*8 eqpf(*),rho(*),rtor,roc

        iplas=neql
        nspl=na1+1
         do i=1,iplas
          rosn(i)=(i-1.d0)/(iplas-1.d0)
         enddo

!! exstrapolation on magn. axis
       
             rh0=0.d0
		              
             rh1=rho(1)           
             rh2=rho(2)           
             rh3=rho(3)           

             ppx1=eqpf(1)
             ppx2=eqpf(2)
             ppx3=eqpf(3)

             !ffx1=eqff(1)
             !ffx2=eqff(2)
             !ffx3=eqff(3)

       call EXTRP2(rh0,ppx0, rh1,rh2,rh3, ppx1,ppx2,ppx3)
       !call EXTRP2(rh0,ffx0, rh1,rh2,rh3, ffx1,ffx2,ffx3)

           ppr(1)=ppr(2)
           !ppr(1)=ppx0
           !ffpr(1)=ffx0
           rhow(1)=0.d0

          do i=2,na1+1
           ppr(i)=eqpf(i-1)
           !ffpr(i)=eqff(i-1)
           rhow(i)=rho(i-1)/roc
          enddo		

           dpdpsi(1)=ppr(1)
           !dfdpsi(1)=ffpr(1)

           dpdpsi(iplas)=ppr(na1+1)
           !dfdpsi(iplas)=ffpr(na1+1)


!linear interpolation

           do i=2,iplas-1
              zrho=rosn(i)

            do ic=1,nspl-1
             if(zrho.le.rhow(ic+1) .AnD. zrho.gt.rhow(ic)) then
              dpdpsi(i)=( ppr(ic)*(rhow(ic+1)-zrho)
     *                   +ppr(ic+1)*(zrho-rhow(ic)) )
     *                      /(rhow(ic+1)-rhow(ic))
              exit             
             endif             
            enddo


           enddo






!           CALL E01BAF(nspl,rhow,ppr,RRK,CCK,
!     *                        nspl+4,WRK,6*nspl+16,IFAIL)

!           do i=2,iplas-1
!
!              zrho=rosn(i)
!              CALL E02BCF(nspl+4,RRK,CCK,zrho,0,CWk,IFAIL)
!              dpdpsi(i)=cwk(1)
!
!           enddo


!           CALL E01BAF(nspl,rhow,ffpr,RRK,CCK,
!     *                        nspl+4,WRK,6*nspl+16,IFAIL)
!
!           do i=2,iplas-1
!
!              zrho=rosn(i)
!              CALL E02BCF(nspl+4,RRK,CCK,zrho,0,CWk,IFAIL)
!              dfdpsi(i)=cwk(1)
!
!           enddo

              do i=1,iplas
                dpdpsi(i)=dpdpsi(i)/rtor
                !dfdpsi(i)=dfdpsi(i)*rtor
              enddo

            return
            end
		 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       subroutine profro_sbc(na1,neql,rtor,btor,cc,Te,cubs,cd,rho,roc)
c	03 JAN 2008
           include 'double.inc'
          parameter(nursp=1000,nursp4=nursp+4,nursp6=nursp4*6)
           include 'dim.inc'
           include 'compol.inc'
	    common /com_sigcd/ C_sig(nrp),T_el(nrp),C_bts(nrp),C_driv(nrp)

          real*8 rosn(nrp),funw(nursp),rhow(nursp)
          real*8 rrk(nursp4),cck(nursp4),wrk(nursp6)
          real*8 cwk(4)

        real*8 rho(*),rtor,btor,roc,cc(*),Te(*),cubs(*),cd(*)

        iplas=neql
        nspl=na1+1
         do i=1,iplas
          rosn(i)=(i-1.d0)/(iplas-1.d0)
         enddo

!! exstrapolation on magn. axis
       
             rh0=0.d0
		              
             rh1=rho(1)           
             rh2=rho(2)           
             rh3=rho(3)           

             ccx1=cc(1)
             ccx2=cc(2)
             ccx3=cc(3)

             cdx1=cd(1)
             cdx2=cd(2)
             cdx3=cd(3)

             cbx1=cubs(1)
             cbx2=cubs(2)
             cbx3=cubs(3)

             cex1=te(1)
             cex2=te(2)
             cex3=te(3)

       call EXTRP2(rh0,ccx0, rh1,rh2,rh3, ccx1,ccx2,ccx3)
       call EXTRP2(rh0,cdx0, rh1,rh2,rh3, cdx1,cdx2,cdx3)
       call EXTRP2(rh0,cbx0, rh1,rh2,rh3, cbx1,cbx2,cbx3)
       call EXTRP2(rh0,cex0, rh1,rh2,rh3, cex1,cex2,cex3)

!!!!!!!!!!!!!!!!!!!cc!!!!!!!!!!!!!!!!!!!!!!!

           funw(1)=ccx0
           rhow(1)=0.d0

          do i=2,na1+1
           funw(i)=cc(i-1)
           rhow(i)=rho(i-1)/roc
          enddo		

           C_sig(1)=funw(1)
           C_sig(iplas)=funw(na1+1)

           CALL E01BAF(nspl,rhow,funw,RRK,CCK,
     *                        nspl+4,WRK,6*nspl+16,IFAIL)

           do i=2,iplas-1

              zrho=rosn(i)
              CALL E02BCF(nspl+4,RRK,CCK,zrho,0,CWk,IFAIL)
              C_sig(i)=cwk(1)

           enddo

 !!!!!!!!!!!!!!!!!!!cc!!!!!!!!!!!!!!!!!!!!!!!
  
!!!!!!!!!!!!!!!!!!!cd!!!!!!!!!!!!!!!!!!!!!!!

           funw(1)=cdx0
           rhow(1)=0.d0

          do i=2,na1+1
           funw(i)=cd(i-1)
           rhow(i)=rho(i-1)/roc
          enddo		

           C_driv(1)=funw(1)
           C_driv(iplas)=funw(na1+1)

           CALL E01BAF(nspl,rhow,funw,RRK,CCK,
     *                        nspl+4,WRK,6*nspl+16,IFAIL)

           do i=2,iplas-1

              zrho=rosn(i)
              CALL E02BCF(nspl+4,RRK,CCK,zrho,0,CWk,IFAIL)
              C_driv(i)=cwk(1)

           enddo

 !!!!!!!!!!!!!!!!!!!cd!!!!!!!!!!!!!!!!!!!!!!!
  
 !!!!!!!!!!!!!!!!!!!cubs!!!!!!!!!!!!!!!!!!!!!!!

           funw(1)=cbx0
           rhow(1)=0.d0

          do i=2,na1+1
           funw(i)=cubs(i-1)
           rhow(i)=rho(i-1)/roc
          enddo		

           C_bts(1)=funw(1)
           C_bts(iplas)=funw(na1+1)

           CALL E01BAF(nspl,rhow,funw,RRK,CCK,
     *                        nspl+4,WRK,6*nspl+16,IFAIL)

           do i=2,iplas-1

              zrho=rosn(i)
              CALL E02BCF(nspl+4,RRK,CCK,zrho,0,CWk,IFAIL)
              C_bts(i)=cwk(1)

           enddo

 !!!!!!!!!!!!!!!!!!!cubs!!!!!!!!!!!!!!!!!!!!!!!
 
 !!!!!!!!!!!!!!!!!!!te!!!!!!!!!!!!!!!!!!!!!!!

           funw(1)=cex0
           rhow(1)=0.d0

          do i=2,na1+1
           funw(i)=te(i-1)
           rhow(i)=rho(i-1)/roc
          enddo		

           T_el(1)=funw(1)
           T_el(iplas)=funw(na1+1)

           CALL E01BAF(nspl,rhow,funw,RRK,CCK,
     *                        nspl+4,WRK,6*nspl+16,IFAIL)

           do i=2,iplas-1

              zrho=rosn(i)
              CALL E02BCF(nspl+4,RRK,CCK,zrho,0,CWk,IFAIL)
              T_el(i)=cwk(1)

           enddo

 !!!!!!!!!!!!!!!!!!!te!!!!!!!!!!!!!!!!!!!!!!!
  
             do i=1,iplas
                C_bts(i)=C_bts(i)*btor
                C_driv(i)=C_driv(i)*btor
                T_el(i)=T_el(i)*1.d3
              enddo

            return
            end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       subroutine S_profro_sbc(equil_in,parameters_spider)
c	 APR 2013

         use imas_ids       
         use parameters_a2spider , only: type_parameters, TWOPI

           include 'double.inc'
          parameter(nursp=1000,nursp4=nursp+4,nursp6=nursp4*6)
           include 'dim.inc'
           include 'compol.inc'
	    common /com_sigcd/ C_sig(nrp),T_el(nrp),C_bts(nrp),C_driv(nrp)

          real*8 rosn(nrp),funw(nursp),rhow(nursp)
          real*8 rrk(nursp4),cck(nursp4),wrk(nursp6)
          real*8 cwk(4)

        !real*8 rho(*),rtor,btor,roc,cc(*),Te(*),cubs(*),cd(*)
       type(type_equilibrium) equil_in
       type(type_parameters) parameters_spider

        iplas=parameters_spider%neql
        nspl=size(equil_in%profiles_1d%psi)
        
         do i=1,iplas
          rosn(i)=(i-1.d0)/(iplas-1.d0)
         enddo

!!!!!!!!!!!!!!!!!!!cc!!!!!!!!!!!!!!!!!!!!!!!
            roc = equil_in%profiles_1d%rho_tor(nspl)

          do i=1,nspl
           rhow(i)=equil_in%profiles_1d%rho_tor(i)/roc
           funw(i)=equil_in%profiles_1d%sigmapar%value(i)
          enddo	
          	

           C_sig(1)=funw(1)
           C_sig(iplas)=funw(nspl)

           CALL E01BAF(nspl,rhow,funw,RRK,CCK,
     *                        nspl+4,WRK,6*nspl+16,IFAIL)

           do i=2,iplas-1

              zrho=rosn(i)
              CALL E02BCF(nspl+4,RRK,CCK,zrho,0,CWk,IFAIL)
              C_sig(i)=cwk(1)

           enddo

 !!!!!!!!!!!!!!!!!!!cc!!!!!!!!!!!!!!!!!!!!!!!!!
  
 !!!!!!!!!!!!!!!!!!!cubs!!!!!!!!!!!!!!!!!!!!!!!

          do i=1,nspl
           funw(i) = equil_in%profiles_1d%jni%value(i)
          enddo	

           C_bts(1)=funw(1)
           C_bts(iplas)=funw(nspl)

           CALL E01BAF(nspl,rhow,funw,RRK,CCK,
     *                        nspl+4,WRK,6*nspl+16,IFAIL)

           do i=2,iplas-1

              zrho=rosn(i)
              CALL E02BCF(nspl+4,RRK,CCK,zrho,0,CWk,IFAIL)
              C_bts(i)=cwk(1)

           enddo

 !!!!!!!!!!!!!!!!!!!cubs!!!!!!!!!!!!!!!!!!!!!!!
 
 !!!!!!!!!!!!!!!!!!!te!!!!!!!!!!!!!!!!!!!!!!!

          do i=1,nspl
           funw(i) = equil_in%profiles_1d%te%value(i)
          enddo	

           T_el(1)=(1)
           T_el(iplas)=funw(nspl)

           CALL E01BAF(nspl,rhow,funw,RRK,CCK,
     *                        nspl+4,WRK,6*nspl+16,IFAIL)

           do i=2,iplas-1

              zrho=rosn(i)
              CALL E02BCF(nspl+4,RRK,CCK,zrho,0,CWk,IFAIL)
              T_el(i)=cwk(1)

           enddo

 !!!!!!!!!!!!!!!!!!!te!!!!!!!!!!!!!!!!!!!!!!!
  
             do i=1,iplas
                C_bts(i)=C_bts(i)*btor
                C_driv(i)=0.d0
                T_el(i)=T_el(i)*1.d3
              enddo

            return
            end
		 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine profro_sbc_L(na1,neql,rtor,btor,cc,Te,cubs,cd,rho,roc)
 
            include 'double.inc'
          parameter(nursp=1000,nursp4=nursp+4,nursp6=nursp4*6)
           include 'dim.inc'
           include 'compol.inc'
	    common /com_sigcd/ C_sig(nrp),T_el(nrp),C_bts(nrp),C_driv(nrp)

          real*8 rosn(nrp),funw(nursp),rhow(nursp)
          real*8 rrk(nursp4),cck(nursp4),wrk(nursp6)
          real*8 cwk(4)

        real*8 rho(*),rtor,btor,roc,cc(*),Te(*),cubs(*),cd(*)

        iplas=neql
        nspl=na1+1
         do i=1,iplas
          rosn(i)=(i-1.d0)/(iplas-1.d0)
         enddo

!! exstrapolation on magn. axis
       
             rh0=0.d0
		              
             rh1=rho(1)           
             rh2=rho(2)           
             rh3=rho(3)           

             ccx1=cc(1)
             ccx2=cc(2)
             ccx3=cc(3)

             cdx1=cd(1)
             cdx2=cd(2)
             cdx3=cd(3)

             cbx1=cubs(1)
             cbx2=cubs(2)
             cbx3=cubs(3)

             cex1=te(1)
             cex2=te(2)
             cex3=te(3)

       call EXTRP2(rh0,ccx0, rh1,rh2,rh3, ccx1,ccx2,ccx3)
       call EXTRP2(rh0,cdx0, rh1,rh2,rh3, cdx1,cdx2,cdx3)
       call EXTRP2(rh0,cbx0, rh1,rh2,rh3, cbx1,cbx2,cbx3)
       call EXTRP2(rh0,cex0, rh1,rh2,rh3, cex1,cex2,cex3)

!!!!!!!!!!!!!!!!!!!cc!!!!!!!!!!!!!!!!!!!!!!!

           !funw(1)=ccx0
           funw(1)=ccx1
           rhow(1)=0.d0

          do i=2,na1+1
           funw(i)=cc(i-1)
           rhow(i)=rho(i-1)/roc
          enddo		

           C_sig(1)=funw(1)
           C_sig(iplas)=funw(na1+1)

!           CALL E01BAF(nspl,rhow,funw,RRK,CCK,
!     *                        nspl+4,WRK,6*nspl+16,IFAIL)
!
!           do i=2,iplas-1
!
!              zrho=rosn(i)
!              CALL E02BCF(nspl+4,RRK,CCK,zrho,0,CWk,IFAIL)
!              C_sig(i)=cwk(1)
!
!           enddo

            nspl=na1+1
           call Linsplin(nspl,funw,rhow, iplas,C_sig,rosn)

 !!!!!!!!!!!!!!!!!!!cc!!!!!!!!!!!!!!!!!!!!!!!
  
!!!!!!!!!!!!!!!!!!!cd!!!!!!!!!!!!!!!!!!!!!!!

           !funw(1)=cdx0
           funw(1)=cdx1
           rhow(1)=0.d0

          do i=2,na1+1
           funw(i)=cd(i-1)
           rhow(i)=rho(i-1)/roc
          enddo		

           C_driv(1)=funw(1)
           C_driv(iplas)=funw(na1+1)

!           CALL E01BAF(nspl,rhow,funw,RRK,CCK,
!     *                        nspl+4,WRK,6*nspl+16,IFAIL)
!
!           do i=2,iplas-1
!
!              zrho=rosn(i)
!              CALL E02BCF(nspl+4,RRK,CCK,zrho,0,CWk,IFAIL)
!              C_driv(i)=cwk(1)
!
!           enddo

            nspl=na1+1
           call Linsplin(nspl,funw,rhow, iplas,C_driv,rosn)

 !!!!!!!!!!!!!!!!!!!cd!!!!!!!!!!!!!!!!!!!!!!!
  
 !!!!!!!!!!!!!!!!!!!cubs!!!!!!!!!!!!!!!!!!!!!!!

           !funw(1)=cbx0
           funw(1)=cbx1
           rhow(1)=0.d0

          do i=2,na1+1
           funw(i)=cubs(i-1)
           rhow(i)=rho(i-1)/roc
          enddo		

           C_bts(1)=funw(1)
           C_bts(iplas)=funw(na1+1)

!           CALL E01BAF(nspl,rhow,funw,RRK,CCK,
!     *                        nspl+4,WRK,6*nspl+16,IFAIL)
!
!           do i=2,iplas-1
!
!              zrho=rosn(i)
!              CALL E02BCF(nspl+4,RRK,CCK,zrho,0,CWk,IFAIL)
!              C_bts(i)=cwk(1)
!
!           enddo

            nspl=na1+1
           call Linsplin(nspl,funw,rhow, iplas,C_bts,rosn)

 !!!!!!!!!!!!!!!!!!!cubs!!!!!!!!!!!!!!!!!!!!!!!
 
 !!!!!!!!!!!!!!!!!!!te!!!!!!!!!!!!!!!!!!!!!!!

           !funw(1)=cex0
           funw(1)=cex1
           rhow(1)=0.d0

          do i=2,na1+1
           funw(i)=te(i-1)
           rhow(i)=rho(i-1)/roc
          enddo		

           T_el(1)=funw(1)
           T_el(iplas)=funw(na1+1)

!           CALL E01BAF(nspl,rhow,funw,RRK,CCK,
!     *                        nspl+4,WRK,6*nspl+16,IFAIL)
!
!           do i=2,iplas-1
!
!              zrho=rosn(i)
!              CALL E02BCF(nspl+4,RRK,CCK,zrho,0,CWk,IFAIL)
!              T_el(i)=cwk(1)
!
!           enddo

            nspl=na1+1
           call Linsplin(nspl,funw,rhow, iplas,T_el,rosn)

 !!!!!!!!!!!!!!!!!!!te!!!!!!!!!!!!!!!!!!!!!!!
  
             do i=1,iplas
                C_bts(i)=C_bts(i)*btor
                C_driv(i)=C_driv(i)*btor
                T_el(i)=T_el(i)*1.d3
              enddo

            return
            end
		 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       subroutine profro_cu_L(na1,neql,rtor,btor,cu,rho,roc)
 
           include 'double.inc'
          parameter(nursp=1000,nursp4=nursp+4,nursp6=nursp4*6)
           include 'dim.inc'
           include 'compol.inc'

         common /com_jb/ BJ_av(nrp),curfi_av(nrp)
          real*8 rosn(nrp),funw(nursp),rhow(nursp)
          real*8 rrk(nursp4),cck(nursp4),wrk(nursp6)
          real*8 cwk(4)

        real*8 rho(*),rtor,btor,roc,cu(*)
        
        iplas=neql
        nspl=na1+1
         do i=1,iplas
          rosn(i)=(i-1.d0)/(iplas-1.d0)
         enddo

!! exstrapolation on magn. axis
       
             rh0=0.d0
		              
             rh1=rho(1)           
             rh2=rho(2)           
             rh3=rho(3)           

             cux1=cu(1)
             cux2=cu(2)
             cux3=cu(3)


       call EXTRP2(rh0,cux0, rh1,rh2,rh3, cux1,cux2,cux3)

!!!!!!!!!!!!!!!!!!!cu!!!!!!!!!!!!!!!!!!!!!!!

           !funw(1)=cux0
           funw(1)=cux1
           rhow(1)=0.d0

          do i=2,na1+1
           funw(i)=cu(i-1)
           rhow(i)=rho(i-1)/roc
          enddo		

           BJ_av(1)=funw(1)
           BJ_av(iplas)=funw(na1+1)

!           CALL E01BAF(nspl,rhow,funw,RRK,CCK,
!     *                        nspl+4,WRK,6*nspl+16,IFAIL)
!
!           do i=2,iplas-1
!
!              zrho=rosn(i)
!              CALL E02BCF(nspl+4,RRK,CCK,zrho,0,CWk,IFAIL)
!              BJ_av(i)=cwk(1)
!
!           enddo

            nspl=na1+1
           call Linsplin(nspl,funw,rhow, iplas,BJ_av,rosn)


 !!!!!!!!!!!!!!!!!!!cu!!!!!!!!!!!!!!!!!!!!!!!
  
             do i=1,iplas
                BJ_av(i)=BJ_av(i)*btor*amu0
              enddo

            return
            end
		 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       subroutine profro_cu(na1,neql,rtor,btor,cu,rho,roc)
 
           include 'double.inc'
          parameter(nursp=1000,nursp4=nursp+4,nursp6=nursp4*6)
           include 'dim.inc'
           include 'compol.inc'

         common /com_jb/ BJ_av(nrp),curfi_av(nrp)
          real*8 rosn(nrp),funw(nursp),rhow(nursp)
          real*8 rrk(nursp4),cck(nursp4),wrk(nursp6)
          real*8 cwk(4)

        real*8 rho(*),rtor,btor,roc,cu(*)
        
        iplas=neql
        nspl=na1+1
         do i=1,iplas
          rosn(i)=(i-1.d0)/(iplas-1.d0)
         enddo

!! exstrapolation on magn. axis
       
             rh0=0.d0
		              
             rh1=rho(1)           
             rh2=rho(2)           
             rh3=rho(3)           

             cux1=cu(1)
             cux2=cu(2)
             cux3=cu(3)


       call EXTRP2(rh0,cux0, rh1,rh2,rh3, cux1,cux2,cux3)

!!!!!!!!!!!!!!!!!!!cu!!!!!!!!!!!!!!!!!!!!!!!

           funw(1)=cux0
           rhow(1)=0.d0

          do i=2,na1+1
           funw(i)=cu(i-1)
           rhow(i)=rho(i-1)/roc
          enddo		

           BJ_av(1)=funw(1)
           BJ_av(iplas)=funw(na1+1)

           CALL E01BAF(nspl,rhow,funw,RRK,CCK,
     *                        nspl+4,WRK,6*nspl+16,IFAIL)

           do i=2,iplas-1

              zrho=rosn(i)
              CALL E02BCF(nspl+4,RRK,CCK,zrho,0,CWk,IFAIL)
              BJ_av(i)=cwk(1)

           enddo

 !!!!!!!!!!!!!!!!!!!cu!!!!!!!!!!!!!!!!!!!!!!!
  
             do i=1,iplas
                BJ_av(i)=BJ_av(i)*btor*amu0
              enddo

            return
            end
		 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       subroutine profro_pres(na1,neql,pres,rho,roc)
       
           include 'double.inc'
          parameter(nursp=1000,nursp4=nursp+4,nursp6=nursp4*6)
           include 'dim.inc'
           include 'compol.inc'
	    common /com_pres/ P_rho(0:nrp),dpdro(nrp),dPdFi(1:nrp),romin

          real*8 rocn(nrp),funw(nursp),rhow(nursp)
          real*8 rrk(nursp4),cck(nursp4),wrk(nursp6)
          real*8 cwk(4)

        real*8 rho(*),roc,pres(*)
        real*8 rokn(nrp)

        iplas=neql
        nspl=na1+1

         do i=1,iplas-1
          rocn(i)=(i-0.5d0)/(iplas-1.d0)
         enddo
         
          rocn(iplas)=1.d0

         do i=1,iplas
          rokn(i)=((i-1.d0)/(iplas-1.d0))
         enddo
         
          rokn(iplas)=1.d0
          rokn(1)=0.d0
          
!! exstrapolation on magn. axis
       
             rh0=0.d0
		              
             rh1=rho(1)           
             rh2=rho(2)           
             rh3=rho(3)           

             px1=pres(1)
             px2=pres(2)
             px3=pres(3)


       call EXTRP2(rh0,px0, rh1,rh2,rh3, px1,px2,px3)

!!!!!!!!!!!!!!!!!!!p!!!!!!!!!!!!!!!!!!!!!!!

           funw(1)=px0
           !funw(1)=pres(1)
           rhow(1)=0.d0

          do i=2,na1+1
           funw(i)=pres(i-1)
           rhow(i)=rho(i-1)/roc
          enddo		

           P_rho(0)=funw(1)
           P_rho(iplas)=funw(na1+1)

           CALL E01BAF(nspl,rhow,funw,RRK,CCK,
     *                        nspl+4,WRK,6*nspl+16,IFAIL)

           do i=1,iplas-1

              zrho=rocn(i)
              CALL E02BCF(nspl+4,RRK,CCK,zrho,0,CWk,IFAIL)
              P_rho(i)=cwk(1)

           enddo
           
           do i=2,na1+1

              zrho=rhow(i)
              CALL E02BCF(nspl+4,RRK,CCK,zrho,0,CWk,IFAIL)
              dPdFi(i)=cwk(2)/zrho

           enddo
           do i=1,iplas

              zrho=rokn(i)
            if(zrho.lt.rhow(2)) zrho=rhow(2)
              CALL E02BCF(nspl+4,RRK,CCK,zrho,0,CWk,IFAIL)
              dPdro(i)=0.5d0*cwk(2)/zrho

           enddo
           
            romin=rhow(2)

!!!!!!!!!!!!!!!!!!!p!!!!!!!!!!!!!!!!!!!!!!!
  
  
!             do i=1,iplas
!                C_bts(i)=C_bts(i)*btor
!                C_driv(i)=C_driv(i)*btor
!                T_el(i)=T_el(i)*1.d3
!              enddo

            return
            end
!!!!!!!!!!!!!!!!!!!p!!!!!!!!!!!!!!!!!!!!!!!
       subroutine proffi_pres(na1,neql,pres,rho,roc)
       
           include 'double.inc'
          !parameter(nursp=1000,nursp4=nursp+4,nursp6=nursp4*6)
           include 'dim.inc'
           include 'compol.inc'
	    common /com_pres/ P_rho(0:nrp),dpdro(nrp),dPdFi(1:nrp),romin

           real*8 rocn(nrp) !,funw(nursp),rhow(nursp)
          !real*8 rrk(nursp4),cck(nursp4),wrk(nursp6)
          !real*8 cwk(4)

        real*8 rho(*),roc,pres(*)
        real*8 rokn(nrp)
        
        real(8), allocatable :: PXIN(:), PYIN(:), PXOUT(:), PYOUT(:),
     +  PYOUTP(:), PYOUTPP(:), PY2(:), PWORK(:),
     +  PAMAT(:,:), PYINNEW(:)

        real(8), allocatable :: p05_prim(:),rho05(:),d2pf(:)

           allocate(rho05(na1+1))
           allocate(p05_prim(na1+1))
           allocate(d2pf(neql))
           
         do i=2,na1
          rho05(i)=0.5d0*( rho(i-1)**2+rho(i)**2 )/rho(na1)**2
         enddo
          rho05(1)=0.d0
          rho05(na1+1)=1.d0
          
         do i=2,na1
          p05_prim(i)=( pres(i)-pres(i-1) )
     &               /(rho(i)**2-rho(i-1)**2)*rho(na1)**2
         enddo

!! exstrapolation on magn. axis
       
             rh0=0.d0
		              
             rh1=rho05(2)           
             rh2=rho05(3)           
             rh3=rho05(4)           

             px1=p05_prim(2)
             px2=p05_prim(3)
             px3=p05_prim(4)

       call EXTRP2(rh0,px0, rh1,rh2,rh3, px1,px2,px3) 
       
             p05_prim(1)=px0

!! exstrapolation to boundry

             rh0=1.d0
		              
             rh1=rho05(na1)           
             rh2=rho05(na1-1)           
             rh3=rho05(na1-2)           

             px1=p05_prim(na1)  
             px2=p05_prim(na1-1)
             px3=p05_prim(na1-2)

       call EXTRP2(rh0,px0, rh1,rh2,rh3, px1,px2,px3) 
       
             p05_prim(na1+1)=px0

        iplas=neql
        nspl=na1+1

         do i=1,iplas-1
          rocn(i)=((i-0.5d0)/(iplas-1.d0))**2
         enddo
         
          rocn(iplas)=1.d0
          
         do i=1,iplas
          rokn(i)=((i-1.d0)/(iplas-1.d0))**2
         enddo
         
          rokn(iplas)=1.d0
          rokn(1)=0.d0



           KNIN = nspl
           KNOUT = iplas
           KOPT = 1
           PTAUS = 1d-6
           allocate(PYINNEW(KNIN))
           allocate(PY2(KNIN))
           allocate(PWORK(KNIN))
           MDAMAT = 7
           PBCLFT = 2.d0
           PBCRGT = 2.d0
           allocate(PAMAT(MDAMAT,KNIN))
         CALL INTRPTAU(rho05,p05_prim,PYINNEW,PY2,KNIN,rokn,dPdFi,d2pf,
     +        PYOUTPP,KNOUT,KOPT,PTAUS,PWORK,PAMAT,MDAMAT,PBCLFT,PBCRGT)

!!!!!!!!!!!!!!!!!!!p!!!!!!!!!!!!!!!!!!!!!!!

            return
            end
		 
!!!!!!!!!!!!!!!!!!!p!!!!!!!!!!!!!!!!!!!!!!!
       subroutine S_proffi_pres(equil_in,parameters_spider)
       use imas_ids       
       use parameters_a2spider , only: type_parameters, TWOPI

       include 'double.inc'
       include 'dim.inc'
       include 'compol.inc'
       
       type(type_equilibrium) equil_in
       type(type_parameters) parameters_spider
       
	 common /com_pres/ P_rho(0:nrp),dpdro(nrp),dPdFi(1:nrp),romin

!       real*8 rocn(nrp)
       real*8 rokn(nrp)
        
        real(8), allocatable :: PXIN(:), PYIN(:), PXOUT(:), PYOUT(:),
     +  PYOUTP(:), PYOUTPP(:), PY2(:), PWORK(:),
     +  PAMAT(:,:), PYINNEW(:)

        real(8), allocatable :: p05_prim(:),rho05(:),d2pf(:)

        iplas=parameters_spider%neql
        na1=size(equil_in%profiles_1d%psi)
        nspl=na1+1

        allocate(rho05(na1+1))
        allocate(p05_prim(na1+1))
        allocate(d2pf(iplas))
           
         do i=2,na1
          rho05(i)=0.5d0*( 
     +      equil_in%profiles_1d%rho_tor(i-1)**2 
     +    + equil_in%profiles_1d%rho_tor(i)**2 
     +             )/equil_in%profiles_1d%rho_tor(na1)**2
         enddo
          rho05(1)=0.d0
          rho05(na1+1)=1.d0
          
         do i=2,na1
          p05_prim(i)=(equil_in%profiles_1d%pressure(i)
     +                 -equil_in%profiles_1d%pressure(i-1))
     +               /(equil_in%profiles_1d%rho_tor(i)**2
     +                 -equil_in%profiles_1d%rho_tor(i-1)**2)
     +               *equil_in%profiles_1d%rho_tor(na1)**2
         enddo
         prescale = 1.d-6
         do i=2,na1
          p05_prim(i)=prescale*p05_prim(i)
         enddo


!! exstrapolation on magn. axis
       
             rh0=0.d0
		              
             rh1=rho05(2)           
             rh2=rho05(3)           
             rh3=rho05(4)           

             px1=p05_prim(2)
             px2=p05_prim(3)
             px3=p05_prim(4)

       call EXTRP2(rh0,px0, rh1,rh2,rh3, px1,px2,px3) 
       
             p05_prim(1)=px0

!! exstrapolation to boundry

             rh0=1.d0
		              
             rh1=rho05(na1)           
             rh2=rho05(na1-1)           
             rh3=rho05(na1-2)           

             px1=p05_prim(na1)  
             px2=p05_prim(na1-1)
             px3=p05_prim(na1-2)

       call EXTRP2(rh0,px0, rh1,rh2,rh3, px1,px2,px3) 
       
             p05_prim(na1+1)=px0

!         do i=1,iplas-1
!          rocn(i)=((i-0.5d0)/(iplas-1.d0))**2
!         enddo
!         
!          rocn(iplas)=1.d0
          
         do i=1,iplas
          rokn(i)=((i-1.d0)/(iplas-1.d0))**2
         enddo
         
          rokn(iplas)=1.d0
          rokn(1)=0.d0

           KNIN = nspl
           KNOUT = iplas
           KOPT = 1
           PTAUS = 1d-6
           allocate(PYINNEW(KNIN))
           allocate(PY2(KNIN))
           !allocate(PYOUTPP(KNOUT))
           allocate(PWORK(KNIN))
           MDAMAT = 7
           PBCLFT = 2.d0
           PBCRGT = 2.d0
           allocate(PAMAT(MDAMAT,KNIN))
         CALL INTRPTAU(rho05,p05_prim,PYINNEW,PY2,KNIN,rokn,dPdFi,d2pf,
     +        PYOUTPP,KNOUT,KOPT,PTAUS,PWORK,PAMAT,MDAMAT,PBCLFT,PBCRGT)

!!!!!!!!!!!!!!!!!!!p!!!!!!!!!!!!!!!!!!!!!!!

            return
            end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       subroutine proffi_pres_(na1,neql,pres,rho,roc)
       
           include 'double.inc'
          parameter(nursp=1000,nursp4=nursp+4,nursp6=nursp4*6)
           include 'dim.inc'
           include 'compol.inc'
	    common /com_pres/ P_rho(0:nrp),dpdro(nrp),dPdFi(1:nrp),romin

          real*8 rocn(nrp),funw(nursp),rhow(nursp)
          real*8 rrk(nursp4),cck(nursp4),wrk(nursp6)
          real*8 cwk(4)

        real*8 rho(*),roc,pres(*)
        real*8 rokn(nrp)

        iplas=neql
        nspl=na1+1

         do i=1,iplas-1
          rocn(i)=((i-0.5d0)/(iplas-1.d0))**2
         enddo
         
          rocn(iplas)=1.d0
          
         do i=1,iplas
          rokn(i)=((i-1.d0)/(iplas-1.d0))**2
         enddo
         
          rokn(iplas)=1.d0
          rokn(1)=0.d0

!! exstrapolation on magn. axis
       
             rh0=0.d0
		              
             rh1=rho(1)           
             rh2=rho(2)           
             rh3=rho(3)           

             px1=pres(1)
             px2=pres(2)
             px3=pres(3)


       call EXTRP2(rh0,px0, rh1,rh2,rh3, px1,px2,px3)

!!!!!!!!!!!!!!!!!!!p!!!!!!!!!!!!!!!!!!!!!!!

           funw(1)=px0
           rhow(1)=0.d0

          do i=2,na1+1
           funw(i)=pres(i-1)
           rhow(i)=(rho(i-1)/roc)**2
          enddo		

           P_rho(0)=funw(1)
           P_rho(iplas)=funw(na1+1)

           CALL E01BAF(nspl,rhow,funw,RRK,CCK,
     *                        nspl+4,WRK,6*nspl+16,IFAIL)

           do i=1,iplas-1

              zrho=rocn(i)
              CALL E02BCF(nspl+4,RRK,CCK,zrho,0,CWk,IFAIL)
              P_rho(i)=cwk(1)

           enddo
           
           do i=1,iplas

              zrho=rokn(i)
              CALL E02BCF(nspl+4,RRK,CCK,zrho,0,CWk,IFAIL)
              dPdFi(i)=cwk(2)

           enddo

!!!!!!!!!!!!!!!!!!!p!!!!!!!!!!!!!!!!!!!!!!!

            return
            end
		 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       subroutine profro_pres_L(na1,neql,pres,rho,roc)
       
           include 'double.inc'
          parameter(nursp=1000)
           include 'dim.inc'
           include 'compol.inc'
	    common /com_pres/ P_rho(0:nrp),dpdro(nrp),dPdFi(1:nrp),romin

          real*8 rocn(nrp),funw(nursp),rhow(nursp)
          real*8 rho(*),roc,pres(*)

        iplas=neql
        nspl=na1+1

         do i=1,iplas-1
          rocn(i)=(i-0.5d0)/(iplas-1.d0)
         enddo
         
          rocn(iplas)=1.d0

!! exstrapolation to magn. axis
       
             rh0=0.d0
		              
             rh1=rho(1)           
             rh2=rho(2)           
             rh3=rho(3)           

             px1=pres(1)
             px2=pres(2)
             px3=pres(3)


       call EXTRP2(rh0,px0, rh1,rh2,rh3, px1,px2,px3)

!!!!!!!!!!!!!!!!!!!p!!!!!!!!!!!!!!!!!!!!!!!

           funw(1)=px0
           rhow(1)=0.d0

          do i=2,na1+1
           funw(i)=pres(i-1)
           rhow(i)=rho(i-1)/roc
          enddo		

           P_rho(0)=funw(1)
           P_rho(iplas)=funw(na1+1)

           do i=1,iplas-1

              zrho=rocn(i)
            do j=1,na1
             if( zrho.gt.rhow(j) .AnD. zrho.le.rhow(j+1) ) then
              Pres_rho=( funw(j)*(rhow(j+1)-zrho) +
     &                   funw(j+1)*(zrho-rhow(j)) ) 
     &                  /( rhow(j+1)-rhow(j) ) 
              go to 100             
             endif
            enddo

 100        P_rho(i)=Pres_rho

           enddo

!!!!!!!!!!!!!!!!!!!p!!!!!!!!!!!!!!!!!!!!!!!
  
            return
            end
		 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       subroutine pres_d_psi
       
           include 'double.inc'
           include 'dim.inc'
          parameter(nrsp=nrp+1,nrsp4=nrsp+4,nrsp6=nrsp4*6)
           include 'compol.inc'
	    common /com_pres/ P_rho(0:nrp),dpdro(nrp),dPdFi(1:nrp),romin

          real*8 rhocn(nrp+1),funw(nrsp),rhokn(nrp)
          real*8 rrk(nrsp4),cck(nrsp4),wrk(nrsp6)
          real*8 cwk(4)


         do i=1,iplas-1
          rhocn(i+1)=(i-0.5d0)/(iplas-1.d0)
         enddo
         
          rhocn(iplas+1)=1.d0
          rhocn(1)=0.d0
          
         do i=1,iplas
          rhokn(i)=(i-1.d0)/(iplas-1.d0)
         enddo

!! exstrapolation q to magn. axis
       
             rh0=0.d0
		              
             rh1=0.5d0           
             rh2=1.5d0           
             rh3=2.5d0           

             qx1=q(1)
             qx2=q(2)
             qx3=q(3)

       call EXTRP2(rh0,qx0, rh1,rh2,rh3, qx1,qx2,qx3)

          !dpdpsi(1)=-amu0*qx0*dpdro(1)/flucfm
          dpdpsi(1)=-amu0*qx0*dpdfi(1)/flucfm
          funw(1)=qx0

         do i=1,iplas
          funw(i+1)=q(i)
         enddo

!!!!!!!!!!!!!!!!!!!!!!! q !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         n3spl=iplas+1
        CALL E01BAF(n3spl,rhocn,funw,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=2,iplas
		       zrho=rhokn(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zrho,0,CWk,IFAIL)
               qkn=CWk(1)
               !dpdpsi(i)=-amu0*qkn*dpdro(i)/flucfm
               dpdpsi(i)=-amu0*qkn*dpdfi(i)/flucfm
            !if(zrho.lt.romin) dpdpsi(i)=dpdpsi(1)
           enddo
               !dpdpsi(1)=dpdpsi(2)


            return
            end
		 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!********************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         subroutine get_rz(rout,zout,ni_p,nj_p)

           include 'double.inc'
           include 'dim.inc'
           include 'compol.inc'

          real*8 rout(ni_p,nj_p),zout(ni_p,nj_p)

        do i=1,iplas
        do j=1,nt
         rout(i,j)=r(i,j)
         zout(i,j)=z(i,j)
        enddo
        enddo

      return
      end


	subroutine	NEWGRD_ss( ROC, HRO, NB1, NA1, RHO,btor, HROA  )
C----------------------------------------------------------------------|
C Define size of the edge cell to be:   .6 <= HROA/HRO < 1.8
C	with a hysteresis of 0.2*HRO, so that
C		HROA<=0.6*HRO  jumps to  HROA<=1.6*HRO
C		HROA>=1.8*HRO  jumps to  HROA>=0.8*HRO
C	write(*,*)NA1,TIME,RHO(NA1)," -> ",ROC
C----------------------------------------------------------------------|
         include 'double.inc'
         include 'dim.inc'
         include 'compol.inc'
         dimension RHO(*)

           na=na1-1
          roc=sqrt(flx_fi(iplas)/(pi*btor))

C	write(*,*)NA1,ROC,HRO*(NA+0.1),HRO*(NA1+0.3),HRO*(NB1-0.49)
C	write(*,'(1P,6E13.5)')(RHO(j),j=NA-2,NA1+1)
C	write(*,'(1P,6E13.5)')(AMETR(j),j=NA-2,NA1+1)
	if (ROC .gt. HRO*(NB1-0.49))	then
C rho_edge gets out of the grid
	   write(*,*)'>>> WARNING: the allocated grid is too small'
	   write(*,*)'Time =',TIME,'   Rho_b =',ROC,
     >		'   Rho_max = ',HRO*(NB1-0.5)
	   write(*,*)'Try to increase AWALL'
	   write(*,*)'If this does not help inspect equilibrium input'
	   NA = NB1-1
	elseif (ROC .le. HRO*(NA+0.1))	then	! 0.1 <=> 0.6-0.5
!	    if (HROA < 0.6*HRO) then   reduce NA
C	   write(*,*)" <- ",ROC,HRO*NA,RHO(NA1),HRO*NA1,NA1
	   do	j=NA-1,1,-1
		NA = j
C	   	write(*,*)j,NA*HRO,ROC,(NA+1)*HRO,ROC/HRO-(NA-0.5)
		if (ROC .gt. HRO*(j+0.1))	goto	10
	   enddo
	elseif (ROC .gt. HRO*(NA1+0.3))	then	! 0.3 <=> 1.8-1.5
!	    if (HROA > 1.8*HRO) then   increase NA
C	   write(*,*)" -> ",ROC,HRO*NA,RHO(NA1),HRO*NA1,NA1
	   do	j=NA1,NB1
		NA = j
C	   	write(*,*)j,NA*HRO,ROC,(NA+1)*HRO
		if (ROC .le. HRO*(j+1.3))	goto	10
	   enddo
	endif
 10	continue

C HRO*(NA+0.6) < ROC <= HRO*(NA+1.8)
C	write(*,*)AB,ABC,TIME,TSTART
C	write(*,*)NA+1,NA+1-NA1,HRO*(NA+0.6),ROC,HRO*(NA+1.8)

	NA1 = NA+1
	do	J=1,NB1
	   RHO(J)=(J-0.5)*HRO
	enddo
	HROA = ROC-RHO(NA)
	RHO(NA1) = ROC
	
      return
	end
C======================================================================|

	subroutine	get_roc(ROC,btor)

c	implicit none
         include 'double.inc'
         include 'dim.inc'
         include 'compol.inc'

	 real*8 ROC,btor
          roc=sqrt(flx_fi(iplas)/(pi*btor))
 
	return
	end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine	get_phi(phi)

c	implicit none
         include 'double.inc'
         include 'dim.inc'
         include 'compol.inc'

	 real*8 phi(*)
	 integer i
	 do i=1,iplas
          phi(i)=flx_fi(i)
       enddo
 
	return
	end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         subroutine pla_volt(dt)
        include 'double.inc'
        include 'dim.inc'
        include 'compol.inc'
          common /com_volt/ upls(nrp)
          common /com_psisave/ psi_0(nrp)
 
          call get_psibon(psi_bn1)
          
        do i=1,iplas
         !upls(i)=(psia(i)*psim-psi_0(i))/dt
         upls(i)=(psi(i,2)+psi_bn1-psi_0(i))/dt
         enddo
C	write(*,*) 'uplss ' , upls(iplas),psi(iplas,2)
C	write(*,*) psi_bn1,psi_0(iplas),dt          
        return
        end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         subroutine savepsi
        include 'double.inc'
        include 'dim.inc'
        include 'compol.inc'
          common /com_psisave/ psi_0(nrp)
          
         call get_psibon(psi_bn1)
          
         do i=1,iplas
          !psi_0(i)=psim*psia(i)
          psi_0(i)=psi(i,2)+psi_bn1
         enddo
          
        return
        end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         subroutine get_psibon(psi_bnd)
        include 'double.inc'
        include 'dim.inc'
      common/selcon/ psi_d(nrp),fi_d(nrp),f_d(nrp),ri_d(nrp),
     *               ps_pnt(nrp),del_psb,psi_bn1
          
          psi_bnd=psi_bn1
          
        return
        end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         subroutine flux_state(flx_st)
          include 'double.inc'
          include 'param.inc'
          include 'comblc.inc'
          

            Sj=0.d0
            Spsj=0.d0

          do j=1,nj
          do i=1,ni
          
           if(ipr(i,j).eq.1) then
            Sj=Sj+curf(i,j)
            Spsj=Spsj+curf(i,j)*ue(i,j)
           endif
          
          enddo
          enddo
          
            flx_st=Spsj/Sj
          
        return
        end
	
         subroutine Linsplin(nspl,usp,xsp, n,u,x)
        include 'double.inc'
        real*8 usp(nspl),xsp(nspl)
        real*8 u(n),x(n)
        
         u(1)=usp(1)
         u(n)=usp(nspl)
        
        
           do i=2,n-1
              zx=x(i)

            do ic=1,nspl-1
             if(zx.le.xsp(ic+1) .AnD. zx.gt.xsp(ic)) then
                   u(i)=( usp(ic)*(xsp(ic+1)-zx)
     *                   + usp(ic+1)*(zx-xsp(ic)) )
     *                      /(xsp(ic+1)-xsp(ic))
              exit             
             endif             
            enddo

           enddo
                 
        return
        end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       subroutine prof_AST_SP_L(na1,neql,ArrAS,ArrSP,rho,roc)
 
!!!!!ASTRA - SPIDER arrays transformation 
!!!!!ArrAS - ASTRA array(input) 
!!!!!ArrSP - SPIDER array(output)
           include 'double.inc'
           include 'dim.inc'
           include 'compol.inc'

          parameter(nursp=1000,nursp4=nursp+4,nursp6=nursp4*6)

          real*8 rosn(nrp),funw(nursp),rhow(nursp)
          real*8 rrk(nursp4),cck(nursp4),wrk(nursp6)
          real*8 cwk(4)

        real*8 rho(*),ArrAS(*),ArrSP(*),roc
        
        iplas=neql
        nspl=na1+1
         do i=1,iplas
          rosn(i)=(i-1.d0)/(iplas-1.d0)
         enddo

!! exstrapolation on magn. axis
       
             rh0=0.d0
		              
             rh1=rho(1)           
             rh2=rho(2)           
             rh3=rho(3)           

             cux1=ArrAS(1)
             cux2=ArrAS(2)
             cux3=ArrAS(3)


       call EXTRP2(rh0,cux0, rh1,rh2,rh3, cux1,cux2,cux3)

!!!!!!!!!!!!!!!!!!!cu!!!!!!!!!!!!!!!!!!!!!!!

           funw(1)=cux0
           !funw(1)=cux1
           rhow(1)=0.d0

          do i=2,na1+1
           funw(i)=ArrAS(i-1)
           rhow(i)=rho(i-1)/roc
          enddo		

           ArrSP(1)=funw(1)
           ArrSP(iplas)=funw(na1+1)

!           CALL E01BAF(nspl,rhow,funw,RRK,CCK,
!     *                        nspl+4,WRK,6*nspl+16,IFAIL)
!
!           do i=2,iplas-1
!
!              zrho=rosn(i)
!              CALL E02BCF(nspl+4,RRK,CCK,zrho,0,CWk,IFAIL)
!              BJ_av(i)=cwk(1)
!
!           enddo

            nspl=na1+1
           call Linsplin(nspl,funw,rhow, iplas,ArrSP,rosn)


            return
            end
		 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine input2spider(nbnd,rzbnd,na1,eqpf,eqff,fp,ipl,
     *                    rtor,btor,rho,pres,cu,
     *                    nstep, key_equil,             
     *                    equil_in,params)
       
      use imas_ids       
      use parameters_a2spider, only: type_parameters, TWOPI
 
      implicit none
 
      
      integer nbnd, key_equil, na1, nstep
      real(8) dt,time,dpsdt,yreler
      real(8) ipl, rtor, btor,roc
      real(8) rzbnd(1:2*nbnd), 
     *        eqpf(na1), eqff(na1), fp(na1),
     *        rho(na1), pres(na1), cu(na1)
     
     
      type(type_equilibrium) equil_in
      type(type_parameters) params
      
      !locals
      integer neql,nteta
      integer kpr,k_con,kname
      integer k_fixfree,key_ini,key_start,key_0stp,key_pres,key_dmf
      integer k_grid,k_auto
      real(8), allocatable :: mu(:), cc(:),Te(:),cubc(:),cd(:)
      character*80 eqdfn
      character*80 prename
      integer i
      
      allocate( mu(na1), cc(na1), Te(na1), cubc(na1), cd(na1) )
      cc = equil_in%profiles_1d%sigmapar%value
      cubc = equil_in%profiles_1d%jni%value
      Te = equil_in%profiles_1d%te%value
      cc = equil_in%profiles_1d%sigmapar%value
      cd = 0
      do i =1,na1
        mu(i) = 1.d0/equil_in%profiles_1d%q(i)
      enddo

      !setup parameters for equilibrium problem to solve
      params%nstep = nstep
      if(nstep .ne. 0) then
        params%key_dmf = key_equil
      else
        if(key_equil .eq. 0) then
            params%key_0stp = key_equil
        elseif(key_equil .eq. -10) then  
            params%key_0stp = 1
        else
            params%key_0stp = key_equil
        endif
      endif
      !p' input for key_equil=0 (otherwise the default key_pres is used)
      if(key_equil .eq. 0) then
         params%key_pres = 0
      endif

      kpr=params%kpr     !print in spider
      call  kpr_calc(kpr)
      k_con=0      
      call put_key_con(k_con)
      prename=params%prename
      kname=params%kname      
      call  put_name(prename,kname)

      nstep = params%nstep
      time = params%time
      dt = params%dt
      key_dmf = params%key_dmf
      k_grid = params%k_grid  ! k_grid= 0   rect. grid
                              ! k_grid= 1   adap. grid
      k_auto = params%k_auto  ! k_auto= 1->   full initialization,
                              ! k_auto= 0-> preinitialization is assumed to be done
      
      k_fixfree = params%k_fixfree   !=0->only fixed boundary spider 
      dpsdt = params%dpsdt
      key_ini = params%key_ini       ! =1 astra profiles, =0 start from EQDSK and SPIDER profiles
      key_start=params%key_start     ! 
      eqdfn=params%eqdfn
      key_0stp=params%key_0stp   
      key_pres=params%key_pres   
      neql=params%neql   
      nteta=params%nteta  


        call aspid_flag(1)
        !call put_key_fix(k_fixfree)
        roc=rho(na1)
             
        call astra2spider(neql,nteta,nbnd,rzbnd,key_dmf,
     *                    na1,eqpf,eqff,fp,ipl,
     *                    rtor,btor,rho,roc,nstep,yreler,mu,
     *                    cc,Te,cubc,cd,key_ini,eqdfn,pres,cu,
     *                    key_0stp,key_pres                    )

      
      return
      end   
             
      subroutine spider2output(
     *                       equil,params,
     *                       rtor,btor,rho,roc,na1,
     *                       g11,g22,g33,vr,vrs,slat,gradro,rocnew,
     *                       mu,ipol,bmaxt,bmint,bdb02,b0db2,bdb0,droda,
     *				         yreler,yli3,platok,cu_out,fp,
     *                       pres_s )

      use imas_ids       
      use parameters_a2spider, only: type_parameters, TWOPI
      
      implicit none
      
      type(type_equilibrium) equil
      type(type_parameters) params
      
      integer na1,ni_p,nj_p
      real(8) rtor,btor, platok,yli3,roc,rocnew,yreler
      real(8)  
     *        g11(na1), g22(na1), g33(na1),
     *        rho(na1), pres_s(na1), cu_out(na1),fp(na1),
     *        vr(na1),vrs(na1),slat(na1),
     *        gradro(na1),droda(na1),mu(na1),     
     *        ipol(na1),bmaxt(na1),bmint(na1),      
     *        bdb02(na1),b0db2(na1),bdb0(na1)      
      
      !locals
      real(8), allocatable :: W_Dj(:), yFOFB(:)
      integer i, j, neql, nteta
      real(8), allocatable :: rout(:,:), zout(:,:)
      
      allocate( W_Dj(na1), yFOFB(na1) )
      ni_p = params%neql
      nj_p = params%nteta
      allocate( rout(ni_p,nj_p), zout(ni_p,nj_p) )
      
      call spider2astra(rout,zout,rtor,btor,rho,roc,na1,
     *                    g11,g22,g33,vr,vrs,slat,gradro,rocnew,
     *                    mu,ipol,bmaxt,bmint,bdb02,b0db2,bdb0,droda,
     *				      yreler,yli3,ni_p,nj_p,platok,cu_out,fp,
     *                    pres_s,W_Dj,yFOFB )
 
!        !equil_out
!	  neql = params%neql
!	  nteta = params%nteta
!        !grid coordinates into coord_sys
!        allocate( equil%coord_sys%position%r(neql,nteta-2) )
!        allocate( equil%coord_sys%position%z(neql,nteta-2) )
!        do i=1,neql
!            do j=1,nteta-2
!                equil%coord_sys%position%r(i,j)
!     &            = rout(i,j)
!                equil%coord_sys%position%z(i,j)
!     &            = zout(i,j)
!            enddo
!        enddo
      
      return
      end
      
	subroutine	get_eq(equil_out, psi_boundary,press0,btor)
      use imas_ids
			use parameters_a2spider, only: type_parameters, TWOPI
   	!implicit none
      include 'double.inc'
      include 'dim.inc'
      include 'compol.inc'
      type(type_equilibrium) equil_out
      real*8 psi_boundary,btor
      real*8 pscale, fscale
      real*8 um, up,press0
      real*8, allocatable :: wrk(:)
      real*8, pointer :: rho_p(:)

	common /com_sigcd/ C_sig(nrp),T_el(nrp),C_bts(nrp),C_driv(nrp)
      common /com_jb/ BJ_av(nrp),curfi_av(nrp)
     
      pscale = 2*TWOPI*1.d-7 * TWOPI
      fscale = TWOPI
      
      if(.NOT. associated(equil_out%profiles_1d%psi)) then
          allocate(equil_out%profiles_1d%psi(iplas))
          allocate(equil_out%profiles_1d%pressure(iplas))
          allocate(equil_out%profiles_1d%phi(iplas))
          allocate(equil_out%profiles_1d%pprime(iplas))
          allocate(equil_out%profiles_1d%ffprime(iplas))
          allocate(equil_out%profiles_1d%F_dia(iplas))
          allocate(equil_out%profiles_1d%q(iplas))
       endif
          
      if(.NOT. associated(equil_out%coord_sys%position%r)) then
          allocate(equil_out%coord_sys%position%r(iplas,nt1))
          allocate(equil_out%coord_sys%position%z(iplas,nt1))    
          allocate(equil_out%coord_sys%position%teta2d(nt1))    
          allocate(equil_out%coord_sys%position%rmin(iplas,nt1))    
          allocate(equil_out%coord_sys%position%psirz(iplas,nt1))    
      endif
      
      do i=1,iplas
      
        ! psim, psip for fixed boundary
        !equil_out%profiles_1d%psi(i) =((psim-psip)*psia(i)+psip)*2.d0*PI
        equil_out%profiles_1d%psi(i) = 
     &                       ((psim-psip)*psia(i)+psi_boundary)*TWOPI
     
        equil_out%profiles_1d%phi(i) = flx_fi(i)
        equil_out%profiles_1d%pprime(i) = -dpdpsi(i)/pscale
        equil_out%profiles_1d%ffprime(i) = -dfdpsi(i)/fscale
        !equil_out%profiles_1d%F_dia(i) = -f(i)
        !equil_out%profiles_1d%q(i) = q(i)/TWOPI
        do j=1,nt1
          equil_out%coord_sys%position%r(i,j) = r(i,j+1)
          equil_out%coord_sys%position%z(i,j) = z(i,j+1)
          equil_out%coord_sys%position%teta2d(j) = teta(j+1)
          equil_out%coord_sys%position%rmin(i,j) = ro(i,j+1)
          equil_out%coord_sys%position%psirz(i,j) = -psi(i,j+1)
        enddo     
      enddo

      equil_out%global_param%toroid_field%r0 = r0ax
      equil_out%global_param%toroid_field%b0 = b0ax
      equil_out%global_param%i_plasma = tokp*1e6
         
      if(.NOT. associated(equil_out%eqgeometry%boundary%r)) then
          allocate(equil_out%eqgeometry%boundary%r(nt1-1))
          allocate(equil_out%eqgeometry%boundary%z(nt1-1))
      endif
      equil_out%eqgeometry%boundary%npoints = nt1-1    
      equil_out%eqgeometry%boundary%r = 
     &  equil_out%coord_sys%position%r(iplas,1:nt1-1)
      equil_out%eqgeometry%boundary%z = 
     &  equil_out%coord_sys%position%z(iplas,1:nt1-1)

 !pressure     
      
          allocate(wrk(0:iplas+1))
          
       wrk(iplas)=0.d0
       equil_out%profiles_1d%pressure(iplas) =
     & wrk(iplas)/2.d0/TWOPI*1.d7
      do i=iplas-1,1,-1 
       dpsi=(psim-psip)*(psia(i+1)-psia(i))
       wrk(i)=wrk(i+1)-0.5d0*(dpdpsi(i+1)+dpdpsi(i))*dpsi
       equil_out%profiles_1d%pressure(i) = 
     & wrk(i)/2.d0/TWOPI*1.d7+press0
      enddo
	deallocate(wrk)
      
      if(.NOT. associated(equil_out%profiles_1d%rho_tor)) then
         allocate( equil_out%profiles_1d%rho_tor(iplas) )
         allocate( equil_out%profiles_1d%jparallel(iplas) )
         allocate( equil_out%profiles_1d%sigmapar%value(iplas) )
         allocate( equil_out%profiles_1d%jni%value(iplas) )
         allocate( equil_out%profiles_1d%te%value(iplas) )
      endif
      do i=1,iplas
         equil_out%profiles_1d%rho_tor(i) = 
     &                         sqrt( flx_fi(i)/(0.5*TWOPI*b0ax) )
         equil_out%profiles_1d%jparallel(i) = BJ_av(i)/b0ax/(.2d0*TWOPI)
         equil_out%profiles_1d%sigmapar%value(i) = C_sig(i)
         equil_out%profiles_1d%jni%value(i) = C_bts(i)/b0ax
         equil_out%profiles_1d%te%value(i) = T_el(i)*1d-3
      enddo

           rho_p=>equil_out%profiles_1d%rho_tor

          allocate(wrk(iplas))
          
	    call cell2node(q,wrk,rho_p,iplas)
      do i=1,iplas
        equil_out%profiles_1d%q(i) = wrk(i)/TWOPI
      enddo
	    call cell2node(f,wrk,rho_p,iplas)
      do i=1,iplas
        equil_out%profiles_1d%F_dia(i) = -wrk(i)
      enddo
         
	    deallocate(wrk)

	return
	end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	subroutine	cell2node(arr_c,arr_n,rho,iplas)
 	include 'double.inc'
	real*8 arr_c(*),arr_n(*),rho(*)
      real*8, allocatable :: x(:)
      real*8, allocatable :: x05(:)
      real*8, allocatable :: RRK(:),CCK(:),WRK(:),CWK(:)
      
         allocate(RRK(iplas+4),CCK(iplas+4),WRK((iplas+4)*6))
         allocate(CWK(4))
	
         allocate( x(iplas) )
         allocate( x05(iplas-1) )
	
	   do i=1,iplas
	    !x(i)=(i-1.d0)/(iplas-1.d0)
	    x(i)=rho(i)/rho(iplas)
         enddo
         
	   do i=1,iplas-1
	    !x05(i)=(i-0.5d0)/(iplas-1.d0)
	    x05(i)=0.5d0*(x(i)+x(i+1))
         enddo
         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          x0=x(1)
          x1=x05(1)
          x2=x05(2)
          x3=x05(3)
          
          u1=arr_c(1)
          u2=arr_c(2)
          u3=arr_c(3)
 
	call EXTRP2(x0,u0, X1,X2,X3, u1,u2,u3)
          
          arr_n(1)=u0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          x0=x(iplas)
          x1=x05(iplas-1)
          x2=x05(iplas-2)
          x3=x05(iplas-3)
          
          u1=arr_c(iplas-1)
          u2=arr_c(iplas-2)
          u3=arr_c(iplas-3)
 
	call EXTRP2(x0,u0, X1,X2,X3, u1,u2,u3)
          
          arr_n(iplas)=u0

 
        n3spl=iplas-1

        CALL E01BAF(n3spl,x05,arr_c,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=2,iplas-1
               zspl=x(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               arr_n(i)=CWk(1)
           enddo
	
         deallocate(RRK,CCK,WRK)
         deallocate(CWK)
	
         deallocate( x )
         deallocate( x05 )
	
 	return
	end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	subroutine	get_umup(umt,upt)
	!implicit none
	include 'double.inc'
      include 'param.inc'
      include 'comblc.inc'
	real*8 umt, upt
      umt = um
      upt = up
	return
	end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	subroutine	f_get_umup(umt,upt)
	!implicit none
	include 'double.inc'
      include 'param.inc'
      include 'dimpl1.inc'
      include 'compol.inc'
	real*8 umt, upt
      umt = psim
      upt = psip
	return
	end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine vol2d_eq(equil_out,i1,i2,fcell2d,fvav)
      use vol2d_eq_face
      use imas_ids
	implicit none
		integer i1,i2
      type(type_equilibrium) equil_out
      real(8) :: fcell2d(i1,i2)
      real(8) :: fvav(i1+1)
          
      real*8, allocatable :: fvavs(:), rhowr(:), rhowrh(:)
      integer iplas, nt1
      real*8 r1,r2,r3,r4,r0, z1,z2,z3,z4, sss
      real*8 dvoli, fvavi
      real*8, allocatable :: RRK(:),CCK(:),WRK(:)
      real*8 CWK(4)
      real*8 zspl
      integer n3spl,IFAIL
      
      integer i,j,jerr
      integer ip0,ip1,ip2,ip3
      real*8 x0,x1,x2,x3
      real*8 funsq
      
      iplas = size(equil_out%coord_sys%position%r,1)
      nt1 = size(equil_out%coord_sys%position%r,2)

      allocate(fvavs(iplas+1), rhowr(iplas), rhowrh(iplas+1))
      
      do i=1,iplas-1
        dvoli = 0
        fvavi = 0
        do j=1,nt1-1
          r1=equil_out%coord_sys%position%r(i,j)
          r2=equil_out%coord_sys%position%r(i+1,j)
          r3=equil_out%coord_sys%position%r(i+1,j+1)
          r4=equil_out%coord_sys%position%r(i,j+1)
          r0=(r1+r2+r3+r4)*0.25d0
          z1=equil_out%coord_sys%position%z(i,j)
          z2=equil_out%coord_sys%position%z(i+1,j)
          z3=equil_out%coord_sys%position%z(i+1,j+1)
          z4=equil_out%coord_sys%position%z(i,j+1)
          sss=funsq(r1,r2,r3,r4,z1,z2,z3,z4)
          dvoli = dvoli + r0*sss
          fvavi = fvavi + fcell2d(i,j)*r0*sss
        enddo
        fvavs(i+1) = fvavi/dvoli
      enddo
      do i=1,iplas
        rhowr(i)=
     &   sqrt(equil_out%profiles_1d%phi(i)
     &   /equil_out%profiles_1d%phi(iplas))
      enddo
      do i=2,iplas
        rhowrh(i)=(rhowr(i)+rhowr(i-1))*0.5d0
      enddo
      rhowrh(iplas+1)=rhowr(iplas)
      rhowrh(1)=0.d0
!!!!!!!!!!!!!!!!extrapolation to the axis!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ip0=1
      ip1=2
      ip2=3
      ip3=4
      x3=rhowrh(ip3)
      x2=rhowrh(ip2)
      x1=rhowrh(ip1)
      x0=rhowrh(ip0)
      call	b_extrp(x0,x1,x2,x3, ip1,ip2,ip3,ip0,fvavs,jerr)
!!!!!!!!!!!!!!!!extrapolation to the boundary!!!!!!!!!!!!!!!!!!!!!!!!!       
      ip0=iplas+1
      ip1=iplas-2
      ip2=iplas-1
      ip3=iplas
      x1=rhowrh(ip3)
      x2=rhowrh(ip2)
      x3=rhowrh(ip1)
      x0=rhowrh(ip0)
      call	b_extrp(x0,x1,x2,x3, ip3,ip2,ip1,ip0,fvavs,jerr)  
!!!!!!!!!!!!!!!!!!!!!!!!!!cubic splines!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      n3spl=iplas+1
      allocate(RRK(n3spl+4),CCK(n3spl+4),WRK(6*n3spl+16))
      CALL E01BAF(n3spl,rhowrh,fvavs,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)
      do i=1,iplas
        zspl=rhowr(i)
        CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWK,IFAIL)
        fvav(i)=CWK(1)
      enddo
      fvav(iplas)=fvavs(iplas+1)
      
      return
      end
          









!!!!!!!!!!! new routine
	subroutine	put_eq(equil_out,parameters_spider)
	   use imas_ids
         use parameters_a2spider , only: type_parameters, TWOPI
       	 !implicit none
         include 'double.inc'
         include 'dim.inc'
         include 'compol.inc'
          common /com_volt/ upls(nrp)
         type(type_equilibrium) equil_out
       type(type_parameters) parameters_spider
	integer nfour_maxx,jindx
         parameter(nfour_maxx=5) 

         parameter(nrpl=nrp+1) 
         parameter(nrpl4=nrpl+4,nrpl6=nrpl4*6)
          
          real*8 RRK(nrpl4),CCK(nrpl4),WRK(nrpl6)
          real*8 CWK(4),dsqi,Zcurr,totcurr
          real*8 Rcurr,dumm1,dumm2
          
          real(8), allocatable :: rhowr(:), rhowrh(:),rhos(:)
          real*8 xxxx1(3),yyyy1(3),pppp1(3)

          real(8), allocatable :: g11(:), g22(:), g33(:), g2int(:)
          real(8), allocatable :: g41(:)
          real(8), allocatable :: gradro(:),droda(:),dpsidv(:),rbp_b2(:)
          real(8), allocatable :: bplfs(:),dvdpsi(:)
          real(8), allocatable :: surface_s(:)
          real(8), allocatable :: surface_x(:)

          real(8), allocatable :: acosB2a(:,:)
          real(8), allocatable :: asinB2a(:,:)
          real(8), allocatable :: acosBlnBa(:,:)
          real(8), allocatable :: asinBlnBa(:,:)
          
          real(8), allocatable :: g11s(:), g22s(:), g33s(:), g2ints(:)
          real(8), allocatable :: g41s(:)
          real(8), allocatable :: gradros(:),drodas(:),dpsidvs(:)
          real(8), allocatable :: bplfss(:),dvdpsis(:),sa(:)

          real(8), allocatable :: acosB2as(:,:)
          real(8), allocatable :: asinB2as(:,:)
          real(8), allocatable :: acosBlnBas(:,:)
          real(8), allocatable :: asinBlnBas(:,:)
          
          real(8), allocatable :: volum(:),dvdro(:)
										real(8) :: theta_poloidal
          real(8), allocatable :: arr2(:,:)
          
          
           real*8, allocatable :: b_maxt(:),b_mint(:),
     *                            b_db02(:),b_db0(:),b_0db2(:)
         
           real*8, allocatable :: bmaxt(:),bmint(:),
     *                            bdb02(:),bdb0(:),b0db2(:)
          
           real*8, allocatable :: Btot(:,:),fcell2d1(:,:),fcell2d2(:,:)
           real*8, allocatable :: fcell2d21(:,:,:),fcell2d22(:,:,:)
           real*8, allocatable :: fcell2d23(:,:,:),fcell2d24(:,:,:)
           real*8, allocatable :: traps(:),dvoliz(:)
           real*8, allocatable :: yFOFB(:)
           real*8, allocatable :: areats(:)
           real*8, allocatable :: areat(:)
           real*8, allocatable :: perims(:)
           real*8, allocatable :: perim(:)
	real*8 :: dum1           
           !geometry
           real*8, allocatable :: ya(:),yra(:),yri(:),yshif(:),yshiv(:),
     *                        yelon(:),ytria_u(:),ytria_l(:)
          
										real*8 psplexs,dpsplexs,qedge
										real*8 dum_area,dum_perim
										real*8 ggreen,ddllt,psplex,pressvoli

           !2d in cells
            allocate(equil_out%coord_sys%gradvcell(iplas-1,nt1-1))
            allocate(equil_out%coord_sys%bpcell(iplas-1,nt1-1))
            allocate(equil_out%coord_sys%bcell(iplas-1,nt1-1))
            allocate(equil_out%coord_sys%rcell(iplas-1,nt1-1))
! for new computations
            allocate(fcell2d1(iplas-1,nt1-1))
            allocate(fcell2d2(iplas-1,nt1-1))
            allocate(fcell2d21(iplas-1,nt1-1,nfour_maxx))
            allocate(fcell2d22(iplas-1,nt1-1,nfour_maxx))
            allocate(fcell2d23(iplas-1,nt1-1,nfour_maxx))
            allocate(fcell2d24(iplas-1,nt1-1,nfour_maxx))

          
          allocate( g11(iplas), g22(iplas), g33(iplas) )
          allocate( g41(iplas) )
          allocate( rbp_b2(iplas), bplfs(iplas))
          allocate( g2int(iplas) , dpsidv(iplas))
          allocate( gradro(iplas), droda(iplas) )
          allocate( volum(iplas) )
          allocate( dvdpsi(iplas) )
          allocate( dvdro(iplas) )
          allocate( surface_x(iplas) )
          allocate( surface_s(iplas+1) )
          allocate( areat(iplas) )
          allocate( areats(iplas+1) )
          allocate( perim(iplas) )
          allocate( perims(iplas+1) )

          allocate( acosB2a(iplas,nfour_maxx) )
          allocate( asinB2a(iplas,nfour_maxx) )
          allocate( acosBlnBa(iplas,nfour_maxx) )
          allocate( asinBlnBa(iplas,nfour_maxx) )

          allocate( acosB2as(iplas+1,nfour_maxx) )
          allocate( asinB2as(iplas+1,nfour_maxx) )
          allocate( acosBlnBas(iplas+1,nfour_maxx) )
          allocate( asinBlnBas(iplas+1,nfour_maxx) )

          allocate( g11s(iplas+1), g22s(iplas+1), g33s(iplas+1) )
          allocate( g41s(iplas+1) )
          allocate( g2ints(iplas+1),dpsidvs(iplas+1) )
          allocate( gradros(iplas+1), drodas(iplas+1) )
          allocate( bplfss(iplas+1))
          allocate( dvoliz(iplas+1))
          allocate( dvdpsis(iplas+1) )
          allocate( sa(iplas+1) )

          allocate( arr2(iplas,nt) )
          
          allocate( rhowr(iplas), rhowrh(iplas+1) )
          allocate( rhos(iplas) )
          allocate ( Btot(iplas,nt) )

          allocate (b_maxt(iplas+1),b_mint(iplas+1),
     *               b_db02(iplas+1),b_db0(iplas+1),b_0db2(iplas+1))
         
          allocate (bmaxt(iplas),bmint(iplas),
     *               bdb02(iplas),bdb0(iplas),b0db2(iplas))
           
          allocate (traps(iplas+1))
          allocate (yFOFB(iplas))
          
          btor=b0ax
          rtor=r0ax

!
          

           volum(1)=0.d0
           pressvoli=0.d0
         do i=2,iplas
           voli=0.d0
          do j=2,nt1
           voli=voli+vol(i-1,j)
          enddo
           volum(i)=voli+volum(i-1)
	pressvoli=pressvoli+
     &   equil_out%profiles_1d%pressure(i-1)*voli
         enddo
!	write(*,*) volum(1:iplas)

!	pressvoli=pressvoli/volum(iplas)

         do i=1,iplas
          rhowr(i)=sqrt(flx_fi(i)/flx_fi(iplas))
										rhos(i)=sqrt(flx_fi(i)/(pi*btor))
         enddo
!         do i=1,iplas
!          rhowr(i)=sqrt(volum(i)/volum(iplas))
!										rhos(i)=sqrt(flx_fi(i)/(pi*btor))
!         enddo

!note that SPIDER natural grid is equispaced rhos=sqrt(Phi/(pi*btor)) !

         do i=2,iplas
          rhowrh(i)=(rhowr(i)+rhowr(i-1))*0.5d0
         enddo
          rhowrh(iplas+1)=rhowr(iplas)
          rhowrh(1)=0.d0


C		A=GREENI(R,Z,RP,ZP)
	psplexs=0.
	psplex=0.
	dpsplexs=0.
	Zcurr=0.
	Rcurr=0.
	totcurr=0.
         do i=1,iplas-1

           dsqi=0.d0

          dvoli=volum(i+1)-volum(i)
          u1=volum(i)
          u2=volum(i+1)
          u3=volum(i+1)
          u4=volum(i)
												dvdro(i+1)=2.d0*pi*dvoli/(rhos(i+1)-rhos(i))

          do j=2,nt1
          
          cur1=cur(i,j)
          cur2=cur(i+1,j)
          cur3=cur(i+1,j+1)
          cur4=cur(i,j+1)
          cur0=(cur1+cur2+cur3+cur4)*0.25d0
          
	!Zcurr=Zcurr+z(i,j)*cur(i,j)*vol(i,j)/r(i,j)
	!Rcurr=Rcurr+r(i,j)*cur(i,j)*vol(i,j)/r(i,j)
	!totcurr=totcurr+cur(i,j)*vol(i,j)/r(i,j)

          r1=r(i,j)
          r2=r(i+1,j)
          r3=r(i+1,j+1)
          r4=r(i,j+1)
          r0=(r1+r2+r3+r4)*0.25d0

          z1=z(i,j)
          z2=z(i+1,j)
          z3=z(i+1,j+1)
          z4=z(i,j+1)
          z0=(z1+z2+z3+z4)*0.25d0
          
	Zcurr=Zcurr+z0*cur0*s(i,j)
	Rcurr=Rcurr+r0*cur0*s(i,j)
	totcurr=totcurr+cur0*s(i,j)

           sss=s(i,j)
           
           dvdr= 0.5d0*((u1+u2)*(z2-z1)+(u2+u3)*(z3-z2)+
     *                   (u3+u4)*(z4-z3)+(u4+u1)*(z1-z4))/sss

           dvdz=-0.5d0*((u1+u2)*(r2-r1)+(u2+u3)*(r3-r2)+
     *                   (u3+u4)*(r4-r3)+(u4+u1)*(r1-r4))/sss

           gr2=dvdr**2+dvdz**2

            arr2(i,j)=gr2
            
            equil_out%coord_sys%gradvcell(i,j-1) = sqrt(gr2)*2*PI

          enddo
         enddo



	Zcurr=Zcurr/totcurr
	Rcurr=Rcurr/totcurr

        call avr2_c(arr2,iplas,nt,g11s(2))
	areats(1)=0.
	perims(1)=0.
	g2ints=0.
	
         do i=1,iplas-1

           dsqi=0.d0

C          dvoli=volum(i+1)-volum(i)
          dvoli=0.
          dvoliz(i)=dvoli
                     
           av_r2=0.d0
           avr2=0.d0
           av_gr1=0.d0
           av_grr2=0.d0
	dum_area=0.
	dum_perim=0.
          do j=2,nt1

           dsqi=dsqi+sr(i+1,j)
           dvoli=dvoli+vol(i,j)


          r1=r(i,j)
          r2=r(i+1,j)
          r3=r(i+1,j+1)
          r4=r(i,j+1)
          r0=(r1+r2+r3+r4)*0.25d0
          
            gr2=arr2(i,j)
            
            av_gr1=av_gr1+vol(i,j)*dsqrt(gr2)
            av_grr2=av_grr2+vol(i,j)*gr2/r0**2
            av_r2=av_r2+vol(i,j)/r0**2
            avr2=avr2+r0**2.*vol(i,j)
            
            equil_out%coord_sys%rcell(i,j-1) = r0
	!dum_area=dum_area+vol(i,j)/r0         
	dum_area=dum_area+s(i,j)         
	!dum_perim=dum_perim+sr(i+1,j)/r0         
	dum_perim=dum_perim+dlt(i+1,j)         
          enddo

          sa(i+1)=dsqi*2.d0*pi
          g22s(i+1)=av_grr2/dvoli
          g2ints(i+1)=g2ints(i)+av_grr2
          g33s(i+1)=av_r2/dvoli
          g41s(i+1)=avr2/dvoli
          gradros(i+1)=av_gr1/dvoli
          drodas(i+1)=dvoli/(rhowr(i+1)-rhowr(i))
										surface_s(i+1)=sa(i+1)
										areats(i+1)=areats(i)+dum_area
										perims(i+1)=dum_perim
         enddo


cccc      magnetic field

           bp2av=0.d0
		   bp2=0.d0
cc          do i=2,iplas-1 
          sa(1)=0.
          do i=1,iplas-1 

           dsqi=0.d0
           dvoli=0.d0
           bmx=0.d0
           bmn=1.d9*btor
           bavr=0.d0
           bavr2=0.d0
           bavrd2=0.d0
           theta_poloidal=0.d0

           do j=2,nt1

           dum1=dvoli
           dvoli=dvoli+vol(i,j)
C	write(*,*) theta_poloidal
          r1=r(i,j)
          r2=r(i+1,j)
          r3=r(i+1,j+1)
          r4=r(i,j+1)
          r0=(r1+r2+r3+r4)*0.25d0

            bm_pol=psim*(psia(i+1)-psia(i))/st(i,j)
            bp_pol=psim*(psia(i+1)-psia(i))/st(i,j+1)

         if(i.ne.1) then
        bpol2v=bm_pol**2*( vol1(i,j)/sin1(i,j)+vol2(i,j)/sin2(i,j))+ 
     +         bp_pol**2*( vol3(i,j)/sin3(i,j)+vol4(i,j)/sin4(i,j)) 
         else
        bpol2v=bm_pol**2*(                     vol2(i,j)/sin2(i,j))+ 
     +         bp_pol**2*( vol3(i,j)/sin3(i,j)                    ) 
         endif
         
            bpol2=bpol2v /vol(i,j)
!	write(*,*) 'bpol',i,j,bpol2
						b_fi2=(f(i)/r0)**2
         
            equil_out%coord_sys%bpcell(i,j-1) = sqrt(bpol2)
           
            b_tot2=bpol2+b_fi2
            b_tot=dsqrt(b_tot2)
            
            equil_out%coord_sys%bcell(i,j-1) = b_tot
            
            fcell2d1(i,j-1) = bpol2/b_tot2*r0**2.0
C	write(*,*) b_tot2
           dsqi=dsqi+dlt(i+1,j)/sqrt(bpol2)

!!   aai 06/02/10{{{{{{{{{{{
           Btot(i,j)=b_tot
!!   }}}}}}}}}}}

            bmx=dmax1(b_tot,bmx)
            bmn=dmin1(b_tot,bmn)

           bavr2=bavr2+b_tot2*vol(i,j)
           bavr =bavr +b_tot *vol(i,j)
	   bavrd2=bavrd2+vol(i,j)/(b_tot2+1.d-8)
           bp2av=bp2av+bpol2v
		bp2=bp2+bpol2v
          enddo
										dvoliz(i)=dvoli         
	dvdpsis(i+1)=dsqi
	dpsidvs(i+1)=1./dvdpsis(i+1)
           b_maxt(i+1)=bmx
           b_mint(i+1)=bmn
           b_db02(i+1)=bavr2/(dvoli*btor**2)
           b_0db2(i+1)=(bavrd2*btor**2)/(dvoli)
           b_db0(i+1)=bavr/(dvoli*btor)

        enddo 
C	write(*,*) 'li3 = ', tokp,bp2,rtor,4.d0*pi*bp2/(rtor*(.4d0*pi*tokp)**2)    
C	write(*,*) 'betpol = ',bp2,pressvoli,
C     &   (pressvoli/bp2)*2.*(2.*TWOPI*1.E-7) 

	equil_out%global_param%li3	= 4.d0*pi*bp2/
     &  (rtor*(.4d0*pi*tokp)**2)
	equil_out%global_param%betpol	= 
     & (pressvoli/bp2)*2.*(2.*TWOPI*1.E-7)
	equil_out%global_param%wkin	= TWOPI*pressvoli
	equil_out%global_param%bpkin	= TWOPI*bp2/(2.*(2.*TWOPI*1.E-7))

        do i=1,iplas-1
           Bmax=0.d0
           Svol=0.d0         
          do j=2,nt-1
           B_ij=Btot(i,j)
           Bmax=dmax1(Bmax,B_ij)
           Svol=Svol+vol(i,j)
          enddo

           avr_v=0.d0
          do j=2,nt-1
           Bnor=Btot(i,j)/Bmax
           avr_v=avr_v+
     &          (b0ax/Btot(i,j))**2
     &         *( 1.d0-dsqrt(1.d0-Bnor)*(1.d0+.5d0*Bnor) )*vol(i,j)
          enddo
          
           !trap(i)=avr_v/Svol
           traps(i+1)=avr_v/Svol
           
        enddo


! stuff for psi as b.c.
CEfable
	!if (parameters_spider%k_fixfree.eq.1) then

      qedge = equil_out%profiles_1d%q(iplas)
	if (parameters_spider%k_grid.eq.0) then
!	call psib_pla_V(psplexs,btor)   ! my g_integral
!	psplexs = psplexs ! my G integral
	call psib_pla(psplexs)    ! spider G integral
	psplexs = psplexs/(btor*rhos(iplas)**2)*qedge  ! spider G integral
	endif

	if (parameters_spider%k_grid.eq.1) then
!	call psib_pla_V(psplexs,btor)   ! my g_integral
!	psplexs = psplexs ! my G integral
	!call f_psib_pla(psplexs)    ! spider G integral
	call psib_pla(psplexs)    ! spider G integral
	psplexs = psplexs/(btor*rhos(iplas)**2)*qedge  ! spider G integral
	endif

	!endif







cccc      cosinus sinus fourier modes of theta_gen

		ip0 = max(2,minloc(abs(teta(1:nt1)),1))
          do i=1,iplas-1 

           dsqi=0.d0
           dvoli=0.d0
           theta_poloidal=0.d0

           do j=ip0,nt1
C           dvoli=dvoli+(vol(i,j)+vol(i,j-1))/2.

           dvoli=dvoli+equil_out%coord_sys%bcell(i,j-1)*vol(i,j)

	theta_poloidal= TWOPI*dvoli/dvoliz(i)*1./(btor*b_db0(i+1))

C	write(*,*) theta_poloidal

	do jindx=1,nfour_maxx
	fcell2d21(i,j-1,jindx) = equil_out%coord_sys%bcell(i,j-1)**2.*
     &  cos(jindx*theta_poloidal)
	fcell2d22(i,j-1,jindx) = equil_out%coord_sys%bcell(i,j-1)*
     &  log(equil_out%coord_sys%bcell(i,j-1))*cos(jindx*theta_poloidal)
	fcell2d23(i,j-1,jindx) = equil_out%coord_sys%bcell(i,j-1)**2.*
     &  sin(jindx*theta_poloidal)
	fcell2d24(i,j-1,jindx) = equil_out%coord_sys%bcell(i,j-1)*
     &  log(equil_out%coord_sys%bcell(i,j-1))*sin(jindx*theta_poloidal)
	enddo
C	write(*,*) equil_out%coord_sys%bcell(i,j-1)**2.*cos(1*theta_poloidal)
	enddo
           do j=2,ip0-1
C           dvoli=dvoli+(vol(i,j)+vol(i,j-1))/2.
           dvoli=dvoli+equil_out%coord_sys%bcell(i,j-1)*vol(i,j)

	theta_poloidal= TWOPI*dvoli/dvoliz(i)*1./(btor*b_db0(i+1))

C	write(*,*) theta_poloidal

	do jindx=1,nfour_maxx
	fcell2d21(i,j-1,jindx) = equil_out%coord_sys%bcell(i,j-1)**2.*
     &  cos(jindx*theta_poloidal)
	fcell2d22(i,j-1,jindx) = equil_out%coord_sys%bcell(i,j-1)*
     &  log(equil_out%coord_sys%bcell(i,j-1))*cos(jindx*theta_poloidal)
	fcell2d23(i,j-1,jindx) = equil_out%coord_sys%bcell(i,j-1)**2.*
     &  sin(jindx*theta_poloidal)
	fcell2d24(i,j-1,jindx) = equil_out%coord_sys%bcell(i,j-1)*
     &  log(equil_out%coord_sys%bcell(i,j-1))*sin(jindx*theta_poloidal)
	enddo

C	write(*,*) equil_out%coord_sys%bcell(i,j-1)**2.*cos(1*theta_poloidal)
		enddo            
	enddo 


! COmpute <R^2 Bp^2 / B^2>
	call vol2d_eq(equil_out,iplas-1,nt1-1,fcell2d1,rbp_b2)
! dpsidv
!	call vol2d_eq(equil_out,iplas-1,nt1-1,fcell2d2,dpsidv)

	do jindx=1,nfour_maxx
	call vol2d_eq(equil_out,iplas-1,nt1-1,
     &     fcell2d21(:,:,jindx),acosB2a(:,jindx))
	call vol2d_eq(equil_out,iplas-1,nt1-1,
     &     fcell2d22(:,:,jindx),acosBlnBa(:,jindx))
	call vol2d_eq(equil_out,iplas-1,nt1-1,
     &     fcell2d23(:,:,jindx),asinB2a(:,jindx))
	call vol2d_eq(equil_out,iplas-1,nt1-1,
     &     fcell2d24(:,:,jindx),asinBlnBa(:,jindx))
	enddo

!Compute bp at low field side
C	write(*,*) teta(1:nt1)
	dumm1=100.
	do jindx=2,nt1
	dumm2=abs(z(iplas,jindx)-z(1,2))
	if (dumm2.lt.dumm1 
     & .and. r(iplas,jindx).gt.r(1,2)) then
		j=jindx
		dumm1=dumm2
	endif
	enddo

	j = minloc(abs(teta(2:nt1)),1)
	do i=1,iplas-1 
		bplfss(i+1)= equil_out%coord_sys%bpcell(i,j)
	enddo
!!!!!!!!!!!!!!!!extrapolation to the axis!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	  ip0=1
	  ip1=2
	  ip2=3
	  ip3=4

          x3=rhowrh(ip3)
          x2=rhowrh(ip2)
          x1=rhowrh(ip1)
          x0=rhowrh(ip0)
          
C 	  call	b_extrp(x0,X1,X2,X3, ip1,ip2,ip3,ip0,g11s,jerr)
         g11s(1)=0.d0
C 	  call	b_extrp(x0,X1,X2,X3, ip1,ip2,ip3,ip0,g22s,jerr)
         g22s(1)=0.d0
C 	  call	b_extrp(x0,X1,X2,X3, ip1,ip2,ip3,ip0,g2ints,jerr)
         g2ints(1)=0.d0
C 	  call	b_extrp(x0,X1,X2,X3, ip1,ip2,ip3,ip0,g33s,jerr)
         g33s(1)=1.d0/rm**2
         g41s(1)=rm**2
 	  call	b_extrp(x0,X1,X2,X3, ip1,ip2,ip3,ip0,dpsidvs,jerr)
 	  call	b_extrp(x0,X1,X2,X3, ip1,ip2,ip3,ip0,gradros,jerr)
C 	  call	b_extrp(x0,X1,X2,X3, ip1,ip2,ip3,ip0,drodas,jerr)
         drodas(1)=0.d0
         areats(1)=0.d0
         perims(1)=0.d0
C 	  call	b_extrp(x0,X1,X2,X3, ip1,ip2,ip3,ip0,bplfss,jerr)
         bplfss(1)=0.d0
         surface_s(1)=0.d0
!!!!!!!!!!!!b_db02(1)
	call	b_extrp(x0, X1,X2,X3, ip1,ip2,ip3,ip0,b_db02,jerr)

!!!!!!!!!!!!b_0db2(1)
	call	b_extrp(x0, X1,X2,X3, ip1,ip2,ip3,ip0,b_0db2,jerr)
 
!!!!!!!!!!!!b_db0(1)
	call	b_extrp(x0, X1,X2,X3, ip1,ip2,ip3,ip0,b_db0,jerr)
!!!!!!!!!!!!tpaps(1)
	call	b_extrp(x0,X1,X2,X3, ip1,ip2,ip3,ip0,traps,jerr)

	b_maxt(1)=b_db0(1)*btor
	b_mint(1)=b_maxt(1)




C	acosB2as(1,:)=0.
C	asinB2as(1,:)=0.
C	acosBlnBas(1,:)=0.
C	asinBlnBas(1,:)=0.
!!!!!!!!!!!!!!!!extrapolation to the boundary!!!!!!!!!!!!!!!!!!!!!!!!!
          


  	ip0=iplas+1
	ip1=iplas-2
	ip2=iplas-1
	ip3=iplas

          x1=rhowrh(ip3)
          x2=rhowrh(ip2)
          x3=rhowrh(ip1)
          x0=rhowrh(ip0)


C  			  call EXTRAP_EF(rhowrh(1:iplas-1),
C     &      g22s(1:iplas-1),x1,
C     &      iplas-1,g22s(iplas),1,iplas-1)			

  			  call EXTRAP_EFSP(rhowrh(1:iplas),
     &      g22s(1:iplas),x0,
     &      iplas,g22s(ip0),1,iplas)			

  			  call EXTRAP_EFSP(rhowrh(1:iplas),
     &      dpsidvs(1:iplas),x0,
     &      iplas,dpsidvs(ip0),2,iplas)			
C	do i=2,iplas+1
C	if (g22s(i).lt.g22s(i-1)) g22s(i)=g22s(i-1)
C	enddo

C  			  call EXTRAP_EFSP(rhowrh(iplas-2:iplas),
C     &      psplexs(iplas-2:iplas),x0,
C     &      3,psplex(iplas),1,3)			
	psplex=psplexs
C	write(777,*)  psplexs,psplex

C  			  call EXTRAP_EF(rhowrh(1:iplas-1),
C     &      g33s(1:iplas-1),x1,
C     &      iplas-1,g33s(iplas),1,iplas-1)			

C  			  call EXTRAP_EF(rhowrh(1:iplas),
C     &      g33s(1:iplas),x0,
C     &      iplas,g33s(ip0),1,iplas)			


 	  call	b_extrp(x0,X1,X2,X3, ip3,ip2,ip1,ip0,g11s,jerr)
C 	  call	b_extrp(x0,X1,X2,X3, ip3,ip2,ip1,ip0,g22s,jerr)
 	  call	b_extrp(x0,X1,X2,X3, ip3,ip2,ip1,ip0,g2ints,jerr)
 	  call	b_extrp(x0,X1,X2,X3, ip3,ip2,ip1,ip0,g33s,jerr)
 	  call	b_extrp(x0,X1,X2,X3, ip3,ip2,ip1,ip0,g41s,jerr)
C 	  call	b_extrp(x0,X1,X2,X3, ip3,ip2,ip1,ip0,dpsidvs,jerr)
 	  call	b_extrp(x0,X1,X2,X3, ip3,ip2,ip1,ip0,gradros,jerr)
 	  call	b_extrp(x0,X1,X2,X3, ip3,ip2,ip1,ip0,drodas,jerr)
 	  call	b_extrp(x0,X1,X2,X3, ip3,ip2,ip1,ip0,areats,jerr)
 	  call	b_extrp(x0,X1,X2,X3, ip3,ip2,ip1,ip0,perims,jerr)
 	  call	b_extrp(x0,X1,X2,X3, ip3,ip2,ip1,ip0,bplfss,jerr)
 	  call	b_extrp(x0,X1,X2,X3, ip3,ip2,ip1,ip0,surface_s,jerr)
       
!!!!!!!!!!!!b_db02(iplas+1)
	call	b_extrp(x0, X1,X2,X3, ip3,ip2,ip1,ip0,b_db02,jerr)
!!!!!!!!!!!!b_0db2(iplas+1)
	call	b_extrp(x0, X1,X2,X3, ip3,ip2,ip1,ip0,b_0db2,jerr)
!!!!!!!!!!!!b_db0(iplas+1)
	call	b_extrp(x0, X1,X2,X3, ip3,ip2,ip1,ip0,b_db0,jerr)
!!!!!!!!!!!!tpaps(iplas+1)
	call	b_extrp(x0,X1,X2,X3, ip3,ip2,ip1,ip0,traps,jerr)
!!!!!!!!!!!!b_maxt(1)
	call	b_extrp(x0, X1,X2,X3, ip3,ip2,ip1,ip0,b_maxt,jerr)
!!!!!!!!!!!!b_mint(1)
	call	b_extrp(x0, X1,X2,X3, ip3,ip2,ip1,ip0,b_mint,jerr)

C	do jindx=1,nfour_maxx
C	call	b_extrp(x0, X1,X2,X3, ip3,ip2,ip1,ip0,acosB2as(:,jindx),jerr)
C	call	b_extrp(x0, X1,X2,X3, ip3,ip2,ip1,ip0,asinB2as(:,jindx),jerr)
C	call	b_extrp(x0, X1,X2,X3, ip3,ip2,ip1,ip0,acosBlnBas(:,jindx),jerr)
C	call	b_extrp(x0, X1,X2,X3, ip3,ip2,ip1,ip0,asinBlnBas(:,jindx),jerr)
C	enddo


!!!!!!!!!!!!!!!!extrapolation to the boundary!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          
        n3spl=iplas+1

!!!!!!!!!!!!!!!!!!!!!!!!!!g11!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CALL E01BAF(n3spl,rhowrh,g11s,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,iplas
               zspl=rhowr(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               g11(i)=CWk(1)
           enddo
               g11(1)=g11s(1)
               g11(iplas)=g11s(iplas+1)
!!!!!!!!!!!!!!!!!!!!!!!!!!g22!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CALL E01BAF(n3spl,rhowrh,g22s,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,iplas
               zspl=rhowr(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               g22(i)=CWk(1)
           enddo
               g22(1)=g22s(1)
               g22(iplas)=g22s(iplas+1)


!!!!!!!!!!!!!!!!!!!!!!!!!!g33!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CALL E01BAF(n3spl,rhowrh,g33s,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,iplas
               zspl=rhowr(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               g33(i)=CWk(1)
           enddo
               g33(1)=g33s(1)
               g33(iplas)=g33s(iplas+1)
!!!!!!!!!!!!!!!!!!!!!!!!!!g41!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CALL E01BAF(n3spl,rhowrh,g41s,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,iplas
               zspl=rhowr(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               g41(i)=CWk(1)
           enddo
               g41(1)=g41s(1)
               g41(iplas)=g41s(iplas+1)

!!!!!!!!!!!!!!!!!!!!!!!!!!dpsidv!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CALL E01BAF(n3spl,rhowrh,dpsidvs,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,iplas
               zspl=rhowr(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               dpsidv(i)=CWk(1)
           enddo
               dpsidv(1)=dpsidvs(1)
               dpsidv(iplas)=dpsidvs(iplas+1)



!!!!!!!!!!!!!!!!!!!!!!!!!!g33!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C	do jindx=1,nfour_maxx
C        CALL E01BAF(n3spl,rhowrh,acosB2as(:,jindx),RRK,CCK,
C     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

C          do i=1,iplas
C               zspl=rhowr(i)
C              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
C               acosB2a(i,jindx)=CWk(1)
C           enddo
C               acosB2a(1,jindx)=acosB2as(1,jindx)
C               acosB2a(iplas,jindx)=acosB2as(iplas+1,jindx)


C        CALL E01BAF(n3spl,rhowrh,asinB2as(:,jindx),RRK,CCK,
C     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)
C          do i=1,iplas
C               zspl=rhowr(i)
C             CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
C              asinB2a(i,jindx)=CWk(1)
C           enddo
C               asinB2a(1,jindx)=asinB2as(1,jindx)
C               asinB2a(iplas,jindx)=asinB2as(iplas+1,jindx)


C        CALL E01BAF(n3spl,rhowrh,acosBlnBas(:,jindx),RRK,CCK,
C     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)
C          do i=1,iplas
C               zspl=rhowr(i)
C              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
C               acosBlnBa(i,jindx)=CWk(1)
C           enddo
C               acosBlnBa(1,jindx)=acosBlnBas(1,jindx)
C               acosBlnBa(iplas,jindx)=acosBlnBas(iplas+1,jindx)



C        CALL E01BAF(n3spl,rhowrh,asinBlnBas(:,jindx),RRK,CCK,
C     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)
C          do i=1,iplas
C               zspl=rhowr(i)
C              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
C               asinBlnBa(i,jindx)=CWk(1)
C           enddo
C               asinBlnBa(1,jindx)=asinBlnBas(1,jindx)
C               asinBlnBa(iplas,jindx)=asinBlnBas(iplas+1,jindx)
C	enddo





!!!!!!!!!!!!!!!!!!!!!!!!!!g22!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              
C	call qinterp_ef(rhowrh(1:iplas+1),g22s(1:iplas+1),
C     *  iplas+1,rhowr(1:iplas),g22(1:iplas),iplas)														
C														 g22(1)=g22s(1)
C               g22(iplas)=g22s(iplas+1)


C	call qinterp_ef(rhowrh(1:iplas+1),g33s(1:iplas+1),
C     *  iplas+1,rhowr(1:iplas),g33(1:iplas),iplas)														
C														 g33(1)=g33s(1)
C               g33(iplas)=g33s(iplas+1)



!!!!!!!!!!!!!!!!!!!!!!!!!!g2int!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CALL E01BAF(n3spl,rhowrh,g2ints,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,iplas
               zspl=rhowr(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               g2int(i)=CWk(1)
           enddo
               g2int(1)=g2ints(1)
               g2int(iplas)=g2ints(iplas+1)
!!!!!!!!!!!!!!!!!!!!!!!!!!gradro!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CALL E01BAF(n3spl,rhowrh,gradros,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,iplas
               zspl=rhowr(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               gradro(i)=CWk(1)
           enddo
               gradro(1)=gradros(1)
               gradro(iplas)=gradros(iplas+1)
!!!!!!!!!!!!!!!!!!!!!!!!!!droda!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CALL E01BAF(n3spl,rhowrh,drodas,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,iplas
               zspl=rhowr(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               droda(i)=CWk(1)
           enddo
               droda(1)=drodas(1)
               droda(iplas)=drodas(iplas+1)
               
!!!!!!!!!!!!!!!!!!!!!!!!!!bplfs!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CALL E01BAF(n3spl,rhowrh,bplfss,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,iplas
               zspl=rhowr(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               bplfs(i)=CWk(1)
           enddo
               bplfs(1)=bplfss(1)
               bplfs(iplas)=bplfss(iplas+1)
               
!!!!!!!!!!!!!!!!!!!!! bdb02 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CALL E01BAF(n3spl,rhowrh,b_db02,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,iplas
               !zspl=(rhowr(i)+rhowr(i+1))*0.5d0
               zspl=rhowr(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               bdb02(i)=CWk(1)
           enddo
               bdb02(1)=b_db02(1)
               bdb02(iplas)=b_db02(iplas)

!!!!!!!!!!!!!!!!!!!!! b0db2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CALL E01BAF(n3spl,rhowrh,b_0db2,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,iplas
               !zspl=(rhowr(i)+rhowr(i+1))*0.5d0
               zspl=rhowr(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               b0db2(i)=CWk(1)
           enddo
               b0db2(1)=b_0db2(1)
               b0db2(iplas)=b_0db2(iplas)

cw	write(*,*) 'b0db2'

!!!!!!!!!!!!!!!!!!!!! bdb0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CALL E01BAF(n3spl,rhowrh,b_db0,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,iplas
               !zspl=(rhowr(i)+rhowr(i+1))*0.5d0
               zspl=rhowr(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               bdb0(i)=CWk(1)
          enddo
              bdb0(1)=b_db0(1)
              bdb0(iplas)=b_db0(iplas)
!!!!!!!!!!!!!!!!!!!!!!! yFOFB !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CALL E01BAF(n3spl,rhowrh,traps,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,iplas
               !zspl=(rhowr(i)+rhowr(i+1))*0.5d0
               zspl=rhowrh(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               yFOFB(i)=CWk(1)
           enddo
               yFOFB(1)=traps(1)
               yFOFB(iplas)=traps(iplas+1)
!!!!!!!!!!!!!!!!!!!!!!!!!!surface!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CALL E01BAF(n3spl,rhowrh,surface_s,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,iplas
               zspl=rhowr(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               surface_x(i)=CWk(1)
           enddo
               surface_x(1)=surface_s(1)
               surface_x(iplas)=surface_s(iplas+1)
!!!!!!!!!!!!!!!!!!!!!!!!!!area!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CALL E01BAF(n3spl,rhowrh,areats,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,iplas
               zspl=rhowr(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               areat(i)=CWk(1)
           enddo
               areat(1)=areats(1)
               areat(iplas)=areats(iplas+1)
!!!!!!!!!!!!!!!!!!!!!!!!!!perim!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CALL E01BAF(n3spl,rhowrh,perims,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,iplas
               zspl=rhowr(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               perim(i)=CWk(1)
           enddo
               perim(1)=perims(1)
               perim(iplas)=perims(iplas+1)

!!!!!!!!!!!!!!!!!!!!!!!!!!bmint!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CALL E01BAF(n3spl,rhowrh,b_mint,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,iplas
               !zspl=(rhowr(i)+rhowr(i+1))*0.5d0
               zspl=rhowrh(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               bmint(i)=CWk(1)
           enddo
               bmint(iplas)=b_mint(iplas)

cw	write(*,*) 'bmint'

!!!!!!!!!!!!!!!!!!!!! bmaxt !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CALL E01BAF(n3spl,rhowrh,b_maxt,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,iplas
               !zspl=(rhowr(i)+rhowr(i+1))*0.5d0
               zspl=rhowrh(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               bmaxt(i)=CWk(1)
           enddo
               bmaxt(iplas)=b_maxt(iplas)

       if(.NOT. associated(equil_out%profiles_1d%g1)) then
           allocate(equil_out%profiles_1d%gm1(iplas))
           allocate(equil_out%profiles_1d%gm4(iplas))
           allocate(equil_out%profiles_1d%gm5(iplas))
           allocate(equil_out%profiles_1d%gm41(iplas))

           allocate(equil_out%profiles_1d%rbp_b2(iplas))
           allocate(equil_out%profiles_1d%bplfs(iplas))

           allocate(equil_out%profiles_1d%acosB2a(iplas,nfour_maxx))
           allocate(equil_out%profiles_1d%asinB2a(iplas,nfour_maxx))
           allocate(equil_out%profiles_1d%acosBlnBa(iplas,nfour_maxx))
           allocate(equil_out%profiles_1d%asinBlnBa(iplas,nfour_maxx))
           
           allocate(equil_out%profiles_1d%g1(iplas))
           allocate(equil_out%profiles_1d%g2(iplas))
           allocate(equil_out%profiles_1d%g2int(iplas))
           allocate(equil_out%profiles_1d%fofb(iplas))
           allocate(equil_out%profiles_1d%areat(iplas))
           allocate(equil_out%profiles_1d%perim(iplas))
           allocate(equil_out%profiles_1d%ggradro(iplas))
           allocate(equil_out%profiles_1d%gdroda(iplas))
           allocate(equil_out%profiles_1d%bmaxt(iplas))
           allocate(equil_out%profiles_1d%bmint(iplas))
           allocate(equil_out%profiles_1d%bdb0(iplas))
           allocate(equil_out%profiles_1d%dPSIdV(iplas))
           allocate(equil_out%profiles_1d%surface(iplas))
       endif

! compute dpsidv

!        equil_out%profiles_1d%psi(i) = 
!     &                       ((psim-psip)*psia(i)+psi_boundary)*TWOPI

C	i=1
C       equil_out%profiles_1d%dPSIdV(i) =-(psim-psip)* 
C     *  (psia(2)-psia(1))/(volum(2)-volum(1))
C	do i=2,iplas-1
C       equil_out%profiles_1d%dPSIdV(i) = -(psim-psip)*
C     *  (psia(i+1)-psia(i-1))/(volum(i+1)-volum(i-1))
C	enddo
C	i=iplas
C  			  call EXTRAP_EF(volum(1:iplas-1),
C     &  equil_out%profiles_1d%dPSIdV(1:iplas),
C     *  volum(iplas),
C     &  iplas-1,equil_out%profiles_1d%dPSIdV(iplas),
C     *  1,iplas-1)			

	equil_out%global_param%Zcurr	= Zcurr
	equil_out%global_param%Rcurr	= Rcurr
	equil_out%global_param%Vloop	= upls(iplas)
C	write(*,*) 'upls ',upls(iplas)
							equil_out%profiles_1d%dPSIdV = dpsidv
							equil_out%profiles_1d%gm1 = g33
       equil_out%profiles_1d%gm4 = bdb02 * btor**2
       equil_out%profiles_1d%gm41 = g41 
       equil_out%profiles_1d%gm5 = b0db2 / btor**2
       equil_out%profiles_1d%surface = surface_x

       equil_out%profiles_1d%rbp_b2 = rbp_b2
       equil_out%profiles_1d%rbp_b2(1) = 0.
       equil_out%profiles_1d%bplfs = bplfs

       equil_out%profiles_1d%acosB2a = acosB2a
       equil_out%profiles_1d%asinB2a = asinB2a
       equil_out%profiles_1d%acosBlnBa = acosBlnBa
       equil_out%profiles_1d%asinBlnBa = asinBlnBa
       
       equil_out%profiles_1d%g1 = g11*4*PI**2
       equil_out%profiles_1d%g2 = g22*4*PI**2
       equil_out%profiles_1d%g2int = g2int*4*PI**2
       equil_out%profiles_1d%fofb = yFOFB
       equil_out%profiles_1d%areat = areat
       equil_out%profiles_1d%perim = perim
       equil_out%profiles_1d%ggradro = gradro*2*PI
       equil_out%profiles_1d%gdroda = droda
       equil_out%profiles_1d%bmaxt = bmaxt
       equil_out%profiles_1d%bmint = bmint
       equil_out%profiles_1d%bdb0 = bdb0 * btor
       
       !geometry
C================== calculation of 3M from SPIDER output geometry
       allocate(ya(iplas),yra(iplas),yri(iplas),
     *  yshif(iplas),yshiv(iplas),
     *  yelon(iplas),ytria_u(iplas),ytria_l(iplas))

	 do ji=2,iplas
		    yrmax=-99999.d0
		    yrmin=99999.d0
		    yzmax=-99999.d0
		    yzmin=99999.d0

	i=minloc(z(ji,2:nt1),1)  
	i=i+1
	i1=i-1
	i2=i+1
	if (i.eq.1) i1=nt1
	if (i.eq.nt1) i2=2
	xxxx1(1)=r(ji,i1)
	xxxx1(2)=r(ji,i)
	xxxx1(3)=r(ji,i2)
	yyyy1(1)=z(ji,i1)
	yyyy1(2)=z(ji,i)
	yyyy1(3)=z(ji,i2)
	call polyfitcc_spid(xxxx1,yyyy1,pppp1)
	yrzmin=-pppp1(2)/(2.*pppp1(1))
	yzmin=pppp1(1)*yrzmin**2.0+pppp1(2)*yrzmin+pppp1(3)
	
	i=maxloc(z(ji,2:nt1),1)  
	i=i+1
	i1=i-1
	i2=i+1
	if (i.eq.1) i1=nt1
	if (i.eq.nt1) i2=2
	xxxx1(1)=r(ji,i1)
	xxxx1(2)=r(ji,i)
	xxxx1(3)=r(ji,i2)
	yyyy1(1)=z(ji,i1)
	yyyy1(2)=z(ji,i)
	yyyy1(3)=z(ji,i2)
	call polyfitcc_spid(xxxx1,yyyy1,pppp1)
	yrzmax=-pppp1(2)/(2.*pppp1(1))
	yzmax=pppp1(1)*yrzmax**2.0+pppp1(2)*yrzmax+pppp1(3)
!	            if(yzmax.le.z(ji,jj)) then 
!		            yzmax=z(ji,jj)
!		            yrzmax=r(ji,jj)
!	            endif
	i=minloc(r(ji,2:nt1),1)  
	yrmin=r(ji,i+1)
	i=maxloc(r(ji,2:nt1),1)  
	yrmax=r(ji,i+1)
!	            if(yrmin.ge.r(ji,jj)) then 
!		            yrmin=r(ji,jj)
!	            endif
!	            if(yrmax.le.r(ji,jj)) then 
!		            yrmax=r(ji,jj)
!	            endif
            yrr=.5d0*(yrmax+yrmin)
            ya(ji)=.5*(yrmax-yrmin)
            yra(ji) = yrmax
            yri(ji) = yrmin
		    yshif(ji)=yrr-rtor
		    yshiv(ji)=.5d0*(yzmin+yzmax)
		    yelon(ji)=(yzmax-yzmin)/(yrmax-yrmin)
		    !ytria(ji)=(yrr-0.5d0*(yrzmin+yrzmax))/ya(ji)
		    ytria_u(ji)=(yrr-yrzmax)/ya(ji)
		    ytria_l(ji)=(yrr-yrzmin)/ya(ji)
	  enddo
	  ya(1)=0.d0
	  yra(1)=r(1,2)
	  yri(1)=r(1,2)
		yelon(1)=yelon(2)
		ytria_u(1)=ytria_u(2)
		ytria_l(1)=ytria_l(2)
		yshif(1)=r(1,2)-rtor
		yshiv(1)=z(1,2)


       if(.NOT. associated(equil_out%profiles_1d%volume)) then
           allocate(equil_out%profiles_1d%volume(iplas))
           allocate(equil_out%profiles_1d%r_inboard(iplas))
           allocate(equil_out%profiles_1d%r_outboard(iplas))
           allocate(equil_out%profiles_1d%elongation(iplas))
           allocate(equil_out%profiles_1d%tria_upper(iplas))
           allocate(equil_out%profiles_1d%tria_lower(iplas))
           
           allocate(equil_out%profiles_1d%shif(iplas))
           allocate(equil_out%profiles_1d%shiv(iplas))
       endif
       equil_out%profiles_1d%volume = volum*2*PI  
       equil_out%profiles_1d%r_inboard = yri
       equil_out%profiles_1d%r_outboard = yra
       equil_out%profiles_1d%elongation = yelon
       equil_out%profiles_1d%tria_upper = ytria_u
       equil_out%profiles_1d%tria_lower = ytria_l 
             
       equil_out%profiles_1d%shif = yshif
       equil_out%profiles_1d%shiv = yshiv       

!get rectangular psi
	call get2d_eq_ef(equil_out)

!get psplex
	equil_out%global_param%psplex=psplex

!Deallocation     
          deallocate( g11, g22, g33 )
          deallocate( g41 )
          deallocate( rbp_b2, bplfs)
          deallocate( g2int , dpsidv)
          deallocate( gradro, droda )
          deallocate( volum )
          deallocate( dvdpsi )
          deallocate( dvdro )
          deallocate( surface_s )
          deallocate( surface_x )
          deallocate( areat )
          deallocate( areats )
          deallocate( perim )
          deallocate( perims )

          deallocate( acosB2a )
          deallocate( asinB2a )
          deallocate( acosBlnBa )
          deallocate( asinBlnBa )

          deallocate( acosB2as )
          deallocate( asinB2as )
          deallocate( acosBlnBas )
          deallocate( asinBlnBas )

          deallocate( g11s, g22s, g33s )
          deallocate( g41s)
          deallocate( g2ints,dpsidvs )
          deallocate( gradros, drodas )
          deallocate( bplfss)
          deallocate( dvoliz)
          deallocate( dvdpsis )
          deallocate( sa )

          deallocate( arr2 )
          
          deallocate( rhowr, rhowrh )
          deallocate( rhos )
          deallocate ( Btot )

          deallocate (b_maxt,b_mint,
     *               b_db02,b_db0,b_0db2)
         
          deallocate (bmaxt,bmint,
     *               bdb02,bdb0,b0db2)
           
          deallocate (traps)
          deallocate (yFOFB)
       deallocate(ya,yra,yri,
     *  yshif,yshiv,
     *  yelon,ytria_u,ytria_l)

       return
       end



!get rectangular part
	subroutine	get2d_eq_ef(equil_out)
	   use imas_ids
		 use parameters_a2spider, only: TWOPI
       	 !implicit none
        INCLUDE 'double.inc'
        INCLUDE 'param.inc'
        INCLUDE 'comblc.inc'
         type(type_equilibrium) equil_out

       if(.NOT. associated(equil_out%eqgeometry%rectgrid%r2d)) then
           allocate(equil_out%eqgeometry%rectgrid%r2d(ni1))
           allocate(equil_out%eqgeometry%rectgrid%z2d(nj1))
           allocate(equil_out%eqgeometry%rectgrid%psirz2d(ni1,nj1))
       endif

			 equil_out%eqgeometry%rectgrid%npointsr=ni1
			 equil_out%eqgeometry%rectgrid%npointsz=nj1
	do j=1,nj1
		do i=1,ni1
			 equil_out%eqgeometry%rectgrid%r2d(i)=r(i)
			 equil_out%eqgeometry%rectgrid%z2d(j)=z(j)
			 equil_out%eqgeometry%rectgrid%psirz2d(i,j)=u(i,j)
		enddo
	enddo

			 equil_out%eqgeometry%rectgrid%psi_axis=um
			 equil_out%eqgeometry%rectgrid%psi_boundary=up

      equil_out%global_param%psibound=-TWOPI*up
      equil_out%global_param%psiaxis=-TWOPI*um


	return
	end

	subroutine polyfitcc_spid(x,y,P)
	implicit none
	double precision x(3),y(3),P(3)
	double precision y21,y32,x21,x32,h21,h32

	y32=y(3)-y(2)				 
	y21=y(2)-y(1)				 
	x32=x(3)-x(2)				 
	x21=x(2)-x(1)				 
	h32=x(3)+x(2)				 
	h21=x(2)+x(1)				 

	P(1)=(x21*y32-x32*y21)/(x21*x32*(h32-h21))
	P(2)=y21/x21-P(1)*h21
	P(3)=y(3)-P(1)*x(3)**2.-P(2)*x(3)


	end

	subroutine polyfitcc_spid_1(x,y,P)
	implicit none
	double precision x(2),y(2),P(2)
	double precision y21,x21,h21

	y21=y(2)-y(1)				 
	x21=x(2)-x(1)				 

	P(1)=y21/x21
	P(2)=y(1)-x(1)*y21/x21


	end

	subroutine polyfitcc_spid_0(x,y,P)
	implicit none
	double precision x(1),y(1),P(1)

	P(1)=y(1)


	end



C======================================================================|
C EXTRAP_EF computes extrapolations 
C======================================================================|
C Assume that x is of r-type, i.e. interpolation in 0 has zero odd derivatives

	subroutine	EXTRAP_EFSP(x_input,y_input,x_extrap,
     &  j_extrap,y_extrap,ex_order,nagrid)
	implicit none
	integer	i,j,k1,k2,k3,ex_order,j_extrap,nagrid
	integer	jsign,k4
	double precision x_input(nagrid),x_extrap
	double precision y_input(nagrid),y_extrap
	double precision y1tmp,y2tmp,y3tmp,drho,y4tmp
	double precision x1tmp,x2tmp,x3tmp,x4tmp
	double precision Acoef,Bcoef,Ccoef,Dcoef
	double precision G,H,F,P(3)
         
				 if (j_extrap.eq.nagrid)  jsign=-1
				 if (j_extrap.eq.1)  jsign=1
						
C Constant interpolation
      if (ex_order.eq.0) then
				   k1=j_extrap
							call polyfitcc_spid_0(x_input(k1),y_input(k1),P(1))
       y_extrap=P(1)
						endif

C Linear interpolation
      if (ex_order.eq.1) then
			   if (jsign.lt.0) then
			   k1=j_extrap+jsign*1
			   k2=j_extrap
				 	call polyfitcc_spid_1(x_input(k1:k2),y_input(k1:k2),P(1:2))
	        y_extrap=P(1)*x_extrap+P(2)
				 endif
			   if (jsign.gt.0) then
			   k1=j_extrap
							call polyfitcc_spid_0(x_input(k1),y_input(k1),P(1))
       y_extrap=P(1)
				 endif
			endif

C Quadratic interpolation
      if (ex_order.eq.2) then
			   if (jsign.lt.0) then
			   k1=j_extrap+jsign*2
			   k2=j_extrap+jsign*1
			   k3=j_extrap
	call polyfitcc_spid(x_input(k1:k3),y_input(k1:k3),P)
				y_extrap=P(1)*x_extrap**2.0+P(2)*x_extrap+P(3)
			  endif
			   if (jsign.gt.0) then
			   k1=j_extrap
			   k2=j_extrap+jsign*1
			   k3=j_extrap+jsign*2
	call polyfitcc_spid(x_input(k1:k3),y_input(k1:k3),P)
				y_extrap=P(1)*x_extrap**2.0+P(2)*x_extrap+P(3)
				 endif
				
				
			endif


	end
C======================================================================|
CEfable
