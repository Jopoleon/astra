          include 'double.inc'
          include 'dim.inc'
          include 'parevo.inc'
            !PARAMETER ( NJLIM=1550 )
            !PARAMETER ( NPFC0 = 40 )
            !PARAMETER ( NpLIM = 400 )
          include 'dimcp.inc'
          parameter(ni_p=nrp,nj_p=ntp)
          parameter(np=1000,nb1=np+1) 
          parameter(nbtabp=1000,nbtabp2=nbtabp*2)
          parameter(nursp=1000,nursp4=nursp+4,nursp6=nursp4*6)

        real*8 rbtab(nbtabp),zbtab(nbtabp),nbtab
        real*8 pstab(nursp),pptab(nursp),fptab(nursp),nutab
        real*8 Cc(np),CUbs(np),Te(np),Cd(np),pres(np)
        real*8 Cu(np)
        real*8 Cu_out(np)
        real*8 pres_s(np)

        real*8 rzbnd(nbtabp2),eqpf(np),eqff(np),fp(np),ipl,ybetpl,yli3,
     *       rtor,btor,rho(np),roc,
     *       g11(np),g22(np),g33(np),                          
     *       vr(np),vrs(np),slat(np),gradro(np),rocnew, 
     *       mu(np),ipol(np),bmaxt(np),bmint(np),bdb02(np),bdb0(np),
     *	     b0db2(np),droda(np),W_Dj(np),yFOFB(np)                      

        real*8 rout(ni_p,nj_p),zout(ni_p,nj_p)
              character*40 prename
              character*40 eqdfn
       dimension currpf_ref(npfc0),Dgap_ref(ngapp)
       dimension arsp(nrp)
       
!       save R_ref,Z_ref,tok_ref
!       save Dgap_ref,currpf_ref
!       save ipl,key_dmf,rtor
!       save btor
           
! !!!!!!!!!!!!+++++test++++++++++!!!!!!!!!!!!!!!!!!!!!!
!
!          call blic_d(r0,z0, r1,r2,r3,r4, z1,z2,z3,z4,
!     *                             u1,u2,u3,u4, u0, dudr,dudz)
!
!
! !!!!!!!!!!!!+++++test++++++++++!!!!!!!!!!!!!!!!!!!!!!

 !!!!!!!!!!!!+++++temporary++++++++++!!!!!!!!!!!!!!!!!!!!!!
         open(1,file='spidat.da1.NBI_wns02')  
         !open(1,file='spidat.dat')  
         !open(1,file='spidat_0.dat')  !asdex

        ! open(1,file='spidat.dat.fail')
         !open(1,file='spidat.dat-1')
         !open(1,file='spidat.dat.1500s')
         !open(1,file='spidat.dat.500.4%-peaked')
         !open(1,file='spidat.dat.500.4%-flat')
         !open(1,file='spidatp.dat')
         !open(1,file='spidat150509.dat')
         !open(1,file='spidat290509.dat')
         !open(1,file='spidat.dat.flat.5%')
         !open(1,file='spidat.dat-DT-flat-3%-1')
         !open(1,file='spidat.075.iter2008')
         !open(1,file='spidat.083.iter2008')
         !open(1,file='spidat_D.dat')
         !open(1,file='spidat.da1.iter2008')
         !open(1,file='spidat.dat.last-SS')
         !open(1,file='spidat.dat.70')
         !open(1,file='spidat.dat.150')
         !open(1,file='spidat_after_spider.dat')
         !open(1,file='spidat.dat.1480')
         !open(1,file='spidat.dat.ss-k2-0.1440')
         !open(1,file='spidat52.dat')
        ! !open(1,file='spidat.dat.ss-k2-0.1440')
         !open(1,file='spidat100709.dat')
         !open(1,file='spidat.dat.L-n3-p41-t243')
         !open(1,file='spidat.dat.n3.p31.5.t230s')
         !open(1,file='spidat.3.5keV-flat')
         !open(1,file='spidat.4.5keV')
         !open(1,file='spidat.3.5keV')
         !open(1,file='spidat.6keV')
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

       !call prof_AST_SP_L(na1,neql,Te,ArSP,rho,roc)


           hro=rho(2)-rho(1)

!!!!!!!!!!!!+++++temporary++++++++++++++++++!!!!!!!!!!!!!!!!!!!!!!
!            cnor=0.523107
!           do i=1,na1
!            eqpf(i)=eqpf(i)*cnor           
!            eqff(i)=eqff(i)*cnor           
!           enddo
!!!!!!!!!!!!+++++++++++++++++++++++++++++++!!!!!!!!!!!!!!!!!!!!!!
        na=na1-1
        nutab=na1

        eqff(na1)=eqff(na)
                  !neql=200
                  !nteta=110
!!!!!!!!!!!!+++++++++++++++++++++++++++++++!!!!!!!!!!!!!!!!!!!!!!

	    kpr=1   !print in spider
          call  kpr_calc(kpr)
          k_con=0          
          call put_key_con(k_con)
	    prename=''
          kname=4
      call  put_name(prename,kname)

          k_grid= 0   ! k_grid= 0   rect. grid
                      ! k_grid= 1   adap. grid

             k_auto= 1   ! k_auto= 1->   full initialization, 
                         ! k_auto= 0-> preinitialization is assumed to be done
             nstep=0     !nstep=0->initialization,#0->time stepping
             key_dmf=-3  !=1->diff.mag.field, =0->without
             k_fixfree=0   !=0->only fixed boundary spider 
             key_ini=1   ! =1 astra profiles, =0 start from EQDSK and SPIDER profiles
             key_start=1   ! =1 astra profiles, =0 start from EQDSK and SPIDER profiles
             eqdfn='Scen_4.txt'
             key_0stp=1   
             key_pres=1   
             
             dt=1.d-2
             time=0.d0
             !dpsdt=-1.d0
             dpsdt=0.d0

        call aspid_flag(1)
        !call put_key_fix(k_fixfree)

        call astra2spider(neql,nteta,nbnd,rzbnd,key_dmf,
     *                    na1,eqpf,eqff,fp,ipl,
     *                    rtor,btor,rho,roc,nstep,yreler,mu,
     *                    cc,Te,cubc,cd,key_ini,eqdfn,pres,cu,
     *                    key_0stp,key_pres                    )

        !call aspid_flag(0) !test

        call spider(nstep,time,dt,key_dmf,k_grid,k_auto,k_fixfree,
     *              dpsdt,key_start,voltpf)  

        !call cur_avg
        call  wrb
        call wr_spik
        !call flux_state(flx_st)
          ! call field_c       
                    !stop 
!!!!!!!!!!!!!!!!!!!!!!gap test!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!         call get_psix(r_xp,z_xp,psi_xp)
!         rg=5.575d0             !8.2806d0
!         zg=-3.891d0            !0.4665d0
!         rt=rg
!         zt=zg
!         ut=psi_xp
!        
!         call gap(ut,rt,zt,rg,zg)
!
!        D_gap=dsqrt((rg-rt)**2+(zg-zt)**2)
!
!!!!!!!!!!!!!!!!!!!!!!gap test!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        call spider2astra(rout,zout,rtor,btor,rho,roc,na1,
     *                    g11,g22,g33,vr,vrs,slat,gradro,rocnew,
     *                    mu,ipol,bmaxt,bmint,bdb02,b0db2,bdb0,droda,
     *	    yreler,yli3,ni_p,nj_p,platok,cu_out,fp,pres_s,W_Dj,
     *                  yFOFB )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!time !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!! steping!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        do nstep=1,1500
  
         !if(nstep.ge.4) then
         ! key_dmf=-3
         ! k_fixfree=0
         !else 
         ! key_dmf=-2
         !endif  
              !call put_key_fix(k_fixfree)
 
              time=time+dt
              call put_tim(dt,time)
              call savepsi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              do i=1,na1
!                eqpf(i)=eqpf(i)/amu0
!                eqff(i)=eqff(i)/amu0
!              enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if(nstep.eq.3) then
              do i=1,na1
                !pres(i)=pres(i)*1.2
              enddo
         endif
         
        call astra2spider(neql,nteta,nbnd,rzbnd,key_dmf,
     *                    na1,eqpf,eqff,fp,ipl,
     *                    rtor,btor,rho,roc,nstep,yreler,mu,
     *                    cc,Te,cubs,cd,key_ini,eqdfn,pres,cu,
     *                    key_0stp,key_pres                    )


        call spider(nstep,time,dt,key_dmf,k_grid,k_auto,k_fixfree,
     *              dpsdt,key_ini)  
          call  pla_volt(dt)
          call cur_avg
          call  wrb

	  !call get_conp(R_ref,Z_ref,tok_ref,Dgap_ref,currpf_ref)
 	  !call put_conp(R_ref,Z_ref,tok_ref,Dgap_ref,currpf_ref)
        
         !write(*,*) 'before NEWGRD_s'
         !write(*,*) 'roc hro nb1 na1 btor hroa '
         !write(*,*) roc, hro, nb1, na1, btor, hroa 

	  !call	NEWGRD_ss( ROC, HRO, NB1, NA1, RHO,btor, HROA  )
	  !call get_roc(ROC,btor)
	  !call	NEWGRD
         
         !write(*,*) 'after NEWGRD_s'
         !write(*,*) 'roc hro nb1 na1 btor hroa '
         !write(*,*) roc, hro, nb1, na1, btor, hroa 
         !write(*,*) 'rho'
         !write(*,*) (rho(i),i=1,na1)
 
        call spider2astra(rout,zout,rtor,btor,rho,roc,na1,
     *                    g11,g22,g33,vr,vrs,slat,gradro,rocnew,
     *                    mu,ipol,bmaxt,bmint,bdb02,b0db2,bdb0,droda,
     *				    yreler,yli3,ni_p,nj_p,platok,cu_out,fp,pres_s,W_Dj,
     *                  yFOFB )




        enddo


        stop
        end
