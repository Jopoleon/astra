	subroutine feqis_main(nucoils,ucoils,parameters_spider,ifplasma & 
	& ,equil_in, &
	&  equil_out)

      use imas_ids       
      use parameters_a2spider
      use circuit, only: ncoils, nrho, nteta, nr2, nz2, &
          voltage, ucoils, &
          psiplasmatoconduc, psirz, psiextrz, psplex, psibndp, &
          psi_cur_old, psiaxisp, &
          psi_external_calc_ef
      implicit none
 
      type(type_equilibrium) equil_in, equil_out
       type(type_parameters) parameters_spider
      integer ifplasma,nucoils
			integer nrplasma
	double precision ucoils(nucoils)
	integer j_call,j_vacplas
	data j_call/0/	
	data j_vacplas/0/	
	integer j_init
	data j_init/0/	
	save j_call,j_init,j_vacplas
	
!allocates output
! defintiions various
	write(*,*) ifplasma
	call definitions_ef_equil(equil_in,parameters_spider,j_init,ifplasma)
	ncoils=nucoils
	voltage(1:ncoils)=ucoils(1:ncoils) ! voltage inputs for active conductors
	write(*,*) ifplasma
nrplasma=nrho
	if (ifplasma.eq.1) then
          allocate(equil_out%profiles_1d%psi(nrplasma))
          allocate(equil_out%profiles_1d%pressure(nrplasma))
          allocate(equil_out%profiles_1d%phi(nrplasma))
          allocate(equil_out%profiles_1d%pprime(nrplasma))
          allocate(equil_out%profiles_1d%ffprime(nrplasma))
          allocate(equil_out%profiles_1d%F_dia(nrplasma))
          allocate(equil_out%profiles_1d%q(nrplasma))
          allocate(equil_out%coord_sys%position%r(nrplasma,nteta))
          allocate(equil_out%coord_sys%position%z(nrplasma,nteta))    
          allocate(equil_out%coord_sys%position%teta2d(nteta))    
          allocate(equil_out%coord_sys%position%rmin(nrplasma,nteta))    
          allocate(equil_out%coord_sys%position%psirz(nrplasma,nteta))    
          allocate(equil_out%eqgeometry%boundary%r(nteta))
          allocate(equil_out%eqgeometry%boundary%z(nteta))
         allocate( equil_out%profiles_1d%rho_tor(nrplasma) )
         allocate( equil_out%profiles_1d%jparallel(nrplasma) )
         allocate( equil_out%profiles_1d%sigmapar%value(nrplasma) )
         allocate( equil_out%profiles_1d%jni%value(nrplasma) )
         allocate( equil_out%profiles_1d%te%value(nrplasma) )
            allocate(equil_out%coord_sys%gradvcell(nrplasma,nteta))
            allocate(equil_out%coord_sys%bpcell(nrplasma,nteta))
            allocate(equil_out%coord_sys%bcell(nrplasma,nteta))
            allocate(equil_out%coord_sys%rcell(nrplasma,nteta))
           allocate(equil_out%profiles_1d%gm1(nrplasma))
           allocate(equil_out%profiles_1d%gm4(nrplasma))
           allocate(equil_out%profiles_1d%gm5(nrplasma))
           allocate(equil_out%profiles_1d%gm41(nrplasma))

           allocate(equil_out%profiles_1d%rbp_b2(nrplasma))
           allocate(equil_out%profiles_1d%bplfs(nrplasma))

           allocate(equil_out%profiles_1d%acosB2a(nrplasma,5))
           allocate(equil_out%profiles_1d%asinB2a(nrplasma,5))
           allocate(equil_out%profiles_1d%acosBlnBa(nrplasma,5))
           allocate(equil_out%profiles_1d%asinBlnBa(nrplasma,5))
           
           allocate(equil_out%profiles_1d%g1(nrplasma))
           allocate(equil_out%profiles_1d%g2(nrplasma))
           allocate(equil_out%profiles_1d%g2int(nrplasma))
           allocate(equil_out%profiles_1d%fofb(nrplasma))
           allocate(equil_out%profiles_1d%areat(nrplasma))
           allocate(equil_out%profiles_1d%perim(nrplasma))
           allocate(equil_out%profiles_1d%ggradro(nrplasma))
           allocate(equil_out%profiles_1d%gdroda(nrplasma))
           allocate(equil_out%profiles_1d%bmaxt(nrplasma))
           allocate(equil_out%profiles_1d%bmint(nrplasma))
           allocate(equil_out%profiles_1d%bdb0(nrplasma))
           allocate(equil_out%profiles_1d%dPSIdV(nrplasma))
           allocate(equil_out%profiles_1d%surface(nrplasma))
           allocate(equil_out%profiles_1d%volume(nrplasma))
           allocate(equil_out%profiles_1d%r_inboard(nrplasma))
           allocate(equil_out%profiles_1d%r_outboard(nrplasma))
           allocate(equil_out%profiles_1d%elongation(nrplasma))
           allocate(equil_out%profiles_1d%tria_upper(nrplasma))
           allocate(equil_out%profiles_1d%tria_lower(nrplasma))
           
           allocate(equil_out%profiles_1d%shif(nrplasma))
           allocate(equil_out%profiles_1d%shiv(nrplasma))
           allocate(equil_out%eqgeometry%rectgrid%r2d(nr2))
           allocate(equil_out%eqgeometry%rectgrid%z2d(nz2))
           allocate(equil_out%eqgeometry%rectgrid%psirz2d(nr2,nz2))
	endif

	write(*,*) 'fix and nstep',parameters_spider%k_fixfree, & 
	& parameters_spider%nstep,nrho,nteta

!init coils and grid
	if (j_call.eq.0) then
	if (parameters_spider%k_fixfree.eq.1) then
	call equil_ef_init_circ
	endif
	endif

	write(*,*) parameters_spider%k_fixfree,j_init
! call fix boundary code, also here boundary comes from experiment
! only circuit equations solved!!!!!!!!!!!!!!!!!!!!!!
	if (parameters_spider%k_fixfree.eq.1) then

	if (ifplasma.eq.0) then
	write(*,*) 'vacuum'
	 psi_cur_old=0.
	 psiplasmatoconduc=0.
 		call circuit_eq_advance_ef(j_call)	
 		call psi_external_calc_ef
		psirz=psiextrz
		j_vacplas=0
	endif
! only circuit equations solved!!!!!!!!!!!!!!!!!!!!!!


! full plasma solved!!!!!!!!!!!!!!!!!!!!!!
	if (ifplasma.eq.1) then
		if (j_vacplas.eq.0) j_call=0
		call full_system_advance_ef(j_call)
		if (j_call.eq.-1) then
			call convert_boundary_to_pbe
			call fix_boundary_ef(1)
		endif
		j_vacplas=1		
		endif
! full plasma!!!!!!!!!!!!!!!!!!!!!!
	endif !kfixfree

	if (parameters_spider%k_fixfree.eq.0) then
			write(*,*) 'call fix equil code'
			call fix_boundary_ef(j_init)
				equil_out%global_param%psplex=psplex
	equil_out%global_param%psibound=psibndp
	equil_out%global_param%psiaxis=psiaxisp
	call assignment_of_equilout_stuff(equil_out)
			write(*,*) 'end fix equil code'
	j_init=1
	write(*,*) 'end equil code'
	return
			endif

	write(*,*) 'end equil code'
	j_call=1
	j_init=1

	if (ifplasma.eq.1) then
		call assignment_of_equilout_stuff(equil_out)
	endif

	return
	end

	subroutine assignment_of_equilout_stuff(equil_out)

      use imas_ids       
      use parameters_a2spider
      use circuit
			use transfer_functions
   use pi_vars, only: GPI2

	implicit none

      type(type_equilibrium) equil_out



!global params
	equil_out%global_param%psplex=psplex
	equil_out%global_param%psibound=-GPI2*psibnd
	equil_out%global_param%psiaxis=-GPI2*psiaxis
	equil_out%global_param%li3=li3
	equil_out%global_param%betpol=betapol	
	equil_out%global_param%i_plasma=iplasma*1.e6

	equil_out%coord_sys%position%psirz(1:nrho,1:nteta)=psirhoteta(1:nrho,1:nteta)/GPI2
	equil_out%coord_sys%position%r(1:nrho,1:nteta)= rpbez(1:nrho,1:nteta)
	equil_out%coord_sys%position%z(1:nrho,1:nteta)= zpbez(1:nrho,1:nteta)
	equil_out%coord_sys%position%teta2d(1:nteta)=t2dbez(1:nteta)
	equil_out%coord_sys%position%rmin(1:nrho,1:nteta)=rmin2dbez(1:nrho,1:nteta)
	equil_out%coord_sys%bpcell(1:nrho,1:nteta)=bpcell2dbez(1:nrho,1:nteta)
	equil_out%coord_sys%bcell(1:nrho,1:nteta)= bcell2dbez(1:nrho,1:nteta)

	equil_out%eqgeometry%rectgrid%npointsr=nr2
	equil_out%eqgeometry%rectgrid%npointsz=nz2
	equil_out%eqgeometry%rectgrid%r2d(1:nr2)=r(1:nr2)
	equil_out%eqgeometry%rectgrid%z2d(1:nz2)=z(1:nz2)	
	equil_out%eqgeometry%rectgrid%psirz2d(1:nr2,1:nz2)=psirz(1:nr2,1:nz2)
	equil_out%eqgeometry%rectgrid%psi_axis=psiaxis
	equil_out%eqgeometry%rectgrid%psi_boundary=psibnd


  equil_out%profiles_1d%ffprime(1:nrho)=0.
  equil_out%profiles_1d%pprime(1:nrho)=0.
  equil_out%profiles_1d%pressure(1:nrho)=0.
  equil_out%profiles_1d%rho_tor(1:nrho)=0.
  equil_out%profiles_1d%F_dia(1:nrho)=0.

	equil_out%profiles_1d%dPSIdV(1:nrho)= dpsidvbez(1:nrho)
	equil_out%profiles_1d%psi(1:nrho)=  psibez(1:nrho) ! psi from spider
	equil_out%profiles_1d%g2(1:nrho)=g2bez(1:nrho)
  equil_out%profiles_1d%g2int(1:nrho)=g2ibez(1:nrho)
  equil_out%profiles_1d%gm1(1:nrho)=gm1bez(1:nrho)
  equil_out%profiles_1d%r_outboard(1:nrho)=routbez(1:nrho)
  equil_out%profiles_1d%r_inboard(1:nrho)=rinbez(1:nrho)
  equil_out%profiles_1d%volume(1:nrho)=vbez(1:nrho)
  equil_out%profiles_1d%g1(1:nrho)=g1bez(1:nrho)
  equil_out%profiles_1d%gm41(1:nrho)=gm41bez(1:nrho)
  equil_out%profiles_1d%ggradro(1:nrho)=ggrhobez(1:nrho)
  equil_out%profiles_1d%bmaxt(1:nrho)=bmaxbez(1:nrho)
  equil_out%profiles_1d%bmint(1:nrho)=bminbez(1:nrho)
  equil_out%profiles_1d%gm4(1:nrho)=gm4bez(1:nrho)
  equil_out%profiles_1d%bdb0(1:nrho)=bdb0bez(1:nrho)
  equil_out%profiles_1d%gm5(1:nrho)=gm5bez(1:nrho)
  equil_out%profiles_1d%fofb(1:nrho)=fofbbez(1:nrho)
  equil_out%profiles_1d%areat(1:nrho)=areatbez(1:nrho)
  equil_out%profiles_1d%perim(1:nrho)=perimbez(1:nrho)
  equil_out%profiles_1d%shif(1:nrho)=shifbez(1:nrho)
  equil_out%profiles_1d%elongation(1:nrho)=kbez(1:nrho)
  equil_out%profiles_1d%surface(1:nrho)= surfbez(1:nrho) ! lateral surface
  equil_out%profiles_1d%tria_upper(1:nrho)=triaubez(1:nrho)
	equil_out%profiles_1d%tria_lower=equil_out%profiles_1d%tria_upper(1:nrho)
  equil_out%profiles_1d%phi(1:nrho)=phibez(1:nrho)
  equil_out%profiles_1d%q(1:nrho)=qbez(1:nrho)
	equil_out%profiles_1d%rbp_b2(1:nrho)=rbp2_b2bez(1:nrho)


	return
	end





	subroutine convert_boundary_to_pbe


	use circuit
	use pi_vars, only: GPI, GPI2

	implicit none	
		double precision teta_fbe(i_dim5)
		double precision x1,x2,x3,x4,t1,t2,t3,t4,z1,z2,z3,z4
		double precision x11,x22,x33,x44,t11,t22,t33,t44,z11,z22,z33,z44
		double precision dx
		integer i,j,k,i1,i2,i3,i4,j1,j2,j3,j4
		
	do i=1,nteta+1
		teta(i)=GPI2*(i-1.)/(nteta+0.)
	enddo
	dteta=teta(2)-teta(1)

	dx =sqrt(dr**2.+dz**2.)

! find boundary
	j=jaxis
	do i=iaxis,nr2
		if (psirz(i,j).ge.psibnd) k=i
		if (psirz(i,j).le.psibnd) goto 10
	enddo	

10 continue

	if (psirz(k,j).eq.psibnd) then
		rbnd(1)=r(k)
	else
		rbnd(1)=r(k)-(psirz(k,j)-psibnd)/(psirz(k,j)-psirz(k-1,j))*dr		
	endif	
		zbnd(1)=z(j)
                teta_fbe(1) = ATAN2(zbnd(1) - zax, rbnd(1) - rax)
			i=1
	write(*,*) i,rbnd(i),zbnd(i),teta_fbe(i)/GPI*180.,x1,t1,t2,t3,z1,z2,z3,psibnd
			
	do i=2,nteta

		dx =sqrt(dr**2.+dz**2.)
		
		x1=sqrt((rbnd(i-1)-rax)**2.+(zbnd(i-1)-zax)**2.)
		teta_fbe(i)=teta_fbe(i-1)+dteta

		t1=rax+x1*cos(teta_fbe(i))
		t2=zax+x1*sin(teta_fbe(i))

		if (t1.ge.raus) then
			x1=(raus-rax)/cos(teta_fbe(i))
		endif
		if (t1.le.rinner) then
			x1=(rinner-rax)/cos(teta_fbe(i))
		endif
		if (t2.ge.ztop) then
			x1=(ztop-zax)/sin(teta_fbe(i))
		endif
		if (t2.le.zbot) then
			x1=(zbot-zax)/sin(teta_fbe(i))
		endif
		t1=rax+x1*cos(teta_fbe(i))
		t2=zax+x1*sin(teta_fbe(i))
		call find_fields_interp_ef_psionly(t1,t2,t3) !give back psi,br,bz at r0,z0
		if (t3.eq.psibnd) then
			rbnd(i)=t1
			zbnd(i)=t2
			goto 2
		endif

		if (t3.lt.psibnd) then
4 continue
			z1=rax+(x1-dx)*cos(teta_fbe(i))
			z2=zax+(x1-dx)*sin(teta_fbe(i))
			call find_fields_interp_ef_psionly(z1,z2,z3) 
			if (z3.lt.psibnd) then
				dx=1.1*dx
				goto 4
			endif
			x11=x1-(t3-psibnd)/(t3-z3)*dx		
			rbnd(i)=rax+x11*cos(teta_fbe(i))
			zbnd(i)=zax+x11*sin(teta_fbe(i))
			goto 2
		endif

		if (t3.gt.psibnd) then
3	continue
			j4=0
			x2=x1+dx
			z1=rax+x2*cos(teta_fbe(i))
			z2=zax+x2*sin(teta_fbe(i))
		if (z1.ge.raus) then
			x2=(raus-rax)/cos(teta_fbe(i))
			j4=1
		endif
		if (z1.le.rinner) then
			x2=(rinner-rax)/cos(teta_fbe(i))
			j4=1
		endif
		if (z2.ge.ztop) then
			x2=(ztop-zax)/sin(teta_fbe(i))
			j4=1
		endif
		if (z2.le.zbot) then
			x2=(zbot-zax)/sin(teta_fbe(i))
			j4=1
		endif
			z1=rax+x2*cos(teta_fbe(i))
			z2=zax+x2*sin(teta_fbe(i))
			call find_fields_interp_ef_psionly(z1,z2,z3) 
			if (z3.gt.psibnd) then
				if (j4.eq.1) then
					rbnd(i)=z1
					zbnd(i)=z2
					goto 2
				endif
				dx=dx*1.1
				goto 3
			endif
			x11=x1+(t3-psibnd)/(t3-z3)*dx		
			rbnd(i)=rax+x11*cos(teta_fbe(i))
			zbnd(i)=zax+x11*sin(teta_fbe(i))
			goto 2
		endif

2	continue

!	write(*,*) i,rbnd(i),zbnd(i),teta_fbe(i)/GPI*180.,x1,t1,t2,t3,z1,z2,z3,psibnd,x11
!	write(8878,*) rbnd(i),zbnd(i)
!	write(88784,*) nteta

	enddo		

	nbnd=nteta

!	open(32,file='fort.607')
!	do i=1,nbnd
!		write(32,'(8E25.11)') rbnd(i),zbnd(i)
!	enddo
!	close(32)


!	if (teta_fbe(1).gt.teta_fbe(nbnd)) then !reorder
!		j=1
!		do k=1,nbnd-1
!			if (teta_fbe(k+1).lt.teta_fbe(k)) j=k+1   ! j is the first teta above 0
!		enddo		
!		if (j.gt.1)	teta_fbe(1:j-1)=teta_fbe(1:j-1)-GPI2
!	endif

	teta(1:nteta)=teta_fbe(1:nteta)
	teta(nteta+1)=teta(nteta)+dteta
	rbndp(1:nteta)=rbnd(1:nteta)
	zbndp(1:nteta)=zbnd(1:nteta)
	rbndp(nteta+1)=rbnd(nteta)
	zbndp(nteta+1)=zbnd(nteta)

	raxp=rax
	zaxp=zax
	psibndp=psibnd
	psiaxisp=psiaxis
	
	return
	end
	

!!!!!!!!!!!!!!!!!!!!!!!!!!! fix boundary

	subroutine fix_boundary_ef(j_init,equil_out)

	use circuit
	use exchange_with_astra
	use imas_ids
	use metric_coefficients_pbe	
	use transfer_functions
        use pi_vars, only: GPI2
	
	implicit none

	integer j_init,i,j
	type(type_equilibrium) equil_out

! initial guess
	if (j_init.eq.0) then
	raxp=raxis_astra
	zaxp=zaxis_astra
	psiaxisp=psi0_astra
	psibndp=psib_astra	
!boundary from experiment
	endif	


	
!boundary from previous time step
!	write(*,*) psiaxisp,psibndp
	call PHI_EQ_2d_PBE( &! Inputs
     &    nrho,nteta,psigrida(1:nrho),iplasma,pressure(1:nrho), &
     &    ffprime(1:nrho),pprime(1:nrho),btor0,rgeom0,&
     &    rbndp(1:nteta),zbndp(1:nteta),Raxp,Zaxp,psiaxisp, &
     &    psibndp,ipol(1:nrho), & 
		 & rpbez(1:nrho,1:nteta), & 
			&	zpbez(1:nrho,1:nteta), & 
			& psirhoteta(1:nrho,1:nteta), & 
			& dpsidvbez(1:nrho), & 
  & psibez(1:nrho), &   ! psi from spider
    & g2bez(1:nrho), &
    & gm1bez(1:nrho), &
   & routbez(1:nrho), &
   & rinbez(1:nrho), &
  & vbez(1:nrho), &
  & g1bez(1:nrho), &
  & gm41bez(1:nrho), &
  & ggrhobez(1:nrho), &
  & bmaxbez(1:nrho), &
  & bminbez(1:nrho), &
  & gm4bez(1:nrho), &
  & bdb0bez(1:nrho), &
  & gm5bez(1:nrho), &
  & fofbbez(1:nrho), &
  & areatbez(1:nrho), &
  & perimbez(1:nrho), &
  & shifbez(1:nrho), &
  & kbez(1:nrho), &
  & surfbez(1:nrho), & ! lateral surface
  & triaubez(1:nrho), &
  & phibez(1:nrho), &
  & qbez(1:nrho), &
	& t2dbez(1:nteta), &
	& rmin2dbez(1:nrho,1:nteta), &
	& bpcell2dbez(1:nrho,1:nteta), &
	& bcell2dbez(1:nrho,1:nteta), & 
	& rbp2_b2bez(1:nrho),jrhoteta(1:nrho,1:nteta), &
	& solve_fix,li3,betapol,dator(1:nrho,1:nteta), &
	& rpol(1:nrho,1:nteta),zpol(1:nrho,1:nteta),psplex)

	rpul(1:nrho,1:nteta)=rpbez(1:nrho,1:nteta)
	zpul(1:nrho,1:nteta)=zpbez(1:nrho,1:nteta)

! additional info from rectangular grid
 
	do j=1,nteta
	do i=1,nrho-1
		rho(i,j)=sqrt((rpol(i,j)-raxp)**2.+(zpol(i,j)-zaxp)**2.)
	enddo
	enddo
	teta(1:nteta)=t2dbez(1:nteta)
	rho(1:nrho-1,nteta+1)=rho(1:nrho-1,1)
	teta(nteta+1)=teta(1)+GPI2

	psia_1d(1:nrho)=psigrida(1:nrho)
	ffp_1d(1:nrho)=ffprime(1:nrho)
	ppp_1d(1:nrho)=pprime(1:nrho)

	return
	end

!---------------------------------------------------------------------
	subroutine solve_gs2d(g)

	use circuit
use pi_vars, only: sintable, costable, mu0

	implicit none
	double precision g(i_dim2,i_dim2)
	double precision gt(i_dim2,i_dim2)
	double precision rhs(i_dim2,i_dim2)
	double precision wrhs(i_dim2,i_dim2)
	double precision trhs(i_dim2,i_dim2)
	double precision A(258),z_fourier(258)
	double precision B(258),x1,x2
	double precision C(258),tin,tup
	integer i,j,k,k_fourier
	integer j_init
	data j_init/0/
	save A,B,C,j_init,z_fourier
	
! g(1,jz), g(nr2,jz), g(jr,1), gr(jr,nz2) contain the b.c.
! solution is g(2:nr1,2:nz1) --> nr*nz points
! solution goes back to g

	do j=1,nz2
	do i=1,nr2
	rhs(i,j)=-mu0*r(i)*jrz(i,j)
	enddo
	enddo
	
	rhs(2,2:nz1)=rhs(2,2:nz1)-g(1,2:nz1)/dr**2.*2./(r(1)+r(2))*r(2)
	rhs(nr1,2:nz1)=rhs(nr1,2:nz1)-g(nr2,2:nz1)/dr**2.*2./(r(nr2)+r(nr1))*r(nr1)
	wrhs=rhs


!	CALL CPU_TIME(tin)
	do i=2,nr1
  	call discrete_sine_transform_ef(nz,wrhs(i,2:nz1))
		do j=2,nz1
			wrhs(i,j)=wrhs(i,j)-g(i,nz2)/dz**2.*sintable(nz,j-1) &
      & -g(i,1)/dz**2.*sintable(1,j-1)
		enddo
	enddo
!	CALL CPU_TIME(tup)
!	write(*,*) 'four transf' ,tup-tin	


	k_fourier = nz
!create inverse matrix for gs2d
	if (j_init.eq.0) then
		A=0.
		B=0.
		C=0.
		do i=2,nr-1		
		x1=0.5*(rcomp(i)+rcomp(i-1))		
		x2=0.5*(rcomp(i+1)+rcomp(i))		
		B(i)=-rcomp(i)/dr**2.*(1./x2+1./x1) 
		A(i)=rcomp(i)/x2/dr**2. 
		C(i)=rcomp(i)/x1/dr**2. 
		enddo
		i=1
		x1=0.5*(rcomp(i)+r(1))		
		x2=0.5*(rcomp(i+1)+rcomp(i))		
		B(i)=-rcomp(i)/dr**2.*(1./x2+1./x1) 
		A(i)=rcomp(i)/x2/dr**2. 
		i=nr
		x1=0.5*(rcomp(i)+rcomp(i-1))		
		x2=0.5*(r(nr2)+rcomp(i))		
		B(i)=-rcomp(i)/dr**2.*(1./x2+1./x1) 
		C(i)=rcomp(i)/x1/dr**2. 
		j_init=1
		do k=2,nz1
			z_fourier(k)=2./dz**2.*(costable(1,k-1)-1.)
		enddo
	endif
	
	gt=0.
!solve matrix
		do k=2,nz1
			call solve_tridiag_fbe_ef(C(1:nr),B(1:nr)+z_fourier(k),A(1:nr),wrhs(2:nr1,k),gt(2:nr1,k),nr)
		enddo
	!invert fourier from gt(1:nr,1:kfourier) to g(2:nr1,2:nz1)
!			gt(i,k)=sum(invMM_gs2d(i-1,1:nr,k-1)*wrhs(2:nr1,k))

	do i=2,nr1
  	call discrete_sine_transform_ef(nz,gt(i,2:nz1))
	enddo
	g(2:nr1,2:nz1)=2./(nz+1)*gt(2:nr1,2:nz1)


!	write(*,*) 'g',g(130,130)
!	pause

!	g(2:nr1,2:nz1)=1./(2*nz+2)*gt(2:nr1,2:nz1)

	return
	end



	subroutine solve_tridiag_fbe_ef(A,B,C,R,f,Ngrid)

! Provides solution of the system:
!
!   Aj fj-1  + Bj fj  + Cj fj+1 = Rj
!
!   where j = 1..Ngrid
!
!   bcbound = 1  -> given f_NA1
!
!  eximp = 2: implicit
!
! by means of this method:

	implicit none
	integer	i,j,k,Ngrid,NgridS,bcbound
	double precision A(Ngrid),B(Ngrid)
	double precision C(Ngrid),R(Ngrid)
	double precision f(Ngrid),f_bound
	double precision alpha(Ngrid),beta(Ngrid),f0
	double precision Bstar,Cstar,Rstar
	integer eximp

	alpha=0.
	beta=0.

         j = 1
			   alpha(j) = C(j)/B(j)
         beta(j) = R(j)/B(j)
			 
			   do j=2,Ngrid-1
			    alpha(j) = C(j)/(B(j)-A(j)*alpha(j-1))
			   enddo
			   do j=2,Ngrid
			    beta(j) = (R(j)-A(j)*beta(j-1))/(B(j)-A(j)*alpha(j-1))
			   enddo

! Note that boundary value is assumed to be on the main last grid point, so 1-dx/2. This has to be
! corrected later on... CEfable
			   j = Ngrid

!C			     f(j-1)=(f_bound-beta(j-1))/alpha(j-1)
!CTHis one should be appropriate with extrapolation... but now go back to real b.c.
!C			     f(j-1)=(2./3.*f_bound-beta(j-1))/(alpha(j-1)-1./3.)
			     f(j)=beta(j)
					 do k=1,Ngrid-1
					  j=Ngrid-1-k+1
					  f(j)=beta(j)-alpha(j)*f(j+1)
					 enddo

	return
	end	

!---------------------------------------------------------------------
	subroutine find_new_axis_part1	
	
	use circuit
	
	implicit none
	
	integer jaold,iaold,i_mode
	integer i,j,k,iax,jax,i_ixpoint
	double precision errtol,tolerr,raxm,zaxm
	double precision br1,bz1,dbrr,dbrz,dbzr,dbzz
	double precision br2,bz2,darax,dazax,bub(100),xub(100),yub(100)
	double precision br3,bz3,icase,psi0,ccc(6)
	double precision br4,bz4,determ,rleft,rright
	double precision x1,x2,x3,x4,x5,x6,x7,x8,x9,psibtmp,x22,x33,x10,x11
	double precision ppx(2)
	integer i1,i2,i3,i4,i5,i6,i7,i8,i9
	integer j1,j2,j3,j4,j5,j6,j7,j8,j9,i0,j0
	!the used function is psirz
	
! 1) find new magnetic axis
!old axis
!	call find_actual_index_ef(rax,zax,iax,jax)
	iax=iaxis
	jax=jaxis
	i1=10000
	j1=10000
	i2=10000
	j2=10000
	errtol=1.e6
	tolerr=10.
	j=1
	raxm=rax
	zaxm=zax
!	write(44125,*) ' '
  do while(errtol.gt.tolerr)
	i1=i2
	j1=j2
	i2=iax
	j2=jax
	if (psirz(iax+1,jax).gt.psirz(iax,jax)) then
		iax=iax+1
		jax=jax
	endif
	if (psirz(iax-1,jax).gt.psirz(iax,jax)) then
		iax=iax-1
		jax=jax
	endif
	if (psirz(iax,jax+1).gt.psirz(iax,jax)) then
		iax=iax
		jax=jax+1
	endif
	if (psirz(iax,jax-1).gt.psirz(iax,jax)) then
		iax=iax
		jax=jax-1
	endif
	if (psirz(iax+1,jax-1).gt.psirz(iax,jax)) then
		iax=iax+1
		jax=jax-1
	endif
	if (psirz(iax-1,jax+1).gt.psirz(iax,jax)) then
		iax=iax-1
		jax=jax+1
	endif
	if (psirz(iax-1,jax-1).gt.psirz(iax,jax)) then
		iax=iax-1
		jax=jax-1
	endif
	if (psirz(iax+1,jax+1).gt.psirz(iax,jax)) then
		iax=iax+1
		jax=jax+1
	endif
	
!	write(44125,*) iax,jax,psirz(iax,jax)

!	call find_actual_index_ef(rax,zax,iax,jax)
	
	j=j+1
!	write(*,*) rax,zax,raxm,zaxm
	if ((iax.eq.i1).and.(jax.eq.j1)) goto 101
	if (j.ge.100000) goto 101		
	enddo
101 continue	

	iaxis=iax
	jaxis=jax	
	write(*,*) 'ax',rax,zax,iax,jax,r(iax),z(jax)
	
	call nine_point_regression(r(iax),z(jax),ppx,derivpsi,psiaxis)

	rax=ppx(1)  !r(iaxis)
	zax=ppx(2) !z(jaxis)
	trax=rax
	tzax=zax


			if (isnan(rax)) then
			write(*,*) 'rax is nan in fbe find axis'
			stop
			endif
	
	return
	end

!---------------------------------------------------------------------
	subroutine find_psi_boundary

	use circuit
use pi_vars, only: GPI
use errors_params, only: err_find_oxpoints_derivs
	
	implicit none
	
	integer jaold,iaold,i_mode,niter
	integer i,j,k,iax,jax,i_ixpoint
	double precision errtol,tolerr,raxm,zaxm
	double precision br1,bz1,dbrr,dbrz,dbzr,dbzz,rstart
	double precision br2,bz2,darax,dazax,bub(9),xub(9),yub(9)
	double precision br3,bz3,icase,psi0,ccc(6),psiloc,polfield
	double precision br4,bz4,determ,rleft,rright,omega_temp(300,300)
	double precision x1,x2,x3,x4,x5,x6,x7,x8,x9,psibtmp,x22,x33,x10,x11
	double precision xbnd(i_dim5),ybnd(i_dim5),r_temp(i_dim2),z_temp(i_dim2),rminz
	double precision pos_xpoint(2),ddipsi(8),tin,tup,frlim_ef
	double precision r_limp,z_limp,psi_limp(500),hard_left,hard_right
	double precision psi_xpoint(max_xpoints),ddpsi(5)
 	integer ipath(5),jpath(5),oldpointnum
	integer i1,i2,i3,i4,i5,i6,i7,i8,i9
	integer j1,j2,j3,j4,j5,j6,j7,j8,j9,i0,j0,i_county
	data i_county/0/
	save i_county,rstart,oldpointnum

	i_plasmatype=0

	if (i_county.eq.1) then

!	go through old x-points and see where they end up
		if (n_of_xpoints.ge.1) then
			iaold=n_of_xpoints
			do i=1,iaold
				niter=0

289	niter=niter+1				
!				call find_actual_index_ef(r_xpoint(i),z_xpoint(i),j,k)
				call nine_point_regression_follow(r_xpoint(i),z_xpoint(i),pos_xpoint,ddipsi,x1)
!				write(5671,*) niter,r_xpoint(i),z_xpoint(i),pos_xpoint
				r_xpoint(i)=pos_xpoint(1)
				z_xpoint(i)=pos_xpoint(2)
				if ( (abs(ddipsi(1))+abs(ddipsi(2))).le.err_find_oxpoints_derivs) goto 290
				if (niter.gt.100000) goto 291 !no point found
				if ((pos_xpoint(1).gt.r(nr2)-dr).or.(pos_xpoint(1).lt.r(1)+dr).or. & 
				 & (pos_xpoint(2).gt.z(nz2)-dz).or.(pos_xpoint(2).lt.z(1)+dz)) then
			! xpoint doesnt exist anymore
			!write(5672,*) 'killed',i
					r_xpoint(i)=1.e6
					z_xpoint(i)=0.
					deriv_x(1:5,i)=1.e6
					goto 317
				endif
				goto 289
291 continue
					r_xpoint(i)=1.e6
					z_xpoint(i)=0.
					deriv_x(1:5,i)=1.e6
					goto 317
290	continue


!			write(5672,*) 'old point new ',i,pos_xpoint(1),pos_xpoint(2)

!check if point outside of domain
				if ((pos_xpoint(1).gt.r(nr2)-dr).or.(pos_xpoint(1).lt.r(1)+dr).or. & 
				 & (pos_xpoint(2).gt.z(nz2)-dz).or.(pos_xpoint(2).lt.z(1)+dz)) then
			! xpoint doesnt exist anymore
			!write(5672,*) 'killed',i
					r_xpoint(i)=1.e6
					z_xpoint(i)=0.
					deriv_x(1:5,i)=1.e6
					goto 317
				endif

				r_xpoint(i)=pos_xpoint(1)
				z_xpoint(i)=pos_xpoint(2)
				deriv_x(1:5,i)=ddpsi(1:5)
!			write(5672,*) 'refound',i,r_xpoint(i),z_xpoint(i)
317	continue
!				write(5671,*) 'old point new ',i,pos_xpoint(1),pos_xpoint(2)
			enddo			
		endif

			j=2
			do i=2,nr1
	call nine_point_regression(r(i),z(j),pos_xpoint,ddipsi,x1)
	
							x5=(ddipsi(5)**2.-ddipsi(3)*ddipsi(4))
	
	if ((pos_xpoint(1).ge.r(i)-dr).and. & 
							& (pos_xpoint(1).le.r(i)+dr).and. & 
				& (pos_xpoint(2).ge.z(j)-dz).and. & 
				& (pos_xpoint(2).le.z(j)+dz).and. & 
				& (x5.ge.0.)) then
	
			n_of_xpoints=min(max_xpoints,n_of_xpoints+1)
			r_xpoint(n_of_xpoints)=pos_xpoint(1)
			z_xpoint(n_of_xpoints)=pos_xpoint(2)
			deriv_x(1:5,n_of_xpoints)=ddpsi(1:5)
	
	endif
			enddo
			j=nz1
			do i=2,nr1
	call nine_point_regression(r(i),z(j),pos_xpoint,ddipsi,x1)
	
							x5=(ddipsi(5)**2.-ddipsi(3)*ddipsi(4))
	
	if ((pos_xpoint(1).ge.r(i)-dr).and. & 
							& (pos_xpoint(1).le.r(i)+dr).and. & 
				& (pos_xpoint(2).ge.z(j)-dz).and. & 
				& (pos_xpoint(2).le.z(j)+dz).and. & 
				& (x5.ge.0.)) then
	
			n_of_xpoints=min(max_xpoints,n_of_xpoints+1)
			r_xpoint(n_of_xpoints)=pos_xpoint(1)
			z_xpoint(n_of_xpoints)=pos_xpoint(2)
			deriv_x(1:5,n_of_xpoints)=ddpsi(1:5)
	
	endif
			enddo
			i=2
			do j=2,nz1
	call nine_point_regression(r(i),z(j),pos_xpoint,ddipsi,x1)
	
							x5=(ddipsi(5)**2.-ddipsi(3)*ddipsi(4))
	
	if ((pos_xpoint(1).ge.r(i)-dr).and. & 
							& (pos_xpoint(1).le.r(i)+dr).and. & 
				& (pos_xpoint(2).ge.z(j)-dz).and. & 
				& (pos_xpoint(2).le.z(j)+dz).and. & 
				& (x5.ge.0.)) then
	
			n_of_xpoints=min(max_xpoints,n_of_xpoints+1)
			r_xpoint(n_of_xpoints)=pos_xpoint(1)
			z_xpoint(n_of_xpoints)=pos_xpoint(2)
			deriv_x(1:5,n_of_xpoints)=ddpsi(1:5)
	
	endif
			enddo
			i=nr1
			do j=2,nz1
	call nine_point_regression(r(i),z(j),pos_xpoint,ddipsi,x1)
	
							x5=(ddipsi(5)**2.-ddipsi(3)*ddipsi(4))
	
	if ((pos_xpoint(1).ge.r(i)-dr).and. & 
							& (pos_xpoint(1).le.r(i)+dr).and. & 
				& (pos_xpoint(2).ge.z(j)-dz).and. & 
				& (pos_xpoint(2).le.z(j)+dz).and. & 
				& (x5.ge.0.)) then
	
			n_of_xpoints=min(max_xpoints,n_of_xpoints+1)
			r_xpoint(n_of_xpoints)=pos_xpoint(1)
			z_xpoint(n_of_xpoints)=pos_xpoint(2)
			deriv_x(1:5,n_of_xpoints)=ddpsi(1:5)
	
	endif
			enddo
!	write(5671,*) 'new points ',n_of_xpoints,r_xpoint(1:n_of_xpoints),z_xpoint(1:n_of_xpoints)
!	write(5672,*) 'new points ',n_of_xpoints,r_xpoint(1:n_of_xpoints),z_xpoint(1:n_of_xpoints)
	endif
!END of tracking of xpoints


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if (i_county.eq.0) then ! do a full pass to find all X-points
		n_of_xpoints=0
		do j=2,nz1
			do i=2,nr1
				call nine_point_regression(r(i),z(j),pos_xpoint,ddipsi,x1)
							x5=(ddipsi(5)**2.-ddipsi(3)*ddipsi(4))
!							write(5671,*) 'xgrid',r(i),z(j),pos_xpoint,ddipsi,x1,dr,dz,x5
							if ((pos_xpoint(1).ge.r(i)-dr).and. & 
							& (pos_xpoint(1).le.r(i)+dr).and. & 
				& (pos_xpoint(2).ge.z(j)-dz).and. & 
				& (pos_xpoint(2).le.z(j)+dz).and. & 
				& (x5.ge.0.)) then

								if (abs(pos_xpoint(1)-rax).le.2.*dr.and.abs(pos_xpoint(2)-zax).le.2.*dz) then
									else
									n_of_xpoints=min(max_xpoints,n_of_xpoints+1)
									r_xpoint(n_of_xpoints)=pos_xpoint(1)
									z_xpoint(n_of_xpoints)=pos_xpoint(2)
									deriv_x(1:5,n_of_xpoints)=ddpsi(1:5)
!									write(5671,*) 'yes point',n_of_xpoints,r_xpoint(n_of_xpoints),z_xpoint(n_of_xpoints)
								endif

							endif
			enddo
		enddo
		i_county=1
		oldpointnum=n_of_xpoints
	endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


		!remove disappeared x-points and double counts
318	i=0
319	i=i+1	
			if (r_xpoint(i).ge.1.e5) then
					if (n_of_xpoints.eq.1) then
						n_of_xpoints=0
						goto 320
					endif
					r_xpoint(i:n_of_xpoints-1)=r_xpoint(i+1:n_of_xpoints)
					z_xpoint(i:n_of_xpoints-1)=z_xpoint(i+1:n_of_xpoints)
					n_of_xpoints=n_of_xpoints-1
					i=i-1
			endif

	do k=1,i-1 ! check if double counted
		if ((abs(r_xpoint(i)-r_xpoint(k)).le.2.*dr).and. & 
		 & (abs(z_xpoint(i)-z_xpoint(k)).le.2.*dz)) then
			r_xpoint(k)=0.5*(r_xpoint(i)+r_xpoint(k))	
			z_xpoint(k)=0.5*(z_xpoint(i)+z_xpoint(k))	
			deriv_x(1:5,k)=0.5*(deriv_x(1:5,i)+deriv_x(1:5,k))
					r_xpoint(i:n_of_xpoints-1)=r_xpoint(i+1:n_of_xpoints)
					z_xpoint(i:n_of_xpoints-1)=z_xpoint(i+1:n_of_xpoints)
					n_of_xpoints=n_of_xpoints-1
					i=i-1
		endif
	enddo


		if (i.lt.n_of_xpoints) goto 319
320	continue

	oldpointnum=n_of_xpoints

!ignore limiter if use_limiter_astra is 0, da trasferirsi in init
	if (use_limiter_yesno.eq.0) then
		limiterR=r(nr1)
		limiterZ=z(nz1)	
	endif	
!calculate limiter flux 
	do i=1,nlimiter
		call find_fields_interp_ef_psionly(limiterR(i),limiterZ(i),psi_limp(i)) !give back psi,br,bz at r0,z0
	enddo


!	write(*,*) limiterr(nlimiter),limiterz(1:nlimiter),psi_limp(1:nlimiter)

! now, remove limiters that are in the shadow of xpoints
	ztop=1.e6
	zbot=-1.e6
	raus=1.e6
	rinner=0.
	if (n_of_xpoints.ge.1) then
		i_plasmatype=1
		do i=1,n_of_xpoints
			call find_actual_index_ef(r_xpoint(i),z_xpoint(i),j,k)
	pos_xpoint(1)=r_xpoint(i)
	pos_xpoint(2)=z_xpoint(i)
	call find_fields_interp_ef_psionly(pos_xpoint(1),pos_xpoint(2),psi_xpoint(i))

!	write(*,*) i,pos_xpoint(1),pos_xpoint(2),psi_xpoint(i)


			if (zlimpotential(j,k).gt.0.) then
                                x1 = ATAN2(z_xpoint(i) - zax, r_xpoint(i) - rax)
				if ((r_xpoint(i).gt.rax).and.(x1.ge.7./4.*GPI.or.x1.le.GPI/4.)) raus=min(raus,r_xpoint(i))
				if ((z_xpoint(i).gt.zax).and.(x1.ge.GPI/4..and.x1.le.3./4.*GPI)) ztop=min(ztop,z_xpoint(i))
				if ((r_xpoint(i).lt.rax).and.(x1.ge.3./4.*GPI.and.x1.le.5./4.*GPI)) rinner=max(rinner,r_xpoint(i))
				if ((z_xpoint(i).lt.zax).and.(x1.ge.5./4.*GPI.and.x1.le.7./4.*GPI)) zbot=max(zbot,z_xpoint(i))
				do j=1,nlimiter
					if (limiterr(j).le.rinner) psi_limp(j)=-1.e6
					if (limiterr(j).ge.raus) psi_limp(j)=-1.e6
					if (limiterz(j).le.zbot) psi_limp(j)=-1.e6
					if (limiterz(j).ge.ztop) psi_limp(j)=-1.e6
!					x1=frlim_ef(deriv_x(1:5,i),limiterz(j),pos_xpoint(1),pos_xpoint(2),rax,zax)					
!					x2=frlim_ef(deriv_x(1:5,i),zax,pos_xpoint(1),pos_xpoint(2),rax,zax)				
!					write(*,*) 'frlim',i,j,x1,x2,limiterr(j),rax	
!					if ((limiterr(j)-x1)*(rax-x2).le.0.) psi_limp(j)=-1.e6

				enddo			
			else
			! xpoint is outside of the limiter, doesnt count
			psi_xpoint(i)=-1.e6
			endif			
		enddo
!psi boundary is the highest of all the values
		x1=maxval(psi_xpoint(1:n_of_xpoints),1)
		x2=maxval(psi_limp(1:nlimiter),1)
		i4=maxloc(psi_xpoint(1:n_of_xpoints),1)
		i5=maxloc(psi_limp(1:nlimiter),1)
		psibnd=max(x1,x2)
		if (x2.gt.x1) i_plasmatype=0
		if (x1.ge.x2) i_plasmatype=1
!	write(*,*) psibnd,x1,x2,i_plasmatype,psi_xpoint(i4),psi_limp(i5),limiterr(i5),limiterz(i5),r_xpoint(i4),z_xpoint(i4)
	else
	!no x-points, take largest limiter flux
		psibnd = maxval(psi_limp(1:nlimiter),1)
		i_plasmatype=0
!	write(*,*) psibnd,i_plasmatype
	endif

!	write(114511,*) psibnd,i_plasmatype,curconduc(12:15),psiplasmatoconduc(12:15),rax,zax,r_xpoint(1:2),z_xpoint(1:2),psi_xpoint(1:2)

!	if (i_plasmatype.eq.1) 
	psibnd=psiaxis+(psibnd-psiaxis)*alpsep


!!!!TEST EFABLE
!			call find_actual_index_ef(-4.2d0,8.0d0,j,k)
!	psibnd=psirz(j,k)
!	write(*,*) 'boundarz ,',j,k,r(j),z(k),psibnd,psirz(j,k) !!
!!!!TEST EFABLE



!normalized flux
	u_n(1:nr2,1:nz2)=(psirz(1:nr2,1:nz2)-psiaxis)/(psibnd-psiaxis)


!	write(*,*) 'raus...',raus,rinner,ztop,zbot


!	open(32,file='fort.44442')
!		write(32,*) r(1:nr2),z(1:nz2),psirz(1:nr2,1:nz2),psibnd
!	close(32)
!	pause

!	write(5671,*) ' a ',psibnd,i_plasmatype,alpsep,n_of_xpoints
!	write(5672,*) ' a ',psibnd,i_plasmatype,alpsep,n_of_xpoints
	do k=1,n_of_xpoints 
!	write(5671,*) 'x',r_xpoint(k),z_xpoint(k),psi_xpoint(k)
!	write(5672,*) 'x',r_xpoint(k),z_xpoint(k),psi_xpoint(k)
	enddo
	
	
!	if (i_plasmatype.eq.0) pause
	
!	write(*,*) i_plasmatype,n_of_xpoints,r_xpoint(1:n_of_xpoints),z_xpoint(1:n_of_xpoints),psibnd,zbot,ztop,raus,rinner


!	open(32,file='fort.44598')
!		write(32,*) psiaxis,psibnd,u_n(1:nr2,1:nz2)
!	close(32)

!	open(32,file='fort.445991')
!		write(32,*) zlimpotential(1:nr2,1:nz2)
!	close(32)

!	open(32,file='fort.45918')
!	do i=1,n_of_xpoints
!					write(32,'(6E25.11)') r_xpoint(i),z_xpoint(i),psi_xpoint(i)
!	enddo
!	close(32)
!	pause
!	close(5671)

	return
	end





	subroutine new_jrz_ef ! calculate new right hand side given new boundary!

	use circuit
	use exchange_with_astra 
use rcurr_zcurr_2def, only: R_curr_2D,Z_curr_2D

	implicit none

	integer i,j,k,i1,i2,i3,i4,i5,j1,j2,j3,j4,j5
	double precision dum1,dum2,dum3,zeta,dumc(i_dim2,i_dim2)
	double precision t1,t2,t3,t4,x,y,alp,bet,gam,det,det0
	integer quadrant
	double precision rpluz,zpluz
	integer ipluz,jpluz,qipluz,qjpluz
	integer ilast,jlast,totpoints,istart,jcallaz
	integer external_griddo_j(90000,2),j_griddo_j
	integer internal_griddo(90000,2),i_griddo_j
	data jcallaz/0/
	save jcallaz

! in entry: rbnd,zbnd,nbnd,psiaxis,psibnd,u_n

	dumc=0.	
	area_eff=dr*dz
	i_griddo_j=0
	j_griddo_j=0

! external_griddo_j(1000,2)
!	internal_griddo(5000,2)

!	open(6712,file='fort.6712')

	totpoints=0
	do quadrant=1,4


! sweep from axis to exterior and fill in the current
!write(*,*) ' quadrant',quadrant

!start from axis position
	i1=floor((rax-rmin)/dr+1.)	
	j1=floor((zax-zmin)/dz+1.)
!write(6712,*)  ' quadrant',quadrant
	if (quadrant.eq.1) then
		i=i1+1
		j=j1+1
		rpluz=1.
		zpluz=1.
		ipluz=1
		jpluz=1
		qipluz=1
		qjpluz=1
	endif
	if (quadrant.eq.2) then
		i=i1
		j=j1+1
		rpluz=-1.
		zpluz=1.
		ipluz=-1
		jpluz=1
		qipluz=0
		qjpluz=1
	endif
	if (quadrant.eq.3) then
		i=i1
		j=j1
		rpluz=-1.
		zpluz=-1.
		ipluz=-1
		jpluz=-1
		qipluz=0
		qjpluz=0
	endif
	if (quadrant.eq.4) then
		i=i1+1
		j=j1
		rpluz=1.
		zpluz=-1.
		ipluz=1
		jpluz=-1
		qipluz=1
		qjpluz=0
	endif
	istart=i


558 continue
	totpoints=totpoints+1
	if (totpoints.gt.nz2*nr2) then
		write(88881,*) 'error in (totpoints.gt.nz2*nr2)'
		stop
	endif
			call fill_in_current(r(i),z(j), & 
					& nrho,psia_2d(1:nrho),ppp_2d(1:nrho),ffp_2d(1:nrho),dumc(i,j),u_n(i,j))
	i_griddo_j=i_griddo_j+1
	internal_griddo(i_griddo_j,1)=i
	internal_griddo(i_griddo_j,2)=j
!	write(6712,*)  r(i),z(j),area_eff(i,j)
		i1=i
		i2=j
		t1=(1.-u_n(i1+ipluz,i2))/(u_n(i1,i2)-u_n(i1+ipluz,i2))
		if (t1.lt.0.or.t1.gt.1.) t1=1.e6
		t2=(1.-u_n(i1,i2+jpluz))/(u_n(i1,i2)-u_n(i1,i2+jpluz))
		if (t2.lt.0.or.t2.gt.1.) t2=1.e6
	if (t2.gt.0..and.t2.le.1.) then
			j_griddo_j=j_griddo_j+1
			external_griddo_j(j_griddo_j,1)=i
			external_griddo_j(j_griddo_j,2)=j+jpluz
	endif

	if (t1.gt.1.e5) then
!move horizontally to the right
		i=i+ipluz
		goto 558
	endif

	ilast=i+ipluz

	j_griddo_j=j_griddo_j+1
	external_griddo_j(j_griddo_j,1)=ilast
	external_griddo_j(j_griddo_j,2)=j

!found boundary, go back, check vertically

	i=istart
59781	continue
		t2=(1.-u_n(i,j+jpluz))/(u_n(i,j)-u_n(i,j+jpluz))
	if (t2.lt.0.or.t2.gt.1.) then
		goto 5581
	endif

	goto 5582

5581	continue
	j=j+jpluz
	goto 558

!found boundary on Z, need to advance 1 more
5582	continue

	j_griddo_j=j_griddo_j+1
	external_griddo_j(j_griddo_j,1)=i
	external_griddo_j(j_griddo_j,2)=j+jpluz


	i=i+ipluz
	istart=i
	if (i.eq.ilast) then
		!terminated quadrant
		goto 559
	endif
	goto 59781



559 continue

	enddo !quadrant cycle
	
	close(6712)

! fill current in external griddo
!	open(4314,file='fort.4314')
!	open(4315,file='fort.4315')
	do i=1,j_griddo_j
!		write(4314,*) external_griddo_j(i,1),external_griddo_j(i,2)
		i1=external_griddo_j(i,1)
		i2=external_griddo_j(i,2)
		t1=(1.-u_n(i1-1,i2))/(u_n(i1,i2)-u_n(i1-1,i2))
		if (t1.lt.0.or.t1.gt.1.) t1=0.
		t2=(1.-u_n(i1,i2+1))/(u_n(i1,i2)-u_n(i1,i2+1))
		if (t2.lt.0.or.t2.gt.1.) t2=0.
		t3=(1.-u_n(i1+1,i2))/(u_n(i1,i2)-u_n(i1+1,i2))
		if (t3.lt.0.or.t3.gt.1.) t3=0.
		t4=(1.-u_n(i1,i2-1))/(u_n(i1,i2)-u_n(i1,i2-1))
		if (t4.lt.0.or.t4.gt.1.) t4=0.
			call fill_in_current(r(i1),z(i2), & 
					& nrho,psia_2d(1:nrho),ppp_2d(1:nrho),ffp_2d(1:nrho),det,0.99999999999999d0)

! defining S1 = dR - d1, S2 = dZ-d2, S3 = dZ-d3, S4 = dR-d4
! Jvacuum = C*Jb
! if 1 point only: C = 1 - S1/dR
! if 2 points: C = 1 - S1*S2/(dR*dZ)
! if 3 points: C = 1 - (S1*S2 + S1*S3)/(2*dR*dZ)
! if 4 points: C = 1 - (S1*S2 + S1*S3 + S3*S4 + S4*S2)/(4*dR*dZ)
! t_j = d_j /(dR or dZ depending on direction)

!	det0=t1**2.+t2**2.+t3**2.+t4**2.-t1*t2-t2*t3-t3*t4-t4*t1
	if (t1+t2+t3+t4.gt.0.) then
		!det0=(t1**2.+t2**2.+t3**2.+t4**2.)/(t1+t2+t3+t4)
		det0=t1+t2+t3+t4-t1*t2-t1*t3-t1*t4-t2*t3-t2*t4-t3*t4
		dumc(i1,i2)=det*det0
	else
		dumc(i1,i2)=0.
	endif

	enddo

		jrz(1:nr2,1:nz2)=dumc(1:nr2,1:nz2)

!rescale current
	dum1=sum(jrz*area_eff) !*area_eff) !*area_eff)*dr*dz
	jrz=jrz/dum1*iplasma

	if (isnan(dum1)) then
	write(88881,*) 'total current is nan in fbe current rescaling'
	stop
	endif
	
	t1=0.
	t2=0.
	t3=0.
	do j=1,nz2
	do i=1,nr1
		t1=t1+r(i)**2.*jrz(i,j)*area_eff(i,j)
		t2=t2+z(j)*jrz(i,j)*area_eff(i,j)
		t3=t3+jrz(i,j)*area_eff(i,j)
	enddo
	enddo

	
	R_curr_2D= sqrt(t1/t3)
	Z_curr_2D=	t2/t3
	
	return
	end





!!


	subroutine fill_in_current(r0,z0,nx,psia_2d,ppp_2d,ffp_2d,dumc,un)

	implicit none
	double precision r0,z0,dr,dz,t1,t2,t3,t4,area_eff,zeta,un
	double precision psibnd,dumc
	integer k,nx
	double precision ppp_2d(nx),ffp_2d(nx),psia_2d(nx)
	
	if (un.le.0.) un=0.

	if (un.le.1.) then
				zeta=un*(nx-1.)+1.
				k = max(1,floor(zeta))
				if (k.eq.nx) then
					dumc=ppp_2d(k)*r0+ &
						& ffp_2d(k)/r0
				else
					dumc=((ppp_2d(k+1)*(zeta-k)+ppp_2d(k)*(k+1.-zeta))*r0+ &
						& (ffp_2d(k+1)*(zeta-k)+ffp_2d(k)*(k+1.-zeta))/r0)
				endif
	else
		dumc=0.0
	endif

	return
	end




















































































































































	
	
























































