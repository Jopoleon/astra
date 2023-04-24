	subroutine definitions_ef_equil(equil_in,params,j_call,ifplasma)

	use parameters_a2spider
	use imas_ids
	use ef_circuit
	use exchange_with_astra
	
	implicit none
	
	integer j_call,i,j,k,ifplasma
	
	type(type_parameters) params
       type(type_equilibrium) equil_in
	
	if (j_call.eq.0) then
		data_dir=params%prename &
    &  (1:params%kname)
	data_dir_k=params%kname
	nteta=equil_in%eqgeometry%boundary%npoints
	nrho=params%neql
	psistabR=0.
	psistabZ=0.
	!normalized psi from 0 axis to 1 edge, equispaced
	do i=1,nrho
	psigrid(i)=(i-1.)/(nrho-1.)
!	psia_2d(i)=(i-1.)/(nrho-1.)
	enddo

	dr_factor_init=dr_factor_init_astra
	dz_factor_init=dz_factor_init_astra

!errors
	open(32,file=data_dir(1:data_dir_k)//'errors_data.dat')	
	read(32,*) err_circ_plasma_iter
	read(32,*) err_find_oxpoints
	read(32,*) err_find_oxpoints_derivs
	read(32,*) err_find_psistab
	read(32,*) err_find_delr
	read(32,*) err_find_biquad
	read(32,*) err_epsilon
	read(32,*) err_gaptolez
	read(32,*) err_fix_boundary
	close(32)

	!constants
!	GPI=3.141592653589793
!	GPI2=2.*GPI
!	GPI4=GPI2**2.0
	mu0=0.4*GPI
		Rgeom0=equil_in%global_param%toroid_field%r0
	tau_circuit_ef=0.001 !default value	
	tau_gseq_ef=0.001	 !default value
	activate_coil_ef=1. ! when 0., coil is forced to 0 current
	current_limit_ef(:,1)=1e6 ! cant be higher than 1e6 MA
	current_limit_ef(:,2)=-1e6 ! cant be lower than -1e6 MA
	max_iter=10000 !hardwired
! teta for polar grid, goes from 0 to 2*pi-dteta, but point nt+1 is the periodic one
	omega_pl=0.
	psi0_astra=equil_in%profiles_1d%psi(1)
	psib_astra=equil_in%profiles_1d%psi(nrho)
	raxp=raxis_astra
	zaxp=zaxis_astra	
	psibnd=-1.e6
	use_limiter_yesno=use_limiter_astra
	voltage_old=0.
	voltage=0.
	endif

	
	if (ifplasma.eq.1) then

	dteta=GPI2/(nteta+1+0.)
	do i=1,nteta+1
		teta(i)=GPI2*(i-1.)/(nteta+0.)
	enddo

		btor0=equil_in%global_param%toroid_field%b0
		iplasma=equil_in%global_param%i_plasma/1.e6
!	write(*,*) nrho
	pressure(1:nrho)=equil_in%profiles_1d%pressure(1:nrho)
	pprime(1:nrho)=equil_in%profiles_1d%pprime(1:nrho)
	ffprime(1:nrho)=equil_in%profiles_1d%ffprime(1:nrho)
	psigrida(1:nrho)=equil_in%profiles_1d%psi(1:nrho) !unnormalized
	psigrida(1:nrho)=(psigrida(1:nrho)-psigrida(1))/(psigrida(nrho)-psigrida(1)) !normalized 0 axis, 1 boundary

	do i=1,nrho
		psia_2d(i)=(i-1.)/(nrho-1.)
	enddo
	call linterp_ef(psigrida(1:nrho),ffprime(1:nrho),nrho,psia_2d(1:nrho),ffp_2d(1:nrho),nrho)
	call linterp_ef(psigrida(1:nrho),pprime(1:nrho),nrho,psia_2d(1:nrho),ppp_2d(1:nrho),nrho)
	ffp_2d=-GPI2/mu0*ffp_2d
	ppp_2d=-GPI2*1.e-6*ppp_2d

!	open(32,file='fort.367891')
!		write(32,*) psia_2d(1:nrho),ppp_2d(1:nrho),ffp_2d(1:nrho)
!	close(32)

	ipol(1:nrho)=equil_in%profiles_1d%F_dia(1:nrho)

	if (params%k_fixfree.eq.0) then !if 1, comes from free boundary
		rexp(1:nteta)=equil_in%eqgeometry%boundary%r(1:nteta)
		zexp(1:nteta)=equil_in%eqgeometry%boundary%z(1:nteta)
		do i=1,nteta
			call find_angle_ef(raxp,zaxp,rexp(i),zexp(i),tetaexp(i))
		enddo	
!	write(*,*) 'psigrid',psigrida(1:nrho),raxp,zaxp
!	write(*,*) 'psigrid',rexp(1:nteta),zexp(1:nteta),tetaexp(1:nteta)
	!reorder
	if (tetaexp(1).gt.tetaexp(nteta)) then !reorder
		j=1
		do k=1,nteta-1
			if (tetaexp(k+1).lt.tetaexp(k)) j=k+1   ! j is the first teta above 0
		enddo		
		if (j.gt.1)	tetaexp(1:j-1)=tetaexp(1:j-1)-GPI2
	endif
	do i=1,nteta
	rexp(nteta+i)=rexp(i)
	zexp(nteta+i)=zexp(i)
		tetaexp(nteta+i)=tetaexp(i)+GPI2
	enddo
		call linterp_ef_equilef(tetaexp(1:nteta*2),rexp(1:nteta*2),nteta*2, &
     &   teta(1:nteta),rbndp(1:nteta),nteta)
		call linterp_ef_equilef(tetaexp(1:nteta*2),zexp(1:nteta*2),nteta*2, &
     &   teta(1:nteta),zbndp(1:nteta),nteta)
!	open(32,file='fort.765')
!	do i=1,nteta
!		write(32,*) rbndp(i),zbndp(i)
!	enddo
!	close(32)

	endif
	endif

	return
	end
	
	
	
	
	



























































!!!!!!!!!!!!!!!!!!!!!!!!!!! init_circ
	subroutine equil_ef_init_circ

	use ef_circuit
	use green_matrix
	use exchange_with_astra
	use fenix_params
		
	implicit none

	integer j_files,i,j,k,ii,jj,kk,iii,jjj,kkk,ielem
	double precision dummy1,k_fourier
	double precision tempcoilr(300),tempcoilz(300), &
     & tempcoilangle(300),tempcoildr(300),tempcoildz(300)
 integer identcoil(5200),tempcoilturns(300),tempcoilelem(300)
	integer numeqcump(100),jjelem(300,300),tempnnc(300)
	double precision r1,r2,z1,z2,r3,z3,r4,z4,area,gtemp
	double precision dr1,dr2,dz1,dz2,dr3,dz3,dr4,dz4,area2,gtemp2
	double precision x1,x2,x3,x4,x5,x6,x7,x8,x9,areactmp(300)
	double precision rcetmp(5200),zcetmp(5200),datmp(5200)
	double precision drcetmp(5200),dzcetmp(5200),areazz(300)
	double precision matrix_gs2d(500,500),z_fourier,greenf,ssfw
	integer equivtmp(5200)
	integer equivforce(5200)
	double precision tatmp(5200)





	open(32,file=data_dir(1:data_dir_k)//'refit_ef.dat')
	read(32,*) nr_of_fit_parameters     ! nr of fit parameters from rmag to zxp_fit
	read(32,*) rmag_fit   ! R mag axis
	read(32,*) zmag_fit ! Z mag axis
	read(32,*) k_fit ! elongation at mag axis
	read(32,*) rxp_fit ! X point R
	read(32,*) zxp_fit ! X point Z
	close(32)


	open(32,file=data_dir(1:data_dir_k)//'data.dat')
	read(32,*) j_files
	read(32,*) nr
	read(32,*) nz
	read(32,*) rmin
	read(32,*) rmax
	read(32,*) zmin
	read(32,*) zmax
	read(32,*) alpsep
	close(32)

	

!compatibility with spider
	nr=nr-1  !because for spider its 65,65 for example, but here its 64,64
	nz=nz-1



!define grid
	nr2=nr+2
	nz2=nz+2
	nr1=nr+1
	nz1=nz+1
 do i=1,nr2
 	r(i)=rmin+(i-1.)*(rmax-rmin)/nr1     ! computational domain is r(2:nr+1), boundaries are r(1) and r(nr+2)
 enddo
 do i=1,nz2
 	z(i)=zmin+(i-1.)*(zmax-zmin)/nz1
 enddo
	rcomp(1:nr)=r(2:nr1)
	zcomp(1:nz)=z(2:nz1)
	dr=r(2)-r(1)
	dz=z(2)-z(1)

		do i=1,nz
		do j=1,nz
			sintable(i,j)=sin(i*j*GPI/(nz+1))
			costable(i,j)=cos(i*j*GPI/(nz+1))
		enddo
		enddo


!load coils
	r_cond=0.
	z_cond=0.
	numeqcump=0
	open(32,file=data_dir(1:data_dir_k)//'coil.dat')
	read(32,*) ncoils
	do i=1,ncoils
	read(32,*) nelemcoil(i)
	read(32,*) rcoil(i),zcoil(i),drcoil(i),dzcoil(i),dummy1,anglecoil(i), & 
	& curcoil(i),mturns(i),mequivalence(i)
		curconduc(mequivalence(i))=curcoil(i) !assign current to conductor
		r_cond(mequivalence(i))=r_cond(mequivalence(i))+rcoil(i) !assign current to conductor
		z_cond(mequivalence(i))=z_cond(mequivalence(i))+zcoil(i) !assign current to conductor
    numeqcump(mequivalence(i))=numeqcump(mequivalence(i))+1
	enddo
	anglecoil=anglecoil/180.*GPI
	close(32)
	nconduc=maxval(mequivalence(1:ncoils))
	r_cond(1:nconduc)=r_cond(1:nconduc)/numeqcump(1:nconduc)
	z_cond(1:nconduc)=z_cond(1:nconduc)/numeqcump(1:nconduc)

	nblocks=0

!load coilres
	open(32,file=data_dir(1:data_dir_k)//'coilres.dat')
	read(32,*) nreseqcoil
	do i=1,nreseqcoil
	read(32,*) j
	read(32,*) resconduc(i,1:nreseqcoil)
	enddo
	close(32)
	nactive=nconduc
	write(*,*) nactive 

!load limiter
	open(32,file=data_dir(1:data_dir_k)//'limpnt.dat')
	read(32,*) nlimiter
	do i=1,nlimiter
	read(32,*) limiterr(i),limiterz(i)
	enddo
	read(32,*) lim_maxR,lim_minR,lim_maxZ,lim_minZ
	close(32)
	
!forces limiter to adapt to grid points !EFable test, but this should be better than leaving it floating
	do i=1,nlimiter
			call find_actual_index_ef(limiterr(i),limiterz(i),j,k)
			limiterr(i)=r(j)
			limiterz(i)=z(k)
	enddo	

	ilim_maxR=1+nint((lim_maxR-r(1))/dr)
	lim_maxR=r(ilim_maxR)

	ilim_minR=1+nint((lim_minR-r(1))/dr)
	lim_minR=r(ilim_minR)

	ilim_maxZ=1+nint((lim_maxZ-z(1))/dz)
	lim_maxZ=z(ilim_maxZ)

	ilim_minZ=1+nint((lim_minZ-z(1))/dz)
	lim_minZ=z(ilim_minZ)
	
	


!load first wall (not working yet)
	open(32,file=data_dir(1:data_dir_k)//'blanfw.dat')
	read(32,*) reswall
	read(32,*) nfirstwall
	if (nfirstwall.ge.1) then
	do i=1,nfirstwall
	read(32,*) j,rwall(i)
	enddo
	endif
	close(32)

!load blanket (works)
	open(32,file=data_dir(1:data_dir_k)//'blanbp.dat')
	read(32,*) resblan,widthblan
	read(32,*) nblanket
	nelemblanket=9  !nelemblanket sub element of a blanket element, for now hardwired width to 10 cm
	if (nblanket.ge.1) then
	ssfw=0.
	do i=1,nblanket
	j=2*i-1
	read(32,*) jjj,x1,x2,x3,x4,x5,x6
	x7=x3-x1
	x8=x4-x2
	rblan(j)=x1+1./4.*x7
	rblan(2*i)=x1+3./4.*x7
	zblan(j)=x2+1./4.*x8
	zblan(2*i)=x2+3./4.*x8
	x9=sqrt(x7**2.+x8**2.)
	areablan(j)=x9/2.*widthblan
	areablan(2*i)=x9/2.*widthblan
	ssfw=ssfw+areablan(j)/rblan(j)+areablan(2*i)/rblan(2*i)/2.
	write(*,*) i,x1,x2,x3,x4,x7,x8,j,areablan(j),areablan(2*i),2*i
	enddo
	write(*,*) resblan,nblanket,x1,x2,x3,x4,x5,x6,x7,x8,rblan(1),zblan(1), &
	& areablan(1)
	nblanket=2*nblanket
	do i=1,nblanket
	nconduc=nconduc+1
	curconduc(nconduc)=0.
	r_cond(nconduc)=rblan(i)
	z_cond(nconduc)=zblan(i)
	resconduc(nconduc,nconduc)=resblan*r_cond(nconduc)/areablan(i)*ssfw
	write(*,*) nconduc,areablan(nconduc),resconduc(nconduc,nconduc),ssfw,resblan
	enddo
	endif
	close(32)

!load passive conduc  (works)
	open(32,file=data_dir(1:data_dir_k)//'blanbpc.dat')
	read(32,*) nblanketpc
	nelemblanketpc=9  !nelemblanketpc sub element of a blanket element
	if (nblanketpc.ge.1) then
	do i=1,nblanketpc
	read(32,*) rblanpc(i),zblanpc(i),resblanpc(i),areablanpc(i),curblanpc(i)
	nconduc=nconduc+1
	curconduc(nconduc)=curblanpc(i)
	resconduc(nconduc,nconduc)=resblanpc(i)
	r_cond(nconduc)=rblanpc(i)
	z_cond(nconduc)=zblanpc(i)
	enddo
	endif
	close(32)

! currents are in MA!
	npassive = nconduc-nactive
	nblocks=ncoils+npassive

	if (j_files.eq.0) then
		
!create circuit structure
	nconduc=0	
	nblocks=0
	ielem=0
	tatmp=0.
! coils first
	identcoil(1:ncoils)=1
	do i =1,ncoils
	nblocks=nblocks+1
	k=0
	if (identcoil(i).eq.1) then
		nconduc=nconduc+1
		k=k+1
		tempcoilr(k)=rcoil(i)		
		tempcoilz(k)=zcoil(i)		
		tempcoilangle(k)=anglecoil(i)		
		tempcoildr(k)=drcoil(i)		
		tempcoildz(k)=dzcoil(i)		
		tempcoilturns(k)=mturns(i)		
		tempcoilelem(k)=nelemcoil(i)		
		tempnnc(k)=i
		do j=1,ncoils
			if ((mequivalence(i).eq.mequivalence(j)).and. & 
			 & (i.ne.j)) then
					identcoil(j)=0
					k=k+1
		tempcoilr(k)=rcoil(j)		
		tempcoilz(k)=zcoil(j)		
		tempcoilangle(k)=anglecoil(j)		
		tempcoildr(k)=drcoil(j)		
		tempcoildz(k)=dzcoil(j)		
		tempcoilturns(k)=mturns(j)		
		tempcoilelem(k)=nelemcoil(j)
		tempnnc(k)=j
			 endif
		enddo
		!analysis coil	
!	double precision rcetmp(1200),zcetmp(1200),datmp(1200),equivtmp(1200)
		do j=1,k
		 iii=nint(sqrt(tempcoilelem(j)+1.d-6))
		 dr1=tempcoildr(j)/iii
		 dz1=tempcoildz(j)/iii
		 x1=tempcoildz(j)/sin(tempcoilangle(j))
		 x2=x1*cos(tempcoilangle(j))
		 r1=tempcoilr(j)-tempcoildr(j)/2.-x2/4.
		 r2=tempcoilr(j)+tempcoildr(j)/2.-x2/4.
		 z1=tempcoilz(j)-tempcoildz(j)/2.
		 z2=tempcoilz(j)-tempcoildz(j)/2.
		 r3=r2+x2
		 r4=r1+x2
		 z3=tempcoilz(j)+tempcoildz(j)/2.
		 z4=tempcoilz(j)+tempcoildz(j)/2.
	   do jjj=1,iii
	   do jj=1,iii
			ielem=ielem+1
			rcetmp(ielem)=r1+(jj-1.)*dr1+(jjj-1.)*dz1*cos(tempcoilangle(j))+dr1/2.	 !center					  
			zcetmp(ielem)=z1+(jjj-1.)*dz1+dz1/2. !center
			drcetmp(ielem)=dr1	 !center					  
			dzcetmp(ielem)=dz1 !center
			datmp(ielem)=dr1*dz1/sin(tempcoilangle(j)) !area including turns
			equivtmp(ielem)=nconduc						  
			equivforce(ielem)=tempnnc(j)						  
	tatmp(ielem)=datmp(ielem)*tempcoilturns(j)/(tempcoildr(j)*tempcoildz(j)/sin(tempcoilangle(j))) !area including turns
		 enddo								  
		 enddo !cycles over 1 single coil element
		enddo !cycle over equivalent coils
	endif
	enddo !cycle over equivalent coils
	write(333,*) rcetmp(1:ielem),zcetmp(1:ielem),0,0

!now blanket
write(*,*) nblanketpc
	if (nblanketpc.ge.1) then
	do j =1,nblanketpc
		nconduc=nconduc+1
		nblocks=nblocks+1
!	double precision rcetmp(1200),zcetmp(1200),datmp(1200),equivtmp(1200)
		 iii=nint(sqrt(nelemblanketpc+1.d-6))
		 dr1=sqrt(areablanpc(j))/iii
		 dz1=dr1
		 x1=dz1
		 x2=0
		 r1=rblanpc(j)-dr1/2.
		 r2=rblanpc(j)+dr1/2.
		 z1=zblanpc(j)-dz1/2.
		 z2=zblanpc(j)-dz1/2.
		 r3=r2
		 r4=r1
		 z3=zblanpc(j)+dz1/2.
		 z4=zblanpc(j)+dz1/2.
	   do jjj=1,iii
	   do jj=1,iii
			ielem=ielem+1
			rcetmp(ielem)=r1+(jj-1./2.)*dr1	 !center					  
			zcetmp(ielem)=z1+(jjj-1./2.)*dz1 !center
			drcetmp(ielem)=dr1	 !center					  
			dzcetmp(ielem)=dz1 !center
			datmp(ielem)=dr1*dz1 !area including turns
			equivtmp(ielem)=nconduc						  
			equivforce(ielem)=nblocks						  
	tatmp(ielem)=dr1*dz1/areablanpc(j) !area including turns
		 enddo								  
		 enddo !cycles over 1 single coil element
		enddo !cycle over equivalent coils
	endif


!now blanket
write(*,*) nblanket
	if (nblanket.ge.1) then
	do j =1,nblanket
		nconduc=nconduc+1
		nblocks=nblocks+1
!	double precision rcetmp(1200),zcetmp(1200),datmp(1200),equivtmp(1200)
		 iii=nint(sqrt(nelemblanket+1.d-6))
		 dr1=sqrt(areablan(j))/iii
		 dz1=dr1
		 x1=dz1
		 x2=0
		 r1=rblan(j)-dr1/2.
		 r2=rblan(j)+dr1/2.
		 z1=zblan(j)-dz1/2.
		 z2=zblan(j)-dz1/2.
		 r3=r2
		 r4=r1
		 z3=zblan(j)+dz1/2.
		 z4=zblan(j)+dz1/2.
	   do jjj=1,iii
	   do jj=1,iii
			ielem=ielem+1
			rcetmp(ielem)=r1+(jj-1./2.)*dr1	 !center					  
			zcetmp(ielem)=z1+(jjj-1./2.)*dz1 !center
			drcetmp(ielem)=dr1	 !center					  
			dzcetmp(ielem)=dz1 !center
			datmp(ielem)=dr1*dz1 !area including turns
			equivtmp(ielem)=nconduc
			equivforce(ielem)=nblocks						  
	tatmp(ielem)=dr1*dz1/areablan(j) !area including turns
			write(*,*) 'ielem',ielem,rcetmp(ielem),zcetmp(ielem),iii,dr1,dz1,x1,x2,r1,r2,z1,z2,z3,z4						  
		 enddo								  
		 enddo !cycles over 1 single coil element
		enddo !cycle over equivalent coils
	endif




!calculate self-inductances
	write(*,*) 'self',nconduc,npassive,nactive,ielem
	indconduc=0.
	areactmp=0.
	identcoil=1
	jjelem=0
	do i=1,ielem
			iii=equivtmp(i)
			do j=1,ielem
				if ((equivtmp(j).eq.iii).and.(i.ne.j)) then			
					write(*,*) i,iii,rcetmp(i),j,rcetmp(j),zcetmp(i),zcetmp(j)
					call green_function(rcetmp(i),zcetmp(i),rcetmp(j),zcetmp(j),gtemp)
					indconduc(iii,iii)=indconduc(iii,iii)+ & 
							& mu0/GPI*gtemp*tatmp(i)*tatmp(j)
				endif
				if ((equivtmp(j).eq.iii).and.(i.eq.j)) then			
					jjelem(iii,iii)=jjelem(iii,iii)+1
					call green_function_identity(rcetmp(i),zcetmp(i),gtemp,drcetmp(i),dzcetmp(i))
					indconduc(iii,iii)=indconduc(iii,iii)+ & 
							& mu0/GPI*gtemp/2*tatmp(i)*tatmp(j) 
				endif
			enddo		
	enddo
	do i=1,nconduc
	indconduc(i,i)=GPI2*indconduc(i,i)
	enddo

	
!calculate mutual-inductances
	write(*,*) 'mutual'
	identcoil=1
	do i=1,ielem
			iii=equivtmp(i)
			do j=1,ielem
				if ((equivtmp(j).ne.iii)) then			
					jjelem(equivtmp(j),iii)=jjelem(equivtmp(j),iii)+1
					call green_function(rcetmp(i),zcetmp(i),rcetmp(j),zcetmp(j),gtemp)
					indconduc(iii,equivtmp(j))=indconduc(iii,equivtmp(j))+ & 
							& mu0/GPI*gtemp*tatmp(i)*tatmp(j)
					indconduc(equivtmp(j),iii)=indconduc(iii,equivtmp(j))
				endif
			enddo		
	enddo
	do j=1,nconduc
	do i=1,nconduc
	if (i.ne.j) then
	indconduc(i,j)=GPI2*indconduc(i,j)/2.
	endif
	enddo
	enddo
	
	

!calculate grid-inductances
	write(*,*) 'grid'
	greeni=0.
	do i=1,ielem
			iii=equivtmp(i)
			do jj=1,nz2
				do ii=1,nr2
					call green_function(rcetmp(i),zcetmp(i),r(ii),z(jj),gtemp)
					greeni(ii,jj,iii)=greeni(ii,jj,iii)+ & 
							& mu0/GPI*gtemp*tatmp(i)
				enddo		
			enddo
	enddo
	
	
!calculate inter-blocks forces
dgreenirpl=0.
dgreenizpl=0.
dgreenirj=0.
dgreenizj=0.

	do i=1,ielem
		iii=equivforce(i)	
			do j=1,ielem
				if ((equivforce(j).ne.iii)) then
				
				call green_function(rcetmp(i),zcetmp(i),rcetmp(j)+dr/2.,zcetmp(j),x1)
				call green_function(rcetmp(i),zcetmp(i),rcetmp(j)-dr/2.,zcetmp(j),x2)
					dgreenirj(iii,equivforce(j))=dgreenirj(iii,equivforce(j)) & 
					    & -mu0/GPI*GPI2*tatmp(i)*tatmp(j)* &
 							& (x1-x2)/dr

				call green_function(rcetmp(i),zcetmp(i),rcetmp(j),zcetmp(j)+dz/2.,x1)
				call green_function(rcetmp(i),zcetmp(i),rcetmp(j),zcetmp(j)-dz/2.,x2)
					dgreenizj(iii,equivforce(j))=dgreenizj(iii,equivforce(j)) & 
					    & -mu0/GPI*GPI2*tatmp(i)*tatmp(j)* &
 							& (x1-x2)/dz

!					dgreenirj(equivforce(j),iii)=dgreenirj(iii,equivforce(j))
!					dgreenizj(equivforce(j),iii)=dgreenizj(iii,equivforce(j))

				endif			
			enddo
	enddo	
	
	write(*,*) 'grid forces'
	do i=1,ielem
			iii=equivforce(i)
			do jj=1,nz2
				do ii=1,nr2

				call green_function(r(ii),z(jj),rcetmp(i)+dr/2.,zcetmp(i),x1)
				call green_function(r(ii),z(jj),rcetmp(i)-dr/2.,zcetmp(i),x2)
					dgreenirpl(ii,jj,iii)=dgreenirpl(ii,jj,iii) & 
					    & -mu0/GPI*GPI2*tatmp(i)* &
 							& (x1-x2)/dr

				call green_function(r(ii),z(jj),rcetmp(i),zcetmp(i)+dz/2.,x1)
				call green_function(r(ii),z(jj),rcetmp(i),zcetmp(i)-dz/2.,x2)
					dgreenizpl(ii,jj,iii)=dgreenizpl(ii,jj,iii) & 
					    & -mu0/GPI*GPI2*tatmp(i)* &
 							& (x1-x2)/dz

				enddo		
			enddo
	enddo

!!!!!!!!!!!!!!!

!generate zlimpotential
	zlimpotential=1.	
	do j=1,nz2
	do i=1,nr2   
!	call generate_zlim_potential(r(i),z(j),nlimiter,limiterR(1:nlimiter), & 
!	& limiterZ(1:nlimiter),zlimpotential(i,j))
		if (i.lt.ilim_minR) zlimpotential(i,j)=0
		if (i.gt.ilim_maxR) zlimpotential(i,j)=0
		if (j.lt.ilim_minZ) zlimpotential(i,j)=0
		if (j.gt.ilim_maxZ) zlimpotential(i,j)=0
		if (use_zlim_pot.eq.0) zlimpotential(i,j)=1.
	enddo
	enddo

! write everything on file
	open(32,file=data_dir(1:data_dir_k)//'induc_matrix.dat')	
	write(32,*) nconduc
	do i=1,nconduc	
		write(32,*) indconduc(i,1:nconduc)
	enddo
	close(32)
	open(32,file=data_dir(1:data_dir_k)//'resistance_matrix.dat')	
	write(32,*) nconduc
	do i=1,nconduc	
		write(32,*) resconduc(i,1:nconduc)
	enddo
	close(32)
	open(32,file=data_dir(1:data_dir_k)//'greeni_matrix.dat')
	do i=1,nconduc	
	do j=1,nr2
		write(32,*) greeni(j,1:nz2,i)	
	enddo	
	enddo
	close(32)

!force matrix
	open(32,file=data_dir(1:data_dir_k)//'force_matrix.dat')
		write(32,*) nblocks	
	do j=1,nblocks	
	do i=1,nblocks
		write(32,*) dgreenirj(i,j),dgreenizj(i,j)	
	enddo	
	enddo
	do ii=1,nblocks
	do j=1,nz2	
	do i=1,nr2
		write(32,*) dgreenirpl(i,j,ii),dgreenizpl(i,j,ii)	
	enddo	
	enddo
	enddo
	close(32)

	open(32,file=data_dir(1:data_dir_k)//'zlim_potential.dat')
	do jj=1,nz2
	do ii=1,nr2
		write(32,*) zlimpotential(ii,jj)	
	enddo	
	enddo
	close(32)




	write(*,*) 'bound'
	open(32,file=data_dir(1:data_dir_k)//'green_boundary.dat')
	write(32,*) nint((2.*nr+2.*nz)*(2.*nr1+2.*nz1))
! lower side
	do i=2,nr1
		do j=1,nr1
			call green_function(r(i),z(1),r(j)+dr/2.,z(1),greenf)
			write(32,*) greenf
		enddo
!right side	
		do j=1,nz1
			call green_function(r(i),z(1),r(nr2),z(j)+dz/2.,greenf)
			write(32,*) greenf
		enddo
! upper side
		do j=1,nr1
			call green_function(r(i),z(1),r(j)+dr/2.,z(nz2),greenf)
			write(32,*) greenf
		enddo
!left side	
		do j=1,nz1
			call green_function(r(i),z(1),r(1),z(j)+dz/2.,greenf)
			write(32,*) greenf
		enddo
	enddo

! right side
	do i=2,nz1
		do j=1,nr1
			call green_function(r(nr2),z(i),r(j)+dr/2.,z(1),greenf)
			write(32,*) greenf
		enddo
!right side	
		do j=1,nz1
			call green_function(r(nr2),z(i),r(nr2),z(j)+dz/2.,greenf)
			write(32,*) greenf
		enddo
! upper side
		do j=1,nr1
			call green_function(r(nr2),z(i),r(j)+dr/2.,z(nz2),greenf)
			write(32,*) greenf
		enddo
!left side	
		do j=1,nz1
			call green_function(r(nr2),z(i),r(1),z(j)+dz/2.,greenf)
			write(32,*) greenf
		enddo
	enddo



! upper side
	do i=2,nr1
		do j=1,nr1
			call green_function(r(i),z(nz2),r(j)+dr/2.,z(1),greenf)
			write(32,*) greenf
		enddo
!right side	
		do j=1,nz1
			call green_function(r(i),z(nz2),r(nr2),z(j)+dz/2.,greenf)
			write(32,*) greenf
		enddo
! upper side
		do j=1,nr1
			call green_function(r(i),z(nz2),r(j)+dr/2.,z(nz2),greenf)
			write(32,*) greenf
		enddo
!left side	
		do j=1,nz1
			call green_function(r(i),z(nz2),r(1),z(j)+dz/2.,greenf)
			write(32,*) greenf
		enddo
	enddo



! left side
	do i=2,nz1
		do j=1,nr1
			call green_function(r(1),z(i),r(j)+dr/2.,z(1),greenf)
			write(32,*) greenf
		enddo
!right side	
		do j=1,nz1
			call green_function(r(1),z(i),r(nr2),z(j)+dz/2.,greenf)
			write(32,*) greenf
		enddo
! upper side
		do j=1,nr1
			call green_function(r(1),z(i),r(j)+dr/2.,z(nz2),greenf)
			write(32,*) greenf
		enddo
!left side	
		do j=1,nz1
			call green_function(r(1),z(i),r(1),z(j)+dz/2.,greenf)
			write(32,*) greenf
		enddo
	enddo
	
	close(32)
	
	
	
	write(*,*) 'files done,please restart'
	call err_catch_a




!
	else
	
! load everything from file
	open(32,file=data_dir(1:data_dir_k)//'induc_matrix.dat')	
	read(32,*) nconduc
	do i=1,nconduc	
		read(32,*) indconduc(i,1:nconduc)	
	enddo
	close(32)
	open(32,file=data_dir(1:data_dir_k)//'resistance_matrix.dat')	
	read(32,*) nconduc
	do i=1,nconduc	
		read(32,*) resconduc(i,1:nconduc)	
	enddo
	close(32)
	open(32,file=data_dir(1:data_dir_k)//'greeni_matrix.dat')
	do i=1,nconduc	
	do j=1,nr2
		read(32,*) greeni(j,1:nz2,i)	
	enddo	
	enddo
	close(32)
	open(32,file=data_dir(1:data_dir_k)//'green_boundary.dat')
	read(32,*) ngbnd
	read(32,*) green_bnd_f(1:ngbnd)
	close(32)


	open(32,file=data_dir(1:data_dir_k)//'zlim_potential.dat')
	do jj=1,nz2
	do ii=1,nr2
		read(32,*) zlimpotential(ii,jj)	
	enddo	
	enddo
	close(32)
	
	
	!force matrix
	open(32,file=data_dir(1:data_dir_k)//'force_matrix.dat')
		read(32,*) nblocks	
	do j=1,nblocks	
	do i=1,nblocks
		read(32,*) dgreenirj(i,j),dgreenizj(i,j)	
	enddo	
	enddo
	do ii=1,nblocks
	do j=1,nz2	
	do i=1,nr2
		read(32,*) dgreenirpl(i,j,ii),dgreenizpl(i,j,ii)	
	enddo	
	enddo
	enddo
	close(32)

	
	

!use spider values !!!!!!!!!! temporary, to improve my calculations
! 	open(32,file=data_dir(1:data_dir_k)//'ppind_mat.wr')	
!	read(32,*) nconduc
!		read(32,*) ((indconduc(i,j),i=1,nconduc),j=1,nconduc)
!	close(32)
!	indconduc=indconduc*GPI2
!use spider values
! 	open(32,file=data_dir(1:data_dir_k)//'res_mat.wr')	
!		read(32,*) ((resconduc(i,j),i=1,nconduc),j=1,nconduc)
!	close(32)
!	resconduc=resconduc*GPI2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	





	endif



	
		return
	end


























































