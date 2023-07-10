!!!!!!!!!!!!!!!!!!!!!!!!!!! init_circ
	subroutine generate_files_feqis(data_dir2)

	use ef_circuit
	use green_matrix
		
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
	double precision dummatrix(52,52) !Efable test
	integer nctype(5200) !Efable test
	double precision dhoriz(200),dvert(200) !blanket elements lengths
	integer equivtmp(5200)
	integer equivforce(5200)
	double precision tatmp(5200)
	character*80 data_dir2



!Efable test, remove this later
	data_dir=data_dir2
	mu0=0.4*GPI
!!!!!!!!!!!!

	open(32,file=trim(data_dir)//'refit_ef.dat')
	read(32,*) nr_of_fit_parameters     ! nr of fit parameters from rmag to zxp_fit
	read(32,*) rmag_fit   ! R mag axis
	read(32,*) zmag_fit ! Z mag axis
	read(32,*) k_fit ! elongation at mag axis
	read(32,*) rxp_fit ! X point R
	read(32,*) zxp_fit ! X point Z
	close(32)


	open(32,file=trim(data_dir)//'data.dat')
	read(32,*) j_files
	read(32,*) nr
	read(32,*) nz
	read(32,*) rmin
	read(32,*) rmax
	read(32,*) zmin
	read(32,*) zmax
	read(32,*) alpsep
	close(32)
	j_files=0

!compatibility with spider
! in spider ni=65, ni1=64
	nr=nr-2  !because for spider its 65,65 for example, but here its 63,63
	nz=nz-2  ! this was -1 before!!!!!



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
	open(32,file=trim(data_dir)//'coil.dat')
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
	open(32,file=trim(data_dir)//'coilres.dat')
	read(32,*) nreseqcoil
	do i=1,nreseqcoil
	read(32,*) j
	read(32,*) resconduc(i,1:nreseqcoil)
	enddo
	close(32)
	nactive=nconduc
	write(*,*) nactive 

!load limiter
	open(32,file=trim(data_dir)//'limpnt.dat')
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
	open(32,file=trim(data_dir)//'blanfw.dat')
	read(32,*) reswall
	read(32,*) nfirstwall
	if (nfirstwall.ge.1) then
	do i=1,nfirstwall
	read(32,*) j,rwall(i)
	enddo
	endif
	close(32)

!load blanket (works)
	open(32,file=trim(data_dir)//'blanbp.dat')
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
	dhoriz(j)=x7/2.
	dvert(j)=x8/2.
	dhoriz(2*i)=x7/2.
	dvert(2*i)=x8/2.
	areablan(j)=x9*widthblan
	areablan(2*i)=x9*widthblan
	ssfw=ssfw+areablan(j)/rblan(j)+areablan(2*i)/rblan(2*i)
!	write(*,*) i,x1,x2,x3,x4,x7,x8,j,areablan(j),areablan(2*i),2*i
	enddo
!	write(*,*) resblan,nblanket,x1,x2,x3,x4,x5,x6,x7,x8,rblan(1),zblan(1), &
!	& areablan(1)
	nblanket=2*nblanket
	do i=1,nblanket
	nconduc=nconduc+1
	curconduc(nconduc)=0.
	r_cond(nconduc)=rblan(i)
	z_cond(nconduc)=zblan(i)
	resconduc(nconduc,nconduc)=resblan*r_cond(nconduc)/areablan(i)*ssfw
!	write(*,*) nconduc,areablan(nconduc),resconduc(nconduc,nconduc),ssfw,resblan
	enddo
	endif
	close(32)

!load passive conduc  (works)
	open(32,file=trim(data_dir)//'blanbpc.dat')
	read(32,*) nblanketpc
	nelemblanketpc=9  !nelemblanketpc sub element of a blanket element, hardwired to 9 for now
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

!create circuit structure
	nconduc=0	
	nblocks=0
	ielem=0
	tatmp=0.
! coils first
	identcoil(1:ncoils)=1
	do i=1,ncoils
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
		do j=i+1,ncoils
			if ((mequivalence(i).eq.mequivalence(j))) then
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
		x1=tempcoilr(j)-0.5*(tempcoildr(j)+tempcoildz(j)*cos(tempcoilangle(j))/sin(tempcoilangle(j)))
		x2=tempcoilz(j)-0.5*(tempcoildz(j))
		r1=x1
		z1=x2
		
		x1=tempcoildr(j)
		x2=tempcoildz(j)/sin(tempcoilangle(j))
		
		x3=tempcoildr(j)
		x4=0.
		x5=tempcoildz(j)*cos(tempcoilangle(j))/sin(tempcoilangle(j))
		x6=tempcoildz(j)
		
		 iii=nint(sqrt(tempcoilelem(j)*x1/x2)+0.5)
		 jjj=nint(sqrt(tempcoilelem(j)*x2/x1)+0.5)

			x3=x3/iii
			x4=x4/iii
			x5=x5/jjj
			x6=x6/jjj	

		 dr1=x3 !tempcoildr(j)/iii
		 dz1=x6 !tempcoildz(j)/jjj
!		 x1=tempcoildz(j)/sin(tempcoilangle(j))
!		 x2=x1*cos(tempcoilangle(j))
!		 r1=tempcoilr(j)-tempcoildr(j)/2.-x2/4.
!		 z1=tempcoilz(j)-tempcoildz(j)/2.
			r1=r1+0.5*(x3+x5)
			z1=z1+0.5*(x4+x6)
	   do jj=1,jjj
	   do ii=1,iii
			ielem=ielem+1
!			rcetmp(ielem)=r1+(ii-1.)*dr1+(jj-1.)*dz1*cos(tempcoilangle(j))+dr1/2.	 !center					  
!			zcetmp(ielem)=z1+(jj-1.)*dz1+dz1/2. !center
			rcetmp(ielem)=r1+(ii-1.)*x3+(jj-1.)*x5	 !center					  
			zcetmp(ielem)=z1+(ii-1.)*x4+(jj-1.)*x6 !center
			drcetmp(ielem)=dr1	 !center					  
			dzcetmp(ielem)=dz1 !center
			datmp(ielem)=dr1*dz1/sin(tempcoilangle(j)) !area including turns
			equivtmp(ielem)=nconduc						  
			equivforce(ielem)=tempnnc(j)						  
			tatmp(ielem)=datmp(ielem)*tempcoilturns(j)/(tempcoildr(j)*tempcoildz(j)/sin(tempcoilangle(j))) !area including turns
	write(*,'(6E25.11)') ielem+0.,rcetmp(ielem),zcetmp(ielem),0.+equivtmp(ielem),0.+equivforce(ielem)
			nctype(ielem)=2 !rectangular coil block
		 enddo								  
		 enddo !cycles over 1 single coil element
		enddo !cycle over equivalent coils
	endif
	enddo !cycle over equivalent coils
!	stop

!now blanket
!write(*,*) nblanketpc
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
			nctype(ielem)=2 !rectangular passive element
		 enddo								  
		 enddo !cycles over 1 single coil element
		enddo !cycle over equivalent coils
	endif


!now blanket
!write(*,*) nblanket
	if (nblanket.ge.1) then
	do j =1,nblanket
		nconduc=nconduc+1
		nblocks=nblocks+1
!	double precision rcetmp(1200),zcetmp(1200),datmp(1200),equivtmp(1200)
			ielem=ielem+1
			rcetmp(ielem)=rblan(j)	 !center					  
			zcetmp(ielem)=zblan(j) !center
			drcetmp(ielem)=dhoriz(j)	 !center					  
			dzcetmp(ielem)=dvert(j) !center
			datmp(ielem)=dhoriz(j)*dvert(j) !area including turns
			equivtmp(ielem)=nconduc
			equivforce(ielem)=nblocks						  
			tatmp(ielem)=1. !area including turns
			nctype(ielem)=1 !segment block
		enddo !cycle over equivalent coils
	endif




!calculate self-inductances ! Fable, trying to match SPIDER inductances....
	write(*,*) 'self',nconduc,npassive,nactive,ielem
	indconduc=0.
	areactmp=0.
	identcoil=1
	jjelem=0
	do i=1,ielem
				iii=equivtmp(i)
			do j=1,ielem
				if (iii.eq.equivtmp(j)) then
					if ((equivforce(j).eq.equivforce(i)).and.(i.ne.j)) then			
						call green_function_non_identity(rcetmp(i),zcetmp(i),rcetmp(j),zcetmp(j),gtemp & 
					& ,drcetmp(i),dzcetmp(i),drcetmp(j),dzcetmp(j),nctype(i),nctype(j))
						indconduc(iii,iii)=indconduc(iii,iii)+ & 
								& mu0/GPI*gtemp*tatmp(i)*tatmp(j)
					endif
					if ((equivforce(j).ne.equivforce(i))) then			
						call green_function_non_identity(rcetmp(i),zcetmp(i),rcetmp(j),zcetmp(j),gtemp & 
					& ,drcetmp(i),dzcetmp(i),drcetmp(j),dzcetmp(j),nctype(i),nctype(j))
						indconduc(iii,iii)=indconduc(iii,iii)+ & 
								& mu0/GPI*gtemp*tatmp(i)*tatmp(j)
					endif
					if (i.eq.j) then			
						jjelem(iii,iii)=jjelem(iii,iii)+1
						call green_function_identity(rcetmp(i),zcetmp(i),gtemp,drcetmp(i),dzcetmp(i),nctype(i))
						indconduc(iii,iii)=indconduc(iii,iii)+ & 
								& mu0/GPI*gtemp*tatmp(i)*tatmp(j)/2.
					endif
				endif
			enddo		
	enddo
	do i=1,nconduc
	indconduc(i,i)=GPI2*indconduc(i,i)
	enddo
	write(*,*) indconduc(1,1)/GPI2
	
!calculate mutual-inductances
	write(*,*) 'mutual'
	identcoil=1
	do i=1,ielem
			iii=equivtmp(i)
			do j=1,ielem
				if ((equivtmp(j).ne.iii)) then			
					jjelem(equivtmp(j),iii)=jjelem(equivtmp(j),iii)+1
						call green_function_non_identity(rcetmp(i),zcetmp(i),rcetmp(j),zcetmp(j),gtemp & 
					& ,drcetmp(i),dzcetmp(i),drcetmp(j),dzcetmp(j),nctype(i),nctype(j))
!					call green_function(rcetmp(i),zcetmp(i),rcetmp(j),zcetmp(j),gtemp)
					indconduc(iii,equivtmp(j))=indconduc(iii,equivtmp(j))+ & 
							& mu0/GPI*gtemp*tatmp(i)*tatmp(j)
		!			indconduc(equivtmp(j),iii)=indconduc(iii,equivtmp(j))
				endif
			enddo		
	enddo
	do j=1,nconduc
	do i=1,nconduc
	if (i.ne.j) then
	indconduc(i,j)=GPI2*indconduc(i,j)
	endif
	enddo
	enddo
	
	write(*,*) indconduc(1,2)/GPI2
	

!calculate grid-inductances
	write(*,*) 'grid'
	greeni=0.
	do i=1,ielem
			iii=equivtmp(i)
			do jj=1,nz2
				do ii=1,nr2
						call green_function_non_identity(rcetmp(i),zcetmp(i),r(ii),z(jj),gtemp & 
					& ,drcetmp(i),dzcetmp(i),dr,dz,nctype(i),2)
!					call green_function(rcetmp(i),zcetmp(i),r(ii),z(jj),gtemp)
					greeni(ii,jj,iii)=greeni(ii,jj,iii)+ & 
							& mu0/GPI*gtemp*tatmp(i)
				enddo		
			enddo
	enddo
	
	
!calculate inter-blocks forces, to check if these are right btw
dgreenirpl=0.
dgreenizpl=0.
dgreenirj=0.
dgreenizj=0.

	do i=1,ielem
		iii=equivforce(i)	
			do j=1,ielem
				if ((equivforce(j).ne.iii)) then
				
				call green_function(rcetmp(i)+dr/2.,zcetmp(i),rcetmp(j),zcetmp(j),x1)
				call green_function(rcetmp(i)-dr/2.,zcetmp(i),rcetmp(j),zcetmp(j),x2)
					dgreenirj(iii,equivforce(j))=dgreenirj(iii,equivforce(j)) & 
					    & -mu0*2.*tatmp(i)*tatmp(j)* &
 							& (x1-x2)/dr

				call green_function(rcetmp(i),zcetmp(i)+dz/2.,rcetmp(j),zcetmp(j),x1)
				call green_function(rcetmp(i),zcetmp(i)-dz/2.,rcetmp(j),zcetmp(j),x2)
					dgreenizj(iii,equivforce(j))=dgreenizj(iii,equivforce(j)) & 
					    & -mu0*2.*tatmp(i)*tatmp(j)* &
 							& (x1-x2)/dz

				endif			
			enddo
	enddo	
	
	write(*,*) 'grid forces from plasma to coil'
	do i=1,ielem
			iii=equivforce(i)
			do jj=1,nz2
				do ii=1,nr2

				call green_function(rcetmp(i)+dr/2.,zcetmp(i),r(ii),z(jj),x1)
				call green_function(rcetmp(i)-dr/2.,zcetmp(i),r(ii),z(jj),x2)
					dgreenirpl(ii,jj,iii)=dgreenirpl(ii,jj,iii) & 
					    & -mu0*2.*tatmp(i)* &
 							& (x1-x2)/dr

				call green_function(rcetmp(i),zcetmp(i)+dz/2.,r(ii),z(jj),x1)
				call green_function(rcetmp(i),zcetmp(i)-dz/2.,r(ii),z(jj),x2)
					dgreenizpl(ii,jj,iii)=dgreenizpl(ii,jj,iii) & 
					    & -mu0*2.*tatmp(i)* &
 							& (x1-x2)/dz

				enddo		
			enddo
	enddo

!!!!!!!!!!!!!!!

!generate zlimpotential
	zlimpotential=1.	
	do j=1,nz2
	do i=1,nr2   
		if (i.lt.ilim_minR) zlimpotential(i,j)=0
		if (i.gt.ilim_maxR) zlimpotential(i,j)=0
		if (j.lt.ilim_minZ) zlimpotential(i,j)=0
		if (j.gt.ilim_maxZ) zlimpotential(i,j)=0
	enddo
	enddo

! write everything on file
	open(32,file=trim(data_dir)//'induc_matrix.dat')	
	write(32,*) nconduc
	do i=1,nconduc	
		write(32,*) indconduc(i,1:nconduc)
	enddo
	close(32)
	open(32,file=trim(data_dir)//'resistance_matrix.dat')	
	write(32,*) nconduc
	do i=1,nconduc	
		write(32,*) resconduc(i,1:nconduc)
	enddo
	close(32)


	open(32,file=trim(data_dir)//'greeni_matrix.dat')
	do i=1,nconduc	
	do j=1,nr2
		write(32,*) greeni(j,1:nz2,i)	
	enddo	
	enddo
	close(32)

!force matrix
	open(32,file=trim(data_dir)//'force_matrix.dat')
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

	open(32,file=trim(data_dir)//'zlim_potential.dat')
	do jj=1,nz2
	do ii=1,nr2
		write(32,*) zlimpotential(ii,jj)	
	enddo	
	enddo
	close(32)




	write(*,*) 'bound'
	open(32,file=trim(data_dir)//'green_boundary.dat')
	write(32,*) nint((2.*nr+2.*nz)*(2.*nr+2.*nz))
! lower side
	do i=2,nr1
		do j=2,nr1
			call	green_function_includingsamepoint(r(i),z(1),r(j),z(1),dr,r(i),greenf)
			write(32,*) greenf
		enddo
!right side	
		do j=2,nz1
			call	green_function_includingsamepoint(r(i),z(1),r(nr2),z(j),dz,r(i),greenf)
			write(32,*) greenf
		enddo
! upper side
		do j=2,nr1
			call	green_function_includingsamepoint(r(i),z(1),r(j),z(nz2),dr,r(i),greenf)
			write(32,*) greenf
		enddo
!left side	
		do j=2,nz1
			call	green_function_includingsamepoint(r(i),z(1),r(1),z(j),dz,r(i),greenf)
			write(32,*) greenf
		enddo
	enddo

! right side
	do i=2,nz1
		do j=2,nr1
			call	green_function_includingsamepoint(r(nr2),z(i),r(j),z(1),dr,r(nr2),greenf)
			write(32,*) greenf
		enddo
!right side	
		do j=2,nz1
			call	green_function_includingsamepoint(r(nr2),z(i),r(nr2),z(j),dz,r(nr2),greenf)
			write(32,*) greenf
		enddo
! upper side
		do j=2,nr1
			call	green_function_includingsamepoint(r(nr2),z(i),r(j),z(nz2),dr,r(nr2),greenf)
			write(32,*) greenf
		enddo
!left side	
		do j=2,nz1
			call	green_function_includingsamepoint(r(nr2),z(i),r(1),z(j),dz,r(nr2),greenf)
			write(32,*) greenf
		enddo
	enddo



! upper side
	do i=2,nr1
		do j=2,nr1
			call	green_function_includingsamepoint(r(i),z(nz2),r(j),z(1),dr,r(i),greenf)
			write(32,*) greenf
		enddo
!right side	
		do j=2,nz1
			call	green_function_includingsamepoint(r(i),z(nz2),r(nr2),z(j),dz,r(i),greenf)
			write(32,*) greenf
		enddo
! upper side
		do j=2,nr1
			call	green_function_includingsamepoint(r(i),z(nz2),r(j),z(nz2),dr,r(i),greenf)
			write(32,*) greenf
		enddo
!left side	
		do j=2,nz1
			call	green_function_includingsamepoint(r(i),z(nz2),r(1),z(j),dz,r(i),greenf)
			write(32,*) greenf
		enddo
	enddo



! left side
	do i=2,nz1
		do j=2,nr1
			call	green_function_includingsamepoint(r(1),z(i),r(j),z(1),dr,r(1),greenf)
			write(32,*) greenf
		enddo
!right side	
		do j=2,nz1
			call	green_function_includingsamepoint(r(1),z(i),r(nr2),z(j),dz,r(1),greenf)
			write(32,*) greenf
		enddo
! upper side
		do j=2,nr1
			call	green_function_includingsamepoint(r(1),z(i),r(j),z(nz2),dr,r(1),greenf)
			write(32,*) greenf
		enddo
!left side	
		do j=2,nz1
			call	green_function_includingsamepoint(r(1),z(i),r(1),z(j),dz,r(1),greenf)
			write(32,*) greenf
		enddo
	enddo
	
	close(32)
	
	
	
	write(*,*) 'files done,please restart'
!	stop !Efable test



	
		return
	end
















