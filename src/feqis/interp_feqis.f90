subroutine find_angle_ef(rt,zt,r,z,anglr)


	use pi_grec_vars
! angle from (rt,zt) to (r,z)
	
	implicit none
	double precision rt,zt,r,z,anglr
	integer j

	if (-(zt-z).ge.0 &
     & .and. &
     & -(rt-r).ge.0) then
	anglr=atan((zt-z) &
     & /(rt-r))	
	endif
	if (-(zt-z).ge.0 &
     & .and.  &
     & -(rt-r).lt.0) then
	anglr=GPI+atan((zt-z)  &
     & /(rt-r))	
	endif
	if (-(zt-z).lt.0  &
     & .and.  &
     & -(rt-r).lt.0) then
	anglr=GPI+atan((zt-z)  &
     & /(rt-r))	
	endif
	if (-(zt-z).lt.0  &
     & .and.  &
     & -(rt-r).ge.0) then
	anglr=GPI2+atan((zt-z)  &
     & /(rt-r))	
	endif

	if (anglr.gt.GPI2) anglr=anglr-GPI2

return
end



			
		subroutine inverse_matrix_equilef(a,c,n)
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! c(n,n) - inverse matrix of A
! Based on Doolittle LU factorization for Ax = B
! Alex G. December 2009. www2.odu.eud/~agodunov/computing/programs/book2/Ch06/Inverse.f90
	implicit none 

	integer n

	double precision a(n,n), c(n,n)
	double precision L(n,n), U(n,n), b(n), d(n), x(n)
	double precision coeff

	integer i, j, k

!	write(*,*) 'inverse '


!	write(*,*) 'inverse ',n

	L=0.0
	U=0.0
	b=0.0

!	write(*,*) 'inverse ',n


! step 1: forward elimination
	do k=1, n-1
	   do i=k+1,n
	      coeff=a(i,k)/a(k,k)
	      L(i,k) = coeff
	      do j=k+1,n
	         a(i,j) = a(i,j)-coeff*a(k,j)
	      end do
	   end do
	end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
	do i=1,n
	  L(i,i) = 1.0
	end do
! U matrix is the upper triangular part of A
	do j=1,n
	  do i=1,j
	    U(i,j) = a(i,j)
	  end do
	end do

! Step 3: compute columns of the inverse matrix C
	do k=1,n
	  b(k)=1.0
	  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
	  do i=2,n
	    d(i)=b(i)
	    do j=1,i-1
	      d(i) = d(i) - L(i,j)*d(j)
	    end do
	  end do
! Step 3b: Solve Ux=d using the back substitution
	  x(n)=d(n)/U(n,n)
	  do i = n-1,1,-1
	    x(i) = d(i)
	    do j=n,i+1,-1
	      x(i)=x(i)-U(i,j)*x(j)
	    end do
	    x(i) = x(i)/u(i,i)
	  end do
! Step 3c: fill the solutions x(n) into column k of C
	  do i=1,n
	    c(i,k) = x(i)
	  end do
	  b(k)=0.0
	end do
!	write(*,*) 'end inverse'
	end subroutine inverse_matrix_equilef












	subroutine interp_j_fromrhotorz

	use ef_circuit
	
	implicit none
	
	integer i,j,k,k1,k2
	double precision 	t1,t2,t3,t4
	
!	write(*,*) jrhoteta(1:91,1:90)
!	stop

! go from jrhoteta to jrz
	jrz=0.
	jrhoteta(1:nrho,nteta+1)=jrhoteta(1:nrho,1)
	do j=1,nz2
	do i=1,nr2
		call curinterp_ef(r(i),z(j),jrhoteta(1:nrho,1:nteta+1), & 
		& rho(1:nrho,1:nteta+1),teta(1:nteta+1),raxp,zaxp,nrho,nteta+1,jrz(i,j))
	enddo
	enddo

	return
	end


	subroutine t_find_u_n(i1,j1,i2,j2,t1)		
	
	use ef_circuit, only: u_n
	
	implicit none
	integer, intent(in):: i1,i2,j1,j2
	double precision,intent(out):: t1
		t1=(1.-u_n(i1,j1))/(u_n(i2,j2)-u_n(i1,j1))
	return
	end



	subroutine curinterp_ef(r,z,jrho,rho,teta,rax,zax,nrho,nteta,j)

	use pi_grec_vars
	use numerical_tools, only: linterp
	
	implicit none
	integer, intent(in) :: nrho,nteta
	double precision, intent(in) :: r,z,jrho(nrho,nteta),rho(nrho,nteta),teta(nteta), &
		rax,zax
	double precision,intent(out) :: j

	integer i,k,k1,k2,k3,k4,j1,j2
	double precision anglr,rho0,x1(2),y1(2),jt(2),dum1(1)
	double precision z1,z2,z3,z4,a1,a2,a3,a4

	double precision r1,r2,r3,r4,jj1(4),matrix(4,4),imatrix(4,4)
	double precision coef(4),d1,d2,d3,d4

			
	call find_angle_ef(rax,zax,r,z,anglr)
  rho0=sqrt((r-rax)**2.0+(z-zax)**2.0)
  
	if (anglr.lt.teta(1)) anglr=anglr+GPI2
	
	j1=1
	do i=1,nteta
		if (anglr.ge.teta(i)) j1=i
	enddo
		j2=j1+1

	z1=rho(nrho,j1)
	z2=rho(nrho,j2)
	if (rho0.gt.z1.or.rho0.gt.z2) then
		j=0.
		return
	endif

k1=1
	do i=1,nrho-1
		if (rho0.ge.rho(i,j1)) k1=i
	enddo	
		k2=k1+1

k3=1
	do i=1,nrho-1
		if (rho0.ge.rho(i,j2)) k3=i
	enddo	
		k4=k3+1
	
r1=rax+rho(k1,j1)*cos(teta(j1))		
r2=rax+rho(k2,j1)*cos(teta(j1))			
r3=rax+rho(k4,j2)*cos(teta(j2))
r4=rax+rho(k3,j2)*cos(teta(j2))	
z1=zax+rho(k1,j1)*sin(teta(j1))		
z2=zax+rho(k2,j1)*sin(teta(j1))			
z3=zax+rho(k4,j2)*sin(teta(j2))
z4=zax+rho(k3,j2)*sin(teta(j2))	
jj1(1)=jrho(k1,j1)
jj1(2)=jrho(k2,j1)
jj1(3)=jrho(k4,j2)
jj1(4)=jrho(k3,j2)
d1=sqrt((r1-r)**2.+(z1-z)**2.)
d2=sqrt((r2-r)**2.+(z2-z)**2.)
d3=sqrt((r3-r)**2.+(z3-z)**2.)
d4=sqrt((r4-r)**2.+(z4-z)**2.)

j=(jj1(1)*d2*d3*d4+jj1(2)*d1*d3*d4+jj1(3)*d1*d2*d4+jj1(4)*d1*d2*d3) &
	/(d1*d2*d3+d1*d3*d4+d2*d3*d4+d1*d2*d4)
	

	
	
	
	
	return
	end



























	subroutine discrete_sine_transform_ef(n,y)
	
	use pi_grec_vars
	use fft_mod_eff

	implicit none
	
	integer n,i,j,imethod1,icall,k
	double precision, intent(inout) :: y(n)
!	complex :: z(2*n)	
	double precision z(n)
	complex(kind=dp) :: d(2*(n+1))

	imethod1=1

	if (imethod1.eq.1) then


!	write(*,*) 'pre',y(1:n)
	z=0.
	do i=1,n
		do j=1,n
			z(i)=z(i)+y(j)*sintable(i,j)
		enddo
	enddo
	y=z

!	z=0.
!	do i=1,n
!		do j=1,n
!			z(i)=z(i)+y(j)*sintable(i,j)
!		enddo
!	enddo
!	y=z
	!
!
!	write(*,*) 'post',y(1:n)

!	pause

	return
	
	endif




! fast sine transform 
! this problem is equivalent to DST-I with N = n+1 

	if (imethod1.eq.2) then

	k=2*(n+1)

	z=y
!	write(*,*) 'pre',y(1:n),k
	d(1)=cmplx(0.,0.)
	do i=1,n
		d(i+1)=cmplx(z(i),0)
		d(k-i+1)=cmplx(-z(i),0)
	enddo
	d(k-n)=cmplx(0.,0.)

	call fft_eff(d)
	
!		write(111,*) ' '  	
!	do i=1,2*k
!		write(111,'(3E25.11)') i+0.,real(d(i)),aimag(d(i)) 	
!	enddo

	d(1:k-1)=d(2:k)
	do i=1,n
		z(i)=0.5*aimag(d(i)-d(k-i))	
	enddo

!	d(1)=cmplx(0.,0.)
!	do i=1,n
!		d(i+1)=cmplx(z(i),0)
!		d(k-i+1)=cmplx(-z(i),0)
!	enddo
!	d(k-n)=cmplx(0.,0.)!

!	call fft_eff(d)

!	d(1:k-1)=d(2:k)
!	do i=1,n
!		z(i)=0.5*aimag(d(i)-d(k-i))	
!	enddo

	y=z/2.
!	write(*,*) 'post',y(1:n)
!	pause
	endif


!	pause

	return
	end
	

















	subroutine coil_forces_feqis(ncoilz,force_R,force_Z,plasma_state) !

	use ef_circuit
	use green_matrix
	
	implicit none
	integer i,j,k,ncoilz
	double precision force_R(ncoilz),force_Z(ncoilz)
	double precision x1
	integer nblock_a
	integer plasma_state
	
	force_R=0.
	force_Z=0.
	nblock_a=nblocks-npassive
	
	write(*,*) 'coil forces',nblock_a
	
	if (plasma_state.eq.1) then !not sure about the plasma response...
		do i=1,nblock_a
			x1=sum(jrz(1:nr2,1:nz2)* &
      & area_eff(1:nr2,1:nz2)*dgreenirpl(1:nr2,1:nz2,i))
			force_R(i)=force_R(i)+curconduc(mequivalence(i))*x1
			x1=-sum(jrz(1:nr2,1:nz2)* &
      & area_eff(1:nr2,1:nz2)*dgreenizpl(1:nr2,1:nz2,i))
			force_Z(i)=force_Z(i)+curconduc(mequivalence(i))*x1
		enddo
	endif
	
!block-to-block
	do i=1,nblock_a
		do j=1,nblock_a
				if (i.ne.j) then
					force_R(i)=force_R(i)+ & 
					& curconduc(mequivalence(j))*curconduc(mequivalence(i))*dgreenirj(i,j)
					force_Z(i)=force_Z(i)- & 
					& curconduc(mequivalence(j))*curconduc(mequivalence(i))*dgreenizj(i,j)
				endif
		enddo		
	enddo
	

	
	
	return
	end






	subroutine plasma_psi_to_coils_ef !

	use ef_circuit
	use green_matrix
	
	implicit none
	integer i,j,k
	
	do i=1,nconduc
		psiplasmatoconduc(i)=sum(jrz(1:nr2,1:nz2)* &
      & area_eff(1:nr2,1:nz2)*greeni(1:nr2,1:nz2,i))
	enddo


	return
	end














      subroutine get_zccurb_efff(rc_cur,zc_cur,z2c_cur, &
     & rgeoc,zgeoc,ahorc,bpcell)

	use   ef_circuit  
	use   metric_coefficients_pbe  

	implicit none
	real*8 rc_cur,zc_cur,z2c_cur
	real*8 perimz,ahorc2,ahorc,rgeoc,zgeoc
	real*8 bpcell(i_dim2,i_dim2),avgelem
	integer i,j
	real*8 	sum_z2c,sum_zc,sum_rc,sum_Ipla
	
         sum_z2c=0.d0
         sum_zc=0.d0
         sum_rc=0.d0
         sum_Ipla=0.d0
			  	perimz=0.         
	rgeoc=0.
	zgeoc=0.
	ahorc=0.
	ahorc2=0.
		
         zc_cur=Z_curr_0D
         z2c_cur=zc_cur
         rc_cur=R_curr_0D

!         zc_cur=Z_curr_2D
!         z2c_cur=zc_cur
!         rc_cur=R_curr_2D


!	open(32,file='fort.344')
!		write(32,*) rpol(1:nrho-1,1:nteta),zpol(1:nrho-1,1:nteta), &
!     & jrhoteta(1:nrho-1,1:nteta),dator(1:nrho-1,1:nteta)
!	close(32)

	do i=1,nrho-1
          do j=1,nteta
 avgelem=dator(i,j) !/bpcell(i,j-1)
	perimz=perimz+avgelem
	rgeoc=rgeoc+rpol(i,j)*avgelem
	zgeoc=zgeoc+zpol(i,j)*avgelem
          enddo
	enddo
	rgeoc=rgeoc/perimz
	zgeoc=zgeoc/perimz

	do i=1,nrho-1
          do j=1,nteta
	avgelem=dator(i,j) !/bpcell(i,j-1)
	ahorc2=ahorc2+(rpol(i,j)-rgeoc)**2.*avgelem
          enddo
	enddo
	ahorc2=ahorc2/perimz
	ahorc=2.*sqrt(ahorc2)


	return
	end








	subroutine psib_ext_efff(psiext_out)  !gives back external flux on plasma boundary

	use ef_circuit
	use green_matrix
	
	implicit none
	integer i,j,k,m
	integer i1,j1,k1,m1
	double precision psiext_out,dlt,dllt,dum1,dum2	
	
!cycle over boundary
	psiext_out=0.
	dllt=0.
	do i=1,nbnd-1
		call find_fields_interp_ef_psiext(rbnd(i),zbnd(i),dum1) !give back psi,br,bz at r0,z0
		call find_fields_interp_ef_psiext(rbnd(i+1),zbnd(i+1),dum2) !give back psi,br,bz at r0,z0
		dlt=sqrt((rbnd(i+1)-rbnd(i))**2.+(zbnd(i+1)-zbnd(i))**2.)
		psiext_out=psiext_out+0.5*(dum1+dum2)*dlt
		dllt=dllt+dlt
	enddo	
		call find_fields_interp_ef_psiext(rbnd(nbnd),zbnd(nbnd),dum1) !give back psi,br,bz at r0,z0
		call find_fields_interp_ef_psiext(rbnd(1),zbnd(1),dum2) !give back psi,br,bz at r0,z0
		dlt=sqrt((rbnd(1)-rbnd(nbnd))**2.+(zbnd(1)-zbnd(nbnd))**2.)
		psiext_out=psiext_out+0.5*(dum1+dum2)*dlt
		dllt=dllt+dlt

	psiext_out=psiext_out/dllt


	write(*,*) 'psibbb',psibnd,psiext_out


	return
	end




	subroutine psiplex_calc_ef(dumz)

	use astra2fbe
	use ef_circuit
	implicit none
	double precision dumz
	integer i,j,k
	double precision t1,t2,t3,t4,z1,z2,z3,z4
	double precision x1,x2,x3,x4,y1,y2,y3,y4
	double precision arc1,arc2
	
	write(*,*) 'spid par',psplex_from_fbe
		
	if (psplex_from_fbe.eq.1) then
	else
		dumz=dumz
		return
	endif	
!calculate 
!psplexs = psplexs/(btor*rhos(iplas)**2)*qedge	
!cycle over boundary
!calculate psplex 
!	dydr(1)=0.
!	dydr(2)=0.
!	do j=1,nt
!	do i=1,nt
!		call green_function(X(nx-1,i),Y(nx-1,i),
!     & X_i1(nx-1,j),Y_i1(nx-1,j),greenf)
!		dydr(1)=dydr(1)+greenf/Rmaj2(nx-1,i)*
!     & gradpsia(nx-1,i)*dl_arc1(i)*dl_arc2(j)	
!	enddo
!	enddo
!	psplex=dydr(1)/sum(dl_arc2)
!	psplex=psplex/(GP2*btor*rhoedge**2)*qedge
!	psiext_out=0.
!	dllt=0.
	z4=0.
	z2=0.
	do j=1,nbnd-1
		x1=rbnd(j)
		y1=zbnd(j)
		x3=rbnd(j+1)
		y3=zbnd(j+1)
		x2=0.5*(rbnd(j+1)+rbnd(j))
		y2=0.5*(zbnd(j+1)+zbnd(j))
		arc2=sqrt((x3-x1)**2.+(y3-y1)**2.)
	do i=1,nbnd-1
		x1=rbnd(i)
		y1=zbnd(i)
		x3=rbnd(i+1)
		y3=zbnd(i+1)
		arc1=sqrt((x3-x1)**2.+(y3-y1)**2.)
		call green_function(x1,y1,x2,y2,z3)
		call find_fields_interp_ef_psionly(rbnd(i)+dr/2,zbnd(i),t1)
		call find_fields_interp_ef_psionly(rbnd(i),zbnd(i)+dz/2,t2)
		call find_fields_interp_ef_psionly(rbnd(i)-dr/2,zbnd(i),t3)
		call find_fields_interp_ef_psionly(rbnd(i),zbnd(i)-dz/2,t4)
		z1=sqrt(((t3-t1)/dr)**2.+((t4-t2)/dz)**2.)
		z2=z2+z3/x1*z1*arc1*arc2
	enddo
		z4=z4+arc2
	enddo	
!	psiext_out=psiext_out/dllt

!	psplex=psplex
	dumz=z2/z4*GPI2

	return
	end





	subroutine wrd_equilef  !gives back external flux on plasma boundary

	use ef_circuit
	
	implicit none

!		open(32,file='fort.4444')
!	write(32,*) r(1:nr2),z(1:nz2),jrz(1:nr2,1:nz2)
!	close(32)
!	open(32,file='fort.4445')
!	write(32,*) psiextrz(1:nr2,1:nz2)
!	close(32)
!	open(32,file='fort.4447')
!		write(32,*) psirz(1:nr2,1:nz2)
!	close(32)

	return
	end


	subroutine wrd_equilef_pbe  !gives back external flux on plasma boundary

	use ef_circuit
	
	implicit none

!		open(32,file='fort.4444')
!	write(32,*) rpol(1:nrho,1:nteta),zpol(1:nrho,1:nteta),jrhoteta(1:nrho,1:nteta)
!	close(32)
!	open(32,file='fort.4445')
!	write(32,*) rpul(1:nrho,1:nteta),zpul(1:nrho,1:nteta),psirhoteta(1:nrho,1:nteta)
!	close(32)
!	open(32,file='fort.4447')
!		write(32,*) psia_1d(1:nrho),ffp_1d(1:nrho),ppp_1d(1:nrho)
!	close(32)

	return
	end



!obsolete
	subroutine get_zccur_efff(Rcurr,Zcurr,Zsquad)
	implicit none
	real*8 rcurr,zcurr,zsquad
	return
	end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine find_demo_gaps_efff(ngaps,demo_gaps,geom1d)

	use ef_circuit
	
	implicit none
	integer ngaps,i,j
	integer n_iterz,j1,j2,j3,j4
	double precision demo_gaps(ngaps,4),geom1d(ngaps)
	double precision dumx,dumy,dumz,l_ref,l_gap
	double precision dumx2,dumy2,dumu1,dumu2,u00,up
	double precision dumx3,dumy3,l_gap2,l_gap3,x003,y003,u002
	double precision dumx0,dumy0,dumxx,dumyy,d_step
	double precision gapmin,gapmax,x00,y00,x002,y002
	double precision tolez,bolez,dur1,dur2
	integer onlypos
	tolez=err_gaptolez

	d_step=0.1 !advance in 1 cm steps	
	n_iterz=100
	gapmin=-1.5
	gapmax=2.5
	up=psibnd
	
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
	call find_fields_interp_ef_psionly(dumx,dumy,u00) 
	call find_fields_interp_ef_psionly(dumx2,dumy2,u002) 

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
	call find_fields_interp_ef_psionly(dumx,dumy,u00) 
	call find_fields_interp_ef_psionly(dumx2,dumy2,u002) 
		

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
	if (abs(l_gap).lt.abs(l_gap2)) l_gap3=l_gap
	
	if (isnan(l_gap)) l_gap3=l_gap2
	if (isnan(l_gap2)) l_gap3=l_gap

	write(*,*) 'gaps ',i,dumx0,dumy0,dumx,dumy, &
     & dumx2,dumy2,dumxx,dumyy,l_gap2,l_gap,l_gap3

	geom1d(i)=min(gapmax,max(gapmin,l_gap3))
!	write(999,'(3E25.11)') l_gap,l_gap2,geom1d(i)
	else
		geom1d(i)=min(gapmax,max(gapmin,l_gap))
	endif
	enddo
	
	return
	end



	subroutine find_demo_meas2021_equilef(nbexp,nfexp,br_in,br_out, &
     & flux_in,flux_out,ncoilzzz,pjk,dddc)

	use ef_circuit
	implicit none
	integer nbexp,nfexp,i,j,k,i0,j0,ncoilzzz
	double precision br_in(nbexp,3),br_out(nbexp),pjk(ncoilzzz)
	double precision flux_in(nfexp,2),flux_out(nfexp)
	double precision br1,bz1,brr,brz,bzr,bzz,dum0,dddc(200)

	pjk(1:nconduc)=curconduc(1:nconduc)
	dddc(1:nconduc)=dpc(1:nconduc)

	do i=1,nbexp
		call find_fields_interp_ef(br_in(i,1),br_in(i,2),dum0,br1,bz1,brr,brz,bzr,bzz)
		br_out(i)=br1*cos(br_in(i,3))+bz1*sin(br_in(i,3))	
	enddo

	do i=1,nfexp
		call find_fields_interp_ef_psionly(flux_in(i,1),flux_in(i,2),flux_out(i)) 
	enddo

	return
	end




	subroutine psi_external_calc_ef

	use ef_circuit
	use green_matrix
	
	implicit none
	integer i,j,k
	
	do j=1,nz2	
	do i=1,nr2	
	psiextrz(i,j)=sum(curconduc(1:nconduc)*greeni(i,j,1:nconduc))
	enddo
	enddo


	return
	end






	
	!generate zlimpotential
	subroutine generate_zlim_potential(r,z,nl,rl, & 
	& zl,zlimp)

	use	astra2fbe
	use pi_grec_vars
	use fenix_params
	
	implicit none
	integer nr,nz,nlim,i,j,k,nt,nl
	double precision r,z,rl(nl),zl(nl),zlimp
	double precision anglr(nl)
	double precision discrim,anglf(nl),anglt(nl)	
	double precision diffa(nl),rt,zt	


	zlimp=-1.

	rt=r
	zt=z	

!check if sol current filament is inside plasma
	do j=1,nl

	call find_angle_ef(rt,zt,rl(j),zl(j),anglr(j))
	
	enddo

	nt=nl
	do j=1,nt
	if (j.gt.1)	then
			diffa(j)=abs(anglr(j)-anglr(j-1))
			if (diffa(j).gt.1.) diffa(j)=abs(diffa(j)-GPI2)
		endif
	if (j.eq.1)	then
	diffa(j)=abs(anglr(j)-anglr(nt))
			if (diffa(j).gt.1.) diffa(j)=abs(diffa(j)-GPI2)
	endif
	enddo

	anglt(1:nt)=sin(anglr(1:nt))*diffa(1:nt)
	anglf(1:nt)=cos(anglr(1:nt))*diffa(1:nt)	

	discrim=abs(sum(anglt(1:nt)))/sum(diffa(1:nt))+ &
     & 	abs(sum(anglf(1:nt)))/sum(diffa(1:nt))


	if (discrim.lt.0.6) zlimp=1.

	if (use_zlim_pot.eq.0) zlimp=1.


	return
	end



	subroutine compute_fields_fbe_ef(iaold,jaold,br,bz,brr,brz,bzr,bzz)

	use ef_circuit

	implicit none
	integer iaold,jaold
	double precision br,bz,brr,brz,bzr,bzz
	double precision br2,bz2,br3,bz3
	
		br=1./r(iaold)*(psirz(iaold,jaold+1)-psirz(iaold,jaold))/dz	
		bz=-2./(r(iaold)+r(iaold+1))*(psirz(iaold+1,jaold)-psirz(iaold,jaold))/dr	

		br2=1./r(iaold)*(psirz(iaold,jaold)-psirz(iaold,jaold-1))/dz	
		bz2=-2./(r(iaold)+r(iaold-1))*(psirz(iaold,jaold)-psirz(iaold-1,jaold))/dr	
		br3=1./r(iaold+1)*(psirz(iaold+1,jaold+1)-psirz(iaold+1,jaold))/dz	
		bz3=-2./(r(iaold)+r(iaold+1))*(psirz(iaold+1,jaold+1)-psirz(iaold,jaold+1))/dr	

		brr=1/dr*(br3-br)
		brz=1/dz*(br-br2)
		bzr=1/dr*(bz-bz2)
		bzz=1/dz*(bz3-bz)
	
	return
	end





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE least_square_biquad_ef(r,z,u,n,c,rax,zax,uax,derivs)

	implicit none
	integer n,k
	double precision r(n),z(n),u(n),derivs(8)
	double precision A(6,6),B(6),cc(6),c(6),Ainv(6,6)
	double precision  rax,zax,uax
	double precision sums(21),det,det_r,det_z

	sums=0.

	sums(1)=sum(r**4.)
	sums(2)=sum(z**4.)
	sums(3)=sum(r**2.*z**2.)
	sums(4)=sum(r**2.)
	sums(5)=sum(z**2.)
	sums(6)=n+0.
	sums(7)=sum(r**3.*z)
	sums(8)=sum(r**3.)
	sums(9)=sum(r**2.*z)
	sums(10)=sum(r*z**3.)
	sums(11)=sum(r*z**2.)
	sums(12)=sum(z**3.)
	sums(13)=sum(r*z)
	sums(14)=sum(r)
	sums(15)=sum(z)
	sums(16)=sum(u*r**2.)
	sums(17)=sum(u*z**2.)
	sums(18)=sum(u*r*z)
	sums(19)=sum(u*r)
	sums(20)=sum(u*z)
	sums(21)=sum(u)

	B(1)=-2*sums(16)
	B(2)=-2*sums(17)
	B(3)=-2*sums(18)
	B(4)=-2*sums(19)
	B(5)=-2*sums(20)
	B(6)=-2*sums(21)

	A(1,1)=4*sums(1)
	A(1,2)=4*sums(3)
	A(1,3)=4*sums(7)
	A(1,4)=4*sums(8)
	A(1,5)=4*sums(9)
	A(1,6)=4*sums(4)

	A(2,1)=A(1,2)
	A(2,2)=4*sums(2)
	A(2,3)=4*sums(10)
	A(2,4)=4*sums(11)
	A(2,5)=4*sums(12)
	A(2,6)=4*sums(5)

	A(3,1)=A(1,3)
	A(3,2)=A(2,3)
	A(3,3)=4*sums(3)
	A(3,4)=4*sums(9)
	A(3,5)=4*sums(11)
	A(3,6)=4*sums(13)

	A(4,1)=A(1,4)
	A(4,2)=A(2,4)
	A(4,3)=A(3,4)
	A(4,4)=4*sums(4)
	A(4,5)=4*sums(13)
	A(4,6)=4*sums(14)

	A(5,1)=A(1,5)
	A(5,2)=A(2,5)
	A(5,3)=A(3,5)
	A(5,4)=A(4,5)
	A(5,5)=4*sums(5)
	A(5,6)=4*sums(15)

	A(6,1)=A(1,6)
	A(6,2)=A(2,6)
	A(6,3)=A(3,6)
	A(6,4)=A(4,6)
	A(6,5)=A(5,6)
	A(6,6)=4*sums(6)


!find coefficients
	call inverse_matrix_equilef(A,Ainv,6)

		do k=1,6
			cc(k)=-2*sum(Ainv(k,1:6)*B(1:6))
		enddo
			c(4)=cc(1)
			c(5)=cc(2)
			c(6)=cc(3)
			c(2)=cc(4)
			c(3)=cc(5)
			c(1)=cc(6)


!c   magnetic axis

         det=4.d0*cc(1)*cc(2)-cc(3)**2.

         det_r=-2.d0*cc(2)*cc(4)+cc(3)*cc(5)
         det_z=-2.d0*cc(1)*cc(5)+cc(3)*cc(4)

         rax=det_r/det
         zax=det_z/det

	uax=cc(1)*rax**2.+cc(2)*zax**2.+cc(3)*rax*zax+cc(4)*rax+cc(5)*zax+cc(6)


         derivs(1)=cc(4)
         derivs(2)=cc(5)
         derivs(3)=2.*cc(1)
         derivs(4)=2.*cc(2)
         derivs(5)=cc(3) 
				 
!	write(88881,*) ' cc',c(1:6),B(1:6),det,det_r,det_z,rax,zax,uax

         RETURN
         END
				 
				 
				 
				 
				 
				 
	SUBROUTINE exact_biquad_ef(r,z,u,n,ccc,rax,zax,uax,derivs,dr,dz)
	use errors_params
	implicit none
	integer n,k,i,j
	double precision x(9),y(9),r(9),z(9),u(9),derivs(8)
	double precision A(9,9),B(9),ccc(6),Ainv(9,9)
	double precision  rax,zax,uax,c(9),dr,dz
	double precision s_r,s_z,s_r2,s_z2,s_rz
	double precision s_r3,s_rz2,s_r2z,s_z3
	double precision s_r4,s_r2z2,s_r3z,s_z4,s_rz3
	double precision s_u,s_ur,s_uz,s_ur2,s_uz2,s_urz	
	double precision det, det_r,det_z,tolez

	integer niter,j_success
	
	tolez=err_find_biquad
	
!transformation
	x=(r-r(5))/dr		
	y=(z-z(5))/dz		
	

!	write(*,*) r,z,u


!find coefficients
	call ainv_matrix_def(Ainv)

		do k=1,9
			c(k)=sum(Ainv(k,1:9)*u(1:9))
		enddo

	rax=x(5)
	zax=y(5)

!	write(667,*) Ainv
!	write(667,*) c

!now find axis
	niter=0
	s_r=100000.
	s_z=100000.
1234	continue
	niter=niter+1
s_r2=2*C(1)*rax*zax**2 + 2*C(2)*rax*zax+C(3)*zax**2+C(4)*zax+2*C(5)*rax+C(7);
s_z2=2*C(1)*rax**2*zax**1 + C(2)*rax**2 +2*C(3)*zax**1*rax+C(4)*rax+2*C(6)*zax+C(8);
A(1,1)=2*C(1)*zax**2+2*C(2)*zax+2*C(5);
A(1,2)=4*C(1)*rax**1*zax+2*C(2)*rax+2*C(3)*zax+C(4);
A(2,2)=2*C(1)*rax**2+2*C(3)*rax+2*c(6);
A(2,1)=4*C(1)*rax**1*zax+2*C(2)*rax+2*C(3)*zax+C(4);
B(1)=s_r2;
B(2)=s_z2;
det=(A(1,1)*A(2,2))-(A(1,2)*A(2,1))
s_r2=1/det*(A(2,2)*B(1)-A(1,2)*B(2))
s_z2=1/det*(A(1,1)*B(2)-A(2,1)*B(1))
rax=rax-s_r2
zax=zax-s_z2

s_r2=s_r
s_z2=s_z
s_r=2*C(1)*rax*zax**2 + 2*C(2)*rax*zax+C(3)*zax**2+C(4)*zax+2*C(5)*rax+C(7)
s_z=2*C(1)*rax**2*zax**1 + C(2)*rax**2 +2*C(3)*zax**1*rax+C(4)*rax+2*C(6)*zax+C(8)
!	write(5671,*) niter,rax,zax,abs(s_r2)+abs(s_z2),abs(s_r)+abs(s_z),r(5),z(5)

	if (abs(s_r).lt.tolez.and.abs(s_z).lt.tolez) then
			j_success=1
		goto 1235
	endif
	if (abs(rax).gt.1.) then
		j_success=0.
		goto 1235
	endif
	if (abs(zax).gt.1) then
		j_success=0.
		goto 1235
	endif
	if (niter.gt.100000) then
		j_success=0.
		goto 1235
	endif
	
	goto 1234
	
1235	continue

!pause
	if (j_success.eq.0) then
		rax=1.e6
		zax=1.e6
		uax=-1.e6
		derivs=1.e6
		return	
	endif

!c   magnetic axis

         uax=c(1)*rax**2.*zax**2. + & 
				  & c(2)*rax**2.*zax**1. + &
				  & c(3)*rax**1.*zax**2. + &
				  & c(4)*rax**1.*zax**1. + &
				  & c(5)*rax**2.*zax**0. + &
				  & c(6)*rax**0.*zax**2. + &
				  & c(7)*rax**1.*zax**0. + &
				  & c(8)*rax**0.*zax**1. + &
				  & c(9)*rax**0.*zax**0. 
					

	derivs(1)=1/dr*(2*C(1)*rax*zax**2 + 2*C(2)*rax*zax+C(3)*zax**2+C(4)*zax+2*C(5)*rax+C(7))
	derivs(2)=1/dz*(2*C(1)*rax**2*zax**1 + C(2)*rax**2 +2*C(3)*zax**1*rax+C(4)*rax+2*C(6)*zax+C(8))
	derivs(3)=1/dr**2.*(2.*c(1)*zax**2.+2*c(2)*zax+2*c(5))
	derivs(4)=1/dz**2.*(2.*c(1)*rax**2.+2*c(3)*rax+2*c(6))
	derivs(5)=1/dr/dz*(4*c(1)*rax*zax+2*c(2)*rax+2*c(3)*zax+c(4)) 
				 
	rax=rax*dr+r(5)
	zax=zax*dz+z(5)


         RETURN
         END
				 
				 
				 
				 
				 
				 
				 
	SUBROUTINE exact_biquad_regress_ef(r,z,u,n,ccc,rax,zax,uax,derivs,dr,dz,rx,zx)

	implicit none
	integer n,k,i,j
	double precision x(9),y(9),r(9),z(9),u(9),derivs(8)
	double precision A(9,9),B(9),ccc(6),Ainv(9,9)
	double precision  rax,zax,uax,c(9),dr,dz,rx,zx
	double precision s_r,s_z,s_r2,s_z2,s_rz,xx,yy
	double precision s_r3,s_rz2,s_r2z,s_z3
	double precision s_r4,s_r2z2,s_r3z,s_z4,s_rz3
	double precision s_u,s_ur,s_uz,s_ur2,s_uz2,s_urz	
	double precision det, det_r,det_z,tolez

	integer niter,j_success
	
!transformation
	x=(r-r(5))/dr		
	y=(z-z(5))/dz		
	xx=(rx-r(5))/dr
	yy=(zx-z(5))/dz	
	

!find coefficients
	call ainv_matrix_def(Ainv)

		do k=1,9
			c(k)=sum(Ainv(k,1:9)*u(1:9))
		enddo

	rax=xx
	zax=yy

!	write(667,*) Ainv
!	write(667,*) c

!now find axis
s_r2=2*C(1)*rax*zax**2 + 2*C(2)*rax*zax+C(3)*zax**2+C(4)*zax+2*C(5)*rax+C(7);
s_z2=2*C(1)*rax**2*zax**1 + C(2)*rax**2 +2*C(3)*zax**1*rax+C(4)*rax+2*C(6)*zax+C(8);
A(1,1)=2*C(1)*zax**2+2*C(2)*zax+2*C(5);
A(1,2)=4*C(1)*rax**1*zax+2*C(2)*rax+2*C(3)*zax+C(4);
A(2,2)=2*C(1)*rax**2+2*C(3)*rax+2*c(6);
A(2,1)=4*C(1)*rax**1*zax+2*C(2)*rax+2*C(3)*zax+C(4);
B(1)=s_r2;
B(2)=s_z2;
det=(A(1,1)*A(2,2))-(A(1,2)*A(2,1))
s_r2=1/det*(A(2,2)*B(1)-A(1,2)*B(2))
s_z2=1/det*(A(1,1)*B(2)-A(2,1)*B(1))
rax=rax-s_r2
zax=zax-s_z2

s_r=2*C(1)*rax*zax**2 + 2*C(2)*rax*zax+C(3)*zax**2+C(4)*zax+2*C(5)*rax+C(7)
s_z=2*C(1)*rax**2*zax**1 + C(2)*rax**2 +2*C(3)*zax**1*rax+C(4)*rax+2*C(6)*zax+C(8)

         uax=c(1)*rax**2.*zax**2. + & 
				  & c(2)*rax**2.*zax**1. + &
				  & c(3)*rax**1.*zax**2. + &
				  & c(4)*rax**1.*zax**1. + &
				  & c(5)*rax**2.*zax**0. + &
				  & c(6)*rax**0.*zax**2. + &
				  & c(7)*rax**1.*zax**0. + &
				  & c(8)*rax**0.*zax**1. + &
				  & c(9)*rax**0.*zax**0. 
					

	derivs(1)=1/dr*(2*C(1)*rax*zax**2 + 2*C(2)*rax*zax+C(3)*zax**2+C(4)*zax+2*C(5)*rax+C(7))
	derivs(2)=1/dz*(2*C(1)*rax**2*zax**1 + C(2)*rax**2 +2*C(3)*zax**1*rax+C(4)*rax+2*C(6)*zax+C(8))
	derivs(3)=1/dr**2.*(2.*c(1)*zax**2.+2*c(2)*zax+2*c(5))
	derivs(4)=1/dz**2.*(2.*c(1)*rax**2.+2*c(3)*rax+2*c(6))
	derivs(5)=1/dr/dz*(4*c(1)*rax*zax+2*c(2)*rax+2*c(3)*zax+c(4)) 
				 
	rax=rax*dr+r(5)
	zax=zax*dz+z(5)
!	write(5671,*) rax,zax,abs(s_r)+abs(s_z),s_r2,s_z2,r(5),z(5)


         RETURN
         END
				 
				 
				 
				 
	real*8 function frlim_ef(dp,ylim,rx,zx,rm,zm)



	implicit none

	double precision dp(5),dxx,dxy,dyy,ylim,rx,zx,rm,zm
	double precision disc,cc,cdpls,cdmns,ang1,ang2,c1,dl2x,c2,dl2y

             Dxx=dp(3)
             Dxy=dp(4)
             Dyy=dp(5)

             disc=(Dxy/Dyy)**2 - Dxx/Dyy

          if(Disc.lt.0.) then

          !  write(6,*) ' FRLIM: ***ERROR disc lt.0'
          !  write(6,*) ' rx zx ',rx,zx

             cc=(zm-zx)/(rm-rx)

             cc=-1.d0/cc
             go to 100

          endif

             cdpls = -Dxy/Dyy + dsqrt(disc)
             cdmns = -Dxy/Dyy - dsqrt(disc)

             ang1 = 0.5d0*(datan(cdpls)+datan(cdmns))
             ang2 =-0.5d0*(datan(1.d0/cdpls)+datan(1.d0/cdmns))


!c.........calculation D2u/Dl2(direction ang1    )

             c1=dtan(ang1)
             Dl2x=Dyy*c1*c1+2.d0*Dxy*c1+Dxx

!c.........calculation D2u/Dl2(direction ang2    )

             c2=dtan(ang2)
             Dl2y=Dyy*c2*c2+2.*Dxy*c2+Dxx

          if(Dl2x.lt.0.) then

             cc=c1

          elseif(Dl2y.lt.0.) then

             cc=c2

          else

          !  write(6,*) ' FRLIM: ***ERROR both deriv. ge.0'
          !  write(6,*) ' rx zx ',rx,zx
          !  call wrd

             cc=(zm-zx)/(rm-rx)

             cc=-1.d0/cc

          endif

 100      continue

          frlim_ef=rx+(ylim-zx)/cc

          return
          end

				 
				 
				 
				 
	subroutine find_x_point_type(dp,rx,zx,rm,zm,k)

	implicit none
	
	integer k
	double precision dp(5),dxx,dxy,dyy,ylim,rx,zx,rm,zm,frlim_ef
	double precision disc,cc,cdpls,cdmns,ang1,ang2,c1,dl2x,c2,dl2y

             Dxx=dp(3)
             Dxy=dp(4)
             Dyy=dp(5)

             disc=(Dxy/Dyy)**2. - Dxx/Dyy

             cdpls = -Dxy/Dyy + dsqrt(disc)
             cdmns = -Dxy/Dyy - dsqrt(disc)

             ang1 = 0.5d0*(datan(cdpls)+datan(cdmns))
             ang2 =-0.5d0*(datan(1.d0/cdpls)+datan(1.d0/cdmns))


!c.........calculation D2u/Dl2(direction ang1    )

             c1=dtan(ang1)
             Dl2x=Dyy*c1*c1+2.d0*Dxy*c1+Dxx

!c.........calculation D2u/Dl2(direction ang2    )

             c2=dtan(ang2)
             Dl2y=Dyy*c2*c2+2.*Dxy*c2+Dxx

          if(Dl2x.lt.0.) then

             cc=c1

          elseif(Dl2y.lt.0.) then

             cc=c2

          else

          !  write(6,*) ' FRLIM: ***ERROR both deriv. ge.0'
          !  write(6,*) ' rx zx ',rx,zx
          !  call wrd

             cc=(zm-zx)/(rm-rx)

             cc=-1.d0/cc

          endif

 100      continue

          frlim_ef=rx+(ylim-zx)/cc

          return
          end

				 
				 
				 
				 
				 
	subroutine vector_gap_find(r0,z0,angle,l_gap,dr, & 
	& gapmax,rmin,rmax,zmin,zmax,psib)
! this routine is valid only from the interior of the plasma	
	use errors_params	
	implicit none
	integer i,j,k
	double precision r0,z0,angle
	double precision tolez,d_step,psib
	double precision up,l_ref,dr,rmin,rmax,zmin,zmax
	double precision dumx0,dumy0,dumz
	double precision dumx,dumy,dumx2,dumy2,dur1,dur2
	double precision u00,u002,l_gap,gapmax
		
		tolez=err_gaptolez

	d_step=dr !advance in dr steps	
	up=psib
	
!first positive gap!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	l_ref=d_step

			dumx0=r0 ! R0
			dumy0=z0 ! Z0
			dumz=angle ! angle	 
	dur1=0.
	dur2=l_ref			
		
1781	continue

	dumx=dumx0+dur1*cos(dumz)
	dumy=dumy0+dur1*sin(dumz)
	dumx2=dumx0+dur2*cos(dumz)
	dumy2=dumy0+dur2*sin(dumz)
	if (dumx2.lt.rmin) then
		dumx2=rmin
		dur2=sqrt((dumx2-r0)**2.+(dumy2-z0)**2.)
		l_ref=dur2-dur1
	endif
	if (dumx2.gt.rmax)  then
		dumx2=rmax
		dur2=sqrt((dumx2-r0)**2.+(dumy2-z0)**2.)
		l_ref=dur2-dur1
	endif
	if (dumy2.lt.zmin) then
		 dumy2=zmin
		dur2=sqrt((dumx2-r0)**2.+(dumy2-z0)**2.)
		l_ref=dur2-dur1
	endif
	if (dumy2.gt.zmax)  then
		dumy2=zmax
		dur2=sqrt((dumx2-r0)**2.+(dumy2-z0)**2.)
		l_ref=dur2-dur1
	endif

	call find_fields_interp_ef_psionly(dumx,dumy,u00) 
	call find_fields_interp_ef_psionly(dumx2,dumy2,u002) 

	if (abs(l_ref).lt.tolez) goto 131
		
	if (u00.eq.up) goto 118
	if (u002.eq.up) goto 117
	if (u002.gt.up.and.u00.lt.up) goto 115
	if (u002.lt.up.and.u00.gt.up) goto 116

	
!case with no intersection: larger 1. m
	l_gap=0.5*(dur1+dur2)
	if (l_gap.gt.gapmax)	goto 119

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
	l_gap=5000.
	goto 211
	
131	continue
	l_gap=0.5*(dur1+dur2)

211	continue
	
	return
	end









				 

	subroutine nine_point_regression(r0,z0,pos_xpoint,ddpsi,f00)
	use ef_circuit

	implicit none
	
	integer iax,jax,i,j,k,d
	double precision pos_xpoint(2)
	double precision r0,z0,f00,xub(90),bub(90),yub(90)
	double precision f0,fr0,fz0,frr0,fzz0,frz0,fr2z0,frz20,fr2z20
	double precision f(9),drad,dzad,t1,t2,det,matrix(2,2)
	double precision c1,c2,c3,c4,c5,c6,c7,c8,c9,ddpsi(8),c(6)
	integer i1,i2,i3,i4

	call find_actual_index_ef(r0,z0,iax,jax)
!find true axis
	k=0
	i3=-1
	i1=-1
	i4=1
	i2=1
	if (iax.eq.2) i3=-1
	if (jax.eq.2) i1=-1
	if (iax.eq.nr1) i4=1
	if (jax.eq.nz1) i2=1
		
	d=(i4-i3+1)*(i2-i1+1)


	do j=i3,i4
	do i=i1,i2
		k=k+1
		xub(k)=r(iax+i)
		yub(k)=z(jax+j)
		bub(k)=psirz(iax+i,jax+j)
	enddo
	enddo
!	call least_square_biquad_ef(xub(1:d),yub(1:d),bub(1:d),d, &
!		& c,pos_xpoint(1),pos_xpoint(2),f00,ddpsi)
	call exact_biquad_ef(xub(1:d),yub(1:d),bub(1:d),d, &
		& c,pos_xpoint(1),pos_xpoint(2),f00,ddpsi,dr,dz)

	return
	end



	subroutine nine_point_regression_follow(rx,zx,pos_xpoint,ddpsi,f00)
	use ef_circuit

	implicit none
	
	integer iax,jax,i,j,k,d
	double precision pos_xpoint(2),rx,zx
	double precision r0,z0,f00,xub(90),bub(90),yub(90)
	double precision f0,fr0,fz0,frr0,fzz0,frz0,fr2z0,frz20,fr2z20
	double precision f(9),drad,dzad,t1,t2,det,matrix(2,2)
	double precision c1,c2,c3,c4,c5,c6,c7,c8,c9,ddpsi(8),c(9)
	integer i1,i2,i3,i4

	call find_fields_interp_ef_psionly(rx-dr,zx-dz,bub(1))
	call find_fields_interp_ef_psionly(rx,zx-dz,bub(2))
	call find_fields_interp_ef_psionly(rx+dr,zx-dz,bub(3))
	call find_fields_interp_ef_psionly(rx-dr,zx,bub(4))
	call find_fields_interp_ef_psionly(rx,zx,bub(5))
	call find_fields_interp_ef_psionly(rx+dr,zx,bub(6))
	call find_fields_interp_ef_psionly(rx-dr,zx+dz,bub(7))
	call find_fields_interp_ef_psionly(rx,zx+dz,bub(8))
	call find_fields_interp_ef_psionly(rx+dr,zx+dz,bub(9))

		xub(1)=rx-dr
		yub(1)=zx-dz
		xub(2)=rx
		yub(2)=zx-dz
		xub(3)=rx+dr
		yub(3)=zx-dz
		xub(4)=rx-dr
		yub(4)=zx
		xub(5)=rx
		yub(5)=zx
		xub(6)=rx+dr
		yub(6)=zx
		xub(7)=rx-dr
		yub(7)=zx+dz
		xub(8)=rx
		yub(8)=zx+dz
		xub(9)=rx+dr
		yub(9)=zx+dz

	d=9
	call exact_biquad_regress_ef(xub(1:d),yub(1:d),bub(1:d),d, &
		& c,pos_xpoint(1),pos_xpoint(2),f00,ddpsi,dr,dz,rx,zx)

	return
	end















	subroutine find_actual_index_ef(r0,z0,i,j)

	use ef_circuit

	implicit none
	integer i,j
	double precision r0,z0
	i=nint((r0-rmin)/dr+1.)	
	j=nint((z0-zmin)/dz+1.)	

	
	return
	end



	subroutine find_floor_index_ef(r0,z0,i,j)

	use ef_circuit

	implicit none
	integer i,j
	double precision r0,z0
	

	i=floor((r0-rmin)/dr+1.)	
	j=floor((z0-zmin)/dz+1.)	

	
	return
	end


	subroutine find_fields_interp_ef(r0,z0,psi0,br0,bz0,brr,brz,bzr,bzz) !give back psi,br,bz at r0,z0

	use ef_circuit

	implicit none
	integer i,j
	double precision r0,z0,psi0,br0,bz0
	double precision x1,x2,x3,x4,x5,x6,x7,x8,y1,y2,y3,y4,y5,y6
	double precision z1,z2,z3,z4,z5,z6,z7,z8,t1,t2,t3,t4,t5,t6
	double precision brr,brz,bzr,bzz
	
	i=floor((r0-rmin)/dr+1.)	
	j=floor((z0-zmin)/dz+1.)	
	x1=r(i)
	x2=r(i+1)
	y1=z(j)
	y2=z(j+1)
	z1=psirz(i,j)
	z2=psirz(i+1,j)
	z3=psirz(i,j+1)
	z4=psirz(i+1,j+1)

!bilinear interpolation
	call bilinear_average_ef(x1,x2,y1,y2,r0,z0,z1,z2,z3,z4,psi0)

!br
  t1=-1./x1*(psirz(i,j+1)-psirz(i,j-1))/dz/2.
  t2=-1./x2*(psirz(i+1,j+1)-psirz(i+1,j-1))/dz/2.
  t3=-1./x1*(psirz(i,j+2)-psirz(i,j))/dz/2.
  t4=-1./x2*(psirz(i+1,j+2)-psirz(i+1,j))/dz/2.
	call bilinear_average_ef(x1,x2,y1,y2,r0,z0,t1,t2,t3,t4,br0)
	br0=-br0

!bz
  t1=1./x1*(psirz(i+1,j)-psirz(i-1,j))/dr/2.
  t2=1./x2*(psirz(i+2,j)-psirz(i,j))/dr/2.
  t3=1./x1*(psirz(i+1,j+1)-psirz(i-1,j+1))/dr/2.
  t4=1./x2*(psirz(i+2,j+1)-psirz(i,j+1))/dr/2.
	call bilinear_average_ef(x1,x2,y1,y2,r0,z0,t1,t2,t3,t4,bz0)
	bz0=-bz0

!brr
  t1=(-1./r(i+1)*(psirz(i+1,j+1)-psirz(i+1,j-1))/dz/2.+1./r(i-1)*(psirz(i-1,j+1)-psirz(i-1,j-1))/dz/2.)/dr/2.
  t2=(-1./r(i+2)*(psirz(i+2,j+1)-psirz(i+2,j-1))/dz/2.+1./r(i)*(psirz(i,j+1)-psirz(i,j-1))/dz/2.)/dr/2.
  t3=(-1./r(i+1)*(psirz(i+1,j+2)-psirz(i+1,j))/dz/2.+1./r(i-1)*(psirz(i-1,j+2)-psirz(i-1,j))/dz/2.)/dr/2.
  t4=(-1./r(i+2)*(psirz(i+2,j+2)-psirz(i+2,j))/dz/2.+1./r(i)*(psirz(i,j+2)-psirz(i,j))/dz/2.)/dr/2.
	call bilinear_average_ef(x1,x2,y1,y2,r0,z0,t1,t2,t3,t4,brr)
	brr=-brr

!brz
  t1=(-1./r(i)*(psirz(i,j+1)-psirz(i,j))/dz+1./r(i)*(psirz(i,j)-psirz(i,j-1))/dz)/dz
  t2=(-1./r(i+1)*(psirz(i+1,j+1)-psirz(i+1,j))/dz+1./r(i+1)*(psirz(i+1,j)-psirz(i+1,j-1))/dz)/dz
  t3=(-1./r(i)*(psirz(i,j+2)-psirz(i,j+1))/dz+1./r(i)*(psirz(i,j+1)-psirz(i,j))/dz)/dz
  t4=(-1./r(i+1)*(psirz(i+1,j+2)-psirz(i+1,j+1))/dz+1./r(i+1)*(psirz(i+1,j+1)-psirz(i+1,j))/dz)/dz
	call bilinear_average_ef(x1,x2,y1,y2,r0,z0,t1,t2,t3,t4,brz)
	brz=-brz

!bzr
  t1=(+2./(r(i)+r(i+1))*(psirz(i+1,j)-psirz(i,j))/dr-2./(r(i)+r(i-1))*(psirz(i,j)-psirz(i-1,j))/dr)/dr
  t2=(+2./(r(i+2)+r(i+1))*(psirz(i+2,j)-psirz(i+1,j))/dr-2./(r(i)+r(i+1))*(psirz(i+1,j)-psirz(i,j))/dr)/dr
  t3=(+2./(r(i)+r(i+1))*(psirz(i+1,j+1)-psirz(i,j+1))/dr-2./(r(i)+r(i-1))*(psirz(i,j+1)-psirz(i,j+1))/dr)/dr
  t4=(+2./(r(i+2)+r(i+1))*(psirz(i+2,j+1)-psirz(i+1,j+1))/dr-2./(r(i)+r(i+1))*(psirz(i+1,j+1)-psirz(i,j+1))/dr)/dr
	call bilinear_average_ef(x1,x2,y1,y2,r0,z0,t1,t2,t3,t4,bzr)
	bzr=-bzr

!bzz
  t1=(+1./(r(i))*(psirz(i+1,j+1)-psirz(i-1,j+1))/dr/2.-1./(r(i))*(psirz(i+1,j-1)-psirz(i-1,j-1))/dr/2.)/dz/2.
  t2=(+1./(r(i+1))*(psirz(i+2,j+1)-psirz(i,j+1))/dr/2.-1./(r(i+1))*(psirz(i+2,j-1)-psirz(i,j-1))/dr/2.)/dz/2.
  t3=(+1./(r(i))*(psirz(i+1,j+2)-psirz(i-1,j+2))/dr/2.-1./(r(i))*(psirz(i+1,j)-psirz(i-1,j))/dr/2.)/dz/2.
  t4=(+1./(r(i+1))*(psirz(i+2,j+2)-psirz(i,j+2))/dr/2.-1./(r(i+1))*(psirz(i+2,j)-psirz(i,j))/dr/2.)/dz/2.
	call bilinear_average_ef(x1,x2,y1,y2,r0,z0,t1,t2,t3,t4,bzz)
	bzz=-bzz
	
	return
	end



	subroutine find_fields_interp_ef_psionly(r0,z0,psi0) !give back psi,br,bz at r0,z0

	use ef_circuit

	implicit none
	integer i,j
	double precision r0,z0,psi0,br0,bz0
	double precision x1,x2,x3,x4,x5,x6,x7,x8,y1,y2,y3,y4,y5,y6
	double precision z1,z2,z3,z4,z5,z6,z7,z8,t1,t2,t3,t4,t5,t6
	double precision brr,brz,bzr,bzz
	
	i=min(nr1,max(1,floor((r0-rmin)/dr+1.)))	
	j=min(nr1,max(1,floor((z0-zmin)/dz+1.)))	
	x1=r(i)
	x2=r(i+1)
	y1=z(j)
	y2=z(j+1)
	z1=psirz(i,j)
	z2=psirz(i+1,j)
	z3=psirz(i,j+1)
	z4=psirz(i+1,j+1)

!bilinear interpolation
	call bilinear_average_ef(x1,x2,y1,y2,r0,z0,z1,z2,z3,z4,psi0)

	
	return
	end



	subroutine find_fields_interp_ef_green(r0,z0,psi0,iconduc) !give back psi,br,bz at r0,z0

	use ef_circuit
	use green_matrix       ! declaration of minimal CPOs

	implicit none
	integer i,j,iconduc
	double precision r0,z0,psi0,br0,bz0
	double precision x1,x2,x3,x4,x5,x6,x7,x8,y1,y2,y3,y4,y5,y6
	double precision z1,z2,z3,z4,z5,z6,z7,z8,t1,t2,t3,t4,t5,t6
	double precision brr,brz,bzr,bzz
	
	i=floor((r0-rmin)/dr+1.)	
	j=floor((z0-zmin)/dz+1.)	
	x1=r(i)
	x2=r(i+1)
	y1=z(j)
	y2=z(j+1)
	z1=greeni(i,j,iconduc)
	z2=greeni(i+1,j,iconduc)
	z3=greeni(i,j+1,iconduc)
	z4=greeni(i+1,j+1,iconduc)

!bilinear interpolation
	call bilinear_average_ef(x1,x2,y1,y2,r0,z0,z1,z2,z3,z4,psi0)

	
	return
	end


	subroutine find_fields_interp_ef_psiext(r0,z0,psi0) !give back psi,br,bz at r0,z0

	use ef_circuit

	implicit none
	integer i,j
	double precision r0,z0,psi0,br0,bz0
	double precision x1,x2,x3,x4,x5,x6,x7,x8,y1,y2,y3,y4,y5,y6
	double precision z1,z2,z3,z4,z5,z6,z7,z8,t1,t2,t3,t4,t5,t6
	double precision brr,brz,bzr,bzz
	
	i=floor((r0-rmin)/dr+1.)	
	j=floor((z0-zmin)/dz+1.)	
	x1=r(i)
	x2=r(i+1)
	y1=z(j)
	y2=z(j+1)
	z1=psiextrz(i,j)
	z2=psiextrz(i+1,j)
	z3=psiextrz(i,j+1)
	z4=psiextrz(i+1,j+1)

!bilinear interpolation
	call bilinear_average_ef(x1,x2,y1,y2,r0,z0,z1,z2,z3,z4,psi0)

	
	return
	end



	subroutine find_fields_interp_ef_psionly_neg(r0,z0,psi0) !give back psi,br,bz at r0,z0

	use ef_circuit

	implicit none
	integer i,j
	double precision r0,z0,psi0,br0,bz0
	double precision x1,x2,x3,x4,x5,x6,x7,x8,y1,y2,y3,y4,y5,y6
	double precision z1,z2,z3,z4,z5,z6,z7,z8,t1,t2,t3,t4,t5,t6
	double precision brr,brz,bzr,bzz
	
	

	
	i=floor((r0-rmin)/dr+1.)	
	j=floor((z0-zmin)/dz+1.)	



	x1=r(i)
	x2=r(i+1)
	y1=z(j)
	y2=z(j+1)
	z1=-psirz(i,j)
	z2=-psirz(i+1,j)
	z3=-psirz(i,j+1)
	z4=-psirz(i+1,j+1)

!bilinear interpolation
	call bilinear_average_ef(x1,x2,y1,y2,r0,z0,z1,z2,z3,z4,psi0)

	
	return
	end



	subroutine bilinear_average_ef(x1,x2,y1,y2,x,y,f11,f21,f12,f22,f0)

	implicit none
	
	double precision x1,x2,y1,y2,x,y,f11,f21,f12,f22,f0
	
	f0=1./((x2-x1)*(y2-y1))*( &
	&  f11*(x2-x)*(y2-y)+ & 
	&  f21*(x-x1)*(y2-y)+ & 
	&  f12*(x2-x)*(y-y1)+ &
	&  f22*(x-x1)*(y-y1) )
	
	return
	end

!-----------------------------------------------------------------------------------
    function ellE_green(X, DL) ! gives back the first kind elliptic integral, from K. Lackner, T. Lunt, IPP - Garching 2022

    double precision, intent(in) :: X, DL
    double precision :: ellE_green
    ellE_green = (((.01736506451D0 *X + .04757383546D0)*X + .06260601220D0)*X + .44325141463D0)*X + 1.0D0 - &
           (((.00526449639D0 *X + .04069697526D0)*X + .09200180037D0)*X + .24998368310D0)*X*DL

    end function ellE_green

!-----------------------------------------------------------------------------------
    function ellK_green(X, DL) ! gives back the second kind elliptic integral, from K. Lackner, T. Lunt, IPP - Garching 2022

   double precision, intent(in) :: X, DL
    double precision :: ellK_green
    ellK_green= ((( .01451196212D0*X + .03742563713D0)*X + .03590092383D0)*X + .09666344259D0)*X + 1.38629436112D0 - &
          ((((.00441787012D0*X + .03328355346D0)*X + .06880248576D0)*X + .12498593597D0)*X + .5D0)*DL
 
    end function ellK_green

!----------------------------------------------------------------------------------- 

!-----------------------------------------------------------------------------------
    function ellE_green_dx(X,DL) ! d/dx

    double precision, intent(in) :: X,DL
    double precision :: ellE_green_dx
    ellE_green_dx = (((4.*.01736506451D0 *X + 3.*.04757383546D0)*X + 2.*.06260601220D0)*X + .44325141463D0) - &
           (((4.*.00526449639D0 *X + 3.*.04069697526D0)*X + 2.*.09200180037D0)*X + .24998368310D0)*DL + &
           (((.00526449639D0 *X + .04069697526D0)*X + .09200180037D0)*X + .24998368310D0) 
    end function ellE_green_dx

!-----------------------------------------------------------------------------------
    function ellK_green_dx(X,DL) ! d/dx

   double precision, intent(in) :: X,DL
    double precision :: ellK_green_dx
    ellK_green_dx= ((( 4.*.01451196212D0*X + 3.*.03742563713D0)*X + 2.*.03590092383D0)*X + .09666344259D0)  - &
          ((((4.*.00441787012D0*X + 3.*.03328355346D0)*X + 2.*.06880248576D0)*X + .12498593597D0))*DL+ &
          ((((.00441787012D0*X + .03328355346D0)*X + .06880248576D0)*X + .12498593597D0)*X + .5D0)/X
 
    end function ellK_green_dx

!----------------------------------------------------------------------------------- 

!-----------------------------------------------------------------------------------
    function ellE_green_asy(X, DL) ! for x --> 0

    double precision, intent(in) :: X, DL
    double precision :: ellE_green_asy
    ellE_green_asy = .44325141463D0*X + 1.0D0 - .24998368310D0*X*DL

    end function ellE_green_asy

!-----------------------------------------------------------------------------------
    function ellK_green_asy(X, DL) ! for x-->0

   double precision, intent(in) :: X, DL
    double precision :: ellK_green_asy
    ellK_green_asy= .09666344259D0*X + 1.38629436112D0 - (.12498593597D0*X + .5D0)*DL
 
    end function ellK_green_asy

!----------------------------------------------------------------------------------- 

!-----------------------------------------------------------------------------------
    function dk_dr1(r1,r2,z1,z2)  !identical do dk_dr2

   double precision, intent(in) :: r1,r2,z1,z2
    double precision :: dk_dr1,f
		f=((r2+r1)**2.+(z2-z1)**2.)
		dk_dr1 = sqrt(r2/r1)*sqrt(f)-2*sqrt(r1*r2)*(r1+r2)*f**(-1.5)
 
    end function dk_dr1
    
		function dk_dz1(r1,r2,z1,z2)  !identical to -dk_dz2

   double precision, intent(in) :: r1,r2,z1,z2
    double precision :: dk_dz1,f
		f=((r2+r1)**2.+(z2-z1)**2.)
		dk_dz1 = -2*sqrt(r1*r2)*(z1-z2)*f**(-1.5)
 
    end function dk_dz1

!----------------------------------------------------------------------------------- 
	
	subroutine green_function(r1,z1,r2,z2,greenf)
	
	implicit none
	
	double precision r1,z1,r2,z2,greenf
	double precision TT,K,ELCK,ELCE,ellk_green,elle_green
	double precision s21bbf,s21bcf,acl,alg
	integer ifailk,ifaile
	



	K=sqrt(4.*r1*r2/ & 
	 & ((r2+r1)**2.+(z2-z1)**2.))

	TT = 1.-K**2.

	acl=tt
	alg=dlog(acl)

	ELCK=ellK_green(acl,alg) !              S21BBF(0.D0,TT,1.D0,IFAILK)
	ELCE=ellE_green(acl,alg) ! ELCK-K**2./3.D0*S21BCF(0.D0,TT,1.D0,IFAILE)
	
	greenf = ( (1.D0-K**2./2.)*ELCK-ELCE )*( SQRT(r1*r2)/K )



	return
	end

!----------------------------------------------------------------------------------- 
	
	subroutine green_function_includingsamepoint(r1,z1,r2,z2,dl,R0,greenf)
	
	implicit none
	
	double precision, intent(in) :: r1,r2,z1,z2,dl,r0
	double precision, intent(out) :: greenf
	double precision TT,K,ELCK,ELCE,ellk_green,elle_green
	double precision s21bbf,s21bcf,acl,alg
	integer ifailk,ifaile
	


	if (abs(r1-r2).lt.1.e-6.and.abs(z1-z2).lt.1.e-6) then
	
		tt=dl/(4*R0)
!		greenf=-tt*(log(tt**2.)-4*log(2.)+1.)/4.*8.*r0**2./dl/2.	 
		greenf=-tt*(log(tt**2.)-4*log(2.)+2.)/4.*8.*r0**2./dl/2.	  !as in lackner code !this one works

!  tt = dl/(4*R0)
!		greenf ~ -0.25*(log(tt^2)-4*log(2.) + 2)/4.*8.*r0**2./dl/2.	  lackner version

	else
		K=sqrt(4.*r1*r2/ & 
	 & ((r2+r1)**2.+(z2-z1)**2.))

		TT = 1.-K**2.

		acl=tt
		alg=dlog(acl)

		ELCK=ellK_green(acl,alg) !              S21BBF(0.D0,TT,1.D0,IFAILK)
		ELCE=ellE_green(acl,alg) ! ELCK-K**2./3.D0*S21BCF(0.D0,TT,1.D0,IFAILE)
	
		greenf = ( (1.D0-K**2./2.)*ELCK-ELCE )*( SQRT(r1*r2)/K )
	endif


	return
	end
	
	subroutine green_function_identity(r1,z1,greenf,dr,dz,ntype)
	
	implicit none
	
	double precision r1,z1,r2,z2,greenf,dr,dz
	double precision TT,K,ELCK,ELCE
	double precision s21bbf,s21bcf
	integer ifailk,ifaile,ntype
	integer i,j

	if (ntype == 2) then
		greenf = r1*(log(8.*r1/(0.2236*(dr+dz)))-2.0) ! From A. Kavin, used in SPIDER (A. A. Ivanov and S. Yu. Medvedev), self inductance of a rectangular coil in toroidal direction
	endif
	
	if (ntype == 1) then
          elce=SQRT(dr**2+dz**2)
          greenf=0.
          do i=-1,1
          do j=-1,i
      	    if(i .eq. j) then
     	    	 greenf = greenf+ (r1+dr*i/3.d0) &
     	   &  * (log( 8.*(r1+dr*i/3.d0)/(elce/3.d0) ) - 0.5D0)
     	     else
							call green_function(r1+dr*i/3.,dz*i/3.,r1+dr*j/3.,dz*j/3.,elck)
		          greenf=greenf+4.*elck
    	      endif
          enddo
          enddo
          greenf = greenf/9.d0
	endif	


	return
	end
	

	subroutine green_function_non_identity(r1,z1,r2,z2,greenf,dr,dz,dr2,dz2,ntype1,ntype2)
	
	implicit none
	
	double precision r1,z1,r2,z2,greenf,dr,dz,dr2,dz2
	double precision TT,K,ELCK,ELCE
	double precision s21bbf,s21bcf
	integer ifailk,ifaile,ntype1,ntype2
	integer i,j

	greenf=0.

	if ( ntype1.eq.1.and.ntype2.eq.1) then
          do i=-1,1
          do j=-1,1
						call green_function(r1+dr*i/3.,z1+dz*i/3.,r2+dr2*j/3.,z2+dz2*j/3.,elck)
	          greenf=greenf+elck
          enddo
          enddo
          greenf = greenf/9.d0
	endif
	if ( ntype1.eq.1.and.ntype2.eq.2) then
          do i=-1,1
						call green_function(r1+dr*i/3.,z1+dz*i/3.,r2,z2,elck)
	          greenf=greenf+elck
          enddo
          greenf = greenf/3.d0
	endif
	if ( ntype1.eq.2.and.ntype2.eq.1) then
          do j=-1,1
						call green_function(r1,z1,r2+dr2*j/3.,z2+dz2*j/3.,elck)
	          greenf=greenf+elck
          enddo
          greenf = greenf/3.d0
	endif
	if ( ntype1.eq.2.and.ntype2.eq.2) then
						call green_function(r1,z1,r2,z2,elck)
	          greenf=greenf+elck
	endif


	return
	end





	subroutine boundary_ef(g)

! new bc is integral_over_boundary of -Green * dg/dn * dl
	use ef_circuit

	implicit none
	integer i,j,k
	double precision integr(i_dim2,4),g(i_dim2,i_dim2),dgdn(i_dim2),greenf
	integer jcounty
	
	integr=0.
	jcounty=0
! lower side
	do i=2,nr1
		call bgint_ef(r(i),z(1),integr(i,1),g,jcounty)
	enddo
! right side
	do i=2,nz1
		call bgint_ef(r(nr2),z(i),integr(i,2),g,jcounty)
	enddo
! upper side
	do i=2,nr1
		call bgint_ef(r(i),z(nz2),integr(i,3),g,jcounty)
	enddo
! left side
	do i=2,nz1
		call bgint_ef(r(1),z(i),integr(i,4),g,jcounty)
	enddo
	g(2:nr1,1)=integr(2:nr1,1)/GPI
	g(nr2,2:nz1)=integr(2:nz1,2)/GPI
	g(2:nr1,nz2)=integr(2:nr1,3)/GPI
	g(1,2:nz1)=integr(2:nz1,4)/GPI
	
	return
	end


	subroutine bgint_ef(r0,z0,bgintsol,g,jcounty)

! calculates  integral_over_boundary of -Green * dg/dn * dl for point r0,z0
	use ef_circuit

	implicit none

	integer i,j,k,jcounty
	double precision g(i_dim2,i_dim2),dgdn(i_dim2),greenf
	double precision r0,z0,bgintsol
	
	bgintsol=0.
! lower side
		do j=2,nr1
			jcounty=jcounty+1
			greenf=green_bnd_f(jcounty)
			dgdn(j)=-g(j,2)/dz*greenf*dr/r(j)
		enddo
		bgintsol=bgintsol-sum(dgdn(2:nr1))
!right side	
		do j=2,nz1
			jcounty=jcounty+1
			greenf=green_bnd_f(jcounty)
			dgdn(j)=-g(nr1,j)/dr*greenf*dz/(r(nr2)+r(nr1))*2.
		enddo
		bgintsol=bgintsol-sum(dgdn(2:nz1))
! upper side
		do j=2,nr1
			jcounty=jcounty+1
			greenf=green_bnd_f(jcounty)
			dgdn(j)=-g(j,nz1)/dz*greenf*dr/r(j)
		enddo
		bgintsol=bgintsol-sum(dgdn(2:nr1))
!left side	
		do j=2,nz1
			jcounty=jcounty+1
			greenf=green_bnd_f(jcounty)
			dgdn(j)=-g(2,j)/dr*greenf*dz/(r(1)+r(2))*2.
		enddo
		bgintsol=bgintsol-sum(dgdn(2:nz1))
	
	
	return
	end



	subroutine check_xpoint_connection_axis(rx,zx,rax,zax,dr,dz,icheck) !!
	
	!this routine checks that going from axis to x point, the directed gradient of psi never changes sign
	!(otherwise it means the x point is not connected to the plasma
	
	implicit none
	
	double precision, intent(in) :: rx,zx,rax,zax,dr,dz
	integer, intent(out) :: icheck
	
	double precision :: angl,dbl,dd,t1,t2,t3,t4,t5
	double precision :: z1,z2,psiold,z3
	integer nsteps,i

	call find_angle_ef(rax,zax,rx,zx,angl)
	
	dbl=sqrt((rx-rax)**2.+(zx-zax)**2.)
	dd=sqrt(dr**2.+dz**2.)
	nsteps=nint(dbl/dd)
	dd=dbl/nsteps !perfect ratio	
	
	
	icheck=1
	psiold=0.	
	do i=2,nsteps
		t1=rax+dd*(i-1)*cos(angl)
		t2=zax+dd*(i-1)*sin(angl)
		t3=rax+dd*i*cos(angl)
		t4=zax+dd*i*sin(angl)
		call find_fields_interp_ef_psionly(t1,t2,z1)
		call find_fields_interp_ef_psionly(t3,t4,z2)
		z3=(z2-z1)*psiold
		psiold=(z2-z1)
		if (z3<0) then
			icheck=0
			goto 300		
		endif
		if (dd*i.ge.dbl) then
			goto 300		
		endif
	enddo

300	continue	

	return
	end



	subroutine green_integral_element(x0,y0,x1,y1,x2,y2,greeni) !compute accurately the line integral of the green function.

! compute integral of green function from x1,y1 to x2,y2, when x0,y0 is at the center of the element

	implicit none
	
	
	double precision, intent(in) :: x0,y0,x1,y1,x2,y2
	double precision, intent(out) :: greeni
	
	double precision :: t1,t2,t3,t4,t5,t6,t7,t8,t9,t10
	double precision :: dl,d0,xavg,yavg
	
! computes the expression: int_(x1,y1)^(x2,y2)dl (line integral) of  G(x'-->x) 

	xavg=0.5*(x1+x2)
	yavg=0.5*(y1+y2)

	dl = sqrt((x2-x1)**2.+(y2-y1)**2.) !element length
	d0 = sqrt((xavg-x0)**2.+(yavg-y0)**2.) !distance of x0 from center of element
	
	
	greeni=0.

	return
	end
	
	
	
	subroutine ainv_matrix_def(A)
	
	implicit none
	double precision A(9,9)
	
	A(1,1:9)= (/  0.2500  , -0.5000  ,  0.2500 ,  -0.5000 ,   1.0000 ,  -0.5000  ,  0.2500 ,  -0.5000  ,  0.2500 /)
	A(2,1:9)= (/  -0.2500 ,   0.5000 ,  -0.2500   ,      0.  ,       0.   ,      0.  ,  0.2500 ,  -0.5000 ,   0.2500  /)
	A(3,1:9)= (/  -0.2500  ,       0. ,   0.2500  ,  0.5000  ,       0.  , -0.5000 ,  -0.2500  ,       0. ,   0.2500  /)
	A(4,1:9)= (/  0.2500  ,       0. ,  -0.2500  ,       0.  ,       0.  ,       0. ,  -0.2500  ,       0.  ,  0.2500  /)
	A(5,1:9)= (/  0.   ,      0.    ,     0. ,   0.5000 ,  -1.0000  ,  0.5000   ,      0.    ,     0.   ,      0.  /)
	A(6,1:9)= (/  0.  ,  0.5000    ,     0.   ,      0. ,  -1.0000  ,       0.  ,       0. ,   0.5000   ,      0.  /)
	A(7,1:9)= (/  0.  ,       0.   ,      0.  , -0.5000 ,        0.  ,  0.5000  ,       0.  ,       0.  ,       0.  /)
	A(8,1:9)= (/  0. ,  -0.5000   ,      0.   ,     0.   ,      0.  ,       0.  ,       0. ,   0.5000    ,     0.  /)
	A(9,1:9)= (/  0.  ,       0.   ,      0.   ,      0. ,   1.0000  ,       0. ,        0.  ,       0.  ,       0. /)
	
	
	
	return
	end






! fast sine transform from Zhongde Wang, Signal Processing 19 (1990) 91-102 Elsevier


	subroutine fst1(f,g,d,l,m,n)
	
	dimension f(n-1),g(n-1),d(n-1),l(n-1)
	m1=m-1
	n1=n/2
	n2=n1-1
	do 10 i=1,n2
	g(i)=f(i)+f(n-i)
10	f(n-i)=f(i)-f(n-i)
	g(n1)=f(n1)	
	call fst3(g,d,m1,n1)
	n0=0
20	n0=n1+n0
n1=n1/2
m1=m1-1
n2=n1-1
n3=n0+n1
n4=n3+n3
do 30 i=n0+1,n3-1
g(i)=f(i)+f(n4-i)
30	f(n4-i)=f(n4-i)-f(i)
g(n3)=f(n3)
	call fst3(g(n0+1),d,m1,n1)
	if (n1.ge.2) goto 20
	g(n-1)=f(n-1)
	do 40 i=1,n-1
40	f(l(i))=g(i)
return
end


subroutine fst2(f,d,m,n)
dimension f(n),d(1)
n1=n/2
n2=2
j3=1
do 10 i=1,n1
j3=-j3
i2=i+i
i1=i2-1
t=f(i1)
f(i1)=(t+f(i2))*d(n1+i-j3)
10	f(i2)=t-f(i2)
do 50 i=1,m-1
n1=n1/2
n7=n2
n2=n2+n2
n6=0
if (n1.eq.1) j3=0
do 40 j=1,n1
j3=-j3
j2=n1+j
b=d(j2)+d(j2)
n3=n6+1
n6=n6+n2
n8=n6-n7
n5=n8-1
n4=n8+n8
t=f(n8)
f(n8)=(t+f(n6))*d(j2-j3)
f(n6)=t-f(n6)
do 20 k=n3,n5
k1=n7+k
t=f(k)
f(k)=t+f(k1)
20 f(k1)=b*(t-f(k1))
do 30 k=n3,n5
30 f(n4-k)=f(n4-k)+f(k)
40 continue
50 continue
f(n)=d(1)*f(n)
return
end



subroutine fst3(f,d,m,n)
dimension f(n),d(1)
j1=1
n1=n
j3=-1
f(n)=d(1)*f(n)
do 40 i=1,m-1
n2=n1
n1=n1/2
n6=0
do 30 j=1,j1
j2=j1+j
j3=-j3
a=d(j2-j3)
b=d(j2)+d(j2)
n3=n6+1
n6=n2+n6
n5=n6-n1
n4=n5+n5
do 10 k=n3,n5-1
k1=n4-k
f(k)=f(k)+f(k1)
10 f(k1)=f(k1)*b
f(n5)=f(n5)*a
do 20 k=n3,n5
t=f(k+n1)
f(k+n1)=f(k)-t
20 f(k)=f(k)+t
30 continue
40 j1=j1+j1
50 do 60 i=1,j1
i1=i+i
i2=i1-1
j3=-j3
t=d(j1+i-j3)*f(i2)
f(i2)=t+f(i1)
60 f(i1)=t-f(i1)
return
end

subroutine fst4(f,d,m,n)
dimension f(n),d(1)
j1=1
n1=n
do 40 i =1,m
j2=j1
n2=n1
n1=n1/2
n7=0
do 30 j=1,j1
a=d(j2+j)
a=a+a
n6=n7
n7=n7+n2
n3=n6+1
n4=n3+n7
n5=n6+n1
do 10 k=n3,n5
k1=n4-k
f(k)=f(k)+f(k1)
10 f(k1)=a*f(k1)
do 20 k=n3,n5
t=f(k+n1)
f(k+n1)=f(k)-t
20 f(k)=f(k)+t
30 continue
40 j1=j1+j1
do 50 i=1,j1
i1=i+i
i2=i1-1
f(i1)=d(j1+j2)*f(i1)
50 f(i2)=D(j1+i1)*f(i2)
return
end


subroutine coefs(d,n1)
dimension d(n1),c(15)
double precision c,t
c(1)=sqrt(0.5)
d(1)=c(1)
m=2
m0=1
10 do 20 i=m0,m-1
i1=i+i
i2=i1+1
c(i1)=sqrt(0.5*(1.+c(i)))
c(i2)=sqrt(0.5*(1.-c(i)))
d(i1)=c(i1)
20 d(i2)=c(i2)
m0=m0+m0
m=m0+m0
if (m.ge.n1) goto 70
if (m.ge.16) goto 30
goto 10
30 do 40 i=m0,m-1
i1=i+i
i2=i1+1
t=sqrt(0.5*(1.+c(i)))
d(i+i)=t
t=sqrt(0.5*(1.-c(i)))
40 d(i+i+1)=t
m0=m0+m0
m=m0+m0
if (m.ge.n1) goto 70
50 do 60 i=m0,m-1
d(i+i)=sqrt(.5*(1.+d(i)))
60 d(i+i+1)=sqrt(.5*(1.-d(i)))
m0=m0+m0
m=m0+m0
if (m.ge.n1) goto 70
goto 50
70 do 80 i=n1,2,-1
80 d(i)=d(i-1)
return
end


subroutine hdmod(f,g,n)
dimension f(n),g(1)
if (n.le.2) return
n1=n/2
n0=n1
n2=1
10 i1=n
i2=-n2
do 20 i=1,n1
i1=i1-n2
i2=i2+n2
do 20 j=1,n2
20 g(i2+j)=f(i1+j)
n3=0
n4=n2+n2
i1=n
do 30 i=1,n1
n3=n3+n2
i1=i1-n4
i2=n0-n3
do 30 j=1,n2
30 f(i1+j)=f(i2+j)
i1=-n2
i2=-n2
do 40 i=1,n1
i1=i1+n4
i2=i2+n2
do 40 j=1,n2
40 f(i1+j)=g(i2+j)
n2=n2+n2
n1=n1/2
if (n1.ge.2) goto 10
return
end


subroutine ivhdm(f,g,n)
dimension f(n),g(1)
if (n.le.2) return
n1=n/4
i1=1
10 n2=n1+n1
n3=n2+n2
i2=0
do 40 i=1,i1
i3=i2
i2=i2+n3
i4=i3+n2
do 20 j=1,n2
20 g(j)=f(i3+j+j)
do 30 j=2,n2
30 f(i3+j)=f(i3+j+j-1)
do 40 j=1,n1
j1=j+j
j2=j1-1
f(i4+j)=G(j2)
40 f(i4+j2)=g(j1)
n1=n1/2
i1=i1+i1
if (n1.ge.1) goto 10
return
end

subroutine reord(l,n)
dimension l(n)
l(1)=1
l(2)=3
n2=2
10 n1=n2
n2=n1+n1
do 20 i=n1,1,-1
i1=i+i
l(i1)=n2-l(i)
20 l(i1-1)=l(i)
if (n1.lt.n/2) goto 10
n0=0
n2=n/2
n3=n2
30 n2=n2/2
do 40 i=1,n2
i1=n0+i+i-1
40 l(n3+i)=l(i1)+l(i1)
if (n2.eq.1) return
n0=n3
n3=n3+n2
goto 30
return
end








subroutine re_psi_feqis(psi_new,psi,nr,nt,rb,zb,r0o,z0o,r0,z0, &
				& teta,lambda,dldr,dtdr,dldz,dtdz,Rmaj,ZZ,psiax)  ! transfer PSI from old axis to new axis
				

implicit none

integer, intent(IN) :: nr,nt
double precision, intent(IN) :: r0,z0,r0o,z0o,psiax
double precision, intent(in) , dimension(nt):: rb,zb,teta
double precision, intent(in) , dimension(nr,nt):: psi,lambda,dldr,dtdr,dldz,dtdz,Rmaj,ZZ
double precision, intent(out) :: psi_new(nr,nt)

double precision :: dr,dz,dz0,dr0, & 
	& dr1,dz1,rl1(nr,nt),zl1(nr,nt), & 
	& gpsidpsi(nr,nt),tetaprim(nt),db1(nt),drl(nr,nt),dzl(nr,nt)


double precision :: dpsidl(nr,nt),dpsidteta(nr,nt)

integer i,j,k,ii,jj,kk

!calculate new psi
psi_new=psi
psi_new(1,:)=psiax

dr=r0-r0o
dz=z0-z0o

!write(*,*) 'repsi',r0o,z0o,r0,z0,psiax

!calculate thetaprime
do j=1,nt
	dz0=Zb(j)-Z0o
	dr0=Rb(j)-R0o
		call find_angle_ef(r0,z0,rb(j),zb(j),tetaprim(j))
!	tetaprim(j)=teta(j)+(dz0*dr-dr0*dz)/(dr0**2.+dz0**2.)
	db1(j)=sqrt((Rb(j)-r0)**2.+(Zb(j)-z0)**2.)
	do i=2,nr-1
		rl1(i,j)=r0+lambda(i,j)*db1(j)*cos(tetaprim(j))
		zl1(i,j)=z0+lambda(i,j)*db1(j)*sin(tetaprim(j))
		drl(i,j)=rl1(i,j)-rmaj(i,j)
		dzl(i,j)=zl1(i,j)-zz(i,j)
	enddo
enddo

	if (tetaprim(1).gt.tetaprim(nt)) then !reorder
		j=1
		do k=1,nt-1
			if (tetaprim(k+1).lt.tetaprim(k)) j=k+1   ! j is the first teta above 0
		enddo		
		if (j.gt.1)	tetaprim(1:j-1)=tetaprim(1:j-1)-3.141592*2.
	endif

!write(*,*) tetaprim(1:nt)

! calculate gradpsi
do i=2,nr-1
	do j=2,nt-1
		dpsidl(i,j)=(psi(i+1,j)-psi(i-1,j))/(lambda(i+1,j)-lambda(i-1,j))
		dpsidteta(i,j)=(psi(i,j+1)-psi(i,j-1))/(teta(j+1)-teta(j-1))	
	enddo
	j=1
		dpsidl(i,j)=(psi(i+1,j)-psi(i-1,j))/(lambda(i+1,j)-lambda(i-1,j))
		dpsidteta(i,j)=(psi(i,j+1)-psi(i,nt))/(teta(j+1)-(teta(nt)-3.141592*2.))	
	j=nt
		dpsidl(i,j)=(psi(i+1,j)-psi(i-1,j))/(lambda(i+1,j)-lambda(i-1,j))
		dpsidteta(i,j)=(psi(i,1)-psi(i,j-1))/(teta(1)+3.141592*2.-teta(j-1))			
enddo
	
	
!calculate new psi
do i=2,nr-1
	do j=1,nt
		psi_new(i,j)=psi(i,j)+ & 
		& ((dpsidl(i,j)*dldr(i,j)+dpsidteta(i,j)*dtdr(i,j))*drl(i,j)+ & 
		& (dpsidl(i,j)*dldz(i,j)+dpsidteta(i,j)*dtdz(i,j))*dzl(i,j))
	enddo
enddo	


psi_new=0.49*psi_new+0.51*psi

!write(*,*) sum(abs(psi-psi_new))

return
end







subroutine find_new_X0Y0(Nr, Nt, PSI, X, Y, X0, Y0, & 
& psiax, iax, jax,derivs)

implicit none
integer, intent(in) :: Nr, Nt
double precision, intent(in), dimension(Nr, Nt) :: PSI, X, Y
integer, intent(out) :: iax, jax
double precision, intent(out) :: X0, Y0, psiax
integer :: jt, info, jmin(2),jj,j,ii(nt)
double precision :: g3(1*Nt+1)
double precision :: matrix(1*Nt+1, 6)
double precision :: derivs(8),c(6),yy(nt)
double precision :: xs(1000),ys(1000),fun(1000)
integer nsh
double precision :: det,dp(5), work(2*(1*nt+1)*6)


	do j=1,nt
		ii(j)=minloc(psi(1:nr,j),1)
		yy(j)=psi(ii(j),j)
	enddo
	jj=minloc(yy(1:nt),1)
	psiax=psi(ii(jj),jj)
	x0=x(ii(jj),jj)
	y0=y(ii(jj),jj)
	iax=ii(jj)
	jax=jj

	if (iax.ne.1) return
	if (jax.ne.1) return


!refine axis 
	g3(2:nt+1)=psi(2,1:nt)
	g3(1)=psi(1,1)
	matrix(1,1)=x(1,1)**2.		 
	matrix(1,2)=x(1,1)		 
	matrix(:,3)=1.		 
	matrix(1,4)=y(1,1)**2.		 
	matrix(1,5)=y(1,1)
	matrix(1,6)=x(1,1)*y(1,1)
	
	nsh=1
	xs(1)=x(1,1)
	ys(1)=y(1,1)
	fun(1)=psi(1,1)
	do jj=1,nt
	matrix(jj+1,1)=x(2,jj)**2.		 
	matrix(jj+1,2)=x(2,jj)		 
	matrix(jj+1,4)=y(2,jj)**2.		 
	matrix(jj+1,5)=y(2,jj)		 
	matrix(jj+1,6)=x(2,jj)*y(2,jj)		 
	nsh=nsh+1
	xs(nsh)=x(2,jj)
	ys(nsh)=y(2,jj)
	fun(nsh)=psi(2,jj)
	enddo 

!        call dgels('N', 1*nt + 1, 6, 1, matrix, 1*nt + 1, g3, 1*nt + 1, WORK, 2*(1*nt + 1)*6, INFO)

        det = 4*g3(1)*g3(4) - g3(6)**2
        x0 = (g3(6)*g3(5) - 2*g3(4)*g3(2))/det
        y0 = (g3(6)*g3(2) - 2*g3(1)*g3(5))/det
        psiax = g3(1)*x0**2 + g3(2)*x0 + g3(3) +  &
                g3(4)*y0**2 + g3(5)*y0 + g3(6)*x0*y0




	 

return
end subroutine find_new_X0Y0









subroutine find_closest_xpoints(rx,zx,ierr,n_add)

use ef_circuit, only: rbnd,zbnd,nteta,r,z,nr1,nz1, & 
	& lim_maxR,lim_minR,lim_minZ,lim_maxZ,dr,dz

!this routine finds the x-points close to the plasma boundary, irrespective of other x-points

implicit none

integer i,j,i0,j0,i1,j1,iinc,inow,jnow,nx,ierr,n_add,jinc
double precision posx(2),x1,tolez,toleb,rx(20),zx(20)
double precision ddipsi(5),bx0,bx1,bx2
integer jcycl(250),istart


bx1=sqrt(dr**2.+dz**2.)
tolez=1000
rx=1000.
zx=1000.
nx=0
ierr=1
n_add=0
jinc=0

!write(*,*) 'find clos',nteta,rbnd(1:nteta),zbnd(1:nteta)
if (rbnd(1).le.r(1)) return ! boundary doesnt exist yet
ierr=0

do i=1,nteta !cycle over boundary points
!around each boundary point, do a 3-layer X-point search (25 point search x boundary point)
!	write(*,*) 'find closest x point',i,rbnd(i),zbnd(i)
	call nine_point_regression_follow(rbnd(i),zbnd(i),posx,ddipsi,x1)
!	write(*,*) 'find closest x point',i,rbnd(i),zbnd(i),posx,x1
	if (isnan(x1)) then
		else
		if (posx(1).le.lim_maxR.and.posx(1).ge.lim_minR.and.posx(2).ge.lim_minZ.and.posx(2).le.lim_maxZ) then
			jinc=jinc+1
			jcycl(jinc)=i
		endif
	endif
enddo


if (jinc.eq.0) then
	ierr=1
	return
endif

!write(*,*) lim_maxR,lim_minR,lim_minZ,lim_maxZ

call nine_point_regression_follow(rbnd(jcycl(1)),zbnd(jcycl(1)),posx,ddipsi,x1)
n_add=1
rx(1)=posx(1)
zx(1)=posx(2)

if (jinc.eq.1) then
	return
endif

if (jinc.ge.2) then
	!remove double counts
	do i=2,jinc
		call nine_point_regression_follow(rbnd(jcycl(i)),zbnd(jcycl(i)),posx,ddipsi,x1)
		bx0=sqrt((rx(i-1)-posx(1))**2.+(zx(i-1)-posx(2))**2.)
		if (bx0.le.bx1) then
		else
			n_add=n_add+1
			rx(n_add)=posx(1)
			zx(n_add)=posx(2)		
		endif		
	enddo
endif







!write(*,*) 'nadd',n_add,rx(1:n_add),zx(1:n_add)

if (n_add.eq.0) ierr=1


return
end subroutine find_closest_xpoints
















subroutine sort_PSI_values(y1, &
     &   Nx)

	implicit none
	
	integer j,i,Nx
	double precision y1(Nx),z(nx)
	double precision dum1(Nx),dum2(Nx),dum3(nx)
	double precision z1,z2,z3,dum4(nx)
	
	dum2=y1
	dum3=z
	do i=1,Nx
 	j=minloc(dum2,1)
		dum1(i)=dum2(j)
		dum4(i)=dum3(j)
		dum2(j)=1.E+08
	enddo
	y1=dum1
	z=dum4
	end




	subroutine port_from_fbe_ef(Nx,nt,xx,yy,PSI)
	use ef_circuit
	
	implicit none
	
	integer nx,nt,i,j,i0,j0
	double precision xx(nx,nt),yy(nx,nt),psi(nx,nt)
	double precision dum1,x1,x2,y1,y2,z1,z2,z3,z4
	double precision r0,z0
	
	do j=1,nt
	do i=1,nx
		r0=xx(i,j)
		z0=yy(i,j)
		i0=floor((r0-rmin)/dr+1.)	
		j0=floor((z0-zmin)/dz+1.)	
	x1=r(i0)
	x2=r(i0+1)
	y1=z(j0)
	y2=z(j0+1)
	z1=-psirz(i0,j0)
	z2=-psirz(i0+1,j0)
	z3=-psirz(i0,j0+1)
	z4=-psirz(i0+1,j0+1)
!bilinear interpolation
	call bilinear_average_ef(x1,x2,y1,y2,r0,z0,z1,z2,z3,z4,psi(i,j))

	enddo
	enddo
	
	

	return
	end
	
	



















