subroutine qinterp_feqis(x1, y1, Nx1, x2, y2, Nx2)

implicit none

integer, intent(in) :: Nx1, Nx2
double precision, intent(in) , dimension(Nx1) :: x1, y1
double precision, intent(in) , dimension(Nx2) :: x2
double precision, intent(out), dimension(Nx2) :: y2

integer :: i, j
double precision :: A, B, C, z1, z2, z3, t1, t2, t3, t4

do i=1, Nx2
    t4 = x2(i)

    do j=2, Nx1-1
        z1 = x1(j-1)
        z2 = x1(j)
        z3 = x1(j+1)

        if (t4 == z1) then
            y2(i) = y1(j-1)
            EXIT
        elseif (t4 == z2) then
            y2(i) = y1(j)
            EXIT
        elseif (t4 == z3) then
            y2(i) = y1(j+1)
            EXIT
        elseif (t4 > z1 .and. t4 < z3) then
            t1 = y1(j-1)
            t2 = y1(j)
            t3 = y1(j+1)
            A = (t3 - t2 - (z3 - z2)*(t1 - t2)/(z1 - z2)) / ((z3 - z2)*(z3 - z1))
            B = (t1 - t2)/(z1 - z2) - A*(z1 + z2)
            C = t2 - A*z2**2 - B*z2
            y2(i) = A*t4**2 + B*t4 + C
            EXIT
        elseif (t4 < z1 .and. j == 2) then
            z1 = x1(j-1)**2
            z2 = x1(j)**2
            t1 = y1(j-1)
            t2 = y1(j)
            A = 0.0
            B = (t1 - t2)/(z1 - z2) - A*(z1 + z2)
            C = t2 - A * z2**2 - B*z2
            y2(i) = A * t4**4 + B * t4**2 + C
            EXIT
        elseif (t4 > z3 .and. j == Nx1-1) then
            t1 = y1(j-1)
            t2 = y1(j)
            t3 = y1(j+1)
            A = (t3 - t2 - (z3 - z2)*(t1 - t2)/(z1 - z2)) / ((z3 - z2)*(z3 - z1))
            B = (t1 - t2)/(z1 - z2) - A*(z1 + z2)
            C = t2 - A * z2**2 - B*z2
            y2(i) = A * t4**2 + B*t4 + C
            EXIT
        endif
    enddo
enddo

return
end subroutine qinterp_feqis

!------------------------------------------------------------
subroutine linterp_feqis(x1, y1, Nx1, x2, y2, Nx2)

implicit none

integer, intent(in) :: Nx1, Nx2
double precision, intent(in) , dimension(Nx1) :: x1, y1
double precision, intent(in) , dimension(Nx2) :: x2
double precision, intent(out), dimension(Nx2) :: y2

integer :: i, j
double precision :: A, C, z1, z2, t1, t2, t4

do i=1, Nx2
    t4 = x2(i)
    do j=2, Nx1
        z1 = x1(j-1)
        z2 = x1(j)
        if (t4 == z1) then
            y2(i)=y1(j-1)
            EXIT
        elseif (t4 == z2) then
            y2(i)=y1(j)
            EXIT
        elseif (t4 > z1 .and. t4 < z2) then
            t1 = y1(j-1)
            t2 = y1(j)
            A = (t2 - t1)/(z2 - z1)
            C = t2 - A*z2
            y2(i) = A*t4 + C
            EXIT
        elseif (t4 < z1 .and.j == 2) then
            y2(i)=y1(j-1)
            EXIT
        elseif (t4 > z2 .and. j == nx1) then
            y2(i) = y1(j)
            EXIT
        endif
    enddo
enddo

return
end subroutine linterp_feqis

!------------------------------------------------------------
subroutine polyfitcc_feqis(x, y, P)

implicit none

double precision, intent(in) , dimension(3) :: x, y
double precision, intent(out), dimension(3) :: P

double precision :: y21, y32, x21, x32, h21, h32

y32 = y(3) - y(2)
y21 = y(2) - y(1)
x32 = x(3) - x(2)
x21 = x(2) - x(1)
h32 = x(3) + x(2)
h21 = x(2) + x(1)

P(1) = (x21*y32 - x32*y21)/(x21*x32*(h32 - h21))
P(2) = y21/x21 - P(1)*h21
P(3) = y(3) - P(1)*x(3)**2 - P(2)*x(3)

return
end subroutine polyfitcc_feqis










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








	subroutine qinterp_ef_equilef(x1,y1,Nx1,x2,y2,Nx2)
		 
		 	implicit none

	integer Nx1,Nx2,i,j,k,jdone

	double precision x1(Nx1)
	double precision y1(Nx1)
	double precision x2(Nx2)
	double precision y2(Nx2),A,B,C
	double precision z1,z2,z3,z4,z5,z11
	double precision t1,t2,t3,t4,t5,t11


	do i=1,Nx2
		jdone=0
		
		t4=x2(i)

		do j=2,Nx1-1
		z1=x1(j-1)
		z2=x1(j)
		z3=x1(j+1)
		
	if (t4.eq.z1 .and. jdone.eq.0) then
		y2(i)=y1(j-1)
		jdone=1
	endif

	if (t4.eq.z2 .and. jdone.eq.0) then
		y2(i)=y1(j)
		jdone=1
	endif

	if (t4.eq.z3 .and. jdone.eq.0) then
		y2(i)=y1(j+1)
		jdone=1
	endif
	
	if (t4.gt.z1 .and. t4.lt.z3 .and. jdone.eq.0) then
		t1=y1(j-1)
		t2=y1(j)
		t3=y1(j+1)
		
				A=(t3-t2-(z3-z2)*(t1-t2)/(z1-z2)) &
     &     /((z3-z2)*(z3-z1))
				B=(t1-t2)/(z1-z2) &
     &     -A*(z1+z2)
				C=t2-A*(z2**2.0)-B*z2

		y2(i)=A*(t4**2.0)+B*t4+C

		jdone=1
	endif

	if (t4.lt.z1 .and. jdone.eq.0  &
     &   .and. j.eq.2) then
		z1=x1(j-1)**2.0
		z2=x1(j)**2.0
!		z3=x1(j+1)**2.0
		t1=y1(j-1)
		t2=y1(j)
!		t3=y1(j+1)
		

!				A=(t3-t2-(z3-z2)*(t1-t2)/(z1-z2))
!     &     /((z3-z2)*(z3-z1))
!				B=(t1-t2)/(z1-z2)
!     &     -A*(z1+z2)
!				C=t2-A*(z2**2.0)-B*z2
				A=0.0
				B=(t1-t2)/(z1-z2) &
     &     -A*(z1+z2)
				C=t2-A*(z2**2.0)-B*z2


		y2(i)=A*(t4**4.0)+B*(t4**2.0)+C
!		y2(i)=B*t4**2.0+C
		jdone=1
	endif


	if (t4.gt.z3 .and. jdone.eq.0  &
     &   .and. j.eq.Nx1-1) then

		t1=y1(j-1)
		t2=y1(j)
		t3=y1(j+1)
		

				A=(t3-t2-(z3-z2)*(t1-t2)/(z1-z2)) &
     &     /((z3-z2)*(z3-z1))
				B=(t1-t2)/(z1-z2) &
     &     -A*(z1+z2)
				C=t2-A*(z2**2.0)-B*z2


		y2(i)=A*(t4**2.0)+B*t4+C
		jdone=1
	endif

		enddo
	enddo
	end 




	subroutine linterp_ef_equilef(x1,y1,Nx1,x2,y2,Nx2)
		 
		 	implicit none

	integer Nx1,Nx2,i,j,k,jdone

	double precision x1(Nx1)
	double precision y1(Nx1)
	double precision x2(Nx2)
	double precision y2(Nx2),A,B,C
	double precision z1,z2,z3,z4,z5,z11
	double precision t1,t2,t3,t4,t5,t11


	do i=1,Nx2
		jdone=0
		
		t4=x2(i)

		do j=2,Nx1
		z1=x1(j-1)
		z2=x1(j)
		
	if (t4.eq.z1 .and. jdone.eq.0) then
		y2(i)=y1(j-1)
		jdone=1
	endif

	if (t4.eq.z2 .and. jdone.eq.0) then
		y2(i)=y1(j)
		jdone=1
	endif
	
	if (t4.gt.z1 .and. t4.lt.z2 .and. jdone.eq.0) then
		t1=y1(j-1)
		t2=y1(j)
		
				A=(t2-t1)/(z2-z1)
				C=t2-A*z2

		y2(i)=A*t4+C

		jdone=1
	endif

	if (t4.lt.z1 .and. jdone.eq.0.and.j.eq.2) then
		y2(i)=y1(j-1)
		jdone=1
	endif


	if (t4.gt.z2 .and. jdone.eq.0.and.j.eq.nx1) then

		y2(i)=y1(j)
		jdone=1
	endif

		enddo
	enddo
	end 



	subroutine integrcc_equilef(nx,x,y,sy)

	implicit none
	integer i,j,nx
	double precision x(nx),y(nx),drho,y1tmp
	double precision sy(nx)

	sy=0.
	do i=2,nx
		drho=(x(i)-x(i-1))
		y1tmp=(y(i)+y(i-1))/2.
	 sy(i)=sy(i-1)+y1tmp*drho
	enddo

	end


	subroutine derivcc_equilef(nx,x,y,dy,gga)

	implicit none
	integer i,j,nx,gga
	double precision x(nx),y(nx),dy(nx)
	double precision P(3)
	
	dy=0.
	do i=2,nx-1
		dy(i)=(y(i+1)-y(i-1))/(x(i+1)-x(i-1))
	enddo
	if (gga.eq.1)	dy(1)=0.
	if (gga.eq.2) dy(1)=(y(2)-y(1))/(x(2)-x(1))

	call polyfitcc(x(nx-2:nx),y(nx-2:nx),P)

	dy(nx)=2.*P(1)*x(nx)+P(2)


	end

	subroutine derivcc2_equilef(nx,x,y,dy,itype)

	implicit none
	integer i,j,nx,itype
	double precision x(nx),y(nx),dy(nx)
	double precision P(3)

	dy=0.
	i=1
	if (itype.eq.1) then
		dy(i)=2.*(y(2)-y(1))/(x(2)**2.)
	else
	endif
	

	do i=2,nx-1
	dy(i)=2./(x(i+1)-x(i-1))*( &
     & (y(i+1)-y(i))/(x(i+1)-x(i)) - &
     & (y(i)-y(i-1))/(x(i)-x(i-1)) &
     &  )
	enddo


	call polyfitcc(x(nx-2:nx),y(nx-2:nx),P)
	dy(nx)=2.*P(1)

	if (itype.eq.2) then
	call polyfitcc(x(1:3),y(1:3),P)
	dy(1)=2.*P(1)
	endif

	end


	subroutine polyfitcc_equilef(x,y,P)
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
	
	integer i,j,k

! go from jrhoteta to jrz
	jrz=0.
	do j=1,nz2
	do i=1,nr2
		call curinterp_ef(r(i),z(j),jrhoteta(1:nrho-1,1:nteta), & 
		& rho(1:nrho-1,1:nteta),teta(1:nteta),raxp,zaxp,nrho-1,nteta,jrz(i,j))
	enddo
	enddo

	return
	end






	subroutine curinterp_ef(r,z,jrho,rho,teta,rax,zax,nrho,nteta,j)

	use pi_grec_vars
	implicit none
	integer nrho,nteta,i,k,j1,j2,j3,j4,k1,k2,k3,k4
	double precision r,z,jrho(nrho,nteta),j,rho(nrho,nteta),teta(nteta)
	double precision anglr,rho0,x1(2),y1(2),jt(2),rax,zax
		
	call find_angle_ef(rax,zax,r,z,anglr)
  rho0=sqrt((r-rax)**2.0+(z-zax)**2.0)
  
!	j1=nteta-1
!	write(*,*) rax,zax,r,z,anglr,teta(1:nteta)
	if (anglr.lt.teta(1)) anglr=anglr+GPI2
	do i=1,nteta-1
		if (anglr.ge.teta(i)) j1=i
	enddo
		j2=j1+1
!	write(3333,*) r,z,anglr,rho0,rho(nrho,j1),rho(nrho,j2),jrho(1,1),rax,zax


	if (rho0.gt.rho(nrho,j1).or. &
	& rho0.gt.rho(nrho,j2)) then
	j=0.
	return
	endif

	do i=1,nrho-1
		if (rho0.gt.rho(i,j1)) k1=i
	enddo	
		k2=k1+1
	do i=1,nrho-1
		if (rho0.gt.rho(i,j2)) k3=i
	enddo	
		k4=k3+1
	
	x1(1)=rho(k1,j1)
	x1(2)=rho(k2,j1)
	y1(1)=jrho(k1,j1)
	y1(2)=jrho(k2,j1)
!	write(3333,*) 'k',k1,k2,k3,k4,x1,y1
	call linterp_ef_equilef(x1,y1,2,rho0,jt(1),1)
	x1(1)=rho(k3,j2)
	x1(2)=rho(k4,j2)
	y1(1)=jrho(k3,j2)
	y1(2)=jrho(k4,j2)
!	write(3333,*) 'k',k1,k2,k3,k4,jt(1),x1,y1
	call linterp_ef_equilef(x1,y1,2,rho0,jt(2),1)
	x1(1)=teta(j1)
	x1(2)=teta(j2)
	y1(1)=jt(1)
	y1(2)=jt(2)
	call linterp_ef_equilef(x1,y1,2,anglr,j,1)
!	write(3333,*) 'k',k1,k2,k3,k4,jt(1),jt(2),x1,y1,j

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
	

















	subroutine coil_forces_feqis(ncoilz,force_R,force_Z) !

	use ef_circuit
	use green_matrix
	
	implicit none
	integer i,j,k,ncoilz
	double precision force_R(ncoilz),force_Z(ncoilz)
	double precision x1,nblock_a
	
	force_R=0.
	force_Z=0.
	nblock_a=nblocks-npassive
	do i=1,nblock_a
			x1=sum(jrz(1:nr2,1:nz2)* &
      & area_eff(1:nr2,1:nz2)*dgreenirpl(1:nr2,1:nz2,i))
			force_R(i)=force_R(i)+curconduc(mequivalence(i))*x1
			x1=sum(jrz(1:nr2,1:nz2)* &
      & area_eff(1:nr2,1:nz2)*dgreenizpl(1:nr2,1:nz2,i))
			force_Z(i)=force_Z(i)+curconduc(mequivalence(i))*x1
	enddo
	do i=nblock_a+1,nblocks
			x1=sum(jrz(1:nr2,1:nz2)* &
      & area_eff(1:nr2,1:nz2)*dgreenirpl(1:nr2,1:nz2,i))
			force_R(i)=force_R(i)+curconduc(i-nblock_a+nactive)*x1
			x1=sum(jrz(1:nr2,1:nz2)* &
      & area_eff(1:nr2,1:nz2)*dgreenizpl(1:nr2,1:nz2,i))
			force_Z(i)=force_Z(i)+curconduc(i-nblock_a+nactive)*x1
	enddo

!block-to-block
	do i=1,nblock_a
		do j=1,nblock_a
				force_R(i)=force_R(i)+ & 
				& curconduc(mequivalence(j))*curconduc(mequivalence(i))*dgreenirj(i,j)
				force_Z(i)=force_Z(i)+ & 
				& curconduc(mequivalence(j))*curconduc(mequivalence(i))*dgreenizj(i,j)
		enddo		
		do j=1+nblock_a,nblocks
				force_R(i)=force_R(i)+ & 
				& curconduc(j-nblock_a+nactive)*curconduc(mequivalence(i))*dgreenirj(i,j)
				force_Z(i)=force_Z(i)+ & 
				& curconduc(j-nblock_a+nactive)*curconduc(mequivalence(i))*dgreenizj(i,j)
		enddo		
	enddo
	
	
	
	do i=1+nblock_a,nblocks
		do j=1,nblock_a
				force_R(i)=force_R(i)+ & 
				& curconduc(mequivalence(j))*curconduc(i-nblock_a+nactive)*dgreenirj(i,j)
				force_Z(i)=force_Z(i)+ & 
				& curconduc(mequivalence(j))*curconduc(i-nblock_a+nactive)*dgreenizj(i,j)
		enddo		
		do j=1+nblock_a,nblocks
					force_R(i)=force_R(i)+ & 
					& curconduc(j-nblock_a+nactive)*curconduc(i-nblock_a+nactive)*dgreenirj(i,j)
					force_Z(i)=force_Z(i)+ & 
					& curconduc(j-nblock_a+nactive)*curconduc(i-nblock_a+nactive)*dgreenizj(i,j)
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
	psiext_out=psiext_out/dllt


	write(*,*) 'psibbb',psibnd,psiext_out


	return
	end




	subroutine psiplex_calc_ef(dumz)

	use exchange_with_astra
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

	use	exchange_with_astra
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
	
	subroutine green_function_identity(r1,z1,greenf,dr,dz)
	
	implicit none
	
	double precision r1,z1,r2,z2,greenf,dr,dz
	double precision TT,K,ELCK,ELCE
	double precision s21bbf,s21bcf
	integer ifailk,ifaile
	

	greenf = r1*(log(8.*r1/(0.2236*(dr+dz)))-2.0) ! From A. Kavin, used in spider, self inductance of a rectangular coil in toroidal direction



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
		do j=1,nr1
			jcounty=jcounty+1
			greenf=green_bnd_f(jcounty)
			dgdn(j)=-(g(j,2)-g(j,1)+g(j+1,2)-g(j+1,1))/2./dz*greenf*dr/(r(j)+dr/2.)
		enddo
		bgintsol=bgintsol-sum(dgdn(1:nr1))
!right side	
		do j=1,nz1
			jcounty=jcounty+1
			greenf=green_bnd_f(jcounty)
			dgdn(j)=(g(nr2,j)-g(nr1,j)+g(nr2,j+1)-g(nr1,j+1))/2./dr*greenf*dz/(r(nr2)+r(nr1))*2.
		enddo
		bgintsol=bgintsol-sum(dgdn(1:nz1))
! upper side
		do j=1,nr1
			jcounty=jcounty+1
			greenf=green_bnd_f(jcounty)
			dgdn(j)=(g(j,nz2)-g(j,nz1)+g(j+1,nz2)-g(j+1,nz1))/2./dz*greenf*dr/(r(j)+dr/2.)
		enddo
		bgintsol=bgintsol-sum(dgdn(1:nr1))
!left side	
		do j=1,nz1
			jcounty=jcounty+1
			greenf=green_bnd_f(jcounty)
			dgdn(j)=-(g(2,j)-g(1,j)+g(2,j+1)-g(1,j+1))/2./dr*greenf*dz/(r(1)+r(2))*2.
		enddo
		bgintsol=bgintsol-sum(dgdn(1:nz1))
	
	
	return
	end




