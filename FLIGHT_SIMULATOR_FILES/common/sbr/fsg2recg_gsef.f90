!C----------------------------------------------------------------------|
	subroutine	fsg2recg_gsef(Nr,Nt,PSIf,PSIb, &
     &  amin,IPOL, &
     &  imethod,X,Y,r_a,thetap, &
     &  R2,Z2,NR2,NZ2,PSIr, &
     &  B_Rrect,B_Zrect,B_Trect,amp_fac)

	use parameters_a2equil
	implicit none

	integer Nt,Nr,i1,j1,i2,j2,k,NR2,NZ2,j,i,ipos
	integer p1,p2,p3,p4,p5,use_kernel,II(Nr),JJ
	integer imethod,Nrt,Ntt,nfour
	parameter(Nrt=500,Ntt=500)
	double precision PSIf(Nr,Nt),PSIr(NR2,NZ2),amp_fac
	double precision X(Nr,Nt),R2(NR2)
	double precision Y(Nr,Nt),Z2(NZ2)
	double precision B_R(Nr,Nt),IPOL(Nr)
	double precision B_Z(Nr,Nt),ipolp(NR2,NZ2)
	double precision r_a(Nr,Nt)
	double precision x_a(Nr,Nt),r_a2(Nr,Nt+1)
	double precision PSI_test(Nr,Nt)
	double precision thetap(Nt),thetap2(Nt+1)
	double precision B_T(Nr,Nt),amin
	double precision distance_rz(Nr,Nt)
	double precision distance_tet(Nt)
	double precision kernel_int(Nr,Nt)
	double precision kernel_in1(3,3)
	double precision kernel_in3(2,3)
	double precision kernel_in2(2,Nt)
	double precision g1,g2,g3,g4
	double precision t1,t2,t3,t4
	double precision rr1,zz1,rr2,zz2,PSIb
	double precision B_Rrect(NR2,NZ2)
	double precision B_Zrect(NR2,NZ2)
	double precision dumm(3,3),dum2(2,Nt)
	double precision psimm(3,3),psim2(2,Nt)
	double precision brmm(3,3),brm2(2,Nt)
	double precision bzmm(3,3),bzm2(2,Nt)
	double precision btmm(3,3),btm2(2,Nt)
	double precision dum3(2,3),tum1(3)
	complex CHRMTp(cheb_degree,1+2*four_degree)
	complex CHRMTbr(cheb_degree,1+2*four_degree)
	complex CHRMTbz(cheb_degree,1+2*four_degree)
	complex CHRMTbt(cheb_degree,1+2*four_degree)
	double precision psim3(2,3)
	double precision dumrr(NR2),dumrz(NZ2)
	double precision brm3(2,3)
	double precision bzm3(2,3),dXb0
	double precision tht(NR2,NZ2),rhr(NR2,NZ2)
	double precision btm3(2,3),rmax
	integer mfour(1+2*four_degree)
	double precision B_Trect(NR2,NZ2),yy(Nr)
	
	!write(*,*) 'ddd1'

	thetap2(1:Nt)=thetap(1:Nt)
	thetap2(Nt+1)=GP2+thetap(1)
	r_a2(:,1:Nt)=r_a(:,1:Nt)
	r_a2(:,Nt+1)=r_a(:,1)

	g1 = maxval(maxval(X,2),1)*(1.+amp_fac/100.)
	g2 = maxval(maxval(Y,2),1)*(1.+amp_fac/100.)
	g3 = minval(minval(X,2),1)*(1.-amp_fac/100.)
	g4 = minval(minval(Y,2),1)*(1.+amp_fac/100.)
	
	do j2=1,NZ2
	Z2(j2)=g4+(g2-g4)*(j2-1)/(NZ2-1) 	
	enddo
		do i2=1,NR2
	R2(i2)=g3+(g1-g3)*(i2-1)/(NR2-1) 	
	enddo

!!	write(*,*) 'ddd2'

!C Biquadratic interpolation
	if (imethod.eq.4) then
	
	do j=1,NZ2
		do i=1,NR2
			dXb0=((R2(i)-X(1,1))**2.0+ &
     &  (Z2(j)-Y(1,1))**2.0)**0.5
			if ((Z2(j)-Y(1,1)).ge.0 .and. &
     &   (R2(i)-X(1,1)).ge.0) then 
				tht(i,j)=atan((Z2(j)-Y(1,1))/ &
     &  (R2(i)-X(1,1)))
			endif
			if ((Z2(j)-Y(1,1)).ge.0 .and. &
     &   (R2(i)-X(1,1)).lt.0) then 
				tht(i,j)=GP+atan((Z2(j) &
     &  -Y(1,1))/(R2(i)-X(1,1)))
			endif
			if ((Z2(j)-Y(1,1)).lt.0 .and. &
     &   (R2(i)-X(1,1)).lt.0) then 
				tht(i,j)=GP+atan((Z2(j)- &
     &  Y(1,1))/(R2(i)-X(1,1)))
			endif
			if ((Z2(j)-Y(1,1)).lt.0 .and.  &
     &  (R2(i)-X(1,1)).ge.0) then 
				tht(i,j)=GP2+atan((Z2(j)- &
     &  Y(1,1))/(R2(i)-X(1,1)))
			endif
			rhr(i,j)=dXb0
			rr2 = R2(i)
			zz2 = Z2(j)

!	write(*,*) 'ddd3',r2(i),z2(j)

			distance_rz = ((X-rr2)**2.0+(Y-zz2)**2.0)**0.5
		do j1=1,Nr
		II(j1)=minloc(distance_rz(j1,:),1)
		yy(j1)=minval(distance_rz(j1,:),1)
	enddo
		JJ=minloc(yy,1)
		t1=minval(yy,1)
		k = II(JJ)

	distance_tet=((thetap-tht(i,j))**2.0)**0.5
		k = minloc(distance_tet,1)
!		write(*,*) 'ddd4',JJ,t1,k,distance_tet


		if (JJ.gt.1 .and. JJ.lt.Nr) then
			if (k.gt.1 .and. k.lt.Nt) then
		dumm=r_a(JJ-1:JJ+1,k-1:k+1)
		tum1=thetap(k-1:k+1)
		psimm=PSIf(JJ-1:JJ+1,k-1:k+1)
		brmm(:,1)=IPOL(JJ-1:JJ+1)
		brmm(:,2)=IPOL(JJ-1:JJ+1)
		brmm(:,3)=IPOL(JJ-1:JJ+1)
			endif
			if (k.eq.1) then
		dumm(:,1) = r_a(JJ-1:JJ+1,Nt)
		dumm(:,2:3)=r_a(JJ-1:JJ+1,1:2)
	if (tht(i,j).lt.GP) then
		tum1(1)=thetap(Nt)-GP2
		tum1(2:3)=thetap(1:2)
	else
		tum1(1)=thetap(Nt)
		tum1(2:3)=thetap(1:2)+GP2
	endif
		psimm(:,1)=PSIf(JJ-1:JJ+1,Nt)
		psimm(:,2:3)=PSIf(JJ-1:JJ+1,k:k+1)
		brmm(:,1)=IPOL(JJ-1:JJ+1)
		brmm(:,2)=IPOL(JJ-1:JJ+1)
		brmm(:,3)=IPOL(JJ-1:JJ+1)
			endif
			if (k.eq.Nt) then
		dumm=r_a2(JJ-1:JJ+1,k-1:k+1)
	if (tht(i,j).gt.GP) then
		tum1=thetap2(k-1:k+1)
	else
		tum1=thetap2(k-1:k+1)-GP2
	endif
		psimm(:,1:2)=PSIf(JJ-1:JJ+1,k-1:k)
		psimm(:,3)=PSIf(JJ-1:JJ+1,1)
		brmm(:,1)=IPOL(JJ-1:JJ+1)
		brmm(:,2)=IPOL(JJ-1:JJ+1)
		brmm(:,3)=IPOL(JJ-1:JJ+1)
			endif
	ipos=1
			endif

		if (JJ.eq.Nr) then
			if (k.gt.1 .and. k.lt.Nt) then
		dumm=r_a(JJ-2:JJ,k-1:k+1)
	tum1=thetap(k-1:k+1)
	psimm=PSIf(JJ-2:JJ,k-1:k+1)
			endif
			if (k.eq.1) then
		dumm(:,1) = r_a(JJ-2:JJ,Nt)
		dumm(:,2:3)=r_a(JJ-2:JJ,k:k+1)
	if (tht(i,j).lt.GP) then
		tum1(1)=thetap(Nt)-GP2
		tum1(2:3)=thetap(1:2)
	else
		tum1(1)=thetap(Nt)
		tum1(2:3)=thetap(1:2)+GP2
	endif
		psimm(:,1)=PSIf(JJ-2:JJ,Nt)
		psimm(:,2:3)=PSIf(JJ-2:JJ,k:k+1)
			endif
			if (k.eq.Nt) then
		dumm=r_a2(JJ-2:JJ,k-1:k+1)
	if (tht(i,j).gt.GP) then
		tum1=thetap2(k-1:k+1)
	else
		tum1=thetap2(k-1:k+1)-GP2
	endif
		psimm(:,1:2)=PSIf(JJ-2:JJ,k-1:k)
		psimm(:,3)=PSIf(JJ-2:JJ,1)
			endif
		brmm(:,1)=IPOL(Nr)*(/1.,1.,1./)
		brmm(:,2)=IPOL(Nr)*(/1.,1.,1./)
		brmm(:,3)=IPOL(Nr)*(/1.,1.,1./)
	ipos=1
			endif

		if (JJ.eq.1) then
				dum2=r_a(1:2,:)
				psim2=PSIf(1:2,:)
	do j2=1,Nt
		brm2(:,j2)=IPOL(1:2)
	enddo
	ipos=0
			endif

	if (ipos.eq.1) then
	call biquadratic_interp_gsef(dumm,(tum1),psimm, &
     &  rhr(i,j),(tht(i,j)), &
     &  PSIr(i,j))
	call biquadratic_interp_gsef(dumm,(tum1),brmm, &
     &  rhr(i,j),(tht(i,j)), &
     &  ipolp(i,j))
	endif
	if (ipos.eq.0) then
	call biquadratic_interp_gsef0(dum2,(thetap),psim2, &
     &  rhr(i,j),(tht(i,j)),Nt, &
     &  PSIr(i,j))
	call biquadratic_interp_gsef0(dum2,(thetap),brm2, &
     &  rhr(i,j),(tht(i,j)),Nt, &
     &  ipolp(i,j))
	endif
	if (ipos.eq.2) then
	call liquadratic_interp_gsef(dum3,(tum1),psim3, &
     &  rhr(i,j),(tht(i,j)), &
     &  PSIr(i,j),i,j)
	call liquadratic_interp_gsef(dum3,(tum1),brm3, &
     &  rhr(i,j),(tht(i,j)), &
     &  ipolp(i,j),0,0) 
	endif



			enddo
	enddo
!		write(*,*) 'ddd5'

B_zrect=0.
B_trect=0.
b_rrect=0.
dumrr=0.
dumrz=0.

	do j=1,NZ2
			dumrr(2:nr2-1)=(psir(3:nr2,j)-psir(1:nr2-2,j))/(R2(3:nr2)-R2(1:nr2-2))
			B_Zrect(:,j)=1./(GP2*R2)*dumrr
			B_Trect(:,j)=ipolp(:,j)/R2
	enddo
	do i=1,NR2
			dumrz(2:nz2-1)=(psir(i,3:nz2)-psir(i,1:nz2-2))/(z2(3:nz2)-z2(1:nz2-2))
			B_Rrect(i,:)=-1./(GP2*R2(i))*dumrz
	enddo
	
!write(*,*) 'dd6'
!	open(1,file=name_gsefdir(1:kname_gsef)//'diag_rect.dat');
! 	write(1,*) Nr2,Nz2
!		write(1,*) ((tht(i,j),i=1,Nr2),j=1,Nz2)
!		write(1,*) ((rhr(i,j),i=1,Nr2),j=1,Nz2)
!	close(1)

	endif

return

	end
!C end of GS SOLVER EF



		subroutine biquadratic_interp_gsef(r1,thet1,y1, &
     &  r2,thet2, &
     &  y2)

	implicit none
	integer i,j,k
	double precision r2,y2,A,B,C,yy(3)
	double precision r1(3,3),thet1(3),y1(3,3)
	double precision r0(3),p0(3),t1,thet2,t2,t3,z1,z2,z3

	do j=1,3
		r0=r1(:,j)
		p0=y1(:,j)
		t1=p0(1)
		t2=p0(2)
		t3=p0(3)
		z1=r0(1)
		z2=r0(2)
		z3=r0(3)
				A=(t3-t2-(z3-z2)*(t1-t2)/(z1-z2)) &
     &     /((z3-z2)*(z3-z1))
				B=(t1-t2)/(z1-z2) &
     &     -A*(z1+z2)
				C=t2-A*(z2**2.0)-B*z2

		yy(j)=A*(r2**2.0)+B*r2+C
	enddo

		t1=yy(1)
		t2=yy(2)
		t3=yy(3)
		z1=thet1(1)
		z2=thet1(2)
		z3=thet1(3)
				A=(t3-t2-(z3-z2)*(t1-t2)/(z1-z2)) &
     &     /((z3-z2)*(z3-z1))
				B=(t1-t2)/(z1-z2) &
     &     -A*(z1+z2) 
				C=t2-A*(z2**2.0)-B*z2

		y2=A*(thet2**2.0)+B*thet2+C

	
	end


	subroutine biquadratic_interp_gsef0(r1,t1,y1, &
     &  r2,t2,Nt, &
     &  y2)

	implicit none
	integer i,j,k,Nt,jok,jplus,jdone,GP2
	double precision r2,t2,y2,a,b,c,d(Nt)
	double precision r1(2,Nt),t1(Nt),y1(2,Nt)
	
	GP2=6.283185
	jok=0
	jplus=Nt
	jdone=0

	do i=1,Nt
		d(i)=y1(1,i)+(r2-r1(1,i))*(y1(2,i)-y1(1,i))/(r1(2,i)-r1(1,i))
	if (t2.eq.t1(i) .and. jdone.eq.0) then
	jok=i
	jdone=1
	endif
	if (t2.lt.t1(i) .and. jdone.eq.0) then
	jok=0
	jplus=i
	jdone=1
	endif
	enddo		

	if (jok.eq.1) then
	y2=d(jok)	
	return
	endif

	if (jplus.eq.1) then
	a=t1(Nt)-GP2
	b=t1(1)
	y2=d(Nt)+(t2-a)*(d(1)-d(Nt))/(b-a)
	return
	endif

	if (jplus.eq.Nt) then
	a=t1(Nt)
	b=t1(1)+GP2
	y2=d(Nt)+(t2-a)*(d(1)-d(Nt))/(b-a)
	return
	endif

	a=t1(jplus)
	b=t1(jplus-1)
	y2=d(jplus)+(t2-a)*(d(jplus-1)-d(jplus))/(b-a)
	
	end
	
	
	
	subroutine liquadratic_interp_gsef(r1,thet1,y1, &
     &  r2,thet2, &
     &  y2,j1,j2)

	implicit none
	integer i,j,k,j1,j2
	double precision r2,y2,A,B,C,yy(3)
	double precision r1(2,3),thet1(3),y1(2,3)
	double precision r0(2),p0(2),t1,thet2,t2,t3,z1,z2,z3

	do j=1,3
		r0=r1(:,j)
		p0=y1(:,j)
		t1=p0(1)
		t2=p0(2)
		z1=r0(1)
		z2=r0(2)
		yy(j)=(t2-t1)/(z2-z1)*(r2-z1)+t1
	enddo

		t1=yy(1)
		t2=yy(2)
		t3=yy(3)
		z1=thet1(1)
		z2=thet1(2)
		z3=thet1(3)
				A=(t3-t2-(z3-z2)*(t1-t2)/(z1-z2)) &
     &     /((z3-z2)*(z3-z1))
				B=(t1-t2)/(z1-z2) &
     &     -A*(z1+z2)
				C=t2-A*(z2**2.0)-B*z2

		y2=A*(thet2**2.0)+B*thet2+C
	
	end






	subroutine rule_for_interp_gsef(g,d,k,n1,n2,n10,n20)
	
	implicit none
	integer n1,n2,i,j,n10,n20
	double precision g,d(n1,n2),k(n1,n2)	

	k=exp(-10.*(1.-g/d))
	k(n10,n20)=1.0
	
	
		
	end

		

	subroutine CHRMT_comp_gsef(nfour,mfour, &
     &  PSIf,x_a,thetap,Nr,Nt,CHRMT)

	use parameters_a2equil
	implicit none

	integer Nt,Nr,i,j,k,nfour,mfour(nfour)
	double precision PSIf(Nr,Nt)
	complex CHRMT(cheb_degree,nfour),kernel(Nr-1,Nt)
	double precision x_a(Nr,Nt),thetap(Nt),gradxa(Nr-1)
	double precision cheb_fac,xmid(Nr-1),psimid(Nr-1)
	double precision chebpol(Nr-1),gtheta
! testing orthogonality
	double precision chebtest(cheb_degree,cheb_degree)
	double precision chebpol2(Nr-1),kernelt(Nr-1)



	cheb_fac=4.0/GP
		do j=1,cheb_degree
		do i=1,cheb_degree
		 if (i.eq.1) cheb_fac=2.0/GP
	if (i.gt.1) cheb_fac=4.0/GP

			k=1
				xmid=(x_a(2:Nr,k)+x_a(1:Nr-1,k))/2.0
				gradxa=x_a(2:Nr,k)-x_a(1:Nr-1,k)

	chebpol=(cheb_coefs(i,1)+ &
     &   cheb_coefs(i,2)*xmid**2.0+ &
     &   cheb_coefs(i,3)*xmid**4.0+ &
     &   cheb_coefs(i,4)*xmid**6.0+ &
     &   cheb_coefs(i,5)*xmid**8.0)
	chebpol2=(cheb_coefs(j,1)+ &
     &   cheb_coefs(j,2)*xmid**2.0+ &
     &   cheb_coefs(j,3)*xmid**4.0+ &
     &   cheb_coefs(j,4)*xmid**6.0+ &
     &   cheb_coefs(j,5)*xmid**8.0)
					
				kernelt=gradxa/(1.-xmid**2.0)**0.5
				kernelt=kernelt*chebpol*chebpol2				
				
				chebtest(i,j)=cheb_fac*sum(kernelt,1)
	enddo
	enddo
	do  i=1,cheb_degree

	enddo

	
	

	do j=1,nfour
		do i=1,cheb_degree
		 if (i.eq.1) cheb_fac=2.0/GP
	if (i.gt.1) cheb_fac=4.0/GP


			do k=1,Nt
				xmid=(x_a(2:Nr,k)+x_a(1:Nr-1,k))
				psimid=PSIf(2:Nr,k)+PSIf(1:Nr-1,k)
				gradxa=x_a(2:Nr,k)-x_a(1:Nr-1,k)
	if (k.lt.Nt) then
				xmid=(xmid+x_a(2:Nr,k+1)+x_a(1:Nr-1,k+1))/4.0
				psimid=(psimid+PSIf(2:Nr,k+1)+PSIf(1:Nr-1,k+1))/4.0
				gtheta=(thetap(k+1)-thetap(k))	
	endif
	if (k.eq.Nt) then
				xmid=(xmid+x_a(2:Nr,1)+x_a(1:Nr-1,1))/4.0
				psimid=(psimid+PSIf(2:Nr,1)+PSIf(1:Nr-1,1))/4.0
				gtheta=(GP2+thetap(1)-thetap(Nt))	
	endif
		
	chebpol=cheb_fac*(cheb_coefs(i,1)+ &
     &   cheb_coefs(i,2)*xmid**2.0+ &
     &   cheb_coefs(i,3)*xmid**4.0+ &
     &   cheb_coefs(i,4)*xmid**6.0+ &
     &   cheb_coefs(i,5)*xmid**8.0)
					
				kernel(:,k)=gradxa/(1.-xmid**2.0)**0.5
				kernel(:,k)=kernel(:,k)*chebpol* &
     & exp(-CMPLX(0,1)*mfour(j)*thetap(k))*gtheta
				
				kernel(:,k)=psimid*kernel(:,k)
	enddo

				CHRMT(i,j)=sum(sum(kernel,2),1)
				CHRMT(i,j)=sum(sum(kernel,2),1)

	enddo
	enddo
	

	end
	
	
	
	subroutine CHRMT_to_PSI_gsef(CHRMT,rhr,tht,nf,mf,rmax, &
     &  yy)

	use parameters_a2equil
	implicit none

	integer i,j,k,nf,mf(nf)
	double precision PSIr
	complex CHRMT(cheb_degree,nf)
	double precision yy,dum1,dum2,dum3
	double precision rhr,tht,chebpol,rmax

	dum3=rhr/rmax

	if (dum3.le.1.0) then
	dum1=0.0
	do j=1,nf
	do i=1,cheb_degree
	chebpol=cheb_coefs(i,1)+ &
     &  cheb_coefs(i,2)*dum3**2.0+ &
     &  cheb_coefs(i,3)*dum3**4.0+ &
     &  cheb_coefs(i,4)*dum3**6.0+ &
     &  cheb_coefs(i,5)*dum3**8.0
	dum2=real(CHRMT(i,j)*chebpol*exp(CMPLX(0,1)*mf(j)*tht))
	dum1=dum1+dum2

	enddo
	enddo
	yy=1/GP2*dum1
	endif

	if (dum3.gt.1.0) then
	dum3=1.0
	dum1=0.0
	do j=1,nf
	do i=1,cheb_degree
	chebpol=cheb_coefs(i,1)+ &
     &  cheb_coefs(i,2)*dum3**2.0+ &
     &  cheb_coefs(i,3)*dum3**4.0+ &
     &  cheb_coefs(i,4)*dum3**6.0+ &
     &  cheb_coefs(i,5)*dum3**8.0
	dum2=real(CHRMT(i,j)*chebpol*exp(CMPLX(0,1)*mf(j)*tht))
	dum1=dum1+dum2

	enddo
	enddo
	yy=1/GP2*dum1
	endif
		
	end	
	
	
	
