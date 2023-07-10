	subroutine PBMITER

	use parameter_inc
	use const_inc
	use status_inc

	implicit none

	integer	 J,i,istep
	real*8 ped_width,exbsh(NA1),chi_edge,mmm
	real*8 dum1,dum2(NA1),nnn,c1,c2,c3,c4,P_LH
	real*8 dum4,dum5,dum6,dum7,machinez
	real*8 pfavg,sp1avg,sp2avg,sp3avg
	integer nspec
	parameter (nspec=7)

	real*8 D_j(nspec),enrich(nspec)
	real*8 Dw_j(nspec)
	real*8 a_j(nspec)
	real*8 R_j(nspec),qdtr
	real*8 V_sol,V_div,V_pump
	real*8 T_pump
	real*8 dens_Sj(nspec)
	real*8 dens_Dj(nspec)
	real*8 dens_Sjo(nspec)
	real*8 dens_Djo(nspec)
	real*8 dens_Djoo(nspec,25000)
	real*8 Gplasma(nspec)
	real*8 delayz,bumb1,bumb2,bumb3
	real*8 bumb4
	real*8 losspflow
	real*8 Gvmp(nspec),Gvdiv(nspec)
	real*8 Gpell(nspec),Gneuty(nspec)
	real*8 A_1(nspec)
	real*8 B_1(nspec),n_elS,n_elD
	real*8 C_1(nspec)
	real*8 A_2(nspec)
	real*8 B_2(nspec)
	real*8 C_2(nspec),timeyy,veloc
	real*8 c_eros(nspec)
	integer jgipp

	data timeyy /0./
	data j_countm /0/
	save dens_Sjo,dens_Djo
	save dens_Djoo
	save timeyy,j_countm
	save pfavg,sp1avg,sp2avg,sp3avg,jgipp
	data bumb1 /0./
	data bumb2 /0./
	data bumb3 /0./
	save bumb1,bumb2,bumb3	

	integer imodel,j_countm
	double precision c_model

!	write(*,*) TIME,CDWM1,CSCL2




!	return

!This model specifically works for:
! 1 D , 2 T, 3 H, 4 He, 5 Be , 6 W,  7 Ne
	
!Definition of constants
!volumes
	V_sol=0.8*GP2*RTOR*GP2*&
     & ABC*(1.+ELON(NA1)**2./2.)**0.5  !here there is a length missing...
	V_div=(1.-0.8)*GP2*RTOR*GP2*&
     & ABC*(1.+ELON(NA1)**2./2.)**0.5  !here there is a length missing...
	V_pump=0.3*V_div
!diffus ioniz and recycl coeff
	D_j=RTOR/9.	
	Dw_j=0.1*RTOR/9.	
	a_j=1.	
	R_j=0.01	
	T_pump=CV13
	enrich=1.d0
	enrich(3)=1.2d0
	enrich(4)=1.2d0
	enrich(5)=1.d0
	enrich(6)=20.d0
	enrich(7)=20.d0

	bumb1=0.01
	c_eros(1)=bumb1*2. !D
	c_eros(2)=bumb1*3. !T
	c_eros(3)=bumb1*1. !H
	c_eros(4)=bumb1*4. !He
	c_eros(5)=bumb1*10. !Be
	c_eros(6)=bumb1*180. !W
	c_eros(7)=bumb1*20. !Ne
	

!time scheme
 !initial condition
	if (TIME.le.TAU) then
		dens_Sjo(1)=F1(na1)*V_sol	
		dens_Sjo(2)=F2(na1)*V_sol	
		dens_Sjo(3)=F3(na1)*V_sol	
		dens_Sjo(4)=F4(na1)*V_sol	
		dens_Sjo(5)=F6(na1)*V_sol	
		dens_Sjo(6)=F7(na1)*V_sol	!to check
		dens_Sjo(7)=F8(na1)*V_sol	!to check
	pfavg=0.
	sp1avg=0. 
	sp2avg=0.
	sp3avg=0.
	jgipp=0
	
	
	dens_Djo=dens_Sjo

	endif






! 1 D , 2 T, 3 H, 4 He, 5 Be, 6 W, 7 Ne
!solve per each particle species
!	do j=1,nspec
	
!equation: dy1/dt = A1 y1 + B1 y2 + C1	
!equation: dy2/dt = A2 y1 + B2 y2 + C2	

!valve flow
	Gvmp(1)=	(1.-CDYM3)*CV4 ! D valve  midplane p/s + 2*water
	Gvmp(2)=	CDYM3*CV4     !T valve midplane p/s
	Gvmp(3)=	0. !reinjected He by mistake  ! He
	Gvmp(4)=	0. !reinjected Xe by mistake !Xe
	Gvmp(5)=	10. !reinjected Ar by mistake !Ar
	Gvmp(6)=	0.; !+CDJM8 !reinjected W by mistake !W!+ wall source 0.1*teped
	Gvmp(7)=	0.; !reinjected Kr by mistake + Oxi puff mp


	Gvdiv(1)=	0.
	Gvdiv(2)=	0.
	Gvdiv(3)=	0.
	Gvdiv(4)=	0.
	Gvdiv(5)=	0.  !Ar puff
	Gvdiv(6)=	0.  ! W puff due to erosion, computed above
	Gvdiv(7)=	0.  !Kr puff 

!pellet flow + fusion products
	Gpell(1)=	CPEL1-1./5.632*QDTR(ROC)
	Gpell(2)=	CIMP3-1./5.632*QDTR(ROC)
	Gpell(4)=	1./5.632*QDTR(ROC)
	Gpell(3)=	0.
	Gpell(5)=	0.
	Gpell(6)=	0.
	Gpell(7)=	0.

!neutral plasma flow
	do j=1,nspec
	Gneuty(j)=(1.-a_j(j))*(Gvmp(j)+ &
     & (1.-R_j(j))*Dw_j(j)*dens_Sjo(j))	
	enddo
	
!plasma flow
	Gplasma(1)=QF1(na1)
	Gplasma(2)=QF2(na1)
	Gplasma(3)=QF3(na1)
	Gplasma(4)=QF4(Na1)
	Gplasma(5)=QF6(NA1)
	Gplasma(6)=QF7(NA1)
	Gplasma(7)=QF8(na1)


	do j=1,nspec
		Gplasma(j)=Gplasma(j)
	enddo

	NNWM = Dw_j(1)/(Dw_j(1)+D_j(1))*(Gneuty(1)+Gneuty(2))
	NNCL = D_j(1)/(Dw_j(1)+D_j(1))*(Gneuty(1)+Gneuty(2))

	ENCL= 0.01
	ENWM= 10.

	veloc=sqrt(ENCL*1.e3*9.5756e+07)
	NNCL=NNCL/veloc/SLAT(NA1)
	veloc=sqrt(ENWM*1.e3*9.5756e+07)
	NNWM=NNWM/veloc/SLAT(NA1)

	ZRD30=Gneuty(1)/(Gneuty(1)+Gneuty(2)+1.e-16)

! 1 is sol, 2 is div
!A1 = Dw_j*(-1+aj*(1-Rj))-Dj
!B1 = Dj/enrichj
!C1 = Gplasmaj + aj*Gvmpj

!A2 = Dj
!B2 = -Dj/enrich-Vpump/Vdiv^2 * Tpump
!C2 = Gvdivj

!matrix form: y1(1-A1*dt) + y2(-B1*dt) = y10+C1*dt	
!matrix form: y1(-A2*dt) + y2(1.-B2*dt) = y20+C2*dt	

!solution:
! y1 = 1/((1.-A1*dt)*(1.-B2*dt)-B1*A2*dt*dt)*(((1.-B2*dt))*(y10+C1*dt)+(B1*dt)*(y20+C2*dt))		
! y2 = 1/((1.-A1*dt)*(1.-B2*dt)-B1*A2*dt*dt)*((A2*dt)*(y10+C1*dt)+(1.-A1*dt)*(y20+C2*dt))		
	
	do j=1,nspec
		A_1(j)=Dw_j(j)*(-1.+a_j(j)*(1.-r_j(j)))-D_j(j)
		A_2(j)=D_j(j)
		B_1(j)=D_j(j)/enrich(j)
		B_2(j)=-D_j(j)/enrich(j) &
     & -V_pump/V_div**2.*T_pump
		C_1(j)=Gplasma(j)+a_j(j)*Gvmp(j)
		C_2(j)=Gvdiv(j)


	dum7=1./((1.-A_1(j)*tau)*(1.-B_2(j)*tau) &
     & -B_1(j)*A_2(j)*tau*tau)

	dens_Sj(j)=dum7*(((1.-B_2(j)*tau))*(dens_Sjo(j)+ &
     & C_1(j)*tau)+(B_1(j)*tau)*(dens_Djo(j)+C_2(j)*tau))			

	dens_Dj(j)=dum7*((A_2(j)*tau)*(dens_Sjo(j)+C_1(j)*tau)+ &
     & (1.-A_1(j)*tau)*(dens_Djo(j)+C_2(j)*tau))		

	enddo	
	



!assign final
	dens_Sjo=dens_Sj
	dens_Djo=dens_Dj


!assign back to boundary conditions
	CDBC1 = dens_Sjo(1)/V_sol !f1
	CDBC2 = dens_Sjo(2)/V_sol !f2
	CDBC4 = dens_Sjo(3)/V_sol !f3
	CDBC3 = dens_Sjo(4)/V_sol !f4
	CDBC5 = dens_Sjo(5)/V_sol !f6
	CDBC6 = dens_Sjo(6)/V_sol !f7
	ZRD74 = dens_Sjo(7)/V_sol !f8

	
      end
