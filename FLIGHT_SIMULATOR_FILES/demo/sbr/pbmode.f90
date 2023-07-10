!C======================================================================|
	subroutine PBMODE(machinez)
!C
!C----------------------------------------------------------------------|
	use parameter_inc
	use const_inc
	use status_inc
	implicit none
!	include	'for/parameter.inc'
!	include 'for/const.inc'
!	include 'for/status.inc'
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
! 1 D , 2 T, 3 He, 4 imp1core, 5 imp2div, 6 imp3unwanted, 7 Impsol_Kr
	
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
	enrich(4)=6.d0
	enrich(5)=20.d0
	enrich(6)=6.d0
	enrich(7)=6.d0
!	enrich(7)=1.d0

	bumb1=0.01
	c_eros(1)=bumb1*2. !D
	c_eros(2)=bumb1*2. !T
	c_eros(3)=bumb1*4. !He
	c_eros(4)=bumb1*140. !Xe
	c_eros(5)=bumb1*40. !Ar
	c_eros(6)=bumb1*180. !W
!	c_eros(7)=bumb1*84. !Kr
	c_eros(7)=bumb1*16. !Oxi
	

!time scheme
 !initial condition
	if (TIME.le.TAU) then
		dens_Sjo(1)=F1(na1)*V_sol	
		dens_Sjo(2)=F2(na1)*V_sol	
		dens_Sjo(3)=F4(na1)*V_sol	
		dens_Sjo(4)=F3(na1)*V_sol	
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






! 1 D , 2 T, 3 He, 4 imp1core, 5 imp2sol, 6 imp3unwanted, 7 impsolkr
!solve per each particle species
!	do j=1,nspec
	
!equation: dy1/dt = A1 y1 + B1 y2 + C1	
!equation: dy2/dt = A2 y1 + B2 y2 + C2	

!	CV4=0.
!valve flow
!	Gvmp(1)=	CDYM3*CV4 + CDWM6*CDWM5   ! D valve  midplane p/s 
	Gvmp(1)=	(1.-CDYM3)*CV4 + CDWM6*CDWM5+0.*2.*ZRD73   ! D valve  midplane p/s + 2*water
	Gvmp(2)=	CDYM3*CV4 + (1.-CDWM6)*CDWM5    !T valve midplane p/s
	Gvmp(3)=	0.*CDJM5 !reinjected He by mistake  ! He
	Gvmp(4)=	CBND3+0.*CDJM6 !reinjected Xe by mistake !Xe
	Gvmp(5)=	0.*CDJM7 !reinjected Ar by mistake !Ar
	Gvmp(6)=	0.08; !+CDJM8 !reinjected W by mistake !W!+ wall source 0.1*teped
!	Gvmp(7)=	ZRD70+ZRD73 !reinjected Kr by mistake + Kr puff mp
	Gvmp(7)=	0.*(ZRD70+ZRD73) !reinjected Kr by mistake + Oxi puff mp


!tungsten erosion
	bumb3=0.95 !prompt redeposition factor
	bumb1=sqrt(max(0.,CDWM1-5.))
	bumb2=sum(dens_djo*c_eros)
	CSCL2=0.01*bumb2*bumb1*(1.-bumb3)
!	write(*,*) TIME,CDWM1,CSCL2

	Gvdiv(1)=	0.
	Gvdiv(2)=	0.
	Gvdiv(3)=	0.
	Gvdiv(4)=	0.
	Gvdiv(5)=	CSCL1  !Ar puff
	Gvdiv(6)=	CSCL2  ! W puff due to erosion, computed above
	Gvdiv(7)=	ZRD71  !Kr puff 

!pellet flow + fusion products
	Gpell(1)=	CPEL1-1./5.632*QDTR(ROC)
	Gpell(2)=	CIMP3-1./5.632*QDTR(ROC)
	Gpell(3)=	1./5.632*QDTR(ROC)
	Gpell(4)=	CPEL1*CDJM9  !Xe pellet
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
	Gplasma(3)=QF4(na1)
	Gplasma(4)=QF3(Na1)
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
	

!	write(*,*) ' '
!	write(*,*) time
!	write(*,*) ' '
!	write(*,*) A_1
!	write(*,*) ' '
!	write(*,*) B_1
!	write(*,*) ' '
!	write(*,*) C_1
!	write(*,*) ' '
!	write(*,*) ' '
!	write(*,*) A_2
!	write(*,*) ' '
!	write(*,*) B_2
!	write(*,*) ' '
!	write(*,*) C_2
!	write(*,*) ' '
!	write(*,*) Gplasma
!	write(*,*) ' '
!	write(*,*) dens_Sj/v_sol
!	write(*,*) ' '
!	write(*,*) dens_Dj/v_div
	
!	enddo


!assign final
	dens_Sjo=dens_Sj
	dens_Djo=dens_Dj


!assign back to boundary conditions
	CDBC1 = dens_Sjo(1)/V_sol !f1
	CDBC2 = dens_Sjo(2)/V_sol !f2
	CDBC3 = dens_Sjo(3)/V_sol !f4
	CDBC4 = dens_Sjo(4)/V_sol !f3
	CDBC5 = dens_Sjo(5)/V_sol !f6
	CDBC6 = dens_Sjo(6)/V_sol !f7
	ZRD74 = dens_Sjo(7)/V_sol !f8

   
	CDYM5 = dens_Sjo(4)/v_sol
	CDYM6 = dens_Sjo(5)/v_sol
	CDYM7 = dens_Djo(4)/v_div
	CDYM8 = dens_Djo(5)/v_div
	ZRD75 = dens_Sjo(7)/v_sol
	ZRD76 = dens_Djo(7)/v_div

!	write(*,*) CDBC1,CDBC2,CDBC3,CDBC4,CDBC5,CDBC6,ZRD74
!	write(*,*) CDYM6,CDYM8,CDYM5,CDYM7

!valid only for equidistant time step
!	if (time.gt.delayz) then
!		do j=1,2
!			dens_Djoo(j,1:j_countm-1)=dens_Djoo(j,2:j_countm)
!			dens_Djoo(j,j_countm)=dens_Djo(j)
!		enddo
!	endif
!	write(*,*) ne(1:na1)

!	write(9991,'(111E25.11)') TIME,
!     & CDBC1,CDBC2,CDBC3,CDBC4,CDBC5,CDBC6,
!     & dens_Sjo/v_sol,dens_Djo/v_div,
!     & CDYM5,CDYM6,CDYM7,CDYM8

	CDWM2 = dens_Djo(1)*V_pump/V_div**2. &
     & *T_pump   ! D pumped in p/s
	CDWM7 = dens_Djo(2)*V_pump/V_div**2. &
     & *T_pump   ! T pumped in p/s

	CDJM1 = dens_Djo(3)*V_pump/V_div**2. &
     & *T_pump   ! He pumped in p/s

	CDJM2 = dens_Djo(4)*V_pump/V_div**2. &
     & *T_pump   ! Xe pumped in p/s

	CDJM3 = dens_Djo(5)*V_pump/V_div**2. &
     & *T_pump   ! Ar pumped in p/s

	CDJM4 = dens_Djo(6)*V_pump/V_div**2. &
     & *T_pump   ! W pumped in p/s
					
	ZRD77 = dens_Djo(7)*V_pump/V_div**2. &
     & *T_pump   ! Kr pumped in p/s

	ZRD25=(CDJM1/CDWM2) ! He/D conc in pumping
	ZRD26=(CDJM2/CDWM2) ! Xe/D conc in pumping
	ZRD27=(CDJM3/CDWM2) ! Ar/D conc in pumping

!	write(*,*) 'fases',
!     & CDBC3,CDBC4,CDBC5,CDBC6,ZRD74






	if (TIME.gt.100) then
	
	jgipp=jgipp+1
	pfavg=pfavg+qdtr(roc)*5
	sp1avg=sp1avg+zrd25/2
	sp2avg=sp2avg+zrd26/2
	sp3avg=sp3avg+zrd27/2
	
	endif
	
	
	write(*,*) 'gas ',CDBC1,CDBC2, &
     & sp1avg/jgipp,sp2avg/jgipp,sp3avg/jgipp,pfavg/jgipp

!      1-purity  He       Xe          Ar 
!	fac: 0.    ,   5.95%, 4.51 promil, 2.57 %, pfus=1971 MW
!      0.015 ,  6.% , 4.05 promil, 2.51%, Pfus = 1937 MW
!      0.025,   6%,  3.7 promil,  2.53% , pfus = 1901 
!      0.05 , 5.97% ,  2.8 permil , 2.69%    . Pfus = 1804 
!      0.075 ,  5.96% , 2.153 promil, 2.48%, Pfus = 1759 MW
!      0.1 ,  6.% , 1.65 promil, 2.25%, Pfus = 1710 MW

!	write(1331,'(111E25.11)') time,dens_Djo,dens_Sjo
!global particle balance
!	bumb1=sum(f1(1:na1)*vr(1:na1))*hro
!     & +dens_Sjo(1)+dens_Djo(1)
!	bumb2=Gpell(1)+Gvmp(1)+Gvdiv(1)	
!	write(*,*) 'D lifetime : ',bumb1/bumb2

!	bumb1=sum(f2(1:na1)*vr(1:na1))*hro
!     & +dens_Sjo(2)+dens_Djo(2)
!	bumb2=Gpell(2)+Gvmp(2)+Gvdiv(2)	
!	write(*,*) 'T lifetime : ',bumb1/bumb2


!	bumb1=sum(f4(1:na1)*vr(1:na1))*hro
!     & +dens_Sjo(3)+dens_Djo(3)
!	bumb2=Gpell(3)+Gvmp(3)+Gvdiv(3)	
!	write(*,*) 'He lifetime : ',bumb1/bumb2
!	write(*,*) 'He taustar : ',bumb1/bumb2/TAUER(ROC)

!	bumb1=sum(f3(1:na1)*vr(1:na1))*hro
!     & +dens_Sjo(4)+dens_Djo(4)
!	bumb2=Gpell(4)+Gvmp(4)+Gvdiv(4)	
!	write(*,*) 'Xe lifetime : ',bumb1/bumb2

!	bumb1=sum(f6(1:na1)*vr(1:na1))*hro
!     & +dens_Sjo(5)+dens_Djo(5)
!	bumb2=Gpell(5)+Gvmp(5)+Gvdiv(5)	
!	write(*,*) 'Ar lifetime : ',bumb1/bumb2

!	bumb1=sum(f7(1:na1)*vr(1:na1))*hro
!     & +dens_Sjo(6)+dens_Djo(6)
!	bumb2=Gpell(6)+Gvmp(6)+Gvdiv(6)	
!	write(*,*) 'W lifetime : ',bumb1/bumb2

! partilcle inventory D
!	bumb1=bumb1+tau*(Gpell(1) + 
!     & 1./5.632*QDTR(ROC) + CDYM3*CV4) !fresh D request p / s
!	bumb2=bumb2+tau*(Gpell(2) + 
!     & 1./5.632*QDTR(ROC) + (1.-CDYM3)*CV4) !fresh T request p / s
!T available

!	bumb4=QDTR(ROC)*4./5.
!	bumb3=bumb3+tau*(-(Gpell(2) + 
!     & 1./5.632*QDTR(ROC) + (1.-CDYM3)*CV4)+
!     & bumb4+CDWM7-(1.-CDWM6)*CDWM5)

!	write(*,*) 'D needed,T needed, T inventory',bumb1,bumb2,bumb3

	

















	
      end
