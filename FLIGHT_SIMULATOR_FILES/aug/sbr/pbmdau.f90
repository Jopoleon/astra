!C======================================================================|
subroutine PBMDAU
!C
!C----------------------------------------------------------------------|
	use parameter_inc
	use const_inc
	use status_inc
	use outcmn_inc
      use fenix_params
	use flight_sim_geometrics
	use plasma_state
	
	implicit none

	integer	 J,i,istep,jj
	real*8 exbsh(NA1),chi_edge,mmm
	real*8 dum1,dum2(NA1),nnn,c1,c2,c3,c4,P_LH
	real*8 dum4,dum5,dum6,dum7,sviet
	integer nspec
	parameter (nspec=9)

	real*8 D_j(nspec),enrich(nspec)
	real*8 Dw_j(nspec)
	real*8 Dw2_j(nspec)
	real*8 a_j(nspec)
	real*8 R_j(nspec)
	real*8 V_sol,V_div,V_pump
	real*8 T_pump
	real*8 dens_Sj(nspec)
	real*8 dens_Dj(nspec)
	real*8 dens_Sjo(nspec)
	real*8 dens_Djo(nspec)
	real*8 dens_Djoo(nspec,25000)
	real*8 Gplasma(nspec)
	real*8 delayz,bumb1,bumb2,bumb3
	real*8 bumb4,svrect
	real*8 losspflow
	real*8 Gvmp(nspec),Gvdiv(nspec)
	real*8 Gpell(nspec),Gneuty(nspec)
	real*8 A_1
	real*8 B_1,n_elS,n_elD
	real*8 C_1
	real*8 A_2,dumar1,dumar2
	real*8 B_2,hydspecies,psepdiv
	real*8 C_2,timeyy,veloc
	real*8 c_eros(nspec),w_wallsput,d_sol


	data timeyy /0./
	data j_countm /0/
	save dens_Sjo,dens_Djo
	save dens_Djoo,dumar1,dumar2
	save timeyy,j_countm
	data bumb1 /0./
	data bumb2 /0./
	data bumb3 /0./
	data dumar1 /0./
	data dumar2 /0./
	save bumb1,bumb2,bumb3	

	integer imodel,j_countm
	double precision c_model,dumnn(na1),machz,lparsep



!	write(*,*) 'pbmdau',time,zrd78
!

!equations that are solved per species, for ionized ions:

! cold neutrals at the separatrix: nn_sep = (Gvmp+(1.-R_j)*Dw_j*dens_Sjo)*(1-a_j)/(vel*Surf) * D_j / (Dwj + Dj)
! energy neutrals at the separatrix: nn_sep = (Gvmp+(1.-R_j)*Dw_j*dens_Sjo)*(1-a_j)/(vel*Surf) * Dwj / (Dwj+Dj)

! d nsol / dt = Gvmp*a_j + G_plasma + Dj (n_div / enrich - n_sol) + Dw_j*(-1+aj*(1-Rj)) * n_sol

! d ndiv / dt = Gdiv - Gpump +  + Dj (n_sol - n_div/enrich )  


!	write(*,*) 'cpel1 ',CPEL1
!This model specifically works for:
! 1 D , 2 H, 3 He, 4 B, 5 N,  6 W,  7 Ne   ,  8 Kr,  , 9 Ar
!   F1   F7   F2   F3    F5    F4    F6	       F8      F9
!Definition of constants
!volumes
	V_sol=vsoldiv*GP2*RTOR*GP2* &
     & ABC*(1.+ELON(NA1)**2./2.)**0.5*solwidth
	V_div=(1.-vsoldiv)*GP2*RTOR*GP2*&
     & ABC*(1.+ELON(NA1)**2./2.)**0.5*solwidth

!	V_pump=0.3*V_div
!diffus ioniz and recycl coeff
	machz=sqrt(1.e3*1.6022e-19/1.27e-27/2.d0*TI(na1))
	Lparsep=1./MU(NA1)*GP2*RTOR/2. !	Lparsep=120
	D_j=D_j_m*machz/Lparsep	! time exchange between sol and div

	D_j(8)=5*D_j_m*machz/Lparsep !Kr


	Dw_j=Dw_j_m/(wallpos-geom1d(58))**2.d0	 !time exchange between sol and wall
	Dw2_j=Dw2_j_m/(wallpos-geom1d(58))**2.d0	 !time exchange between sol and wall


	SVIEt	=.0136/TE(na1)
	IF(TE(na1).GT..01)	THEN
	SVIEt=9.7E5*EXP(-SVIEt)*SQRT(SVIEt/(1.+SVIEt))/(SVIEt+.73)
				ELSE
	SVIEt=2.958E5*EXP(-SVIEt)*SQRT(SVIEt)
				ENDIF

	veloc=sqrt(ENCL*9.5756e+07)

	a_j=1. !for all species
	a_j(1)=1.-exp(-0.254*ne(na1)-2.6*te(na1)-1.4*(wallpos-geom1d(58))) !just for D & H & He
	a_j(2)=1.-exp(-0.254*ne(na1)-2.6*te(na1)-1.4*(wallpos-geom1d(58))) !just for D & H & He
	a_j(3)=1.-exp(-0.254*ne(na1)-2.6*te(na1)-1.4*(wallpos-geom1d(58))) !just for D & H & He

	R_j=recycl_wall	!fraction that stays on the wall

	T_pump=CV13

	hydspecies=cdym8+max(0.001,dens_Djo(2)/V_div)
!midplane
	enrich=1.0


!	enrich=V_sol/V_div*(hydspecies)**(2./3.)/16. ! all the same enrichment, to be corrected later on
!	enrich(1)=V_sol/V_div*(hydspecies)**(2./3.)/16. ! all the same enrichment, to be corrected later on
!	enrich(2)=V_sol/V_div*(hydspecies)**(2./3.)/16. ! all the same enrichment, to be corrected later on

!divertor
	enrich(5)=40.
	enrich(8)=10. !Krypton compression factor
	enrich(1)=max(1.,min(20.,0.2*hydspecies**0.67))
	enrich(2)=enrich(1)
!	write(*,*) 'hyd',enrich(1)

	bumb1=0.01
	c_eros(1)=bumb1*2. !D
	c_eros(2)=bumb1*1. !H
	c_eros(3)=bumb1*4. !He
	c_eros(4)=bumb1*10. !B
	c_eros(5)=bumb1*14. !N
	c_eros(6)=bumb1*180. !W
	c_eros(7)=bumb1*28. !Ne
	c_eros(8)=bumb1*84. !Kr
	c_eros(9)=bumb1*50. !Ar
	

!time scheme
 !initial condition
	if (TIME-TSTART.le.TAU.or.plasma_up.eq.0) then
		dens_Sjo(1)=F1(na1)*V_sol	
		dens_Sjo(2)=F7(na1)*V_sol	
		dens_Sjo(3)=F2(na1)*V_sol	
		dens_Sjo(4)=F3(na1)*V_sol	
		dens_Sjo(5)=F5(na1)*V_sol	
		dens_Sjo(6)=F4(na1)*V_sol	!to check
		dens_Sjo(7)=F6(na1)*V_sol	!to check
		dens_Sjo(8)=F8(na1)*V_sol	!to check
		dens_Sjo(9)=F9(na1)*V_sol	!to check
	
		dens_Djo=dens_Sjo
	
		if (IBKDW.lt.0..or.TIME.eq.0.) then
			CDBC1=F1(NA1)
			CDBC5=F7(NA1)
			CDBC2=F2(NA1)
			CDBC3=F3(NA1)
			CDBC6=F5(NA1)
			CDBC4=F4(NA1)
			CDBC7=F6(NA1)
			CDBC8=F8(NA1)
			ZRD74=F9(NA1)
			NNCL=1.e-4
			NNWM=0.
		endif
	
	endif

!	write(*,*) 'b c. pbmdau'
!		write(*,*) CDBC1,CDBC3,NNCL,NNWM
!	write(*,*)'D sd, Dsw, Dws, Dped, nsol,ndiv, ndiv/nsol'
!	write(*,*) D_j(1),Dw_j(1),Dw2_j(1),DN(40),
!     & dens_sjo(1)/V_sol,dens_djo(1)/V_div,
!     & dens_djo(1)/dens_sjo(1)*V_sol/V_div
!	write(*,*) ne(na1),te(na1),wallpos-geom1d(58),a_j(1)
!	write(*,*) 'NEB ',ne(na1),cv4/1e3

	if (plasma_up.eq.0.or.TIME.eq.0.) return

	if (IBKDW.lt.0.) then
!		CDBC1=CDBC1+TAU*(CV4**0.3-CDBC1)
!		CDBC3=0.01
!		NNCL=1.e-4
!		NNWM=0.
!		return
	endif
! 1 D , 2 H, 3 He, 4 B, 5 N,  6 W,  7 Ne   ,  8 Kr,  , 9 Ar
!   F1   F7   F2   F3    F5    F4    F6	       F8      F9
!solve per each particle species
!	do j=1,nspec
	
!equation: dy1/dt = A1 y1 + B1 y2 + C1	
!equation: dy2/dt = A2 y1 + B2 y2 + C2	
!	write(*,*) 'valves ',car33(1:24)


!	CV4=0.
!valve flow from car33 from simulink
	Gvmp(1)=	car33(2)  ! D valve  midplane p/s
	Gvmp(2)=	car33(1)     !H valve midplane p/s
	Gvmp(3)=	car33(4) !CDJM5 ! He
	Gvmp(4)=	car33(20)+boron_wall ! B
	Gvmp(5)=	car33(3) ! N
	Gvmp(6)=	zrd72+1*car33(21)+3.*0.0005*sum(TE(1:NA1))**1.* & !zrd72 from nbispl demo emulation w control for qrad
     & exp(-2.2+(RTOR+ABC+SHIFT)) ! W
	Gvmp(7)=	car33(6) !CDWM5 ! Ne
	Gvmp(8)=	car33(8) 
	!or Psep EFABLECOMMENT
!	if (TIME.lt.0.9+5.4) then
!		psepdiv=(TIME-5.4)/0.9*1
!	endif
!	if (TIME.ge.0.9+5.4) psepdiv=1.
!	if (TIME.lt.3.7+5.4.and.TIME.ge.2.3+5.4) then
!	psepdiv=3.-(TIME-2.3-5.4)/1.4*2.
!	endif
!	if (TIME.ge.3.7+5.4) psepdiv=1.
!	psepdiv=1
!	dumar1=max(0.,dumar1+0.03*(CRAD3-QRADR(ROC)-psepdiv))
!	Gvmp(9)=	max(0.,dumar1) ! Ar
!	write(*,*) 'Ar puff ',CRAD3,QRADR(ROC),psepdiv,dumar1
	Gvmp(9)=car33(7)
	
!tungsten erosion
	bumb3=predep_w !prompt redeposition factor
	bumb1=sqrt(max(0.,CDWM1-5.))
	bumb2=sum(dens_djo*c_eros)
	CSCL2=bumb2*bumb1*(1.-bumb3)

	Gvdiv(1)=	car33(11) !CV4 !D
	Gvdiv(2)=	car33(10) !ZRD73  ! H
	Gvdiv(3)=	car33(13)  ! He
	Gvdiv(4)=	car33(23)  !B
	Gvdiv(5)=	car33(12) !CDWM6  !N puff
	Gvdiv(6)=	car33(24)  !W puff due to erosion, computed above
	Gvdiv(7)=	car33(15)  !Ne puff 
	Gvdiv(8)=	car33(17)  !Kr puff 
	Gvdiv(9)=	car33(16)  !Ar puff 

!	write(*,*) 'valves ',gvmp,gvdiv

	Gpell = 0.

!neutral plasma flow
	Gneuty=0.
	j=1
	Gneuty(j)=(1.-a_j(j))*(Gvmp(j)+ &
     & (1.-R_j(j))*Dw2_j(j)*dens_Sjo(j))	
	
!plasma flow
!	Gplasma(1)=-(sum(f1(1:na1)*vr(1:na1))*hro-
!     & sum(f1o(1:na1)*vr(1:na1))*hro)/tau
!	Gplasma(2)=-(sum(f2(1:na1)*vr(1:na1))*hro-
!     & sum(f2o(1:na1)*vr(1:na1))*hro)/tau
!	Gplasma(3)=-(sum(f4(1:na1)*vr(1:na1))*hro-
!     & sum(f4o(1:na1)*vr(1:na1))*hro)/tau
!	Gplasma(4)=-(sum(f3(1:na1)*vr(1:na1))*hro-
!     & sum(f3o(1:na1)*vr(1:na1))*hro)/tau
!	Gplasma(5)=-(sum(f6(1:na1)*vr(1:na1))*hro-
!     & sum(f6o(1:na1)*vr(1:na1))*hro)/tau
!	Gplasma(6)=-(sum(f7(1:na1)*vr(1:na1))*hro-
!     & sum(f7o(1:na1)*vr(1:na1))*hro)/tau
!	Gplasma(7)=-(sum(f8(1:na1)*vr(1:na1))*hro-
!     & sum(f8o(1:na1)*vr(1:na1))*hro)/tau
!	Gplasma(8)=-(sum(f8(1:na1)*vr(1:na1))*hro-
!     & sum(f8o(1:na1)*vr(1:na1))*hro)/tau
!	Gplasma(9)=-(sum(f8(1:na1)*vr(1:na1))*hro-
!     & sum(f8o(1:na1)*vr(1:na1))*hro)/tau


!	do j=1,nspec
	j=1
!		Gplasma(j)=Gplasma(j)+Gpell(j)+Gneuty(j)+snebtr(roc)

	  Gplasma(1)=QF1(na1)
	  Gplasma(2)=QF7(na1)
	  Gplasma(3)=QF2(na1)
	  Gplasma(4)=QF3(na1)
	  Gplasma(5)=QF5(na1)
	  Gplasma(6)=QF4(na1)
	  Gplasma(7)=QF6(na1)
	  Gplasma(8)=QF8(na1)
	  Gplasma(9)=QF9(na1)

!	write(*,*) 'gplasma ',Gplasma(1)
!		Gplasma(j)=Gplasma(j)+Gpell(j)+Gneuty(j)
!	enddo

	j=1
	NNWM = (1.-a_j(j))*(1.-R_j(j))*Dw2_j(j)*dens_Sjo(j)
	NNCL = (1.-a_j(j))*Gvmp(j)
	j=2
	NNWM = NNWM+(1.-a_j(j))*(1.-R_j(j))*Dw2_j(j)*dens_Sjo(j)
	NNCL = NNCL+(1.-a_j(j))*Gvmp(j)
	j=3
	NNWM = NNWM+(1.-a_j(j))*(1.-R_j(j))*Dw2_j(j)*dens_Sjo(j)
	NNCL = NNCL+(1.-a_j(j))*Gvmp(j)
!	NNWM = 1.e1
!	NNCL = 0.

!	NNWM=0.
!	ENCL= 0.01
!	ENWM= 10.

	veloc=sqrt(ENCL*9.5756e+07)
	NNCL=NNCL/veloc/SLAT(NA1)
	veloc=sqrt(ENWM*9.5756e+07)
	NNWM=NNWM/veloc/SLAT(NA1)

!	write(*,*) 'nnwm,nncl ',NNWM,NNCL
	
	NNWM=min(3.e-4,NNWM) !test
	NNCL=min(3.e-4,NNCL) !test

	veloc=sqrt(ENCL*9.5756e+07)

!calculate NN
!	NN=NNCL*exp((xrho-1.)*ne*te**2.)	
!	if (TIME.eq.0.) NN=0.
!	NN(NA1)=NNCL+NNWM
!	dumnn=nn(1:na1)
!	do jj=2,na1
!	j=na1-jj
!	include 'fml/svie'
!	include 'fml/svii'
!	include 'fml/svrec'
!	NN(j)=max(0.,dumnn(j)+TAU*(
!     & -SVIE*dumnn(j)*ne(j)+
!     & SVREC*NI(j)*NE(j)- 
!     & 1./VR(j)/HRO*(
!     & (VRS(j)*dumnn(j)*veloc-VRS(j-1)*dumnn(j-1)*veloc)))) 
!	write(*,*) 'neut ',NN(j),j,SVIE*dumnn(j)*ne(j),SVREC*NI(j)*NE(j)
!	write(*,*) 'neut ',1./VR(j)/HRO*(
!     & (VRS(j+1)*dumnn(j+1)*veloc-VRS(j)*dumnn(j)*veloc))
!	enddo
!	write(*,*) 'neut ',CSOL4,veloc
!	write(*,*) 'neut ',NNWM,NNCL,veloc*slat(na1),gneuty(1)

!	ZRD30=Gneuty(1)/(Gneuty(1))

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
		A_1=1/tau-Dw2_j(j)*(a_j(j)*(1.-r_j(j)))+D_j(j)+Dw_j(j)
		B_1=-D_j(j)/enrich(j)
		C_1=Gplasma(j)+a_j(j)*Gvmp(j)+dens_Sjo(j)/tau

		A_2=-D_j(j)
		B_2=1/tau+D_j(j)/enrich(j)+T_pump/V_div
		C_2=Gvdiv(j)+dens_Djo(j)/tau


	dum7=1./(A_1*B_2-A_2*B_1)

	dens_Sj(j)=dum7*(B_2*C_1-B_1*C_2)			

	dens_Dj(j)=dum7*(-A_2*C_1+A_1*C_2)	

	enddo	
	

!assign final
	CHE3=(dens_Sj(6)-dens_Sjo(6))/tau ! W influx into plasma
	dens_Sjo=dens_Sj
	dens_Djo=dens_Dj


!assign back to boundary conditions
	CDBC1 = max(0.001,dens_Sjo(1)/V_sol) !d sol
	CDYM8 = max(0.001,dens_Djo(1)/v_div) !d div

	CDBC5 = max(0.001,dens_Sjo(2)/V_sol) !h

	CDBC2 = max(0.000001,dens_Sjo(3)/V_sol) !he

	CDBC3 = max(0.000001,dens_Sjo(4)/V_sol) !b

	CDBC6 = max(0.000001,dens_Sjo(5)/V_sol) !n
	CDYM6 = max(0.000001,dens_Djo(5)/V_div) !n

	CDBC4 = max(0.000001,dens_Sjo(6)/V_sol) !w

	CDBC7 = max(0.000001,dens_Sjo(7)/V_sol) !ne

	CDBC8 = max(0.000001,dens_Sjo(8)/V_sol) !kr

	ZRD74 = max(0.000001,dens_Sjo(9)/V_sol) !ar

	CDJM1= max(0.000001,dens_djo(6)/V_div) !w
	CDJM2= max(0.000001,dens_Sjo(6)/V_sol) !w

!		write(*,*) CDBC1,CDBC2,CDBC3,CDBC4,CDBC5,CDBC6,CDBC7
!		write(*,*) CDBC8,ZRD74
!		write(*,*) ' ' 

!	CDYM5 = dens_Sjo(4)/v_sol
!	CDYM6 = dens_Sjo(5)/v_sol
!	CDYM7 = dens_Djo(4)/v_div
!	ZRD75 = dens_Sjo(7)/v_sol
!	ZRD76 = dens_Djo(7)/v_div

!valve flow from car33 from simulink
!	Gvmp(1)=	car33(2)  ! D valve  midplane p/s
!	Gvmp(2)=	car33(1)     !H valve midplane p/s
!	Gvmp(3)=	car33(4) !CDJM5 ! He
!	Gvmp(4)=	car33(11)+boron_wall ! B
!	Gvmp(5)=	car33(3) ! N
!	Gvmp(6)=	car33(12)+0.00015*sum(TE(1:NA1))**2.*
!     & exp(-2.2+(RTOR+ABC+SHIFT)) ! W
!	Gvmp(7)=	car33(6) !CDWM5 ! Ne
!	Gvmp(8)=	car33(8) 
!	Gvmp(9)=car33(7)

	car34(1:12)=0
	car34(2) = dens_Djo(1)*V_pump/V_div**2. &
     & *T_pump   ! D pumped in p/s
	car34(1) = dens_Djo(2)*V_pump/V_div**2. &
     & *T_pump   ! D pumped in p/s
	car34(4) = dens_Djo(3)*V_pump/V_div**2. &
     & *T_pump   ! D pumped in p/s
	car34(11) = dens_Djo(4)*V_pump/V_div**2. &
     & *T_pump   ! D pumped in p/s
	car34(3) = dens_Djo(5)*V_pump/V_div**2. &
     & *T_pump   ! D pumped in p/s
	car34(12) = dens_Djo(6)*V_pump/V_div**2. &
     & *T_pump   ! D pumped in p/s
	car34(6) = dens_Djo(7)*V_pump/V_div**2. &
     & *T_pump   ! D pumped in p/s
	car34(8) = dens_Djo(8)*V_pump/V_div**2. &
     & *T_pump   ! D pumped in p/s
	car34(7) = dens_Djo(9)*V_pump/V_div**2. &
     & *T_pump   ! D pumped in p/s


	car34(13:24) = 0. !neutral pressures in divertor of gases
	CDJM4=sum(car34(13:24)) !total neutral pressure in divertor of gases
	
      end
