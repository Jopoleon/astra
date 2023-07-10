	subroutine SOLMODINP(zimp1,zimp2,zimp3,qplateyy,&
     & T_sep)

	use parameter_inc
	use const_inc
	use status_inc
	use outcmn_inc

	
	!uses pbmodel
		implicit none
!	include	'for/parameter.inc'
!	include	'for/const.inc'
!	include	'for/status.inc'
!	include	'for/outcmn.inc'

	double precision cz_sol,cz_div,Lparsep,qpar,fx,gamma_up 
	double precision cz_sol2,cz_div2
	double precision cz_sol3,cz_div3,machinez
	double precision lambda_p,Dr,n_up,q_plate,Tup_new,Tsh_new,fm_SE,Lr
	double precision lambda_q,zimp1,zimp2,zimp3, &
     & cz_core1,cz_core2,cz_core3,qplatey,qplateyy,dummy1

	double precision n_sep, T_sep,tg1,tau_relax,tg2
	double precision qeiadr,qtokr,qedwtr,qidwtr,qntotr
	double precision qetotr,qitotr

	double precision psepe,psepi,dum12,comp_facsepsol, comp_facsoldiv
	double precision comp_facsepsol2, comp_facsoldiv2
	double precision comp_facsepsol3, comp_facsoldiv3
	double precision tg_u,fm_g,dqoq_g,ldiv,t_divg,bbolp,t_tung
	integer i_call,itypez,j1,j2,j3,j4
	data i_call /0/
	data qplatey /10./
	data lambda_q /0.01/
	data t_tung /273./
	data dum12/1./
	save i_call,qplatey
	save tg_u,fm_g,dqoq_g	,t_divg
	save lambda_q
	save t_tung,dum12
	
!q_up = P_SOL/(2 *pi*R*lambda_q),  dove R e' il raggio maggiore dell'outer midplane e lambda_q la e-folding length (per DEMO e' dalle parti del mm secondo il famoso Eich scaling, ma ci si puo' giocare). La routine vuole q_up in W/m2.

!Per gamma_up invece tu hai particelle/secondo. Credo che la cosa piu' sensata sia usare particelle/secondo e dividerle per la superficie che corrisponde alla last close flux surface (tipo 2*pi*R * 2*pi*r). Cioe' suppongo che le particelle escano dalla separatrice in modo omogeneo. Poi magari il conto si puo' raffinare...
!	write(*,*) 'solmod'


	if (TIME.lt.2*TAU) then
	else
		
!	write(*,*) 'sol model 2'
	dummy1=0.73e-3*BTOR**(-0.78)*&
     &   (abc**2.*btor/(0.2*rtor*ipl))**1.02*&
     &   (QTOKR(ROC)-QEDWTR(ROC)-QIDWTR(ROC))**0.1*&
     &   (RTOR)**0.02*&
     &    3.
	dummy1=dummy1*3.  ! ~ 2 cm    , nominal *1 is 1 mm for ITER, 1.3 mm for DEMO1, 3* can be feaasible

	if (isnan(dummy1)) then
	else
	lambda_q=dummy1   !gives 1 mm for ITER 15 MA scenario
	endif					!put Eich scaling
!	cpel3=lambda_q			
!	write(*,*) 'sol model 3'
		
	cz_sol=comp_facsepsol*cz_core1 ! Impurity concentration - Conductive Region !Xe
	cz_div=comp_facsoldiv*cz_sol ! Impurity concentration - Convective Region

	cz_sol2=comp_facsepsol2*cz_core2 ! Impurity concentration - Conductive Region
	cz_div2=comp_facsoldiv2*cz_sol2 ! Impurity concentration - Convective Region

	cz_sol3=comp_facsepsol3*cz_core3 ! Impurity concentration - Conductive Region
	cz_div3=comp_facsoldiv3*cz_sol3 ! Impurity concentration - Convective Region

	Lparsep=1./MU(NA1)*GP2*RTOR/2. !	Lparsep=120
	ldiv=lparsep*0.55 !50 m Ldiv in DEMO

!	qpar=max(1.,(QTOKR(ROC)-QEDWTR(ROC)-QIDWTR(ROC))*1e6/
!     & (GP2*(RTOR+ABC+SHIFT)*lambda_q*ABC/RTOR*MU(NA1))/2.)
!	write(*,*) 'sol model 4'

! better calculation with cylindrical q:
	qpar=max(1.,(QTOKR(ROC)-QEDWTR(ROC)-QIDWTR(ROC))*1e6/ &
     & (GP2*(RTOR+ABC+SHIFT)*lambda_q*IPL/(GP*ABC*BTOR))/2.)

	psepe=QETOTR(ROC)-QEIADR(ROC)-QEDWTR(ROC)
	psepi=QITOTR(ROC)+QEIADR(ROC)-QIDWTR(ROC)


!	write(*,*) QTOKR(ROC),QEDWTR(ROC),QIDWTR(ROC)
!	write(*,*) lparsep,ldiv,lambda_q,rtor

!	write(*,*) 'dummy1',dummy1,lambda_q,qpar,
!     &  qtokr(roc),QEDWTR(ROC),QIDWTR(ROC)

!	Dr = QTOKR(ROC)-QEDWTR(ROC)-QIDWTR(ROC)
!	if (Dr.lt.140) pause
	
!	write(*,*) 'fdjsk : ',QTOKR(ROC),qpar/1.e6,zimp1,zimp2,lambda_q

!Choose impurities
! 1 = H, 2 = He, 3 = Be, 4 = C, 5 = N, 6 = O, 7 =  Ne, 8 = Si, 9 = Ar, 10 = Fe, 11 = Ni, 12 = Kr, 13 = Xe , 14 = W
!		write(*,*) TIME
!	write(*,*) 'sol model 5'
	fx=10
		gamma_up=1.e19*(QNTOTR(ROC))/SLAT(NA1)
	lambda_p=0.01 ! [m] Decay length for the radial density profile @ OM
	Dr=0.1 ! [m2/sec] Radial diffusion coefficient
		n_sep=NE(NA1)
		T_sep = TE(NA1)

	gamma_up = 1.e19*NE(NA1)*Dr/lambda_p
	tg1=1.e3*TE(NA1)
!	tg2=Tsh_new
	if (TIME.le.TAU) tg2=10.
	tau_relax=0.1
	if (TIME.ge.TAU) then
	n_up=NE(na1)*1.e19
	if (i_call.eq.0) then
	tg2=30.
	tg_u=100.
	fm_g=0.
	dqoq_g=0.
	i_call=1
	t_divg=40.
	endif
	dqoq_g=min(0.99,dqoq_g)


!	write(*,*) 'sol model 6'

!	j4=10
!	do j3=1,j4
!	do j2=1,j4
!	do j1=1,j4
!	NE(NA1)=2.+(j3-1)/(j4-1.)*6.
!	NE(NA1)=4.4
!	NE(NA1)=4.
!	n_up=NE(na1)*1.e19
!	if (j2.le.20) then
!	CDYM6=(j2-1.)/(j4-1.)*0.25
!	else
!	CDYM6=0.25-(j2-21.)/(j4-1.)*0.25
!	endif
!	CDYM6=0.085
!	CDYM8=CDYM6*4.
!	CDYM5=0.
!	CDYM7=0.
!	gamma_up = 1.e19*NE(NA1)*Dr/lambda_p
!	qpar=(j1+1.e-16)/(j4+1.e-16)*300.
!	qpar=150.
!	qpar=qpar*1.e6/
!     & (GP2*(RTOR+ABC+SHIFT)*lambda_q*ABC/RTOR*MU(NA1))/2.
!	if (j1*j2.eq.1) then
	tg2=700.
	tg_u=1000.
	fm_g=0.
	dqoq_g=0.
	t_divg=800.
!	endif
	tg2=1.
	tg_u=100.
	fm_g=1.
	dqoq_g=1.
	t_divg=10.


!	write(*,*) 'sol model 7',cz_sol,cz_div,Lparsep,qpar,fx,gamma_up,
!     & lambda_p,Dr,n_up,q_plate,Tup_new,
!     & Tsh_new,fm_SE,Lr,zimp1,zimp2,zimp3,
!     & cz_sol2,cz_div2,cz_sol3,cz_div3,tg1,
!     & tau_relax,TAU,tg2,tg_u,fm_g,dqoq_g,Ldiv,t_divg,RTOR,
!!     & QTOKR(ROC),CDYM6,CDYM8,CDYM5,CDYM7,ZRD75,ZRD76


!	call MS_SOL_model_sub(cz_sol,cz_div,Lparsep,qpar,fx,gamma_up,
!     & lambda_p,Dr,n_up,q_plate,Tup_new,
!     & Tsh_new,fm_SE,Lr,zimp1,zimp2,zimp3,
!     & cz_sol2,cz_div2,cz_sol3,cz_div3,tg1,
!     & tau_relax,TAU,tg2,tg_u,fm_g,dqoq_g,Ldiv,t_divg,RTOR,
!     & QTOKR(ROC),CDYM6,CDYM8,CDYM5,CDYM7,ZRD75,ZRD76)

!	write(*,*) 'sol model 8'
	
!	write(*,*) 'sol model info'
!	write(*,*) 'qpar,q_plate',qpar/1e6,q_plate/1e6
!!	write(125,'(111E25.11)') TIME,qpar/1e6,q_plate/1e6,
!     & ne(na1),CDYM6,CDYM8,CDYM5,CDYM7

!for testing fast other things
	if (nint(zrd82).eq.0) then
	tg2=max(0.01,qpar)/1.e7*atan(1./(CDYM8/ &
     & (max(0.01,qpar)/1.e6*0.00033*(4./NE(NA1))))**16.)
	q_plate = tg2*1.e6*2.d0/4.5d0
	tg_u=max(170.,tg2)
	else
	tg2=max(0.01,qpar)/1.e7*atan(1./(CDYM8/ &
     & (max(0.01,qpar)/1.e6*0.00011*(4./NE(NA1))))**16.)
	q_plate = tg2*1.e6*2.d0/4.5d0
	tg_u=max(170.,tg2)
	endif
!
	if (tg2.le.5.) zrd82=1.d0
	if (tg2.gt.5.) zrd82=0.d0

!	write(126,'(111E25.11)') TIME,qpar/1e6,q_plate/1e6,
!     & ne(na1),CDYM6,CDYM8,CDYM5,CDYM7



!	write(*,*) TIME,q_plate/1.d6,tg2,nint(ZRD82)

!	write(1551,'(222E25.11)') NE(NA1),CDYM8,qpar,
!     & q_plate,tg2,tg_u 	


	if (isnan(q_plate)) then
	qplatey=qplatey
	else
		qplatey=q_plate/1.e6
	endif
	qplateyy=qplatey
		CDWM1=tg2  ! T plate
		n_sep=n_up/1e19
		T_sep = min(0.3,tg_u/1.e3)
!	endif
	endif
	if (isnan(T_sep).or.T_sep.lt.1.e-3) then
	 i_call=0
!	write(*,*) 'crap'
	endif
	if (qplatey.gt.1.e6) then
	 i_call=0
!	write(*,*) 'crap'
	endif
!	write(*,*) ' t sep'
!	write(*,*) T_sep

	dum12=max(1.,dum12+0.01*(psepi-psepe))

	T_sep=TE(NA1)*dum12
	write(*,*) 'ti sep',T_sep,psepi,psepe,dum12
!	write(1237,'(111E25.11)') TIME,qpar,CDYM6,NE(NA1),
!     & q_plate,CDYM8,F6(1),F6(1)/NE(1)

!	write(6751,'(13E25.11)') TIME,Tup_new,cz_core1,cz_core2,
!     & cz_sol,cz_sol2,cz_div,cz_div2,
!     & Tsh_new,qpar/1.e6,QDTR(ROC),QTOKR(ROC),CV7
!	bbolp=1.3
!	CSCL3=max(0.,CSCL3+min(bbolp,tau)*1.e-4*
!     & (CV8-5.)/(1.+min(bbolp,tau))) 
!	cscl3=cscl3*zrd42x
	endif

	if (TIME.le.TAU) T_sep =0.05


	end
