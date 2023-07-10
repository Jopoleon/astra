!C======================================================================|
	subroutine TRMODEL(model_type,smear_c,machinez)
!C
!C----------------------------------------------------------------------|

	use parameter_inc
	use const_inc
	use status_inc
	use outcmn_inc
	use fenix_params

	implicit none

!	include	'for/parameter.inc'
!	include 'for/const.inc'
!	include 'for/status.inc'
!	include 'for/outcmn.inc'

	integer	 JJ,J,i,istep
	integer	jk
	real*8 model_type,machinez,tpfft,grne
	real*8 q_pos,p_pos,yhatli,smear_c
	real*8 d_1(nrd),ngree,hatly,lte,lti
	real*8 d_2(nrd),grnee(nrd),hreg
	real*8 d_3(nrd),d_5(nrd)
	real*8 dum5(nrd),dumhh(5),dummh(na1,5)
	real*8 d_4(nrd),delayz,timnow
	real*8 dum1(nrd),wtotr,qtokr,qradr,lne,coulg,nuee,nues
	real*8 dum2(nrd),qedwtr,ftllmr,tpf
	real*8 dum3(nrd),wped,qidwtr,ftav
	real*8 dum4(nrd)
	real*8 dum6(nrd),bumb(na1,5)
	real*8 dum7(nrd)

!	include 'dat/fenixparams.var'

	real*8 dum8(nrd),chi_scal
	real*8 dum9(nrd)
	real*8 dum10(nrd)
	real*8 C_T(nrd)
	real*8 C_P(nrd)
	real*8 C_NU(nrd)
	real*8 calib_f,rnd_nmb
	real*8 nueih
	JJ=nint(model_type)
	
	!	subroutine	SMEARR(ALFA,FO,FN)
!	write(*,*) JJ

!	d_1(1:na1)=min(20.,max(0.001,abs(work(1:na1,21))))
!	d_2(1:na1)=min(20.,max(0.001,abs(work(1:na1,22))))
!	d_3(1:na1)=min(20.,max(0.001,abs(work(1:na1,23))))
!	d_4(1:na1)=min(1.,max(-1.,work(1:na1,24)))
!	d_1(90:94)=d_1(90)
!	d_2(90:94)=d_2(90)
!	d_3(90:94)=d_3(90)
!	d_4(90:94)=d_4(90)
!	call smearr(smear_c,d_1,car41)
!	call smearr(smear_c,d_2,car42)
!	call smearr(smear_c,d_3,car43)
!	call smearr(smear_c,d_4,car44)

!	return	
	
	if (TIME.le.0.) then
	work(1:na1,21)=0.3+xrho(1:na1)**2.	
	work(1:na1,22)=0.3+xrho(1:na1)**2.	
	work(1:na1,23)=0.3+xrho(1:na1)**2.	
	work(1:na1,24)=0.*(0.3+xrho(1:na1)**2.)	

	return	
	endif	
	transp_variance=0.05  ! 5% transport coefficient variance and randomly generated at each time step
		
! hatli model
	if (JJ.eq.1) 		then
	
		do i=1,na1
	J=i
			include 'fml/hatly'
!			write(*,*) hatly
!	call random_number(rnd_nmb)
!	write(*,*) 'rnd',rnd_nmb
	work(i,21)=2.3*abs(hatly)*sqrt(te(i)) &
     & /btor*1.5/(elong/1.737)  !oscillating transport +- 10% of value
			
		enddo			!
!	write(*,*) work(1:na1,21)
!	pause
!	call smearr(0.1*abc/abc,work(:,21),work(:,22))
	d_1=CSOL1*work(:,21) !*exp(2.*(1.-NE/NI)) !chii

	d_2=csol3*csol1*work(:,21) !chie

	d_3=csol4*d_2 !dn



	call smearr(smear_c,d_1,work(:,21))
	call smearr(smear_c,d_2,work(:,22))
	call smearr(smear_c,d_3,work(:,23))
	do i=1,na1
!	if (car53(i).lt.-2) work(i,23)=work(i,23)/1.e4
!	if (car53(i).lt.-2) work(i,21)=work(i,21)/1.e4
!	if (car53(i).lt.-2) work(i,22)=work(i,22)/1.e4
	enddo

		do i=1,na1
		J=i
			include 'fml/lte'


	nueih=rtor*ne(i)/te(i)**2.d0
			work(i,24)=-abs(work(i,23))* &
     & (0.11*rtor/lte+0.11*shear(i) &
     & -0.1*nueih-0.01*ne(i)*te(i)/btor**2.)/rtor !cn !test
	enddo

!	write(*,*) work(1:na1,21)
!	write(*,*) work(1:na1,22)

	endif

! hatli model with F Palermos regression
	if (JJ.eq.111) 		then
	
		do i=1,na1
	J=i
			include 'fml/hatly'
!			write(*,*) hatly
!	call random_number(rnd_nmb)
!	write(*,*) 'rnd',rnd_nmb
	work(i,21)=1.3*abs(hatly)*sqrt(te(i)) &
     & /btor*1.5
		enddo			
		q_pos=0.8
		p_pos=0.91
	work(1:na1,21)=work(1:na1,21)/work(nint(na1/2.),21)

	ngree=10.*ipl/GP/abc**2.
	jk=nint(p_pos*na1)
	wped = wtotr(roc)-wtotr(p_pos*roc)+ &
     & .0024*(NE(JK)*TE(JK)+NI(JK)*TI(JK))*volum(jk)

	hreg=0.48*(ne(jk)/ngree)**0.3563* &
     & te(jk)**0.4901* &
     & rtor**-0.7644* &
     & btor**-0.7164* &
     & exp(-0.0145/50.*qtokr(roc))* &
     & mu(na1-1)**-0.8029* &
     & (rtor/abc)**1.4712

	chi_scal=0.7*((QTOKr(q_pos*ROC)-QIDWTr(q_pos*ROC)- &
     & QEDWTr(q_pos*ROC)+QRADr(q_pos*roc)))*abc**2. &
     & /(hreg*CMHD1*((QTOKr(ROC)-QIDWTr(ROC)- &
     & QEDWTr(ROC)+QRADr(roc)))-wped)
		 			

	write(777,*) hreg,cmhd2

!	write(*,*) work(1:na1,21)
!	pause
!	call smearr(0.1*abc/abc,work(:,21),work(:,22))
	d_1=chi_scal*csol1*work(:,21) !*exp(2.*(1.-NE/NI)) !chii

	d_2=chi_scal*csol3*csol1*work(:,21) !chie

	d_3=csol4*0.8*d_2 !dn



	call smearr(smear_c,d_1,work(:,21))
	call smearr(smear_c,d_2,work(:,22))
	call smearr(smear_c,d_3,work(:,23))
	do i=1,na1
!	if (car53(i).lt.-2) work(i,23)=work(i,23)/1.e4
!	if (car53(i).lt.-2) work(i,21)=work(i,21)/1.e4
!	if (car53(i).lt.-2) work(i,22)=work(i,22)/1.e4
	enddo

		do i=1,na1
		J=i
			include 'fml/lte'


	nueih=rtor*ne(i)/te(i)**2.d0
			work(i,24)=-abs(work(i,23))* &
     & (0.2*rtor/lte+0.22*shear(i))/rtor !cn
	enddo

!	write(*,*) work(1:na1,21)
!	write(*,*) work(1:na1,22)

	endif

!based on some physics
	if (JJ.eq.2) 		then
	calib_f = 2.*3*CSOL1	

!gb factor for m^2/s
	dum1=0.0032*sqrt(te*amj)/zmj/btor/rtor  !rhostar
	dum2=3.0944e+05*sqrt(te/amj)
	dum1=dum1**2.*dum2*rtor !chi_GB

		do i=1,na1
	J=i
			include 'fml/lte'
			include 'fml/lti'
			include 'fml/lne'
			include 'fml/nues'
			include 'fml/tpf'

!"""average mode frequency""" positive for TEM
	dum6(i)=(lti-lte*tpf-lne*tpf*3.)
	if (dum6(i).ne.0.) dum6(i)=dum6(i)/abs(dum6(i))

	dum6(i)=dum6(i)*(0.15)
	
	dum7(i)=TPF

	dum8(i)=NUES

!"""average growth rate"""
	if (dum6(i).ge.0.) then
		dum4(i)=(rtor/lte+rtor/lne*3.)*0.3
!chi_i
			d_1(i)=dum4(i)/3.
!chi_e
			d_2(i)=dum4(i)*dum7(i)
	endif
	if (dum6(i).lt.0.) then
		dum4(i)=sqrt(rtor/lti)*0.2
!chi_i
			d_1(i)=dum4(i)
!chi_e
			d_2(i)=dum4(i)*dum7(i)
	endif

!d_n
			d_3(i)=dum7(i)*dum4(i)*0.3


			C_T(i)=-dum7(i)*dum6(i)

			C_P(i)=-dum7(i)*0.3

			C_NU(i)=dum7(i)*NUES*dum6(i)
			
		enddo			

	call smearr(smear_c,d_1,work(:,21))
	call smearr(smear_c,d_2,work(:,22))
	call smearr(smear_c,d_3,work(:,23))

		do i=1,na1
		J=i
			include 'fml/lte'
			work(i,24)=-abs(work(i,23))* &
     & (C_T(i)*rtor/lte+C_P(i)*shear(i)+C_NU(i))/rtor !cn
	enddo


	work(:,21)=calib_f*dum1*work(:,21)*dum1/MU**1.2
	work(:,22)=calib_f*dum1*work(:,22)*dum1/MU**1.2
	work(:,23)=calib_f*dum1*work(:,23)*dum1/MU**1.2
	work(:,24)=calib_f*dum1*work(:,24)*dum1/MU**1.2

!	write(*,*) work(1:na1,21)
!	write(*,*) work(1:na1,22)

	endif







	if (JJ.eq.3) 		then

!	call glf161d

	d_1=0.1*work(:,1)
	d_2=0.1*work(:,2)
	d_3=0.1*work(:,3)
	d_4=0.1*work(:,4)

	call smearr(smear_c,d_1,work(:,21))
	call smearr(smear_c,d_2,work(:,22))
	call smearr(smear_c,d_3,work(:,23))
	call smearr(smear_c,d_4,work(:,24))

	endif

	
! TGLF 
	if (JJ.eq.324) 		then
		do i=1,na1
	J=i
			include 'fml/hatly'
!			write(*,*) hatly
!	call random_number(rnd_nmb)
!	write(*,*) 'rnd',rnd_nmb
	work(i,21)=0.7*abs(hatly)*sqrt(te(i)) &
     & /btor*1.5 !*max(grti/5.,1.)**2.                   ! added mu to make q dependence shallower
		enddo			
!	write(*,*) work(1:na1,21)
!	pause
!	call smearr(0.1*abc/abc,work(:,21),work(:,22))
	d_1=CSOL1*work(:,21) !*exp(2.*(1.-NE/NI)) !chii

	d_2=csol1*csol3*work(:,21) !chie

	d_3=csol4*d_2 !dn

		do i=1,na1
		J=i
	tpfft=FTLLMR(RHO(J))
			include 'fml/lte'


	nueih=rtor*ne(i)/te(i)**2.d0
!			work(i,24)=-abs(work(i,23))/rtor*
!     & (
!     & tpfft/0.7*(0.2*rtor/lte+0.22*shear(i)-0.1*nueih) 
!     & -(1.-tpfft)*0.1*te(i)*ne(i)/btor**2.
!     & ) !cn
	grne=1./(rho(i+1)-rho(i))* &
     & (ne(i+1)-ne(i))/ne(i)

!	work(i,24)=work(i,24)+d_3(i)*grne
			work(i,24)=-abs(d_3(i))* &
     & (0.2*rtor/lte+0.22*shear(i))/rtor !cn
	
	enddo
	d_4=work(:,24)

	if (time.ge.timnow) then
		dummh=0.
		do i=20,90
!			call nnetwroud(1,i,dummh(i,1:5))
		enddo			
		timnow=time+delayz
	write(*,*) TIME,'writing nnetwork inputs'
	endif

		dummh=0.
		do i=20,90
!			call nnetwroud(0,i,dummh(i,1:5))
			work(i,81)=min(15.,max(0.001,work(i,41)))
			work(i,82)=min(15.,max(0.001,work(i,42)))
			work(i,84)=min(5.,max(-5.,work(i,44)))
			work(i,91)=FTAV(work(i,81),0.5*abc/abc);
			work(i,92)=FTAV(work(i,82),0.5*abc/abc);
			work(i,94)=FTAV(work(i,84),0.5*abc/abc);
		enddo			
	grnee=0.
		do i=20,90
	work(i,21)=abs(work(i,91))               
	work(i,22)=abs(work(i,92))
	work(i,24)=(work(i,94))
	grnee(i)=1./(rho(i+1)-rho(i))* &
     & (ne(i+1)-ne(i))/ne(i)
!	work(i,24)=work(i,44)
		enddo			
	write(*,*) work(50,44)
	
	d_1(20:90)=work(20:90,21) !*exp(2.*(1.-NE/NI)) !chii

	d_2(20:90)=work(20:90,22) !chie

	d_3(20:90)=csol4*0.5*(d_2(20:90)+d_1(20:90)) !dn

	d_4(30:70)=work(30:70,24)+d_3(30:70)*grnee(30:70)
	write(*,*) work(50,24),d_4(50),d_3(50),grnee(50)
		
	call smearr(smear_c,d_1,work(:,21))
	call smearr(smear_c,d_2,work(:,22))
	call smearr(smear_c,d_3,work(:,23))
	call smearr(smear_c,d_4,work(:,24))
!	write(*,*) 'chuis',TIME,work(50,21),d_2(50)



	endif

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
      end
