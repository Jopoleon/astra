	subroutine TRMAUG(model_type,smear_c,machinez)

	use parameter_inc
	use const_inc
	use status_inc
	use outcmn_inc
      use fenix_params
			
	implicit none

	integer	 JJ,J,i,istep
	integer	jk
	real*8 model_type,machinez
	real*8 q_pos,p_pos,yhatli,smear_c
	real*8 d_1(nrd),ngree
	real*8 d_2(nrd)
	real*8 d_3(nrd),d_5(nrd)
	real*8 d_4(nrd),grnee(nrd),hreg
	real*8 dum1(nrd)
	real*8 lte,lne,nuee,nues,tpf,lti,coulg
	real*8 dum2(nrd),delayz,timnow
	real*8 dum3(nrd),wped
	real*8 dum4(nrd),hatly
	real*8 dum5(nrd),dumhh(5),dummh(na1,5)
	real*8 dum6(nrd),bumb(na1,5)
	real*8 dum7(nrd),grte,grti,grne
	real*8 grteth,grtith,grneth,grnethte
	real*8 chiitg,chietg,chitem
	real*8 stifitg,stifetg,stiftem
	real*8 tpfft,gbfac,tauzzz



	real*8 dum8(nrd),chi_scal
	real*8 dum9(nrd)
	real*8 dum10(nrd)
	real*8 C_T(nrd)
	real*8 C_P(nrd)
	real*8 C_NU(nrd)
	real*8 calib_f,rnd_nmb,rnd_nmb2,rnd_nmb3
	real*8 nueih,scalfac
	save timnow !for tglf neuran network
	save rnd_nmb,rnd_nmb2,rnd_nmb3 !for tglf neuran network
	data timnow/0.11/ !for tglf neuran network
	delayz=3. !for tglf neuran network
	JJ=nint(model_type)
	

	

!	subroutine	SMEARR(ALFA,FO,FN)
!	write(*,*) 'trmaug',time
	
	if (TIME.le.0.) then
	work(1:na1,21)=0.3+xrho(1:na1)**2.	
	work(1:na1,22)=0.3+xrho(1:na1)**2.	
	work(1:na1,23)=0.3+xrho(1:na1)**2.	
	work(1:na1,24)=0.*(0.3+xrho(1:na1)**2.)	

	return	
	endif	
		
		
		
		
		
		
		
		
		
		
		
! hatli model
	if (JJ.eq.1) 		then
	
	call random_number(rnd_nmb)
		do i=1,na1
	J=i
			include 'fml/hatly'

	work(i,21)=abs(hatly)*sqrt(te(i))*&
     & (1.+transp_variance*2.*(0.5-rnd_nmb))&
     & /btor*1.5
			
		enddo			
		do i=1,na1
!			call nnetwrout(i,bumb(i,1:5),0)
			if (i.ge.nint(0.95*na1).or.i.le.nint(0.2*na1)) then
				bumb(i,1:5)=0.
			endif
		enddo			

	work(1:na1,502)=min(50.,abs(bumb(1:na1,1)))
	work(1:na1,501)=min(50.,abs(bumb(1:na1,2)))
	work(1:na1,503)=max(-30.,min(30.,bumb(1:na1,3)))
	work(1:na1,504)=max(-30.,min(30.,bumb(1:na1,4)))
	work(1:na1,505)=max(-30.,min(30.,bumb(1:na1,5)))

!	write(*,*) 'bumbe ', work(1:na1,502)

	d_1=CSOL1*work(:,21) !chii
	d_2=CSOL3*CSOL1*work(:,21) !chie
!	d_1 = work(:,501)
!	d_2 = work(:,502)

	d_3=CSOL4*CSOL3*CSOL1*work(:,21) !dn

	call smearr(smear_c,d_1,work(:,21))
	call smearr(smear_c,d_2,work(:,22))
	call smearr(smear_c,d_3,work(:,23))

		do i=1,na1
		J=i
			include 'fml/lte'


	nueih=rtor*ne(i)/te(i)**2.d0
			work(i,24)=-abs(work(i,23))*&
     & max(0.,(0.2*rtor/lte+0.15*shear(i)-nueih/15.d0)/rtor) !cn

!	call nnetwrout(i,dumhh(1:5),1)

	enddo



	endif




		
! TGLF neural network EF version 1
	if (JJ.eq.1001) 		then
	
		do i=1,na1
!			call nnetwrout(i,work(i,21:25),0)
		enddo			

	d_1=work(:,22) !chii
	d_2=work(:,21) !chie
	d_3=work(:,23) !ve
	d_4=work(j,24) !vz
	d_5=work(j,25) !vd

	call smearr(smear_c,d_1,work(:,21))
	call smearr(smear_c,d_2,work(:,22))
	call smearr(smear_c,d_3,work(:,23))
	call smearr(smear_c,d_4,work(:,24))
	call smearr(smear_c,d_5,work(:,25))

!	work(:,21)=work(:,21)*




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
			work(i,24)=-abs(work(i,23))*&
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
	
	
	
      end
