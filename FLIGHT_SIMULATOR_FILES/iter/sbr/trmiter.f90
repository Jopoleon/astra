!C======================================================================|
	subroutine TRMITER
!C
!C----------------------------------------------------------------------|
	
	use parameter_inc
	use status_inc
	use const_inc
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
	
	if (TIME.le.0.) then
	work(1:na1,21)=0.3+xrho(1:na1)**2.	
	work(1:na1,22)=0.3+xrho(1:na1)**2.	
	work(1:na1,23)=0.3+xrho(1:na1)**2.	
	work(1:na1,24)=0.*(0.3+xrho(1:na1)**2.)	

	return	
	endif	
	
	smear_c = 1.e-3
	transp_variance=0.05  ! 5% transport coefficient variance and randomly generated at each time step
		
		do i=1,na1
	J=i
			include 'fml/hatly'

	work(i,21)=2.2*abs(hatly)*sqrt(te(i)) &
     & /btor*1.5  
						
		enddo			!

	d_1=CSOL1*work(:,21) !*exp(2.*(1.-NE/NI)) !chii

	d_2=csol3*csol1*work(:,21) !chie

	d_3=csol4*d_2 !dn



	call smearr(smear_c,d_1,car41)
	call smearr(smear_c,d_2,car42)
	call smearr(smear_c,d_3,car43)

		do i=1,na1
		J=i
			include 'fml/lte'
	nueih=rtor*ne(i)/te(i)**2.d0
			car44(i)=-abs(CAR43(i))* &
     & (0.15*rtor/lte+0.15*shear(i) &
     & -0.*nueih-0.05*ne(i)*te(i)/btor**2.)/rtor !cn !test
	enddo


return
	
	
	
	
	
	
	
	
      end
