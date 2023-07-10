!C======================================================================|
	subroutine EDGMITER
!C
!C----------------------------------------------------------------------|
	use fenix_params
	use const_inc
	use status_inc
	use parameter_inc
	use outcmn_inc
	
	implicit none
!	include	'for/parameter.inc'
!	include 'for/const.inc'
!	include 'for/status.inc'
!	include 'for/outcmn.inc'

!	include 'dat/fenixparams.var'
	integer ix1,jx1,ix2,jx2
	integer	 J,i,istep,j_init
	real*8 ped_width2,exbsh(NA1),chi_edge,mmm
	real*8 dum1,dum2(NA1),nnn,c1,c2,c3,c4,P_LH
	real*8 dum4,dum5,dum6,dum7,c5,machinez,dum8,dum9
	real*8 pres_clamp,plh_multip,psep_minimum,qfluxprev
	real*8 VDIA,ROTSH,CS,qtokr,ftav,qedwtr,qidwtr,qitotr
	data c5 /1./
	save c5
	save c1,c2,c3,c4,qfluxprev
	save	j_init
	integer imodel
	double precision c_model
	data j_init /0/



	P_LH=cdvm2

	dum4=QITOTR(ROC)-QIDWTR(ROC);   ! this is probably not so good when Ti boundary is different than Te boundary, rather this should be at the pedestal top or soemthing so that the region close to the sep does not suck all the energy out of the ions when tib > teb

	dum7=FTAV(dum4,cf7)


        if (TIME-TSTART.lt.5.) dum7=P_LH*1.2 !EF this is to avoid crashing before reaching steady-state...

	dum5=0.+(1.-tanh(20.*(dum7/(P_LH/2.)-1.)))/0.5;

	write(*,*) 'lh ',dum4,dum7,P_LH/2.,cdhj5,cv3
	write(*,*) 'lh ',QITOTR(RHO(90))-QIDWTR(RHO(90))

	dum8=(NE(nint((1.-CV8)*na1))/6.5)**0.95
	dum8=dum8*abs((QTOKR(ROC)-QEDWTR(ROC)-QIDWTR(ROC))/130)**0.47

	dum4=1.
	write(*,*) 'presfac',dum4,P_LH,abs(QTOKR(ROC)-QEDWTR(ROC)-QIDWTR(ROC)),zrd21
	
	i=NA1-nint(ped_width*NA1)

		dum1=1.e3*1.602*(NE(i)*TE(i)+NI(i)*TI(i))
		dum1=dum1*4*GP*1.e-7*2./BTOR**2.
	
	dum1=(TE(nint((1.-CV8)*na1))*NE(nint((1.-CV8)*na1)))*1.1

	if (isnan(ZRD21)) then
	else
	c5=ZRD21
	endif

	c3=0.25


	chi_edge=dum5	+c3*(dum1/c5/dum4 )**4.
	
	
	if (isnan(chi_edge)) chi_edge=2.	

	CV14=chi_edge
	


	return
	
      end
