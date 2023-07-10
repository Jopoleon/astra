!C======================================================================|
	subroutine EDGMDEM(c_model,P_LH,ped_width2,chi_edge, &
     & machinez,pres_clamp)
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


	if (TIME-TSTART.lt.2*TAU) then
	else
	
	imodel=nint(c_model)

	!edge model for edge transport

!	open(32,file='dat/exbdata.dat')
!	read(32,*) c1,c2,c3,c4
!	close(32)

!	write(*,*) DN(NA1-10:na1)
	if (imodel.eq.1) then !exb shearing model
!compute ExB shearing


	do j=1,NA1
!	include 'fml/vdia'
	ER(j)=BTOR*VDIA
!	include 'fml/rotsh'
!	include 'fml/cs'
!		exbsh(j)=ROTSH/CS*RTOR
	enddo

	do j=1,NA1-1
!	dum1=NE(j+1)*TE(j+1)+NI(j+1)*TI(j+1)-
!     &  (NE(j)*TE(j)+NI(j)*TI(j))
	dum1=TI(j+1)-TI(j)
!	write(*,*) dum1,j,NE(j+1)*TE(j+1)+NI(j+1)*TI(j+1),
!     & NE(j)*TE(j)+NI(j)*TI(j)
		dum2(j)=-dum1/HRO/(TI(j))*RTOR
	enddo

	i=NA1-nint(ped_width*NA1)

		dum1=1.e3*1.602*(NE(i)*TE(i)+NI(i)*TI(i))
		dum1=dum1*4*GP*1.e-7*2./BTOR**2.

	
!	write(*,*) exbsh(NA1-nint(ped_width*NA1):NA1),
!     & dum2(NA1-nint(ped_width*NA1):NA1),c1,c2,c3,c4

	mmm=sum(exbsh(NA1-nint(ped_width*NA1):NA1))/nint(ped_width*NA1)
	nnn=sum(dum2(NA1-nint(ped_width*NA1):NA1))/nint(ped_width*NA1)

!average elm effect

	dum1=(TE(nint(0.94*na1))*NE(nint(0.94*na1)))*1.1
	c4=ZRD21
	c3=0.25


	chi_edge = abs(c1*nnn/(1.+c2*mmm**2.)+c3*(dum1/c4 )**4.)


!	write(*,*) mmm,chi_edge,dum1,nnn,
!     & c1*nnn/(1.+c2*mmm**2.)

!	write(115,*) CBND3,mmm,nnn
	endif
	
	
!	write(*,*) 'edgempd'	
	if (imodel.eq.2) then !psep model ala Loarte 2016   + pedestal top pressure constraint
	
	dum4=QTOKR(ROC)-QEDWTR(ROC)-QIDWTR(ROC);
	dum7=max(0.01,FTAV(dum4,CF7));
!	write(*,*) dum7
!        if (iresta.eq.0..and.TIME-TSTART.lt.50.) dum7=P_LH*1.2 !EF this is to avoid crashing before reaching steady-state...
	dum5=0.+(1.-tanh(20.*(dum7/P_LH-1.)))/0.5;

	i=NA1-nint(ped_width*NA1)

		dum1=1.e3*1.602*(NE(i)*TE(i)+NI(i)*TI(i))
		dum1=dum1*4*GP*1.e-7*2./BTOR**2.
	
	dum1=(TE(nint((1.-CV8)*na1))*NE(nint((1.-CV8)*na1)))*1.1

	if (isnan(ZRD21)) then
	else
	c5=ZRD21
	endif

	c3=0.25

	if (nint(pres_clamp).eq.0)	chi_edge=dum5	
	if (nint(pres_clamp).eq.1)	chi_edge=dum5	+c3*(dum1/c5 )**4.
	
!		write(*,*) mmm,chi_edge,dum1,nnn,
!     & c1*nnn/(1.+c2*mmm**2.)
!	write(*,*) dum4,dum7,dum5,dum1,c5,c3,chi_edge
 
!	write(*,*) 'edgempd2'	

	endif	

	if (imodel.eq.21) then !psep model: Ti criterion, Ti critical is 3. keV
	
	dum4=QTOKR(ROC)-QEDWTR(ROC)-QIDWTR(ROC);
	dum7=max(0.01,FTAV(dum4,CF7));

	dum5=0.+(1.-tanh(20.*(TI(nint((1.-CV8)*na1))/3.d0-1.)))/0.5;

	i=NA1-nint(ped_width*NA1)

		dum1=1.e3*1.602*(NE(i)*TE(i)+NI(i)*TI(i))
		dum1=dum1*4*GP*1.e-7*2./BTOR**2.
	
	dum1=(TE(nint((1.-CV8)*na1))*NE(nint((1.-CV8)*na1)))*1.1

	if (isnan(ZRD21)) then
	else
	c5=ZRD21
	endif

	c3=0.25

	if (nint(pres_clamp).eq.0)	chi_edge=dum5	
	if (nint(pres_clamp).eq.1)	chi_edge=dum5	+c3*(dum1/c5 )**4.
	
!		write(*,*) mmm,chi_edge,dum1,nnn,
!     & c1*nnn/(1.+c2*mmm**2.)
!	write(*,*) dum4,dum7,dum5,dum1,c5,c3,chi_edge
 

	endif	

!	write(*,*) 'edgempd'	
	if (imodel.eq.245) then !psep model ala Loarte 2016   + pedestal top pressure constraint, for qi
	
!	dum4=QITOTR(ROC)+QEIADR(ROC)-QIDWTR(ROC);
	dum4=QITOTR(ROC)-QIDWTR(ROC);   ! this is probably not so good when Ti boundary is different than Te boundary, rather this should be at the pedestal top or soemthing so that the region close to the sep does not suck all the energy out of the ions when tib > teb
	dum7=FTAV(dum4,cf7)

!	dum7=(dum4+qfluxprev)/2;
!	qfluxprev=dum4
!        if (iresta.eq.0..and.TIME-TSTART.lt.50.) dum7=P_LH*1.2 !EF this is to avoid crashing before reaching steady-state...
!	write(*,*) dum7,P_LH/2.,
!     & QTOKR(ROC)-QEDWTR(ROC)-QIDWTR(ROC),
!     & QETOTR(ROC)-QEIADR(ROC)-QEDWTR(ROC),
!     & QEIADR(ROC)
!	write(*,*) dum7,P_LH/2.,
!     & QTOKR(ROC)-QEDWTR(ROC)-QIDWTR(ROC),
!     & QETOTR(ROC)-QEDWTR(ROC),
!     & QEIADR(ROC)
        if (TIME-TSTART.lt.50.) dum7=P_LH*1.2 !EF this is to avoid crashing before reaching steady-state...
!	dum7=1000.001
	dum5=0.+(1.-tanh(20.*(dum7/(P_LH/2.)-1.)))/0.5;

	write(*,*) 'lh ',dum4,dum7,P_LH/2.,cdhj5,cv3
	write(*,*) 'lh ',QITOTR(RHO(90))-QIDWTR(RHO(90))

	dum8=(NE(nint((1.-CV8)*na1))/6.5)**0.95
	dum8=dum8*abs((QTOKR(ROC)-QEDWTR(ROC)-QIDWTR(ROC))/130)**0.47
!	dum4=1.1*5.5*6.5*min(1.,max(0.01,dum8))
!	dum4=1.1*3.*6.5*min(1.,max(0.01,dum8))
	dum4=1.
	write(*,*) 'presfac',dum4,P_LH, &
     & abs(QTOKR(ROC)-QEDWTR(ROC)-QIDWTR(ROC)),zrd21
	
	i=NA1-nint(ped_width*NA1)

		dum1=1.e3*1.602*(NE(i)*TE(i)+NI(i)*TI(i))
		dum1=dum1*4*GP*1.e-7*2./BTOR**2.
	
	dum1=(TE(nint((1.-CV8)*na1))*NE(nint((1.-CV8)*na1)))*1.1

	if (isnan(ZRD21)) then
	else
	c5=ZRD21
	endif

	c3=0.25

!	c5=1.0 !testefable

	if (nint(pres_clamp).eq.0)	chi_edge=dum5	
	if (nint(pres_clamp).eq.1)	chi_edge=dum5	+c3*(dum1/c5/dum4 )**4.
	
!		write(*,*) mmm,chi_edge,dum1,nnn,
!     & c1*nnn/(1.+c2*mmm**2.)
!	write(*,*) dum4,dum7,dum5,dum1,c5,c3,chi_edge
 
!	write(*,*) 'edgempd2'	

	endif	

	if (imodel.eq.21) then !psep model: Ti criterion, Ti critical is 3. keV
	
	dum4=QTOKR(ROC)-QEDWTR(ROC)-QIDWTR(ROC);
	dum7=max(0.01,FTAV(dum4,CF7));

	dum5=0.+(1.-tanh(20.*(TI(nint((1.-CV8)*na1))/3.d0-1.)))/0.5;

	i=NA1-nint(ped_width*NA1)

		dum1=1.e3*1.602*(NE(i)*TE(i)+NI(i)*TI(i))
		dum1=dum1*4*GP*1.e-7*2./BTOR**2.
	
	dum1=(TE(nint((1.-CV8)*na1))*NE(nint((1.-CV8)*na1)))*1.1

	if (isnan(ZRD21)) then
	else
	c5=ZRD21
	endif

	c3=0.25

	if (nint(pres_clamp).eq.0)	chi_edge=dum5	
	if (nint(pres_clamp).eq.1)	chi_edge=dum5	+c3*(dum1/c5 )**4.
	
!		write(*,*) mmm,chi_edge,dum1,nnn,
!     & c1*nnn/(1.+c2*mmm**2.)
!	write(*,*) dum4,dum7,dum5,dum1,c5,c3,chi_edge
 

	endif	

	
	if (imodel.eq.3) then !rotsh > gamma_ITG
	
	
	do j=1,NA1
!	include 'fml/vdia'
	ER(j)=BTOR*VDIA
!	include 'fml/rotsh'
!	include 'fml/cs'
!		exbsh(j)=ROTSH/CS*RTOR
!		exbsh(j)=-ER(j)/BTOR/CS*c2
	enddo

	do j=1,NA1-2
!	dum1=NE(j+1)*TE(j+1)+NI(j+1)*TI(j+1)-
!     &  (NE(j)*TE(j)+NI(j)*TI(j))
	dum1=TI(j+1)-TI(j)
!	write(*,*) dum1,j,NE(j+1)*TE(j+1)+NI(j+1)*TI(j+1),
!     & NE(j)*TE(j)+NI(j)*TI(j)
		dum2(j)=-dum1/HRO/(TI(j))*RTOR
	enddo
	dum2=sqrt(max(0.,dum2))*c1


	nnn=sum(dum2(NA1-nint(ped_width*NA1):NA1-1))/nint(ped_width*NA1)
	mmm=sum(exbsh(NA1-nint(ped_width*NA1):NA1-1))/nint(ped_width*NA1)

	write(*,*) 'exb',exbsh(NA1-nint(ped_width*NA1):NA1-1), &
     &     dum2(NA1-nint(ped_width*NA1):NA1-1),mmm,nnn,&
     &     TI(NA1-nint(ped_width*NA1)-1:NA1)

	dum4=mmm/nnn;
	dum5=0.+(1.-tanh(20.*(dum4-1.)))/0.5;
	write(*,*) 'dums',dum4,dum5	
	dum7=FTAV(dum5,CF7);

	i=NA1-nint(ped_width*NA1)

		dum1=1.e3*1.602*(NE(i)*TE(i)+NI(i)*TI(i))
		dum1=dum1*4*GP*1.e-7*2./BTOR**2.

	
!	write(*,*) exbsh(NA1-nint(ped_width*NA1):NA1),
!     & dum2(NA1-nint(ped_width*NA1):NA1),c1,c2,c3,c4


!average elm effect

	dum1=(TE(nint(0.94*na1))*NE(nint(0.94*na1)))*1.1
	c4=ZRD21
	c3=0.25


	chi_edge=dum7	+c3*(dum1/c4 )**4.
	
	write(*,*) 'chi',mmm,chi_edge,dum1,nnn,&
     & c3*(dum1/c4 )**4.,dum7
	
	endif	
	
	
	if (isnan(chi_edge)) chi_edge=2.	

	endif




      end
