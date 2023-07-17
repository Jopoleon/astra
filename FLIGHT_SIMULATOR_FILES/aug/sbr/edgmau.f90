	subroutine EDGMAU(c_model,P_LH,ped_width2,chi_edge,&
     & machinez,pres_clamp)


	use parameter_inc
	use const_inc
	use status_inc
	use outcmn_inc
      use fenix_params
	use  astra2fbe, only: x_point_save, plasma_config       ! declaration of minimal CPOs

	implicit none

	integer ix1,jx1,ix2,jx2
	integer	 J,i,istep,j_init
	real*8 ped_width2,exbsh(NA1),chi_edge,mmm
	real*8 dum1,dum2(NA1),nnn,c1,c2,c3,c4,P_LH
	real*8 dum4,dum5,dum6,dum7,c5,machinez,vdia,rotsh,cs
	real*8 pres_clamp,plh_multip,psep_minimum
	data c5 /1./
	save c5
	save c1,c2,c3,c4
	save	j_init
	integer imodel
	double precision c_model,qtokr,betanr,necar
	double precision qedwtr,qidwtr,ftav,qeiaur,qitotr
	data j_init /0/

	
!	write(*,*) 'edgmau',time,zrd78

	CDVM4=1.

	if (IPEQL.eq.5) then
		alp0=0.995
		if (plasma_config.eq.0) alpnew0=0.
		if (plasma_config.eq.1) alpnew0=alp0
		rx00=x_point_save(20,1)
		zx00=x_point_save(20,2)
	endif	
	

!to redo alpnew stuff since variables name changed...	
	!default
!	if (btipdirec.eq.1)	plh_multip=1. !fav bgradb drift
!	if (btipdirec.eq.-1)	plh_multip=2. !unfav bgradbdrift
	plh_multip=1.
	write(*,*) alp0,alpnew0,rx00,zx00
	if (alpnew0.lt.alp0-1.e-6)	then
		write(*,*) 'LIM'
		if (btipdirec.eq.1)	plh_multip=2. !fav bgradb drift
		if (btipdirec.eq.-1)	plh_multip=2. !unfav bgradbdrift
	endif

	if (alpnew0.ge.alp0-1.e-6.and.zx00.ge.0.) then !USN
		write(*,*) 'USN'
		if (btipdirec.eq.1)	plh_multip=2. !fav bgradb drift
		if (btipdirec.eq.-1)	plh_multip=1. !unfav bgradbdrift
	endif

	if (alpnew0.ge.alp0-1.e-6.and.zx00.lt.0.) then !LSN
		write(*,*) 'LSN'
		if (btipdirec.eq.1)	plh_multip=1. !fav bgradb drift
		if (btipdirec.eq.-1)	plh_multip=2. !unfav bgradbdrift
	endif
	write(*,*) 'plh multip',plh_multip

	if (TRIAN.ge.0.1) then
		ZRD39=5.3*RTOR**(-0.38)*TRIAN**0.83*IPL**1.25*&    !Efable test, increase pedestal height
     & ELONG**0.62*BETANR(ROC)**0.43; !8.3 before
		else
		ZRD39=3.3*RTOR**(-0.38)*exp(5.*TRIAN)/1.5*0.1**0.83&
     & *IPL**1.25*&
     & ELONG**0.62*BETANR(ROC)**0.43; !8.3 before	
	endif

	if (TIME-TSTART.lt.2*TAU) then
	else
	
	imodel=lhmodel

	!edge model for edge transport

	if (j_init.eq.0) then	
!		open(32,file='dat/exbdata.dat')
!		read(32,*) c1,c2,c3,c4
!		close(32)
		j_init=1
	endif
	

!	write(*,*) 'edgemod ',imodel


	if (imodel.eq.2) then !psep model ala Loarte 2016   + pedestal top pressure constraint
	
	dum4=QTOKR(ROC)-QEDWTR(ROC)-QIDWTR(ROC);
	dum7=max(0.01,FTAV(dum4,CF7));


!	write(*,*) 'edge ',dum7,P_LH
!	write(*,*) qedwtr(roc),qidwtr(roc)

!        if (iresta.eq.0..and.TIME-TSTART.lt.0.001) dum7=P_LH*1.2 !EF this is to avoid crashing before reaching steady-state...
	dum5=0.+(1.-tanh(20.*(dum7/P_LH-1.)))/0.5;

	i=NA1-nint(ped_width*NA1)

		dum1=1.e3*1.602*(NE(i)*TE(i)+NI(i)*TI(i))
		dum1=dum1*4*GP*1.e-7*2./BTOR**2.
	
	dum1=(TE(nint((1.-CV8)*na1))*NE(nint((1.-CV8)*na1)))*1.1

	if (isnan(zrd39)) then
	c5=100.
	else
	c5=zrd39
	endif

	c3=0.25
	dum5=dum5/mu(na1)**2./25.
	if (pr_clamp.eq.0)	chi_edge=dum5	
	if (pr_clamp.eq.1)	chi_edge=dum5	+c3*(dum1/c5 )**4.
	if (isnan(chi_edge)) chi_edge = 0.01
!	if (chi_edge.gt.100.) chi_edge = 0.01
	
!	write(*,*) 'chi_edge ',chi_edge,dum1,c5,dum5,zrd39
!	write(*,*) ifdfvx(51),betanr(roc),csol1
!	write(*,*) elong,rtor,trian,ipl,ndeut(na1)

	if (dum7/P_LH.ge.1.) CSOL1=hmodetransp ! H mode core transport
	if (dum7/P_LH.lt.1.) CSOL1=lmodetransp ! L mode core transport

!	write(*,*) 'edgemod ',dum7/P_LH,dum7,P_LH
	endif	






	if (imodel.eq.211) then !psepions = plh/2

!EFABLE TEST, to change this afther alpnew alp e zx0 exist in FEQIS
	if (IPEQL.eq.5) then
	plh_multip=1 !force this for IPEQL =5 Feqis
	endif
!!!!
	
	P_LH=0.0029*(sum(ne(1:na1))/na1)**0.72&
     & *BTOR**0.78*G11(nint(0.95*NA1))**0.93 !H mode qi threshold in MW
	CDVM2=0.5*1.53*btor**.78*rtor**1.75*&
     & (abc/rtor)**0.57*(necar(roc)/10.)**0.66
	CDVM2=0.0021*(sum(ne(1:na1))/na1)**1.07&
     & *BTOR**0.76*G11(nint(0.95*NA1))
	P_LH=0.7*CDVM2 !Efable test
	CDVM3=0.5*0.0139*(sum(ne(na1-8:na1))/9)**1.&
     & *BTOR**0.2*G11(nint(0.95*NA1)) !i mode qi threshold in MW

	
	dum4=QITOTR(ROC)+QEIAUR(ROC)-QIDWTR(ROC); !for implicit 
!	dum4=QITOTR(ROC)+0*QEIAUR(ROC)-QIDWTR(ROC); !for explicit 
	dum7=max(0.01,FTAV(dum4,CF7));
!	write(*,*) 'plh ,',dum4,cf7,dum7,QTOTR(ROC),QITOTR(ROC),
!     & QEIAUR(ROC),QIDWTR(ROC)
	psep_minimum=min(plh_multip*P_LH,CDVM3)
!	psep_minimum=plh_multip*P_LH !,CDVM3)

!	write(9431,'(222E25.11)') TIME,CDVM2,CDVM3,plh_multip,dum7

	if (TIME-ZRD78.lt.0.1) psep_minimum=1e6
	if (TIME-ZRD78.lt.0.1) CDVM3=1e6
	if (TIME-ZRD78.lt.0.1) P_LH=1e6
&
!	write(*,*) 'threshold ',zx0,alp,alpnew,&
!     & plh_multip*P_LH,CDVM3,slat(na1),psep_minimum

	if (dum7/psep_minimum.lt.1.) ZRD15=0. ! L mode
	if (dum7/CDVM3.ge.1.) ZRD15=2. ! I mode 
	if (dum7/(plh_multip*P_LH).ge.1.) ZRD15=1. ! H mode 

	dum5=0.+(1.-tanh(20.*(dum7/psep_minimum-1.)))/0.5;

	i=NA1-nint(ped_width*NA1)

		dum1=1.e3*1.602*(NE(i)*TE(i)+NI(i)*TI(i))
		dum1=dum1*4*GP*1.e-7*2./BTOR**2.
	
	dum1=(TE(nint((1.-CV8)*na1))*NE(nint((1.-CV8)*na1)))*1.1

	if (isnan(zrd39)) then
	c5=100.
		else
	if (ZRD15.eq.1)	c5=zrd39
	if (ZRD15.eq.2)	c5=zrd39/2.d0   !Imode has half of the Hmode pedestal pressure
	endif

	dum5=dum5/mu(na1)**2./25.
	c3=0.025*(dum1/c5 )**2.

	if (pr_clamp.eq.0)	chi_edge=dum5	
	if (pr_clamp.eq.1)	chi_edge=dum5	+c3
	if (isnan(chi_edge)) chi_edge = 0.01
!	if (chi_edge.gt.100.) chi_edge = 0.01
	
!	write(*,*) 'chi_edge ',chi_edge,dum1,c5,dum5,zrd39
!	write(*,*) ifdfvx(51),betanr(roc),csol1
!	write(*,*) elong,rtor,trian,ipl,ndeut(na1)

	ZRD64=0. !no ELMs
	if (ZRD15.ge.1.) CSOL1=hmodetransp ! H mode core transport
	if (ZRD15.eq.1.) ZRD64=1. ! ELMs
	if (ZRD15.lt.1.) CSOL1=lmodetransp ! L mode core transport

	if (ZRD15.lt.1.) CDVM4=1. ! no density pedestal, L mode
	if (ZRD15.eq.2.) CDVM4=1./chi_edge ! no density pedestal, I mode
	if (ZRD15.eq.1.) CDVM4=CBND3 !Hmode


!	write(*,*) 'mode type',CDVM4,ZRD15,psep_minimum,dum7,&
!     & cdvm3,P_LH,plh_multip*P_LH,plh_multip

	CDJM6=plh_multip*P_LH
	
!	write(*,*) 'edgemod ',dum7/P_LH*2.d0,dum7,dum1/c5,&
!     & P_LH/2.d0,qitotr(roc),qidwtr(roc),csol1,chi_edge,dum5

!	write(*,*) CSOL1,ZRD15,ZRD64,CDVM4

	endif	




!write(*,*) 'edgmau diagz eq',time,psplex,psiext,ne(1),te(1), &
! & ipl,iplfbe,iplx,fp(1),fp(na1),ni(1),ti(1),&
! & elong,trian,updwn,abc,shift,ccoil(1:12),vcoil(1:12)

















	if (imodel.eq.21) then !psep model: Ti criterion, Ti critical is 3. keV
	
	dum4=QTOKR(ROC)-QEDWTR(ROC)-QIDWTR(ROC);
	dum7=max(0.01,FTAV(dum4,CF7));

	dum5=0.+(1.-tanh(20.*(TI(nint((1.-CV8)*na1))/3.d0-1.)))/0.5;

	i=NA1-nint(ped_width*NA1)

		dum1=1.e3*1.602*(NE(i)*TE(i)+NI(i)*TI(i))
		dum1=dum1*4*GP*1.e-7*2./BTOR**2.
	
	dum1=(TE(nint((1.-CV8)*na1))*NE(nint((1.-CV8)*na1)))*1.1

	if (isnan(ZRD39)) then
	else
	c5=zrd39
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
	include 'fml/vdia'
	ER(j)=BTOR*VDIA
	include 'fml/rotsh'
	include 'fml/cs'
!		exbsh(j)=ROTSH/CS*RTOR
		exbsh(j)=-ER(j)/BTOR/CS*c2
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

	write(*,*) 'exb',exbsh(NA1-nint(ped_width*NA1):NA1-1),&
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
	c4=zrd39
	c3=0.25


	chi_edge=dum7	+c3*(dum1/c4 )**4.
	
	write(*,*) 'chi',mmm,chi_edge,dum1,nnn,&
     & c3*(dum1/c4 )**4.,dum7
	
	endif	
	
	
	if (isnan(chi_edge)) chi_edge=2.	

	endif


!	write(*,*) DN(NA1-10:na1)
	if (imodel.eq.1) then !exb shearing model
!compute ExB shearing


	do j=1,NA1
	include 'fml/vdia'
	ER(j)=BTOR*VDIA
	include 'fml/rotsh'
	include 'fml/cs'
		exbsh(j)=ROTSH/CS*RTOR
	enddo

	do j=1,NA1-1
!	dum1=NE(j+1)*TE(j+1)+NI(j+1)*TI(j+1)-
!     &  (NE(j)*TE(j)+NI(j)*TI(j))
	dum1=TI(j+1)-  TI(j)
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
	c4=zrd39
	c3=0.25


	chi_edge = abs(c1*nnn/(1.+c2*mmm**2.)+c3*(dum1/c4 )**4.)


!	write(*,*) mmm,chi_edge,dum1,nnn,
!     & c1*nnn/(1.+c2*mmm**2.)

!	write(115,*) CBND3,mmm,nnn
	endif
	
	








      end
