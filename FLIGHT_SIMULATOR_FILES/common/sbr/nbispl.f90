	subroutine	NBISPL

	use parameter_inc
	use const_inc
	use status_inc
	use outcmn_inc

	implicit none

	integer		jj,length,JINOUT,nbimodel
	integer icmodel, ecmodel
	double precision dum1,vint,dum8(na1),icdr
	double precision dumpower
	double precision QETOTR,QITOTR,time0000,qradr,qeiclr
	double precision qefrac_target,qtot_target,qohr
!C----------------------------------------------------------------------|

	nbimodel=nint(ZRD13) ! 2 is rabbit
	ecmodel=nint(ZRD18) ! if not 1, call torba from the model file
	icmodel=nint(ZRD19) ! 

!for hartmut
!	time0000=7.54
!	if (TIME.ge.time0000.and.
!     & time.le.time0000+4.) then
!		car32(25:26)=0.
!		car32(17:24)=0.
!		car32(1:8)=0.
!		car32(21)=5.*(time-time0000)/4.
!		car32(1)=5.-5.*(time-time0000)/4.
!	endif
!		call evaluate_qeqicontrol   !full duality case
!		call evaluate_qeqicontrol2  !ech nb scan for te/ti check
!		call evaluate_qeqicontrol3  !reproduce palpha simple case

!!!!!!!!!!!!!!!!!!



	QNBI=sum(car32(17:24))
	QECR=sum(car32(1:8))
	QICR=sum(car32(25:26))



!NBI model	done by rabbit.f
!	write(*,*) 'qnbi ',qnbi
	if (nbimodel.eq.1) then
	if (QNBI.eq.0.) then
	pebm=0.
	pibm=0.
	NIBM=0.0
	PBLON=0.0
	PBPER=0.0
	else
	PEBM=0.4*QNBI*exp(-XRHO**2./0.3)
	PIBM=0.6*QNBI*exp(-XRHO**2./0.3)
	dum1=vint(pebm+pibm,roc)
	pebm=pebm/dum1*QNBI
	pibm=pibm/dum1*QNBI
	NIBM=0.01*QNBI*exp(-XRHO**2./0.3)
	PBLON=NIBM*70.*TE/TE(1)
	PBPER=NIBM*70.*TE/TE(1)
	SNEBM=(PEBM+PIBM)*625./50.
	endif
	endif
	if (nbimodel.eq.2) call rabbit !rabbit or whatever
!	write(*,*) 'nbi',QNBI


!ICRF model:
	PEICR=0.
	PIICR=0.
	if (icmodel.eq.1) then
	if (QICR.eq.0.) then
	else
	dum8(1:na1)=exp(-(XRHO(1:na1)-0.1)**2./0.2**2.)
	call trapz_gg(VOLUM(1:na1),1, &
     & dum8(1:na1),&
     &    dum1,NA1)!		dum1=vint(dum8,roc)
	PEICR(1:na1)=PEICR(1:na1)+dum8(1:na1)/dum1*QICR
	PIICR(1:na1)=PIICR(1:na1)+dum8(1:na1)/dum1*QICR
	PEICR=0.7*PEICR*0.6    ! 40% of icrh power is lost
	PIICR=0.3*PIICR*0.6
!	NIBM=NIBM+0.01*QICR*exp(-XRHO**2./0.3)
!	PBLON=PBLON
!	PBPER=PBPER 
	endif
	endif

!	write(*,*) 'icr',QICR


	
!ECRH model:	

!test from Hartmut Zohm EFable
!		dumpower=sum(car32(1:8))
!		if (dumpower.le.0.) QECR=0.
!		if (dumpower.gt.0.) QECR=1.5
!		QECR=QECR*(0.5+0.5*cos(time*GP2*100.))
!!!!!
	
	if (ecmodel.eq.1) then
		PEECR = 0.
		CUECR = 0.
		if (QECR.eq.0.) then
		PEECR=0.
		else

	dum8(1:na1)=exp(-(XRHO(1:na1)-0.1)**2./0.02**2.)
	call trapz_gg(VOLUM(1:na1),1,&
     & dum8(1:na1),&
     &    dum1,NA1)!		dum1=vint(dum8,roc)
	PEECR(1:na1)=PEECR(1:na1)+dum8(1:na1)/dum1*QECR


		dum8(1:na1)=exp(-(XRHO(1:na1)-0.3)**2./0.02**2.)
	call trapz_gg(VOLUM(1:na1),1,&
     & dum8(1:na1),&
     &    dum1,NA1)!		dum1=vint(dum8,roc)
		PEECR(1:na1)=PEECR(1:na1)+dum8(1:na1)/dum1*0*car32(2)

	CUECR=0.*0.06*GP2*RTOR*PEECR*TE/NE*6./4.; !current drive from some scaling


!		dum8=exp(-(XRHO-0.15)**2./0.02**2.)
!		dum1=vint(dum8,roc)
!		PEECR=PEECR+dum8/dum1*car32(3)

!		dum8=exp(-(XRHO-0.6)**2./0.02**2.)
!		dum1=vint(dum8,roc)
!		PEECR=PEECR+dum8/dum1*car32(4)

		endif
	endif	
	if (ecmodel.eq.2) call torba
		
!	write(*,*) 'ech',QECR

!!test from Hartmut Zohm EFable, voltages modified btw
!	write(834214,'(6E25.11)') TIME,QECR,
!     & CU(1),CU(21),ICDR(ROC),IPL
	write(834214,'(61111E25.11)') TIME,QECR,&
     & QNBI,QICR,TE(1),TI(1),NE(1),&
     & QETOTR(ROC),QITOTR(ROC),QRADR(ROC),&
     & QEICLR(ROC),TE(20),TI(20),NE(20),&
     & QETOTR(0.5*ROC),QITOTR(0.5*ROC),&
     & QEICLR(0.5*ROC),QRADR(0.5*ROC),&
     & QOHR(ROC),qohr(0.5*ROC),car33(1:24)
!!!!!!!

	end


!C======================================================================|
!C INTEGR_EF computes integrals over x
!C  x_input: x_variable
!C  x_type for GRP style-grid is: 1 main grid -> shifted grid (deriv), 2 shifted !grid -> main grid (deriv)
!C  y_input: y_variable
!C  ys_output: integral
!C  nagrid: number of grid points
!C======================================================================|
	subroutine	trapz_gg(x_input,x_type,y_input,&
     &    ys_output,nagrid)
		 use numerical_tools, only: polyfitcc
	implicit none
	integer	i,j,k,order_d,x_type,nagrid
	double precision x_input(nagrid)
	double precision y_input(nagrid),ys_output
	double precision y1tmp,y2tmp,y3tmp,drho
	double precision x1tmp,x2tmp,x3tmp
	double precision Acoef,Bcoef,Ccoef,y0,P(3)

!C Normalized grid , GRP style
				    if (x_type.eq.1) then
							j=1
							drho=x_input(1)
							ys_output=y_input(1)*drho
							do j=2,nagrid
          			drho=x_input(j)-x_input(j-1)
								y1tmp=y_input(j)
								ys_output=ys_output+y1tmp*drho
							enddo
						endif
				    if (x_type.eq.2) then
!C First interpolate shifted variable to zero
               j=1
	call polyfitcc(x_input(1:3),y_input(1:3),P)
	y0=P(3)
          			drho=x_input(1)
								y1tmp=y0
								ys_output=y1tmp*drho
 						do j=2,nagrid
          			drho=x_input(j)-x_input(j-1)
								y1tmp=y_input(j-1)
								ys_output=ys_output+y1tmp*drho
							enddo
						endif
!CCCCCCCCCCCCCCCCCCCCCCCCCCCC



	end
!C======================================================================|
	subroutine evaluate_qeqicontrol
		 
	use parameter_inc
	use const_inc
	use status_inc
	use outcmn_inc
	implicit none
		 
	integer		jj,length,JINOUT,nbimodel,imidrad
	integer icmodel, ecmodel,icall,imethod
	double precision vint,dum8(na1),icdr
	double precision dumpower
	double precision QETOTR,QITOTR,time0000,qradr,qeiclr
	double precision qefrac_target,qtot_target
	double precision qeqimid,qohr
	double precision qeecmid,qenbimid,qradmid
	double precision qinbimid
	double precision qenet, qtotnet,qtotsep,qenetb
	double precision qauxsep,wtotr,wer
	double precision nbidutyperiod
	double precision echdutyperiod
	double precision nbirisetime
	double precision echrisetime
	double precision nbistepwidth
	double precision echstepwidth
	double precision midrad
	double precision zum1(100)
	double precision qenez(100)
	double precision qenezb(100)
	double precision qtnez(100)
	double precision qsnez(100)
	double precision C_1,C_2
	integer n_sourcesnbi
	integer n_sourcesech,iccc
	integer idum1(100),icount,i_case
	
	double precision ndemo,palpharef,pauxref
	double precision Tfus,Tdemo,Telecdemo
	double precision phmode,waugref
	double precision demoqeqiscal
	double precision radcorefrac,psepscal,naugref
	
	 
	double precision dum1(100)
	data icall/0/
	save dum1,icall,idum1,icount
	save qenez,qsnez,qtnez,qenezb

	ndemo=8.d0
	naugref=6.d0
	palpharef=8.d0
	pauxref=palpharef/10.*0. !ignited plasma if pauxref=0., but ohmic ramp-up.... Interestingly, with the palpharef/10., it will go into H-mode naturally before more power is applied. alpha power helps there.
	Tfus=25.
	Tdemo=13.
	waugref=0.5
	Telecdemo=35.
	phmode=3.
	demoqeqiscal=10.
	radcorefrac=0.5
	psepscal=4.

	car32=0.
	qtot_target=0.        ! separatrix qtot target
	qefrac_target=0.    ! mid radius of net qe/qtot
	time0000=7.54

	 imidrad=floor((na1+0.)/2.)
	 midrad=0.5*roc

	i_case=2
	imethod=2

!evaluate equipartition power
	qeqimid=0.
	if (imethod.eq.1) then
		qeqimid=qeiclr(midrad)
	endif		 
	if (imethod.eq.2) then
		qeqimid=ne(imidrad)**(2.5)* &
     & (2*wer(roc)-wtotr(roc))/ &
     & wer(roc)**(1.5)/105.
	endif		 

	if (i_case.eq.1) then
	if (TIME.ge.time0000.and.&
     & time.le.time0000+1.) then
		qtot_target=5.        ! separatrix qtot target
		qefrac_target=0.3    ! mid radius of net qe/qtot
	endif
	time0000=8.54
	if (TIME.ge.time0000.and.&
     & time.le.time0000+1.) then
		qtot_target=5.        ! separatrix qtot target
		qefrac_target=0.6    ! mid radius of net qe/qtot
	endif
		time0000=9.54
	if (TIME.ge.time0000.and.&
     & time.le.time0000+1.) then
		qtot_target=8.        ! separatrix qtot target
		qefrac_target=0.3    ! mid radius of net qe/qtot
	endif
	time0000=10.54
	if (TIME.ge.time0000.and.&
     & time.le.time0000+1.) then
		qtot_target=8.        ! separatrix qtot target
		qefrac_target=0.6    ! mid radius of net qe/qtot
	endif
	endif
	
	zum1=0.
	C_1=0.
	C_2=0.
	
	if (i_case.eq.2) then !alpha power emulation but radiated power is just imposed to get some Psep, it is not controlled with gas puff of impurities
		if (TIME.ge.5.54.and.&
     & time.le.time0000+4.) then
			if (TIME.gt.time0000.and.TIME.lt.time0000+1.) then
!		if (TIME.ge.time0000-1.and.
!     & time.le.time0000+4.) then
!			if (TIME.ge.7.54.and.time.lt.8.54) then
					zum1(20)=phmode !to enter into H-mode, auxiliary power
				else
					zum1(20)=0.
			endif

			zum1(1)=(wtotr(roc)/waugref) ! normalized energy. reference energy is waugref
!			zum1(1)=(2.*(wtotr(roc)-wer(roc))/waugref) ! normalized energy from ions. reference energy is waugref/2.


			call svdt_alphaemul(zum1(2),zum1(1)*Tfus)
			call svdt_alphaemul(zum1(31),zum1(1)*Tfus+0.1)
			call svdt_alphaemul(zum1(21),Tfus)
			zum1(3)=palpharef*zum1(2)/zum1(21) !alpha power 
!			zum1(3)=palpharef*zum1(2)/zum1(21)*
!     & (sum(ne(1:na1))/na1/naugref)**2.d0 !alpha power including density dependence
			
			zum1(2)=zum1(1)*Tfus/zum1(2)*&
     & (zum1(31)-zum1(2))/0.1 !store cross section exponent

			zum1(4)=pauxref +&
     & zum1(20) !auxiliary power&

	call paion_emul(zum1(22),zum1(1)*Telecdemo,&
     & ndemo,0.01d0,2.5d0,&
     & ndemo*0.9)
			zum1(5)=zum1(3)*(1.-zum1(22)) !alpha power to electrons
			zum1(6)=zum1(3)*zum1(22) !alpha power to ions
			
			zum1(7)=zum1(3)+zum1(4) !total heating


			zum1(8)=max(0.,zum1(7)-psepscal)  !radiated power to get psepscal MW at separatrix, to sub with qradr(roc) and put controller from impurities
!			zum1(8)=qradr(roc) !qradr(roc) has to be controlled to get psep = psepscal


			zum1(9)=qeqimid*demoqeqiscal !equipartition power up to mid-radius,the factor demoqeqiscal is scaled to DEMO, but to be improved
	
!below are prob wrong			
!			zum1(99)=zum1(5)+zum1(4)-zum1(8)-zum1(9) !net electron power at separatrix
!			zum1(89)=zum1(5)+zum1(4)-zum1(8)*radcorefrac-zum1(9) !net electron power at mid radius. Its the same as sep but radiation is cut by radcorefrac. The "core" radiation fraction will ultimately determine the electron flux to be matched. Important quantity
!			zum1(100)=zum1(6)+zum1(9) !net ion power
!!!!!
!obsolete			
!	C_1 = zum1(3)+zum1(4)-zum1(8)+qradr(roc)
!	C_2 = zum1(5)+zum1(4)-zum1(9)+qeqimid+
!     & qradr(midrad)-radcorefrac*zum1(8)
!!!!!!

		 zum1(89)=zum1(5)+zum1(4)-zum1(9)-radcorefrac*zum1(8) !elec midrad power scal. alpha power at mid radius not obviously 100% of total, needs some rescaling to be checked on demo
		 zum1(99)=zum1(6)+zum1(9) !ion midrad power scal
		 zum1(100)=zum1(7)-zum1(8) !total sep power scal
		 
		qtot_target=zum1(100)        ! separatrix qtot target
		qefrac_target=&
     & max(0.,zum1(89)/abs(zum1(89)+zum1(99)))    ! mid radius of net qe/qtot
		
		endif
	endif
		
	if (icall.eq.0) dum1=0. 
	if (icall.eq.0) idum1=0
	if (icall.eq.0) icount=0
	if (icall.eq.0) qenez=0.
	if (icall.eq.0) qenezb=0.
	if (icall.eq.0) qtnez=0.
	if (icall.eq.0) qsnez=0.
	if (icall.eq.0) icall=1
	
	
	qeecmid=vint(peecr,midrad)
	qenbimid=vint(pebm,midrad)
	qinbimid=vint(pibm,midrad)
	qauxsep=vint(peecr,roc)+&
     & vint(pebm,roc)+vint(pibm,roc)+&
     & qohr(roc)
	qradmid=radcorefrac*qradr(roc)	 


	icount=icount+1
	qenez(icount)=qeecmid+qenbimid-qradmid-qeqimid+qohr(midrad)
	qtnez(icount)=qeecmid+qenbimid+qinbimid-qradmid+qohr(midrad)
	qsnez(icount)=qauxsep-qradr(roc)
	qenezb(icount)=qenbimid

!radiation done outside
!	qenez(icount)=qeecmid+qenbimid-qeqimid
!	qtnez(icount)=qeecmid+qenbimid+qinbimid
!	qsnez(icount)=qauxsep

	iccc=20
	if (icount.eq.iccc) icount=0
	
	qenet=sum(qenez(1:iccc))/iccc
	qtotnet=sum(qtnez(1:iccc))/iccc
	qtotsep=sum(qsnez(1:iccc))/iccc
	qenetb=sum(qenezb(1:iccc))/iccc
	
	! car33(21) !W midplane puff to control radiation
!	zrd72=max(0.,zrd72+
!     & 1.1*tau*(qtotsep-psepscal)) !controls Psep and zum1(8) is equal to that
	zrd72=0.
!	zrd72=max(0.,zrd72+
!     & 0.4*tau*(zum1(8)-qradr(roc))) 	! this just sets qrad = pradscal, but doesnt control pradscal actually, in reality this should control qradr with pradscal=qradr, and psep goes to target value
!		 write(11123,*) 'W influx',zum1(8),qradr(roc),zrd72,qtotsep,psepscal

!calculate required electron heating
	dum1(1)=max(0.,dum1(1)+&
     & 250*tau*(qefrac_target-qenet/qtotnet)) 
	dum1(2)=25*(qefrac_target-qenet/qtotnet)

!calculate total heating
	dum1(3)=max(0.,dum1(3)+&
     & 50*tau*(qtot_target-qtotsep)) 
	dum1(4)=1.*(qtot_target-qtotsep)

	dum1(6) = max(0.,dum1(3)+dum1(4)) !total heating in MW


!ech heating
	dum1(5) = max(0.,dum1(1)+dum1(2)) !electron heating in MW
	dum1(5) = min(dum1(6),dum1(5)) !cut ech to total power

! nbi heating
	dum1(7)=max(0.,dum1(6)-dum1(5))

	write(4444,*) time,dum1(1),dum1(2),dum1(3),dum1(4),dum1(6),&
     & dum1(5),dum1(7)
	
	
!!! to check the above for the effect of tot heating vs qrad... not sure this is correct.	
	
!	dum1(7)=max(0.,C_1-C_2+qenetb)
!	dum1(5)=max(0.,C_2-qenetb)
	
	write(89543,'(5555E25.11)') TIME,&
     & qeecmid,qenbimid,qinbimid,&
     & qauxsep,qradmid,qenet,qtotnet,qtotsep,&
     & qeqimid,ne(imidrad)**(2.5)*&
     & (2*wer(roc)-wtotr(roc))/&
     & wer(roc)**(1.5),zum1(1:9),&
     & zum1(20),zum1(99:100),zum1(89),&
     & qradr(roc),zum1(8),radcorefrac,&
     & C_1,C_2,dum1(7),dum1(5),qenetb,&
     & qefrac_target,qtot_target
!	write(834214,'(61111E25.11)') TIME,QECR,
!     & QNBI,QICR,TE(1),TI(1),NE(1),
!     & QETOTR(ROC),QITOTR(ROC),QRADR(ROC),
!     & QEICLR(ROC),TE(20),TI(20),NE(20),
!     & QETOTR(0.5*ROC),QITOTR(0.5*ROC),
!     & QEICLR(0.5*ROC),QRADR(0.5*ROC),QOHR(ROC),qohr(0.5*roc)
	
!gyros 1:8, nbi 17:24 car32
!assume each gyro does 0.7 MW, each NBI 2.5 MW
!NBI modulation cycle is 
	nbidutyperiod=0.02
	echdutyperiod=0.02
	nbirisetime=0.003
	echrisetime=0.001
	dum1(50)=2.5 !nbi power
	dum1(51)=0.7 !ech power
	
	n_sourcesnbi=min(7,max(0,floor(dum1(7)/dum1(50))))
	if (n_sourcesnbi.ge.1) then
		car32(17:17+n_sourcesnbi-1)=dum1(50)
	endif
	
	n_sourcesech=min(7,max(0,floor(dum1(5)/dum1(51))))
	if (n_sourcesech.ge.1) then
		car32(1:1+n_sourcesech-1)=dum1(51)
	endif
	
	dum1(10)=n_sourcesnbi*dum1(50)
	dum1(11)=n_sourcesech*dum1(51)
	
	dum1(12)=min(dum1(50),dum1(7)-dum1(10)) !residual nbi power
	dum1(13)=min(dum1(51),dum1(5)-dum1(11)) !residual ech power
	
	nbistepwidth=dum1(12)/dum1(50)*nbidutyperiod
	echstepwidth=dum1(13)/dum1(51)*echdutyperiod
	

!idum1(1) NBI off/on 0/1
!idum1(2) ECH off/on 0/1
!dum1(80) is time of start NBI duty cycle
!dum1(81) is time of start ECH duty cycle

	if (nbistepwidth.gt.0.) then
		if (idum1(1).eq.0) dum1(80)=time	
		if (idum1(1).eq.0) idum1(1)=1		
	endif
	
	if (echstepwidth.gt.0.) then
		if (idum1(2).eq.0) dum1(81)=time	
		if (idum1(2).eq.0) idum1(2)=1		
	endif
	 
	if (idum1(1).eq.1) then
		dum1(90)=time-dum1(80)
		if (nbistepwidth.ge.dum1(90)) then
			n_sourcesnbi=n_sourcesnbi+1
			car32(16+n_sourcesnbi)=dum1(50)*&
     & (1.-exp(-dum1(90)/nbirisetime))
		endif				
		if (dum1(90).gt.nbidutyperiod) idum1(1)=0
	endif	 	 
	
	if (idum1(2).eq.1) then
		dum1(91)=time-dum1(81)
		if (echstepwidth.ge.dum1(91)) then
			n_sourcesech=n_sourcesech+1
			car32(n_sourcesech)=dum1(51)*&
     & (1.-exp(-dum1(91)/echrisetime))
		endif				
		if (dum1(91).gt.echdutyperiod) idum1(2)=0
	endif	 
	 
	 
	 
	return
	end	 
	 
	 
	subroutine svdt_alphaemul(svdt,ti)

	implicit none
	
	double precision svdt,ti
	
	SVDT = TI ** (-0.33333333)
      SVDT = 8.972*EXP(-19.9826*SVDT)*SVDT*SVDT*&
    & ((TI+1.0134)/(1.+6.386E-3*(TI+1.0134)**2)+&
    & 1.877*EXP(-.16176*TI*SQRT(TI)))
	 
	return
	end
	
		 
	subroutine paion_emul(paion1,te,ne,nalf,amain,&
     & ni)
		 implicit none
		 double precision paion1
		 double precision te,amain,ne,ni
		 double precision nalf
		 double precision yvalp,yllame
		 double precision yy6,yllama,yvc,yy7
		 double precision yllami,yeps
      YVALP=1.2960e+07
      YLLAME=23.9+LOG(1.e3*TE/SQRT(1.e19*NE))
      YY6=SQRT(1.e3*TE/1.e19*NE)
      YY6=YY6*(4.*AMAIN*YVALP)/(4.+AMAIN)      
		if(YY6.lt..1)YY6=0.1
       YLLAMI=14.2+LOG(YY6)
       YY6=SQRT(1.0e3*TE/1.0e19*NE)*2.0*YVALP
		if(YY6.lt..01)YY6=.01
      YLLAMA=14.2+LOG(YY6)
	YY6=YLLAMI*NI/(AMAIN*NE)
	YY6=7.3e-4/YLLAME*(YY6+YLLAMA*NALF/NE)
	
		if(YY6.lt.0.0001)YY6=.0001
	YVC=YY6**0.33*SQRT(2.0*TE*1.7564e+14)
      YEPS=YVALP/(YVC+0.0001)
      YY6=atan(0.577*(2.*YEPS-1.))       
      YY7=log((1.+YEPS)**2/(1.-YEPS+YEPS**2))
      PAION1=2./YEPS**2*(0.577*YY6-0.167*YY7+0.3)
	return
	end	 
	 
	 
	 
	 
	 
	 
		 
!	C======================================================================|
	subroutine evaluate_qeqicontrol2
		 
		 
	use parameter_inc
	use const_inc
	use status_inc
	use outcmn_inc
	
	implicit none
	
	integer		jj,length,JINOUT,nbimodel,imidrad
	integer icmodel, ecmodel,icall,imethod
	double precision vint,dum8(na1),icdr
	double precision dumpower
	double precision QETOTR,QITOTR,time0000,qradr,qeiclr
	double precision qefrac_target,qtot_target
	double precision qeqimid,qohr
	double precision qeecmid,qenbimid,qradmid
	double precision qinbimid
	double precision qenet, qtotnet,qtotsep,qenetb
	double precision qauxsep,wtotr,wer
	double precision nbidutyperiod
	double precision echdutyperiod
	double precision nbirisetime
	double precision echrisetime
	double precision nbistepwidth
	double precision echstepwidth
	double precision midrad
	double precision zum1(100)
	double precision qenez(100)
	double precision qenezb(100)
	double precision qtnez(100)
	double precision qsnez(100)
	double precision C_1,C_2
	integer n_sourcesnbi
	integer n_sourcesech,iccc
	integer idum1(100),icount,i_case
	
	double precision ndemo,palpharef,pauxref
	double precision Tfus,Tdemo,Telecdemo
	double precision phmode,waugref
	double precision demoqeqiscal
	double precision radcorefrac,psepscal,naugref
	
	 
	double precision dum1(100)
	data icall/0/
	save dum1,icall,idum1,icount
	save qenez,qsnez,qtnez,qenezb



	ndemo=8.d0
	naugref=6.d0
	palpharef=8.d0
	pauxref=palpharef/10.*0. !ignited plasma if pauxref=0., but ohmic ramp-up.... Interestingly, with the palpharef/10., it will go into H-mode naturally before more power is applied. alpha power helps there.
	Tfus=25.
	Tdemo=13.
	waugref=0.5
	Telecdemo=35.
	phmode=3.
	demoqeqiscal=3.
	radcorefrac=0.5
	psepscal=4.

	car32=0.
	qtot_target=0.        ! separatrix qtot target
	qefrac_target=0.    ! mid radius of net qe/qtot
	time0000=7.

	 imidrad=floor((na1+0.)/2.)
	 midrad=0.5*roc


!evaluate equipartition power
	qeqimid=0.
		qeqimid=qeiclr(midrad)

	car32=0.

	if (TIME.ge.time0000.and.&
     & time.le.time0000+5.) then
		CAR32(1)=(TIME-time0000)/5.*5.        ! separatrix qtot target
		CAR32(19)=5.-(TIME-time0000)/5.*5.    ! mid radius of net qe/qtot
	endif


	

	 
	write(89543,'(5555E25.11)') TIME, &
     & qeecmid,qenbimid,qinbimid, & 
     & qauxsep,qradmid,qenet,qtotnet,qtotsep, &
     & qeqimid,ne(imidrad)**(2.5)* &
     & (2*wer(roc)-wtotr(roc))/ &
     & wer(roc)**(1.5),zum1(1:9), &
     & zum1(20),zum1(99:100),zum1(89), & 
     & qradr(roc),zum1(8),radcorefrac, &
     & C_1,C_2,dum1(7),dum1(5),qenetb, &
     & qefrac_target,qtot_target
!	write(834214,'(61111E25.11)') TIME,QECR,
!     & QNBI,QICR,TE(1),TI(1),NE(1),
!     & QETOTR(ROC),QITOTR(ROC),QRADR(ROC),
!     & QEICLR(ROC),TE(20),TI(20),NE(20),
!     & QETOTR(0.5*ROC),QITOTR(0.5*ROC),
!     & QEICLR(0.5*ROC),QRADR(0.5*ROC),QOHR(ROC),qohr(0.5*roc)
	
!gyros 1:8, nbi 17:24 car32
!assume each gyro does 0.7 MW, each NBI 2.5 MW
!NBI modulation cycle is 
	 
	 
	 
	return
	end	 
	 
		 
		 
!	C======================================================================|
	subroutine evaluate_qeqicontrol3
		 
	
		use parameter_inc
	use const_inc
	use status_inc
	use outcmn_inc

	
	implicit none
		 
	
	integer		jj,length,JINOUT,nbimodel,imidrad
	integer icmodel, ecmodel,icall,imethod
	double precision vint,dum8(na1),icdr
	double precision dumpower
	double precision QETOTR,QITOTR,time0000,qradr,qeiclr
	double precision qefrac_target,qtot_target
	double precision qeqimid,qohr
	double precision qeecmid,qenbimid,qradmid
	double precision qinbimid
	double precision qenet, qtotnet,qtotsep,qenetb
	double precision qauxsep,wtotr,wer
	double precision nbidutyperiod
	double precision echdutyperiod
	double precision nbirisetime
	double precision echrisetime
	double precision nbistepwidth
	double precision echstepwidth
	double precision midrad
	double precision zum1(100)
	double precision qenez(100)
	double precision qenezb(100)
	double precision qtnez(100)
	double precision qsnez(100)
	double precision C_1,C_2
	integer n_sourcesnbi
	integer n_sourcesech,iccc
	integer idum1(100),icount,i_case
	
	double precision ndemo,palpharef,pauxref
	double precision Tfus,Tdemo,Telecdemo
	double precision phmode,waugref
	double precision demoqeqiscal
	double precision radcorefrac,psepscal,naugref
	
	 
	double precision dum1(100)
	data icall/0/
	save dum1,icall,idum1,icount
	save qenez,qsnez,qtnez,qenezb

	ndemo=8.d0
	naugref=6.d0
	palpharef=8.d0
	pauxref=palpharef/10.*0. !ignited plasma if pauxref=0., but ohmic ramp-up.... Interestingly, with the palpharef/10., it will go into H-mode naturally before more power is applied. alpha power helps there.
	Tfus=25.
	Tdemo=13.
	waugref=0.5
	Telecdemo=35.
	phmode=3.
	demoqeqiscal=3.
	radcorefrac=0.5
	psepscal=4.

	car32=0.
	qtot_target=0.        ! separatrix qtot target
	qefrac_target=0.    ! mid radius of net qe/qtot
	time0000=7.54

	 imidrad=floor((na1+0.)/2.)
	 midrad=0.5*roc

	i_case=1
	imethod=1

!evaluate equipartition power
	qeqimid=0.
		qeqimid=qeiclr(midrad)

	if (i_case.eq.1) then
		if (TIME.ge.time0000.and.&
     & time.le.time0000+5.) then
			qtot_target=(wtotr(roc)/0.5)*6.        ! separatrix qtot target
			qefrac_target=0.0    ! mid radius of net qe/qtot
		endif
	endif
	
	zum1=0.
	
	C_1=qtot_target
	C_2=0.
	
	if (icall.eq.0) dum1=0. 
	if (icall.eq.0) idum1=0
	if (icall.eq.0) icount=0
	if (icall.eq.0) qenez=0.
	if (icall.eq.0) qenezb=0.
	if (icall.eq.0) qtnez=0.
	if (icall.eq.0) qsnez=0.
	if (icall.eq.0) icall=1
	
	
	qeecmid=vint(peecr,midrad)
	qenbimid=vint(pebm,midrad)
	qinbimid=vint(pibm,midrad)
	qauxsep=vint(pebm,roc)+vint(pibm,roc)
	qradmid=radcorefrac*qradr(roc)	 


	icount=icount+1
	qenez(icount)=qeecmid+qenbimid-qradmid-qeqimid+qohr(midrad)
	qtnez(icount)=qeecmid+qenbimid+qinbimid-qradmid+qohr(midrad)
	qsnez(icount)=qauxsep
	qenezb(icount)=qenbimid

!radiation done outside
!	qenez(icount)=qeecmid+qenbimid-qeqimid
!	qtnez(icount)=qeecmid+qenbimid+qinbimid
!	qsnez(icount)=qauxsep

	iccc=20
	if (icount.eq.iccc) icount=0
	
	qenet=sum(qenez(1:iccc))/iccc
	qtotnet=sum(qtnez(1:iccc))/iccc
	qtotsep=sum(qsnez(1:iccc))/iccc
	qenetb=sum(qenezb(1:iccc))/iccc
	zum1(3)=qtotsep
	
	! car33(21) !W midplane puff to control radiation
!	zrd72=max(0.,zrd72+
!     & 1.1*tau*(qtotsep-psepscal)) !controls Psep and zum1(8) is equal to that
	zrd72=0.
!	zrd72=max(0.,zrd72+
!     & 0.4*tau*(zum1(8)-qradr(roc))) 	! this just sets qrad = pradscal, but doesnt control pradscal actually, in reality this should control qradr with pradscal=qradr, and psep goes to target value
!		 write(11123,*) 'W influx',zum1(8),qradr(roc),zrd72,qtotsep,psepscal

!calculate required electron heating
	dum1(1)=max(0.,dum1(1)+&
     & 250*tau*(qefrac_target-qenet/qtotnet)) 
	dum1(2)=25*(qefrac_target-qenet/qtotnet)

!calculate total heating
	dum1(3)=max(0.,dum1(3)+&
     & 50*tau*(qtot_target-qtotsep)) 
	dum1(4)=1.*(qtot_target-qtotsep)

	dum1(6) = max(0.,dum1(3)+dum1(4)) !total heating in MW


!ech heating
	dum1(5) = max(0.,dum1(1)+dum1(2)) !electron heating in MW
	dum1(5) = min(dum1(6),dum1(5)) !cut ech to total power

		if (TIME.ge.time0000.and.&
     & time.le.time0000+1.) then
			dum1(5)=3.
		endif
		if (TIME.ge.time0000+1.) then
			dum1(5)=1.2
		endif






! nbi heating
	dum1(7)=max(0.,dum1(6))

!	write(4444,*) time,dum1(1),dum1(2),dum1(3),dum1(4),dum1(6),
!     & dum1(5),dum1(7)
	
	
!!! to check the above for the effect of tot heating vs qrad... not sure this is correct.	
	
!	dum1(7)=max(0.,C_1-C_2+qenetb)
!	dum1(5)=max(0.,C_2-qenetb)
	
	write(89543,'(5555E25.11)') TIME, &
     & qeecmid,qenbimid,qinbimid, &
     & qauxsep,qradmid,qenet,qtotnet,qtotsep, &
     & qeqimid,ne(imidrad)**(2.5)* &
     & (2*wer(roc)-wtotr(roc))/ &
     & wer(roc)**(1.5),zum1(1:9), &
     & zum1(20),zum1(99:100),zum1(89), &
     & qradr(roc),zum1(8),radcorefrac, &
     & C_1,C_2,dum1(7),dum1(5),qenetb, &
     & qefrac_target,qtot_target
!	write(834214,'(61111E25.11)') TIME,QECR,
!     & QNBI,QICR,TE(1),TI(1),NE(1),
!     & QETOTR(ROC),QITOTR(ROC),QRADR(ROC),
!     & QEICLR(ROC),TE(20),TI(20),NE(20),
!     & QETOTR(0.5*ROC),QITOTR(0.5*ROC),
!     & QEICLR(0.5*ROC),QRADR(0.5*ROC),QOHR(ROC),qohr(0.5*roc)
	
!gyros 1:8, nbi 17:24 car32
!assume each gyro does 0.7 MW, each NBI 2.5 MW
!NBI modulation cycle is 
	nbidutyperiod=0.02
	echdutyperiod=0.02
	nbirisetime=0.003
	echrisetime=0.001
	dum1(50)=2.5 !nbi power
	dum1(51)=0.7 !ech power
	
	n_sourcesnbi=min(7,max(0,floor(dum1(7)/dum1(50))))
	if (n_sourcesnbi.ge.1) then
		car32(17:17+n_sourcesnbi-1)=dum1(50)
	endif
	
	n_sourcesech=min(7,max(0,floor(dum1(5)/dum1(51))))
	if (n_sourcesech.ge.1) then
		car32(1:1+n_sourcesech-1)=dum1(51)
	endif
	
	dum1(10)=n_sourcesnbi*dum1(50)
	dum1(11)=n_sourcesech*dum1(51)
	
	dum1(12)=min(dum1(50),dum1(7)-dum1(10)) !residual nbi power
	dum1(13)=min(dum1(51),dum1(5)-dum1(11)) !residual ech power
	
	nbistepwidth=dum1(12)/dum1(50)*nbidutyperiod
	echstepwidth=dum1(13)/dum1(51)*echdutyperiod
	

!idum1(1) NBI off/on 0/1
!idum1(2) ECH off/on 0/1
!dum1(80) is time of start NBI duty cycle
!dum1(81) is time of start ECH duty cycle

	if (nbistepwidth.gt.0.) then
		if (idum1(1).eq.0) dum1(80)=time	
		if (idum1(1).eq.0) idum1(1)=1		
	endif
	
	if (echstepwidth.gt.0.) then
		if (idum1(2).eq.0) dum1(81)=time	
		if (idum1(2).eq.0) idum1(2)=1		
	endif
	 
	if (idum1(1).eq.1) then
		dum1(90)=time-dum1(80)
		if (nbistepwidth.ge.dum1(90)) then
			n_sourcesnbi=n_sourcesnbi+1
			car32(16+n_sourcesnbi)=dum1(50)*&
     & (1.-exp(-dum1(90)/nbirisetime))
		endif				
		if (dum1(90).gt.nbidutyperiod) idum1(1)=0
	endif	 	 
	
	if (idum1(2).eq.1) then
		dum1(91)=time-dum1(81)
		if (echstepwidth.ge.dum1(91)) then
			n_sourcesech=n_sourcesech+1
			car32(n_sourcesech)=dum1(51)*&
     & (1.-exp(-dum1(91)/echrisetime))
		endif				
		if (dum1(91).gt.echdutyperiod) idum1(2)=0
	endif	 
	 
	 
	 
	return
	end	 
	 
		
