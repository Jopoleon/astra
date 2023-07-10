	SUBROUTINE	MIXINS(RECOND,sawyes)

	use parameter_inc
	use const_inc
	use status_inc
	use outcmn_inc

	implicit none
	double precision RHOS(3),PSIM(3),FPSTAR(NRD),WBPOLR,ICDR,&
     & QRES,taures,shearz(na1)
	double precision sawyes,QPRF(NA1),PHII(NA1),volvol
	integer	IS(3),IOPT,KOPT,JNRES,j,jj,IMIX,JNNRES,BBB(NA1)
	double precision TMIX,RECOND,NEAVR,YMURES,YMUC,YXS,YRS,YWP,YX,&
     &		HFP,YSIGN,YSIGNO,YFPC,YCM,YCM1,D1PS,D2PS,RHOMIX,YV0,DPS,&
     &		YV2,YV3,YV4,YV5,YV6,YNE,YNI,YTE,YTI,YX2,YCU,YM,YM1,YT,&
     &		DEMAX,DIMAX,TEMAX,TIMAX,DEMIN,DIMIN,TEMIN,TIMIN,DDE,DDI,&
     &		DECONT,DICONT,WECONT,WICONT,YTE0,YTI0,YNE0,YDTE,YDTI,a_res, &
     &		yf1(9)
	character*132	STRI
	save	TMIX
	data	TMIX/-99999./



	QPRF=1./MU(1:NA1)
	QRES = 1.
	sawyes=0.
	YMURES = 1./QRES
	JNRES=NA1
	JNNRES=NA1
	PHII=RHO(1:NA1)**2.*BTOR*GP
	FPSTAR(1:NA1)=FP(1:NA1)-FP(1)-PHII+PHII(1)
	do   j=1,NA1
		if (QPRF(j).le.1.) JNRES=j
	enddo

	do   j=2,NA1
		BBB(j)=QPRF(j+1)-QPRF(j)
		if (j.le.jnres.and.BBB(j).le.0.&
     & .and.QPRF(j).gt.QRES) return
	enddo
	do   j=1,NA1
		if (QPRF(j).gt.1.) JNNRES=j
	enddo

	if (JNNRES.le.JNRES) return
	SHEARz(1)=0.
	do   j=2,NA1
		if (FPSTAR(j).gt.0.) IMIX=j
			SHEARz(j)=ametr(j)*mu(j)*(1./mu(j)-1./mu(j-1))&
     & /(ametr(j)-ametr(j-1))
	enddo

!	write(*,*) 'shear',jnres,SHEAR(jnres),recond
	if (SHEARz(jnres).lt.RECOND) return

!	taures=0.4*GP*CC(jnres)/RHO(jnres)**2.


!save sawtooth info
!	write(5345,'(333E25.11)') TIME,XRHO(jnres),
!     & SHEARz(jnres),ICDR(ROC),CDVM6,CDHJ5,
!     & taures 
	

	
!sawteeth yes	

	a_res=AMETR(jnres)	

	sawyes=1.

!clamp the q profile to 1
	QPRF(1:jnres)=1.
	MU(1:na1)=1./QPRF(1:na1)
	volvol=VOLUM(IMIX)


!calculate FP
	do j=1,jnres-1

		FP(jnres-j)=FP(jnres-j+1)+&
     & (PHII(jnres-j)-PHII(jnres-j+1))*&
     & MU(jnres-j+1)
	enddo				

	yne=0.
	yf1=0.
	yni=0.
	yte=0.
	yti=0.
	
	do j=1,IMIX
	yne=yne+NE(j)*VR(j)*HRO
	yni=yni+NI(j)*VR(j)*HRO
	yf1(1)=yf1(1)+F1(j)*VR(j)*HRO
	yf1(2)=yf1(2)+F2(j)*VR(j)*HRO
	yf1(3)=yf1(3)+F3(j)*VR(j)*HRO
	yf1(4)=yf1(4)+F4(j)*VR(j)*HRO
	yf1(5)=yf1(5)+F5(j)*VR(j)*HRO
	yf1(6)=yf1(6)+F6(j)*VR(j)*HRO
	yf1(7)=yf1(7)+F7(j)*VR(j)*HRO
	yf1(8)=yf1(8)+F8(j)*VR(j)*HRO
	yf1(9)=yf1(9)+F9(j)*VR(j)*HRO
	yte=yte+NE(j)*TE(j)*VR(j)*HRO
	yti=yti+NI(j)*TI(j)*VR(j)*HRO
	enddo

	TE(1:IMIX)=yte/yne
	TI(1:IMIX)=yti/yni
	NE(1:IMIX)=yne/volvol
	F1(1:IMIX)=yf1(1)/volvol
	F2(1:IMIX)=yf1(2)/volvol
	F3(1:IMIX)=yf1(3)/volvol
	F4(1:IMIX)=yf1(4)/volvol
	F5(1:IMIX)=yf1(5)/volvol
	F6(1:IMIX)=yf1(6)/volvol
	F7(1:IMIX)=yf1(7)/volvol
	F8(1:IMIX)=yf1(8)/volvol
	F9(1:IMIX)=yf1(9)/volvol


	end
