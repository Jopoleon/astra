C	18-JUL-2013						nbco_2gc.inc
	include 'nbco_1gc.inc'
 		YQBP	=YAQBP(JE)*YDEX
		YF0	=(ARD(JE,JN)-AZ(JR))	!*YBTDB(JN)
		YFA	=YF0-ARD(JE,N1)
		YDDD	=(YAQBP(JE)-YQBP)
	if(ARD(JE,JI).LE.YF0.AND.ARD(JE,N1).GE.YF0)	THEN
!=================================================================== Trapped
		JTRAP	=0
!======================================================================
 !!  	  if(YFA.LE.0..AND.YFA.GE.YR2)	THEN			!Lost
	 if(YFA.LE.0..AND.YFA.GE.YR2.or.rcr.ge.yriplr(ji))THEN 	!======	lost
c...Lost
   		YD	=0.
   		YDCOS	=0.
		YTCOS	=YC2+YCT2(JE)
		YCT2(JE)=YC2
   		JNL	=JN1	! JN?
   		JNR	=JN1	! JN?
   		YJSURF	=1.d0
		JLOSS=1		!lost
   	  ELSE	!============================================ kept

	call nborbco(JE,JN,RCR,RJ,RJT,YC2,YF0,
     ,	JNR,JNL,JNRC,JNLC,JTRAP,JLOSS,JI)

	IF(JLOSS.EQ.0)	THEN	! ========================== kept

		YJSERF	=1./(JNR-JNL+1.d0)
		YD	=YDDD*YJSERF
		YDCOS	=YD*2.d0	
		YTCOS	=YC2+YCT2(JE)
		YCT2(JE)=YC2
	    		ELSE	! ========================== lost
C	WRITE(1,*) 'trap_lost',JN,JR,JH,ixyu
		YD	=0.
		YDCOS	=0.
		YTCOS	=YC2+YCT2(JE)
		YCT2(JE)=YC2
		JNL	=JN1
		JNR	=JN1
		YJSURF	=1.
			ENDIF	!======================== lost/kept
   	  ENDIF !========================================= lost/kept
	ELSE !=================================================== Passing
		JTRAP	=1
!======================================================================

	  if(YFA.LE.YR1.AND.YFA.GE.YR2)	THEN	!================= lost
   		YD	=0.
  		YDCOS	=0.
		YTCOS	=YC2+YCT2(JE)
		YCT2(JE)=YC2
   		JNL	=JN1
   		JNR	=JN1
   		YJSURF	=1.
		JLOSS	=0

            ELSE !================================================= kept
	call nborbco(JE,JN,RCR,RJ,RJT,YC2,YF0,
     ,	JNR,JNL,JNRC,JNLC,JTRAP,JLOSS,JI)
	IF(JLOSS.EQ.0)	THEN !============================ kept
		YJSERF	=1./(JNR-JNL+1)
ctest	WRITE(*,*) 'PASS',YJSERF,JNR,JNL
		YD	=YDDD*YJSERF
		YDCOS	=YD*2.d0	
		YTCOS	=YC2+YCT2(JE)
		YCT2(JE)=YC2
	    		ELSE !============================ lost
c...Lost
C	WRITE(1,*) 'pass_lost',JN,JR,JH,ixyu
		YD	=0.
		YDCOS	=0.
		YTCOS	=YC2+YCT2(JE)
		YCT2(JE)=YC2
		JNL	=JN1
		JNR	=JN1
		YJSURF	=1.
			ENDIF !======================== end kept/lost
!!
	  ENDIF	!============================================ end kept/lost
   	ENDIF	!========================================= end passing/trappped

!     		YATBA(JE,JN1)=YATBA(JE,JN1)+YDDD*YTCOS
!    		YATBA(JE,JN1)=YATBA(JE,JN1)-CONTR*YDDD*AZ(JR)!*2 in SCUBM
!
		YDDD=-2.d0*CONTR*YDDD*AZ(JR)*YJSERF	! orb. aver.
!!	if(JLOSS.eq.1.or.JTRAP.eq.0) then
!!		YDDD=-CONTR*YDDD*YTCOS*YJSERF
!!	else
!!		YDDD=-CONTR*YDDD*2.d0*YJSERF
!!	endif
	do	J_JN=JNL,JNR
!!	if(JLOSS.eq.1.or.JTRAP.eq.0) then
     		YATBA(JE,J_JN)=YATBA(JE,J_JN)+YDDD	
!!	else
!!     		YATBA(JE,J_JN)=YATBA(JE,J_JN)+YDDD*YCOS(J_JN)
!!	endif
     		YANBA(JE,J_JN)=YANBA(JE,J_JN)+YD ! ion source
		YANBA1(JE,JN1)=YANBA1(JE,JN1)+YD ! elecrton source
     		YAQBA(JE,J_JN)=YAQBA(JE,J_JN)+YD
	if(YTRAP(J_JN).le.YCOS(J_JN))
     &     		YACBA(JE,J_JN)=YACBA(JE,J_JN)+YDCOS*YCOS(J_JN)
     		YAQBP(JE)=YQBP
		YDYDT	=YD*YDT
	if(YCOS(J_JN).gt.0.d0) then	! co-motion
		JT	=NTET*YCOS(J_JN)*YBTDB(JN)	!+1
		if(JT.gt.NTET)	JT=NTET
		if(JT.lt.1)	JT=1
		YASBA(JE,J_JN,JT)	=YASBA(JE,J_JN,JT)+YD
	else		! contr-motion
		JT	=(-NTET*YCOS(J_JN)*YBTDB(JN))	!+1
		if(JT.gt.NTET)	JT=NTET
		if(JT.lt.1)	JT=1
		YASBA1(JE,J_JN,JT)	=YASBA1(JE,J_JN,JT)+YD
	endif
	enddo
	enddo			! end JE
