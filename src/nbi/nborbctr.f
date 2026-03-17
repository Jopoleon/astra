	subroutine nborbctr(JE,JN,RCR,RJ,RJT,YC2,YF0,
     ,	JNR,JNL,JNRC,JNLC,ITRAP,ILOSS,JI)
!,RBJN)
!,YCOS)
C================================================== 20-JAN-2013
c....	First orbit analysis taking account  
c	the input source redistribution 
c	due to the first orbit deviation (guide centre apprx.)
c	co_injection at birth
C input:
c	JE,JN,YC2,YF0	starting point of 
c	velocity index,surface index,YC2=<v.B>/B 
c	RCR,RJ,RJT 	= critical, external,internal boundary radii 
C	YF0 		= F value at the birth position
C	RBJN		= major radius of ion birth
C	REJ(j)/RIJ(j)	= major radius of (j) surface with mid plane LFS/HFS
C	RC(j)		= major radius of mag. axis
c output:
c	JNL,JNR	= minimum/maximum surface index of the orbit with finit RLM
c	JNLC,JNRC= minimum/maximum surface index of the orbit g.c.
C	ILOSS		= 0/1 if particle is kept/lost
c	ITRAP 	= 0/1 for banana/passing orbits
c	YCOS(JJN) = <v.B>/vB (x(JJN))
c==================================================================

        use status, only: NRD

	implicit none

	include 'nbi/nbicom.inc' 
!ARD(3,NRD),N1,YCOS(NRD),REJ(NRD),RIJ(NRD),RC(NRD)
	double precision	Y,YY,YRJ,RCR,RJ,RJT,YC2,YF0,RBJN,YY1
!,YCOS(*)
	integer	JE,JN,ITRAP,ILOSS,JNR,JNL,JNRC,JNLC,N,J,J1,JI
	logical trapped  
C====initial values
	
		JNR	=JN
		JNL	=JN
		JNRC	=JN
		JNLC	=JN
!
	ITRAP	=1 	!(def.passing)
	ILOSS	=0	!(def. kept)
 		YCOS(JN)=YC2
		J	=JN
	if(JN.gt.n1) then 
		write(*,*) 'error in nbrbco: JN > N1: ',JN,'>',N1
		return
	endif

	J1	=1	! sign of index increment
C======================================================================			
 1	if(J.gt.n1)	goto 3			! lost

		Y	=ARD(JE,J)-YF0		! motion to LFS (R->RJ)	
		yy	=RCR/2.d0		!*YBTDB(J)
		yy	=yy+dsqrt(yy*yy+y*y)	! R=Rc/2+sqrt((Rc/2)^2+y^2)
	if(j.eq.jn) yy1=yy		
!!		if(yy.le.REJ(j).and.yy.ge.RIJ(j)) then	!Ri(Z=0)<= R<=Re(Z=0) real root
		if(yy.le.REJ(j)) then	!Ri(Z=0)<= R<=Re(Z=0) real root
			if(yy.lt.RJT.or.yy.gt.RJ) goto 3 !lost
				ycos(j)	=y/yy
			if(ycos(j).lt.-1.d0)	ycos(j)=-1.d0
			if(ycos(j).gt.1.d0) 	ycos(j)=1.d0
				JNRC	=J
				JNR	=J
				J	=J+J1	! motion to LFS (R->RJ)	
			goto 1	
		endif	! exit: end of cycle to LFS: fake root
C=======================================================================
		J=JN
		J1=-1

 2	if(J.lt.1) then ! reaches plasma ceter
		ILOSS=0	! kept in plasma
		return
	endif
		Y	=ARD(JE,J)-YF0		! motion to HFS (R->RJT)	
		yy	=RCR/2.d0		!*YBTDB(J)
		yy	=yy+dsqrt(yy*yy+y*y)	! R=Rc/2+sqrt((Rc/2)^2+y^2)
		if(yy.le.REJ(j).and.yy.ge.RIJ(j)) then	!Ri(Z=0)<= R<=Re(Z=0) real root
	if(y.lt.0.) 	itrap=0	!trapped to well

			if(yy.lt.RJT.or.yy.gt.RJ) goto 3  !lost
			ycos(j)	=y/yy
			if(ycos(j).lt.-1.d0)	ycos(j)=-1.d0
			if(ycos(j).gt.1.d0) 	ycos(j)=1.d0
			JNLC	=J
			JNL	=J
			J	=J+J1	! motion to HFS (R->RJT)	
			yy1	=yy
			goto 2	
		endif	! exit: end of cycle to HFS: fake root	
	if(itrap.eq.0) then
!!	write(*,*) '=trapped,kept jrc,jn=',ji,jn,j
          ELSE					!	kept
!!	write(*,*) '=passing,kept jrc,jn=',ji,jn,j
	endif
	return
 3		ILOSS=1	! lost from plasma
!!		write(*,*) 'error 3: ILOSS = 0 => 1 '
	if(itrap.eq.0) then
!!	write(*,*) '=trapped,lost jrc,jn=',ji,jn,yy,rcr
          ELSE					!	kept
!!	write(*,*) '=passing,lost jrc,jn=',ji,jn
	endif
	return
	end	
		

