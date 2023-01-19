	subroutine nbstatus(JINOUT,jNB1,jNA1,
     1	yNE,yNHYDR,yNDEUT,yNTRIT,yNHE3,
     2	yNALF,yNI,yNIZ1,yNIZ2,yNIZ3,
     3	yZIM1,yZIM2,yZIM3,yTE,yTI,
     4	yVR,ySHIF,ySHIV,yELON,yTRIA,
     5	yAMETR,yRHO,yFP,yMU,yAMAIN,yNN,yTN,yZEF,yG33,yIPOL,
     6	yNIBM,yPIBM,yPEBM,yPBLON,yPBPER,
     7	yPBEAM,ySNEBM,ySNNBM,yCUFI,yCUBM,
     8	ySCUBM,ySNIBM1,ySNIBM2,ySNIBM3,
     9	yNNBM1,yNNBM2,yNNBM3)
    
! subroutine fills nbi/nbstatus.inc
        use nbstatus_inc

	implicit none

        double precision	
     1	YNE(*),YNHYDR(*),YNDEUT(*),YNTRIT(*),YNHE3(*),
     2	YNALF(*),YNI(*),YNIZ1(*),YNIZ2(*),YNIZ3(*),
     3	YZIM1(*),YZIM2(*),YZIM3(*),YTE(*),YTI(*),
     4	YVR(*),YSHIF(*),YSHIV(*),YELON(*),YTRIA(*),
     5  YAMETR(*),YRHO(*),YFP(*),YMU(*),YAMAIN(*),
     6	YNIBM(*),YPIBM(*),YPEBM(*),YPBLON(*),YPBPER(*),
     7	YPBEAM(*),YSNEBM(*),YSNNBM(*),YCUFI(*),YCUBM(*),
     8	YSCUBM(*),YSNIBM1(*),YSNIBM2(*),YSNIBM3(*),YG33(*),YIPOL(*),
     9	YNNBM1(*),YNNBM2(*),YNNBM3(*),YNN(*),YTN(*),YZEF(*)
	integer JNB1,J,JNA1,JINOUT
!	write(*,*) 'JINOUT',JINOUT

	if(JINOUT.eq.0) then	!in
	do J=1,JNB1
      		NE(j)	=YNE(j)
		NHYDR(j)=YNHYDR(j)
		NDEUT(j)=YNDEUT(j)
		NTRIT(j)=YNTRIT(j)
		NHE3(j)=YNHE3(j)
      		NALF(j)=YNALF(j)
		NI(j)=YNI(j)
		NIZ1(j)=YNIZ1(j)
		NIZ2(j)=YNIZ2(j)
		NIZ3(j)=YNIZ3(j)
      		ZIM1(j)=YZIM1(j)
		ZIM2(j)=YZIM2(j)
		ZIM3(j)=YZIM3(j)
		TE(j)=YTE(j)
		TI(j)=YTI(j)
      		VR(j)=YVR(j)
		SHIF(j)=YSHIF(j)
		SHIV(j)=YSHIV(j)
		ELON(j)=YELON(j)
		TRIA(j)=YTRIA(j)
        	AMETR(j)=YAMETR(j)
		RHO(j)=YRHO(j)
		FP(j)=YFP(j)
		MU(j)=YMU(j)
		AMAIN(j)=YAMAIN(j)
      		NIBM(j)=YNIBM(j)
		PIBM(j)=YPIBM(j)
		PEBM(j)=YPEBM(j)
		PBLON(j)=YPBLON(j)
		PBPER(j)=YPBPER(j)
      		PBEAM(j)=YPBEAM(j)
		SNEBM(j)=YSNEBM(j)
		SNNBM(j)=YSNNBM(j)
		CUFI(j)=YCUFI(j)
		CUBM(j)=YCUBM(j)
      		SCUBM(j)=YSCUBM(j)
		SNIBM1(j)=YSNIBM1(j)
		SNIBM2(j)=YSNIBM2(j)
		SNIBM3(j)=YSNIBM3(j)
      		NNBM1(j)=YNNBM1(j)
		NNBM2(j)=YNNBM2(j)
		NNBM3(j)=YNNBM3(j)
		NN(j)	=YNN(j)
		TN(j)	=YTN(j)
		ZEF(j)	=YZEF(j)
		G33(j)	=YG33(j)
		IPOL(j) =YIPOL(j)
	enddo
	else	!out
	do J=1,JNB1
      		YNIBM(j)=NIBM(j)
		YPIBM(j)=PIBM(j)
		YPEBM(j)=PEBM(j)
		YPBLON(j)=PBLON(j)
		YPBPER(j)=PBPER(j)
      		YPBEAM(j)=PBEAM(j)
		YSNEBM(j)=SNEBM(j)
		YSNNBM(j)=SNNBM(j)
		YCUFI(j)=CUFI(j)
		YCUBM(j)=CUBM(j)
      		YSCUBM(j)=SCUBM(j)
		YSNIBM1(j)=SNIBM1(j)
		YSNIBM2(j)=SNIBM2(j)
		YSNIBM3(j)=SNIBM3(j)
      		YNNBM1(j)=NNBM1(j)
		YNNBM2(j)=NNBM2(j)
		YNNBM3(j)=NNBM3(j)
	ENDDO
	endif
	return
	end
