	subroutine	build_2dgrid( &
     &      Nr,Nt, &
     &      Rb,Zb, &
     &      X0,Y0,lambda2d,lambda2dp, &
     &      psig,psiax,psib,j_ok, &
     &      psigp,PSI, &
     &      dArea,Rmaj2,dArea2, &
     &      dArc_rp1, &
     &      dArc_rm1, &
     &      dArc_rpt1, &
     &      dArc_rmt1, &
     &      dArc_tp1, &
     &      dArc_tm1, &
     &      dArc_tpr1, &
     &      dArc_tmr1, &
     &      ddr, &
     &      ddr_i, &
     &      dtp, &
     &      dtm, &
     &      dt_i,Rmaj,YY,Y2, &
     &      r,thetap, &
     &      thetap_i, &
     &      Jcbn, &
     &      gradr2,Jcbn2,r_i1,x_i1,y_i1)
		 
	use pi_grec_vars, only: GPI2
	implicit none

	integer Nt,Nr,iax,jax,j_ok
	double precision Rb(Nt),Zb(Nt),X0,Y0,derivs(8)
	double precision lambda2d(Nr,Nt),psiax,psib
	double precision lambda2dold(Nr,Nt)
	double precision lambda2doldp(Nr,Nt)
	double precision lambda2dp(Nr,Nt)
	double precision lambda2di(Nr,Nt)
	double precision lambda2dpi(Nr,Nt)
	double precision psig(Nr),psi(nr,nt)
	double precision psigp(Nr),psin(nr,nt)
	double precision dArea(Nr,nt)
	double precision dArc_rp1(Nr,Nt)
	double precision dArc_rm1(Nr,Nt)
	double precision dArc_rpt1(Nr,Nt)
	double precision dArc_rmt1(Nr,Nt)
	double precision dArc_tp1(Nr,Nt)
	double precision dArc_tm1(Nr,Nt)
	double precision dArc_tpr1(Nr,Nt)
	double precision dArc_tmr1(Nr,Nt)
	double precision dArea2(Nr,nt)
	double precision ddr(Nr,Nt)
	double precision  ddr_i(Nr,Nt)
	double precision  dtp(Nr,Nt)
	double precision  dtm(Nr,Nt)
	double precision  dt_i(Nr,Nt),dlambda
	double precision  Xbnew(Nt),Ybnew(Nt)
	double precision  Xbold(Nt),Ybold(Nt)
	double precision  Xb(Nt),Yb(Nt)
	double precision  dXb0(Nt),thetap(Nt+1)
	double precision  dXb0_i(Nt),thetap_i(Nt+1)
	double precision  r(Nr,nt),rr2(Nr,nt)
	double precision  X(Nr,nt),X2(Nr,nt)
	double precision  XX(Nr,Nt),YY(Nr,Nt),dpsia
	double precision  Y(Nr,nt),Y2(Nr,nt)
	double precision  r_i(Nr,nt),r_i1(Nr,nt)
	double precision  X_i(Nr,nt),X_i1(Nr,nt)
	double precision  Y_i(Nr,nt),Y_i1(Nr,nt)
	double precision  dXdr2(Nr,nt),dYdr2(Nr,nt)
	double precision  dXdr(Nr,nt),dYdr(Nr,nt)
	double precision  dXdri(Nr,nt),dYdri(Nr,nt)
	double precision  dXdri1(Nr,nt),dYdri1(Nr,nt)
	double precision  dXdh2(Nr,nt),dYdh2(Nr,nt)
	double precision  dXdh(Nr,nt),dYdh(Nr,nt)
	double precision  dXdhi(Nr,nt),dYdhi(Nr,nt)
	double precision  dXdhi1(Nr,nt),dYdhi1(Nr,nt)
	double precision  Jcbn(Nr,nt),Jcbni1(Nr,nt)
	double precision  Jcbn2(Nr,nt),Jcbni(Nr,nt)
	double precision  drdX,dhdX,drdY,dhdY
	double precision  gradr(Nr,nt),gradr2(Nr,nt)
	double precision  dphdR(Nr,nt),dphdZ(Nr,nt)
	double precision  gradh(Nr,nt),gradh2(Nr,nt)
	double precision  gradri(Nr,nt),gradri1(Nr,nt)
	double precision  gradhi(Nr,nt),gradhi1(Nr,nt)
	double precision  grt(Nr,nt),t1,t2,t3,t4
	double precision grti1(Nr,nt),Mdet
	double precision grt2(Nr,nt)
	double precision ghht2(Nr,nt)
	double precision grti(Nr,nt),Rmaj(Nr,nt)
	double precision Rmaj2(Nr,nt)
	double precision Rmaji(Nr,nt)
	double precision Rmaji1(Nr,nt)
	double precision dr_el,dt_el,Area
	integer  jok,jfirst,jr,jt,jjt,j,k






	do jt=1,Nt
			dXb0(jt)=sqrt((rb(jt)-X0)**2.0+(zb(jt)-Y0)**2.0)
		call find_angle_ef(X0,Y0,rb(jt),zb(jt),thetap(jt))
	enddo

	if (thetap(1).gt.thetap(nt)) then !reorder
		j=1
		do k=1,nt-1
			if (thetap(k+1).lt.thetap(k)) j=k+1   ! j is the first teta above 0
		enddo		
		if (j.gt.1)	thetap(1:j-1)=thetap(1:j-1)-GPI2
	endif
!	write(*,*) thetap(1:nt)

	
        jt=Nt
        thetap(jt+1)=GPI2+thetap(1)

!C Now find the intermediate theta grid
	do jt=1,Nt
		thetap_i(jt)=(thetap(jt)+thetap(jt+1))/2.0
	enddo
	jt=Nt
 		thetap_i(jt+1)=thetap_i(1)+GPI2
		
	do jt=1,Nt-1
		dXb0_i(jt)=(dXb0(jt)+dXb0(jt+1))/2.0
	enddo
	jt=Nt
		dXb0_i(jt)=dXb0_i(1)


!	write(*,*) 'ax',psiax,psi(1,1),j_ok,iax,jax

	if (j_ok.eq.1) then !relambda
!relambda
		lambda2dold(1:nr,1:nt)=lambda2d(1:nr,1:nt)
	do jt=1,Nt
		psin(1,jt)=0.
		psin(nr,jt)=1.
		do jr=2,Nr
!			ddr(jr,jt)=r(jr,jt)-r(jr-1,jt)
			psin(jr,jt)=(psi(jr,jt)-psiax)/(psib-psiax)
!		  ddr_i(jr,jt)=(psi(jr,jt)-psi(jr-1,jt))/(psib-psiax)
!			dphdR(jr,jt)=r(jr,jt)
!			dphdZ(jr,jt)=rr2(jr,jt)
		enddo
	enddo
	do jt=1,Nt
!		r(2,jt)=
!     & sqrt(psig(2)*(psib-psiax)/
!     & (0.5*derivs(3)*cos(thetap(jt))**2.+
!     & 0.5*derivs(4)*sin(thetap(jt))**2.+
!     & derivs(5)*cos(thetap(jt))*sin(thetap(jt))
!     & ))
!		rr2(1,jt)=
!     & sqrt(psigp(1)*(psib-psiax)/
!     & (0.5*derivs(3)*cos(thetap(jt))**2.+
!     & 0.5*derivs(4)*sin(thetap(jt))**2.+
!     & derivs(5)*cos(thetap(jt))*sin(thetap(jt))
!     & ))
!		lambda2d(2,jt)=r(2,jt)/dxb0(jt)
!		lambda2dp(1,jt)=rr2(1,jt)/dxb0(jt)
!			call linterp_ef_equilef(psin(:,jt),lambda2dold(:,jt),
!     & nr,psig(2:nr-1),lambda2d(2:nr-1,jt),nr-2)
!			call linterp_ef_equilef(psig,lambda2d(:,jt),
!     & nr,psigp(1:nr-1),lambda2dp(1:nr-1,jt),nr-1)
			call linterp_feqis(psin(:,jt),lambda2dold(:,jt), &
     & nr,psig(2:nr-1),lambda2d(2:nr-1,jt),nr-2)
			call linterp_feqis(psig,lambda2d(:,jt), &
     & nr,psigp(1:nr-1),lambda2dp(1:nr-1,jt),nr-1)

		do jr=2,Nr-1
			lambda2d(jr,jt)=min(1.,lambda2d(jr,jt))
			lambda2dp(jr,jt)=min(1.,lambda2dp(jr,jt))
		enddo
	
	
	enddo
!	write(*,*) 'relambda'
!	write(*,*) lambda2d(1:nr,1),lambda2dp(1:nr-1,1),psig(1:nr)
!	pause	

		else

!	write(*,*) 'no relambda'
	
	
	endif




	do jt=1,Nt
		do jr=1,Nr
			r(jr,jt)=lambda2d(jr,jt)*dXb0(jt)
			rr2(jr,jt)=lambda2dp(jr,jt)*dXb0(jt)
			X(jr,jt)=X0+r(jr,jt)*cos(thetap(jt))
			Y(jr,jt)=Y0+r(jr,jt)*sin(thetap(jt))
			X2(jr,jt)=X0+rr2(jr,jt)*cos(thetap(jt))
			Y2(jr,jt)=Y0+rr2(jr,jt)*sin(thetap(jt))
			enddo
	enddo







	do jt=1,Nt-1
		lambda2di(:,jt)=(lambda2d(:,jt)+lambda2d(:,jt+1))/2.0
		lambda2dpi(:,jt)=(lambda2dp(:,jt)+lambda2dp(:,jt+1))/2.0
	enddo
			jt=Nt
		lambda2di(:,jt)=(lambda2d(:,jt)+lambda2d(:,1))/2.0
		lambda2dpi(:,jt)=(lambda2dp(:,jt)+lambda2dp(:,1))/2.0
		
		
		
	do jt=1,Nt
		do jr = 1,Nr
			r_i1(jr,jt)=lambda2di(jr,jt)*dXb0_i(jt)
			r_i(jr,jt)=lambda2dpi(jr,jt)*dXb0_i(jt)
			X_i1(jr,jt)=X0+r_i1(jr,jt)*cos(thetap_i(jt))
			Y_i1(jr,jt)=Y0+r_i1(jr,jt)*sin(thetap_i(jt))
			X_i(jr,jt)=X0+r_i(jr,jt)*cos(thetap_i(jt))
			Y_i(jr,jt)=Y0+r_i(jr,jt)*sin(thetap_i(jt))
	enddo
	enddo
				
!C Now computes jacobian and so on
!C Jacobian is defined as:
!C
!C J = det( {dX/dr dX/dtheta}  {dY/dr dY/dtheta} ) = dX/dr*dY/dtheta - dX/dtheta*dY/dr= J(r,theta)
!C

	jt=1
	do jr = 1,Nr-1
		dXdr2(jr,jt)=(X(jr+1,jt)-X(jr,jt))/   &  !half radial grid, full angle grid
     &   (psig(jr+1)-psig(jr))	
		dYdr2(jr,jt)=(Y(jr+1,jt)-Y(jr,jt))/  & 
     &   (psig(jr+1)-psig(jr))	
		dXdri(jr,jt)=(X_i1(jr+1,jt)-X_i1(jr,jt))/  &     !  hald radial grid, half angle grid
     &   (psig(jr+1)-psig(jr))	
		dYdri(jr,jt)=(Y_i1(jr+1,jt)-Y_i1(jr,jt))/  & 
     &  (psig(jr+1)-psig(jr))	
		dXdh2(jr,jt)=(X_i(jr,jt)-X_i(jr,Nt))/   &  !half radial grid, full angle grid
     &   (thetap_i(Nt+1)-thetap_i(Nt))		
		dYdh2(jr,jt)=(Y_i(jr,jt)-Y_i(jr,Nt))/  & 
     &   (thetap_i(Nt+1)-thetap_i(Nt))
		dXdhi(jr,jt)=(X2(jr,jt+1)-X2(jr,jt))/    &  !  hald radial grid, half angle grid
     &   (thetap(jt+1)-thetap(jt))		
		dYdhi(jr,jt)=(Y2(jr,jt+1)-Y2(jr,jt))/  & 
     &   (thetap(jt+1)-thetap(jt))
	enddo
	do jt=2,Nt-1
	do jr = 1,Nr-1
		dXdr2(jr,jt)=(X(jr+1,jt)-X(jr,jt))/  & 
     &   (psig(jr+1)-psig(jr))		
		dYdr2(jr,jt)=(Y(jr+1,jt)-Y(jr,jt))/  & 
     &   (psig(jr+1)-psig(jr))		
		dXdri(jr,jt)=(X_i1(jr+1,jt)-X_i1(jr,jt))/  & 
     &   (psig(jr+1)-psig(jr))	
		dYdri(jr,jt)=(Y_i1(jr+1,jt)-Y_i1(jr,jt))/  & 
     &  (psig(jr+1)-psig(jr))	
		dXdh2(jr,jt)=(X_i(jr,jt)-X_i(jr,jt-1))/  & 
     &   (thetap_i(jt)-thetap_i(jt-1))		
		dYdh2(jr,jt)=(Y_i(jr,jt)-Y_i(jr,jt-1))/  & 
     &   (thetap_i(jt)-thetap_i(jt-1))
		dXdhi(jr,jt)=(X2(jr,jt+1)-X2(jr,jt))/  & 
     &   (thetap(jt+1)-thetap(jt))		
		dYdhi(jr,jt)=(Y2(jr,jt+1)-Y2(jr,jt))/  & 
     &   (thetap(jt+1)-thetap(jt))
	enddo
	enddo
	jt=Nt
	do jr = 1,Nr-1
		dXdr2(jr,jt)=(X(jr+1,jt)-X(jr,jt))/  & 
     &   (psig(jr+1)-psig(jr))		
		dYdr2(jr,jt)=(Y(jr+1,jt)-Y(jr,jt))/  & 
     &   (psig(jr+1)-psig(jr))		
		dXdri(jr,jt)=(X_i1(jr+1,jt)-X_i1(jr,jt))/  & 
     &   (psig(jr+1)-psig(jr))	
		dYdri(jr,jt)=(Y_i1(jr+1,jt)-Y_i1(jr,jt))/   &   
     &  (psig(jr+1)-psig(jr))	
		dXdh2(jr,jt)=(X_i(jr,jt)-X_i(jr,jt-1))/  & 
     &   (thetap_i(jt)-thetap_i(jt-1))		
		dYdh2(jr,jt)=(Y_i(jr,jt)-Y_i(jr,jt-1))/  & 
     &   (thetap_i(jt)-thetap_i(jt-1))
		dXdhi(jr,jt)=(X2(jr,1)-X2(jr,jt))/  & 
     &   (thetap(jt+1)-thetap(jt))		
		dYdhi(jr,jt)=(Y2(jr,1)-Y2(jr,jt))/  & 
     &   (thetap(jt+1)-thetap(jt))
	enddo

		Jcbn2 = dXdr2*dYdh2-dXdh2*dYdr2	  ! i+1/2, j	
	Jcbni = dXdri*dYdhi-dXdhi*dYdri		! i+1/2, j+1/2 
	
	
	jt = 1
	do jr = 2,Nr-1
		dXdr(jr,jt) = (X2(jr,jt) - X2(jr-1,jt))/  &   ! full grid
     &   (psigp(jr)-psigp(jr-1))		
		dYdr(jr,jt) = (Y2(jr,jt) - Y2(jr-1,jt))/  & 
     &   (psigp(jr)-psigp(jr-1))		
		dXdh(jr,jt) = (X_i1(jr,jt) - X_i1(jr,Nt))/   & 
     &   (thetap_i(Nt+1)-thetap_i(Nt))		
		dYdh(jr,jt) = (Y_i1(jr,jt) - Y_i1(jr,Nt))/  & 
     &   (thetap_i(Nt+1)-thetap_i(Nt))
		dXdri1(jr,jt) = (X_i(jr,jt) - X_i(jr-1,jt))/   &   ! full radial grid, half angle grid
     &   (psigp(jr)-psigp(jr-1))
		dYdri1(jr,jt) = (Y_i(jr,jt) - Y_i(jr-1,jt))/  & 
     &   (psigp(jr)-psigp(jr-1))
		dXdhi1(jr,jt) = (X(jr,jt+1) - X(jr,jt))/  & 
     &   (thetap(jt+1)-thetap(jt))		
		dYdhi1(jr,jt) = (Y(jr,jt+1) - Y(jr,jt))/  & 
     &   (thetap(jt+1)-thetap(jt))
	enddo
	do jt=2,Nt-1
	do jr = 2,Nr-1
		dXdr(jr,jt) = (X2(jr,jt) - X2(jr-1,jt))/  & 
     &   (psigp(jr)-psigp(jr-1))
		dYdr(jr,jt) = (Y2(jr,jt) - Y2(jr-1,jt))/    & 
     &   (psigp(jr)-psigp(jr-1))	
		dXdri1(jr,jt) = (X_i(jr,jt) - X_i(jr-1,jt))/  & 
     &   (psigp(jr)-psigp(jr-1))	   
		dYdri1(jr,jt) = (Y_i(jr,jt) - Y_i(jr-1,jt))/ & 
     &   (psigp(jr)-psigp(jr-1))	
		dXdh(jr,jt) = (X_i1(jr,jt) - X_i1(jr,jt-1))/ & 
     &   (thetap_i(jt)-thetap_i(jt-1))		
		dYdh(jr,jt) = (Y_i1(jr,jt) - Y_i1(jr,jt-1))/ & 
     &   (thetap_i(jt)-thetap_i(jt-1))
		dXdhi1(jr,jt) = (X(jr,jt+1) - X(jr,jt))/ & 
     &   (thetap(jt+1)-thetap(jt))		
		dYdhi1(jr,jt) = (Y(jr,jt+1) - Y(jr,jt))/ & 
     &   (thetap(jt+1)-thetap(jt))
	enddo
	enddo
	jt=Nt
	do jr = 2,Nr-1
		dXdr(jr,jt) = (X2(jr,jt) - X2(jr-1,jt))/ & 
     &   (psigp(jr)-psigp(jr-1))
		dYdr(jr,jt) = (Y2(jr,jt) - Y2(jr-1,jt))/ & 
     &   (psigp(jr)-psigp(jr-1))	
		dXdri1(jr,jt) = (X_i(jr,jt) - X_i(jr-1,jt))/ & 
     &   (psigp(jr)-psigp(jr-1))	
		dYdri1(jr,jt) = (Y_i(jr,jt) - Y_i(jr-1,jt))/ & 
     &   (psigp(jr)-psigp(jr-1))	
		dXdh(jr,jt) = (X_i1(jr,jt) - X_i1(jr,jt-1))/ & 
     &   (thetap_i(jt)-thetap_i(jt-1))		
		dYdh(jr,jt) = (Y_i1(jr,jt) - Y_i1(jr,jt-1))/ & 
     &   (thetap_i(jt)-thetap_i(jt-1))
		dXdhi1(jr,jt) = (X(jr,1) - X(jr,jt))/ & 
     &   (thetap(jt+1)-thetap(jt))		
		dYdhi1(jr,jt) = (Y(jr,1) - Y(jr,jt))/ & 
     &   (thetap(jt+1)-thetap(jt))
	enddo

	jt = 1
	jr = 1
		dXdr(jr,jt) = (X2(jr,jt) - X(jr,jt))/  &  ! full grid
     &   (psigp(jr))		
		dYdr(jr,jt) = (Y2(jr,jt) - Y(jr,jt))/ & 
     &   (psigp(jr))		
		dXdh(jr,jt) = (X_i(jr,jt) - X_i(jr,Nt))/ & 
     &   (thetap_i(Nt+1)-thetap_i(Nt))		
		dYdh(jr,jt) = (Y_i(jr,jt) - Y_i(jr,Nt))/ & 
     &   (thetap_i(Nt+1)-thetap_i(Nt))
		dXdri1(jr,jt) = (X_i(jr,jt) - X(jr,jt))/  &   ! full radial grid, half angle grid
     &   (psigp(jr))
		dYdri1(jr,jt) = (Y_i(jr,jt) - Y(jr,jt))/ & 
     &   (psigp(jr))
		dXdhi1(jr,jt) = (X2(jr,jt+1) - X2(jr,jt))/ & 
     &   (thetap(jt+1)-thetap(jt))		
		dYdhi1(jr,jt) = (Y2(jr,jt+1) - Y2(jr,jt))/ & 
     &   (thetap(jt+1)-thetap(jt))
	do jt=2,Nt-1
		dXdr(jr,jt) = (X2(jr,jt) - X(jr,jt))/ & 
     &   (psigp(jr))
		dYdr(jr,jt) = (Y2(jr,jt) - Y(jr,jt))/ & 
     &   (psigp(jr))	
		dXdri1(jr,jt) = (X_i(jr,jt) - X(jr,jt))/ & 
     &   (psigp(jr))	
		dYdri1(jr,jt) = (Y_i(jr,jt) - Y(jr,jt))/ & 
     &   (psigp(jr))	
		dXdh(jr,jt) = (X_i(jr,jt) - X_i(jr,jt-1))/ & 
     &   (thetap_i(jt)-thetap_i(jt-1))		
		dYdh(jr,jt) = (Y_i(jr,jt) - Y_i(jr,jt-1))/ & 
     &   (thetap_i(jt)-thetap_i(jt-1))
		dXdhi1(jr,jt) = (X2(jr,jt+1) - X2(jr,jt))/ & 
     &   (thetap(jt+1)-thetap(jt))		
		dYdhi1(jr,jt) = (Y2(jr,jt+1) - Y2(jr,jt))/ & 
     &   (thetap(jt+1)-thetap(jt))
	enddo
	jt=Nt
		dXdr(jr,jt) = (X2(jr,jt) - X(jr,jt))/ & 
     &   (psigp(jr))
		dYdr(jr,jt) = (Y2(jr,jt) - Y(jr,jt))/ & 
     &   (psigp(jr))	
		dXdri1(jr,jt) = (X_i(jr,jt) - X(jr,jt))/ & 
     &   (psigp(jr))	
		dYdri1(jr,jt) = (Y_i(jr,jt) - Y(jr,jt))/ & 
     &   (psigp(jr))	
		dXdh(jr,jt) = (X_i(jr,jt) - X_i(jr,jt-1))/ & 
     &   (thetap_i(jt)-thetap_i(jt-1))		
		dYdh(jr,jt) = (Y_i(jr,jt) - Y_i(jr,jt-1))/ & 
     &   (thetap_i(jt)-thetap_i(jt-1))
		dXdhi1(jr,jt) = (X2(jr,1) - X2(jr,jt))/ & 
     &   (thetap(jt+1)-thetap(jt))		
		dYdhi1(jr,jt) = (Y2(jr,1) - Y2(jr,jt))/ & 
     &   (thetap(jt+1)-thetap(jt))

	jt = 1
	jr =Nr
		dXdr(jr,jt) = (X(jr,jt) - X2(jr-1,jt))/ &   ! full grid
     &   (psig(jr)-psigp(jr-1))		
		dYdr(jr,jt) = (Y(jr,jt) - Y2(jr-1,jt))/ & 
     &   (psig(jr)-psigp(jr-1))		
		dXdh(jr,jt) = (X_i1(jr,jt) - X_i1(jr,Nt))/ & 
     &   (thetap_i(Nt+1)-thetap_i(Nt))		
		dYdh(jr,jt) = (Y_i1(jr,jt) - Y_i1(jr,Nt))/ & 
     &   (thetap_i(Nt+1)-thetap_i(Nt))
		dXdri1(jr,jt) = (X_i1(jr,jt) - X_i(jr-1,jt))/  &   ! full radial grid, half angle grid
     &   (psig(jr)-psigp(jr-1))
		dYdri1(jr,jt) = (Y_i1(jr,jt) - Y_i(jr-1,jt))/ & 
     &   (psig(jr)-psigp(jr-1))
		dXdhi1(jr,jt) = (X(jr,jt+1) - X(jr,jt))/ & 
     &   (thetap(jt+1)-thetap(jt))		
		dYdhi1(jr,jt) = (Y(jr,jt+1) - Y(jr,jt))/ & 
     &   (thetap(jt+1)-thetap(jt))
	do jt=2,Nt-1
		dXdr(jr,jt) = (X(jr,jt) - X2(jr-1,jt))/ & 
     &   (psig(jr)-psigp(jr-1))
		dYdr(jr,jt) = (Y(jr,jt) - Y2(jr-1,jt))/ & 
     &   (psig(jr)-psigp(jr-1))	
		dXdri1(jr,jt) = (X_i1(jr,jt) - X_i(jr-1,jt))/ & 
     &   (psig(jr)-psigp(jr-1))	
		dYdri1(jr,jt) = (Y_i1(jr,jt) - Y_i(jr-1,jt))/ & 
     &   (psig(jr)-psigp(jr-1))	
		dXdh(jr,jt) = (X_i1(jr,jt) - X_i1(jr,jt-1))/  &  
     &   (thetap_i(jt)-thetap_i(jt-1))		
		dYdh(jr,jt) = (Y_i1(jr,jt) - Y_i1(jr,jt-1))/ & 
     &   (thetap_i(jt)-thetap_i(jt-1))
		dXdhi1(jr,jt) = (X(jr,jt+1) - X(jr,jt))/  &  
     &   (thetap(jt+1)-thetap(jt))		
		dYdhi1(jr,jt) = (Y(jr,jt+1) - Y(jr,jt))/  &      
     &   (thetap(jt+1)-thetap(jt))
	enddo
	jt=Nt
		dXdr(jr,jt) = (X(jr,jt) - X2(jr-1,jt))/ &  
     &   (psig(jr)-psigp(jr-1))
		dYdr(jr,jt) = (Y(jr,jt) - Y2(jr-1,jt))/ &  
     &   (psig(jr)-psigp(jr-1))	
		dXdri1(jr,jt) = (X_i1(jr,jt) - X_i(jr-1,jt))/ &  
     &   (psig(jr)-psigp(jr-1))	
		dYdri1(jr,jt) = (Y_i1(jr,jt) - Y_i(jr-1,jt))/ &  
     &   (psig(jr)-psigp(jr-1))	
		dXdh(jr,jt) = (X_i1(jr,jt) - X_i1(jr,jt-1))/ &  
     &   (thetap_i(jt)-thetap_i(jt-1))		
		dYdh(jr,jt) = (Y_i1(jr,jt) - Y_i1(jr,jt-1))/ &  
     &   (thetap_i(jt)-thetap_i(jt-1))
		dXdhi1(jr,jt) = (X(jr,1) - X(jr,jt))/ &  
     &   (thetap(jt+1)-thetap(jt))		
		dYdhi1(jr,jt) = (Y(jr,1) - Y(jr,jt))/ &  
     &   (thetap(jt+1)-thetap(jt))


		Jcbn = dXdr*dYdh-dXdh*dYdr		! i,j
		Jcbni1 = dXdri1*dYdhi1-dXdhi1*dYdri1  ! i, j+1/2		



				
	do jt=1,Nt
		do jr = 1,Nr-1
			Mdet = dXdr(jr,jt)*dYdh(jr,jt)- &  
     &   dXdh(jr,jt)*dYdr(jr,jt)
			drdX = dYdh(jr,jt)/Mdet     ! dlambda/dR
			dhdX = -dYdr(jr,jt)/Mdet    ! dtheta/dR
			drdY = -dXdh(jr,jt)/Mdet    ! dlambda/dZ
			dhdY = dXdr(jr,jt)/Mdet     ! dtheta/dZ 
			dphdR(jr,jt)=drdX
			dphdZ(jr,jt)=drdY
			gradr(jr,jt) = drdX**2.0+drdY**2.0		 ! ]grad(lambda)]^2
			gradh(jr,jt) = dhdX**2.0+dhdY**2.0		 ! grad(theta)^2
			grt(jr,jt)=drdX*dhdX+drdY*dhdY        !   grad(theta)*grad(lambda)  ! i,j

			Mdet = dXdr2(jr,jt)*dYdh2(jr,jt)- &  
     &   dXdh2(jr,jt)*dYdr2(jr,jt)
			drdX = dYdh2(jr,jt)/Mdet
			dhdX = -dYdr2(jr,jt)/Mdet
			drdY = -dXdh2(jr,jt)/Mdet
			dhdY = dXdr2(jr,jt)/Mdet
		gradr2(jr,jt) = drdX**2.0+drdY**2.0	
		gradh2(jr,jt) = dhdX**2.0+dhdY**2.0
		grt2(jr,jt) = drdX*dhdX+drdY*dhdY            ! i+1/2, j
			
	Mdet = gradr2(jr,jt)*gradh2(jr,jt)- &     
     &   grt2(jr,jt)**2.0
	ghht2(jr,jt) = gradr2(jr,jt)/Mdet
	
			Mdet = dXdri(jr,jt)*dYdhi(jr,jt)- & 
     &   dXdhi(jr,jt)*dYdri(jr,jt)
			drdX = dYdhi(jr,jt)/Mdet
			dhdX = -dYdri(jr,jt)/Mdet
			drdY = -dXdhi(jr,jt)/Mdet
			dhdY = dXdri(jr,jt)/Mdet
			gradri(jr,jt) = drdX**2.0+drdY**2.0		
			gradhi(jr,jt) = dhdX**2.0+dhdY**2.0		
			grti(jr,jt) = drdX*dhdX+drdY*dhdY        ! i+1/2,j+1/2

			Mdet = dXdri1(jr,jt)*dYdhi1(jr,jt)- & 
     &   dXdhi1(jr,jt)*dYdri1(jr,jt)
			drdX = dYdhi1(jr,jt)/Mdet
			dhdX = -dYdri1(jr,jt)/Mdet
			drdY = -dXdhi1(jr,jt)/Mdet
			dhdY = dXdri1(jr,jt)/Mdet
			gradri1(jr,jt) = drdX**2.0+drdY**2.0		
			gradhi1(jr,jt) = dhdX**2.0+dhdY**2.0		
			grti1(jr,jt) = drdX*dhdX+drdY*dhdY       ! i, j+1/2

	enddo
	enddo
    

	do jt=1,Nt				
	gradr(1,jt) = 0		
	dphdR(1,jt) = 0		
	dphdZ(1,jt) = 0		
	gradri1(1,jt) = 0		
	gradh(1,jt) = 0		
	gradhi1(1,jt) = 0		
	grt(1,jt) = 0		
	grti1(1,jt) = 0		
	enddo
				
!CDefine major radius
!C	if lrgaspect==1 then
!C	Rmaj = Xgeo*ones(Nr,Nt+1)
!C	else
	do jt=1,Nt
	do jr=1,Nr
	Rmaj(jr,jt)=X(jr,jt) 
	Rmaj2(jr,jt)=X2(jr,jt) 
	Rmaji(jr,jt)=X_i(jr,jt) 
	Rmaji1(jr,jt)=X_i1(jr,jt) 
	enddo
	enddo










!Area at grid points
	jr=1
	jt=1
	dt_el = (thetap_i(Nt+1)-thetap_i(Nt))
	dr_el = (psigp(jr))/2.
	dArea(jr,jt)=Jcbn(jr,jt)*dr_el*dt_el
	do jt=2,Nt
		dt_el = (thetap_i(jt)-thetap_i(jt-1))
		dr_el = (psigp(jr))/2.
		dArea(jr,jt)=Jcbn(jr,jt)*dr_el*dt_el
	enddo

	do jr=2,Nr-1
	jt=1
	dt_el = (thetap_i(Nt+1)-thetap_i(Nt))
	dr_el = (psigp(jr)-psigp(jr-1))
	dArea(jr,jt)=Jcbn(jr,jt)*dr_el*dt_el
	do jt=2,Nt
		dt_el = (thetap_i(jt)-thetap_i(jt-1))
		dr_el =(psigp(jr)-psigp(jr-1))
		dArea(jr,jt)=Jcbn(jr,jt)*dr_el*dt_el
	enddo
	enddo

	jr=Nr
	jt=1
	dt_el = (thetap_i(Nt+1)-thetap_i(Nt))
	dr_el = (psig(jr)-psigp(jr-1))
	dArea(jr,jt)=Jcbn(jr,jt)*dr_el*dt_el
	do jt=2,Nt
		dt_el = (thetap_i(jt)-thetap_i(jt-1))
		dr_el =(psig(jr)-psigp(jr-1))
		dArea(jr,jt)=Jcbn(jr,jt)*dr_el*dt_el
	enddo


!hald grid points
	do jr=1,Nr-1
	jt=1
	dt_el = (thetap_i(Nt+1)-thetap_i(Nt))
	dr_el = (psig(jr+1)-psig(jr))
	dArea2(jr,jt)=Jcbn2(jr,jt)*dr_el*dt_el
	do jt=2,Nt
		dt_el = (thetap_i(jt)-thetap_i(jt-1))
		dr_el =(psig(jr+1)-psig(jr))
		dArea2(jr,jt)=Jcbn2(jr,jt)*dr_el*dt_el
	enddo
	enddo

!	do jr=1,Nr
!	dArea(jr,Nt+1)=dArea(jr,1)
!	enddo

!CdAreas are ok

!CCompute total area, note that this area, in the limit of large Nr,Nt, coincides with pi r^2 if it is
!Ccircle, so is ok!

	Area=0.0
	do jt=1,Nt
	do jr=1,Nr
	Area = Area+dArea(jr,jt)
	enddo
	enddo


!CNow compute arc lengths do fluxes
   
	jr=1
	jt=1
		dArc_rp1(jr,jt)=Jcbn2(jr,jt)*gradr2(jr,jt)*   &    ! i+1/2,j dl^2
     &   (thetap_i(Nt+1)-thetap_i(Nt))/Rmaj2(jr,jt)
		dArc_rpt1(jr,jt)=Jcbn2(jr,jt)*grt2(jr,jt)*  & 
     &   (thetap_i(Nt+1)-thetap_i(Nt))/Rmaj2(jr,jt) ! i+1/2,j dl*dt
	dArc_rm1(jr,jt)=0         ! i-1/2,j
	dArc_tp1(jr,jt)=0  ! i,j+1/2
	dArc_tm1(jr,jt)=0  ! i,j-1/2
	do jt=2,Nt
		dArc_rp1(jr,jt)=Jcbn2(jr,jt)*gradr2(jr,jt)*  & 
     &   (thetap_i(jt)-thetap_i(jt-1))/Rmaj2(jr,jt)
		dArc_rpt1(jr,jt)=Jcbn2(jr,jt)*grt2(jr,jt)*  & 
     &   (thetap_i(jt)-thetap_i(jt-1))/Rmaj2(jr,jt)
		dArc_rm1(jr,jt)=0
		dArc_tp1(jr,jt)=0
		dArc_tm1(jr,jt)=0
	enddo

	do jr=2,Nr
	jt=1
		dArc_rp1(jr,jt)=Jcbn2(jr,jt)*gradr2(jr,jt)* &
     &   (thetap_i(Nt+1)-thetap_i(Nt))/Rmaj2(jr,jt)
		dArc_rm1(jr,jt)=Jcbn2(jr-1,jt)*gradr2(jr-1,jt)* &
     &   (thetap_i(Nt+1)-thetap_i(Nt))/Rmaj2(jr-1,jt)
		dArc_rpt1(jr,jt)=Jcbn2(jr,jt)*grt2(jr,jt)* &
     &   (thetap_i(Nt+1)-thetap_i(Nt))/Rmaj2(jr,jt)
		dArc_rmt1(jr,jt)=Jcbn2(jr-1,jt)*grt2(jr-1,jt)* &
     &   (thetap_i(Nt+1)-thetap_i(Nt))/Rmaj2(jr-1,jt)
		dArc_tp1(jr,jt)=Jcbni1(jr,jt)*gradhi1(jr,jt)* &
     &   (psigp(jr)-psigp(jr-1))/Rmaji1(jr,jt)
		dArc_tm1(jr,jt)=Jcbni1(jr,Nt)*gradhi1(jr,Nt)* &
     &   (psigp(jr)-psigp(jr-1))/Rmaji1(jr,Nt)
		dArc_tpr1(jr,jt)=Jcbni1(jr,jt)*grti1(jr,jt)* &
     &   (psigp(jr)-psigp(jr-1))/Rmaji1(jr,jt)
		dArc_tmr1(jr,jt)=Jcbni1(jr,Nt)*grti1(jr,Nt)* &
     &  (psigp(jr)-psigp(jr-1))/Rmaji1(jr,Nt)
	do jt=2,Nt
		dArc_rp1(jr,jt)=Jcbn2(jr,jt)*gradr2(jr,jt)* &
     &   (thetap_i(jt)-thetap_i(jt-1))/Rmaj2(jr,jt)
		dArc_rm1(jr,jt)=Jcbn2(jr-1,jt)*gradr2(jr-1,jt)* &
     &   (thetap_i(jt)-thetap_i(jt-1))/Rmaj2(jr-1,jt)
		dArc_rpt1(jr,jt)=Jcbn2(jr,jt)*grt2(jr,jt)* &
     &   (thetap_i(jt)-thetap_i(jt-1))/Rmaj2(jr,jt)
		dArc_rmt1(jr,jt)=Jcbn2(jr-1,jt)*grt2(jr-1,jt)* &
     &   (thetap_i(jt)-thetap_i(jt-1))/Rmaj2(jr-1,jt)
		dArc_tp1(jr,jt)=Jcbni1(jr,jt)*gradhi1(jr,jt)* &
     &   (psigp(jr)-psigp(jr-1))/Rmaji1(jr,jt)
		dArc_tm1(jr,jt)=Jcbni1(jr,jt-1)*gradhi1(jr,jt-1)* &
     &   (psigp(jr)-psigp(jr-1))/Rmaji1(jr,jt-1)
		dArc_tpr1(jr,jt)=Jcbni1(jr,jt)*grti1(jr,jt)* &
     &   (psigp(jr)-psigp(jr-1))/Rmaji1(jr,jt)
		dArc_tmr1(jr,jt)=Jcbni1(jr,jt-1)*grti1(jr,jt-1)* &
     &   (psigp(jr)-psigp(jr-1))/Rmaji1(jr,jt-1)
	enddo
	enddo


!CCompute differentials

	do jr=1,Nr-1
	do jt=1,Nt
		ddr(jr,jt) = (psig(jr+1)-psig(jr))
	enddo
	enddo

	do jr=2,Nr-1
	do jt=1,Nt
		ddr_i(jr,jt) =  (psigp(jr)-psigp(jr-1))
	enddo
	enddo



	do jr=1,Nr-1

	jt=1
		dtp(jr,jt) = (thetap(jt+1)-thetap(jt))
		dtm(jr,jt) = (thetap(Nt+1)-thetap(Nt))
		dt_i(jr,jt) = (thetap_i(Nt+1)-thetap_i(Nt))

	do jt=2,Nt-1
		dtp(jr,jt) = (thetap(jt+1)-thetap(jt))
		dtm(jr,jt) = (thetap(jt)-thetap(jt-1))
		dt_i(jr,jt) = (thetap_i(jt)-thetap_i(jt-1))
	enddo

	jt=Nt
	dtp(jr,jt) = (thetap(Nt+1)-thetap(Nt))
	dtm(jr,jt) = (thetap(jt)-thetap(jt-1))
	dt_i(jr,jt) = (thetap_i(jt)-thetap_i(jt-1))

	enddo

	XX=X(1:Nr,1:Nt)
	YY=Y(1:Nr,1:Nt)



	end 






























