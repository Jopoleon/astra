subroutine PHI_EQ_2d_PBE(Nr, Nt, psin_grid, iplasma, pressure, &
			 ffprimp, pprimp, btor, r0, & 
				Rb, Zb, Rax, Zax, &
    psiax,PSIb, IPOL, XX, YY, PSI, & 
				dpsidv, psin_grid_out, &
    g2, G3, & 
				r_out, r_in, volum, G1, G41, & 
				GRADRO, &
    BMAXT, BMINT, BDB02, BDB0, B0DB2, &
				& FOFB, &
    areat, perim, shif, elon, slat, tria, & 
				& phi,qqsg,thetap_out, &
    rmin, bpcell,bcell,rbp2,jrhoteta, & 
				& solve_fix, li3, betapol,dator,rpol,zpol,psplex)

use pi_grec_vars, only: GPI, GPI2, GPI4
use metric_coefficients_pbe, only: lambda2d,lambda2dp,Z_curr_0D,R_curr_0D
use errors_params

implicit none

integer :: max_iter
double precision, parameter :: muvac=4.e-7*GPI

integer, intent(in) :: Nr, Nt, solve_fix
double precision, intent(in) :: iplasma, R0, btor, li3, betapol
double precision, intent(in) , dimension(Nr) :: pressure, ipol, &
    g1, g2, g3, volum, gradro, bmaxt, bmint, bdb02, bdb0, b0db2, &
    fofb, areat, perim, slat
double precision, intent(in) , dimension(Nt) :: Rb, Zb
double precision, intent(out), dimension(Nr) :: psin_grid_out, r_out, r_in, &
    g41, shif, elon, tria
double precision, intent(out), dimension(Nt) :: thetap_out
double precision, intent(out), dimension(Nr, Nt) :: Psi, rmin, jrhoteta, XX, YY
double precision, intent(out), dimension(Nr, Nt) :: dator,rpol,zpol, bpcell,bcell
double precision, intent(inout) :: PSIb, rax, zax
double precision, intent(inout), dimension(Nr) :: psin_grid, &
    ffprimp, pprimp,dpsidv,PHI,qqsg,rbp2
double precision, intent(inout) :: psiax,psplex

integer :: ip0,ip1,ip2,ip3
double precision :: rhoedge,t1,t2,t3,t4,qedge

integer :: i, i1, i2, j, jr, jt, ji, j_ok, iax, jax, &
    jiter, nan_count
double precision :: X0, Y0, X0o, Y0o, cnorm, &
    UPDWN, yrr, ya, &
    yrmax, yrmin, yzmax, yzmin, yrzmax, yrzmin
double precision, dimension(3) :: xxxx1, yyyy1, pppp1
double precision, dimension(Nr) :: PSIn_gridp, effprimp, epprimp, r,rhot,rhoa
double precision, dimension(Nt+1) :: thetap, thetap_i
double precision, dimension(Nr, Nt) :: dArea, Rmaj2, &
    known_term,  dt_i, psio, &
    dArc_rp1, dArc_rm1, dArc_rpt1, dArc_rmt1, &
    dArc_tp1, dArc_tm1, dArc_tpr1, dArc_tmr1, &
    ddr, ddr_i, dtp, dtm, r_a, &
    gradh, gradr2, gradh2, dArea2, &
    B_R, B_Z, B_T,YY2, r_i1, X_i1, Y_i1
double precision :: PSIsave(256,256),derivs(8)
integer :: i_call_save,k
data i_call_save /0/
	save i_call_save
	save PSIsave


do j=1, Nr
    epprimp(j)  = -GPI4*muvac*pprimp(j)
    effprimp(j) = -GPI4*ffprimp(j)
enddo

if (i_call_save.eq.0) then
	do jr=1, Nr
    do jt=1, Nt
        lambda2d(jr, jt) = (jr - 1)/(Nr - 1.)
        lambda2dp(jr, jt) = (jr - 1+0.5)/(Nr - 1.)
    enddo
	enddo 
!	lambda2dp = lambda2d + 0.5/(Nr - 1.) 
	do j=1, Nt
	    PSI(1: Nr, j) = psin_grid(1: Nr)
	enddo
else
	psi(1:nr,1:nt)=psisave(1:nr,1:nt)
endif


do jr=1, Nr-1
    psin_gridp(jr) = 0.5*(psin_grid(jr+1) + psin_grid(jr))
enddo



	do i=1,nr
		rhoa(i)=(i-1.)/(nr-1.) !full grid
		rhot(i)=(i-1.+0.5)/(nr-1.) !half grid
	enddo


psio = psi  
X0  = Rax
Y0  = Zax
X0o = X0
Y0o = Y0

iax   = 1
jax   = 1
j_ok  = 0


	if (solve_fix.gt.0)	max_iter=solve_fix
	if (solve_fix.le.0)	max_iter=250
	if (solve_fix.eq.-2)	max_iter=250  !from free boundary
		
		


!External iterations

iter_loop: do jiter=1, max_iter

! recalculate psin_grid based on ffprime
    if (jiter >= 2 .and. (iax == 1 .and. jax == 1)) then
        gradh(Nr, 1) = btor*r0
        gradh(Nr-1, 1) = sqrt((btor*r0)**2. - ffprimp(Nr) * (psin_grid(Nr) - psin_grid(Nr-1)) * (psib - psiax))
        do j=Nr-2, 1, -1
            gradh(j, 1) = sqrt(gradh(j+1, 1)**2. - ffprimp(j+1) * &
                (psin_grid(j+2) - psin_grid(j)) * (psib - psiax))
        enddo

        tria(1) = 0.
        do jr=2, Nr ! toroidal flux on full grid
            tria(jr) = tria(jr-1) + gradh(jr-1, 1) * &
                sum(dArea2(jr-1, 1: Nt)/Rmaj2(jr-1, 1: Nt))
        enddo

	rhoedge=sqrt(tria(nr)/(GPI*btor))

        do jr=1, Nr-1 ! safety factor at half grid
            elon(jr) = (tria(jr+1) - tria(jr))/(psin_grid(jr+1) - psin_grid(jr)) / &
                (psib - psiax)
        enddo
	ip0=nr
	ip1=nr-3
	ip2=nr-2
	ip3=nr-1

          t1=rhot(ip3)
          t2=rhot(ip2)
          t3=rhot(ip1)
          t4=rhoa(ip0)

 	  call	b_extrp_ef(t4,t1,t2,t3, ip3,ip2,ip1,ip0,nr,elon,k)
				qedge=elon(nr)

        gradh2(1, 1) = 0.
        do jr=2, Nr ! new psin grid
            gradh2(jr, 1) = gradh2(jr-1, 1) + (2.*jr - 3.)/elon(jr-1)/(Nr - 1.)**2.
        enddo
        psin_grid = 0.9*psin_grid + 0.1*gradh2(:, 1)/gradh2(Nr, 1)

        do jr=1, Nr-1
            psin_gridp(jr) = 0.5*(psin_grid(jr+1) + psin_grid(jr))
        enddo
    endif

    call build_2dgrid(Nr, Nt, Rb, Zb, X0, Y0, &
        lambda2d(1: Nr, 1: Nt), lambda2dp(1: Nr, 1: Nt), psin_grid, &
        psiax, psib, j_ok, psin_gridp, PSI, &
        dArea, Rmaj2, dArea2, &
        dArc_rp1, dArc_rm1, dArc_rpt1, dArc_rmt1, &
        dArc_tp1, dArc_tm1, dArc_tpr1, dArc_tmr1, &
        ddr, ddr_i, dtp, dtm, dt_i, XX, YY,YY2, &
        r_a, thetap, thetap_i, &
        gradh, gradr2, gradh2, & 
				& r_i1,x_i1,y_i1)




!	write(88882,*) 'b',nr,nt,x0,y0,lambda2d(1:nr,1),lambda2dp(1:nr,1), & 
!	& psiax,psib,j_ok,darea(1:nr,1),rmaj2(1:nr,1),darea2(1:nr,1), & 
!	& darc_rp1(1:nr,1),ddr(1:nr,1)






    do jt=1, Nt
        do jr=1, Nr
            known_term(jr, jt) = (effprimp(jr) * dArea(jr, jt)/XX(jr, jt) + &
                XX(jr, jt)*epprimp(jr)*dArea(jr, jt))
        enddo
    enddo

    cnorm = sum(known_term)/iplasma/GPI2/0.4/GPI
    known_term = known_term/cnorm

!	write(88882,*) 'c',cnorm,known_term(1:nr,1)
	
if (solve_fix.eq.-2) then
		call port_from_fbe_ef(Nr,nt,xx,yy,PSI)
		X0=xx(1,1)
		y0=yy(1,1)
		X0o=xx(1,1)
		y0o=yy(1,1)
		psib=psi(nr,1)
		psiax=psi(1,1)
		iax=1
		jax=1	
		j_ok=1


	else

    do jt=1, Nt
        PSI(Nr, jt) = PSIb
    enddo

			call solver_dphi2(PSIb,Nr,Nt, &
     &      known_term,&
     &      dArc_rp1,&
     &      dArc_rm1,&
     &      dArc_rpt1,&
     &      dArc_rmt1,&
     &      dArc_tp1,&
     &      dArc_tm1,&
     &      dArc_tpr1,&
     &      dArc_tmr1,&
     &      ddr,&
     &      ddr_i,&
     &      dtp,&
     &      dtm,dt_i,&
     &      PSI)

!	write(88882,*) PSIb,nr,Nt, &
!     &      sum(known_term), &
!     &      sum(dArea), &
!     &      sum(dArc_rp1), &
!     &      sum(dArc_rm1), &
!     &      sum(dArc_rpt1), &
!     &      sum(dArc_rmt1), &
!     &      sum(dArc_tp1), & 
!     &      sum(dArc_tm1), &
!     &      sum(dArc_tpr1), &
!     &      sum(dArc_tmr1), &
!     &      sum(ddr), &
!     &      sum(ddr_i), &
!     &      sum(dtp), &
!     &      sum(dtm),sum(dt_i), &
!     &      sum(PSI),0., &
!     &      0.


!	write(88882,*) 'psi',psi(1:nr,1),psi(1:nr,10)

	endif

    call find_new_X0Y0(Nr, Nt, PSI, XX, YY, X0, Y0, PSIax, iax, jax,derivs)

    if (iax == 1 .and. jax == 1) then
        j_ok = 1   
    else
        j_ok = 0
    endif

    nan_count = 0
    do jr=1, Nr
        do jt=2, Nt
            if (psi(jr, jt) /= psi(jr, jt)) nan_count = nan_count + 1
        enddo
    enddo
				t1 = sum(abs(psi-psio))/nr/nt


!	write(88882,*) jiter,x0o,x0,y0o,y0,psiax,iax,jax,j_ok,t1


    psio = psi

    if ((abs(X0) >= 20.) .and. (abs(y0) >= 20.)) call err_catch_a

!Check convergence
!    t1 = abs(x0 - x0o)/x0 + abs(y0 - y0o)/x0
    X0o = X0
    Y0o = Y0
    if (t1 < err_fix_boundary .and.jiter>1) then
        EXIT ! Convergence
    endif
enddo iter_loop

	psisave(1:nr,1:nt)=psi(1:nr,1:nt) ! save for next iteration

call build_2dgrid(Nr, Nt,  Rb, Zb,  X0, Y0, &
    lambda2d(1: Nr, 1: Nt), lambda2dp(1: Nr, 1: Nt), &
    psin_grid, psiax, psib, j_ok,  psin_gridp, PSI, &
    dArea, Rmaj2, dArea2, &
    dArc_rp1, dArc_rm1, dArc_rpt1, dArc_rmt1, &
    dArc_tp1, dArc_tm1, dArc_tpr1, dArc_tmr1, &
    ddr, ddr_i, dtp, dtm, dt_i, XX, YY, YY2, &
    r_a, thetap, thetap_i, &
    gradh, gradr2, gradh2, & 
				& r_i1,x_i1,y_i1)

!regrid
call build_2dgrid2(Nr, Nt, psin_grid, &
    XX,YY, Rmaj2, r_a, thetap, thetap_i, gradr2, gradh2, &
    PSI, r0, pressure, btor, ipol, iplasma, &
    G2, G3, areat, perim, volum, G1, GRADRO, &
    BMAXT, BMINT, BDB02, BDB0, B0DB2, FOFB, &
    slat, li3, betapol, &
				& psplex,qedge,rhoedge,bpcell,bcell, & 
				& r_i1,x_i1,y_i1)

dator=darea2
rpol=rmaj2
zpol=yy2


rmin(1: Nr, 1: Nt) = r_a(1: Nr, 1: Nt)

!do jt=1, Nt
!    do jr=1, Nr
!        jrhoteta(jr, jt) = -GPI2*(ffprimp(jr) * 1./XX(jr, jt)/0.4/GPI + &
!            XX(jr, jt) * 1.e-6*pprimp(jr))/cnorm
!    enddo
!enddo
	do jt=1,nt
	do jr=1,nr-1
		 jrhoteta(jr,jt)=( &
     & 0.5*(effprimp(jr)+effprimp(jr+1))* &
     &  1./Rmaj2(jr,jt) &
     &  +Rmaj2(jr,jt) &
     &  * &
     & 0.5*(epprimp(jr+1)+epprimp(jr)))/GPI2/0.4/GPI/cnorm
	enddo
	enddo

	jr=nr-1 !nint(1.*nx)
	Z_curr_0D=sum(jrhoteta(1:jr,:)*YY2(1:jr,:)*darea2(1:jr,:))/ &
     & sum(jrhoteta(1:jr,:)*darea2(1:jr,:))
	R_curr_0D=sqrt(sum(jrhoteta(1:jr,:)*Rmaj2(1:jr,:)**2.*darea2(1:jr,:))/ &
     & sum(jrhoteta(1:jr,:)*darea2(1:jr,:)))



UPDWN = YY(1, 1)

do ji=1, Nr
    i = minloc(yy(ji, 1: Nt), 1)  
    i1 = i - 1
    i2 = i + 1
    if (i == 1 ) i1 = Nt
    if (i == Nt) i2 = 1
    xxxx1(1) = xx(ji, i1)
    xxxx1(2) = xx(ji, i)
    xxxx1(3) = xx(ji, i2)
    yyyy1(1) = yy(ji, i1)
    yyyy1(2) = yy(ji, i)
    yyyy1(3) = yy(ji, i2)
    call polyfitcc_feqis(xxxx1, yyyy1, pppp1)
    yrzmin = -pppp1(2)/(2.*pppp1(1))
    yzmin = pppp1(1)*yrzmin**2.0 + pppp1(2)*yrzmin + pppp1(3)

    i = maxloc(yy(ji, 1: Nt), 1)  
    i1 = i - 1
    i2 = i + 1
    if (i ==  1) i1 = Nt
    if (i == Nt) i2 = 1
    xxxx1(1) = xx(ji, i1)
    xxxx1(2) = xx(ji, i)
    xxxx1(3) = xx(ji, i2)
    yyyy1(1) = yy(ji, i1)
    yyyy1(2) = yy(ji, i)
    yyyy1(3) = yy(ji, i2)
    call polyfitcc_feqis(xxxx1, yyyy1, pppp1)
    yrzmax = -pppp1(2)/(2.*pppp1(1))
    yzmax = pppp1(1)*yrzmax**2.0 + pppp1(2)*yrzmax + pppp1(3)

	i=minloc(xx(ji,1:nt),1)  
	yrmin=xx(ji,i)
	i=maxloc(xx(ji,1:nt),1)  
	yrmax=xx(ji,i)
!    yrmin = MINVAL(xx(ji, :))
!    yrmax = MAXVAL(xx(ji, :))

    yrr = .5*(yrmax + yrmin)
    ya  = .5*(yrmax - yrmin)

    SHIF (ji) = yrr - R0
    ELON(ji) = (yzmax - yzmin)/(yrmax - yrmin)
    TRIA(ji) = (yrr - 0.5d0*(yrzmin + yrzmax))/ya
    r_out(ji) = yrmax
    r_in(ji) = yrmin
enddo

r_out(1) = xx(1, 1)
r_in(1) = xx(1, 1)
ELON(1) = ELON(2)
TRIA(1) = 0.d0
SHIF (1) = XX(1, 1) - R0

G41 = G1 ! to be fixed

rax = X0
zax = Y0

psin_grid_out(1: Nr) = PSI(1: Nr, 1) 

thetap_out(1: Nt) = thetap(1: Nt)

	i_call_save=i_call_save+1

	if (i_call_save.gt.1) i_call_save=1

return
end subroutine PHI_EQ_2d_PBE




subroutine sort_PSI_values(y1, &
     &   Nx)

	implicit none
	
	integer j,i,Nx
	double precision y1(Nx),z(nx)
	double precision dum1(Nx),dum2(Nx),dum3(nx)
	double precision z1,z2,z3,dum4(nx)
	
	dum2=y1
	dum3=z
	do i=1,Nx
 	j=minloc(dum2,1)
		dum1(i)=dum2(j)
		dum4(i)=dum3(j)
		dum2(j)=1.E+08
	enddo
	y1=dum1
	z=dum4
	end




	subroutine port_from_fbe_ef(Nx,nt,xx,yy,PSI)
	use ef_circuit
	
	implicit none
	
	integer nx,nt,i,j,i0,j0
	double precision xx(nx,nt),yy(nx,nt),psi(nx,nt)
	double precision dum1,x1,x2,y1,y2,z1,z2,z3,z4
	double precision r0,z0
	
	do j=1,nt
	do i=1,nx
		r0=xx(i,j)
		z0=yy(i,j)
		i0=floor((r0-rmin)/dr+1.)	
		j0=floor((z0-zmin)/dz+1.)	
	x1=r(i0)
	x2=r(i0+1)
	y1=z(j0)
	y2=z(j0+1)
	z1=-psirz(i0,j0)
	z2=-psirz(i0+1,j0)
	z3=-psirz(i0,j0+1)
	z4=-psirz(i0+1,j0+1)
!bilinear interpolation
	call bilinear_average_ef(x1,x2,y1,y2,r0,z0,z1,z2,z3,z4,psi(i,j))

	enddo
	enddo
	
	

	return
	end
	
	
