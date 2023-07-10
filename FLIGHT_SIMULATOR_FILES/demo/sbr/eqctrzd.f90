!C======================================================================|
	subroutine EQCTRZD !for demo
!C
!C----------------------------------------------------------------------|
	use    exchange_with_astra       ! declaration of minimal CPOs
	use fenix_params
	use parameter_inc
	use status_inc
	use const_inc
	use outcmn_inc
	use plasma_state
	use flight_sim_geometrics

	use parameters_a2equil, only: equil_now
	use debugger, only: flightsim


	implicit none
!	include	'for/parameter.inc'
!	include 'for/const.inc'
!	include 'for/status.inc'
!	include 'for/outcmn.inc'
	include 'tmp/declar.fnc'

	integer ictrl,i,j
	double precision r_mag,z_mag,icur,voltaz(15)
	double precision r00,r01,r10,r11
	double precision z00,z01,z10,z11,dteqz
	double precision i00,i11,yfircs(5),dt_eqnew
	double precision dtmaxe,ztt0,ztt1,dztt1
	double precision rtt0,rtt1,drtt1,err_given
	double precision rtt2,ztt2,drtt2,dztt2,err_checka
	double precision ioh1,ioh2,ip_ref,F_dia(nrd)
	double precision icc1(11),icc2(11),dccdt(11)
	double precision psiprima,pressure(nrd)
	integer diagzz
	integer j_change,j_coil
!
	double precision a_min_b

	double precision a_ratio,r_pl,z_pl
	double precision tbkdw,LL,RR
	double precision ipl_threshold
	double precision conduc_cur(130),dumm1
	character*80 fname
	integer ii_pl,jj_pl
	double precision dumm2,dumm3,dumm9,dumm8
	double precision dumm9r,dumm9z,dumm71,dumbu(10)
	double precision r_lim(300),z_lim(300)
	double precision r_p_grid(30,30),z_p_grid(30,30)
	double precision ds_p_grid(30,30)
	double precision t_p_grid(30),rho_p_grid(30)
	integer n_p_grid
	integer n_limz
	common/comst1/ PJK
	double precision PJK(1550)
	double precision PJK2(1550)
!	common /psi_tocoils_bkdw/ psibkdwcoils,nequizzz,
!     & use_zlim_poz,time_astra2
		 integer nequizzz
		 double precision psibkdwcoils(500,2),time_astra2
		 integer use_zlim_poz


!	include 'dat/fenixparams.var'
!!!

	real*8 m_c_e(257,257,130),r_c_e(257),z_c_e(257)
	integer seed
	double precision hoop_force(257,257)
	double precision dhb_f(257,257),err_checkrzp
	integer	ni_p
	integer	nj_p

	double precision psi_c_e(257,257),b_r_e(257,257)
	DOUBLE PRECISION b_z_e(257,257),psplexold,ccuscita
	double precision ccuscita2,timez_0
	DOUBLE PRECISION dist_lim(300),dist_lim_min
	double precision a_cd(8,8),b_cd(8,39)
	double precision c_cd(12,8),d_cd(12,39)
	double precision x_cd(8,2),v_cd(12)
	double precision u_cd(39),gaps0(52)
	double precision br_in(258,3)
	double precision br_out(258),dpc(200)
	double precision flux_in(129,2),dt_afazt
	double precision flux_out(129)
	double precision KCURR(11,11),xk1,xk2
	double precision I_MEAS(11),prezz
	double precision I_reff(11),circoils(15)
	double precision t_diag,t_diagref,demo_gaps(52,3)
	double precision rand
	integer gapnum(6)
	
	integer N_95_POS
	double precision V_95_POS

	!!!!!!!!!!!!!!
	integer i_lim_min,j_counta
	integer j_currr,j_rampdown,j_readgrid !
	integer n_a1,n_b1,n_c1,n_d1,contr_type
	integer n_a2,n_b2,n_c2,n_d2,j_init
	integer neqlp,ntetap
	double precision yrout(256,256),yzout(256,256)
	logical file_existence
	data r00/0./
	data j_init/0/
	data t_diag/0./
	data xk1/0./
	data xk2/0./
	data r01/0./
	data dt_eqnew/1.e-6/
	data ioh2/0./
	data psiprima/0./
	data j_rampdown/0/
	data j_readgrid/0/
	data j_counta/0/
	data ztt0/0./
	data j_change/1/
	data j_coil/1/
	data rtt0/0./
	data timez_0/20./
!	data r_pl/1.65/
!	data z_pl/0.1d0/
!	data dumm9r/1.65/
!	data dumm9z/0.1d0/
!	data dumm71/0.5d0/!
!	data a_min_b/0.5/
	save r00,r01,r10,r11,a_min_b,dumbu,j_init
	save z10,dt_eqnew,dumm71,dt_afazt,I_reff
	save i00,i11,dumm9r,dumm9z,psplexold
	save ztt0,ztt1,dztt1,ztt2,dztt2,timez_0
	save j_change,j_coil,rtt0,rtt1,drtt1,rtt2,drtt2
	save ioh1,ioh2,icc1,icc2,dccdt
	save psiprima,r_pl,z_pl,xk1,xk2
	save tbkdw,j_currr,j_rampdown,j_readgrid
	save m_c_e,r_c_e,z_c_e,r_lim,z_lim
	save n_limz,j_counta,x_cd,u_cd
	save a_cd,b_cd,c_cd,d_cd,gaps0,kcurr
	save br_in,br_out,flux_in,flux_out	
	save n_a1,n_b1,n_c1,n_d1
	save n_a2,n_b2,n_c2,n_d2,contr_type
      namelist / fenix / yesfitcc,dt_adapt,dt_fazt, &
     & kastr2,wallpos,d_j_m,dw_j_m,dw2_j_m,dt_afazt,&
     & vsoldiv,solwidth,transp_variance,&
     & recycl_wall,boron_wall,predep_w,&
     & hmodetransp,lmodetransp,lhmodel,&
     & solmod_dt,tctr_dt,equi_dt,&
     & saves_dt,savep_dt,neocl_dt,trmod_ty,nbieqmix_dt,&
     & torba_dt,simdtmultip,ped_width,chie_chii,&
     & D_chie,D_ped_mult,gs2d_tmin_multip,btipdirec,pr_clamp


	tau_circuit_ef=tau	
	tau_gseq_ef=tau	


	time_astra2=time
	time_astra=time
	
	max_max_iterb=10 ! fast iteration limits in Beqb
	cmnd_dioh2s=1
	cmnd_dioh2u=1
	ispid_contour=0

	neqlp=NEQUIL
	ntetap=MEQUIL+1


	yrout(1:neqlp,1:ntetap-1)=equil_now%coord_sys%position%r(1:neqlp,1:ntetap-1) !yrout(neqlp,j)
	yzout(1:neqlp,1:ntetap-1)=equil_now%coord_sys%position%z(1:neqlp,1:ntetap-1) !yrout(neqlp,j)
	yrout(1:neqlp,ntetap)=yrout(1:neqlp,1)
	yzout(1:neqlp,ntetap)=yzout(1:neqlp,1)
	yrout(1:neqlp,ntetap+1)=yrout(1:neqlp,2)
	yzout(1:neqlp,ntetap+1)=yzout(1:neqlp,2)





	activate_coil_ef=1
	sign_coil=1. ! sign of coils currents w.r.t. plasma current
	current_limit_ef(:,1)=1.e6 ! 1 is upper, 2 is lower
	current_limit_ef(:,2)=-1.e6 ! 1 is upper, 2 is lower
	force_coil=0. ! where it is 1, forces coil i,i to current of i,j
!	force_coil(2,1)=1.
!	force_coil(3,2)=1.
	new_equivalence=0

!OH circuit up until entrance of dioh2s and dioh2u later
	use_reduce_circuit=0
	reconnect_circuits=0		
!	if (TIME.le.ITFBE+2*TAU) reconnect_circuits=1		
!	n_equivalence=1
!	new_equivalence(1,1)=1
!	new_equivalence(2,1)=1
!	new_equivalence(3,1)=1
	if (j_init.eq.0) then
	raxis_astra=RTOR+SHIFT
	zaxis_astra=UPDWN
	psi0_astra=FP(1)
	psib_astra=FP(NA1)
	else
	raxis_astra=RTOR+SHIF(1)
	zaxis_astra=UPDWN
	psi0_astra=FP(1)
	psib_astra=FP(NA1)
	endif
	j_init=1


!	if (IPART.eq.1) nonegcurr = 1
!	if (IPART.eq.2) nonegcurr = 0

!	nonegcurr=0

!	write(*,*) raxis_astra,zaxis_astra,ipart,nonegcurr

!
	n_fourier_restab_boundary=5
	refit_mode=0   ! if -1 - 1 turn only, 0 - stab method, if 1 - restab with prescribed axis 
	solve_fix=0   ! if 0 - solve full fix boundary problem, >0 - N pass only, -2 - uses fbe solution 
!	execute_plasma=1   ! if 0 - only circuit equations, if 1 - solve plasma gseq too 
	use_zlim_pot=1
	n_of_newton_iterations=8

!	write(*,*) 'plasma ',plasma_up,IFBEY,ITFBE,TIME,IPART
!	write(*,*) 'plasma params ',PSIEXT,PSPLEX
	res_trigts06=0
	psitok=0.
	nequiz=1
!CREATE TEST
!      call	INTVAR             ! Set exp scalars
			zibkdw=IBKDW
			zifbey=IFBEY

	plasma_up=1.  !for real fbe with plasma

	if (TIME.eq.0) then
		open(32,file='exp/equ/dem_/kcurr.dat')
			do j=1,11
				read(32,*) KCURR(j,1:11)
			enddo
		close(32)
	I_meas=0.
	I_reff=0.
	VCOIL=0.
	voltaz=0.
	circoils=0.
	endif
!	write(*,*) plasma_up,ifbey

	
	
!		write(321,*) TIME,UPDWN
	

	use_zlim_pot=0

	voltaz=0. !for fbe evolution with plasma test
!Position control
! Input Zcur, output V_coil12
! x(k+1)=0.8635*x(k)+1*Zcur(k)
! V12(k)=3.0371e4*x(k)-2.8571e5*zcur(k)
!	solve_fix=0
!	if (TIME.ge.0.01) solve_fix=1	
		psplex_from_fbe=0

	if (TIME.ge.0.05) then
		solve_fix=15
	endif
	if (TIME.ge.0.524) then
		solve_fix=0
		psplex_from_fbe=1
	endif
		fast_mode=0
		s_fazt=0		
	if (TIME.ge.0.55) then
!		solve_fix=-2
!		fast_mode=1
!		s_fazt=1		
	endif
!	if (TIME.ge.ITFBE+dt_fazt) then 
!		s_fazt=1		
!	endif
		s_fazt=0
	
!	rhoedge=roc
!	qedge=1./mu(na1)
	
	cmnd_dioh2s=-1
	cmnd_dioh2u=-1
	
	dt_eqnew=TAU
	ipsmk2=flightsim
!		eq_cmd=1		
!	write(*,*) 'currents',TIME, VCOIL(1:15),CCOIL(1:15)
	circoils(1:11)=CCOIL(1:11)
	circoils(12)=sum(CCOIL(7:8))-sum(CCOIL(9:10))
	circoils(13)=CCOIL(12)
	circoils(14)=CCOIL(14)
	circoils(15)=CCOIL(15)
!	write(444,'(146E25.12)') TIME,pjk(1:115)*1e3,
!     & vcoil(1:15),circoils(1:15)  !voltages and currents


!	write(448,'(3E25.11)') TIME,UPDWN,geom1d(98) !vertical position

	if (TIME.le.TAU) then
		open(32,file='exp/equ/dem_/bvec_in.txt')
			read(32,*) i	
		do i=1,258
			read(32,*) br_in(i,1:3)		
		enddo
		close(32)
		open(32,file='exp/equ/dem_/flux_in.txt')
			read(32,*) i		
		do i=1,129
			read(32,*) flux_in(i,1:2)		
		enddo
	endif

	if (nint(IPEQL).ne.4) then
		if (TIME.gt.ITFBE+TAU) then
		if (TIME.lt.ITFBE+5*TAU) then
!	if (plasma_up.eq.0)	call &
!     & find_demo_meas2021vacuum(258,129,br_in,br_out, &
!     & flux_in,flux_out,pjk(1:115),115,1)
!	if (plasma_up.eq.1)	call find_demo_meas2021(258,129, &
!     & br_in,br_out, &
!     & flux_in,flux_out)
!			else
!	if (plasma_up.eq.0)	call &
!     & find_demo_meas2021vacuum(258,129,br_in,br_out, &
!     & flux_in,flux_out,pjk(1:115),115,0)
!	if (plasma_up.eq.1)	call find_demo_meas2021 &
!     & (258,129,br_in,br_out, &
!     & flux_in,flux_out)
		endif
		endif
	else
		if (TIME.gt.ITFBE+TAU) then
!			call find_demo_meas2021_equilef(258,129,br_in,br_out, &
!     & flux_in,flux_out,115,pjk(1:115),dpc(1:115))
!		 geom1d(98)=z_curr_2d
		endif
	endif	
!	write(445,'(388E25.12)') TIME,br_out,flux_out !measurements
!	write(449,'(45E25.12)') TIME,geom1d(51:94) !gaps

	DTEQL=0.
	if (TIME.le.ITFBE-0.1) DTEQL=0.1


	write(6611,'(5555E25.11)') TIME,ELONG,TE(1),NE(1),TE(90),NE(90),TI(1), &
     & QDTR(ROC),QRADR(ROC),CRAD3,SHIF(1),ZEF(1),VOLUME, &
     & BETANR(ROC),ZRD21,he(80),xi(80),cmhd2,1./mu(50), &
		 LI3R(ROC),V_95_POS(1./MU(1:na1))
	DTEQL=0.5
!	if (TIME.le.ITFBE-2.) DTEQL=1.







!save results
	t_diagref=0.1
	
!		open(32,file='exp/equ/dem_/demo_gaps.data')
!	read(32,*) i
!!	write(*,*) i_gaps
!		do i=1,i
!			read(32,*) demo_gaps(i,1),
!     & demo_gaps(i,2),demo_gaps(i,3)
!		enddo
!		close(32)		

!	call new_gaps_calc(52,demo_gaps,75,geom2d(51,1:75,3),
!     & geom2d(51,1:75,4),
!     & geom1d(94-52+1:94))

	pressure(1:na1)=1602.* &
     & (ne(1:na1)*te(1:na1)+ &
     & ni(1:na1)*ti(1:na1)+pfast(1:na1))
	f_dia(1:na1)=ipol(1:na1)*btor*rtor;
		 
	if (t_diag.ge.t_diagref.and.TIME.ge.0.02) then
		if (IPEQL.eq.4) then
			write(9845,'(34443E25.11)') TIME, &
     & circoils(1:15),pjk(1:15)*1e3, &
     & pjk(16:115)*1e3,vcoil(1:15),UPDWN,geom1d(98) &
     & ,br_out,flux_out,geom1d(94-52+1:94), &
     & geom1d(299:300) &
     & ,1./MU(NA1-1),yrout(neqlp,1:ntetap), &
     & yzout(neqlp,1:ntetap),ipl,geom1d(97), &
     & RTOR+SHIF(1),RTOR,ABC,ELON(NA1),TRIA(NA1), &
     & VOLUM(NA1),FP(NA1),FP(1),PSPLEX,PSIEXT,dpc(1:115), &
     & FP(1:NA1),pressure(1:na1),f_dia(1:na1)
	write(5353,'(40E25.11)') TIME,u_cd(1:39)
!			call wrd_equilef
!			call wrd_equilef_pbe
		else
			write(98461,'(34443E25.11)') TIME, &
     & circoils(1:15),pjk(1:15)*1e3, &
     & pjk(16:115)*1e3,vcoil(1:15),UPDWN,geom1d(98) &
     & ,br_out,flux_out,geom1d(94-52+1:94), &
     & geom1d(299:300) &
     & ,1./MU(NA1-1),yrout(neqlp,1:ntetap), &
     & yzout(neqlp,1:ntetap),ipl,geom1d(97), &
     & RTOR+SHIF(1),RTOR,ABC,ELON(NA1),TRIA(NA1), &
     & VOLUM(NA1),FP(NA1),FP(1),PSPLEX,PSIEXT, &
     & (psibkdwcoils(1:115,2)- &
     & psibkdwcoils(1:115,1))/tau*GP2, &
     & FP(1:NA1),pressure(1:na1),f_dia(1:na1)
	write(53531,'(40E25.11)') TIME,u_cd(1:39)
!			call wrd
!			call wrd_pbe
		endif
	t_diag=0.
	endif
	t_diag=t_diag+tau
	
!	write(97543,'(3E25.11)') time,r_curr_2d,  &
!     & z_curr_2d
	
!	TPAUSE=ITFBE-0.01
! plasma disturbances
!1) transport fluctuations + beta drop over taue
!  
	CSCL4=1.
!	if (time.ge.0.01) then  !random transport fluctuations
!		call random_number(CSCL4)
!		CSCL4=1.+0.1*(-0.5+CSCL4)
!	endif

!	if (time.ge.20.) then  !beta drop over tau_E
!		CSCL4=0.5+0.5*exp(-(time-20.)/4.)
!	endif
!


	
	
!	open(32,file='fort.12345')
!		write(32,*) geom2d(1:51,1:75,3), &
!     & geom2d(1:51,1:75,4), &
!     & geom2d(1:51,1:75,1)
!	close(32)
	
	
!	write(*,*) 'li3, betapol : ',geom1d(299:300),fp(na1)-fp(1),perim(na1)
!	open(32,file='fort.z')
!	do i=1,NA1
!	prezz=NE(i)*TE(i)+NI(i)*TI(i)
!	write(32,'(4E25.11)') FP(i),1./MU(i),prezz,IPOL(i)*RTOR*BTOR
!	enddo



!	close(32)

	I_meas(1)=CCOIL(1)
	I_meas(2)=CCOIL(2)
	I_meas(3)=CCOIL(3)
	I_meas(4)=CCOIL(4)
	I_meas(5)=CCOIL(5)
	I_meas(6)=CCOIL(6)
	I_meas(7)=CCOIL(7)
	I_meas(8)=CCOIL(8)
	I_meas(9)=CCOIL(9)
	I_meas(10)=CCOIL(10)
	I_meas(11)=CCOIL(11)

	if (TIME.le.ITFBE) then
		open(32,file='ssmA_demo.dat')
		do j=1,8
		do i=1,8
		read(32,*) a_cd(i,j)
		enddo
		enddo
		close(32)
		open(32,file='ssmB_demo.dat')
		do j=1,39
		do i=1,8
		read(32,*) b_cd(i,j)
		enddo
		enddo
		close(32)
		open(32,file='ssmC_demo.dat')
		do j=1,8
		do i=1,12
		read(32,*) c_cd(i,j)
		enddo
		enddo
		close(32)
		open(32,file='ssmD_demo.dat')
		do j=1,39
		do i=1,12
		read(32,*) d_cd(i,j)
		enddo
		enddo
		close(32)
		x_cd=0.
	endif	

	gapnum(1)=16
	gapnum(2)=23
	gapnum(3)=30
	gapnum(4)=36
	gapnum(5)=46
	gapnum(6)=47
!	if (TIME.le.ITFBE+0.002) then
!		IPL=IPLX
!	else
	if (TIME.le.ITFBE+1.01) then
!		IPLX=IPL
		IPL=IPLX
		IPLFBE=IPL
	endif
		IPL=IPLFBE  !b.c. for current diffusion	
!	endif

	if (TIME.ge.ITFBE+100000.) s_fazt = 1  ! try fast mode


	if (TIME.ge.ITFBE+0.035) then  !when control starts

!references






	if (TIME.lt.ITFBE+5.) then
		do i=1,6
		u_cd(i)=geom1d(94-52+1+gapnum(i)-1)
		enddo
	endif
!		u_cd(2)=geom1d(94-52+1+gapnum(2)-1)
!		u_cd(3)=geom1d(94-52+1+gapnum(3)-1)
!		u_cd(5)=geom1d(94-52+1+gapnum(5)-1)
!		u_cd(6)=geom1d(94-52+1+gapnum(6)-1)



!feedbacks
		do i=1,6
		u_cd(19+i)=geom1d(94-52+1+gapnum(i)-1)
		enddo
		u_cd(26:36)=I_meas(1:11)*1e3
		u_cd(37)=geom1d(98)
		u_cd(38)=CCOIL(12)*1e3
		u_cd(39)=IPL*1.e6

!control algorithm
	do i=1,8
		x_cd(i,2)=sum(a_cd(i,:)*x_cd(:,1))+ &
     & sum(b_cd(i,:)*u_cd)
	enddo
	do i=1,12
		v_cd(i)=sum(1e-0*c_cd(i,:)*x_cd(:,1))+ &
     & sum(d_cd(i,:)*u_cd)
	enddo
		x_cd(:,1)=x_cd(:,2)
		voltaz(1:11)=1e0*v_cd(1:11)
		voltaz(12:12)=v_cd(12)
		
		do i=1,12
			voltaz(i)=max(-1.e6,voltaz(i))
			voltaz(i)=min(1.e6,voltaz(i))
		enddo
		
	else
!references at 0 time
		do i=1,6
		u_cd(i)=geom1d(94-52+1+gapnum(i)-1)
		enddo
		I_reff(1:11)=CCOIL(1:11)   !reference
		u_cd(7:17)=I_reff(1:11)*1e3
		u_cd(18)=geom1d(98)
		ip_ref=IPLX
		u_cd(19)=ip_ref*1.e6
		x_cd=0.
	endif


	voltaz(13:15)=0.


!	VCOIL(1:5)=voltaz(1:5)
!	VCOIL(6)=voltaz(6)
!	VCOIL(11)=voltaz(11)
!	VCOIL(7)=voltaz(7)+voltaz(15)
!	VCOIL(8)=voltaz(8)+voltaz(15)
!	VCOIL(9)=voltaz(9)-voltaz(15)
!	VCOIL(10)=voltaz(10)-voltaz(15)
!	VCOIL(12)=voltaz(12)
!	VCOIL(13)=-voltaz(12)
!	VCOIL(14:15)=voltaz(13:14)










	return
	end












	subroutine new_gaps_calc(ngap,npol,demo_gaps,R, &
     & Z,gap)

	implicit none
	integer ngap,npol
	double precision demo_gaps(ngap,3)
	double precision R(npol),Z(npol)
	double precision gap(ngap)
	
	integer i,j,k
	double precision x1,x2,x3,x4,z1,z2,z3,z4
	double precision t1,t2,t3,t4,d1(npol),d2(npol),d3(npol)

	do i=1,ngap
		do j=1,npol	 
			call find_angle_ef(demo_gaps(i,1),demo_gaps(i,2), &
     & r(j),z(j),d1(j))
		enddo
	
	
	
	enddo

	return
	end




!C======================================================================|
	subroutine alphapow_est_demo
!C
!C----------------------------------------------------------------------|

	use parameter_inc
	use const_inc
	use status_inc
	use outcmn_inc
!	use declar_fnc
!	use declar_fml

	implicit none
!	include	'for/parameter.inc'
!	include 'for/const.inc'
!	include 'for/status.inc'
!	include 'for/outcmn.inc'
	include 'tmp/declar.fnc'
	include 'tmp/declar.fml'

	integer ictrl,i

	double precision npedreac,tpedreac
	double precision tslowd,ssdt(na1)
	double precision titemp(na1),netemp(na1)
	double precision agenerate(na1),paaa(na1)
	double precision naaa(na1),patot
		
	npedreac=7.
	tpedreac=5.5

	tslowd=1.*tauer(roc)  !check for DEMO the factor
	titemp=ti(1:na1)
	ictrl=nint(0.9*na1)
	
	do j=1,na1
		ti(j)=te(j)/te(ictrl)*tpedreac
		include 'fml/svdt'
		ssdt(j)=svdt
	enddo
		ti(1:na1)=titemp
	netemp=ne(1:na1)/ne(ictrl)*npedreac

	agenerate=1.*netemp**2.*ssdt   !scale to have nalph ~ 2% at steady state

!time advance

	do j=1,na1
		naaa(j)=1./(1./TAU+1./tslowd)* &
     & (agenerate(j)+naaa(j)/TAU)
	enddo


!sawtooth effect
	if (ZRD12.eq.1) then
	
	
	endif

	write(*,*) 'calphas(0): ',naaa(1)/ne(1)

	paaa=naaa/tslowd*5.3*1.  !scale
	patot=0.
	do j=2,na1
		patot=patot+paaa(j-1)*VR(j)*HRO
	enddo

	write(*,*) 'alpha power: ',paaa(1),patot

	return
	end



!		ztt1=0.7778*ztt0+1.*(ztt2-geom1d(98))
!		ccuscita=min(6e3,max(-6e3,-3.1941e6*ztt0+
!     & 1.4403e7*(ztt2-geom1d(98))))
!		ztt0=ztt1
!		ccuscita2=min(3e3,max(-3e3,5e4*(rtt2-geom1d(97))))
!		VCOIL(1:11)=0.
!		VCOIL(7:8)=ccuscita+1.*ccuscita2
!		VCOIL(9:10)=-ccuscita+1.*ccuscita2	

