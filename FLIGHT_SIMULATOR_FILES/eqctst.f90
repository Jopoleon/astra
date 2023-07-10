!C======================================================================|
	subroutine EQCTST !copy this to sbr/
!C
!C----------------------------------------------------------------------|
	use  exchange_with_astra       ! declaration of minimal CPOs
	use parameter_inc
	use const_inc
	use status_inc
	use outcmn_inc
	use plasma_state
	use fenix_params
	use fs_coupling_variables, only: fs_dt_tctrl
	use parameters_a2equil, only: equil_now
	use debugger, only: flightsim
	
	implicit none

	integer ictrl,i,j
	double precision r_mag,z_mag,icur
	double precision r00,r01,r10,r11
	double precision z00,z01,z10,z11,dteqz
	double precision i00,i11,yfircs(5),dt_eqnew
	double precision dtmaxe,ztt0,ztt1,dztt1
	double precision rtt0,rtt1,drtt1,err_given
	double precision rtt2,ztt2,drtt2,dztt2,err_checka
	double precision ioh1,ioh2,li3r
	double precision icc1(12),icc2(12),dccdt(12)
	double precision psiprima
	integer  diagzz
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
	double precision r_lim(200),z_lim(200)
	double precision r_p_grid(30,30),z_p_grid(30,30)
	double precision ds_p_grid(30,30)
	double precision t_p_grid(30),rho_p_grid(30)
	integer n_p_grid
	integer n_limz

!!!

	real*8 m_c_e(257,257,130),r_c_e(257),z_c_e(257)

	double precision hoop_force(257,257),dt_afazt
	double precision dhb_f(257,257),err_checkrzp
	integer	ni_p
	integer	nj_p

	double precision psi_c_e(257,257),b_r_e(257,257)
	DOUBLE PRECISION b_z_e(257,257),psplexold,jumbo(5)
	DOUBLE PRECISION dist_lim(200),dist_lim_min
	double precision v_95_pos

	!!!!!!!!!!!!!!
	integer i_lim_min,j_counta,j_init
	integer j_currr,j_rampdown,j_readgrid !
	logical file_existence
	data r00/0./
	data j_init/0/
	data r01/0./
	data dt_eqnew/1.e-6/
	data ioh2/0./
	data psiprima/0./
	data j_rampdown/0/
	data j_readgrid/0/
	data j_counta/0/
!	data r_pl/1.65/
!	data z_pl/0.1d0/
!	data dumm9r/1.65/
!	data dumm9z/0.1d0/
!	data dumm71/0.5d0/
!	data a_min_b/0.5/
	save r00,r01,r10,r11,a_min_b,dumbu,j_init
	save z10,dt_eqnew,dumm71,dt_afazt
	save i00,i11,dumm9r,dumm9z,psplexold
	save ztt0,ztt1,dztt1,ztt2,dztt2
	save rtt0,rtt1,drtt1,rtt2,drtt2
	save ioh1,ioh2,icc1,icc2,dccdt
	save psiprima,r_pl,z_pl
	save tbkdw,j_currr,j_rampdown,j_readgrid
	save m_c_e,r_c_e,z_c_e,r_lim,z_lim
	save n_limz,j_counta
	


	use_zlim_pot=1 ! for asdex

	tau_circuit_ef=tau	
	tau_gseq_ef=tau	
	time_astra=time

	activate_coil_ef=1 ! if 0, coil is forced to 0 current
	sign_coil=1. ! sign of coils currents w.r.t. plasma current
	current_limit_ef(:,1)=1.e6 ! 1 is upper, 2 is lower
	current_limit_ef(:,2)=-1.e6 ! 1 is upper, 2 is lower
	force_coil=0. ! where it is 1, forces coil i,i to current of i,j
!	force_coil(2,1)=1. ! coil 2  has to have same current of coil 1
!	force_coil(3,2)=1. ! coil 3 has to have same current of coil 2


	use_reduce_circuit=0
	reconnect_circuits=0		
	raxis_astra=RTOR+SHIFT
	zaxis_astra=UPDWN
	psi0_astra=FP(1)
	psib_astra=FP(NA1)
	n_fourier_restab_boundary=5
	use_limiter_astra=1   ! do not use limiter for DEMO
	refit_mode=0   ! if -1 - 1 turn only, 0 - stab method, if 1 - restab with prescribed axis , 3 - full fit like spider but only for eddy currents, 101 - only Z stab
	solve_fix=0   ! if 0 - solve full fix boundary problem, >0 - N pass only, -2 - uses fbe solution 
	execute_plasma=1   ! if 0 - only circuit equations, if 1 - solve plasma gseq too 
	use_zlim_pot=1
	n_of_newton_iterations=150

!factors of dr and dz for initial iterations
	dr_factor_init_astra=1.
	dz_factor_init_astra=1.
	fast_mode=0
	psplex_from_fbe=0
	zibkdw=ibkdw !!!
	zifbey=ifbey
	ipl_threshold=IPLX !Ip at which astra starts to work
	
	res_trigts06=0

		dt_eqnew=taumin
	dteqz=dt_eqnew ! spider GS solver time step
	dteqz2=dteqz

	eq_cmd=1
		execute_plasma=eq_cmd
		tau_gseq_ef=dteqz2	


	if (MACHINE(1:3).eq.'aug') then
		cur_init(1:12)=CCOIL(1:12)/1.e3
		cur_init(13:52)=0.
	endif


	write(*,*) 'eqtime',time,fast_mode,execute_plasma,tau_gseq_ef,ncnb,ncnbt,cur_init(1:12)
	write(*,*) V_95_POS(1./mu(1:na1))
	



	return
	end











