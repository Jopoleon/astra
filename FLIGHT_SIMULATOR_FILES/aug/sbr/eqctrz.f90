!C======================================================================|
	subroutine EQCTRZ
!C
!C----------------------------------------------------------------------|
	use  astra2fbe       ! declaration of minimal CPOs
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
	
      namelist / fenix / yesfitcc,dt_adapt,dt_fazt, &
     & kastr2,wallpos,d_j_m,dw_j_m,dw2_j_m,dt_afazt, &
     & vsoldiv,solwidth,transp_variance, &
     & recycl_wall,boron_wall,predep_w, &
     & hmodetransp,lmodetransp,lhmodel, &
     & solmod_dt,tctr_dt,equi_dt, &
     & saves_dt,savep_dt,neocl_dt,trmod_ty,nbieqmix_dt, &
     & torba_dt,simdtmultip,ped_width,chie_chii, &
     & D_chie,D_ped_mult,gs2d_tmin_multip,btipdirec, &
     & CBND3, &
     & CF7, &
     & CFUS3, &
     & CFUS4, & !
     & CNEUT1, &
     & CNEUT2, &
     & CNEUT4, &
     & CDMJ2,ZRD16,ZRD17, &
     & CDMJ3,ZRD13,ZRD18,ZRD19,ZRD20, &
     & CDMJ4,ZRD11, &
     & CDMJ8, &
     & DROUT, &
     & DTOUT, &
     & DPOUT , &
     & TSCALE ,ZRD22, &
     & TPAUSE,ZRD33, &
     & SGNIP,ZRD92X, &
     & SGNBT,pr_clamp,ZRD61,ZRD62,ZRD63, &
     & ABC,IPL,ELONG,TRIAN,SHIFT,UPDWN,BTOR,IBKDW

!	if (TIME.ge.12.) call err_catch_a  !TEST 

!	if (j_counta>1) then
!		write(*,*) 'ne',g22(na1)
!		call nelig(5,jumbo(1:5))
!		write(*,*) 'ne',jumbo
!	stop
!	endif
!	j_counta=j_counta+1


!	write(9913,'(42E25.11)') TIME,x_point_save(1:20,1),x_point_save(1:20,2),plasma_config+0.

	diagzz=1
	ni_p=257
	nj_p=257
	n_p_grid=30
	use_zlim_pot=1 ! for asdex

	if (j_currr.eq.2) then !rampdown at low current!!
		plasma_up=0
	
		include 'dat/bkdwcurrenteq.dat'
	
		if (diagzz.eq.1) then
!			include 'dat/savestuff.dat'
		endif
		return	
	endif

!
write(*,*) 'eqctrz0'

	if (TIME.le.0) then
	ipl_bf_bkdw=0.0001
!initialized namelists
		yesfitcc=1 !if 1 
		dt_adapt=0.4 !time after breakdown to switch to adaptive grid
		dt_fazt=0.001 !time after breakdown to switch to faster solver
		dt_afazt=0.2 !time after breakdown to switch to faster solver
		kastr2=0. !levae 0
		wallpos=2.2 !levae 0
		solmod_dt=0.002
		tctr_dt=0.0001
		equi_dt=1.e-9
		saves_dt=10.
		savep_dt=10.
		neocl_dt=0.002
		trmod_ty=1.
		chie_chii=0.6
		nbieqmix_dt=0.002
		torba_dt=0.005
		simdtmultip=1.
		btipdirec=1   ! 1 is favourable bgradb drift, -1 is unfavourable bgradbdrift
		pr_clamp = 1 ! 1 - clamp pressure to eped scaling, 0 do not do that
		CBND3  =  1. ! D ped multiplier H mode
		CF7    =  0.003 !time average for psep in edgemod
		CFUS3  = 0.000001 !ELM frequency in Hz
		CFUS4  = 1. !Elm amplitude per frequency
		CNEUT1 =  0.4 !pellet pos
		CNEUT2 =  0.3 !pellet width!
		CNEUT4 =  0.1 ! critical shear for sawtooth triggering
		CDMJ2 = 2.  !m - poloidal number
		CDMJ3 = 1.  !n - toroidal number
		CDMJ4 = 0.  !wseed - width/a starting of island width
		CDMJ8  =  1. ! multiplier of taumin for 2D equilibrium calls
		DROUT  =  0.0002 ! radial profile plotting
		DTOUT  =  0.0002 ! time plotting
		DPOUT  =  0.0002 ! savings
		TSCALE =  20. !scale of simulation
		TPAUSE = 150. ! when to pause
		SGNIP  = -1.  ! sign of IP w.r.t. astra convention 
		SGNBT  = 1.   ! sign of BT w.r.t. astra convention
		ZRD61  = 0.1; !ELM area affected in rho toroidal normalized
		ZRD62  = 10.; !ELM resistivity enhancement
		ZRD63  = 0.1; !bootstrap current reduction due to elms
		ZRD92X = -100.; ! diohdt threshold in kA/s for breakdown (Ip = 10 kA)
	ZRD22=1. ! plasma species: 1 D, 2 H, 3 He
	ZRD16 = 0.2 ! IPL threshold in rampdown in MA
	ZRD17 = 1. !error threshold for speeding up equil
		ZRD33 = 300. !pellet speed in m/s
	ZRD13=1. !nbi model
	ZRD18=1. !EC model
	ZRD19=1. !IC model
	ZRD20=1. !pellet model
	
		com_solver=0   ! 0 uses standard spider, 1 uses matrix inversion but computes the matrix, 2 read the matrix from file dat/matrixgigantspid.dat, now done only for 63 points!

	ABC      =  0.6	
	ELONG    =  1.        
	TRIAN    =  0.        
	SHIFT    =  -0.08      
	UPDWN    =  0.0  
	IPL      =   0.1
	BTOR      =   2.5
	IBKDW     =   0.0
!
	
		write(fname,'(a)') 'exp/nml/'//trim(exp_file)
		INQUIRE( FILE=trim(fname), EXIST= file_existence) 
		if (file_existence .eqv. .true.) then		
         open(53, FILE=trim(fname))
         read(53, nml=fenix)
         close(53)
		endif		


	IPLX      =   IPL

!
		!some parameters
			ZRD53 = solmod_dt; !call of solmodinp and neut dt   solmod_dt
			fs_dt_tctrl = tctr_dt; !time step for time control     tctr_dt
			ZRD35  = equi_dt; !Equilibrium call dt            equi_dt 
			ZRD36  = saves_dt;! call of savectrls dt           saves_dt
			ZRD37  = savep_dt;! call of savectrlp dt           savep_dt
			ZRD47  = neocl_dt; ! call of nclass dt             neocl_dt
			ZRD44  = trmod_ty; !transport model type              trmod_ty
			ZRD67  = nbieqmix_dt; !NBI dt eqctrl nelig and mixins nbieqmix_dt
			ZRD90  = torba_dt; !torbeam                        torba_dt
			CV8    = ped_width !pedestal width in normalized rho_tor
			CSOL3  = chie_chii ! chie / chii ratio
			CSOL4  =  D_chie  ! D/chie ratio
			CSOL2  =  D_ped_mult ! D/chie in pedestal
			CDMJ8  =  gs2d_tmin_multip ! taumin multiplier for 2d gs solver calls max
	com_solver=nint(ZRD11)   ! 0 uses standard spider, 1 uses matrix inversion but computes the matrix, 2 read the matrix from file dat/matrixgigantspid.dat, now done only for 63 points!
			max_max_iterb=3 !iterations in eqb when kpr = -2
			max_max_iterj=1 ! in eqa_ax how many time steps to keep current 			
			max_max_iteri=-2 !iterations in sstepon_r  when kpr=-2  constants when kpr=-2 
			ispid_contour=0 !SPIDER does contouring instead of gs choose 1!
	if (com_solver.eq.10) then
			com_solver=0
			ispid_contour=1 !SPIDER does contouring instead of gs choose 1
	endif
	if (com_solver.eq.11) then
			com_solver=1
			ispid_contour=1 !SPIDER does contouring instead of gs choose 1
	endif
	if (com_solver.eq.12) then
			com_solver=2
			ispid_contour=1 !SPIDER does contouring instead of gs choose 1
	endif
	!load limiter ponits!
		open(32,file='exp/equ/'//trim(MACHINE) &
     &  //'/limpnt.dat')	
			read(32,*) n_limz
			do i=1,n_limz
				read(32,*) r_lim(i),z_lim(i)
			enddo
		close(32)
!		write(*,*) n_limz,r_lim,z_lim
!		call err_catch_a
		a_min_b=ABC



	write(*,*) 'ENOB',ELONG

	endif
write(*,*) 'eqctrz1'


!write(*,*) 'eqctrz',plasma_up

!for feqis
	tau_circuit_ef=tau	
	tau_gseq_ef=tau	
	time_astra=time

	activate_coil_ef=1 ! if 0, coil is forced to 0 current
	sign_coil=1. ! sign of coils currents w.r.t. plasma current
	current_limit_ef(:,1)=1.e6 ! 1 is upper, 2 is lower
	current_limit_ef(:,2)=-1.e6 ! 1 is upper, 2 is lower
	force_coil=0. ! where it is 1, forces coil i,i to current of i,j
	force_coil(2,1)=1. ! coil 2  has to have same current of coil 1
	force_coil(3,2)=1. ! coil 3 has to have same current of coil 2
	if (TIME.gt.tbkdw+ZRD95X) then
		force_coil(2,1)=0.
		force_coil(3,2)=1.
	endif
	if (TIME.gt.tbkdw+ZRD94X) then   
		force_coil(2,1)=0.
		force_coil(3,2)=0.
	endif


	use_reduce_circuit=0
	reconnect_circuits=0		
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
	write(*,*) 'eqctrz',raxis_astra,zaxis_astra,rtor,shift
	n_fourier_restab_boundary=5
	use_limiter_astra=1   ! do not use limiter for DEMO
	refit_mode=3   ! if -1 - 1 turn only, 0 - stab method, if 1 - restab with prescribed axis , 3 - full fit like spider but only for eddy currents, 101 - only Z stab
	solve_fix=0   ! if 0 - solve full fix boundary problem, >0 - N pass only, -2 - uses fbe solution 
	execute_plasma=1   ! if 0 - only circuit equations, if 1 - solve plasma gseq too 
	use_zlim_pot=1
	n_of_newton_iterations=150

!factors of dr and dz for initial iterations
	dr_factor_init_astra=0.01
	dz_factor_init_astra=0.01
	fast_mode=0
	psplex_from_fbe=0
	if (TIME.ge.tbkdw+dt_fazt) then 
		fast_mode=1
		solve_fix=5	
	endif









!plasma species
!bypass input and use actual plasma composition
	if (F1(1).ge.F7(1).and.F1(1).ge.F2(1)) ZRD22=1.
	if (F7(1).ge.F1(1).and.F7(1).ge.F2(1)) ZRD22=2.
	if (F2(1).ge.F1(1).and.F2(1).ge.F7(1)) ZRD22=3.
!plasma species
	if (ZRD22.eq.1) then
		AMJ=2.
		ZMJ=1.
	endif
	if (ZRD22.eq.2) then
		AMJ=1.
		ZMJ=1.
	endif
	if (ZRD22.eq.3) then
		AMJ=4.
		ZMJ=2.
	endif
		AMAIN=AMJ
		ZMAIN=ZMJ	
	if (ZRD22.eq.1) NMAIN=F1  !D
	if (ZRD22.eq.2) NMAIN=F7  !H
	if (ZRD22.eq.3) NMAIN=F2  !He

!	write(*,*) 'zm',zrd22,nmain(1),amain(1),zmain(1),ne(1),ni(1)
!	write(*,*) 'qrad',QRADR(ROC)
!	write(*,*) prad(1:na1),pet(1:na1),pe(1:na1)
	
	zibkdw=ibkdw !!!
	zifbey=ifbey
	ipl_threshold=IPLX !Ip at which astra starts to work
	
	if (RTOR+SHIF(1).gt.2.) then
		write(*,*) 'mag axis major radius > 2 m'
		call err_catch_a	
	endif

!!!
	diohdtthreshold=ZRD92X
	ioh1=ioh2 !
	ioh2=ccoil(1)
	diohdt=(ioh2-ioh1)/tau
	plasma_trig=0 !triggers the plasma initiation. 1 just for 1 time step
	icc1(1:12)=icc2(1:12)
	icc2(1:12)=ccoil(1:12)
	dccdt=(icc2-icc1)/tau


	if (TIME.eq.0.) then
		plasma_up=0
		zrd78 = 1e12
	s_adapt=0
	s_fazt=0
	isafazt=0
	tbkdw=1e12
	j_currr=0
	ipl_bf_bkdw=0.001

	endif

	res_trigts06=0
!	resres_oh6=370000./GP2
!	if (TIME.ge.5.423.and.TIME.lt.5.423+0.042) res_trigts06=1
!	write(*,*) 'dioh',diohdt
	if (flightsim.ge.1.and.ZRD91X.lt.0..and.plasma_up.eq.0) then

		if (diohdt.lt.diohdtthreshold &
     & .and.TIME.ge.ZRD88X) j_currr=1

		if (ipl_bf_bkdw.ge.ipl_threshold) plasma_trig=1

		if (plasma_trig.eq.1 &
     & .and.TIME.ge.ZRD88X) then
			write(*,*) 'plasma triggered'
			write(*,*) diohdt,diohdtthreshold
			plasma_up=1
			tbkdw = TIME
			ZRD78 = TIME
!        <SIGNAL_TRAJECTORY name="rts:Act/EV/TSE.dp" rts_id="0x4f810100">
!          <ENTRY_RULE is="None"/>
!          <EXIT_RULE is="Last"/>
!          <EXECUTION_RULE is="Step"/>
!          <REFERENCE>
!            <POINT time="0.968">      <!--  ! here is where the zero of the discharge is defined , ts06 --> 
!              <VALUE  dimension="1" type="uint32_t">0x4f8a0000</VALUE>
!            </POINT>
!            <POINT time="1.01">  <!--  ! here is where the ts07 of the discharge is defined  --> 
!              <VALUE  dimension="1" type="uint32_t">0x4f8e0000</VALUE>
!            </POINT>
!          </REFERENCE>
!        </SIGNAL_TRAJECTORY>

	write(*,*) SHIFT,TRIAN,ELONG,UPDWN,ABC
!assign new geometrical quantities
!	SHIFT=r_pl-RTOR;
!	TRIAN=0.;
!	ELONG=1.;
!	UPDWN=z_pl;
!	ABC=a_min_b;
	write(*,*) SHIFT,TRIAN,ELONG,UPDWN,ABC
!	 		call err_catch_a
		endif
!!
!!!!

	else

		plasma_up=1
	if (ZRD91X.ge.0.) zrd78 = 0.
	if (ZRD91X.ge.0.) tbkdw = 0.


	endif
!	simdtmultip=5. ! 1 ms astra



!	write(*,*) 'zrd78',time,zrd78

!	if (IPEQL.ne.4.) then
!		write(8891,'(33E25.11)') TIME,CCOIL(1:12), &
!     & UPDWN,RTOR+SHIFT,abc
!		else
!		write(8892,'(33E25.11)') TIME,CCOIL(1:12), &
!     & UPDWN,RTOR+SHIFT,abc
!	endif





	if (TIME.ge.tbkdw+dt_adapt-0.001) then 
		isafazt=2 !rect grid slow
	endif


!	write(*,*) 'plasma current',TIME,IPL,IPLFBE


	if (TIME.ge.tbkdw+dt_adapt) then 

!		isafazt=2 !rect grid slow

		IFBEG=1
		s_adapt=1		
		isafazt=0

!		simdtmultip=1. 
	j_counta=j_counta+1		
		ztt0=ztt1
		ztt1=ztt2
		ztt2=UPDWN
		rtt0=rtt1
		rtt1=rtt2
		rtt2=SHIF(1)
		
	!apply time scheme do cdmj8
		dztt1=2*abs(ztt2-ztt0)/TAU/abc*0.001			
		drtt1=2*abs(rtt2-rtt0)/TAU/abc*0.001			
		dztt2=abs(ztt2-2*ztt1+ztt0)/TAU**2./abc*0.001**2.			
		drtt2=abs(rtt2-2*rtt1+rtt0)/TAU**2./abc*0.001**2.			
		if (dztt1+drtt1+drtt2+dztt2.eq.0.) then
			 err_checka=5.
			else
				err_checka=(dztt1+dztt2+drtt1+drtt2)*200.
				err_checka=max(0.1,min(5.,5./err_checka))
				err_checka=1./10.*nint(10.*err_checka)
				if (err_checka.lt.0.16) err_given=0.1
				if (err_checka.ge.0.16) err_given=0.2
				if (err_checka.ge.0.23) err_given=0.25
				if (err_checka.ge.0.4) err_given=0.5
				if (err_checka.ge.0.8) err_given=1.
				if (err_checka.ge.1.8) err_given=2.
				if (err_checka.ge.2.5) err_given=3.
				if (err_checka.ge.3.5) err_given=4.
				if (err_checka.ge.4.5) err_given=5.
		endif	
	
!	write(*,*) 'simdt ',j_counta,drtt1,
!     & dztt1,drtt2,dztt2,err_given
!	write(443,'(33E25.11)') TIME,drtt1,
!     & dztt1,drtt2,dztt2,err_given,UPDWN
!	simdtmultip=err_given
!	simdtmultip=1.

	endif

!	isafazt=0

	if (TIME.ge.tbkdw+dt_adapt+dt_afazt) then 
!		isafazt=3 !for rect grid go back to fast		
		isafazt=1 !for rect grid go back to fast		
!		psplexavg = 0.1 ! 0.01
	err_checkrzp		=abs(PSPLEX-psplexold)
	if (100.*err_checkrzp.ge.ZRD17) then
			simdtmultip=0.5
			else
			simdtmultip=min(5.,2.*simdtmultip)
		endif
	endif
!	psplexold=PSPLEX


!!!!!
	if (TIME.ge.tbkdw+dt_fazt) then 
		s_fazt=1		
	endif
!	write(*,*) isafazt,dt_afazt+dt_adapt+tbkdw,time

!write(*,*) 'taus eqctrz',cdvm7,simdtmultip

	taumin=CDVM7*simdtmultip ! astra trasnport time step, spider circuit eq time step
	taumax=CDVM7*simdtmultip ! astra transport time step, spider circuit eq time step



	if (TIME.le.ITFBE+2*TAU) then
	UEXT=1. !
		psiprima=PSIEXT
	LEXT=1.
	else
		UEXT=(PSIEXT-psiprima)/TAU
		psiprima=PSIEXT
	LEXT=1. !PSPLEX*ROC

		endif
!	write(*,*) 'current ',LL,RR,VV,IPL,PSPLEX,PSIEXT
	if (plasma_up.eq.1) then
		IPL=IPLFBE  !b.c. for current diffusion	


!	call sawmod !sawtooth model for critical shear

!	RR=VINT(1./CC,ROC)/(GP*ABC**2.*SQRT(1.+ELONG**2./2.))
!	VV=UEXT
!	if (TIME.le.ZRD78+3.*TAU) VV=-0.1*diohdt
!	LL=0.2*GP*RTOR*li3r(roc)+LEXT

!	IPL=(LL*IPL+VV*TAU)/(LL+RR*TAU)

	endif
!	write(*,*) 'current ',IPL,PSPLEX,PSIEXT,G22(NA1),ROC,UPDWN

!		write(357,'(555E25.11)') time,diohdt,ipl_bf_bkdw,
!     & psitok(12,1:2),tau,ccoil(11:12),vcoil(1),ccoil(1),
!     & (psitok(1,2)-psitok(1,1))/tau
	if (plasma_up.eq.0.and.j_currr.eq.1) then
!!
!!!!!!
	!calculate psi to coils	
	if (j_readgrid.eq.0) then
		call grid_spid
		nequiz=52
!		call f_rdexf(nequiz)
!		r_pl=rtor
!		z_pl=updwn
		r_pl=rtor+shift
		z_pl=updwn
		a_min_b=abc
		dumm71=abc
		dumm9r=r_pl
		dumm9z=z_pl
	endif
!!!!
!	psitok(1:nequiz,1)=psitok(1:nequiz,2) !!

!	call flux_r_efable(psitok(1:nequiz,2),52,
!     & ipl_bf_bkdw,r_pl,z_pl,j_readgrid,
!     & m_c_e,r_c_e,z_c_e,conduc_cur)

!	if (j_readgrid.eq.0) then
!		psitok(1:nequiz,1)=psitok(1:nequiz,2) !
!	endif 

		include 'dat/bkdwcurrenteq.dat' !!!

!		r_pl=r_c_e(ii_pl)
!		z_pl=z_c_e(jj_pl)
!
!	geom1d(81)=r_pl!
!	geom1d(82)=z_pl

!	if (j_readgrid.eq.0) then
!		call flux_r_efable(psitok(1:nequiz,2),52,
!     & ipl_bf_bkdw,r_pl,z_pl,j_readgrid,
!     & m_c_e,r_c_e,z_c_e,conduc_cur)
!	endif 

	psitok=0. !


	j_readgrid=1
	endif


!for ramp-down, suppose flattop starts at 1. s
!	if (plasma_up.eq.1
!     & .and.TIME-tbkdw.ge.1.
!     & .and.IPL.lt.0.3) then
!		IFBEG=0.
!		s_adapt=0		
!	endif
	if (plasma_up.eq.1 &
     & .and.TIME-tbkdw.ge.1. &
     & .and.IPL.le.ZRD16) then
		j_currr=2	
		ipl_bf_bkdw=IPL
	ZRD53=1.e6
	endif




	cmnd_dioh2s=0
	cmnd_dioh2u=0
	if (TIME.gt.tbkdw+ZRD95X) cmnd_dioh2s=1
	if (TIME.gt.tbkdw+ZRD94X) cmnd_dioh2u=1   !repaired ! this is for broken oh switch of 2022

	if (diagzz.eq.1) then
!			include 'dat/savestuff.dat'
	endif


!




	if (TIME.le.0.)	tau=taumin

!	ztt0=ztt1
!	ztt1=UPDWN
!	rtt0=rtt1
!	rtt1=RTOR+SHIF(1)
	!apply time scheme do cdmj8
!	dztt1=abs(ztt1-ztt0)/abc/dt_eqnew			
!	drtt1=abs(rtt1-rtt0)/abc/dt_eqnew			
!	if (dztt1+drtt1.eq.0.) then
!	else
!	dt_eqnew=max(taumin,min(dtmaxe,dtmaxe/(drtt1+dztt1)))
!	endif	

!	dt_eqnew=taumin

	if (s_fazt.eq.0) then
		dt_eqnew=taumin
	else

		err_checkrzp		=abs(PSPLEX-dumbu(3))+ &
     & abs(UPDWN-dumbu(1))/abc +  &
     & abs(SHIF(1)-dumbu(2))/abc +  &
     & abs(PSIEXT-dumbu(4))
		write(*,*) 'err cjecl',err_checkrzp,zrd17,psplex,updwn,shif(1),psiext,dumbu(1:4)


		dtmaxe=taumin*CDMJ8
		if (100.*err_checkrzp.ge.ZRD17) then
				dt_eqnew=max(TAUMIN,dt_eqnew/1.2)
				else
				dt_eqnew=min(dtmaxe,1.1*dt_eqnew)
			endif

	write(*,*) dtmaxe, taumin, cdmj8,dt_eqnew

	endif

	dteqz=dt_eqnew ! spider GS solver time step
	dteqz2=dteqz

	ipsmk2=flightsim
	eq_cmd=0
	r00=TIME-r01

!	write(*,*) 'time eq',CDMJ8,dtmaxe,dteqz2,eq_cmd,r00,time,r01

	if (r00.lt.dteqz-1.e-6) then
	else
		dteqz2=r00
		eq_cmd=1
		r00=0.
		r01=TIME
		dumbu(1)=UPDWN
		dumbu(2)=SHIF(1)
		dumbu(3)=PSPLEX
		dumbu(4)=PSIEXT
	endif
!	write(*,*) 'time eq',eq_cmd
		
	if (plasma_up.eq.0) eq_cmd=0		
	if (plasma_trig.eq.1) then
		eq_cmd=1		
		r00=0.
		r01=TIME
		dteqz2=dteqz
		dumbu(1)=UPDWN
		dumbu(2)=SHIF(1)
		dumbu(3)=PSPLEX
		dumbu(4)=PSIEXT
	endif

!equivalent alpha power calculation for discharge scenario design
!	call	alphapow_est_demo
!	call	alphapow_est_iter
!

!	write(77771,*) ' '
!	write(77771,*) time,fast_mode,eq_cmd,dteqz2

	if (fast_mode.eq.1) then
		execute_plasma=eq_cmd
		tau_gseq_ef=dteqz2	
	endif

!	write(77771,*) time,execute_plasma,tau_gseq_ef


!write(*,*) 'eqctrz',plasma_up



!write(*,*) 'eqctrz diag',time,psplex,psiext,ne(1),te(1), &
! & ipl,iplfbe,iplx,fp(1),fp(na1),ni(1),ti(1),&
! & elong,trian,updwn,abc,shift,ccoil(1:12),vcoil(1:12)

	write(*,*) 'eqtime',time,fast_mode,execute_plasma,tau_gseq_ef



	return
	end














!C======================================================================|
	subroutine alphapow_est_demo2
!C
!C----------------------------------------------------------------------|
	use parameter_inc
	use const_inc
	use status_inc
	use outcmn_inc

	implicit none

	integer ictrl,i,j

	double precision npedreac,tpedreac
	double precision tslowd,ssdt(na1)
	double precision titemp(na1),netemp(na1)
	double precision agenerate(na1),paaa(na1)
	double precision naaa(na1),patot
	double precision tauer,svdt
		
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
