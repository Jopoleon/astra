!C======================================================================|
	subroutine SIMULSIM !for iter
!C
!C----------------------------------------------------------------------|
	use parameter_inc
	use status_inc
	use const_inc
	use outcmn_inc
	use fenix_params
	use debugger, only: flightsim
	use fs_coupling_variables, only: fs_dt_smlk

	implicit none

	double precision time_ext
	
	if (flightsim.eq.1) return
	
	
	    call set_fen_inputs( &
        CPEL1, CIMP3, CV4,  CBND3, &
        CSCL1, ZRD71, ZRD73, CDYM3, &
        CDWM5, CDWM6, &
        CDJM5, CDJM6, CDJM7,  CDJM8, ZRD70, &
        CV13, CSOL1, &
        CHE1, CDHJ1, CDHJ2, &
        CHE3, CDHJ3, CDHJ4, &
        CV3, CDHJ5, &
        CDMJ1, CDMJ2, CDMJ3, CDMJ4, VCOIL(1:12), &
        time_ext, CV6, fs_dt_smlk)
        CV13 = MAX(1., CV13)  ! finite pump speed to avoid NaN

	
	
	return
	end



	subroutine set_fen_inputs(dpel,tpel,dtmidppuff,xemidpuff,ardivpuff, & 
		krdivpuff,krmidpuff,tdtrat,midrecirc,tdtrecpuffrat,herec,xerec,herec2,& 
		wrec,wrec2,pump,transp,rech,rec2,rec3,pec1,pec2,pec3,pnb1,pnb2, & 
		ntm1,ntm2,ntm3,ntm4,vcoil,time,tend,dt)
	
	implicit none
	
	double precision dpel,tpel,dtmidppuff,xemidpuff,ardivpuff, & 
		krdivpuff,krmidpuff,tdtrat,midrecirc,tdtrecpuffrat,herec,xerec,herec2,& 
		wrec,wrec2,pump,transp,rech,rec2,rec3,pec1,pec2,pec3,pnb1,pnb2, & 
		ntm1,ntm2,ntm3,ntm4,vcoil(12),time,tend,dt
	
	
		dpel=300.
		tpel=dpel
		dtmidppuff=0.
		xemidpuff=0.
		ardivpuff=  0.
		krdivpuff=0.
		krmidpuff=0.
		tdtrat=0.5
		midrecirc=0.
		tdtrecpuffrat=0.
		herec=0.
		xerec=0.
		herec2= 0.
		wrec=0.
		wrec2=0.
		pump=100.
		transp=0.6
		rech=0.05
		rec2=0.1
		rec3=0.5
		pec1=25.
		pec2=0.
		pec3=0.
		pnb1=25.
		pnb2=0.
		ntm1=0.
		ntm2=0.
		ntm3=0.
		ntm4=0.
		vcoil=0.
		time=0.
		tend=10000.
		dt=0.005
	
	
	
	
	return
	end
	
