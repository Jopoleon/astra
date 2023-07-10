!C======================================================================|
	subroutine EQCTRL
!C
!C----------------------------------------------------------------------|

	use parameter_inc
	use status_inc
	use const_inc
	use outcmn_inc
	use fenix_params
	use plasma_state

	implicit none
!	include	'for/parameter.inc'
!	include 'for/const.inc'
!	include 'for/status.inc'
!	include 'for/outcmn.inc'

	integer ictrl
	double precision r_mag,z_mag,icur
	double precision r00,r01,r10,r11
	double precision z00,z01,z10,z11
	double precision dtmaxe,ztt0,ztt1,dztt1
	double precision i00,i11,yfircs(5)
	double precision t000,t001,cpelfreq
	double precision t_eq1,t_eq2
	integer jdonepel
	data t000/0./
	data t001/0./
	data t_eq1/0./
	data t_eq2/0./
	data jdonepel/0/
	data cpelfreq/0./
	save r00,r01,r10,r11
	save z00,z01,z10,z11
	save i00,i11,t_eq1,t_eq2
	save ztt0,ztt1,t001,t000

!	include 'dat/fenixparams.var'

	save r_mag,z_mag,icur
	
	call slowdn
	call tstopp(cv6)
!	call NTMMOD(ZRD64) !Efable not yet in f90
	plasma_up=1
!		write(*,*) 'tend',TEND,CV6,DN(40)
!		if (IPART.ne.1.) DTEQL = ZRD35
!	write(6661,'(222E25.11)') TIME,NE(90),
!     & betp3r(roc),li3r(roc),qdtr(roc),
!     & qtokr(roc),ndeut(90),ntrit(90),
!     & QTOKR(ROC)-QIDWTR(ROC)-QEDWTR(ROC),
!     & TE(1),TE(90),CPEL1+CIMP3,
!     & ZRD73,CDJM8,NE(NA1),F7(na1),F8(na1),
!     & F7(1),F8(1),F7(50),F8(50),qradr(roc),
!     & NALF(1),DF5(10),DF5(30),DF5(50),
!     & IPL,NE(1),CV3,F4(1)


	call smearr(5.e-4,CAR61,work(:,301))
	call smearr(5.e-4,CAR62,work(:,302))
	call smearr(5.e-4,CAR63,work(:,303))



!calculate pellet frequency
	if (CPEL1+CIMP3.gt.0.001) then
		if (jdonepel.eq.0) then
			cpelfreq=1./(TIME-t000)
			t000=TIME
			jdonepel=1
		endif
		else
		jdonepel=0
	endif	

	VPOL=CV7; !TGLF model
	VPOL(NA1-1)=1./cpelfreq
	VPOL(NA1-2)=t000


!	write(878,*) 'exb parameters: ',TIME,
!     & TTRQ(1),UPAR(90:100),ER(90:100)/BTOR,
!     & ER(90:100)/BTOR/AMETR(90:100)/MU(90:100),
!     & 6.9e6*BTOR/sqrt(MRHO(90:100)),XUPAR(90:100),
!     & SLAT(96),NE(96) 



!	write(*,*) 'source parameter'
!	write(*,*) 

!calculate cost functions
	call cost_functions_actuators
	call cost_functions_diags
	call cost_functions_ramps




!	NCNB=11
!	NCTP=1
!	NCNBT=0

	

	return


	end      


	subroutine cost_functions_actuators
	implicit none
	
	end
	
	subroutine cost_functions_diags
	implicit none
	end
	
	subroutine cost_functions_ramps
	implicit none


	end
