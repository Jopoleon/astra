	subroutine EQCTRA


	use parameter_inc
	use const_inc
	use status_inc
	use outcmn_inc
      use fenix_params
	use plasma_state
	use flight_sim_geometrics

	use debugger, only: flightsim
	use exchange_with_astra, only: cur_init

	implicit none

	integer ictrl
	double precision r_mag,z_mag,icur
	double precision r00,r01,r10,r11
	double precision z00,z01,z10,z11
	double precision dtmaxe,ztt0,ztt1,dztt1
	double precision i00,i11,yfircs(5),li3r
	double precision t000,t001,cpelfreq
	integer jdonepel
	data t000/0./
	data t001/0./
	data jdonepel/0/
	data cpelfreq/0./
	save r00,r01,r10,r11
	save z00,z01,z10,z11
	save i00,i11
	save ztt0,ztt1,t001,t000

	save r_mag,z_mag,icur
	
!	call generate_files_feqis('exp/equ/aug_/',13)

	if (IPART.ne.1.) DTEQL = ZRD35

	if (IPART.eq.1.) call read_coils_dat

! diagnostics
	if (plasma_up.eq.1)	call nelig(5,geom1d(72:76))  ! interferometer h1 h5
	if (IPART.eq.1)	call nelig(5,geom1d(72:76))  ! interferometer h1 h5
	
	if (ITFBP.eq.0.) IPLFBE=IPLX
	if (ITFBP.ne.0.) IPLFBE=IPLFBE

!	if (TIME-TSTART.ge.TAU) DTEQL = ZRD35

!	write(*,*) 'zger',geom1d(66:67)
!	write(*,*) 'cv4 ',CV4,cv13

	if (time.lt.itfbe) then
	r00=0.
	r01=0.
	r10=0.
	r11=0.
	z00=0.
	z01=0.
	z10=0.
	z11=0.
	i00=0.
	i11=0.
	endif
	
	ictrl=1

!	write(*,*) 'ctrl , ',icur,r_mag,z_mag,updwn
!
	!write(*,*) 'ifbey eqctr',IFBEY,2


!	write(895,'(333E25.11)') TIME,CCOIL(1:12),VCOIL(1:12)

!	write(555,*) ' '
!	write(555,*) TIME,CPEL1,
!     & CV4,
!     & ZRD73,
!     & CDJM5,
!     & CBND3,
!     & CDJM8,
!     & CDWM6,
!     & CDWM5,
!     & CSCL1,
!     & ZRD71,
!     & CV13,
!     & CSOL1,
!     & CDMJ1,
!     & CDMJ2,
!     & CDMJ3,
!     & CDMJ4,
!     & 'powers ',CAR32(1:42),'voltages ',vcoil(1:10),
!     & 'currs ',geom1d(1:10),'geom ',geom1d(51:77)


!	write(*,*) 'nn', nn(na1),nncl+nnwm

	!edge model for edge transport
!	if (TIME-TSTART.gt.TAU) DTEQL=ZRD35

	NCNB=12
	NCNBT=0

	cur_init=0. !vacuum initial condition
	

!	write(*,*) 'coils currs', geom1d(1:NCNB)
!	write(*,*) 'res', geom1d(51:51+NCTP-1)

	return

      end
      
!======================================================================|
	subroutine read_coils_dat
!
!C----------------------------------------------------------------------|

	use parameter_inc
	use const_inc
	use status_inc
	use outcmn_inc


	implicit none


	character*80 dir_dat
	character*80 dum1
	double precision ddata(7),pjk(12),dum2
	integer dum3
	
	dir_dat='exp/equ/'//trim(MACHINE)//&
     &  '/coil.dat'
     
     	open(32,file=dir_dat)
	read(32,*) dum3,dum3
	read(32,*) dum3
	read(32,*) ddata
	pjk(1)=ddata(7)	
	read(32,*) dum3
	read(32,*) ddata
	pjk(2)=ddata(7)	
	read(32,*) dum3
	read(32,*) ddata
	pjk(3)=ddata(7)	
	read(32,*) dum3
	read(32,*) ddata
	pjk(4)=ddata(7)	
	read(32,*) dum3
	read(32,*) ddata
	read(32,*) dum3
	read(32,*) ddata
	read(32,*) dum3
	read(32,*) ddata
	pjk(5)=ddata(7)	
	read(32,*) dum3
	read(32,*) ddata
	read(32,*) dum3
	read(32,*) ddata
	read(32,*) dum3
	read(32,*) ddata
	pjk(6)=ddata(7)	
	read(32,*) dum3
	read(32,*) ddata
	read(32,*) dum3
	read(32,*) ddata
	pjk(7)=ddata(7)	
	read(32,*) dum3
	read(32,*) ddata
	read(32,*) dum3
	read(32,*) ddata
	pjk(8)=ddata(7)	
	read(32,*) dum3
	read(32,*) ddata
	read(32,*) dum3
	read(32,*) ddata
	read(32,*) dum3
	read(32,*) ddata
	read(32,*) dum3
	read(32,*) ddata
	pjk(9)=ddata(7)	
	read(32,*) dum3
	read(32,*) ddata
	pjk(10)=ddata(7)	
	read(32,*) dum3
	read(32,*) ddata
	pjk(11)=ddata(7)	
	read(32,*) dum3
	read(32,*) ddata
	pjk(12)=ddata(7)	



	close(32)

!	stop
	
	ccoil(1)=pjk(1)
	ccoil(2)=pjk(7)
	ccoil(3)=pjk(8)
	ccoil(4)=pjk(6)
	ccoil(5)=pjk(5)
	ccoil(6)=pjk(4)
	ccoil(7)=pjk(2)-pjk(1)
	ccoil(8)=pjk(3)-pjk(2)
	ccoil(9)=pjk(9)
	ccoil(10)=pjk(10)
	ccoil(11)=pjk(11)
	ccoil(12)=pjk(12)

	ccoil(1:12)=ccoil(1:12)*1.e3

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
