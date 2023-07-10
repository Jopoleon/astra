!C======================================================================|
	subroutine EQCITEL
!C
!C----------------------------------------------------------------------|
	use const_inc
	use status_inc
	use parameter_inc
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

	if (TIME.lt.ITFBE-0.1) DTEQL = 1.
	if (TIME.ge.ITFBE-0.1) DTEQL = 0.

	

	

	return


	end      


