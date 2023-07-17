!C======================================================================|
	subroutine EQCTRZT !for demo
!C
!C----------------------------------------------------------------------|

	use fenix_params
	use const_inc
	use status_inc
	use outcmn_inc
	use    astra2fbe , only:  	dr_factor_init_astra,dz_factor_init_astra
	

	implicit none
!	include	'for/parameter.inc'
!	include 'for/const.inc'
!	include 'for/status.inc'
!	include 'for/outcmn.inc'
!	include 'tmp/declar.fnc'

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
	double precision t_eq1,t_eq2

	double precision a_ratio,r_pl,z_pl
	double precision tbkdw,LL,RR
	double precision ipl_threshold
	double precision conduc_cur(130),dumm1
	character*80 fname
	integer ii_pl,jj_pl
	data t_eq1/0./
	data t_eq2/0./
	save t_eq1,t_eq2
	save r00,r01,r10,r11
	save z00,z01,z10,z11
	save i00,i11
	save ztt0,ztt1


!	include 'dat/fenixparams.var'
!!!
!factors of dr and dz for initial iterations
	dr_factor_init_astra=0.1
	dz_factor_init_astra=0.1

	eq_cmd=0		
	dteqz2=TAU
	z00=z01
	z01=UPDWN
	r00=r01
	r01=SHIFT
	z10=abs(z01-z00)/tau/abc
	r10=abs(r01-r00)/tau/abc
	
	if (s_fazt.eq.0) t_eq1=time
	if (s_fazt.eq.0) eq_cmd=1
	if (s_fazt.eq.0) return
	if (TIME.le.ITFBE) return

!gs solver command
	t_eq2=time
	dteqz2=tau*(1.+8.*exp(-z10/0.01-r10/0.01))	
	write(1131,*) 'times',t_eq1,t_eq2
	if ((t_eq2-t_eq1).ge.dteqz2-1.e-6) then
		dteqz2=t_eq2-t_eq1
	write(1131,*) 'cmd',t_eq1,t_eq2,dteqz2,tau
		eq_cmd=1
		t_eq1=time
	endif
	write(1131,*) 'times',eq_cmd,t_eq1,t_eq2,dteqz2

	return
	end
