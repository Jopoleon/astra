	subroutine DIAGZ

	use parameter_inc
	use const_inc
	use status_inc
	use outcmn_inc
      use fenix_params
			use plasma_state
			use flight_sim_geometrics

	implicit none

	integer ictrl
	double precision r_mag,z_mag,icur
	double precision r00,r01,r10,r11
	double precision z00,z01,z10,z11,dteqz
	double precision i00,i11,yfircs(5),dt_eqnew
	double precision dtmaxe,ztt0,ztt1,dztt1
	double precision rtt0,rtt1,drtt1
	double precision ioh1,ioh2,qohr
	double precision icc1(12),icc2(12),dccdt(12)
	double precision psiprima,wtozr,li3r,qradr,lintr
	integer  diagzz

	double precision a_ratio,r_pl,z_pl
	double precision tbkdw,LL,RR
	character*80 fname



	integer j_currr,j_rampdown,j_readgrid !



	if (btipdirec.eq.1) then
			if (plasma_up.eq.0) then

				write(765,'(222E25.11)') TIME,geom1d(72)/1.d19,&
     & maxval(te,1),maxval(ne,1),&
     & ccoil(1:12),geom1d(60),geom1d(59),geom1d(53),&
     & wtozr(roc),ipl_bf_bkdw,li3r(roc),1./mu(na1-1),vcoil(1:10),&
     &  psplex,psiext,G22(NA),(FP(NA1)-FP(NA))/HRO,IPOL(NA1),&
     & geom1d(57:58),geom1d(56),geom1d(52),geom1d(81:82),&
     & geom1d(54),geom1d(66),geom1d(67),abc,1./mu(1),&
     & 1./mu(na1),VV,qradr(roc),TRIAN,ZRD15,maxval(ti,1),&
     & CRAD3,car54(1:42),btor,FP(1),FP(NA1),&
     & volume,lintr(roc),car32(551),qohr(roc),car32(552:553)

			else

				write(765,'(222E25.11)') TIME,geom1d(72)/1.d19,&
     & maxval(te,1),maxval(ne,1),&
     & ccoil(1:12),geom1d(60),geom1d(59),geom1d(53),&
     & wtozr(roc),iplfbe,li3r(roc),1./mu(na1-1),vcoil(1:10),&
     &  psplex,psiext,G22(NA),(FP(NA1)-FP(NA))/HRO,IPOL(NA1),&
     & geom1d(57:58),geom1d(56),geom1d(52),geom1d(81:82),&
     & geom1d(54),geom1d(66),geom1d(67),abc,1./mu(1),&
     & 1./mu(na1),UPL(NA1),qradr(roc),TRIAN,ZRD15,maxval(ti,1),&
     & CRAD3,car54(1:42),btor,FP(1),FP(NA1),volume,lintr(roc),&
     & car32(551),qohr(roc),car32(552:553)
!				write(767,'(2222E25.11)') AMETR(1:na1),te(1:na1),
!     & ti(1:na1),ne(1:na1),zef(1:na1),qe(1:na1),
!     & qi(1:na1),1./mu(1:na1),FP(1:na1),peecr(1:na1),pebm(1:na1),
!     & pibm(1:na1),prad(1:na1),peicr(1:na1),piicr(1:na1),
!     & pe(1:na1)-pebm(1:na1)-peecr(1:na1)-peicr(1:na1),
!     & volum(1:na1),rho(1:na1)
!				write(768,'(222E25.11)') TIME,
!     & yrout(neqlp,1:ntetap),yzout(neqlp,1:ntetap)

			endif
	else
			if (plasma_up.eq.0) then

				write(765,'(222E25.11)') TIME,geom1d(72)/1.d19,&
     & maxval(te,1),maxval(ne,1),&
     & -ccoil(1:12),geom1d(60),geom1d(59),geom1d(53),&
     & wtozr(roc),-ipl_bf_bkdw,li3r(roc),1./mu(na1-1),-vcoil(1:10),&
     &  -psplex,-psiext,G22(NA),(FP(NA1)-FP(NA))/HRO,IPOL(NA1),&
     & geom1d(57:58),geom1d(56),geom1d(52),geom1d(81:82),&
     & geom1d(54),geom1d(66),geom1d(67),abc,-1./mu(1),&
     & -1./mu(na1),-VV,qradr(roc),TRIAN,ZRD15,maxval(ti,1),&
     & CRAD3,car54(1:42),btor,FP(1),FP(NA1),volume,lintr(roc),&
     & car32(551),qohr(roc),car32(552:553)

			else

				write(765,'(222E25.11)') TIME,geom1d(72)/1.d19,&
     & maxval(te,1),maxval(ne,1),&
     & -ccoil(1:12),geom1d(60),geom1d(59),geom1d(53),&
     & wtozr(roc),-iplfbe,li3r(roc),1./mu(na1-1),-vcoil(1:10),&
     &  -psplex,-psiext,G22(NA),(FP(NA1)-FP(NA))/HRO,IPOL(NA1),&
     & geom1d(57:58),geom1d(56),geom1d(52),geom1d(81:82),&
     & geom1d(54),geom1d(66),geom1d(67),abc,-1./mu(1),&
     & -1./mu(na1),-UPL(NA1),qradr(roc),TRIAN,ZRD15,maxval(ti,1),&
     & CRAD3,car54(1:42),btor,FP(1),FP(NA1),volume,lintr(roc),&
     & car32(551),qohr(roc),car32(552:553)
!				write(767,'(2222E25.11)') AMETR(1:na1),te(1:na1),
!     & ti(1:na1),ne(1:na1),zef(1:na1),qe(1:na1),
!     & qi(1:na1),1./mu(1:na1),FP(1:na1),peecr(1:na1),pebm(1:na1),
!     & pibm(1:na1),prad(1:na1),peicr(1:na1),piicr(1:na1),
!     & pe(1:na1)-pebm(1:na1)-peecr(1:na1)-peicr(1:na1),
!     & volum(1:na1),rho(1:na1)
!				write(768,'(222E25.11)') TIME,
!     & yrout(neqlp,1:ntetap),yzout(neqlp,1:ntetap)

			endif
	
	
	endif
	

	return
	end
