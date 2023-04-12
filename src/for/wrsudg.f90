	subroutine wrsudg

	use parameter_inc
	use const_inc
	use status_inc
	use outcmn_inc
      use fenix_params
			use flight_sim_geometrics
	use plasma_state
use parameters_a2spider, only: equil_now

	implicit none

	double precision pfusss,pradss,neavss
	double precision coil_forcesR(na1)
	double precision coil_forcesZ(na1)

	double precision betpoledge
	double precision yroutfull(660)
	integer mmequi,j,ntetap
	double precision yrout(256,256),yzout(256,256)
	double precision QRADR,QTOKR,QDTR,QEDWTR,QIDWTR
	double precision WTOZR,LI3R,BETP3R,wer,wtotr
	double precision ne_avvg,qezzo,qizzo,wmhdooo
	double precision wemhdooo
	integer gapnum(6),neqlp
	
	double precision coil_forces(na1,2)
	double precision betp3r_0, LI3R_0, IPL_0, upl_0, &
     & QRADR_0, QTOKR_0, &
     & CRAD4_0, CMHD2_0, CSCL4_0, CDWM1_0, &
     & CDWM2_0, CDJM6_0, CDJM1_0, &
     & CDJM2_0, CDJM4_0, &
     & CDJM7_0, CHE3_0, CDJM3_0, ZRD77_0, &
     & CDMJ5_0, CDMJ6_0, CDMJ7_0, &
     & CDHJ7_0, &
     & te_0(100),ne_0(100), &
     & ccoil_0(10), &
     & geom1d_0(32),time_0
	double precision betp3r_1, LI3R_1, IPL_1, upl_1, &
     & QRADR_1, QTOKR_1, &
     & CRAD4_1, CMHD2_1, CSCL4_1, CDWM1_1, &
     & CDWM2_1, CDJM6_1, CDJM1_1, &
     & CDJM2_1, CDJM4_1, &
     & CDJM7_1, CHE3_1, CDJM3_1, ZRD77_1, &
     & CDMJ5_1, CDMJ6_1, CDMJ7_1, &
     & CDHJ7_1, &
     & te_1(100),ne_1(100), &
     & ccoil_1(10), &
     & geom1d_1(32),time_1,dum1,dum2
	double precision betp3r_now, LI3R_now, IPL_now, upl_now, &
     & QRADR_now, QTOKR_now, &
     & CRAD4_now, CMHD2_now, CSCL4_now, CDWM1_now, &
     & CDWM2_now, CDJM6_now, CDJM1_now, &
     & CDJM2_now, CDJM4_now, &
     & CDJM7_now, CHE3_now, CDJM3_now, ZRD77_now, &
     & CDMJ5_now, CDMJ6_now, CDMJ7_now, &
     & CDHJ7_now, &
     & te_now(100),ne_now(100), &
     & ccoil_now(10), &
     & geom1d_now(32),time_now,time_ext,dt_smlk

	save betp3r_0, LI3R_0, IPL_0, upl_0, &
     & QRADR_0, QTOKR_0, &
     & CRAD4_0, CMHD2_0, CSCL4_0, CDWM1_0, &
     & CDWM2_0, CDJM6_0, CDJM1_0, &
     & CDJM2_0, CDJM4_0, &
     & CDJM7_0, CHE3_0, CDJM3_0, ZRD77_0, &
     & CDMJ5_0, CDMJ6_0, CDMJ7_0, &
     & CDHJ7_0, &
     & te_0,ne_0, &
     & ccoil_0, &
     & geom1d_0,time_0,betp3r_1, LI3R_1, IPL_1, upl_1, &
     & QRADR_1, QTOKR_1, &
     & CRAD4_1, CMHD2_1, CSCL4_1, CDWM1_1, &
     & CDWM2_1, CDJM6_1, CDJM1_1, &
     & CDJM2_1, CDJM4_1, &
     & CDJM7_1, CHE3_1, CDJM3_1, ZRD77_1, &
     & CDMJ5_1, CDMJ6_1, CDMJ7_1, &
     & CDHJ7_1, &
     & te_1,ne_1, &
     & ccoil_1, &
     & geom1d_1,time_1,betp3r_now, LI3R_now, IPL_now, upl_now, &
     & QRADR_now, QTOKR_now, &
     & CRAD4_now, CMHD2_now, CSCL4_now, CDWM1_now, &
     & CDWM2_now, CDJM6_now, CDJM1_now, &
     & CDJM2_now, CDJM4_now, &
     & CDJM7_now, CHE3_now, CDJM3_now, ZRD77_now, &
     & CDMJ5_now, CDMJ6_now, CDMJ7_now, &
     & CDHJ7_now, &
     & te_now,ne_now, &
     & ccoil_now, &
     & geom1d_now,time_now




	double precision qalp1,qrad1,ne1(na1)
	double precision qalp2,qrad2,ne2(na1),magnetics(500)
	integer i_error












      if (MACHINE.eq.'dem_') then
          i_error = 0 ! whether add noise latencies errors to diagnostics
          qalp1=QDTR(ROC)
          qrad1=QRADR(ROC)
          NE1=NE(1:NA1)
          if (i_error .ne. 0) then
              call add_errors(qalp1,qrad1,ne1,qalp2,qrad2,ne2)
              ! QDTR: fusion power
              ! QRAD: radiated power
              ! NE: electron density profile
              qalp1=qalp2
              qrad1=qrad2
              ne1=ne2
          endif

	gapnum(1)=16
	gapnum(2)=23
	gapnum(3)=30
	gapnum(4)=36
	gapnum(5)=46
	gapnum(6)=47
				do j=1,6
					magnetics(j)=geom1d(94-52+1+gapnum(j)-1)
				enddo
					magnetics(7)=geom1d(97) !Rcurr
					magnetics(8)=geom1d(98) !Zcurr
					magnetics(9)=ipl*1.e6   !Ipl A

          call shmw( &
     &        TIME, &
     &        qalp1, qrad1, &
     &        IPL, QTOKR(ROC) - qedwtr(ROC) - qidwtr(ROC), &
     &        CRAD4, CMHD2, CSCL4, CDWM1, &
     &        CDWM2, CDWM7, CDJM1, CDJM2, CDJM3, CDJM4, ZRD77, &
     &        CDMJ5, CDMJ6, CDMJ7, &
     &        CNEUT1, CNEUT2, TE(1:NA1), NE1(1:NA1), &
     &        1d0/MU(1:NA1), &
     &        CU(1:NA1), &
     &     CCOIL(1:15)*1.e3,magnetics(1:439))
      endif










	if (MACHINE.eq.'aug_') then

!       write(*,*) 'ccoil', ccoil(8) 


	geom1d(64)=WTOZR(ROC)*1e6 !total Wmhd including fast ions  in MJ
!	geom1d(72)=sum(NE(1:NA1))/NA1*1e19 !H-1 1019 m-3 !calculated in eqctrl in nelig



	time_1=time-tstart
	time_0=time-tstart-TAU

!	write(*,*) 'wrsudg ',time_0,time_1,TIME-tstart,TAU,ZRD93
	!get coil forces
	if (TIME.ge.ZRD78) then
!		if (nint(IPEQL).ne.4) then
			call coil_force2(coil_forces(1:100, 1:2))
!		else
!			call coil_forces_feqis(100,coil_forces(1:100,1), &
!     & coil_forces(1:100,2)) !
!			coil_forces=-coil_forces
!		endif
	endif

!	write(*,*) 'forcesR',coil_forces(1:21,1)
!	write(*,*) 'forcesZ',coil_forces(1:21,2)
	car54(1:21) = coil_forces(1:21,1)
	car54(22:42) = coil_forces(1:21,2)

	neqlp=abs(nint(NEQUIL))
	ntetap=abs(nint(MEQUIL))+1
	yrout(1:neqlp,1:ntetap-1)=equil_now%coord_sys%position%r(1:neqlp,1:ntetap) !yrout(neqlp,j)
	yzout(1:neqlp,1:ntetap-1)=equil_now%coord_sys%position%z(1:neqlp,1:ntetap) !yrout(neqlp,j)
	yrout(1:neqlp,ntetap)=yrout(1:neqlp,1)
	yzout(1:neqlp,ntetap)=yzout(1:neqlp,1)
	yrout(1:neqlp,ntetap+1)=yrout(1:neqlp,2)
	yzout(1:neqlp,ntetap+1)=yzout(1:neqlp,2)
	


	mmequi=2*ntetap !ntetap
	do 	j=1,ntetap !ntetap
		yroutfull(2*j-1)=yrout(neqlp,j) !yrout(neqlp,j)
		yroutfull(2*j)=yzout(neqlp,j) !yzout(neqlp,j)
	enddo
	if (TIME-TSTART.le.ZRD93+1.e-8) then
!	write(*,*) 'synchronized',TIME,ZRD93,plasma_up,IPL,IPL_0,IPL_1
	if (plasma_up.eq.0) then
	call shmw(TIME-TSTART, &
     & betp3r(roc), LI3R(ROC), ipl_bf_bkdw, upl(na1), &
     & QRADR(ROC), QTOKR(ROC), &
     & CRAD4, CMHD2, CSCL4, CDWM1, &
     & ZRD15, CDJM6, CDJM1, &
     &  CDJM4, &
     & CDJM7, CHE3,CDMJ5, CDMJ6, CDMJ7, &
     & te(1:na1),ne(1:na1), &
     & ZEF(1:NA1),&
     & F4(1:NA1), &
     & coil_forces(1:na1,1),&
     & coil_forces(1:na1,2),&
     & TI(1:na1),&
     & XI(1:na1),&
     & HE(1:na1),&
     & PE(1:na1),&
     & PI(1:na1),&
     & NIBM(1:na1),&
     & PRAD(1:na1),&
     & ccoil(1:10),&
     & geom1d(51:82),&
     & car34(1:12),&
     & yroutfull(1:mmequi),&
     & CDHJ7) !na1+1.d-16  !& geom1d(1:10),
	else
	call shmw(TIME-TSTART,&
     & betp3r(roc), LI3R(ROC), IPLFBE, upl(na1),&
     & QRADR(ROC), QTOKR(ROC),&
     & CRAD4, CMHD2, CSCL4, CDWM1,&
     & ZRD15, CDJM6, CDJM1,&
     & CDJM4,&
     & CDJM7, CHE3,& 
     & CDMJ5, CDMJ6, CDMJ7,&
     & te(1:na1),ne(1:na1),&
     & ZEF(1:NA1),&
     & F4(1:NA1),&
     & coil_forces(1:na1,1),&
     & coil_forces(1:na1,2),&
     & TI(1:na1),&
     & XI(1:na1),&
     & HE(1:na1),&
     & PE(1:na1),&
     & PI(1:na1),&
     & NIBM(1:na1),&
     & PRAD(1:na1),&
     & ccoil(1:10),&
     & geom1d(51:82),&
     & car34(1:12),&
     & yroutfull(1:mmequi),&
     & cdhj7) !na1+1.d-16  !& geom1d(1:10),

	endif

	betp3r_0=betp3r(roc)
	LI3R_0=LI3R(ROC)
	if (plasma_up.eq.0) then
	IPL_0=ipl_bf_bkdw
	else
	IPL_0=IPL
	endif
	upl_0=upl(na1)
	QRADR_0=QRADR(ROC)
	QTOKR_0=QTOKR(ROC)
	CRAD4_0=CRAD4
	CMHD2_0=CMHD2
	CSCL4_0=CSCL4
	CDWM1_0=CDWM1
	CDWM2_0=CDWM2
	CDJM6_0=CDJM6
	CDJM1_0=CDJM1
	CDJM2_0=CDJM2
	CDJM4_0=CDJM4
	CDJM7_0=CDJM7
	CHE3_0=CHE3
	CDJM3_0=CDJM3
	ZRD77_0=ZRD77
	CDMJ5_0=CDMJ5
	CDMJ6_0=CDMJ6
	CDMJ7_0=CDMJ7
	CDHJ7_0=CDHJ7
	te_0(1:na1)=te(1:na1)
	ne_0(1:na1)=ne(1:na1)
	ccoil_0=ccoil(1:10)
	geom1d_0=geom1d(51:82)
!	endif

	else

200	continue
!	write(*,*) 'not synchronized ',time_0,time_1,ZRD93
!	call err_catch_a
	betp3r_1=betp3r(roc)
	LI3R_1=LI3R(ROC)
	if (plasma_up.eq.0) then
	IPL_1=ipl_bf_bkdw
	else
	IPL_1=IPL
	endif
	upl_1=upl(na1)
	QRADR_1=QRADR(ROC)
	QTOKR_1=QTOKR(ROC)
	CRAD4_1=CRAD4
	CMHD2_1=CMHD2
	CSCL4_1=CSCL4
	CDWM1_1=CDWM1
	CDWM2_1=CDWM2
	CDJM6_1=CDJM6
	CDJM1_1=CDJM1
	CDJM2_1=CDJM2
	CDJM4_1=CDJM4
	CDJM7_1=CDJM7
	CHE3_1=CHE3
	CDJM3_1=CDJM3
	ZRD77_1=ZRD77
	CDMJ5_1=CDMJ5
	CDMJ6_1=CDMJ6
	CDMJ7_1=CDMJ7
	CDHJ7_1=CDHJ7
	te_1(1:na1)=te(1:na1)
	ne_1(1:na1)=ne(1:na1)
	ccoil_1=ccoil(1:10)
	geom1d_1=geom1d(51:82)

	time_now=ZRD93
	dum1=(time_1-time_0)**0.5
	dum2=(time_now-time_0)**0.5
	betp3r_now=betp3R_0+(betp3R_1-betp3R_0)/dum1*dum2
	LI3R_now=LI3R_0+(LI3R_1-LI3R_0)/dum1*dum2
	IPL_now=IPL_0 + (IPL_1-IPL_0)/dum1*dum2
	upl_now=upl_0 + (upl_1-upl_0)/dum1*dum2
	QRADR_now=QRADR_0 + (qradr_1-qradr_0)/dum1*dum2
	QTOKR_now=QTOKR_0 + (QTOKR_1-qtokr_0)/dum1*dum2
	CRAD4_now=CRAD4_0 + (CRAD4_1-crad4_0)/dum1*dum2
	CMHD2_now=CMHD2_0 + (CMHD2_1-CMHD2_0)/dum1*dum2
	CSCL4_now=CSCL4_0 + (CSCL4_1-CSCL4_0)/dum1*dum2
	CDWM1_now=CDWM1_0 + (CDWM1_1-CDWM1_0)/dum1*dum2
	CDWM2_now=CDWM2_0 + (CDWM2_1-CDWM2_0)/dum1*dum2
	CDJM6_now=CDJM6_0 + (CDJM6_1-CDJM6_0)/dum1*dum2
	CDJM1_now=CDJM1_0 + (CDJM1_1-CDJM1_0)/dum1*dum2
	CDJM2_now=CDJM2_0 + (CDJM2_1-CDJM2_0)/dum1*dum2
	CDJM4_now=CDJM4_0 + (CDJM4_1-CDJM4_0)/dum1*dum2
	CDJM7_now=CDJM7_0 + (CDJM7_1-CDJM7_0)/dum1*dum2
	CHE3_now=CHE3_0 + (CHE3_1-CHE3_0)/dum1*dum2
	CDJM3_now=CDJM3_0 + (CDJM3_1-CDJM3_0)/dum1*dum2
	ZRD77_now=ZRD77_0 + (ZRD77_1-ZRD77_0)/dum1*dum2
	CDMJ5_now=CDMJ5_0 + (CDMJ5_1-CDMJ5_0)/dum1*dum2
	CDMJ6_now=CDMJ6_0 + (CDMJ6_1-CDMJ6_0)/dum1*dum2
	CDMJ7_now=CDMJ7_0 + (CDMJ7_1-CDMJ7_0)/dum1*dum2
	CDHJ7_now=CDHJ7_0 + (CDHJ7_1-CDHJ7_0)/dum1*dum2
	te_now(1:na1)=te_0(1:na1)+(te_1(1:na1)-te_0(1:na1))/dum1*dum2
	ne_now(1:na1)=ne_0(1:na1)+(ne_1(1:na1)-ne_0(1:na1))/dum1*dum2
	ccoil_now=ccoil_0+(ccoil_1-ccoil_0)/dum1*dum2
	geom1d_now=geom1d_0+(geom1d_1-geom1d_0)/dum1*dum2

!	write(*,*) 'now ', time_now,time_0,time_1,IPL_0,IPL_1,IPL_now
!	if (btipdirec.eq.-1) then
!	call shmw_aug(time_now,
!     & betp3r_now, LI3R_now, -IPL_now*1e6, upl_now,
!     & QRADR_now, QTOKR_now,
!     & CRAD4_now, CMHD2_now, CSCL4_now, CDWM1_now,
!     & CDWM2_now, CDJM6_now, CDJM1_now,
!     & CDJM2_now, CDJM4_now,
!     & CDJM7_now, CHE3_now, CDJM3_now, ZRD77_now,
!     & CDMJ5_now, CDMJ6_now, CDMJ7_now,
!     & CDHJ7_now,
!     & te_now(1:na1),ne_now(1:na1)*1.e19,
!     & -ccoil_now(1:10)*1e3,
!     & geom1d_now(1:32),
!     & na1) !na1+1.d-16  !& geom1d(1:10),
!	else
	call shmw(time_now,&
     & betp3r_now, LI3R_now, IPL_now, upl_now,&
     & QRADR_now, QTOKR_now,&
     & CRAD4_now, CMHD2_now, CSCL4_now, CDWM1_now,&
     & ZRD15, CDJM6_now, CDJM1_now,&
     & CDJM4_now,&
     & CDJM7_now, CHE3_now,&
     & CDMJ5_now, CDMJ6_now, CDMJ7_now,&
     & te_now(1:na1),ne_now(1:na1),&
     & ZEF(1:NA1),&
     & F4(1:NA1),&
     & coil_forces(1:na1,1),&
     & coil_forces(1:na1,2),&
     & TI(1:na1),&
     & XI(1:na1),&
     & HE(1:na1),&
     & PE(1:na1),&
     & PI(1:na1),&
     & NIBM(1:na1),&
     & PRAD(1:na1),&
     & ccoil(1:10),&
     & geom1d(51:82),&
     & car34(1:12),&
     & yroutfull(1:mmequi),&
     & cdhj7) !na1+1.d-16  !& geom1d(1:10),

!	write(9871,*) 'wrsud',time_now,&
!     & betp3r_now, LI3R_now, IPL_now, upl_now,&
!     & QRADR_now, QTOKR_now,&
!     & CRAD4_now, CMHD2_now, CSCL4_now, CDWM1_now,&
!     & ZRD15, CDJM6_now, CDJM1_now,&
!     & CDJM4_now,&
!     & CDJM7_now, CHE3_now,&
!     & CDMJ5_now, CDMJ6_now, CDMJ7_now,&
!     & te_now(1:na1),ne_now(1:na1),&
!     & ZEF(1:NA1),&
!     & F4(1:NA1),&
!    & coil_forces(1:na1,1),&
!     & coil_forces(1:na1,2),&
!     & TI(1:na1),&
!     & XI(1:na1),&
!     & HE(1:na1),&
!     & PE(1:na1),&
!     & PI(1:na1),&
!     & NIBM(1:na1),&
!     & PRAD(1:na1),&
!     & ccoil(1:10),&
!     & geom1d(51:82),&
!     & car34(1:12),&
!     & yroutfull(1:mmequi),&
!     & cdhj7



!	endif
!	call read_input_constant_file(time_ext,dt_smlk)
!	write(*,*) time_now,time_0,time_1,ZRD93
	if (TIME-TSTART.le.ZRD93+1.e-8) goto 300
!	write(*,*) 'redo'
	call read_input_constant_file(time_ext,dt_smlk)
	ZRD93=time_ext+dt_smlk
	goto 200
	endif


300	continue
!	write(*,*) 'end'
	betp3r_0=betp3r(roc)
	LI3R_0=LI3R(ROC)
	if (plasma_up.eq.0) then
	IPL_0=ipl_bf_bkdw
	else
	IPL_0=IPL
	endif
	upl_0=upl(na1)
	QRADR_0=QRADR(ROC)
	QTOKR_0=QTOKR(ROC)
	CRAD4_0=CRAD4
	CMHD2_0=CMHD2
	CSCL4_0=CSCL4
	CDWM1_0=CDWM1
	CDWM2_0=CDWM2
	CDJM6_0=CDJM6
	CDJM1_0=CDJM1
	CDJM2_0=CDJM2
	CDJM4_0=CDJM4
	CDJM7_0=CDJM7
	CHE3_0=CHE3
	CDJM3_0=CDJM3
	ZRD77_0=ZRD77
	CDMJ5_0=CDMJ5
	CDMJ6_0=CDMJ6
	CDMJ7_0=CDMJ7
	CDHJ7_0=CDHJ7
	te_0(1:na1)=te(1:na1)
	ne_0(1:na1)=ne(1:na1)
	ccoil_0=ccoil(1:10)
	geom1d_0=geom1d(51:82)

	endif





	end subroutine wrsudg











      subroutine add_errors(qalp1,qrad1,ne1,qalp2,qrad2,ne2)  !for demo

	use parameter_inc
	use const_inc
	use status_inc
	use outcmn_inc
      implicit none


      integer i
      double precision qalp1,qrad1,ne1(na1)
      double precision qalp2,qrad2,ne2(na1)
      double precision dum1,dum2,dum3,dum4,dum5
      double precision dum11,dum21,dum31,dum41,dum51
      double precision dum12,dum22,dum32(na1),dum42(na1),dum52
      save dum4,dum5

      qalp2=qalp1
      qrad2=qrad1
      ne2=ne1

!add error on alpha power!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      dum1 = 0.01 !error on alpha power reconstruction
      dum2 = 0.4 ! integration time in s
      if (TIME-dum5.ge.dum2) then
          dum5=TIME
          dum4=qalp1
      endif
          qalp2=dum4
      call random_number(dum3)
      qalp2=qalp2*(1.+dum1*(1.-2.*dum3))
!add error on radiated power!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      dum11 = 0.1 !error on alpha power reconstruction
      dum21 = 0.5 ! integration time in s
      if (TIME-dum51.ge.dum21) then
          dum51=TIME
          dum41=qrad1
      endif
          qrad2=dum41
      call random_number(dum31)
      qrad2=qrad2*(1.+dum11*(1.-2.*dum31))
!add error on density power!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      dum12 = 0.02 !error on alpha power reconstruction
      dum22 = 0.001 ! integration time in s
      if (TIME-dum52.ge.dum22) then
          dum52=TIME
          dum42=ne(1:na1)
      endif
          ne2=dum42
      do i=1,na1
          call random_number(dum32(i))
      enddo
      ne2=ne2*(1.+dum12*(1.-2.*dum32))

      end
