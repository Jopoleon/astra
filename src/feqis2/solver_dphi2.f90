	subroutine	solver_dphi2(PSIb,Nr,Nt, &
     &      known_term,&
     &      dArc_rp1,&
     &      dArc_rm1,&
     &      dArc_rpt1,&
     &      dArc_rmt1,&
     &      dArc_tp1,&
     &      dArc_tm1,&
     &      dArc_tpr1,&
     &      dArc_tmr1,&
     &      ddr,&
     &      ddr_i,&
     &      dtp,&
     &      dtm,dt_i,&
     &      PSI)

	implicit none

	integer Nt,Nr,jt,jr,mm,Ndims,LDAB
	double precision known_term(Nr,Nt)
	double precision dArc_rp1(Nr,Nt)
	double precision dArc_rm1(Nr,Nt)
	double precision dArc_rpt1(Nr,Nt)
	double precision dArc_rmt1(Nr,Nt)
	double precision dArc_tp1(Nr,Nt)
	double precision dArc_tm1(Nr,Nt)
	double precision dArc_tpr1(Nr,Nt)
	double precision dArc_tmr1(Nr,Nt)
	double precision PSI(Nr,Nt)
	double precision ddr(Nr,Nt)
	double precision ddr_i(Nr,Nt)
	double precision dtp(Nr,Nt)
	double precision dtm(Nr,Nt),jt_2
	double precision dt_i(Nr,Nt),tt0
	double precision PSIb


	Ndims=1+(Nr-2)*Nt
	LDAB=2*(2*Nt)+(2*Nt)+1

	call	solver_inversion_matrix_gsef(PSIb,Nr,Nt,&
     &      known_term,Ndims,LDAB,&
     &      dArc_rp1,&
     &      dArc_rm1,&
     &      dArc_rpt1,&
     &      dArc_rmt1,&
     &      dArc_tp1,&
     &      dArc_tm1,&
     &      dArc_tpr1,&
     &      dArc_tmr1,&
     &      ddr,&
     &      ddr_i,&
     &      dtp,&
     &      dtm,dt_i,&
     &      PSI)






	end 
