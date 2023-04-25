subroutine solver_inversion_matrix_gsef(PSIb, Nr, Nt, &
   &  known_term, Ndims, LDAB, &
   &  dArc_rp1, dArc_rm1, dArc_rpt1, dArc_rmt1, &
   &  dArc_tp1, dArc_tm1, dArc_tpr1, dArc_tmr1, &
   &  ddr, ddr_i, dtp, dtm, dt_i, PSI)

implicit none

integer Nt,Nr,jt,jr,Ndims,NRHS,Nrp
	integer KL,KU,LDAB,LDB,IPIV(Ndims)
	integer INFO,i,iilow,iiup
	integer jrm,jrp,jdim00,jdim0p,jdimp0,jdim0m
	integer jdimm0,jdimmm,jdimmp,jdimpm,jdimpp
	integer inctype,jtm,jtp,jrawm
	double precision rm_arc,tp_arc,tm_arc,rp_arc
	double precision known_term(Nr,Nt)
	double precision dArc_rp1(Nr,Nt),denom
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
	double precision MATRIX(LDAB,Ndims)
	double precision BB(Ndims)
	double precision PSIb,PSI_0
	integer jrr,jjt

	Nrp=Nr-1
	NRHS=1
	KL=2*Nt
	KU=2*Nt
	
	MATRIX = 0.0
	BB = 0.0
	PSI = 0.0


	do jr = 1,Nrp

	if (jr.eq.1) then


	do jt = 1,Nt
						jtm=jt-1
						jtp=jt+1
	if (jt.eq.1) then
						jtm=Nt
	endif
	if (jt.eq.Nt) then
						jtp=1
	endif						
						jrp=2
						jrm=1

		jdim00=1

		jdimp0=1+jt+(jrp-2)*Nt			
		jdimpp=1+jtp+(jrp-2)*Nt			
		jdimpm=1+jtm+(jrp-2)*Nt			


		jrawm=KL+KU+1+jdim00
	

	BB(1)=BB(1)+known_term(1,jt)


	iilow=max(1,1-KU)
	iiup=min(Ndims,1+KL)
	if (jdim00.ge.iilow .and. jdim00.le.iiup) then
	MATRIX(jrawm-1,1)=MATRIX(jrawm-1,1)- &
     &   dArc_rp1(1,jt)/ddr(1,jt)
	endif

	iilow=max(1,jdimp0-KU)
	iiup=min(Ndims,jdimp0+KL)
	if (jdim00.ge.iilow .and. jdim00.le.iiup) then
	MATRIX(jrawm-jdimp0,jdimp0)=MATRIX(jrawm-jdimp0,jdimp0)+ &
     &   dArc_rp1(jr,jt)/ddr(jr,jt)
	endif

	iilow=max(1,jdimpp-KU)
	iiup=min(Ndims,jdimpp+KL)
	if (jdim00.ge.iilow .and. jdim00.le.iiup) then
	MATRIX(jrawm-jdimpp,jdimpp)=MATRIX(jrawm-jdimpp,jdimpp)+ &
     &   1./4.*dArc_rpt1(jr,jt)/dt_i(jr,jt)
	endif

	iilow=max(1,jdimpm-KU)
	iiup=min(Ndims,jdimpm+KL)
	if (jdim00.ge.iilow .and. jdim00.le.iiup) then
	MATRIX(jrawm-jdimpm,jdimpm)=MATRIX(jrawm-jdimpm,jdimpm)- &
     &   1./4.*dArc_rpt1(jr,jt)/dt_i(jr,jt)
	endif

	enddo


	endif


	if (jr.eq.2) then

	do jt = 1,Nt
									
						jtm=jt-1
						jtp=jt+1
	if (jt.eq.1) then
						jtm=Nt
	endif
	if (jt.eq.Nt) then
						jtp=1
	endif						
						jrp=jr+1
						jrm=1
							
							
		jdim00=1+jt+(jr-2)*Nt			
		jdimp0=1+jt+(jrp-2)*Nt			
		jdimm0=1+jt+(jrm-2)*Nt			
		jdimmm=1+jtm+(jrm-2)*Nt			
		jdimmp=1+jtp+(jrm-2)*Nt			
		jdim0p=1+jtp+(jr-2)*Nt			
		jdim0m=1+jtm+(jr-2)*Nt			
		jdimpp=1+jtp+(jrp-2)*Nt			
		jdimpm=1+jtm+(jrp-2)*Nt			


		jrawm=KL+KU+1+jdim00

	BB(jdim00)=known_term(jr,jt)


	rm_arc=-dArc_tmr1(jr,jt)/ddr_i(jr,jtm)
	rp_arc=dArc_tpr1(jr,jt)/ddr_i(jr,jt)
	tm_arc=-dArc_rmt1(jr,jt)/dt_i(jrm,jt)
	tp_arc=dArc_rpt1(jr,jt)/dt_i(jr,jt)
	denom=-(dArc_rp1(jr,jt)/ddr(jr,jt)+ &
     & dArc_rm1(jr,jt)/ddr(jrm,jt)+ &
     & dArc_tp1(jr,jt)/dtp(jr,jt)+ &
     & dArc_tm1(jr,jt)/dtm(jr,jt))




	iilow=max(1,jdim00-KU)
	iiup=min(Ndims,jdim00+KL)
	if (jdim00.ge.iilow .and. jdim00.le.iiup) then
	MATRIX(jrawm-jdim00,jdim00)=denom
	endif

	iilow=max(1,jdimp0-KU)
	iiup=min(Ndims,jdimp0+KL)
	if (jdim00.ge.iilow .and. jdim00.le.iiup) then
	MATRIX(jrawm-jdimp0,jdimp0)=dArc_rp1(jr,jt)/ddr(jr,jt) &
     &   +1./4.*(rp_arc+rm_arc)
	endif

	iilow=max(1,jdim0p-KU)
	iiup=min(Ndims,jdim0p+KL)
	if (jdim00.ge.iilow .and. jdim00.le.iiup) then
	MATRIX(jrawm-jdim0p,jdim0p)=dArc_tp1(jr,jt)/dtp(jr,jt) &
     &   +1./4.*(tp_arc+tm_arc)
	endif

	iilow=max(1,jdim0m-KU)
	iiup=min(Ndims,jdim0m+KL)
	if (jdim00.ge.iilow .and. jdim00.le.iiup) then
	MATRIX(jrawm-jdim0m,jdim0m)=dArc_tm1(jr,jt)/dtm(jr,jt) &
     &   +1./4.*(-tp_arc-tm_arc)
	endif

	iilow=max(1,jdimpp-KU)
	iiup=min(Ndims,jdimpp+KL)
	if (jdim00.ge.iilow .and. jdim00.le.iiup) then
	MATRIX(jrawm-jdimpp,jdimpp)=1./4.*(tp_arc+rp_arc)
	endif

	iilow=max(1,jdimpm-KU)
	iiup=min(Ndims,jdimpm+KL)
	if (jdim00.ge.iilow .and. jdim00.le.iiup) then
	MATRIX(jrawm-jdimpm,jdimpm)=1./4.*(-tp_arc+rm_arc)
	endif

	iilow=max(1,1-KU)
	iiup=min(Ndims,1+KL)
	if (jdim00.ge.iilow .and. jdim00.le.iiup) then
	MATRIX(jrawm-1,1)=dArc_rm1(jr,jt)/ddr(jrm,jt) &
     &   +1./2.*(-rp_arc-rm_arc) 
	endif


	enddo

	endif

	if (jr.gt.2 .and. jr.lt.Nrp) then

	do jt = 1,Nt
									
						jtm=jt-1
						jtp=jt+1
	if (jt.eq.1) then
						jtm=Nt
	endif
	if (jt.eq.Nt) then
						jtp=1
	endif						
						jrp=jr+1
						jrm=jr-1
							
							
		jdim00=1+jt+(jr-2)*Nt			
		jdimp0=1+jt+(jrp-2)*Nt			
		jdimm0=1+jt+(jrm-2)*Nt			
		jdimmm=1+jtm+(jrm-2)*Nt			
		jdimmp=1+jtp+(jrm-2)*Nt			
		jdim0p=1+jtp+(jr-2)*Nt			
		jdim0m=1+jtm+(jr-2)*Nt			
		jdimpp=1+jtp+(jrp-2)*Nt			
		jdimpm=1+jtm+(jrp-2)*Nt			


		jrawm=KL+KU+1+jdim00
					
					
	BB(jdim00)=known_term(jr,jt)

	rm_arc=-dArc_tmr1(jr,jt)/ddr_i(jr,jtm)
	rp_arc=dArc_tpr1(jr,jt)/ddr_i(jr,jt)
	tm_arc=-dArc_rmt1(jr,jt)/dt_i(jrm,jt)
	tp_arc=dArc_rpt1(jr,jt)/dt_i(jr,jt)
	denom=-(dArc_rp1(jr,jt)/ddr(jr,jt)+ &
     & dArc_rm1(jr,jt)/ddr(jrm,jt)+ &
     & dArc_tp1(jr,jt)/dtp(jr,jt)+ &
     & dArc_tm1(jr,jt)/dtm(jr,jt))



	iilow=max(1,jdim00-KU)
	iiup=min(Ndims,jdim00+KL)
	if (jdim00.ge.iilow .and. jdim00.le.iiup) then
	MATRIX(jrawm-jdim00,jdim00)=denom
	endif

	iilow=max(1,jdim0p-KU)
	iiup=min(Ndims,jdim0p+KL)
	if (jdim00.ge.iilow .and. jdim00.le.iiup) then
	MATRIX(jrawm-jdim0p,jdim0p)=dArc_tp1(jr,jt)/dtp(jr,jt) &
     &   +1./4.*(tp_arc+tm_arc)
	endif

	iilow=max(1,jdim0m-KU)
	iiup=min(Ndims,jdim0m+KL)
	if (jdim00.ge.iilow .and. jdim00.le.iiup) then
	MATRIX(jrawm-jdim0m,jdim0m)=dArc_tm1(jr,jt)/dtm(jr,jt) &
     &   +1./4.*(-tp_arc-tm_arc)
	endif

	iilow=max(1,jdimp0-KU)
	iiup=min(Ndims,jdimp0+KL)
	if (jdim00.ge.iilow .and. jdim00.le.iiup) then
	MATRIX(jrawm-jdimp0,jdimp0)=dArc_rp1(jr,jt)/ddr(jr,jt) &
     &   +1./4.*(rp_arc+rm_arc)
	endif

	iilow=max(1,jdimpp-KU)
	iiup=min(Ndims,jdimpp+KL)
	if (jdim00.ge.iilow .and. jdim00.le.iiup) then
	MATRIX(jrawm-jdimpp,jdimpp)=+1./4.*(rp_arc+tp_arc)
	endif

	iilow=max(1,jdimpm-KU)
	iiup=min(Ndims,jdimpm+KL)
	if (jdim00.ge.iilow .and. jdim00.le.iiup) then
	MATRIX(jrawm-jdimpm,jdimpm)=+1./4.*(-tp_arc+rm_arc)
	endif

	iilow=max(1,jdimm0-KU)
	iiup=min(Ndims,jdimm0+KL)
	if (jdim00.ge.iilow .and. jdim00.le.iiup) then
	MATRIX(jrawm-jdimm0,jdimm0)=dArc_rm1(jr,jt)/ddr(jrm,jt) &
     &   +1./4.*(-rp_arc-rm_arc)
	endif

	iilow=max(1,jdimmp-KU)
	iiup=min(Ndims,jdimmp+KL)
	if (jdim00.ge.iilow .and. jdim00.le.iiup) then
	MATRIX(jrawm-jdimmp,jdimmp)=+1./4.*(tm_arc-rp_arc)
	endif

	iilow=max(1,jdimmm-KU)
	iiup=min(Ndims,jdimmm+KL)
	if (jdim00.ge.iilow .and. jdim00.le.iiup) then
	MATRIX(jrawm-jdimmm,jdimmm)=+1./4.*(-rm_arc-tm_arc)
	endif

	enddo





	
	endif

	if (jr.eq.Nrp) then

	do jt = 1,Nt
									
						jtm=jt-1
						jtp=jt+1
	if (jt.eq.1) then
						jtm=Nt
	endif
	if (jt.eq.Nt) then
						jtp=1
	endif						
						jrp=jr+1
						jrm=jr-1
							
							
		jdim00=1+jt+(jr-2)*Nt			
		jdimp0=1+jt+(jrp-2)*Nt			
		jdimm0=1+jt+(jrm-2)*Nt			
		jdimmm=1+jtm+(jrm-2)*Nt			
		jdimmp=1+jtp+(jrm-2)*Nt			
		jdim0p=1+jtp+(jr-2)*Nt			
		jdim0m=1+jtm+(jr-2)*Nt			
		jdimpp=1+jtp+(jrp-2)*Nt			
		jdimpm=1+jtm+(jrp-2)*Nt			


		jrawm=KL+KU+1+jdim00

	rm_arc=-dArc_tmr1(jr,jt)/ddr_i(jr,jtm)
	rp_arc=dArc_tpr1(jr,jt)/ddr_i(jr,jt)
	tm_arc=-dArc_rmt1(jr,jt)/dt_i(jrm,jt)
	tp_arc=dArc_rpt1(jr,jt)/dt_i(jr,jt)
	denom=-(dArc_rp1(jr,jt)/ddr(jr,jt)+ &
     & dArc_rm1(jr,jt)/ddr(jrm,jt)+ &
     & dArc_tp1(jr,jt)/dtp(jr,jt)+ &
     & dArc_tm1(jr,jt)/dtm(jr,jt))



		BB(jdim00)=known_term(jr,jt)- &
     &   dArc_rp1(jr,jt)/ddr(jr,jt)*PSIb- &
     &   1./2.*(rp_arc+rm_arc)*PSIb



	iilow=max(1,jdim00-KU)
	iiup=min(Ndims,jdim00+KL)
	if (jdim00.ge.iilow .and. jdim00.le.iiup) then
	MATRIX(jrawm-jdim00,jdim00)=denom
	endif

	iilow=max(1,jdim0p-KU)
	iiup=min(Ndims,jdim0p+KL)
	if (jdim00.ge.iilow .and. jdim00.le.iiup) then
	MATRIX(jrawm-jdim0p,jdim0p)=dArc_tp1(jr,jt)/dtp(jr,jt) &
     &   +1./4.*(tp_arc+tm_arc)
	endif

	iilow=max(1,jdim0m-KU)
	iiup=min(Ndims,jdim0m+KL)
	if (jdim00.ge.iilow .and. jdim00.le.iiup) then
	MATRIX(jrawm-jdim0m,jdim0m)=dArc_tm1(jr,jt)/dtm(jr,jt) &
     &   +1./4.*(-tp_arc-tm_arc)
	endif

	iilow=max(1,jdimm0-KU)
	iiup=min(Ndims,jdimm0+KL)
	if (jdim00.ge.iilow .and. jdim00.le.iiup) then
	MATRIX(jrawm-jdimm0,jdimm0)=dArc_rm1(jr,jt)/ddr(jrm,jt) &
     &   +1./4.*(-rp_arc-rm_arc)
	endif

	iilow=max(1,jdimmp-KU)
	iiup=min(Ndims,jdimmp+KL)
	if (jdim00.ge.iilow .and. jdim00.le.iiup) then
	MATRIX(jrawm-jdimmp,jdimmp)=+1./4.*(tm_arc-rp_arc)
	endif

	iilow=max(1,jdimmm-KU)
	iiup=min(Ndims,jdimmm+KL)
	if (jdim00.ge.iilow .and. jdim00.le.iiup) then
	MATRIX(jrawm-jdimmm,jdimmm)=+1./4.*(-rm_arc-tm_arc)
	endif


	enddo

	endif

	enddo



call DGBSV(Ndims, KL, KU, NRHS, matrix, LDAB, IPIV, BB, Ndims, INFO)



	do jt = 1,Nt
	jr=1
	jdim00=1			
	PSI(jr,jt)=BB(jdim00)
		do jr = 2,Nrp
		jdim00=1+jt+(jr-2)*Nt			
	
			PSI(jr,jt)=BB(jdim00)
		enddo
	enddo
			PSI(Nr,:)=PSIb
	

return
end subroutine solver_inversion_matrix_gsef
