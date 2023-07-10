      subroutine PELAUG
!C
!C	Elms routine
!C
!C-------------------------------
	use parameter_inc
	use const_inc
	use status_inc
	use outcmn_inc
!      use fenix_params

      implicit none
      integer		i,j,k,idone
      double precision pelpos,pelwid, &
     & chi_elm,scal_fac,scal_fac2,pelmass
	double precision t_elm,ts_elm,tau_elm,dw_elm,dt_elm
	double precision ngw,xx(8,2),const1,const2,qqq,avgcpel
	double precision xfan,dt_pelz,dtt_pelz,tpel1,tpel2
	integer pedpos,imodel		,solll	,iki			
	save xx,const1,const2,tpel1,tpel2,avgcpel
	data tpel1 /0./
	data tpel2 /0./
	data avgcpel /0./
!	save pelpos, pelwid				
!%ipvec(jscan,1)=ipscan;
!%tepedvec(jscan,1)=tepedtopscan;
!%tesepvec(jscan,1)=tesep;
!%fgwvec(jscan,1)=fgw;
!%te0vec(jscan,1)=te0;
!%nesepfgwvec(jscan,1)=nesepfgw/fgw;
!%nepedtopfgwvec(jscan,1)=nepedtopfgw/fgw;
!%nepikvec(jscan,1)=ne0./nepedtopfgw;
!%pmassvec(jscan,1)=pelletmassscan;


!	return

	iki=2
	
	if (iki.eq.1) then
		dt_pelz=0.8
		dtt_pelz=0.015	

		cpel1=0.
		cimp3=0.	
	
		if (TIME-tpel1.lt.dt_pelz) tpel2=0.
	
		if (TIME-tpel1.ge.dt_pelz) then
			if (tpel2.le.dtt_pelz) then	
				cpel1=30000.
				cimp3=30000.
				tpel2=tpel2+tau
			endif
			if (tpel2.gt.dtt_pelz) then	
				cpel1=0.
				cimp3=0.
				tpel1=TIME
			endif
		endif
		else
	endif

!	write(454,'(955E25.11)') TIME,CPEL1+CIMP3,NE(1:NA1),
!     & TE(1:NA1),HE(1:NA1),XI(1:NA1),CN(1:NA1)


	pedpos=nint(0.91*na1)
	ngw=10.*ipl/(GP*abc**2.)



	imodel=nint(ZRD20)
	solll=1
	

	
	if (imodel.eq.1) then

	cneut1=0.6
	cneut2=0.2
	else

	if (solll.eq.1.and.TIME.le.ZRD78+3*TAU) then ! solution 1 of P. Lang, not penetrating into blanket 
	open(32,file='dat/pelparams_AUG.dat')
	read(32,*) xx(1,1),xx(1,2),const1,const2
	read(32,*) xx(2,1),xx(2,2),const1,const2
	read(32,*) xx(3,1),xx(3,2),const1,const2
	read(32,*) xx(4,1),xx(4,2),const1,const2
	read(32,*) xx(5,1),xx(5,2),const1,const2
	read(32,*) xx(6,1),xx(6,2),const1,const2
	read(32,*) xx(7,1),xx(7,2),const1,const2
	read(32,*) xx(8,1),xx(8,2),const1,const2
	read(32,*) xx(9,1),xx(9,2),const1,const2
	read(32,*) xx(10,1),xx(10,2),const1,const2
	close(32)
	endif
	
	qqq=1./mu(pedpos)

	if (CPEL1.gt.1.e-6) then	
	pelmass=24.

	cneut1=(1.- &
     & (const1*IPL**xx(1,1) * te(pedpos)**xx(2,1) &
     &  * te(na1)**xx(3,1) * te(1)**xx(4,1) &
     &  * (ne(na1)/ngw)**xx(5,1) * &
     &  (ne(pedpos)/ngw)**xx(6,1) *  &
     & (ne(1)/ne(pedpos))**xx(7,1) *  &
     & (pelmass)**xx(8,1) *  &
     & (ZRD33/1.e3)**xx(9,1)*BTOR**xx(10,1)) &
     &  )

	cneut2=min(0.5,max(0.001,const2* &
     & IPL**xx(1,2) * te(pedpos)**xx(2,2) &
     &  * te(na1)**xx(3,2) * te(1)**xx(4,2) &
     &  * (ne(na1)/ngw)**xx(5,2) * &
     &  (ne(pedpos)/ngw)**xx(6,2) *  &
     & (ne(1)/ne(pedpos))**xx(7,2) *  &
     & (pelmass)**xx(8,2) *  &
     & (ZRD33/1.e3)**xx(9,2)*BTOR**xx(10,2))) 
!	write(*,*) 'pelparms a',cneut1,cneut2

	cneut1=xfan(cneut1)
	cneut2=xfan(cneut2)


	endif





	endif	

!	write(*,*) const1,IPL,xx(1,1),te(pedpos),xx(2,1)
!     &  , te(na1),xx(3,1),te(1),xx(4,1)
!     &  , (ne(na1)/ngw),xx(5,1),
!     &  (ne(pedpos)/ngw),xx(6,1),(ne(1)/ne(pedpos)),xx(7,1), 
!     & (cpel1*0.001),xx(8,1) ,
!     & (ZRD33/1.e3),xx(9,1),BTOR,xx(10,1)


!	write(*,*) 'pelparms rho',cneut1,cneut2,
!     & CPEL1
               
                         
	end




