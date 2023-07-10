!C======================================================================|
      subroutine PELscal
!C----------------------------------------------------------------------|
!C
!C	Elms routine
!C
!C----------------------------------------------------------------------|
	use parameter_inc
	use status_inc
	use const_inc
	use outcmn_inc
	
      implicit none
!      include	'for/parameter.inc'
!      include 'for/const.inc'
!      include 'for/outcmn.inc'
!      include 'for/status.inc'
      integer		i,j,k,idone
      double precision pelpos,pelwid,&
     & chi_elm,scal_fac,scal_fac2
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



	imodel=2
	solll=1
	if (TIME.gt.30..and.CPEL1+CIMP3.gt.0.1) then
!	write(*,*) 'pellet',TIME,NE(90)
	
	endif


	
	if (imodel.eq.1) then

	cneut1=0.9
	cneut2=0.1
	else

	if (solll.eq.1.and.TIME.le.TAU) then ! solution 1 of P. Lang, not penetrating into blanket 
	open(32,file='dat/pelparams.dat')
	read(32,*) xx(1,1),xx(1,2),const1,const2
	read(32,*) xx(2,1),xx(2,2),const1,const2
	read(32,*) xx(3,1),xx(3,2),const1,const2
	read(32,*) xx(4,1),xx(4,2),const1,const2
	read(32,*) xx(5,1),xx(5,2),const1,const2
	read(32,*) xx(6,1),xx(6,2),const1,const2
	read(32,*) xx(7,1),xx(7,2),const1,const2
	read(32,*) xx(8,1),xx(8,2),const1,const2
	close(32)
	endif
	if (solll.eq.2.and.TIME.le.TAU) then ! solution 2 of P. Lang, penetrating into blanket
	open(32,file='dat/pelparams2.dat')
	read(32,*) xx(1,1),xx(1,2),const1,const2
	read(32,*) xx(2,1),xx(2,2),const1,const2
	read(32,*) xx(3,1),xx(3,2),const1,const2
	read(32,*) xx(4,1),xx(4,2),const1,const2
	read(32,*) xx(5,1),xx(5,2),const1,const2
	read(32,*) xx(6,1),xx(6,2),const1,const2
	read(32,*) xx(7,1),xx(7,2),const1,const2
	read(32,*) xx(8,1),xx(8,2),const1,const2
	close(32)
	endif
	if (solll.eq.3.and.TIME.le.TAU) then ! solution 2 of P. Lang, penetrating into blanket
	open(32,file='dat/pelparams3.dat')
	read(32,*) xx(1,1),xx(1,2),const1,const2
	read(32,*) xx(2,1),xx(2,2),const1,const2
	read(32,*) xx(3,1),xx(3,2),const1,const2
	read(32,*) xx(4,1),xx(4,2),const1,const2
	read(32,*) xx(5,1),xx(5,2),const1,const2
	read(32,*) xx(6,1),xx(6,2),const1,const2
	read(32,*) xx(7,1),xx(7,2),const1,const2
	read(32,*) xx(8,1),xx(8,2),const1,const2
	close(32)
	endif
	if (solll.eq.4.and.TIME.le.TAU) then ! solution 2 of P. Lang, penetrating into blanket
	open(32,file='dat/pelparams4.dat')
	read(32,*) xx(1,1),xx(1,2),const1,const2
	read(32,*) xx(2,1),xx(2,2),const1,const2
	read(32,*) xx(3,1),xx(3,2),const1,const2
	read(32,*) xx(4,1),xx(4,2),const1,const2
	read(32,*) xx(5,1),xx(5,2),const1,const2
	read(32,*) xx(6,1),xx(6,2),const1,const2
	read(32,*) xx(7,1),xx(7,2),const1,const2
	read(32,*) xx(8,1),xx(8,2),const1,const2
	close(32)
	endif
	
	qqq=61.*mu(pedpos)

	if (CPEL1+CIMP3.gt.1.e-6) then	


	cneut1=(1.-const1*qqq**xx(1,1) * te(pedpos)**xx(2,1) &
     &  * te(na1)**xx(3,1) * te(1)**xx(4,1) &
     &  * (ne(na1)/ngw)**xx(5,1) * &
     &  (ne(pedpos)/ngw)**xx(6,1) * (ne(1)/ne(pedpos))**xx(7,1) * &
     & (1.*(cpel1+cimp3)*0.01)**xx(8,1))

	cneut2=max(0.001,const2*qqq**xx(1,2) * &
     &  te(pedpos)**xx(2,2) * &
     &  te(na1)**xx(3,2) * te(1)**xx(4,2)  &
     & * (ne(na1)/ngw)**xx(5,2) *  &
     & (ne(pedpos)/ngw)**xx(6,2) * &
     &  (ne(1)/ne(pedpos))**xx(7,2) * & 
     & (1.*(cpel1+cimp3)*0.01)**xx(8,2)) 
!	write(*,*) 'pelparms a',cneut1,cneut2

	cneut1=xfan(cneut1)
!	write(*,*) 'pelparms rho',cneut1,cneut2
	write(*,*) 'pelellet stuf ',pedpos,cneut1,cneut2,cpel1,cimp3
	endif
	endif	

	if (TIME.ge.150.and.time.lt.181) then
!	write(967,'(101E25.11)') TIME,NE(1:NA1)	
!	write(968,'(101E25.11)') TIME,TE(1:NA1)	
!	write(965,'(101E25.11)') TIME,DN(1:NA1)	
!	write(966,'(101E25.11)') TIME,car53(1:NA1)	
	
	endif

!	write(*,*) 'pelellet stuf ',qqq,te(pedpos),
!     &   te(na1),te(1),
!     &  (ne(na1)/ngw),
!     &  (ne(pedpos)/ngw),(ne(1)/ne(pedpos)),
!     & (2.*cpel1*0.01)
!	write(*,*) 'pelellet stuf ',xrho(pedpos),ametr(pedpos)/abc

	write(*,*) 'pelellet stuf ',pedpos,cneut1,cneut2,cpel1,cimp3
	write(*,*) 'ne',time,ne(1)
	if (NE(1).ge.30) stop


               
                         
	end



!C======================================================================|

