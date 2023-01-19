	subroutine CTAUG_STD(NCNB,NCTP,CTRLM,
     &       yvcoilx,yccoilx,t_starts,CCOIL,
     &       TIME,TAU,DUMCT,DUMCTP,VCOIL)

!	subroutine CTAUG_STD(
!
! NCNB       number of coils
! NCTP       number of quantities to be controlled
! CTRLM    (3)  (Ip,Pos,Shape) control mode
! yvcoilx  (ncoils)  feed forward voltages
! yccoilx  (ncoils) feed forward coil currents
! t_starts  time when control should start
! TIME     present time
! TAU     dt
! DUMCT  (nquant) target quantities values
! DUMCTP   (nquant) actual quantities values
! VCOIL   (ncoils)  output voltages
!

! values can be: 0 -> feedforward
!																> 0 -> feedback control mode
!																-1 -> freeze actual coil current values

!use		call GETCOILS(VCOIL,CCOIL) to get ff coil currents to be used
!
!DUMCT contains the target values of controlled quantities
!DUMCTP should contain the actual values of controlled quantities
! Note that currents per winding are the real cable currents!
!Quantities:
! 1 Ip    
! 2 Raus
! 3 Rcurr
! 4 Rgeo
! 5 Zcurr
! 6 Zsquad
! 7 Zoben
! 8 Zgera2b
! 9 Zgeri2b
! 10 delRoben
! 11 Rin
! (12) -9 CoIo I
! (13) -10 CoIu I
! controller subroutine
	
	implicit none
C	include	'for/parameter.inc'
C	include 'for/const.inc'
C	include 'for/outcmn.inc'
C	include 'for/status.inc'

	integer t_ip,t_pos,t_shape
	integer NCNB,NCTP
	double precision z_ip,z_pos,z_shape
	double precision CTRLM(3)
	double precision VCOIL(NCNB)
	double precision CCOIL(NCNB)
	double precision DUMCT(NCTP)
	double precision DUMCTP(NCTP)

! AUG coil description
	integer i,j,k,i_quantities,ncoilv
	integer npassive
	parameter(npassive=2)
	integer use_coil_I(NCNB)
	integer use_coil_P(NCNB)
	integer use_coil_S(NCNB)
	integer P_scalep(NCNB)
	double precision Q_scalep(NCNB),dumy0,dumy1(1),dumy2(1,1)
	double precision t_starts
	integer use_q_I(NCTP),ismooths
	integer use_q_P(NCTP)
	integer use_q_S(NCTP),i_usedc(NCNB)
	integer i_nousedc(NCNB)
	integer ncoils,nCqq_c,nAqq_c,i_control_i
	integer ioffseti,icommandx,dumi0,dumi1(1),dumi2(1,1)
	integer ni_sets,np_sets,ns_sets,nq_shape,nc00

	double precision ic_sets(10),ip_sets(10),is_sets(10),TIME,TAU
	double precision Vlimits(NCNB),ydumpp(NCTP),ydump(NCTP)
		
	character*12  :: name_gsefdir='exp/equ/aug/' ! working directory path
	character*29  :: name_cgsefdir='exp/equ/aug/controller_files/' ! working directory path
	character*13  :: Ictr_fb_mode_file='I_ctrl_fbmode' ! time evolution of control modes file
	character*13 ::  Pctr_fb_mode_file='P_ctrl_fbmode' ! time evolution of control modes file
	character*13 :: Sctr_fb_mode_file='S_ctrl_fbmode' ! time evolution of control modes file
	character*13 :: PV_control_file='pv_params.dat' ! time evolution of control modes file
	character*80  Ifile,Pfile,Sfile
	character*11 :: currents_now_file = 'currents.wr' ! controller file
	double precision  yvcoilx(NCNB),yccoilx(NCNB)
	double precision  yccoili(NCNB)
	double precision  yccoiln(NCNB)
	double precision  yccoilnt(NCNB-npassive)
	double precision  yccoilit(NCNB-npassive)
	double precision  yccoilitt(NCNB-npassive)
	double precision  yccoilot(NCNB-npassive)
	double precision  yvcoilt(NCNB-npassive)
	double precision  yvcoilo(NCNB-npassive)
	double precision I_K_P(NCNB-npassive,NCTP)
	double precision I_K_I(NCNB-npassive,NCTP)
	double precision P_K_P(NCNB-npassive,NCTP)
	double precision P_K_I(NCNB-npassive,NCTP)
	double precision S_K_P(NCNB-npassive,NCTP)
	double precision S_K_I(NCNB-npassive,NCTP)
	double precision combci(NCNB-npassive,NCNB-npassive)
	double precision V_K_P(NCNB-npassive,NCNB-npassive)
	double precision V_K_I(NCNB-npassive,NCNB-npassive)
	double precision combco(NCNB-npassive,NCNB-npassive)
	double precision combvi(NCNB-npassive,NCNB-npassive)
	double precision combvo(NCNB-npassive,NCNB-npassive)
	double precision scalfac
	character*1 dumcc
	integer imodenowc,imodenowp,imodenows,istepctrl
	data imodenowc /0/
	data imodenowp /0/
	data imodenows /0/
	data istepctrl /0/
	save imodenowc,imodenowp,imodenows,istepctrl
	
	z_ip=CTRLM(1)
	z_pos=CTRLM(2)
	z_shape=CTRLM(3)

	ncoils=NCNB
	ncoilv=NCNB-npassive
	nq_shape=0
	i_quantities=NCTP
	
	V_K_I=0.
	V_K_P=0.
	combvi=0.
	combvo=0.
	P_scalep=0
	Q_scalep=0.

	nc00=1
	dumi0=0
	dumi1=0
	dumi2=0
	dumy0=0.
	dumy1=0.
	dumy2=0.
	
	use_coil_I=0.
	use_coil_P=0.
	use_coil_S=0.
	use_q_I=0.
	use_q_P=0.
	use_q_S=0.
	I_K_P=0.
	I_K_I=0.
	P_K_P=0.
	P_K_I=0.
	S_K_P=0.
	S_K_I=0.

!Set coil values to experimental ones for ff purposes
C	call GETCOILS(yvcoilx,yccoilx)

! yccoilx is the experimental coil currents
! yccoili is the actual controlled coil currents CCOIL
! yccoiln is the actual simulation currents values from SPIDER
! CCOIL is the controlled currents after the controller

!Note that controllore values assume A for the coil currents, while the voltage controller accepts kA

	yccoili(1:NCNB)=CCOIL(1:NCNB)
	yccoilit=0.
	yccoilot=0.

!Do only if in free boundary mode



!read spider currents (simulated ones)
	open(32,file=name_gsefdir//currents_now_file)
	read(32,*) i,i
	read(32,*) yccoiln(1:ncoils)
	close(32)

C	IPCTRL = 1.  !Set IPCTRL variable to 1 


	ic_sets(1)=1. 
	ic_sets(2)=2.
	ic_sets(3)=3. 

	ip_sets(1)=1. 
	ip_sets(2)=2. 
	ip_sets(3)=3. 
	ip_sets(4)=4. 
	ip_sets(5)=5. 
	ip_sets(6)=6. 

	is_sets(1)=1.
	is_sets(2)=2.
	is_sets(3)=51.
	is_sets(4)=101.
	is_sets(5)=121.
	is_sets(6)=511.
	is_sets(7)=601.

	ni_sets=3
	np_sets=6
	ns_sets=7

	t_ip=nint(z_ip)
	t_pos=nint(z_pos)
	t_shape=nint(z_shape)


	if (istepctrl.eq.0) then
	imodenowc=t_ip
	imodenowp=t_pos
	imodenows=t_shape
		istepctrl=1
	endif




!Plasma current control, at the moment limits are ignored
	i_control_i=0
	if (t_ip.gt.0) then
		write(Ifile,'(A,A,i0,A,F7.4,A)') Ictr_fb_mode_file,'_',t_ip,
     @    '.dat'
	open(32,file=name_cgsefdir//Ifile)
	read(32,*) use_coil_I(1:ncoilv)
	read(32,*) use_q_I(1)
	read(32,*) I_K_P(1,1),I_K_I(1,1)
	close(32)
	nCqq_c=1
	nAqq_c=1
!Call current control
	scalfac=1.E+03
	ioffseti=1
	P_scalep(1)=0
	Q_scalep(1)=0.
	icommandx=1
	ismooths=1
	
!Check if mode has changed
	if (t_ip.lt.imodenowc .or. t_ip.gt.imodenowc) then
	if (ismooths.eq.1) then
!Initialize controller for smooth beginning
	call CONTROLLER_AUG(dumy0,dumy0,nc00,nc00,dumy2,dumy2,
     & dumy1,dumy1,dumy1,dumy1,i_control_i,dumi0,
     & dumi1,dumy1,dumy0,0,dumi0)
	imodenowc	= t_ip
	endif
	endif
	
	
	CCOIL(1)=yccoiln(1)*1.E+03 !initial condition
	call CONTROLLER_AUG(TIME,TAU,nCqq_c,nAqq_c,I_K_P,I_K_I,
     & DUMCT(1),DUMCTP(1),CCOIL(1),yccoilot(1),i_control_i,ioffseti,
     & P_scalep(1),Q_scalep(1),scalfac,icommandx,ismooths)
C	write(*,*) 'call current'
	CCOIL(1)=yccoilot(1)
	endif
	if (t_ip.eq.0) then
 !set OH current to ff value
	icommandx=0
				CCOIL(1)=yccoilx(1)  
	call CONTROLLER_AUG(dumy0,dumy0,nc00,nc00,dumy2,dumy2,
     & dumy1,dumy1,dumy1,dumy1,i_control_i,dumi0,
     & dumi1,dumy1,dumy0,icommandx,dumi0)
C	write(*,*) 'reinit current'
	endif
	if (t_ip.lt.0) then
		!freeze CCOIL value, do nothing
	CCOIL(1)=yccoili(1)
	icommandx=0
	call CONTROLLER_AUG(dumy0,dumy0,nc00,nc00,dumy2,dumy2,
     & dumy1,dumy1,dumy1,dumy1,i_control_i,dumi0,
     & dumi1,dumy1,dumy0,icommandx,dumi0)
C	write(*,*) 'reinit current'
	endif
	












!Plasma position control
	i_control_i=1
	if (t_pos.gt.0) then
		write(Pfile,'(A,A,i0,A,F7.4,A)') Pctr_fb_mode_file,'_',t_pos,
     @    '.dat'
	open(32,file=name_cgsefdir//Pfile)
	read(32,*) (use_coil_P(i),i=1,ncoilv)
	read(32,*) (use_q_P(i),i=1,sum(use_coil_P))
	read(32,*) (P_scalep(i),i=1,sum(use_coil_P))

	Q_scalep(1)=DUMCT(1)
	Q_scalep(2)=DUMCT(1)

			read(32,*) ((P_K_P(j,i),i=1,sum(use_coil_P)),
     & j=1,sum(use_coil_P))
			read(32,*) ((P_K_I(j,i),i=1,sum(use_coil_P)),
     & j=1,sum(use_coil_P))
	close(32)
	nCqq_c=2
	nAqq_c=2

!Call position control
	scalfac=1.E+03 ! from kA to A
	ioffseti=1
	icommandx=1
	ismooths=1

	CCOIL(9:10)=yccoiln(9:10)*1.E+03 !initial condition

	yccoilit(1:2)=CCOIL(9:10)
	ydump(1:2)=DUMCT(use_q_P(1:sum(use_coil_P)))
	ydumpp(1:2)=DUMCTP(use_q_P(1:sum(use_coil_P)))

!Check if mode has changed
	if (t_pos.lt.imodenowp .or. t_pos.gt.imodenowp) then
	if (ismooths.eq.1) then
!Initialize controller for smooth beginning
	call CONTROLLER_AUG(dumy0,dumy0,nc00,nc00,dumy2,dumy2,
     & dumy1,dumy1,dumy1,dumy1,i_control_i,dumi0,
     & dumi1,dumy1,dumy0,0,dumi0)
	imodenowp	= t_pos
	endif
	endif


	call CONTROLLER_AUG(TIME,TAU,nCqq_c,nAqq_c,
     & P_K_P(1:2,1:2),
     & P_K_I(1:2,1:2),
     & ydump(1:2),
     & ydumpp(1:2),
     & yccoilit(1:2),yccoilot(1:2),
     & i_control_i,
     & ioffseti,P_scalep(1:2),
     & Q_scalep(1:2),scalfac,icommandx,ismooths)

	CCOIL(9:10)=yccoilot(1:2)
	endif
	if (t_pos.eq.0) then
				CCOIL(9)=yccoilx(9)   !set CoIo current to ff value
				CCOIL(10)=yccoilx(10)   !set CoIu current to ff value
	icommandx=0
	call CONTROLLER_AUG(dumy0,dumy0,nc00,nc00,dumy2,dumy2,
     & dumy1,dumy1,dumy1,dumy1,i_control_i,dumi0,
     & dumi1,dumy1,dumy0,icommandx,dumi0)
C	write(*,*) 'reinit pos'
	endif
	if (t_pos.lt.0) then
				CCOIL(9)=yccoili(9)   !set CoIo current to freeze value
				CCOIL(10)=yccoili(10)   !set CoIu current to freeze value
		!freeze CCOIL value, do nothing
	icommandx=0
	call CONTROLLER_AUG(dumy0,dumy0,nc00,nc00,dumy2,dumy2,
     & dumy1,dumy1,dumy1,dumy1,i_control_i,dumi0,
     & dumi1,dumy1,dumy0,icommandx,dumi0)
C	write(*,*) 'reinit pos'
	endif
!Check discrepancies










!Plasma shape control
	i_control_i=2
	if (t_shape.gt.0) then
	write(Sfile,'(A,A,i0,A,F7.4,A)') Sctr_fb_mode_file,
     @   '_',t_shape,
     @    '.dat'
	open(32,file=name_cgsefdir//Sfile)
	read(32,*) use_coil_S(1:ncoilv)
	read(32,*) nq_shape
	read(32,*) use_q_S(1:nq_shape)
	read(32,*) (P_scalep(i),i=1,nq_shape)

	Q_scalep(1:nq_shape)=DUMCT(1)

	j=0
	do i=1,ncoilv
	if (use_coil_S(i).eq.1) then
		j=j+1
		i_usedc(j)=i
	endif
	enddo
	
	
	do i=1,nq_shape
	if (use_q_S(i).eq.-9) 	use_q_S(i)=12
	if (use_q_S(i).eq.-10) 	use_q_S(i)=13
	enddo


			read(32,*) ((combci(j,i),i=1,ncoilv),j=1,ncoilv)
			read(32,*) ((combco(j,i),i=1,ncoilv),j=1,ncoilv)

			read(32,*) ((S_K_P(j,i),i=1,nq_shape),
     @   j=1,sum(use_coil_S))
			read(32,*) ((S_K_I(j,i),i=1,nq_shape),
     @   j=1,sum(use_coil_S))

	close(32)
	nCqq_c=nq_shape
	nAqq_c=sum(use_coil_S)
!Call shape control
!Combine currents in input
		yccoilnt=0.
		yccoilit=0.
		yccoilitt=0.
	CCOIL(2:8)=yccoilx(2:8) !ff condition
	do i=1,ncoilv
		yccoilit(i)=sum(combci(i,:)*yccoiln(1:ncoilv)*1.E+03) !initial condition
	enddo

		i=1	
		yccoilitt(i)=CCOIL(1) !use CS current as it is

		i=2	
		yccoilitt(i)=yccoilx(2)-yccoilx(1) !ff condition

		i=3	
		yccoilitt(i)=yccoilx(3)-yccoilx(2) !ff condition

	do i=4,ncoilv
		yccoilitt(i)=sum(combci(i,1:ncoilv)*yccoilx(1:ncoilv)) !ff condition
	enddo


	scalfac=1.E+03
	ioffseti=1
	icommandx=1
	ismooths=1

!Check if mode has changed
	if (t_shape.lt.imodenows .or. t_shape.gt.imodenows) then
	if (ismooths.eq.1) then
!Initialize controller for smooth beginning
	call CONTROLLER_AUG(dumy0,dumy0,nc00,nc00,dumy2,dumy2,
     & dumy1,dumy1,dumy1,dumy1,i_control_i,dumi0,
     & dumi1,dumy1,dumy0,0,dumi0)
	imodenows	= t_shape
	endif
	endif

	call CONTROLLER_AUG(TIME,TAU,nCqq_c,nAqq_c,
     & S_K_P(1:sum(use_coil_S),1:nq_shape),
     & S_K_I(1:sum(use_coil_S),1:nq_shape),
     & DUMCT(use_q_S(1:nq_shape)),
     & DUMCTP(use_q_S(1:nq_shape)),
     & yccoilit(i_usedc(1:sum(use_coil_S))),
     & yccoilnt(1:sum(use_coil_S)),
     & i_control_i,ioffseti,
     & P_scalep(1:nq_shape),
     & Q_scalep(1:nq_shape),
     & scalfac,icommandx,ismooths)
C	write(*,*)  'fdfsfds ', yccoilit(i_usedc(1:sum(use_coil_S)))
C	write(*,*)  'fdfsfds ', yccoilnt(1:sum(use_coil_S))
C	write(*,*)  'fdfsfds ', i_control_i,ioffseti,
C     & i_usedc(1:sum(use_coil_S)),
C     & P_scalep(1:sum(use_coil_S)),
C     & Q_scalep(1:sum(use_coil_S)),
C     & scalfac,icommandx,ismooths
C	write(*,*) 'call shape'
!Combine currents in output
C	write(*,*) i_usedc(1:sum(use_coil_S))
	yccoilot(i_usedc(1:sum(use_coil_S)))=yccoilnt(1:sum(use_coil_S))
!Now look for coils that did no nothing
	j=0
	i_nousedc=0
	do i=1,8
	if (use_coil_S(i).eq.0) then
		j=j+1
		i_nousedc(j)=i
	endif
	enddo
!set them to ff value
	yccoilot(i_nousedc(1:j))=yccoilitt(i_nousedc(1:j))  
C	write(*,*) i_usedc
C	write(*,*) i_nousedc
C	write(*,*) yccoilot(1:8)
C	write(*,*) yccoilitt(1:8)
C	write(*,*) yccoilit(1:8)
C	write(*,*) yccoilx(6:7)
C	write(*,*) yccoiln(6:7)*1.E+03


	do i=1,sum(use_coil_S)
		CCOIL(i_usedc(i))=sum(combco(i_usedc(i),
     @  1:8)*
     @  yccoilot(1:8))
C	write(*,*) combco(i_usedc(i),
C     @  i_usedc(1:sum(use_coil_S)))
	enddo

	do i=1,j
		CCOIL(i_nousedc(i))=sum(combco(i_nousedc(i),
     @  1:8)*
     @  yccoilot(1:8))
	enddo
C	write(*,*) 'biog dcoils ',yccoilot(1:8)
C	write(*,*) 'biog coils ',CCOIL(1:8)


C	pause
	endif
	if (t_shape.eq.0) then
				CCOIL(2:8)=yccoilx(2:8)   !set shape control currents to ff value
				CCOIL(2)=yccoilx(2)-yccoilx(1)+CCOIL(1)   !for dIOH2o is equal to CS + ff
				CCOIL(3)=yccoilx(3)-yccoilx(2)+CCOIL(2)   !  dIOH2u is ok
	icommandx=0
	call CONTROLLER_AUG(dumy0,dumy0,nc00,nc00,dumy2,dumy2,
     & dumy1,dumy1,dumy1,dumy1,i_control_i,dumi0,
     & dumi1,dumy1,dumy0,icommandx,dumi0)
C	write(*,*) 'reinit shape'
		endif
	if (t_shape.lt.0) then
	CCOIL(2:8)=yccoili(2:8)  !freeze currents
		!freeze CCOIL value, do nothing
	icommandx=0
	call CONTROLLER_AUG(dumy0,dumy0,nc00,nc00,dumy2,dumy2,
     & dumy1,dumy1,dumy1,dumy1,i_control_i,dumi0,
     & dumi1,dumy1,dumy0,icommandx,dumi0)
C	write(*,*) 'reinit shape'
	endif









!Now do voltage control	
	i_control_i=3
	open(32,file=name_gsefdir//PV_control_file)
	read(32,*) (V_K_P(i,i),i=1,ncoilv)
	read(32,*) ((combvi(i,j),j=1,ncoilv),i=1,ncoilv)
	read(32,*) ((combvo(i,j),j=1,ncoilv),i=1,ncoilv)
	close(32)
	VCOIL(1:ncoils)=0.
	nAqq_c=ncoilv

C	if (yvcoilx(2).gt.1.E+04) then
C				CCOIL(2)=CCOIL(1)
C				yccoiln(2)=yccoiln(1)
C	endif
C	if (yvcoilx(3).gt.1.E+04) then
C				CCOIL(3)=CCOIL(2)
C				yccoiln(3)=yccoiln(2)
C	endif

!Combine currents in input
	do i=1,ncoilv
		yccoilit(i)=sum(combvi(i,:)*CCOIL(1:ncoilv))  ! reference controller currents
		yccoilnt(i)=sum(combvi(i,:)*yccoiln(1:ncoilv))*1.E+03  !from MA to KA, actual coil currents
	enddo

	scalfac=1.
	ioffseti=0
	P_scalep(1:ncoilv)=0
	Q_scalep(1:ncoilv)=0.
	icommandx=1
	ismooths=0
	call CONTROLLER_AUG(TIME,TAU,nAqq_c,nAqq_c,V_K_P,V_K_I,
     & yccoilit(1:ncoilv),yccoilnt(1:ncoilv),
     & yvcoilt(1:ncoilv),yccoilot(1:ncoilv),i_control_i,ioffseti,
     & P_scalep(1:ncoilv),
     & Q_scalep(1:ncoilv),scalfac,icommandx,ismooths)
	
	yvcoilt(1:ncoilv)=yccoilot(1:ncoilv)

	
!Combine voltages in output
	do i=1,ncoilv
		yvcoilo(i)=sum(combvo(i,:)*yvcoilt(1:ncoilv))
	enddo

!Voltage limits
	Vlimits(1)=5000.
	Vlimits(2:3)=5000.
	Vlimits(4)=1100.
	Vlimits(5)=1100.
	Vlimits(6)=2000.
	Vlimits(7)=2000.
	Vlimits(8)=2000.
	Vlimits(9:10)=500.

	do i=1,ncoilv
		if (yvcoilo(i).ge.Vlimits(i)) yvcoilo(i)=Vlimits(i)
		if (yvcoilo(i).le.-Vlimits(i)) yvcoilo(i)=-Vlimits(i)
	enddo



!Update voltages
	VCOIL(1:ncoilv)=yvcoilo(1:ncoilv)	



	if (TIME.le.t_starts) VCOIL(1:NCNB)=0.

	VCOIL(11:11+npassive-1)=0.  !At the moment set zero voltage in passive loops

C	stop




!Call controller
!		call CONTROLLER_AUG(TIME,TAU,ncoils)
		
	
	
	
	end
