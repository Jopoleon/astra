!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
       SUBROUTINE  FEQISUPDATE(machine,coilzzz,time_nowz,nccc)

	use ef_circuit
	use exchange_with_astra


	implicit none


	integer nccc
		character*4 machine
	double precision coilzzz(nccc),time_nowz


	
	if (machine.eq.'aug '.or.machine.eq.'aug_' &
     & .or.machine.eq.'aug5') then 
	coilzzz(1)=curconduc(1)
	coilzzz(2)=curconduc(2)-curconduc(1)
	coilzzz(3)=curconduc(3)-curconduc(2)
	coilzzz(4)=curconduc(4)
	coilzzz(5)=curconduc(5)
	coilzzz(6)=curconduc(6)
	coilzzz(7)=curconduc(7)
	coilzzz(8)=curconduc(8)
	coilzzz(9)=curconduc(9)
	coilzzz(10)=curconduc(10)
	coilzzz(11)=curconduc(11)
	coilzzz(12)=curconduc(12)

	coilzzz(1:12)=coilzzz(1:12)*1.e3
	else
	
	coilzzz(1:nccc)=curconduc(1:nccc)
	coilzzz(1:nccc)=coilzzz(1:nccc)*1.e3
	

	endif

	cur_con_old(1:nconduc)=curconduc(1:nconduc)
!	voltage_old(1:nconduc)=voltage(1:nconduc)

	if (fast_mode.eq.0) then
			psi_cur_old(1:nconduc)=psiplasmatoconduc(1:nconduc)
	endif
		



      RETURN
      END






!circuit eq advance
	subroutine circuit_eq_advance_ef(j_init)

	use ef_circuit
	use exchange_with_astra       ! declaration of minimal CPOs

	implicit none
	integer j_init,i,ic,j,k,iii,jjj,invertcommand,i_equivalence	
	double precision restemp(200,200),indtemp(200,200)
	double precision dpctemp(i_dim1),vtemp(i_dim1)
	double precision curotemp(i_dim1),curtemp(i_dim1)
	integer rem_coils(200),i_cnew,firstcall
	data firstcall/0/
	save restemp,indtemp,rem_coils,i_cnew,firstcall

	tau_new=tau_circuit_ef
	invertcommand=0
	
	if (tau_new.ne.tau_old) then
		invertcommand=1
		tau_old=tau_new
	endif
	if (j_init.eq.0.or.firstcall.eq.0) then
		invertcommand=1
		tau_old=tau_new
		rem_coils=0
	endif

	ic=nconduc

!	call plasma_psi_to_coils_ef
	
	if (j_init.eq.0) cur_con_old=curconduc
	if (j_init.eq.0) i_cnew=ic
	if (j_init.eq.0) return

	firstcall=1

	dpc(1:ic)=GPI2*(psiplasmatoconduc(1:ic)-psi_cur_old(1:ic))/tau_gseq_ef !plasma contribution

	do j=1,nconduc
		if (activate_coil_ef(j).eq.0) cur_con_old(j)=0.
		if (activate_coil_ef(j).eq.0) dpc(j)=0.
	enddo





	if (reconnect_circuits.eq.1) then
	invertcommand=1
	i_cnew=ic
	indtemp(1:ic,1:ic)=indconduc(1:ic,1:ic)
	restemp(1:ic,1:ic)=resconduc(1:ic,1:ic)
	do k=1,n_equivalence
	i_equivalence=1
		do i=1,nconduc
			if (new_equivalence(i,k).gt.0.and.i_equivalence.eq.1) then
				i_equivalence=0

				do j=i+1,nconduc
					if (new_equivalence(j,k).eq.new_equivalence(i,k)) then
						indtemp(i,i)=indconduc(i,i)+indconduc(j,j)+indconduc(i,j)+indconduc(j,i)
						restemp(i,i)=resconduc(i,i)+resconduc(j,j)+resconduc(i,j)+resconduc(j,i)
						do iii=1,nconduc
							if (iii.ne.i.and.iii.ne.j) then
								indtemp(i,iii)=indconduc(i,iii)+indconduc(j,iii)
								restemp(i,iii)=resconduc(i,iii)+resconduc(j,iii)
								indtemp(iii,i)=indconduc(iii,i)+indconduc(iii,j)
								restemp(iii,i)=resconduc(iii,i)+resconduc(iii,j)
							endif
						enddo
						rem_coils(j)=new_equivalence(j,k)
					ic=ic-1
					endif			
				enddo	
!reove obsolete columns and rows
	do jjj=1,nconduc
		j=nconduc-jjj+1
		if (rem_coils(j).gt.0) then
						indtemp(1:nconduc,j:nconduc)=indtemp(1:nconduc,j+1:nconduc+1)
						restemp(1:nconduc,j:nconduc)=restemp(1:nconduc,j+1:nconduc+1)
						indtemp(j:nconduc,1:nconduc)=indtemp(j+1:nconduc+1,1:nconduc)
						restemp(j:nconduc,1:nconduc)=restemp(j+1:nconduc+1,1:nconduc)
		endif
	enddo

			endif
		enddo	
	enddo
		i_cnew=ic
	endif



	if (use_reduce_circuit.eq.1) then
	ic=nconduc
	dpctemp(1:ic)=dpc(1:ic)
	vtemp(1:ic)=voltage(1:ic)
	curotemp(1:ic)=cur_con_old(1:ic)
	do k=1,n_equivalence
	i_equivalence=1
		do i=1,nconduc
			if (new_equivalence(i,k).gt.0.and.i_equivalence.eq.1) then
				i_equivalence=0
				do jjj=i+1,nconduc
					j=nconduc-jjj+i+1
					if (new_equivalence(j,k).eq.new_equivalence(i,k)) then
						dpctemp(i)=dpctemp(i)+dpctemp(j)
					ic=ic-1
					endif			
				enddo	
			endif
		enddo	
	enddo
	do jjj=1,nconduc
		j=nconduc-jjj+1
		if (rem_coils(j).gt.0) then
						dpctemp(j:nconduc)=dpctemp(j+1:nconduc+1)	
						vtemp(j:nconduc)=vtemp(j+1:nconduc+1)	
						curotemp(j:nconduc)=curotemp(j+1:nconduc+1)
		endif
	enddo
	endif




	i=i_cnew
!	write(*,*) 'voltages ',voltage(1:15),cur_con_old(1:15)

!	write(*,*) 'mmm',use_reduce_circuit,i,indconduc(5,5), & 
!	& resconduc(5,5),cur_con_old(5),curconduc(5),voltage(5),dpc(5),tau_new,invertcommand

	if (use_reduce_circuit.eq.0) then
		call solve_circuit_equations(i,& 
		& indconduc(1:i,1:i),resconduc(1:i,1:i), & 
		& cur_con_old(1:i),curconduc(1:i),& 
		& 	voltage(1:i),dpc(1:i),tau_new,invertcommand)
	else
!	write(*,*) 'ic',i_cnew,indtemp(1,1),restemp(1,1),vtemp(1)
		call solve_circuit_equations(i,& 
		& indtemp(1:i,1:i),restemp(1:i,1:i), & 
		& curotemp(1:i),curtemp(1:i),& 
		& 	vtemp(1:i),dpctemp(1:i),tau_new,invertcommand)	
!	write(*,*) curotemp(1),curtemp(1)
!readapt currents
		do i=1,nconduc
			if (rem_coils(i).eq.0) then
					curconduc(i)=curtemp(i)
			endif
			if (rem_coils(i).gt.0) then
				curconduc(i)=curtemp(rem_coils(i))
				curtemp(i+1:i_cnew+1)=curtemp(i:i_cnew)
			endif
		enddo
	endif
!	write(*,*) 'cur ',curconduc(1:15)

	ic=nconduc
	i=ic
	do j=1,ic
		if (activate_coil_ef(j).eq.0) curconduc(j)=0.
		curconduc(j)=max(curconduc(j),current_limit_ef(j,2))
		curconduc(j)=min(curconduc(j),current_limit_ef(j,1))
		if (sum(force_coil(j,1:i)).gt.0.5) then
			do k=1,i
				if (force_coil(j,k).eq.1.) then
					curconduc(j)=curconduc(k)
				endif
			enddo
		endif
	enddo

	reconnect_circuits=0
	
	return
	end
	


	





	subroutine solve_circuit_equations(nc,im,rm,I0,I1,& 
			& 	V,dpc,tau,invertcommand)
	
	implicit none
	
	integer i,j,k,nc,invertcommand
	double precision im(nc,nc),rm(nc,nc),i0(nc), &
    & i1(nc),v(nc),dpc(nc),tau,matrix(nc,nc)	
	double precision b(nc),invmatrix(200,200)
	save invmatrix

!equation is im*(i1-i0)/tau + rm*i1 = v-dpc

	do i=1,nc
	b(i) = v(i)-dpc(i)+sum(im(i,1:nc)*i0(1:nc))/tau
	enddo
	
	
	if (invertcommand.eq.1) then
	matrix(1:nc,1:nc) = im(1:nc,1:nc)/tau+rm(1:nc,1:nc)
	call	inverse_matrix_equilef(matrix,invmatrix(1:nc,1:nc),nc)
	else
	endif

!	write(*,*) 'ii',i0(12:15)
	do i=1,nc
	i1(i) = sum(invmatrix(i,1:nc)*b(1:nc))
	enddo
!	write(*,*) 'ii',i1(12:15)


!	open(32,file='fort.33872')
!	write(32,*) i0(1),i1(1),b(1),v(1),dpc(1),sum(im(1,1:nc)*i0(1:nc))/tau
!	close(32)



	return
	end
