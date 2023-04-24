!full system advance
	subroutine full_system_advance_ef(j_init)

	use ef_circuit
	use exchange_with_astra       ! declaration of minimal CPOs
	use parameters_a2spider


	implicit none
	integer j_init,j_iter
	double precision error_temp	
	double precision cur_temp(300)	

! time stepping
! at iteration 0, dpsidt = 0
	if (j_init.eq.0) then
	write(*,*) 'init full system'
	!first do full equilibrium solution at time t=0
 		call psi_external_calc_ef
		call solve_gse2d_fbe_full_ef(0)	
		call plasma_psi_to_coils_ef
		psi_cur_old(1:nconduc)=psiplasmatoconduc(1:nconduc)

!		open(32,file='fort.4444')
!	write(32,*) r(1:nr2),z(1:nz2),jrz(1:nr2,1:nz2)
!	close(32)

!	open(32,file='fort.4445')
!	write(32,*) psiextrz(1:nr2,1:nz2)
!	close(32)
!	open(32,file='fort.4447')
!		write(32,*) psirz(1:nr2,1:nz2)
!	close(32)
!	open(32,file='fort.4448')
!		write(32,*) psiplasmatoconduc(1:nconduc)
!	close(32)
!	open(32,file='fort.4449')
!		write(32,*) curconduc(1:nconduc)
!	close(32)
	write(*,*) 'init done'
			call circuit_eq_advance_ef(0)	
			j_init=-1
		return
	endif

! full iterations
!	write(*,*) 'iter full system'



!		open(32,file='fort.44987')
!		write(32,*) psiextrz(1:nr2,1:nz2)
!		close(32)

!		if (execute_plasma.eq.1) then
!			psi_cur_old(1:nconduc)=psiplasmatoconduc(1:nconduc)
!			call solve_gse2d_fbe_full_ef_1turn(1,0,0.,0.)	
!			call plasma_psi_to_coils_ef
!		endif
		

	if (fast_mode.eq.1.and.execute_plasma.eq.1) then
		psi_cur_old(1:nconduc)=psiplasmatoconduc(1:nconduc)
 		call psi_external_calc_ef
		call solve_gse2d_fbe_full_ef_1turn(1,0,0.,0.)	
		call plasma_psi_to_coils_ef
	endif

!	write(12112,*) ' '
!	write(12112,*) time_astra,cur_con_old(1),2*3.141592*(psiplasmatoconduc(1)-psi_cur_old(1)),voltage(1), &
!     & tau_gseq_ef 

	do j_iter=1,2*max_iter

	cur_temp(1:nconduc)=curconduc(1:nconduc)
		call circuit_eq_advance_ef(1)	

!	curconduc(1:nconduc)=0.75*curconduc(1:nconduc)+0.25*cur_temp(1:nconduc)

	if (fast_mode.eq.0) then
 		call psi_external_calc_ef
		!open(32,file='fort.44987')
		!write(32,*) psiextrz(1:nr2,1:nz2)
		!close(32)
		call solve_gse2d_fbe_full_ef_1turn(1,0,0.,0.)	
		call plasma_psi_to_coils_ef
	endif

	error_temp=sum(abs(cur_temp(1:nconduc)-curconduc(1:nconduc)))/(nconduc+err_epsilon)/iplasma
	if (j_iter.gt.100) then
		write(*,*) 'oscillating solution'
		error_temp=0.99*err_circ_plasma_iter
	endif

	if (j_iter.eq.1.and.fast_mode.eq.0) then
		error_temp=100.
	endif

		
!	write(*,*) 'iterations current ',j_iter,rax,zax,error_temp

!	write(8898,*) rax,zax

!	pause
!	if (j_iter.ge.90) pause


	

	write(*,*) 'error',error_temp
	if (error_temp.le.err_circ_plasma_iter) then
		if (execute_plasma.eq.1.or.fast_mode.eq.0) then
			j_init=-1
		endif
!		open(32,file='fort.4444')
!	write(32,*) r(1:nr2),z(1:nz2),jrz(1:nr2,1:nz2)
!	close(32)
!	open(32,file='fort.4445')
!	write(32,*) psiextrz(1:nr2,1:nz2)
!	close(32)
!	open(32,file='fort.4447')
!		write(32,*) psirz(1:nr2,1:nz2)
!	close(32)
!	open(32,file='fort.4448')
!		write(32,*) psiplasmatoconduc(1:nconduc)
!	close(32)
!	open(32,file='fort.4449')
!		write(32,*) curconduc(1:nconduc)
!	close(32)
	write(*,*) 'init done'
!	write(12112,*) curconduc(1),2*3.141592*(psiplasmatoconduc(1)-psi_cur_old(1)),voltage(1), &
!     & tau_gseq_ef 
		return	
	endif



	if (j_iter.gt.max_iter) then
		write(*,*) 'circuit equations not converging, max number of iterations override... stopping',error_temp
		stop
	endif

!	if (fast_mode.eq.1) then
!		if (execute_plasma.eq.1) then
!			j_init=-1
!		endif
!		return
!	endif

	enddo



return
end























	subroutine solve_gse2d_fbe_full_ef(j_init)	

	use ef_circuit
	use exchange_with_astra
	
	implicit none
	
	integer j_init,i,j,k,ii,jj,kk,iii,jjj,kkk
	integer j_iter,j_iter2,iax,jax,j_cyclo
	double precision g(300,300),temp_err,raxold,zaxold
	double precision g2(300,300)
	double precision temp_err2,raxoldo,zaxoldo
	double precision g000(300,300),br1,bz1,dum6,dum7,dum8,dum9,dum10
	double precision dum1,dum2,raxtmp,zaxtmp,dum3,det
	double precision psistab1o,psistab2o,dbrr,dbrz,dbzr,dbzz
	double precision curpastmp(npassive),correction
	double precision curpastmpo(npassive),psicdtmp(300)
	double precision voltdz(300),curref(300),rdtl,psro,pszo,dist1,dist2
	double precision przsav(6),psistabrr,psistabzz,cibapr,cibazr
	double precision rleft,rright,zup,zdown,dcrdr,dcrdz,dczdr,dczdz
	integer jeppa
!first, initialized initial guess coming from prescribed boundary current density: jrhoteta

!	call err_catch_a







	if (refit_mode.eq.-1) then ! 1 turn only
		call solve_gse2d_fbe_full_ef_1turn(j_init,0,0.,0.)	
		return
	endif
























	if (refit_mode.eq.0) then
	!start iterations to find self-consistent solution
	g000(1:nr2,1:nz2)=psiextrz(1:nr2,1:nz2)
		call find_actual_index_ef(raxp,zaxp,iaxis,jaxis)
		rax=r(iaxis)
		zax=z(jaxis)
	raxold=rax
	zaxold=zax
	raxoldo=rax
	zaxoldo=zax
	psistabR=0.
	psistabZ=0.
	psistab1o=0.
	psistab2o=0.
	temp_err=100.
	temp_err2=100.
	rleft=raxold
	rright=raxold
	zup=zaxold
	zdown=zaxold
! outer cycle, calculate new axis
	j_cyclo=0

!	raxold=9.2466
!	zaxold=0.1285
	jeppa=1
	dist1=dr*dr_factor_init
	dist2=dz*dz_factor_init

	do j_iter2=1,10000000

	write(22291,'(4E25.11)') raxold,zaxold,psistabr,psistabz


	if (j_cyclo.eq.1) then
		raxold=raxoldo+dist1
		zaxold=zaxoldo
		psistab1o=psistabr
		psistab2o=psistabz
	endif
	if (j_cyclo.eq.2) then
		raxold=raxoldo
		zaxold=zaxoldo+dist2
		dcrdr=(psistabr-psistab1o)/dist1
		dczdr=(psistabz-psistab2o)/dist1
	endif

	if (j_cyclo.eq.3) then      ! inverse of (dcrdr dcrdz  ; dczdr dczdz ) = ( dczdz -dcrdz ; -dczdr dcrdr)
		dcrdz=(psistabr-psistab1o)/dist2
		dczdz=(psistabz-psistab2o)/dist2
		det=(dcrdr*dczdz-dcrdz*dczdr)
		cibapr=(dczdz*psistab1o-dcrdz*psistab2o)/det
		cibazr=(-dczdr*psistab1o+dcrdr*psistab2o)/det
!		cibapr=cibapr/(abs(cibapr)+1.e-16)*min(abs(cibapr),dr)
!		cibazr=cibazr/(abs(cibazr)+1.e-16)*min(abs(cibazr),dz)
		raxold=raxoldo-cibapr
		zaxold=zaxoldo-cibazr
!		dist1=max(dr/10.,min(dr,abs(cibapr)))
!		dist2=max(dz/10.,min(dz,abs(cibazr)))
		raxoldo=raxold
		zaxoldo=zaxold
	endif
	
	write(2222,*) ' '
	write(2222,*) j_iter2,j_cyclo,raxoldo,zaxoldo,raxold,zaxold,psistabr,psistabz,dcrdr,dcrdz,dczdr,dczdz,cibapr,cibazr,det,dist1,dist2

!	pause
!	raxoldo=raxold
!	zaxoldo=zaxold	

!inner cycle, calculate psi1 and psi2

	psro=1000.
	pszo=1000.
	redo_bnd=1


	do j_iter=1,100000000

	raxtmp=trax
	zaxtmp=tzax

!	open(5671,file='fort.5671')
			call solve_gse2d_fbe_full_ef_1turn(j_iter-1+j_iter2-1,1,raxold,zaxold)	
!	close(5671)
!	pause

	temp_err=(abs(psro-psistabr)+abs(pszo-psistabz))

!	write(*,*)  'psistab',j_iter,raxold,zaxold,rax,zax,temp_err,psistabR,psistabZ
!	write(444,'(8E25.11)')  trax,tzax,rax,zax,psistabR,psistabZ
	
	if (temp_err.le.err_find_psistab)	goto 941

	if (j_iter.ge.400) then	
		write(*,*) 'oscillating solution'
		pause	
	endif
	
	psro=psistabr
	pszo=psistabz
	
	enddo
	write(*,*) 'total iterations passed!'
	stop
	return


941 continue

			write(65722,'(333E25.11)') rax,zax,psistabr,psistabz


	if (j_cyclo.eq.3) then	
		temp_err2=abs(psistabr)+abs(psistabZ)
		write(3217,'(4E25.11)')  rax,zax,psistabr,psistabz
		j_cyclo=1
		write(*,*) n_of_newton_iterations,nint((0.+j_iter2)/3.)
!		pause
		if (nint((0.+j_iter2)/3.).ge.n_of_newton_iterations) goto 609
	else
		j_cyclo=j_cyclo+1
	endif

	if (j_iter2.ge.400000)	goto 609
	if (temp_err2.le.err_find_psistab)		goto 609



	enddo


609 continue

!		open(32,file='fort.4444')
!	write(32,*) r(1:nr2),z(1:nz2),jrz(1:nr2,1:nz2)
!	close(32)
!	open(32,file='fort.4445')
!	write(32,*) psiextrz(1:nr2,1:nz2)
!	close(32)
!	open(32,file='fort.4447')
!		write(32,*) psirz(1:nr2,1:nz2)
!	close(32)!

!	write(*,*) 'code converged',j_iter2
!	pause
	



	write(*,*) 'total iterations passed!'

	return

	
	
	endif	









	if (refit_mode.eq.101) then !only vertical stab
	!start iterations to find self-consistent solution
	g000(1:nr2,1:nz2)=psiextrz(1:nr2,1:nz2)
		call find_actual_index_ef(raxp,zaxp,iaxis,jaxis)
		rax=r(iaxis)
		zax=z(jaxis)
	raxold=rax
	zaxold=zax
	raxoldo=rax
	zaxoldo=zax
	psistabR=0.
	psistabZ=0.
	psistab1o=0.
	psistab2o=0.
	temp_err=100.
	temp_err2=100.
	rleft=raxold
	rright=raxold
	zup=zaxold
	zdown=zaxold
! outer cycle, calculate new axis
	j_cyclo=0

!	raxold=9.2466
!	zaxold=0.1285
	jeppa=1
	dist1=dr
	dist2=dz

	do j_iter2=1,10000000
	if (j_iter2.gt.250) stop

	write(22291,'(4E25.11)') raxold,zaxold,psistabr,psistabz


	if (j_cyclo.eq.1) then
		zaxold=zaxoldo+dist2
		psistab2o=psistabz
	endif
	if (j_cyclo.eq.2) then      ! inverse of (dcrdr dcrdz  ; dczdr dczdz ) = ( dczdz -dcrdz ; -dczdr dcrdr)
		dczdz=(psistabz-psistab2o)/dist2
		cibazr=psistab2o/dczdz
		zaxold=zaxoldo-cibazr
		zaxoldo=zaxold
	endif
	
	write(2222,*) ' '
	write(2222,*) j_iter2,j_cyclo,raxoldo,zaxoldo,raxold,zaxold,psistabr,psistabz,dcrdr,dcrdz,dczdr,dczdz,cibapr,cibazr,det,dist1,dist2

!	pause
!	raxoldo=raxold
!	zaxoldo=zaxold	

!inner cycle, calculate psi1 and psi2

	psro=1000.
	pszo=1000.
	redo_bnd=1


	do j_iter=1,100000000

	raxtmp=trax
	zaxtmp=tzax

!	open(5671,file='fort.5671')
			call solve_gse2d_fbe_full_ef_1turn(j_iter-1+j_iter2-1,1,-1.d6,zaxold)	
!	close(5671)
!	pause

	temp_err=(abs(pszo-psistabz))

!	write(*,*)  'psistab',j_iter,raxold,zaxold,rax,zax,temp_err,psistabR,psistabZ
!	write(444,'(8E25.11)')  trax,tzax,rax,zax,psistabR,psistabZ
	
	if (temp_err.le.err_find_psistab)	goto 9411

	if (j_iter.ge.400) then	
		write(*,*) 'oscillating solution'
		pause	
	endif
	
	psro=psistabr
	pszo=psistabz
	
	enddo
	write(*,*) 'total iterations passed!'
	stop
	return


9411 continue

			write(65722,'(333E25.11)') rax,zax,psistabr,psistabz


	if (j_cyclo.eq.2) then	
		temp_err2=abs(psistabZ)
		write(3217,'(4E25.11)')  rax,zax,psistabr,psistabz
		j_cyclo=1
		write(*,*) n_of_newton_iterations,nint((0.+j_iter2)/2.)
!		pause
		if (nint((0.+j_iter2)/2.).ge.n_of_newton_iterations) goto 6091
	else
		j_cyclo=j_cyclo+1
	endif

	if (j_iter2.ge.400000)	goto 6091
	if (temp_err2.le.err_find_psistab)		goto 6091



	enddo


6091 continue

!		open(32,file='fort.4444')
!	write(32,*) r(1:nr2),z(1:nz2),jrz(1:nr2,1:nz2)
!	close(32)
!	open(32,file='fort.4445')
!	write(32,*) psiextrz(1:nr2,1:nz2)
!	close(32)
!	open(32,file='fort.4447')
!		write(32,*) psirz(1:nr2,1:nz2)
!	close(32)!

!	write(*,*) 'code converged',j_iter2
!	pause
	



	write(*,*) 'total iterations passed!'

	return

	
	
	endif	
















	if (refit_mode.eq.1) then
		call restab_axis_with_furier_wall
		return	
	endif	


	if (refit_mode.eq.2) then
		call restab_boundary_with_furier_wall !doesnt work well
		return	
	endif	


	if (refit_mode.eq.3) then
		call restab_F_function_full_fonfit
		return	
	endif	




	return
	end

	
	
	
	subroutine restab_F_function_full_fonfit
	
	!refits all currents
	
	use ef_circuit
	use exchange_with_astra
	use green_matrix

	implicit none
	
	integer j_init,i,j,k,ii,jj,kk,iii,jjj,kkk
	integer j_iter,j_iter2,iax,jax
	double precision g(300,300),temp_err,raxold,zaxold
	double precision g000(300,300),br1,bz1,dum6,dum7,dum8,dum9,dum10
	double precision psistab1,psistab2,dum1,dum2,zum1,zum2
	double precision psistab1o,psistab2o,dbrr,dbrz,dbzr,dbzz
	double precision voltdz(300),rdtl,tin,tup,cum1,cum2
	double precision bub(9),xub(9),yub(9),ccc(6),ddipsi(8)
	double precision G_00(nconduc,nteta)
	double precision G_00c(nconduc) !average
	double precision G_00r(nconduc)
	double precision G_00z(nconduc)
	double precision f_correction,x1,x2,x3,x4
	double precision psibt0,psibt1,psicorr(nteta)
	double precision ggep(nteta)
	double precision matrix(nconduc,nconduc)
	double precision invmatrix(nconduc,nconduc)
	double precision Fderiv(nconduc),Ffunc,Ffunc_old
	integer info,whichcoil
	double precision sigma_B, sigma_axis,sigma_coils(nconduc)
	double precision curref(nconduc),curnow(nconduc),curdiff(nconduc)
	double precision raxref,zaxref,rbref(500),zbref(500)
	data cum1/0./
	data cum2/0./
	save cum1,cum2,dum8,dum9

	
	
	psicorr=0.
	
	Ffunc_old=1.e6

!first, initialized initial guess coming from prescribed boundary current density: jrhoteta
	write(*,*) 'reinterp curr, restab'
		call interp_j_fromrhotorz
!rescale plasma
		dum1=0.
		do j=1,nz2
		do i=1,nr2
			dum1=dum1+jrz(i,j)*dr*dz
		enddo
		enddo
		jrz=jrz/dum1*iplasma

		rax=raxp
		zax=zaxp
		write(*,*) raxp,zaxp,dum1




		call find_actual_index_ef(rax,zax,iaxis,jaxis)
		iax=iaxis
		jax=jaxis

!axis is given by raxp, zaxp, boundary by rbndp, zbndp, coilref by curconduc(passive)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	

!read efonfit.dat
	open(32,file=data_dir(1:data_dir_k)//'efonfit.dat')
	!read sigma_B, sigma_axis, sigma_coils
	read(32,*) sigma_B, sigma_axis
	read(32,*) sigma_coils(1:nactive)  !active conductors
	read(32,*) sigma_coils(nactive+1)  !passive conductors
	close(32) 
	sigma_coils(nactive+1:nconduc)=sigma_coils(nactive+1)

	curref(1:nconduc)=curconduc(1:nconduc)
	curnow(1:nconduc)=curconduc(1:nconduc)
	curdiff=0.
	
	raxref=raxp
	zaxref=zaxp
	rbref(1:nteta)=rbndp(1:nteta)
	zbref(1:nteta)=zbndp(1:nteta)


! calculate the matrix F_li of the F function, including the green function terms
	matrix=0.
	invmatrix=0.
	do j=1,nconduc
		do k=1,nteta
			call find_fields_interp_ef_green(rbref(k),zbref(k),G_00(j,k),j) !give back psi,br,bz at r0,z0
		enddo
			G_00c(j)=sum(G_00(j,1:nteta))/(0.+nteta)
			call find_fields_interp_ef_green(raxref-dr/2.,zaxref,bub(1),j) !give back psi,br,bz at r0,z0
			call find_fields_interp_ef_green(raxref+dr/2.,zaxref,bub(2),j) !give back psi,br,bz at r0,z0
			call find_fields_interp_ef_green(raxref,zaxref-dz/2.,bub(3),j) !give back psi,br,bz at r0,z0
			call find_fields_interp_ef_green(raxref,zaxref+dz/2.,bub(4),j) !give back psi,br,bz at r0,z0
			G_00r(j)=(bub(2)-bub(1))/dr
			G_00z(j)=(bub(4)-bub(3))/dz
	enddo

	do j=1,nconduc
		do i=1,nconduc
			if (i.eq.j) matrix(i,j)=matrix(i,j)+2.*sigma_coils(i)
			matrix(i,j)=matrix(i,j)+ & 
			& 2.*sigma_B*sum((G_00(i,1:nteta)-G_00c(i))*(G_00(j,1:nteta)-G_00c(j)))+ &
		& 2.*sigma_axis*(G_00r(i)*G_00r(j)+G_00z(i)*G_00z(j))
		enddo
	enddo
		
!calculate inverse
	call inverse_matrix_equilef(matrix,invmatrix,nconduc)



	do j_iter=1,300000 !iterations to find currents

	if (j_iter.gt.50) stop

!		open(32,file='fort.44490')
!		write(32,*) jrz(1:nr2,1:nz2)
!	close(32)
	CALL CPU_TIME(tin)	
	g=0.
	call solve_gs2d(g) !jrz as right hand side
	CALL CPU_TIME(tup)
	!write(*,*) 'solve g0',tup-tin	
	!	open(32,file='fort.44491')
	!	write(32,*) g(1:nr2,1:nz2)
	!close(32)

!	write(*,*) 'stop here'
!	call err_catch_a
	write(*,*) 'stop here2'
	!find boundary condition using g
	CALL CPU_TIME(tin)	
	call boundary_ef(g)  ! gbound = integral (Green*dg/dn) over the boundary
	write(*,*) 'stop here3'
	CALL CPU_TIME(tup)
	write(*,*) tup-tin	

	CALL CPU_TIME(tin)	
	call solve_gs2d(g) ! again jrz as right hand side
	psiplasrz(1:nr2,1:nz2)=g(1:nr2,1:nz2)
	write(*,*) 'stop here4'
	CALL CPU_TIME(tup)
	write(*,*) tup-tin	
		!open(32,file='fort.4446')
		!write(32,*) g(1:nr2,1:nz2)
	!close(32)
	


	

!	write(*,*) 'stop here5'
!	open(32,file='fort.4447')
!		write(32,*) psirz(1:nr2,1:nz2)
!	close(32)

	


	CALL CPU_TIME(tin)	


	!construct correction
	do j=1,nz2
	do i=1,nr2
		f_correction=0.
		do k=1,nconduc
			f_correction=f_correction+curdiff(k)*greeni(i,j,k) 
		enddo
		psirz(i,j)=psiplasrz(i,j)+psiextrz(i,j)+f_correction !total flux
	enddo
	enddo

	
	CALL CPU_TIME(tup)
	write(*,*) tup-tin	
	

	
	
			do j=1,nteta
			xub(1)=rbref(j)
			yub(1)=zbref(j)
			call find_fields_interp_ef_psionly(xub(1),yub(1),psicorr(j))  !psi on the boundary
		enddo

	x1=sum(psicorr)/(nteta+0.) !average psi on the boundary
	
	!derivative at ref axis
			call find_fields_interp_ef_psionly(raxref-dr/2.,zaxref,bub(1)) !give back psi,br,bz at r0,z0
			call find_fields_interp_ef_psionly(raxref+dr/2.,zaxref,bub(2)) !give back psi,br,bz at r0,z0
			call find_fields_interp_ef_psionly(raxref,zaxref-dz/2.,bub(3)) !give back psi,br,bz at r0,z0
			call find_fields_interp_ef_psionly(raxref,zaxref+dz/2.,bub(4)) !give back psi,br,bz at r0,z0
	x2 = (bub(2)-bub(1))/dr ! psir
	x3 = (bub(4)-bub(3))/dz ! psiz


	Ffunc=sigma_B*sum((psicorr-x1)**2.)+sum(sigma_coils*curdiff**2.)+sigma_axis*(x2**2.+x3**2.)

!calculate F derivative
	do i=1,nconduc
		Fderiv(i)=2.*(sigma_coils(i)*curdiff(i)+ & 
		& sigma_B*sum((psicorr-x1)*(G_00(i,1:nteta)-G_00c(i)))+ &
		& sigma_axis*(x2*G_00r(i)+x3*G_00z(i)) & 
		& )
	enddo
	
!calculate new currents

	do i=1,nconduc
		curnow(i)=curnow(i)-sum(invmatrix(i,1:nconduc)*Fderiv)
	enddo	
		curdiff=curnow-curref

!	write(16761,*) j_iter
!	write(16761,*) curdiff(1:nconduc)


	!construct correction
	do j=1,nz2
	do i=1,nr2
		f_correction=0.
		do k=1,nconduc
			f_correction=f_correction+curdiff(k)*greeni(i,j,k) 
		enddo
		psirz(i,j)=psiplasrz(i,j)+psiextrz(i,j)+f_correction !total flux
	enddo
	enddo

	
	CALL CPU_TIME(tup)
	write(*,*) tup-tin	
	

	write(*,*) 'stop here41'
	CALL CPU_TIME(tin)	
	call find_new_axis_part1	
	CALL CPU_TIME(tup)
	write(*,*) tup-tin	

	write(*,*) j_iter,rax,zax,raxp,zaxp,psistabr,psistabz

	CALL CPU_TIME(tin)	
	write(*,*) 'stop here55'
	call find_psi_boundary
	CALL CPU_TIME(tup)
	write(*,*) tup-tin	

	CALL CPU_TIME(tin)	
	call new_jrz_ef  ! calculate new right hand side
	write(*,*) 'stop here7'
	CALL CPU_TIME(tup)
	write(*,*) 'current',tup-tin	
!	write(*,*) psiaxis,psibnd
!	open(32,file='fort.44442')
!		write(32,*) jrz(1:nr2,1:nz2)
!	close(32)

	temp_err=abs(Ffunc-Ffunc_old)

	Ffunc_old=Ffunc

	write(*,*) 'temp err',temp_err,Ffunc
!	write(16762,'(335553E25.11)') rax,zax,(Ffunc),raxp,zaxp,rbref(1:nteta),zbref(1:nteta),curref(1:nconduc),curdiff(1:nconduc), & 
!	 & temp_err,sum(abs(Fderiv))

	
	if (temp_err.le.err_find_psistab) goto 1032

	enddo

	
1032 continue

	write(*,*) raxp,zaxp,rax,zax

!check
	do i=1,nconduc
		curconduc(i)= curnow(i)
	enddo

	call psi_external_calc_ef
	psirz(1:nr2,1:nz2)=psiplasrz(1:nr2,1:nz2)+psiextrz(1:nr2,1:nz2)
	call find_new_axis_part1	
	call find_psi_boundary
	call new_jrz_ef  ! calculate new right hand side

	write(*,*) curconduc(1:nconduc),rax,zax	
	write(*,*) 'full fonfit eddy currents converged'




!	pause(0.1)



	return
	end











	subroutine restab_boundary_with_furier_wall !not working well
	
	use ef_circuit
	use exchange_with_astra
	use green_matrix

	implicit none
	
	integer j_init,i,j,k,ii,jj,kk,iii,jjj,kkk
	integer j_iter,j_iter2,iax,jax
	double precision g(300,300),temp_err,raxold,zaxold
	double precision g000(300,300),br1,bz1,dum6,dum7,dum8,dum9,dum10
	double precision psistab1,psistab2,dum1,dum2,zum1,zum2
	double precision psistab1o,psistab2o,dbrr,dbrz,dbzr,dbzz
	double precision curpastmp(npassive),correction,delr,delz
	double precision curpastmpo(npassive),psicdtmp(300),anglr(npassive)
	double precision voltdz(300),curref(300),rdtl,tin,tup,cum1,cum2
	double precision bub(9),xub(9),yub(9),ccc(6),ddipsi(8)
	double precision g_0(npassive),g0_r(npassive),g0_z(npassive)
	double precision S_00(n_fourier_restab_boundary),C_00(n_fourier_restab_boundary)
	double precision G_00c(nteta,n_fourier_restab_boundary)
	double precision G_00s(nteta,n_fourier_restab_boundary)
	double precision f_correction,x1,x2,x3,x4
	double precision psibt0,psibt1,psicorr(nteta)
	double precision matrix(nteta,2*n_fourier_restab_boundary)
	integer info
	double precision work(4*nteta*n_fourier_restab_boundary)
	data cum1/0./
	data cum2/0./
	save cum1,cum2,dum8,dum9




!first, initialized initial guess coming from prescribed boundary current density: jrhoteta
	write(*,*) 'reinterp curr, restab'
		call interp_j_fromrhotorz
!rescale plasma
		dum1=0.
		do j=1,nz2
		do i=1,nr2
			dum1=dum1+jrz(i,j)*dr*dz
		enddo
		enddo
		jrz=jrz/dum1*iplasma

		rax=raxp
		zax=zaxp
		write(*,*) raxp,zaxp,dum1


		call find_actual_index_ef(rax,zax,iaxis,jaxis)
		iax=iaxis
		jax=jaxis
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	G_00c=0.
	G_00s=0.
	
! evaluate coils things
	do i=1,npassive
		call find_angle_ef(rax,zax,r_cond(nactive+i),z_cond(nactive+i),anglr(i))
	!find true axis
	do j=1,nteta
		xub(1)=rbndp(j)
		yub(1)=zbndp(j)
		call find_fields_interp_ef_green(xub(1),yub(1),bub(2),nactive+i) !give back psi,br,bz at r0,z0
		do k=1,n_fourier_restab_boundary
			G_00c(j,k)=G_00c(j,k)+cos(k*anglr(i))*bub(2)
			G_00s(j,k)=G_00s(j,k)+sin(k*anglr(i))*bub(2)
		enddo
	enddo
	enddo
	
	do j=1,nteta
				do k=1,n_fourier_restab_boundary
		matrix(j,k)=G_00c(j,k)
		matrix(j,n_fourier_restab_boundary+k)=G_00s(j,k)
	enddo
	enddo
	
	psibt0=1000.
	psibt1=1000.
		
	do j_iter=1,30

!		open(32,file='fort.44490')
!		write(32,*) jrz(1:nr2,1:nz2)
!	close(32)
	CALL CPU_TIME(tin)	
	g=0.
	call solve_gs2d(g) !jrz as right hand side
	CALL CPU_TIME(tup)
	write(*,*) 'solve g0',tup-tin	
		!open(32,file='fort.44491')
	!	write(32,*) g(1:nr2,1:nz2)
	!close(32)

!	write(*,*) 'stop here'
!	call err_catch_a
	write(*,*) 'stop here2'
	!find boundary condition using g
	CALL CPU_TIME(tin)	
	call boundary_ef(g)  ! gbound = integral (Green*dg/dn) over the boundary
	write(*,*) 'stop here3'
	CALL CPU_TIME(tup)
	write(*,*) tup-tin	

	CALL CPU_TIME(tin)	
	call solve_gs2d(g) ! again jrz as right hand side
	psiplasrz(1:nr2,1:nz2)=g(1:nr2,1:nz2)
	write(*,*) 'stop here4'
	CALL CPU_TIME(tup)
	write(*,*) tup-tin	
		!open(32,file='fort.4446')
	!	write(32,*) g(1:nr2,1:nz2)
	!close(32)
	


	
	
	

!	write(*,*) 'stop here5'
!	open(32,file='fort.4447')
!		write(32,*) psirz(1:nr2,1:nz2)
!	close(32)

	


	CALL CPU_TIME(tin)	


			psirz(1:nr2,1:nz2)=psiplasrz(1:nr2,1:nz2)+psiextrz(1:nr2,1:nz2) !total flux

		do j=1,nteta
			xub(1)=rbndp(j)
			yub(1)=zbndp(j)
			call find_fields_interp_ef_psionly(xub(1),yub(1),psicorr(j)) 
		enddo

	x1=sum(psicorr)/(nteta+0.)
	psibt1=x1
	psicorr=x1-psicorr
	
!least square fit solution

	call dgels('N', &
     &		nteta,&
     &  	n_fourier_restab_boundary*2,&
     &  	1,&
     &  		matrix,&
     &  	nteta,&
     &  	psicorr,&
     &  	nteta,&
     &  	WORK,&
     &  	2*(nteta)*n_fourier_restab_boundary*2,&
     &  	INFO ) 		 

	!construct correction
	do j=1,nz2
	do i=1,nr2
		f_correction=0.
		do k=1,n_fourier_restab_boundary
			f_correction=f_correction+ & 
			& psicorr(k)*sum(greeni(i,j,nactive+1:nactive+npassive)*cos(k*anglr(1:npassive)))	+ & 
			& psicorr(n_fourier_restab_boundary+k)*sum(greeni(i,j,nactive+1:nactive+npassive)*sin(k*anglr(1:npassive)))
		enddo
		psirz(i,j)=psiplasrz(i,j)+psiextrz(i,j)+f_correction !total flux
	enddo
	enddo




	
	CALL CPU_TIME(tup)
	write(*,*) tup-tin	
	

	write(*,*) 'stop here41'
	CALL CPU_TIME(tin)	
	call find_new_axis_part1	
	CALL CPU_TIME(tup)
	write(*,*) tup-tin	

	write(*,*) j_iter,rax,zax,raxp,zaxp,psistabr,psistabz

	CALL CPU_TIME(tin)	
	write(*,*) 'stop here55'
	call find_psi_boundary
	CALL CPU_TIME(tup)
	write(*,*) tup-tin	

	CALL CPU_TIME(tin)	
	call new_jrz_ef  ! calculate new right hand side
	write(*,*) 'stop here7'
	CALL CPU_TIME(tup)
	write(*,*) 'current',tup-tin	
!	write(*,*) psiaxis,psibnd
!	open(32,file='fort.44442')
!		write(32,*) jrz(1:nr2,1:nz2)
!	close(32)

	temp_err=abs(psibt0-psibt1)

	write(*,*) 'temp err',temp_err,psibt0,psibt1

	psibt0=psibt1


	
	if (temp_err.le.err_find_psistab) goto 1032

	enddo

	
1032 continue

	write(*,*) raxp,zaxp,rax,zax,psistabr,psistabz	

!check
	do i=1,npassive
		do j=1,n_fourier_restab_boundary
		curconduc(nactive+i)= curconduc(nactive+i)+& 
		 & C_00(j)*cos(j*anglr(i))+S_00(j)*sin(j*anglr(i))
		enddo
	enddo

	call psi_external_calc_ef
	psirz(1:nr2,1:nz2)=psiplasrz(1:nr2,1:nz2)+psiextrz(1:nr2,1:nz2)
	call find_new_axis_part1	
	call find_psi_boundary
	call new_jrz_ef  ! calculate new right hand side

	write(*,*) curconduc(1:nconduc),rax,zax	

		stop
	return
	end






	subroutine restab_axis_with_furier_wall
	
	use ef_circuit
	use exchange_with_astra
	use green_matrix       ! declaration of minimal CPOs

	implicit none
	
	integer j_init,i,j,k,ii,jj,kk,iii,jjj,kkk
	integer j_iter,j_iter2,iax,jax
	double precision g(300,300),temp_err,raxold,zaxold
	double precision g000(300,300),br1,bz1,dum6,dum7,dum8,dum9,dum10
	double precision psistab1,psistab2,dum1,dum2,zum1,zum2
	double precision psistab1o,psistab2o,dbrr,dbrz,dbzr,dbzz
	double precision curpastmp(npassive),correction,delr,delz
	double precision curpastmpo(npassive),psicdtmp(300),anglr(npassive)
	double precision voltdz(300),curref(300),rdtl,tin,tup,cum1,cum2
	double precision bub(9),xub(9),yub(9),ccc(6),ddipsi(8)
	double precision g_0(npassive),g0_r(npassive),g0_z(npassive)
	double precision S_00r,C_00r,S_00z,C_00z,C_00(258,258),S_00(258,258)
	data cum1/0./
	data cum2/0./
	save cum1,cum2,dum8,dum9




!first, initialized initial guess coming from prescribed boundary current density: jrhoteta
	write(*,*) 'reinterp curr, restab'
		call interp_j_fromrhotorz
!rescale plasma
		dum1=0.
		do j=1,nz2
		do i=1,nr2
			dum1=dum1+jrz(i,j)*dr*dz
		enddo
		enddo
		jrz=jrz/dum1*iplasma

		rax=raxp
		zax=zaxp
		write(*,*) raxp,zaxp,dum1


		call find_actual_index_ef(rax,zax,iaxis,jaxis)
		iax=iaxis
		jax=jaxis
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
! evaluate coils things
	do i=1,npassive
		call find_angle_ef(rax,zax,r_cond(nactive+i),z_cond(nactive+i),anglr(i))
!find true axis
	xub(1)=r(iax-1)
	xub(2)=r(iax)
	xub(3)=r(iax+1)
	xub(4)=r(iax)
	xub(5)=r(iax)
	xub(6)=r(iax-1)
	xub(7)=r(iax-1)
	xub(8)=r(iax+1)
	xub(9)=r(iax+1)
	yub(1)=z(jax)
	yub(2)=z(jax)
	yub(3)=z(jax)
	yub(4)=z(jax-1)
	yub(5)=z(jax+1)
	yub(6)=z(jax-1)
	yub(7)=z(jax+1)
	yub(8)=z(jax-1)
	yub(9)=z(jax+1)
	bub(1)=greeni(iax-1,jax,nactive+i)
	bub(2)=greeni(iax,jax,nactive+i)
	bub(3)=greeni(iax+1,jax,nactive+i)
	bub(4)=greeni(iax,jax-1,nactive+i)
	bub(5)=greeni(iax,jax+1,nactive+i)
	bub(6)=greeni(iax-1,jax-1,nactive+i)
	bub(7)=greeni(iax-1,jax+1,nactive+i)
	bub(8)=greeni(iax+1,jax-1,nactive+i)
	bub(9)=greeni(iax+1,jax+1,nactive+i)
		call least_square_biquad_ef(xub,yub,bub,9,ccc,dum1,dum2,zum1,ddipsi)
		g0_r(i)=ddipsi(1)
		g0_z(i)=ddipsi(2)
	enddo
	
	do j=1,nz2
	do i=1,nr2
	C_00(i,j)=sum(greeni(i,j,nactive+1:nactive+npassive)*cos(anglr))
	S_00(i,j)=sum(greeni(i,j,nactive+1:nactive+npassive)*sin(anglr))
	enddo
	enddo
	C_00r=sum(g0_r*cos(anglr))
	S_00r=sum(g0_r*sin(anglr))
	C_00z=sum(g0_z*cos(anglr))
	S_00z=sum(g0_z*sin(anglr))
	
	write(*,*) 'coils structure',c_00(iax,jax),s_00(iax,jax),c_00r,s_00r,c_00z,s_00z,iplasma


	psistab1o=1000.
	psistab2o=1000.
	
	do j_iter=1,30000

!		open(32,file='fort.44490')
!		write(32,*) jrz(1:nr2,1:nz2)
!	close(32)
	CALL CPU_TIME(tin)	
	g=0.
	call solve_gs2d(g) !jrz as right hand side
	CALL CPU_TIME(tup)
	write(*,*) 'solve g0',tup-tin	
		!open(32,file='fort.44491')
		!write(32,*) g(1:nr2,1:nz2)
	!close(32)

!	write(*,*) 'stop here'
!	call err_catch_a
	write(*,*) 'stop here2'
	!find boundary condition using g
	CALL CPU_TIME(tin)	
	call boundary_ef(g)  ! gbound = integral (Green*dg/dn) over the boundary
	write(*,*) 'stop here3'
	CALL CPU_TIME(tup)
	write(*,*) tup-tin	

	CALL CPU_TIME(tin)	
	call solve_gs2d(g) ! again jrz as right hand side
	psiplasrz(1:nr2,1:nz2)=g(1:nr2,1:nz2)
	write(*,*) 'stop here4'
	CALL CPU_TIME(tup)
	write(*,*) tup-tin	
		!open(32,file='fort.4446')
	!	write(32,*) g(1:nr2,1:nz2)
!	close(32)
	


	
	
	

!	write(*,*) 'stop here5'
!	open(32,file='fort.4447')
!		write(32,*) psirz(1:nr2,1:nz2)
!	close(32)

	


	CALL CPU_TIME(tin)	
	psistabr=0.
	psistabz=0.
	delr=0.
	delz=0.
			psirz(1:nr2,1:nz2)=psiplasrz(1:nr2,1:nz2)+psiextrz(1:nr2,1:nz2) !total flux
			call find_new_axis_part1	

	dum1=C_00r*S_00z-C_00z*S_00r	
	call find_fields_interp_ef_psionly(raxp+dr,zaxp,bub(1)) !give back psi,br,bz at r0,z0
	call find_fields_interp_ef_psionly(raxp-dr,zaxp,bub(2)) !give back psi,br,bz at r0,z0
	call find_fields_interp_ef_psionly(raxp,zaxp+dz,bub(3)) !give back psi,br,bz at r0,z0
	call find_fields_interp_ef_psionly(raxp,zaxp-dz,bub(4)) !give back psi,br,bz at r0,z0

	xub(1)=(bub(1)-bub(2))/(2.*dr)
	yub(1)=(bub(3)-bub(4))/(2.*dz)

	delr=(S_00z*(-xub(1))+(-C_00z)*(-yub(1)))/dum1
	delz=(-S_00r*(-xub(1))+C_00r*(-yub(1)))/dum1


		psistabr=delr
		psistabz=delz




	write(*,*) j_iter,rax,zax,raxp,zaxp,psistabr,psistabz


			psirz(1:nr2,1:nz2)=psiplasrz(1:nr2,1:nz2)+psiextrz(1:nr2,1:nz2)+ &
       & psistabr*C_00(1:nr2,1:nz2)+psistabz*S_00(1:nr2,1:nz2) !total flux



	
	CALL CPU_TIME(tup)
	write(*,*) tup-tin	
	

	write(*,*) 'stop here41'
	CALL CPU_TIME(tin)	
	call find_new_axis_part1	
	CALL CPU_TIME(tup)
	write(*,*) tup-tin	

	write(*,*) j_iter,rax,zax,raxp,zaxp,psistabr,psistabz

	CALL CPU_TIME(tin)	
	write(*,*) 'stop here55'
	call find_psi_boundary
	CALL CPU_TIME(tup)
	write(*,*) tup-tin	

	CALL CPU_TIME(tin)	
	call new_jrz_ef  ! calculate new right hand side
	write(*,*) 'stop here7'
	CALL CPU_TIME(tup)
	write(*,*) 'current',tup-tin	
!	write(*,*) psiaxis,psibnd
!	open(32,file='fort.44442')
!		write(32,*) jrz(1:nr2,1:nz2)
!	close(32)

	temp_err=(abs(psistab1o-psistabr)+abs(psistab2o-psistabz))

	write(*,*) 'temp err',temp_err,psistab1o,psistab2o,psistabr,psistabz

	psistab1o=psistabr
	psistab2o=psistabz


	
	if (temp_err.le.err_find_psistab) goto 1032

	enddo

	
1032 continue

	write(*,*) raxp,zaxp,rax,zax,psistabr,psistabz	

!check
	do i=1,npassive
		curconduc(nactive+i)= curconduc(nactive+i)+& 
		 & psistabr*cos(anglr(i))+psistabz*sin(anglr(i))
	enddo

	call psi_external_calc_ef
	psirz(1:nr2,1:nz2)=psiplasrz(1:nr2,1:nz2)+psiextrz(1:nr2,1:nz2)
	call find_new_axis_part1	
	call find_psi_boundary
	call new_jrz_ef  ! calculate new right hand side

	write(*,*) curconduc(1:nconduc),rax,zax	

	write(*,*) 'code restab axis with fourier wall converged'
!	pause
	
	return
	end


















	subroutine solve_gse2d_fbe_full_ef_1turn(j_init,j_stab,raxold,zaxold)	

	use ef_circuit
	use exchange_with_astra
	
	implicit none
	
	integer j_init,i,j,k,ii,jj,kk,iii,jjj,kkk
	integer j_iter,j_iter2,iax,jax,j_stab,j_count
	double precision g(300,300),temp_err,raxold,zaxold
	double precision g000(300,300),br1,bz1,dum6,dum7,dum8,dum9,dum10
	double precision psistab1,psistab2,dum1,dum2,zum1,zum2
	double precision psistab1o,psistab2o,dbrr,dbrz,dbzr,dbzz
	double precision curpastmp(npassive),correction,delr,delz
	double precision curpastmpo(npassive),psicdtmp(300)
	double precision voltdz(300),curref(300),rdtl,tin,tup,cum1,cum2
	data cum1/0./
	data cum2/0./
	data j_count/0/
	save cum1,cum2,dum8,dum9,j_count




!first, initialized initial guess coming from prescribed boundary current density: jrhoteta
	if (j_init.eq.0) then
	write(*,*) 'reinterp curr'
		call interp_j_fromrhotorz
!rescale plasma
		dum1=0.
		do j=1,nz2
		do i=1,nr2
			dum1=dum1+jrz(i,j)*dr*dz
		enddo
		enddo
		jrz=jrz/dum1*iplasma

		rax=raxp
		zax=zaxp
		call find_actual_index_ef(rax,zax,iaxis,jaxis)
		rax=r(iaxis)
		zax=z(jaxis)
		write(*,*) raxp,zaxp,rax,zax
	endif

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	
	
	
	
	


		!open(32,file='fort.44490')
		!write(32,*) jrz(1:nr2,1:nz2)
	!close(32)
	CALL CPU_TIME(tin)	
	g=0.
	call solve_gs2d(g) !jrz as right hand side
	CALL CPU_TIME(tup)
	write(*,*) 'solve g0',tup-tin	
		!open(32,file='fort.44491')
		!	write(32,*) g(1:nr2,1:nz2)
	!close(32)

!	write(*,*) 'stop here'
!	call err_catch_a
	write(*,*) 'stop here2'
	!find boundary condition using g
	CALL CPU_TIME(tin)	
	call boundary_ef(g)  ! gbound = integral (Green*dg/dn) over the boundary
	write(*,*) 'stop here3'
	CALL CPU_TIME(tup)
	write(*,*) tup-tin	

	CALL CPU_TIME(tin)	
	call solve_gs2d(g) ! again jrz as right hand side
	psiplasrz(1:nr2,1:nz2)=g(1:nr2,1:nz2)
	write(*,*) 'stop here4'
	CALL CPU_TIME(tup)
	write(*,*) tup-tin	
		!open(32,file='fort.4446')
		!write(32,*) g(1:nr2,1:nz2)
	!close(32)
	


	
	
	

	write(*,*) 'stop here5'
	!open(32,file='fort.4447')
	!	write(32,*) psiextrz(1:nr2,1:nz2)
	!close(32)

!	write(3321,*) curconduc(1:nconduc)

!	stop
!	call err_catch_a
	


	CALL CPU_TIME(tin)	
	if (j_stab.eq.1) then
	psistabr=0.
	psistabz=0.
	delr=0.
	delz=0.
			psirz(1:nr2,1:nz2)=psiplasrz(1:nr2,1:nz2)+psiextrz(1:nr2,1:nz2) !total flux
	do ii=1,100000
		do i=1,nr2
		do j=1,nz2
			psirz(i,j)=psirz(i,j)+delr*r(i)**2.+delz*z(j) !total flux
		enddo
		enddo
	call find_new_axis_part1	
	if (raxold.lt.-1.d5) then
		delr=-(derivpsi(3)*(rax-rax)+derivpsi(5)*(zaxold-zax))*0.5/rax
		delz=-(derivpsi(4)*(zaxold-zax)+derivpsi(5)*(rax-rax))
	else
		delr=-(derivpsi(3)*(raxold-rax)+derivpsi(5)*(zaxold-zax))*0.5/raxold
		delz=-(derivpsi(4)*(zaxold-zax)+derivpsi(5)*(raxold-rax))
	endif
		psistabr=psistabr+delr
		psistabz=psistabz+delz
	!write(5671,*) 'psi',delr,delz,psistabr,psistabz,raxold,rax,zaxold,zax,derivpsi(3:5)
	if (abs(delr)+abs(delz).gt.100.) then
		write(*,*) 'delr,delz dont converge'
		stop !test, to put this better
	endif
	if (abs(delr)+abs(delz).lt.err_find_delr) goto 2010

	enddo

	write(*,*) 'axis not converged, stop'
	stop

2010 continue

	write(*,*) 'raxx',rax,zax,trax,tzax,raxold,zaxold,psistabr,psistabz
	

	else

		psirz(1:nr2,1:nz2)=psiplasrz(1:nr2,1:nz2)+psiextrz(1:nr2,1:nz2) !total flux

	endif
	CALL CPU_TIME(tup)
	write(*,*) tup-tin	
	
	

	write(*,*) 'stop here41'
	CALL CPU_TIME(tin)	
	call find_new_axis_part1	
	CALL CPU_TIME(tup)
	write(*,*) tup-tin	

	CALL CPU_TIME(tin)	
	write(*,*) 'stop here55'
	call find_psi_boundary
	CALL CPU_TIME(tup)
	write(*,*) tup-tin	

	CALL CPU_TIME(tin)	
	call new_jrz_ef  ! calculate new right hand side
	write(*,*) 'stop here7'
	CALL CPU_TIME(tup)
	write(*,*) 'current',tup-tin,nr2,nz2	
!	write(*,*) psiaxis,psibnd
	!open(32,file='fort.444421')
	!	write(32,*) jrz(1:nr2,1:nz2)
	!close(32)

	



	
	return
	end












