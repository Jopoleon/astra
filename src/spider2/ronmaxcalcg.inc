!calculate rhonmaxg
!	write(*,*) 'calc ronmax'
			do j=1,nlimadap
				if ((zlimadap(j)-zm).ge.0.and.
     &    (rlimadap(j)-rm).ge.0.) then
				th_rzt(j)=
     & atan((zlimadap(j)-zm)/
     &      (rlimadap(j)-rm))
				endif
				if ((zlimadap(j)-zm).ge.0.and.
     &    (rlimadap(j)-rm).lt.0.) then
				th_rzt(j)=3.141592+
     & atan((zlimadap(j)-zm)/
     &      (rlimadap(j)-rm))
				endif
				if ((zlimadap(j)-zm).lt.0.and.
     &    (rlimadap(j)-rm).lt.0.) then
				th_rzt(j)=3.141592+
     & atan((zlimadap(j)-zm)/
     &      (rlimadap(j)-rm))
				endif
				if ((zlimadap(j)-zm).lt.0.and.
     &    (rlimadap(j)-rm).ge.0.) then
				th_rzt(j)=3.141592*2.+
     & atan((zlimadap(j)-zm)/
     &      (rlimadap(j)-rm))
				endif
			enddo
!	write(*,*) 'calc ronmax2'
	jjj=1
		do j=2,nt1

	if (teta(j).gt.2.*3.141592) th_temp=teta(j)-2.*3.141592
	if (teta(j).lt.-2.*3.141592) th_temp=teta(j)+2.*3.141592
!	write(*,*) 'calc ronmax',th_temp,th_rzt(1:nlimadap)

	j3=minloc(abs(th_rzt(1:nlimadap)-th_temp),1)
	jjj=j3		
		
	ron_max_g(j)=((rlimadap(jjj)-rm)**2.+
     & 		(zlimadap(jjj)-zm)**2.)**0.5

			ron_max_g(j)=ron_max_g(j)/ro(i_bc,j)
		enddo


!		write(*,*) th_rzt(1:nlimadap),ron_max_g(2:nt1)
!		call err_catch_a
