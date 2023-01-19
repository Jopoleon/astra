	jjj=1

	if (teta(j).gt.2.*3.141592) th_temp=teta(j)-2.*3.141592
	if (teta(j).lt.-2.*3.141592) th_temp=teta(j)+2.*3.141592


	j3=minloc(abs(th_rzt(1:nlimadap)-th_temp),1)
	jjj=j3		
!	write(*,*) 'calc ronmax',th_temp,th_rzt(jjj)
		
	rob(j)=((rlimadap(jjj)-rm)**2.+
     & 		(zlimadap(jjj)-zm)**2.)**0.5

!			rob(j)=max(ron_max*ro(iplas,j),rob(j))


!		write(*,*) rob(j)
!		call err_catch_a
