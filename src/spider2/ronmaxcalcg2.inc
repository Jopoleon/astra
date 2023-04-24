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
