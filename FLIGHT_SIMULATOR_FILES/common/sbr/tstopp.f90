	subroutine tstopp(tennd)
	
	
	use parameter_inc
	use const_inc
	use status_inc
	use outcmn_inc

		implicit none

	double precision tennd

!	write(*,*) 'tend',TIME,tennd


	if (TIME.ge.tennd) then
	call	CPUSE(6)
	call a_stop
	endif
	if (TIME.ge.tennd-2*tau) then
	write(*,*) 'stopping sim by end time', &
     & time,tennd,tennd-2*tau,tau
	call	err_catch_a
	endif
!	write(*,*) 'tend2',tennd
	

!	write(888,'(111E25.11)') time,updwn,ccoil(9)/100,
!     & vcoil(9)/1e4,ccoil(1)

	


	return




	end
