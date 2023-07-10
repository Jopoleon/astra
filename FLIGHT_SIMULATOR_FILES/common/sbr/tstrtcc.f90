	subroutine TSTRTCC(t_startpp,itp_opt)

	use parameter_inc
	use const_inc
	
	implicit none
	double precision t_startpp,itp_opt            

	if (TIME.ge.t_startpp) then
	ITFBP  =  itp_opt
	else
	ITFBP  = 0.
	endif


!	write(*,*) 'tstrtcc ',time,t_startpp,itp_opt,itfbp,ipl,iplfbe,iplx


	return






	if (TIME.ge.t_startpp) then
!solve the plasma current	

	endif

	end
