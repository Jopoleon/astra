	subroutine write_output_diag_file

	use parameter_inc
	use const_inc

	implicit none
	
		if (nint(IPSMK).eq.1) call wrsudg

	end subroutine write_output_diag_file
