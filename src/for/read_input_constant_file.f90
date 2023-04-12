	subroutine read_input_constant_file(time_ext,dt_smlk)

	use parameter_inc
	use status_inc
	use const_inc
	use outcmn_inc
	
	implicit none

	real*8  time_ext,dt_smlk
	
		if (nint(IPSMK).eq.1) call rrsudg(time_ext,dt_smlk)

	end subroutine read_input_constant_file
