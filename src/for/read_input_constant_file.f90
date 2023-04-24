subroutine read_input_constant_file(time_ext, dt_smlk)

use const_inc, only: IPSMK

implicit none

real*8 :: time_ext, dt_smlk
	
if (nint(IPSMK) == 1) call rrsudg(time_ext, dt_smlk)

return
end subroutine read_input_constant_file
