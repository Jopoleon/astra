subroutine write_output_diag_file

use const_inc, only: IPSMK

implicit none
	
if (nint(IPSMK) == 1) call wrsudg

return
end subroutine write_output_diag_file
