subroutine write_output_diag_file

use debugger, only: flightsim

implicit none

if (flightsim == 1) call wrsudg

return
end subroutine write_output_diag_file
