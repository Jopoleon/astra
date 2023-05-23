subroutine read_input_constant_file(time_ext, dt_smlk)

use debugger, only: flightsim

implicit none

real*8, intent(out) :: time_ext, dt_smlk

if (flightsim == 1) call rrsudg(time_ext, dt_smlk)

return
end subroutine read_input_constant_file
