subroutine pid_control(Kp_in, Ki_in, Kd_in, error_in, control_signal_in, control_signal_out)
  
use const_inc, only: NA1, TIME, TAU
use standard_functions, only: TIMDER, TIMINT
implicit none

double precision, intent(in):: Kp_in, Ki_in, Kd_in, error_in, control_signal_in
double precision, intent(out) :: control_signal_out

double precision :: error_integral, error_derivative

error_integral = TIMINT(error_in)
error_derivative = TIMDER(error_in)
control_signal_out = control_signal_in + Kp_in*error_in + Kd_in*error_derivative + Ki_in*error_integral
  
return  
end subroutine pid_control
