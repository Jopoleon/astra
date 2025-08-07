module ipc_mod

implicit none

integer, parameter :: n_sbp_arr_in=33, n_sbp_arr_out=15
double precision, allocatable, dimension(:, :) :: mem_tglf, mem_qlkz, mem_neo

end module ipc_mod
