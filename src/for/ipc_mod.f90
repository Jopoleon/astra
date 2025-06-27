module ipc_mod

use parameter_inc, only: NRD

implicit none

!double precision, allocatable, dimension(:, :) :: mem_tglf, mem_qlkz, mem_neo
double precision, dimension(NRD, 25) :: mem_tglf, mem_qlkz, mem_neo

end module ipc_mod
