subroutine torbeam_rabbit

use status_inc, only: PEECR, CUECR

implicit none
! If using ctr2rz, one needs time consistency between ctr2rz and TORBEAM, RABBIT

call ctr2rz
call A2RABBIT
CALL A2TORBEAM

return
end subroutine torbeam_rabbit
