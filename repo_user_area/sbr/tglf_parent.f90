subroutine tglf_parent

use mpi

use parameter_inc, only: NRD
use const_inc, only: BTOR, RTOR, GP, GP2, ABC, &
    AMJ, AIM1, AIM2, AIM3, ZMJ, PSIAX, PSIBO, &
    NA1, NA1N, NA1E, NA1I
use status_inc, only: NE, TE, NI, TI, &
    ZEF, ZIM1, ZIM2, ZIM3, PBLON, PBPER, &
    PFAST, NIZ3, AMAIN, ER, MU, FP_NORM, &
    RHO, AMETR, SHIF, ELON, &
    NDEUT, NIZ1, NTRIT, NIZ2, NHE3, &
    TRIA, VTOR, NIBM, G11, VPOL, VRS, SHEAR
use parameters_a2equil, only: equil_now

implicit none

integer, parameter :: n_inputs=25, n_outputs=15
double precision, parameter :: &
   k0   = 1.6022E-12, &       ! erg/ev
   e0   = 4.8032E-10, &       ! elementary charge (statcoulombs)
   e00  = 1.6020e-19, &       ! elementary charge (C)
   c0   = 2.9979E+10, &       ! speed of light (cm/sec)
   mp   = 1.6726E-24, &       ! proton mass (g)
   mpp  = 1.6726E-27, &       ! proton mass (kg)
   pi   = 3.141592653589793

integer :: nrho_tg, ierr, intercomm, errcodes(100), status(MPI_STATUS_SIZE)
integer :: i, i1, i2, chunk, nprocs, nworkers, dims(3)
double precision, allocatable :: profs_tg(:,:), out_tg(:,:)
character(len=256) :: worker_exe

worker_exe = "/shares/departments/AUG/users/git/a8/xpr/tglf.x"

call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

nrho_tg = 61
nworkers = 10  ! SET THIS CLEARLY
chunk = nrho_tg / nworkers + 1
nrho_tg = chunk * nworkers  ! Ensure exact division

call MPI_Comm_spawn(worker_exe, MPI_ARGV_NULL, nworkers, MPI_INFO_NULL, 0, MPI_COMM_SELF, intercomm, errcodes, ierr)
print *, "MPI workers = ", nworkers, nrho_tg

allocate(profs_tg(nrho_tg, n_inputs), out_tg(nrho_tg, n_outputs))
profs_tg = 0.5  ! Dummy data

dims(1) = nrho_tg
dims(2) = n_inputs
dims(3) = n_outputs

! Send dimensions and data to workers
do i=0, nworkers-1
    call MPI_Send(dims, 3, MPI_INTEGER, i, 0, intercomm, ierr)
enddo
do i=0, nworkers-1
    i1 = i * chunk + 1
    i2 = (i + 1) * chunk
    call MPI_Send(profs_tg(i1:i2, :), chunk * n_inputs, MPI_DOUBLE_PRECISION, i, 0, intercomm, ierr)
enddo

! Receive results from each worker
do i=0, nworkers-1
    i1 = i * chunk + 1
    i2 = (i + 1) * chunk
    call MPI_Recv(out_tg(i1:i2, :), chunk * n_outputs, MPI_DOUBLE_PRECISION, i, 1, intercomm, status, ierr)
enddo

write(*, *) 'Out', out_tg(:, 1)

return
end subroutine tglf_parent
