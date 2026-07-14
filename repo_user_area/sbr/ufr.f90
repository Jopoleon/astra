!------------------------------------------------------------------
! UFNAME  - U-file name to read
! TIME    - prescrbed time for output
! arr1d   - Output array at given time
!------------------------------------------------------------------

subroutine UF1DR(ufnam, tim_in, val)

use read_input, only: ufheader, ufrd

implicit none

double precision, intent(in) :: tim_in
character(len=120), intent(in) :: ufnam
double precision, intent(out) :: val

logical :: uf_exists
integer :: nt_u, nx_u, nscal_u, ndim_u, j, jt, jch, jtprev
double precision :: tim1, tim2, tmp1, tmp2
double precision, dimension(:), allocatable :: t_u, x_u, var_u
character(len=30) :: rholbl

INQUIRE(file=TRIM(ufnam), exist=uf_exists)
if (.not. uf_exists) then
    write(*, '(3A)') 'File ' // TRIM(ufnam) // ' not found, returning 0'
    val = 0.
    return
endif

call ufheader(TRIM(ufnam), nscal_u, ndim_u, nt_u, nx_u, rholbl)
call ufrd(TRIM(ufnam), nscal_u, ndim_u, nt_u, nx_u, t_u, x_u, var_u)

! Interpolate between time points
if (tim_in <= t_u(1)) then
    val = var_u(1)
else if (tim_in >= t_u(nt_u)) then
    val = var_u(nt_u)
else
    jt = 1
    do while(t_u(jt) < tim_in)
        jt = jt + 1
    enddo
    jtprev = jt - 1
    tim1 = t_u(jtprev)
    tim2 = t_u(jtprev + 1)
    tmp1 = var_u(jtprev)
    tmp2 = var_u(jtprev + 1)
    val = tmp1 + (tim_in - tim1) * (tmp2 - tmp1)/(tim2 - tim1)
endif

return
end subroutine UF1DR

!------------------------------------------------------------------
! UFNAME  - U-file name to read
! TIME    - prescrbed time for output
! arr1d   - Output array at given time
!------------------------------------------------------------------

subroutine UF2DR(ufnam, tim_in, arr1d)

use read_input, only: ufheader, ufrd

implicit none

double precision, intent(in) :: tim_in
character(len=120), intent(in) :: ufnam
double precision, intent(out) :: arr1d(*)

logical :: uf_exists
integer :: nt_u, nx_u, nscal_u, ndim_u, j, jt, jch, jtprev
double precision :: tim1, tim2
double precision, dimension(:), allocatable :: t_u, x_u, var_u, tmp1, tmp2
character(len=30) :: rholbl

INQUIRE(file=TRIM(ufnam), exist=uf_exists)
if (.not. uf_exists) then
    write(*, '(3A)') 'File ' // TRIM(ufnam) // ' not found, returning 0'
    arr1d(1: 20) = 0.
    return
endif

call ufheader(TRIM(ufnam), nscal_u, ndim_u, nt_u, nx_u, rholbl)
allocate(tmp1(nx_u))
allocate(tmp2(nx_u))
call ufrd(TRIM(ufnam), nscal_u, ndim_u, nt_u, nx_u, t_u, x_u, var_u)

! Interpolate between time points
if (tim_in <= t_u(1)) then
    arr1d(1: nx_u) = (/ (var_u((jch-1)*nt_u + 1 ), jch=1, nx_u) /)
else if (tim_in >= t_u(nt_u)) then
    arr1d(1: nx_u) = (/ (var_u((jch-1)*nt_u + nt_u), jch=1, nx_u) /)
else
    jt = 1
    do while(t_u(jt) < tim_in)
        jt = jt + 1
    enddo
    jtprev = jt - 1
    tim1 = t_u(jtprev)
    tim2 = t_u(jtprev + 1)
    tmp1 = (/ (var_u((jch-1)*nt_u + jtprev)    , jch=1, nx_u) /)
    tmp2 = (/ (var_u((jch-1)*nt_u + jtprev + 1), jch=1, nx_u) /)
    arr1d(1:nx_u) = tmp1(1: nx_u) + (tim_in - tim1) * &
                   (tmp2(1: nx_u) - tmp1(1: nx_u))/(tim2 - tim1)
endif

deallocate(tmp1)
deallocate(tmp2)

end subroutine UF2DR
