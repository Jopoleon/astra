!quadratic inerpolation
!----------------------------------------------------------------------|
subroutine qinterp(x_in, y_in, n_in, x_out, y_out, n_out)
!----------------------------------------------------------------------|
implicit none

integer, intent(in) :: n_in, n_out
integer :: i, j

real, intent(in) , dimension(n_in)  :: x_in, y_in
real, intent(in) , dimension(n_out) :: x_out
real, intent(out), dimension(n_out) :: y_out
real A, B
real z1, z2, z3, z4
real t1, t3, t4

do i=1, n_out

    t1 = x_out(i)

    do j=2, n_in
        z1 = x_in(j-1)
        z2 = x_in(j)
 
        if (t1 .eq. z1) then
            y_out(i) = y_in(j-1)
            EXIT
        endif

        if (t1 .eq. z2) then
            y_out(i) = y_in(j)
            EXIT
        endif

        if (t1 .gt. z1 .and. t1 .lt. z2) then
            t3 = x_in(j-1)
            t4 = x_in(j)
 
            z3 = y_in(j-1)
            z4 = y_in(j)
            A = (z4 - z3)/(t4 - t3)
            B = z3 - t3*A
            y_out(i) = A*t1 + B

            EXIT
        endif

        if (t1 .lt. z1 .and. j .eq. 2) then
            t3 = x_in(j-1)**2.0
            t4 = x_in(j)**2.0
 
            z3 = y_in(j-1)
            z4 = y_in(j)
            A = (z4 - z3)/(t4 - t3)
            B = z3 - t3*A
            y_out(i) = A*t1**2.0 + B

            EXIT
        endif

        if (t1 .gt. z2 .and. j .eq. n_in) then
            t3 = x_in(j-1)
            t4 = x_in(j)
 
            z3 = y_in(j-1)
            z4 = y_in(j)
            A = (z4 - z3)/(t4 - t3)
            B = z3 - t3*A
            y_out(i) = A*t1 + B

            EXIT
        endif

    enddo
enddo

end subroutine qinterp
