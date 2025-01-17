subroutine qinterp(x_in, y_in, Nx_in, x_out, y_out, Nx_out)

! Quadratic inerpolation

integer, intent(in) :: Nx_in, Nx_out
double precision, intent(in), dimension(Nx_in) :: x_in, y_in
double precision, intent(in)  :: x_out(Nx_out)
double precision, intent(out) :: y_out(Nx_out)

integer :: i, j
double precision :: A, B, C, x1, x2, x3, y1, y2, y3, xloc_out

do i=1, Nx_out
    xloc_out = x_out(i)
    do j=2, Nx_in-1
        x1 = x_in(j-1)
        x2 = x_in(j)
        x3 = x_in(j+1)
        y1 = y_in(j-1)
        y2 = y_in(j)
        y3 = y_in(j+1)
        if (xloc_out == x1) then
            y_out(i) = y1
            EXIT
        else if (xloc_out == x2) then
            y_out(i) = y2
            EXIT
        else if (xloc_out == x3) then
            y_out(i) = y3
            EXIT
        else if (xloc_out < x1 .and. j == 2) then
            x1 = x_in(j-1)**2
            x2 = x_in(j)**2
            B = (y1 - y2)/(x1 - x2)
            C = y2 - B*x2
            y_out(i) = B * xloc_out**2 + C
            EXIT
        else if ( (xloc_out > x1 .and. xloc_out < x3) .or. &
                  (xloc_out > x3 .and. j == Nx_in-1) ) then
            A = (y3 - y2 - (x3 - x2)*(y1 - y2)/(x1 - x2)) / ((x3 - x2)*(x3 - x1))
            B = (y1 - y2)/(x1 - x2) - A*(x1 + x2)
            C = y2 - A* x2**2 - B*x2
            y_out(i) = A * xloc_out**2 + B*xloc_out + C
            EXIT
        endif
    enddo
enddo

return
end subroutine qinterp
