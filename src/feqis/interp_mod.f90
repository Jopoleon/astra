module interp_mod

implicit none

contains

!---------------------------------------------------------------------
    subroutine qinterp(x1, y1, Nx1, x2, y2, Nx2)

    integer, intent(in) :: Nx1, Nx2
    double precision, intent(in) , dimension(Nx1) :: x1, y1
    double precision, intent(in) , dimension(Nx2) :: x2
    double precision, intent(out), dimension(Nx2) :: y2

    integer :: i, j
    double precision :: A, B, C, z1, z2, z3, t1, t2, t3, t4

    do i=1, Nx2
        t4 = x2(i)

        do j=2, Nx1-1
            z1 = x1(j-1)
            z2 = x1(j)
            z3 = x1(j+1)

            if (t4 == z1) then
                y2(i) = y1(j-1)
                EXIT
            elseif (t4 == z2) then
                y2(i) = y1(j)
                EXIT
            elseif (t4 == z3) then
                y2(i) = y1(j+1)
                EXIT
            elseif (t4 > z1 .and. t4 < z3) then
                t1 = y1(j-1)
                t2 = y1(j)
                t3 = y1(j+1)
                A = (t3 - t2 - (z3 - z2)*(t1 - t2)/(z1 - z2)) / ((z3 - z2)*(z3 - z1))
                B = (t1 - t2)/(z1 - z2) - A*(z1 + z2)
                C = t2 - A*z2**2 - B*z2
                y2(i) = A*t4**2 + B*t4 + C
                EXIT
            elseif (t4 < z1 .and. j == 2) then
                z1 = x1(j-1)**2
                z2 = x1(j)**2
                t1 = y1(j-1)
                t2 = y1(j)
                A = 0.0
                B = (t1 - t2)/(z1 - z2) - A*(z1 + z2)
                C = t2 - A * z2**2 - B*z2
                y2(i) = A * t4**4 + B * t4**2 + C
                EXIT
            elseif (t4 > z3 .and. j == Nx1-1) then
                t1 = y1(j-1)
                t2 = y1(j)
                t3 = y1(j+1)
                A = (t3 - t2 - (z3 - z2)*(t1 - t2)/(z1 - z2)) / ((z3 - z2)*(z3 - z1))
                B = (t1 - t2)/(z1 - z2) - A*(z1 + z2)
                C = t2 - A * z2**2 - B*z2
                y2(i) = A * t4**2 + B*t4 + C
                EXIT
            endif
        enddo
    enddo

    return
    end subroutine qinterp

!------------------------------------------------------------
    subroutine linterp(x1, y1, Nx1, x2, y2, Nx2)

    integer, intent(in) :: Nx1, Nx2
    double precision, intent(in) , dimension(Nx1) :: x1, y1
    double precision, intent(in) , dimension(Nx2) :: x2
    double precision, intent(out), dimension(Nx2) :: y2

    integer :: i, j
    double precision :: A, C, z1, z2, t1, t2, t4

    do i=1, Nx2
        t4 = x2(i)
        do j=2, Nx1
            z1 = x1(j-1)
            z2 = x1(j)
            if (t4 == z1) then
                y2(i)=y1(j-1)
                EXIT
            elseif (t4 == z2) then
                y2(i)=y1(j)
                EXIT
            elseif (t4 > z1 .and. t4 < z2) then
                t1 = y1(j-1)
                t2 = y1(j)
                A = (t2 - t1)/(z2 - z1)
                C = t2 - A*z2
                y2(i) = A*t4 + C
                EXIT
            elseif (t4 < z1 .and.j == 2) then
                y2(i)=y1(j-1)
                EXIT
            elseif (t4 > z2 .and. j == nx1) then
                y2(i) = y1(j)
                EXIT
            endif
        enddo
    enddo

    return
    end subroutine linterp

!------------------------------------------------------------
    subroutine polyfitcc(x, y, P)

    double precision, intent(in) , dimension(3) :: x, y
    double precision, intent(out), dimension(3) :: P

    double precision :: y21, y32, x21, x32, h21, h32

    y32 = y(3) - y(2)
    y21 = y(2) - y(1)
    x32 = x(3) - x(2)
    x21 = x(2) - x(1)
    h32 = x(3) + x(2)
    h21 = x(2) + x(1)
    
    P(1) = (x21*y32 - x32*y21)/(x21*x32*(h32 - h21))
    P(2) = y21/x21 - P(1)*h21
    P(3) = y(3) - P(1)*x(3)**2 - P(2)*x(3)

    return
    end subroutine polyfitcc

!-------------------------------------------------------------------
    function EXTRAPOLATE(x_interp, j1, j2, j3, narr, rho, y_in) result(f_out)
    
    integer, intent(in) :: j1, j2, j3, narr
    double precision, intent(in) :: x_interp
    double precision, intent(in), dimension(narr) :: rho, y_in
    double precision :: f_out
    double precision :: x1, x2, x3, f1, f2, f3, dfdx, d2fdx

    f1 = y_in(j1)
    f2 = y_in(j2)
    f3 = y_in(j3)
    x1 = rho(j1)
    x2 = rho(j2)
    x3 = rho(j3)
    dfdx = (f3 - f1)/(x3 - x1)
    d2fdx = ((f3 - f2)/(x3 - x2) - (f2 - f1)/(x2 - x1))/(x3 - x1)

    f_out = f2 + (x_interp - x2)*(dfdx + d2fdx*(x_interp - x2))
      
    return
    end function EXTRAPOLATE

end module interp_mod
