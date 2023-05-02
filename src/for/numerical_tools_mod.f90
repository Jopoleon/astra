module numerical_tools

implicit none

contains

!---------------------------------------------------------------------
    subroutine reinterp_back(x_in, y_in, Nx_in, x_out, y_out, Nx_out, interp_switch)

    integer, intent(in) :: Nx_in, Nx_out, interp_switch
    double precision, intent(in), dimension(Nx_in) :: x_in, y_in
    double precision, intent(in) :: x_out(Nx_out)
    double precision, intent(out) :: y_out(Nx_out)

    integer :: i
    double precision :: zspl
    double precision, dimension(Nx_in) :: b, c, d

    if (interp_switch == 1) then
        call qinterp(x_in, y_in, Nx_in, x_out, y_out, Nx_out)
    endif

    if (interp_switch == 2) then
        call spline_numerict(x_in, y_in, b, c, d, Nx_in)
        do i=2, Nx_out-1
            zspl = x_out(i)
            y_out(i) = ispline_nt(zspl, x_in, y_in, b, c, d, Nx_in)
        enddo

        if (x_out(Nx_out) > x_in(Nx_in)) then
            call EXTRAP(x_out(1: Nx_out-1), y_out(1: Nx_out-1), x_out(Nx_out), &
                Nx_out-1, y_out(Nx_out), 2, Nx_out-1)
        else
            zspl = x_out(Nx_out)
            y_out(Nx_out) = ispline_nt(zspl, x_in, y_in, b, c, d, Nx_in)
        endif

        if (x_out(1) < x_in(1)) then
            call EXTRAP(x_out(2: Nx_out), y_out(2: Nx_out), x_out(1), 1, y_out(1), 2, Nx_out-1)
        else
            zspl = x_out(1)
            y_out(1) = ispline_nt(zspl, x_in, y_in, b, c, d, Nx_in)
        endif
    endif

    return
    end subroutine reinterp_back

!---------------------------------------------------------------------
    subroutine reinterp_back_quad(x_in, y_in, Nx_in, x_out, y_out, Nx_out)

    integer, intent(in) :: Nx_in, Nx_out
    double precision, intent(in), dimension(Nx_in) :: x_in, y_in
    double precision, intent(in) :: x_out(Nx_out)
    double precision, intent(out) :: y_out(Nx_out)

    integer :: i, j
    double precision :: A, B, C, D, dy1, dy2, dy0, &
        x1, x2, x3, x4, xloc_out, y1, y2, y3, y4, s1, s2

    do i=1, Nx_out
        xloc_out = x_out(i)
        do j=2, Nx_in
            x1 = x_in(j-2)
            x2 = x_in(j-1)
            x3 = x_in(j) 
            x4 = x_in(j+1)
            y1 = y_in(j-2)
            y2 = y_in(j-1)
            y3 = y_in(j)
            y4 = y_in(j+1)
            if (xloc_out == x2) then
                y_out(i) = y_in(j-1)
                EXIT
            else if (xloc_out == x3) then
                y_out(i) = y_in(j)
                EXIT
            else if (xloc_out > x2 .and. xloc_out < x3) then
                if (j == 2) then
                    dy2 = 0.5*((y4 - y3)/(x4 - x3) + (y3 - y2)/(x3 - x2))
                    s2 = x3
                    dy1 = (y3 - y2)/(x3 - x2)
                    s1 = 0.5*(x2 + x3)
                    dy0 = x2*(dy2 - dy1)/(s2 - s1) +  (dy1 - s1*((dy2 - dy1)/(s2 - s1)))
                    dy1 = max(0., dy0)
                else if (j == Nx_in) then
                    dy1 = 0.5*((y2 - y1)/(x2 - x1) + (y3 - y2)/(x3 - x2))
                    dy2 = (y3 - y2)/(x3 - x2)
                else
                    dy1 = 0.5*((y2 - y1)/(x2 - x1) + (y3 - y2)/(x3 - x2))
                    dy2 = 0.5*((y4 - y3)/(x4 - x3) + (y3 - y2)/(x3 - x2))
                endif
                A = (dy1*x2 - dy1*x3  +  dy2*x2 - dy2*x3 - 2*y2  +  2*y3)/ &
                     (x2**3 - 3*x2**2*x3  +  3*x2*x3**2 - x3**3)
                B = -(dy1*x2**2  +  dy1*x2*x3 - 2*dy1*x3**2  + &
                     2*dy2*x2**2 - dy2*x2*x3 - dy2*x3**2 - 3*x2*y2  + &
                     3*x2*y3 - 3*x3*y2  +  3*x3*y3)/ &
                     ((x2 - x3)*(x2**2 - 2*x2*x3  +  x3**2))
                C = (2*dy1*x2**2*x3 - dy1*x2*x3**2 - dy1*x3**3  + &
                     dy2*x2**3  +  dy2*x2**2*x3 - 2*dy2*x2*x3**2 - &
                     6*x2*x3*y2  +  6*x2*x3*y3)/ &
                     ((x2 - x3)*(x2**2 - 2*x2*x3  +  x3**2))
                D = -(dy1*x2**2*x3**2 - dy1*x2*x3**3  +  dy2*x2**3*x3 - &
                      dy2*x2**2*x3**2 - x2**3*y3  +  3*x2**2*x3*y3 - &
                      3*x2*x3**2*y2  +  x3**3*y2)/ &
                      ((x2 - x3)*(x2**2 - 2*x2*x3  +  x3**2))
                y_out(i) = A*xloc_out**3 + B*xloc_out**2 + C*xloc_out + D
                EXIT
            endif
        enddo
    enddo

    if (x_out(Nx_out) > x_in(Nx_in)) then
        call EXTRAP(x_out(1: Nx_out-1), y_out(1: Nx_out-1), x_out(Nx_out), &
            Nx_out-1, y_out(Nx_out), 2, Nx_out-1)
    endif

    return
    end subroutine reinterp_back_quad

!---------------------------------------------------------------------
    subroutine spline_numerict(x, y, b, c, d, n)

!  Calculate the coefficients b(i), c(i), and d(i), i = 1, 2, ..., n
!  for cubic spline interpolation
!  s(x) = y(i) + b(i)*(x - x(i)) + c(i)*(x - x(i))**2 + d(i)*(x - x(i))**3
!  for  x(i) <=  x <=  x(i+1)
!  Alex G: January 2010
!---------------------------------------------------------------------
!  input..
!  x = the arrays of data abscissas (in strictly increasing order)
!  y = the arrays of data ordinates
!  n = size of the arrays xi() and yi() (n>= 2)
!  output..
!  b, c, d  = arrays of spline coefficients
!  comments ...
!  spline.f90 program is based on fortran version of program spline.f
!  the accompanying function fspline can be used for interpolation

    integer, intent(in) :: n
    double precision, intent(in) , dimension(n) :: x, y
    double precision, intent(out), dimension(n) :: b, c, d

    integer :: i, j, gap
    double precision :: h

    gap = n - 1
! check input
    if (n < 2) return
    if (n < 3) then
        b(1) = (y(2) - y(1))/(x(2) - x(1))    ! linear interpolation
        c(1) = 0.
        d(1) = 0.
        b(2) = b(1)
        c(2) = 0.
        d(2) = 0.
        return
    endif

! step 1: preparation

    d(1) = x(2) - x(1)
    c(2) = (y(2) - y(1))/d(1)
    do i=2, gap
        d(i) = x(i+1) - x(i)
        b(i) = 2.0*(d(i-1) + d(i))
        c(i+1) = (y(i+1) - y(i))/d(i)
        c(i) = c(i+1) - c(i)
    enddo

! step 2: end conditions 
    
    b(1) = -d(1)
    b(n) = -d(n-1)
    c(1) = 0.0
    c(n) = 0.0
    if(n /= 3) then
        c(1) = c(3)/(x(4) - x(2)) - c(2)/(x(3) - x(1))
        c(n) = c(n-1)/(x(n) - x(n-2)) - c(n-2)/(x(n-1) - x(n-3))
        c(1) = c(1)*d(1)**2/(x(4) - x(1))
        c(n) = -c(n)*d(n-1)**2/(x(n) - x(n-3))
    endif

! step 3: forward elimination 

    do i=2, n
        h = d(i-1)/b(i-1)
        b(i) = b(i) - h*d(i-1)
        c(i) = c(i) - h*c(i-1)
    enddo

! step 4: back substitution

    c(n) = c(n)/b(n)
    do j=1, gap
        i = n - j
        c(i) = (c(i) - d(i)*c(i+1))/b(i)
    enddo

! step 5: compute spline coefficients

    b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2.0*c(n))
    do i=1, gap
        b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0*c(i))
        d(i) = (c(i+1) - c(i))/d(i)
        c(i) = 3.*c(i)
    enddo
    c(n) = 3.0*c(n)
    d(n) = d(n-1)

    return
    end subroutine spline_numerict

!---------------------------------------------------------------------
    function ispline_nt(u, x, y, b, c, d, n) result(f_out)

! function ispline evaluates the cubic spline interpolation at point z
! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
! where  x(i) <=  u <=  x(i+1)
!---------------------------------------------------------------------
! input..
! u         = the abscissa at which the spline is to be evaluated
! x, y     = the arrays of given data points
! b, c, d = arrays of spline coefficients computed by spline
! n         = the number of data points
! output:
! ispline = interpolated value at point u

    integer, intent(in) :: n
    double precision, intent(in) :: u
    double precision, intent(in), dimension(n) :: x, y, b, c, d
    double precision :: f_out

    integer :: i, j, k
    double precision :: dx

! if u is ouside the x() interval take a boundary value (left or right)
    if (u <= x(1)) then
        f_out = y(1)
        return
    else if (u >= x(n)) then
        f_out = y(n)
        return
    endif

!  binary search for for i, such that x(i) <=  u <=  x(i+1)

    i = 1
    j = n + 1
    do while (j > i+1)
        k = (i + j)/2
        if(u < x(k)) then
            j = k
        else
            i = k 
        endif
    enddo

! evaluate spline interpolation

    dx = u - x(i)
    f_out = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))

    return
    end function ispline_nt

!---------------------------------------------------------------------
    subroutine DERIV_CDE(x_in, x_type, y_in, yd_out, nagrid)

! DERIV computes first and second derivative over x
!
!  x_in: x_variable
!  x_type for GRP style-grid is: 1 main grid -> shifted grid (deriv), 2 shifted grid -> main grid (deriv)
!  y_in: y_variable
!  yd_out: derivative
!  order_d: 1st or 2nd derivative
!  nagrid: number of grid points
! iextrap: 1 if yes interpolate last grid point, 0 do not interpolate last grid point

    integer, intent(in) :: x_type, nagrid
    double precision, intent(in), dimension(nagrid) :: x_in, y_in
    double precision, intent(out), dimension(nagrid) :: yd_out

    integer :: j, iextrap, j_end, OEXTRAP
    double precision :: x1tmp, y1tmp, y2tmp, drho

    OEXTRAP = 1
    iextrap = 1

    if (iextrap == 1) j_end = 1
    if (iextrap == 0) j_end = 0

! Normalized grid , GRP style
    j = 1
    if (x_type == 1) then
        yd_out(j) = 0.0
    endif
    if (x_type == 2) then
        x1tmp = x_in(2) - x_in(1)
        y1tmp = y_in(2) - y_in(1)
        yd_out(j) = y1tmp/x1tmp
    endif

    if (x_type == 3) then
        x1tmp = x_in(2) - x_in(1)
        y1tmp = y_in(2) - y_in(1)
        yd_out(j) = y1tmp/x1tmp
    endif

    do j=2, nagrid-j_end
        drho  = (x_in(j+1) - x_in(j-1))
        y1tmp = (y_in(j+1) - y_in(j-1))
        y2tmp = drho
        yd_out(j) = y1tmp/y2tmp
    enddo

    ! Interpolate to last grid point
    j = nagrid
    drho  = (x_in(j) - x_in(j-1))
    y1tmp = (y_in(j) - y_in(j-1))
    y2tmp = drho
    yd_out(j) = y1tmp/y2tmp

! Interpolate to last grid point
    j = nagrid
    call EXTRAP(x_in(1: j-1), yd_out(1: j-1), x_in(j), & 
        j-1, yd_out(j), OEXTRAP, j-1)    

    return
    end subroutine DERIV_CDE

!---------------------------------------------------------------------
    subroutine linterp(x_in, y_in, Nx_in, x_out, y_out, Nx_out)

! Linear inerpolation

    integer, intent(in) :: Nx_in, Nx_out
    double precision, intent(in), dimension(Nx_in) :: x_in, y_in
    double precision, intent(in)  :: x_out(Nx_out)
    double precision, intent(out) :: y_out(Nx_out)

    integer :: i, j
    double precision :: A, B, x1, x2, y1, y2, xloc_out

    do i=1, Nx_out
        xloc_out = x_out(i)
        do j=2, Nx_in
            x1 = x_in(j-1)
            x2 = x_in(j)
            y1 = y_in(j-1)
            y2 = y_in(j)
            if (xloc_out == x1) then
                y_out(i) = y1
                EXIT
            else if (xloc_out == x2) then
                y_out(i) = y2
                EXIT
            else if (xloc_out < x1 .and. j == 2) then
                x1 = x_in(j-1)**2
                x2 = x_in(j)**2
                A = (y2 - y1)/(x2 - x1)
                B = y1 - x1*A
                y_out(i) = A * xloc_out**2 + B
                EXIT
            else if ( (xloc_out > x1 .and. xloc_out < x2) .or. &
                      (xloc_out > x2 .and. j == Nx_in) ) then
                A = (y2 - y1)/(x2 - x1)
                B = y1 - x1*A
                y_out(i) = A*xloc_out + B
                EXIT
            endif
        enddo
    enddo

    return
    end subroutine linterp

!---------------------------------------------------------------------
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
                y_out(i) = y_in(j-1)
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
                      (xloc_out > x3 .and. j == Nx_out-1) ) then
                A = (y3 - y2 - (x3 - x2)*(y1 - y2)/(x1 - x2)) / ((x3 - x2)*(x3 - x1))
                B = (y1 - y2)/(x1 - x2) - A*(x1 + x2)
                C = y2 - A*(x2**2) - B*x2
                y_out(i) = A * xloc_out**2 + B*xloc_out + C
                EXIT
            endif
        enddo
    enddo

    return
    end subroutine qinterp

!---------------------------------------------------------------------
    subroutine integrcc(nx, x, y, sy)

    integer, intent(in) :: nx
    double precision, intent(in), dimension(nx) :: x, y
    double precision, intent(out) :: sy(nx)

    integer :: i
    double precision :: drho, y1tmp

    sy = 0.
    do i=2, nx
        drho = (x(i) - x(i-1))
        y1tmp = (y(i) + y(i-1))/2.
        sy(i) = sy(i-1) + y1tmp*drho
    enddo

    return
    end subroutine integrcc

!---------------------------------------------------------------------
    subroutine derivcc(nx, x, y, dy, gga)

    integer, intent(in) :: nx, gga
    double precision, intent(in), dimension(nx) :: x, y
    double precision, intent(out) :: dy(nx)

    integer :: i
    double precision :: P(3)

    dy = 0.
    do i=2, nx-1
        dy(i) = (y(i+1) - y(i-1))/(x(i+1) - x(i-1))
    enddo
    if (gga == 1) dy(1) = 0.
    if (gga == 2) dy(1) = (y(2) - y(1))/(x(2) - x(1))

    call polyfitcc(x(nx-2: nx), y(nx-2: nx), P)

    dy(nx) = 2.*P(1)*x(nx) + P(2)

    return
    end subroutine derivcc

!---------------------------------------------------------------------
    subroutine polyfitcc(x, y, P)

    double precision, intent(in), dimension(3) :: x, y
    double precision, intent(out) :: P(3)
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

!---------------------------------------------------------------------
    subroutine EXTRAP(x_in, y_in, x_extrap, j_extrap, y_extrap, ex_order, nagrid)

! Assume that x is of r-type, i.e. interpolation in 0 has zero odd derivatives

    integer, intent(in) :: ex_order, j_extrap, nagrid
    double precision, intent(in)  :: x_in(nagrid), x_extrap, y_in(nagrid)
    double precision, intent(out) :: y_extrap

    integer :: k1, k2, k3, jsign
    double precision :: P(3)
            
    if (j_extrap == nagrid) jsign = -1
    if (j_extrap == 1) jsign = 1
      
! Constant interpolation
    if (ex_order == 0) then
        k1 = j_extrap
        P(1) = y_in(k1)
        y_extrap = P(1)
    endif

    if (ex_order == 1) then ! Linear extrapolation
        if (jsign < 0) then
            k1 = j_extrap + jsign
            k2 = j_extrap
            call polyfitcc_1(x_in(k1: k2), y_in(k1: k2), P(1: 2))
            y_extrap = P(1)*x_extrap + P(2)
        endif
        if (jsign > 0) then
            k1 = j_extrap
            k2 = j_extrap + 1
            call polyfitcc_1(x_in(k1: k2), y_in(k1: k2), P(1: 2))
            y_extrap = P(1)*x_extrap + P(2)
        endif
    else if (ex_order == 2) then! Quadratic extrapolation
        if (jsign < 0) then
            k1 = j_extrap + jsign*2
            k2 = j_extrap + jsign
            k3 = j_extrap
            call polyfitcc(x_in(k1: k3), y_in(k1: k3), P)
            y_extrap = P(1) * x_extrap**2 + P(2)*x_extrap + P(3)
        else if (jsign > 0) then
            k1 = j_extrap
            k2 = j_extrap + jsign
            k3 = j_extrap + jsign*2
            call polyfitcc(x_in(k1: k3), y_in(k1: k3), P)
            y_extrap = P(1) * x_extrap**2 + P(2)*x_extrap + P(3)
        endif
    endif

    return
    end subroutine EXTRAP

!---------------------------------------------------------------------
    subroutine polyfitcc_1(x, y, P)

    double precision, intent(in) , dimension(2) :: x, y
    double precision, intent(out), dimension(2) :: P

    double precision :: y21, x21

    y21 = y(2) - y(1)      
    x21 = x(2) - x(1)      

    P(1) = y21/x21
    P(2) = y(1) - x(1)*y21/x21

    return
    end subroutine polyfitcc_1

!---------------------------------------------------------------------
    subroutine DERIV(x_in, x_out, x_type, y_in, yd_out, order_d, nagrid, iextrap)

! DERIV computes first or second derivative over x
!
!  x_in: x_variable
!  x_type for GRP style-grid is: 1 main grid -> shifted grid (deriv), 2 shifted grid -> main grid (deriv)
!  y_in: y_variable
!  yd_out: derivative
!  order_d: 1st or 2nd derivative
!  nagrid: number of grid points
! iextrap: 1 if yes interpolate last grid point, 0 do not interpolate last grid point

    integer, intent(in) :: order_d, x_type, nagrid, iextrap
    double precision, intent(in) , dimension(nagrid) :: x_in, x_out, y_in
    double precision, intent(out), dimension(nagrid) :: yd_out

    integer :: j, j_end, OEXTRAP
    double precision dx, dy, y0, P(3)
     
    OEXTRAP = 1

    if (iextrap == 1) j_end=1
    if (iextrap == 0) j_end=0
  
! dy/dx
    if (order_d == 1) then

        if (x_type == 1) then
            do j=1, nagrid - j_end
                dx = x_in(j+1) - x_in(j)
                dy = y_in(j+1) - y_in(j)
                yd_out(j) = dy/dx
            enddo
        endif

        if (x_type == 2) then
            j = 1
            call polyfitcc(x_in(1: 3), y_in(1: 3), P)
            y0 = P(3)
            dx = x_in(1)
            dy = (y_in(1) - y0)
            yd_out(j) = dy/dx
            do j=2, nagrid
                dx = x_in(j) - x_in(j-1)
                dy = y_in(j) - y_in(j-1)
                yd_out(j) = dy/dx
            enddo
        endif

    else if (order_d == 2) then ! d2y/dx^2
        if (x_type == 1) then
            j = 1
            dx = (x_in(j+1) - x_in(j))**2
            dy = (y_in(j+1) - y_in(j))
            yd_out(j) = dy/dx
    
            do j=2, nagrid - j_end
                dx = (x_in(j+1) - x_in(j))**2
                dy = (y_in(j+1) - 2.*y_in(j) + y_in(j-1))
                yd_out(j) = dy/dx
            enddo
        endif
        if (x_type == 2) then
! First interpolate shifted variable to zero
            j = 1
            call polyfitcc(x_in(1:3), y_in(1:3), P)
            y0 = P(3)
            dx = x_in(j)**2
            dy = (y_in(j+1) - 2.0*y_in(j) + y0)
            yd_out(j) = dy/dx
            do j=2, nagrid - j_end
                dx = (x_in(j+1) - x_in(j))**2
                dy = (y_in(j+1) - 2.0*y_in(j) + y_in(j-1))
                yd_out(j) = dy/dx
            enddo
        endif
    endif

! Interpolate to last grid point
    if (iextrap == 1 .and. order_d == 1) then
        j = nagrid
        call polyfitcc(x_in(j-2: j), y_in(j-2: j), P)
        yd_out(j) = 2.*P(1)*x_out(j) + P(2)
    endif

    if (iextrap == 1 .and. order_d == 2) then
        j = nagrid
        call polyfitcc(x_in(j-2: j), y_in(j-2: j), P)
        yd_out(j) = 2.*P(1)
    endif

    return
    end subroutine DERIV

!---------------------------------------------------------------------
    subroutine INTEGR(x_in, x_type, y_in, ys_out, nagrid)

! INTEGR computes integrals over x
!  x_in: x_variable
!  x_type for GRP style-grid is: 1 main grid -> shifted grid (deriv), 2 shifted grid -> main grid (deriv)
!  y_in: y_variable
!  ys_out: integral
!  nagrid: number of grid points

    integer, intent(in) :: x_type, nagrid
    double precision, intent(in) , dimension(nagrid) :: x_in, y_in
    double precision, intent(out), dimension(nagrid) :: ys_out

    integer :: j
    double precision :: P(3)

! Normalized grid , GRP style
    if (x_type == 1) then
        ys_out(1) = y_in(1)*x_in(1)
        do j=2, nagrid
            ys_out(j) = ys_out(j-1) + y_in(j)*(x_in(j) - x_in(j-1))
        enddo
    else if (x_type == 2) then ! First interpolate shifted variable to zero
        call polyfitcc(x_in(1: 3), y_in(1: 3), P)
        ys_out(1) = P(3)*x_in(1)
        do j=2, nagrid
            ys_out(j) = ys_out(j-1) + y_in(j-1)*(x_in(j) - x_in(j-1))
        enddo
    endif

    return
    end subroutine INTEGR

!---------------------------------------------------------------------
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


end module numerical_tools
