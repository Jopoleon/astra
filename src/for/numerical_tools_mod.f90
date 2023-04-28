module numerical_tools

implicit none

contains

subroutine reinterp_back(x1, y1, Nx1, x2, y2, Nx2, interp_switch)

integer, intent(in) :: Nx1, Nx2, interp_switch
double precision, intent(in), dimension(Nx1) :: x1, y1
double precision, intent(in) :: x2(Nx2)
double precision, intent(out) :: y2(Nx2)

integer :: i
double precision :: zspl
double precision, dimension(Nx1) :: b, c, d

if (interp_switch == 1) then
    call qinterp(x1, y1, Nx1, x2, y2, Nx2)
endif

if (interp_switch == 2) then
    call spline_numerict(x1, y1, b, c, d, Nx1)
    do i=2, Nx2-1
        zspl = x2(i)
        y2(i) = ispline_nt(zspl, x1, y1, b, c, d, Nx1)
    enddo

    if (x2(Nx2) > x1(Nx1)) then
        call EXTRAP(x2(1: Nx2-1), y2(1: Nx2-1), x2(Nx2), &
            Nx2-1, y2(Nx2), 2, Nx2-1)
    else
        zspl = x2(Nx2)
        y2(Nx2) = ispline_nt(zspl, x1, y1, b, c, d, Nx1)
    endif

    if (x2(1) < x1(1)) then
        call EXTRAP(x2(2: Nx2), y2(2: Nx2), x2(1), 1, y2(1), 2, Nx2-1)
    else
        zspl = x2(1)
        y2(1) = ispline_nt(zspl, x1, y1, b, c, d, Nx1)
    endif
endif

end subroutine reinterp_back

!--------------------------------------------------
subroutine reinterp_back_quad(x1, y1, Nx1, x2, y2, Nx2)

integer, intent(in) :: Nx1, Nx2
double precision, intent(in), dimension(Nx1) :: x1, y1
double precision, intent(in) :: x2(Nx2)
double precision, intent(out) :: y2(Nx2)

integer :: i, j

double precision :: A, B, C, D, dy1, dy2, dy0, B1, D1, C1, C2, &
    z1, z2, z3, z4, z0, t1, t2, t3, t4, s1, s2

do i = 1, Nx2
    z0 = x2(i) 
    do j = 2, Nx1
        z2 = x1(j-1)
        z3 = x1(j)
 
        if (z0 == z2) then
            y2(i) = y1(j-1)
            EXIT
        endif

        if (z0 == z3) then
            y2(i) = y1(j)
            EXIT
        endif

        if (z0 > z2 .and. z0 < z3) then
            if (j == 2) then
                z4 = x1(j+1)
                t2 = y1(j-1)
                t3 = y1(j)
                t4 = y1(j+1)
                dy2 = 0.5*((t4 - t3)/(z4 - z3) + (t3 - t2)/(z3 - z2))
                s2 = z3
                dy1 = (t3 - t2)/(z3 - z2)
                s1 = 0.5*(z2 + z3)
                dy0 = z2*(dy2 - dy1)/(s2 - s1) +  (dy1 - s1*((dy2 - dy1)/(s2 - s1)))
                dy1 = max(0., dy0)
            endif
            if (j == Nx1) then
                z1 = x1(j-2)
                t1 = y1(j-2)
                t2 = y1(j-1)
                t3 = y1(j)
                dy1 = 0.5*((t2 - t1)/(z2 - z1) + (t3 - t2)/(z3 - z2))
                dy2 = 1.0*((t3 - t2)/(z3 - z2))
            endif
            if (j > 2.and.j < Nx1) then
                z1 = x1(j-2)
                z4 = x1(j+1)
                t1 = y1(j-2)
                t2 = y1(j-1)
                t3 = y1(j)
                t4 = y1(j + 1)
                dy1 = 0.5*((t2 - t1)/(z2 - z1) + (t3 - t2)/(z3 - z2))
                dy2 = 0.5*((t4 - t3)/(z4 - z3) + (t3 - t2)/(z3 - z2))
            endif

            D1 = (dy2 - dy1)/(2.*(z3 - z2)) 
            B1 = -3./2.*(z3**2. - z2**2.)/(z3 - z2)
            C1 = dy1 - 2.*D1*z2
            C2 = -(3.*z2**2. + 2.*B1*z2)

            A = (z3 - z2 - B1*(z3**2. - z2**2.) - C1*(z3 - z2))/ &
                 (z3**3. + B1*z3**2. - z2**3. - B1*z2**2. + C2*(z3 - z2))
            B = D1 + B1*A
            C = C1 + C2*A
            D = t2 - A*z2**3. - B*z2**2. - C*z2

            A = (dy1*z2 - dy1*z3  +  dy2*z2 - dy2*z3 - 2*t2  +  2*t3)/ &
                 (z2**3 - 3*z2**2*z3  +  3*z2*z3**2 - z3**3)
            B = -(dy1*z2**2  +  dy1*z2*z3 - 2*dy1*z3**2  + &
                 2*dy2*z2**2 - dy2*z2*z3 - dy2*z3**2 - 3*z2*t2  + &
                 3*z2*t3 - 3*z3*t2  +  3*z3*t3)/ &
                 ((z2 - z3)*(z2**2 - 2*z2*z3  +  z3**2))
            C = (2*dy1*z2**2*z3 - dy1*z2*z3**2 - dy1*z3**3  + &
                 dy2*z2**3  +  dy2*z2**2*z3 - 2*dy2*z2*z3**2 - &
                 6*z2*z3*t2  +  6*z2*z3*t3)/ &
                 ((z2 - z3)*(z2**2 - 2*z2*z3  +  z3**2))
            D = -(dy1*z2**2*z3**2 - dy1*z2*z3**3  +  dy2*z2**3*z3 - &
                  dy2*z2**2*z3**2 - z2**3*t3  +  3*z2**2*z3*t3 - &
                  3*z2*z3**2*t2  +  z3**3*t2)/ &
                  ((z2 - z3)*(z2**2 - 2*z2*z3  +  z3**2))

            y2(i) = A*z0**3.0 + B*z0**2. + C*z0 + D
            EXIT
        endif

    enddo
enddo

if (x2(Nx2) > x1(Nx1)) then
    call EXTRAP(x2(1: Nx2-1), y2(1: Nx2-1), x2(Nx2), &
        Nx2-1, y2(Nx2), 2, Nx2-1)
endif

end subroutine reinterp_back_quad

!------------------------------------------------
subroutine spline_numerict (x, y, b, c, d, n)

!  Calculate the coefficients b(i), c(i), and d(i), i = 1, 2, ..., n
!  for cubic spline interpolation
!  s(x) = y(i) + b(i)*(x - x(i)) + c(i)*(x - x(i))**2 + d(i)*(x - x(i))**3
!  for  x(i) <=  x <=  x(i+1)
!  Alex G: January 2010
!----------------------------------------------------------------------
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
double precision, intent(in), dimension(n) :: x, y
double precision, intent(out), dimension(n) :: b, c, d

integer :: i, j, gap
double precision :: h

gap = n - 1
! check input
if ( n < 2 ) return
if ( n < 3 ) then
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
do i = 2, gap
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

do i = 2, n
    h = d(i-1)/b(i-1)
    b(i) = b(i) - h*d(i-1)
    c(i) = c(i) - h*c(i-1)
enddo

! step 4: back substitution

c(n) = c(n)/b(n)
do j = 1, gap
    i = n - j
    c(i) = (c(i) - d(i)*c(i+1))/b(i)
enddo

! step 5: compute spline coefficients

b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2.0*c(n))
do i = 1, gap
    b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0*c(i))
    d(i) = (c(i+1) - c(i))/d(i)
    c(i) = 3.*c(i)
enddo
c(n) = 3.0*c(n)
d(n) = d(n-1)

end subroutine spline_numerict

!-----------------------------------------------
function ispline_nt(u, x, y, b, c, d, n) result(f_out)

! function ispline evaluates the cubic spline interpolation at point z
! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
! where  x(i) <=  u <=  x(i+1)
!----------------------------------------------------------------------
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
if(u <= x(1)) then
    f_out = y(1)
    return
endif
if(u >= x(n)) then
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

!-------------------------------------------------------------
!Efable DERIV computes first and second derivative over x
!
!  x_input: x_variable
!  x_type for GRP style-grid is: 1 main grid -> shifted grid (deriv), 2 shifted grid -> main grid (deriv)
!  y_input: y_variable
!  yd_output: derivative
!  order_d: 1st or 2nd derivative
!  nagrid: number of grid points
! iextrap: 1 if yes interpolate last grid point, 0 do not interpolate last grid point
!-------------------------------------------------------------
subroutine DERIV_CDE(x_input, x_type, y_input, yd_output, nagrid)

integer, intent(in) :: x_type, nagrid
double precision, intent(in), dimension(nagrid) :: x_input, y_input
double precision, intent(out), dimension(nagrid) :: yd_output

integer :: j, iextrap, j_end, OEXTRAP
double precision :: x1tmp, y1tmp, y2tmp, drho

OEXTRAP = 1
iextrap = 1

if (iextrap == 1) j_end = 1
if (iextrap == 0) j_end = 0

! Normalized grid , GRP style
j = 1
if (x_type == 1) then
    yd_output(j) = 0.0
endif
if (x_type == 2) then
    x1tmp = x_input(2) - x_input(1)
    y1tmp = y_input(2) - y_input(1)
    yd_output(j) = y1tmp/x1tmp
endif

if (x_type == 3) then
    x1tmp = x_input(2) - x_input(1)
    y1tmp = y_input(2) - y_input(1)
    yd_output(j) = y1tmp/x1tmp
endif

do j = 2, nagrid-j_end
    drho  = (x_input(j+1) - x_input(j-1))
    y1tmp = (y_input(j+1) - y_input(j-1))
    y2tmp = drho
    yd_output(j) = y1tmp/y2tmp
enddo

! Interpolate to last grid point
j = nagrid
drho  = (x_input(j) - x_input(j-1))
y1tmp = (y_input(j) - y_input(j-1))
y2tmp = drho
yd_output(j) = y1tmp/y2tmp

! Interpolate to last grid point
j = nagrid
call EXTRAP(x_input(1: j-1), yd_output(1: j-1), x_input(j), & 
    j-1, yd_output(j), OEXTRAP, j-1)    

end subroutine DERIV_CDE

!linear inerpolation
!----------------------------------------------------------------------|
subroutine linterp(x1, y1, Nx1, x2, y2, Nx2)

integer, intent(in) :: Nx1, Nx2
double precision, intent(in), dimension(Nx1) :: x1, y1
double precision, intent(in) :: x2(Nx2)
double precision, intent(out) :: y2(Nx2)

integer :: i, j
double precision :: A, B, z1, z2, z3, z4, t1, t3, t4

do i = 1, Nx2
    t1 = x2(i)

    do j = 2, Nx1
        z1 = x1(j-1)
        z2 = x1(j)
        if (t1 == z1) then
            y2(i) = y1(j-1)
            EXIT
        endif
        if (t1 == z2) then
            y2(i) = y1(j)
            EXIT
        endif

        if (t1 > z1 .and. t1 < z2) then
            z3 = y1(j-1)
            z4 = y1(j)
            t3 = x1(j-1)
            t4 = x1(j)
 
            A = (z4 - z3)/(t4 - t3)
            B = z3 - t3*A
            y2(i) = A*t1 + B
            EXIT
        endif

        if (t1 < z1 .and. j == 2) then
            z3 = y1(j-1)
            z4 = y1(j)
            t3 = x1(j-1)**2.0
            t4 = x1(j)**2.0
            A = (z4 - z3)/(t4 - t3)
            B = z3 - t3*A
            y2(i) = A*t1**2.0 + B
            EXIT
        endif

        if (t1 > z2 .and. j == Nx1) then
            z3 = y1(j-1)
            z4 = y1(j)
            t3 = x1(j-1)
            t4 = x1(j)
            A = (z4 - z3)/(t4 - t3)
            B = z3 - t3*A
            y2(i) = A*t1 + B
            EXIT
        endif
    enddo
enddo

end subroutine linterp

!quadratic inerpolation
!----------------------------------------------------------------------|
subroutine qinterp(x1, y1, Nx1, x2, y2, Nx2)

integer, intent(in) :: Nx1, Nx2
double precision, intent(in), dimension(Nx1) :: x1, y1
double precision, intent(in) :: x2(Nx2)
double precision, intent(out) :: y2(Nx2)

integer :: i, j
double precision :: A, B, C, z1, z2, z3, t1, t2, t3, t4

do i = 1, Nx2
    t4 = x2(i)

    do j = 2, Nx1-1
        z1 = x1(j-1)
        z2 = x1(j)
        z3 = x1(j+1)
 
        if (t4 == z1) then
            y2(i) = y1(j-1)
            EXIT
        endif

        if (t4 == z2) then
            y2(i) = y1(j)
            EXIT
        endif

        if (t4 == z3) then
            y2(i) = y1(j+1)
            EXIT
        endif

        if (t4 > z1 .and. t4 < z3) then
            t1 = y1(j-1)
            t2 = y1(j)
            t3 = y1(j+1)
 
            A = (t3 - t2 - (z3 - z2)*(t1 - t2)/(z1 - z2)) / ((z3 - z2)*(z3 - z1))
            B = (t1 - t2)/(z1 - z2) - A*(z1 + z2)
            C = t2 - A*(z2**2.0) - B*z2

            y2(i) = A*(t4**2.0) + B*t4 + C
            EXIT
        endif

        if (t4 < z1 .and. j == 2) then
            z1 = x1(j-1)**2.0
            z2 = x1(j)**2.0
            t1 = y1(j-1)
            t2 = y1(j)

            A = 0.0
            B = (t1 - t2)/(z1 - z2) - A*(z1 + z2)
            C = t2 - A*(z2**2.0) - B*z2

            y2(i) = A*(t4**4.0) + B*(t4**2.0) + C
            EXIT
        endif

        if (t4 > z3 .and. j == Nx1-1) then
            t1 = y1(j-1)
            t2 = y1(j)
            t3 = y1(j+1)
 
            A = (t3 - t2 - (z3 - z2)*(t1 - t2)/(z1 - z2)) / ((z3 - z2)*(z3 - z1))
            B = (t1 - t2)/(z1 - z2) - A*(z1 + z2)
            C = t2 - A*(z2**2.0) - B*z2

            y2(i) = A*(t4**2.0) + B*t4 + C
            EXIT
        endif

    enddo
enddo

end subroutine qinterp

!-----------------------------------------------------
subroutine integrcc(nx, x, y, sy)

integer, intent(in) :: nx
double precision, intent(in), dimension(nx) :: x, y
double precision, intent(out) :: sy(nx)

integer :: i
double precision :: drho, y1tmp

sy = 0.
do i = 2, nx
    drho = (x(i) - x(i-1))
    y1tmp = (y(i) + y(i-1))/2.
    sy(i) = sy(i-1) + y1tmp*drho
enddo

end subroutine integrcc

!------------------------------------------------------
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

end subroutine derivcc

!------------------------------------------------------
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
P(3) = y(3) - P(1)*x(3)**2. - P(2)*x(3)

end subroutine polyfitcc

!---------------------------------------------------------------------=|
! Assume that x is of r-type, i.e. interpolation in 0 has zero odd derivatives
subroutine EXTRAP(x_input, y_input, x_extrap, j_extrap, y_extrap, ex_order, nagrid)

integer, intent(in) :: ex_order, j_extrap, nagrid
double precision, intent(in)  :: x_input(nagrid), x_extrap, y_input(nagrid)
double precision, intent(out) :: y_extrap

integer :: k1, k2, k3, jsign
double precision :: P(3)
            
if (j_extrap == nagrid) jsign = -1
if (j_extrap == 1) jsign = 1
      
! Constant interpolation
if (ex_order == 0) then
    k1 = j_extrap
    P(1) = y_input(k1)
    y_extrap = P(1)
endif

if (ex_order == 1) then ! Linear extrapolation
    if (jsign < 0) then
        k1 = j_extrap + jsign
        k2 = j_extrap
        call polyfitcc_1(x_input(k1: k2), y_input(k1: k2), P(1: 2))
        y_extrap = P(1)*x_extrap + P(2)
    endif
    if (jsign > 0) then
        k1 = j_extrap
        k2 = j_extrap + 1
        call polyfitcc_1(x_input(k1: k2), y_input(k1: k2), P(1: 2))
        y_extrap = P(1)*x_extrap + P(2)
    endif
else if (ex_order == 2) then! Quadratic extrapolation
    if (jsign < 0) then
        k1 = j_extrap + jsign*2
        k2 = j_extrap + jsign
        k3 = j_extrap
        call polyfitcc(x_input(k1: k3), y_input(k1: k3), P)
        y_extrap = P(1) * x_extrap**2.0 + P(2)*x_extrap + P(3)
    else if (jsign > 0) then
        k1 = j_extrap
        k2 = j_extrap + jsign
        k3 = j_extrap + jsign*2
        call polyfitcc(x_input(k1: k3), y_input(k1: k3), P)
        y_extrap = P(1) * x_extrap**2.0 + P(2)*x_extrap + P(3)
    endif
endif

return
end subroutine EXTRAP

!---------------------------------------------------------------------=|
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
!---------------------------------------------------------------------=|
!Efable DERIV computes first or second derivative over x
!
!  x_input: x_variable
!  x_type for GRP style-grid is: 1 main grid -> shifted grid (deriv), 2 shifted grid -> main grid (deriv)
!  y_input: y_variable
!  yd_output: derivative
!  order_d: 1st or 2nd derivative
!  nagrid: number of grid points
! iextrap: 1 if yes interpolate last grid point, 0 do not interpolate last grid point
!---------------------------------------------------------------------=|
subroutine DERIV(x_input, x_output, x_type, y_input, yd_output, order_d, nagrid, iextrap)

integer, intent(in) :: order_d, x_type, nagrid, iextrap
double precision, intent(in) , dimension(nagrid) :: x_input, x_output, y_input
double precision, intent(out), dimension(nagrid) :: yd_output

integer :: j, j_end, OEXTRAP
double precision dx, dy, y0, P(3)
     
OEXTRAP = 1

if (iextrap == 1) j_end=1
if (iextrap == 0) j_end=0
  
! dy/dx
if (order_d == 1) then

    if (x_type == 1) then
        do j=1, nagrid - j_end
            dx = x_input(j+1) - x_input(j)
            dy = y_input(j+1) - y_input(j)
            yd_output(j) = dy/dx
        enddo
    endif

    if (x_type == 2) then
        j = 1
        call polyfitcc(x_input(1: 3), y_input(1: 3), P)
        y0 = P(3)
        dx = x_input(1)
        dy = (y_input(1) - y0)
        yd_output(j) = dy/dx
        do j=2, nagrid
            dx = x_input(j) - x_input(j-1)
            dy = y_input(j) - y_input(j-1)
            yd_output(j) = dy/dx
        enddo
    endif

else if (order_d == 2) then ! d2y/dx^2
    if (x_type == 1) then
        j = 1
        dx = (x_input(j+1) - x_input(j))**2.0
        dy = (y_input(j+1) - y_input(j))
        yd_output(j) = dy/dx

        do j=2, nagrid - j_end
            dx = (x_input(j+1) - x_input(j))**2.0
            dy = (y_input(j+1) - 2.*y_input(j) + y_input(j-1))
            yd_output(j) = dy/dx
        enddo
    endif
    if (x_type == 2) then
! First interpolate shifted variable to zero
        j = 1
        call polyfitcc(x_input(1:3), y_input(1:3), P)
        y0 = P(3)
        dx = x_input(j)**2.0
        dy = (y_input(j+1) - 2.0*y_input(j) + y0)
        yd_output(j) = dy/dx
        do j=2, nagrid - j_end
            dx = (x_input(j+1) - x_input(j))**2.0
            dy = (y_input(j+1) - 2.0*y_input(j) + y_input(j-1))
            yd_output(j) = dy/dx
        enddo
    endif
endif

! Interpolate to last grid point
if (iextrap == 1 .and. order_d == 1) then
    j = nagrid
    call polyfitcc(x_input(j-2: j), y_input(j-2: j), P)
    yd_output(j) = 2.*P(1)*x_output(j) + P(2)
endif

if (iextrap == 1 .and. order_d == 2) then
    j = nagrid
    call polyfitcc(x_input(j-2: j), y_input(j-2: j), P)
    yd_output(j) = 2.*P(1)
endif

return
end subroutine DERIV

!---------------------------------------------------------------------=|
! INTEGR_EF computes integrals over x
!  x_input: x_variable
!  x_type for GRP style-grid is: 1 main grid -> shifted grid (deriv), 2 shifted grid -> main grid (deriv)
!  y_input: y_variable
!  ys_output: integral
!  nagrid: number of grid points
!---------------------------------------------------------------------=|
subroutine INTEGR(x_input, x_type, y_input, ys_output, nagrid)

integer, intent(in) :: x_type, nagrid
double precision, intent(in) , dimension(nagrid) :: x_input, y_input
double precision, intent(out), dimension(nagrid) :: ys_output

integer :: j
double precision :: P(3)

! Normalized grid , GRP style
if (x_type == 1) then
    ys_output(1) = y_input(1)*x_input(1)
    do j=2, nagrid
        ys_output(j) = ys_output(j-1) + y_input(j)*(x_input(j) - x_input(j-1))
    enddo
else if (x_type == 2) then ! First interpolate shifted variable to zero
    call polyfitcc(x_input(1: 3), y_input(1: 3), P)
    ys_output(1) = P(3)*x_input(1)
    do j=2, nagrid
        ys_output(j) = ys_output(j-1) + y_input(j-1)*(x_input(j) - x_input(j-1))
    enddo
endif

return
end subroutine INTEGR

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


end module numerical_tools
