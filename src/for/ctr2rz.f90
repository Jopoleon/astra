SUBROUTINE quad_int(x_out, x_in, y_in, y_out)

implicit none

double precision, intent(in) :: x_out, x_in(3), y_in(3)
double precision, intent(out) :: y_out

double precision :: ratio, A, B, C

ratio = (y_in(3) - y_in(2))/(x_in(3) - x_in(2))
A     = (y_in(1) - y_in(2))/(x_in(1) - x_in(2)) - ratio
A = A/(x_in(1) - x_in(3))
B = ratio - A*(x_in(3) + x_in(2))
C = y_in(1) - A*x_in(1)**2 - B*x_in(1)

y_out =  A*x_out**2 + B*x_out + C

return
end SUBROUTINE quad_int

SUBROUTINE lin_int(x_out, x_in, y_in, y_out)

implicit none

double precision, intent(in) :: x_out, x_in(2), y_in(2)
double precision, intent(out) :: y_out

double precision m

m = (y_in(2)-y_in(1)) / (x_in(2)-x_in(1))
y_out =  m*(x_out-x_in(2)) + y_in(2)

return
end SUBROUTINE lin_int

!----------------------------------------------------------
SUBROUTINE interp_norm(n_the, r_in, z_in, th_ref, norm_out)

use const_inc, only: GP2

implicit none

integer, intent(in) :: n_the
double precision, intent(in) :: th_ref
double precision, dimension(n_the), intent(in) :: r_in, z_in
double precision, intent(out) :: norm_out

integer :: jmin(1), jmid
double precision, dimension(3) :: theta3, norm3
double precision, dimension(n_the) :: theta

theta = ATAN2(z_in, r_in) - th_ref
jmin = MINLOC(ABS(theta))
jmid = jmin(1)
theta3(2) = theta(jmid)
norm3(2) = r_in(jmid)**2 + z_in(jmid)**2
if (jmid == 1) then
    norm3(1) = r_in(n_the)**2 + z_in(n_the)**2
    norm3(3) = r_in(2    )**2 + z_in(2    )**2
    theta3(1) = theta(n_the) - GP2
    theta3(3) = theta(2)
else if (jmid == n_the) then
    norm3(1) = r_in(n_the-1)**2 + z_in(n_the-1)**2
    norm3(3) = r_in(1      )**2 + z_in(1      )**2
    theta3(1) = theta(n_the-1)
    theta3(3) = theta(1) + GP2
else
    norm3(1) = r_in(jmid-1)**2 + z_in(jmid-1)**2
    norm3(3) = r_in(jmid+1)**2 + z_in(jmid+1)**2
    theta3(1) = theta(jmid-1)
    theta3(3) = theta(jmid+1)
endif

norm3 = SQRT(norm3)
CALL quad_int(0., theta3, norm3, norm_out)

return
end SUBROUTINE interp_norm

!----------------------------------------------------------
SUBROUTINE ctr2rz_b(n_rho, n_the, pf1d, ipol, X, Y, Nrrect, Nzrect, Rgrid, Zgrid, &
    pfm, B_R, B_Z, B_T)

use const_inc, only: GP2
use numerical_tools, only: deriv_cde

implicit none

integer, intent(in) :: n_rho, n_the, Nrrect, Nzrect
double precision, intent(in) :: Rgrid(Nrrect), Zgrid(Nzrect)
double precision, intent(in), dimension(n_rho) :: ipol, pf1d
double precision, intent(in), dimension(n_rho, n_the) :: X, Y

double precision, intent(out), dimension(Nrrect, Nzrect) :: pfm, B_R, B_Z, B_T

integer :: jr, jz, irho, jmin(2), jleft, jrho
double precision :: norm, tht, Rpos, Zpos
double precision, dimension(Nrrect) :: Rrect, dum_r, pfmr
double precision, dimension(Nzrect) :: Zrect, dum_z, pfmz
double precision, dimension(Nrrect, Nzrect) :: ipolp
double precision, dimension(n_the-1) :: Rctr1, Zctr1
double precision, dimension(n_rho-1, n_the-1) :: Rctr, Zctr, rdist, zdist
double precision, dimension(3) :: norm3, pf3, rb3

! Use reference arounf Rmag, Zmag
do jrho = 1, n_rho-1
    Rctr(jrho, 1:n_the-1) = X(jrho+1, 1:n_the-1) - X(1, 1) 
    Zctr(jrho, 1:n_the-1) = Y(jrho+1, 1:n_the-1) - Y(1, 1) 
enddo
Rrect = Rgrid - X(1, 1)
Zrect = Zgrid - Y(1, 1)

! Biquadratic interpolation

do jr=1, Nrrect
    Rpos = Rrect(jr)
    rdist = (Rctr - Rpos)**2
    do jz=1, Nzrect
        Zpos = Zrect(jz)
        zdist = (Zctr - Zpos)**2
        tht  = ATAN2(Zpos, Rpos)
        norm = SQRT(Rpos**2 + Zpos**2)

        jmin = MINLOC(rdist + zdist)
        irho = jmin(1)
        if (irho == 1) then
            norm3(1) = 0.
            Rctr1 = Rctr(1, :)
            Zctr1 = Zctr(1, :)
            CALL interp_norm(n_the - 1, Rctr1, Zctr1, tht, norm3(2))
            Rctr1 = Rctr(2, :)
            Zctr1 = Zctr(2, :)
            CALL interp_norm(n_the - 1, Rctr1, Zctr1, tht, norm3(3))
            pf3 = pf1d(1:3)
            rb3 = ipol(1:3)
        else
            if (irho == n_rho-1) then
                jleft = n_rho - 3
            else
                jleft = irho - 1
            endif
            do jrho = 1, 3
                Rctr1 = Rctr(jleft + jrho - 1, :)
                Zctr1 = Zctr(jleft + jrho - 1, :)
                CALL interp_norm(n_the - 1, Rctr1, Zctr1, tht, norm3(jrho))
                pf3(jrho) = pf1d(jleft + jrho)
                rb3(jrho) = ipol(jleft + jrho)
            enddo
        endif
! Extrapolate flat
        if ( (irho == n_rho-1) .and. (norm > norm3(3)) ) then
            pfm  (jr, jz) = 2.*pf3(3) - pf3(2)
            ipolp(jr, jz) = 2.*rb3(3) - rb3(2)
        else
            CALL quad_int(norm, norm3, pf3, pfm(jr, jz))
            CALL quad_int(norm, norm3, rb3, ipolp(jr, jz))
        endif
    enddo
enddo

! Magnetic field components
do jz=1, Nzrect
    pfmr = pfm(:, jz)
    call DERIV_CDE(Rgrid, 3, pfmr, dum_r, Nrrect)
    B_Z(:, jz) = -dum_r/(GP2*Rgrid)
    B_T(:, jz) = ipolp(:, jz)/Rgrid
enddo
do jr=1, Nrrect
    pfmz = pfm(jr, :)
    call DERIV_CDE(Zgrid, 3, pfmz, dum_z, Nzrect)
    B_R(jr, :) = dum_z/(GP2*Rgrid(jr))
enddo

return
end SUBROUTINE ctr2rz_b
!----------------------------------------------------------
SUBROUTINE ctr2rz_fun3(n_rho, n_the, f1d, X, Y, Nrrect, Nzrect, Rgrid, Zgrid, f2d)

use const_inc, only: GP2

implicit none

integer, intent(in) :: n_rho, n_the, Nrrect, Nzrect
double precision, intent(in) :: Rgrid(Nrrect), Zgrid(Nzrect)
double precision, intent(in), dimension(n_rho) :: f1d
double precision, intent(in), dimension(n_rho, n_the) :: X, Y

double precision, intent(out), dimension(Nrrect, Nzrect) :: f2d

integer :: jr, jz, jmin, imin, jstart, jend, jwhere, jscale, jcount
double precision :: norm, tht, Rpos, Zpos, Rgeo, Zgeo
double precision, dimension(Nrrect) :: Rrect
double precision, dimension(Nzrect) :: Zrect
double precision, dimension(n_rho, n_the+1) :: Rctr, Zctr, rdist
double precision, dimension(n_the+1) :: thet
double precision, dimension(3) :: x3, y3, yout

f2d = 1.e8

Rgeo = X(1, 1)
Zgeo = Y(1, 1)

Rctr(:, 1:n_the) = X - Rgeo 
Zctr(:, 1:n_the) = Y - Zgeo 
Rctr(:, n_the+1) = Rctr(:, 1)
Zctr(:, n_the+1) = Zctr(:, 1)

rdist = sqrt(Rctr**2 + Zctr**2)
thet = atan2(Zctr(2, :), Rctr(2, :))
do jz=1, n_the+1 
    if (thet(jz) < 0.) thet(jz) = thet(jz) + GP2
enddo
thet(n_the+1) = thet(1) + GP2

Rrect = Rgrid - Rgeo
Zrect = Zgrid - Zgeo

! Biquadratic interpolation
do jz=1, Nzrect
    do jr=1, Nrrect
  
! grid point under investigation
        Rpos = Rrect(jr)
        Zpos = Zrect(jz)
        norm = SQRT(Rpos**2 + Zpos**2)
        tht  = ATAN2(Zpos, Rpos)
        if (tht < 0.) tht = tht + GP2
        jmin = minloc(abs(tht - thet), 1)            ! nearest polar point theta
        imin = minloc(abs(norm - rdist(:, jmin)), 1) ! nearest polar point radius
        if (imin < n_rho) then ! internal points
            jstart = -1
            jend = 1
            jscale = 3
            if (imin == 1) then
                jstart = 0
                jend = 2
                jscale = 3
            endif
            jcount = 0
            do jwhere=imin+jstart, imin+jend
                jcount = jcount + 1
                if (jmin == n_the+1) then
                    x3(1) = thet(jmin-1)
                    x3(2) = thet(jmin)
                    x3(3) = thet(2) + GP2
                    y3(1) = rdist(jwhere, jmin-1)
                    y3(2) = rdist(jwhere, jmin)
                    y3(3) = rdist(jwhere, 2)		
                endif
                if (jmin == 1) then
                    x3(1) = thet(n_the) - GP2
                    x3(2) = thet(1)
                    x3(3) = thet(2)
                    y3(1) = rdist(jwhere, n_the)
                    y3(2) = rdist(jwhere, 1)
                    y3(3) = rdist(jwhere, 2)		
                endif
                if (jmin > 1 .and. jmin < n_the+1) then
                    x3(1) = thet(jmin-1)
                    x3(2) = thet(jmin)
                    x3(3) = thet(jmin+1)
                    y3(1) = rdist(jwhere, jmin-1)
                    y3(2) = rdist(jwhere, jmin)
                    y3(3) = rdist(jwhere, jmin+1)
                endif
                call quad_int(tht, x3(1:3), y3(1:3), yout(jcount)) !rdist a tht, imin
            enddo
            if (imin == 1) then
                x3(1) = yout(1)
                y3(1) = f1d(imin)
                x3(2) = yout(2)
                y3(2) = f1d(imin+1)
                x3(3) = yout(3)
                y3(3) = f1d(imin+2)
                call quad_int(norm, x3(1:3), y3(1:3), f2d(jr, jz))
            else
                x3(1) = yout(1)
                y3(1) = f1d(imin-1)
                x3(2) = yout(2)
                y3(2) = f1d(imin)
                x3(3) = yout(3)
                y3(3) = f1d(imin+1)
                call quad_int(norm, x3(1:3), y3(1:3), f2d(jr, jz))
            endif
        else ! boundary and external points
            jcount = 0
            do jwhere=imin-2, imin
                jcount = jcount + 1
                if (jmin == n_the+1) then
                    x3(1) = thet(jmin-1)
                    x3(2) = thet(jmin)
                    x3(3) = thet(2) + GP2		
                    y3(1) = rdist(jwhere, jmin-1)
                    y3(2) = rdist(jwhere, jmin)
                    y3(3) = rdist(jwhere, 2)		
                endif
                if (jmin == 1) then
                    x3(1) = thet(n_the) - GP2
                    x3(2) = thet(1)
                    x3(3) = thet(2)		
                    y3(1) = rdist(jwhere, n_the)
                    y3(2) = rdist(jwhere, 1)
                    y3(3) = rdist(jwhere, 2)		
                endif	 
                if (jmin > 1 .and. jmin < n_the+1) then
                    x3(1) = thet(jmin-1)
                    x3(2) = thet(jmin)
                    x3(3) = thet(jmin+1)		
                    y3(1) = rdist(jwhere, jmin-1)
                    y3(2) = rdist(jwhere, jmin)
                    y3(3) = rdist(jwhere, jmin+1)		
                endif	 
                call quad_int(tht, x3(1:3), y3(1:3), yout(jcount))
            enddo
            x3(1) = yout(1)
            y3(1) = f1d(imin-2)
            x3(2) = yout(2)
            y3(2) = f1d(imin-1)
            x3(3) = yout(3)
            y3(3) = f1d(imin)
            call lin_int(norm, x3(2:3), y3(2:3), f2d(jr, jz))
        endif
    enddo
enddo

return
end SUBROUTINE ctr2rz_fun3


!----------------------------------------------------------
SUBROUTINE ctr2rz_fun(n_rho, n_the, f1d, X, Y, Nrrect, Nzrect, Rgrid, Zgrid, f2d)

implicit none

integer, intent(in) :: n_rho, n_the, Nrrect, Nzrect
double precision, intent(in) :: Rgrid(Nrrect), Zgrid(Nzrect)
double precision, intent(in), dimension(n_rho) :: f1d
double precision, intent(in), dimension(n_rho, n_the) :: X, Y

double precision, intent(out), dimension(Nrrect, Nzrect) :: f2d

integer :: jr, jz, irho, jmin(2), jleft, jrho

double precision :: norm, tht, Rpos, Zpos
double precision, dimension(Nrrect) :: Rrect
double precision, dimension(Nzrect) :: Zrect
double precision, dimension(n_the-1) :: Rctr1, Zctr1
double precision, dimension(n_rho-1, n_the-1) :: Rctr, Zctr, rdist, zdist
double precision, dimension(3) :: norm3, f3

f2d = 1.e8

! Use reference arounf Rmag, Zmag
DO jrho = 1, n_rho-1
    Rctr(jrho, 1:n_the-1) = X(jrho+1, 1:n_the-1) - X(1, 1) 
    Zctr(jrho, 1:n_the-1) = Y(jrho+1, 1:n_the-1) - Y(1, 1) 
ENDDO
Rrect = Rgrid - X(1, 1)
Zrect = Zgrid - Y(1, 1)

! Biquadratic interpolation
do jr=1, Nrrect
    Rpos = Rrect(jr)
    rdist = (Rctr - Rpos)**2
    do jz=1, Nzrect
        Zpos = Zrect(jz)
        zdist = (Zctr - Zpos)**2
        tht  = ATAN2(Zpos, Rpos)
        norm = SQRT(Rpos**2 + Zpos**2)

        jmin = MINLOC(rdist + zdist)
        irho = jmin(1)
        if (irho == 1) then
            norm3(1) = 0.
            Rctr1 = Rctr(1, :)
            Zctr1 = Zctr(1, :)
            CALL interp_norm(n_the - 1, Rctr1, Zctr1, tht, norm3(2))
            Rctr1 = Rctr(2, :)
            Zctr1 = Zctr(2, :)
            CALL interp_norm(n_the - 1, Rctr1, Zctr1, tht, norm3(3))
            f3 = f1d(1:3)
        else
            if (irho == n_rho-1) then
                jleft = n_rho - 3
            else
                jleft = irho - 1
            endif
            do jrho = 1, 3
                Rctr1(:) = Rctr(jleft + jrho - 1, :)
                Zctr1(:) = Zctr(jleft + jrho - 1, :)
                CALL interp_norm(n_the - 1, Rctr1, Zctr1, tht, norm3(jrho))
                f3(jrho) = f1d(jleft + jrho)
            enddo
        endif
        if ( (irho == n_rho-1) .and. (norm > norm3(3)) ) then
!            f2d(jr, jz) = 2.*f3(3) - f3(2)
         !linear extrapolation (for Rabbit, such that orbits outside of the last closed flux surface can be calculated.)
            CALL lin_int(norm, norm3(2:3), f3(2:3), f2d(jr, jz))
        else
            CALL quad_int(norm, norm3, f3, f2d(jr, jz))
        endif
    enddo
enddo

return
end SUBROUTINE ctr2rz_fun

!----------------------------------------------------------
SUBROUTINE ctr2rz
! 2D interpolation of Psi, Fdia from contours(rho, theta) to Cartesian R, z 2D-grid
  
use parameters_a2equil, only: equil_now
use const_inc, only: IFBEY

implicit none

integer :: jr, jz, irho, jleft, jrho, n_rho, n_the, nr_rect, nz_rect
integer, dimension(2) :: jmin
double precision :: norm, tht, Rpos, Zpos, Rmag, Zmag
double precision, allocatable, dimension(:) :: Rctr1, Zctr1
double precision, allocatable, dimension(:, :) :: Rctr, Zctr, rdist, zdist
double precision, dimension(3) :: norm3, pf3, rb3

! Use reference arounf Rmag, Zmag
n_rho = SIZE(equil_now%coord_sys%position%r, 1)
n_the = SIZE(equil_now%coord_sys%position%r, 2)
nr_rect = SIZE(equil_now%eqgeometry%rectgrid%r2d)
nz_rect = SIZE(equil_now%eqgeometry%rectgrid%z2d)
if (.not. allocated(Rctr1)) then
    allocate(Rctr1(n_the-1))
    allocate(Zctr1(n_the-1))
    allocate(Rctr(n_rho-1, n_the-1), rdist(n_rho-1, n_the-1), &
             Zctr(n_rho-1, n_the-1), zdist(n_rho-1, n_the-1))
endif
 
Rmag = equil_now%coord_sys%position%r(1, 1)
Zmag = equil_now%coord_sys%position%z(1, 1)
Rctr = equil_now%coord_sys%position%r(2:, :n_the-1) - Rmag
Zctr = equil_now%coord_sys%position%z(2:, :n_the-1) - Zmag

! Biquadratic interpolation
do jr=1, nr_rect
    Rpos = equil_now%eqgeometry%rectgrid%r2d(jr) - Rmag
    rdist = (Rctr - Rpos)**2
    do jz=1, nz_rect
        Zpos = equil_now%eqgeometry%rectgrid%z2d(jz) - Zmag
        zdist = (Zctr - Zpos)**2
        tht  = ATAN2(Zpos, Rpos)
        norm = SQRT(Rpos**2 + Zpos**2)

        jmin = MINLOC(rdist + zdist)
        irho = jmin(1)
        if (irho == 1) then
            norm3(1) = 0.
            Rctr1 = Rctr(1, :)
            Zctr1 = Zctr(1, :)
            CALL interp_norm(n_the - 1, Rctr1, Zctr1, tht, norm3(2))
            Rctr1 = Rctr(2, :)
            Zctr1 = Zctr(2, :)
            CALL interp_norm(n_the - 1, Rctr1, Zctr1, tht, norm3(3))
            pf3 = equil_now%profiles_1d%psi(1:3)
            rb3 = equil_now%profiles_1d%F_dia(1:3)
        else
            if (irho == n_rho-1) then
                jleft = n_rho - 3
            else
                jleft = irho - 1
            endif
            do jrho = 1, 3
                Rctr1(:) = Rctr(jleft + jrho - 1, :)
                Zctr1(:) = Zctr(jleft + jrho - 1, :)
                CALL interp_norm(n_the - 1, Rctr1, Zctr1, tht, norm3(jrho))
                pf3(jrho) = equil_now%profiles_1d%psi(jleft + jrho)
                rb3(jrho) = equil_now%profiles_1d%F_dia(jleft + jrho)
            enddo
        endif
        if ( (irho == n_rho-1) .and. (norm > norm3(3)) ) then
!            f2d(jr, jz) = 2.*pf3(3) - pf3(2)
         !linear extrapolation (for Rabbit, such that orbits outside of the last closed flux surface can be calculated.)
            if (IFBEY == 0) CALL lin_int(norm, norm3(2:3), pf3(2:3), equil_now%eqgeometry%rectgrid%psirz2d(jr, jz))
            CALL lin_int(norm, norm3(2:3), rb3(2:3), equil_now%eqgeometry%rectgrid%fdia2d(jr, jz))
        else
            if (IFBEY == 0) CALL quad_int(norm, norm3, pf3, equil_now%eqgeometry%rectgrid%psirz2d(jr, jz))
            CALL quad_int(norm, norm3, rb3, equil_now%eqgeometry%rectgrid%fdia2d(jr, jz))
        endif
    enddo
enddo

return
end SUBROUTINE ctr2rz
