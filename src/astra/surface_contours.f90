module surface_contours

implicit none
contains

!---------------------------------------------------------------------
    subroutine quad_int(x_out, x_in, y_in, y_out)

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
    end subroutine quad_int

!---------------------------------------------------------------------
    subroutine lin_int(x_out, x_in, y_in, y_out)

    double precision, intent(in) :: x_out, x_in(2), y_in(2)
    double precision, intent(out) :: y_out

    double precision m

    m = (y_in(2)-y_in(1)) / (x_in(2)-x_in(1))
    y_out =  m*(x_out-x_in(2)) + y_in(2)

    return
    end subroutine lin_int

!---------------------------------------------------------------------
    subroutine interp_norm(n_the, r_in, z_in, th_ref, norm_out)

    use pi_const, only: GP2

    integer, intent(in) :: n_the
    double precision, intent(in) :: th_ref
    double precision, dimension(n_the), intent(in) :: r_in, z_in
    double precision, intent(out) :: norm_out

    integer :: jmin(1), jmid, ind3(3)
    double precision, dimension(3) :: theta3, norm3
    double precision, dimension(n_the) :: theta

    theta = ATAN2(z_in, r_in) - th_ref
    jmin = MINLOC(ABS(theta))
    jmid = jmin(1)
    ind3 = (/ MOD(jmid-2 + n_the, n_the) + 1, jmid, MOD(jmid, n_the) + 1 /)

    theta3 = theta(ind3)
    if (jmid == 1) then
        theta3(1) = theta(n_the) - GP2
    else if (jmid == n_the) then
        theta3(3) = theta(1) + GP2
    endif

    norm3 = SQRT(r_in(ind3)**2 + z_in(ind3)**2)

    CALL quad_int(0., theta3, norm3, norm_out)

    return
    end subroutine interp_norm

!---------------------------------------------------------------------
    subroutine ctr2rz_b(Rgrid, Zgrid, pfm, fdiam, pf, fdia, B_R, B_Z, B_T)

    use pi_const, only: GP2
    use numerical_tools, only: deriv_cde
    use parameters_a2equil, only: equil_now

    double precision, intent(in), dimension(:) :: Rgrid, Zgrid
    double precision, intent(in), optional, dimension(*) :: pf, fdia
    double precision, intent(out), dimension(:, :), allocatable :: pfm
    double precision, intent(out), dimension(:, :), allocatable, optional :: &
         fdiam, B_R, B_Z, B_T

    logical :: b_flag, fd_flag, fd_b_flag
    integer :: NrRect, NzRect, n_rho, n_the
    integer :: jr, jz, irho, jmin(2), jleft, jrho
    double precision :: norm, tht, Rpos, Zpos, Rmag, Zmag
    double precision, dimension(:), allocatable :: Rrect, Zrect, &
        dpf_dr, dpf_dz, Rctr1, Zctr1, pf1d, fdia1d
    double precision, dimension(:, :), allocatable :: fdia2d, Rctr, Zctr, rdist, zdist
    double precision, dimension(3) :: norm3, pf3, rb3

    b_flag  = (present(B_r))
    fd_flag = (present(fdiam)) 
    fd_b_flag = (b_flag .or. fd_flag) 

! Allocate variables

    n_rho = SIZE(equil_now%coord_sys%position%r, 1)
    n_the = SIZE(equil_now%coord_sys%position%r, 2)
    NrRect = SIZE(Rgrid)
    NzRect = SIZE(Zgrid)
    allocate(Rctr1(n_the-1), Zctr1(n_the-1))
    allocate(Rctr(n_rho-1, n_the-1), rdist(n_rho-1, n_the-1))
    allocate(Zctr(n_rho-1, n_the-1), zdist(n_rho-1, n_the-1))
    allocate(pf1d(n_rho), fdia1d(n_rho))
    allocate(Rrect(NrRect), dpf_dr(NrRect))
    allocate(Zrect(NzRect), dpf_dz(NzRect))
    allocate(pfm(NrRect, NzRect), fdia2d(NrRect, NzRect))
    if (present(pf)) then
        pf1d = pf(1: n_rho)
    else
        pf1d = equil_now%profiles_1d%psi
    endif
    if (present(fdia)) then
        fdia1d = fdia(1: n_rho)
    else
        fdia1d = equil_now%profiles_1d%F_dia
    endif

! Shift by Rmag, Zmag

    Rmag = equil_now%coord_sys%position%r(1, 1)
    Zmag = equil_now%coord_sys%position%z(1, 1)
    Rctr = equil_now%coord_sys%position%r(2:, :n_the-1) - Rmag
    Zctr = equil_now%coord_sys%position%z(2:, :n_the-1) - Zmag
    Rrect = Rgrid(1: NrRect) - Rmag
    Zrect = Zgrid(1: NzRect) - Zmag

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
                if (fd_b_flag) rb3 = fdia1d(1:3)
            else
                if (irho == n_rho-1) then
                    jleft = n_rho - 3
                else
                    jleft = irho - 1
                endif
                do jrho=1, 3
                    Rctr1 = Rctr(jleft + jrho - 1, :)
                    Zctr1 = Zctr(jleft + jrho - 1, :)
                    CALL interp_norm(n_the - 1, Rctr1, Zctr1, tht, norm3(jrho))
                    pf3(jrho) = pf1d(jleft + jrho)
                    if (fd_b_flag) rb3(jrho) = fdia1d(jleft + jrho)
                enddo
            endif
! Outside sep: linear extrapolation (for Rabbit, such that orbits outside of the last closed flux surface can be calculated.)
            if ( (irho == n_rho-1) .and. (norm > norm3(3)) ) then
                CALL lin_int(norm, norm3(2:3), pf3(2:3), pfm(jr, jz))
                if (fd_b_flag) CALL lin_int(norm, norm3(2:3), rb3(2:3), fdia2d(jr, jz))
            else
                CALL quad_int(norm, norm3, pf3, pfm(jr, jz))
                if (fd_b_flag) CALL quad_int(norm, norm3, rb3, fdia2d(jr, jz))
            endif
        enddo
    enddo

! Magnetic field components

    if (b_flag) then
        allocate(B_R(NrRect, NzRect), B_Z(NrRect, NzRect), B_T(NrRect, NzRect))
        do jz=1, Nzrect
            call DERIV_CDE(Rgrid(1: NrRect), 3, pfm(:, jz), dpf_dr, Nrrect)
            B_Z(:, jz) = -dpf_dr/(GP2*Rgrid(1: NrRect))
            B_T(:, jz) = fdia2d(:, jz)/Rgrid(1: NrRect)
        enddo
        do jr=1, Nrrect
            call DERIV_CDE(Zgrid(1: NzRect), 3, pfm(jr, :), dpf_dz, Nzrect)
            B_R(jr, :) = dpf_dz/(GP2*Rgrid(jr))
        enddo
    endif

    if (fd_flag) then
        allocate(fdiam(NrRect, NzRect))
        fdiam = fdia2d
    endif

    deallocate(Rctr1, Zctr1, Rctr, Zctr, rdist, zdist, pf1d, fdia1d, fdia2d)
    deallocate(Rrect, Zrect, dpf_dr, dpf_dz)

    return
    end subroutine ctr2rz_b

!---------------------------------------------------------------------
    subroutine ctr2rz_fun3(n_rho, n_the, f1d, X, Y, Nrrect, Nzrect, Rgrid, Zgrid, f2d)

    use pi_const, only: GP2

    integer, intent(in) :: n_rho, n_the, Nrrect, Nzrect
    double precision, intent(in) :: Rgrid(Nrrect), Zgrid(Nzrect)
    double precision, intent(in), dimension(n_rho) :: f1d
    double precision, intent(in), dimension(n_rho, n_the) :: X, Y

    double precision, intent(out), dimension(Nrrect, Nzrect) :: f2d

    integer :: jr, jz, jmin, imin, jstart, jend, jwhere, jcount
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
                if (imin == 1) then
                    jstart = 0
                    jend = 2
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
                    else if (jmin == 1) then
                        x3(1) = thet(n_the) - GP2
                        x3(2) = thet(1)
                        x3(3) = thet(2)
                        y3(1) = rdist(jwhere, n_the)
                        y3(2) = rdist(jwhere, 1)
                        y3(3) = rdist(jwhere, 2)
                    else
                        x3 = thet(jmin-1: jmin+1)
                        y3 = rdist(jwhere, jmin-1: jmin+1)
                    endif
                    call quad_int(tht, x3, y3, yout(jcount)) !rdist a tht, imin
                enddo
                x3 = yout
                if (imin == 1) then
                    y3 = f1d(imin: imin+2)
                else
                    y3 = f1d(imin-1: imin+1)
                endif
                call quad_int(norm, x3, y3, f2d(jr, jz))
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
                    else if (jmin == 1) then
                        x3(1) = thet(n_the) - GP2
                        x3(2) = thet(1)
                        x3(3) = thet(2)
                        y3(1) = rdist(jwhere, n_the)
                        y3(2) = rdist(jwhere, 1)
                        y3(3) = rdist(jwhere, 2)
                    else
                        x3 = thet(jmin-1: jmin+1)
                        y3 = rdist(jwhere, jmin-1: jmin+1)
                    endif
                    call quad_int(tht, x3, y3, yout(jcount))
                enddo
                x3 = yout
                y3(2: 3) = f1d(imin-1: imin)
                call lin_int(norm, x3(2:3), y3(2:3), f2d(jr, jz))
            endif
        enddo
    enddo

    return
    end subroutine ctr2rz_fun3

!---------------------------------------------------------------------
    subroutine ctr2rz()

    use parameters_a2equil, only: equil_now
    use scalars, only: IFBEY

    double precision, dimension(:, :), allocatable :: pfm, fdiam

    call ctr2rz_b(equil_now%eqgeometry%rectgrid%r2d, equil_now%eqgeometry%rectgrid%z2d, pfm, fdiam=fdiam)

    equil_now%eqgeometry%rectgrid%fdia2d = fdiam
    if (IFBEY == 0) equil_now%eqgeometry%rectgrid%psirz2d = pfm

    return
    end subroutine ctr2rz

end module surface_contours
