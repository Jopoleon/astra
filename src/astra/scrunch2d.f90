module scrunch2d

use numerical_tools, only: argsort, sort_by_index, mom2rz

implicit none

contains

!---------------------------------------------------------------------
    subroutine surf2surf(mpol, ntheta_out, r_in, z_in, r_out, z_out)

    integer, intent(in) :: mpol, ntheta_out
    double precision, intent(in), dimension(:) :: r_in, z_in
    double precision, intent(out), dimension(ntheta_out) :: r_out, z_out

    double precision, dimension(mpol) :: rcos, rsin, zcos, zsin

    call scrunch(mpol, r_in, z_in, ntheta_out, rcos, rsin, zcos, zsin)
    call mom2rz(mpol, ntheta_out, rcos, rsin, zcos, zsin, r_out, z_out)

    end subroutine surf2surf

!---------------------------------------------------------------------
    subroutine scrunch(mpol, r_in, z_in, ntheta_out, rcos, rsin, zcos, zsin)

    integer, parameter :: niter=1000, max_iterate=10, n_resets=100, pexp=4
    double precision, parameter:: ftol=1d-8, bmax=0.15d0, pi=3.14159265359d0

    integer, intent(in) :: mpol, ntheta_out
    double precision, intent(in), dimension(:) :: r_in, z_in
    double precision, intent(out), dimension(mpol) :: rcos, rsin, zcos, zsin

    logical :: reset
    integer :: jmom, jthe, iterate, jiter, nresets, ntheta
    integer, dimension(:), allocatable :: ind_sort
    double precision :: dnorm, raxis, zaxis, r_cos, r_sin, z_cos, z_sin, arg, xi, yi, &
        t1fac, r10, r10sq, delt, dtau, gmin, gnorm, gnorm_old, xm_ratio, phiangle, &
        elongate, tcon
    double precision, dimension(mpol) :: dm1, faccon, xmpq3, xmpq4, fac, ygc, ygs, &
        grc, grs, gzc, gzs
    double precision, dimension(:), allocatable :: xc, yc, xangle, r1, z1, rt1, zt1, &
        rcon, zcon, gcon, gcon2, gangle, gr, gz, pol_angle, rin, zin
    double precision, dimension(:, :), allocatable :: cosa, sina
    double precision, dimension(:), allocatable :: xvec, xstore, xdot, gvec

    ntheta = SIZE(r_in)

    allocate(ind_sort(ntheta))
    allocate(xc(ntheta), yc(ntheta), xangle(ntheta), pol_angle(ntheta))
    allocate(r1(ntheta), z1(ntheta), rt1(ntheta), zt1(ntheta), rcon(ntheta), rin(ntheta), zin(ntheta))
    allocate(zcon(ntheta), gcon(ntheta), gcon2(ntheta), gangle(ntheta), gr(ntheta), gz(ntheta))
    allocate(cosa(ntheta, mpol), sina(ntheta, mpol))
    allocate(xvec(4*mpol + ntheta), xstore(4*mpol + ntheta), xdot(4*mpol + ntheta), gvec(4*mpol + ntheta))

! Initialise, it's crucial
    rcos = 0.d0
    rsin = 0.d0
    zcos = 0.d0
    zsin = 0.d0

    dnorm = dble(2)/dble(ntheta)
    do jmom=1, mpol
        dm1(jmom) = dble(jmom-1)
    enddo
    xmpq4 = dm1**pexp
    xmpq4(1:2) = 0
    xmpq3 = SQRT(xmpq4)
    faccon = 0.125*dnorm/(1. + dm1)**pexp

    raxis = SUM(r_in)/dble(ntheta)
    zaxis = SUM(z_in)/dble(ntheta)
    xc = r_in - raxis
    yc = z_in - zaxis

    pol_angle = ATAN2(yc, xc)

    do jthe=1, ntheta
        if (pol_angle(jthe) < 0) pol_angle(jthe) = pol_angle(jthe) + 2.d0*pi
        xangle(jthe) = dble(jthe-1)
    enddo
    xangle = xangle*pi*dnorm

    call argsort(ntheta, pol_angle, ind_sort)
    call sort_by_index(ntheta, ind_sort, r_in, rin)
    call sort_by_index(ntheta, ind_sort, z_in, zin)
    xc = rin - raxis
    yc = zin - zaxis

    do iterate=1, max_iterate
        r_cos = SUM(xc*COS(xangle))
        r_sin = SUM(xc*SIN(xangle))
        z_cos = SUM(yc*COS(xangle))
        z_sin = SUM(yc*SIN(xangle))
        elongate = z_sin/r_cos
        phiangle = ATAN2(elongate*z_cos - r_sin, elongate*z_sin + r_cos)
        xangle = xangle + phiangle
        if (ABS(phiangle) < 0.02) then
            EXIT
        endif
    enddo

    do jthe=1, ntheta
        arg = xangle(jthe)
        xi = dnorm*(rin(jthe) - raxis)
        yi = dnorm*(zin(jthe) - zaxis)
        do jmom=2, mpol
            rcos(jmom) = rcos(jmom) + cos((jmom-1)*arg)*xi
            rsin(jmom) = rsin(jmom) + sin((jmom-1)*arg)*xi
            zcos(jmom) = zcos(jmom) + cos((jmom-1)*arg)*yi
            zsin(jmom) = zsin(jmom) + sin((jmom-1)*arg)*yi
        enddo
    enddo
    rcos(1) = raxis
    zcos(1) = zaxis
    rsin(2) = 0.5*(rsin(2) + zcos(2))
    zcos(2) = rsin(2)
    xvec(       1: mpol  ) = rcos
    xvec(  mpol+1: 2*mpol) = rsin
    xvec(2*mpol+1: 3*mpol) = zcos
    xvec(3*mpol+1: 4*mpol) = zsin
    xvec(4*mpol+1:       ) = xangle

    r10 = 0.5*( ABS(rcos(2)) + ABS(rsin(2)) + ABS(zcos(2)) + ABS(zsin(2)) )
    r10sq = 1./r10**2
    delt = 1.
    xm_ratio = 1./SQRT(1. + xmpq3(mpol))
    nresets = 0
    gmin = 1.d8
    reset = .True.
    xdot = 0.d0
    xstore = 0.d0
    gnorm_old = 1.d0

    iter_loop: do jiter=1, niter
        rt1  = 0.d0
        zt1  = 0.d0
        rcon = 0.d0
        zcon = 0.d0
        r1 = -rin
        z1 = -zin
        do jmom=1, mpol
            do jthe=1, ntheta
                arg = dm1(jmom) * xangle(jthe)
                cosa(jthe, jmom) = cos(arg)
                sina(jthe, jmom) = sin(arg)
                gr(jthe) = rcos(jmom)*cosa(jthe, jmom) + rsin(jmom)*sina(jthe, jmom)
                gz(jthe) = zcos(jmom)*cosa(jthe, jmom) + zsin(jmom)*sina(jthe, jmom)
                r1(jthe)   = r1(jthe)  + gr(jthe)
                z1(jthe)   = z1(jthe)  + gz(jthe)
                rt1(jthe)  = rt1(jthe) + dm1(jmom)*(rsin(jmom)*cosa(jthe, jmom) - rcos(jmom)*sina(jthe, jmom))
                zt1(jthe)  = zt1(jthe) + dm1(jmom)*(zsin(jmom)*cosa(jthe, jmom) - zcos(jmom)*sina(jthe, jmom))
                rcon(jthe) = rcon(jthe) + xmpq4(jmom)*gr(jthe)
                zcon(jthe) = zcon(jthe) + xmpq4(jmom)*gz(jthe)
            enddo
        enddo

        t1fac = 1./MAXVAL(rt1**2 + zt1**2)
        gangle = (r1*rt1 + z1 *zt1)*t1fac
        gcon = rcon*rt1 + zcon*zt1
        tcon = t1fac*xm_ratio ! scalar
        ygc = 0.d0
        ygs = 0.d0
        do jmom=2, mpol-1
            ygc(jmom) = dot_product (cosa(:, jmom), gcon(:))*faccon(jmom)
            ygs(jmom) = dot_product (sina(:, jmom), gcon(:))*faccon(jmom)
        enddo
        gcon2 = 0.d0
        do jmom=2, mpol-1
            gcon2 = gcon2 + ygc(jmom)*cosa(:, jmom) + ygs(jmom)*sina(:, jmom)
        enddo
        gcon2 = tcon*xmpq3(mpol)*gcon2
        rcon = r1 + gcon2*rt1
        zcon = z1 + gcon2*zt1
        do jmom=1, mpol
            grc(jmom) = dnorm * dot_product(cosa(:, jmom), rcon)
            grs(jmom) = dnorm * dot_product(sina(:, jmom), rcon)
            gzc(jmom) = dnorm * dot_product(cosa(:, jmom), zcon)
            gzs(jmom) = dnorm * dot_product(sina(:, jmom), zcon)
        enddo
        gzc(2) = 0.5*(gzc(2) + grs(2))
        grs(2) = gzc(2)
        grc(1) = 0.5*grc(1)
        gnorm = SUM(grc**2 + grs**2 + gzc**2 + gzs**2) * r10sq
        gnorm = gnorm + dnorm*SUM(gangle**2)

        if (jiter > 1) then
            fac = 1./(1. + tcon*xmpq3)
            grc = grc*fac
            grs = grs*fac
            gzc = gzc*fac
            gzs = gzs*fac
            gvec(       1:   mpol) = grc
            gvec(  mpol+1: 2*mpol) = grs
            gvec(2*mpol+1: 3*mpol) = gzc
            gvec(3*mpol+1: 4*mpol) = gzs
            gvec(4*mpol+1:       ) = gangle

            dtau = 0.5*min(bmax, ABS(1. - gnorm/gnorm_old) + 0.001)
            gmin = min(gmin, gnorm)
            xdot = ((1. - dtau)*xdot - delt*gvec)/(1. + dtau)
            xvec = xvec + xdot*delt

            if (gnorm/gmin > 1.e6) then
                reset = .False.
            endif
            if ((.not. reset) .or. (gmin == gnorm)) then
                if (reset) then
                    xstore = xvec
                else
                    xdot = 0.
                    xvec = xstore
                    delt = delt*0.95
                    reset = .True.
                    nresets = nresets + 1
                    if (nresets >= 100) then
                        write(6, *) 'Time step reduced 100 times without convergence', gnorm
                        rcos = 0.
                        rsin = 0.
                        zcos = 0.
                        zsin = 0.
                        return
                    endif
                endif
            endif
        endif

        rcos   = xvec(       1:   mpol)
        rsin   = xvec(  mpol+1: 2*mpol)
        zcos   = xvec(2*mpol+1: 3*mpol)
        zsin   = xvec(3*mpol+1: 4*mpol)
        xangle = xvec(4*mpol+1:)

        if (gnorm < ftol**2) EXIT   ! convergence

        gnorm_old = gnorm
        do jmom=1, mpol
            do jthe=1,ntheta
                arg = dm1(jmom) * xangle(jthe)
                cosa(jthe, jmom) = cos(arg)
                sina(jthe, jmom) = sin(arg)
            enddo
        enddo

    enddo iter_loop

    deallocate(ind_sort)
    deallocate(rin, zin, xc, yc, xangle, r1, z1, rt1, zt1, &
        rcon, zcon, gcon, gcon2, gangle, gr, gz, pol_angle)
    deallocate(cosa, sina)
    deallocate(xvec, xstore, xdot, gvec)

    end subroutine scrunch

end module scrunch2d
