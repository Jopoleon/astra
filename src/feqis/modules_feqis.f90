module feqis_dimensions

    integer, parameter :: i_dim1=300 !coil stuff
    integer, parameter :: i_dim2=300 !plasma grids
    integer, parameter :: i_dim3=300
    integer, parameter :: i_dim4=300
    integer, parameter :: i_dim5=5200
    integer, parameter :: nrho2d=3999

end module feqis_dimensions

!---------------------------------------------------------------------
module errors_params

    double precision :: err_circ_plasma_iter, err_find_oxpoints, &
        err_find_oxpoints_derivs, err_find_psistab, err_find_delr, &
        err_find_biquad, err_epsilon, err_gaptolez, err_fix_boundary

end module errors_params

!---------------------------------------------------------------------
module rcurr_zcurr_2def

    double precision :: R_curr_2D, Z_curr_2D

end module rcurr_zcurr_2def

!---------------------------------------------------------------------
module feqis_circuit

    use feqis_dimensions, only: i_dim1, i_dim2, i_dim5, nrho2d

    implicit none 

    integer :: use_limiter_yesno

!time stepping
    double precision :: tau_old, tau_new
    double precision, dimension(i_dim1) :: psi_cur_old, dpc

!circuits
    integer :: ncoils, nreseqcoil, nlimiter
    integer, dimension(i_dim1) :: mequivalence
    double precision, dimension(500) :: limiterR, limiterZ
    double precision :: lim_maxR, lim_minR, lim_maxZ, lim_minZ
    integer :: ilim_maxR, ilim_minR, ilim_maxZ, ilim_minZ
    double precision, dimension(i_dim1) :: Rcoil, Zcoil, drcoil, dzcoil, &
        anglecoil

    integer :: nconduc, nblocks, npassive, nactive
    double precision, dimension(i_dim1) :: curconduc, voltage, voltage_old, &
        cur_con_old, r_cond, z_cond
    double precision, dimension(i_dim1, i_dim1) :: resconduc, indconduc
    double precision :: psiplasmatoconduc(i_dim1) !plasma --> conduc at t

! coordinates:
! r, z --> rectangular grid in meters
! rho --> distance of point from magnetic axis in meters
! teta --> angle from LFS midplane which is 0, coincides with Zmag plane
! psigrid --> in poloidal flux equispaced

!grids
    integer :: nr, nz, nrho, nteta, nr2, nz2, nr1, nz1, nbnd, ngbnd, &
        redo_bnd

    double precision :: rmin, rmax, zmin, zmax, dr, dz, dteta, &
        zbot, ztop, raus, rinner
    double precision, dimension(i_dim2) :: r, z, teta, rcomp, zcomp, psigrid  
    double precision, dimension(i_dim2, i_dim2) :: rho, area_eff, &
        rpol, zpol, rpul, zpul, u_n, omega_pl, zlimpotential, &
        psirz, psirhoteta, psiextrz, psiplasrz
! r(z)pol: R, Z in polar coordinates half radial grid
! r(z)pul: R, Z in polar coordinates full radial grid
    double precision, dimension(nrho2d) :: psia_2d, ffp_2d, ppp_2d
    double precision, dimension(i_dim2) :: psia_1d, ffp_1d, ppp_1d
    double precision, dimension(8) :: derivpsi

! boundary and axis FBE, PBE
    integer, parameter :: max_xpoints=500
    integer :: i_plasmatype !(0-limited, 1-single null, 2-double null)
    integer :: iaxis, jaxis, n_of_xpoints
    double precision :: psibnd, psiaxis, psibndp, psiaxisp, &
        rax, zax, trax, tzax, raxp, zaxp, &
        alpsep, psistabR, psistabZ, dr_factor_init, dz_factor_init
    double precision, dimension(i_dim5) :: rbnd, zbnd
    double precision, dimension(max_xpoints) :: r_xpoint, z_xpoint
    double precision :: deriv_x(5, max_xpoints)
    double precision :: green_bnd_f(16*i_dim2**2)
    double precision, dimension(i_dim2) :: rbndp, zbndp, rexp, zexp, tetaexp

! plasma parameters
    double precision :: iplasma, btor0, rgeom0, psplex, li3, betapol
    double precision, dimension(i_dim2) :: pprime, ffprime, pressure, psigrida, ipol
    double precision, dimension(i_dim2, i_dim2) :: jrz, jrhoteta

end module feqis_circuit

!---------------------------------------------------------------------
module transfer_functions

    use feqis_dimensions, only: i_dim2

    double precision, dimension(i_dim2) ::tetabez, dpsidvbez, psibez, &
        g2bez, g2ibez, gm1bez, routbez, rinbez, vbez, g1bez, gm41bez, &
        ggrhobez, bmaxbez, bminbez, gm4bez, bdb0bez, gm5bez, fofbbez, &
        areatbez, perimbez, shifbez, kbez, surfbez, triaubez, phibez, &
        qbez, t2dbez, rbp2_b2bez, &
        ffprimebez, pprimebez, pressbez, ipolbez

    double precision, dimension(i_dim2, i_dim2) :: rpbez, zpbez, &
        rminbez, bpcellbez, bcellbez, rmin2dbez, &
        bpcell2dbez, bcell2dbez

end module transfer_functions


!---------------------------------------------------------------------
module metric_coefficients_pbe

    use feqis_dimensions, only: i_dim2

! metric coefficients in polar coordinates 
    double precision :: R_curr_0D, Z_curr_0D
    double precision, dimension(i_dim2) :: vol, sator, gg1, gg2, gg3
    double precision, dimension(i_dim2, i_dim2) :: dl, dator, dalat, dvol, &
        bpoloidal, bphi, lambda2d
    
end module metric_coefficients_pbe

!---------------------------------------------------------------------
module green_matrix

    use feqis_dimensions, only: i_dim1, i_dim2

    double precision, dimension(i_dim2, i_dim2, i_dim1) :: greeni, &
        dgreenirpl, dgreenizpl
    double precision, dimension(i_dim1, i_dim1) :: dgreenirj, dgreenizj

end module green_matrix

!---------------------------------------------------------------------
module fft_mod_eff

    use pi_vars, only: GPI2

    implicit none

    integer, parameter :: dp=selected_real_kind(15, 300)
    double precision, dimension(256, 256) :: sintable, costable

    contains
 
  ! In place Cooley-Tukey FFT
        recursive subroutine fft_eff(x)

        complex(kind=dp), dimension(:), intent(inout)  :: x
        complex(kind=dp) :: t
        integer :: N, i
        complex(kind=dp), dimension(:), allocatable :: even, odd
 
        N = size(x)
 
        if(N .le. 1) return
 
        allocate(odd((N+1)/2))
        allocate(even(N/2))

! divide
        odd  = x(1:N:2)
        even = x(2:N:2)
 
! conquer
        call fft_eff(odd)
        call fft_eff(even)
 
! combine
        do i=1, N/2
            t = exp(cmplx(0.0, -GPI2*(i - 1.)/(N + 0.)))*even(i)
            x(i)     = odd(i) + t
            x(i+N/2) = odd(i) - t
        enddo
 
        deallocate(odd)
        deallocate(even)
 
        end subroutine fft_eff
 
end module fft_mod_eff
