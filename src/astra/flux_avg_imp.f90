module flux_avg_imp

implicit none

contains

!---------------------------------------------------------------------
    subroutine lfs2fsa_impDV(geom_type, Zimp_in, Aimp_in, e0imp_out, FVimp_out)
!---------------------------------------------------------------------
! Calculate the terms that transform the low field side (LFS) diffusive and convective
! transport coefficients of an impurity species into flux surface averaged (FSA) coefficients
! following the appendix of Angioni 2014 Nucl. Fusion 54 083028
!
! --- D. Fajardo, Dec 2025
!
! * INPUTS
! --------
! - geom_type --> 0: circular, 1: Miller (STILL TO DO), 2: full flux surface
! - Zimp_in ----> impurity charge (profile) [-]
! - Aimp_in ----> impurity mass [-]
!
! * OUTPUTS
! ---------
! - e0imp_out -> transformation coefficient of the diffusion [-]
! - FVimp_out -> transformation coefficient of the convection [1/m]
!---------------------------------------------------------------------

    use scalars, only: NA1, RTOR, ZMJ, AMJ
    use status, only: TE, TI, ZEF, VTOR, NRD

    integer, intent(in) :: geom_type
    double precision, intent(in) :: Aimp_in
    double precision, dimension(NRD), intent(in) :: Zimp_in
    double precision, dimension(NRD), intent(out) :: e0imp_out, FVimp_out

! physical constants
    integer, parameter :: n_theta=101 !301
    double precision, parameter :: q_e=1.602176634e-19   ! electron charge [C]
    double precision, parameter :: m_p=1.67262192369e-27 ! proton mass [kg]

    integer :: i
    double precision, dimension(NA1) :: Z_imp, T_e, T_i, Z_eff, Mach_i_sq, Mz_star_sq, &
         r_min, R_LFS, e0imp, FVimp, d_e0imp_dr
    double precision, dimension(n_theta) :: theta
    double precision, dimension(NA1, n_theta) :: R_2D, Z_2D, Jacobian, exp_m_Eimp

! Plasma profiles
    T_e = 1.0e3*max(TE(1:NA1), 1.e-5)                                ! electron temperature in eV
    T_i = 1.0e3*max(TI(1:NA1),1.e-5)                                 ! main ion temperature in eV  
    Z_eff     = ZEF(1:NA1)                                           ! effective charge
    Z_imp     = Zimp_in(1:NA1)                                       ! impurity charge
    Mach_i_sq = abs(VTOR(1:NA1))**2/(2.*q_e*T_i/(AMJ*m_p))           ! main ion Mach number squared
    Mz_star_sq = Mach_i_sq*(Aimp_in/AMJ - (Z_imp/ZMJ)*Z_eff/(Z_eff + T_i/T_e)) ! effective impurity Mach number squared
 
! Magnetic geometry  
    call flux_surf_geom(geom_type, n_theta, r_min, theta, R_2D, Z_2D, Jacobian, R_LFS)

! Normalized energy exponential
    do i=1, NA1
        exp_m_Eimp(i, :)  = exp(Mz_star_sq(i)*((R_2D(i, :)**2 - R_LFS(i)**2)/RTOR**2))
    enddo

! Surface averaging
    do i=1, NA1
        call flux_surf_avg(n_theta, theta , exp_m_Eimp(i, :), Jacobian(i, :), e0imp(i))
    enddo

    do i=2, NA1-1
        d_e0imp_dr(i) = (e0imp(i+1)-e0imp(i-1))/(r_min(i+1)-r_min(i-1))
    enddo

    d_e0imp_dr(1)  = (e0imp(2)-e0imp(1))/(r_min(2)-r_min(1))
    d_e0imp_dr(NA1) = (e0imp(NA1)-e0imp(NA1-1))/(r_min(NA1)-r_min(NA1-1))
    FVimp = d_e0imp_dr/e0imp

! Outputs
    e0imp_out(1: NA1) = e0imp(1: NA1)
    FVimp_out(1: NA1) = FVimp(1: NA1)
  
    end subroutine lfs2fsa_impDV

!---------------------------------------------------------------------
    subroutine flux_surf_avg(nth, thetay, AF, JJ, Aavg)
!---------------------------------------------------------------------
! Flux surface average of a function AF
!
! * INPUTS
! --------
! - nth ----> size of theta-dependent arrays [-] {int} 
! - thetay -> field-line angle [-] {arr, nth}
! - AF -----> function to average {arr, nth}
! - JJ -----> jacobian of coordinate system
!
! * OUTPUTS
! ---------
! - Aavg -----> flux surface average of AF function {float}
!-----------------------------------------------------------------------------------

    double precision, dimension(nth), intent(in) :: AF, JJ, thetay
    double precision, intent(out) :: Aavg
    integer :: nth, ith
    double precision :: denom

    Aavg = 0.
    denom = 0.

    do ith=2,nth
        denom = denom + 0.5*(thetay(ith)-thetay(ith-1))*(JJ(ith) + JJ(ith-1))
    enddo

    do ith=2,nth
        Aavg = Aavg+ 0.5*(thetay(ith)-thetay(ith-1))*(AF(ith)*JJ(ith) + AF(ith-1)*JJ(ith-1))
    enddo

    if (abs(denom) > 0.) then
        Aavg = Aavg/denom
    else
        Aavg = sum(AF)/nth
    endif

    end subroutine flux_surf_avg

!---------------------------------------------------------------------
    subroutine flux_surf_geom(geom_type, ntheta_in, rmin_out, theta_out, R_out, Z_out, Jacobian_out, R_LFS_out)
!---------------------------------------------------------------------
! Calculate flux surface contours R(r,theta), Z(r,theta) and interpolate them to minor
! radius and poloidal grids. The poloidal grid goes from -pi to pi, with theta = 0 at the LFS
! The Jacobian of (R,Z) -> (r,theta) and the low field side major radius are also calculated
!
! --- D. Fajardo, May 2023, updated December 2023
!
! * INPUTS
! --------
! - geom_type ----> 0: circular, 1: Miller (STILL TO DO), 2: full flux surface
! - ntheta_in ----> discretization of the poloidal grid [-]
!
! * OUTPUTS
! ---------
! - rmin_out -----> minor radius coordinate [m]
! - theta_out ----> poloidal coordinate [-]
! - R_out --------> major radius coordinate [m]
! - Z_out --------> vertical coordinate [m]
! - Jacobian_out -> jacobian [m]
! - R_LFS_out ----> low field side major radius [m]
!---------------------------------------------------------------------

    use scalars, only: NA1, RTOR, time, tau, tstart
    use status, only: AMETR, SHIF, SHIV
    use parameters_a2equil, only: equil_now
    use numerical_tools, only: qinterp

    integer, intent(in) :: geom_type, ntheta_in
    double precision, dimension(NA1), intent(out) :: rmin_out, R_LFS_out
    double precision, dimension(ntheta_in), intent(out) :: theta_out
    double precision, dimension(NA1, ntheta_in), intent(out) :: R_out, Z_out, Jacobian_out

! mathematical constants
    double precision, parameter :: pi=4.d0*ATAN(1.d0)

    integer :: i, j
    integer :: idxmpi
    integer :: nrho_surf, nthe_surf

! intermediate variables
    double precision, allocatable, dimension(:) :: rmin_equ, th0, th1, th2, pf_eq, rho_eq
    double precision, allocatable, dimension(:, :) :: R_1, R_2, R_3, Z_1, Z_2, Z_3
    double precision, dimension(NA1, ntheta_in) :: dRdr, dRdth, dZdr, dZdth, grr, grt, gtt
 
!-----------------------------------------------------------------------------------------------!
! magnetic geometry

    rmin_out = AMETR(1:NA1)   ! minor radius [m]

! flux surface contours from equilibrium
    nrho_surf = SIZE(equil_now%coord_sys%position%r, 1)
    nthe_surf = SIZE(equil_now%coord_sys%position%r, 2)
    allocate(pf_eq(nrho_surf), rho_eq(nrho_surf))

! allocate other quantities
    allocate(rmin_equ(nrho_surf))
    allocate(th0(nthe_surf), th1(nthe_surf))
    allocate(R_1(nrho_surf, nthe_surf), Z_1(nrho_surf, nthe_surf))
    allocate(th2(nthe_surf + 2))
    allocate(R_2(nrho_surf, nthe_surf+2), Z_2(nrho_surf, nthe_surf+2))
    allocate(R_3(nrho_surf, ntheta_in), Z_3(nrho_surf, ntheta_in))

! initial poloidal coordinate
    do j=1, nthe_surf
        th0(j) = atan2(equil_now%coord_sys%position%z(nrho_surf, j) - SHIV(NA1), equil_now%coord_sys%position%r(nrho_surf, j)-(RTOR+SHIF(NA1)))
    enddo

! find location where th0 is closest to -pi
    idxmpi = minloc(abs(th0-(-pi)), 1)

! roll coordinates over this 
    th1(1:nthe_surf-idxmpi+1) = th0(idxmpi:)
    th1(nthe_surf-idxmpi+2:)  = th0(1:idxmpi-1)

! now that th1 is correct, roll R_1, Z_1
    if (time-tstart > tau) then
        R_1(:, 1: nthe_surf-idxmpi+1) = equil_now%coord_sys%position%r(:, idxmpi:)
        R_1(:, nthe_surf-idxmpi+2:)   = equil_now%coord_sys%position%r(:, 1: idxmpi-1)
        Z_1(:, 1: nthe_surf-idxmpi+1) = equil_now%coord_sys%position%z(:, idxmpi:)
        Z_1(:, nthe_surf-idxmpi+2:)   = equil_now%coord_sys%position%z(:, 1: idxmpi-1)
    endif
! finish the poloidal turn:
    th2(2:nthe_surf+1) = th1
    th2(1) = -pi
    th2(nthe_surf+2) = pi

    R_2(:, 2: nthe_surf+1) = R_1
    Z_2(:, 2: nthe_surf+1) = Z_1

    R_2(:, 1) = 0.5*(R_1(:, 1) + R_1(:, nthe_surf))
    R_2(:, nthe_surf+2) = 0.5*(R_1(:, 1) + R_1(:, nthe_surf))

    Z_2(:, 1) = 0.5*(Z_1(:, 1) + Z_1(:, nthe_surf))
    Z_2(:, nthe_surf+2) = 0.5*(Z_1(:, 1) + Z_1(:, nthe_surf))
    
! interpolate

    do j=1, ntheta_in-1
        theta_out(j) = -pi + 2.*pi*float(j-1)/float(ntheta_in-1)
    enddo
    theta_out(ntheta_in) = pi

    do i=1, nrho_surf
        rmin_equ(i) = max(0.5*(maxval(equil_now%coord_sys%position%r(i, :)) - minval(equil_now%coord_sys%position%r(i, :))), 1.e-3)
    enddo

    if (geom_type == 2) then ! full flux surface geometry, interpolate in theta
        do i=1, nrho_surf
            call qinterp(th2, R_2(i, :), nthe_surf + 2, theta_out, R_3(i, :), ntheta_in)
            call qinterp(th2, Z_2(i, :), nthe_surf + 2, theta_out, Z_3(i, :), ntheta_in)
        enddo

! interpolate in r
        do j=1, ntheta_in
            call qinterp(rmin_equ, R_3(:, j), nrho_surf, rmin_out, R_out(:, j), NA1)
            call qinterp(rmin_equ, Z_3(:, j), nrho_surf, rmin_out, Z_out(:, j), NA1)
        enddo
    elseif (geom_type == 0) then ! circular geometry
        do j=1, ntheta_in
            R_out(:, j) = RTOR + SHIF(1:NA1) + rmin_out*cos(theta_out(j))
            Z_out(:, j) = SHIV(1:NA1) + rmin_out*sin(theta_out(j))
        enddo 
    elseif (geom_type == 1) then
        print *, "Miller still to do"
    else
        print *, "select a geom_type = 0, 1 or 2"
    endif

! Jacobian

! dR/dr
    do i=2, NA1-1
        dRdr(i, :) = (R_out(i+1, :) - R_out(i-1, :))/(rmin_out(i+1) - rmin_out(i-1))
    enddo
    dRdr(  1, :) = (R_out(  2, :) - R_out(    1, :))/(rmin_out(  2) - rmin_out(1))
    dRdr(NA1, :) = (R_out(NA1, :) - R_out(NA1-1, :))/(rmin_out(NA1) - rmin_out(NA1-1))

! dZ/dr
    do i=2, NA1-1
        dZdr(i, :) = (Z_out(i+1, :) - Z_out(i-1, :))/(rmin_out(i+1) - rmin_out(i-1))
    enddo
    dZdr(  1, :) = (Z_out(  2, :) - Z_out(    1, :))/(rmin_out(  2) - rmin_out(1))
    dZdr(NA1, :) = (Z_out(NA1, :) - Z_out(NA1-1, :))/(rmin_out(NA1) - rmin_out(NA1-1))

! dR/dtheta
    do j=2, ntheta_in-1
        dRdth(:, j) = (R_out(:, j+1) - R_out(:, j-1))/(theta_out(j+1) - theta_out(j-1))
    enddo

    dRdth(:,       1)   = (R_out(:,         2) - R_out(:,           1))/(theta_out(        2) - theta_out(1))
    dRdth(:, ntheta_in) = (R_out(:, ntheta_in) - R_out(:, ntheta_in-1))/(theta_out(ntheta_in) - theta_out(ntheta_in-1))

! dZ/dtheta
    do j=2, ntheta_in-1
        dZdth(:, j) = (Z_out(:, j+1) - Z_out(:, j-1))/(theta_out(j+1) - theta_out(j-1))
    enddo

    dZdth(:,         1) = (Z_out(:,         2) - Z_out(:,           1))/(theta_out(        2) - theta_out(1))
    dZdth(:, ntheta_in) = (Z_out(:, ntheta_in) - Z_out(:, ntheta_in-1))/(theta_out(ntheta_in) - theta_out(ntheta_in-1)) 


    grr = dRdr**2 + dZdr**2
    grt = dRdr*dRdth + dZdr*dZdth
    gtt = dRdth**2 + dZdth**2
    Jacobian_out = R_out*sqrt(grr*gtt - grt**2)

! finally, LFS major radius

    R_LFS_out = RTOR + SHIF(1:NA1) + rmin_out

    deallocate(pf_eq, rho_eq, rmin_equ, th0, th1, R_1, Z_1, th2, R_2, Z_2, R_3, Z_3)

    end subroutine flux_surf_geom

end module flux_avg_imp
