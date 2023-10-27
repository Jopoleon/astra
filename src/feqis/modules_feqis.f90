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
