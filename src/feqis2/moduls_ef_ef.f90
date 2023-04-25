module dimensions_ef_parameters

	integer, parameter :: i_dim1=300 !coil stuff
	integer, parameter :: i_dim2=300 !plasma grids
	integer, parameter :: i_dim3=300
	integer, parameter :: i_dim4=300
	integer, parameter :: i_dim5=5200

end module dimensions_ef_parameters


module errors_params
	double precision :: err_circ_plasma_iter
	double precision :: err_find_oxpoints
	double precision :: err_find_oxpoints_derivs
	double precision :: err_find_psistab
	double precision :: err_find_delr
	double precision :: err_find_biquad
	double precision :: err_epsilon
	double precision :: err_gaptolez
	double precision :: err_fix_boundary
end module errors_params


module rcurr_zcurr_2def

	double precision :: R_curr_2D,Z_curr_2D

end module rcurr_zcurr_2def

module pi_vars

	double precision, parameter :: GPI=3.1415926,GPI2=2.*GPI,GPI4=GPI2*GPI2,mu0=0.4*GPI
	double precision :: sintable(256,256)
	double precision :: costable(256,256)

end module pi_vars


module transfer_functions
	use dimensions_ef_parameters

	double precision :: rpbez(i_dim2,i_dim2)
	double precision :: zpbez(i_dim2,i_dim2)
	double precision :: tetabez(i_dim2)
	double precision :: rminbez(i_dim2,i_dim2)
	double precision :: bpcellbez(i_dim2,i_dim2)
	double precision :: bcellbez(i_dim2,i_dim2)
	double precision :: dpsidvbez(i_dim2)
	double precision :: psibez(i_dim2)
	double precision :: g2bez(i_dim2)
	double precision :: g2ibez(i_dim2)
	double precision :: gm1bez(i_dim2)
	double precision :: routbez(i_dim2)
	double precision :: rinbez(i_dim2)
	double precision :: vbez(i_dim2)
	double precision :: g1bez(i_dim2)
	double precision :: gm41bez(i_dim2)
	double precision :: ggrhobez(i_dim2)
	double precision :: bmaxbez(i_dim2)
	double precision :: bminbez(i_dim2)
	double precision :: gm4bez(i_dim2)
	double precision :: bdb0bez(i_dim2)
	double precision :: gm5bez(i_dim2)
	double precision :: fofbbez(i_dim2)
	double precision :: areatbez(i_dim2)
	double precision :: perimbez(i_dim2)
	double precision :: shifbez(i_dim2)
	double precision :: kbez(i_dim2)
	double precision :: surfbez(i_dim2)
	double precision :: triaubez(i_dim2)
	double precision :: phibez(i_dim2)
	double precision :: qbez(i_dim2)
	double precision :: t2dbez(i_dim2)
	double precision :: rmin2dbez(i_dim2,i_dim2)
	double precision :: bpcell2dbez(i_dim2,i_dim2)
	double precision :: bcell2dbez(i_dim2,i_dim2)
	double precision :: rbp2_b2bez(i_dim2)

end module transfer_functions






module metric_coefficients_pbe
	use dimensions_ef_parameters

! metric coefficients in polar coordinates 
	double precision :: dl(i_dim2,i_dim2),dator(i_dim2,i_dim2),dalat(i_dim2,i_dim2),dvol(i_dim2,i_dim2)
	double precision :: vol(i_dim2),sator(i_dim2)
	double precision :: gg1(i_dim2),gg2(i_dim2),gg3(i_dim2) !metrics
	double precision :: bpoloidal(i_dim2,i_dim2),bphi(i_dim2,i_dim2)	 !fields
	double precision :: lambda2d(i_dim2,i_dim2)
	double precision :: lambda2dp(i_dim2,i_dim2),R_curr_0D,Z_curr_0D
	
end module metric_coefficients_pbe









module    green_matrix       ! declaration of minimal CPOs
	use dimensions_ef_parameters

	double precision greeni(i_dim2,i_dim2,i_dim1)
	double precision dgreenirpl(i_dim2,i_dim2,i_dim1)
	double precision dgreenizpl(i_dim2,i_dim2,i_dim1)
	double precision dgreenirj(i_dim1,i_dim1)
	double precision dgreenizj(i_dim1,i_dim1)
!	double precision dgreenir(i_dim2,i_dim2,i_dim1)
!	double precision dgreeniz(i_dim2,i_dim2,i_dim1)
!	double precision dgreenirr(i_dim2,i_dim2,i_dim1)
!	double precision dgreenizz(i_dim2,i_dim2,i_dim1)

end module green_matrix




module    exchange_with_astra       ! declaration of minimal CPOs
	use dimensions_ef_parameters
	
	double precision :: tau_circuit_ef	
	double precision :: tau_gseq_ef
	double precision :: time_astra

	double precision :: activate_coil_ef(i_dim1)
	double precision :: sign_coil(i_dim1)
	double precision :: current_limit_ef(i_dim1,2) ! 1 is upper, 2 is lower
	double precision :: force_coil(i_dim1,i_dim1) ! where it is 1, forces coil i,i to current of i,j


	integer :: use_limiter_astra  ! 1-uses limiter, 0-ignore limiter

	integer :: refit_mode   ! if -1 - 1 pass only , 0 - self-consistent solution, if 1 - stab axis using passive wall currents fourier modes cos and sin, if 2 - same as 1 but uses boundary points using 5 fourier modes, 3-uses full currents fit using analytic F function and fit file efonfit.dat
	integer :: solve_fix   ! if 0 - solve full fix boundary problem, if 1 - 1 iteration only , 2 - only contouring
	integer :: execute_plasma   ! if 0 - only circuit equations, if 1 - solve plasma gseq too 
	integer :: nonegcurr ! nonegcurr = 0 --> no negative current allowed in plasma
	integer :: fast_mode   ! 0 - normale, 1 - domnt do iterationsin fbe gse
	integer :: reconnect_circuits   ! 0-nothing, 1-recompute matrix with new circuits
	integer :: new_equivalence(i_dim1,10)   ! if reconnect, says what is the new_equivalence, for example (1,1,1,0,0,0,0,..) means coil 1,2,3 become 1,1,1
	integer :: n_equivalence   ! number of equivalences
	integer :: use_reduce_circuit   ! 0-all coils solved. 1 - some coils not solved
	integer :: n_of_newton_iterations   ! to find actual mag axis. recommended between 5 - 10 
	double precision :: raxis_astra,zaxis_astra,psi0_astra,psib_astra
	integer :: n_fourier_restab_boundary   ! nr of fourier modes for boundary restab, default = 5
	integer :: psplex_from_fbe ! put 1 to get psplex fromfree boundary
!	double precision :: rhoedge,qedge
	double precision :: dr_factor_init_astra, dz_factor_init_astra !factors of dr and dz for initial iterations

end module exchange_with_astra




module fft_mod_eff
	use pi_vars
  implicit none
  integer,       parameter :: dp=selected_real_kind(15,300)
contains
 
  ! In place Cooley-Tukey FFT
  recursive subroutine fft_eff(x)
    complex(kind=dp), dimension(:), intent(inout)  :: x
    complex(kind=dp)                               :: t
    integer                                        :: N
    integer                                        :: i
    complex(kind=dp), dimension(:), allocatable    :: even, odd
 
    N=size(x)
 
    if(N .le. 1) return
 
    allocate(odd((N+1)/2))
    allocate(even(N/2))
 
    ! divide
    odd =x(1:N:2)
    even=x(2:N:2)
 
    ! conquer
    call fft_eff(odd)
    call fft_eff(even)
 
    ! combine
    do i=1,N/2
       t=exp(cmplx(0.0,-2.0*GPI*(i-1.)/(N+0.)))*even(i)
       x(i)     = odd(i) + t
       x(i+N/2) = odd(i) - t
    end do
 
    deallocate(odd)
    deallocate(even)
 
  end subroutine fft_eff
 
end module fft_mod_eff


