module dimensions_ef_parameters

	integer, parameter :: i_dim1=300 !coil stuff
	integer, parameter :: i_dim2=300 !plasma grids
	integer, parameter :: i_dim3=300
	integer, parameter :: i_dim4=300
	integer, parameter :: i_dim5=5200
	integer, parameter :: nrho2d=3999

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

module pi_grec_vars
	use pi_vars, only: GPI,GPI2,GPI4

!	double precision, parameter :: GPI=3.1415926,GPI2=2.*GPI,GPI4=GPI2*GPI2
	double precision :: sintable(256,256)
	double precision :: costable(256,256)

end module pi_grec_vars


module    ef_circuit       ! declaration of minimal CPOs

	use errors_params
	use pi_grec_vars
	use dimensions_ef_parameters
	use rcurr_zcurr_2def
	implicit none 

!generic
	character*80 :: data_dir
	
	
	integer :: iteration_step
	
	integer :: use_limiter_yesno

	double precision :: mu0 !=0.4*GP*1.E-06

!fit
	integer :: nr_of_fit_parameters     ! nr of fit parameters from rmag to zxp_fit
	double precision :: rmag_fit   ! R mag axis
	double precision :: zmag_fit ! Z mag axis
	double precision :: k_fit ! elongation at mag axis
	double precision :: rxp_fit ! X point R
	double precision :: zxp_fit ! X point Z

!time stepping
	double precision psi_cur_old(i_dim1),dpc(i_dim1)
	double precision tau_old,tau_new



!circuits
	integer :: ncoils, nsubcoils,nreseqcoil
	integer :: nfirstwall
	integer :: nlimiter
	integer :: nblanket,nblanketpc,nelemblanketpc,nelemblanket
	integer :: nelemcoil(i_dim1),mturns(i_dim1),mequivalence(i_dim1)
	double precision :: limiterR(500),limiterZ(500)
	double precision :: lim_maxR,lim_minR,lim_maxZ,lim_minZ
	integer :: ilim_maxR,ilim_minR,ilim_maxZ,ilim_minZ
	double precision :: Rcoil(i_dim1),Zcoil(i_dim1),drcoil(i_dim1),dzcoil(i_dim1)
	double precision :: anglecoil(i_dim1),rescoil(i_dim1,i_dim1),indcoil(i_dim1) !self induct
	double precision :: curcoil(i_dim1)

!if areas are 0 --> treated as filaments, otherwise they are rectangle areas
	double precision :: Rwall(i_dim1),Zwall(i_dim1),reswall,areawall(i_dim1),curwall(i_dim1), & 
	& indwall(300)	
	double precision :: Rblan(i_dim1),Zblan(i_dim1),resblan,widthblan
	double precision :: areablan(i_dim1),curblan(i_dim1), &
	& indblan(300)	
	double precision :: Rblanpc(i_dim1),Zblanpc(i_dim1),resblanpc(i_dim1),areablanpc(i_dim1),curblanpc(i_dim1), &
	& indblanpc(i_dim1)

	integer :: nconduc,nblocks,npassive,nactive
	double precision :: curconduc(i_dim1),resconduc(i_dim1,i_dim1),indconduc(i_dim1,i_dim1),voltage(i_dim1),voltage_old(i_dim1) !self and mutual induc
	double precision :: psiplasmatoconduc(i_dim1) !plasma --> conduc at t
	double precision ::  psiconductoplasma !conduc --> plasma boundary
	double precision :: cur_con_old(i_dim1)
	double precision :: r_cond(i_dim1),z_cond(i_dim1)

! coordinates:
! r,z --> rectangular grid in meters
! rho --> distance of point from magnetic axis in meters
! teta --> angle from LFS midplane which is 0, coincides with Zmag plane
! psigrid --> in poloidal flux equispaced


!grids
	integer :: nr, nz,nrho,nteta,nr2,nz2,nr1,nz1   !nteta+1 is the periodic point. nrho is the plasma boundary

	double precision :: rmin,rmax,zmin,zmax ! rho is defined as rho_toroidal as in astra
	
	double precision :: r(i_dim2),z(i_dim2),rho(i_dim2,i_dim2),teta(i_dim2) ! rho is defined as actual distance in meters as in astra
	double precision :: rcomp(i_dim2),zcomp(i_dim2) ! rho is defined as distance
	double precision :: dr,dz,drho(i_dim2,i_dim2),dteta,psigrid(i_dim2) ! teta and psigrid are equispaced
	double precision :: rpol(i_dim2,i_dim2),zpol(i_dim2,i_dim2) ! R,Z in polar coordinates half radial grid
	double precision :: rpul(i_dim2,i_dim2),zpul(i_dim2,i_dim2) ! R,Z in polar coordinates full radial grid
	double precision :: psia_2d(nrho2d),ffp_2d(nrho2d),ppp_2d(nrho2d),ipol_2d(i_dim2),pres_2d(i_dim2)
	double precision :: psia_1d(i_dim2),ffp_1d(i_dim2),ppp_1d(i_dim2)
	double precision :: area_eff(i_dim2,i_dim2) !

!inversion matrix for R,Z solution 
	double precision :: invMM_gs2d(i_dim2,i_dim2,i_dim2)	
	
! potential
	double precision :: psirz(i_dim2,i_dim2),psirhoteta(i_dim2,i_dim2)
	double precision :: u_n(i_dim2,i_dim2)
	double precision :: omega_pl(i_dim2,i_dim2)
	double precision :: psiextrz(i_dim2,i_dim2),psiplasrz(i_dim2,i_dim2)
	double precision :: phirhoteta(i_dim2,i_dim2) !toroidal flux
	double precision :: zlimpotential(i_dim2,i_dim2)
	double precision :: 	derivpsi(8),zbot,ztop,raus,rinner

! boundary and axis FBE
	integer :: nbnd,ngbnd,redo_bnd !redo_bnd is temporary
	double precision :: rbnd(i_dim5),zbnd(i_dim5),psibnd,psiaxis
	integer ibnd(i_dim5),jbnd(i_dim5)
	double precision :: rax,zax,alpsep
	integer :: i_plasmatype !(0-limited, 1-single null, 2-double null)
	integer :: iaxis,jaxis
	integer :: n_of_xpoints
	double precision :: r_xpoint(500),z_xpoint(500),deriv_x(5,500)
	integer :: max_xpoints=500
	double precision :: trax,tzax !true axis for more precision
	double precision :: psistabR,psistabZ !stab terms
	double precision :: green_bnd_f(16*i_dim2**2)
	double precision :: dr_factor_init, dz_factor_init

! boundary and axis PBE
	double precision :: rbndp(i_dim2),zbndp(i_dim2),psibndp,psiaxisp
	double precision :: raxp,zaxp
	double precision :: rexp(i_dim2),zexp(i_dim2),tetaexp(i_dim2)

	
	
	
! plasma parameters
	double precision :: iplasma,btor0,rgeom0,psplex,li3,betapol
	double precision :: pprime(i_dim2),ffprime(i_dim2),pressure(i_dim2),psigrida(i_dim2)	 !these 3 come from astra, psi is FP of astra
	double precision :: jrz(i_dim2,i_dim2),jrhoteta(i_dim2,i_dim2),ipol(i_dim2)	 !current density
	
	
end module ef_circuit


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
	double precision :: ffprimebez(i_dim2)
	double precision :: pprimebez(i_dim2)
	double precision :: pressbez(i_dim2)
	double precision :: ipolbez(i_dim2)

end module transfer_functions






module metric_coefficients_pbe
	use dimensions_ef_parameters

! metric coefficients in polar coordinates 
	double precision :: dl(i_dim2,i_dim2),dator(i_dim2,i_dim2),dalat(i_dim2,i_dim2),dvol(i_dim2,i_dim2)
	double precision :: vol(i_dim2),sator(i_dim2)
	double precision :: gg1(i_dim2),gg2(i_dim2),gg3(i_dim2) !metrics
	double precision :: bpoloidal(i_dim2,i_dim2),bphi(i_dim2,i_dim2)	 !fields
	double precision :: lambda2d(i_dim2,i_dim2)
	double precision :: R_curr_0D,Z_curr_0D
	
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





module fft_mod_eff
	use pi_grec_vars
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


