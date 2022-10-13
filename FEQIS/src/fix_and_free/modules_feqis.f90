module dimensions_ef_parameters

implicit none

integer, parameter :: i_dim1=300, i_dim2=300, i_dim5=5400

end module dimensions_ef_parameters

!-------------------------------------------
module pi_grec_vars

implicit none

double precision, parameter :: GPI=3.141592653589793, &
    GPI2=2.*GPI, GPI4=GPI2**2.0, mu0=0.4*GPI, muvac=4.e-7*GPI

end module pi_grec_vars

!-------------------------------------------
module ef_circuit       ! declaration of minimal CPOs

use dimensions_ef_parameters, only: i_dim1, i_dim2, i_dim5

implicit none 

!generic
character*120 :: data_dir
integer :: data_dir_k, iteration_step, max_iter, nr_of_fit_parameters

double precision :: errgs2dfull,  tau_old, tau_new, &
    rmag_fit, zmag_fit, k_fit, & ! R, z, elongation at mag axis
    rxp_fit, zxp_fit ! X point R, Z

!time stepping
double precision psi_cur_old(i_dim1), dpc(i_dim1)

!circuits
integer :: ncoils, nsubcoils, nreseqcoil
integer :: nfirstwall
integer :: nlimiter
integer :: nblanket, nblanketpc, nelemblanketpc, nelemblanket
integer :: nelemcoil(i_dim1), mturns(i_dim1), mequivalence(i_dim1)
double precision :: limiterR(500), limiterZ(500)
double precision :: Rcoil(i_dim1), Zcoil(i_dim1), drcoil(i_dim1), dzcoil(i_dim1)
double precision :: anglecoil(i_dim1), indcoil(i_dim1) !self induct
double precision :: curcoil(i_dim1)

!if areas are 0 --> treated as filaments, otherwise they are rectangle areas
double precision :: Rwall(i_dim1), Zwall(i_dim1), reswall, areawall(i_dim1), curwall(i_dim1)    
double precision :: Rblan(i_dim1), Zblan(i_dim1), resblan, areablan(i_dim1), curblan(i_dim1)    
double precision :: Rblanpc(i_dim1), Zblanpc(i_dim1), resblanpc(i_dim1), areablanpc(i_dim1), curblanpc(i_dim1)

integer :: nconduc, npassive, nactive
double precision :: curconduc(i_dim1), resconduc(i_dim1, i_dim1), indconduc(i_dim1, i_dim1), voltage(i_dim1) !self and mutual induc
double precision :: psiplasmatoconduc(i_dim1) !plasma --> conduc at t
double precision ::  psiconductoplasma !conduc --> plasma boundary
double precision :: cur_con_old(i_dim1)
double precision :: r_cond(i_dim1), z_cond(i_dim1)

! coordinates:
! r, z --> rectangular grid in meters
! rho --> distance of point from magnetic axis in meters
! teta --> angle from LFS midplane which is 0, coincides with Zmag plane
! psigrid --> in poloidal flux equispaced

!grids
integer :: nr, nz, nrho, nteta, nr2, nz2, nr1, nz1   !nteta+1 is the periodic point. nrho is the plasma boundary

double precision :: rmin, rmax, zmin, zmax ! rho is defined as rho_toroidal as in astra

double precision :: r(i_dim2), z(i_dim2), rho(i_dim2, i_dim2), teta(i_dim2) ! rho is defined as actual distance in meters as in astra
double precision :: rcomp(i_dim2), zcomp(i_dim2) ! rho is defined as distance
double precision :: dr, dz ! teta and psigrid are equispaced
double precision :: rpol(i_dim2, i_dim2), zpol(i_dim2, i_dim2) ! R, Z in polar coordinates
double precision :: psia_2d(i_dim2), ffp_2d(i_dim2), ppp_2d(i_dim2), ipol_2d(i_dim2), pres_2d(i_dim2)
double precision :: area_eff(i_dim2, i_dim2) !

! potential
double precision :: psirz(i_dim2, i_dim2), psirhoteta(i_dim2, i_dim2)
double precision :: u_n(i_dim2, i_dim2)
double precision :: psiextrz(i_dim2, i_dim2), psiplasrz(i_dim2, i_dim2)
double precision :: zlimpotential(i_dim2, i_dim2)
double precision :: derivpsi(8), zbot, ztop, raus, rinner

! boundary and axis FBE
integer :: nbnd, ngbnd, redo_bnd !redo_bnd is temporary
double precision :: rbnd(i_dim5), zbnd(i_dim5), psibnd, psiaxis
integer ibnd(i_dim5), jbnd(i_dim5)
double precision :: rax, zax, alpsep
integer :: i_plasmatype !(0-limited, 1-single null, 2-double null)
integer :: iaxis, jaxis
integer :: n_of_xpoints
integer, parameter :: max_xpoints=80
double precision :: r_xpoint(max_xpoints), z_xpoint(max_xpoints)
double precision :: trax, tzax !true axis for more precision
double precision :: psistabR, psistabZ !stab terms
double precision :: green_bnd_f(16*i_dim2**2)

! boundary and axis PBE
double precision :: rbndp(i_dim2), zbndp(i_dim2), psibndp, psiaxisp
double precision :: raxp, zaxp, alpsepp
double precision :: rexp(i_dim2), zexp(i_dim2), tetaexp(i_dim2)

! plasma parameters
double precision :: iplasma, btor0, rgeom0, psplex, li3, betapol
double precision :: pprime(i_dim2), ffprime(i_dim2), pressure(i_dim2), psigrida(i_dim2)     !these 3 come from astra, psi is FP of astra
! current density
double precision :: jrz(i_dim2, i_dim2), jrhoteta(i_dim2, i_dim2), ipol(i_dim2)

end module ef_circuit

!------------------------------------
module metric_coefficients_pbe

use dimensions_ef_parameters, only: i_dim2

implicit none

! metric coefficients in polar coordinates 
double precision, dimension(i_dim2, i_dim2) :: lambda2d, lambda2dp

end module metric_coefficients_pbe

!------------------------------------
module green_matrix

use dimensions_ef_parameters, only: i_dim1, i_dim2

implicit none

double precision, dimension(i_dim2, i_dim2, i_dim1) :: greeni

end module green_matrix

!------------------------------------
module exchange_with_astra       ! declaration of minimal CPOs

use dimensions_ef_parameters, only: i_dim1

double precision :: tau_circuit_ef    
double precision :: tau_gseq_ef    

double precision :: activate_coil_ef(i_dim1)
double precision :: sign_coil(i_dim1)
double precision :: current_limit_ef(i_dim1, 2) ! 1 is upper, 2 is lower
double precision :: force_coil(i_dim1, i_dim1) ! where it is 1, forces coil i, i to current of i, j

integer :: refit_mode   ! if -1 - 1 pass only , 0 - self-consistent solution, if 1 - stab axis using passive wall currents fourier modes cos and sin, if 2 - ... 
integer :: execute_plasma   ! if 0 - only circuit equations, if 1 - solve plasma gseq too 
integer :: use_zlim_pot, nonegcurr ! nonegcurr = 0 --> no negative current allowed in plasma
integer :: fast_mode   ! 0 - normale, 1 - domnt do iterationsin fbe gse
integer :: reconnect_circuits   ! 0-nothing, 1-recompute matrix with new circuits
integer :: new_equivalence(i_dim1, 10)   ! if reconnect, says what is the new_equivalence, for example (1, 1, 1, 0, 0, 0, 0, ..) means coil 1, 2, 3 become 1, 1, 1
integer :: n_equivalence   ! number of equivalences
integer :: use_reduce_circuit   ! 0-all coils solved. 1 - some coils not solved
integer :: n_of_newton_iterations   ! to find actual mag axis. recommended between 5 - 10 
double precision :: raxis_astra, zaxis_astra, psi0_astra, psib_astra

end module exchange_with_astra

!----------------------------------------
module parameters_gsef

implicit none

integer, parameter :: Nrrect=60, Nzrect=61, diagnostic_gsef=0
double precision, parameter :: sorparam=1.58, lambdatol=1.E-3

character(len=10), parameter :: file_eqoutr='outpr.gsef', file_fields='bfields.wr'
character(len=9), parameter :: file_eqout='outp.gsef'

integer :: max_iter, save_gsef=1
character(len=120) :: name_gsefdir = 'exp/equ/aug/' ! working directory path

end module

!----------------------------------------
module fft_mod_feqis

use pi_grec_vars, only: gpi
   
implicit none

contains
 
! In place Cooley-Tukey FFT
    recursive subroutine fft_feqis(x)

    complex, dimension(:), intent(inout)  :: x
    complex :: t
    integer :: N, i
    complex, dimension(:), allocatable :: even, odd
 
    N = size(x)
 
    if(N .le. 1) return
 
    allocate(odd((N+1)/2))
    allocate(even(N/2))
 
! divide
    odd  = x(1: N: 2)
    even = x(2: N: 2)
 
! conquer
    call fft_feqis(odd)
    call fft_feqis(even)
 
! combine
    do i=1, N/2
        t = exp(cmplx(0.0, -2.0*gpi*(i-1.)/N))*even(i)
        x(i)     = odd(i) + t
        x(i+N/2) = odd(i) - t
    enddo
 
    deallocate(odd)
    deallocate(even)
  
    end subroutine fft_feqis
 
end module fft_mod_feqis

!-----------------
module debugger_ef

implicit none

integer :: debug
character(len=132) :: last_mark, sec_last_mark

contains

   subroutine markloc_ef(str_in, debug_lev)

   character(len=*), intent(in) :: str_in
   integer, optional, intent(in) :: debug_lev

   integer :: verbose

   if (PRESENT(debug_lev)) then
      verbose = debug_lev
   else
      verbose = debug
   endif

   SELECT CASE(verbose)
   CASE(1)
      write(*, '(/A)') 'EQEF-Trackback:'
      write(*, '(4X, A)') TRIM(str_in)
   CASE(2)
      write(*, '(/A)') 'EQEF-Trackback:'
      write(*, '(4X, A)') TRIM(last_mark)
      write(*, '(8X, A)') TRIM(str_in)
   CASE(3)
      write(*, '(/A)') 'EQEF-Trackback:'
      write(*, '( 4X, A)') TRIM(sec_last_mark)
      write(*, '( 8X, A)') TRIM(last_mark)
      write(*, '(12X, A)') TRIM(str_in)
   END SELECT

   sec_last_mark = TRIM(last_mark)
   last_mark = TRIM(str_in)

   return
   end subroutine markloc_ef

end module debugger_ef
