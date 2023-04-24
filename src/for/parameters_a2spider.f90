module parameters_a2spider       ! declaration of code parameters

use imas_ids, only: type_equilibrium

integer :: fix_adapgrid
integer, parameter, private :: DP=kind(1.0D0)
double precision, parameter :: GP=3.14159265359, GP2=2.*GP, GP4=GP2**2.0, muvac=0.4*GP*1.E-06
double precision :: epsf_tol, epss_tol, ydiff, ydiff2, epsv_tol, sorparam, epsg_tol, &
    epstol, itertol, lambdatol, urelax, urelax2, time_fix_eqpff, murelax2, cheb_coefs(5, 5)
character(len=80) :: name_gsefdir = 'exp/equ/' ! working directory path

integer :: max_iter, miter_ext, diagnostic_gsef, interp_routine, &
   do_adcmp, interp_method_rect, cheb_degree, four_degree, &
   advanced_methods, iter_one_only_fbe, spidat_yes, i3method, &
   fix_eqpf_eqff
integer :: key_no_startz  ! if 0, if key_start=1, skip iterations
integer :: key_no_refits  ! if 0, does not overwrite coil.dat with new refit currents
type(type_equilibrium) :: equil_now

type type_parameters  
    integer :: kpr=0  ! print in spider (0 - no print, -1 - no write)
    integer :: k_grid= 0   ! k_grid= 0   rect. grid
                           ! k_grid= 1   adap. grid
    integer :: k_auto= 1   ! k_auto= 1->   full initialization
    integer :: key_dmf=0 !=1->diff.mag.field, =0->without

    integer :: nstep = 0   ! nstep=0 - initial eq., nstep>0 using computed eq.

    character(len=80) :: prename = 'exp/equ/' ! working directory path
    integer :: kname = 8   ! path name length >=1
    character(len=80) :: eqdfn = ''
    
    integer ::  k_fixfree=1  !=0->only fixed boundary spider 
    integer ::  k_filesss=1  !=1->use files, 0 use memory 
    integer ::  key_ini=1   ! =1 astra profiles, =0 start from EQDSK and SPIDER profiles
    integer ::  key_start=0 ! =1 reconstruction, =0 direct for free-boundary equilibrium
    integer ::  key_0stp=0   !initial eq. only =0 - p',ff'; =1 - p,cu
    integer ::  key_pres=0   !=1 - pressure profile, =0 - p' profile  
    integer ::  i_eqdsk=0   !=1 - eqdsk file as input, =0 - other user input 
    integer ::  key_plc=1   !=1 - precribed Ip, =0 - no G-S rhs renormalization , ipl is an output anyway 
    integer ::  key_out=0   !=1 - circuit equations with currents and inductive voltages update, =0 - no update, 0 is to do iterations, last one has to have key_out=1

    integer ::  key_psibcf=1   !=1 - compute psiext and psipl for b.c. apt for current control
        
    real(DP) :: dt=1.d-3
    real(DP) :: time=0.d0
    real(DP) :: dpsdt=0.d0   ! boundary vloop in input

    real(DP) :: epsro=1.0d-7 ! fixed boundary equilibrium accuracy
    real(DP) :: enels=1.0d-6 ! circuit equation accuracy
    
    integer :: neql = 100   ! number of nodes in radial
    integer :: nteta = 90  ! number of intervals in poloidal + 2

    integer :: n_dmf = 3  ! number of iterations of cde in SPIDER with rectangular grid

endtype

end module parameters_a2spider
