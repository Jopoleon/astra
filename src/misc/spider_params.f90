module spider_params

integer, parameter, private :: DP=kind(1.0D0)

type type_parameters

    integer :: kpr=0  ! print in equil (0 - no print, -1 - no write)
    integer :: k_grid= 0  ! k_grid= 0   rect. grid-
                          ! k_grid= 1   adap. grid
    integer :: k_auto= 1  ! k_auto= 1->   full initialization
    integer :: key_dmf=0  !=1->diff.mag.field, =0->without
    integer :: nstep = 0  ! nstep=0 - initial eq., nstep>0 using computed eq.

    character(len=80) :: prename = 'exp/equ/' ! working directory path
    character(len=80) :: eqdfn = ''

    integer ::  k_fixfree=1 !=0->only fixed boundary spider 
    integer ::  k_filesss=1 !=1->use files, 0 use memory 
    integer ::  key_ini=1   ! =1 astra profiles, =0 start from EQDSK and SPIDER profiles
    integer ::  key_start=0 ! =1 reconstruction, =0 direct for free-boundary equilibrium
    integer ::  key_0stp=0  !initial eq. only =0 - p',ff'; =1 - p,cu
    integer ::  key_pres=0  !=1 - pressure profile, =0 - p' profile  
    integer ::  i_eqdsk=0   !=1 - eqdsk file as input, =0 - other user input 
    integer ::  key_plc=1   !=1 - precribed Ip, =0 - no G-S rhs renormalization , ipl is an output anyway 
    integer ::  key_out=0   !=1 - circuit equations with currents and inductive voltages update, =0 - no update, 0 is to do iterations, last one has to have key_out=1

    integer ::  key_psibcf=1   !=1 - compute psiext and psipl for b.c. apt for current control

    real(DP) :: dt=1.d-3
    real(DP) :: time=0.d0
    real(DP) :: dpsdt=0.d0   ! boundary vloop in input

    real(DP) :: epsro=1.0d-7 ! fixed boundary equilibrium accuracy
    real(DP) :: enels=1.0d-6 ! circuit equation accuracy

    integer :: neql = 100 ! number of nodes in radial
    integer :: ntheta = 90 ! number of intervals in poloidal + 2
    integer :: n_dmf = 3  ! number of iterations of cde in SPIDER with rectangular grid
    integer :: no_circuit_eq = 0 ! if 1, doesnt do circuit equations

endtype

end module spider_params
