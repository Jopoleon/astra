subroutine facit(Z_imp_in, A_imp_in, N_imp_in, rot_mod, Dz_out, Vz_out)

  !============================================================================================!
  !
  ! Calculate flux surface averaged neoclassical impurity transport coefficients with FACIT
  !
  ! --- D. Fajardo, July 2024
  !
  ! * INPUTS
  ! --------
  ! - Z_imp_in -> impurity charge [-] (profile, radial array NA1)
  ! - A_imp_in -> impurity mass [-] (float)
  ! - N_imp_in -> impurity density [10^19/m^3] (profile, radial array NA1)
  ! - rot_mod --> 
  !
  ! * OUTPUTS
  ! ---------
  ! - Dz_out ---> FSA impurity diffusion coefficient [m^2/s] (profile, radial array NA1)
  ! - Vz_out ---> FSA impurity convective velocity [m/s] (profile, radial array NA1)
  !
  !============================================================================================!

  use parameter_inc, only: NRD
  use const_inc,     only: NA1, BTOR, RTOR, AMJ
  use status_inc,    only: TE, TI, NE, IPOL, MU, SQEPS, AMETR, VTOR, ZEF, ZMAIN, NMAIN, VRS, G11, RHO

  implicit none

  ! inputs and outputs
  double precision, dimension(NRD), intent(in) :: Z_imp_in, N_imp_in ! input impurity charge and density profiles
  double precision, intent(in) :: A_imp_in, rot_mod                  ! input impurity mass and rotation model
  double precision, dimension(NRD), intent(out) :: Dz_out, Vz_out    ! diffusion coefficient and convection velocity for output

  ! physical constants
  double precision, parameter :: mp         = 1.67262192369e-27 ! proton mass [kg]
  double precision, parameter :: q_e        = 1.602176634e-19   ! electron charge [C]

  integer :: i,j,jmax ! dummy for loops

  ! INPUTS

  integer :: nx                                         ! grid parameters
  double precision, dimension(NA1)    :: r_min, epsilon ! grid parameters

  double precision, dimension(NA1)  :: Zimp, Zi ! impurity and main ion charge
  double precision :: Aimp, Ai          ! impurity and main ion mass

  double precision, dimension(NA1)  :: T_e, T_i, N_e, N_i, N_z, Mach, Zeff ! plasma profiles
  double precision, dimension(NA1)  :: grad_Ti, grad_Ni, grad_Nz            ! gradients

  double precision :: B0, R0                     ! equilibrium
  double precision, dimension(NA1)  :: qmag, FF  ! equilibrium

  integer :: rotation_model

  double precision, dimension(NA1) :: qstar, torflux_ov2pi, dtorflux_dr  ! modified q profile for flux surface geometry
  double precision, dimension(NA1) :: grad_rho_sq, r_tor, drtor_drmin !gradrhosq_exp, drhodr, drho, drmin ! to convert to ASTRA grid
  
  ! OUTPUTS

  double precision, dimension(NA1) :: Dz_lfs, Vz_lfs, Dz_fsa, Vz_fsa

  ! OTHER:
  !double precision:: cimp, fH, bC, sigH, nsigH, TperpsTpar_axis
  !double precision, dimension(NA1) :: TperpsTpar, Zeff

  nx  = int(NA1)    ! number of radial points

  !-----------------------------------------------------------------------------------------------!
  ! model options
  
  rotation_model = int(rot_mod) ! control rotation model from equ file input to facit sbr

  !-----------------------------------------------------------------------------------------------!
  ! plasma profiles
  
  T_e   = 1.0e3*max(TE(1:NA1),1.e-5)               ! electron temperature in eV
  T_i   = 1.0e3*max(TI(1:NA1),1.e-5)               ! main ion temperature in eV
  N_e   = 1.0e19*max(NE(1:NA1),1.e-5)              ! electron density in 1/m^3
  N_i   = 1.0e19*max(NMAIN(1:NA1),1.e-5)           ! main ion density in 1/m^3
  N_z   = 1.0e19*max(N_imp_in(1:NA1),1.e-8)        ! impurity density in 1/m^3 from equ file input to impflux sbr
  
  Zeff  = ZEF(1:NA1)                               ! effective charge 
  Zimp  = Z_imp_in(1:NA1)                          ! impurity charge profile from equ file input to impflux sbr
  Aimp  = A_imp_in                                 ! impurity mass from equ file input to impflux sbr
  Zi    = ZMAIN(1:NA1)                             ! main ion charge number
  Ai    = AMJ                                      ! main ion mass number
  Mach = abs(VTOR(1:NA1))/sqrt(2*q_e*T_i/(Ai*mp))  ! main ion Mach number
  
  !-----------------------------------------------------------------------------------------------!
  ! equilibrium

  r_min   = AMETR(1:NA1)     ! local minor radius [m]
  epsilon = SQEPS(1:NA1)**2. ! local inverse aspect ratio

  B0       = abs(BTOR)       ! magnetic field at magnetic axis
  R0       = RTOR            ! major radius at magnetic axis
  qmag     = 1.0/MU(1:NA1)   ! safety factor

  FF = IPOL(1:NA1)*R0*B0

  ! modified q profile for flux surface geometry
  r_tor = RHO(1:NA1)
  torflux_ov2pi = 0.5 * B0 * r_tor**2

  do i=2,nx-1
     dtorflux_dr(i) = (torflux_ov2pi(i+1)-torflux_ov2pi(i-1))/(r_min(i+1)-r_min(i-1))
  enddo
  dtorflux_dr(nx) = (torflux_ov2pi(nx)-torflux_ov2pi(nx-1))/(r_min(nx)-r_min(nx-1))
  dtorflux_dr(1) = dtorflux_dr(2)
  
  qstar = qmag * (epsilon * FF / dtorflux_dr)
  qstar(1) = qstar(2)
  
  !-----------------------------------------------------------------------------------------------!
  ! gradients
  
  do i=2,nx-1
     grad_Ni(i) = (N_i(i+1)-N_i(i-1))/(r_min(i+1)-r_min(i-1))
     grad_Ti(i) = (T_i(i+1)-T_i(i-1))/(r_min(i+1)-r_min(i-1))
     grad_Nz(i) = (N_z(i+1)-N_z(i-1))/(r_min(i+1)-r_min(i-1))
  enddo

  grad_Ni(nx) = (N_i(nx)-N_i(nx-1))/(r_min(nx)-r_min(nx-1))
  grad_Ti(nx) = (T_i(nx)-T_i(nx-1))/(r_min(nx)-r_min(nx-1))
  grad_Nz(nx) = (N_z(nx)-N_z(nx-1))/(r_min(nx)-r_min(nx-1))

  ! force gradients to be zero at rho = 0.0
  grad_Ni(1) = 0.0
  grad_Ti(1) = 0.0
  grad_Nz(1) = 0.0
  
  !-----------------------------------------------------------------------------------------------!
  ! FACIT call for LFS coefficients

  !print *, "Calling FACIT, rotation model", rotation_model

  call FACIT_LFS(nx, epsilon, &                          ! grid parameters
                 Zimp, Aimp, Zi, Ai, &                   ! impurity and main ion charge and mass
                 T_e, T_i, N_e, N_i, N_z, Mach, Zeff, &  ! plasma profiles
                 grad_Ti, grad_Ni, grad_Nz, &            ! gradients
                 B0, R0, qstar, FF, &                    ! equilibrium
                 rotation_model, &                       ! model options
                 Dz_lfs, Vz_lfs)                         ! output coefficients

  !print *, "Done with FACIT call"

  !-----------------------------------------------------------------------------------------------!

  ! Transform to ASTRA grid (r_min -> r_tor)
  grad_rho_sq = G11(1:NA1)/VRS(1:NA1)

  do i=2,nx-1
     drtor_drmin(i) = (r_tor(i+1)-r_tor(i-1))/(r_min(i+1)-r_min(i-1))
  enddo

  drtor_drmin(nx) = (r_tor(nx)-r_tor(nx-1))/(r_min(nx)-r_min(nx-1))
  drtor_drmin(1) = (r_tor(2)-r_tor(1))/(r_min(2)-r_min(1))
  
  ! In case rotation model = 2, transform LFS to FSA coefficients, otherwise already output is FSA

  ! extract output for arguments of equ file call to FACIT sbr
  if (rotation_model .eq. 0 .or. rotation_model .eq. 1) then
     
     Dz_out(1:NA1) = Dz_lfs*(drtor_drmin**2/grad_rho_sq)
     Vz_out(1:NA1) = Vz_lfs*(drtor_drmin/grad_rho_sq)

  elseif (rotation_model .eq. 2) then

     call lfs2fsa_impDV(2, 1, Zimp, Aimp, Dz_lfs, Vz_lfs, Dz_fsa, Vz_fsa)

     Dz_out(1:NA1) = Dz_fsa*(drtor_drmin**2/grad_rho_sq)
     Vz_out(1:NA1) = Vz_fsa*(drtor_drmin/grad_rho_sq)

  else
     print *, "Please select correct rotation model = 0, 1 or 2. It's the 4th argument of the FACIT call"

  endif

end subroutine facit

!*******************************************************************************
!**************************** Main subroutine **********************************
!*******************************************************************************


subroutine FACIT_LFS(nx, eps, &                                 ! grid parameters
                     Zz, Az, Zi, Ai, &                          ! impurity and main ion charge and mass
                     T_e, T_i, N_e, N_i, N_z, Mach_ii, Z_eff, & ! plasma profiles
                     gradTi, gradNi, gradNz, &                  ! gradients
                     B0, R0, qmag, FV, &                        ! equilibrium
                     rotation_model, &                          ! model options
                     Dz, Vz)                                    ! output coefficients

!*******************************************************************************
! Self-consistent calculation of collisional impurity transport and poloidal
! distribution of the impurity density
!------------
! References:
! Maget et al 2020 Plasma Phys. Control. Fusion 62 105001
! Fajardo et al 2022 Plasma Phys. Control. Fusion 64 055017
! Fajardo et al 2023 Plasma Phys. Control. Fusion 65 035021
!-------------------
! License statement:
! Software name : FACIT
! Authors : P. Maget and D. Fajardo, C. Angioni, P. Manas
! Copyright holders : Commissariat à l’Energie Atomique et aux Energies Alternatives
!                     (CEA), France; Max-Planck Institut für Plasmaphysik, Germany
! CEA and IPP authorize the use of the FACIT software under the CeCILL-C open
! source license https://cecill.info/licences/Licence_CeCILL-C_V1-en.html  
! The terms and conditions of the CeCILL-C license are deemed to be accepted upon downloading 
! the software and/or exercising any of the rights granted under the CeCILL-C license.
!*******************************************************************************
! INPUTS:
! -----------> description [unit] {variable type, shape}
! ------------
! nx --------> size of radial arrays [-] {int}
! nth -------> size of poloidal arrays [-] {int}
! xn --------> radial coordinate [-] {arr, nx}
! theta -----> poloidal coordinate [-] {arr, nth}
! Za --------> impurity charge number [-] {arr, nx}
! Aa --------> impurity mass number [-] {float}
! Zi --------> main ion charge number [-] {float}
! Ai --------> main ion mass number [-] {float}
! Te --------> electron temperature [eV] {arr, nx}
! Ti --------> main ion temperature [eV] {arr, nx}
! Ne --------> electron density [1/m^3] {arr, nx}
! Ni --------> main ion denisty [1/m^3] {arr, nx}
! Na --------> impurity density [1/m^3] {arr, nx}
! Machi -----> Mach number of main ion [-] {arr, nx}
! Zeff ------> Effective charge [-] {arr, nx}
! gradTi ----> main ion temperature gradient [eV/m] {arr, nx}
! gradNi ----> main ion density gradient [1/m^3/m] {arr, nx}
! gradNz ----> impurity density gradient [1/m^3/m] {arr, nx}
! invaspct --> inverse aspect ratio [-] {float}
! B0 --------> magnetic field at magnetic axis [T] {float}
! R0 --------> major radius at magnetic axis [m] {float}
! qmag ------> safety factor [-] {arr, nx}
! dpsidx ----> radial derivative of poloidal flux [V*s/-] {arr, nx}
! FV --------> poloidal current flux function [T*m] {arr, nx}
! BV --------> magnetic field [T] {arr, (nx,nth)}
! RV --------> major radius [m] {arr, (nx,nth)}
! jacob -----> Jacobian of the coordinate system [m/T] {arr, (nx,nth)}
! AsymPhi ---> poloidal asymmetry of electrostatic potential [-] {arr, (nx,2)}
! AsymN -----> poloidal asymmetry of main ion density [-] {arr, (nx,2)}
! pol_asym --> poloidally symmetric (False) or asymmetric (True) system {boolean}
! full_geom -> analytical (False) or iterative (True) geometry {boolean}
! rotation --> consider toroidal rotation {boolean}
! regulopt --> options for iterative calculations [-] {arr, 4}
!*******************************************************************************
! OUTPUTS:
! --------
! Da --------> diffusion coefficient [m^2/s] {arr, nx}
! Vconv -----> convective velocity [m/s] {arr, nx}
!*******************************************************************************
  
  implicit none

  ! mathematical constants
  double precision, parameter :: pi    = 4.d0*ATAN(1.d0)
  double precision, parameter :: sqrt2 = SQRT(2.d0)

  ! physical constants
  double precision, parameter :: mp         = 1.67262192369e-27           ! proton mass [kg]
  double precision, parameter :: me         = 9.1093837015e-31            ! electron mass [kg]
  double precision, parameter :: q_e        = 1.602176634e-19             ! electron charge [C]
  double precision, parameter :: eps0       = 8.8541878128e-12            ! vacuum permittivity [F/m]
  double precision, parameter :: eps_pi_fac = 3.0*eps0**2*(2.0*pi/q_e)**1.5/q_e ! to be used in Coulomb logarithms [F^2/(m^2 C^(5/2))]

  !------------------------------ declarations --------------------------------

  ! INPUTS
  integer, intent(in) :: nx, rotation_model
  double precision, dimension(nx), intent(in)  :: T_e, T_i, N_e, N_i, N_z, gradTi, gradNi, gradNz
  double precision, dimension(nx), intent(in)  :: qmag, FV, Zz, Mach_ii, Z_eff, eps
  double precision, intent(in) :: B0, R0, Ai, Az, Zi

  ! OUTPUTS

  double precision, dimension(nx), intent(out) :: Dz, Vz

  
  double precision, dimension(nx) :: Flux_imp, dmin, dmaj
  double precision, dimension(nx) :: Da_BP, Da_PS, Da_CL, Ka_BP, Ka_PS, Ka_CL, Ha_BP, Ha_PS, Ha_CL, Va_BP, Va_PS, Va_CL

  ! OTHER
  integer :: i, p, n, it, ix!, info, ierr, ierrmax
  double precision :: amin, mi, mz, ftrap_lfs, ki_Redl_lfs, ki_rot_lfs, C2_lfs
  double precision, dimension(nx) :: grad_ln_ni, grad_ln_Ti, grad_ln_nz, Ti15, Te15, Zeff
  double precision, dimension(nx) :: epsk, eps15, eps2, ft, C0a, wca, deltaM, ki, g, dD2
  double precision, dimension(nx) :: f1, f2, f3, fG, fU, adps, ahbp, fvv
  double precision, dimension(nx) :: y11zb, y11zp, y11zps, y11ib, y11ip, y11ips
  double precision, dimension(nx) :: y12zb, y12zp, y12zps, y12ib, y12ip, y12ips
  !double precision, dimension(nx) :: LneeNRL, LneiNRL, LneimpNRL
  double precision, dimension(nx) :: LniiNRL, LniimpNRL, LnimpimpNRL!, LnimpeNRL
  !double precision, dimension(nx) :: Tauee, Tauei, Taueimp
  double precision, dimension(nx) :: Tauii, Tauiimp!, Tauie
  double precision, dimension(nx) :: Tauimpi, Tauimpimp!, Tauimpe
  double precision, dimension(nx) :: wii, wimpimp, nuistar, rhoLimp2!, wee, nuestar, nuimpstar
  double precision, dimension(nx) :: L11impi, nuswca, mu_ie, alpha!, Zeff
  double precision, dimension(nx) :: Cgeo_G, Cgeo_U, Cgeo_Gcl
  double precision, dimension(nx) :: K11a, K12a, K22a, K11i, K12i, K22i
  double precision, dimension(nx) :: Vra_BP, Vra_PS, Vra_CL
  double precision, dimension(nx) :: Ka, Ha
  double precision, dimension(nx) :: Mzstar, Mach_ion

  
  !---------------------------------------------------------------------------

  !print *, Mach_ii

  mi = Ai*mp ! main ion mass
  mz = Az*mp ! impurity mass

  eps15 = eps**1.5
  eps2  = eps**2.
  

  if (maxval(Z_eff).eq.1.0) then
     Zeff = (Zz**2*N_z + Zi**2*N_i)/(N_e) ! effective charge
  else
     Zeff = Z_eff
  endif

  do i = 1, nx
    ! logarithmic gradients
    grad_ln_ni(i) = gradNi(i)/(N_i(i) + 1.e-33)
    grad_ln_Ti(i) = gradTi(i)/(T_i(i) + 1.e-33)
    grad_ln_nz(i) = gradNz(i)/(N_z(i) + 1.e-33)
    ! trapped particle fraction
    ft(i) = ftrap_lfs(eps(i))
  enddo


  if (rotation_model.eq.0) then
     Mach_ion = 0.0
  else
     Mach_ion = Mach_ii
  endif

  Mzstar = Mach_ion*sqrt(Az/Ai - (Zz/Zi)*Zeff/(Zeff + T_i/T_e)) ! effective impurity Mach number

  ! Coulomb Logarithms (from NRL formulary)

  LniiNRL     = 23. - log(Zi*Zi*sqrt(2*(N_i/1.e6)*Zi**2)) + 1.5*log(T_i)
  LniimpNRL   = 23. - log(Zi*Zz*sqrt((N_i/1.e6)*Zi**2+(N_z/1.e6)*Zz**2)) + 1.5*log(T_i)
  LnimpimpNRL = 23. - log(Zz*Zz*sqrt((N_z/1.e6)*Zz**2+(N_z/1.e6)*Zz**2)) + 1.5*log(T_i)


  ! Collision times (Braginskii)
  Ti15 = T_i**1.5
  Te15 = T_e**1.5

  Tauii     = (eps_pi_fac*sqrt(mi)*Ti15)/(Zi**4*N_i*LniiNRL)
  Tauiimp   = (eps_pi_fac*sqrt(mi)*Ti15)/(Zi**2*Zz**2*N_z*LniimpNRL)
  Tauimpi   = (eps_pi_fac*sqrt(mz)*Ti15)/(Zi**2*Zz**2*N_i*LniimpNRL)
  Tauimpimp = (eps_pi_fac*sqrt(mz)*Ti15)/(Zz**4*N_z*LnimpimpNRL)


  ! impurity collision frequency
  L11impi = 1.0/(sqrt(1.0 + Az/Ai)*Tauimpi)

  ! transit frequencies
  wii     = (2.0*q_e*T_i/mi)**0.5/(R0*qmag)
  wimpimp = (2.0*q_e*T_i/mz)**0.5/(R0*qmag)

  ! collisionalities
  nuistar   = 1.0/((eps15 + 1.e-33)*wii*Tauii)
  wca    = q_e*Zz*B0/mz ! impurity cyclotron frequency
  

  ! fitted factors
  call facs_lfs(nx, Zz, Az, ft, Mzstar, rotation_model, &
            f1, f2, f3, fG, fU, adps, ahbp, fvv, &
            y11zb, y11zp, y11zps, y11ib, y11ip, y11ips, &
            y12zb, y12zp, y12zps, y12ib, y12ip, y12ips)


  g     = nuistar*eps15 ! collisionality parameter
  mu_ie = (96.0*sqrt2/125.0)*(1/Zi**2)*sqrt(me/mi)*(Ti15/Te15) ! ion-electron heat exchange term (Fülöp-Helander PoP '01)
  alpha = N_z*Zz**2/(N_i*Zi**2) ! impurity strength parameter

  do i = 1, nx
     C0a(i) = C2_lfs(alpha(i), g(i), f1(i), f2(i), Az, Ai)/(1 + f3(i)*mu_ie(i)*g(i)**2) ! coefficient of ion heat flux in impurity-ion friction

     if (rotation_model.eq.2) then
        ki(i) = ki_rot_lfs(nuistar(i), ft(i), Zeff(i), Mach_ion(i)) ! neoclassical ion flow coefficient
     else
        ki(i) = ki_Redl_lfs(nuistar(i), ft(i), Zeff(i)) ! neoclassical ion flow coefficient
     endif
     
  enddo

  rhoLimp2 = (2.0*q_e*T_i/mz)/wca**2 ! Impurity Larmor radius (squared)

  !---------------------------------------------------------------------------
  !-----------------------  Poloidal asymmetry  ------------------------------
  !---------------------------------------------------------------------------


  if (rotation_model.eq.0) then

     dmin = 0.0
     dmaj = 0.0

     Cgeo_G = 2*eps2*(0.96*(1 - 0.54*ft**4.5))
     Cgeo_U = 0.0
     Cgeo_Gcl = 1.0 + 2*eps2

  elseif (rotation_model.eq.1) then

     dmin = 0.0
     dmaj = 0.0

     Cgeo_G = 2*eps2*(0.96*(1 - 0.54*ft**4.5))
     Cgeo_U = 0.0
     Cgeo_Gcl = 1.0 + 2*eps2

  elseif (rotation_model.eq.2) then

      Cgeo_G = 2*eps2*(0.96*(1 - 0.54*ft**4.5))*fG
      Cgeo_U = -2*eps2*(0.96*(1 - 0.54*ft**4.5))*fU
      Cgeo_Gcl = 1.0 + 2*eps2

      dmin = 0.0
      dmaj = 0.0

  endif


  !---------------------------------------------------------------------------
  !---------------------------  Impurity flux  -------------------------------
  !---------------------------------------------------------------------------

  ! General form: Va = -D*grad_ln_na + (K*grad_ln_ni + H*grad_ln_Ti + Vrot)
  !                  = -D*grad_ln_na + Vconv


  ! Pfirsch-Schlüter flux

  Da_PS   = adps*qmag**2*rhoLimp2*L11impi*(Cgeo_G/(2.0*eps2))
  Ka_PS   = (Zz/Zi)*Da_PS
  Ha_PS   = -((1.0 + (Zz/Zi)*(C0a - 1.0)) + (Cgeo_U/Cgeo_G)*(Zz/Zi)*(C0a + ki))*Da_PS

  Vra_PS  = -Da_PS*grad_ln_nz + Ka_PS*grad_ln_ni + Ha_PS*grad_ln_Ti
  Va_PS  = Ka_PS*grad_ln_ni + Ha_PS*grad_ln_Ti


  ! Classical flux

  Da_CL   = (2.0*eps2*Cgeo_Gcl/Cgeo_G)*Da_PS/(2*qmag**2)
  Ka_CL   = (Zz/Zi)*Da_CL
  Ha_CL   = -(1.0 + (Zz/Zi)*(C0a - 1.0))*Da_CL

  Vra_CL  = -Da_CL*grad_ln_nz + Ka_CL*grad_ln_ni + Ha_CL*grad_ln_Ti
  Va_CL = Ka_CL*grad_ln_ni + Ha_CL*grad_ln_Ti


  ! Banana-Plateau flux

  call K_VISC_lfs(nx, N_i, N_z, T_i, wii, wimpimp, Zi, Zz, Ai, Az, Tauii, &
              Tauimpi, Tauiimp, Tauimpimp, eps2, R0, qmag, ft, &
              y11zb, y11zp, y11zps, y11ib, y11ip, y11ips, &
              y12zb, y12zp, y12zps, y12ib, y12ip, y12ips, &
              K11a, K12a, K22a, K11i, K12i, K22i)


  Da_BP = (1.5*q_e*T_i*(1.0/(1.0/K11a + 1.0/K11i))/(Zz**2*q_e**2*FV**2*N_z))
  Ka_BP = (Zz/Zi)*Da_BP
  Ha_BP = ahbp*((Zz/Zi)*(K12i/K11i - fvv) - (K12a/K11a - fvv))*Da_BP

  Vra_BP = -Da_BP*grad_ln_nz + Ka_BP*grad_ln_ni + Ha_BP*grad_ln_Ti
  Va_BP  = Ka_BP*grad_ln_ni + Ha_BP*grad_ln_Ti
 

  ! Total transport coefficients

  Dz = Da_PS + Da_BP + Da_CL ! total diffusion coefficient
  Ka = Ka_PS + Ka_BP + Ka_CL
  Ha = Ha_PS + Ha_BP + Ha_CL

  !Vra = Vra_PS + Vra_BP + Vra_CL

  ! total convective velocity
  Vz = Ka*grad_ln_ni + Ha*grad_ln_Ti

  ! Total surface-averaged flux
  Flux_imp = -Dz*gradNz + N_z*Vz


  return
end subroutine FACIT_LFS



!*******************************************************************************
!*********************** Complementary subroutines *****************************
!*******************************************************************************


function ftrap_lfs(epsK)
!*******************************************************************************
! Trapped particle fraction as a function of the local inverse aspect ratio
!*******************************************************************************
! INPUT:
! ------
! epsK --> local inverse aspect ratio [-] {float}
!*******************************************************************************
! OUTPUT:
! ------
! ftrap -> trapped particle fraction [-] {float}
!*******************************************************************************
  implicit none

  double precision :: ftrap_lfs
  double precision, intent(inout) :: epsK

  ftrap_lfs = 1. - (1. - epsK)**1.5/(sqrt(1. + epsK)*(1 + 1.46*sqrt(epsK)))

  return
end function ftrap_lfs




function ki_Redl_lfs(nuistar, ft, Zeff)

!*******************************************************************************
! Neoclassical ion flow coefficient, fitted w.r.t. NEO in Redl PoP (2021)
!*******************************************************************************
! INPUTS:
! -------
! nuistar -> main ion collisionality [-]
! ft ------> trapped particle fraction [-]
! Zeff ----> effective charge [-]
!*******************************************************************************
! OUTPUT:
! ------
! ki ------> main ion flow coefficient [-]
!*******************************************************************************

  implicit none

  double precision, intent(in) :: nuistar, ft, Zeff
  double precision :: alpha0
  double precision :: ki_Redl_lfs

  alpha0 = -(0.62 + 0.055*(Zeff - 1.0))*(1.0 - ft)/((0.53 + 0.17*(Zeff - 1.0))*&
            (1.0 - (0.31 - 0.065*(Zeff - 1.0))*ft - 0.25*ft**2))

  ki_Redl_lfs = ((alpha0 + 0.7*Zeff*sqrt(ft*nuistar))/(1.0 + 0.18*sqrt(nuistar))-&
              0.002*nuistar**2*ft**6)/(1.0 + 0.004*nuistar**2*ft**6)

  return
end function ki_Redl_lfs



function ki_rot_lfs(nuistar, ft, Zeff, Mach_i)

!*******************************************************************************
! Neoclassical ion flow coefficient, fitted w.r.t. NEO with rotation
! in Fajardo PPCF (202X)
!*******************************************************************************
! INPUTS:
! -------
! nuistar -> main ion collisionality [-]
! ft ------> trapped particle fraction [-]
! Zeff ----> effective charge [-]
! Machi ---> Main ion Mach number [-]
!*******************************************************************************
! OUTPUT:
! ------
! ki ------> main ion flow coefficient [-]
!*******************************************************************************

  implicit none

  double precision, intent(in) :: nuistar, ft, Zeff, Mach_i
  double precision :: ki0, c01, c02, c03, l1k, l2k, l3k, l4k, l5k, l6k
  double precision :: ki_rot_lfs


  c01 = 0.53*(1.0 + 0.646*Mach_i**1.5)
  c02 = 1.158*(1.0 - 0.968*Mach_i**1.56)
  c03 = -0.98*(1.0 - 1.228*Mach_i**1.7)
        
  l1k = 5.7*(1.0-ft)**6.7 + 0.38
  l2k = (-1.52 + 38.4*(1.0-ft)**3.02*ft**2.07) + (-1.0 + 2.6*ft)*Mach_i**(2.5*(1.0 - 0.6*ft))
  l3k = 0.25 + 1.2*(1.0-ft)**3.65
  l5k = 0.1*(1.0-ft)**1.46*ft**4.33 + 0.051*(1.0 - 0.82*ft)*Mach_i**2.5
  l4k = 0.8 + (1.25*ft + 0.585)*Mach_i
  l6k = (-0.05 + 1.95*ft**2.5)/(1.0 + 2.55*ft**17) + -((0.217 + 14.57*ft**6.3)/(1.0 + 5.62*ft**5.72))*Mach_i**1.5
    
  ki0 = -(c01 + 0.055*(Zeff-1.0))*(1.0-ft)/((0.53 + 0.17*(Zeff-1.0))*&
         (1.0-(c02 - 0.065*(Zeff-1.0))*ft - c03*ft**2))
        
  ki_rot_lfs = ((ki0 + (l1k)*Zeff*sqrt(ft*nuistar) + l2k*nuistar**(0.25))/(1+l3k*sqrt(nuistar))-&
           l4k*l5k*nuistar**2*ft**6 + l6k*nuistar**(0.25))/(1.0 + l5k*nuistar**2*ft**6)

  return
end function ki_rot_lfs


function C2_lfs(alpha, g, f1, f2, Aimp, Ai)

!*******************************************************************************
! new parametrized version of the C2 function
! see Hirshman-Sigmar NF (1981), page 61 (1138)
!*******************************************************************************
! INPUTS:
! -------
! alpha -> impurity strength parameter [-]
! g -----> main ion collisionality parameter [-]
! f1 ----> f1(Zimp) factor [-]
! f2 ----> f2(Zimp) factor [-]
!*******************************************************************************
! OUTPUT:
! ------
! C2 ----> term in impurity-ion friction C0a [-]
!*******************************************************************************

  implicit none

  double precision, intent(in) :: alpha, g, f1, f2, Aimp, Ai
  double precision :: C2_lfs

  C2_lfs = 1.5/(1.0 + (2.0/Aimp)*f1) - (0.29 + 0.68*alpha)/(0.59 + alpha + (1.34 + f2)/g**2)

  return
end function C2_lfs



subroutine facs_lfs(nx, Zimp, Aimp, ft, Mzstar, rotation_model, &
                f1, f2, f3, fG, fU, adps, ahbp, fvv, &
                y11zb, y11zp, y11zps, y11ib, y11ip, y11ips, &
                y12zb, y12zp, y12zps, y12ib, y12ip, y12ips)

!*******************************************************************************
! set of factors fitted with respect to NEO
!*******************************************************************************
! INPUTS:
! -------
! nx --------------> size of radial arrays [-]
! Zimp ------------> impurity charge [-]
! Aimp ------------> impurity mass [-]
! ft --------------> trapped particle fraction [-]
! Mzstar-----------> effective impurity Mach number [-]
! rotattion_model -> int: 0,1,2
!*******************************************************************************
! OUTPUTS:
! --------
! f1 --------> low collisionality saturation of PS C0a coefficient [-]
! f2 --------> correction to low-Z impurity-ion friction in C0a coefficient [-]
! f3 --------> factor in i-e heat exchange contribution to C0a [-]
! fG,fU -----> factors in Cgeo_G, Cgeo_U [-] (only if rotation_model = 2)
! adps,ahbp -> factors multpliying Da_PS, Ha_BP [-]
! fvv -------> factor in Ha_BP structure [-]
! yjksc -----> factors of c = B, P, PS component of the j,k viscosity
!              coefficient of species s = i,z [-]
!*******************************************************************************

  implicit none

  integer :: nx, rotation_model
  double precision, dimension(nx) :: Zimp, ft, Mzstar
  double precision :: Aimp

  double precision, dimension(nx) :: f1, f2, f3, fG, fU
  double precision, dimension(nx) :: adps, ahbp, fvv
  double precision, dimension(nx) :: y11zb, y11zp, y11zps, y11ib, y11ip, y11ips
  double precision, dimension(nx) :: y12zb, y12zp, y12zps, y12ib, y12ip, y12ips

  double precision, dimension(nx) :: c1, c2, c3, c4, c5, c6     ! auxiliary
  double precision, dimension(nx) :: l1, l2, l3, l4, l5, l6, l7 ! auxiliary


  ! f's:
  f1 = (1.74*(1.0 - 0.028*Aimp) + 10.25/(1.0 + Aimp/3.0)**2.8) - 0.423*(1.0 - 0.136*Aimp)*ft**(5.0/4.0)
  f2 = (88.28389935 + 10.50852772*Zimp)/(1.0 + 0.2157175*Zimp**2.57338463)
  f3 = (-4.45719179e+06 + 2.72813878e+06*Zimp)/(1.0 + 5.26716920e+06*Zimp**8.33610108e-01)

  adps = ((0.711 + 2.08e-3*Zimp**1.26)/(1.0 + 1.06e-11*Zimp**5.78))

  fvv = 1.5

  ! y's
  y11zb  = (8.0*(1.0-ft)**20 - 0.66*ft**4.9 + 0.94)*((1.085e4 + 9.3e3*14.0/Zimp**(5.0/3.0))/(1.0 + 14.0/Zimp**(5.0/3.0)))*ft**4.6*(1.0-ft)/&
           (1.0 + 9.44e3*(1.0 - 2e-3*Zimp)*ft**4.16)
  y11zp  = 1.0
  y11zps = 1.0
        
  y12zb  = (1.8e3 + 7.54*Zimp**1.8)*ft**(4.05*(1.0 + 0.0039*Zimp))*(1.0-ft)**(0.7*(1.0 + 0.015*Zimp))/(1.0 + 1276.0*(1.0 + 0.053*Zimp)*ft**3.6)
  y12zp  = 1.0
  y12zps = 1.0
        
        
  y11ib  = ((5.91e-5*Zimp + 0.812 + 0.806/Zimp**0.44) + (0.013*Zimp + 0.098 - 7.03/Zimp**1.04)*ft + &
            (-0.047*Zimp + 0.79 + 13.8/Zimp**1.26)*ft**2 + (0.04*Zimp - 0.575 - 9.64/Zimp**1.31)*ft**3)*&
            (571.6*(1.0 + 0.84*Zimp + 7.8e-7*Zimp**5.15)*ft**(3.43*(1.0 + 0.012*Zimp))*&
            (1.0-ft)**(1.2*(1.0 + 1.6e-8*Zimp**5.1))/(1.0 + 500*ft**(8.0/3.0)) + 1.0e-3)
  y11ip  = 1.0
  y11ips = 1.0/((1+(99.0/44.0**6)*Zimp**6))
        
  y12ib  = (571.6*(1.0 + 0.84*Zimp + 7.8e-7*Zimp**5.15)*ft**(3.43*(1.0 + 0.012*Zimp))*&
            (1.0-ft)**(1.2*(1.0 + 1.6e-8*Zimp**5.1))/(1.0 + 500.0*ft**(8.0/3.0)) + 1.0e-3)
  y12ip  = 1.0
  y12ips = 1.0


  if (rotation_model.eq.2) then

     c1 = -1.4*ft**7 + 2.23*(1.0-0.31*ft)
     c2 = 2.8*(1.0-0.63*ft)
     c3 = 3.5*(1.0 - ft)/c2
     c4 = 4.0*ft
     c5 = 0.38*ft**4
     c6 = 3.95*(1.0 + 0.424*ft*(1.0 - 0.65*ft))
     fG = (1.0 + c1*Mzstar**c2)**(c3)*(1.0 + 0.2*Mzstar**c4)/(1.0 + c5*Mzstar**c6)
        
        
     c1 = 2.72*(1.0-0.91*ft)
     c2 = 2.98*(1.0-0.471*ft)
     c3 = 0.95*(1.0-ft)**4
     c4 = 4.0*ft
     c5 = 0.1314*ft**2.84 + 3.178*ft**11.4
     c6 = -9.38*(ft-0.5)**2 + 4.64
     fU = (c1*Mzstar**c2)*(1.0 + c3*Mzstar**c4)/(1.0 + c5*Mzstar**c6)

     ahbp = (0.135 + 2.647e-3*Zimp**1.464 + 3.478e-10*Zimp**5.347)*&
            ((1.0 + (3.0/(1.0 + 1.0e-7*Zimp**6))*(1.0/(1.0 + 1.2e5*ft**12))*Mzstar)/&
             (1.0 + (3.0/(1.0 + 1.0e-7*Zimp**6))*(1.208 - 4.46*ft + 4.394*ft**2)*Mzstar**2))

     c1 = 1.0 + 14.86*(1 - ft)**16.45 + 15.27*ft**7.4
     c2 = 0.77*(1.0 + 4.11*ft)
     c3 = 0.01*(1.0 + 359*ft**2.5 + 1078*ft**12)
     l1 = (1.0 + c1*Mzstar**c2)*exp(-c3*Mzstar**2)
            
     c1 = 6.39*(1.0 -ft)**15.8 + 0.1
     c2 = 0.943*(1.0 + 3.5*ft)
     l2 = (1.0 + c1*Mzstar**0.5)/(1.0 + c2*Mzstar**(10.0/3.0))
            
     l3 = 1.0/(1.0 + 2.0*ft*Mzstar)
            
     c1 = 6.13 + 28.18*ft**2.13 + 336.25*(1.0-ft)**11.65
     c2 = 0.5 + 9.55*ft**1.14*(1.0-ft)**1.42
     c3 = (0.0087 + 4.49*ft**3.48)/(1.0 + 0.873*ft**3.48)
     c4 = 3.6*(1.0 - 0.36*ft)
     l4 =  (1.0 + c1*Mzstar**c2)/(1.0 + c3*Mzstar**(c4))
            
     c1 = (1.0-ft)**8
     c2 = 113.5*ft**8.46
     c3 = 11*(1.0-ft)     
     l5 = (1.0 + c1*c2*Mzstar**(c3))/(1.0 + c2*Mzstar**(c3))
            
     l6 = (1.0 + 0.035*10.0*Mzstar**4)/(1.0 + 10.0*Mzstar**4)
            
     l7 = exp(-10.0*Mzstar**2)

     f2 = f2*l7
     f3 = f3*(1.0 + (1.0 + 1.86e6*ft**11.07*(1.0-ft)**7.36)*Mzstar**4)*exp(-0.8*Mzstar**2)

     fvv = fvv*l7
            
     y11zb  = y11zb*l1
     y11zp  = y11zp*l1/l2
     y11zps = y11zps*l1/(l2*l3)
            
     y12zb  = y12zb*l4
     y12zp  = y12zp*l4/l5
     y12zps = y12zps*l4/(l5*l6)
            
     y11ib  = y11ib*l1
     y11ip  = y11ip*l1
     y11ips = y11ips*l1
            
     y12ib  = y12ib*l7
     y12ip  = y12ip*l7
     y12ips = y12ips*l7


  else
     ahbp = (1.01579172 + -1.78923911e-3*Zimp)/(1.0 + 6.60170647e-13*Zimp**6.66398825)

  endif

        
end subroutine facs_lfs




subroutine K_VISC_lfs(nx, ni, nimp, Ti, wii, wimpimp, Zi, Zimp, Ai, Aimp, Tauii, &
                  Tauimpi, Tauiimp, Tauimpimp, eps2, R0, qmag, ft, &
                  y11zb, y11zp, y11zps, y11ib, y11ip, y11ips, &
                  y12zb, y12zp, y12zps, y12ib, y12ip, y12ips, &
                  K11a, K12a, K22a, K11i, K12i, K22i)

!*******************************************************************************
! Calculates the positive definite matrix of neoclassical viscosity coefficients
! in a fully analytical way by solving for the coefficients in the individual
! banana, plateau and Pfirsch-Schlüter collisionality regimes and using a
! rational approximation to interpolate between them, as well as evaluating
! kinetic integrals for a Maxwellian distribution
!*******************************************************************************
! INPUTS:
! -------
! nx --------> size of radial arrays [-]
! ni --------> main ion density [1/m^3]
! nimp ------> impurity density [1/m^3]
! Ti --------> main ion temperature [eV]
! wii -------> ion transit frequency [1/s]
! wimpimp ---> impurity transit frequecy [1/s]
! Zi --------> main ion charge [-]
! Zimp ------> impurity charge [-]
! Ai --------> main ion mass [-]
! Aimp ------> impurity mass [-]
! Tauii -----> ion-ion collision time [s]
! Tauimpi ---> impurity-ion collision time [s]
! Tauiimp ---> ion-impurity collision time [s]
! Tauimpimp -> impurity impurity collision time [s]
! eps2 ------> local inverse aspect ratio squared [-]
! R0 --------> major radius at magnetic axis [m]
! qmag ------> safety factor [-]
! ft --------> trapped particle fraction [-]
! yjksc -----> fitted factors of c = B, P, PS component of the j,k viscosity
!              coefficient of species s = i,z, from facs subroutine [-]
!*******************************************************************************
! OUTPUTS:
! --------
! K11a ------> (1,1) impurity viscosity coefficient [kg/(m*s)]
! K12a ------> (1,2)=(2,1) impurity viscosity coefficient [kg/(m*s)]
! K22a ------> (2,2) impurity viscosity coefficient [kg/(m*s)]
! K11i ------> (1,1) main ion viscosity coefficient [kg/(m*s)]
! K12i ------> (1,2)=(2,1) main ion viscosity coefficient [kg/(m*s)]
! K22i ------> (2,2) main ion viscosity coefficient [kg/(m*s)]
!*******************************************************************************

  implicit none

  ! mathematical constants
  double precision, parameter :: pi    = 4.d0*ATAN(1.d0)
  double precision, parameter :: sqrt2 = SQRT(2.d0)

  ! physical constants
  double precision, parameter :: mp         = 1.67262192369e-27          ! proton mass [kg]
  double precision, parameter :: q_e        = 1.602176634e-19            ! electron charge [C]

  integer :: nx
  double precision, dimension(nx) :: ni, nimp, Ti, wii, wimpimp, Tauii, Tauimpi, Tauiimp, Tauimpimp, eps2, qmag, ft
  double precision :: Zi, Ai, Aimp, R0, mimp, mi
  double precision, dimension(nx) :: Zimp
  double precision, dimension(nx) :: K11a, K12a, K22a, K11i, K12i, K22i
  double precision, dimension(nx) :: y11zb, y11zp, y11zps, y11ib, y11ip, y11ips
  double precision, dimension(nx) :: y12zb, y12zp, y12zps, y12ib, y12ip, y12ips

  double precision, dimension(nx) :: fac_a_P, fac_i_P, K11aP, K12aP, K22aP, K11iP, K12iP, K22iP
  double precision :: r00, r01, r11, xai, xia, x2ai, x2ia, xfac_ai, xfac_ia
  double precision :: qaa00, qaa01, qaa11, qai00, qia00, qai01, qia01, qai11, qia11
  double precision, dimension(nx) :: fac_qai_PS, fac_qia_PS
  double precision, dimension(nx) :: qa00, qi00, qa01, qi01, qa11, qi11, Qa, Qi
  double precision, dimension(nx) :: la11, la12, la22, li11, li12, li22
  double precision, dimension(nx) :: fac_imp_PS, fac_ion_PS
  double precision, dimension(nx) :: K11aPS, K12aPS, K22aPS, K11iPS, K12iPS, K22iPS
  double precision, dimension(nx) :: fac_B, nuDai, nuD2ai, nuD4ai
  double precision, dimension(nx) :: nuDia, nuD2ia, nuD4ia
  double precision, dimension(nx) :: nuDaa, nuD2aa, nuD4aa
  double precision, dimension(nx) :: nuDii, nuD2ii, nuD4ii
  double precision, dimension(nx) :: K11aB, K12aB, K22aB, K11iB, K12iB, K22iB


  mimp = Aimp*mp
  mi = Ai*mp


  ! Plateau regime

  fac_a_P = nimp*(q_e*Ti)*sqrt(pi)/(3.0*wimpimp)
  fac_i_P = ni*(q_e*Ti)*sqrt(pi)/(3.0*wii)

  K11aP = fac_a_P*2.0
  K12aP = fac_a_P*2.0*3.0
  K22aP = fac_a_P*2.0*3.0*4.0

  K11iP = fac_i_P*2.0
  K12iP = fac_i_P*2.0*3.0
  K22iP = fac_i_P*2.0*3.0*4.0

  ! Pfirsch-Schlüter regime

  r00 = 1.0/sqrt2
  r01 = 1.5/sqrt2
  r11 = 3.75/sqrt2

  xai = sqrt(Aimp/Ai)
  xia = 1.0/xai

  x2ai = xai**2
  x2ia = xia**2

  xfac_ai = (1+x2ai)**0.5
  xfac_ia = (1+x2ia)**0.5

  qaa00 = 8.0/2**1.5
  qaa01 = 15.0/2**2.5
  qaa11 = 132.5/2**3.5

  qai00 = (3.0+5.0*x2ai)/xfac_ai**3
  qia00 = (3.0+5.0*x2ia)/xfac_ia**3

  qai01 = 1.5*(3.0+7.0*x2ai)/xfac_ai**5
  qia01 = 1.5*(3.0+7.0*x2ia)/xfac_ia**5

  qai11 = (35.0*x2ai**3 + 38.5*x2ai**2 + 46.25*x2ai + 12.75)/xfac_ai**7
  qia11 = (35.0*x2ia**3 + 38.5*x2ia**2 + 46.25*x2ia + 12.75)/xfac_ia**7

  fac_qai_PS = (ni*Zi**2/(nimp*Zimp**2))
  fac_qia_PS = (nimp*Zimp**2/(ni*Zi**2))


  qa00 = fac_qai_PS*qai00 +  qaa00 - r00
  qi00 = fac_qia_PS*qia00 + qaa00 - r00

  qa01 = fac_qai_PS*qai01 + qaa01 - r01
  qi01 = fac_qia_PS*qia01 + qaa01 - r01

  qa11 = fac_qai_PS*qai11 + qaa11 - r11
  qi11 = fac_qia_PS*qia11 + qaa11 - r11


  Qa = 0.4*(qa00*qa11-qa01*qa01)
  Qi = 0.4*(qi00*qi11-qi01*qi01)


  la11 = qa11/Qa
  la12 = 3.5*(qa11+qa01)/Qa
  la22 = 12.25*(qa11+qa00+2*qa01)/Qa

  li11 = qi11/Qi
  li12 = 3.5*(qi11 + qi01)/Qi
  li22 = 12.25*(qi11+qi00 + 2.0*qi01)/Qi

  fac_imp_PS = nimp*(q_e*Ti)*Tauimpimp
  fac_ion_PS = ni*(q_e*Ti)*Tauii


  K11aPS = fac_imp_PS*la11
  K12aPS = fac_imp_PS*la12
  K22aPS = fac_imp_PS*la22

  K11iPS = fac_ion_PS*li11
  K12iPS = fac_ion_PS*li12
  K22iPS = fac_ion_PS*li22



  ! Banana regime

  fac_B = (ft/(1.0 - ft))*(2.0*R0**2*qmag**2/(3.0*eps2))

  ! Maxwellian integrals

  nuDai  = (xfac_ai + x2ai*log(xai/(1 + xfac_ai)))/Tauimpi
  nuD2ai = 1.0/(xfac_ai*Tauimpi)
  nuD4ai = 2.0*(1.0 + 1.25*x2ai)/(xfac_ai**3*Tauimpi)

  nuDia  = (xfac_ia + x2ia*log(xia/(1 + xfac_ia)))/Tauiimp
  nuD2ia = 1.0/(xfac_ia*Tauiimp)
  nuD4ia = 2.0*(1.0 + 1.25*x2ia)/(xfac_ia**3*Tauiimp)

  nuDaa  = (sqrt2 + log(1/(1 + sqrt2)))/Tauimpimp
  nuD2aa = 1.0/(sqrt2*Tauimpimp)
  nuD4aa = 4.5/(2.0**1.5*Tauimpimp)

  nuDii  = (sqrt2 + log(1/(1 + sqrt2)))/Tauii
  nuD2ii = 1.0/(sqrt2*Tauii)
  nuD4ii = 4.5/(2.0**1.5*Tauii)

  ! Banana regime viscosity coefficient

  K11aB = fac_B*nimp*mimp*(nuDai + nuDaa)
  K12aB = fac_B*nimp*mimp*(nuD2ai + nuD2aa)
  K22aB = fac_B*nimp*mimp*(nuD4ai + nuD4aa)

  K11iB = fac_B*ni*mi*(nuDia + nuDii)
  K12iB = fac_B*ni*mi*(nuD2ia + nuD2ii)
  K22iB = fac_B*ni*mi*(nuD4ia + nuD4ii)

  ! total viscosity coefficients, rational approximation for interpolation

  K11a = y11zb*K11aB/((1 + y11zb*K11aB/(y11zp*K11aP))*(1 + y11zp*K11aP/(y11zps*K11aPS)))
  K12a = y12zb*K12aB/((1 + y12zb*K12aB/(y12zp*K12aP))*(1 + y12zp*K12aP/(y12zps*K12aPS)))
  K22a = K22aB/((1.0 + K22aB/(K22aP))*(1.0 + K22aP/(K22aPS)))
        
  K11i = y11ib*K11iB/((1 + y11ib*K11iB/(y11ip*K11iP))*(1 + y11ip*K11iP/(y11ips*K11iPS)))
  K12i = y12ib*K12iB/((1 + y12ib*K12iB/(y12ip*K12iP))*(1 + y12ip*K12iP/(y12ips*K12iPS)))
  K22i = K22iB/((1.0 + K22iB/(K22iP))*(1.0 + K22iP/(K22iPS)))



  return
end subroutine K_VISC_lfs


subroutine flux_surf_geom(geom_type, ntheta_in, rmin_out, theta_out, R_out, Z_out, Jacobian_out, R_LFS_out)
    
  !============================================================================================!
  !
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
  !
  !============================================================================================!

  use const_inc, only: NA1, RTOR, time, tau, tstart
  use status_inc, only: AMETR, SHIF, SHIV
  use parameters_a2equil, only: equil_now
  use numerical_tools, only: qinterp

  implicit none

  
  ! inputs and outputs
  integer, intent(in) :: geom_type
  integer, intent(in) :: ntheta_in
  double precision, dimension(NA1), intent(out) :: rmin_out, R_LFS_out
  double precision, dimension(ntheta_in), intent(out) :: theta_out
  double precision, dimension(NA1, ntheta_in), intent(out) :: R_out, Z_out, Jacobian_out

  ! mathematical constants
  double precision, parameter :: pi = 4.d0*ATAN(1.d0)
  

  integer :: i,j ! dummy for loops
  integer :: idxmpi
  integer :: nrho_surf, nthe_surf, jmax

  ! intermediate variables
  double precision, allocatable, dimension(:)   :: rmin_equ 
  double precision, allocatable, dimension(:)   :: th0, th1
  double precision, allocatable, dimension(:,:) :: R_1, Z_1
  double precision, allocatable, dimension(:)   :: th2
  double precision, allocatable, dimension(:,:) :: R_2, Z_2
  double precision, allocatable, dimension(:,:) :: R_3, Z_3
  double precision, allocatable, dimension(:)   :: pf_eq, rho_eq
  double precision, dimension(NA1, ntheta_in) :: dRdr, dRdth, dZdr, dZdth, grr, grt, gtt
 
  !-----------------------------------------------------------------------------------------------!
  ! magnetic geometry

  rmin_out  = AMETR(1:NA1)     ! minor radius [m]

  ! flux surface contours from equilibrium
  nrho_surf = SIZE(equil_now%coord_sys%position%r, 1)
  nthe_surf = SIZE(equil_now%coord_sys%position%r, 2)
  allocate(pf_eq(nrho_surf), rho_eq(nrho_surf))

  ! allocate other quantities
  allocate(rmin_equ(nrho_surf))
  allocate(th0(nthe_surf), th1(nthe_surf))
  allocate(R_1(nrho_surf,nthe_surf), Z_1(nrho_surf,nthe_surf))
  allocate(th2(nthe_surf + 2))
  allocate(R_2(nrho_surf,nthe_surf+2), Z_2(nrho_surf,nthe_surf+2))
  allocate(R_3(nrho_surf,ntheta_in), Z_3(nrho_surf,ntheta_in))

  ! initial poloidal coordinate
  do j=1, nthe_surf
     th0(j) = atan2(equil_now%coord_sys%position%z(nrho_surf, j) - SHIV(NA1), equil_now%coord_sys%position%r(nrho_surf,j)-(RTOR+SHIF(NA1)))
  enddo

  ! find location where th0 is closest to -pi
  idxmpi = minloc(abs(th0-(-pi)),1)

  ! roll coordinates over this 
  th1(1:nthe_surf-idxmpi+1) = th0(idxmpi:)
  th1(nthe_surf-idxmpi+2:)  = th0(1:idxmpi-1)

  ! now that th1 is correct, roll R_1, Z_1
  if (time-tstart > tau) then
    R_1(:, 1:nthe_surf-idxmpi+1) = equil_now%coord_sys%position%r(:, idxmpi:)
    R_1(:, nthe_surf-idxmpi+2:)  = equil_now%coord_sys%position%r(:, 1:idxmpi-1)
    Z_1(:, 1:nthe_surf-idxmpi+1) = equil_now%coord_sys%position%z(:, idxmpi:)
    Z_1(:, nthe_surf-idxmpi+2:)  = equil_now%coord_sys%position%z(:, 1:idxmpi-1)
  endif
  ! finish the poloidal turn:
  th2(2:nthe_surf+1) = th1
  th2(1) = -pi
  th2(nthe_surf+2) = pi

  R_2(:,2:nthe_surf+1) = R_1
  Z_2(:,2:nthe_surf+1) = Z_1

  R_2(:,1)=0.5*(R_1(:,1) + R_1(:,nthe_surf))
  R_2(:,nthe_surf+2)=0.5*(R_1(:,1) + R_1(:,nthe_surf))

  Z_2(:,1)=0.5*(Z_1(:,1) + Z_1(:,nthe_surf))
  Z_2(:,nthe_surf+2)=0.5*(Z_1(:,1) + Z_1(:,nthe_surf))
  
  ! interpolate

  do j=1, ntheta_in-1
     theta_out(j) = -pi + 2.*pi*float(j-1)/float(ntheta_in-1)
  enddo
  theta_out(ntheta_in) = pi

  do i=1, nrho_surf
     rmin_equ(i) = max(0.5*(maxval(equil_now%coord_sys%position%r(i,:))-minval(equil_now%coord_sys%position%r(i,:))), 1.e-3)
  enddo

  if (geom_type.eq.2) then
     ! full flux surface geometry
     ! interpolate in theta
     do i=1, nrho_surf
        call qinterp(th2, R_2(i,:), nthe_surf + 2, theta_out, R_3(i,:), ntheta_in)
        call qinterp(th2, Z_2(i,:), nthe_surf + 2, theta_out, Z_3(i,:), ntheta_in)
     enddo
     
     ! interpolate in r
     do j=1, ntheta_in
        call qinterp(rmin_equ, R_3(:,j), nrho_surf, rmin_out, R_out(:,j), int(NA1))
        call qinterp(rmin_equ, Z_3(:,j), nrho_surf, rmin_out, Z_out(:,j), int(NA1))
     enddo

  elseif (geom_type .eq. 0) then
     ! circular geometry
     do j=1, ntheta_in
        R_out(:,j) = RTOR + SHIF(1:NA1) + rmin_out*cos(theta_out(j))
        Z_out(:,j) = SHIV(1:NA1) + rmin_out*sin(theta_out(j))
     enddo
     
  elseif (geom_type .eq. 1) then
     print *, "Miller still to do"
  else
     print *, "select a geom_type = 0, 1 or 2"
  endif

  ! Jacobian

  ! dR/dr
  do i=2,int(NA1)-1
     dRdr(i,:) = (R_out(i+1,:)-R_out(i-1,:))/(rmin_out(i+1)-rmin_out(i-1))
  enddo

  dRdr(1,:)  = (R_out(2,:)-R_out(1,:))/(rmin_out(2)-rmin_out(1))
  dRdr(int(NA1),:) = (R_out(int(NA1),:)-R_out(int(NA1)-1,:))/(rmin_out(int(NA1))-rmin_out(int(NA1)-1))

  ! dZ/dr
  do i=2,int(NA1)-1
     dZdr(i,:) = (Z_out(i+1,:)-Z_out(i-1,:))/(rmin_out(i+1)-rmin_out(i-1))
  enddo

  dZdr(1,:)  = (Z_out(2,:)-Z_out(1,:))/(rmin_out(2)-rmin_out(1))
  dZdr(int(NA1),:) = (Z_out(int(NA1),:)-Z_out(int(NA1)-1,:))/(rmin_out(int(NA1))-rmin_out(int(NA1)-1))

  ! dR/dtheta
  do j=2,ntheta_in-1
     dRdth(:,j) = (R_out(:,j+1)-R_out(:,j-1))/(theta_out(j+1)-theta_out(j-1))
  enddo

  dRdth(:,1)   = (R_out(:,2)-R_out(:,1))/(theta_out(2)-theta_out(1))
  dRdth(:,ntheta_in) = (R_out(:,ntheta_in)-R_out(:,ntheta_in-1))/(theta_out(ntheta_in)-theta_out(ntheta_in-1))

  ! dZ/dtheta
  do j=2,ntheta_in-1
     dZdth(:,j) = (Z_out(:,j+1)-Z_out(:,j-1))/(theta_out(j+1)-theta_out(j-1))
  enddo

  dZdth(:,1)   = (Z_out(:,2)-Z_out(:,1))/(theta_out(2)-theta_out(1))
  dZdth(:,ntheta_in) = (Z_out(:,ntheta_in)-Z_out(:,ntheta_in-1))/(theta_out(ntheta_in)-theta_out(ntheta_in-1)) 


  grr = dRdr**2 + dZdr**2
  grt = dRdr*dRdth + dZdr*dZdth
  gtt = dRdth**2 + dZdth**2

  Jacobian_out = R_out*sqrt(grr*gtt - grt**2)

  ! finally, LFS major radius

  R_LFS_out = RTOR + SHIF(1:NA1) + rmin_out

  deallocate(pf_eq, rho_eq, rmin_equ, th0, th1, R_1, Z_1, th2, R_2, Z_2, R_3, Z_3)

  
  return
end subroutine flux_surf_geom


subroutine lfs2fsa_impDV(geom_type, output_op, Zimp_in, Aimp_in, Dz_lfs_in, Vz_lfs_in, Dz_fsa_out, Vz_fsa_out)
    
  !============================================================================================!
  !
  ! Calculate the terms that transform the low field side (LFS) diffusive and convective
  ! transport coefficients of an impurity species into flux surface averaged (FSA) coefficients
  ! following the appendix of Angioni 2014 Nucl. Fusion 54 083028
  !
  ! --- D. Fajardo, May 2023
  !
  ! * INPUTS
  ! --------
  ! - geom_type --> 0: circular, 1: Miller (STILL TO DO), 2: full flux surface
  ! - output_op --> 0: only e0imp on D and V, 1: full e0imp and FV
  ! - Zimp_in ----> impurity charge (profile) [-]
  ! - Aimp_in ----> impurity mass [-]
  ! - Dz_lfs_in --> LFS impurity diffusion coefficient [m^2/s]
  ! - Vz_lfs_in --> LFS impurity convective velocity [m/s]
  !
  ! * OUTPUTS
  ! ---------
  ! - Dz_fsa_out -> FSA impurity diffusion coefficient [m^2/s]
  ! - Vz_fsa_out -> FSA impurity convective velocity [m/s]
  !
  !============================================================================================!

  use parameter_inc, only: NRD
  use const_inc,     only: NA1, RTOR, ZMJ, AMJ
  use status_inc,    only: TE, TI, ZEF, VTOR, SHIF

  implicit none

  
  ! inputs and outputs
  integer, intent(in) :: geom_type, output_op
  double precision, intent(in) :: Aimp_in
  double precision, dimension(NRD), intent(in) :: Zimp_in, Dz_lfs_in, Vz_lfs_in
  double precision, dimension(NRD), intent(out) :: Dz_fsa_out, Vz_fsa_out

  ! physical constants
  double precision, parameter :: q_e = 1.602176634e-19   ! electron charge [C]
  double precision, parameter :: m_p = 1.67262192369e-27 ! proton mass [kg]

  ! intermediate variables
  integer :: i,j,nr ! dummy for loops
  integer, parameter :: n_theta = 101 !301

  double precision, dimension(NA1) :: Z_imp, T_e, T_i, Z_eff, Mach_i_sq, Mz_star_sq
  double precision, dimension(NA1) :: r_min, R_LFS
  double precision, dimension(n_theta) :: theta
  double precision, dimension(NA1, n_theta) :: R_2D, Z_2D, Jacobian

  double precision, dimension(NA1, n_theta) :: exp_m_Eimp
  double precision, dimension(NA1) :: e0imp, FVimp, d_e0imp_dr

  !-----------------------------------------------------------------------------------------------!
  ! plasma profiles

  T_e        = 1.0e3*max(TE(1:NA1),1.e-5)                                    ! electron temperature in eV
  T_i        = 1.0e3*max(TI(1:NA1),1.e-5)                                    ! main ion temperature in eV  
  Z_eff      = ZEF(1:NA1)                                                    ! effective charge
  Z_imp      = Zimp_in(1:NA1)                                                ! impurity charge
  Mach_i_sq  = abs(VTOR(1:NA1))**2./(2.*q_e*T_i/(AMJ*m_p))                   ! main ion Mach number squared
  Mz_star_sq = Mach_i_sq*(Aimp_in/AMJ - (Z_imp/ZMJ)*Z_eff/(Z_eff + T_i/T_e)) ! effective impurity Mach number squared
 
  !-----------------------------------------------------------------------------------------------!
  ! magnetic geometry
  
  call flux_surf_geom(geom_type, n_theta, r_min, theta, R_2D, Z_2D, Jacobian, R_LFS)

  !-----------------------------------------------------------------------------------------------!

  nr = int(NA1)

  ! normalized energy exponential
  do i=1, nr
     exp_m_Eimp(i,:)  = exp(Mz_star_sq(i)*((R_2D(i,:)**2. - R_LFS(i)**2.)/RTOR**2.))
  enddo

  !-----------------------------------------------------------------------------------------------!
  ! take averages

  do i=1,nr
     call flux_surf_avg(n_theta, theta , exp_m_Eimp(i,:), Jacobian(i,:), e0imp(i))
  enddo

  do i=2,nr-1
     d_e0imp_dr(i) = (e0imp(i+1)-e0imp(i-1))/(r_min(i+1)-r_min(i-1))
  enddo

  d_e0imp_dr(1)  = (e0imp(2)-e0imp(1))/(r_min(2)-r_min(1))
  d_e0imp_dr(nr) = (e0imp(nr)-e0imp(nr-1))/(r_min(nr)-r_min(nr-1))

  FVimp = d_e0imp_dr/e0imp

  !-----------------------------------------------------------------------------------------------!
  ! outputs
  Dz_fsa_out(1:NA1) = Dz_lfs_in(1:NA1)/e0imp

  if (output_op.eq.0) then

     Vz_fsa_out(1:NA1) = Vz_lfs_in(1:NA1)/e0imp
     
  elseif (output_op.eq.1) then
     
     Vz_fsa_out(1:NA1) = Vz_lfs_in(1:NA1)/e0imp + Dz_fsa_out(1:NA1)*FVimp

  else
     print *, "select a correct output_op"
  endif
  
  return
end subroutine lfs2fsa_impDV



subroutine flux_surf_avg(nth,thetay,AF,JJ,Aavg)

  !============================================================================================!
  !
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
  !
  !============================================================================================!

  implicit none

  double precision, dimension(nth) :: AF, JJ
  double precision, dimension(nth) :: thetay
  double precision :: denom, Aavg
  integer   :: nth, ith

  Aavg = 0.
  denom = 0.

  do ith=2,nth
    denom = denom + 0.5*(thetay(ith)-thetay(ith-1))*(JJ(ith) + JJ(ith-1))
  enddo

  do ith=2,nth
    Aavg = Aavg+ 0.5*(thetay(ith)-thetay(ith-1))*(AF(ith)*JJ(ith) + AF(ith-1)*JJ(ith-1))
  enddo

  if (abs(denom).gt.0.) then
    Aavg = Aavg/denom
  else
    Aavg = sum(AF)/nth
  endif

  return
end subroutine flux_surf_avg
