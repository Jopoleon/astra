subroutine NEOCL4

!    Astra interface to Houlberg's code NCLASS
!
!             May-2004, A.Zolotukhin
!             Apr-2021 G. Tardini -> f90
!             Feb-2023 G. Tardini: removed work, using module nclass_mod instead
!---------------------------------------------------------------------
! NCLASS calculates the neoclassical transport properties of a multiple
!   species axisymmetric plasma using k_order parallel and radial force
!   balance equations for each species
! References:
!   Houlberg, Shaing, Hirshman, Zarnstorff, Phys Plasmas 4 (1997) 3230
!   Hirshman, Sigmar, Nucl Fusion 21 (1981) 1079
!   W.A.Houlberg 3/99
!---------------------------------------------------------------------
! Usage example in an ASTRA model:
!  - for electron density
!    DN = ... + "dn_e_nc(j)";
!    CN = ... + "cn_e_nc(j)";
!
!  - for electron heat transport
!    HE = ... + "xe_nc(j)";
!    CE = ... + "ce_nc(j)";
!
!  - for a density of impurities Z1 
!                                (equation F1 , for example)
!    DF1 = ... + "dn_imp1_nc(j)";
!    VF1 = ... + "cn_imp1_nc(j)";
!
!  - for a bootstrap current
!      DC=0.;  HC=0.;  XC=0.;
!                   CD=... + "jbs_nc(j)";  
!---------------------------------------------------------------------
! Below the net fluxes are understood as a flux G_i of the 
! particular species through the entire flux surface according to
!                 dn_i    d   
!                 ---- = ----[G_i] + Source_i
!                  dt     dV  
! It coincides with the Astra total fluxes QN, QF1, etc., as
!---------------------------------------------------------------------
! Shared, for each species:
!   - electrons, main ions, proton, D, T, He3, He4, Imp1, Imp2, Imp3
!  1) - particle radial net flux Gamma [1.e19/s]
!  2) - particle diffusion coefficient Dn [m**2/s]
!  3) - convective velocity Vn [m/s]
!  4) - radial net heat conduction flux q_cond [MW] 
!  5) - heat conductivity Chi [m**2/s]
!  6) - heat convective velocity V_heat [m/s]
!  7) - radial energy (conduction+convection) flux [MW] 
!  8) - bootstrap current on (p'/p), [MA/m**2]
!  9) - bootstrap current on (T'/T), [MA/m**2]
! 10) - poloidal flow velocity on outside midplane, [m/s]
!-----------------------------------------------
!               - Miscellanious parameters
! bootstrap current [MA/m**2]
! external current [MA/m**2]
! current conductivity [MS/m=1/(microOhm*m)]
! density of main ions as they are defined above
!---------------------------------------------------------------------

use const_inc, only: GP2, ABC, ROC, BTOR, RTOR, HRO, NA, NA1, &
   AMJ, AIM1, AIM2, AIM3, ZMJ
use status_inc, only: BDB0, B0DB2, BDB02, BMAXT, FOFB, IPOL, &
   ULON, ER, VRS, G11, &
   MU, ELON, SHIF, TE, TI, &
   NE, NHYDR, NDEUT, NTRIT, NHE3, NALF, ZIM1, ZIM2, ZIM3, NIZ1, NIZ2, NIZ3, NMAIN
use nclass_mod

implicit none

!  Physical and conversion constants
real, parameter :: z_coulomb=1.6022e-19, z_electronmass=9.1095e-31, &
   z_protonmass=1.6726e-27, z_j7kv=1.6022e-16, z_mu0=1.2566e-06

double precision :: YGRRdB2(1000),YNGRTHETA(1000),YFM(3,1000)

double precision :: y_grrho2, grti, y_den, yh
real, dimension(mx_ms) :: dq_s, vq_s
real p_eps
integer :: k_electron, k_mainion, k_proton, k_deuteron, k_triton, &
   k_he3, k_alpha, k_impZ1, k_impZ2, k_impZ3, narray
real :: ybbmax, ybbmax2, yftupper, yftlower
real rdum(8)
!Declaration of input to NCLASS
integer :: k_order, k_potato, m_i, m_z
real :: c_den, c_potb, c_potl
real :: p_b2, p_bm2, p_eb, p_fhat, p_fm(3), p_ft, p_grbm2, p_grphi, &
   p_gr2phi, p_ngrth
real, dimension(mx_mi) :: amu_i, grt_i, temp_i
real :: den_iz(mx_mi,mx_mz), fex_iz(3,mx_mi,mx_mz), grp_iz(mx_mi,mx_mz)
real :: YNMAIN, YNMAIN1
!Declaration of output from NCLASS
integer :: iflag, m_s
integer, dimension(mx_ms) :: jm_s, jz_s
real :: p_bsjb, p_etap, p_exjb
real :: calm_i(3,3,mx_mi)
real, dimension(3,3,mx_mi,mx_mi) :: caln_ii, capm_ii, capn_ii
real, dimension(mx_ms) :: bsjbp_s, bsjbt_s, dn_s, sqz_s, vn_s, veb_s, qeb_s, xi_s
real, dimension(5,mx_ms) :: gfl_s, qfl_s  
real, dimension(3, 3, mx_ms) :: upar_s, utheta_s, ymu_s
real, dimension(mx_ms, mx_ms) :: chip_ss, chit_ss, dp_ss, dt_ss

!Declaration of local variables
character(len=120) :: label
integer :: i, iza, j, jj, j1, k, im, nout, idum(8)

! Set control and radially independent data
!  k_order-order of v moments to be solved [-]
!         =2 u and p_q
!         =3 u, p_q, and u2
!         =else error
!  k_potato-option to include potato orbits [-]
!          =0 off
!          =else on
k_order  = 2
k_potato = 1
!  c_den-density cutoff below which species is ignored (/m**3)
!  c_potb-kappa(0)*Bt(0)/[2*q(0)**2] (T)
!  c_potl-q(0)*R(0) (m)
c_den  = 1.0e10
y_den  = 1.e-19*c_den                  ! c_den in ASTRA units
c_potb = -0.5*ELON(1)*BTOR*MU(1)**2
! Parameters:
!  mx_mi=9  - max number of isotopes
!  mx_mz=18 - max charge
!  mx_ms=40 - max number of species
!  m_i-number of isotopes (1<mi<mx_mi+1)
!  m_z-highest charge state of all species (0<mz<mx_mz+1)
!  grt_i(i)-temperature gradient of i (keV/rho)
!  grp_iz(i,z)-pressure gradient of i,z (keV/m**3/rho)

call ZBFAUX(yGRRdB2, yNGRTheta, YFM)

do j=1,NA
! Set radially dependent data
!  p_eps-inverse aspect ratio [-]
!  p_grphi-radial electric field Phi' (V/rho)
!  p_gr2phi-radial electric field gradient Psi'(Phi'/Psi')' (V/rho**2)
!  p_q-safety factor [-] (-1./MU(j))
!  p_eb-<E.B> (V*T/m)
!  p_b2-<B**2> (T**2)
!  p_bm2-<1/B**2> (/T**2)
!  p_fhat-mu_0*F/(dPsi/dr) (rho/m)
!  p_ft-trapped fraction [-]
!  p_grbm2-<grad(rho)**2/B**2> (rho**2/m**2/T**2)
!  p_ngrth-<n.grad(Theta)> (1/m)
!  c_potl   = -(RTOR+SHIF(j))/MU(1)
   k_electron = 0
   k_mainion  = 0
   k_proton   = 0
   k_deuteron = 0
   k_triton   = 0
   k_he3      = 0
   k_alpha    = 0
   k_impZ1    = 0
   k_impZ2    = 0
   k_impZ3    = 0

   YH = HRO

   c_potl   = -(RTOR+SHIF(j))/MU(1)
   y_grrho2 = G11(j)/VRS(j)

   p_eps    = HRO*j/ROC*ABC/(RTOR+SHIF(j))
   p_grphi  = ER(j)
   p_gr2phi = MU(j)*j*(ER(j+1)/MU(j+1)/((j+1)*YH) - ER(j)/MU(j)/(j*YH))
! Warning: The NCLASS version 1.2 returns NaN resistivity
!          if p_eb = 0.
!         p_eb     = max(.00001,-ULON(j)*BTOR/(GP2*RTOR))
   p_eb     = -ULON(j)*BTOR/(GP2*RTOR)
   if (abs(p_eb) .lt. 1.d-3) p_eb = sign(0.001,p_eb)
   p_b2     = BTOR**2*BDB02(j)
   p_bm2    = B0DB2(j)/BTOR**2
   p_fhat   = -IPOL(j)*RTOR/(MU(j)*j*YH)
! Trapped particle fraction
   ybbmax   = BTOR/BMAXT(j)*BDB0(j)      !->  <h>=<B/B_max>
   ybbmax2  = (BTOR/BMAXT(j))**2*BDB02(j)!->  <h2>=<(B/B_max)2>
   yftupper = 1.- ybbmax2/ybbmax**2* &
              (1.- sqrt(1. - ybbmax)*(1. + 0.5*ybbmax))
   yftlower = 1. - ybbmax2*(BMAXT(j)/BTOR)**2*FOFB(j)
   p_ft     = 0.75*yftupper + 0.25*yftlower
   p_grbm2  = yGRRdB2(j)   ! <grad(rho)**2/B**2>
   p_ngrth  = yNGRTheta(j) ! <n.grad(Theta)> 

!  p_fm(3)-poloidal moments of geometric factor for PS viscosity [-]
   do i=1,3
      p_fm(i)=0.0
   enddo
   do i=1,3
      p_fm(i) = YFM(i,j)
   enddo

   do j1=1,mx_mi
      do jj=1,mx_mz
         den_iz(j1,jj) = 0.
         grp_iz(j1,jj) = 0.
         do i=1,3
            fex_iz(i,j1,jj) = 0.
         enddo
      enddo
   enddo

   m_i = 1                ! Reserved for electrons
   k_electron = 1
   m_z = 1
   amu_i(m_i)  = 5.4463e-4
   temp_i(m_i) = TE(j)
   grt_i(m_i)  = (TE(j+1)-TE(j))/YH
   den_iz(m_i,m_z) = 1.e19*NE(j)
   grp_iz(m_i,m_z) = 1.e19/YH*(TE(j+1)*NE(j+1)-TE(j)*NE(j)) 
   grti = (TI(j+1)-TI(j))/YH

   if (NHYDR(j) .gt. y_den) then
      m_i = m_i+1
      k_proton = m_i
      m_z = 1
      amu_i(m_i) = 1.
      grt_i(m_i) = grti
      temp_i(m_i) = TI(j)
      den_iz(m_i,m_z)=1.e19*NHYDR(j)
      grp_iz(m_i,m_z) = 1.e19/YH*(TI(j+1)*NHYDR(j+1)-TI(j)*NHYDR(j))
   endif

   if (NDEUT(j) .gt. y_den) then
      m_i = m_i+1
      k_deuteron = m_i
      amu_i(m_i) = 2.
      grt_i(m_i) = grti
      temp_i(m_i) = TI(j)
      den_iz(m_i,m_z) = 1.e19*NDEUT(j)
      grp_iz(m_i,m_z) = 1.e19/YH*(TI(j+1)*NDEUT(j+1)-TI(j)*NDEUT(j)) 
   endif
   if (NTRIT(j) .gt. y_den) then
      m_i = m_i+1
      k_triton = m_i
      amu_i(m_i) = 3.
      grt_i(m_i) = grti
      temp_i(m_i) = TI(j)
      den_iz(m_i,m_z) = 1.e19*NTRIT(j)
      grp_iz(m_i,m_z) = 1.e19/YH*(TI(j+1)*NTRIT(j+1)-TI(j)*NTRIT(j))
   endif

   if (m_i .eq. 1 .and. mx_mz .lt. 2 .and. ZMJ .ge. 2 ) then
      write(*,*)">>> NEOCL: ion composition error:"
      write(*,*)" mx_mz < 2 and no hydrogen. Call ignored."
      return
   endif
! It is assumed here that the ion species are not specified explicitly
! in the model
! or/and
! is used to satisfy the quasineutrality condition
   YNMAIN = NMAIN(J)
   ni_nc(j) = YNMAIN
   YNMAIN1 = NMAIN(J+1)
   if (YNMAIN .gt. y_den) then
      m_i = m_i + 1
      k_mainion = m_i
      m_z = ZMJ
      amu_i(m_i) = AMJ
      temp_i(m_i) = TI(j)
      grt_i(m_i) = grti
      den_iz(m_i,m_z) = 1.e19*YNMAIN
      grp_iz(m_i,m_z) = 1.e19/YH*(TI(j+1)*YNMAIN1-TI(j)*YNMAIN) 
   endif
   if (NHE3(j) .gt. y_den) then
      m_i = m_i+1
      k_he3 = m_i
      amu_i(m_i) = 3.
      grt_i(m_i) = grti
      temp_i(m_i) = TI(j)
      m_z = 2
      den_iz(m_i,m_z) = 1.e19*NHE3(j)
      grp_iz(m_i,m_z) = 1.e19/YH*(TI(j+1)*NHE3(j+1)-TI(j)*NHE3(j)) 
   endif
   
   if (m_i .eq. mx_mi) goto 5

   if (NALF(j) .gt. y_den) then
      m_i = m_i+1
      k_alpha = m_i
      amu_i(m_i) = 4.
      grt_i(m_i) = grti
      temp_i(m_i) = TI(j)
      m_z = 2
      den_iz(m_i,m_z) = 1.e19*NALF(j)
      grp_iz(m_i,m_z) = 1.e19/YH*(TI(j+1)*NALF(j+1)-TI(j)*NALF(j)) 
   endif

! consider up to three impurity species
   jj = nint(ZIM1(j))
   if (jj <= mx_mz) then
      if (NIZ1(j) .gt. y_den) then
         m_i = m_i+1
         k_impZ1 = m_i
         amu_i(m_i) = AIM1
         grt_i(m_i) = grti
         temp_i(m_i) = TI(j)
         if (jj .gt. m_z) m_z = jj
         den_iz(m_i,jj) = 1.e19*NIZ1(j)
         grp_iz(m_i,jj) = 1.e19/YH*(TI(j+1)*NIZ1(j+1)-TI(j)*NIZ1(j)) 
      endif
      if (m_i .eq. mx_mi) goto 5
   endif

   jj = nint(ZIM2(j))
   if (jj <= mx_mz) then
      if (NIZ2(j) .gt. y_den) then
         m_i = m_i+1
         k_impZ2 = m_i
         amu_i(m_i) = AIM2
         grt_i(m_i) = grti
         temp_i(m_i) = TI(j)
         if (jj .gt. m_z) m_z = jj
         den_iz(m_i,jj) = 1.e19*NIZ2(j)
         grp_iz(m_i,jj) = 1.e19/YH*(TI(j+1)*NIZ2(j+1)-TI(j)*NIZ2(j)) 
      endif
      if (m_i .eq. mx_mi) goto 5
   endif

   jj = nint(ZIM3(j))
   if (jj <= mx_mz) then
      if (NIZ3(j) .gt. y_den) then
         m_i = m_i+1
         k_impZ3 = m_i
         amu_i(m_i) = AIM3
         grt_i(m_i) = grti
         temp_i(m_i) = TI(j)
         if (jj .gt. m_z) m_z = jj
         den_iz(m_i,jj) = 1.e19*NIZ3(j)
         grp_iz(m_i,jj) = 1.e19/YH*(TI(j+1)*NIZ3(j+1)-TI(j)*NIZ3(j)) 
      endif
   endif
  
 5 continue

   call NCLASS(k_order,k_potato,m_i,m_z,c_den,c_potb,c_potl,p_b2, &
               p_bm2,p_eb,p_fhat,p_fm,p_ft,p_grbm2,p_grphi,p_gr2phi, &
               p_ngrth,amu_i,grt_i,temp_i,den_iz,fex_iz,grp_iz, &
! output:
               m_s,jm_s,jz_s,p_bsjb,p_etap,p_exjb,calm_i,caln_ii,capm_ii, &
               capn_ii,bsjbp_s,bsjbt_s,dn_s,gfl_s,qfl_s,sqz_s,upar_s, &
               utheta_s,vn_s,veb_s,qeb_s,xi_s,ymu_s,chip_ss,chit_ss, &
               dp_ss,dt_ss,iflag)

   if (iflag.gt.0) then
      nout=6
      open(unit=nout, file='user')     
      if(iflag.eq.1) THEN
         label='ERROR:NCLASS-k_order must be 2 or 3, k_order='
         idum(1)=k_order
         call WRITE_LINE_IR(nout,label,1,idum,0,rdum,0)
      elseif(iflag.eq.2) THEN
         label='ERROR:NCLASS-require 1<m_i<mx_mi, m_i='
         idum(1)=m_i
         call WRITE_LINE_IR(nout,label,1,idum,0,rdum,0)
      elseif(iflag.eq.3) THEN
         label='ERROR:NCLASS-require 0<m_z<mx_mz, m_z='
         idum(1)=m_z
         call WRITE_LINE_IR(nout,label,1,idum,0,rdum,0)
      elseif(iflag.eq.4) THEN
         label='ERROR:NCLASS-require 0<m_s<mx_ms, m_s='
         idum(1)=m_s
         call WRITE_LINE_IR(nout,label,1,idum,0,rdum,0)
      elseif(iflag.eq.5) THEN
         label='ERROR:NCLASS-inversion of flow matrix failed'
         call WRITE_LINE(nout,label,0,0)
      endif
      return
   endif

! Output
! Electrons
   if (k_electron .gt. 0) then
! Total radial particle flux (electrons)
      call RARRAY_COPY(5,gfl_s(1,k_electron),1,rdum,1)
      rdum(6) = RARRAY_SUM(5,rdum,1)
      gamma_e_nc(j) = rdum(6)*VRS(j)*1.e-19
! Particle diffusion and velocity (electrons)
      im  = jm_s(k_electron)
      iza = IABS(jz_s(k_electron))
      dn_e_nc(j) = dn_s(k_electron)/y_grrho2
      cn_e_nc(j) = (vn_s(k_electron) + veb_s(k_electron) + &
                    gfl_s(5,k_electron)/den_iz(im,iza)) / y_grrho2
! Radial conduction flux (electrons)
      call RARRAY_COPY(5,qfl_s(1,k_electron),1,rdum,1)
      rdum(6) = RARRAY_SUM(5,rdum,1)
      qcond_e_nc(j) = rdum(6)*VRS(j)*1.e-6
! Heat conduction and velocity (electrons)
! Conduction total is sum of components
      rdum(1)=RARRAY_SUM(5,qfl_s(1,k_electron),1)
! Diagonal conductivity plus convective velocity
      dq_s(k_electron)=chit_ss(k_electron,k_electron) + &
                       chip_ss(k_electron,k_electron)
      vq_s(k_electron)=rdum(1)/(den_iz(im,iza)*z_j7kv*temp_i(im)) + &
                       dq_s(k_electron)*grt_i(im)/temp_i(im)
      xe_nc(j) = dq_s(k_electron)/y_grrho2
      ce_nc(j) = vq_s(k_electron)/y_grrho2
! Total radial energy flux (electrons)
      do k=1,5
         rdum(k)=qfl_s(k,k_electron) + &
                 2.5*gfl_s(k,k_electron)*temp_i(im)*z_j7kv
      enddo
      rdum(6)=RARRAY_SUM(5,rdum,1)
      qen_e_nc(j) = rdum(6)*VRS(j)*1.e-6
! Bootstrap current on p'/p (electrons)
      rdum(1)=bsjbp_s(k_electron)
      rdum(2)=bsjbt_s(k_electron)
      bs_pe_nc(j) = -1.e-6*rdum(1)/BTOR
      bs_te_nc(j) = -1.e-6*rdum(2)/BTOR
! Poloidal flow velociry on outside midplane (electrons)
      polflow_e_nc(j) = ( utheta_s(1,1,k_electron) +  &
                      utheta_s(1,2,k_electron) + &
                      utheta_s(1,3,k_electron)  ) * &
                      BTOR/(1.0 + p_eps)/p_fhat
   endif
!---  Main Ions  ---
   if (k_mainion .gt. 0) then
! Total radial particle flux (main ions)
      call RARRAY_COPY(5,gfl_s(1,k_mainion),1,rdum,1)
      rdum(6) = RARRAY_SUM(5,rdum,1)
      gamma_i_nc(j) = rdum(6)*VRS(j)*1.e-19
! Particle diffusion and velocity (main ions)
      im  = jm_s(k_mainion)
      iza = IABS(jz_s(k_mainion))
      dn_i_nc(j) = dn_s(k_mainion)/y_grrho2
      cn_i_nc(j) = (vn_s(k_mainion) + veb_s(k_mainion) + &
                    gfl_s(5,k_mainion)/den_iz(im,iza)) / y_grrho2
! Radial conduction flux (main ions)
      call RARRAY_COPY(5,qfl_s(1,k_mainion),1,rdum,1)
      rdum(6) = RARRAY_SUM(5,rdum,1)
      qcond_i_nc(j) = rdum(6)*VRS(j)*1.e-6
! Heat conduction and velocity (main ions)
! Conduction total is sum of components
      rdum(1)=RARRAY_SUM(5,qfl_s(1,k_mainion),1)
! Diagonal conductivity plus convective velocity
      dq_s(k_mainion)=chit_ss(k_mainion,k_mainion) + &
                      chip_ss(k_mainion,k_mainion)
      vq_s(k_mainion)=rdum(1)/(den_iz(im,iza)* z_j7kv*temp_i(im)) + &
                      dq_s(k_mainion)*grt_i(im)/temp_i(im)
      xi_nc(j) = dq_s(k_mainion)/y_grrho2
      ci_nc(j) = vq_s(k_mainion)/y_grrho2
! Total radial energy flux (main ions)
      do k=1,5
         rdum(k)=qfl_s(k,k_mainion) + &
                 2.5*gfl_s(k,k_mainion)*temp_i(im)*z_j7kv
      enddo
      rdum(6)=RARRAY_SUM(5,rdum,1)
      qen_i_nc(j) = rdum(6)*VRS(j)*1.e-6
! Bootstrap current on p'/p (main ions)
      rdum(1)=bsjbp_s(k_mainion)
      rdum(2)=bsjbt_s(k_mainion)
      bs_pi_nc(j) = -1.e-6*rdum(1)/BTOR
      bs_ti_nc(j) = -1.e-6*rdum(2)/BTOR
! Poloidal flow velociry on outside midplane (main ions)
      polflow_i_nc(j) = ( utheta_s(1,1,k_mainion) +  &
                      utheta_s(1,2,k_mainion) + &
                      utheta_s(1,3,k_mainion)  ) * &
                      BTOR/(1.0 + p_eps)/p_fhat
   endif

!--- Protons ---
   if (k_proton .gt. 0) then
! Total radial particle flux (protons  )
      call RARRAY_COPY(5,gfl_s(1,k_proton),1,rdum,1)
      rdum(6) = RARRAY_SUM(5,rdum,1)
      gamma_p_nc(j) = rdum(6)*VRS(j)*1.e-19
! Particle diffusion and velocity (protons  )
      im  = jm_s(k_proton)
      iza = IABS(jz_s(k_proton))
      dn_p_nc(j) = dn_s(k_proton)/y_grrho2
      cn_p_nc(j) = (vn_s(k_proton) + veb_s(k_proton) + &
                     gfl_s(5,k_proton)/den_iz(im,iza))  / y_grrho2
! Radial conduction flux (protons  )
      call RARRAY_COPY(5,qfl_s(1,k_proton),1,rdum,1)
      rdum(6) = RARRAY_SUM(5,rdum,1)
      qcond_p_nc(j) = rdum(6)*VRS(j)*1.e-6
! Heat conduction and velocity (protons  )
! Conduction total is sum of components
      rdum(1)=RARRAY_SUM(5,qfl_s(1,k_proton),1)
! Diagonal conductivity plus convective velocity
      dq_s(k_proton)=chit_ss(k_proton,k_proton) + chip_ss(k_proton,k_proton)
      vq_s(k_proton)=rdum(1)/(den_iz(im,iza)* z_j7kv*temp_i(im)) + &
                     dq_s(k_proton)*grt_i(im)/temp_i(im)
      xp_nc(j) = dq_s(k_proton)/y_grrho2
      cp_nc(j) = vq_s(k_proton)/y_grrho2
! Total radial energy flux (protons  )
      do k=1,5
         rdum(k)=qfl_s(k,k_proton) + &
                 2.5*gfl_s(k,k_proton)*temp_i(im)*z_j7kv
      enddo
      rdum(6)=RARRAY_SUM(5,rdum,1)
      qen_p_nc(j) = rdum(6)*VRS(j)*1.e-6
! Bootstrap current on p'/p (protons  )
      rdum(1)=bsjbp_s(k_proton)
      rdum(2)=bsjbt_s(k_proton)
      bs_pp_nc(j) = -1.e-6*rdum(1)/BTOR
      bs_tp_nc(j) = -1.e-6*rdum(2)/BTOR
! Poloidal flow velociry on outside midplane (protons)
      polflow_p_nc(j) = ( utheta_s(1,1,k_proton) +  &
                      utheta_s(1,2,k_proton) + &
                      utheta_s(1,3,k_proton)  ) * &
                      BTOR/(1.0 + p_eps)/p_fhat
   endif

!--- Deuterons ---
   if (k_deuteron .gt. 0) then
! Total radial particle flux (deuterons  )
      call RARRAY_COPY(5,gfl_s(1,k_deuteron),1,rdum,1)
      rdum(6) = RARRAY_SUM(5,rdum,1)
      gamma_d_nc(j) = rdum(6)*VRS(j)*1.e-19
! Particle diffusion and velocity (deuterons  )
      im  = jm_s(k_deuteron)
      iza = IABS(jz_s(k_deuteron))
      dn_d_nc(j) = dn_s(k_deuteron)/y_grrho2
      cn_d_nc(j) = (vn_s(k_deuteron) + veb_s(k_deuteron) + &
                     gfl_s(5,k_deuteron)/den_iz(im,iza)) / y_grrho2
! Radial conduction flux (deuterons  )
      call RARRAY_COPY(5,qfl_s(1,k_deuteron),1,rdum,1)
      rdum(6) = RARRAY_SUM(5,rdum,1)
      qcond_d_nc(j) = rdum(6)*VRS(j)*1.e-6
! Heat conduction and velocity (deuterons  )
! Conduction total is sum of components
      rdum(1)=RARRAY_SUM(5,qfl_s(1,k_deuteron),1)
! Diagonal conductivity plus convective velocity
      dq_s(k_deuteron)=chit_ss(k_deuteron,k_deuteron) + &
                       chip_ss(k_deuteron,k_deuteron)
      vq_s(k_deuteron)=rdum(1)/(den_iz(im,iza)* z_j7kv*temp_i(im)) + &
                       dq_s(k_deuteron)*grt_i(im)/temp_i(im)
      xd_nc(j) = dq_s(k_deuteron)/y_grrho2
      cd_nc(j) = vq_s(k_deuteron)/y_grrho2
! Total radial energy flux (deuterons  )
      do k=1,5
         rdum(k)=qfl_s(k,k_deuteron) + &
                 2.5*gfl_s(k,k_deuteron)*temp_i(im)*z_j7kv
      enddo
      rdum(6)=RARRAY_SUM(5,rdum,1)
      qen_d_nc(j) = rdum(6)*VRS(j)*1.e-6
! Bootstrap current on p'/p (deuterons  )
      rdum(1)=bsjbp_s(k_deuteron)
      rdum(2)=bsjbt_s(k_deuteron)
      bs_pd_nc(j) = -1.e-6*rdum(1)/BTOR
      bs_td_nc(j) = -1.e-6*rdum(2)/BTOR
! Poloidal flow velociry on outside midplane (deuterons)
      polflow_d_nc(j) = ( utheta_s(1,1,k_deuteron) +  &
                      utheta_s(1,2,k_deuteron) + &
                      utheta_s(1,3,k_deuteron)  ) * &
                      BTOR/(1.0 + p_eps)/p_fhat
   endif

!--- Tritons ---
   if (k_triton .gt. 0) then
! Total radial particle flux (tritons  )
      call RARRAY_COPY(5,gfl_s(1,k_triton),1,rdum,1)
      rdum(6) = RARRAY_SUM(5,rdum,1)
      gamma_t_nc(j) = rdum(6)*VRS(j)*1.e-19
! Particle diffusion and velocity (tritons  )
      im  = jm_s(k_triton)
      iza = IABS(jz_s(k_triton))
      dn_t_nc(j) = dn_s(k_triton)/y_grrho2
      cn_t_nc(j) = (vn_s(k_triton) + veb_s(k_triton) + &
                     gfl_s(5,k_triton)/den_iz(im,iza)) / y_grrho2
! Radial conduction flux (tritons  )
      call RARRAY_COPY(5,qfl_s(1,k_triton),1,rdum,1)
      rdum(6) = RARRAY_SUM(5,rdum,1)
      qcond_t_nc(j) = rdum(6)*VRS(j)*1.e-6
! Heat conduction and velocity (tritons  )
! Conduction total is sum of components
      rdum(1)=RARRAY_SUM(5,qfl_s(1,k_triton),1)
! Diagonal conductivity plus convective velocity
      dq_s(k_triton)=chit_ss(k_triton,k_triton) + &
                     chip_ss(k_triton,k_triton)
      vq_s(k_triton)=rdum(1)/(den_iz(im,iza)* z_j7kv*temp_i(im)) + &
                     dq_s(k_triton)*grt_i(im)/temp_i(im)
      xt_nc(j) = dq_s(k_triton)/y_grrho2
      ct_nc(j) = vq_s(k_triton)/y_grrho2
! Total radial energy flux (tritons  )
      do k=1,5
         rdum(k)=qfl_s(k,k_triton) + &
                 2.5*gfl_s(k,k_triton)*temp_i(im)*z_j7kv
      enddo
      rdum(6)=RARRAY_SUM(5,rdum,1)
      qen_t_nc(j) = rdum(6)*VRS(j)*1.e-6
! Bootstrap current on p'/p (tritons  )
      rdum(1)=bsjbp_s(k_triton)
      rdum(2)=bsjbt_s(k_triton)
      bs_pt_nc(j) = -1.e-6*rdum(1)/BTOR
      bs_tt_nc(j) = -1.e-6*rdum(2)/BTOR
! Poloidal flow velociry on outside midplane (tritons)
      polflow_t_nc(j) = ( utheta_s(1,1,k_triton) +  &
                      utheta_s(1,2,k_triton) + &
                      utheta_s(1,3,k_triton)  ) * &
                      BTOR/(1.0 + p_eps)/p_fhat
   endif

!--- He3 Particles ---
   if (k_he3 .gt. 0) then
! Total radial particle flux (He3 particles)
      call RARRAY_COPY(5,gfl_s(1,k_he3),1,rdum,1)
      rdum(6) = RARRAY_SUM(5,rdum,1)
      gamma_he3_nc = rdum(6)*VRS(j)*1.e-19
! Particle diffusion and velocity (He3 particles)
      im  = jm_s(k_he3)
      iza = IABS(jz_s(k_he3))
      dn_he3_nc(j) = dn_s(k_he3)/y_grrho2
      cn_he3_nc(j) = (vn_s(k_he3)  + veb_s(k_he3) + &
                     gfl_s(5,k_he3)/den_iz(im,iza))  / y_grrho2
! Radial conduction flux (He3 particles)
      call RARRAY_COPY(5,qfl_s(1,k_he3),1,rdum,1)
      rdum(6) = RARRAY_SUM(5,rdum,1)
      qcond_he3_nc(j) = rdum(6)*VRS(j)*1.e-6
! Heat conduction and velocity (He3 particles)
! Conduction total is sum of components
      rdum(1)=RARRAY_SUM(5,qfl_s(1,k_he3),1)
! Diagonal conductivity plus convective velocity
      dq_s(k_he3)=chit_ss(k_he3,k_he3) + chip_ss(k_he3,k_he3)
      vq_s(k_he3)=rdum(1)/(den_iz(im,iza)* z_j7kv*temp_i(im)) + &
                  dq_s(k_he3)*grt_i(im)/temp_i(im)
      xhe3_nc(j) = dq_s(k_he3)/y_grrho2
      che3_nc(j) = vq_s(k_he3)/y_grrho2
! Total radial energy flux (He3 particles)
      do k=1,5
         rdum(k)=qfl_s(k,k_he3) + &
                 2.5*gfl_s(k,k_he3)*temp_i(im)*z_j7kv
      enddo
      rdum(6)=RARRAY_SUM(5,rdum,1)
      qen_he3_nc(j) = rdum(6)*VRS(j)*1.e-6
! Bootstrap current on p'/p (He3 particles)
      rdum(1)=bsjbp_s(k_he3)
      rdum(2)=bsjbt_s(k_he3)
      bs_phe3_nc(j) = -1.e-6*rdum(1)/BTOR
      bs_the3_nc(j) = -1.e-6*rdum(2)/BTOR
! Poloidal flow velociry on outside midplane (He3 particles)
      polflow_he3_nc(j) = ( utheta_s(1,1,k_he3) +  &
                      utheta_s(1,2,k_he3) + &
                      utheta_s(1,3,k_he3)  ) * &
                      BTOR/(1.0 + p_eps)/p_fhat
   endif
!--- Alpha Particles ---
   if (k_alpha .gt. 0) then
! Total radial particle flux (alpha particles)
      call RARRAY_COPY(5,gfl_s(1,k_alpha),1,rdum,1)
      rdum(6) = RARRAY_SUM(5,rdum,1)
      gamma_he4_nc(j) = rdum(6)*VRS(j)*1.e-19
! Particle diffusion and velocity (alpha particles)
      im  = jm_s(k_alpha)
      iza = IABS(jz_s(k_alpha))
      dn_he4_nc(j) = dn_s(k_alpha)/y_grrho2
      cn_he4_nc(j) = (vn_s(k_alpha) + veb_s(k_alpha) + &
                     gfl_s(5,k_alpha)/den_iz(im,iza)) / y_grrho2
! Radial conduction flux (alpha particles)
      call RARRAY_COPY(5,qfl_s(1,k_alpha),1,rdum,1)
      rdum(6) = RARRAY_SUM(5,rdum,1)
      qcond_he4_nc(j) = rdum(6)*VRS(j)*1.e-6
! Heat conduction and velocity (alpha particles)
! Conduction total is sum of components
      rdum(1)=RARRAY_SUM(5,qfl_s(1,k_alpha),1)
! Diagonal conductivity plus convective velocity
      dq_s(k_alpha)=chit_ss(k_alpha,k_alpha) + chip_ss(k_alpha,k_alpha)
      vq_s(k_alpha)=rdum(1)/(den_iz(im,iza)* z_j7kv*temp_i(im)) + &
                    dq_s(k_alpha)*grt_i(im)/temp_i(im)
      xhe4_nc(j) = dq_s(k_alpha)/y_grrho2
      che4_nc(j) = vq_s(k_alpha)/y_grrho2
! Total radial energy flux (alpha particles)
      do k=1,5
         rdum(k)=qfl_s(k,k_alpha) + &
                 2.5*gfl_s(k,k_alpha)*temp_i(im)*z_j7kv
      enddo
      rdum(6)=RARRAY_SUM(5,rdum,1)
      qen_he4_nc(j) = rdum(6)*VRS(j)*1.e-6
! Bootstrap current on p'/p (alpha particles)
      rdum(1)=bsjbp_s(k_alpha)
      rdum(2)=bsjbt_s(k_alpha)
      bs_phe4_nc(j) = -1.e-6*rdum(1)/BTOR
      bs_the4_nc(j) = -1.e-6*rdum(2)/BTOR
! Poloidal flow velociry on outside midplane (alpha particles)
      polflow_he4_nc(j) = ( utheta_s(1,1,k_alpha) +  &
                      utheta_s(1,2,k_alpha) + &
                      utheta_s(1,3,k_alpha)  ) * &
                      BTOR/(1.0 + p_eps)/p_fhat
   endif
!--- Impurity Z1 ---
   if (k_impZ1 .gt. 0) then
! Total radial particle flux (impZ1 particles)
      call RARRAY_COPY(5,gfl_s(1,k_impZ1),1,rdum,1)
      rdum(6) = RARRAY_SUM(5,rdum,1)
      gamma_imp1_nc(j) = rdum(6)*VRS(j)*1.e-19
! Particle diffusion and velocity (impZ1 particles)
      im  = jm_s(k_impZ1)
      iza = IABS(jz_s(k_impZ1))
      dn_imp1_nc(j) = dn_s(k_impZ1)/y_grrho2
      cn_imp1_nc(j) = (vn_s(k_impZ1) + veb_s(k_impZ1) + &
                     gfl_s(5,k_impZ1)/den_iz(im,iza)) / y_grrho2
! Radial conduction flux (impZ1 particles)
      call RARRAY_COPY(5,qfl_s(1,k_impZ1),1,rdum,1)
      rdum(6) = RARRAY_SUM(5,rdum,1)
      qcond_imp1_nc(j) = rdum(6)*VRS(j)*1.e-6
! Heat conduction and velocity (impZ1 particles)
! Conduction total is sum of components
      rdum(1)=RARRAY_SUM(5,qfl_s(1,k_impZ1),1)
! Diagonal conductivity plus convective velocity
      dq_s(k_impZ1)=chit_ss(k_impZ1,k_impZ1) + chip_ss(k_impZ1,k_impZ1)
      vq_s(k_impZ1)=rdum(1)/(den_iz(im,iza)* z_j7kv*temp_i(im)) + &
                    dq_s(k_impZ1)*grt_i(im)/temp_i(im)
      ximp1_nc(j) = dq_s(k_impZ1)/y_grrho2
      cimp1_nc(j) = vq_s(k_impZ1)/y_grrho2
! Total radial energy flux (impZ1 particles)
      do k=1,5
         rdum(k)=qfl_s(k,k_impZ1) + &
                 2.5*gfl_s(k,k_impZ1)*temp_i(im)*z_j7kv
      enddo
      rdum(6)=RARRAY_SUM(5,rdum,1)
      qen_imp1_nc(j) = rdum(6)*VRS(j)*1.e-6
! Bootstrap current on p'/p (impZ1 particles)
      rdum(1)=bsjbp_s(k_impZ1)
      rdum(2)=bsjbt_s(k_impZ1)
      bs_pimp1_nc(j) = -1.e-6*rdum(1)/BTOR
      bs_timp1_nc(j) = -1.e-6*rdum(2)/BTOR
! Poloidal flow velociry on outside midplane (impZ1 particles)
      polflow_imp1_nc(j) = ( utheta_s(1,1,k_impZ1) +  &
                      utheta_s(1,2,k_impZ1) + &
                      utheta_s(1,3,k_impZ1)  ) * &
                      BTOR/(1.0 + p_eps)/p_fhat
   endif
!--- Impurity Z2 ---
   if (k_impZ2 .gt. 0) then
! Total radial particle flux (impZ2 particles)
      call RARRAY_COPY(5,gfl_s(1,k_impZ2),1,rdum,1)
      rdum(6) = RARRAY_SUM(5,rdum,1)
      gamma_imp2_nc(j) = rdum(6)*VRS(j)*1.e-19
! Particle diffusion and velocity (impZ2 particles)
      im  = jm_s(k_impZ2)
      iza = IABS(jz_s(k_impZ2))
      dn_imp2_nc(j) = dn_s(k_impZ2)/y_grrho2
      cn_imp2_nc(j) = (vn_s(k_impZ2) + veb_s(k_impZ2) + &
                     gfl_s(5,k_impZ2)/den_iz(im,iza))  / y_grrho2
! Radial conduction flux (impZ2 particles)
      call RARRAY_COPY(5,qfl_s(1,k_impZ2),1,rdum,1)
      rdum(6) = RARRAY_SUM(5,rdum,1)
      qcond_imp2_nc(j) = rdum(6)*VRS(j)*1.e-6
! Heat conduction and velocity (impZ2 particles)
! Conduction total is sum of components
      rdum(1)=RARRAY_SUM(5,qfl_s(1,k_impZ2),1)
! Diagonal conductivity plus convective velocity
      dq_s(k_impZ2)=chit_ss(k_impZ2,k_impZ2) + chip_ss(k_impZ2,k_impZ2)
      vq_s(k_impZ2)=rdum(1)/(den_iz(im,iza)* z_j7kv*temp_i(im)) + &
                    dq_s(k_impZ2)*grt_i(im)/temp_i(im)
      ximp2_nc(j) = dq_s(k_impZ2)/y_grrho2
      cimp2_nc(j) = vq_s(k_impZ2)/y_grrho2
! Total radial energy flux (impZ2 particles)
      do k=1,5
         rdum(k)=qfl_s(k,k_impZ2) + &
                 2.5*gfl_s(k,k_impZ2)*temp_i(im)*z_j7kv
      enddo
      rdum(6)=RARRAY_SUM(5,rdum,1)
      qen_imp2_nc(j) = rdum(6)*VRS(j)*1.e-6
! Bootstrap current on p'/p (impZ2 particles)
      rdum(1)=bsjbp_s(k_impZ2)
      rdum(2)=bsjbt_s(k_impZ2)
      bs_pimp2_nc(j) = -1.e-6*rdum(1)/BTOR
      bs_timp2_nc(j) = -1.e-6*rdum(2)/BTOR
! Poloidal flow velociry on outside midplane (impZ2 particles)
      polflow_imp2_nc(j) = ( utheta_s(1,1,k_impZ2) +  &
                      utheta_s(1,2,k_impZ2) + &
                      utheta_s(1,3,k_impZ2)  ) * &
                      BTOR/(1.0 + p_eps)/p_fhat
   endif
!--- Impurity Z3 ---
   if (k_impZ3 .gt. 0) then
! Total radial particle flux (impZ3 particles)
      call RARRAY_COPY(5,gfl_s(1,k_impZ3),1,rdum,1)
      rdum(6) = RARRAY_SUM(5,rdum,1)
      gamma_imp3_nc(j) = rdum(6)*VRS(j)*1.e-19
! Particle diffusion and velocity (impZ3 particles)
      im  = jm_s(k_impZ3)
      iza = IABS(jz_s(k_impZ3))
      dn_imp3_nc(j) = dn_s(k_impZ3)/y_grrho2
      cn_imp3_nc(j) = (vn_s(k_impZ3) + veb_s(k_impZ3) + &
                     gfl_s(5,k_impZ3)/den_iz(im,iza)) / y_grrho2
! Radial conduction flux (impZ3 particles)
      call RARRAY_COPY(5,qfl_s(1,k_impZ3),1,rdum,1)
      rdum(6) = RARRAY_SUM(5,rdum,1)
      qcond_imp3_nc(j) = rdum(6)*VRS(j)*1.e-6
! Heat conduction and velocity (impZ3 particles)
! Conduction total is sum of components
      rdum(1)=RARRAY_SUM(5,qfl_s(1,k_impZ3),1)
! Diagonal conductivity plus convective velocity
      dq_s(k_impZ3)=chit_ss(k_impZ3,k_impZ3) + chip_ss(k_impZ3,k_impZ3)
      vq_s(k_impZ3)=rdum(1)/(den_iz(im,iza)* z_j7kv*temp_i(im)) + &
                    dq_s(k_impZ3)*grt_i(im)/temp_i(im)
      ximp3_nc(j) = dq_s(k_impZ3)/y_grrho2
      cimp3_nc(j) = vq_s(k_impZ3)/y_grrho2
! Total radial energy flux (impZ3 particles)
      do k=1,5
         rdum(k)=qfl_s(k,k_impZ3) + &
                 2.5*gfl_s(k,k_impZ3)*temp_i(im)*z_j7kv
      enddo
      rdum(6)=RARRAY_SUM(5,rdum,1)
      qen_imp3_nc(j) = rdum(6)*VRS(j)*1.e-6
! Bootstrap current on p'/p (impZ3 particles)
      rdum(1)=bsjbp_s(k_impZ3)
      rdum(2)=bsjbt_s(k_impZ3)
      bs_pimp3_nc(j) = -1.e-6*rdum(1)/BTOR
      bs_timp3_nc(j) = -1.e-6*rdum(2)/BTOR
! Poloidal flow velociry on outside midplane (impZ3 particles)
      polflow_imp3_nc(j) = ( utheta_s(1,1,k_impZ3) +  &
                      utheta_s(1,2,k_impZ3) + &
                      utheta_s(1,3,k_impZ3)  ) * &
                      BTOR/(1.0 + p_eps)/p_fhat
   endif

!--- Miscellaneous Parameters ---

! Bootstrap current
   rdum(1)=p_bsjb
   jbs_nc(j)  = -1.e-6*rdum(1)/BTOR
! External current
   rdum(2)=p_exjb
   jext_nc(j)  = -1.e-6*rdum(2)/BTOR
! Current conductivity
   cc_nc(j)  = 1.e-6/p_etap

   if (isnan(cn_imp1_nc(j))) cn_imp1_nc(j) = 0.
   if (isnan(cn_imp2_nc(j))) cn_imp2_nc(j) = 0.
   if (isnan(cn_imp3_nc(j))) cn_imp3_nc(j) = 0.
   if (isnan(dn_imp1_nc(j))) dn_imp1_nc(j) = 0.
   if (isnan(dn_imp2_nc(j))) dn_imp2_nc(j) = 0.
   if (isnan(dn_imp3_nc(j))) dn_imp3_nc(j) = 0.

enddo ! Main radial loop

! Boundary value

gamma_e_nc(NA1) = gamma_e_nc(NA)
dn_e_nc(NA1)    = dn_e_nc(NA)
cn_e_nc(NA1)    = cn_e_nc(NA)
qcond_e_nc(NA1) = qcond_e_nc(NA)
xe_nc(NA1)      = xe_nc(NA)
ce_nc(NA1)      = ce_nc(NA)
qen_e_nc(NA1)   = qen_e_nc(NA)
bs_pe_nc(NA1)   = bs_pe_nc(NA)
bs_te_nc(NA1)   = bs_te_nc(NA)
polflow_e_nc(NA1) = polflow_e_nc(NA)

gamma_i_nc(NA1) = gamma_i_nc(NA)
dn_i_nc(NA1)    = dn_i_nc(NA)
cn_i_nc(NA1)    = cn_i_nc(NA)
qcond_i_nc(NA1) = qcond_i_nc(NA)
xi_nc(NA1)      = xi_nc(NA)
ci_nc(NA1)      = ci_nc(NA)
qen_i_nc(NA1)   = qen_i_nc(NA)
bs_pi_nc(NA1)   = bs_pi_nc(NA)
bs_ti_nc(NA1)   = bs_ti_nc(NA)
polflow_i_nc(NA1) = polflow_i_nc(NA)

gamma_p_nc(NA1) = gamma_p_nc(NA)
dn_p_nc(NA1)    = dn_p_nc(NA)
cn_p_nc(NA1)    = cn_p_nc(NA)
qcond_p_nc(NA1) = qcond_p_nc(NA)
xp_nc(NA1)      = xp_nc(NA)
cp_nc(NA1)      = cp_nc(NA)
qen_p_nc(NA1)   = qen_p_nc(NA)
bs_pp_nc(NA1)   = bs_pp_nc(NA)
bs_tp_nc(NA1)   = bs_tp_nc(NA)
polflow_p_nc(NA1) = polflow_p_nc(NA)

gamma_d_nc(NA1) = gamma_d_nc(NA)
dn_d_nc(NA1)    = dn_d_nc(NA)
cn_d_nc(NA1)    = cn_d_nc(NA)
qcond_d_nc(NA1) = qcond_d_nc(NA)
xd_nc(NA1)      = xd_nc(NA)
cd_nc(NA1)      = cd_nc(NA)
qen_d_nc(NA1)   = qen_d_nc(NA)
bs_pd_nc(NA1)   = bs_pd_nc(NA)
bs_td_nc(NA1)   = bs_td_nc(NA)
polflow_d_nc(NA1) = polflow_d_nc(NA)

gamma_t_nc(NA1) = gamma_t_nc(NA)
dn_t_nc(NA1)    = dn_t_nc(NA)
cn_t_nc(NA1)    = cn_t_nc(NA)
qcond_t_nc(NA1) = qcond_t_nc(NA)
xt_nc(NA1)      = xt_nc(NA)
ct_nc(NA1)      = ct_nc(NA)
qen_t_nc(NA1)   = qen_t_nc(NA)
bs_pt_nc(NA1)   = bs_pt_nc(NA)
bs_tt_nc(NA1)   = bs_tt_nc(NA)
polflow_t_nc(NA1) = polflow_t_nc(NA)

gamma_he3_nc(NA1) = gamma_he3_nc(NA)
dn_he3_nc(NA1)    = dn_he3_nc(NA)
cn_he3_nc(NA1)    = cn_he3_nc(NA)
qcond_he3_nc(NA1) = qcond_he3_nc(NA)
xhe3_nc(NA1)      = xhe3_nc(NA)
che3_nc(NA1)      = che3_nc(NA)
qen_he3_nc(NA1)   = qen_he3_nc(NA)
bs_phe3_nc(NA1)   = bs_phe3_nc(NA)
bs_the3_nc(NA1)   = bs_the3_nc(NA)
polflow_he3_nc(NA1) = polflow_he3_nc(NA)

gamma_he4_nc(NA1) = gamma_he4_nc(NA)
dn_he4_nc(NA1)    = dn_he4_nc(NA)
cn_he4_nc(NA1)    = cn_he4_nc(NA)
qcond_he4_nc(NA1) = qcond_he4_nc(NA)
xhe4_nc(NA1)      = xhe4_nc(NA)
che4_nc(NA1)      = che4_nc(NA)
qen_he4_nc(NA1)   = qen_he4_nc(NA)
bs_phe4_nc(NA1)   = bs_phe4_nc(NA)
bs_the4_nc(NA1)   = bs_the4_nc(NA)
polflow_he4_nc(NA1) = polflow_he4_nc(NA)

gamma_imp1_nc(NA1) = gamma_imp1_nc(NA)
dn_imp1_nc(NA1)    = dn_imp1_nc(NA)
cn_imp1_nc(NA1)    = cn_imp1_nc(NA)
qcond_imp1_nc(NA1) = qcond_imp1_nc(NA)
ximp1_nc(NA1)      = ximp1_nc(NA)
cimp1_nc(NA1)      = cimp1_nc(NA)
qen_imp1_nc(NA1)   = qen_imp1_nc(NA)
bs_pimp1_nc(NA1)   = bs_pimp1_nc(NA)
bs_timp1_nc(NA1)   = bs_timp1_nc(NA)
polflow_imp1_nc(NA1) = polflow_imp1_nc(NA)

gamma_imp2_nc(NA1) = gamma_imp2_nc(NA)
dn_imp2_nc(NA1)    = dn_imp2_nc(NA)
cn_imp2_nc(NA1)    = cn_imp2_nc(NA)
qcond_imp2_nc(NA1) = qcond_imp2_nc(NA)
ximp2_nc(NA1)      = ximp2_nc(NA)
cimp2_nc(NA1)      = cimp2_nc(NA)
qen_imp2_nc(NA1)   = qen_imp2_nc(NA)
bs_pimp2_nc(NA1)   = bs_pimp2_nc(NA)
bs_timp2_nc(NA1)   = bs_timp2_nc(NA)
polflow_imp2_nc(NA1) = polflow_imp2_nc(NA)

gamma_imp3_nc(NA1) = gamma_imp3_nc(NA)
dn_imp3_nc(NA1)    = dn_imp3_nc(NA)
cn_imp3_nc(NA1)    = cn_imp3_nc(NA)
qcond_imp3_nc(NA1) = qcond_imp3_nc(NA)
ximp3_nc(NA1)      = ximp3_nc(NA)
cimp3_nc(NA1)      = cimp3_nc(NA)
qen_imp3_nc(NA1)   = qen_imp3_nc(NA)
bs_pimp3_nc(NA1)   = bs_pimp3_nc(NA)
bs_timp3_nc(NA1)   = bs_timp3_nc(NA)
polflow_imp3_nc(NA1) = polflow_imp3_nc(NA)

jbs_nc(NA1)  = jbs_nc(NA)
jext_nc(NA1) = jext_nc(NA)
cc_nc(NA1)   = cc_nc(NA)
ni_nc(NA1)   = ni_nc(NA)

return
end subroutine neocl4

!----------------------------------------------------------------
subroutine ZBFAUX(GRRdB2, NGRTheta, YFM)
! Returns <(grad(rho)/B)2>
!         <n.grad(theta)>
!         and F_m - poloidal moments of geometric factor for PS viscosity
!                                  A.V.Zolotukhin 09-09-2003
!                  corrected by I.Yu. Senichenkov (IYS) February 2008
!-----------------------------------------------------------------------

use const_inc, only: ROC, HRO, NA, NA1, GP2, BTOR, RTOR
use status_inc, only: BDB0, BDB02, IPOL, MU, TRIA, ELON, SHIF, &
   AMETR, RHO, DRODA

implicit none

integer, parameter :: ntheta=32
integer :: j_theta, j_theta1, j_rho, j_mvisc
double precision :: YH, yametr, ydroda, yrho, I_str, THETA_cap,  &
   y_gamma, yr_theta, theta_s, yLambda, yLambdap, yDelta, yDeltap, &
   ytriang, ytriangp, COS_th, COS2_th, yd, y_sqg, y_sqg1, grr2dB2, &
   yBm1, Ynablarho2, yBmod, yBtheta, THETA_gen, theta1_s, &
   COS_th1, COS2_th1, yr_theta1, yd1, yBmod1, yBtheta1, y_sqg11, y_gamma1
double precision, dimension(1000) :: GRRdB2, NGRTheta, yVRR, yBm
double precision :: YFM(3,1000)
double precision, dimension(3) :: yBlnBCOS, yBlnBSIN, yB2COS, yB2SIN

do j_rho=1, NA
   YH = AMETR(j_rho+1) - AMETR(j_rho)
   yametr    = AMETR(j_rho)
   yrho      = RHO(j_rho)
   ydroda    = DRODA(j_rho)
   I_str     = IPOL(j_rho)*RTOR*BTOR
   THETA_cap = BTOR*yrho*ydroda*MU(j_rho)/I_str

! Magnetic surface ellipticity and its derivative
   yLambda  = ELON(j_rho)
   yLambdap = (ELON(j_rho+1)-ELON(j_rho))/YH

! Magnetic axis shift and its derivative
   yDelta   = SHIF(j_rho)
   yDeltap  = (SHIF(j_rho+1)-SHIF(j_rho))/YH

! Magnetic surface triangularity and its derivative
   ytriang  = TRIA(j_rho)                     !!**!!
   ytriangp = (TRIA(j_rho+1)-TRIA(j_rho))/YH  !!**!!


   GRRdB2(j_rho)   = 0.
   NGRTheta(j_rho) = 0.
   yVRR(j_rho)     = 0.
   yBm(j_rho)      = 0.
   do j_mvisc=1,3
      YFM(j_mvisc,j_rho) = 0.
      yBlnBCOS(j_mvisc) = 0.
      yBlnBSIN(j_mvisc) = 0.
      yB2COS(j_mvisc)   = 0.
      yB2SIN(j_mvisc)   = 0.
   enddo

   do j_theta=NTHETA,1,-1
      theta_s = GP2*j_theta/NTHETA
      COS_th  = COS(theta_s)
      COS2_th = COS_th*COS_th
      yr_theta = RTOR + yDelta + yametr*(COS_th - ytriang * (1.-COS2_th))
      yd = yLambda * (1.0 + yDeltap * Cos_th) + (1.-COS2_th) * &
          (yametr * yLambdap * (1.0 + 2.0 * ytriang * Cos_th) +  &
           yLambda * (ytriang - yametr * ytriangp) * Cos_th)

      y_sqg = yr_theta*yametr*yd ! Square root of the determinant g
      y_sqg1 = yr_theta*yd ! sqrt(g)/rho

! B - magnetic field

      yBmod = I_str/(yr_theta*yd)*sqrt(yd**2 + THETA_cap**2 * &
              ((1.0 + 2.0 * ytriang * Cos_th)**2 * (1.-COS2_th) + &
              yLambda**2 * Cos2_th))
! 1/B
      yBm1 = 1./yBmod
! B_theta - contravariant component of the magnetic field
      yBtheta = BTOR*MU(j_rho)*yrho*ydroda/y_sqg
! nabla(rho)2
      ynablarho2 = ydroda**2 * ((1.0 + 2.0 * ytriang * Cos_th)**2 *  &
                  (1.-COS2_th) + yLambda**2 * Cos2_th) / yd**2
! nabla(rho)2/B2
      grr2dB2 = ynablarho2/(yBmod*yBmod)
! Integral[sqrt(g)/rho, 0, 2pi]
      yVRR(j_rho) = yVRR(j_rho) + y_sqg1*GP2/NTHETA
! Integral[(nabla(rho)/B)2*SQRT(g)/rho, 0, 2pi]
      GRRdB2(j_rho)   = GRRdB2(j_rho)+ grr2dB2*y_sqg1*GP2/NTHETA
! Integral[1/B, 0, 2pi]
      THETA_gen = 0.
      do j_theta1=NTHETA,1,-1
         theta1_s   = theta_s*j_theta1/NTHETA
         COS_th1    = COS(theta1_s)
         COS2_th1   = COS_th1*COS_th1
         yr_theta1 = RTOR + yDelta + yametr*(COS_th1 - ytriang * (1.-COS2_th1))
         yd1 = yLambda * (1.0 + yDeltap * Cos_th1) + (1.-COS2_th1) * &
              (yametr * yLambdap * (1.0 + 2.0 * ytriang * Cos_th) +  &
               yLambda * (ytriang - yametr * ytriangp) * Cos_th1)
         y_sqg11    = yr_theta1*yd1

! B1 - magnetic field
         yBmod1 = I_str/(yr_theta1*yd1)*sqrt(yd1**2 + THETA_cap**2 * &
                 ((1.0 + 2.0 * ytriang * Cos_th1)**2 * (1.-COS2_th1) + &
                  yLambda**2 * Cos2_th1))
! B_theta1 - contravariant component of the magnetic field
         yBtheta1 = BTOR*MU(j_rho)*yrho*ydroda/y_sqg11/yametr
         THETA_gen = THETA_gen + yBmod1/yBtheta1*theta_s/NTHETA
         if (j_theta.eq.NTHETA) y_gamma1 = THETA_gen
      enddo  ! j_theta1

      if (j_theta.eq.NTHETA) y_gamma = GP2/y_gamma1
      THETA_gen = y_gamma*THETA_gen

      do j_mvisc=1,3
         yBlnBCOS(j_mvisc) = yBlnBCOS(j_mvisc) + &
                             yBmod*LOG(yBmod)*COS(j_mvisc*THETA_gen) * &
                             y_sqg1*GP2/NTHETA
         yBlnBSIN(j_mvisc) = yBlnBSIN(j_mvisc) + &
                             yBmod*LOG(yBmod)*SIN(j_mvisc*THETA_gen) * &
                             y_sqg1*GP2/NTHETA
         yB2COS(j_mvisc)   = yB2COS(j_mvisc) + &
                             yBmod*yBmod*COS(j_mvisc*THETA_gen) * &
                             y_sqg1*GP2/NTHETA
         yB2SIN(j_mvisc)   = yB2SIN(j_mvisc) + &
                             yBmod*yBmod*SIN(j_mvisc*THETA_gen) * &
                             y_sqg1*GP2/NTHETA
      enddo
      yBm(j_rho) = yBm(j_rho) + yBmod*y_sqg1*GP2/NTHETA
   enddo

   GRRdB2(j_rho)   = GRRdB2(j_rho)/yVRR(j_rho)
   NGRTheta(j_rho) = y_gamma
   yBm(j_rho) = yBm(j_rho)/yVRR(j_rho)

   do j_mvisc=1,3
      yBlnBCOS(j_mvisc) = yBlnBCOS(j_mvisc)/yVRR(j_rho)
      yBlnBSIN(j_mvisc) = yBlnBSIN(j_mvisc)/yVRR(j_rho)
      yB2COS(j_mvisc)   = yB2COS(j_mvisc)/yVRR(j_rho)
      yB2SIN(j_mvisc)   = yB2SIN(j_mvisc)/yVRR(j_rho)
      YFM(j_mvisc,j_rho) = 2.*(j_mvisc*y_gamma)**2/ &
                           (BTOR**3*BDB02(j_rho)*BDB0(j_rho))* &
                           (  yBlnBCOS(j_mvisc)*yB2COS(j_mvisc)  + &
                              yBlnBSIN(j_mvisc)*yB2SIN(j_mvisc)   )
   enddo

enddo

GRRdB2(NA1) = GRRdB2(NA)*(ROC/HRO-(NA-1)) + GRRdB2(NA-1)*(NA-ROC/HRO)
NGRTheta(NA1) = NGRTheta(NA)*(ROC/HRO-(NA-1)) + NGRTheta(NA-1)*(NA-ROC/HRO)
do j_mvisc=1,3
   YFM(j_mvisc,NA1) = YFM(j_mvisc,NA)*(ROC/HRO-(NA-1))+ &
                      YFM(j_mvisc,NA-1)*(NA-ROC/HRO)
enddo

return
end subroutine zbfaux
