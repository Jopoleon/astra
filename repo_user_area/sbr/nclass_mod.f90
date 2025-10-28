module nclass_mod

use parameter_inc, only: NRD

implicit none

double precision, dimension(NRD) :: &
    gamma_e_nc, dn_e_nc, cn_e_nc, qcond_e_nc, xe_nc, ce_nc, qen_e_nc, &
    bs_pe_nc, bs_te_nc, polflow_e_nc, &
    gamma_i_nc, dn_i_nc, cn_i_nc, qcond_i_nc, xi_nc, ci_nc, qen_i_nc, &
    bs_pi_nc, bs_ti_nc, polflow_i_nc, &
    gamma_p_nc, dn_p_nc, cn_p_nc, qcond_p_nc, xp_nc, cp_nc, qen_p_nc, &
    bs_pp_nc, bs_tp_nc, polflow_p_nc, &
    gamma_d_nc, dn_d_nc, cn_d_nc, qcond_d_nc, xd_nc, cd_nc, qen_d_nc, &
    bs_pd_nc, bs_td_nc, polflow_d_nc, &
    gamma_t_nc, dn_t_nc, cn_t_nc, qcond_t_nc, xt_nc, ct_nc, qen_t_nc, &
    bs_pt_nc, bs_tt_nc, polflow_t_nc, &
    gamma_he3_nc, dn_he3_nc, cn_he3_nc, qcond_he3_nc, xhe3_nc, che3_nc, qen_he3_nc, &
    bs_phe3_nc, bs_the3_nc, polflow_he3_nc, &
    gamma_he4_nc, dn_he4_nc, cn_he4_nc, qcond_he4_nc, xhe4_nc, che4_nc, qen_he4_nc, &
    bs_phe4_nc, bs_the4_nc, polflow_he4_nc, &
    gamma_imp1_nc, dn_imp1_nc, cn_imp1_nc, qcond_imp1_nc, ximp1_nc, cimp1_nc, qen_imp1_nc, &
    bs_pimp1_nc, bs_timp1_nc, polflow_imp1_nc, &
    gamma_imp2_nc, dn_imp2_nc, cn_imp2_nc, qcond_imp2_nc, ximp2_nc, cimp2_nc, qen_imp2_nc, &
    bs_pimp2_nc, bs_timp2_nc, polflow_imp2_nc, &
    gamma_imp3_nc, dn_imp3_nc, cn_imp3_nc, qcond_imp3_nc, ximp3_nc, cimp3_nc, qen_imp3_nc, &
    bs_pimp3_nc, bs_timp3_nc, polflow_imp3_nc, &
    jbs_nc, jext_nc, cc_nc, ni_nc

integer, parameter :: mx_mi=9, mx_ms=40, mx_mz=100

contains

!----------------------------------------------------------------
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
    ULON, ER, VRS, G11, RHO, AMETR, &
    MU, ELON, SHIF, TE, TI, &
    NE, NHYDR, NDEUT, NTRIT, NHE3, NALF, ZIM1, ZIM2, ZIM3, NIZ1, NIZ2, NIZ3, NMAIN

!  Physical and conversion constants
real, parameter :: z_coulomb=1.6022e-19, z_electronmass=9.1095e-31, &
    z_protonmass=1.6726e-27, z_j7kv=1.6022e-16, z_mu0=1.2566e-06

double precision :: YGRRdB2(1000), YNGRTHETA(1000), YFM(3, 1000)

double precision :: y_grrho2, grti, y_den, yh
real, dimension(mx_ms) :: dq_s, vq_s
real :: eps_shear_fac
integer :: k_electron, k_mainion, k_proton, k_deuteron, k_triton, &
    k_he3, k_alpha, k_impZ1, k_impZ2, k_impZ3, narray
real :: ybbmax, ybbmax2, yftupper, yftlower
real rdum(8)
!Declaration of input to NCLASS
integer :: k_order, k_potato, m_i, m_z
real :: c_den, c_potb, c_potl
real :: p_b2, p_bm2, p_eb, p_fhat, p_fm(3), p_ft, p_grbm2, p_grphi, p_gr2phi, p_ngrth
real, dimension(mx_mi) :: amu_i, grt_i, temp_i
real :: den_iz(mx_mi, mx_mz), fex_iz(3, mx_mi, mx_mz), grp_iz(mx_mi, mx_mz)
real :: YNMAIN, YNMAIN1
!Declaration of output from NCLASS
integer :: iflag, m_s
integer, dimension(mx_ms) :: jm_s, jz_s
real :: p_bsjb, p_etap, p_exjb
real :: calm_i(3, 3, mx_mi)
real, dimension(3, 3, mx_mi, mx_mi) :: caln_ii, capm_ii, capn_ii
real, dimension(mx_ms) :: bsjbp_s, bsjbt_s, dn_s, sqz_s, vn_s, veb_s, qeb_s, xi_s
real, dimension(5, mx_ms) :: gfl_s, qfl_s  
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
!  grp_iz(i, z)-pressure gradient of i, z (keV/m**3/rho)

call ZBFAUX(yGRRdB2, yNGRTheta, YFM)

do j=1, NA
! Set radially dependent data
!  eps_shear_fac-geometrical factor
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

    c_potl   = -(RTOR + SHIF(j))/MU(1)
    y_grrho2 = G11(j)/VRS(j)

    if (j == 1) then
        eps_shear_fac = BTOR*RHO(j)*MU(j)/RTOR/ &
                ((SHIF(j+1) + AMETR(j+1) - SHIF(j) - AMETR(j))/ &
                (RHO(j) - RHO(j+1)))
    else
        eps_shear_fac = BTOR*RHO(j)*MU(j)/RTOR/ &
                ((SHIF(j+1) + AMETR(j+1) - SHIF(j-1) - AMETR(j-1))/ &
                (RHO(j-1) - RHO(j+1)))
    endif
    p_grphi  = ER(j)
    p_gr2phi = MU(j)*j*(ER(j+1)/MU(j+1)/((j+1)*YH) - ER(j)/MU(j)/(j*YH))
! Warning: The NCLASS version 1.2 returns NaN resistivity
!          if p_eb = 0.
!         p_eb     = max(.00001, -ULON(j)*BTOR/(GP2*RTOR))
    p_eb     = -ULON(j)*BTOR/(GP2*RTOR)
    if (abs(p_eb) < 1.d-3) p_eb = sign(0.001, p_eb)
    p_b2     = BTOR**2*BDB02(j)
    p_bm2    = B0DB2(j)/BTOR**2
    p_fhat   = -IPOL(j)*RTOR/(MU(j)*j*YH)
! Trapped particle fraction
    ybbmax   = BTOR/BMAXT(j)*BDB0(j)      !->  <h>=<B/B_max>
    ybbmax2  = (BTOR/BMAXT(j))**2*BDB02(j)!->  <h2>=<(B/B_max)2>
    yftupper = 1.- ybbmax2/ybbmax**2* (1.- sqrt(1. - ybbmax)*(1. + 0.5*ybbmax))
    yftlower = 1. - ybbmax2*(BMAXT(j)/BTOR)**2*FOFB(j)
    p_ft     = 0.75*yftupper + 0.25*yftlower
    p_grbm2  = yGRRdB2(j)   ! <grad(rho)**2/B**2>
    p_ngrth  = yNGRTheta(j) ! <n.grad(Theta)> 

!  p_fm(3)-poloidal moments of geometric factor for PS viscosity [-]
    do i=1, 3
        p_fm(i) = 0.0
    enddo
    do i=1, 3
        p_fm(i) = YFM(i, j)
    enddo

    do j1=1, mx_mi
        do jj=1, mx_mz
            den_iz(j1, jj) = 0.
            grp_iz(j1, jj) = 0.
            do i=1, 3
                fex_iz(i, j1, jj) = 0.
            enddo
        enddo
    enddo

    m_i = 1                ! Reserved for electrons
    k_electron = 1
    m_z = 1
    amu_i(m_i)  = 5.4463e-4
    temp_i(m_i) = TE(j)
    grt_i(m_i)  = (TE(j+1)-TE(j))/YH
    den_iz(m_i, m_z) = 1.e19*NE(j)
    grp_iz(m_i, m_z) = 1.e19/YH*(TE(j+1)*NE(j+1) - TE(j)*NE(j)) 
    grti = (TI(j+1) - TI(j))/YH

    ni_nc(j) = 0. !density of all thermal ions

    if (NHYDR(j) > y_den) then
        ni_nc(j) = ni_nc(j) + nhydr(j)
        m_i = m_i + 1
        k_proton = m_i
        m_z = 1
        amu_i(m_i) = 1.
        grt_i(m_i) = grti
        temp_i(m_i) = TI(j)
        den_iz(m_i, m_z) = 1.e19*NHYDR(j)
        grp_iz(m_i, m_z) = 1.e19/YH*(TI(j+1)*NHYDR(j+1) - TI(j)*NHYDR(j))
    endif

    if (NDEUT(j) > y_den) then
        ni_nc(j) = ni_nc(j) + ndeut(j)
        m_i = m_i + 1
        k_deuteron = m_i
        amu_i(m_i) = 2.
        grt_i(m_i) = grti
        temp_i(m_i) = TI(j)
        den_iz(m_i, m_z) = 1.e19*NDEUT(j)
        grp_iz(m_i, m_z) = 1.e19/YH*(TI(j+1)*NDEUT(j+1) - TI(j)*NDEUT(j)) 
    endif

    if (NTRIT(j) > y_den) then
        ni_nc(j) = ni_nc(j) + ntrit(j)
        m_i = m_i + 1
        k_triton = m_i
        amu_i(m_i) = 3.
        grt_i(m_i) = grti
        temp_i(m_i) = TI(j)
        den_iz(m_i, m_z) = 1.e19*NTRIT(j)
        grp_iz(m_i, m_z) = 1.e19/YH*(TI(j+1)*NTRIT(j+1) - TI(j)*NTRIT(j))
    endif

    if (m_i == 1 .and. mx_mz < 2 .and. ZMJ >= 2 ) then
        write(*, *)">>> NEOCL: ion composition error:"
        write(*, *)" mx_mz < 2 and no hydrogen. Call ignored."
        return
    endif

! It is assumed here that the ion species are not specified explicitly
! in the model
! or/and
! is used to satisfy the quasineutrality condition
    YNMAIN = 1./ZMJ*(NE(J) - ZIM1(J)*NIZ1(J) - ZIM2(J)*NIZ2(J) - ZIM3(J)*NIZ3(J) - &
             NHYDR(J) - NDEUT(J) - NTRIT(J) - 2.*NHE3(J) - 2.*NALF(J))
    YNMAIN1 = 1./ZMJ*(NE(J+1) - ZIM1(J+1)*NIZ1(J+1) - ZIM2(J+1)*NIZ2(J+1) - ZIM3(J+1)*NIZ3(J+1) - &
              NHYDR(J+1) - NDEUT(J+1) - NTRIT(J+1) - 2.*NHE3(J+1) - 2.*NALF(J+1))
    if (YNMAIN > y_den) then
        ni_nc(j) = ni_nc(j) + ynmain
        m_i = m_i + 1
        k_mainion = m_i
        m_z = ZMJ
        amu_i(m_i) = AMJ
        temp_i(m_i) = TI(j)
        grt_i(m_i) = grti
        den_iz(m_i, m_z) = 1.e19*YNMAIN
        grp_iz(m_i, m_z) = 1.e19/YH*(TI(j+1)*YNMAIN1 - TI(j)*YNMAIN) 
    endif

    if (NHE3(j) > y_den) then
        ni_nc(j) = ni_nc(j) + nhe3(j)
        m_i = m_i + 1
        k_he3 = m_i
        amu_i(m_i) = 3.
        grt_i(m_i) = grti
        temp_i(m_i) = TI(j)
        m_z = 2
        den_iz(m_i, m_z) = 1.e19*NHE3(j)
        grp_iz(m_i, m_z) = 1.e19/YH*(TI(j+1)*NHE3(j+1) - TI(j)*NHE3(j)) 
    endif

    if (m_i == mx_mi) goto 5

    if (NALF(j) > y_den) then
        m_i = m_i + 1
        k_alpha = m_i
        amu_i(m_i) = 4.
        grt_i(m_i) = grti
        temp_i(m_i) = TI(j)
        m_z = 2
        den_iz(m_i, m_z) = 1.e19*NALF(j)
        grp_iz(m_i, m_z) = 1.e19/YH*(TI(j+1)*NALF(j+1) - TI(j)*NALF(j)) 
    endif

! consider up to three impurity species
    jj = nint(ZIM1(j))
    if (jj <= mx_mz) then
        if (NIZ1(j) > y_den) then
            ni_nc(j) = ni_nc(j) + niz1(j)
            m_i = m_i + 1
            k_impZ1 = m_i
            amu_i(m_i) = AIM1
            grt_i(m_i) = grti
            temp_i(m_i) = TI(j)
            if (jj > m_z) m_z = jj
            den_iz(m_i, jj) = 1.e19*NIZ1(j)
            grp_iz(m_i, jj) = 1.e19/YH*(TI(j+1)*NIZ1(j+1) - TI(j)*NIZ1(j)) 
        endif
        if (m_i == mx_mi) goto 5
    endif

    jj = nint(ZIM2(j))
    if (jj <= mx_mz) then
        if (NIZ2(j) > y_den) then
            ni_nc(j) = ni_nc(j) + niz2(j)
            m_i = m_i + 1
            k_impZ2 = m_i
            amu_i(m_i) = AIM2
            grt_i(m_i) = grti
            temp_i(m_i) = TI(j)
            if (jj > m_z) m_z = jj
            den_iz(m_i, jj) = 1.e19*NIZ2(j)
            grp_iz(m_i, jj) = 1.e19/YH*(TI(j+1)*NIZ2(j+1) - TI(j)*NIZ2(j)) 
        endif
        if (m_i == mx_mi) goto 5
    endif

    jj = nint(ZIM3(j))
    if (jj <= mx_mz) then
        if (NIZ3(j) > y_den) then
            ni_nc(j) = ni_nc(j) + niz3(j)
            m_i = m_i + 1
            k_impZ3 = m_i
            amu_i(m_i) = AIM3
            grt_i(m_i) = grti
            temp_i(m_i) = TI(j)
            if (jj > m_z) m_z = jj
            den_iz(m_i, jj) = 1.e19*NIZ3(j)
            grp_iz(m_i, jj) = 1.e19/YH*(TI(j+1)*NIZ3(j+1) - TI(j)*NIZ3(j)) 
        endif
    endif
  
 5 continue

    call NCLASS(k_order, k_potato, m_i, m_z, c_den, c_potb, c_potl, p_b2, &
                p_bm2, p_eb, p_fhat, p_fm, p_ft, p_grbm2, p_grphi, p_gr2phi, &
                p_ngrth, amu_i, grt_i, temp_i, den_iz, fex_iz, grp_iz, &
! output:
                m_s, jm_s, jz_s, p_bsjb, p_etap, p_exjb, calm_i, caln_ii, capm_ii, &
                capn_ii, bsjbp_s, bsjbt_s, dn_s, gfl_s, qfl_s, sqz_s, upar_s, &
                utheta_s, vn_s, veb_s, qeb_s, xi_s, ymu_s, chip_ss, chit_ss, &
                dp_ss, dt_ss, iflag)

    if (iflag > 0) then
        nout = 6
        open(unit=nout, file='user')
        SELECT CASE(iflag)
        CASE(1)
            label = 'ERROR:NCLASS-k_order must be 2 or 3, k_order='
            idum(1) = k_order
            call WRITE_LINE_IR(nout, label, 1, idum, 0, rdum, 0)
        CASE(2)
            label='ERROR:NCLASS-require 1<m_i<mx_mi, m_i='
            idum(1) = m_i
            call WRITE_LINE_IR(nout, label, 1, idum, 0, rdum, 0)
        CASE(3)
            label = 'ERROR:NCLASS-require 0<m_z<mx_mz, m_z='
            idum(1) = m_z
            call WRITE_LINE_IR(nout, label, 1, idum, 0, rdum, 0)
        CASE(4)
            label = 'ERROR:NCLASS-require 0<m_s<mx_ms, m_s='
            idum(1) = m_s
            call WRITE_LINE_IR(nout, label, 1, idum, 0, rdum, 0)
        CASE(5)
            label = 'ERROR:NCLASS-inversion of flow matrix failed'
            call WRITE_LINE(nout, label, 0, 0)
        CASE(6)
            label = 'ERROR:NCLASS-Trapped fraction not between 0 and 1'
            call WRITE_LINE(nout, label, 0, 0)
        END SELECT
        return
    endif

! Output
! Electrons
    if (k_electron > 0) then
! Total radial particle flux (electrons)
        call RARRAY_COPY(5, gfl_s(1, k_electron), 1, rdum, 1)
        rdum(6) = RARRAY_SUM(5, rdum, 1)
        gamma_e_nc(j) = rdum(6)*VRS(j)*1.e-19
! Particle diffusion and velocity (electrons)
        im  = jm_s(k_electron)
        iza = IABS(jz_s(k_electron))
        dn_e_nc(j) = dn_s(k_electron)/y_grrho2
        cn_e_nc(j) = (vn_s(k_electron) + veb_s(k_electron) + &
                      gfl_s(5, k_electron)/den_iz(im, iza)) / y_grrho2
! Radial conduction flux (electrons)
        call RARRAY_COPY(5, qfl_s(1, k_electron), 1, rdum, 1)
        rdum(6) = RARRAY_SUM(5, rdum, 1)
        qcond_e_nc(j) = rdum(6)*VRS(j)*1.e-6
! Heat conduction and velocity (electrons)
! Conduction total is sum of components
        rdum(1) = RARRAY_SUM(5, qfl_s(1, k_electron), 1)
! Diagonal conductivity plus convective velocity
        dq_s(k_electron) = chit_ss(k_electron, k_electron) + &
                           chip_ss(k_electron, k_electron)
        vq_s(k_electron) = rdum(1)/(den_iz(im, iza)*z_j7kv*temp_i(im)) + &
                           dq_s(k_electron)*grt_i(im)/temp_i(im)
        xe_nc(j) = dq_s(k_electron)/y_grrho2
        ce_nc(j) = vq_s(k_electron)/y_grrho2
! Total radial energy flux (electrons)
        do k=1, 5
            rdum(k) = qfl_s(k, k_electron) + 2.5*gfl_s(k, k_electron)*temp_i(im)*z_j7kv
        enddo
        rdum(6) = RARRAY_SUM(5, rdum, 1)
        qen_e_nc(j) = rdum(6)*VRS(j)*1.e-6
! Bootstrap current on p'/p (electrons)
        rdum(1) = bsjbp_s(k_electron)
        rdum(2) = bsjbt_s(k_electron)
        bs_pe_nc(j) = -1.e-6*rdum(1)/BTOR
        bs_te_nc(j) = -1.e-6*rdum(2)/BTOR
! Poloidal flow velociry on outside midplane (electrons)
        polflow_e_nc(j) = eps_shear_fac * ( &
            utheta_s(1, 1, k_electron) + &
            utheta_s(1, 2, k_electron) + &
            utheta_s(1, 3, k_electron) ) 
    endif

!---  Main Ions  ---
    if (k_mainion > 0) then
! Total radial particle flux (main ions)
        call RARRAY_COPY(5, gfl_s(1, k_mainion), 1, rdum, 1)
        rdum(6) = RARRAY_SUM(5, rdum, 1)
        gamma_i_nc(j) = rdum(6)*VRS(j)*1.e-19
! Particle diffusion and velocity (main ions)
        im  = jm_s(k_mainion)
        iza = IABS(jz_s(k_mainion))
        dn_i_nc(j) = dn_s(k_mainion)/y_grrho2
        cn_i_nc(j) = (vn_s(k_mainion) + veb_s(k_mainion) + &
                      gfl_s(5, k_mainion)/den_iz(im, iza)) / y_grrho2
! Radial conduction flux (main ions)
        call RARRAY_COPY(5, qfl_s(1, k_mainion), 1, rdum, 1)
        rdum(6) = RARRAY_SUM(5, rdum, 1)
        qcond_i_nc(j) = rdum(6)*VRS(j)*1.e-6
! Heat conduction and velocity (main ions)
! Conduction total is sum of components
        rdum(1) = RARRAY_SUM(5, qfl_s(1, k_mainion), 1)
! Diagonal conductivity plus convective velocity
        dq_s(k_mainion) = chit_ss(k_mainion, k_mainion) + &
                          chip_ss(k_mainion, k_mainion)
        vq_s(k_mainion) = rdum(1)/(den_iz(im, iza)* z_j7kv*temp_i(im)) + &
                          dq_s(k_mainion)*grt_i(im)/temp_i(im)
        xi_nc(j) = dq_s(k_mainion)/y_grrho2
        ci_nc(j) = vq_s(k_mainion)/y_grrho2
! Total radial energy flux (main ions)
        do k=1, 5
            rdum(k) = qfl_s(k, k_mainion) + 2.5*gfl_s(k, k_mainion)*temp_i(im)*z_j7kv
        enddo
        rdum(6) = RARRAY_SUM(5, rdum, 1)
        qen_i_nc(j) = rdum(6)*VRS(j)*1.e-6
! Bootstrap current on p'/p (main ions)
        rdum(1) = bsjbp_s(k_mainion)
        rdum(2) = bsjbt_s(k_mainion)
        bs_pi_nc(j) = -1.e-6*rdum(1)/BTOR
        bs_ti_nc(j) = -1.e-6*rdum(2)/BTOR
! Poloidal flow velociry on outside midplane (main ions)
        polflow_i_nc(j) = eps_shear_fac * ( &
            utheta_s(1, 1, k_mainion) + &
            utheta_s(1, 2, k_mainion) + &
            utheta_s(1, 3, k_mainion) )
    endif

!--- Protons ---
    if (k_proton > 0) then
! Total radial particle flux (protons  )
        call RARRAY_COPY(5, gfl_s(1, k_proton), 1, rdum, 1)
        rdum(6) = RARRAY_SUM(5, rdum, 1)
        gamma_p_nc(j) = rdum(6)*VRS(j)*1.e-19
! Particle diffusion and velocity (protons  )
        im  = jm_s(k_proton)
        iza = IABS(jz_s(k_proton))
        dn_p_nc(j) = dn_s(k_proton)/y_grrho2
        cn_p_nc(j) = (vn_s(k_proton) + veb_s(k_proton) + &
                       gfl_s(5, k_proton)/den_iz(im, iza))  / y_grrho2
! Radial conduction flux (protons  )
        call RARRAY_COPY(5, qfl_s(1, k_proton), 1, rdum, 1)
        rdum(6) = RARRAY_SUM(5, rdum, 1)
        qcond_p_nc(j) = rdum(6)*VRS(j)*1.e-6
! Heat conduction and velocity (protons  )
! Conduction total is sum of components
        rdum(1)=RARRAY_SUM(5, qfl_s(1, k_proton), 1)
! Diagonal conductivity plus convective velocity
        dq_s(k_proton)=chit_ss(k_proton, k_proton) + chip_ss(k_proton, k_proton)
        vq_s(k_proton)=rdum(1)/(den_iz(im, iza)* z_j7kv*temp_i(im)) + &
                        dq_s(k_proton)*grt_i(im)/temp_i(im)
        xp_nc(j) = dq_s(k_proton)/y_grrho2
        cp_nc(j) = vq_s(k_proton)/y_grrho2
! Total radial energy flux (protons  )
        do k=1, 5
            rdum(k) = qfl_s(k, k_proton) + 2.5*gfl_s(k, k_proton)*temp_i(im)*z_j7kv
        enddo
        rdum(6)=RARRAY_SUM(5, rdum, 1)
        qen_p_nc(j) = rdum(6)*VRS(j)*1.e-6
! Bootstrap current on p'/p (protons  )
        rdum(1)=bsjbp_s(k_proton)
        rdum(2)=bsjbt_s(k_proton)
        bs_pp_nc(j) = -1.e-6*rdum(1)/BTOR
        bs_tp_nc(j) = -1.e-6*rdum(2)/BTOR
! Poloidal flow velociry on outside midplane (protons)
        polflow_p_nc(j) = eps_shear_fac * ( &
            utheta_s(1, 1, k_proton) +  &
            utheta_s(1, 2, k_proton) + &
            utheta_s(1, 3, k_proton) )
    endif

!--- Deuterons ---
    if (k_deuteron > 0) then
! Total radial particle flux (deuterons  )
        call RARRAY_COPY(5, gfl_s(1, k_deuteron), 1, rdum, 1)
        rdum(6) = RARRAY_SUM(5, rdum, 1)
        gamma_d_nc(j) = rdum(6)*VRS(j)*1.e-19
! Particle diffusion and velocity (deuterons  )
        im  = jm_s(k_deuteron)
        iza = IABS(jz_s(k_deuteron))
        dn_d_nc(j) = dn_s(k_deuteron)/y_grrho2
        cn_d_nc(j) = (vn_s(k_deuteron) + veb_s(k_deuteron) + &
                       gfl_s(5, k_deuteron)/den_iz(im, iza)) / y_grrho2
! Radial conduction flux (deuterons  )
        call RARRAY_COPY(5, qfl_s(1, k_deuteron), 1, rdum, 1)
        rdum(6) = RARRAY_SUM(5, rdum, 1)
        qcond_d_nc(j) = rdum(6)*VRS(j)*1.e-6
! Heat conduction and velocity (deuterons  )
! Conduction total is sum of components
        rdum(1)=RARRAY_SUM(5, qfl_s(1, k_deuteron), 1)
! Diagonal conductivity plus convective velocity
        dq_s(k_deuteron) = chit_ss(k_deuteron, k_deuteron) + &
                           chip_ss(k_deuteron, k_deuteron)
        vq_s(k_deuteron) = rdum(1)/(den_iz(im, iza)* z_j7kv*temp_i(im)) + &
                           dq_s(k_deuteron)*grt_i(im)/temp_i(im)
        xd_nc(j) = dq_s(k_deuteron)/y_grrho2
        cd_nc(j) = vq_s(k_deuteron)/y_grrho2
! Total radial energy flux (deuterons  )
        do k=1, 5
            rdum(k) = qfl_s(k, k_deuteron) + 2.5*gfl_s(k, k_deuteron)*temp_i(im)*z_j7kv
        enddo
        rdum(6) = RARRAY_SUM(5, rdum, 1)
        qen_d_nc(j) = rdum(6)*VRS(j)*1.e-6
! Bootstrap current on p'/p (deuterons  )
        rdum(1) = bsjbp_s(k_deuteron)
        rdum(2) = bsjbt_s(k_deuteron)
        bs_pd_nc(j) = -1.e-6*rdum(1)/BTOR
        bs_td_nc(j) = -1.e-6*rdum(2)/BTOR
! Poloidal flow velociry on outside midplane (deuterons)
        polflow_d_nc(j) = eps_shear_fac * ( &
            utheta_s(1, 1, k_deuteron) + &
            utheta_s(1, 2, k_deuteron) + &
            utheta_s(1, 3, k_deuteron) )
    endif
     
!--- Tritons ---
    if (k_triton > 0) then
! Total radial particle flux (tritons  )
        call RARRAY_COPY(5, gfl_s(1, k_triton), 1, rdum, 1)
        rdum(6) = RARRAY_SUM(5, rdum, 1)
        gamma_t_nc(j) = rdum(6)*VRS(j)*1.e-19
! Particle diffusion and velocity (tritons  )
        im  = jm_s(k_triton)
        iza = IABS(jz_s(k_triton))
        dn_t_nc(j) = dn_s(k_triton)/y_grrho2
        cn_t_nc(j) = (vn_s(k_triton) + veb_s(k_triton) + &
                      gfl_s(5, k_triton)/den_iz(im, iza)) / y_grrho2
! Radial conduction flux (tritons  )
        call RARRAY_COPY(5, qfl_s(1, k_triton), 1, rdum, 1)
        rdum(6) = RARRAY_SUM(5, rdum, 1)
        qcond_t_nc(j) = rdum(6)*VRS(j)*1.e-6
! Heat conduction and velocity (tritons  )
! Conduction total is sum of components
        rdum(1) = RARRAY_SUM(5, qfl_s(1, k_triton), 1)
! Diagonal conductivity plus convective velocity
        dq_s(k_triton) = chit_ss(k_triton, k_triton) + &
                         chip_ss(k_triton, k_triton)
        vq_s(k_triton) = rdum(1)/(den_iz(im, iza)* z_j7kv*temp_i(im)) + &
                         dq_s(k_triton)*grt_i(im)/temp_i(im)
        xt_nc(j) = dq_s(k_triton)/y_grrho2
        ct_nc(j) = vq_s(k_triton)/y_grrho2
! Total radial energy flux (tritons  )
        do k=1, 5
            rdum(k) = qfl_s(k, k_triton) + 2.5*gfl_s(k, k_triton)*temp_i(im)*z_j7kv
        enddo
        rdum(6) = RARRAY_SUM(5, rdum, 1)
        qen_t_nc(j) = rdum(6)*VRS(j)*1.e-6
! Bootstrap current on p'/p (tritons  )
        rdum(1) = bsjbp_s(k_triton)
        rdum(2) = bsjbt_s(k_triton)
        bs_pt_nc(j) = -1.e-6*rdum(1)/BTOR
        bs_tt_nc(j) = -1.e-6*rdum(2)/BTOR
! Poloidal flow velociry on outside midplane (tritons)
        polflow_t_nc(j) = eps_shear_fac * ( &
            utheta_s(1, 1, k_triton) + &
            utheta_s(1, 2, k_triton) + &
            utheta_s(1, 3, k_triton) )
    endif

!--- He3 Particles ---
    if (k_he3 > 0) then
! Total radial particle flux (He3 particles)
        call RARRAY_COPY(5, gfl_s(1, k_he3), 1, rdum, 1)
        rdum(6) = RARRAY_SUM(5, rdum, 1)
        gamma_he3_nc = rdum(6)*VRS(j)*1.e-19
! Particle diffusion and velocity (He3 particles)
        im  = jm_s(k_he3)
        iza = IABS(jz_s(k_he3))
        dn_he3_nc(j) = dn_s(k_he3)/y_grrho2
        cn_he3_nc(j) = (vn_s(k_he3)  + veb_s(k_he3) + &
                       gfl_s(5, k_he3)/den_iz(im, iza))  / y_grrho2
! Radial conduction flux (He3 particles)
        call RARRAY_COPY(5, qfl_s(1, k_he3), 1, rdum, 1)
        rdum(6) = RARRAY_SUM(5, rdum, 1)
        qcond_he3_nc(j) = rdum(6)*VRS(j)*1.e-6
! Heat conduction and velocity (He3 particles)
! Conduction total is sum of components
        rdum(1) = RARRAY_SUM(5, qfl_s(1, k_he3), 1)
! Diagonal conductivity plus convective velocity
        dq_s(k_he3) = chit_ss(k_he3, k_he3) + chip_ss(k_he3, k_he3)
        vq_s(k_he3) = rdum(1)/(den_iz(im, iza)* z_j7kv*temp_i(im)) + &
                      dq_s(k_he3)*grt_i(im)/temp_i(im)
        xhe3_nc(j) = dq_s(k_he3)/y_grrho2
        che3_nc(j) = vq_s(k_he3)/y_grrho2
! Total radial energy flux (He3 particles)
        do k=1, 5
            rdum(k) = qfl_s(k, k_he3) + 2.5*gfl_s(k, k_he3)*temp_i(im)*z_j7kv
        enddo
        rdum(6) = RARRAY_SUM(5, rdum, 1)
        qen_he3_nc(j) = rdum(6)*VRS(j)*1.e-6
! Bootstrap current on p'/p (He3 particles)
        rdum(1) = bsjbp_s(k_he3)
        rdum(2) = bsjbt_s(k_he3)
        bs_phe3_nc(j) = -1.e-6*rdum(1)/BTOR
        bs_the3_nc(j) = -1.e-6*rdum(2)/BTOR
! Poloidal flow velociry on outside midplane (He3 particles)
        polflow_he3_nc(j) = eps_shear_fac * ( &
            utheta_s(1, 1, k_he3) + &
            utheta_s(1, 2, k_he3) + &
            utheta_s(1, 3, k_he3) )
    endif

!--- Alpha Particles ---
    if (k_alpha > 0) then
! Total radial particle flux (alpha particles)
        call RARRAY_COPY(5, gfl_s(1, k_alpha), 1, rdum, 1)
        rdum(6) = RARRAY_SUM(5, rdum, 1)
        gamma_he4_nc(j) = rdum(6)*VRS(j)*1.e-19
! Particle diffusion and velocity (alpha particles)
        im  = jm_s(k_alpha)
        iza = IABS(jz_s(k_alpha))
        dn_he4_nc(j) = dn_s(k_alpha)/y_grrho2
        cn_he4_nc(j) = (vn_s(k_alpha) + veb_s(k_alpha) + &
                        gfl_s(5, k_alpha)/den_iz(im, iza)) / y_grrho2
! Radial conduction flux (alpha particles)
        call RARRAY_COPY(5, qfl_s(1, k_alpha), 1, rdum, 1)
        rdum(6) = RARRAY_SUM(5, rdum, 1)
        qcond_he4_nc(j) = rdum(6)*VRS(j)*1.e-6
! Heat conduction and velocity (alpha particles)
! Conduction total is sum of components
        rdum(1) = RARRAY_SUM(5, qfl_s(1, k_alpha), 1)
! Diagonal conductivity plus convective velocity
        dq_s(k_alpha) = chit_ss(k_alpha, k_alpha) + chip_ss(k_alpha, k_alpha)
        vq_s(k_alpha) = rdum(1)/(den_iz(im, iza)* z_j7kv*temp_i(im)) + &
                        dq_s(k_alpha)*grt_i(im)/temp_i(im)
        xhe4_nc(j) = dq_s(k_alpha)/y_grrho2
        che4_nc(j) = vq_s(k_alpha)/y_grrho2
! Total radial energy flux (alpha particles)
        do k=1, 5
            rdum(k) = qfl_s(k, k_alpha) + 2.5*gfl_s(k, k_alpha)*temp_i(im)*z_j7kv
        enddo
        rdum(6) = RARRAY_SUM(5, rdum, 1)
        qen_he4_nc(j) = rdum(6)*VRS(j)*1.e-6
! Bootstrap current on p'/p (alpha particles)
        rdum(1) = bsjbp_s(k_alpha)
        rdum(2) = bsjbt_s(k_alpha)
        bs_phe4_nc(j) = -1.e-6*rdum(1)/BTOR
        bs_the4_nc(j) = -1.e-6*rdum(2)/BTOR
! Poloidal flow velociry on outside midplane (alpha particles)
        polflow_he4_nc(j) = eps_shear_fac * ( &
            utheta_s(1, 1, k_alpha) + &
            utheta_s(1, 2, k_alpha) + &
            utheta_s(1, 3, k_alpha) )
    endif

!--- Impurity Z1 ---
    if (k_impZ1 > 0) then
! Total radial particle flux (impZ1 particles)
        call RARRAY_COPY(5, gfl_s(1, k_impZ1), 1, rdum, 1)
        rdum(6) = RARRAY_SUM(5, rdum, 1)
        gamma_imp1_nc(j) = rdum(6)*VRS(j)*1.e-19
! Particle diffusion and velocity (impZ1 particles)
        im  = jm_s(k_impZ1)
        iza = IABS(jz_s(k_impZ1))
        dn_imp1_nc(j) = dn_s(k_impZ1)/y_grrho2
        cn_imp1_nc(j) = (vn_s(k_impZ1) + veb_s(k_impZ1) + &
                        gfl_s(5, k_impZ1)/den_iz(im, iza)) / y_grrho2
! Radial conduction flux (impZ1 particles)
        call RARRAY_COPY(5, qfl_s(1, k_impZ1), 1, rdum, 1)
        rdum(6) = RARRAY_SUM(5, rdum, 1)
        qcond_imp1_nc(j) = rdum(6)*VRS(j)*1.e-6
! Heat conduction and velocity (impZ1 particles)
! Conduction total is sum of components
        rdum(1) = RARRAY_SUM(5, qfl_s(1, k_impZ1), 1)
! Diagonal conductivity plus convective velocity
        dq_s(k_impZ1) = chit_ss(k_impZ1, k_impZ1) + chip_ss(k_impZ1, k_impZ1)
        vq_s(k_impZ1) = rdum(1)/(den_iz(im, iza)* z_j7kv*temp_i(im)) + &
                        dq_s(k_impZ1)*grt_i(im)/temp_i(im)
        ximp1_nc(j) = dq_s(k_impZ1)/y_grrho2
        cimp1_nc(j) = vq_s(k_impZ1)/y_grrho2
! Total radial energy flux (impZ1 particles)
        do k=1, 5
            rdum(k) = qfl_s(k, k_impZ1) + 2.5*gfl_s(k, k_impZ1)*temp_i(im)*z_j7kv
        enddo
        rdum(6) = RARRAY_SUM(5, rdum, 1)
        qen_imp1_nc(j) = rdum(6)*VRS(j)*1.e-6
! Bootstrap current on p'/p (impZ1 particles)
        rdum(1) = bsjbp_s(k_impZ1)
        rdum(2) = bsjbt_s(k_impZ1)
        bs_pimp1_nc(j) = -1.e-6*rdum(1)/BTOR
        bs_timp1_nc(j) = -1.e-6*rdum(2)/BTOR
! Poloidal flow velociry on outside midplane (impZ1 particles)
        polflow_imp1_nc(j) = eps_shear_fac * ( &
            utheta_s(1, 1, k_impZ1) + &
            utheta_s(1, 2, k_impZ1) + &
            utheta_s(1, 3, k_impZ1) )
    endif

!--- Impurity Z2 ---
    if (k_impZ2 > 0) then
! Total radial particle flux (impZ2 particles)
        call RARRAY_COPY(5, gfl_s(1, k_impZ2), 1, rdum, 1)
        rdum(6) = RARRAY_SUM(5, rdum, 1)
        gamma_imp2_nc(j) = rdum(6)*VRS(j)*1.e-19
! Particle diffusion and velocity (impZ2 particles)
        im  = jm_s(k_impZ2)
        iza = IABS(jz_s(k_impZ2))
        dn_imp2_nc(j) = dn_s(k_impZ2)/y_grrho2
        cn_imp2_nc(j) = (vn_s(k_impZ2) + veb_s(k_impZ2) + &
                         gfl_s(5, k_impZ2)/den_iz(im, iza))  / y_grrho2
! Radial conduction flux (impZ2 particles)
        call RARRAY_COPY(5, qfl_s(1, k_impZ2), 1, rdum, 1)
        rdum(6) = RARRAY_SUM(5, rdum, 1)
        qcond_imp2_nc(j) = rdum(6)*VRS(j)*1.e-6
! Heat conduction and velocity (impZ2 particles)
! Conduction total is sum of components
        rdum(1) = RARRAY_SUM(5, qfl_s(1, k_impZ2), 1)
! Diagonal conductivity plus convective velocity
        dq_s(k_impZ2) = chit_ss(k_impZ2, k_impZ2) + chip_ss(k_impZ2, k_impZ2)
        vq_s(k_impZ2) = rdum(1)/(den_iz(im, iza)* z_j7kv*temp_i(im)) + &
                        dq_s(k_impZ2)*grt_i(im)/temp_i(im)
        ximp2_nc(j) = dq_s(k_impZ2)/y_grrho2
        cimp2_nc(j) = vq_s(k_impZ2)/y_grrho2
! Total radial energy flux (impZ2 particles)
        do k=1, 5
            rdum(k) = qfl_s(k, k_impZ2) + 2.5*gfl_s(k, k_impZ2)*temp_i(im)*z_j7kv
        enddo
        rdum(6) = RARRAY_SUM(5, rdum, 1)
        qen_imp2_nc(j) = rdum(6)*VRS(j)*1.e-6
! Bootstrap current on p'/p (impZ2 particles)
        rdum(1) = bsjbp_s(k_impZ2)
        rdum(2) = bsjbt_s(k_impZ2)
        bs_pimp2_nc(j) = -1.e-6*rdum(1)/BTOR
        bs_timp2_nc(j) = -1.e-6*rdum(2)/BTOR
! Poloidal flow velociry on outside midplane (impZ2 particles)
        polflow_imp2_nc(j) = eps_shear_fac * ( &
            utheta_s(1, 1, k_impZ2) + &
            utheta_s(1, 2, k_impZ2) + &
            utheta_s(1, 3, k_impZ2) )
    endif

!--- Impurity Z3 ---
    if (k_impZ3 > 0) then
! Total radial particle flux (impZ3 particles)
        call RARRAY_COPY(5, gfl_s(1, k_impZ3), 1, rdum, 1)
        rdum(6) = RARRAY_SUM(5, rdum, 1)
        gamma_imp3_nc(j) = rdum(6)*VRS(j)*1.e-19
! Particle diffusion and velocity (impZ3 particles)
        im  = jm_s(k_impZ3)
        iza = IABS(jz_s(k_impZ3))
        dn_imp3_nc(j) = dn_s(k_impZ3)/y_grrho2
        cn_imp3_nc(j) = (vn_s(k_impZ3) + veb_s(k_impZ3) + &
                         gfl_s(5, k_impZ3)/den_iz(im, iza)) / y_grrho2
! Radial conduction flux (impZ3 particles)
        call RARRAY_COPY(5, qfl_s(1, k_impZ3), 1, rdum, 1)
        rdum(6) = RARRAY_SUM(5, rdum, 1)
        qcond_imp3_nc(j) = rdum(6)*VRS(j)*1.e-6
! Heat conduction and velocity (impZ3 particles)
! Conduction total is sum of components
        rdum(1) = RARRAY_SUM(5, qfl_s(1, k_impZ3), 1)
! Diagonal conductivity plus convective velocity
        dq_s(k_impZ3) = chit_ss(k_impZ3, k_impZ3) + chip_ss(k_impZ3, k_impZ3)
        vq_s(k_impZ3) = rdum(1)/(den_iz(im, iza)* z_j7kv*temp_i(im)) + &
                        dq_s(k_impZ3)*grt_i(im)/temp_i(im)
        ximp3_nc(j) = dq_s(k_impZ3)/y_grrho2
        cimp3_nc(j) = vq_s(k_impZ3)/y_grrho2
! Total radial energy flux (impZ3 particles)
        do k=1, 5
            rdum(k) = qfl_s(k, k_impZ3) + 2.5*gfl_s(k, k_impZ3)*temp_i(im)*z_j7kv
        enddo
        rdum(6) = RARRAY_SUM(5, rdum, 1)
        qen_imp3_nc(j) = rdum(6)*VRS(j)*1.e-6
! Bootstrap current on p'/p (impZ3 particles)
        rdum(1) = bsjbp_s(k_impZ3)
        rdum(2) = bsjbt_s(k_impZ3)
        bs_pimp3_nc(j) = -1.e-6*rdum(1)/BTOR
        bs_timp3_nc(j) = -1.e-6*rdum(2)/BTOR
! Poloidal flow velociry on outside midplane (impZ3 particles)
        polflow_imp3_nc(j) = eps_shear_fac * ( &
            utheta_s(1, 1, k_impZ3) + &
            utheta_s(1, 2, k_impZ3) + &
            utheta_s(1, 3, k_impZ3) )
    endif

!--- Miscellaneous Parameters ---

! Bootstrap current
    rdum(1) = p_bsjb
    jbs_nc(j)  = -1.e-6*rdum(1)/BTOR
! External current
    rdum(2) = p_exjb
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

!---------------------------------------------------------------------
subroutine ZBFAUX(GRRdB2, NGRTheta, YFM)
! Returns <(grad(rho)/B)2>
!         <n.grad(theta)>
!         and F_m - poloidal moments of geometric factor for PS viscosity
!                                  A.V.Zolotukhin 09-09-2003
!                  corrected by I.Yu. Senichenkov (IYS) February 2008
!---------------------------------------------------------------------

use const_inc, only: ROC, HRO, NA, NA1, GP2, BTOR, RTOR
use status_inc, only: BDB0, BDB02, IPOL, MU, TRIA, ELON, SHIF, &
    AMETR, RHO, DRODA

integer, parameter :: ntheta=32
integer :: j_theta, j_theta1, j_rho, j_mvisc
double precision :: YH, yametr, ydroda, yrho, I_str, THETA_cap,  &
    y_gamma, yr_theta, theta_s, yLambda, yLambdap, yDelta, yDeltap, &
    ytriang, ytriangp, COS_th, COS2_th, yd, y_sqg, y_sqg1, grr2dB2, &
    yBm1, Ynablarho2, yBmod, yBtheta, THETA_gen, theta1_s, &
    COS_th1, COS2_th1, yr_theta1, yd1, yBmod1, yBtheta1, y_sqg11, y_gamma1
double precision, dimension(1000) :: GRRdB2, NGRTheta, yVRR, yBm
double precision :: YFM(3, 1000)
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
    yLambdap = (ELON(j_rho+1) - ELON(j_rho))/YH

! Magnetic axis shift and its derivative
    yDelta   = SHIF(j_rho)
    yDeltap  = (SHIF(j_rho+1) - SHIF(j_rho))/YH

! Magnetic surface triangularity and its derivative
    ytriang  = TRIA(j_rho)                     !!**!!
    ytriangp = (TRIA(j_rho+1) - TRIA(j_rho))/YH  !!**!!

    GRRdB2(j_rho)   = 0.
    NGRTheta(j_rho) = 0.
    yVRR(j_rho)     = 0.
    yBm(j_rho)      = 0.
    do j_mvisc=1, 3
        YFM(j_mvisc, j_rho) = 0.
        yBlnBCOS(j_mvisc) = 0.
        yBlnBSIN(j_mvisc) = 0.
        yB2COS(j_mvisc)   = 0.
        yB2SIN(j_mvisc)   = 0.
    enddo

    do j_theta=NTHETA, 1, -1
        theta_s = GP2*j_theta/NTHETA
        COS_th  = COS(theta_s)
        COS2_th = COS_th*COS_th
        yr_theta = RTOR + yDelta + yametr*(COS_th - ytriang * (1. - COS2_th))
        yd = yLambda * (1.0 + yDeltap * Cos_th) + (1. - COS2_th) * &
            (yametr * yLambdap * (1.0 + 2.0 * ytriang * Cos_th) +  &
             yLambda * (ytriang - yametr * ytriangp) * Cos_th)

        y_sqg = yr_theta*yametr*yd ! Square root of the determinant g
        y_sqg1 = yr_theta*yd ! sqrt(g)/rho

! B - magnetic field

        yBmod = I_str/(yr_theta*yd)*sqrt(yd**2 + THETA_cap**2 * &
                ((1.0 + 2.0 * ytriang * Cos_th)**2 * (1. - COS2_th) + &
                yLambda**2 * Cos2_th))
! 1/B
        yBm1 = 1./yBmod
! B_theta - contravariant component of the magnetic field
        yBtheta = BTOR*MU(j_rho)*yrho*ydroda/y_sqg
! nabla(rho)2
        ynablarho2 = ydroda**2 * ((1.0 + 2.0 * ytriang * Cos_th)**2 *  &
                    (1. - COS2_th) + yLambda**2 * Cos2_th) / yd**2
! nabla(rho)2/B2
        grr2dB2 = ynablarho2/(yBmod*yBmod)
! Integral[sqrt(g)/rho, 0, 2pi]
        yVRR(j_rho) = yVRR(j_rho) + y_sqg1*GP2/NTHETA
! Integral[(nabla(rho)/B)2*SQRT(g)/rho, 0, 2pi]
        GRRdB2(j_rho)   = GRRdB2(j_rho)+ grr2dB2*y_sqg1*GP2/NTHETA
! Integral[1/B, 0, 2pi]
        THETA_gen = 0.
        do j_theta1=NTHETA, 1, -1
            theta1_s   = theta_s*j_theta1/NTHETA
            COS_th1    = COS(theta1_s)
            COS2_th1   = COS_th1*COS_th1
            yr_theta1 = RTOR + yDelta + yametr*(COS_th1 - ytriang * (1. - COS2_th1))
            yd1 = yLambda * (1.0 + yDeltap * Cos_th1) + (1. - COS2_th1) * &
                 (yametr * yLambdap * (1.0 + 2.0 * ytriang * Cos_th) +  &
                  yLambda * (ytriang - yametr * ytriangp) * Cos_th1)
            y_sqg11    = yr_theta1*yd1

! B1 - magnetic field
            yBmod1 = I_str/(yr_theta1*yd1)*sqrt(yd1**2 + THETA_cap**2 * &
                     ((1.0 + 2.0 * ytriang * Cos_th1)**2 * (1. - COS2_th1) + &
                     yLambda**2 * Cos2_th1))
! B_theta1 - contravariant component of the magnetic field
            yBtheta1 = BTOR*MU(j_rho)*yrho*ydroda/y_sqg11/yametr
            THETA_gen = THETA_gen + yBmod1/yBtheta1*theta_s/NTHETA
            if (j_theta == NTHETA) y_gamma1 = THETA_gen
        enddo  ! j_theta1

        if (j_theta == NTHETA) y_gamma = GP2/y_gamma1
        THETA_gen = y_gamma*THETA_gen

        do j_mvisc=1, 3
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

    do j_mvisc=1, 3
        yBlnBCOS(j_mvisc) = yBlnBCOS(j_mvisc)/yVRR(j_rho)
        yBlnBSIN(j_mvisc) = yBlnBSIN(j_mvisc)/yVRR(j_rho)
        yB2COS(j_mvisc)   = yB2COS(j_mvisc)/yVRR(j_rho)
        yB2SIN(j_mvisc)   = yB2SIN(j_mvisc)/yVRR(j_rho)
        YFM(j_mvisc, j_rho) = 2.*(j_mvisc*y_gamma)**2/ &
                             (BTOR**3*BDB02(j_rho)*BDB0(j_rho)) * &
                             ( yBlnBCOS(j_mvisc)*yB2COS(j_mvisc) + &
                               yBlnBSIN(j_mvisc)*yB2SIN(j_mvisc) )
    enddo

enddo

GRRdB2(NA1) = GRRdB2(NA)*(ROC/HRO - (NA-1)) + GRRdB2(NA-1)*(NA - ROC/HRO)
NGRTheta(NA1) = NGRTheta(NA)*(ROC/HRO - (NA-1)) + NGRTheta(NA-1)*(NA - ROC/HRO)
do j_mvisc=1, 3
    YFM(j_mvisc, NA1) = YFM(j_mvisc, NA)*(ROC/HRO - (NA-1)) + &
                        YFM(j_mvisc, NA-1)*(NA - ROC/HRO)
enddo

return
end subroutine zbfaux

!----------------------------------------------------------------
subroutine NCLASS(k_order,k_potato,m_i,m_z,c_den,c_potb,c_potl, &
                  p_b2,p_bm2,p_eb,p_fhat,p_fm,p_ft,p_grbm2, &
                  p_grphi,p_gr2phi,p_ngrth,amu_i,grt_i,temp_i, &
                  den_iz,fex_iz,grp_iz,m_s,jm_s,jz_s,p_bsjb, &
                  p_etap,p_exjb,calm_i,caln_ii,capm_ii,capn_ii, &
                  bsjbp_s,bsjbt_s,dn_s,gfl_s,qfl_s,sqz_s,upar_s, &
                  utheta_s,vn_s,veb_s,qeb_s,xi_s,ymu_s,chip_ss, &
                  chit_ss,dp_ss,dt_ss,iflag)
!----------------------------------------------------------------
!NCLASS calculates the neoclassical transport properties of a multiple
!  species axisymmetric plasma using k_order parallel and radial force
!  balance equations for each species
!References:
!  Houlberg, Shaing, Hirshman, Zarnstorff, Phys Plasmas 4 (1997) 3230
!  Hirshman, Sigmar, Nucl Fusion 21 (1981) 1079
!  W.A.Houlberg 6/99
!Input:
!  k_order-order of v moments to be solved (-)
!       =2 u and q
!       =3 u, q, and u2
!       =else error
!  k_potato-option to include potato orbits (-)
!          =0 off
!          =else on
!  m_i-number of isotopes (1<m_i<mx_mi+1)
!  m_z-highest charge state of all species (0<m_z<mx_mz+1)
!  c_den-density cutoff below which species is ignored (/m**3)
!  c_potb-kappa(0)*Bt(0)/[2*q(0)**2] (T)
!  c_potl-q(0)*R(0) (m)
!  p_b2-<B**2> (T**2)
!  p_bm2-<1/B**2> (/T**2)
!  p_eb-<E.B> (V*T/m)
!  p_fhat-mu_0*F/(dPsi/dr) (rho/m)
!  p_fm(3)-poloidal moments of geometric factor for PS viscosity (-)
!  p_ft-trapped fraction (-)
!  p_grbm2-<grad(rho)**2/B**2> (rho**2/m**2/T**2)
!  p_grphi-radial electric field Phi' (V/rho)
!  p_gr2phi-radial electric field gradient Psi'(Phi'/Psi')' (V/rho**2)
!  p_ngrth-<n.grad(Theta)> (1/m)
!  amu_i(i)-atomic mass number of i (-)
!  grt_i(i)-temperature gradient of i (keV/rho)
!  temp_i(i)-temperature of i (keV)
!  den_iz(i,z)-density of i,z (/m**3)
!  fex_iz(3,i,z)-moments of external parallel force on i,z (T*j/m**3)
!  grp_iz(i,z)-pressure gradient of i,z (keV/m**3/rho)
!Output:
!  m_s-number of species (1<ms<mx_ms+1)
!  jm_s(s)-isotope number of s (-)
!  jz_s(s)-charge state of s (-)
!  p_bsjb-<J_bs.B> (A*T/m**2)
!  p_etap-parallel electrical resistivity (Ohm*m)
!  p_exjb-<J_ex.B> current response to fex_iz (A*T/m**2)
!  calm_i(3,3,i)-tp eff friction matrix for i (kg/m**3/s)
!  caln_ii(3,3,i1,i2)-fp eff friction matrix for i1 on i2 (kg/m**3/s)
!  capm_ii(3,3,i1,i2)-test part (tp) friction matrix for i1 on i2 (-)
!  capn_ii(3,3,i1,i2)-field part (fp) friction matrix for i1 on i2 (-)
!  bsjbp_s(s)-<J_bs.B> driven by unit p'/p of s (A*T*rho/m**3)
!  bsjbt_s(s)-<J_bs.B> driven by unit T'/T of s (A*T*rho/m**3)
!  dn_s(s)-diffusion coefficient (diag comp) of s (rho**2/s)
!  gfl_s(m,s)-radial particle flux comps of s (rho/m**3/s)
!             m=1, banana-plateau, p' and T'
!             m=2, Pfirsch-Schluter
!             m=3, classical
!             m=4, banana-plateau, <E.B>
!             m=5, banana-plateau, external parallel force fex_iz
!  qfl_s(m,s)-radial heat conduction flux comps of s (W*rho/m**3)
!             m=1, banana-plateau, p' and T'
!             m=2, Pfirsch-Schluter
!             m=3, classical
!             m=4, banana-plateau, <E.B>
!             m=5, banana-plateau, external parallel force fex_iz
!  sqz_s(s)-orbit squeezing factor for s (-)
!  upar_s(3,m,s)-parallel flow of s from force m (T*m/s)
!             m=1, p', T', Phi'
!             m=2, <E.B>
!             m=3, fex_iz
!  utheta_s(3,m,s)-poloidal flow of s from force m (m/s/T)
!             m=1, p', T'
!             m=2, <E.B>
!             m=3, fex_iz
!  vn_s(s)-convection velocity (off diag comps-p', T') of s (rho/s)
!  veb_s(s)-<E.B> particle convection velocity of s (rho/s)
!  qeb_s(s)-<E.B> heat convection velocity of s (rho/s)
!  xi_s(s)-charge weighted density factor of s (-)
!  ymu_s(s)-normalized viscosity for s (kg/m**3/s)
!  chip_ss(s1,s2)-heat cond coefficient of s2 on p'/p of s1 (rho**2/s)
!  chit_ss(s1,s2)-heat cond coefficient of s2 on T'/T of s1 (rho**2/s)
!  dp_ss(s1,s2)-diffusion coefficient of s2 on p'/p of s1 (rho**2/s)
!  dt_ss(s1,s2)-diffusion coefficient of s2 on T'/T of s1 (rho**2/s)
!  iflag-warning and error flag
!       =-4 warning: no viscosity
!       =-3 warning: no banana viscosity
!       =-2 warning: no Pfirsch-Schluter viscosity
!       =-1 warning: no potato orbit viscosity
!       =0 no warnings or errors
!       =1 error: order of v moments to be solved must be 2 or 3
!       =2 error: number of species must be 1<m_i<mx_mi+1
!       =3 error: number of species must be 0<m_z<mx_mz+1
!       =4 error: number of species must be 1<m_s<mx_ms+1
!       =5 error: inversion of flow matrix failed
!       =6 error: trapped fraction must be 0.0.le.p_ft.le.1.0

!Declaration of input variables
integer        k_order,                 k_potato
integer        m_i,                     m_z
real           c_den,                   c_potb, &
            c_potl
real           p_b2,                    p_bm2, &
            p_eb,                    p_fhat, &
            p_fm(3),                 p_ft, &
            p_grbm2,                 p_grphi, &
            p_gr2phi,                p_ngrth
real           amu_i(mx_mi),            grt_i(mx_mi), &
            temp_i(mx_mi)
real           den_iz(mx_mi,mx_mz),     fex_iz(3,mx_mi,mx_mz), &
            grp_iz(mx_mi,mx_mz)
!Declaration of output variables
integer        iflag,                   m_s
integer        jm_s(mx_ms),             jz_s(mx_ms)
real           p_bsjb,                  p_etap, &
            p_exjb
real           calm_i(3,3,mx_mi)
real           caln_ii(3,3,mx_mi,mx_mi),capm_ii(3,3,mx_mi,mx_mi), &
            capn_ii(3,3,mx_mi,mx_mi)
real           bsjbp_s(mx_ms),          bsjbt_s(mx_ms), &
            dn_s(mx_ms),             gfl_s(5,mx_ms), &
            qfl_s(5,mx_ms),          sqz_s(mx_ms), &
            upar_s(3,3,mx_ms),       utheta_s(3,3,mx_ms), &
            vn_s(mx_ms),             veb_s(mx_ms), &
            qeb_s(mx_ms),            xi_s(mx_ms), &
            ymu_s(3,3,mx_ms)
real           chip_ss(mx_ms,mx_ms),    chit_ss(mx_ms,mx_ms), &
            dp_ss(mx_ms,mx_ms),      dt_ss(mx_ms,mx_ms)
!Declaration of local variables
integer        k_banana,                k_pfirsch
integer        i,                       im, &
            iz,                      iza, &
            jflag,                   jm, &
            k,                       l
real           dent
real           z_coulomb,               z_electronmass, &
            z_j7kv,                  z_mu0, &
            z_pi,                    z_protonmass
real           denz2(mx_mi),            vt_i(mx_mi)
real           pgrp_iz(mx_mi,mx_mz)
real           amnt_ii(mx_mi,mx_mi)
real           tau_ss(mx_ms,mx_ms)
!Initialization
!  Error flag
iflag=0
!  Consistency checks
!     Order must be two or three moments
if(k_order.lt.2.or.k_order.gt.3) then
  iflag=1
  return
endif
!     At least two but not greater than mx_mi species
if(m_i.lt.2.or.m_i.gt.mx_mi) then
  iflag=2
  return
endif
!     Highest charge state at least 1 but not greater than mx_mz
if(m_z.lt.1.or.m_i.gt.mx_mz) then
  iflag=3
  return
endif
!     Trapped fraction between 0 and 1 inclusive
if((p_ft.lt.0.0).or.(p_ft.gt.1.0)) then
  iflag=6
  return
endif
!     Potato orbit contribution to viscosity
if((ABS(c_potb).gt.0.0).and.(ABS(c_potl).gt.0.0) &
.and.(k_potato.ne.0)) then
  k_potato=1
else
  k_potato=0
  iflag=-1
endif
!     Pfirsch-Schluter contribution to viscosity
if(ABS(p_fm(1)+p_fm(2)+p_fm(3)).gt.0.0) then
  k_pfirsch=1
else
  k_pfirsch=0
  iflag=-2
endif
!     Banana contribution to viscosity
if(ABS(p_ft).gt.0.0) then
  k_banana=1
else
  k_banana=0
  iflag=-3
endif
!     No viscsoity
if((k_banana.eq.0).and.(k_pfirsch.eq.0)) then
  k_potato=0
  iflag=-4
endif
!  Physical and conversion constants
z_coulomb=1.6022e-19
z_electronmass=9.1095e-31
z_j7kv=1.6022e-16
z_mu0=1.2566e-06
z_pi=ACOS(-1.0)
z_protonmass=1.6726e-27
!Find significant charge states and mapping
m_s=0
do im=1,m_i
  do iza=1,m_z
    if(den_iz(im,iza).gt.c_den) then
      m_s=m_s+1
!           Set isotope number and charge state for this species
      jm_s(m_s)=im
      if(amu_i(im).lt.0.5) then
        jz_s(m_s)=-iza
      else
        jz_s(m_s)=iza
      endif
    endif
  enddo
enddo
if(m_s.lt.2.or.m_s.gt.mx_ms) then
  iflag=4
  return
endif
!Get friction coefficients
call NCLASS_MN(k_order,m_i,amu_i,temp_i,capm_ii,capn_ii)
!Calculate thermal velocity
do im=1,m_i
  vt_i(im)=SQRT(2.0*z_j7kv*temp_i(im)/amu_i(im)/z_protonmass)
enddo
!Get collision times
call NCLASS_TAU(m_i,m_s,jm_s,jz_s,amu_i,temp_i,vt_i,den_iz, &
             amnt_ii,tau_ss)
!Calculate reduced friction coefficients
call RARRAY_ZERO(9*m_i,calm_i)
do im=1,m_i
  do jm=1,m_i
    do k=1,k_order
      do l=1,k_order
!             Sum over isotopes b for test particle component
        calm_i(k,l,im)=calm_i(k,l,im) &
                   +amnt_ii(im,jm)*capm_ii(k,l,im,jm)
!             Field particle component.
        caln_ii(k,l,im,jm)=amnt_ii(im,jm)*capn_ii(k,l,im,jm)            
      enddo   
    enddo    
  enddo   
enddo   
!Calculate species charge state density factor, total nT, squeezing
dent=0.0
do im=1,m_i
  denz2(im)=0.0
  do iza=1,m_z
    if(den_iz(im,iza).gt.c_den) then
      denz2(im)=denz2(im)+den_iz(im,iza)*iza**2
      dent=dent+den_iz(im,iza)*temp_i(im)
    endif
  enddo   
enddo   
do i=1,m_s                                                            
  im=jm_s(i)
  iz=jz_s(i)
  iza=IABS(iz)
  xi_s(i)=den_iz(im,iza)*iz**2/denz2(im)
  sqz_s(i)=1.0+p_fhat**2/p_b2*amu_i(im)*z_protonmass &
            *ABS(p_gr2phi/(z_coulomb*iz))
enddo
!Get normalized viscosities
call NCLASS_MU(k_order,k_banana,k_pfirsch,k_potato,m_s,jm_s,jz_s, &
            c_potb,c_potl,p_fm,p_ft,p_ngrth,amu_i,temp_i,vt_i, &
            den_iz,sqz_s,ymu_s,tau_ss)
!Add potential gradient to pressure gradient
do i=1,m_s
  im=jm_s(i)
  iz=jz_s(i)
  iza=IABS(iz)
  pgrp_iz(im,iza)=grp_iz(im,iza) &
               +p_grphi*den_iz(im,iza)*iz*z_coulomb/z_j7kv
enddo
!Get normalized parallel flows within a surface
jflag=0
call NCLASS_FLOW(k_order,m_i,m_s,jm_s,jz_s,p_b2,p_bm2,p_eb,p_fhat, &
              p_grbm2,grt_i,temp_i,calm_i,caln_ii,den_iz, &
              fex_iz,pgrp_iz,xi_s,ymu_s,p_bsjb,p_etap,p_exjb, &
              bsjbp_s,bsjbt_s,gfl_s,qfl_s,upar_s,chip_ss, &
              chit_ss,dp_ss,dt_ss,jflag)
if(jflag.gt.0) then
  iflag=5
  return
endif
!Calculate poloidal velocity from parallel velocity
do i=1,m_s
  im=jm_s(i)
  iz=jz_s(i)
  iza=IABS(iz)
  do k=1,k_order
    do l=1,3
      utheta_s(k,l,i)=upar_s(k,l,i)/p_b2
      if(l.eq.1) then
        if(k.eq.1) then
          utheta_s(k,l,i)=utheta_s(k,l,i) &
                       +p_fhat*pgrp_iz(im,iza)*z_j7kv &
                       /(z_coulomb*iz*den_iz(im,iza))/p_b2
        elseif(k.eq.2) then
          utheta_s(k,l,i)=utheta_s(k,l,i) &
                       +p_fhat*z_j7kv*grt_i(im) &
                       /(iz*z_coulomb*p_b2)
        endif
      endif
    enddo
  enddo
enddo
!Convert to diffusivities and conductivities
!  Full coefficient matrices 
do i=1,m_s
  im=jm_s(i)
  iza=IABS(jz_s(i))
!  Diagonal diffusivity
  dn_s(i)=dp_ss(i,i)
!  Off-diagonal expressed as radial velocity
  vn_s(i)=(gfl_s(1,i)+gfl_s(2,i)+gfl_s(3,i) &
       +dn_s(i)*(grp_iz(im,iza)-den_iz(im,iza)*grt_i(im)) &
       /temp_i(im))/den_iz(im,iza)
!  <E.B> particle and heat convection velocities
  veb_s(i)=gfl_s(4,i)/den_iz(im,iza)
  qeb_s(i)=qfl_s(4,i)/den_iz(im,iza)/temp_i(im)/z_j7kv
enddo

return
end subroutine nclass

!----------------------------------------------------------------
subroutine NCLASS_FLOW(k_order,m_i,m_s,jm_s,jz_s,p_b2,p_bm2,p_eb, &
                       p_fhat,p_grbm2,grt_i,temp_i,calm_i,caln_ii, &
                       den_iz,fex_iz,grp_iz,xi_s,ymu_s,p_bsjb, &
                       p_etap,p_exjb,bsjbp_s,bsjbt_s,gfl_s,qfl_s, &
                       upar_s,chip_ss,chit_ss,dp_ss,dt_ss,iflag)
!----------------------------------------------------------------
!NCLASS_FLOW calculates the k_order neoclassical flows u and q/p (and
!  u2) plus other transport properties
!References:                                                     
!  Hirshman, Sigmar, Nucl Fusion 21 (1981) 1079
!  Houlberg, Shaing, Hirshman, Zarnstorff, Phys Plasmas 4 (1997) 3230
!  W.A.Houlberg 6/99
!Input:
!  k_order-order of v moments to be solved (-)
!       =2 u and q
!       =3 u, q, and u2
!       =else error
!  m_i-number of isotopes (1<m_i<mx_mi+1)
!  m_s-number of species (1<ms<mx_ms+1)
!  jm_s(s)-isotope number of s (-)
!  jz_s(s)-charge state of s (-)
!  p_b2-<B**2> (T**2)
!  p_bm2-<1/B**2> (/T**2)
!  p_eb-<E.B> (V*T/m)
!  p_fhat-mu_0*F/(dPsi/dr) (rho/m)
!  p_grbm2-<grad(rho)**2/B**2> (rho**2/m**2/T**2)
!  grt_i(i)-temperature gradient of i (keV/rho)
!  temp_i(i)-temperature of i (keV)
!  calm_i(3,3,i)-tp eff friction matrix for i (kg/m**3/s)
!  caln_ii(3,3,i1,i2)-fp eff friction matrix for i1 on i2 (kg/m**3/s)
!  den_iz(i,z)-density of i,z (/m**3)
!  fex_iz(3,i,z)-moments of external parallel force on i,z (T*j/m**3)
!  grp_iz(i,z)-pressure gradient of i,z (keV/m**3/rho)
!  xi_s(s)-charge weighted density factor of s (-)
!  ymu_s(s)-normalized viscosity for s (kg/m**3/s)
!Output:
!  p_bsjb-<J_bs.B> (A*T/m**2)
!  p_etap-parallel electrical resistivity (Ohm*m)
!  p_exjb-<J_ex.B> current response to fex_iz (A*T/m**2)
!  bsjbp_s(s)-<J_bs.B> driven by unit p'/p of s (A*T*rho/m**3)
!  bsjbt_s(s)-<J_bs.B> driven by unit T'/T of s (A*T*rho/m**3)
!  gfl_s(m,s)-radial particle flux comps of s (rho/m**3/s)
!             m=1, banana-plateau, p' and T'
!             m=2, Pfirsch-Schluter
!             m=3, classical
!             m=4, banana-plateau, <E.B>
!             m=5, banana-plateau, external parallel force fex_iz
!  qfl_s(m,s)-radial heat conduction flux comps of s (W*rho/m**3)
!             m=1, banana-plateau, p' and T'
!             m=2, Pfirsch-Schluter
!             m=3, classical
!             m=4, banana-plateau, <E.B>
!             m=5, banana-plateau, external parallel force fex_iz
!  upar_s(3,m,s)-parallel flow of s from force m (T*m/s)
!                m=1, p' and T'
!                m=2, <E.B>
!                m=3, fex_iz
!  chip_ss(s1,s2)-heat cond coefficient of s2 on p'/p of s1 (rho**2/s)
!  chit_ss(s1,s2)-heat cond coefficient of s2 on T'/T of s1 (rho**2/s)
!  dp_ss(s1,s2)-diffusion coefficient of s2 on p'/p of s1 (rho**2/s)
!  dt_ss(s1,s2)-diffusion coefficient of s2 on T'/T of s1 (rho**2/s)
!  iflag-error flag
!       =0 no errors
!       =1 inversion of flow matrix failed
!----------------------------------------------------------------

!Declaration of input variables
integer        k_order,                 m_i, &
            m_s
integer        jm_s(mx_ms),             jz_s(mx_ms)       
real           p_b2,                    p_bm2, &
            p_eb,                    p_fhat,                     &
            p_grbm2
real           temp_i(mx_mi),           grt_i(mx_mi), &
            calm_i(3,3,mx_mi)
real           caln_ii(3,3,mx_mi,mx_mi)
real           den_iz(mx_mi,mx_mz),     grp_iz(mx_mi,mx_mz), &
            fex_iz(3,mx_mi,mx_mz)
real           xi_s(mx_ms),             ymu_s(3,3,mx_ms)
!Declaration of output variables
integer        iflag
real           p_bsjb,                  p_etap, &
            p_exjb
real           bsjbp_s(mx_ms),          bsjbt_s(mx_ms), &
            gfl_s(5,mx_ms),          qfl_s(5,mx_ms), &
            upar_s(3,3,mx_ms)
real           chip_ss(mx_ms,mx_ms),    chit_ss(mx_ms,mx_ms), &
            dp_ss(mx_ms,mx_ms),      dt_ss(mx_ms,mx_ms)
!Declaration of local variables
integer        i,                       im, &
            iz,                      iza, &
            j,                       jm, &
            k,                       l,                      &
            l1,                      m, &
            m1
integer        indx(3),                 indxb(3*mx_mi)
real           cbp,                     cbpa, &
            cbpaq,                   cc, &
            ccl,                     ccla, &
            cclaq,                   cclb, &
            cclbq,                   ceb, &
            cps,                     cpsa, &
            cpsaq,                   cpsb, &
            cpsbq,                   d, &
            denzc,                   p_ohjb, &
            z_coulomb,               z_j7kv
real           aa(3,3),                 xl(3)
real           crhat(3,6,mx_mi),        rhat(3,6,mx_ms)
real           srcth(3,mx_ms),          srcthp(mx_ms),             &
            srctht(mx_ms)
real           ab(3*mx_mi,3*mx_mi),     xab(3*mx_mi,3)
real           crhatp(3,mx_ms,mx_mi),   crhatt(3,mx_ms,mx_mi)
real           rhatp(3,mx_ms,mx_ms),    rhatt(3,mx_ms,mx_ms), &
            uaip(3,mx_ms,mx_ms),     uait(3,mx_ms,mx_ms)
real           xabp(3*mx_mi,mx_ms),     xabt(3*mx_mi,mx_ms)
!Initialization
!  Error flag
iflag=0
!  Physical and conversion constants
z_coulomb=1.6022e-19
z_j7kv=1.6022e-16
!  Zero out arrays
call RARRAY_ZERO(3*6*m_i,crhat)
call RARRAY_ZERO(3*mx_ms*m_i,crhatp)
call RARRAY_ZERO(3*mx_ms*m_i,crhatt)
call RARRAY_ZERO(3*mx_mi*3*m_i,ab)
call RARRAY_ZERO(5*m_s,gfl_s)
call RARRAY_ZERO(5*m_s,qfl_s)
p_etap=0.0
p_bsjb=0.0
p_exjb=0.0
p_ohjb=0.0
cc=(p_fhat/z_coulomb)*z_j7kv
cbp=p_fhat/p_b2/z_coulomb
cps=(p_fhat/z_coulomb)*(1.0/p_b2-p_bm2)
ccl=(p_grbm2/z_coulomb)/p_fhat
ceb=p_fhat/p_b2*p_eb
!Calculate responses for each species
do i=1,m_s
  im=jm_s(i)
  iz=jz_s(i)
  iza=IABS(iz)
!  Set up response matrix for each charge state 
  do k=1,k_order
    do l=1,k_order
      aa(k,l)=xi_s(i)*calm_i(k,l,im)-ymu_s(k,l,i)
    enddo
  enddo
!  Get lu decomposition of response matrix
  call U_LU_DECOMP(aa,k_order,3,indx,d,iflag)
  if(iflag.ne.0) return
!  Get sources and evaluate responses from back substitution 
!       Lambda terms involving isotopic flows 
  do l=1,k_order
    do k=1,k_order
      if(k.eq.l) then
        rhat(k,l,i)=xi_s(i)
      else
        rhat(k,l,i)=0.0
      endif
    enddo
    call U_LU_BACKSUB(aa,k_order,3,indx,rhat(1,l,i))
    do k=1,k_order
      crhat(k,l,im)=crhat(k,l,im)+xi_s(i)*rhat(k,l,i)
    enddo
  enddo
!       Poloidal source (p' and T') terms 
  srcth(1,i)=(cc/iz)*grp_iz(im,iza)/den_iz(im,iza)
  srcth(2,i)=(cc/iz)*grt_i(im)
  srcth(3,i)=0.0
  do k=1,k_order
    rhat(k,4,i)=0.0
    do l=1,k_order
      rhat(k,4,i)=rhat(k,4,i)+srcth(l,i)*ymu_s(k,l,i)
    enddo
  enddo
  call U_LU_BACKSUB(aa,k_order,3,indx,rhat(1,4,i))
  do k=1,k_order
    crhat(k,4,im)=crhat(k,4,im)+xi_s(i)*rhat(k,4,i)
  enddo
!       Unit p'/p and T'/T terms for decomposition of fluxes 
  srcthp(i)=-(cc/iz)*temp_i(im)
  srctht(i)=-(cc/iz)*temp_i(im)
  do k=1,k_order
    rhatp(k,i,i)=srcthp(i)*ymu_s(k,1,i)
    rhatt(k,i,i)=srctht(i)*ymu_s(k,2,i)
  enddo
  call U_LU_BACKSUB(aa,k_order,3,indx,rhatp(1,i,i))
  call U_LU_BACKSUB(aa,k_order,3,indx,rhatt(1,i,i))
  do k=1,k_order
    crhatp(k,i,im)=crhatp(k,i,im)+xi_s(i)*rhatp(k,i,i)
    crhatt(k,i,im)=crhatt(k,i,im)+xi_s(i)*rhatt(k,i,i)
  enddo
!       Parallel electric field terms for resistivity 
  rhat(1,5,i)=-iz*z_coulomb*den_iz(im,iza)
  rhat(2,5,i)=0.0
  rhat(3,5,i)=0.0
  call U_LU_BACKSUB(aa,k_order,3,indx,rhat(1,5,i))
  do k=1,k_order
    crhat(k,5,im)=crhat(k,5,im)+xi_s(i)*p_eb*rhat(k,5,i)
  enddo
!       External force terms 
  rhat(1,6,i)=-fex_iz(1,im,iza)
  rhat(2,6,i)=-fex_iz(2,im,iza)
  if(k_order.eq.3) then
    rhat(3,6,i)=-fex_iz(3,im,iza)
  else
    rhat(3,6,i)=0.0
  endif
  call U_LU_BACKSUB(aa,k_order,3,indx,rhat(1,6,i))
  do k=1,k_order
    crhat(k,6,im)=crhat(k,6,im)+xi_s(i)*rhat(k,6,i)
  enddo
enddo   
!Load coefficient matrix and source terms for isotopic flows
do im=1,m_i
  do m=1,k_order
    m1=im+(m-1)*m_i
!  Diagonal coefficients       
    ab(m1,m1)=1.0
!  Source terms 
!         p' and T'
    xab(m1,1)=crhat(m,4,im)
!         Unit p'/p and T'/T
    do j=1,m_s
      xabp(m1,j)=crhatp(m,j,im)
      xabt(m1,j)=crhatt(m,j,im)
    enddo
!         <E.B>
    xab(m1,2)=crhat(m,5,im)
!         External source
    xab(m1,3)=crhat(m,6,im)
!  Field particle friction       
    do jm=1,m_i
      do l=1,k_order
        l1=jm+(l-1)*m_i
        do k=1,k_order
          ab(m1,l1)=ab(m1,l1)+caln_ii(k,l,im,jm)*crhat(m,k,im)
        enddo  
      enddo
    enddo
  enddo        
enddo   
!Get lu decomposition of coefficient matrix 
call U_LU_DECOMP(ab,k_order*m_i,3*mx_mi,indxb,d,iflag)
if(iflag.ne.0) return
!Evaluate isotopic flows from back substitution for each source 
!  xab(1,k) to xab(m_i,k) are the isotopic velocities 
!  xab(m_i+1,k) to xab(2*m_i,k) are the isotopic heat flows 
!  xab(2*m_i+1,k) to xab(3*m_i,k) are the u2 flows 
!  Evaluate species flows 
do k=1,3
  call U_LU_BACKSUB(ab,k_order*m_i,3*mx_mi,indxb,xab(1,k))
enddo
do i=1,m_s
  call U_LU_BACKSUB(ab,k_order*m_i,3*mx_mi,indxb,xabp(1,i))
  call U_LU_BACKSUB(ab,k_order*m_i,3*mx_mi,indxb,xabt(1,i))
enddo
do i=1,m_s
  im=jm_s(i)
!  Source contributions 
  do m=1,3        
    do k=1,k_order
      if(m.eq.2) then
        upar_s(k,m,i)=p_eb*rhat(k,5,i)
      else
        upar_s(k,m,i)=rhat(k,m+3,i)
      endif
    enddo
!         Response contributions           
    do jm=1,m_i
      call RARRAY_ZERO(k_order,xl)
      do l=1,k_order
        l1=jm+(l-1)*m_i
        do k=1,k_order
          xl(k)=xl(k)-caln_ii(k,l,im,jm)*xab(l1,m)
        enddo
      enddo
      do l=1,k_order
        do k=1,k_order
          upar_s(k,m,i)=upar_s(k,m,i)+xl(l)*rhat(k,l,i)
        enddo
      enddo
    enddo
  enddo
!  Unit p'/p and T'/T
  do j=1,m_s
    do k=1,k_order
      uaip(k,j,i)=rhatp(k,j,i)
      uait(k,j,i)=rhatt(k,j,i)
    enddo
!         Response contributions           
    do jm=1,m_i
      call RARRAY_ZERO(k_order,xl)
      do l=1,k_order
        l1=jm+(l-1)*m_i
        do k=1,k_order
          xl(k)=xl(k)-caln_ii(k,l,im,jm)*xabp(l1,j)
        enddo
      enddo
      do l=1,k_order
        do k=1,k_order
          uaip(k,j,i)=uaip(k,j,i)+xl(l)*rhat(k,l,i)
        enddo
      enddo
      call RARRAY_ZERO(k_order,xl)
      do l=1,k_order
        l1=jm+(l-1)*m_i
        do k=1,k_order
          xl(k)=xl(k)-caln_ii(k,l,im,jm)*xabt(l1,j)
        enddo
      enddo
      do l=1,k_order
        do k=1,k_order
          uait(k,j,i)=uait(k,j,i)+xl(l)*rhat(k,l,i)
        enddo
      enddo
    enddo
  enddo
enddo
!Currents and fluxes
call RARRAY_ZERO(m_s,bsjbp_s)
call RARRAY_ZERO(m_s,bsjbt_s)
call RARRAY_ZERO(mx_ms*m_s,dp_ss)
call RARRAY_ZERO(mx_ms*m_s,dt_ss)
call RARRAY_ZERO(mx_ms*m_s,chip_ss)
call RARRAY_ZERO(mx_ms*m_s,chit_ss)
do i=1,m_s
  im=jm_s(i)
  iz=jz_s(i)
  iza=IABS(iz)
!<J_bs.B> 
  denzc=den_iz(im,iza)*iz*z_coulomb
  p_bsjb=p_bsjb+denzc*upar_s(1,1,i)
!<J_OH.B>
  p_ohjb=p_ohjb+denzc*upar_s(1,2,i)
!<J_ex.B>
  p_exjb=p_exjb+denzc*upar_s(1,3,i)
!  Unit p'/p and T'/T
  do j=1,m_s
    bsjbp_s(j)=bsjbp_s(j)+denzc*uaip(1,j,i)
    bsjbt_s(j)=bsjbt_s(j)+denzc*uait(1,j,i)
  enddo
!Fluxes
!  Banana-Plateau
  cbpa=cbp/iz
  cbpaq=cbpa*(z_j7kv*temp_i(im))
!       Unit p'/p and T'/T
  dp_ss(i,i)=dp_ss(i,i)-cbpa*ymu_s(1,1,i)*srcthp(i)
  dt_ss(i,i)=dt_ss(i,i)-cbpa*ymu_s(1,2,i)*srctht(i)
  chip_ss(i,i)=chip_ss(i,i)-cbpaq*ymu_s(2,1,i)*srcthp(i)
  chit_ss(i,i)=chit_ss(i,i)-cbpaq*ymu_s(2,2,i)*srctht(i)
  do k=1,k_order
!         p' and T'
    gfl_s(1,i)=gfl_s(1,i)-cbpa*ymu_s(1,k,i) &
                       *(upar_s(k,1,i)+srcth(k,i))
    qfl_s(1,i)=qfl_s(1,i)-cbpaq*ymu_s(2,k,i) &
                       *(upar_s(k,1,i)+srcth(k,i))
!         Unit p'/p and T'/T            
    do j=1,m_s
      dp_ss(j,i)=dp_ss(j,i)-cbpa*ymu_s(1,k,i)*uaip(k,j,i)
      dt_ss(j,i)=dt_ss(j,i)-cbpa*ymu_s(1,k,i)*uait(k,j,i)
      chip_ss(j,i)=chip_ss(j,i)-cbpaq*ymu_s(2,k,i)*uaip(k,j,i)
      chit_ss(j,i)=chit_ss(j,i)-cbpaq*ymu_s(2,k,i)*uait(k,j,i)
    enddo
!         <E.B> 
    gfl_s(4,i)=gfl_s(4,i)-cbpa*ymu_s(1,k,i)*upar_s(k,2,i)
    qfl_s(4,i)=qfl_s(4,i)-cbpaq*ymu_s(2,k,i)*upar_s(k,2,i)
!         External force 
    gfl_s(5,i)=gfl_s(5,i)-cbpa*ymu_s(1,k,i)*upar_s(k,3,i)
    qfl_s(5,i)=qfl_s(5,i)-cbpaq*ymu_s(2,k,i)*upar_s(k,3,i)
  enddo
!  Pfirsch-Schluter and classical 
!       Test particle               
  cpsa=cps*(xi_s(i)/iz)
  cpsaq=cpsa*(z_j7kv*temp_i(im))
  ccla=ccl*(xi_s(i)/iz)
  cclaq=ccla*(z_j7kv*temp_i(im))
  do k=1,k_order
!         Pfirsch-Schluter 
    gfl_s(2,i)=gfl_s(2,i)-cpsa*calm_i(1,k,im)*srcth(k,i)
    qfl_s(2,i)=qfl_s(2,i)-cpsaq*calm_i(2,k,im)*srcth(k,i)
!         Classical 
    gfl_s(3,i)=gfl_s(3,i)+ccla*calm_i(1,k,im)*srcth(k,i)
    qfl_s(3,i)=qfl_s(3,i)+cclaq*calm_i(2,k,im)*srcth(k,i)
  enddo
!       Unit p'/p and T'/T 
  dp_ss(i,i)=dp_ss(i,i)-(cpsa-ccla)*calm_i(1,1,im)*srcthp(i)
  dt_ss(i,i)=dt_ss(i,i)-(cpsa-ccla)*calm_i(1,2,im)*srctht(i)
  chip_ss(i,i)=chip_ss(i,i)-(cpsaq-cclaq)*calm_i(2,1,im)*srcthp(i)
  chit_ss(i,i)=chit_ss(i,i)-(cpsaq-cclaq)*calm_i(2,2,im)*srctht(i)
!       Field particle 
  do j=1,m_s
    jm=jm_s(j)
    cpsb=cpsa*xi_s(j)
    cpsbq=cpsb*(z_j7kv*temp_i(im))
    cclb=ccla*xi_s(j)
    cclbq=cclb*(z_j7kv*temp_i(im))
    do k=1,k_order
!           Pfirsch-Schluter 
      gfl_s(2,i)=gfl_s(2,i)-cpsb*caln_ii(1,k,im,jm)*srcth(k,j)
      qfl_s(2,i)=qfl_s(2,i)-cpsbq*caln_ii(2,k,im,jm)*srcth(k,j)
!           Classical 
      gfl_s(3,i)=gfl_s(3,i)+cclb*caln_ii(1,k,im,jm)*srcth(k,j)
      qfl_s(3,i)=qfl_s(3,i)+cclbq*caln_ii(2,k,im,jm)*srcth(k,j)
    enddo
!         Unit p'/p and T'/T 
    dp_ss(j,i)=dp_ss(j,i)-(cpsb-cclb)*caln_ii(1,1,im,jm)*srcthp(j)
    dt_ss(j,i)=dt_ss(j,i)-(cpsb-cclb)*caln_ii(1,2,im,jm)*srctht(j)
    chip_ss(j,i)=chip_ss(j,i) &
              -(cpsbq-cclbq)*caln_ii(2,1,im,jm)*srcthp(j)
    chit_ss(j,i)=chit_ss(j,i) &
              -(cpsbq-cclbq)*caln_ii(2,2,im,jm)*srctht(j)
  enddo
enddo
!Electrical resistivity
p_etap=p_eb/p_ohjb
!Convert to diffusivities and conductivities
!  Full coefficient matrices 
do i=1,m_s
  im=jm_s(i)
  iza=IABS(jz_s(i))
  do j=1,m_s
    dp_ss(j,i)=dp_ss(j,i)/den_iz(im,iza)
    dt_ss(j,i)=dt_ss(j,i)/den_iz(im,iza)
    chip_ss(j,i)=chip_ss(j,i)/den_iz(im,iza)/temp_i(im)/z_j7kv
    chit_ss(j,i)=chit_ss(j,i)/den_iz(im,iza)/temp_i(im)/z_j7kv
  enddo
enddo

return
end subroutine NCLASS_FLOW

!----------------------------------------------------------------
subroutine NCLASS_K(k_banana,k_pfirsch,k_potato,m_s,jm_s,jz_s, &
                    c_potb,c_potl,p_fm,p_ft,p_ngrth,x,amu_i, &
                    temp_i,vt_i,sqz_s,ykb_s,ykp_s,ykpo_s,ykpop_s, &
                    tau_ss)
!----------------------------------------------------------------
!NCLK calculates the velocity-dependent neoclassical viscosity
!  coefficients, K
!References:                                                      
!  Hirshman, Sigmar, Nucl Fusion 21 (1981) 1079                     
!  Kessel, Nucl Fusion 34 (1994) 1221                              
!  Shaing, Yokoyama, Wakatani, Hsu, Phys Plasmas 3 (1996) 965      
!  Houlberg, Shaing, Hirshman, Zarnstorff, Phys Plasmas 4 (1997) 3230                     
!  W.A.Houlberg 6/99                          
!Input:
!  k_banana-option to include banana viscosity (-)
!          =0 off
!          =else on
!  k_pfirsch-option to include Pfirsch-Schluter viscosity (-)
!          =0 off
!          =else on
!  k_potato-option to include potato orbits (-)
!          =0 off
!          =else on
!  m_s-number of species (1<m_s<mx_ms+1)
!  jm_s(s)-isotope number of s (-)
!  jz_s(s)-charge state of s (-)
!  c_potb-kappa(0)*Bt(0)/[2*q(0)**2] (T)
!  c_potl-q(0)*R(0) (m)
!  p_fm(3)-poloidal moments of geometric factor for PS viscosity (-)
!  p_ft-trapped fraction (-)
!  p_ngrth-<n.grad(Theta)> (1/m)
!  x-velocity normalized to thermal velocity v/(2kT/m)**0.5 (-)
!  amu_i(i)-atomic mass number of i (-)
!  temp_i(i)-temperature of i (keV)
!  vt_i(i)-thermal velocity of i (m/s)
!  sqz_s(s)-orbit squeezing factor for s (-)
!Output:
!  ykb_s(s)-banana viscosity for s (kg/m**3/s)
!  ykp_s(s)-Pfirsch-Schluter viscosity for s (kg/m**3/s)
!  ykpo_s(s)-potato viscosity for s (kg/m**3/s)
!  ykpop_s(s)-potato-plateau viscosity for s (kg/m**3/s)
!  tau_ss(s1,s2)-90 degree scattering time of s1 on s2 (s)
!----------------------------------------------------------------

!Declaration of input variables
integer        k_banana,                k_pfirsch, &
            k_potato,                m_s
integer        jm_s(mx_ms),             jz_s(mx_ms)
real           c_potb,                  c_potl, &
            p_fm(3),                 p_ft, &
            p_ngrth,                 x
real           amu_i(mx_mi),            temp_i(mx_mi), &
            vt_i(mx_mi)
real           sqz_s(mx_ms)
!Declaration of output variables
real           ykb_s(mx_ms),            ykp_s(mx_ms), &
            ykpo_s(mx_ms),           ykpop_s(mx_ms)
real           tau_ss(mx_ms,mx_ms)
!Declaration of local variables
integer        i,                       im, &
            iz,                      m
real           c1,                      c2, &
            c3,                      c4, &
            z_coulomb,               z_pi, &
            z_protonmass
real           ynud_s(mx_ms),           ynut_s(mx_ms), &
            ynutis(3,mx_ms)
!Initialization
!  Physical and conversion constants
z_coulomb=1.6022e-19
z_pi=ACOS(-1.0)
z_protonmass=1.6726e-27
!  Zero out arrays
call RARRAY_ZERO(m_s,ykb_s)
call RARRAY_ZERO(m_s,ykp_s)
call RARRAY_ZERO(m_s,ykpo_s)
call RARRAY_ZERO(m_s,ykpop_s)
!  Get collisional frequencies
call NCLASS_NU(m_s,jm_s,p_ngrth,x,temp_i,vt_i,tau_ss,ynud_s, &
            ynut_s,ynutis)
!  Set velocity dependent viscosities (K's)
c1=1.5*x**2
c2=3.0*2.19/(2.0**1.5)*x**(1.0/3.0)
if(k_potato.ne.0) then
  c3=3.0*z_pi/(64.0*2.0**0.33333)/ABS(c_potl)
endif
do i=1,m_s
  im=jm_s(i)
  iz=jz_s(i)
  if(k_banana.ne.0) then
!         Provide cutoff to eliminate failure at unity trapped fraction
!         At A=>1 viscosity will go over to Pfirsch-Schluter value
    c4=1.0-p_ft
    if(c4.lt.1.0e-3) c4=1.0e-3
    ykb_s(i)=p_ft/c4/sqz_s(i)**1.5*ynud_s(i)
  endif
  if(k_pfirsch.ne.0) then
    do m=1,3
      ykp_s(i)=ykp_s(i)+c1*vt_i(im)**2*p_fm(m) &
                     *(ynutis(m,i)/ynut_s(i))
    enddo
  endif
  if(k_potato.ne.0) then
    c4=ABS(amu_i(im)*z_protonmass*vt_i(im) &
    /(iz*z_coulomb*c_potb*c_potl))
    ykpo_s(i)=c2*c4**(1.0/3.0)*ynud_s(i)/sqz_s(i)**(5.0/3.0)
    ykpop_s(i)=c3*vt_i(im)*c4**(4.0/3.0)
  endif
enddo
         
return
end subroutine NCLASS_K

!----------------------------------------------------------------
subroutine NCLASS_MN(k_order,m_i,amu_i,temp_i,capm_ii,capn_ii)
!----------------------------------------------------------------
!NCLASS_MN calculates the k_order*korder matrix of test particle (M) and
!  field particle (N) coefficients of the collision operator using the
!  Laguerre polynomials of order 3/2 as basis functions for each
!  isotopic species combination
!References:                                         
!  Hirshman, Sigmar, Nucl Fusion 21 (1981) 1079 (HS81)
!  Hirshman, Phys Fluids 20 (1977) 589
!  Houlberg, Shaing, Hirshman, Zarnstorff, Phys Plasmas 4 (1997) 3230
!  W.A.Houlberg 6/99
!Input:
!  k_order-order of v moments to be solved (-)
!       =2 u and q
!       =3 u, q, and u2
!       =else error
!  m_i-number of isotopes (1<m_i<mx_mi+1)
!  amu_i(a)-atomic mass of a (-)
!  temp_i(a)-temperature of a (keV)
!Output:
!  capm_ii(3,3,i1,i2)-test part (tp) friction matrix for i1 on i2 (-)
!  capn_ii(3,3,i1,i2)-field part (fp) friction matrix for i1 on i2 (-)
!Comments:                                                     
!The indices on the M and N matrices are one greater than the notation
!  in the review article so as to avoid 0 as an index
!----------------------------------------------------------------

!Declaration of input variables
integer        k_order,                 m_i
real           amu_i(mx_mi),            temp_i(mx_mi)
!Declaration of output variables
real           capm_ii(3,3,mx_mi,mx_mi), capn_ii(3,3,mx_mi,mx_mi)
!Declaration of local variables
integer        im,                      jm
real           xab,                     xab2,                     &
            xmab,                    xtab, &
            yab32,                   yab52, &
            yab72,                   yab92
!Loop over isotope a
do im=1,m_i
!Loop over isotope b
  do jm=1,m_i
!         Ratio of masses
    xmab=amu_i(im)/amu_i(jm)
!         Ratio of temperatures
    xtab=temp_i(im)/temp_i(jm)
!         Ratio of thermal velocities, vtb/vta
    xab=SQRT(xmab/xtab)
!  Elements of M
    xab2=xab**2
    yab32=(1.0+xab2)*SQRT(1.0+xab2)
    yab52=(1.0+xab2)*yab32
    if(k_order.eq.3) then
      yab72=(1.0+xab2)*yab52
      yab92=(1.0+xab2)*yab72
    endif
!         Eqn 4.11 for M00 (HS81)
    capm_ii(1,1,im,jm)=-(1.0+xmab)/yab32
!         Eqn 4.12 for M01 (HS81)
    capm_ii(1,2,im,jm)=3.0/2.0*(1.0+xmab)/yab52
!         Eqn 4.8 for M10 (HS81)
    capm_ii(2,1,im,jm)=capm_ii(1,2,im,jm)
!         Eqn 4.13 for M11 (HS81)
    capm_ii(2,2,im,jm)=-(13.0/4.0+xab2*(4.0+xab2*15.0/2.0))/yab52
    if(k_order.eq.3) then
!           Eqn 4.15 for M02 (HS81)
      capm_ii(1,3,im,jm)=-15.0/8.0*(1.0+xmab)/yab72          
!           Eqn 4.16 for M12 (HS81)
      capm_ii(2,3,im,jm)=(69.0/16.0+xab2*(6.0+xab2*63.0/4.0)) &
                      /yab72
!           Eqn 4.8 for M20 (HS81)
      capm_ii(3,1,im,jm)=capm_ii(1,3,im,jm)
!           Eqn 4.8 for M21 (HS81)
      capm_ii(3,2,im,jm)=capm_ii(2,3,im,jm)
!           Eqn 5.21 for M22 (HS81)
      capm_ii(3,3,im,jm)=-(433.0/64.0+xab2*(17.0+xab2*(459.0/8.0 &
                      +xab2*(28.0+xab2*175.0/8.0))))/yab92
    endif     
!  Elements of N
!       Momentum conservation, Eqn 4.11 for N00 (HS81)
    capn_ii(1,1,im,jm)=-capm_ii(1,1,im,jm)
!       Eqn 4.9 and 4.12 for N01 (HS81)
    capn_ii(1,2,im,jm)=-xab2*capm_ii(1,2,im,jm)
!       Momentum conservation, Eqn 4.12 for N10 (HS81)
    capn_ii(2,1,im,jm)=-capm_ii(2,1,im,jm)
!       Eqn 4.14 for N11 (HS81) - corrected rhs
    capn_ii(2,2,im,jm)=(27.0/4.0)*SQRT(xtab)*xab2/yab52
    if(k_order.eq.3) then
!       Eqn 4.15 for N02 (HS81) - corrected rhs by Ta/Tb
      capn_ii(1,3,im,jm)=-xab2**2*capm_ii(1,3,im,jm)
!       Eqn 4.17 for N12 (HS81)
      capn_ii(2,3,im,jm)=-225.0/16.0*xtab*xab2**2/yab72
!       Momentum conservation for N20 (HS81)
      capn_ii(3,1,im,jm)=-capm_ii(3,1,im,jm)
!       Eqn 4.9 and 4.17 for N21 (HS81)
      capn_ii(3,2,im,jm)=-225.0/16.0*xab2**2/yab72
!       Eqn 5.22 for N22 (HS81) 
      capn_ii(3,3,im,jm)=2625.0/64.0*xtab*xab2**2/yab92
    endif
  enddo   
enddo

return
end subroutine NCLASS_MN

!----------------------------------------------------------------
subroutine NCLASS_MU(k_order,k_banana,k_pfirsch,k_potato,m_s,jm_s, &
                  jz_s,c_potb,c_potl,p_fm,p_ft,p_ngrth,amu_i, &
                  temp_i,vt_i,den_iz,sqz_s,ymu_s,tau_ss)
!----------------------------------------------------------------
!NCLASS_MU calculates the k_order*k_order matrix of neoclassical
!  viscosities by integrating the velocity-dependent banana and Pfirsch-
!  Schluter contributions
!References:                                             
!  Shaing, Yokoyama, Wakatani, Hsu, Phys Plasmas 3 (1996) 965
!  Houlberg, Shaing, Hirshman, Zarnstorff, Phys Plasmas 4 (1997) 3230
!  W.A.Houlberg 6/99
!Input:
!  k_order-order of v moments to be solved (-)
!       =2 u and q
!       =3 u, q, and u2
!       =else error
!  k_banana-option to include banana viscosity (-)
!          =0 off
!          =else on
!  k_pfirsch-option to include Pfirsch-Schluter viscosity (-)
!          =0 off
!          =else on
!  k_potato-option to include potato orbits (-)
!          =0 off
!          =else on
!  m_s-number of species (1<m_s<mx_ms+1)
!  jm_s(s)-isotope number of species s (-)
!  jz_s(s)-charge state of s (-)
!  c_potb-kappa(0)*Bt(0)/[2*q(0)**2] (T)
!  c_potl-q(0)*R(0) (m)
!  p_fm(3)-poloidal moments of geometric factor for PS viscosity (-)
!  p_ft-trapped fraction (-)
!  p_ngrth-<n.grad(Theta)> (1/m)
!  amu_i(i)-atomic mass number of i (-)
!  temp_i(i)-temperature of i (keV)
!  vt_i(i)-thermal velocity of i (m/s)
!  den_iz(i,z)-density of i,z (/m**3)
!  sqz_s(s)-orbit squeezing factor for s (-)
!Output:
!  ymu_s(s)-normalized viscosity for s (kg/m**3/s)
!  tau_ss(s1,s2)-90 degree scattering time of s1 on s2 (s)
!----------------------------------------------------------------

!Declaration of input variables
integer        k_banana,                k_order, &
            k_pfirsch,               k_potato, &
            m_s
integer        jm_s(mx_ms),             jz_s(mx_ms)
real           c_potb,                  c_potl, &
            p_ft,                    p_ngrth, &
            p_fm(3)
real           amu_i(mx_mi),            temp_i(mx_mi), &
            vt_i(mx_mi)
real           den_iz(mx_mi,mx_mz)
real           sqz_s(mx_ms)
!Declaration of output variables
real           ymu_s(3,3,mx_ms)
real           tau_ss(mx_ms,mx_ms)
!Declaration of local variables
integer        mpnts
PARAMETER     (mpnts=13)
integer        i,                       init, &
            im,                      iza, &
            k,                       l, &
            m
real           bmax,                    c1, &
            c2,                      ewt, &
            expmx2,                  h, &
            x2,                      x4, &
            x6,                      x8, &
            x10,                     x12, &
            xk,                      xx, &
            z_pi,                    z_protonmass
real           x(mpnts),                w(mpnts,5)
real           ykb_s(mx_ms),            ykp_s(mx_ms), &
            ykpo_s(mx_ms),           ykpop_s(mx_ms), &
            ymubs(3,3,mx_ms),        ymubps(3,3,mx_ms), &
            ymupps(3,3,mx_ms),       ymupos(3,3,mx_ms)
real           dum(3)
SAVE           init
SAVE           c1,                      h, &
            x,                       w, &
            z_protonmass
DATA bmax/     3.2/
DATA init/     0/
!Initialization
if(init.eq.0) then
!  Physical and conversion constants
  z_pi=ACOS(-1.0)
  z_protonmass=1.6726e-27
!  Set integration points and weights
  h=bmax/(mpnts-1)
  x(1)=0.0
  w(1,1)=0.0
  w(1,2)=0.0
  w(1,3)=0.0
  w(1,4)=0.0
  w(1,5)=0.0
  do m=2,mpnts
    x(m)=h*(m-1)
    x2=x(m)*x(m)
    expmx2=EXP(-x2)
    x4=x2*x2
    w(m,1)=x4*expmx2
    x6=x4*x2
    w(m,2)=x6*expmx2
    x8=x4*x4
    w(m,3)=x8*expmx2
    x10=x4*x6
    w(m,4)=x10*expmx2
    x12=x6*x6
    w(m,5)=x12*expmx2
  enddo   
  c1=8.0/3.0/SQRT(z_pi)*h
  init=1
endif
call RARRAY_ZERO(9*m_s,ymu_s)
call RARRAY_ZERO(9*m_s,ymubs)
call RARRAY_ZERO(9*m_s,ymubps)
call RARRAY_ZERO(9*m_s,ymupos)
call RARRAY_ZERO(9*m_s,ymupps)
if((k_banana.ne.0).or.(k_pfirsch.ne.0)) then
!Loop over grid (first node has null value)
  do m=2,mpnts
    if(m.eq.mpnts) then
!           Use half weight for end point
      ewt=0.5
    else
!           Use full weight
      ewt=1.0
    endif
    xx=x(m)
!  Get velocity-dependent k values        
    call NCLASS_K(k_banana,k_pfirsch,k_potato,m_s,jm_s,jz_s, &
               c_potb,c_potl,p_fm,p_ft,p_ngrth,xx,amu_i,temp_i, &
               vt_i,sqz_s,ykb_s,ykp_s,ykpo_s,ykpop_s,tau_ss)
!  Loop over species
    do i=1,m_s
      im=jm_s(i)
      iza=IABS(jz_s(i))
      c2=c1*ewt*den_iz(im,iza)*amu_i(im)*z_protonmass
      dum(1)=c2*w(m,1)
      dum(2)=c2*(w(m,2)-5.0/2.0*w(m,1))
      dum(3)=c2*(w(m,3)-5.0*w(m,2)+(25.0/4.0)*w(m,1))
      if(k_banana.eq.0) then
        xk=ykp_s(i)
      elseif(k_pfirsch.eq.0) then
        xk=ykb_s(i)
      else
        xk=ykb_s(i)*ykp_s(i)/(ykb_s(i)+ykp_s(i))
      endif
      ymubs(1,1,i)=ymubs(1,1,i)+ykb_s(i)*dum(1)
      ymubs(1,2,i)=ymubs(1,2,i)+ykb_s(i)*dum(2)
      ymubs(2,2,i)=ymubs(2,2,i)+ykb_s(i)*dum(3)
      ymubps(1,1,i)=ymubps(1,1,i)+xk*dum(1)
      ymubps(1,2,i)=ymubps(1,2,i)+xk*dum(2)
      ymubps(2,2,i)=ymubps(2,2,i)+xk*dum(3)
      if(k_potato.ne.0) then
        ymupos(1,1,i)=ymupos(1,1,i)+ykpo_s(i)*dum(1)
        ymupos(1,2,i)=ymupos(1,2,i)+ykpo_s(i)*dum(2)
        ymupos(2,2,i)=ymupos(2,2,i)+ykpo_s(i)*dum(3)
        xk=ykpo_s(i)*ykpop_s(i)/(ykpo_s(i)+ykpop_s(i))
        ymupps(1,1,i)=ymupps(1,1,i)+xk*dum(1)
        ymupps(1,2,i)=ymupps(1,2,i)+xk*dum(2)
        ymupps(2,2,i)=ymupps(2,2,i)+xk*dum(3)
      endif
      if(k_order.eq.3) then
        dum(1)=c2*((1.0/2.0)*w(m,3)-(7.0/2.0)*w(m,2) &
            +(35.0/8.0)*w(m,1))
        dum(2)=c2*((1.0/2.0)*w(m,4)-(19.0/4.0)*w(m,3) &
            +(105.0/8.0)*w(m,2)-(175.0/16.0)*w(m,1))
        dum(3)=c2*((1.0/4.0)*w(m,5)-(7.0/2.0)*w(m,4) &
            +(133.0/8.0)*w(m,3)-(245.0/8.0)*w(m,2) &
            +(1225.0/64.0)*w(m,1))
        ymubs(1,3,i)=ymubs(1,3,i)+ykb_s(i)*dum(1)
        ymubs(2,3,i)=ymubs(2,3,i)+ykb_s(i)*dum(2)
        ymubs(3,3,i)=ymubs(3,3,i)+ykb_s(i)*dum(3)
        xk=ykb_s(i)*ykp_s(i)/(ykb_s(i)+ykp_s(i))
        ymubps(1,3,i)=ymubps(1,3,i)+xk*dum(1)
        ymubps(2,3,i)=ymubps(2,3,i)+xk*dum(2)
        ymubps(3,3,i)=ymubps(3,3,i)+xk*dum(3)
        if(k_potato.ne.0) then
          ymupos(1,3,i)=ymupos(1,3,i)+ykpo_s(i)*dum(1)
          ymupos(2,3,i)=ymupos(2,3,i)+ykpo_s(i)*dum(2)
          ymupos(3,3,i)=ymupos(3,3,i)+ykpo_s(i)*dum(3)
          xk=ykpo_s(i)*ykpop_s(i)/(ykpo_s(i)+ykpop_s(i))
          ymupps(1,3,i)=ymupps(1,3,i)+xk*dum(1)
          ymupps(2,3,i)=ymupps(2,3,i)+xk*dum(2)
          ymupps(3,3,i)=ymupps(3,3,i)+xk*dum(3)
        endif
      endif
    enddo    
  enddo
!Load net viscosity
  do i=1,m_s
    do l=1,k_order
      do k=1,l
        if(k_potato.eq.0) then
!               Banana Pfirsch-Schluter
          ymu_s(k,l,i)=ymubps(k,l,i)
        else
!               Banana Pfirsch-Schluter plus potato potato-plateau
!orig                ymu_s(k,l,i)=(ymupos(k,l,i)**3*ymupps(k,l,i)
!orig     #                       +ymubs(k,l,i)**3*ymubps(k,l,i))
!orig     #                       /(ymupos(k,l,i)**3+ymubs(k,l,i)**3)
! The 3 above line have been replaced to avoid division by zero
!     when ymupos^3 and ymubs^3 become very small
!     G.Pereverzev 16-10-2008
          ymu_s(k,l,i) = (ymupos(k,l,i)/ymubs(k,l,i))**3
          ymu_s(k,l,i)=(ymu_s(k,l,i)*ymupps(k,l,i)+ymubps(k,l,i)) &
                    /(1.+ymu_s(k,l,i))
        endif
      enddo
    enddo
  enddo
!Fill viscosity matrix using symmetry
  do i=1,m_s
    do l=1,k_order-1
      do k=l+1,k_order
          ymu_s(k,l,i)=ymu_s(l,k,i)
      enddo
    enddo
  enddo
endif

return
end subroutine NCLASS_MU

!----------------------------------------------------------------
subroutine NCLASS_NU(m_s,jm_s,p_ngrth,x,temp_i,vt_i,tau_ss,ynud_s, &
                  ynut_s,ynutis)
!----------------------------------------------------------------
!NCLASS_NU calculates the velocity dependent pitch angle diffusion and
!  anisotropy relaxation rates, nu_D, nu_T, and nu_T*I_Rm
!References:                                                     
!  Hirshman, Sigmar, Phys Fluids 19 (1976) 1532
!  Hirshman, Sigmar, Nucl Fusion 21 (1981) 1079
!  Shaing, Yokoyama, Wakatani, Hsu, Phys Plasmas 3 (1996) 965
!  Houlberg, Shaing, Hirshman, Zarnstorff, Phys Plasmas 4 (1997) 3230
!  W.A.Houlberg 6/99
!Input:
!  m_s-number of species (1<m_s<mx_ms+1)
!  jm_s(s)-isotope number of s (-)
!  p_ngrth-<n.grad(Theta)> (1/m)
!  x-velocity normalized to thermal velocity v/(2kT/m)**0.5 (-)
!  temp_i(i)-temperature of i (keV)
!  vt_i(i)-thermal velocity of i (m/s)
!  tau_ss(s1,s2)-90 degree scattering time of s1 on s2 (s)
!Output:
!  ynud_s(s)-pitch angle diffusion rate for s (/s)
!  ynut_s(s)-anisotropy relaxation rate for s (/s)
!  ynutis(3,s)-PS anisotropy relaxation rates for s (/s)
!----------------------------------------------------------------

integer, parameter :: mx_mi=9, mx_ms=40

!Declaration of input variables
integer        m_s
integer        jm_s(mx_ms)       
real           p_ngrth,                  x
REAl           temp_i(mx_mi),            vt_i(mx_mi)
real           tau_ss(mx_ms,mx_ms)
!Declaration of output variables
real           ynud_s(mx_ms),            ynut_s(mx_ms), &
            ynutis(3,m_s)
!Declaration of local variables
integer        i,                       im, &
            j,                       jm, &
            m
real           c1,                      c2, &
            c3,                      g, &
            phi,                     z_pi
!Declaration of external functions

!Initializaton
!  Physical and conversion constants
z_pi=ACOS(-1.0)
!  Zero out arrays
call RARRAY_ZERO(m_s,ynud_s)
call RARRAY_ZERO(m_s,ynut_s)
call RARRAY_ZERO(3*m_s,ynutis)
do i=1,m_s
  im=jm_s(i)
!Calculate nu_D and nu_T
  do j=1,m_s
    jm=jm_s(j)
    c1=vt_i(jm)/vt_i(im)
    c2=x/c1
    phi=U_ERF(c2)
    g=(phi-c2*(2.0/SQRT(z_pi))*EXP(-c2**2))/(2.0*c2**2)
    ynud_s(i)=ynud_s(i)+(3.0*SQRT(z_pi)/4.0)*(phi-g)/x**3 &
                     /tau_ss(i,j)
    ynut_s(i)=ynut_s(i)+((3.0*SQRT(z_pi)/4.0)*((phi-3.0*g)/x**3 &
                    +4.0*(temp_i(im)/temp_i(jm) &
                    +1.0/c1**2)*g/x))/tau_ss(i,j)
  enddo   
! Calculate nu_T*I_m
  do m=1,3
    if(ABS(p_ngrth).gt.0.0) then
      c1=x*vt_i(im)*m*p_ngrth
      c2=(ynut_s(i)/c1)**2
      if(c2.gt.9.0) then
        ynutis(m,i)=0.4
      else
        c3=ynut_s(i)/c1*ATAN(c1/ynut_s(i))
        ynutis(m,i)=0.5*c3+c2*(3.0*(c3-0.5)+c2*4.5*(c3-1.0))
      endif
    else
      ynutis(m,i)=0.4
    endif
  enddo   
enddo

return
end subroutine NCLASS_NU

!----------------------------------------------------------------
subroutine NCLASS_TAU(m_i,m_s,jm_s,jz_s,amu_i,temp_i,vt_i,den_iz, &
                   amnt_ii,tau_ss)
!----------------------------------------------------------------
!NCLASS_TAU calculates the collision time for 90 degree scattering
!  assuming a common value for the Coulomb logarithm for each isotope
!References:                                      
!  Hirshman, Sigmar, Nucl Fusion 21 (1981) 1079
!  Houlberg, Shaing, Hirshman, Zarnstorff, Phys Plasmas 4 (1997) 3230
!  W.A.Houlberg 6/99
!Input:
!  m_i-number of isotopes (1<m_i<mx_mi+1)
!  m_s-number of (1<m_s<mx_ms+1)
!  jm_s(s)-isotope number of s (-)
!  jz_s(s)-charge state of s (-)
!  amu_i(i)-atomic mass number of i (-)
!  temp_i(i)-temperature of i (keV)
!  vt_i(i)-thermal velocity of i (m/s)
!  den_iz(i,z)-density of i,z (/m**3)
!Output:
!  amnt_ii(s1,s2)-eff relaxation rate for s1 on s2 (kg/m**3/s)
!  tau_ss(s1,s2)-90 degree scattering time of s1 on s2 (s)
!----------------------------------------------------------------

!Declaration of input variables
integer        m_i,                      m_s
integer        jm_s(mx_ms),              jz_s(mx_ms)
real           amu_i(mx_mi),             temp_i(mx_mi), &
            vt_i(mx_mi)
real           den_iz(mx_mi,mx_mz)
!Declaration of output variables
real           amnt_ii(mx_mi,mx_mi)
real           tau_ss(mx_ms,mx_ms)
!Declaration of local variables
integer        i,                       im, &
            iz,                      iza, &
            j,                       jm, &
            jz,                      jza
real           c1,                      c2, &
            clnab,                   z_coulomb, &
            z_epsilon0,              z_pi, &
            z_protonmass
real           xlnab(mx_mi,mx_mi),      xn(mx_mi), &
            xnz(mx_mi),              xz(mx_mi), &
            xnz2(mx_mi)
!Initialization
!  Physical and conversion constants
z_coulomb=1.6022e-19
z_epsilon0=8.8542e-12
z_pi=ACOS(-1.0)
z_protonmass=1.6726e-27
!  Zero out arrays
call RARRAY_ZERO(mx_mi*m_i,amnt_ii)
call RARRAY_ZERO(mx_ms*m_s,tau_ss)
call RARRAY_ZERO(m_i,xn)
call RARRAY_ZERO(m_i,xnz)
call RARRAY_ZERO(m_i,xnz2)
c1=4.0/3.0/SQRT(z_pi)*4.0*z_pi*(z_coulomb &
/(4.0*z_pi*z_epsilon0))**2*(z_coulomb/z_protonmass)**2
!Coulomb logarithm
do i=1,m_s
  im=jm_s(i)
  iz=jz_s(i)
  iza=IABS(iz)
  if(iz.lt.0) then
!         Electrons
    clnab=37.8-ALOG(SQRT(den_iz(im,iza))/temp_i(im))
  endif
  xn(im)=xn(im)+den_iz(im,iza)
  xnz(im)=xnz(im)+iz*den_iz(im,iza)
  xnz2(im)=xnz2(im)+iz**2*den_iz(im,iza)
enddo
do i=1,m_i
  xz(i)=xnz(i)/xn(i)
enddo
do im=1,m_i
  iz=jz_s(im)
  do jm=1,m_i
    jz=jz_s(jm)
    if(iz.lt.0.or.jz.lt.0) then
!  Electrons
      xlnab(im,jm)=clnab
    else
!  Ions
      xlnab(im,jm)=40.3 &
                -ALOG(xz(im)*xz(jm)*(amu_i(im)+amu_i(jm)) &
                /(amu_i(im)*temp_i(jm)+amu_i(jm)*temp_i(im)) &
                *SQRT(xnz2(im)/temp_i(im)+xnz2(jm)/temp_i(jm)))
    endif
  enddo
enddo
!Collision times and mass density weighted collision rates
do i=1,m_s
  im=jm_s(i)
  iz=jz_s(i)
  iza=IABS(iz)
  c2=(vt_i(im)**3)*amu_i(im)**2/c1
  do j=1,m_s
    jm=jm_s(j)
    jz=jz_s(j)
    jza=IABS(jz)
    tau_ss(i,j)=c2/xlnab(im,jm)/iz**2 &
             /(den_iz(jm,jza)*jz**2)
    amnt_ii(im,jm)=amnt_ii(im,jm)+amu_i(im)*z_protonmass &
                *den_iz(im,iza)/tau_ss(i,j)
  enddo          
enddo

return
end subroutine NCLASS_TAU

!----------------------------------------------------------------
subroutine WRITE_C(nout,n_c,c,n_l)
!----------------------------------------------------------------
!WRITE_C writes out n_c character variables in fields of length n_l
!Input:
!  nout-output file unit number (-)
!  n_c-number of character variables
!  c()-character array
!  n_l-length of field
!  W.A. Houlberg 2/99
!----------------------------------------------------------------

!Declaration of input variables
character*(*)  c(*)
integer        n_c,                     n_l, &
            nout
!Declaration of local variables
character*30   char
character*2    c_c,                     c_l
integer        i
WRITE(c_c,'(i2)') n_c
WRITE(c_l,'(i2)') n_l
char='('//c_c//'a'//c_l//')'
WRITE(nout,char) (c(i),i=1,n_c)

return
end subroutine WRITE_C

!----------------------------------------------------------------
subroutine WRITE_IR(nout,n_i,i,n_r,r,k_format)
!----------------------------------------------------------------
!WRITE_IR writes a writes n_i integer variables followed by n_r real
!  variables to the unit nout
!Input:
!  nout-output file unit number (-)
!  n_i-number of integer variables (-)
!  i()-integer array
!  n_r-number of real variables (-)
!  r()-real array
!  k_format-real format option
!          =1 use f
!          =2 use 1pe
!          =else use e
!  W.A. Houlberg 2/99
!----------------------------------------------------------------

!Declaration of input variables
integer        k_format,                n_i, &
            n_r,                     nout
integer        i(*)
real           r(*)
!Declaration of local variables
character*30   char
character*2    c_i,                     c_r
integer        j
WRITE(c_i,'(i2)') n_i
WRITE(c_r,'(i2)') n_r
if(n_i.eq.0) then
!Real data only
  if(k_format.eq.1) then
    char='('//c_r//'(f12.6))'
  elseif(k_format.eq.2) then
    char='('//c_r//'(1pe12.4))'
  else
    char='('//c_r//'(e12.4))'
  endif
  WRITE(nout,char) (r(j),j=1,n_r)
elseif(n_r.eq.0) then
!Integer data only
  char='('//c_i//'(i12))'
  WRITE(nout,char) (i(j),j=1,n_i)
else
!Both integer and real data
  if(k_format.eq.1) then
    char='('//c_i//'(i12),'//c_r//'(f12.6))'
  elseif(k_format.eq.2) then
    char='('//c_i//'(i12),'//c_r//'(1pe12.4))'
  else
    char='('//c_i//'(i12),'//c_r//'(e12.4))'
  endif
  WRITE(nout,char) (i(j),j=1,n_i),(r(j),j=1,n_r)
endif

return
end subroutine WRITE_IR

!----------------------------------------------------------------
subroutine WRITE_LINE(nout,label,nabove,nbelow)
!----------------------------------------------------------------
!WRITE_LINE writes a label preceeded by nabove blank lines and
!  followed by nbelow blank lines to the unit nout
!Input:
!  nout-output file unit number (-)
!  label-label to be printed
!  nabove-number of blank lines above label
!  nbelow-number of blanklines below label
!  W.A. Houlberg 2/99
!----------------------------------------------------------------

!Declaration of input variables
character*(*)  label
integer        nabove,                  nbelow, &
            nout
!Declaration of local variables
integer        j
if(nabove.gt.0) then
  do j=1,nabove
    WRITE(nout,'( )')
  enddo
endif
WRITE(nout,'(a)') label
if(nbelow.gt.0) then
  do j=1,nbelow
    WRITE(nout,'( )')
  enddo
endif

return
end subroutine WRITE_LINE

!----------------------------------------------------------------
subroutine WRITE_LINE_IR(nout,label,n_i,i,n_r,r,k_format)
!----------------------------------------------------------------
!WRITE_LINE_IR writes a writes a label followed by n_i integer numbers
!  and n_r real numbers to the unit nout
!Input:
!  nout-output file unit number
!  label-label to be printed
!  n_i-number of integer variables
!  i()-integer array
!  n_r-number of real variables
!  r()-real array
!  k_format-format option for real variables
!          =1 use f
!          =2 use 1pe
!          =else use e
!  W.A. Houlberg 12/98
!----------------------------------------------------------------

!Declaration of input variables
character*(*)  label
integer        k_format,                n_i, &
            n_r,                     nout
integer        i(*)
real           r(*)
!Declaration of local variables
character*30   char
character*2    c_i,                     c_r
integer        j
WRITE(c_i,'(i2)') n_i
WRITE(c_r,'(i2)') n_r
if(n_i.eq.0) then
!Real data only
  if(k_format.eq.1) then
    char='(a48,'//c_r//'(f12.6))'
  elseif(k_format.eq.2) then
    char='(a48,'//c_r//'(1pe12.4))'
  else
    char='(a48,'//c_r//'(e12.4))'
  endif
  WRITE(nout,char) label(1:48),(r(j),j=1,n_r)
elseif(n_r.eq.0) then
!Integer data only
  char='(a48,'//c_i//'(i12))'
  WRITE(nout,char) label(1:48),(i(j),j=1,n_i)
else
!Both integer and real data
  if(k_format.eq.1) then
    char='(a48,'//c_i//'(i12),'//c_r//'(f12.6))'
  elseif(k_format.eq.2) then
    char='(a48,'//c_i//'(i12),'//c_r//'(1pe12.4))'
  else
    char='(a48,'//c_i//'(i12),'//c_r//'(e12.4))'
  endif
  WRITE(nout,char) label(1:48),(i(j),j=1,n_i),(r(j),j=1,n_r)
endif

return
end subroutine WRITE_LINE_IR

!----------------------------------------------------------------
subroutine U_LU_BACKSUB(a,n,ndim,indx,b)
!----------------------------------------------------------------
!U_LU_BACKSUB solves the matrix equation a x = b, where a is in LU form as
!  generated by a prior call to U_LU_DECOMP
!References:
!  Flannery, Teukolsky, Vetterling, Numerical Recipes
!  W.A.Houlberg 1/99
!Input:
!  a-matrix in lu decomposed form
!  n-dimension of matrix
!  ndim-first dimension of a array
!  indx-vector showing row permutations due to partial pivoting
!  b-right hand side of equation (overwritten on return)
!Output:
!  b-solution vector
!Comments: 
!  The complete procedure is, given the equation a*x = b
!    call U_LU_DECOMP(a,n,ndim,indx,d,iflag)
!    if(iflag.eq.0) call U_LU_BACKSUB(a,n,ndim,indx,b) 
!  and b is now the solution vector
!  To solve with a different right hand side, just reload b as  desired
!  and use the same LU decomposition, call U_LU_BACKSUB(a,n,ndim,indx,b)                                 *
!----------------------------------------------------------------

!Declaration of input variables
integer        indx(*),               n, &
            ndim
real           a(ndim,*),             b(*)
!Declaration of local variables
integer        i,                     ii, &
            j,                     k
real           sum
!Find the index of the first nonzero element of b
ii=0
do i=1,n
  k=indx(i)
  sum=b(k)
  b(k)=b(i)
  if(ii.ne.0) then
    do j=ii,i-1
      sum=sum-a(i,j)*b(j)
    enddo   
  elseif (sum.ne.0.0) then
    ii=i
  endif
  b(i)=sum
enddo    
!  Begin back substitution
do i=n,1,-1
  sum=b(i)
  if(i.lt.n) then
    do j=i+1,n
      sum=sum-a(i,j)*b(j)
    enddo   
  endif
  b(i)=sum/a(i,i)
enddo   
return
end subroutine U_LU_BACKSUB

!----------------------------------------------------------------
subroutine U_LU_DECOMP(a,n,ndim,indx,d,iflag)
!----------------------------------------------------------------
!U_LU_DECOMP performs an LU decomposition of the matrix a ans is called
!   prior to U_LU_BACKSUB to solve linear equations or to invert a
!   matrix
!References:
!  Flannery, Teukolsky, Vetterling, Numerical Recipes
!  W.A.Houlberg 1/99
!Input:
!  a(ndim,ndim)-square matrix, overwritten on return
!  n-number of equations to be solved
!  ndim-first dimension of a array
!Output:
!  indx-vector showing row permutations due to partial pivoting
!  d-flag for number of row exchanges
!   =1.0 even number
!   =-1.0 odd number
!  iflag-error flag.
!       =0 no errors
!       =1 singular matrix
!----------------------------------------------------------------

integer         nmax
PARAMETER      (nmax=100)
!Declaration of input
integer         iflag,                indx(*), &
             n,                    ndim
real            a(ndim,*),            d
!Declaration of local variables
integer         i,                    imax, &
             j,                    k
real            aamax,                dum, &
             sum,                  vv(nmax)
!Initialization
iflag=0
d=1.0
!Loop over rows to get the implicit scaling information
do i=1,n
  aamax=0.0
  do j=1,n
    if(ABS(a(i,j)).gt.aamax) aamax=ABS(a(i,j))
  enddo    
  if(aamax.eq.0.0) then
    iflag=1
    return
  endif
  vv(i)=1.0/aamax
enddo   
!Loop over columns using Crout's method
do j=1,n
  do i=1,j-1
    sum=a(i,j)
    do k=1,i-1
      sum=sum-a(i,k)*a(k,j)
    enddo   
    a(i,j)=sum
  enddo    
!  Search for largest pivot element using dum as a figure of merit
  aamax=0.0
  do i=j,n
    sum=a(i,j)
    do k=1,j-1
      sum=sum-a(i,k)*a(k,j)
    enddo    
    a(i,j)=sum
    dum=vv(i)*ABS(sum)
    if(dum.ge.aamax) then
      imax=i
      aamax=dum
    endif
  enddo   
  if(j.ne.imax) then
!         Interchange rows
    do k=1,n
      dum=a(imax,k)
      a(imax,k)=a(j,k)
      a(j,k)=dum
    enddo    
    d=-d
    vv(imax)=vv(j)
  endif
  indx(j)=imax
  if(a(j,j).eq.0.0) then
    iflag=1
    return
  endif
  if(j.ne.n) then
!         Divide by pivot element
    dum=1.0/a(j,j)
    do i=j+1,n
      a(i,j)=a(i,j)*dum
    enddo   
  endif
enddo   

return
end subroutine U_LU_DECOMP

!----------------------------------------------------------------
real FUNCTION U_ERF(x)
!----------------------------------------------------------------
!U_ERF evaluates the error function
!  W.A.Houlberg 1/99
!Input:
!  x-argument of error function

real           x

integer        i
real           a,                       b, &
            c,                       d, &
            del,                     gln, &
            h,                       x2
gln=5.723649e-1
x2=x**2
U_ERF=0.0
if(x2.lt.1.5) then
  a=0.5
  b=2.0
  del=b
  do i=1,100
    a=a+1.0
    del=del*x2/a
    b=b+del
    if(ABS(del).lt.ABS(b)*3.0e-7) then
      U_ERF=b*EXP(-x2+0.5*LOG(x2)-gln)
      if(x.lt.0.0) U_ERF=-U_ERF
      return
    endif
  enddo
else
  b=x2+0.5
  c=1.0/1.0e-30
  d=1.0/b
  h=d
  do i=1,100
    a=-i*(i-0.5)
    b=b+2.0
    d=a*d+b
    if(ABS(d).lt.1.0e-30) d=1.0e-30
    c=b+a/c
    if(ABS(c).lt.1.0e-30) c=1.0e-30
    d=1.0/d
    del=d*c
    h=h*del
    if(ABS(del-1.0).lt.3.0e-7) then
      U_ERF=1.0-EXP(-x2+0.5*LOG(x2)-gln)*h
      if(x.lt.0.0) U_ERF=-U_ERF
      return
    endif
  enddo
endif

return
end function U_ERF

!----------------------------------------------------------------
subroutine RARRAY_ZERO(n,x)
!----------------------------------------------------------------
!RARRAY_ZERO sets the elements of array x to 0.0
!  W.A.Houlberg 12/98
!Input:
!  n-number of elements to be zeroed
!  x-array to be zeroed
!Output:
!  x-zeroed array
!----------------------------------------------------------------

integer        n
real           x(*)

integer        i
do i=1,n
  x(i)=0.0
enddo   

return
end subroutine RARRAY_ZERO

!----------------------------------------------------------------
real function RARRAY_SUM(n,x,incx)
!----------------------------------------------------------------
!  RARRAY_SUM is the sum of the elements of the array x
!  W.A.Houlberg 12/98
!  Input:
!  n-number of elements to be summed
!  x-array to be summed
!  incx-increment in sx index
!----------------------------------------------------------------

!Declaration of input variables
integer        incx,                    n
real           x(*)
!Declaration of local variables
integer        i,                       ix
RARRAY_SUM=0.0
ix=1
do i=1,n
  RARRAY_SUM=RARRAY_SUM+x(ix)
  ix=ix+incx
enddo   
return
end function RARRAY_SUM

!----------------------------------------------------------------
subroutine RARRAY_COPY(n,x,incx,y,incy)
!----------------------------------------------------------------
!RARRAY_COPY copies elements of the array x into y
!  W.A.Houlberg 12/98
!Input:  
!  n-number of elements to be copied
!  x-array to be copied
!  incx-increment in x index
!  incy-increment in y index
!Output:
!  y-new array
!----------------------------------------------------------------

!Declaration of input variables
integer        incx,                    incy, &
            n
real           x(*)
!Declaration of output variables
real           y(*)
!Declaration of local variables
integer        i,                       ix, &
            iy
ix=1
iy=1
do i=1,n
  y(iy)=x(ix)
  ix=ix+incx
  iy=iy+incy
enddo   

return
end subroutine RARRAY_COPY

end module nclass_mod
