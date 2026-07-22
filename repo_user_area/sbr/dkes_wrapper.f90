SUBROUTINE DKES_INTERFACE_WRAPPER(er_field,no_edge)
!E. Buglione-Ceresa, Master Thesis, 2025 TUM Garching & MPG-IPP Garching
! only used ASTRA variable of interest is er_field = ER/1000. (kV/m)
! E. Fable added er_field as variable

  use status  !dkes astra variables are defined here as (NRD, 50)
  use scalars
  use a2dkes

  IMPLICIT NONE

  double precision er_field(NRD),no_edge
  integer i,i_option,i_edge

  dkes2astra_variables=0.d0
  
  CALL CALC_GRADIENT(NA1, HRO, NE, dkes2astra_variables(:,1))
  CALL CALC_GRADIENT(NA1, HRO, TE, dkes2astra_variables(:,2))
  CALL CALC_GRADIENT(NA1, HRO, TI, dkes2astra_variables(:,3))
  CALL CALC_GRADIENT(NA1, HRO, NI, dkes2astra_variables(:,4))

  IF (MAXVAL(NIZ1) > 1.0E-10_8) CALL CALC_GRADIENT(NA1, HRO, NIZ1, dkes2astra_variables(:,5))
  IF (MAXVAL(NIZ2) > 1.0E-10_8) CALL CALC_GRADIENT(NA1, HRO, NIZ2, dkes2astra_variables(:,6))
  IF (MAXVAL(NIZ3) > 1.0E-10_8) CALL CALC_GRADIENT(NA1, HRO, NIZ3, dkes2astra_variables(:,7))

! if electric field goes banaNaNs, restart
 do i=1,nrd
  if (isnan(er_field(i))) er_field(i)=0.
 enddo

  i_option=1  !  0 - full calculation in the interface. 1 only calculate fluxes at 1 position, then astra does Er
! be careful, er_field or electric field in input here must be kV/m !!!, so ER in astra is er_field*1000.
! dkes2astra_variables(1:NA1,8) output is the radial current as source to Ser_field.
    DO i = 2, NA1
      CALL CALCULATE_FLUXES_FOR_ASTRA(NA1, BTOR, ZMJ, AMJ, & 
        TE(1:NA1),TI(1:NA1),NE(1:NA1),NI(1:NA1), & 
        dkes2astra_variables(1:NA1,1),dkes2astra_variables(1:NA1,2), & 
	dkes2astra_variables(1:NA1,3),dkes2astra_variables(1:NA1,4), & 
        NIZ1(1:NA1), ZIM1(1:NA1), AIM1, dkes2astra_variables(1:NA1,5), &
        NIZ2(1:NA1), ZIM2(1:NA1), AIM2, dkes2astra_variables(1:NA1,6), &
        NIZ3(1:NA1), ZIM3(1:NA1), AIM3, dkes2astra_variables(1:NA1,7), &
        RHO(1:NA1),VR(1:NA1),VRS(1:NA1),er_field(1:NA1), dkes2astra_variables(1:NA1,8), &
        dkes2astra_variables(1:NA1,9),dkes2astra_variables(1:NA1,10),&
        dkes2astra_variables(1:NA1,11),dkes2astra_variables(1:NA1,12),&
        dkes2astra_variables(1:NA1,13),dkes2astra_variables(1:NA1,14),&        
        dkes2astra_variables(1:NA1,15),dkes2astra_variables(1:NA1,16),&        
        dkes2astra_variables(1:NA1,17),dkes2astra_variables(1:NA1,18),&        
        dkes2astra_variables(1:NA1,19),dkes2astra_variables(1:NA1,20),&
        dkes2astra_variables(1:NA1,21),dkes2astra_variables(1:NA1,22),&
        dkes2astra_variables(1:NA1,23),dkes2astra_variables(1:NA1,24),&
        dkes2astra_variables(1:NA1,25),dkes2astra_variables(1:NA1,26),&
        dkes2astra_variables(1:NA1,27),dkes2astra_variables(1:NA1,28),&
        dkes2astra_variables(1:NA1,29),dkes2astra_variables(1:NA1,30),&
        dkes2astra_variables(1:NA1,31),dkes2astra_variables(1:NA1,32),&
        1,i, G11(1:NA1))
    END DO

 do i=1,nrd
  if (isnan(dkes2astra_variables(i,14))) dkes2astra_variables(i,14)=1.
  if (isnan(dkes2astra_variables(i,15))) dkes2astra_variables(i,15)=0.
  if (isnan(dkes2astra_variables(i,8))) dkes2astra_variables(i,8)=0.
 enddo

! put radial current to zero at some edge point... defined by user via no_edge = rho_tor_norm of no er region 
 i_edge = nint(no_edge*NA1)
 dkes2astra_variables(i_edge:NA1,8) = 0.0

! Input: 

! dkv(1) --> dne/drho
! dkv(2) --> dTe/drho
! dkv(3) --> dTi/drho
! dkv(4) --> dni/drho
! dkv(5) --> dniz1/drho
! dkv(6) --> dniz2/drho
! dkv(7) --> dniz3/drho

! Output: 

! dkv(8) --> -Jradial/epsilon_perp (source for Er equation)

! dkv(9) --> elec part flux
! dkv(10) --> elec heat flux
! dkv(11) --> total ion heat flux
! dkv(12) --> Jradial = Z_i*gamma_i-gamma_e
! dkv(13) --> main ion part flux

! dkv(14) = sigma_parallel
! dkv(15) = j_bootstrap_parallel in kA/m2
! dkv(16) = j_bootstrap_e
! dkv(17) = j_bootstrap_i
! dkv(18)    : Impurity 1 heat flux (Q_z1)
! dkv(19)    : Impurity 1 bootstrap current (J_bs_z1)
! dkv(20)    : Impurity 2 heat flux (Q_z2)
! dkv(21)    : Impurity 2 bootstrap current (J_bs_z2)
! dkv(22)    : Impurity 3 heat flux (Q_z3)
! dkv(23)      : Impurity 3 bootstrap current (J_bs_z3)

! dkv(24) = HE_NC
! dkv(25) = DE_NC
! dkv(26) = XI_NC
! dkv(27) = DI_NC
! dkv(28) = CE_NC
! dkv(29) = CI_NC
! dkv(30)    : Impurity 1 particle flux (Gamma_z1)
! dkv(31)    : Impurity 2 particle flux (Gamma_z2)
! dkv(32)    : Impurity 3 particle flux (Gamma_z3)


return



CONTAINS

  SUBROUTINE CALC_GRADIENT(N_PTS, H_STEP, PROF_IN, GRAD_OUT)
    INTEGER, INTENT(IN) :: N_PTS
    REAL(KIND=8), INTENT(IN) :: H_STEP
    REAL(KIND=8), DIMENSION(N_PTS), INTENT(IN) :: PROF_IN
    REAL(KIND=8), DIMENSION(N_PTS), INTENT(OUT) :: GRAD_OUT
    INTEGER :: k

    GRAD_OUT(1) = 0.0_8
    
    k = 2
    GRAD_OUT(k) = (PROF_IN(k+1) - PROF_IN(k-1)) / (2.0_8 * H_STEP)

    DO k = 3, N_PTS-2
        GRAD_OUT(k) = (-PROF_IN(k+2) + 8.0_8*PROF_IN(k+1) - 8.0_8*PROF_IN(k-1) + PROF_IN(k-2)) / (12.0_8 * H_STEP)
    END DO

    k = N_PTS-1
    GRAD_OUT(k) = (PROF_IN(k+1) - PROF_IN(k-1)) / (2.0_8 * H_STEP)

    k = N_PTS
    GRAD_OUT(k) = (3.0_8*PROF_IN(k) - 4.0_8*PROF_IN(k-1) + PROF_IN(k-2)) / (2.0_8 * H_STEP)
  END SUBROUTINE CALC_GRADIENT



END SUBROUTINE DKES_INTERFACE_WRAPPER
