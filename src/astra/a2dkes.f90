MODULE a2dkes

  IMPLICIT NONE

  CHARACTER(LEN=256) :: DKES_FILE_PATH
  CHARACTER(LEN=256) :: VMEC_HEADER_FILE
  CHARACTER(LEN=256) :: B00_PROFILE_FILE
  CHARACTER(LEN=256) :: MINOR_RADIUS_W7AS_FILE
  NAMELIST /ASTRA_DKES_INTERFACE/ DKES_FILE_PATH, VMEC_HEADER_FILE, B00_PROFILE_FILE, MINOR_RADIUS_W7AS_FILE

  INTEGER, SAVE :: nr
  
  LOGICAL, SAVE :: uses_new_gradient_def = .FALSE.

  REAL(KIND=8), ALLOCATABLE, SAVE :: F1_previous_step(:)
  REAL(KIND=8), ALLOCATABLE, SAVE :: gamma_e_prev(:), q_e_prev(:)
  REAL(KIND=8), ALLOCATABLE, SAVE :: gamma_i_prev(:), q_i_prev(:)
  REAL(KIND=8), ALLOCATABLE, SAVE :: B00_PHYSICAL_PROFILE(:) 
  REAL(KIND=8), PARAMETER :: comparison_epsilon = 1.0e-10_8
  REAL(KIND=8), ALLOCATABLE, SAVE :: last_er_used_for_fluxes(:)
  REAL(KIND=8), PARAMETER :: ELECTRON_MASS = 9.1093837E-31
  REAL(KIND=8), PARAMETER :: PROTON_MASS = 1.6726219E-27
  REAL(KIND=8), PARAMETER :: ELEMENTARY_CHARGE_C = 1.6021766E-19
  REAL(KIND=8), PARAMETER :: PERMITTIVITY = 8.8541878E-12
  REAL(KIND=8), PARAMETER :: MY_PI = 3.14159265358979323846
  
  
  TYPE dkes_radial_block
    REAL(KIND=8) :: r_val, R00, B00_norm, iota, kn, ftrap, b2_norm 
    REAL(KIND=8) :: g11_0, EFIELD_upper, g11_Er, ex_Er
    INTEGER :: num_points
    REAL(KIND=8), ALLOCATABLE :: cmul_list(:)
    REAL(KIND=8), ALLOCATABLE :: efield_list(:)
    REAL(KIND=8), ALLOCATABLE :: g11_list(:)
    REAL(KIND=8), ALLOCATABLE :: g31_list(:)
    REAL(KIND=8), ALLOCATABLE :: g33_list(:)
    INTEGER :: num_unique_cmuls
    REAL(KIND=8), ALLOCATABLE :: unique_cmuls(:)
    INTEGER, ALLOCATABLE :: slice_counts(:)
    REAL(KIND=8), ALLOCATABLE :: spline_x(:,:)
    REAL(KIND=8), ALLOCATABLE :: spline_x0(:,:)
    REAL(KIND=8), ALLOCATABLE :: spline_x1(:,:)
    REAL(KIND=8), ALLOCATABLE :: spline_x2(:,:)
    REAL(KIND=8), ALLOCATABLE :: spline_y2_0(:,:)
    REAL(KIND=8), ALLOCATABLE :: spline_y2_1(:,:)
    REAL(KIND=8), ALLOCATABLE :: spline_y2_2(:,:)
  END TYPE dkes_radial_block
  
  TYPE(dkes_radial_block), ALLOCATABLE, SAVE :: dkes_data(:)
   
  INTEGER, PARAMETER :: N_GAUSS_LAGUERRE = 25
  
  REAL(KIND=8), PARAMETER, DIMENSION(N_GAUSS_LAGUERRE) :: gauss_laguerre_x = [ &
    2.24159D-02, 1.18123D-01, 2.90366D-01, 5.39286D-01, 8.65037D-01, &
    1.26781D+00, 1.74786D+00, 2.30546D+00, 2.94096D+00, 3.65475D+00, &
    4.44727D+00, 5.31900D+00, 6.27050D+00, 7.30237D+00, 8.41528D+00, &
    9.60994D+00, 1.08872D+01, 1.22478D+01, 1.36927D+01, 1.52230D+01, &
    1.68397D+01, 1.85439D+01, 2.03370D+01, 2.22202D+01, 2.41951D+01 ]

  REAL(KIND=8), PARAMETER, DIMENSION(N_GAUSS_LAGUERRE) :: gauss_laguerre_w = [ &
    5.62529D-02, 1.19024D-01, 1.57496D-01, 1.67547D-01, 1.53353D-01, &
    1.24221D-01, 9.03423D-02, 5.94778D-02, 3.56276D-02, 1.94804D-02, &
    9.74360D-03, 4.46431D-03, 1.87536D-03, 7.22647D-04, 2.55488D-04, &
    8.28714D-05, 2.46569D-05, 6.72671D-06, 1.68179D-06, 3.85081D-07, &
    8.06883D-08, 1.54573D-08, 2.70450D-09, 4.31679D-10, 6.27784D-11 ]


  LOGICAL, SAVE :: dkes_is_initialized = .FALSE.
  
  REAL(KIND=8), SAVE :: module_ABC
  REAL(KIND=8), SAVE :: module_psi_a
  REAL(KIND=8), SAVE :: minor_radius_W7AS_value

  PUBLIC :: CALCULATE_FLUXES_FOR_ASTRA

CONTAINS

  ! Main subroutine called by ASTRA at each time step
  SUBROUTINE CALCULATE_FLUXES_FOR_ASTRA(NA1,BTOR,ZMJ,AMJ, &
       TE,TI,NE,NI,CAR10,CAR11,CAR12,NI_GRAD, &
       CAR50,ZIM1,AIM1,CAR13, &
       CAR51,ZIM2,AIM2,CAR14, &
       CAR52,ZIM3,AIM3,CAR15, &
       RHO,VPRIME,VPRIME_S,F1,CAR30, &
       CAR21,CAR22,CAR23,CAR24,CAR25,CAR26,CAR27,CAR28,CAR29, &
       CAR31,CAR32,CAR33,CAR34,CAR35,CAR36, &
       CAR41,CAR42,CAR43,CAR44,CAR46,CAR47, &
       CAR37,CAR38,CAR39, &
       i_option,j_pos,G11)
    IMPLICIT NONE

    integer, intent(IN):: NA1, j_pos, i_option
    double precision, intent(IN):: BTOR,ZMJ,AMJ
    double precision, intent(IN):: AIM1, AIM2, AIM3
    double precision, dimension(NA1), intent(IN):: ZIM1, ZIM2, ZIM3
    double precision, dimension(NA1), intent(IN):: TE, TI, NE, NI, VPRIME, VPRIME_S, G11
    double precision, dimension(NA1), intent(IN):: CAR10,CAR11,CAR12, NI_GRAD
    double precision, dimension(NA1), intent(IN):: CAR50, CAR13
    double precision, dimension(NA1), intent(IN):: CAR51, CAR14
    double precision, dimension(NA1), intent(IN):: CAR52, CAR15
    double precision, dimension(NA1) :: RHO
    double precision, dimension(NA1), intent(INOUT):: F1
    double precision, dimension(NA1), intent(OUT):: CAR30, CAR21, CAR22, CAR23,CAR24,CAR25, CAR26, CAR27, CAR28, CAR29
    double precision, dimension(NA1), intent(OUT):: CAR41, CAR42, CAR43, CAR44, CAR46, CAR47
    double precision, dimension(NA1), intent(OUT):: CAR31, CAR32, CAR33, CAR34, CAR35, CAR36
    double precision, dimension(NA1), intent(OUT):: CAR37, CAR38, CAR39

    INTEGER :: J, er_status
    INTEGER :: unit_b00_profile
    REAL(KIND=8) :: dlnndr_rho, dlntedr_rho, dlntidr_rho, dlnidr_rho
    REAL(KIND=8) :: J_radial
    REAL(KIND=8) :: r_dkes_lookup, rho_astra, ROC
    REAL(KIND=8) :: b00_physical_loc, s_astra
    REAL(KIND=8) :: r_physical, r_lookup_scaled
    REAL(KIND=8) :: sigma_parallel
    REAL(KIND=8) :: j_bootstrap_parallel
    REAL(KIND=8) :: j_bootstrap_e
    REAL(KIND=8) :: j_bootstrap_i
    REAL(KIND=8) :: j_bootstrap_imp1, j_bootstrap_imp2, j_bootstrap_imp3
    REAL(KIND=8) :: HE_NC, DE_NC, XI_NC, DI_NC, CE_NC, CI_NC
    REAL(KIND=8) :: epsilon_local
    REAL(KIND=8) :: drho_dr

    REAL(KIND=8) :: b00_physical_scaled, dr_over_drtilde_loc, b2_norm_interpolated

    INTEGER :: nml_unit, io_stat
    CHARACTER(LEN=256) :: namelist_file

    INTEGER :: N_spec, ispec
    REAL(KIND=8) :: Z_spec(5), A_spec(5)
    REAL(KIND=8) :: dens_loc(5), temp_loc(5), grad_n_loc(5), grad_T_loc(5)
    REAL(KIND=8) :: flux_part_loc(5), flux_heat_loc(5)
    REAL(KIND=8) :: total_ion_heat_flux  

    ROC=RHO(NA1)
     
    IF (.NOT. dkes_is_initialized) THEN
      namelist_file = 'vmec_io/stell_files.nml'
      OPEN(newunit=nml_unit, file=TRIM(namelist_file), status='old', action='read', iostat=io_stat)
      IF (io_stat /= 0) STOP 'Namelist Read Error'
      READ(nml_unit, NML=ASTRA_DKES_INTERFACE, iostat=io_stat)
      CLOSE(nml_unit)
      IF (io_stat /= 0) STOP 'Namelist Read Error'
      
      CALL read_minorradiusW7AS(er_status)
      IF (er_status /= 0) STOP 'Minor r W7AS error'

      CALL read_vmec_header_data(er_status)
      IF (er_status /= 0) STOP 'VMEC Header Read Error'
      OPEN(newunit=unit_b00_profile, file=TRIM(B00_PROFILE_FILE), status='old', action='read', iostat=io_stat)
      IF (io_stat /= 0) STOP 'B00 Read Error'

      ALLOCATE(B00_PHYSICAL_PROFILE(NA1))
      READ(unit_b00_profile, *) B00_PHYSICAL_PROFILE
      CLOSE(unit_b00_profile)

      CALL read_dkes_data(DKES_FILE_PATH, er_status)
      IF (er_status == 0) THEN
        dkes_is_initialized = .TRUE.
        ALLOCATE(last_er_used_for_fluxes(NA1))
        ALLOCATE(F1_previous_step(NA1))
        F1_previous_step = 0.0_8
      ELSE
        STOP 'DKES read error'
      END IF
      ALLOCATE(gamma_e_prev(NA1), q_e_prev(NA1), gamma_i_prev(NA1), q_i_prev(NA1))
      gamma_e_prev = 0.0_8; q_e_prev = 0.0_8;
      gamma_i_prev = 0.0_8; q_i_prev = 0.0_8;
    END IF
    
    N_spec = 2
    IF (MAXVAL(CAR50) > 1.0E-10_8) N_spec = 3
    IF (MAXVAL(CAR51) > 1.0E-10_8) N_spec = 4
    IF (MAXVAL(CAR52) > 1.0E-10_8) N_spec = 5

  drho_dr = ROC / module_ABC  

  if (i_option==0) then
    DO J = 1, NA1
      IF (J == 1) THEN
        F1(J) = 0.0_8; CAR30(J) = 0.0_8; CAR21(J) = 0.0_8; CAR22(J) = 0.0_8; CAR23(J) = 0.0_8; CAR24(J) = 0.0_8; CAR25(J) = 0.0_8;
        CAR26(J) = 1.0E-6_8; CAR27(J) = 0.0_8; CAR28(J) = 0.0_8; CAR29(J) = 0.0_8
        CAR31(J) = 0.0_8; CAR32(J) = 0.0_8; CAR33(J) = 0.0_8; CAR34(J) = 0.0_8; CAR35(J) = 0.0_8; CAR36(J) = 0.0_8
        CAR37(J) = 0.0_8; CAR38(J) = 0.0_8; CAR39(J) = 0.0_8
        CAR41(J) = 0.0_8; CAR42(J) = 0.0_8; CAR43(J) = 0.0_8; CAR44(J) = 0.0_8; CAR46(J) = 0.0_8; CAR47(J) = 0.0_8
        last_er_used_for_fluxes(J) = 0.0_8
        CYCLE
      ENDIF
      
      ! Map from ASTRA rho coordinate to the physical minor radius r
      rho_astra = RHO(J)
      s_astra = (rho_astra / ROC)**2
      r_physical = module_ABC * SQRT(s_astra)
      r_dkes_lookup = r_physical / module_ABC * minor_radius_W7AS_value
      !r_dkes_lookup = r_physical
      !WRITE(*,*) 'DEBUG: AMJ = ', AMJ
      !PRINT *, 'rho_astra = ', rho_astra
      !PRINT *, 's_astra = ', s_astra
      !PRINT *, 'r_physical = ',r_physical
      !PRINT *, 'r_dkes_lookup = ', r_dkes_lookup
      
      ! Gradients from ASTRA CAR 
      dlnndr_rho  = CAR10(J) / NE(J)
      dlntedr_rho = CAR11(J) / TE(J)
      dlntidr_rho = CAR12(J) / TI(J)
      dlnidr_rho  = NI_GRAD(J) / NI(J)
      !PRINT *, 'dlntedr_rho = ', dlntedr_rho
      !PRINT *, 'dlntidr_rho = ', dlntidr_rho
      !PRINT *, 'dlnndr_rho = ', dlnndr_rho
      
      ! Geometric factor 
      b00_physical_loc = B00_PHYSICAL_PROFILE(J) 
      b00_physical_scaled = b00_physical_loc * (BTOR / B00_PHYSICAL_PROFILE(1))
      
      IF (uses_new_gradient_def) THEN
          dr_over_drtilde_loc = 1.0_8
      ELSE
          dr_over_drtilde_loc = (module_ABC * minor_radius_W7AS_value * b00_physical_loc) / (2.0_8 * ABS(module_psi_a))
      END IF
      
      b2_norm_interpolated = interpolate_b2_norm_only(r_dkes_lookup)

      Z_spec(1) = -1.0_8; A_spec(1) = ELECTRON_MASS / PROTON_MASS
      dens_loc(1) = NE(J); temp_loc(1) = TE(J)
      grad_n_loc(1) = dlnndr_rho; grad_T_loc(1) = dlntedr_rho

      Z_spec(2) = ZMJ; A_spec(2) = AMJ
      dens_loc(2) = NI(J); temp_loc(2) = TI(J)
      grad_n_loc(2) = dlnidr_rho; grad_T_loc(2) = dlntidr_rho

      IF (N_spec >= 3) THEN
          Z_spec(3) = ZIM1(J); A_spec(3) = AIM1
          dens_loc(3) = CAR50(J); temp_loc(3) = TI(J)
          grad_T_loc(3) = dlntidr_rho
          IF (CAR50(J) > 1.0E-10_8) THEN; grad_n_loc(3) = CAR13(J) / CAR50(J); ELSE; grad_n_loc(3) = 0.0_8; END IF
      END IF

      IF (N_spec >= 4) THEN
          Z_spec(4) = ZIM2(J); A_spec(4) = AIM2
          dens_loc(4) = CAR51(J); temp_loc(4) = TI(J)
          grad_T_loc(4) = dlntidr_rho
          IF (CAR51(J) > 1.0E-10_8) THEN; grad_n_loc(4) = CAR14(J) / CAR51(J); ELSE; grad_n_loc(4) = 0.0_8; END IF
      END IF

      IF (N_spec >= 5) THEN
          Z_spec(5) = ZIM3(J); A_spec(5) = AIM3
          dens_loc(5) = CAR52(J); temp_loc(5) = TI(J)
          grad_T_loc(5) = dlntidr_rho
          IF (CAR52(J) > 1.0E-10_8) THEN; grad_n_loc(5) = CAR15(J) / CAR52(J); ELSE; grad_n_loc(5) = 0.0_8; END IF
      END IF
      
      CALL calculate_all_fluxes_multispecies(J, r_dkes_lookup, N_spec, Z_spec, A_spec, &
                                dens_loc, temp_loc, grad_n_loc, grad_T_loc, &
                                F1_previous_step(J) * 1000.0_8, & 
                                flux_part_loc, flux_heat_loc, sigma_parallel, j_bootstrap_parallel, &
                                j_bootstrap_e, j_bootstrap_i, j_bootstrap_imp1, j_bootstrap_imp2, j_bootstrap_imp3, &
                                b00_physical_loc, b00_physical_scaled, dr_over_drtilde_loc, b2_norm_interpolated, &
                                RHO(J), ROC, G11(J), VPRIME(J), VPRIME_S(J), epsilon_local, &
                                HE_NC, DE_NC, XI_NC, DI_NC, CE_NC, CI_NC)

      J_radial = 0.0_8
      DO ispec = 1, N_spec; J_radial = J_radial + Z_spec(ispec) * flux_part_loc(ispec); END DO
      J_radial = ELEMENTARY_CHARGE_C * J_radial

      total_ion_heat_flux = 0.0_8
      DO ispec = 2, N_spec; total_ion_heat_flux = total_ion_heat_flux + flux_heat_loc(ispec); END DO

      CAR21(J) = flux_part_loc(1) * VPRIME(J) * drho_dr / 1.0D19
      CAR22(J) = flux_heat_loc(1) * VPRIME(J) * drho_dr/ 1.0E6
      CAR23(J) = total_ion_heat_flux * VPRIME(J) * drho_dr / 1.0E6   
      CAR24(J) = J_radial
      CAR25(J) = flux_part_loc(2) * VPRIME(J) * drho_dr / 1.0D19
      CAR30(J) = -(J_radial / epsilon_local) / 1000.0
      CAR26(J) = sigma_parallel
      CAR27(J) = j_bootstrap_parallel
      CAR28(J) = j_bootstrap_e
      CAR29(J) = j_bootstrap_i
      CAR41(J) = HE_NC; CAR42(J) = DE_NC; CAR43(J) = XI_NC; CAR44(J) = DI_NC; CAR46(J) = CE_NC; CAR47(J) = CI_NC
      
      IF (N_spec >= 3) THEN
          CAR37(J) = flux_part_loc(3) * VPRIME(J) * drho_dr / 1.0D19
          CAR31(J) = flux_heat_loc(3) * VPRIME(J) * drho_dr / 1.0E6
          CAR32(J) = j_bootstrap_imp1
      END IF
      IF (N_spec >= 4) THEN
          CAR38(J) = flux_part_loc(4) * VPRIME(J) * drho_dr / 1.0D19
          CAR33(J) = flux_heat_loc(4) * VPRIME(J) * drho_dr / 1.0E6
          CAR34(J) = j_bootstrap_imp2
      END IF
      IF (N_spec >= 5) THEN
          CAR39(J) = flux_part_loc(5) * VPRIME(J) * drho_dr / 1.0D19
          CAR35(J) = flux_heat_loc(5) * VPRIME(J) * drho_dr / 1.0E6
          CAR36(J) = j_bootstrap_imp3
      END IF

    END DO

    ! Current Er profile for the next time step
    F1_previous_step = F1
  endif   

  if (i_option==1) then
    ! Main loop over ASTRA radial grid
    DO J = j_pos, j_pos
        CAR30(J) = 0.0_8; CAR21(J) = 0.0_8; CAR22(J) = 0.0_8; CAR23(J) = 0.0_8; CAR24(J) = 0.0_8; CAR25(J) = 0.0_8;
        CAR26(J) = 1.0E-6_8; CAR27(J) = 0.0_8; CAR28(J) = 0.0_8; CAR29(J) = 0.0_8;
        CAR31(J) = 0.0_8; CAR32(J) = 0.0_8; CAR33(J) = 0.0_8; CAR34(J) = 0.0_8; CAR35(J) = 0.0_8; CAR36(J) = 0.0_8;
        CAR37(J) = 0.0_8; CAR38(J) = 0.0_8; CAR39(J) = 0.0_8;
        CAR41(J) = 0.0_8; CAR42(J) = 0.0_8; CAR43(J) = 0.0_8; CAR44(J) = 0.0_8; CAR46(J) = 0.0_8; CAR47(J) = 0.0_8
       
      rho_astra = RHO(J)
      s_astra = (rho_astra / ROC)**2
      r_physical = module_ABC * SQRT(s_astra)
      r_dkes_lookup = r_physical / module_ABC * minor_radius_W7AS_value
      !r_dkes_lookup = r_physical
      !WRITE(*,*) 'DEBUG: AMJ = ', AMJ
      !PRINT *, 'rho_astra = ', rho_astra
      !PRINT *, 's_astra = ', s_astra
      !PRINT *, 'r_physical = ',r_physical
      !PRINT *, 'r_dkes_lookup = ', r_dkes_lookup
      
      ! Gradients from ASTRA CAR 
      dlnndr_rho  = CAR10(J) / NE(J)
      dlntedr_rho = CAR11(J) / TE(J)
      dlntidr_rho = CAR12(J) / TI(J)
      dlnidr_rho  = NI_GRAD(J) / NI(J)
      !PRINT *, 'dlntedr_rho = ', dlntedr_rho
      !PRINT *, 'dlntidr_rho = ', dlntidr_rho
      !PRINT *, 'dlnndr_rho = ', dlnndr_rho
      
      ! Geometric factor
      b00_physical_loc = B00_PHYSICAL_PROFILE(J) 
      b00_physical_scaled = b00_physical_loc * (BTOR / B00_PHYSICAL_PROFILE(1))
      IF (uses_new_gradient_def) THEN
          dr_over_drtilde_loc = 1.0_8
      ELSE
          dr_over_drtilde_loc = (module_ABC * minor_radius_W7AS_value * b00_physical_loc) / (2.0_8 * ABS(module_psi_a))
      END IF
      b2_norm_interpolated = interpolate_b2_norm_only(r_dkes_lookup)

      Z_spec(1) = -1.0_8; A_spec(1) = ELECTRON_MASS / PROTON_MASS
      dens_loc(1) = NE(J); temp_loc(1) = TE(J)
      grad_n_loc(1) = dlnndr_rho; grad_T_loc(1) = dlntedr_rho

      Z_spec(2) = ZMJ; A_spec(2) = AMJ
      dens_loc(2) = NI(J); temp_loc(2) = TI(J)
      grad_n_loc(2) = dlnidr_rho; grad_T_loc(2) = dlntidr_rho

      IF (N_spec >= 3) THEN
          Z_spec(3) = ZIM1(J); A_spec(3) = AIM1
          dens_loc(3) = CAR50(J); temp_loc(3) = TI(J)
          grad_T_loc(3) = dlntidr_rho
          IF (CAR50(J) > 1.0E-10_8) THEN; grad_n_loc(3) = CAR13(J) / CAR50(J); ELSE; grad_n_loc(3) = 0.0_8; END IF
      END IF

      IF (N_spec >= 4) THEN
          Z_spec(4) = ZIM2(J); A_spec(4) = AIM2
          dens_loc(4) = CAR51(J); temp_loc(4) = TI(J)
          grad_T_loc(4) = dlntidr_rho
          IF (CAR51(J) > 1.0E-10_8) THEN; grad_n_loc(4) = CAR14(J) / CAR51(J); ELSE; grad_n_loc(4) = 0.0_8; END IF
      END IF

      IF (N_spec >= 5) THEN
          Z_spec(5) = ZIM3(J); A_spec(5) = AIM3
          dens_loc(5) = CAR52(J); temp_loc(5) = TI(J)
          grad_T_loc(5) = dlntidr_rho
          IF (CAR52(J) > 1.0E-10_8) THEN; grad_n_loc(5) = CAR15(J) / CAR52(J); ELSE; grad_n_loc(5) = 0.0_8; END IF
      END IF

      CALL calculate_all_fluxes_multispecies(J, r_dkes_lookup, N_spec, Z_spec, A_spec, &
                                dens_loc, temp_loc, grad_n_loc, grad_T_loc, &
                                F1(J) * 1000.0_8, &
                                flux_part_loc, flux_heat_loc, sigma_parallel, j_bootstrap_parallel, &
                                j_bootstrap_e, j_bootstrap_i, j_bootstrap_imp1, j_bootstrap_imp2, j_bootstrap_imp3, &
                                b00_physical_loc, b00_physical_scaled, dr_over_drtilde_loc, b2_norm_interpolated, &
                                RHO(J), ROC, G11(J), VPRIME(J), VPRIME_S(J), epsilon_local, &
                                HE_NC, DE_NC, XI_NC, DI_NC, CE_NC, CI_NC)

      J_radial = 0.0_8
      DO ispec = 1, N_spec; J_radial = J_radial + Z_spec(ispec) * flux_part_loc(ispec); END DO
      J_radial = ELEMENTARY_CHARGE_C * J_radial

      total_ion_heat_flux = 0.0_8
      DO ispec = 2, N_spec; total_ion_heat_flux = total_ion_heat_flux + flux_heat_loc(ispec); END DO

      CAR21(J) = flux_part_loc(1) * VPRIME(J) * drho_dr / 1.0D19
      CAR22(J) = flux_heat_loc(1) * VPRIME(J) * drho_dr / 1.0E6
      CAR23(J) = total_ion_heat_flux * VPRIME(J) * drho_dr / 1.0E6  
      CAR24(J) = J_radial
      CAR25(J) = flux_part_loc(2) * VPRIME(J) * drho_dr / 1.0D19
      CAR30(J) = -(J_radial / epsilon_local) / 1000.0
      CAR26(J) = sigma_parallel
      CAR27(J) = j_bootstrap_parallel
      CAR28(J) = j_bootstrap_e
      CAR29(J) = j_bootstrap_i
      CAR41(J) = HE_NC; CAR42(J) = DE_NC; CAR43(J) = XI_NC; CAR44(J) = DI_NC; CAR46(J) = CE_NC; CAR47(J) = CI_NC
      
      IF (N_spec >= 3) THEN
          CAR37(J) = flux_part_loc(3) * VPRIME(J) * drho_dr / 1.0D19
          CAR31(J) = flux_heat_loc(3) * VPRIME(J) * drho_dr / 1.0E6
          CAR32(J) = j_bootstrap_imp1
      END IF
      IF (N_spec >= 4) THEN
          CAR38(J) = flux_part_loc(4) * VPRIME(J) * drho_dr / 1.0D19
          CAR33(J) = flux_heat_loc(4) * VPRIME(J) * drho_dr / 1.0E6
          CAR34(J) = j_bootstrap_imp2
      END IF
      IF (N_spec >= 5) THEN
          CAR39(J) = flux_part_loc(5) * VPRIME(J) * drho_dr / 1.0D19
          CAR35(J) = flux_heat_loc(5) * VPRIME(J) * drho_dr / 1.0E6
          CAR36(J) = j_bootstrap_imp3
      END IF

    END DO
  endif

  END SUBROUTINE CALCULATE_FLUXES_FOR_ASTRA


  SUBROUTINE calculate_all_fluxes_multispecies(J, r_lookup_in, N_spec, Z_spec, A_spec, &
                                dens_loc, temp_loc, grad_n_loc, grad_T_loc, Er_loc_Vm, &
                                flux_part_out, flux_heat_out, sigma_parallel_out, j_bootstrap_out, &
                                jbs_e_out, jbs_i_out, jbs_imp1_out, jbs_imp2_out, jbs_imp3_out, &
                                b00_physical_at_r, b00_physical_scaled, dr_over_drtilde_loc, b2_norm_interpolated, &
                                RHO, ROC_in, G11_loc, VPRIME_loc, VPRIME_S_loc, epsilon_perp_out, &
                                HE_NC_out, DE_NC_out, XI_NC_out, DI_NC_out, CE_NC_out, CI_NC_out)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: J, N_spec
    REAL(KIND=8), INTENT(IN) :: r_lookup_in, Er_loc_Vm
    REAL(KIND=8), DIMENSION(N_spec), INTENT(IN) :: Z_spec, A_spec, dens_loc, temp_loc, grad_n_loc, grad_T_loc
    REAL(KIND=8), INTENT(IN) :: b00_physical_at_r, b00_physical_scaled, dr_over_drtilde_loc, b2_norm_interpolated
    REAL(KIND=8), INTENT(IN) :: RHO, ROC_in, G11_loc, VPRIME_loc, VPRIME_S_loc
    
    REAL(KIND=8), DIMENSION(N_spec), INTENT(OUT) :: flux_part_out, flux_heat_out
    REAL(KIND=8), INTENT(OUT) :: sigma_parallel_out, j_bootstrap_out, jbs_e_out, jbs_i_out
    REAL(KIND=8), INTENT(OUT) :: jbs_imp1_out, jbs_imp2_out, jbs_imp3_out
    REAL(KIND=8), INTENT(OUT) :: HE_NC_out, DE_NC_out, XI_NC_out, DI_NC_out, CE_NC_out, CI_NC_out, epsilon_perp_out
    
    REAL(KIND=8), DIMENSION(N_spec) :: L11, L12, L22
    REAL(KIND=8), DIMENSION(N_spec) :: M_00, M_01, M_02, M_11, M_12, M_22
    REAL(KIND=8), DIMENSION(N_spec) :: N_00, N_01, N_02, N_11, N_12, N_22
    REAL(KIND=8), DIMENSION(N_spec) :: grad_n_r, grad_T_r, er_term
    REAL(KIND=8), DIMENSION(N_spec) :: A1_DKES, A2_DKES
    REAL(KIND=8), DIMENSION(N_spec) :: X1_SN, X2_SN
    REAL(KIND=8), DIMENSION(N_spec) :: jbs_array_corr

    INTEGER :: ispec
    CHARACTER(LEN=10) :: species_name
    REAL(KIND=8) :: drho_dr_factor
    REAL(KIND=8) :: m_ion_kg, n_ion_SI, ASTRA_map_factor, grad_map, er_map
    REAL(KIND=8), PARAMETER :: MJ_PER_KEV_1E19 = 1.6021766E-3_8

    drho_dr_factor = ROC_in / module_ABC

    DO ispec = 1, N_spec
        IF (Z_spec(ispec) < 0.0_8) THEN; species_name = 'electron'; ELSE; species_name = 'ion'; END IF

        CALL energy_convolution_multispecies(J, r_lookup_in, ispec, N_spec, &
             dens_loc, temp_loc, Z_spec, A_spec, Er_loc_Vm, species_name, &
             L11(ispec), L12(ispec), L22(ispec), &
             M_00(ispec), M_01(ispec), M_02(ispec), M_11(ispec), M_12(ispec), M_22(ispec), &
             N_00(ispec), N_01(ispec), N_02(ispec), N_11(ispec), N_12(ispec), N_22(ispec), &
             b00_physical_scaled, b00_physical_at_r, dr_over_drtilde_loc, b2_norm_interpolated)

        grad_n_r(ispec) = grad_n_loc(ispec) * drho_dr_factor
        grad_T_r(ispec) = grad_T_loc(ispec) * drho_dr_factor

        er_term(ispec) = -(Z_spec(ispec) * ELEMENTARY_CHARGE_C * Er_loc_Vm) / (temp_loc(ispec) * 1000.0_8 * ELEMENTARY_CHARGE_C)
        A1_DKES(ispec) = grad_n_r(ispec) - 1.5_8 * grad_T_r(ispec) + er_term(ispec)
        A2_DKES(ispec) = grad_T_r(ispec)
	
        ! Sugama-Nishimura thermodynamic forces 
        ! X1 = - T * (n'/n + T'/T - q*Er/T)
        X1_SN(ispec) = - (temp_loc(ispec) * 1000.0_8 * ELEMENTARY_CHARGE_C) * (grad_n_r(ispec) + grad_T_r(ispec) + er_term(ispec))
        ! X2 = + T * (T'/T) 
        X2_SN(ispec) =   (temp_loc(ispec) * 1000.0_8 * ELEMENTARY_CHARGE_C) * grad_T_r(ispec)
    END DO
    
    CALL apply_momentum_correction_multispecies(N_spec, dens_loc, temp_loc, Z_spec, A_spec, &
         L11, L12, L22, M_00, M_01, M_02, M_11, M_12, M_22, N_00, N_01, N_02, N_11, N_12, N_22, &
         X1_SN, X2_SN, A1_DKES, A2_DKES, b2_norm_interpolated, Er_loc_Vm, drho_dr_factor, r_lookup_in, &
         b00_physical_scaled, flux_part_out, flux_heat_out, jbs_array_corr, sigma_parallel_out)
    
    j_bootstrap_out = SUM(jbs_array_corr) / 1000.0_8
    jbs_e_out = jbs_array_corr(1) / 1000.0_8
    jbs_i_out = jbs_array_corr(2) / 1000.0_8
    
    IF (N_spec >= 3) THEN; jbs_imp1_out = jbs_array_corr(3) / 1000.0_8; ELSE; jbs_imp1_out = 0.0_8; END IF
    IF (N_spec >= 4) THEN; jbs_imp2_out = jbs_array_corr(4) / 1000.0_8; ELSE; jbs_imp2_out = 0.0_8; END IF
    IF (N_spec >= 5) THEN; jbs_imp3_out = jbs_array_corr(5) / 1000.0_8; ELSE; jbs_imp3_out = 0.0_8; END IF

    m_ion_kg = PROTON_MASS * A_spec(2) 
    n_ion_SI = dens_loc(2) * 1.0D19
    epsilon_perp_out = (n_ion_SI * m_ion_kg) / (b00_physical_at_r**2)
    
    ASTRA_map_factor = VPRIME_S_loc / G11_loc          ! = 1/<(grad rho)^2>, exact
    grad_map = drho_dr_factor**2 * ASTRA_map_factor
    er_map   = drho_dr_factor    * ASTRA_map_factor
    
    HE_NC_out = (L22(1) - 1.5_8*L12(1)) * grad_map
    DE_NC_out =  L12(1)                 * grad_map
    XI_NC_out = (L22(2) - 1.5_8*L12(2)) * grad_map
    DI_NC_out =  L12(2)                 * grad_map
    CE_NC_out =  L12(1) * (-1.0_8*Er_loc_Vm)    / (temp_loc(1)*1000.0_8) * er_map
    CI_NC_out =  L12(2) * (Z_spec(2)*Er_loc_Vm) / (temp_loc(2)*1000.0_8) * er_map
    
  END SUBROUTINE calculate_all_fluxes_multispecies

  ! Energy convolution to get thermal coefficients for N species
  SUBROUTINE energy_convolution_multispecies(J, r_lookup, ispec, N_spec, &
                                dens_loc, temp_loc, Z_spec, A_spec, Er_loc_Vm, species_name, &
                                L11_out, L12_out, L22_out, &
                                M_00_out, M_01_out, M_02_out, M_11_out, M_12_out, M_22_out, &
                                N_00_out, N_01_out, N_02_out, N_11_out, N_12_out, N_22_out, &
                                b00_physical_scaled, b00_physical_at_r, dr_over_drtilde_loc, b2_norm_val)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: J, ispec, N_spec
    REAL(KIND=8), INTENT(IN) :: r_lookup, Er_loc_Vm, b00_physical_scaled, dr_over_drtilde_loc, b2_norm_val, b00_physical_at_r
    REAL(KIND=8), DIMENSION(N_spec), INTENT(IN) :: dens_loc, temp_loc, Z_spec, A_spec
    CHARACTER(LEN=*), INTENT(IN) :: species_name
    REAL(KIND=8), INTENT(OUT) :: L11_out, L12_out, L22_out
    REAL(KIND=8), INTENT(OUT) :: M_00_out, M_01_out, M_02_out, M_11_out, M_12_out, M_22_out 
    REAL(KIND=8), INTENT(OUT) :: N_00_out, N_01_out, N_02_out, N_11_out, N_12_out, N_22_out 
    
    INTEGER :: i, extrapk_py, min_idx
    REAL(KIND=8) :: m_s, q_s, v_th, T_joules, K, w, v_k, efield_k, x
    REAL(KIND=8) :: min_val, current_val, f_curr, f_next
    REAL(KIND=8) :: dMat_11, dMat_12, dMat_22
    REAL(KIND=8) :: dMat_M_00, dMat_M_01, dMat_M_02, dMat_M_11, dMat_M_12, dMat_M_22
    REAL(KIND=8) :: dMat_N_00, dMat_N_01, dMat_N_02, dMat_N_11, dMat_N_12, dMat_N_22
    REAL(KIND=8) :: fak1, bracfac, oneminusbracfac, nu_sik, OmegaoverB_si, W_sik
    REAL(KIND=8) :: D_mesik_1, D_mesik_2, M_esk_core, N_esk_core
    REAL(KIND=8) :: dr_over_drtilde_sq
    REAL(KIND=8) :: Son_0, Son_1, Son_2, weight_factor
    
    REAL(KIND=8) :: cmul_arr(N_GAUSS_LAGUERRE), g11_arr(N_GAUSS_LAGUERRE)
    REAL(KIND=8) :: g31_arr(N_GAUSS_LAGUERRE), g33_arr(N_GAUSS_LAGUERRE)
    REAL(KIND=8) :: Y_arr(N_GAUSS_LAGUERRE)
    REAL(KIND=8), PARAMETER :: g33_extrap_lim = 1.0e-4_8

    REAL(KIND=8) :: py_B00, py_g33_FSAb2, py_FSAB2
    REAL(KIND=8) :: py_D33starCB
    REAL(KIND=8) :: a_ratio
    
    a_ratio = module_ABC / minor_radius_W7AS_value
    
    py_B00 = b00_physical_scaled
    py_g33_FSAb2 = b2_norm_val
    py_FSAB2 = (py_B00**2) * b2_norm_val

    IF (Z_spec(ispec) < 0.0_8) THEN
        m_s = ELECTRON_MASS
    ELSE
        m_s = (A_spec(ispec) * PROTON_MASS) - (Z_spec(ispec) * ELECTRON_MASS)
    END IF
    
    q_s = ELEMENTARY_CHARGE_C * Z_spec(ispec)
      
    dr_over_drtilde_sq = dr_over_drtilde_loc**2

    T_joules = temp_loc(ispec) * 1000.0_8 * ELEMENTARY_CHARGE_C 
    v_th = SQRT(2.0_8 * T_joules / m_s)

    dMat_11 = 0.0_8; dMat_12 = 0.0_8; dMat_22 = 0.0_8
    dMat_M_00 = 0.0_8; dMat_M_01 = 0.0_8; dMat_M_02 = 0.0_8; dMat_M_11 = 0.0_8; dMat_M_12 = 0.0_8; dMat_M_22 = 0.0_8
    dMat_N_00 = 0.0_8; dMat_N_01 = 0.0_8; dMat_N_02 = 0.0_8; dMat_N_11 = 0.0_8; dMat_N_12 = 0.0_8; dMat_N_22 = 0.0_8

    DO i = 1, N_GAUSS_LAGUERRE
        K = gauss_laguerre_x(i); v_k = v_th * SQRT(K)
        cmul_arr(i) = calculate_dkes_cmul_multispecies(ispec, N_spec, K, v_k, m_s, q_s, dens_loc, temp_loc, Z_spec, A_spec)         

        efield_k = normalize_er_for_dkes(Er_loc_Vm, b00_physical_at_r, v_k, dr_over_drtilde_loc)
      
        CALL interpolate_dkes_coeffs(r_lookup, cmul_arr(i), efield_k, g11_arr(i), g31_arr(i), g33_arr(i))
        
        py_D33starCB = -g33_arr(i) * 1.5_8 * cmul_arr(i) / b2_norm_val
        
        Y_arr(i) = 1.0_8 - py_D33starCB
	
    END DO

    min_idx = 0
    IF (Y_arr(1) < 1.0E-4_8) THEN
        DO i = 1, N_GAUSS_LAGUERRE - 1
            IF (Y_arr(i) < 1.0E-4_8 .AND. Y_arr(i+1) >= 1.0E-4_8) THEN
                min_idx = i
                EXIT
            END IF
        END DO
    END IF
        
    
    IF (min_idx > 0) THEN
        DO i = min_idx, 1, -1
            Y_arr(i) = Y_arr(min_idx+1) * (cmul_arr(min_idx+1) / cmul_arr(i))**2
            END DO
    END IF
    

    DO i = 1, N_GAUSS_LAGUERRE
        K = gauss_laguerre_x(i); w = gauss_laguerre_w(i); x = SQRT(K)

        dMat_11 = dMat_11 + w * K**(2.0_8) * g11_arr(i)
        dMat_12 = dMat_12 + w * K**(3.0_8) * g11_arr(i)
        dMat_22 = dMat_22 + w * K**(4.0_8) * g11_arr(i)

        oneminusbracfac = Y_arr(i)
        bracfac = 1.0_8 - oneminusbracfac

        nu_sik = cmul_arr(i) * x * v_th
        OmegaoverB_si = q_s / m_s
        
        W_sik = 1.5_8 * m_s * nu_sik / (T_joules * K * py_FSAB2)

	D_mesik_1 = -g31_arr(i) * (dr_over_drtilde_loc) * K * v_th**2 / (2.0_8 * OmegaoverB_si)
        D_mesik_2 = oneminusbracfac / W_sik

        M_esk_core = (D_mesik_2 / bracfac) * (nu_sik**2) * (m_s**2) / T_joules
        N_esk_core = (D_mesik_1 / bracfac) * nu_sik * m_s / T_joules

        Son_0 = 1.0_8
        Son_1 = 2.5_8 - K
        Son_2 = 0.5_8 * (K**2 - 7.0_8 * K + 8.75_8)
        weight_factor = (2.0_8 / SQRT(MY_PI)) * w * x

        dMat_M_00 = dMat_M_00 + weight_factor * Son_0 * Son_0 * M_esk_core
        dMat_M_01 = dMat_M_01 + weight_factor * Son_0 * Son_1 * M_esk_core
        dMat_M_02 = dMat_M_02 + weight_factor * Son_0 * Son_2 * M_esk_core
        dMat_M_11 = dMat_M_11 + weight_factor * Son_1 * Son_1 * M_esk_core
        dMat_M_12 = dMat_M_12 + weight_factor * Son_1 * Son_2 * M_esk_core
        dMat_M_22 = dMat_M_22 + weight_factor * Son_2 * Son_2 * M_esk_core

        dMat_N_00 = dMat_N_00 + weight_factor * Son_0 * Son_0 * N_esk_core
        dMat_N_01 = dMat_N_01 + weight_factor * Son_0 * Son_1 * N_esk_core
        dMat_N_02 = dMat_N_02 + weight_factor * Son_0 * Son_2 * N_esk_core
        dMat_N_11 = dMat_N_11 + weight_factor * Son_1 * Son_1 * N_esk_core
        dMat_N_12 = dMat_N_12 + weight_factor * Son_1 * Son_2 * N_esk_core
        dMat_N_22 = dMat_N_22 + weight_factor * Son_2 * Son_2 * N_esk_core
    END DO

    fak1  = -(m_s**2 * v_th**3) * dr_over_drtilde_sq / (2.0_8 * q_s**2 * py_B00**2)
    L11_out = fak1  * (2.0_8 / SQRT(MY_PI)) * dMat_11
    L12_out = fak1  * (2.0_8 / SQRT(MY_PI)) * dMat_12
    L22_out = fak1  * (2.0_8 / SQRT(MY_PI)) * dMat_22
    
    M_00_out = dMat_M_00 * dens_loc(ispec) * 1.0D19
    M_01_out = dMat_M_01 * dens_loc(ispec) * 1.0D19
    M_02_out = dMat_M_02 * dens_loc(ispec) * 1.0D19
    M_11_out = dMat_M_11 * dens_loc(ispec) * 1.0D19
    M_12_out = dMat_M_12 * dens_loc(ispec) * 1.0D19
    M_22_out = dMat_M_22 * dens_loc(ispec) * 1.0D19

    N_00_out = dMat_N_00 * dens_loc(ispec) * 1.0D19
    N_01_out = dMat_N_01 * dens_loc(ispec) * 1.0D19
    N_02_out = dMat_N_02 * dens_loc(ispec) * 1.0D19
    N_11_out = dMat_N_11 * dens_loc(ispec) * 1.0D19
    N_12_out = dMat_N_12 * dens_loc(ispec) * 1.0D19
    N_22_out = dMat_N_22 * dens_loc(ispec) * 1.0D19
  END SUBROUTINE energy_convolution_multispecies


  ! Calculates CMUL = nu_total/v for multi-species
  FUNCTION calculate_dkes_cmul_multispecies(ispec, N_spec, K_energy, v_particle, m_s, charge_s, &
                                            dens_loc, temp_loc, Z_spec, A_spec) RESULT(cmul_out)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ispec, N_spec
    REAL(KIND=8), INTENT(IN) :: K_energy, v_particle, m_s, charge_s
    REAL(KIND=8), DIMENSION(N_spec), INTENT(IN) :: dens_loc, temp_loc, Z_spec, A_spec
    REAL(KIND=8) :: cmul_out, nu_total_v
    INTEGER :: jspec
    REAL(KIND=8) :: m_b, Z_a, Z_b

    Z_a = Z_spec(ispec)
    nu_total_v = 0.0_8

    ! Sum collision frequencies over all species in the plasma
    DO jspec = 1, N_spec
        IF (Z_spec(jspec) < 0.0_8) THEN
            m_b = ELECTRON_MASS
        ELSE
            m_b = (A_spec(jspec) * PROTON_MASS) - (Z_spec(jspec) * ELECTRON_MASS)
        END IF
        Z_b = Z_spec(jspec)

        nu_total_v = nu_total_v + nu_ab_v_multispecies(m_s, temp_loc(ispec), v_particle, Z_a, dens_loc(ispec), &
                                                       m_b, temp_loc(jspec), Z_b, dens_loc(jspec), &
                                                       temp_loc(1), dens_loc(1)) ! species(1) is electrons
    END DO

    IF (v_particle > 1.0E-6_8) THEN
      cmul_out = nu_total_v / v_particle
    ELSE
      cmul_out = 0.0_8
    END IF

  CONTAINS
  
    ! nu for collisions between species a and b 
    FUNCTION nu_ab_v_multispecies(m_a, T_a_keV, v_a, Z_a, n_a_1e19, &
                                  m_b, T_b_keV, Z_b, n_b_1e19, &
                                  Te_keV_bg, ne_1e19_bg) RESULT(nu_out)
      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN) :: m_a, T_a_keV, v_a, Z_a, n_a_1e19
      REAL(KIND=8), INTENT(IN) :: m_b, T_b_keV, Z_b, n_b_1e19
      REAL(KIND=8), INTENT(IN) :: Te_keV_bg, ne_1e19_bg
      REAL(KIND=8) :: nu_out, n_b_si, ln_lambda, v_th_b, xb, G_xb, Phi_xb
      REAL(KIND=8) :: prefactor, sum_term

      IF (v_a < 1.0E-6_8) THEN; nu_out = 0.0_8; RETURN; END IF
      n_b_si = n_b_1e19 * 1.0E19_8
      IF (T_b_keV < 1.0E-3_8) THEN; nu_out = 0.0_8; RETURN; END IF
      
      v_th_b = SQRT(2.0_8 * (T_b_keV * 1000.0_8 * ELEMENTARY_CHARGE_C) / m_b)
      
      ! Coulomb logarithm 
      ln_lambda = calculate_ln_lambda_multispecies(m_a, T_a_keV, n_a_1e19, Z_a, &
                                                   m_b, T_b_keV, n_b_1e19, Z_b, &
                                                   Te_keV_bg, ne_1e19_bg)
      
      ! Chandrasekhar function G(x) and error function Phi(x)
      xb = v_a / v_th_b
      IF (xb < 1.0E-4_8) THEN
        Phi_xb = (2.0_8/SQRT(MY_PI))*(xb - xb**3/3.0_8)
        G_xb = (2.0_8/(3.0_8*SQRT(MY_PI)))*xb
      ELSE
        Phi_xb = ERF(xb)
        G_xb = (Phi_xb - xb*(2.0_8/SQRT(MY_PI))*EXP(-xb**2))/(2.0_8*xb**2)
      END IF

      prefactor = (ELEMENTARY_CHARGE_C**4 * Z_a**2) / &
                  ((4.0_8 * MY_PI * PERMITTIVITY**2) * (m_a**2 * v_a**4))
      sum_term = Z_b**2 * n_b_si * (Phi_xb - G_xb) * ln_lambda
      nu_out = v_a * prefactor * sum_term
      
    END FUNCTION nu_ab_v_multispecies

    ! Coulomb logarithm
    FUNCTION calculate_ln_lambda_multispecies(m_a, T_a_keV, n_a_1e19, Z_a, &
                                              m_b, T_b_keV, n_b_1e19, Z_b, &
                                              Te_keV_bg, ne_1e19_bg) RESULT(ln_lambda_out)
      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN) :: m_a, T_a_keV, n_a_1e19, Z_a
      REAL(KIND=8), INTENT(IN) :: m_b, T_b_keV, n_b_1e19, Z_b
      REAL(KIND=8), INTENT(IN) :: Te_keV_bg, ne_1e19_bg
      REAL(KIND=8) :: ln_lambda_out
      REAL(KIND=8) :: Te_eV, ne_m3, Ta_eV, Tb_eV, na_m3, nb_m3
      REAL(KIND=8) :: term1, term_sqrt

      ! Collision involving electrons (Z_a < 0 or Z_b < 0)
      IF (Z_a < 0.0_8 .OR. Z_b < 0.0_8) THEN
        Te_eV = Te_keV_bg * 1000.0_8
        ne_m3 = ne_1e19_bg * 1.0E19_8

        IF (Te_eV > 1.0E-6_8 .AND. ne_m3 > 1.0E-6_8) THEN
          ln_lambda_out = 32.2_8 + 1.15_8 * LOG10(Te_eV**2 / ne_m3)
        ELSE
          ln_lambda_out = 17.0_8
        END IF

      ! ion-ion collisions (including impurity-impurity)
      ELSE
        Ta_eV = T_a_keV * 1000.0_8
        Tb_eV = T_b_keV * 1000.0_8
        na_m3 = n_a_1e19 * 1.0E19_8
        nb_m3 = n_b_1e19 * 1.0E19_8

        IF (Ta_eV > 0.0_8 .AND. Tb_eV > 0.0_8) THEN
          term1 = Z_a * Z_b * (m_a + m_b) / (m_a * Tb_eV + m_b * Ta_eV)
          term_sqrt = na_m3 * Z_a**2 / Ta_eV + nb_m3 * Z_b**2 / Tb_eV
          IF (term_sqrt > 0.0_8 .AND. term1 > 0.0_8) THEN
            ln_lambda_out = 30.3_8 - LOG(term1 * SQRT(term_sqrt))
          ELSE
            ln_lambda_out = 17.0_8 
          END IF
        ELSE
          ln_lambda_out = 17.0_8 
        END IF
      END IF

    END FUNCTION calculate_ln_lambda_multispecies

  END FUNCTION calculate_dkes_cmul_multispecies


  SUBROUTINE apply_momentum_correction_multispecies(N_spec, dens_loc, temp_loc, Z_spec, A_spec, &
       L11, L12, L22, M_00, M_01, M_02, M_11, M_12, M_22, N_00, N_01, N_02, N_11, N_12, N_22, &
       X1_SN, X2_SN, A1_DKES, A2_DKES, b2_norm_val, Er_loc_Vm, drho_dr, r_lookup, &
       b00_in, flux_part_out, flux_heat_out, jbs_array_out, sigma_parallel_out)
       
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N_spec
    REAL(KIND=8), DIMENSION(N_spec), INTENT(IN) :: dens_loc, temp_loc, Z_spec, A_spec
    REAL(KIND=8), DIMENSION(N_spec), INTENT(IN) :: L11, L12, L22
    REAL(KIND=8), DIMENSION(N_spec), INTENT(IN) :: M_00, M_01, M_02, M_11, M_12, M_22
    REAL(KIND=8), DIMENSION(N_spec), INTENT(IN) :: N_00, N_01, N_02, N_11, N_12, N_22
    REAL(KIND=8), DIMENSION(N_spec), INTENT(IN) :: X1_SN, X2_SN, A1_DKES, A2_DKES
    REAL(KIND=8), INTENT(IN) :: b2_norm_val, Er_loc_Vm, drho_dr, r_lookup, b00_in
    
    REAL(KIND=8), DIMENSION(N_spec), INTENT(OUT) :: flux_part_out, flux_heat_out, jbs_array_out
    REAL(KIND=8), INTENT(OUT) :: sigma_parallel_out
    
    REAL(KIND=8), ALLOCATABLE :: Mat(:,:), RHS(:,:), Sol(:,:)
    
    INTEGER :: ierr, a, b
    REAL(KIND=8) :: n_si(N_spec), T_joules(N_spec), m_kg(N_spec), vT(N_spec), q_si(N_spec)
    REAL(KIND=8) :: nu_hat_ab, tau_ab, ln_lambda, x2, x4, x6, x8, yplus1, invsq
    REAL(KIND=8) :: term1, term_sqrt
    REAL(KIND=8) :: x_abi(N_spec, N_spec), TaoverTb_abi(N_spec, N_spec), y_ab(N_spec, N_spec)
    REAL(KIND=8) :: M_mat(N_spec, N_spec, 3, 3), N_mat(N_spec, N_spec, 3, 3)
    REAL(KIND=8) :: l_ab(N_spec, N_spec, 3, 3)
    REAL(KIND=8) :: M_sum(N_spec, 3, 3) 
    
    REAL(KIND=8), PARAMETER :: eps0 = 8.8541878E-12_8
    REAL(KIND=8), PARAMETER :: e_ch = 1.6021766E-19_8
    REAL(KIND=8), PARAMETER :: m_p = 1.6606E-27_8 
    REAL(KIND=8), PARAMETER :: m_e = 9.1093837E-31_8
    
    ALLOCATE(Mat(3*N_spec, 3*N_spec), RHS(3*N_spec, 2), Sol(3*N_spec, 2))
    Mat = 0.0_8; RHS = 0.0_8; Sol = 0.0_8

    ! Uncorr radial fluxes
    DO a = 1, N_spec
        n_si(a) = dens_loc(a) * 1.0E19_8
        T_joules(a) = temp_loc(a) * 1000.0_8 * e_ch
        q_si(a) = Z_spec(a) * e_ch
        IF (Z_spec(a) < 0.0_8) THEN
            m_kg(a) = m_e
        ELSE
            m_kg(a) = (A_spec(a) * m_p) - (Z_spec(a) * m_e)
        END IF
        vT(a) = SQRT(2.0_8 * T_joules(a) / m_kg(a))
        
        flux_part_out(a) = -n_si(a) * (L11(a)*A1_DKES(a) + L12(a)*A2_DKES(a))
        flux_heat_out(a) = -n_si(a) * T_joules(a) * (L12(a)*A1_DKES(a) + L22(a)*A2_DKES(a))
    END DO
    
    DO a = 1, N_spec
        DO b = 1, N_spec
            x_abi(a,b) = vT(b) / vT(a)
            TaoverTb_abi(a,b) = T_joules(a) / T_joules(b)
            y_ab(a,b) = m_kg(a) / m_kg(b)
        END DO
    END DO

    ! Classical friction matrix (Hirshman formulation )
    M_mat = 0.0_8; N_mat = 0.0_8
    DO a = 1, N_spec
        DO b = 1, N_spec
            x2 = x_abi(a,b)**2; x4 = x2**2; x6 = x4*x2; x8 = x4*x4
            yplus1 = 1.0_8 + y_ab(a,b)
            invsq = 1.0_8 / SQRT(1.0_8 + x2)
            
            ! 0-0 terms (particle flow drag)
            M_mat(a,b,1,1) = -yplus1 * invsq**3
            N_mat(a,b,1,1) = -M_mat(a,b,1,1)
            
            ! 0-1 and 1-0 terms (thermal friction)
            M_mat(a,b,1,2) = -1.5_8 * yplus1 * invsq**5
            M_mat(a,b,2,1) = M_mat(a,b,1,2)
            N_mat(a,b,2,1) = 1.5_8 * yplus1 * invsq**5
	    N_mat(a,b,1,2) = (1.5_8 * (1.0_8 + y_ab(b,a)) * &
                 (1.0_8 / SQRT(1.0_8 + x_abi(b,a)**2))**5) / TaoverTb_abi(b,a) * x_abi(b,a)
            
            ! 1-1 terms (heat flow drag)
            M_mat(a,b,2,2) = -(3.25_8 + 4.0_8*x2 + 7.5_8*x4) * invsq**5
            N_mat(a,b,2,2) = 6.75_8 * TaoverTb_abi(a,b) * x2 * invsq**5
            
            ! 0-2 and 2-0 terms (higher order)
            M_mat(a,b,1,3) = -(15.0_8 / 8.0_8) * yplus1 * invsq**7
            M_mat(a,b,3,1) = M_mat(a,b,1,3)
            N_mat(a,b,3,1) = -M_mat(a,b,1,3)
            N_mat(a,b,1,3) = -M_mat(a,b,1,3) * x4 
            
            ! 1-2 and 2-1 terms
            M_mat(a,b,2,3) = -(69.0_8 / 16.0_8 + 6.0_8 * x2 + (63.0_8 / 4.0_8) * x4) * invsq**7
            M_mat(a,b,3,2) = M_mat(a,b,2,3)
            N_mat(a,b,2,3) = (225.0_8 / 16.0_8) * TaoverTb_abi(a,b) * x4 * invsq**7
            N_mat(a,b,3,2) = N_mat(a,b,2,3) / TaoverTb_abi(a,b) / x2
            
            ! 2-2 terms
            M_mat(a,b,3,3) = -(433.0_8 / 64.0_8 + 17.0_8 * x2 + (459.0_8 / 8.0_8) * x4 + 28.0_8 * x6 + (175.0_8 / 8.0_8) * x8) * invsq**9
            N_mat(a,b,3,3) = (2625.0_8 / 64.0_8) * TaoverTb_abi(a,b) * x4 * invsq**9
        END DO
    END DO
    
    l_ab = 0.0_8
    M_sum = 0.0_8
    
    DO a = 1, N_spec
        DO b = 1, N_spec
            
            ! Calculate tau_ab (thermal deflection time)
            IF (Z_spec(a) < 0.0_8 .OR. Z_spec(b) < 0.0_8) THEN
                IF (temp_loc(1) > 1.0E-6_8 .AND. dens_loc(1) > 1.0E-6_8) THEN
                      ln_lambda = 32.2_8 + 1.15_8 * LOG10((temp_loc(1)*1000.0_8)**2 / (dens_loc(1)*1.0D19))
                ELSE
                    ln_lambda = 17.0_8
                END IF
            ELSE
                IF (temp_loc(a) > 0.0_8 .AND. temp_loc(b) > 0.0_8) THEN
                    term1 = Z_spec(a)*Z_spec(b)*(m_kg(a)+m_kg(b)) / (m_kg(a)*temp_loc(b)*1000.0_8 + m_kg(b)*temp_loc(a)*1000.0_8)
                    term_sqrt = n_si(a)*Z_spec(a)**2 / (temp_loc(a)*1000.0_8) + n_si(b)*Z_spec(b)**2 / (temp_loc(b)*1000.0_8)
                      IF (term_sqrt > 0.0_8 .AND. term1 > 0.0_8) THEN
                    ln_lambda = 30.3_8 - LOG(term1 * SQRT(term_sqrt))
                ELSE
                    ln_lambda = 17.0_8
                END IF
                  ELSE
                      ln_lambda = 17.0_8
                  END IF
            END IF
            
            nu_hat_ab = (n_si(b) * Z_spec(a)**2 * Z_spec(b)**2 * e_ch**4 * ln_lambda) / (4.0_8 * MY_PI * eps0**2 * m_kg(a)**2 * vT(a)**3)
            IF (nu_hat_ab > 1.0E-14_8) THEN; tau_ab = 0.75_8 * SQRT(MY_PI) / nu_hat_ab; ELSE; tau_ab = 1.0E14_8; END IF

            l_ab(a,b,1,1) = N_mat(a,b,1,1) * n_si(a)*m_kg(a) / tau_ab
            l_ab(a,b,1,2) = N_mat(a,b,1,2) * n_si(a)*m_kg(a) / tau_ab
            l_ab(a,b,1,3) = N_mat(a,b,1,3) * n_si(a)*m_kg(a) / tau_ab
            l_ab(a,b,2,1) = N_mat(a,b,2,1) * n_si(a)*m_kg(a) / tau_ab
            l_ab(a,b,2,2) = N_mat(a,b,2,2) * n_si(a)*m_kg(a) / tau_ab
            l_ab(a,b,2,3) = N_mat(a,b,2,3) * n_si(a)*m_kg(a) / tau_ab
            l_ab(a,b,3,1) = N_mat(a,b,3,1) * n_si(a)*m_kg(a) / tau_ab
            l_ab(a,b,3,2) = N_mat(a,b,3,2) * n_si(a)*m_kg(a) / tau_ab
            l_ab(a,b,3,3) = N_mat(a,b,3,3) * n_si(a)*m_kg(a) / tau_ab

            M_sum(a,1,1) = M_sum(a,1,1) + M_mat(a,b,1,1) * n_si(a)*m_kg(a) / tau_ab
            M_sum(a,1,2) = M_sum(a,1,2) + M_mat(a,b,1,2) * n_si(a)*m_kg(a) / tau_ab
            M_sum(a,1,3) = M_sum(a,1,3) + M_mat(a,b,1,3) * n_si(a)*m_kg(a) / tau_ab
            M_sum(a,2,1) = M_sum(a,2,1) + M_mat(a,b,2,1) * n_si(a)*m_kg(a) / tau_ab
            M_sum(a,2,2) = M_sum(a,2,2) + M_mat(a,b,2,2) * n_si(a)*m_kg(a) / tau_ab
            M_sum(a,2,3) = M_sum(a,2,3) + M_mat(a,b,2,3) * n_si(a)*m_kg(a) / tau_ab
            M_sum(a,3,1) = M_sum(a,3,1) + M_mat(a,b,3,1) * n_si(a)*m_kg(a) / tau_ab
            M_sum(a,3,2) = M_sum(a,3,2) + M_mat(a,b,3,2) * n_si(a)*m_kg(a) / tau_ab
            M_sum(a,3,3) = M_sum(a,3,3) + M_mat(a,b,3,3) * n_si(a)*m_kg(a) / tau_ab
        END DO
    END DO
    
    DO a = 1, N_spec
        l_ab(a,a,1,1) = l_ab(a,a,1,1) + M_sum(a,1,1)
        l_ab(a,a,1,2) = l_ab(a,a,1,2) + M_sum(a,1,2)
        l_ab(a,a,1,3) = l_ab(a,a,1,3) + M_sum(a,1,3)
        l_ab(a,a,2,1) = l_ab(a,a,2,1) + M_sum(a,2,1)
        l_ab(a,a,2,2) = l_ab(a,a,2,2) + M_sum(a,2,2)
        l_ab(a,a,2,3) = l_ab(a,a,2,3) + M_sum(a,2,3)
        l_ab(a,a,3,1) = l_ab(a,a,3,1) + M_sum(a,3,1)
        l_ab(a,a,3,2) = l_ab(a,a,3,2) + M_sum(a,3,2)
        l_ab(a,a,3,3) = l_ab(a,a,3,3) + M_sum(a,3,3)
    END DO
    
    ! parallel momentum balance: M_NC - L_class * Flows = - N_NC * X_SN
    
    ! left side Mat = M_NC - L_class
    DO a = 1, N_spec
        Mat(3*a-2, 3*a-2) = M_00(a)
        Mat(3*a-2, 3*a-1) = M_01(a)
        Mat(3*a-2, 3*a)   = M_02(a)
        Mat(3*a-1, 3*a-2) = M_01(a)
        Mat(3*a-1, 3*a-1) = M_11(a)
        Mat(3*a-1, 3*a)   = M_12(a)
        Mat(3*a,   3*a-2) = M_02(a)
        Mat(3*a,   3*a-1) = M_12(a)
        Mat(3*a,   3*a)   = M_22(a)
        
        DO b = 1, N_spec
	    
            Mat(3*a-2, 3*b-2) = Mat(3*a-2, 3*b-2) - l_ab(a,b,1,1) * (b00_in**2 * b2_norm_val)
            Mat(3*a-2, 3*b-1) = Mat(3*a-2, 3*b-1) - l_ab(a,b,1,2) * (b00_in**2 * b2_norm_val)
            Mat(3*a-2, 3*b)   = Mat(3*a-2, 3*b)   - l_ab(a,b,1,3) * (b00_in**2 * b2_norm_val)
            Mat(3*a-1, 3*b-2) = Mat(3*a-1, 3*b-2) - l_ab(a,b,2,1) * (b00_in**2 * b2_norm_val)
            Mat(3*a-1, 3*b-1) = Mat(3*a-1, 3*b-1) - l_ab(a,b,2,2) * (b00_in**2 * b2_norm_val)
            Mat(3*a-1, 3*b)   = Mat(3*a-1, 3*b)   - l_ab(a,b,2,3) * (b00_in**2 * b2_norm_val)
            Mat(3*a,   3*b-2) = Mat(3*a,   3*b-2) - l_ab(a,b,3,1) * (b00_in**2 * b2_norm_val)
            Mat(3*a,   3*b-1) = Mat(3*a,   3*b-1) - l_ab(a,b,3,2) * (b00_in**2 * b2_norm_val)
            Mat(3*a,   3*b)   = Mat(3*a,   3*b)   - l_ab(a,b,3,3) * (b00_in**2 * b2_norm_val)
        END DO
        
        ! right side RHS = - N_NC * X_SN
        RHS(3*a-2, 1) = -( N_00(a)*X1_SN(a) + N_01(a)*X2_SN(a) )
        RHS(3*a-1, 1) = -( N_01(a)*X1_SN(a) + N_11(a)*X2_SN(a) )
        RHS(3*a,   1) = -( N_02(a)*X1_SN(a) + N_12(a)*X2_SN(a) )

        ! Parallel conductivity source term 
        RHS(3*a-2, 2) = n_si(a) * q_si(a) * b00_in
        RHS(3*a-1, 2) = 0.0_8
        RHS(3*a,   2) = 0.0_8
    END DO

    CALL solve_NxM(3*N_spec, 2, Mat, RHS, Sol, ierr)
    
    sigma_parallel_out = 0.0_8
    IF (ierr == 0) THEN
        DO a = 1, N_spec
            jbs_array_out(a) = Sol(3*a-2, 1) * b00_in * b2_norm_val * q_si(a) * n_si(a)
            sigma_parallel_out = sigma_parallel_out + (Sol(3*a-2, 2) * b00_in * b2_norm_val * q_si(a) * n_si(a))
        END DO  
    ELSE
        jbs_array_out = 0.0_8
    END IF

    DEALLOCATE(Mat, RHS, Sol)
  END SUBROUTINE apply_momentum_correction_multispecies
  

  SUBROUTINE solve_NxM(N, M_cols, A, B, X, ierr)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, M_cols
    REAL(KIND=8), DIMENSION(N,N), INTENT(IN) :: A
    REAL(KIND=8), DIMENSION(N,M_cols), INTENT(IN) :: B
    REAL(KIND=8), DIMENSION(N,M_cols), INTENT(OUT) :: X
    INTEGER, INTENT(OUT) :: ierr  

    REAL(KIND=8), DIMENSION(N,N) :: M
    REAL(KIND=8), DIMENSION(N) :: tmp_row
    REAL(KIND=8) :: tmp_val, pivot_val, max_val, factor
    INTEGER :: i, j, k, p
    
    M = A; X = B; ierr = 0
    DO k = 1, N
       p = k; max_val = ABS(M(k,k))
       DO i = k+1, N
          IF (ABS(M(i,k)) > max_val) THEN
             max_val = ABS(M(i,k)); p = i
          END IF
       END DO
       
       IF (max_val < 1.0D-20) THEN; ierr = 1; RETURN; END IF
       
       IF (p /= k) THEN
          tmp_row = M(p,:); M(p,:) = M(k,:); M(k,:) = tmp_row
          DO j = 1, M_cols
             tmp_val = X(p,j); X(p,j) = X(k,j); X(k,j) = tmp_val
          END DO
       END IF
       
       pivot_val = M(k,k); M(k,:) = M(k,:) / pivot_val; X(k,:) = X(k,:) / pivot_val
       DO i = 1, N
          IF (i /= k) THEN
             factor = M(i,k); M(i,:) = M(i,:) - factor * M(k,:); X(i,:) = X(i,:) - factor * X(k,:)
          END IF
       END DO
    END DO
  END SUBROUTINE solve_NxM
  
  FUNCTION normalize_er_for_dkes(Er_Vm, b00_physical_at_r, v_particle, dr_dr_tilde) RESULT(efield_out)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: Er_Vm, b00_physical_at_r, v_particle, dr_dr_tilde
    REAL(KIND=8) :: efield_out
    
    efield_out = (Er_Vm * dr_dr_tilde) / (v_particle * b00_physical_at_r)

  END FUNCTION normalize_er_for_dkes
  
  
  ! Interpolate b2_norm radially for sigma parallel
  FUNCTION interpolate_b2_norm_only(r_lookup_target) RESULT(interpolated_value)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: r_lookup_target
    REAL(KIND=8) :: interpolated_value
    INTEGER :: startrind, Nradii, i, current_r_idx
    REAL(KIND=8), DIMENSION(4) :: rs_stencil, vals_at_r
    REAL(KIND=8) :: r0, r1, r2, r3, hr0, hr1, hr2, hr3
    REAL(KIND=8) :: r02, r12, r22, r03, r13, r23, ha, hb, hg
    REAL(KIND=8), PARAMETER :: r_tolerance = 1.0E-3_8
    INTEGER :: exact_match_idx

    IF (nr < 4) THEN
        startrind = 1
        Nradii = nr
        IF (nr == 0) THEN
             interpolated_value = 1.0_8
             RETURN
        END IF
    END IF

    exact_match_idx = 0
    DO i = 1, nr
        IF (ABS(r_lookup_target - dkes_data(i)%r_val) < r_tolerance) THEN
            exact_match_idx = i
            EXIT
        END IF
    END DO

    IF (exact_match_idx > 0) THEN
        startrind = exact_match_idx
        Nradii = 1
    ELSE IF (r_lookup_target < dkes_data(2)%r_val) THEN 
        startrind = 1
        Nradii = 3
    ELSE IF (r_lookup_target >= dkes_data(nr-1)%r_val) THEN 
        startrind = nr - 2
        Nradii = 3
    ELSE 
        startrind = 0
        Nradii = 4
        DO i = 2, nr - 2
            IF (r_lookup_target >= dkes_data(i)%r_val .AND. r_lookup_target < dkes_data(i+1)%r_val) THEN
                startrind = i - 1
                EXIT
            END IF
        END DO
        IF (startrind == 0) STOP 'Interpolation stencil error'
    END IF

    DO i = 1, Nradii
        current_r_idx = startrind + i - 1
        rs_stencil(i) = dkes_data(current_r_idx)%r_val
        vals_at_r(i) = dkes_data(current_r_idx)%b2_norm 
    END DO

    ! Radial interpolation/extrapolation
    IF (Nradii == 4) THEN
        r0 = rs_stencil(1); r1 = rs_stencil(2); r2 = rs_stencil(3); r3 = rs_stencil(4)
        hr0 = (r_lookup_target-r1)*(r_lookup_target-r2)*(r_lookup_target-r3) / ((r0-r1)*(r0-r2)*(r0-r3))
        hr1 = (r_lookup_target-r0)*(r_lookup_target-r2)*(r_lookup_target-r3) / ((r1-r0)*(r1-r2)*(r1-r3))
        hr2 = (r_lookup_target-r0)*(r_lookup_target-r1)*(r_lookup_target-r3) / ((r2-r0)*(r2-r1)*(r2-r3))
        hr3 = (r_lookup_target-r0)*(r_lookup_target-r1)*(r_lookup_target-r2) / ((r3-r0)*(r3-r1)*(r3-r2))
        interpolated_value = hr0*vals_at_r(1) + hr1*vals_at_r(2) + hr2*vals_at_r(3) + hr3*vals_at_r(4)
    ELSE IF (Nradii == 3) THEN
        r0 = rs_stencil(1); r1 = rs_stencil(2); r2 = rs_stencil(3)
        IF (r_lookup_target < dkes_data(2)%r_val .AND. nr >= 2) THEN ! Core: zero-derivative extrapolation
            r02 = r0**2; r12 = r1**2; r22 = r2**2
            r03 = r0**3; r13 = r1**3; r23 = r2**3
            
            ha = (((vals_at_r(3)-vals_at_r(2))/(r23-r13)) - ((vals_at_r(3)-vals_at_r(1))/(r23-r03))) / &
                 (((r22-r12)/(r23-r13)) - ((r22-r02)/(r23-r03)))
            hb = (((vals_at_r(3)-vals_at_r(2))/(r22-r12)) - ((vals_at_r(3)-vals_at_r(1))/(r22-r02))) / &
                 (((r23-r13)/(r22-r12)) - ((r23-r03)/(r22-r02)))
            hg = vals_at_r(1) - r02*ha - r03*hb
            
            interpolated_value = hg + (r_lookup_target**2) * ha + (r_lookup_target**3) * hb
        ELSE ! Edge or fallback
            hr0 = (r_lookup_target-r1)*(r_lookup_target-r2) / ((r0-r1)*(r0-r2))
            hr1 = (r_lookup_target-r0)*(r_lookup_target-r2) / ((r1-r0)*(r1-r2))
            hr2 = (r_lookup_target-r0)*(r_lookup_target-r1) / ((r2-r0)*(r2-r1))
            interpolated_value = hr0*vals_at_r(1) + hr1*vals_at_r(2) + hr2*vals_at_r(3)
        END IF
    ELSE IF (Nradii == 2) THEN
        r0 = rs_stencil(1); r1 = rs_stencil(2)
        IF (ABS(r1 - r0) > comparison_epsilon) THEN
             interpolated_value = vals_at_r(1) + (vals_at_r(2) - vals_at_r(1)) * (r_lookup_target - r0) / (r1 - r0)
        ELSE
             interpolated_value = vals_at_r(1)
        END IF
    ELSE IF (Nradii == 1) THEN 
        interpolated_value = vals_at_r(1)
    ELSE 
        interpolated_value = 1.0_8
    END IF

    IF (interpolated_value < 1.0E-6_8) THEN
!         WRITE(*,*) 'WARNING: Interpolated b2_norm is near zero or negative at r_lookup=', r_lookup_target, ' Value=', interpolated_value, '. Using fallback value 1.0.'
         interpolated_value = 1.0_8
    END IF
  END FUNCTION interpolate_b2_norm_only
  
  SUBROUTINE interpolate_dkes_coeffs(r, cmul, efield, g11_val, g31_val, g33_val)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: r, cmul, efield
    REAL(KIND=8), INTENT(OUT) :: g11_val, g31_val, g33_val
    INTEGER :: startrind, Nradii, i, current_r_idx
    REAL(KIND=8), DIMENSION(4) :: rs, x0_at_r, x1_at_r, x2_at_r
    REAL(KIND=8) :: r0, r1, r2, r3, hr0, hr1, hr2, hr3
    REAL(KIND=8) :: r02, r12, r22, r03, r13, r23, ha, hb, hg
    REAL(KIND=8) :: x0_val, x1_val, x2_val
    REAL(KIND=8) :: a_ratio
    REAL(KIND=8), PARAMETER :: r_tolerance = 1.0E-3_8
    INTEGER :: exact_match_idx

    exact_match_idx = 0
    DO i = 1, nr
        IF (ABS(r - dkes_data(i)%r_val) < r_tolerance) THEN
            exact_match_idx = i
            EXIT
        END IF
    END DO

    IF (exact_match_idx > 0) THEN
        startrind = exact_match_idx
        Nradii = 1
    ELSE IF (r < dkes_data(2)%r_val) THEN 
        startrind = 1
        Nradii = 3
    ELSE IF (r >= dkes_data(nr-1)%r_val) THEN 
        startrind = nr - 2
        Nradii = 3
    ELSE 
        startrind = 0 
        Nradii = 4
        DO i = 2, nr - 2
            IF (r >= dkes_data(i)%r_val .AND. r < dkes_data(i+1)%r_val) THEN
                startrind = i - 1
                EXIT
            END IF
        END DO
        IF (startrind == 0) STOP 'Interpolation stencil error'
    END IF

    ! Get coefficients at the DKES grid points
    DO i = 1, Nradii
        current_r_idx = startrind + i - 1
        rs(i) = dkes_data(current_r_idx)%r_val
        CALL interpolate_on_block(dkes_data(current_r_idx), cmul, efield, r, x0_at_r(i), x1_at_r(i), x2_at_r(i))
    END DO

    ! INTERPOLATION / EXTRAPOLATION 
    IF (Nradii == 4) THEN 
        ! Interior: 4-point Lagrange 
        r0 = rs(1); r1 = rs(2); r2 = rs(3); r3 = rs(4)
        hr0 = (r-r1)*(r-r2)*(r-r3) / ((r0-r1)*(r0-r2)*(r0-r3))
        hr1 = (r-r0)*(r-r2)*(r-r3) / ((r1-r0)*(r1-r2)*(r1-r3))
        hr2 = (r-r0)*(r-r1)*(r-r3) / ((r2-r0)*(r2-r1)*(r2-r3))
        hr3 = (r-r0)*(r-r1)*(r-r2) / ((r3-r0)*(r3-r1)*(r3-r2))
        
        x0_val = hr0*x0_at_r(1) + hr1*x0_at_r(2) + hr2*x0_at_r(3) + hr3*x0_at_r(4)
        x1_val = hr0*x1_at_r(1) + hr1*x1_at_r(2) + hr2*x1_at_r(3) + hr3*x1_at_r(4)
        x2_val = hr0*x2_at_r(1) + hr1*x2_at_r(2) + hr2*x2_at_r(3) + hr3*x2_at_r(4)

    ELSE IF (Nradii == 3) THEN
        r0 = rs(1); r1 = rs(2); r2 = rs(3)
        IF (r < dkes_data(2)%r_val) THEN 
            ! Core: f(r) = hg + ha*r^2 + hb*r^3
            r02 = r0**2; r12 = r1**2; r22 = r2**2
            r03 = r0**3; r13 = r1**3; r23 = r2**3
            
            ! Fit x0
            ha = (((x0_at_r(3)-x0_at_r(2))/(r23-r13)) - ((x0_at_r(3)-x0_at_r(1))/(r23-r03))) / &
                 (((r22-r12)/(r23-r13)) - ((r22-r02)/(r23-r03)))
            hb = (((x0_at_r(3)-x0_at_r(2))/(r22-r12)) - ((x0_at_r(3)-x0_at_r(1))/(r22-r02))) / &
                 (((r23-r13)/(r22-r12)) - ((r23-r03)/(r22-r02)))
            hg = x0_at_r(1) - r02*ha - r03*hb
            x0_val = hg + r**2 * ha + r**3 * hb
        
            ! Fit x1
            ha = (((x1_at_r(3)-x1_at_r(2))/(r23-r13)) - ((x1_at_r(3)-x1_at_r(1))/(r23-r03))) / &
                 (((r22-r12)/(r23-r13)) - ((r22-r02)/(r23-r03)))
            hb = (((x1_at_r(3)-x1_at_r(2))/(r22-r12)) - ((x1_at_r(3)-x1_at_r(1))/(r22-r02))) / &
                 (((r23-r13)/(r22-r12)) - ((r23-r03)/(r22-r02)))
            hg = x1_at_r(1) - r02*ha - r03*hb
            x1_val = hg + r**2 * ha + r**3 * hb
        
            ! Fit x2
            ha = (((x2_at_r(3)-x2_at_r(2))/(r23-r13)) - ((x2_at_r(3)-x2_at_r(1))/(r23-r03))) / &
                 (((r22-r12)/(r23-r13)) - ((r22-r02)/(r23-r03)))
            hb = (((x2_at_r(3)-x2_at_r(2))/(r22-r12)) - ((x2_at_r(3)-x2_at_r(1))/(r22-r02))) / &
                 (((r23-r13)/(r22-r12)) - ((r23-r03)/(r22-r02)))
            hg = x2_at_r(1) - r02*ha - r03*hb
            x2_val = hg + r**2 * ha + r**3 * hb
            
        ELSE 
            ! Edge: 3-point Lagrange extrapolation
            hr0 = (r-r1)*(r-r2) / ((r0-r1)*(r0-r2))
            hr1 = (r-r0)*(r-r2) / ((r1-r0)*(r1-r2))
            hr2 = (r-r0)*(r-r1) / ((r2-r0)*(r2-r1))
            
            x0_val = hr0*x0_at_r(1) + hr1*x0_at_r(2) + hr2*x0_at_r(3)
            x1_val = hr0*x1_at_r(1) + hr1*x1_at_r(2) + hr2*x1_at_r(3)
            x2_val = hr0*x2_at_r(1) + hr1*x2_at_r(2) + hr2*x2_at_r(3)
        END IF
        
    ELSE IF (Nradii == 1) THEN
        x0_val = x0_at_r(1)
        x1_val = x1_at_r(1)
        x2_val = x2_at_r(1)
    END IF

    a_ratio = module_ABC / minor_radius_W7AS_value

    g11_val = -(10.0_8**x0_val) * (a_ratio**2)
    g31_val = x1_val * (1.0_8 / (1.0_8 + cmul**2)) * a_ratio
    g33_val = -x2_val / cmul

  END SUBROUTINE interpolate_dkes_coeffs


  SUBROUTINE interpolate_on_block(block, target_cmul, target_efield, r_lookup, x0_out, x1_out, x2_out)
    IMPLICIT NONE
    TYPE(dkes_radial_block), INTENT(IN) :: block
    REAL(KIND=8), INTENT(IN) :: target_cmul, target_efield, r_lookup 
    REAL(KIND=8), INTENT(OUT) :: x0_out, x1_out, x2_out

    INTEGER, PARAMETER :: NUM_SPLINE_PTS = 4
    INTEGER :: i, n_unique
    INTEGER :: spline_indices(NUM_SPLINE_PTS)
    REAL(KIND=8) :: log_cmul_pts(NUM_SPLINE_PTS)
    REAL(KIND=8) :: x0_pts(NUM_SPLINE_PTS), x1_pts(NUM_SPLINE_PTS), x2_pts(NUM_SPLINE_PTS)
    REAL(KIND=8) :: cmul_val, lgcmul
    REAL(KIND=8) :: tok_g11, tok_g31, tok_g33
    REAL(KIND=8) :: g11_e0, g13_e0, g33_e0, cfit_g11, dum
    REAL(KIND=8) :: g11_corr, g31_corr, g33_corr, g11_nu
    REAL(KIND=8) :: target_efield_eff

    n_unique = block%num_unique_cmuls
    lgcmul = LOG10(target_cmul)
    
    IF (target_cmul < block%unique_cmuls(1)) THEN
    
      ! resonance broadened effective electric field
      target_efield_eff = target_efield * (block%unique_cmuls(1) / target_cmul)**(1.0_8 / 3.0_8)
        
      IF (ABS(target_efield) > block%EFIELD_upper .AND. block%EFIELD_upper > 0.0_8) THEN

          CALL evaluate_cached_spline(block, 1, target_efield_eff, r_lookup, tok_g11, dum, dum)
          g11_nu = 10.0_8**(tok_g11 + LOG10(target_cmul) - LOG10(block%unique_cmuls(1)))

          CALL evaluate_cached_spline(block, 1, target_efield, r_lookup, dum, x1_out, x2_out)

          x0_out = LOG10(MAX(1.0E-20_8, g11_nu))

      ELSE
          CALL calculate_low_coll_fit(block, target_cmul, target_efield, cfit_g11)
          
          CALL evaluate_cached_spline(block, 1, block%EFIELD_upper, r_lookup, tok_g11, dum, dum)
          g11_nu = 10.0_8**(tok_g11 + LOG10(target_cmul) - LOG10(block%unique_cmuls(1)))

          x0_out = LOG10(MAX(1.0E-20_8, g11_nu + ABS(cfit_g11)))
          
          CALL evaluate_cached_spline(block, 1, target_efield, r_lookup, dum, x1_out, x2_out)
    END IF

    ELSE IF (target_cmul > block%unique_cmuls(n_unique)) THEN
        
      CALL calculate_tokamak_fit(block%r_val, target_cmul, dkes_data(nr)%R00, block%iota, block%kn, &
                                 tok_g11, tok_g31, tok_g33)
    
      CALL calculate_high_efield_correction(block, target_efield, target_cmul, &
                                            tok_g11, tok_g31, tok_g33, &
                                            g11_corr, g31_corr, g33_corr)
            
      x0_out = LOG10(MAX(1.0E-20_8, ABS(g11_corr)))
      x1_out = g31_corr / (1.0_8 / (1.0_8 + target_cmul**2))
      x2_out = ABS(g33_corr) * target_cmul
            

        ELSE
        
        CALL find_nearest_spline_indices(block%unique_cmuls, n_unique, target_cmul, NUM_SPLINE_PTS, spline_indices)
    
    DO i = 1, NUM_SPLINE_PTS
        cmul_val = block%unique_cmuls(spline_indices(i))
        log_cmul_pts(i) = LOG10(cmul_val)
            CALL evaluate_cached_spline(block, spline_indices(i), target_efield, r_lookup, x0_pts(i), x1_pts(i), x2_pts(i))
    END DO

        IF (spline_indices(1) == 1) THEN
            CALL hermite_cmul_interp(log_cmul_pts, x0_pts, lgcmul, x0_out, leftBC_dydx_in=-0.5_8)
        ELSE
            CALL hermite_cmul_interp(log_cmul_pts, x0_pts, lgcmul, x0_out)
        END IF
        CALL hermite_cmul_interp(log_cmul_pts, x1_pts, lgcmul, x1_out)
        CALL hermite_cmul_interp(log_cmul_pts, x2_pts, lgcmul, x2_out)
    END IF
  END SUBROUTINE interpolate_on_block


  SUBROUTINE evaluate_cached_spline(block, cmul_idx, efield_val, r_lookup, x0_out, x1_out, x2_out)
    IMPLICIT NONE
    TYPE(dkes_radial_block), INTENT(IN) :: block
    INTEGER, INTENT(IN) :: cmul_idx
    REAL(KIND=8), INTENT(IN) :: efield_val, r_lookup 
    REAL(KIND=8), INTENT(OUT) :: x0_out, x1_out, x2_out
    
    INTEGER :: count
    REAL(KIND=8) :: target_x_mapped, a_ratio, scale_factor

    a_ratio = module_ABC / minor_radius_W7AS_value
    scale_factor = 1.0E-7_8

    count = block%slice_counts(cmul_idx)

    IF (r_lookup > 1.0E-6_8) THEN
        target_x_mapped = ASINH(ABS(efield_val) / (r_lookup * scale_factor))
    ELSE
        target_x_mapped = 0.0_8
    END IF
    
    IF (count >= 3) THEN
        x0_out = spline_interp(target_x_mapped, block%spline_x(cmul_idx, 1:count), &
                               block%spline_x0(cmul_idx, 1:count), block%spline_y2_0(cmul_idx, 1:count), count)
        x1_out = spline_interp(target_x_mapped, block%spline_x(cmul_idx, 1:count), &
                               block%spline_x1(cmul_idx, 1:count), block%spline_y2_1(cmul_idx, 1:count), count)
        x2_out = spline_interp(target_x_mapped, block%spline_x(cmul_idx, 1:count), &
                               block%spline_x2(cmul_idx, 1:count), block%spline_y2_2(cmul_idx, 1:count), count)
    ELSE
        x0_out = block%spline_x0(cmul_idx, 1)
        x1_out = block%spline_x1(cmul_idx, 1)
        x2_out = block%spline_x2(cmul_idx, 1)
    END IF
  END SUBROUTINE evaluate_cached_spline
            
            
  SUBROUTINE hermite_interp_1d_array(x_arr, y_arr, n, x_target, y_out)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    REAL(KIND=8), INTENT(IN) :: x_arr(n), y_arr(n), x_target
    REAL(KIND=8), INTENT(OUT) :: y_out
            
    INTEGER :: i1, i2
    REAL(KIND=8) :: w, dgl, dgu, dx21, t, ha, hb
    REAL(KIND=8) :: xa, xb, xc, ya, yb, yc, c0, c1, c2
    REAL(KIND=8) :: dgmin, dgmax
    REAL(KIND=8), PARAMETER :: mix_param = 0.1_8
    REAL(KIND=8), PARAMETER :: mix_trans = 0.1_8

    IF (x_target <= x_arr(1)) THEN
        y_out = y_arr(1)
        RETURN
    ELSE IF (x_target >= x_arr(n)) THEN
        y_out = y_arr(n)
        RETURN
            END IF
    
    CALL find_indices_and_weight_1d(x_arr, n, x_target, i1, i2, w)

    IF (n < 3) THEN
        y_out = linear_interp_1d(y_arr(i1), y_arr(i2), w)
        RETURN
    END IF

    dx21 = x_arr(i2) - x_arr(i1)
    t = (x_target - x_arr(i1)) / dx21

    IF (i1 == 1) THEN
        dgl = 0.0_8
    ELSE
        xa = x_arr(i1-1); xb = x_arr(i1); xc = x_arr(i1+1)
        ya = y_arr(i1-1); yb = y_arr(i1); yc = y_arr(i1+1)
        c0 = (xb-xc)/((xa-xb)*(xa-xc))
        c1 = (2.0_8*xb-xa-xc)/((xb-xa)*(xb-xc))
        c2 = (xb-xa)/((xc-xa)*(xc-xb))
        dgl = c0*ya + c1*yb + c2*yc
    END IF

    IF (i2 == n) THEN
        dgu = 3.0_8 * (y_arr(n) - y_arr(n-1)) / dx21 - 0.5_8 * dgl
    ELSE
        xa = x_arr(i2-1); xb = x_arr(i2); xc = x_arr(i2+1)
        ya = y_arr(i2-1); yb = y_arr(i2); yc = y_arr(i2+1)
        c0 = (xb-xc)/((xa-xb)*(xa-xc))
        c1 = (2.0_8*xb-xa-xc)/((xb-xa)*(xb-xc))
        c2 = (xb-xa)/((xc-xa)*(xc-xb))
        dgu = c0*ya + c1*yb + c2*yc
    END IF

    IF (dgl * dgu < 0.0_8) THEN
        dgmin = MIN(ABS(dgl), ABS(dgu))
        dgmax = MAX(ABS(dgl), ABS(dgu))
        IF (dgmax > 0.0_8) THEN
            IF (dgmin >= mix_trans * dgmax) THEN
                dgl = mix_param * dgl
                dgu = mix_param * dgu
            ELSE
                dgl = dgl * (1.0_8 - (1.0_8 - mix_param) * dgmin / mix_trans / dgmax)
                dgu = dgu * (1.0_8 - (1.0_8 - mix_param) * dgmin / mix_trans / dgmax)
            END IF
        END IF
    END IF

    ha = 3.0_8*(y_arr(i2)-y_arr(i1)) - (2.0_8*dgl + dgu)*dx21
    hb = -2.0_8*(y_arr(i2)-y_arr(i1)) + (dgl + dgu)*dx21

    y_out = y_arr(i1) + t*(dgl*dx21 + t*(ha + t*hb))
        
  END SUBROUTINE hermite_interp_1d_array


  SUBROUTINE hermite_cmul_interp(x, y, x_target, y_out, leftBC_dydx_in)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: x(4), y(4), x_target
    REAL(KIND=8), INTENT(OUT) :: y_out
    REAL(KIND=8), INTENT(IN), OPTIONAL :: leftBC_dydx_in

    REAL(KIND=8) :: d(4), dx_int, t, ha, hb
    REAL(KIND=8) :: mix_param = 0.2_8, mix_trans = 0.1_8
    REAL(KIND=8) :: Dxla, Dxlb, Dxlmin, Dxlmax
    INTEGER :: i1, i2

    IF (PRESENT(leftBC_dydx_in)) THEN
        d(1) = leftBC_dydx_in
    ELSE
        d(1) = ((2.0_8*x(1)-x(2)-x(3))/((x(1)-x(2))*(x(1)-x(3)))) * y(1) + &
               ((x(1)-x(3))/((x(2)-x(1))*(x(2)-x(3)))) * y(2) + &
               ((x(1)-x(2))/((x(3)-x(1))*(x(3)-x(2)))) * y(3)
    END IF

    d(2) = ((x(2)-x(3))/((x(1)-x(2))*(x(1)-x(3)))) * y(1) + &
           ((2.0_8*x(2)-x(1)-x(3))/((x(2)-x(1))*(x(2)-x(3)))) * y(2) + &
           ((x(2)-x(1))/((x(3)-x(1))*(x(3)-x(2)))) * y(3)
    Dxla = y(2) - y(1); Dxlb = y(2) - y(3)
    IF (Dxla * Dxlb > 0.0_8) THEN
        Dxlmin = MIN(ABS(Dxla), ABS(Dxlb)); Dxlmax = MAX(ABS(Dxla), ABS(Dxlb))
        IF (Dxlmin >= mix_trans * Dxlmax) THEN
            d(2) = mix_param * d(2)
        ELSE
            d(2) = d(2) * (1.0_8 - (1.0_8 - mix_param) * Dxlmin / (mix_trans * Dxlmax))
        END IF
    END IF

    d(3) = ((x(3)-x(4))/((x(2)-x(3))*(x(2)-x(4)))) * y(2) + &
           ((2.0_8*x(3)-x(2)-x(4))/((x(3)-x(2))*(x(3)-x(4)))) * y(3) + &
           ((x(3)-x(2))/((x(4)-x(2))*(x(4)-x(3)))) * y(4)
    Dxla = y(3) - y(2); Dxlb = y(3) - y(4)
    IF (Dxla * Dxlb > 0.0_8) THEN
        Dxlmin = MIN(ABS(Dxla), ABS(Dxlb)); Dxlmax = MAX(ABS(Dxla), ABS(Dxlb))
        IF (Dxlmin >= mix_trans * Dxlmax) THEN
            d(3) = mix_param * d(3)
        ELSE
            d(3) = d(3) * (1.0_8 - (1.0_8 - mix_param) * Dxlmin / (mix_trans * Dxlmax))
        END IF
    END IF

    d(4) = ((x(4)-x(3))/((x(2)-x(4))*(x(2)-x(3)))) * y(2) + ((2.0_8*x(4)-x(2)-x(3))/((x(4)-x(2))*(x(4)-x(3)))) * y(3) + &
           ((x(4)-x(2))/((x(3)-x(4))*(x(3)-x(2)))) * y(4)

    IF (x_target <= x(2)) THEN
        i1 = 1; i2 = 2
    ELSE IF (x_target >= x(3)) THEN
        i1 = 3; i2 = 4
        ELSE
        i1 = 2; i2 = 3
    END IF

    dx_int = x(i2) - x(i1)
    t = (x_target - x(i1)) / dx_int
    ha = 3.0_8*(y(i2)-y(i1)) - (2.0_8*d(i1) + d(i2))*dx_int
    hb = -2.0_8*(y(i2)-y(i1)) + (d(i1) + d(i2))*dx_int

    y_out = y(i1) + t*(d(i1)*dx_int + t*(ha + t*hb))

  END SUBROUTINE hermite_cmul_interp
  
  SUBROUTINE build_spline_cache(block)
    IMPLICIT NONE
    TYPE(dkes_radial_block), INTENT(INOUT) :: block
    INTEGER :: j, i, k, count_raw, valid_count, count
    INTEGER :: ii, jj
    REAL(KIND=8) :: cmul_val, efield_norm, g13_norm
    REAL(KIND=8) :: a_ratio, scale_factor, temp_val
    REAL(KIND=8), ALLOCATABLE :: temp_x(:), temp_x0(:), temp_x1(:), temp_x2(:)
    REAL(KIND=8), ALLOCATABLE :: x_slice(:), x0_slice(:), x1_slice(:), x2_slice(:)
    REAL(KIND=8), ALLOCATABLE :: y2_0(:), y2_1(:), y2_2(:)
    REAL(KIND=8), PARAMETER :: E_NORM_MIN = 1.0E-10_8

    a_ratio = module_ABC / minor_radius_W7AS_value
    scale_factor = 1.0E-7_8 
    
    ALLOCATE(block%slice_counts(block%num_unique_cmuls))
    ALLOCATE(block%spline_x(block%num_unique_cmuls, 200))
    ALLOCATE(block%spline_x0(block%num_unique_cmuls, 200))
    ALLOCATE(block%spline_x1(block%num_unique_cmuls, 200))
    ALLOCATE(block%spline_x2(block%num_unique_cmuls, 200))
    ALLOCATE(block%spline_y2_0(block%num_unique_cmuls, 200))
    ALLOCATE(block%spline_y2_1(block%num_unique_cmuls, 200))
    ALLOCATE(block%spline_y2_2(block%num_unique_cmuls, 200))

    DO j = 1, block%num_unique_cmuls
       cmul_val = block%unique_cmuls(j)
       g13_norm = 1.0_8 / (1.0_8 + cmul_val**2)

       count_raw = 0
       DO i = 1, block%num_points
           IF (ABS(block%cmul_list(i) - cmul_val) < comparison_epsilon) count_raw = count_raw + 1
       END DO

       ALLOCATE(temp_x(count_raw), temp_x0(count_raw), temp_x1(count_raw), temp_x2(count_raw))

       k = 1
       DO i = 1, block%num_points
           IF (ABS(block%cmul_list(i) - cmul_val) < comparison_epsilon) THEN
               efield_norm = ABS(block%efield_list(i)) / (block%r_val * scale_factor)
               efield_norm = MAX(E_NORM_MIN, efield_norm)
               temp_x(k) = ASINH(efield_norm)
               temp_x0(k) = LOG10(MAX(1.0E-20_8, ABS(block%g11_list(i))))
               temp_x1(k) = block%g31_list(i) / g13_norm
               temp_x2(k) = ABS(block%g33_list(i)) * cmul_val
               k = k + 1
           END IF
       END DO

       DO ii = 1, count_raw - 1
           DO jj = 1, count_raw - ii
               IF (temp_x(jj) > temp_x(jj+1)) THEN
                   temp_val = temp_x(jj); temp_x(jj) = temp_x(jj+1); temp_x(jj+1) = temp_val
                   temp_val = temp_x0(jj); temp_x0(jj) = temp_x0(jj+1); temp_x0(jj+1) = temp_val
                   temp_val = temp_x1(jj); temp_x1(jj) = temp_x1(jj+1); temp_x1(jj+1) = temp_val
                   temp_val = temp_x2(jj); temp_x2(jj) = temp_x2(jj+1); temp_x2(jj+1) = temp_val
               END IF
           END DO
       END DO

       valid_count = count_raw
       DO WHILE (valid_count > 1 .AND. temp_x0(valid_count) < -11.0_8)
           valid_count = valid_count - 1
       END DO

       count = valid_count + 6
       block%slice_counts(j) = count

       block%spline_x(j, 4:valid_count+3)  = temp_x(1:valid_count)
       block%spline_x0(j, 4:valid_count+3) = temp_x0(1:valid_count)
       block%spline_x1(j, 4:valid_count+3) = temp_x1(1:valid_count)
       block%spline_x2(j, 4:valid_count+3) = temp_x2(1:valid_count)

       block%spline_x(j, 1) = 0.0_8
       block%spline_x(j, 2) = 0.33_8 * block%spline_x(j, 4)
       block%spline_x(j, 3) = 0.67_8 * block%spline_x(j, 4)
       block%spline_x0(j, 1:3) = block%spline_x0(j, 4)
       block%spline_x1(j, 1:3) = block%spline_x1(j, 4)
       block%spline_x2(j, 1:3) = block%spline_x2(j, 4)

       block%spline_x(j, valid_count+4) = temp_x(valid_count) + 1.0_8
       block%spline_x(j, valid_count+5) = temp_x(valid_count) + 2.0_8
       block%spline_x(j, valid_count+6) = temp_x(valid_count) + 5.0_8

       block%spline_x0(j, valid_count+4) = temp_x0(valid_count) - 0.5_8 * 1.0_8
       block%spline_x0(j, valid_count+5) = temp_x0(valid_count) - 0.5_8 * 2.0_8
       block%spline_x0(j, valid_count+6) = temp_x0(valid_count) - 0.5_8 * 5.0_8

       block%spline_x1(j, valid_count+4:valid_count+6) = temp_x1(valid_count)
       block%spline_x2(j, valid_count+4:valid_count+6) = temp_x2(valid_count)

       ALLOCATE(x_slice(count), x0_slice(count), x1_slice(count), x2_slice(count))
       ALLOCATE(y2_0(count), y2_1(count), y2_2(count))

       x_slice = block%spline_x(j, 1:count)
       x0_slice = block%spline_x0(j, 1:count)
       x1_slice = block%spline_x1(j, 1:count)
       x2_slice = block%spline_x2(j, 1:count)

       CALL spline_calc_derivs_clamped_start(x_slice, x0_slice, count, 0.0_8, y2_0)
       CALL spline_calc_derivs(x_slice, x1_slice, count, y2_1)
       CALL spline_calc_derivs(x_slice, x2_slice, count, y2_2)

       block%spline_y2_0(j, 1:count) = y2_0
       block%spline_y2_1(j, 1:count) = y2_1
       block%spline_y2_2(j, 1:count) = y2_2

       DEALLOCATE(temp_x, temp_x0, temp_x1, temp_x2)
       DEALLOCATE(x_slice, x0_slice, x1_slice, x2_slice, y2_0, y2_1, y2_2)
    END DO
  END SUBROUTINE build_spline_cache
  
  SUBROUTINE find_nearest_spline_indices(grid, n_pts, val, n_out, indices)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n_pts, n_out
    REAL(KIND=8), INTENT(IN) :: val
    REAL(KIND=8), DIMENSION(n_pts), INTENT(IN) :: grid
    INTEGER, DIMENSION(n_out), INTENT(OUT) :: indices
    INTEGER :: i, j1, start_idx

    j1 = 0
    DO i = 1, n_pts - 1
        IF (grid(i) <= val .AND. val < grid(i+1)) THEN
            j1 = i
            EXIT
        END IF
    END DO

    IF (j1 == 0) THEN 
        IF (val >= grid(n_pts)) THEN
            j1 = n_pts - 1
        ELSE
            j1 = 1
        END IF
    END IF
    
    start_idx = j1 - (n_out/2 - 1)
    
    start_idx = MAX(1, start_idx)
    start_idx = MIN(start_idx, n_pts - n_out + 1)
    
    DO i = 1, n_out
        indices(i) = start_idx + i - 1
    END DO
  END SUBROUTINE find_nearest_spline_indices
  
  SUBROUTINE spline_calc_derivs_clamped_start(x, y, n, yp1, y2)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    REAL(KIND=8), INTENT(IN) :: yp1
    REAL(KIND=8), DIMENSION(n), INTENT(IN) :: x, y
    REAL(KIND=8), DIMENSION(n), INTENT(OUT) :: y2
    INTEGER :: i
    REAL(KIND=8) :: p, sig, qn, un
    REAL(KIND=8), ALLOCATABLE :: u(:)
  
    ALLOCATE(u(n))

    y2(1) = -0.5_8
    u(1) = (3.0_8 / (x(2) - x(1))) * ((y(2) - y(1)) / (x(2) - x(1)) - yp1)
  
    DO i = 2, n - 1
      sig = (x(i) - x(i-1)) / (x(i+1) - x(i-1))
      p = sig * y2(i-1) + 2.0_8
      y2(i) = (sig - 1.0_8) / p
      u(i) = (y(i+1) - y(i)) / (x(i+1) - x(i)) - (y(i) - y(i-1)) / (x(i) - x(i-1))
      u(i) = (6.0_8 * u(i) / (x(i+1) - x(i-1)) - sig * u(i-1)) / p
    END DO

    qn = 0.0_8
    un = 0.0_8
    y2(n) = (un - qn * u(n-1)) / (qn * y2(n-1) + 1.0_8)

    DO i = n - 1, 1, -1
      y2(i) = y2(i) * y2(i+1) + u(i)
    END DO
  
    DEALLOCATE(u)
   END SUBROUTINE spline_calc_derivs_clamped_start


  ! High-collisionality (Pfirsch-Schlueter) fit 
  SUBROUTINE calculate_tokamak_fit(r_in, cmul_in, Rmajor_in, iota_in, kn_in, &
                                   g11_out, g31_out, g33_out)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: r_in, cmul_in, Rmajor_in, iota_in, kn_in
    REAL(KIND=8), INTENT(OUT) :: g11_out, g31_out, g33_out
    REAL(KIND=8) :: a1, a2, a3, a4, a5, b2, b3, b4, b6, c1, c2
    REAL(KIND=8) :: d0, d1, d2, d3, d4, d5, d6, d7
    REAL(KIND=8) :: aio, amul, akn, rmul, xkn2, epsilon_tokfit
    REAL(KIND=8) :: g1, g2, g3, g4
    REAL(KIND=8) :: db0, gb, db, dpl, dps, dbp, d31_0
    REAL(KIND=8) :: d11n_ps, g33n, d33n_ps, gc4, d11_0, d33_0
    REAL(KIND=8), PARAMETER :: b0 = 1.0_8

    a1=0.9733; a2=0.3899; a3=0.06778; a4=-1.75; a5=-1.75
    b2=1.0340; b3=-0.6689; b4=0.6666667; b6=0.3333333
    c1=-0.9665; c2=1.75
    d0=4.0/3.0; d1=3.4229; d2=-2.5766; d3=-0.6039; d4=2.0/3.0
    d5=-1.1776; d6=0.6756; d7=1.8436

    aio  = ABS(iota_in)
    amul = ABS(cmul_in)
    akn  = ABS(kn_in)
    rmul = Rmajor_in * amul
    xkn2 = akn**2
    epsilon_tokfit = r_in / Rmajor_in

    g1 = SQRT(akn / epsilon_tokfit) / aio
    g2 = epsilon_tokfit * xkn2
    g3 = epsilon_tokfit * xkn2 * aio
    g4 = epsilon_tokfit * aio

    db0 = a1 * g1
    gb  = 1.0 / (epsilon_tokfit**b4 * akn * aio**b6)
    db  = db0 * (1.0 + b3 * (akn * epsilon_tokfit)**2) / (1.0 + b2 * gb * SQRT(rmul))
    dpl = a2 * g2 / rmul
    dps = a3 * g3 / rmul**2
    dbp = (db**a4 + dpl**a4)**(1.0/a4)
    d31_0 = (dbp**a5 + dps**a5)**(1.0/a5)

    d11n_ps = d0 * (akn / aio)**2 * (1.0 + d1*epsilon_tokfit**3.6 * (1.0 + d2*aio**1.6) + &
                                         d3*epsilon_tokfit**2 * (1.0 - xkn2))
    g33n    = 1.0 + d5*(epsilon_tokfit*akn)**d7 + d6*epsilon_tokfit**3 * aio**2.5
    d33n_ps = d4 * g33n

    gc4   = (1.0 + c1*(epsilon_tokfit*akn)**c2) * g4
    d11_0 = (d11n_ps + d31_0 / gc4) * amul
    d33_0 = (d33n_ps - d31_0 * gc4) / amul

    g11_out = -d11_0 / (b0**2)
    g31_out = d31_0 * iota_in / (aio * b0)
    g33_out = -d33_0 / g33n
  END SUBROUTINE calculate_tokamak_fit
  

  FUNCTION linear_interp_on_slice(x_grid, y_slice, n_pts, x_val)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n_pts
    REAL(KIND=8), INTENT(IN) :: x_val
    REAL(KIND=8), DIMENSION(n_pts), INTENT(IN) :: x_grid, y_slice
    REAL(KIND=8) :: linear_interp_on_slice
    INTEGER :: i1, i2
    REAL(KIND=8) :: w
    CALL find_indices_and_weight_1d(x_grid, n_pts, x_val, i1, i2, w)
    linear_interp_on_slice = linear_interp_1d(y_slice(i1), y_slice(i2), w)
  END FUNCTION linear_interp_on_slice
  
  
  SUBROUTINE calculate_high_efield_correction(block, target_efield, target_cmul, &
                                            g11_e0, g13_e0, g33_e0, &
                                            g11_out, g31_out, g33_out)
    IMPLICIT NONE
    TYPE(dkes_radial_block), INTENT(IN) :: block
    REAL(KIND=8), INTENT(IN) :: target_efield, target_cmul
    REAL(KIND=8), INTENT(IN) :: g11_e0, g13_e0, g33_e0
    REAL(KIND=8), INTENT(OUT) :: g11_out, g31_out, g33_out

    ! variables from neotransp
    REAL(KIND=8) :: a, b, a2, b2, xln, atg, fac, epsilon_tokfit
    REAL(KIND=8) :: val11, val31, val110, val310, d11, d31, d33, eres1
    REAL(KIND=8), PARAMETER :: pa1=3.733333, pa2=1.408160, pa3=0.469388, pa4=0.333333
    REAL(KIND=8), PARAMETER :: pa5=2.333333, pa6=1.523808, pa7=1.600000, pa8=0.380952, pa9=0.177777

    epsilon_tokfit = block%r_val / block%R00
    fac = epsilon_tokfit * block%iota
    
    ! a is normalized e-field, b is normalized collisionality
    a = ABS(target_efield / fac)
    b = ABS(target_cmul * block%R00 / block%iota)
    
    a2 = a**2
    b2 = b**2

    IF (a2 + b2 >= 16.0_8 .AND. b >= 1.0_8) THEN
        val11 = pa1 / (a2 + b2) * (1.0_8 - pa3 / (a2 + b2) * (1.0_8 - 4.0_8 * a2 / (a2 + b2)))
        val31 = -1.0_8 / (a2 + b2) * (pa7 - (pa8 - (pa9 + pa6 * a2) / (a2 + b2)) / (a2 + b2))
    ELSE IF (a >= 4.0_8) THEN
        val11 = pa1 * (1.0_8 - (b2 - pa2) / a2) / a2
        val31 = -(pa7 - (pa8 - b2) / a2) / a2
    ELSE
        xln = LOG(((1.0_8 + a)**2 + b2) / ((1.0_8 - a)**2 + b2))
        IF (b > 1.0E-6_8) THEN
            IF (a > 2.0_8) THEN
                atg = ATAN(2.0_8 * b / (a2 + b2 - 1.0_8)) / b
            ELSE
                atg = (ATAN((1.0_8 - a) / b) + ATAN((1.0_8 + a) / b)) / b
            END IF
        ELSE
            IF (a < 1.0_8) THEN
                atg = MY_PI / b
            ELSE
                atg = 0.0_8
            END IF
        END IF
        val11 = 2.0_8*(pa5+3.0_8*a2-b2-a*(1.0_8+a2-b2)*xln) + (1.0_8+a2*(2.0_8+a2)-b2*(2.0_8-b2+6.0_8*a2))*atg
        val31 = 2.0_8*(pa4+3.0_8*a2-b2-a*(a2-b2)*xln) - (1.0_8-a2**2+6.0_8*a2*b2-b2**2)*atg
    END IF

    IF (b > 4.0_8) THEN
        val110 = pa1 * (1.0_8 - pa3 / b2) / b2
        val310 = -(pa7 - (pa8 - pa9 / b2) / b2) / b2
    ELSE IF (b < 1.0E-3_8) THEN
        val110 = 2.0_8 * pa5 + MY_PI / b
        val310 = 2.0_8 * pa4 - MY_PI / b
    ELSE
        atg = ATAN(1.0_8 / b) / b
        val110 = 2.0_8 * (pa5 - b2 + (1.0_8 - b2 * (2.0_8 - b2)) * atg)
        val310 = 2.0_8 * (pa4 - b2 - (1.0_8 - b2**2) * atg)
    END IF

    ! Corrections to E=0 coefficients 
    d11 = g11_e0 * val11 / val110
    d31 = g13_e0 * val31 / val310
    d33 = g33_e0 ! D33 independent of E-field 

    ! PS correction for D11
    eres1 = pa4 * fac / b
    d11 = d11 / (1.0_8 + (target_efield / eres1)**2)

    g11_out = d11
    g31_out = d31
    g33_out = d33

  END SUBROUTINE calculate_high_efield_correction
  
  SUBROUTINE calculate_low_coll_fit(block, target_cmul, target_efield, g11_out)
    IMPLICIT NONE
    TYPE(dkes_radial_block), INTENT(IN) :: block
    REAL(KIND=8), INTENT(IN) :: target_cmul, target_efield
    REAL(KIND=8), INTENT(OUT) :: g11_out
    
    REAL(KIND=8) :: arg_1overnu, arg_sqrtnu
    REAL(KIND=8) :: efield_eff
    
    efield_eff = MIN(ABS(target_efield), block%EFIELD_upper)

    IF (block%EFIELD_upper <= 0.0_8 .OR. efield_eff < 1.0E-10_8) THEN
        g11_out = block%g11_0 / target_cmul
    ELSE
        arg_1overnu = (target_cmul / block%g11_0)**block%ex_Er
        arg_sqrtnu  = ( (efield_eff**1.5_8) / (block%g11_Er * SQRT(target_cmul)) )**block%ex_Er
        
        g11_out = (arg_1overnu + arg_sqrtnu)**(-1.0_8 / block%ex_Er)
    END IF
    
    g11_out = -ABS(g11_out)
  END SUBROUTINE calculate_low_coll_fit
  

  SUBROUTINE find_indices_and_weight_1d(grid_array,n_pts,val,i1,i2,w)
    IMPLICIT NONE
    REAL(KIND=8),DIMENSION(:),INTENT(IN) :: grid_array
    INTEGER,INTENT(IN) :: n_pts
    REAL(KIND=8),INTENT(IN) :: val
    INTEGER,INTENT(OUT) :: i1,i2
    REAL(KIND=8),INTENT(OUT) :: w
    INTEGER::j

    i1 = 0; i2 = 0; w = 0.0
    IF(n_pts < 1) RETURN
    IF(val <= grid_array(1)) THEN
      i1=1; i2=1; w=0.0; RETURN
    END IF
    IF(val >= grid_array(n_pts)) THEN
      i1=n_pts; i2=n_pts; w=1.0; RETURN
    END IF
    DO j=1,n_pts-1
      IF(val >= grid_array(j) .AND. val <= grid_array(j+1)) THEN
        i1=j; i2=j+1; EXIT
      END IF
    END DO
    IF(i1 > 0) THEN
      IF(ABS(grid_array(i2) - grid_array(i1)) < 1.0E-20) THEN
        w = 0.0
      ELSE
        w = (val - grid_array(i1)) / (grid_array(i2) - grid_array(i1))
      END IF
    END IF
  END SUBROUTINE find_indices_and_weight_1d

  FUNCTION linear_interp_1d(y1, y2, w)
    IMPLICIT NONE
    REAL(KIND=8),INTENT(IN) :: y1,y2,w
    REAL(KIND=8) :: linear_interp_1d

    IF(y1 < -900.0 .OR. y2 < -900.0) THEN
      linear_interp_1d = -999.0
    ELSE
      linear_interp_1d = y1*(1.0-w) + y2*w
    END IF
  END FUNCTION linear_interp_1d
  
  SUBROUTINE spline_calc_derivs(x, y, n, y2)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    REAL(KIND=8), DIMENSION(n), INTENT(IN) :: x, y
    REAL(KIND=8), DIMENSION(n), INTENT(OUT) :: y2
    INTEGER :: i
    REAL(KIND=8) :: p, sig
    REAL(KIND=8), ALLOCATABLE :: u(:)
    ALLOCATE(u(n))
    y2(1) = 0.0_8; u(1) = 0.0_8
    DO i = 2, n - 1
      sig = (x(i) - x(i-1)) / (x(i+1) - x(i-1))
      p = sig * y2(i-1) + 2.0_8
      y2(i) = (sig - 1.0_8) / p
      u(i) = (y(i+1) - y(i)) / (x(i+1) - x(i)) - (y(i) - y(i-1)) / (x(i) - x(i-1))
      u(i) = (6.0_8 * u(i) / (x(i+1) - x(i-1)) - sig * u(i-1)) / p
    END DO
    y2(n) = 0.0_8
    DO i = n - 1, 1, -1
      y2(i) = y2(i) * y2(i+1) + u(i)
    END DO
    DEALLOCATE(u)
  END SUBROUTINE spline_calc_derivs

  FUNCTION spline_interp(xa, x, y, y2, n)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    REAL(KIND=8), INTENT(IN) :: xa
    REAL(KIND=8), DIMENSION(n), INTENT(IN) :: x, y, y2
    REAL(KIND=8) :: spline_interp
    INTEGER :: k, klo, khi
    REAL(KIND=8) :: h, b, a
    klo = 1; khi = n
    DO WHILE (khi - klo > 1)
      k = (khi + klo) / 2
      IF (x(k) > xa) THEN; khi = k; ELSE; klo = k; END IF
    END DO
    h = x(khi) - x(klo)
    IF (h == 0.0_8) THEN; spline_interp = y(klo); RETURN; END IF
    a = (x(khi) - xa) / h
    b = (xa - x(klo)) / h
    spline_interp = a*y(klo) + b*y(khi) + ((a**3-a)*y2(klo) + (b**3-b)*y2(khi))*(h**2)/6.0_8
  END FUNCTION spline_interp
  
  ! minorradiusW7AS value from its text file
  SUBROUTINE read_minorradiusW7AS(status_out)
    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: status_out
    INTEGER           :: file_unit, io_stat_check
    LOGICAL           :: file_exists_flag

    status_out = 0

    IF (LEN_TRIM(MINOR_RADIUS_W7AS_FILE) == 0) THEN
      MINOR_RADIUS_W7AS_FILE = 'vmec_io/minorradiusW7AS.txt'
      WRITE(*,*) 'WARNING: MINOR_RADIUS_W7AS_FILE not in namelist'
    END IF
    
    INQUIRE(file=TRIM(MINOR_RADIUS_W7AS_FILE), exist=file_exists_flag)
    IF (.NOT. file_exists_flag) THEN
      WRITE(*,*) "ERROR: '", TRIM(MINOR_RADIUS_W7AS_FILE), "' not found"
      status_out = 301
      RETURN
    END IF

    OPEN(newunit=file_unit, file=TRIM(MINOR_RADIUS_W7AS_FILE), status='old', action='read', iostat=io_stat_check)

    IF (io_stat_check /= 0) THEN
      status_out = 302
      RETURN
    END IF

    READ(file_unit, *, iostat=io_stat_check) minor_radius_W7AS_value
    IF (io_stat_check /= 0) THEN
      status_out = 303
      CLOSE(file_unit)
      RETURN
    END IF

    CLOSE(file_unit)
    PRINT *, 'Loaded minorradiusW7AS = ', minor_radius_W7AS_value
  END SUBROUTINE read_minorradiusW7AS
  
   ! header data (a, psi_a) from VMEC file
  SUBROUTINE read_vmec_header_data(status_out)
    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: status_out
    INTEGER           :: header_unit, io_stat_check
    LOGICAL           :: file_exists_flag

    status_out = 0

    IF (LEN_TRIM(VMEC_HEADER_FILE) == 0) THEN
      VMEC_HEADER_FILE = 'vmec_io/vmec_header_data.txt'
      WRITE(*,*) 'WARNING: VMEC_HEADER_FILE not in namelist'
    END IF
    
    INQUIRE(file=TRIM(VMEC_HEADER_FILE), exist=file_exists_flag)
    IF (.NOT. file_exists_flag) THEN
      WRITE(*,*) "ERROR: '", TRIM(VMEC_HEADER_FILE), "' not found"
      status_out = 301
      RETURN
    END IF

    OPEN(newunit=header_unit, file=TRIM(VMEC_HEADER_FILE), status='old', action='read', iostat=io_stat_check)

    IF (io_stat_check /= 0) THEN
      status_out = 302
      RETURN
    END IF

    READ(header_unit, *, iostat=io_stat_check) module_ABC
    IF (io_stat_check /= 0) THEN
      status_out = 303
      CLOSE(header_unit)
      RETURN
    END IF

    READ(header_unit, *, iostat=io_stat_check) module_psi_a
    IF (io_stat_check /= 0) THEN
      status_out = 304
      CLOSE(header_unit)
      RETURN
    END IF

    CLOSE(header_unit)

    PRINT *, 'Loaded ABC = ', module_ABC
    PRINT *, 'Loaded psi_a   = ', module_psi_a

  END SUBROUTINE read_vmec_header_data


  SUBROUTINE sort_arrays_by_cmul(n, c_arr, g_arr, idx_arr)
    INTEGER, INTENT(IN) :: n
    REAL(KIND=8), INTENT(INOUT) :: c_arr(n), g_arr(n)
    INTEGER, INTENT(INOUT) :: idx_arr(n)
    INTEGER :: i, j, tmp_idx
    REAL(KIND=8) :: tmp_c, tmp_g
    DO i = 1, n - 1
        DO j = 1, n - i
            IF (c_arr(j) > c_arr(j+1)) THEN
                tmp_c = c_arr(j); c_arr(j) = c_arr(j+1); c_arr(j+1) = tmp_c
                tmp_g = g_arr(j); g_arr(j) = g_arr(j+1); g_arr(j+1) = tmp_g
                tmp_idx = idx_arr(j); idx_arr(j) = idx_arr(j+1); idx_arr(j+1) = tmp_idx
            END IF
        END DO
    END DO
  END SUBROUTINE sort_arrays_by_cmul

  SUBROUTINE sort_array(arr, n)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    REAL(KIND=8), DIMENSION(n), INTENT(INOUT) :: arr
    INTEGER :: i, j
    REAL(KIND=8) :: temp
    IF (n < 2) RETURN
    DO i = 1, n - 1
        DO j = 1, n - i
            IF (arr(j) > arr(j+1)) THEN
                temp = arr(j)
                arr(j) = arr(j+1)
                arr(j+1) = temp
            END IF
        END DO
    END DO
  END SUBROUTINE sort_array

  SUBROUTINE post_process_and_sort_dkes_data()
    IMPLICIT NONE
    INTEGER :: i, j, count_unique, count_unique_efield
    REAL(KIND=8), ALLOCATABLE :: temp_unique_cmuls(:), temp_unique_efields(:)
    REAL(KIND=8), ALLOCATABLE :: sorted_unique_efields(:)
    
    DO i = 1, nr
        IF (dkes_data(i)%num_points > 0) THEN
            ALLOCATE(temp_unique_cmuls(dkes_data(i)%num_points))
            temp_unique_cmuls = 0.0_8
            
            count_unique = 0
            DO j = 1, dkes_data(i)%num_points
                IF (count_unique == 0 .OR. .NOT. ANY(ABS(dkes_data(i)%cmul_list(j) - temp_unique_cmuls(1:count_unique)) < comparison_epsilon)) THEN
                    count_unique = count_unique + 1
                    temp_unique_cmuls(count_unique) = dkes_data(i)%cmul_list(j)
                END IF
            END DO
            
            dkes_data(i)%num_unique_cmuls = count_unique
            ALLOCATE(dkes_data(i)%unique_cmuls(count_unique))
            dkes_data(i)%unique_cmuls = temp_unique_cmuls(1:count_unique)
            DEALLOCATE(temp_unique_cmuls)
            
            CALL sort_array(dkes_data(i)%unique_cmuls, dkes_data(i)%num_unique_cmuls)
            
            CALL build_spline_cache(dkes_data(i))
            
        ELSE
            dkes_data(i)%num_unique_cmuls = 0
        END IF

        WRITE(*,*) 'Diagnostic for r ', i, ' (r = ', dkes_data(i)%r_val, ' m)'

        IF (dkes_data(i)%num_points > 0) THEN
            WRITE(*,*) 'Total data read: ', dkes_data(i)%num_points
            WRITE(*,*) 'Unique CMUL values found: ', dkes_data(i)%num_unique_cmuls
            WRITE(*,*) 'CMUL range: ', &
                dkes_data(i)%unique_cmuls(1), ' to ', dkes_data(i)%unique_cmuls(dkes_data(i)%num_unique_cmuls)
!            PRINT *, 'Unique CMUL values:'
!            WRITE(*,*) (dkes_data(i)%unique_cmuls(j), j=1, dkes_data(i)%num_unique_cmuls)

            ! EFIELD analysis
            ALLOCATE(temp_unique_efields(dkes_data(i)%num_points))
            temp_unique_efields = 0.0_8
            count_unique_efield = 0
            DO j = 1, dkes_data(i)%num_points
                IF (count_unique_efield == 0 .OR. .NOT. ANY(ABS(dkes_data(i)%efield_list(j) - temp_unique_efields(1:count_unique_efield)) < comparison_epsilon)) THEN
                    count_unique_efield = count_unique_efield + 1
                    temp_unique_efields(count_unique_efield) = dkes_data(i)%efield_list(j)
                END IF
            END DO
            
            ALLOCATE(sorted_unique_efields(count_unique_efield))
            sorted_unique_efields = temp_unique_efields(1:count_unique_efield)
            DEALLOCATE(temp_unique_efields)
            CALL sort_array(sorted_unique_efields, count_unique_efield)

            WRITE(*,*) 'Unique EFIELD values found: ', count_unique_efield
            WRITE(*,*) 'EFIELD Range: ', &
                sorted_unique_efields(1), ' to ', sorted_unique_efields(count_unique_efield)
!            PRINT *, 'Unique EFIELD values:'
!            WRITE(*,*) (sorted_unique_efields(j), j=1, count_unique_efield)
            DEALLOCATE(sorted_unique_efields)

        ELSE
            PRINT *, 'No data points found for this r'
        END IF
    END DO
  END SUBROUTINE post_process_and_sort_dkes_data
    
  
  SUBROUTINE read_dkes_data(dkes_filepath_arg, status_out)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: dkes_filepath_arg
    INTEGER, INTENT(OUT) :: status_out
    
    INTEGER :: dkes_unit, io_stat_check, i
    CHARACTER(LEN=1024) :: line_buffer
    INTEGER, ALLOCATABLE :: points_per_radius(:)
    INTEGER :: current_r_idx, current_point_idx
    LOGICAL :: expecting_error_line
    LOGICAL, SAVE :: uses_new_gradient_def = .FALSE.
    
    REAL(KIND=8) :: val1, val2, val3, val4, val5, val6, val7

    PRINT *, 'read_dkes_data subroutine'
    status_out = 0
    
    OPEN(newunit=dkes_unit, file=TRIM(dkes_filepath_arg), status='old', action='read', iostat=io_stat_check)
    IF(io_stat_check /= 0) THEN; status_out = 210; RETURN; END IF
    
    ALLOCATE(points_per_radius(100))
    points_per_radius = 0
    nr = 0
    expecting_error_line = .FALSE.
    
    first_pass_loop: DO
        READ(dkes_unit, '(A)', iostat=io_stat_check) line_buffer
        IF (io_stat_check /= 0) EXIT first_pass_loop
        line_buffer = ADJUSTL(TRIM(line_buffer))
        
        IF (TRIM(line_buffer) == '' .OR. line_buffer(1:1)=='c' .OR. line_buffer(1:1)=='C' .OR. line_buffer(1:1)=='!') CYCLE

        IF (expecting_error_line) THEN
            expecting_error_line = .FALSE.
            CYCLE
        END IF

        READ(line_buffer, *, iostat=io_stat_check) val1, val2, val3, val4, val5, val6, val7
        IF (io_stat_check == 0) THEN
            nr = nr + 1
            CYCLE
        END IF

        READ(line_buffer, *, iostat=io_stat_check) val1, val2, val3, val4, val5
        IF (io_stat_check == 0 .AND. nr > 0) THEN
            points_per_radius(nr) = points_per_radius(nr) + 1
            expecting_error_line = .TRUE. 
        END IF
    END DO first_pass_loop
    
    PRINT *, 'Dimensions found: nr=', nr
    DO i=1, nr
        PRINT *, 'Radius #', i, ' has ', points_per_radius(i), ' data points.'
    END DO

    ! Allocate memory
    IF (ALLOCATED(dkes_data)) DEALLOCATE(dkes_data)
    ALLOCATE(dkes_data(nr))
    DO i = 1, nr
        dkes_data(i)%num_points = points_per_radius(i)
        IF (dkes_data(i)%num_points > 0) THEN
            ALLOCATE( dkes_data(i)%cmul_list(dkes_data(i)%num_points), dkes_data(i)%efield_list(dkes_data(i)%num_points), &
                      dkes_data(i)%g11_list(dkes_data(i)%num_points),  dkes_data(i)%g31_list(dkes_data(i)%num_points), &
                      dkes_data(i)%g33_list(dkes_data(i)%num_points) )
        END IF
    END DO
    
    !Read and fill data
    REWIND(dkes_unit)
    
    current_r_idx = 0
    current_point_idx = 0
    expecting_error_line = .FALSE.
    
    second_pass_loop: DO
        READ(dkes_unit, '(A)', iostat=io_stat_check) line_buffer
        IF (io_stat_check /= 0) EXIT second_pass_loop
        line_buffer = ADJUSTL(TRIM(line_buffer))
        
        IF (TRIM(line_buffer) == '') CYCLE
        
        IF (line_buffer(1:1) == 'c' .OR. line_buffer(1:1) == 'C' .OR. line_buffer(1:1) == '!') THEN
	     IF (INDEX(line_buffer, 'dPsi/dr = 2 psi_a') > 0 .OR. INDEX(line_buffer, 'dPsi/dr =2 psi_a') > 0) THEN
        	     uses_new_gradient_def = .TRUE.
   	     END IF
             IF (line_buffer(1:4) == 'cfit' .AND. current_r_idx > 0) THEN
                 READ(line_buffer(5:), *, iostat=io_stat_check) val1, val2, val3, val4, val5
                 IF (io_stat_check == 0) THEN
                     dkes_data(current_r_idx)%g11_0=val2; dkes_data(current_r_idx)%EFIELD_upper=val3
                     dkes_data(current_r_idx)%g11_Er=val4; dkes_data(current_r_idx)%ex_Er=val5
                 END IF
             END IF
             CYCLE
        END IF
        
        IF (expecting_error_line) THEN
            expecting_error_line = .FALSE.
            CYCLE
        END IF
        
        READ(line_buffer, *, iostat=io_stat_check) val1, val2, val3, val4, val5, val6, val7
        IF (io_stat_check == 0) THEN
            current_r_idx = current_r_idx + 1
            dkes_data(current_r_idx)%r_val=val1; dkes_data(current_r_idx)%R00=val2; dkes_data(current_r_idx)%B00_norm=val3
            dkes_data(current_r_idx)%iota=val4; dkes_data(current_r_idx)%kn=val5; dkes_data(current_r_idx)%ftrap=val6
            dkes_data(current_r_idx)%b2_norm=val7
            current_point_idx = 0
        ELSE
            READ(line_buffer, *, iostat=io_stat_check) val1, val2, val3, val4, val5
            IF(io_stat_check == 0 .AND. current_r_idx > 0) THEN
                current_point_idx = current_point_idx + 1
                dkes_data(current_r_idx)%cmul_list(current_point_idx) = val1
                dkes_data(current_r_idx)%efield_list(current_point_idx) = val2
                dkes_data(current_r_idx)%g11_list(current_point_idx) = val3
                dkes_data(current_r_idx)%g31_list(current_point_idx) = val4
                dkes_data(current_r_idx)%g33_list(current_point_idx) = val5
                expecting_error_line = .TRUE.
            END IF
        END IF
    END DO second_pass_loop
    CLOSE(dkes_unit)
    DEALLOCATE(points_per_radius)
    PRINT *, 'DKES data loaded'

    CALL post_process_and_sort_dkes_data()

  END SUBROUTINE read_dkes_data
 
END MODULE a2dkes
