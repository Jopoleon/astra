module json_vars

use json_module, only : json_file, json_core, json_value, json_ck

implicit none

integer, parameter :: nunit=25
integer :: jid
integer :: n_var, n_varx, n_const, n_intern, n_intern2, n_prof, &
    n_profx, n_equil_sig, n_equil_prof, n_equil_rect, n_equil_coord
character(len=6), allocatable, dimension(:) :: varNames, &
    varxNames, constNames, internNames, intern2Names, &
    profNames, profxNames, equil_sigNames, equil_profNames, &
    equil_rectNames, equil_coordNames
type(json_file) :: astra_vars
type(json_value), pointer :: varPtr, varxPtr, constPtr, &
     internPtr, intern2Ptr, profPtr, profxPtr, &
     equil_sigPtr, equil_profPtr, equil_rectPtr, equil_coordPtr

contains

!---------------------------------------------------------------------
    subroutine get_subdict(label, nvars, names_out, jsonOut)

    character(len=*), intent(in) :: label
    integer, intent(out) :: nvars  
    character(len=6), intent(out), allocatable, dimension(:) :: names_out
    type(json_value), pointer, intent(out) :: jsonOut

    logical :: found
    integer :: j
    character(KIND=JSON_CK, len=:), allocatable :: sname
    type(json_value), pointer :: dictPointer
    type(json_core) :: jCore

    call astra_vars%info(label, n_children=nvars)
    allocate(names_out(nvars))
    call astra_vars%get(label, jsonOut, found)
    do j=1, nvars
        call jCore%get_child(jsonOut, j, dictPointer, found)
        call jCore%info(dictPointer, name=sname)
        names_out(j) = sname
    enddo

    return
    end subroutine get_subdict
  
!---------------------------------------------------------------------
    subroutine read_metadata

    logical :: status_ok
    character(len=:), allocatable :: error_msg
    character(len=240) :: file_in
    file_in = 'astra_variables.json'
    call astra_vars%initialize()
    call astra_vars%load(filename=trim(file_in))
    if (astra_vars%failed()) then
        write(*, *) 'Error: '
        call astra_vars%check_for_errors(status_ok, error_msg)    
        write(*, *) 'Error: ', error_msg
    endif

! Get sub-dictionaries objects & lists
    call get_subdict('variables'  , n_var    ,     varNames,     varPtr)
    call get_subdict('variables_x', n_varx   ,    varxNames,    varxPtr)
    call get_subdict('constants'  , n_const  ,   constNames,   constPtr)
    call get_subdict('internal'   , n_intern ,  internNames,  internPtr)
    call get_subdict('intern2'    , n_intern2, intern2Names, intern2Ptr)
    call get_subdict('profiles'   , n_prof   ,    profNames,    profPtr)
    call get_subdict('profiles_x' , n_profx  ,   profxNames,   profxPtr)
    call get_subdict('equil_signals' , n_equil_sig  , equil_sigNames  , equil_sigPtr)
    call get_subdict('equil_profiles', n_equil_prof , equil_profNames , equil_profPtr)
    call get_subdict('equil_rect'    , n_equil_rect , equil_rectNames , equil_rectPtr)
    call get_subdict('equil_coord'   , n_equil_coord, equil_coordNames, equil_coordPtr)

    return
    end subroutine read_metadata

!---------------------------------------------------------------------
    subroutine write_json

    use parameters_a2equil, only: equil_now
    use parameter_inc
    use const_inc
    use status_inc
    use io_mod, only: AWD, exp_file, equ_file
    use debugger, only: debug

    integer :: j, jrho, ios, j_call=1, nrho_surf, nthe_surf, nR, nZ
    character(len=120) :: json_out
    double precision, dimension(200) :: time_traces, equil_traces
    character(KIND=JSON_CK, len=:), allocatable :: sunit, sdesc, sname
    type(json_core) :: jCore

    save j_call

    nrho_surf = SIZE(equil_now%profiles_1d%rho_tor_norm)
    nthe_surf = SIZE(equil_now%coord_sys%position%teta2d)
    nR = SIZE(equil_now%eqgeometry%rectgrid%r2d)
    nZ = SIZE(equil_now%eqgeometry%rectgrid%z2d)

    if (debug > 0) write(*, '(A, 3i)') 'Starting a2json', j_call, nrho_surf, nthe_surf

    write(json_out, '(5A, i0, A)') TRIM(awd), '/ncdf_out/', TRIM(exp_file), TRIM(equ_file), '-', j_call, '.json'

    open(nunit, file=TRIM(json_out), iostat=ios)
    write(nunit, '(A/)') '{'

!------
! ASTRA
!------
    
    write(nunit, '(A/)') '"astra": {'

! Scalars
    call write_scalar_block(varPtr, DEVAR)
    call write_scalar_block(varxPtr, DEVARX)
    call write_scalar_block(constPtr, constValues)
    call write_scalar_block(internPtr, DELOUT)

! Sparse scalars

    time_traces(1)  = TSTART
    time_traces(2)  = TAU
    time_traces(3)  = TAUPRP
    time_traces(4)  = HRO
    time_traces(5)  = HROX
    time_traces(6)  = ALBPL
    time_traces(7)  = NNCX
    time_traces(8)  = QETB
    time_traces(9)  = QFF0B
    time_traces(10) = QFF1B
    time_traces(11) = QFF2B
    time_traces(12) = QFF3B
    time_traces(13) = QFF4B
    time_traces(14) = QFF5B
    time_traces(15) = QFF6B
    time_traces(16) = QFF7B
    time_traces(17) = QFF8B
    time_traces(18) = QFF9B
    time_traces(19) = QITB
    time_traces(20) = QNNB
    time_traces(21) = FTO
    time_traces(22) = FTN
    time_traces(23) = BTN
    time_traces(24) = IPLN
    time_traces(25) = IFBEY
    time_traces(26) = ROB
    time_traces(27) = ROWALL
    time_traces(28) = ROC
    time_traces(29) = ROCO
    time_traces(30) = RON
    time_traces(31) = ROE
    time_traces(32) = ROI
    time_traces(33) = RO0
    time_traces(34) = RO1
    time_traces(35) = RO2
    time_traces(36) = RO3
    time_traces(37) = RO4
    time_traces(38) = RO5
    time_traces(39) = RO6
    time_traces(40) = RO7
    time_traces(41) = RO8
    time_traces(42) = RO9
    time_traces(43) = ROU
    time_traces(44) = VOLUME
    time_traces(45) = PSIAX
    time_traces(46) = PSIBO
    time_traces(47) = GP
    time_traces(48) = GP2
    time_traces(49) = PSIFB
    time_traces(50) = PSIFBO
    time_traces(51) = PSIEXT
    time_traces(52) = PSPLEX
    time_traces(53) = PSIEXO
    time_traces(54) = PSPLXO
    time_traces(55) = IPLFBE
    time_traces(56) = ATREQ
    time_traces(57) = PTREQ
    time_traces(58) = RBDOT
    time_traces(59) = BBDOT
    time_traces(60) = NA
    time_traces(61) = NA1
    time_traces(62) = NA1N
    time_traces(63) = NA1E
    time_traces(64) = NA1I
    time_traces(65) = NA1U
    time_traces(66) = NA10
    time_traces(67) = NA11
    time_traces(68) = NA12
    time_traces(69) = NA13
    time_traces(70) = NA14
    time_traces(71) = NA15
    time_traces(72) = NA16
    time_traces(73) = NA17
    time_traces(74) = NA18
    time_traces(75) = NA19
    time_traces(76) = NAB
    time_traces(77) = ITREQ
    time_traces(78) = NITOT
    time_traces(79) = NSTEPS
    time_traces(80) = QBEAM
    call write_scalar_block(intern2Ptr, time_traces)

!-----------
! Profiles
!-----------

! Exp profiles
    jid = 0    
! Get sub-dictionaries dimensions
    do j=1, n_profx
        call write_array((/NA1/), EXT(1:NA1, j), profxPtr)
    enddo

! CAR profiles
    jid = 0    
    do j=1, 128 ! CAR*
        call write_array((/NA1/), CAR(1:NA1, j), profPtr)
    enddo

! Sparse profiles
    call write_array((/NA1/),  AIMPT(1:NA1), profPtr)
    call write_array((/NA1/),  AMAIN(1:NA1), profPtr)
    call write_array((/NA1/),  AMETR(1:NA1), profPtr)
    call write_array((/NA1/),  AREAT(1:NA1), profPtr)
    call write_array((/NA1/),  B0DB2(1:NA1), profPtr)
    call write_array((/NA1/),   BDB0(1:NA1), profPtr)
    call write_array((/NA1/),  BDB02(1:NA1), profPtr)
    call write_array((/NA1/),  BMAXT(1:NA1), profPtr)
    call write_array((/NA1/),  BMINT(1:NA1), profPtr)
    call write_array((/NA1/),     CC(1:NA1), profPtr)
    call write_array((/NA1/),     CD(1:NA1), profPtr)
    call write_array((/NA1/),     CE(1:NA1), profPtr)
    call write_array((/NA1/),     CI(1:NA1), profPtr)
    call write_array((/NA1/),     CN(1:NA1), profPtr)
    call write_array((/NA1/),  CNPAD(1:NA1), profPtr)
    call write_array((/NA1/),  CNPAP(1:NA1), profPtr)
    call write_array((/NA1/),  CNPAR(1:NA1), profPtr)
    call write_array((/NA1/),     CU(1:NA1), profPtr)
    call write_array((/NA1/),   CUBM(1:NA1), profPtr)
    call write_array((/NA1/),   CUBS(1:NA1), profPtr)
    call write_array((/NA1/),  CUECR(1:NA1), profPtr)
    call write_array((/NA1/),   CUFI(1:NA1), profPtr)
    call write_array((/NA1/),   CUFW(1:NA1), profPtr)
    call write_array((/NA1/),  CUICR(1:NA1), profPtr)
    call write_array((/NA1/),   CULH(1:NA1), profPtr)
    call write_array((/NA1/),  CUTOR(1:NA1), profPtr)
    call write_array((/NA1/),     CV(1:NA1), profPtr)
    call write_array((/NA1/),     DC(1:NA1), profPtr)
    call write_array((/NA1/),     DN(1:NA1), profPtr)
    call write_array((/NA1/),  DDNEO(1:NA1), profPtr)
    call write_array((/NA1/), DDNEOD(1:NA1), profPtr)
    call write_array((/NA1/),  DIMP1(1:NA1), profPtr)
    call write_array((/NA1/),  DIMP2(1:NA1), profPtr)
    call write_array((/NA1/),  DIMP3(1:NA1), profPtr)
    call write_array((/NA1/),  DLNEO(1:NA1), profPtr)
    call write_array((/NA1/), DLNEOD(1:NA1), profPtr)
    call write_array((/NA1/),  DRODA(1:NA1), profPtr)
    call write_array((/NA1/),   ELON(1:NA1), profPtr)
    call write_array((/NA1/),   EQFF(1:NA1), profPtr)
    call write_array((/NA1/),   EQPF(1:NA1), profPtr)
    call write_array((/NA1/),     ER(1:NA1), profPtr)
    call write_array((/NA1/),     F0(1:NA1), profPtr)
    call write_array((/NA1/),    F0O(1:NA1), profPtr)
    call write_array((/NA1/),     F1(1:NA1), profPtr)
    call write_array((/NA1/),    F1O(1:NA1), profPtr)
    call write_array((/NA1/),     F2(1:NA1), profPtr)
    call write_array((/NA1/),    F2O(1:NA1), profPtr)
    call write_array((/NA1/),     F3(1:NA1), profPtr)
    call write_array((/NA1/),    F3O(1:NA1), profPtr)
    call write_array((/NA1/),     F4(1:NA1), profPtr)
    call write_array((/NA1/),    F4O(1:NA1), profPtr)
    call write_array((/NA1/),     F5(1:NA1), profPtr)
    call write_array((/NA1/),    F5O(1:NA1), profPtr)
    call write_array((/NA1/),     F6(1:NA1), profPtr)
    call write_array((/NA1/),    F6O(1:NA1), profPtr)
    call write_array((/NA1/),     F7(1:NA1), profPtr)
    call write_array((/NA1/),    F7O(1:NA1), profPtr)
    call write_array((/NA1/),     F8(1:NA1), profPtr)
    call write_array((/NA1/),    F8O(1:NA1), profPtr)
    call write_array((/NA1/),     F9(1:NA1), profPtr)
    call write_array((/NA1/),    F9O(1:NA1), profPtr)
    call write_array((/NA1/),   FOFB(1:NA1), profPtr)
    call write_array((/NA1/),     FP(1:NA1), profPtr)
    call write_array((/NA1/),    FPO(1:NA1), profPtr)
    call write_array((/NA1/),FP_NORM(1:NA1), profPtr)
    call write_array((/NA1/),     FV(1:NA1), profPtr)
    call write_array((/NA1/),    G11(1:NA1), profPtr)
    call write_array((/NA1/),    G22(1:NA1), profPtr)
    call write_array((/NA1/),   G22E(1:NA1), profPtr)
    call write_array((/NA1/),    G33(1:NA1), profPtr)
    call write_array((/NA1/),   G33E(1:NA1), profPtr)
    call write_array((/NA1/),    G41(1:NA1), profPtr)
    call write_array((/NA1/),    G42(1:NA1), profPtr)
    call write_array((/NA1/),    G43(1:NA1), profPtr)
    call write_array((/NA1/),    G44(1:NA1), profPtr)
    call write_array((/NA1/),    G45(1:NA1), profPtr)
    call write_array((/NA1/),     GN(1:NA1), profPtr)
    call write_array((/NA1/), GRADRO(1:NA1), profPtr)
    call write_array((/NA1/),     HC(1:NA1), profPtr)
    call write_array((/NA1/),     HE(1:NA1), profPtr)
    call write_array((/NA1/),   IPOL(1:NA1), profPtr)
    call write_array((/NA1/),   MRHO(1:NA1), profPtr)
    call write_array((/NA1/),     MU(1:NA1), profPtr)
    call write_array((/NA1/),     MV(1:NA1), profPtr)
    call write_array((/NA1/),   NALF(1:NA1), profPtr)
    call write_array((/NA1/),  NDEUT(1:NA1), profPtr)
    call write_array((/NA1/),     NE(1:NA1), profPtr)
    call write_array((/NA1/),    NEO(1:NA1), profPtr)
    call write_array((/NA1/),   NHE3(1:NA1), profPtr)
    call write_array((/NA1/),  NHYDR(1:NA1), profPtr)
    call write_array((/NA1/),     NI(1:NA1), profPtr)
    call write_array((/NA1/),   NIBM(1:NA1), profPtr)
    call write_array((/NA1/),  NIMPT(1:NA1), profPtr)
    call write_array((/NA1/),    NIO(1:NA1), profPtr)
    call write_array((/NA1/),   NIZ1(1:NA1), profPtr)
    call write_array((/NA1/),   NIZ2(1:NA1), profPtr)
    call write_array((/NA1/),   NIZ3(1:NA1), profPtr)
    call write_array((/NA1/),  NMAIN(1:NA1), profPtr)
    call write_array((/NA1/),     NN(1:NA1), profPtr)
    call write_array((/NA1/),  NNBM1(1:NA1), profPtr)
    call write_array((/NA1/),  NNBM2(1:NA1), profPtr)
    call write_array((/NA1/),  NNBM3(1:NA1), profPtr)
    call write_array((/NA1/),  NRATE(1:NA1), profPtr)
    call write_array((/NA1/),  NTRIT(1:NA1), profPtr)
    call write_array((/NA1/),  PBEAM(1:NA1), profPtr)
    call write_array((/NA1/),  PBLON(1:NA1), profPtr)
    call write_array((/NA1/),  PBOL1(1:NA1), profPtr)
    call write_array((/NA1/),  PBOL2(1:NA1), profPtr)
    call write_array((/NA1/),  PBOL3(1:NA1), profPtr)
    call write_array((/NA1/),  PBPER(1:NA1), profPtr)
    call write_array((/NA1/),    PDE(1:NA1), profPtr)
    call write_array((/NA1/),    PDI(1:NA1), profPtr)
    call write_array((/NA1/),     PE(1:NA1), profPtr)
    call write_array((/NA1/),   PEBM(1:NA1), profPtr)
    call write_array((/NA1/),  PEECR(1:NA1), profPtr)
    call write_array((/NA1/),   PEFW(1:NA1), profPtr)
    call write_array((/NA1/),  PEICR(1:NA1), profPtr)
    call write_array((/NA1/),  PEIQI(1:NA1), profPtr)
    call write_array((/NA1/),   PELH(1:NA1), profPtr)
    call write_array((/NA1/),  PELON(1:NA1), profPtr)
    call write_array((/NA1/),  PEPER(1:NA1), profPtr)
    call write_array((/NA1/),  PERIM(1:NA1), profPtr)
    call write_array((/NA1/),  PETOT(1:NA1), profPtr)
    call write_array((/NA1/),  PFAST(1:NA1), profPtr)
    call write_array((/NA1/),     PI(1:NA1), profPtr)
    call write_array((/NA1/),   PIBM(1:NA1), profPtr)
    call write_array((/NA1/),   PIFW(1:NA1), profPtr)
    call write_array((/NA1/),  PIICR(1:NA1), profPtr)
    call write_array((/NA1/),  PITOT(1:NA1), profPtr)
    call write_array((/NA1/),   PRAD(1:NA1), profPtr)
    call write_array((/NA1/),   PRES(1:NA1), profPtr)
    call write_array((/NA1/),  PSXR1(1:NA1), profPtr)
    call write_array((/NA1/),  PSXR2(1:NA1), profPtr)
    call write_array((/NA1/),  PSXR3(1:NA1), profPtr)
    call write_array((/NA1/),     QE(1:NA1), profPtr)
    call write_array((/NA1/),    QF0(1:NA1), profPtr)
    call write_array((/NA1/),    QF1(1:NA1), profPtr)
    call write_array((/NA1/),    QF2(1:NA1), profPtr)
    call write_array((/NA1/),    QF3(1:NA1), profPtr)
    call write_array((/NA1/),    QF4(1:NA1), profPtr)
    call write_array((/NA1/),    QF5(1:NA1), profPtr)
    call write_array((/NA1/),    QF6(1:NA1), profPtr)
    call write_array((/NA1/),    QF7(1:NA1), profPtr)
    call write_array((/NA1/),    QF8(1:NA1), profPtr)
    call write_array((/NA1/),    QF9(1:NA1), profPtr)
    call write_array((/NA1/),     QI(1:NA1), profPtr)
    call write_array((/NA1/),     QN(1:NA1), profPtr)
    call write_array((/NA1/),     QU(1:NA1), profPtr)
    call write_array((/NA1/),    RHO(1:NA1), profPtr)
    call write_array((/NA1/),RHO_POL(1:NA1), profPtr)
    call write_array((/NA1/),  RUPAR(1:NA1), profPtr)
    call write_array((/NA1/),  RUPFR(1:NA1), profPtr)
    call write_array((/NA1/),  RUPYR(1:NA1), profPtr)
    call write_array((/NA1/),  SCUBM(1:NA1), profPtr)
    call write_array((/NA1/),    SD0(1:NA1), profPtr)
    call write_array((/NA1/),    SD1(1:NA1), profPtr)
    call write_array((/NA1/),    SD2(1:NA1), profPtr)
    call write_array((/NA1/),    SD3(1:NA1), profPtr)
    call write_array((/NA1/),    SD4(1:NA1), profPtr)
    call write_array((/NA1/),    SD5(1:NA1), profPtr)
    call write_array((/NA1/),    SD6(1:NA1), profPtr)
    call write_array((/NA1/),    SD7(1:NA1), profPtr)
    call write_array((/NA1/),    SD8(1:NA1), profPtr)
    call write_array((/NA1/),    SD9(1:NA1), profPtr)
    call write_array((/NA1/),    SDN(1:NA1), profPtr)
    call write_array((/NA1/), SF0TOT(1:NA1), profPtr)
    call write_array((/NA1/), SF1TOT(1:NA1), profPtr)
    call write_array((/NA1/), SF2TOT(1:NA1), profPtr)
    call write_array((/NA1/), SF3TOT(1:NA1), profPtr)
    call write_array((/NA1/), SF4TOT(1:NA1), profPtr)
    call write_array((/NA1/), SF5TOT(1:NA1), profPtr)
    call write_array((/NA1/), SF6TOT(1:NA1), profPtr)
    call write_array((/NA1/), SF7TOT(1:NA1), profPtr)
    call write_array((/NA1/), SF8TOT(1:NA1), profPtr)
    call write_array((/NA1/), SF9TOT(1:NA1), profPtr)
    call write_array((/NA1/),  SGNEO(1:NA1), profPtr)
    call write_array((/NA1/), SGNEOD(1:NA1), profPtr)
    call write_array((/NA1/),  SHEAR(1:NA1), profPtr)
    call write_array((/NA1/),   SHIF(1:NA1), profPtr)
    call write_array((/NA1/),   SHIV(1:NA1), profPtr)
    call write_array((/NA1/),   SLAT(1:NA1), profPtr)
    call write_array((/NA1/),     SN(1:NA1), profPtr)
    call write_array((/NA1/),  SNEBM(1:NA1), profPtr)
    call write_array((/NA1/), SNIBM1(1:NA1), profPtr)
    call write_array((/NA1/), SNIBM2(1:NA1), profPtr)
    call write_array((/NA1/), SNIBM3(1:NA1), profPtr)
    call write_array((/NA1/),  SNNBM(1:NA1), profPtr)
    call write_array((/NA1/),  SNTOT(1:NA1), profPtr)
    call write_array((/NA1/),  SQEPS(1:NA1), profPtr)
    call write_array((/NA1/), SQUARN(1:NA1), profPtr)
    call write_array((/NA1/),   SRHO(1:NA1), profPtr)
    call write_array((/NA1/),   SXHO(1:NA1), profPtr)
    call write_array((/NA1/),     TE(1:NA1), profPtr)
    call write_array((/NA1/),    TEO(1:NA1), profPtr)
    call write_array((/NA1/),     TI(1:NA1), profPtr)
    call write_array((/NA1/),    TIO(1:NA1), profPtr)
    call write_array((/NA1/),     TN(1:NA1), profPtr)
    call write_array((/NA1/),   TRIA(1:NA1), profPtr)
    call write_array((/NA1/),   TTRQ(1:NA1), profPtr)
    call write_array((/NA1/),  TTRQI(1:NA1), profPtr)
    call write_array((/NA1/),   ULON(1:NA1), profPtr)
    call write_array((/NA1/),   UPAR(1:NA1), profPtr)
    call write_array((/NA1/),  UPARO(1:NA1), profPtr)
    call write_array((/NA1/),    UPL(1:NA1), profPtr)
    call write_array((/NA1/),   UPS0(1:NA1), profPtr)
    call write_array((/NA1/),  UPS0O(1:NA1), profPtr)
    call write_array((/NA1/),   UPS1(1:NA1), profPtr)
    call write_array((/NA1/),  UPS1O(1:NA1), profPtr)
    call write_array((/NA1/),   UPS2(1:NA1), profPtr)
    call write_array((/NA1/),  UPS2O(1:NA1), profPtr)
    call write_array((/NA1/),  VIMP1(1:NA1), profPtr)
    call write_array((/NA1/),  VIMP2(1:NA1), profPtr)
    call write_array((/NA1/),  VIMP3(1:NA1), profPtr)
    call write_array((/NA1/),  VOLUM(1:NA1), profPtr)
    call write_array((/NA1/),     VP(1:NA1), profPtr)
    call write_array((/NA1/),   VPFP(1:NA1), profPtr)
    call write_array((/NA1/),   VPOL(1:NA1), profPtr)
    call write_array((/NA1/),     VR(1:NA1), profPtr)
    call write_array((/NA1/),    VRO(1:NA1), profPtr)
    call write_array((/NA1/),    VRS(1:NA1), profPtr)
    call write_array((/NA1/),   VTOR(1:NA1), profPtr)
    call write_array((/NA1/),     XC(1:NA1), profPtr)
    call write_array((/NA1/),     XI(1:NA1), profPtr)
    call write_array((/NA1/),   XRHO(1:NA1), profPtr)
    call write_array((/NA1/),  XUPAD(1:NA1), profPtr)
    call write_array((/NA1/),  XUPAP(1:NA1), profPtr)
    call write_array((/NA1/),  XUPAR(1:NA1), profPtr)
    call write_array((/NA1/),    ZEF(1:NA1), profPtr)
    call write_array((/NA1/),   ZEF1(1:NA1), profPtr)
    call write_array((/NA1/),   ZEF2(1:NA1), profPtr)
    call write_array((/NA1/),   ZEF3(1:NA1), profPtr)
    call write_array((/NA1/),   ZIM1(1:NA1), profPtr)
    call write_array((/NA1/),   ZIM2(1:NA1), profPtr)
    call write_array((/NA1/),   ZIM3(1:NA1), profPtr)
    call write_array((/NA1/),  ZIMPT(1:NA1), profPtr)
    call write_arr(  (/NA1/),  ZMAIN(1:NA1), profPtr) ! No comma at the end

    write(nunit, '(A/)') '},' ! End of "astra" dictionary

!------------
! Equilibrium
!------------

    write(nunit, '(A/)') '"equil": {'

! Scalars

    equil_traces(1) = equil_now%global_param%toroid_field%b0
    equil_traces(2) = equil_now%global_param%betpol
    equil_traces(3) = equil_now%global_param%i_plasma
    equil_traces(4) = equil_now%global_param%li3
    equil_traces(5) = equil_now%global_param%psibound
    equil_traces(6) = equil_now%global_param%psiaxis
    equil_traces(7) = equil_now%global_param%toroid_field%r0
    equil_traces(8) = equil_now%global_param%Vloop
    call write_scalar_block(equil_sigPtr, equil_traces(1: 8))

! 1 d profiles
    jid = 0
    if (debug > 0) write(*, *) 'Writing equilibrium entries to json file'
    call write_array((/nrho_surf/), equil_now%profiles_1d%areat  , equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%bdb0   , equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%bmaxt  , equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%bmint  , equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%dpsidv , equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%elongation, equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%f_dia  , equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%ffprime, equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%fofb   , equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%g1     , equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%g2     , equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%ggradro, equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%gm1    , equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%gm4    , equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%gm41   , equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%gm5    , equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%perim  , equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%phi    , equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%pprime , equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%pressure, equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%psi    , equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%q      , equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%r_inboard , equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%r_outboard, equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%rho_tor_norm, equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%shif   , equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%surface, equil_profPtr)
    call write_array((/nrho_surf/), equil_now%profiles_1d%volume , equil_profPtr)

    if (debug > 0) write(*, *) 'Writing mag.surf quantities'
    jid = 0
    call write_array((/nrho_surf, nthe_surf/), equil_now%coord_sys%position%r, equil_coordPtr)
    call write_array((/nrho_surf, nthe_surf/), equil_now%coord_sys%position%rmin, equil_coordPtr)
    call write_array((/nrho_surf, nthe_surf/), equil_now%coord_sys%position%psirz, equil_coordPtr)
    call write_array((/nthe_surf/), equil_now%coord_sys%position%teta2d, equil_coordPtr)
    call write_array((/nrho_surf, nthe_surf/), equil_now%coord_sys%position%z, equil_coordPtr)

    jid = 0
    if (debug > 0) write(*, *) 'Writing equil cartesian quantities'
    call write_array((/nR, nZ/), equil_now%eqgeometry%rectgrid%psirz2d, equil_rectPtr)
    call write_array((/nR, nZ/), equil_now%eqgeometry%rectgrid%fdia2d , equil_rectPtr)
    call write_array((/nR/), equil_now%eqgeometry%rectgrid%r2d, equil_rectPtr)
    call write_arr((/nZ/), equil_now%eqgeometry%rectgrid%z2d, equil_rectPtr) ! No comma
    write(nunit, '(A)') '}' ! End of "equil" dictionary

!-----------
! Close json
!-----------

    write(nunit, '(A)') '}'
    close(nunit)

    j_call = j_call + 1

    write(*, '(A)') '   Written file ' // TRIM(json_out)

    return
    end subroutine write_json

!---------------------------------------------------------------
    subroutine write_scalar_block(json_in, scalar_list)

    double precision, intent(in), dimension(*) :: scalar_list
    type(json_value), intent(in), pointer :: json_in

    logical :: found
    integer :: nvars, j
    character(KIND=JSON_CK, len=:), allocatable :: sunit, sdesc, sname
    type(json_value), pointer :: dictPointer
    type(json_core) :: jCore

101 format('    "', A, '": {"units": "', A, '", "long_name": "', A, '", "data": ', es16.8e3, '},')

    call jCore%info(json_in, n_children=nvars)
    do j=1, nvars
        call jCore%get_child(json_in, j, dictPointer, found)
        call jCore%info(dictPointer, name=sname)
        call jCore%get(dictPointer, 'units', sunit, found)
        call jCore%get(dictPointer, 'desc' , sdesc, found)
        write(nunit, 101) sname, sunit, sdesc, scalar_list(j)
    enddo
    
    return
    end subroutine write_scalar_block

!---------------------------------------------------------------
    subroutine write_array(dims, arr_in, json_in)

    integer, intent(in), dimension(:) :: dims
    double precision, intent(in), dimension(*) :: arr_in
    type(json_value), intent(in), pointer :: json_in

    call write_arr(dims, arr_in, json_in)
    write(nunit, '(A)') ','

    return
    end subroutine write_array

!---------------------------------------------------------------
    subroutine ndim_string(dims, ndim, sdim)

    integer, intent(in), dimension(:) :: dims
    integer, intent(out) :: ndim  
    character(len=120), intent(out) :: sdim

    if (SIZE(dims) == 1) then
        write(sdim, '(A,i0,A)') '[', dims(1), ']'
        ndim = dims(1)
    else if (SIZE(dims) == 2) then
        write(sdim, '(A,i0,A,i0,A)') '[', dims(1), ', ', dims(2), ']'
        ndim = dims(1)*dims(2)
    endif
     
    return
    end subroutine ndim_string

!---------------------------------------------------------------
    subroutine write_arr(dims, arr_in, json_in)

    integer, intent(in), dimension(:) :: dims
    double precision, intent(in), dimension(*) :: arr_in
    type(json_value), intent(in), pointer :: json_in

    logical :: found
    integer :: i, j, ndim, ij
    double precision :: abs_val
    double precision, dimension(:), allocatable :: array
    double precision, dimension(:, :), allocatable :: array2
    character(KIND=JSON_CK, len=:), allocatable :: sunit, sdesc, sname
    character(len=120) :: sdim
    type(json_value), pointer :: dictPointer
    type(json_core) :: jCore

102 format('    "', A, '": {"dims": ', A, ', "units": "', A, '", "long_name": "', A, '", "data": [')

    call ndim_string(dims, ndim, sdim)
    jid = jid + 1
    call jCore%get_child(json_in, jid, dictPointer, found)
    call jCore%get(dictPointer, 'units', sunit, found)
    call jCore%get(dictPointer, 'desc', sdesc, found)
    call jCore%info(dictPointer, name=sname)
    write(nunit, 102, advance="no") sname, TRIM(sdim), sunit, sdesc


    if (SIZE(dims) == 1) then
        allocate(array(ndim))
        do i=1, ndim
            abs_val = ABS(arr_in(i))
            if (abs_val < 1e-20 .or. abs_val > 1e20) then
                array(i) = 0.
            else
                array(i) = arr_in(i)
            endif
        enddo
        if (MAXVAL(ABS(array)) > 0.) then
            write(nunit, *)
            write(nunit, '(5(es16.8e3, ","))') (array(i), i=1, ndim-1)
            write(nunit, '(es16.8e3)', advance='no') array(ndim) ! No comma after last array entry
        endif
 
    else if (SIZE(dims) == 2) then
        allocate(array2(dims(1), dims(2)))
        do i=1, dims(1)
            do j=1, dims(2)
                ij = dims(1)*(j-1) + i
                abs_val = ABS(arr_in(ij))
                if (abs_val < 1e-20 .or. abs_val > 1e20) then
                    array2(i, j) = 0.
                else
                    array2(i, j) = arr_in(ij)
                endif
            enddo
        enddo
        if (MAXVAL(ABS(array2)) > 0.) then
            do i=1, dims(1)
                write(nunit, '(A)') '['
                write(nunit, '(5(es16.8e3, ","))') (array2(i, j), j=1, dims(2)-1)
                write(nunit, '(es16.8e3)', advance='no') array2(i, dims(2)) ! No comma after last array entry
                if (i == dims(1)) then
                    write(nunit, '(A)') ']'
                else
                    write(nunit, '(A)') '],'
                endif
            enddo
        endif
    endif

    write(nunit, '(A)') ']}' ! No comma after last array entry

    return
    end subroutine write_arr

end module json_vars
