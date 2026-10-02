! Part of install_astra8.sh (installer/a7modules/spider): glue between ASTRA 8 and
! the SPIDER sources of a local ASTRA 7 copy. Contains no SPIDER code.
!
! spider_run as called by ASTRA 8 (src/astra/gs_solver.F90, A_EQUIL).
! Calls the A7 SPIDER entry (renamed to spider_run_a7 at build time) and fills the
! fields that exist only in the A8 type_equilibrium, so that later A8 code
! (json_rw.f90, metrics.f90) never touches unassociated pointers.
subroutine spider_run(ncoils, vcoils, equil_in, equil_out, params)
    use imas_ids, only: type_equilibrium
    use spider_params, only: type_parameters
    implicit none
    integer, intent(in) :: ncoils
    double precision :: vcoils(ncoils)
    type(type_equilibrium) :: equil_in, equil_out
    type(type_parameters) :: params
    integer :: nrho, nthe, nr, nz, i
    double precision :: psiext
    interface
        subroutine spider_run_a7(ncoils, vcoils, equil_in, equil_out, params)
            use imas_ids, only: type_equilibrium
            use parameters, only: type_parameters
            integer :: ncoils
            double precision :: vcoils(ncoils)
            type(type_equilibrium) :: equil_in, equil_out
            type(type_parameters) :: params
        end subroutine spider_run_a7
    end interface

    params%a7_parameters%nteta = params%ntheta
! no_circuit_eq (A8 static free-boundary switch) has no counterpart in A7 SPIDER: ignored
    call spider_run_a7(ncoils, vcoils, equil_in, equil_out, params%a7_parameters)

    nrho = 0
    if (associated(equil_out%profiles_1d%psi)) nrho = size(equil_out%profiles_1d%psi)
    nthe = 0
    if (associated(equil_out%coord_sys%position%r)) nthe = size(equil_out%coord_sys%position%r, 2)

    if (.not. associated(equil_out%profiles_1d%rho_tor_norm)) then
        allocate(equil_out%profiles_1d%rho_tor_norm(nrho))
        if (nrho > 1) then
            if (associated(equil_out%profiles_1d%rho_tor)) then
                if (equil_out%profiles_1d%rho_tor(nrho) > 0.d0) then
                    equil_out%profiles_1d%rho_tor_norm = equil_out%profiles_1d%rho_tor/equil_out%profiles_1d%rho_tor(nrho)
                else
                    equil_out%profiles_1d%rho_tor_norm = (/ ((i-1.d0)/(nrho-1.d0), i=1, nrho) /)
                endif
            else
                equil_out%profiles_1d%rho_tor_norm = (/ ((i-1.d0)/(nrho-1.d0), i=1, nrho) /)
            endif
        endif
    endif
    if (.not. associated(equil_out%profiles_1d%squareness)) then
        allocate(equil_out%profiles_1d%squareness(nrho))
        equil_out%profiles_1d%squareness = 0.d0
    endif
    if (.not. associated(equil_out%coord_sys%position%theta2d)) then
        allocate(equil_out%coord_sys%position%theta2d(nthe))
        if (nthe > 1) equil_out%coord_sys%position%theta2d = (/ ((i-1.d0)*8.d0*atan(1.d0)/(nthe-1.d0), i=1, nthe) /)
    endif
    if (.not. associated(equil_out%coord_sys%darea)) then
        allocate(equil_out%coord_sys%darea(nrho, nthe))
        equil_out%coord_sys%darea = 0.d0
    endif
    if (.not. associated(equil_out%coord_sys%jphi)) then
        allocate(equil_out%coord_sys%jphi(nrho, nthe))
        equil_out%coord_sys%jphi = 0.d0
    endif
    nr = 0
    nz = 0
    if (associated(equil_out%eqgeometry%rectgrid%psirz2d)) then
        nr = size(equil_out%eqgeometry%rectgrid%psirz2d, 1)
        nz = size(equil_out%eqgeometry%rectgrid%psirz2d, 2)
    endif
    if (.not. associated(equil_out%eqgeometry%rectgrid%fdia2d)) then
        allocate(equil_out%eqgeometry%rectgrid%fdia2d(nr, nz))
        equil_out%eqgeometry%rectgrid%fdia2d = 0.d0
    endif
    if (nrho > 0) then
        equil_out%global_param%psiaxis = equil_out%profiles_1d%psi(1)
    endif
! A7 (EQUIL/COUPLING_SCHEME/a_spider.f) took PSIEXT from psib_ext/f_psib_ext after the run;
! A8 reads it as -2*pi*equil_out%global_param%psiext.
    psiext = 0.d0
    if (params%k_fixfree == 1) then
        if (params%k_grid == 0) call psib_ext(psiext)
        if (params%k_grid == 1) call f_psib_ext(psiext)
    endif
    equil_out%global_param%psiext = psiext
end subroutine spider_run
