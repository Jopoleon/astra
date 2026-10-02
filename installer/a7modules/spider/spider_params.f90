! A8 expects "use spider_params, only: type_parameters" (src/astra/gs_solver.F90)
! with components ntheta and no_circuit_eq. The A7 SPIDER type (module "parameters")
! has nteta and no no_circuit_eq; extend it and map in spider_run (spider_a8_adapter.f90).
module spider_params
    use parameters, only: a7_parameters => type_parameters
    implicit none
    type, extends(a7_parameters) :: type_parameters
        integer :: ntheta = 92
        integer :: no_circuit_eq = 1
    end type type_parameters
end module spider_params
