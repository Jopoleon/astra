! Replacement for the A7 SPIDER module "schemas": the A8 IMAS-like types
! (src/misc/imas_ids.f90) are a superset of the A7 CPO types.
module schemas
    use imas_ids
    implicit none
    integer, parameter :: RKIND = SELECTED_REAL_KIND(10)
    real(RKIND), parameter :: TWOPI = 2*3.141592653589793238462643383279502884197_rkind
end module schemas
