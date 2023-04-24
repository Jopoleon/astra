module vol2d_eq_face
INTERFACE
SUBROUTINE vol2d_eq2(equil_out,fcell2d,fvav)
use imas_ids
 implicit none
 type(type_equilibrium) equil_out
 real(8), dimension(:,:) :: fcell2d
 real(8), dimension(:) :: fvav
END SUBROUTINE vol2d_eq2
END INTERFACE
end module
