module bnd_modul

    implicit none

    integer, parameter, private :: DP=kind(1.0D0)
    integer :: nbtab 
    real(DP), allocatable :: rbtab(:),zbtab(:)
    
end module bnd_modul
