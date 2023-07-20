program generate_grids_equil

implicit none

character(len=80), parameter :: dir='exp/cnf/'

integer :: iaug=1, iiter=0, idemo=0
	
!write(*,*) 'do aug? (0 = no, 1 = yes)'
!read(*,*) iaug
!write(*,*) 'do demo? (0 for no or 1 for yes)'
!read(*,*) idemo
!write(*,*) 'do iter? (0 for no or 1 for yes)'
!read(*,*) iiter

if (iaug == 1) then
    call generate_files_feqis(dir, 'aug')
endif

if (iiter == 1) then
    call generate_files_feqis(dir, 'demo')
endif

if (idemo == 1) then
    call generate_files_feqis(dir, 'iter')
endif

end program generate_grids_equil
