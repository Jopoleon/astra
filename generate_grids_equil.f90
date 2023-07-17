program generate_grids_equil

implicit none

 character*80 directori !

	integer iaug,iiter,idemo
	
	write(*,*) 'do aug? (0 for no or 1 for yes)'
	read(*,*) iaug
	write(*,*) 'do demo? (0 for no or 1 for yes)'
	read(*,*) idemo
	write(*,*) 'do iter? (0 for no or 1 for yes)'
	read(*,*) iiter

		directori='exp/cnf/'

	if (iaug.eq.1) then
		call generate_files_feqis(directori,'aug')
!		call CEPPON(0,0,0,0.01,0., & 
!			0.,0.,0.,0,'exp/equ/aug_/',13)
	endif

	if (iiter.eq.1) then
		call generate_files_feqis(directori,'demo')
!		call CEPPON(0,0,0,0.01,0., & 
!			0.,0.,0.,0,'exp/equ/iter/',13)
	endif

	if (idemo.eq.1) then
		call generate_files_feqis(directori,'iter')
!		call CEPPON(0,0,0,0.01,0., & 
!			0.,0.,0.,0,'exp/equ/dem_/',13)
	endif


end
