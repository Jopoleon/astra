CEfable This subroutine write coil currents using a reference file

      subroutine coil2spider(ccoils,ncoils,parameters_spider)
      
         use parameters 
			implicit none

       type(type_parameters) parameters_spider
      integer i,j,ncoils
      double precision       ccoils(*)
			integer nlinetot,nlines,ntimes
			parameter (nlinetot=100)
      character*560 reffile
      character*560 destfile
      character rlineget(nlinetot)
      integer lline(nlinetot)
			double precision  rline1(nlinetot),
     &   rline2(nlinetot),rline3(nlinetot)
			double precision  rline4(nlinetot),
     &   rline5(nlinetot),rline6(nlinetot)
			double precision  ccoilold(nlinetot)
			integer  rline8(nlinetot),ncoilold(nlinetot)
			integer npieces(nlinetot)

			destfile=parameters_spider%prename
     >  (1:parameters_spider%kname)//'coil.dat'
	if (ncoils.gt.0) then
			reffile=parameters_spider%prename
     >  (1:parameters_spider%kname)//'req2d_input/coil_ref.dat'
C			write(*,*)  reffile
			
			open(32,file=reffile,status='old')
			read(32,*) nlines,ntimes,rlineget(1)
C			write(*,*)  nlines,ntimes,rlineget(1)
			do i=1,nlines
		  read(32,*) npieces(i)
			read(32,*) rline1(i),rline2(i),rline3(i),rline4(i),
     &     rline5(i),rline6(i),ccoilold(i),rline8(i),ncoilold(i),
     &     rlineget(i)
   		enddo
			
      close(32)			
			
	endif			
			
C    	write(*,*) 'coil2spider2',ccoils(1:ncoils),ncoils
	if (ncoils.eq.0) then
		ncoils=0
	endif
C			write(*,*) npieces(1:nlines)
C			write(*,*) ccoilold(1:nlines)
C			write(*,*) ncoilold(1:nlines)
			
			
	if (ncoils.gt.0) then
			! remapping coil currents
			do i=1,ncoils
			do j=1,nlines
			   if (ncoilold(j).eq.i) then
				   ccoilold(j)=ccoils(i)*1.E-3   !from kA to MA
				 endif
			enddo
			enddo
			
			!write file
			open(32,file=destfile,status='unknown')
			write(32,*) nlines,ntimes,rlineget(1)
			do i=1,nlines
		  write(32,*) npieces(i)
			write(32,100) rline1(i),
     &     rline2(i),rline3(i),rline4(i),
     &     rline5(i),rline6(i),ccoilold(i),
     &     rline8(i),ncoilold(i),
     &     '  ',rlineget(i)
   		enddo
			
      close(32)			
			
!			pause
			
	endif			

 100	format(7F15.10,2I8,2A)			

      end
CEfable


      
