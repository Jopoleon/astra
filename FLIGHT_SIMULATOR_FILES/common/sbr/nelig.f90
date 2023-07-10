!C  Interface ASTRA-sawmod.f

!C Emiliano 2005

      SUBROUTINE NELIG(nchords,yfircs)
	use parameter_inc
	use const_inc
	use status_inc
use parameters_a2equil, only: equil_now
	use numerical_tools, only: qinterp

      implicit none

      real*8 yfircs(NRD)
      integer	k,jmax,OPTION,J,i,OPTION2
		integer nchords	
	integer nline
	parameter(nline=30)
	double precision rc1(5),rc2(5)
	double precision zc1(5),zc2(5)
	character*19 filefirdesc
	parameter(filefirdesc='dat/filefirdesc.dat')
	integer Nx,Nt
	double precision 	Rb(700),rhopol(NA1)
	double precision 	PSIn_grid(700)
	double precision 	PSIs_grid(700)
	double precision 	PSIb
	double precision 	dchord,ddchord(5),delta_R,delta_Z
	double precision 	rc_l(nline,5),zc_l(nline,5),th_chord,nc_l(nline)
	integer i_method
	integer i_chord,i1(556),i2(556),i3
	double precision psi(25,25), R(25), z(25),dycl(556,556)
	integer neqlp,ntetap
	double precision ne_zequil(556),dyy(556)
	double precision ne_equil(556,556)
	double precision psiii(556,556)                  
	double precision yrout(256,256),yzout(256,256)
	integer ieqtype ,nrp   ,nwcd ,jcall                
!	common /com_flux2/psiii
	data jcall/0/
	save rc1,zc1,rc2,zc2,rc_l,zc_l,ddchord,jcall
		
	yfircs(1:nchords)=0.
	

	if (jcall.lt.2) then
	jcall=jcall+1
	open(32,file=filefirdesc)
	do i=1,nchords
		read(32,*) rc1(i),zc1(i),rc2(i),zc2(i)			
	enddo	
	close(32)
	
!	write(*,*) 'fir',rc1,zc1,rc2,zc2
	
	
	do i_chord=1,nchords   ! chords cycle
	
!Construct chord points
	delta_R = rc2(i_chord)-rc1(i_chord)
	delta_Z = zc2(i_chord)-zc1(i_chord)
	dchord = sqrt(delta_R**2.+delta_Z**2.)
	ddchord(i_chord)=dchord/nline

			if (delta_Z.ge.0 .and.&
     &   delta_R.ge.0) then 
				th_chord=atan(delta_Z/delta_R)
			endif
			if (delta_Z.ge.0 .and.&
     &   delta_R.lt.0) then 
				th_chord=GP+atan(delta_Z/delta_R)
			endif
			if (delta_Z.lt.0 .and.&
     &   delta_R.lt.0) then 
				th_chord=GP+atan(delta_Z/delta_R)
			endif
			if (delta_Z.lt.0 .and. &
     &  delta_R.ge.0) then 
				th_chord=GP2+atan(delta_Z/delta_R)
			endif

		do j=1,nline
			rc_l(j,i_chord)=rc1(i_chord)+ddchord(i_chord)*j*cos(th_chord)
			zc_l(j,i_chord)=zc1(i_chord)+ddchord(i_chord)*j*sin(th_chord)
		enddo

	enddo !chord number
	return
	endif


	neqlp=abs(nint(NEQUIL))
	ntetap=abs(nint(MEQUIL))+1

!	write(*,*) '1'	
	yrout(1:neqlp,1:ntetap-1)=equil_now%coord_sys%position%r(1:neqlp,1:ntetap-1) !yrout(neqlp,j)
	yzout(1:neqlp,1:ntetap-1)=equil_now%coord_sys%position%z(1:neqlp,1:ntetap-1) !yrout(neqlp,j)
	psiii(1:neqlp,1:ntetap-1)=equil_now%coord_sys%position%psirz(1:neqlp,1:ntetap-1) !yrout(neqlp,j)
	yrout(1:neqlp,ntetap)=yrout(1:neqlp,1)
	yzout(1:neqlp,ntetap)=yzout(1:neqlp,1)
	psiii(1:neqlp,ntetap)=psiii(1:neqlp,1)
	yrout(1:neqlp,ntetap+1)=yrout(1:neqlp,2)
	yzout(1:neqlp,ntetap+1)=yzout(1:neqlp,2)
	psiii(1:neqlp,ntetap+1)=psiii(1:neqlp,2)


!	write(998,*) yrout(1:neqlp,1:ntetap),yzout(1:neqlp,1:ntetap),psiii(1:neqlp,1:ntetap)
	

!	if (time.gt.3.03) call a_stop

!		write(*,*) 'here',nline,nchords
!		stop
!
!	do j=1,NA1
		rhopol(1:NA1)=((FP(1:NA1)-FP(1)))**0.5/(FP(na1)-FP(1))**0.5
!	enddo
	do i=1,Ntetap

	PSIs_grid(1:neqlp)=sqrt((psiii(1:neqlp,i)-psiii(1,i))/&
     & (psiii(neqlp,i)-psiii(1,i)))

	call	qinterp(rhopol(1:NA1),&
     &  NE(1:NA1),NA1,&
     &  PSIs_grid(1:neqlp),ne_equil(1:neqlp,i),neqlp)
	enddo
!	write(*,*) ne_equil(1:neqlp,1:ntetap)

!		write(*,*) 'here3'
!        call omp_set_dynamic(.false.)  ! this is really essential!

!Follow the chord line
	
	do i_chord=1,nchords   ! chords cycle
		do j=1,nline
			dycl(1:neqlp,1:ntetap)=sqrt(&
     & (yrout(1:neqlp,1:ntetap)-rc_l(j,i_chord))**2.d0+&
     & (yzout(1:neqlp,1:ntetap)-zc_l(j,i_chord))**2.d0)
	i1(1:ntetap)=minloc(dycl(1:neqlp,1:ntetap),1)
	dyy(1:ntetap)=minval(dycl(1:neqlp,1:ntetap),1)
	i2(1)=minloc(dyy(1:ntetap),1)
	dyy(1)=dycl(i1(i2(1)),i2(1))
	if (dyy(1).gt.0.1) then
		nc_l(j)=0.
	else
		nc_l(j)=ne_equil(i1(i2(1)),i2(1))	
	endif
!	write(*,*) j,i_chord,rc_l(j,i_chord),zc_l(j,i_chord),
!     & dyy(1),i1(i2(1)),i2(1),
!     & yrout(i1(i2(1)),i2(1)),yzout(i1(i2(1)),i2(1)),
!     & psis_grid(i1(i2(1))),ne(na1)
!	write(*,*) i_chord,j,nc_l(j)


		enddo
! sum up
	yfircs(i_chord)=ddchord(i_chord)*sum(nc_l(1:nline))*1.d19
	enddo !chord number


!	write(*,*) 'interferometer : ',ne(1),ne(20),ne(41),yfircs(1:nchords)
!	write(*,*) 'interferometer : ',fp(1),fp(20),fp(41),psiii(1,1),psiii(neqlp,1)
!	stop






101	format(6F12.4)

      
			end
