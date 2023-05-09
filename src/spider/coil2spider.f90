! This subroutine write coil currents using a reference file

subroutine coil2spider(ccoils,ncoils,parameters_spider)

use parameters_a2spider, only: type_parameters

implicit none

integer, parameter :: nlinetot=100
integer, intent(in) :: ncoils
double precision, intent(in) :: ccoils(*)
type(type_parameters), intent(in) :: parameters_spider

integer i, j, nlines, ntimes
integer, dimension(nlinetot) :: rline8, ncoilold, npieces
double precision, dimension(nlinetot) :: rline1, rline2, rline3, &
    rline4, rline5, rline6, ccoilold
character(len=160) :: reffile, destfile
character rlineget(nlinetot)

destfile = TRIM(parameters_spider%prename) // '/coil.dat'
if (ncoils.gt.0) then
   reffile = TRIM(parameters_spider%prename) // '/req2d_input/coil_ref.dat'

   open(32, file=reffile, status='old')
   read(32, *) nlines, ntimes, rlineget(1)
   do i=1, nlines
      read(32,*) npieces(i)
      read(32,*) rline1(i), rline2(i), rline3(i)  , rline4(i), &
                 rline5(i), rline6(i), ccoilold(i), rline8(i), &
                 ncoilold(i), rlineget(i)
   enddo
   close(32)
endif      

if (ncoils.gt.0) then
! remapping coil currents
   do i=1,ncoils
      do j=1,nlines
         if (ncoilold(j).eq.i) then
           ccoilold(j) = ccoils(i)*1.E-3   !from kA to MA
         endif
      enddo
   enddo

!write file
   open(32, file=destfile, status='unknown')
   write(32, *) nlines, ntimes, rlineget(1)
   do i=1, nlines
      write(32,*) npieces(i)
      write(32,100) rline1(i), rline2(i), rline3(i), &
                    rline4(i), rline5(i), rline6(i), ccoilold(i), &
                    rline8(i), ncoilold(i), '  ', rlineget(i)
   enddo

   close(32)      
   
endif      

 100  format(7F15.10, 2I8, 2A) 

return
end subroutine coil2spider
