subroutine f_ext_fil(pcequi, ncequi)

use sp_parameters, only: njlim, amu0, nip, njp
use comrec, only: ni, nj, ue, zaindk

implicit none

integer, intent(in) :: ncequi
real*8, intent(in) :: pcequi(*)

integer :: i, j, iq
real*8 :: egcurr, vrftfa

do j=1,nj
   do i=1,ni
      ue(i,j)=0.d0
   enddo
enddo

do iq=1,ncequi
   egcurr=PCEQui(iq)*amu0	  
   do i=1,ni
      do j=1,nj
         vrftfa=zaindk(i,j,iq)
         ue(i,j)=ue(i,j)+vrftfa*egcurr
      enddo
   enddo
enddo

return
end subroutine f_ext_fil

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine f_rdexf(ncequi)

use sp_parameters, only: njlim, nip, njp
use iopath, only: path
use comrec
use comevl

implicit none

integer, intent(in) :: ncequi

integer :: i, j, iq, nk
real*8 aindk(nip,njp)
character(len=80) :: fname

write(fname,'(a,a)') TRIM(path),'/exf.wr'
open(1,file=fname)

read(1,*) nk,nequi

do iq=1,ncequi
   read(1,*) ((aindk(i,j),j=1,nj),i=1,ni)
   do j=1,nj
      do i=1,ni
         zaindk(i,j,iq)=aindk(i,j)
      enddo
   enddo
enddo

close(1)

return
end subroutine f_rdexf
