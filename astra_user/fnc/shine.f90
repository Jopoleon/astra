! SDDN [10#19 n/s]:  Integral {0, R} ( SVD1*NDEUT**2 ) dV
! Neutron rate from thermal D-D reaction
!   (Polevoy 25-JUN-97,  29-JUN-2016)
double precision FUNCTION SHINER(YR)

use const_inc, only: CNB1

implicit none

double precision, intent(in) :: YR
integer :: J, JSRNUM
double precision :: Y, Y1
character*12 SHNAME
DATA SHNAME/'dat/shth.dat'/


JSRNUM = ABS(CNB1)
Y = 0.
DO J=1, JSRNUM 
   if(JSRNUM > 99) then
      write(*, *) '>>> Too many NBI sources > 99 (see NBI1TR)'
      return
   endif   
   if(JSRNUM < 10) write(SHNAME(12:12), '(I1)') JSRNUM
   if(JSRNUM > 9 .and. JSRNUM < 100) write(SHNAME(11:12), '(I2)') JSRNUM
   if(JSRNUM > 99) then
      write(*, *) '>>> Too many NBI sources > 99 (see NBI1TR)'
      return
   endif   
   open(38, FILE=SHNAME, ERR=1)
   READ(38, *, ERR=1, END=1) Y1
   Y = Y + Y1
   goto 2
 1 SHINER =Y
 2 close(38)
ENDDO
SHINER =Y

end function shiner
