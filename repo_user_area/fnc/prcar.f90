! PRCARR [MW*m^3]: Carbon Radiation power
!       POST D.E. JENSEN R.V. e.a.,  Atomic Data and Nuclear Tables, 
!        Vol.20 (1977) 397.
!      Pereverzev 29-07-96
! Usage:  PE=...+PRCAR*NE*NIZ1
double precision FUNCTION PRCARR(YR)

use const_inc, only: HRO, NA1
use status_inc, only: TE

implicit none

double precision, intent(in) :: YR
integer :: JK
double precision YT, YZ, Y

JK = int(YR/HRO + .5)
if (JK >= NA1) JK = NA1
if (JK <= 0) JK = 1

YT = TE(JK)*1000.
YZ = LOG10(TE(JK))
if(YT <= 13.)                  Y = -14.4*(YZ + 2.2)**2
if(YT > 13.  .and. YT <= 25. ) Y = -1.3 - 3.67*(YZ + 1.9)
if(YT > 25.  .and. YT <= 63. ) Y = -2.8 + 10.*(YZ + 1.4)**2
if(YT > 63.  .and. YT <= 159.) Y = -2.12 - 7.2*(YZ + 1.)**2
if(YT > 159. .and. YT <= 500.) Y = -2.4 - 1.2*(YZ + 0.8)
if(YT > 500. .and. YT <= 4.e3) Y = -3.3 + 0.957*(YZ - 0.26)**2
if(YT > 4.E3)                  Y = -3.2 + 0.357*(YZ - 0.6)

PRCARR = 10.**(Y + 1.)

end function PRCARR
