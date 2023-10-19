! LI3 []: Internal inductance li(r) [4]
!   Pereverzev 03-APR-08
!   /
!       c*c |
!     Li = ------------*|(Bpol**2)dV   li = Li/(2*pi*R0)
!    4*pi*Ipl*Ipl |
!   /
!   V
!
!   SI: Li[H] = 2*Wi/Ipl**2 ; li[dim.less] = Li[mkHn]/(0.2*pi*R0[m])
!          2*pi*R0[m]*li[dim.less] = 10*Li[mkHn]
!----------------------------------------------------------------------|
! The same as LINT except for calculation of the denominator 
!     as an integral of the current density
!----------------------------------------------------------------------|

double precision function LIG3R(YRO)

use const_inc, only: ROC, HRO, BTOR, NA, GP
use status_inc, only: IPOL, G22, CU, VR, MU
use flight_sim_geometrics ! for flight simulator diagnostics


implicit none

double precision yro


lig3r=geom1d(299)  !li3 from spider or feqis, iter definition

return
end function LIG3R
