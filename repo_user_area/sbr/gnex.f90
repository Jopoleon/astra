subroutine GNEX
! Electron flux GNX due to all neutral sources
! GNX  [10^19 particle/m^2/s]
!    (Pereverzev 23-FEB-98)
!   2-NOV-99 GNX redetermined to give the flux density

use scalars, only: HRO, NA1, NAB, TAU, NNCL, NNWM
use status, only: VR, VRO, NE, NEO, TE, NI, NN, SLAT, GNX, SNEBM

implicit none

integer :: j
double precision :: YSN1, YSN2, YSN3, YH, Y, &
    SNNEU, SVIE, SVII, SVREC, SNNR, SNNI

YSN1 = 0.
YSN2 = 0.
YSN3 = 0.
YH = HRO
do J=1, NA1
    include 'fml/snneu'
    YSN1 = YSN1 + VR(J)*NE(J)
    YSN2 = YSN2 + VRO(J)*NEO(J)
    YSN3 = YSN3 + VR(j)*(NE(j)*SNNEU + SNEBM(j))
! Electron sources due to
!  (i)   ionization of the wall neutrals,
!  (ii)  ionization of the beam neutrals
!  (iii) dNE/dt
    GNX(J) = (YSN3 - (YSN1 - YSN2)/TAU )*YH/SLAT(J)
enddo

if (NA1 < NAB) then
    do J=NA1+1, NAB
        GNX(J) = GNX(NA1)
    enddo
endif

end subroutine gnex
