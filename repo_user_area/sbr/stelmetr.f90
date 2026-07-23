subroutine STELMETR

use scalars
use status
use read_input, only: AWD, exp_file, equ_file
use pi_const

implicit none

integer, parameter :: Nrrect=256, Nzrect=256, eqdsk_unit=11
integer :: i, j, nrho_surf, nthe_surf, n_Rrect, n_Zrect
double precision :: raxis, zaxis, drhot_eq, dpsin_rect, psib, Rmin, Rmax, zmin, zmax, dr, dz


write(*,*) 'fps',FP(1),FPO(1),FP(NA1),FPO(NA1),FP(NA),FPO(NA),sum(FP(1:NA1)-FPO(1:NA1))/TAU/NA1
write(*,*) 'fps',(FP(NA1)-FP(NA))/HRO/(GP2*BTOR*SRHO(NA)),0.2/SRHO(NA)/BTOR/SG11(NA)*IPL+MV(NA)
write(*,*) 'fps',(SG11(NA)*(FP(NA1)-FP(NA))/HRO/(GP2*BTOR*SRHO(NA))+SG12(NA))/SG22(NA)*SRHO(NA)
write(*,*) 'fps',(SG11(1)*(FP(2)-FP(1))/HRO/(GP2*BTOR*SRHO(1))+SG12(1))/SG22(1)*SRHO(1)
write(*,*) 'fps',sum(1./(2/RHO(1:NA1)*BTOR/0.4/CC(1:NA1)*(SG22(1:na1))**2.))*HRO
     
if ((IPEQL) == 6) return

SG12=0. ! half grid
SG21=0. ! half grid
SG11=G22*IPOL/RTOR  ! half grid
SG22=RTOR*IPOL ! half grid

SG12(1:na1) = -(0.25*xrho(1:na1)**5 + 0.1)*SG11(1:na1)  !only example SG12 for MV=(0.25*xrho(1:na1)**5 + 0.1)
SG21(1:na1) = -(0.07*xrho(1:na1)**5 + 0.02)*xrho(1:na1)  !only example SG21

DTEQL  =  0.05
open(32,file='vmec_io/input_metric.dat')
write(32,*) 1
write(32,*) 0.,shif(1:na1),elon(1:na1), &
 tria(1:na1),g33(1:na1),ipol(1:na1),vr(1:na1), &
 slat(1:na1),g11(1:na1),g22(1:na1),droda(1:na1), &
 shiv(1:na1),squarn(1:na1),ametr(1:na1),volum(1:na1),vrs(1:na1), &
 rho(1:na1),sg11(1:na1),sg12(1:na1), &
 sg21(1:na1),sg22(1:na1)
close(32)
return
end subroutine STELMETR
