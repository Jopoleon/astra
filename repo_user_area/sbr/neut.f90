subroutine NEUT
!-----------------------------------------------------------22.01.97---|
! Input: ABC,NA1,NA,AMJ,NAB,NNCX
!        AMAIN(j),ZMAIN(j),TE(j),TI(j),NE(j),NI(j),ZEF(j),SNNBM(j)
!        ENCL,ENWM or wall neutral distribution
!        NNCL,NNWM
! Warning: ENCL > 0.5; ENWM =/= 0 if NNWM =/= 0
! Output: NN, TN, ALBPL
!---------------------------------------------CHANGED BY POLEVOY-------|

use scalars, only: NA, NA1, ABC, NAB, ENCL, ENWM, NNCL, NNWM, AMJ, NNCX, ALBPL
use status, only: NRD, NN, TN, NE, TE, NI, TI, SNNBM, AMAIN

implicit none

integer :: J, JJ, JN
double precision, dimension(NRD) :: SCXNI, TEN, SRCNN, YVI, &
   NIN, NN0, THICKN
double precision, dimension(NRD*NRD) :: YKERN, YKERNF
double precision :: AKDLT, Y, SVIE, SVREC, SVCX, YNN0, YXI1, YXI2, &
      YV1, YV2, YHA, YTI, YN0, YZZN, YZ1, YZ2, YZZ, YZZT

!----------------------------------------------------------------------|
if (ENCL .lt. 0.0001) then
   write(*, '(A, F6.3, A)') &
       ">>> NEUT >>> Too low energy of incoming neutrals ENCL =", 1000.*ENCL, "eV"
   write(*, *) "    Setting ENCL = 2 eV"
   ENCL = 0.002
endif
if (ENWM .lt. 0.002) then
   write(*, '(A, F6.3, A)') &
       ">>> NEUT >>> Too low energy of incoming neutrals ENWM =", 1000.*ENWM, "eV"
   write(*, *) "    Setting ENWM = 2 eV"
   ENWM = 0.002
endif

! Neutral generations
YHA = ABC/(NA1 - .5)
do J=1, NA1
   include 'fml/svcx'
   SCXNI(J) = SVCX*NI(J)
   include 'fml/svie'
   YKERN(NRD*J + 4) = SVIE*NE(J)
   include 'fml/svrec'
   SRCNN(J) = SVREC*NE(J)*NI(J) + SNNBM(J)
! V_i is divided by SQRT(3)
   YTI = max(TI(j), 1.d-3)
   YVI(J) = 2.52E5*SQRT(YTI/AMAIN(J))
enddo

YVI(NA1) = YVI(NA)
!...      THICKN(1)=0.5*(SCXNI(1)+YKERN(1,4))*YHA
THICKN(1) = 0.
! YKERN(NRD*J + 4) or YKERN(NRD*4 + J) ?
do J=1, NA
   THICKN(J+1) = THICKN(J) + (SCXNI(J) + YKERN(NRD*J + 4) + SCXNI(J+1) + YKERN(NRD*(J+1) + 4))*YHA*0.5
enddo
do JJ=1, NA1
   SCXNI(JJ) = 0.5*SCXNI(JJ)/YVI(JJ)
   do J=1, NA1
      YZ1 = (THICKN(J) + THICKN(JJ))/YVI(JJ)
      if (J .GE. JJ) then
         YZ2 = (THICKN(J) - THICKN(JJ))/YVI(JJ)
      else
         YZ2 = (THICKN(JJ) - THICKN(J))/YVI(JJ)
      endif
      YKERN (NRD*J + JJ) = AKDLT(YZ1) + AKDLT(YZ2)
      YKERNF(NRD*J + JJ) = YKERN(NRD*J + JJ)*SCXNI(JJ)
   enddo
enddo
!-------------------------------      Zero generation:
YNN0 = NNCL + NNWM + 1.E-11
YXI1 = (NNCL + 1.E-11)/YNN0
YXI2 = NNWM/YNN0
YV1  = 4.37E5*SQRT(ENCL/AMJ)
YV2  = 4.37E5*SQRT(ENWM/AMJ)
do J=1, NA1
   NIN(J) = 0.5*SRCNN(J)/YNN0/YVI(J)
   YVI(J) = NIN(J)*TI(J)
enddo
do J=1, NA1
   YZZ  = 0.5*(NIN(1)*YKERN(NRD*J +1) + NIN(NA1)*YKERN(NRD*J +NA1))
   YZZT = 0.5*(YVI(1)*YKERN(NRD*J +1) + YVI(NA1)*YKERN(NRD*J +NA1))
   do JJ=2, NA
      YZZ  = YZZ  + NIN(JJ)*YKERN(NRD*J +JJ)
      YZZT = YZZT + YVI(JJ)*YKERN(NRD*J +JJ)
   enddo
   NN0(J) = YHA*YZZ
   TEN(J) = YHA*YZZT
   YZ1    = (THICKN(NA1) + THICKN(J))/YV1
   YZ2    = (THICKN(NA1) - THICKN(J))/YV1
   YN0    = YXI1*(AKDLT(YZ1) + AKDLT(YZ2))
   NN0(J) = NN0(J) + YN0
   TEN(J) = TEN(J) + ENCL*YN0
   if (NNWM .GT. 0.) then
      YZ1 = (THICKN(NA1) + THICKN(J))/YV2
      YZ2 = (THICKN(NA1) - THICKN(J))/YV2
      YN0 = YXI2*(AKDLT(YZ1) + AKDLT(YZ2))
      NN0(J) = NN0(J) + YN0
      TEN(J) = TEN(J) + ENWM*YN0
   endif
enddo

!-----------------------      Zero generation density is ready
!-------------------------------      Iterations:
do J=1, NA1
   NN(J) = NN0(J)
enddo
do JN=1, NNCX
   do J=1, NA1
      NIN(J) = NN0(J)
   enddo
   do J=1, NA1
      YZZ = 0.5*(NIN(1)*YKERNF(NRD*J + 1) + NIN(NA1)*YKERNF(NRD*J + NA1))
      do JJ=2, NA
         YZZ = YZZ + NIN(JJ)*YKERNF(NRD*J + JJ)
      enddo
      NN0(J) = YHA*YZZ
      NN(J)  = NN(J) + NN0(J)
   enddo
enddo

!------------------- End of iterations
do J=1, NA1
   YVI(J) = TI(J)*NN(J)
enddo
do J=1, NA1
   YZZ  = 0.5*(YVI(1)*YKERNF(NRD*J + 1) + YVI(NA1)*YKERNF(NRD*J +NA1))
   YZZN = 0.5*(NN(1) *YKERNF(NRD*J + 1) + NN(NA1) *YKERNF(NRD*J +NA1))
   do JJ=2, NA
      YZZN = YZZN + NN(JJ)*YKERNF(NRD*J + JJ)
      YZZ  = YZZ + YVI(JJ)*YKERNF(NRD*J + JJ)
   enddo
   TN(J) = (YZZ*YHA + TEN(J))/NN(J)
enddo
! albpl: [d/l] Plasma albedo
!            Pereverzev      15-05-95
! (Neutral_outflux)/(Neutral_influx) =
!      = sqrt{(NN-N1-N2)(NN*TN-N1*E1-N2*E2)}/(N1*sqrt(E1)+N2*sqrt(E2))
ALBPL = sqrt((NNCL + NNWM)*(NN(NA1) - 1) * &
            ((NNCL + NNWM)*NN(NA1)*TN(NA1) - NNCL*ENCL - NNWM*ENWM))/ &
             (NNCL*sqrt(ENCL) + NNWM*sqrt(ENWM))
if (NA1 .ge. NAB) return
do J=NA1+1, NAB
   NN(J) = NN(NA1)
   TN(J) = TN(NA1)
enddo
end subroutine NEUT

!======================================================================|
double precision function AKDLT(X)

implicit none

double precision, intent(in) :: X

if (X .GT. 30.) then
   AKDLT = 0.
else
   AKDLT = EXP(-X)
endif
end function AKDLT

