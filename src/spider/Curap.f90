subroutine TREFW( NOUT, NTER, NINFW, NP, KD, RP1, ZP1, RP2, ZP2, &
                  RESSEG, CURSEG, ND, NTYD, RD, ZD, VV, HH, &
                  RESD, CURD, R1, Z1, R2, Z2 )

! VV(I)  - VERTICAL   (Z) SIZE  OF THIN PLATE
! HH(I)  - HORISONTAL (R) SIZE  OF THIN PLATE

use iopath, only: path

implicit none

integer, intent(in) :: NOUT, NTER, NINFW
integer, intent(out) :: NP
integer, intent(out), dimension(*) :: KD, NTYD
real*8, intent(out), dimension(*) :: RP1, ZP1, RP2, ZP2, RESSEG, &
   CURSEG, RD, ZD, VV, HH, RESD, CURD, R1, Z1, R2, Z2

integer :: I, K, L, ND
real*8 :: DELR, DELZ, SSSFW, DDI, fr_cam, RCEFW
character(len=80) :: fname

write(fname,'(a,a)') TRIM(path), '/blanfw.dat'
open(NINFW,file=fname,form='formatted')

read(NINFW,*)  RCEFW
read(NINFW,*)  NP

ND = 0

if( NP.NE.0 )  then
   do L=1,NP
      read(NINFW,*) KD(L), RP1(L), ZP1(L), RP2(L), ZP2(L), &
                    RESSEG(L), CURSEG(L)
   enddo  
   do L=1,NP
      I = ND + 1
      DELR  = (RP2(L) - RP1(L)) / KD(L)
      DELZ  = (ZP2(L) - ZP1(L)) / KD(L)
      RD(I) = RP1(L) + 0.5D0*DELR
      ZD(I) = ZP1(L) + 0.5D0*DELZ
      VV(I) = DELZ
      HH(I) = DELR
      DDI   = dsqrt( HH(I)**2 + VV(I)**2 )
      R1(I) = RP1(L)
      Z1(I) = ZP1(L)
      R2(I) = RP1(L) + DELR
      Z2(I) = ZP1(L) + DELZ
      SSSFW = DDI/RD(I)
      if( KD(L) .GT. 1 ) then
         do K=1,KD(L)-1
            RD(I+K) = RD(I+K-1) + DELR
            ZD(I+K) = ZD(I+K-1) + DELZ
            VV(I+K) = VV(I)
            HH(I+K) = HH(I)
            R1(I+K) = R1(I+K-1) + DELR
            Z1(I+K) = Z1(I+K-1) + DELZ
            R2(I+K) = R2(I+K-1) + DELR
            Z2(I+K) = Z2(I+K-1) + DELZ
            SSSFW   = SSSFW + DDI/RD(I+K)
         enddo
      endif
      SSSFW = SSSFW * RESSEG(L)
      do K=1,KD(L)
         RESD(I+K-1) = SSSFW * RD(I+K-1)/DDI
         NTYD(I+K-1) = 1
         CURD(I+K-1) = CURSEG(L) / KD(L)
      enddo
      ND = ND + KD(L)
   enddo

   fr_cam=0.d0
   do i=1,np
      fr_cam=fr_cam+1.d0/RESSEG(i)
   enddo
   fr_cam=1.d0/fr_cam

   if( RCEFW.GT.0 )  then
      SSSFW = 0.0D0
      do L=1,ND
         SSSFW = SSSFW + dsqrt( HH(L)**2 + VV(L)**2 )/RD(L)
      enddo
      SSSFW = SSSFW * RCEFW
      do L=1,ND
         RESD(L) = SSSFW * RD(L) / dsqrt( HH(L)**2 + VV(L)**2 )
      enddo
   endif

endif

close(NINFW)

return
end subroutine TREFW

!------------------------------------------------------------------
subroutine TREBP( NOUT, NTER, NINFW, NP, KD, RP1, ZP1, RP2, ZP2, &
                  RESSEG, CURSEG, ND, NTYD, RD, ZD, VV, HH, &
                  RESD, CURD, R1, Z1, R2, Z2 )

!  VV(I)  - VERTICAL   (Z) SIZE OF THIN PLATE
!  HH(I)  - HORISONTAL (R) SIZE OF THIN PLATE

use iopath, only: path

implicit none

integer, intent(in) :: NOUT, NTER, NINFW
integer, intent(out) :: NP
integer, intent(out), dimension(*) :: KD, NTYD
real*8, intent(out), dimension(*) :: RP1, ZP1, RP2, ZP2, RESSEG, &
    CURSEG, RD, ZD, VV, HH, RESD, CURD, R1, Z1, R2, Z2

integer :: I, K, L, ND
real*8 :: DELR, DELZ, SSSFW, DDI, fr_cam, RCEBP
character(len=80) :: fname

write(fname,'(a,a)') TRIM(path), '/blanbp.dat'
open(NINFW,file=fname,form='formatted')
read(NINFW,*) RCEBP
read(NINFW,*) NP
ND = 0

if( NP.NE.0 )  then

   do L=1,NP
      read(NINFW,*) KD(L), RP1(L), ZP1(L), RP2(L), ZP2(L), &
                    RESSEG(L), CURSEG(L)
   enddo

   do L=1,NP
      I = ND + 1
      DELR  = (RP2(L) - RP1(L)) / KD(L)
      DELZ  = (ZP2(L) - ZP1(L)) / KD(L)
      RD(I) = RP1(L) + 0.5D0*DELR
      ZD(I) = ZP1(L) + 0.5D0*DELZ
      VV(I) = DELZ
      HH(I) = DELR
      DDI   = dsqrt( DELR**2 + DELZ**2 )
      R1(I) = RP1(L)
      Z1(I) = ZP1(L)
      R2(I) = RP1(L) + DELR
      Z2(I) = ZP1(L) + DELZ
      SSSFW = DDI/RD(I)
      if( KD(L) .GT. 1 ) then
         do K=1,KD(L)-1
            RD(I+K) = RD(I+K-1) + DELR
            ZD(I+K) = ZD(I+K-1) + DELZ
            VV(I+K) = DELZ
            HH(I+K) = DELR
            R1(I+K) = R1(I+K-1) + DELR
            Z1(I+K) = Z1(I+K-1) + DELZ
            R2(I+K) = R2(I+K-1) + DELR
            Z2(I+K) = Z2(I+K-1) + DELZ
            SSSFW   = SSSFW + DDI/RD(I+K)
         enddo
      endif
      SSSFW = SSSFW * RESSEG(L)
      do K=1,KD(L)
         RESD(I+K-1) = SSSFW * RD(I+K-1)/DDI
         NTYD(I+K-1) = 1
         CURD(I+K-1) = CURSEG(L) / KD(L)
      enddo
      ND = ND + KD(L)
   enddo

   if( RCEBP.GT.0 )  then
      SSSFW = 0.0D0
      do L=1,ND
         SSSFW = SSSFW + dsqrt( HH(L)**2 + VV(L)**2 )/RD(L)
      enddo
      SSSFW = SSSFW * RCEBP
      do L=1,ND
         RESD(L) = SSSFW * RD(L) / dsqrt( HH(L)**2 + VV(L)**2 )
      enddo
   endif

endif

close(NINFW)

return
end subroutine TREBP

!-----------------------------------------------------------------
subroutine TREVV( NOUT, NTER, NINFW, NP, KD, RP1, ZP1, RP2, ZP2, &
                  RESSEG, CURSEG, ND, NTYD, RD, ZD, VV, HH, &
                  RESD, CURD, R1, Z1, R2, Z2 )

! VV(I)  - VERTICAL   (Z) SIZE OF THIN PLATE
! HH(I)  - HORISONTAL (R) SIZE OF THIN PLATE

use iopath, only: path

implicit none

integer, intent(in) :: NOUT, NTER, NINFW
integer, intent(out) :: NP
integer, intent(out), dimension(*) :: KD, NTYD
real*8, intent(out), dimension(*) :: RP1, ZP1, RP2, ZP2, RESSEG, &
   CURSEG, RD, ZD, VV, HH, RESD, CURD, R1, Z1, R2, Z2

integer :: I, K, L, ND
real*8 :: DELR, DELZ, SSSFW, DDI, fr_cam, RCEVV
character(len=80) :: fname

write(fname,'(a,a)') TRIM(path), '/vacves.dat'
open(NINFW,file=fname,form='formatted')
read(NINFW,*)  RCEVV
read(NINFW,*)  NP
ND = 0

if( NP.NE.0 )  then

   do L=1,NP
      read(NINFW,*) KD(L), RP1(L), ZP1(L), RP2(L), ZP2(L), &
                    RESSEG(L), CURSEG(L)
   enddo

   do L=1,NP
      I = ND + 1
      DELR  = (RP2(L) - RP1(L)) / KD(L)
      DELZ  = (ZP2(L) - ZP1(L)) / KD(L)
      RD(I) = RP1(L) + 0.5D0*DELR
      ZD(I) = ZP1(L) + 0.5D0*DELZ
      VV(I) = DELZ
      HH(I) = DELR
      DDI   = dsqrt( DELR**2 + DELZ**2 )
      R1(I) = RP1(L)
      Z1(I) = ZP1(L)
      R2(I) = RP1(L) + DELR
      Z2(I) = ZP1(L) + DELZ
      SSSFW = DDI/RD(I)
      if( KD(L) .GT. 1 ) then
         do K=1,KD(L)-1
            RD(I+K) = RD(I+K-1) + DELR
            ZD(I+K) = ZD(I+K-1) + DELZ
            VV(I+K) = DELZ
            HH(I+K) = DELR
            R1(I+K) = R1(I+K-1) + DELR
            Z1(I+K) = Z1(I+K-1) + DELZ
            R2(I+K) = R2(I+K-1) + DELR
            Z2(I+K) = Z2(I+K-1) + DELZ
            SSSFW   = SSSFW + DDI/RD(I+K)
         enddo
      endif
      SSSFW = SSSFW * RESSEG(L)
      do K=1,KD(L)
         RESD(I+K-1) = SSSFW * RD(I+K-1)/DDI
         NTYD(I+K-1) = 1
         CURD(I+K-1) = CURSEG(L) / KD(L)
      enddo
      ND = ND + KD(L)
   enddo

   if( RCEVV.GT.0 )  then
      SSSFW = 0.0D0
      do L=1,ND
         SSSFW = SSSFW + dsqrt( HH(L)**2 + VV(L)**2 )/RD(L)
      enddo
      SSSFW = SSSFW * RCEVV
      do L=1,ND
         RESD(L) = SSSFW * RD(L) / dsqrt( HH(L)**2 + VV(L)**2 )
      enddo
   endif

endif

close(NINFW)

return
end subroutine TREVV
