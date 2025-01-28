      subroutine DPCGRC(IDO, N, X, P, R, Z, RELERR, ITMAX)

      implicit none

      integer, intent(in) :: N, ITMAX
      real*8, intent(in) :: RELERR
      real*8, intent(in) , dimension(N) :: Z
      real*8, intent(out), dimension(N) :: P
      integer, intent(inout) :: IDO
      real*8, intent(inout), dimension(N) :: X, R

      integer :: i
      COMMON/DPCGR1/ITER
      integer ITER
      COMMON/DPCGR2/RSNORM
      real*8 RSNORM
      COMMON/DPCGR3/C
      real*8 C
      real*8 ALFA, BETA, ERROR, C1

      if(IDO == 0) then
         do i=1, N
            P(i) = X(i)
         enddo
         RSNORM = 0.
         do i=1, N
            RSNORM = RSNORM + R(i)*R(i)
         enddo
         RSNORM = SQRT(RSNORM)
         IDO = 1
         ITER = 0
      else if(IDO == 1 .AND. ITER == 0) then
         do i=1, N
            R(i) = R(i) - Z(i)
         enddo
         ERROR = 0.
         do i = 1, N
            ERROR = ERROR + R(i)*R(i)
         enddo
         ERROR = SQRT(ERROR)/RSNORM
         if(ERROR <= RELERR) then
            IDO = 3
            return
         endif
         IDO = 2
         ITER = 0
      else if(IDO == 2 .AND.ITER == 0) then
         do i=1, N
            P(i) = Z(i)
         enddo
         C = 0.
         do i=1, N
            C = C + Z(i)*R(i)
         enddo
         IDO = 1
         ITER = 1
      else if(IDO == 1) then
         ALFA = 0.
         do i=1, N
            ALFA = ALFA + Z(i)*P(i)
         enddo
         ALFA = C/ALFA
         do i=1, N
            X(i) = X(i) + ALFA*P(i)
         enddo
         do i=1, N
            R(i) = R(i) - ALFA*Z(i)
         enddo
         ERROR = 0.
         do i = 1, N
            ERROR = ERROR + R(i)*R(i)
         enddo
         ERROR = SQRT(ERROR)/RSNORM
         if(ERROR <= RELERR .OR. ITER >= ITMAX) then
            IDO = 3
            return
         endif
         IDO = 2
         ITER = ITER + 1
      else if(IDO == 2) then
         C1 = 0.
         do i=1, N
            C1 = C1 + Z(i)*R(i)
         enddo
         BETA = C1/C
         C = C1
         do i=1, N
            P(i) = Z(i) + BETA*P(i)
         enddo
         IDO = 1
      endif

      return
      end subroutine DPCGRC
